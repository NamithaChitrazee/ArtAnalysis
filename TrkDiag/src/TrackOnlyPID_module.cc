//
// Reconstruction-level PID determination from a track without the use of calorimeter information.
//
// Produces a collection of PID scores per input track collection, aligned by index
//
// Original author: Michael MacKenzie
//

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// utilities
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedDtDt.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"

// models
#include "ArtAnalysis/TrkDiag/inc/TrkOnlyPID_v0.hxx"

// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;

namespace TMVA_SOFIE_TrkOnlyPID_v0 {
  class Session;
}

namespace mu2e {
  class TrackOnlyPID : public art::EDProducer {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Sequence<std::string>         kalSeeds    { Name("KalSeeds")    , Comment("KalSeed (Ptr) collection names") };
        fhicl::OptionalSequence<std::string> kalSeedDtDts{ Name("KalSeedDtDts"), Comment("KalSeedDtDt collection names") };
        fhicl::Atom<std::string>             datFilename { Name("DatFilename") , Comment("Filename for the .dat file to use")};
        fhicl::Atom<int>                     MVAVersion  { Name("MVAVersion")  , Comment("MVA version ID to determine the input features")};
        fhicl::Atom<bool>                    printMVA    { Name("PrintMVA")    , Comment("Print the MVA used"), false};
        fhicl::Atom<int>                     debug       { Name("DebugLevel")  , Comment("Debug printout Level"), 0};
      };

      using Parameters = art::EDProducer::Table<Config>;
      TrackOnlyPID(const Parameters& conf);

    private:
      void produce(art::Event& event) override;
      void initializeMVA(std::string xmlfilename);
      float evaluateMVA(const KalSeed& seed, const KalSeedDtDt* dtdt);

      // MVA version-specific calls
      float evaluateMVA_v0(const KalSeed& seed, const KalSeedDtDt* dtdt);

      std::vector<std::string> kalSeeds_;
      std::vector<std::string> kalSeedDtDts_;
      int  MVAVersion_;
      bool printMVA_;
      int debugLevel_;

      std::shared_ptr<TMVA_SOFIE_TrkOnlyPID_v0::Session> mva_v0_;
  };

  //================================================================
  TrackOnlyPID::TrackOnlyPID(const Parameters& conf)
    : art::EDProducer(conf)
    , kalSeeds_(conf().kalSeeds())
    , MVAVersion_(conf().MVAVersion())
    , printMVA_(conf().printMVA())
    , debugLevel_(conf().debug())
  {

    if(MVAVersion_ > 0) throw cet::exception("RECO") << "Unknown MVA feature version " << MVAVersion_;

    if(conf().kalSeedDtDts(kalSeedDtDts_)) {
      if(kalSeeds_.size() != kalSeedDtDts_.size()) throw cet::exception("RECO") << "KalSeed and KalSeedDtDt lists must match";
    } else if(MVAVersion_ == 0) {
      throw cet::exception("RECO") << "KalSeedDtDt is not provided but is required for MVA feature version 0";
    }

    // Produce a PID collection per track collection
    for(const auto& name : kalSeeds_) {
      produces<MVAResultCollection>(name);
    }

    // Define the MVA
    ConfigFileLookupPolicy configFile;
    if     (MVAVersion_ == 0) mva_v0_ = std::make_shared<TMVA_SOFIE_TrkOnlyPID_v0::Session>(configFile(conf().datFilename()));
  }

  //================================================================
  float TrackOnlyPID::evaluateMVA_v0(const KalSeed& seed, const KalSeedDtDt* dtdt) {
    double score = -999.; // < -1 = invalid

    if(!dtdt) return score;

    // Define the features
    const float fitcon  = seed.fitConsistency();
    const float slope   = dtdt->slope();
    const int   nHits   = seed.nHits(false);
    const int   nActive = seed.nHits(true);
    int nNull = 0;
    for(auto hit = seed.hits().begin(); hit != seed.hits().end(); ++hit) {
      if(hit->strawHitState() > WireHitState::inactive && hit->strawHitState()==WireHitState::null) ++nNull;
    }
    const float activeFrac = nActive * 1.f / nHits;
    const float nullFrac   = nNull * 1.f / nActive;


    // Assign the features
    std::array<float,4> features = {activeFrac, nullFrac, fitcon, slope};

    // Evaluate the MVA
    const auto mvaout = mva_v0_->infer(features.data());
    score = mvaout[0];
    if(debugLevel_ > 0) {
      printf("[TrackOnlyPID::%s] Input features: {%.3f, %.3f, %.3f, %.3f} output: %.5f\n", __func__,
             features[0], features[1], features[2], features[3], score);
    }

    // Return the score
    return score;
  }

  //================================================================
  float TrackOnlyPID::evaluateMVA(const KalSeed& seed, const KalSeedDtDt* dtdt) {
    float score = -999.; // < -1 = invalid

    // Check if it's an acceptable track fit
    static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
    if(!seed.status().hasAllProperties(goodfit)) return score;

    // Evaluate the MVA using the proper call for the given version
    switch(MVAVersion_) {
    case  0: score = evaluateMVA_v0(seed, dtdt); break;
    default: throw cet::exception("RECO") << "Unknown MVA feature version " << MVAVersion_;
    }

    // Return the score
    return score;
  }

  void TrackOnlyPID::produce(art::Event& event ) {
      // Loop over all KalSeed collections
    for(size_t index = 0; index < kalSeeds_.size(); ++index) {
      const auto& name = kalSeeds_.at(index);

      // Retrieve the KalSeed collection from the event, checking if it is a collection of KalSeed or KalSeedPtr
      art::Handle<KalSeedCollection> seedHandle;
      event.getByLabel(name, seedHandle);
      art::Handle<KalSeedPtrCollection> seedPtrHandle;
      const bool isSeedCollection = seedHandle.isValid();
      if(!isSeedCollection) {
        event.getByLabel(name, seedPtrHandle);
        if(!seedPtrHandle.isValid()) {
          throw cet::exception("RECO") << "TrackOnlyPID: No KalSeed or KalSeedPtr collection with label " << name << std::endl;
        }
      }
      const auto nseeds = (isSeedCollection) ? seedHandle->size() : seedPtrHandle->size();

      // Retrieve the KalSeedDtDt collection, if requested
      const KalSeedDtDtCollection* dtdtCollection = nullptr;
      if(!kalSeedDtDts_.empty()) {
        auto dtdt_handle = event.getValidHandle<KalSeedDtDtCollection>(kalSeedDtDts_.at(index));
        dtdtCollection = &(*dtdt_handle);
        if(dtdtCollection->size() != nseeds) throw cet::exception("RECO") << "KalSeed and KalSeedDtDt collections must have matching size. N(seeds) = "
                                                                          << nseeds << " N(dtdts) = " << dtdtCollection->size()
                                                                          << " Seeds tag = " << name << " DtDts tag = " << kalSeedDtDts_.at(index);
      }


      // Create the output collection
      std::unique_ptr<MVAResultCollection> mvaCol(new MVAResultCollection());

      // Loop over all seeds in the collection and evaluate the PID score
      if(debugLevel_ > 0) std::cout << "[TrackOnlyPID::" << __func__ << "] Processing " << nseeds << " seeds from collection " << name << std::endl;
      for(size_t iseed = 0; iseed < nseeds; ++iseed) {
        const auto& seed = (isSeedCollection) ? seedHandle->at(iseed) : *seedPtrHandle->at(iseed);
        if(debugLevel_ > 1) std::cout << "[TrackOnlyPID::" << __func__ << "] Processing seed " << iseed << " with " << seed.hits().size() << " hits" << std::endl;
        mvaCol->emplace_back(MVAResult(evaluateMVA(seed, (dtdtCollection) ? &dtdtCollection->at(iseed) : nullptr)));
      }

      // Put the results into the event
      event.put(std::move(mvaCol), name);
    }
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackOnlyPID)
