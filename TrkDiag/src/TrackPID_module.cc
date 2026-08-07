//
// Reconstruction-level PID determination from a track.  The current implementation
// uses the calorimeter matching information (TrkCaloHit) in addition to basic track-level
// information in an MVA.
//
// Produces a collection of PID scores per input track collection, aligned by index
//
// Original author: Dave Brown (LBNL)
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
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

// models
#include "ArtAnalysis/TrkDiag/inc/TrackPID_v0.hxx"
#include "ArtAnalysis/TrkDiag/inc/TrackPID_v1.hxx"

// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;

namespace TMVA_SOFIE_TrackPID_v0 {
  class Session;
}
namespace TMVA_SOFIE_TrackPID_v1 {
  class Session;
}

namespace mu2e {
  class TrackPID : public art::EDProducer {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Atom<float>                   maxDE       { Name("MaxDE")       , Comment("Maximum E(calo) - P(trk)")};
        fhicl::Sequence<std::string>         kalSeeds    { Name("KalSeeds")    , Comment("KalSeed (Ptr) collection names") };
        fhicl::OptionalSequence<std::string> kalSeedDtDts{ Name("KalSeedDtDts"), Comment("KalSeedDtDt collection names") };
        fhicl::Atom<std::string>             datFilename { Name("DatFilename") , Comment("Filename for the .dat file to use")};
        fhicl::Atom<int>                     MVAVersion  { Name("MVAVersion")  , Comment("MVA version ID to determine the input features")};
        fhicl::Atom<bool>                    printMVA    { Name("PrintMVA")    , Comment("Print the MVA used"), false};
        fhicl::Atom<int>                     debug       { Name("DebugLevel")  , Comment("Debug printout Level"), 0};
      };

      using Parameters = art::EDProducer::Table<Config>;
      TrackPID(const Parameters& conf);

    private:
      void produce(art::Event& event) override;
      void initializeMVA(std::string xmlfilename);
      float evaluateMVA(const KalSeed& seed, const KalSeedDtDt* dtdt);
      XYZVectorD momentumAtCalo(const KalSeed& seed);

      // MVA version-specific calls
      float evaluateMVA_v0(const KalSeed& seed);
      float evaluateMVA_v1(const KalSeed& seed, const KalSeedDtDt* dtdt);

      float maxDe_;
      std::vector<std::string> kalSeeds_;
      std::vector<std::string> kalSeedDtDts_;
      int  MVAVersion_;
      bool printMVA_;
      int debugLevel_;

      const mu2e::Calorimeter* calo_;

      std::shared_ptr<TMVA_SOFIE_TrackPID_v0::Session> mva_v0_;
      std::shared_ptr<TMVA_SOFIE_TrackPID_v1::Session> mva_v1_;
  };

  //================================================================
  TrackPID::TrackPID(const Parameters& conf)
    : art::EDProducer(conf)
    , maxDe_(conf().maxDE())
    , kalSeeds_(conf().kalSeeds())
    , MVAVersion_(conf().MVAVersion())
    , printMVA_(conf().printMVA())
    , debugLevel_(conf().debug())
  {

    if(MVAVersion_ > 1) throw cet::exception("RECO") << "Unknown MVA feature version " << MVAVersion_;

    if(conf().kalSeedDtDts(kalSeedDtDts_)) {
      if(kalSeeds_.size() != kalSeedDtDts_.size()) throw cet::exception("RECO") << "KalSeed and KalSeedDtDt lists must match";
    } else if(MVAVersion_ == 1) {
      throw cet::exception("RECO") << "KalSeedDtDt is not provided but is required for MVA feature version 1";
    }

    // Produce a PID collection per track collection
    for(const auto& name : kalSeeds_) {
      produces<MVAResultCollection>(name);
    }

    // Define the MVA
    ConfigFileLookupPolicy configFile;
    if     (MVAVersion_ == 0) mva_v0_ = std::make_shared<TMVA_SOFIE_TrackPID_v0::Session>(configFile(conf().datFilename()));
    else if(MVAVersion_ == 1) mva_v1_ = std::make_shared<TMVA_SOFIE_TrackPID_v1::Session>(configFile(conf().datFilename()));
  }

  //================================================================
  XYZVectorD TrackPID::momentumAtCalo(const KalSeed& seed) {
    if(!seed.hasCaloCluster()) return XYZVectorD(0.,0.,0.); // no cluster

    // Find the segment that is closest to the cluster position
    auto const& tchs = seed.caloHit();
    return seed.nearestSegment(tchs._rptoca)->momentum3();
  }

  //================================================================
  float TrackPID::evaluateMVA_v0(const KalSeed& seed) {
    std::array<float,4> features{-9999,-9999,-9999,-9999}; // features used for training
    double score = -999.; // < -1 = invalid
    constexpr float dtoffset = -1.15; // value used for this version of the MVA

    // Require a calo cluster
    if(!seed.hasCaloCluster() || !seed.caloHit()._flag.hasAllProperties(StrawHitFlag::active)) {
      if(debugLevel_ > 1) printf("[TrackPID::%s] Input track has no calo cluster\n", __func__);
      return score;
    }

    // Fill the features
    auto const& tchs = seed.caloHit();
    auto const& cc = tchs.caloCluster();
    auto trkmom = momentumAtCalo(seed);
    features[0] = cc->energyDep() - sqrt(trkmom.Mag2());
    // move into detector coordinates.  Yikes!!
    XYZVectorF cpos = XYZVectorF(calo_->mu2eToTracker(calo_->diskFFToMu2e( cc->diskID(), cc->cog3Vector())));
    features[1] = sqrt(cpos.Perp2());
    // compute transverse direction WRT position
    cpos.SetZ(0.0);
    trkmom.SetZ(0.0);
    features[2] = cpos.Dot(trkmom)/sqrt(cpos.Mag2()*trkmom.Mag2());
    // the following includes the (Calibrated) light-propagation time delay.  It should eventually be put in the reconstruction FIXME!
    // This velocity should come from conditions FIXME!
    features[3] = tchs.t0().t0()-tchs.time()- std::min((float)200.0,std::max((float)0.0,tchs.hitLen()))*0.005 - dtoffset;
    // hard cut on the energy difference.  This rejects cosmic rays which hit the calo and produce an upstream-going track that is then
    // reconstructed as a downstream particle associated to this cluster
    if(features[0] < maxDe_) {
      // evaluate the MVA
      auto mvaout = mva_v0_->infer(features.data());
      score = mvaout[0];
    } else if(debugLevel_ > 0) {
      printf("[TrackPID::%s] Energy difference at %.3f, above threshold %.3f", __func__, features[0], maxDe_);
    }

    if(debugLevel_ > 0) {
      printf("[TrackPID::%s] Input features: {%.3f, %.3f, %.3f, %.3f} output: %.5f\n", __func__,
             features[0], features[1], features[2], features[3], score);
    }

    return score;
  }

  //================================================================
  float TrackPID::evaluateMVA_v1(const KalSeed& seed, const KalSeedDtDt* dtdt) {
    double score = -999.; // < -1 = invalid

    if(!dtdt) return score;

    // Require a calo cluster
    if(!seed.hasCaloCluster() || !seed.caloHit()._flag.hasAllProperties(StrawHitFlag::active)) {
      if(debugLevel_ > 1) printf("[TrackPID::%s] Input track has no calo cluster\n", __func__);
      return score;
    }

    // Retrieve the calo cluster
    auto const& tchs = seed.caloHit();
    auto const& cc   = tchs.caloCluster();
    const float edep = cc->energyDep();

    // Retrieve the track momentum at the calo
    const auto  momVec = momentumAtCalo(seed);
    const float mom    = std::sqrt(momVec.Mag2());

    // Check if this track should be skipped
    if(edep - mom > maxDe_) return score; // Selection provided by the user
    if(mom == 0.)           return score; // E / P not defined

    // Define the features
    const float ep     = edep / mom;
    const float dt     = tchs._udt; // unbiased time difference
    const float fitcon = seed.fitConsistency();
    const float slope  = dtdt->slope();

    // Assign the features
    std::array<float,4> features = {ep, dt, fitcon, slope};

    // Evaluate the MVA
    const auto mvaout = mva_v1_->infer(features.data());
    score = mvaout[0];
    if(debugLevel_ > 0) {
      printf("[TrackPID::%s] Input features: {%.3f, %.3f, %.3f, %.3f} output: %.5f\n", __func__,
             features[0], features[1], features[2], features[3], score);
    }

    // Return the score
    return score;
  }

  //================================================================
  float TrackPID::evaluateMVA(const KalSeed& seed, const KalSeedDtDt* dtdt) {
    float score = -999.; // < -1 = invalid

    // Check if it's an acceptable track fit
    static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
    if(!seed.status().hasAllProperties(goodfit)) return score;

    // Evaluate the MVA using the proper call for the given version
    switch(MVAVersion_) {
    case  0: score = evaluateMVA_v0(seed); break;
    case  1: score = evaluateMVA_v1(seed, dtdt); break;
    default: throw cet::exception("RECO") << "Unknown MVA feature version " << MVAVersion_;
    }

    // Return the score
    return score;
  }

  void TrackPID::produce(art::Event& event ) {
    mu2e::GeomHandle<mu2e::Calorimeter> calo;
    calo_ = calo.get();

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
          throw cet::exception("RECO") << "TrackPID: No KalSeed or KalSeedPtr collection with label " << name << std::endl;
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
      if(debugLevel_ > 0) std::cout << "[TrackPID::" << __func__ << "] Processing " << nseeds << " seeds from collection " << name << std::endl;
      for(size_t iseed = 0; iseed < nseeds; ++iseed) {
        const auto& seed = (isSeedCollection) ? seedHandle->at(iseed) : *seedPtrHandle->at(iseed);
        if(debugLevel_ > 1) std::cout << "[TrackPID::" << __func__ << "] Processing seed " << iseed << " with " << seed.hits().size() << " hits" << std::endl;
        mvaCol->emplace_back(MVAResult(evaluateMVA(seed, (dtdtCollection) ? &dtdtCollection->at(iseed) : nullptr)));
      }

      // Put the results into the event
      event.put(std::move(mvaCol), name);
    }
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackPID)
