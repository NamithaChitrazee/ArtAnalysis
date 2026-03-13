#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include "TTree.h"

using namespace mu2e;

class ComboHitBkgDiag : public art::EDAnalyzer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> ComboHitCollection { fhicl::Name("ComboHitCollection") };
    fhicl::Atom<art::InputTag> StrawDigiMCCollection { fhicl::Name("StrawDigiMCCollection") };
  };

  explicit ComboHitBkgDiag(const art::EDAnalyzer::Table<Config>& config);

  void beginJob() override;
  void analyze(const art::Event& e) override;
  void endJob() override;

private:
  art::ProductToken<ComboHitCollection> _chToken;
  art::ProductToken<StrawDigiMCCollection> _mcdigiToken;

  const ComboHitCollection* _chcol = nullptr;
  const StrawDigiMCCollection* _mcdigis = nullptr;

  TTree* _tree = nullptr;

  int _iev;
  int _creationCode;
  int _pdg;
  bool _flaggedBkg;
  int _totalce, _totalbkg, _totalsig;
  int _cemis, _bkgtrue, _sigmis;
 
};

ComboHitBkgDiag::ComboHitBkgDiag(const art::EDAnalyzer::Table<Config>& config) :
  art::EDAnalyzer{config},
  _chToken(consumes<ComboHitCollection>(config().ComboHitCollection())),
  _mcdigiToken(consumes<StrawDigiMCCollection>(config().StrawDigiMCCollection()))
{}

void ComboHitBkgDiag::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  _totalce = 0;
  _totalbkg = 0;
  _totalsig = 0;
  _bkgtrue = 0;
  _cemis = 0;
  _sigmis = 0;
  _tree = tfs->make<TTree>("chbkg", "ComboHit background truth check");

  _tree->Branch("iev", &_iev, "iev/I");
  _tree->Branch("creationCode", &_creationCode, "creationCode/I");
  _tree->Branch("pdg", &_pdg, "pdg/I");
  _tree->Branch("flaggedBkg", &_flaggedBkg, "flaggedBkg/O");
}

void ComboHitBkgDiag::analyze(const art::Event& event) {
  _iev = event.id().event();

  _chcol = event.getValidHandle(_chToken).product();
  _mcdigis = event.getValidHandle(_mcdigiToken).product();
  for (size_t ich = 0; ich < _chcol->size(); ++ich) {
    ComboHit const& ch = _chcol->at(ich);

    std::vector<StrawDigiIndex> dids;
    _chcol->fillStrawDigiIndices(ich, dids);

    StrawDigiMC const& mcdigi = _mcdigis->at(dids[0]);
    auto const& sgsp = mcdigi.earlyStrawGasStep();
    art::Ptr<SimParticle> const& sp = sgsp->simParticle();

    _pdg = sp->pdgId();
    _creationCode = sp->creationCode();
    if(_creationCode == 12 or _creationCode == 17 or _creationCode == 40 or _creationCode == 165 or _creationCode == 13
       or _creationCode == 74 or _creationCode == 58 or _creationCode == 114 or _creationCode == 166
       or _creationCode == 21 or _creationCode == 115 or _creationCode == 97 or _creationCode == 101 or _creationCode == 78
       or _creationCode == 113 or _creationCode == 16 or _creationCode == 23 or _creationCode == 35 or _creationCode == 72
       or _creationCode == 56  or _creationCode == 2 or _creationCode == 31 or _creationCode == 59  or _creationCode == 109)
        _totalbkg += ch.nStrawHits();
    else if(_creationCode == 167)
        _totalce += ch.nStrawHits(); 
    else
        _totalsig += ch.nStrawHits();
    _flaggedBkg = ch.flag().hasAllProperties(StrawHitFlag::bkg);
    if(_flaggedBkg){
      if(_creationCode == 12 or _creationCode == 17 or _creationCode == 40 or _creationCode == 165 or _creationCode == 13
       or _creationCode == 74 or _creationCode == 58 or _creationCode == 114 or _creationCode == 166
       or _creationCode == 21 or _creationCode == 115 or _creationCode == 97 or _creationCode == 101 or _creationCode == 78
       or _creationCode == 113 or _creationCode == 16 or _creationCode == 23 or _creationCode == 35 or _creationCode == 72
         or _creationCode == 56  or _creationCode == 2 or _creationCode == 31 or _creationCode == 59  or _creationCode == 109){
        _bkgtrue += ch.nStrawHits();
      }
      else if(_creationCode == 167)
        _cemis += ch.nStrawHits(); 
      else
        _sigmis += ch.nStrawHits();
    }
    _tree->Fill();
  }
}

void ComboHitBkgDiag::endJob() {
  std::cout << "====================================\n";
  std::cout << " FINAL COUNTS OVER ALL EVENTS\n";
  std::cout << " Total background hits = " << _totalbkg <<" True background hits = " << _bkgtrue << "\n";
  std::cout << " Total CE hits = " << _totalce <<" CE mis-tagged hits   = " << _cemis   << "\n";
  std::cout << " Total Signal hits = "<< _totalsig <<" Signal mis-tagged    = " << _sigmis  << "\n";
  std::cout << "====================================\n";
}

DEFINE_ART_MODULE(ComboHitBkgDiag)
