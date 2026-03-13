//
// Low-energy electron background diagnostics.  Split out of FlagBkgHits
//
// Original author D. Brown
//
//
/// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
// data
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "ArtAnalysis/TrkDiag/inc/BkgHitInfo.hh"
// art
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
using std::string;
using namespace ROOT::Math::VectorUtil;
namespace mu2e
{

  class BkgDiag : public art::EDAnalyzer {
    public:

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<int>              diag{                     Name("diagLevel"),                 Comment("Diag Level"),0 };
        fhicl::Atom<int>              debug{                    Name("debugLevel"),                Comment("Debug Level"),0 };
        fhicl::Atom<bool>             mcdiag{                   Name("MCDiag"),                    Comment("MonteCarlo Diag"), true };
        fhicl::Atom<bool>             hdiag{                    Name("HitDiag"),                   Comment("Hit-level Diag tuple"), false };
        fhicl::Atom<float>            maxdt{                    Name("maxTimeDifference"),         Comment("Max Time Difference"), 50.0 };
        fhicl::Atom<float>            maxdrho{                  Name("maxRhoDifference"),          Comment("Max Rho Difference"), 50.0 };
        fhicl::Atom<art::InputTag>    ComboHitCollection{       Name("ComboHitCollection"),        Comment("ComboHit collection name") };
        fhicl::Atom<art::InputTag>    BkgClusterCollection{     Name("BkgClusterCollection"),      Comment("BackgroundCluster collection name") };
        fhicl::Atom<art::InputTag>    BkgClusterHitCollection{  Name("BkgClusterHitCollection"),   Comment("BackgroundClusterHit collection name") };
        fhicl::Atom<art::InputTag>    StrawDigiMCCollection{    Name("StrawDigiMCCollection"),     Comment("StrawDigiMC collection name") };
        fhicl::Atom<art::InputTag>    MCPrimary{                Name("MCPrimary"),                 Comment("MC Primary Particle") };
        fhicl::Atom<std::string>      expectedCHLevel{          Name("ExpectedCHLevel"),           Comment("Level of the input ComboHitCollection") };
     };

      explicit BkgDiag(const art::EDAnalyzer::Table<Config>& config);
      virtual ~BkgDiag();
      virtual void beginJob();
      virtual void analyze(const art::Event& e);
    private:
      // helper functions
      void fillStrawHitInfo(size_t ich, StrawHitInfo& bkghinfo) const;
      void fillStrawHitInfoMC(StrawDigiMC const& mcdigi, art::Ptr<SimParticle>const& mptr, StrawHitInfo& shinfo) const;
      bool findData(const art::Event& e);
      void findMajority(std::vector<StrawDigiIndex>const& dids, art::Ptr<SimParticle>& mptr,XYZVectorF& mmom, int& mostCommonCode, float& fraction) const;
      void findMain(std::vector<StrawDigiIndex>const& dids, art::Ptr<SimParticle>& mptr,XYZVectorF& mmom, std::vector<int>& icontrib) const;

      // control flags
      int _diag,_debug;
      bool _mcdiag, _hdiag;
      float _maxdt, _maxdrho;
      // data tags
      art::ProductToken<ComboHitCollection> _chToken;
      art::ProductToken<BkgClusterCollection> _bkgcToken;
      art::ProductToken<BkgClusterHitCollection> _bkghToken;
      art::ProductToken<StrawDigiMCCollection> _mcdigisToken;
      art::ProductToken<PrimaryParticle> _mcprimaryToken;
      // time offset
      // cache of event objects
      const ComboHitCollection* _chcol;
      const StrawDigiMCCollection *_mcdigis;
      const PrimaryParticle *_mcprimary;
      const BkgClusterCollection *_bkgccol;
      const BkgClusterHitCollection *_bkghitcol;

      //Input Level Check
      StrawIdMask::Level expectedCHLevel_;
      std::string expectedCHLevelName_;

      // background diagnostics
      TTree* _bcdiag = 0;
      TTree* _bhdiag = 0;
      int _iev = 0;
      XYZVectorF _cpos;
      float _rmscposx = 0;
      float _rmscposy = 0;
      float _rmscrho = 0;
      float _rmscphi = 0;
      float _clusterdensity = 0;
      float _cchi2 = 0;
      float _ccons = 0;
      float _ctime = 0;
      float _cedep = 0;
      float _cpitch = 0;
      float _cyaw = 0;
      float _cqual = 0;
      float _csthqual = 0;
      float _ccomqual = 0;
      float _cecc = 0;
      float _rmsctime = 0;
      float _avecedep = 0;
      float _mindt = 0;
      float _mindrho = 0;
      bool _isinit = false;
      bool _isbkg = false;
      bool _isref = false;
      bool _isolated = false;
      bool _stereo = false;
      int _cluIdx, _nactive, _nch, _nsh, _nsth, _nsha, _nbkg;
      int _cndof;
      float _crho;
      float _zmin;
      float _zmax;
      float _zgap;
      float _deltaPhi;
    //float _slope;
      float _zdiff;
      float _phidiff;
      float _zgap2;
      float _pfrac;
      float _kQ;
      int _np, _fp, _lp, _pgap;

      // MC truth variables
      int _mpdg, _mproc, _ncontrib, _icontrib[512];
    //int _truebkghits, _truesighits, _totalbkghits, _totalsighits, _totalCEhits, _trueCEhits, _misCEhits;
      int _prel;
      int _code;
      float _frac;
      XYZVectorF _mmom, _mopos;
      int _nconv, _ndelta, _ncompt, _ngconv, _nebkg, _nprot, _nmain, _nsibling, _nrel;
      std::vector<BkgHitInfo> _bkghinfo;

      int   _nhits,_hitPdg[8192],_hitproc[8192],_hitnch[8192], _hitnsh[8192];
      std::vector<XYZVectorF> _hitPos;
      float _hitTime[8192];
  };

  BkgDiag::BkgDiag(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    _diag( config().diag() ),
    _debug( config().debug() ),
    _mcdiag( config().mcdiag() ),
    _hdiag( config().hdiag() ),
    _maxdt( config().maxdt() ),
    _maxdrho( config().maxdrho() ),
    _chToken{ consumes<ComboHitCollection>(config().ComboHitCollection() ) },
    _bkgcToken{ consumes<BkgClusterCollection>(config().BkgClusterCollection() ) },
    _bkghToken{ consumes<BkgClusterHitCollection>(config().BkgClusterHitCollection() ) },
    _mcdigisToken{ consumes<StrawDigiMCCollection>(config().StrawDigiMCCollection() ) },
    _mcprimaryToken{ consumes<PrimaryParticle>(config().MCPrimary() ) }
  {
    StrawIdMask mask(config().expectedCHLevel());
    expectedCHLevel_ = mask.level();
    expectedCHLevelName_ = mask.levelName();
  }

  BkgDiag::~BkgDiag(){}

  void BkgDiag::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    _iev=0;
    // detailed delta diagnostics
    _bcdiag=tfs->make<TTree>("bkgcdiag","background cluster diagnostics");
    // general branches
    _bcdiag->Branch("iev",&_iev,"iev/I");
    // cluster info branches
    _bcdiag->Branch("cpos",&_cpos);
    _bcdiag->Branch("rmscposx",&_rmscposx,"rmscposx/F");
    _bcdiag->Branch("rmscposy",&_rmscposy,"rmscposy/F");
    _bcdiag->Branch("rmscrho",&_rmscrho,"rmscrho/F");
    _bcdiag->Branch("rmscphi",&_rmscphi,"rmscphi/F");
    _bcdiag->Branch("clusterdensity",&_clusterdensity,"clusterdensity/F");
    _bcdiag->Branch("cpitch",&_cpitch,"cpitch/F");
    _bcdiag->Branch("cyaw",&_cyaw,"cyaw/F");
    _bcdiag->Branch("cqual",&_cqual,"cqual/F");
    _bcdiag->Branch("csthqual",&_csthqual,"csthqual/F");
    _bcdiag->Branch("ccomqual",&_ccomqual,"ccomqual/F");
    _bcdiag->Branch("cecc",&_cecc,"cecc/F");
    _bcdiag->Branch("cchi2",&_cchi2,"cchi2/F");
    _bcdiag->Branch("ccons",&_ccons,"ccons/F");
    _bcdiag->Branch("ctime",&_ctime,"ctime/F");
    _bcdiag->Branch("cedep",&_cedep,"cedep/F");
    _bcdiag->Branch("rmsctime",&_rmsctime,"rmsctime/F");
    _bcdiag->Branch("avecedep",&_avecedep,"avecedep/F");
    _bcdiag->Branch("isinit",&_isinit,"isinit/B");
    _bcdiag->Branch("isbkg",&_isbkg,"isbkg/B");
    _bcdiag->Branch("isref",&_isref,"isref/B");
    _bcdiag->Branch("isolated",&_isolated,"isolated/B");
    _bcdiag->Branch("stereo",&_stereo,"stereo/B");
    _bcdiag->Branch("mindt",&_mindt,"mindt/F");
    _bcdiag->Branch("mindrho",&_mindrho,"mindrho/F");
    _bcdiag->Branch("nch",&_nch,"nch/I");
    _bcdiag->Branch("nsh",&_nsh,"nsh/I");
    _bcdiag->Branch("nsth",&_nsth,"nsth/I");
    _bcdiag->Branch("nactive",&_nactive,"nactive/I");
    _bcdiag->Branch("nsha",&_nsha,"nsha/I");
    _bcdiag->Branch("nbkg",&_nbkg,"nbkg/I");
    _bcdiag->Branch("cndof",&_cndof,"cndof/I");
    _bcdiag->Branch("cluIdx",&_cluIdx,"cluIdx/I");
    // cluster hit info branch
    if(_diag > 0)
      _bcdiag->Branch("bkghinfo",&_bkghinfo);
    _bcdiag->Branch("crho",&_crho,"crho/F");
    _bcdiag->Branch("zmin",&_zmin,"zmin/F");
    _bcdiag->Branch("zmax",&_zmax,"zmax/F");
    _bcdiag->Branch("zgap",&_zgap,"zgap/F");
    _bcdiag->Branch("deltaPhi",&_deltaPhi,"deltaPhi/F");
    // _bcdiag->Branch("slope",&_slope,"slope/F");
    _bcdiag->Branch("zdiff",&_zdiff,"zdiff/F");
    _bcdiag->Branch("phidiff",&_phidiff,"phidiff/F");
    _bcdiag->Branch("zgap2",&_zgap2,"zgap2/F");
    _bcdiag->Branch("np",&_np,"np/I");
    _bcdiag->Branch("pfrac",&_pfrac,"pfrac/F");
    _bcdiag->Branch("kQ",&_kQ,"kQ/F");
    _bcdiag->Branch("fp",&_fp,"fp/I");
    _bcdiag->Branch("lp",&_lp,"lp/I");
    _bcdiag->Branch("pgap",&_pgap,"pgap/I");
    _bcdiag->Branch("nhits",&_nhits,"nhits/I");
    // mc truth branches
    if(_mcdiag){
      _bcdiag->Branch("mmom",&_mmom);
      _bcdiag->Branch("mopos",&_mopos);
      _bcdiag->Branch("mpdg",&_mpdg,"mpdg/I");
      _bcdiag->Branch("mproc",&_mproc,"mproc/I");
      _bcdiag->Branch("prel",&_prel,"prel/I");
      _bcdiag->Branch("code",&_code,"code/I");
      _bcdiag->Branch("frac",&_frac,"frac/F");
      _bcdiag->Branch("nmain",&_nmain,"nmain/I");
      _bcdiag->Branch("nsibling",&_nsibling,"nsibling/I");
      _bcdiag->Branch("nrel",&_nrel,"nrel/I");
      _bcdiag->Branch("nconv",&_nconv,"nconv/I");
      _bcdiag->Branch("ndelta",&_ndelta,"ndelta/I");
      _bcdiag->Branch("ncompt",&_ncompt,"ncompt/I");
      _bcdiag->Branch("ngconv",&_ngconv,"ngconv/I");
      _bcdiag->Branch("nebkg",&_nebkg,"nebkg/I");
      _bcdiag->Branch("nprot",&_nprot,"nprot/I");
      _bcdiag->Branch("ncontrib",&_ncontrib,"ncontrib/I");
      _bcdiag->Branch("icontrib",&_icontrib,"icontrib[ncontrib]/I");
      //_bcdiag->Branch("truebkghits",&_truebkghits);
    }
    if(_hdiag){
      _bhdiag = tfs->make<TTree>("bkghdiag","background hit diagnostics");
      _bhdiag->Branch("iev",        &_iev,          "iev/I");
      _bhdiag->Branch("nhits",      &_nhits,        "nhits/I");
      _bhdiag->Branch("pos",     &_hitPos);
      _bhdiag->Branch("time",    &_hitTime,      "hitTime[nhits]/F");
      _bhdiag->Branch("nch",  &_hitnch,    "hitnch[nhits]/I");
      _bhdiag->Branch("nsh",  &_hitnsh,    "hitnsh[nhits]/I");
      if(_mcdiag){
        _bhdiag->Branch("mcpdg",   &_hitPdg,       "hitPdg[nhits]/I");
        _bhdiag->Branch("mcproc",&_hitproc,    "hitproc[nhits]/I");
      }
    }
  }

  void BkgDiag::analyze(const art::Event& event ) {
    if(!findData(event))
      throw cet::exception("RECO")<<"mu2e::BkgDiag: data missing or incomplete"<< std::endl;
    if( !(_chcol->level() == expectedCHLevel_) ){
      throw cet::exception("RECO")<< "mu2e::BkgDiag: inconsistent outputlevel with input combo hits.\n"
                                  << "FlagBkgHits outputlevel must be "<< expectedCHLevelName_ <<" for training purposes."<< std::endl;
    }
    // loop over background clusters

    _nhits=0;
    _hitPos.clear();
    _hitPos.reserve(_chcol->size());
    for(size_t ich=0;ich<_chcol->size();++ich){
      _hitPos.push_back(_chcol->at(ich).pos());
      _hitTime[_nhits]   = _chcol->at(ich).time();
      _hitnch[_nhits] = _chcol->at(ich).nCombo();
      _hitnsh[_nhits] = _chcol->at(ich).nStrawHits();
      art::Ptr<SimParticle> spp; // main SimParticle for this
      if(_mcdiag){
        std::vector<StrawDigiIndex> dids;
        _chcol->fillStrawDigiIndices(ich,dids);
        StrawDigiMC const& mcdigi = _mcdigis->at(dids[0]);// taking 1st digi: is there a better idea??
        art::Ptr<SimParticle> const& spp = mcdigi.earlyStrawGasStep()->simParticle();
        _hitPdg[_nhits] = spp->pdgId();
        _hitproc[_nhits] = spp->creationCode();
      }
      ++_nhits;
    }

    if(_hdiag)_bhdiag->Fill();

    _cluIdx=0;
    //std::cout<<"bkgccol size = "<<_bkgccol->size()<<std::endl;
    for (size_t ibkg=0;ibkg<_bkgccol->size();++ibkg){
      BkgCluster const& cluster = _bkgccol->at(ibkg);
      // fill cluster info
      _kQ = cluster.getKerasQ();
      _crho = sqrtf(cluster.pos().perp2());
      _cpos = cluster.pos();
      _ctime = cluster.time();
      _cedep = cluster.edep();
      if(cluster.getDistMethod() == BkgCluster::chi2){
        _cchi2 = cluster.points().chisquared();
        _ccons = cluster.points().consistency();
        _cndof = cluster.points().nDOF();
      }
      _isinit = cluster.flag().hasAllProperties(BkgClusterFlag::init);
      _isbkg = cluster.flag().hasAllProperties(BkgClusterFlag::bkg);
      _isref = cluster.flag().hasAllProperties(BkgClusterFlag::refined);
      _isolated = cluster.flag().hasAllProperties(BkgClusterFlag::iso);
      _stereo = cluster.flag().hasAllProperties(BkgClusterFlag::stereo);
      // info on nearest cluster
      _mindt = _mindrho = 1.0e3;
      for(size_t jbkg = 0; jbkg < _bkgccol->size(); ++jbkg){
        if(ibkg != jbkg){
          BkgCluster const& ocluster = _bkgccol->at(jbkg);
          double dt = fabs(ocluster.time() - cluster.time());
          double drho = sqrt((ocluster.pos()-cluster.pos()).Perp2());
          if(dt < _mindt) _mindt = dt;
          if(drho < _mindrho) _mindrho = drho;
        }
      }
      // fill mc info
      art::Ptr<SimParticle> mptr;
      // loop over hits in this cluster and classify them
      _nconv = 0;
      _nprot = 0;
      _nebkg = 0;
      _nmain = 0;
      _mmom = XYZVectorF();
      _mopos = XYZVectorF();
      _mpdg = _mproc = 0;
      _ncontrib = 0;
      _prel=-1;
      if(_mcdiag){
        // fill vector of indices to all digis used in this cluster's hits
        // this goes recursively through the ComboHit chain
        std::vector<StrawDigiIndex> cdids;
        //std::cout<<"Cluster hits size = "<<cluster.hits().size()<<std::endl;
        for(auto const& ich : cluster.hits()){
          // get the list of StrawHit indices associated with this ComboHit
          _chcol->fillStrawDigiIndices(ich,cdids);
        }
        std::vector<int> icontrib;
        findMain(cdids,mptr,_mmom,icontrib);
        findMajority(cdids,mptr,_mmom,_code, _frac);
        _mpdg = mptr->pdgId();
        _mproc = mptr->creationCode();
        _mopos = mptr->startPosXYZ();
        for (int ic : icontrib) {_icontrib[_ncontrib]=ic; ++_ncontrib;}
        if(_code == 12 or _code == 17 or _code == 40 or _code == 165 or _code == 13
           or _code == 74 or _code == 58 or _code == 114 or _code == 166
           or _code == 21 or _code == 115 or _code == 97 or _code == 101 or _code == 78
           or _code == 113 or _code == 16 or _code == 23 or _code == 35 or _code == 72
           or _code == 56  or _code == 2 or _code == 31 or _code == 59  or _code == 109)
          _prel = -2; // This is for a cluster, so not an ACCURATE check
        else if(_code == 167){
          _prel = 1;
          //std::cout<<"CE cluster"<<" mom = "<<_mmom.R()<<" PDG = "<<_mpdg<<std::endl;
        }
        else{
          _prel = 0;
          std::cout<<"Cluster hits size = "<<cluster.hits().size()<<" Background prel = "
                   <<_prel<<" isbkg = "<<_isbkg<<" code = "<<_code<<" pdg = "<<_mpdg<<std::endl;
        }
        /*if(_prel == -2) std::cout<<"No. of true background clusters = "<<_prel<<std::endl;
        if(_prel == 0) std::cout<<"No. of true signal clusters = "<<_prel<<std::endl;
        if(_prel == -2 and _isbkg) std::cout<<"Background classified correctly = "<<_prel<<" "<<_isbkg<<std::endl;
        if(_prel == 0 and !_isbkg) std::cout<<"Signal classified correctly = "<<_prel<<" "<<_isbkg<<std::endl;
        if(_prel == 0 and _isbkg) std::cout<<"Signal classified INCORRECTLY = "<<_prel<<" "<<_isbkg<<std::endl;*/
        //if(_prel == 0) std::cout<<"Analyze:: After find majority = "<<_code<<" frac = "<<_frac<<" pdg = "<<_mpdg<<" mom = "<<_mmom<<" prel = "<<_prel<<std::endl;
        /*if(mptr.isNonnull()){
          _mpdg = mptr->pdgId();
          _mproc = mptr->creationCode();
          _mopos = mptr->startPosXYZ();
          std::cout<<"Analyze:: mpdg = "<<_mpdg<<" creation code = "<<_mproc<<std::endl;
          for(auto const& mcmptr : _mcprimary->primarySimParticles()){
            MCRelationship rel(mcmptr,mptr);
            if(rel.relationship() > MCRelationship::none){
              if(_prel > MCRelationship::none)
                _prel = std::min(_prel,(int)rel.relationship());
              else
                _prel = rel.relationship();
              std::cout<<"prel = "<<_prel<<" mPDG ID = "<<_mpdg<<" creation code = "<<_mproc<<std::endl;
            }
          }
        }*/
      }
      // fill cluster hit info
      _bkghinfo.clear();
      _bkghinfo.reserve(cluster.hits().size());
      _nch = cluster.hits().size();
      /*if(_nch > 0) {
        std::cout<<"Cluster size = "<<_nch<<" code = "<<_code<<" pdg = "<<_mpdg<<std::endl;
        }*/
      _nsh = _nsth = _nactive = _nsha = _nbkg = _nrel = 0;
      //int misCE(0);
      float sumEdep(0.);
      float sumEcc(0.);
      float sqrSumDeltaTime(0.);
      float sqrSumDeltaX(0.);
      float sqrSumDeltaY(0.);
      float sqrSumDeltaPhi(0.);
      float sqrSumQual(0.);
      float sumPitch(0.);
      float sumYaw(0.);
      float sumwPitch(0.);
      float sumwYaw(0.);
      float sumwEcc(0.);
      float phihit(0.), phiclust(0.), phidiff(0.);
      float phimin = std::numeric_limits<float>::max();
      float phimax = -std::numeric_limits<float>::max();
      std::vector<int> strawIds;
      std::vector<int> panelIds;
      std::vector<float> hz;
      std::vector<float> hphi;
      std::array<bool,StrawId::_nplanes> hp{false};
      phiclust = _cpos.phi();
      if(phiclust < -M_PI) phiclust +=2*M_PI;
      if(phiclust > M_PI) phiclust -=2*M_PI;
      for(auto const& ich : cluster.hits()){
        ComboHit const& ch = _chcol->at(ich);
        hphi.push_back(ch.pos().phi());
        hz.push_back(ch.pos().Z());
        hp[ch.strawId().plane()] = true;
        BkgClusterHit const& bhit = _bkghitcol->at(ich);
        sumEdep +=  ch.energyDep()/ch.nStrawHits();
        sqrSumDeltaX += std::pow(ch.pos().X() - _cpos.X(),2);
        sqrSumDeltaY += std::pow(ch.pos().Y() - _cpos.Y(),2);
        phihit = ch.pos().phi();
        phidiff = phihit- phiclust;
        if(phidiff > M_PI) phidiff -= 2*M_PI;
        if(phidiff < -M_PI) phidiff += 2*M_PI;
        if(phidiff < phimin) phimin = phidiff;
        if(phidiff > phimax) phimax = phidiff;
        sqrSumDeltaPhi += std::pow(phidiff,2);
        sqrSumDeltaTime += std::pow(ch.time() - _ctime,2);
        auto hdir = ch.hDir();
        auto wecc = ch.nStrawHits();
        sumEcc += std::sqrt(1-(ch.vVar()/ch.uVar()))*wecc;
        sumwEcc += wecc;
        if(strawIds.size() == 0 && panelIds.size() == 0){
          strawIds.push_back(ch.strawId().straw());
          panelIds.push_back(ch.strawId().uniquePanel());
        }
        else{
          int found = 0;
          for(size_t i = 0;i<strawIds.size();i++){
            if (strawIds.at(i) == ch.strawId().straw() && panelIds.at(i) == ch.strawId().uniquePanel()){
              found = 1;
            }
          }
          if(found == 0){
            strawIds.push_back(ch.strawId().straw());
            panelIds.push_back(ch.strawId().uniquePanel());
          }
        }
        if(ch.flag().hasAllProperties(StrawHitFlag::sline)){
          //quality of SLine fit
          sqrSumQual += std::pow(ch.qual(),2);

          //angle with Mu2e-Y
          float varPitch = std::pow(TMath::ACos(std::sqrt(ch.hcostVar())),2);
          float wPitch = 1/varPitch;
          float signPitch = hdir.Y()/std::abs(hdir.Y());
          sumPitch += signPitch*wPitch*hdir.theta();
          sumwPitch += wPitch;

          ROOT::Math::XYZVectorF z = {0,0,1};
          ROOT::Math::XYZVectorF dxdz = {hdir.X(),0,hdir.Z()};
          float magdxdz = std::sqrt(dxdz.Mag2());

          // angle with Mu2e-Z
          float varYaw = std::sqrt(ch.hphiVar() + varPitch);
          float wYaw = 1/varYaw;
          float signYaw = hdir.X()/std::abs(hdir.X());
          sumYaw += signYaw*wYaw*TMath::ACos(dxdz.Dot(z)/magdxdz);
          sumwYaw += wYaw;

          // # of stereo hits with SLine
          _nsth++;
        }
        _nsh += ch.nStrawHits();
        StrawHitFlag const& shf = bhit.flag();
        //std::cout<<"BkgDiag strawhit flag = "<<shf<<std::endl;
        //std::cout<<"BkgDiag ch flag = "<<ch.flag()<<std::endl;
        if(shf.hasAllProperties(StrawHitFlag::active)){
          _nactive += ch.nStrawHits();
          if(StrawIdMask::station == ch._mask) _nsha+= ch.nStrawHits();
        }
        if(shf.hasAllProperties(StrawHitFlag::bkg))_nbkg+= ch.nStrawHits();
        // fill hit-specific information
        BkgHitInfo bkghinfo;
        // fill basic straw hit info
        fillStrawHitInfo(ich,bkghinfo);
        if(_mcdiag){
          std::vector<StrawDigiIndex> dids;
          _chcol->fillStrawDigiIndices(ich,dids);
          StrawDigiMC const& mcdigi = _mcdigis->at(dids[0]);// taking 1st digi: is there a better idea??
          fillStrawHitInfoMC(mcdigi,mptr,bkghinfo);
          //global counting for the cluster: count signal hits only, but background from background is OK
          if(bkghinfo._mrel==MCRelationship::same){
            _nmain += ch.nStrawHits(); // same as main particle
          } else if(bkghinfo._mrel>MCRelationship::sibling)
            _nsibling += ch.nStrawHits();
          if(bkghinfo._mrel>MCRelationship::none)_nrel += ch.nStrawHits();
          if(bkghinfo._mcpdg == PDGCode::proton)_nprot += ch.nStrawHits();
          /*if(bkghinfo._prel == -2 and _prel == -2) _truebkghits += ch.nStrawHits();
          if((bkghinfo._prel == 0 or bkghinfo._prel == 1) and _prel == 0) _truesighits += ch.nStrawHits();
          if(bkghinfo._prel == -2) _totalbkghits += ch.nStrawHits();
          if(bkghinfo._prel == 0 or bkghinfo._prel == 1) _totalsighits +=ch.nStrawHits();
          if(bkghinfo._prel == 1 and _prel == 1) _trueCEhits +=ch.nStrawHits();
          if(bkghinfo._prel == 1) _totalCEhits +=ch.nStrawHits();
          if(bkghinfo._prel == 1 and _prel == -2) {
            _misCEhits +=ch.nStrawHits();
            misCE +=ch.nStrawHits();
          }*/
        }
        bkghinfo._active = shf.hasAllProperties(StrawHitFlag::active);
        bkghinfo._cbkg = shf.hasAllProperties(StrawHitFlag::bkg);
        //std::cout<<"Bkg hit info = "<<bkghinfo._cbkg<<" hit prel = "<<bkghinfo._prel<<" prel = "<<_prel<<" nbkg = "<<_nbkg<<std::endl;
        bkghinfo._gdist = bhit.distance();
        bkghinfo._index = ich;
        // calculate separation to cluster
        auto psep = ch.pos()-cluster.pos();
        auto pdir = PerpVector(psep,GenVector::ZDir()).Unit();
        bkghinfo._rpos = psep;
        bkghinfo._rerr = std::max(float(2.5),ch.posRes(ComboHit::wire)*fabs(pdir.Dot(ch.uDir())));
        _bkghinfo.push_back(bkghinfo);
      }
      _avecedep = sumEdep/_nch;
      _cecc = sumEcc/sumwEcc;
      _rmscposx = std::sqrt(sqrSumDeltaX/_nch);
      _rmscposy = std::sqrt(sqrSumDeltaY/_nch);
      _rmscrho = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/_nch);
      _rmscphi = std::sqrt(sqrSumDeltaPhi/_nch);
      _rmsctime = std::sqrt(sqrSumDeltaTime/_nch);
      float areaXY = M_PI*_rmscposx*_rmscposy;
      if(areaXY > 0)
        _clusterdensity = static_cast<float>(_nch)/areaXY;
      else
        _clusterdensity = 0.0;
      _cpitch = _nsth > 0 ? sumPitch/sumwPitch : 0.;
      _cyaw = _nsth > 0 ? sumYaw/sumwYaw : 0.;
      _cqual = std::sqrt(sqrSumQual/_nch);//average SLine fit quality of the cluster
      _csthqual = float(_nsth)/float(_nch);//SLine ratio the cluster
      _ccomqual = _cqual + _csthqual;//combined quality metric
      std::sort(hz.begin(),hz.end());
      _zgap = 0.0;
      _deltaPhi = 0.0;
      //_slope = 0.0;
      _zdiff = 0.0;
      _phidiff = 0.0;
      _zgap2 = 0.0;
      for (unsigned iz=1;iz<hz.size();++iz){
        _zgap = std::max(_zgap,hz[iz]-hz[iz-1]);
        _deltaPhi = _deltaPhi + (hphi[iz]-hphi[iz-1]); 
        _zgap2 = _zgap2 + (hz[iz] - hz[iz-1]);
      }
      _zgap2 = _zgap2/hz.size();
      _deltaPhi = _deltaPhi/hz.size();
      //std::cout<<"Delta phi = "<<_deltaPhi<<std::endl;
      _zmin = hz.front();
      _zmax = hz.back();
      //if(_zgap == 0 or _zgap2 == 0 or _zdiff == 0) std::cout<<_zgap<<" "<<_zgap2<<" "<<_zdiff<<std::endl;
      if(_nch > 2) _zdiff = _zmax - _zmin;
      if(_nch > 2) _phidiff = phimax - phimin; //hphi.back() - hphi.front();
      //std::cout<<"zdiff = "<<_zdiff<<" phidiff = "<<_phidiff<<std::endl;
      /*if(misCE > 0) std::cout<<"Mis CE = "<<misCE<<" Cluster rho = "<<_crho<<" zdiff = "<<_zdiff<<" phidiff = "
               <<_phidiff<<" nhits = "<<_nsh<<" RMS rho = "<<_rmscrho<<" RMS time = "
               <<_rmsctime<<" RMS phi = "<<_rmscphi<<" keras Q = "<<_kQ
               <<" pos = "<<_cpos.X()<<" "<<_cpos.Y()<<" "<<_cpos.Z()<<" event = "<<_iev<<" code = "<<_code<<std::endl;*/
      _lp = -1; // last plane in cluster
      _fp = StrawId::_nplanes; // first plane in cluster
      _np = 0;//# of planes
      _pgap = 0; // largest plane gap
      int lp(-1); // last plane seen
      for(int ip=0;ip < StrawId::_nplanes; ++ip){
        if(hp[ip]){
          _np++;
          if(lp > 0 && ip - lp -1 > _pgap)_pgap = ip - lp -1;
          if(ip > _lp)_lp = ip;
          if(ip < _fp)_fp = ip;
          lp = ip;
        }
      }
      _pfrac = static_cast<float>(_np)/static_cast<float>(_lp - _fp);
      _bcdiag->Fill();
      ++_cluIdx;
    }
    //std::cout<<"Truebkghitinfo = "<<_truebkghits<<" sighits = "<<_truesighits<<" CE hits = "
    //         <<_trueCEhits<<" Total bkg hits = "<<_totalbkghits<<" total sig hits = "<<_totalsighits
    //         <<" Total CE hits = "<<_totalCEhits<<" mis CE hits = "<<_misCEhits<<std::endl;   
    ++_iev;
  }

  bool BkgDiag::findData(const art::Event& evt){
    _chcol = 0; _bkgccol = 0; _mcdigis = 0;
    // nb: getValidHandle does the protection (exception) on handle validity so I don't have to
    auto chH = evt.getValidHandle(_chToken);
    _chcol = chH.product();
    auto bkgcH = evt.getValidHandle(_bkgcToken);
    _bkgccol = bkgcH.product();
    auto bkghH = evt.getValidHandle(_bkghToken);
    _bkghitcol = bkghH.product();
    if(_mcdiag){
      auto mcdH = evt.getValidHandle(_mcdigisToken);
      _mcdigis = mcdH.product();
      auto mcpH = evt.getValidHandle(_mcprimaryToken);
      _mcprimary = mcpH.product();
    }
    return _chcol != 0 && _bkgccol != 0
      && ( (_mcdigis != 0 && _mcprimary != 0)  || !_mcdiag);
  }

  void BkgDiag::findMajority(std::vector<uint16_t>const& dids, art::Ptr<SimParticle>& mptr,XYZVectorF& mmom, int& mostCommonCode, float& fraction)
    const {
    // find the unique simparticles which produced these hits
    mostCommonCode = -1;
    fraction = 0.0;
    std::map<int, unsigned> codeCount;
    for(auto id : dids) {
      StrawDigiMC const& mcdigi = _mcdigis->at(id);
      art::Ptr<SimParticle> const& spp = mcdigi.earlyStrawGasStep()->simParticle();
      if(spp.isNonnull()){
        int creationCode = spp->creationCode();
        //if(dids.size() < 2) std::cout<<"cluster has "<<dids.size()<<" creation code = "<<creationCode<<" pdg ID = "<<spp->pdgId()<<std::endl;
        ++codeCount[creationCode];
        //std::cout<<"findMajority::MCDigi ID = "<<id<<" simParticle = "<<spp->pdgId()<<" creation code = "
        //         <<spp->creationCode()<<std::endl;
      }
    }
    // Find most common code
    unsigned maxCount = 0;
    for(auto const& kv : codeCount) {
      if(kv.second > maxCount) {
        maxCount = kv.second;
        mostCommonCode = kv.first;
      }
    }
    fraction = static_cast<float>(maxCount)/static_cast<float>(dids.size());
    // Second pass: find momentum of first matching digi
    for(auto id : dids) {
      StrawDigiMC const& mcdigi = _mcdigis->at(id);
      auto const& sgsp = mcdigi.earlyStrawGasStep();
      art::Ptr<SimParticle> const& spp = sgsp->simParticle();;
      if(spp.isNull()) continue;
      mptr = spp;
      if(spp->creationCode() == mostCommonCode) {
        mmom = sgsp->momentum();
        break;
      }
    }
    //std::cout<<"findMajority::most common code = "<<mostCommonCode<<" fraction = "<<fraction<<" mom = "<<mmom<<" pdg = "<<mptr->pdgId()<<std::endl;
  }

  void BkgDiag::findMain(std::vector<uint16_t>const& dids, art::Ptr<SimParticle>& mptr,XYZVectorF& mmom, std::vector<int>& icontrib) const {
    // find the unique simparticles which produced these hits
    std::set<art::Ptr<SimParticle> > pp;
    for(auto id : dids) {
      StrawDigiMC const& mcdigi = _mcdigis->at(id);
      art::Ptr<SimParticle> const& spp = mcdigi.earlyStrawGasStep()->simParticle();
      if(spp.isNonnull()){
        pp.insert(spp);
        //std::cout<<"findMain::MCDigi ID = "<<id<<" simParticle = "<<spp->pdgId()<<" creation code = "
        //         <<spp->creationCode()<<std::endl;
      }
    }
    // map these particles back to each other, to compress out particles generated inside the cluster
    std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> > spmap;
    // look for particles produced at the same point, like conversions.  It's not enough to look for the same parent,
    // as that parent could produce multiple daughters at different times.  Regardless of mechanism or genealogy, call these 'the same'
    // as they will contribute equally to the spiral
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> sppi = *ipp;
      spmap[sppi] = sppi;
    }
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> sppi = *ipp;
      if(sppi->genParticle().isNull()){
        std::set<art::Ptr<SimParticle> >::iterator jpp=ipp;++jpp;
        for(;jpp!=pp.end();++jpp){
          art::Ptr<SimParticle> sppj = *jpp;
          if(sppj->genParticle().isNull()){
            // call the particles 'the same' if they are related and were produced near each other
            MCRelationship rel(sppi,sppj);
            if(rel==MCRelationship::daughter || rel == MCRelationship::udaughter){
              spmap[sppi] = sppj;
              break;
            } else if(rel == MCRelationship::mother || rel == MCRelationship::umother){
              spmap[sppj] = sppi;
            } else if(rel == MCRelationship::sibling || rel == MCRelationship::usibling){
              double dist = (sppj->startPosition() - sppi->startPosition()).mag();
              if(dist < 10.0){
                if(sppi->id().asInt() > sppj->id().asInt())
                  spmap[sppi] = sppj;
                else
                  spmap[sppj] = sppi;
              }
            }
          }
        }
      }
    }
    // check for remapping
    bool changed(true);
    while(changed){
      changed = false;
      for(std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> >::iterator im = spmap.begin();im!=spmap.end();++im){
        std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> >::iterator ifnd = spmap.find(im->second);
        if( !(ifnd->second == ifnd->first)){
          changed = true;
          spmap[im->first] = ifnd->second;
        }
      }
    }
    // find the most likely ultimate parent for this cluster.  Also fill general info
    std::map<int,int> mode;
    for(std::set<art::Ptr<SimParticle> >::iterator ipp=pp.begin();ipp!=pp.end();++ipp){
      art::Ptr<SimParticle> spp = *ipp;
      int mcid(-1);
      // map back to the ultimate parent
      spp = spmap[spp];
      mcid = spp->id().asInt();
      std::map<int,int>::iterator ifnd = mode.find(mcid);
      if(ifnd != mode.end())
        ++(ifnd->second);
      else
        mode[mcid] = 1;
    }
    int max(0);
    std::map<int,int>::iterator imax = mode.end();
    for(std::map<int,int>::iterator im=mode.begin();im!=mode.end();++im){
      icontrib.push_back(im->first);
      if(im->second>max){
        imax=im;
        max = im->second;
      }
    }
    unsigned pid(0);
    if(imax != mode.end())
      pid=imax->first;
    for(std::map<art::Ptr<SimParticle>,art::Ptr<SimParticle> >::iterator im = spmap.begin();im!=spmap.end();++im){
      if(im->first->id().asInt() == pid){
        mptr = im->first;
        break;
      }
    }
    // find the momentum for the first step point from the primary particle in this delta
    for(auto id : dids) {
      StrawDigiMC const& mcdigi = _mcdigis->at(id);
      auto const& sgsp = mcdigi.earlyStrawGasStep();
      art::Ptr<SimParticle> const& spp = sgsp->simParticle();
      if(spp == mptr){
        mmom = sgsp->momentum();
        //std::cout<<"mmom = "<<mmom<<" pdgID = "<<spp->pdgId()<<" creation code = "<<spp->creationCode()<<std::endl;
        break;
      }
    }
  }

  void BkgDiag::fillStrawHitInfoMC(StrawDigiMC const& mcdigi, art::Ptr<SimParticle>const& mptr, StrawHitInfo& shinfo) const {
    // use early step to define the MC match
    auto const& sgsp = mcdigi.earlyStrawGasStep();
    art::Ptr<SimParticle> const& spp = sgsp->simParticle();
    shinfo._mctime = sgsp->time();
    shinfo._mcht = mcdigi.wireEndTime(mcdigi.earlyEnd());
    shinfo._mcpdg = spp->pdgId();
    shinfo._mcproc = spp->creationCode();
    shinfo._mcedep = mcdigi.energySum();
    shinfo._mcpos = sgsp->position();
    shinfo._mcmom = sqrt(sgsp->momentum().mag2());
    double cosd = cos(sgsp->momentum().Theta());
    shinfo._mctd = cosd/sqrt(1.0-cosd*cosd);
    // relationship to main particle
    shinfo._mrel=MCRelationship::none;
    if(sgsp.isNonnull() && mptr.isNonnull()){
      art::Ptr<SimParticle> const& spp = sgsp->simParticle();
      if(spp.isNonnull()){
        MCRelationship rel(spp,mptr);
        int code = spp->creationCode();
        //int momcreationCode = mptr->creationCode();
        //std::cout<<"fillStrawHitInfoMC creationCode = "<<code<<" mother = "<<momcreationCode<<std::endl;
        if(code == 12 or code == 17 or code == 40 or code == 165 or code == 13
           or code == 74 or code == 58 or code == 114 or code == 166
           or code == 21 or code == 115 or code == 97 or code == 101 or code == 78
           or code == 113 or code == 16 or code == 23 or code == 35 or code == 72
           or code == 56  or code == 2 or code == 31 or code == 59  or code == 109)
          shinfo._prel = -2; //Background hit
        else if(code == 167)
          shinfo._prel = 1;  // CE hit
        else{
          shinfo._prel = 0; //CE like hit

        }
        shinfo._mrel = rel.relationship();
        /*for(auto const& mcmptr : _mcprimary->primarySimParticles()){
          MCRelationship rel(spp,mcmptr);
          if(rel.relationship() > MCRelationship::none){
            if(shinfo._prel > MCRelationship::none)
              shinfo._prel = std::min(shinfo._prel,(int)rel.relationship());
            else
              shinfo._prel = rel.relationship();
          }
          }*/
        
      }
    }
  }

  void BkgDiag::fillStrawHitInfo(size_t ich, StrawHitInfo& shinfo) const {
    ComboHit const& ch = _chcol->at(ich);
    auto const& shf = ch.flag();

    shinfo._stereo = shf.hasAllProperties(StrawHitFlag::stereo);
    shinfo._tdiv = shf.hasAllProperties(StrawHitFlag::tdiv);
    shinfo._esel = shf.hasAllProperties(StrawHitFlag::energysel);
    shinfo._rsel = shf.hasAllProperties(StrawHitFlag::radsel);
    shinfo._tsel = shf.hasAllProperties(StrawHitFlag::timesel);
    shinfo._strawxtalk = shf.hasAllProperties(StrawHitFlag::strawxtalk);
    shinfo._elecxtalk = shf.hasAllProperties(StrawHitFlag::elecxtalk);
    shinfo._isolated = shf.hasAllProperties(StrawHitFlag::isolated);
    shinfo._bkg = shf.hasAllProperties(StrawHitFlag::bkg);
    shinfo._bkgc = shf.hasAllProperties(StrawHitFlag::bkgclust);

    shinfo._pos = ch.pos();
    shinfo._time = ch.correctedTime();
    shinfo._wdist = ch.wireDist();
    shinfo._wres = ch.posRes(ComboHit::wire);
    shinfo._tres = ch.posRes(ComboHit::trans);
    // info depending on stereo hits
    shinfo._chisq = ch.qual();
    shinfo._edep = ch.energyDep();
    StrawId const& sid = ch.strawId();
    shinfo._plane = sid.plane();
    shinfo._panel = sid.panel();
    shinfo._layer = sid.layer();
    shinfo._straw = sid.straw();
    shinfo._stereo = ch.flag().hasAllProperties(StrawHitFlag::stereo);
    shinfo._tdiv = ch.flag().hasAllProperties(StrawHitFlag::tdiv);
  }
} // mu2e namespace

// Part of the magic that makes this class a module.
using mu2e::BkgDiag;
DEFINE_ART_MODULE(BkgDiag)
