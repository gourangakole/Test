// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/CommonPlots_lt.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/TransverseMass.h"
#include "Math/GenVector/VectorUtil.h"

#include "TDirectory.h"
#include "TMatrixD.h"
#include "TMatrixT.h"

class TauFakeRate_mt_MET: public BaseSelector {
public:
  explicit TauFakeRate_mt_MET(const ParameterSet& config, const TH1* skimCounters);
  virtual ~TauFakeRate_mt_MET() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;

  std::pair<double,double> METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, int year, bool isMC, int npv);
private:
  // Common plots
  CommonPlots_lt fCommonPlots;
  //CommonPlots fCommonPlots;

  // Event selection classes and event counters (in same order like they are applied)
  Count cAllEvents;
  Count cTrigger;
  METFilterSelection fMETFilterSelection;
  Count cVertexSelection;
  ElectronSelection fElectronSelection;
  MuonSelection fMuonSelection;
  Count cLeptonOSCounter;
  Count cLeptonMassCounter;
  TauSelection fTauSelection;
  TauSelection fLooseTauSelection;
  Count cTauNCounter;
  Count cTauSFCounter;
  Count cFakeTauSFCounter;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  // FatJetSelection fFatJetSelection;
  Count cSelected;
    
  void doLooseTaus(const Event& event, const TauSelection::Data& tauData);
  void doTightTaus(const Event& event, const TauSelection::Data& looseTauData);
  double CollinearMass( const MuonSelection::Data& muData, const TauSelection::Data& looseTauData, const METSelection::Data& metData);

  // Run
  WrappedTH1* hNSelectedVsRunNumber_Cuts; // For data only

  // Non-common histograms
  WrappedTH1 *hTauPt_num_dm0;
  WrappedTH1 *hTauPt_num_dm1;
  WrappedTH1 *hTauPt_num_dm10;
  WrappedTH1 *hTauPt_den_dm0;
  WrappedTH1 *hTauPt_den_dm1;
  WrappedTH1 *hTauPt_den_dm10;
  WrappedTH1 *hTauPt_num_g_dm0;
  WrappedTH1 *hTauPt_num_g_dm1;
  WrappedTH1 *hTauPt_num_g_dm10;
  WrappedTH1 *hTauPt_den_g_dm0;
  WrappedTH1 *hTauPt_den_g_dm1;
  WrappedTH1 *hTauPt_den_g_dm10;
  WrappedTH1 *hTauEta_num_dm0;
  WrappedTH1 *hTauEta_num_dm1;
  WrappedTH1 *hTauEta_num_dm10;
  WrappedTH1 *hTauEta_den_dm0;
  WrappedTH1 *hTauEta_den_dm1;
  WrappedTH1 *hTauEta_den_dm10;

  // Barrel ( |eta| < 1.5)
  WrappedTH1 *hTauPt_num_dm0_barrel;
  WrappedTH1 *hTauPt_num_dm1_barrel;
  WrappedTH1 *hTauPt_num_dm10_barrel;
  WrappedTH1 *hTauPt_den_dm0_barrel;
  WrappedTH1 *hTauPt_den_dm1_barrel;
  WrappedTH1 *hTauPt_den_dm10_barrel;
  WrappedTH1 *hTauPt_num_g_dm0_barrel;
  WrappedTH1 *hTauPt_num_g_dm1_barrel;
  WrappedTH1 *hTauPt_num_g_dm10_barrel;
  WrappedTH1 *hTauPt_den_g_dm0_barrel;
  WrappedTH1 *hTauPt_den_g_dm1_barrel;
  WrappedTH1 *hTauPt_den_g_dm10_barrel;
  WrappedTH1 *hTauEta_num_dm0_barrel;
  WrappedTH1 *hTauEta_num_dm1_barrel;
  WrappedTH1 *hTauEta_num_dm10_barrel;
  WrappedTH1 *hTauEta_den_dm0_barrel;
  WrappedTH1 *hTauEta_den_dm1_barrel;
  WrappedTH1 *hTauEta_den_dm10_barrel;

  // Endcap ( |eta| >= 1.5)
  WrappedTH1 *hTauPt_num_dm0_endcap;
  WrappedTH1 *hTauPt_num_dm1_endcap;
  WrappedTH1 *hTauPt_num_dm10_endcap;
  WrappedTH1 *hTauPt_den_dm0_endcap;
  WrappedTH1 *hTauPt_den_dm1_endcap;
  WrappedTH1 *hTauPt_den_dm10_endcap;
  WrappedTH1 *hTauPt_num_g_dm0_endcap;
  WrappedTH1 *hTauPt_num_g_dm1_endcap;
  WrappedTH1 *hTauPt_num_g_dm10_endcap;
  WrappedTH1 *hTauPt_den_g_dm0_endcap;
  WrappedTH1 *hTauPt_den_g_dm1_endcap;
  WrappedTH1 *hTauPt_den_g_dm10_endcap;
  WrappedTH1 *hTauEta_num_dm0_endcap;
  WrappedTH1 *hTauEta_num_dm1_endcap;
  WrappedTH1 *hTauEta_num_dm10_endcap;
  WrappedTH1 *hTauEta_den_dm0_endcap;
  WrappedTH1 *hTauEta_den_dm1_endcap;
  WrappedTH1 *hTauEta_den_dm10_endcap;

  WrappedTH1 *hLepLooseTauMass_BeforeOnZSelection;
  WrappedTH1 *hLepLooseTauMass_AfterDRSelection;
  WrappedTH1 *hLepLooseTauMass_METProjZSelection;
  WrappedTH1 *hLepLooseTauCollMass_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMassFlag_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMass40_80_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMassFlag_mass_40_80_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMassDiobjPt_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_AfterMetSelection;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_AfterMetSelection;

  // dm =0 ; barrel and endcap
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_endcap;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_dm0_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_dm0_endcap;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_dm0_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_dm0_endcap;
  // dm =1 ; barrel and endcap
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_endcap;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_dm1_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_dm1_endcap;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_dm1_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_dm1_endcap;
  // dm =10 ; barrel and endcap
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_endcap;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_dm10_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinSS_dm10_endcap;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_dm10_barrel;
  WrappedTH1 *hLepLooseTauCollMassDiobjPtCollmassBinOS_dm10_endcap;

  WrappedTH1 *hLepLooseTauCollMassGenuTau_AfterMetSelection;
  WrappedTH2 *hLepLooseTauPtVsLepLooseTauMass_AfterMetSelection;
  WrappedTH2 *hLepLooseTauPtVsLepLooseTauCollMass_AfterMetSelection;
  WrappedTH1 *hLepLooseTauMassGenuTau_AfterMetSelection;
  WrappedTH1 *hMet_AfterMetSelection;
  WrappedTH1 *hMetPhi_AfterMetSelection;
  WrappedTH1 *hCorrMet_AfterMetSelection;
  WrappedTH1 *hCorrMetPhi_AfterMetSelection;
  WrappedTH1 *hLepLooseTauPt_AfterMetSelection;

  // MET after lepton
  WrappedTH1 *hMet_AfterLepton;
  WrappedTH1 *hMetPhi_AfterLepton;
  WrappedTH1 *hCorrMet_AfterLepton;
  WrappedTH1 *hCorrMetPhi_AfterLepton;

  WrappedTH1 *hLeptonN_AfterLeptonSelection;
  WrappedTH1 *hLeptonPt_AfterLeptonSelection;
  WrappedTH1 *hLeptonEta_AfterLeptonSelection;
  WrappedTH1 *hLepLooseTauPt_AfterLeptonSelection;
  WrappedTH1 *hLepLooseTauEta_AfterLeptonSelection;
  WrappedTH1 *hLepLooseTauMass_AfterLeptonSelection;
  
  WrappedTH1 *hNTau_AfterTauSelection;
  WrappedTH1 *hNTau_AfterJetSelection;
  WrappedTH1 *hNTau_AfterBJetSelection;
  WrappedTH1 *hNTau_AfterMetSelection;

  WrappedTH1 *hLepLooseTauMass_AfterTauSelection;
  WrappedTH1 *hLepLooseTauMass_AfterJetSelection;
  WrappedTH1 *hLepLooseTauMass_AfterBJetSelection;
  WrappedTH1 *hLepLooseTauMass_AfterMetSelection;
  WrappedTH1 *hLepLooseTauMass40_80_AfterMetSelection;

  WrappedTH1 *hDPhiTauMet_AfterMetSelection;
  WrappedTH1 *hDPhiTauLep_AfterMetSelection;
  WrappedTH1 *hDPhiMetLep_AfterMetSelection;

  // new sets of plots for MET corrections
  WrappedTH1 *hdeltaPhiMuMu_Cuts;
  WrappedTH1 *hptRatioMuMu_Cuts;
  WrappedTH1 *hMet_Cuts;
  WrappedTH1 *hMetPhi_Cuts;
  WrappedTH1 *hCorrMet_Cuts;
  WrappedTH1 *hCorrMetPhi_Cuts;
  WrappedTH1 *hMuMuMass_Cuts;

  //plot after medium Tau
  WrappedTH1 *hLepMedTauMass_AfterMetSelection;
  WrappedTH1 *hLepMedTauPt_AfterMetSelection;
  WrappedTH1 *hLepMedTauCollMass_AfterMetSelection;
  WrappedTH1 *hLepMedTauCollMassDiobjPt_AfterMetSelection;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_AfterMetSelection;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_AfterMetSelection;
  WrappedTH2 *hLepMedTauPtVsLepMedTauCollMass_AfterMetSelection;
  // dm=0; barrel and endcap
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_dm0_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_dm0_endcap;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_dm0_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_dm0_endcap;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_dm0_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_dm0_endcap;
  // dm=1; barrel and endcap
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_dm1_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_dm1_endcap;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_dm1_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_dm1_endcap;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_dm1_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_dm1_endcap;
  // dm=10; barrel and endcap
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_dm10_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinIncl_dm10_endcap;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_dm10_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinSS_dm10_endcap;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_dm10_barrel;
  WrappedTH1 *hLepMedTauCollMassDiobjPtCollmassBinOS_dm10_endcap;

  
  WrappedTH1 *hNTau_AfterAllSelections;
  WrappedTH1 *hTauSrc_AfterAllSelections;
  WrappedTH1 *hTauSrcDM0_AfterAllSelections;
  WrappedTH1 *hTauSrcDM1_AfterAllSelections;
  WrappedTH1 *hTauSrcDM10_AfterAllSelections;
  WrappedTH1 *hNJet_AfterAllSelections;
  WrappedTH1 *hNBjet_AfterAllSelections;
  WrappedTH1 *hMet_AfterAllSelections;
  WrappedTH1 *hLepLooseTauPt_AfterAllSelections;
  WrappedTH1 *hLepLooseTauEta_AfterAllSelections;
  WrappedTH1 *hLepLooseTauMass_AfterAllSelections;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(TauFakeRate_mt_MET);

TauFakeRate_mt_MET::TauFakeRate_mt_MET(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    //fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2hwAnalysisWithTop, fHistoWrapper),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots_lt::kTauFakeRateMeasurement, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cLeptonOSCounter(fEventCounter.addCounter("#mu OS")),
    cLeptonMassCounter(fEventCounter.addCounter("m_{#mu#mu} window")),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    // fTauSelection(config.getParameter<ParameterSet>("TauSelection")), // Fixes "An object with name tauSelection_ exists already"
    fLooseTauSelection(config.getParameter<ParameterSet>("LooseTauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cTauNCounter(fEventCounter.addCounter("#geq 1 loose #tau")),
    cTauSFCounter(fEventCounter.addCounter("#tau SF")),
    cFakeTauSFCounter(fEventCounter.addCounter("Fake #tau SF")),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection")),
    //fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b-tag SF")),
    // fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fMETSelection(config.getParameter<ParameterSet>("METSelection")),
    // fFatJetSelection(config.getParameter<ParameterSet>("FatJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    cSelected(fEventCounter.addCounter("Selected Events"))
{ }


void TauFakeRate_mt_MET::book(TDirectory *dir) {

  if (0) std::cout << "=== TauFakeRate_mt_MET::book()" << std::endl;
  // Book common plots histograms
  fCommonPlots.book(dir, isData());

  // Book histograms in event selection classes
  fMETFilterSelection.bookHistograms(dir);
  fElectronSelection.bookHistograms(dir);
  fMuonSelection.bookHistograms(dir);
  fTauSelection.bookHistograms(dir);
  fLooseTauSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);
  fBJetSelection.bookHistograms(dir);
  fMETSelection.bookHistograms(dir);
  // fFatJetSelection.bookHistograms(dir);

  // Get binning from cfg file
  const int   ptN    = fCommonPlots.getPtBinSettings().bins();
  const float ptMin  = fCommonPlots.getPtBinSettings().min();
  const float ptMax  = fCommonPlots.getPtBinSettings().max();
  const int   etaN   = fCommonPlots.getEtaBinSettings().bins();
  const float etaMin = fCommonPlots.getEtaBinSettings().min();
  const float etaMax = fCommonPlots.getEtaBinSettings().max();
  const int   mN     = fCommonPlots.getInvMassBinSettings().bins();
  const float mMin   = fCommonPlots.getInvMassBinSettings().min();
  const float mMax   = fCommonPlots.getInvMassBinSettings().max();
  const int   nN     = fCommonPlots.getNjetsBinSettings().bins();
  const float nMin   = fCommonPlots.getNjetsBinSettings().min();
  const float nMax   = fCommonPlots.getNjetsBinSettings().max();
  const int   metN   = fCommonPlots.getMetBinSettings().bins();
  const float metMin = fCommonPlots.getMetBinSettings().min();
  const float metMax = fCommonPlots.getMetBinSettings().max();
  
  if (isData()) 
    {
      hNSelectedVsRunNumber_Cuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, 
							      "NSelectedVsRunNumber_Cuts", "NSelectedVsRunNumber;Run number;N_{events}", 6000, 270000, 276000); // runnb >=272007 &&runnb<=275376 (B)
    }

  // Book non-common histograms 
  double bin[8] = {20,25,30,35,40,50,60,120}; // iro-fixme

  hTauPt_num_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10", ";p_{T} (GeV)", 7, bin);
  hTauEta_num_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10", ";#eta", etaN, etaMin, etaMax);

  // Barrel ( |eta| < 1.5)
  hTauPt_num_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10_barrel", ";p_{T} (GeV)", 7, bin);
  hTauEta_num_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10_barrel", ";#eta", etaN, etaMin, etaMax);

  // Endcap ( |eta| >= 1.5)
  hTauPt_num_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10_endcap", ";p_{T} (GeV)", 7, bin);
  hTauEta_num_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10_endcap", ";#eta", etaN, etaMin, etaMax);

  hLepLooseTauMass_BeforeOnZSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_BeforeOnZSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauMass_AfterDRSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterDRSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauMass_METProjZSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_METProjZSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMass_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMass_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassFlag_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassFlag_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMass40_80_AfterMetSelection =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMass40_80_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassFlag_mass_40_80_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassFlag_mass_40_80_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPt_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPt_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  // dm=0; barrel and endcap
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm0_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_dm0_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm0_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_dm0_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_dm0_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_dm0_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_dm0_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_dm0_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  // dm=1; barrel and endcap
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm1_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_dm1_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm1_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_dm1_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_dm1_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_dm1_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_dm1_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_dm1_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  // dm=10; barrel and endcap
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm10_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_dm10_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm10_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinSS_dm10_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_dm10_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_dm10_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauCollMassDiobjPtCollmassBinOS_dm10_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassDiobjPtCollmassBinOS_dm10_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);

  hLepLooseTauCollMassGenuTau_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauCollMassGenuTau_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);

  hLepLooseTauPtVsLepLooseTauMass_AfterMetSelection = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, dir, "LepLooseTauPtVsLepLooseTauMass_AfterMetSelection", "LepLooseTauPtVsLepLooseTauMass", ptN , ptMin , ptMax , mN, mMin, mMax);
  hLepLooseTauPtVsLepLooseTauCollMass_AfterMetSelection = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, dir, "LepLooseTauPtVsLepLooseTauCollMass_AfterMetSelection", "LepLooseTauPtVsLepLooseTauCollMass", ptN , ptMin , ptMax , mN, mMin, mMax);

  hLeptonN_AfterLeptonSelection      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonN_AfterLeptonSelection"     , ";lepton multiplicity", nN  , nMin  , nMax );
  hLeptonPt_AfterLeptonSelection     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonPt_AfterLeptonSelection"    , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hLeptonEta_AfterLeptonSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonEta_AfterLeptonSelection"   , ";#eta"       , etaN, etaMin, etaMax);
  hLepLooseTauPt_AfterLeptonSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauPt_AfterLeptonSelection"  , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hLepLooseTauEta_AfterLeptonSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauEta_AfterLeptonSelection" , ";#eta"       , etaN, etaMin, etaMax);
  hLepLooseTauMass_AfterLeptonSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterLeptonSelection", ";m_{l#tau} (GeV)", mN , mMin  , mMax  );

  hLepLooseTauMass_AfterTauSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterTauSelection"  , ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauMass_AfterJetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterJetSelection"  , ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauMass_AfterBJetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterBJetSelection" , ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauMass_AfterMetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterMetSelection"  , ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepLooseTauMassGenuTau_AfterMetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMassGenuTau_AfterMetSelection"  , ";m_{l#tau} (GeV)", mN, mMin, mMax);

  hLepLooseTauMass40_80_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass40_80_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hMet_AfterMetSelection            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_AfterMetSelection", ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hMetPhi_AfterMetSelection         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "MetPhi_AfterMetSelection",";E_{T}^{miss} #phi (rads);Events", 64, -3.2, 3.2);
  hCorrMet_AfterMetSelection            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "CorrMet_AfterMetSelection", ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hCorrMetPhi_AfterMetSelection         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "CorrMetPhi_AfterMetSelection",";E_{T}^{miss} #phi (rads);Events", 64, -3.2, 3.2);
  hLepLooseTauPt_AfterMetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauPt_AfterMetSelection"  , ";p_{T} (GeV)", ptN , ptMin , ptMax );

  // MET after lepton
  hMet_AfterLepton            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_AfterLepton", ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hMetPhi_AfterLepton         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "MetPhi_AfterLepton",";E_{T}^{miss} #phi (rads);Events", 64, -3.2, 3.2);
  hCorrMet_AfterLepton            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "CorrMet_AfterLepton", ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hCorrMetPhi_AfterLepton         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "CorrMetPhi_AfterLepton",";E_{T}^{miss} #phi (rads);Events", 64, -3.2, 3.2);
  

  hDPhiTauMet_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DPhiTauMet_AfterMetSelection", ";#Delta#phi(#tau jet,MET)",180, 0., 180);
  hDPhiTauLep_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DPhiTauLep_AfterMetSelection", ";#Delta#phi(#tau jet,l)",180, 0., 180);
  hDPhiMetLep_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DPhiMetLep_AfterMetSelection", ";#Delta#phi(MET,l)",180, 0., 180);

  hNTau_AfterTauSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterTauSelection"  , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterJetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterJetSelection"  , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterBJetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterBJetSelection" , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterMetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterMetSelection"  , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );

  // plots for MET correction
  hdeltaPhiMuMu_Cuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "deltaPhiMuMu_Cuts", ";#Delta#phi( #mu_{1},#mu_{2})",64, -3.2, 3.2);
  hptRatioMuMu_Cuts =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "ptRatioMuMu_Cuts", ";pt_ratio",10,0.5,1.5);
  hMet_Cuts            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_Cuts", ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hMetPhi_Cuts         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "MetPhi_Cuts",";E_{T}^{miss} #phi (rads);Events", 64, -3.2, 3.2);
  hCorrMet_Cuts            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "CorrMet_Cuts", ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hCorrMetPhi_Cuts         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "CorrMetPhi_Cuts",";E_{T}^{miss} #phi (rads);Events", 64, -3.2, 3.2);
  hMuMuMass_Cuts = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "MuMuMass_Cuts", ";m_{#mu #mu} (GeV)", mN, mMin, mMax);

  // plots after medium taus
  hLepMedTauMass_AfterMetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauMass_AfterMetSelection"  , ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauPt_AfterMetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauPt_AfterMetSelection"  , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hLepMedTauCollMass_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMass_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPt_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPt_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_AfterMetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_AfterMetSelection", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  // dm=0; barrel and endcap
  hLepMedTauCollMassDiobjPtCollmassBinIncl_dm0_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_dm0_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinIncl_dm0_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_dm0_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_dm0_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_dm0_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_dm0_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_dm0_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_dm0_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_dm0_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_dm0_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_dm0_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  // dm=1; barrel and endcap
  hLepMedTauCollMassDiobjPtCollmassBinIncl_dm1_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_dm1_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinIncl_dm1_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_dm1_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_dm1_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_dm1_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_dm1_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_dm1_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_dm1_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_dm1_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_dm1_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_dm1_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  // dm=10; barrel and endcap
  hLepMedTauCollMassDiobjPtCollmassBinIncl_dm10_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_dm10_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinIncl_dm10_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinIncl_dm10_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_dm10_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_dm10_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinSS_dm10_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinSS_dm10_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_dm10_barrel = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_dm10_barrel", ";m_{l#tau} (GeV)", mN, mMin, mMax);
  hLepMedTauCollMassDiobjPtCollmassBinOS_dm10_endcap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepMedTauCollMassDiobjPtCollmassBinOS_dm10_endcap", ";m_{l#tau} (GeV)", mN, mMin, mMax);

  hLepMedTauPtVsLepMedTauCollMass_AfterMetSelection = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, dir, "LepMedTauPtVsLepMedTauCollMass_AfterMetSelection", "LepMedTauPtVsLepMedTauCollMass", ptN , ptMin , ptMax , mN, mMin, mMax);

  hTauSrc_AfterAllSelections       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrc_AfterAllSelections"    , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM0_AfterAllSelections    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM0_AfterAllSelections" , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM1_AfterAllSelections    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM1_AfterAllSelections" , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM10_AfterAllSelections   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM10_AfterAllSelections", ";#tau source", 80, 0.0, 80.0);

  hNTau_AfterAllSelections         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterAllSelections"        , ";loose #tau-jet multiplicity", nN  , nMin  , nMax  );
  hNJet_AfterAllSelections         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NJet_AfterAllSelections"        , ";jet multiplicity"  , nN  , nMin  , nMax  );
  hNBjet_AfterAllSelections        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NBjet_AfterAllSelections"       , ";b-jet multiplicity", nN  , nMin  , nMax  );
  hMet_AfterAllSelections          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_AfterAllSelections"         , ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hLepLooseTauPt_AfterAllSelections   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauPt_AfterAllSelections"  , ";p_{T} (GeV)" , ptN , ptMin , ptMax );
  hLepLooseTauEta_AfterAllSelections  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauEta_AfterAllSelections" , ";#eta"        , etaN, etaMin, etaMax);
  hLepLooseTauMass_AfterAllSelections = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LepLooseTauMass_AfterAllSelections", ";m_{l#tau} (GeV)", mN  , mMin  , mMax  );

  return;
}


void TauFakeRate_mt_MET::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void TauFakeRate_mt_MET::process(Long64_t entry) {

  //====== Initialize
  if (0) cout << "Welcome to muon+tau analyzer" << endl;
  fCommonPlots.initialize();
  //  fCommonPlots.setFactorisationBinForEvent(std::vector<float> {});
  cAllEvents.increment();

  //================================================================================================   
  // Apply trigger 
  //================================================================================================   
  if (0) std::cout << "=== Trigger" << std::endl;
  if ( !(fEvent.passTriggerDecision()) ) return;  
  cTrigger.increment();
  int nVertices = fEvent.vertexInfo().value(); // it returns "nGoodOfflineVertices" but for data "nPUvertices" is not filled
  fCommonPlots.setNvertices(nVertices);

  // Fill histos
  fCommonPlots.fillControlPlotsAfterTrigger(fEvent);

  //================================================================================================   
  // MET filters (to remove events with spurious sources of fake MET)
  //================================================================================================   
  if (0) std::cout << "=== MET Filter" << std::endl;
  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  //if (!metFilterData.passedSelection()) return; // 19Mar

  // Fill histos
  //fCommonPlots.fillControlPlotsAfterMETFilter(fEvent);  // 19Mar


  //================================================================================================   
  // Primarty Vertex (Check that a PV exists)
  //================================================================================================   
  if (0) std::cout << "=== Vertices" << std::endl;
  //  if (nVertices < 1) return; // 19Mar
  //cVertexSelection.increment(); // 19Mar

  // Fill histos
  //fCommonPlots.fillControlPlotsAtVertexSelection(fEvent); // 19Mar
  

  //================================================================================================   
  // Electron Selection
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  //if (eData.hasIdentifiedElectrons()) return; // 19Mar

  // Fill histos
  //fCommonPlots.fillControlPlotsAtElectronSelection(fEvent, eData); // 19Mar

  //================================================================================================
  // Muon Selection
  //================================================================================================
  if (0) std::cout << "=== Muon Selection" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if(!(muData.hasIdentifiedMuons())) return;

  // Require exactly 2 muons
  if (muData.getSelectedMuons().size() != 2) return; // note: remember to disable trigger-matching option if using a single muon trigger
  
  // Apply muon trigger scale factors
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/python/parameters/scaleFactors.py
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/data/TriggerEfficiency/muonPAGEff.json
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(muData.getMuonTriggerSF());
    }

  // Apply muon ID scale factors
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(muData.getMuonIDSF());
    }
  // control plots at muon selection
  fCommonPlots.fillControlPlotsAtMuonSelection(fEvent, muData);

  //================================================================================================
  // MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  // const METSelection::Data metData = fMETSelection.silentAnalyze(fEvent, nVertices);
  const METSelection::Data metData = fMETSelection.analyze(fEvent, nVertices);
  if (!metData.passedSelection()) return;
  //  fCommonPlots.fillControlPlotsAtMETSelection(fEvent, metData);
  
  hMet_AfterLepton->Fill(metData.getMET().R());
  hMetPhi_AfterLepton->Fill(metData.getMET().phi());
  // MET correction
  std::pair<double,double> abcd1 = METXYCorr_Met_MetPhi(metData.getMET().R(), metData.getMET().phi(), fEvent.eventID().run(), 2016, fEvent.isMC(), nVertices);

  hCorrMet_AfterLepton->Fill(abcd1.first);
  hCorrMetPhi_AfterLepton->Fill(abcd1.second);
  

  // recoomended cuts
  double mu1_pt = muData.getSelectedMuons()[0].p4().pt();
  double mu2_pt = muData.getSelectedMuons()[1].p4().pt();
  //double mu1_phi =  muData.getSelectedMuons()[0].p4().phi();
  //double mu2_phi =  muData.getSelectedMuons()[1].p4().phi();

  double deltaPhiMuMu = ROOT::Math::VectorUtil::DeltaPhi(muData.getSelectedMuons()[0].p4(),muData.getSelectedMuons()[1].p4());
  double ptRatioMuMu = mu1_pt/mu2_pt;
  
  if((0.8 < ptRatioMuMu && ptRatioMuMu < 1.2) && (fabs(deltaPhiMuMu) > 2.8 )){
    hdeltaPhiMuMu_Cuts->Fill(deltaPhiMuMu);
    hptRatioMuMu_Cuts->Fill(ptRatioMuMu);
    hMet_Cuts->Fill(metData.getMET().R());
    hMetPhi_Cuts->Fill(metData.getMET().phi());
    // MET correction
    std::pair<double,double> abcd2 = METXYCorr_Met_MetPhi(metData.getMET().R(), metData.getMET().phi(), fEvent.eventID().run(), 2016, fEvent.isMC(), nVertices);
    hCorrMet_Cuts->Fill(abcd2.first);
    hCorrMetPhi_Cuts->Fill(abcd2.second);
    double dimuon_invMass = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).M();
    hMuMuMass_Cuts->Fill(dimuon_invMass);
    if (!fEvent.isMC()) {
      hNSelectedVsRunNumber_Cuts->Fill(fEvent.eventID().run());
    }

  }
  //================================================================================================   
  // Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  const TauSelection::Data looseTauData = fLooseTauSelection.analyzeLoose(fEvent);

  if (!looseTauData.hasIdentifiedTaus() ) return;
  if(looseTauData.getSelectedTaus().size() != 1) return; // not the same as line above!

  if (0) std::cout << "=== Tau Selection: Has at least 1 loose tau(s)" << std::endl;
  cTauNCounter.increment(); 

  // Determine if genuine or fake tau. If more than one tau present use the first one in the list (pt-sorted)
  bool bEleToTau     = false; // ele  -> tau fakes only
  bool bMuonToTau    = false; // muon -> tau fakes only
  // bool bJetToTau     = false; // jet  -> tau fakes only
  bool bLightQToTau  = false; // light quark-> tau fakes only (u,d,s)
  bool bHeavyQToTau  = false; // heavy quark-> tau fakes only (c,b)
  bool bGluonToTau   = false; // jet  -> tau fakes only
  bool bUnknownToTau = false; // unknown -> tau fakes only
  bool bGenuineTau   = false;
  unsigned int tau_dm = looseTauData.getSelectedTau().decayMode();
  if ( fEvent.isMC() )
    {
      bEleToTau     = looseTauData.getSelectedTau().isElectronToTau();
      bMuonToTau    = looseTauData.getSelectedTau().isMuonToTau();
      // bJetToTau     = looseTauData.getSelectedTau().isJetToTau();
      bLightQToTau  = looseTauData.getSelectedTau().isQuarkToTau(1) || looseTauData.getSelectedTau().isQuarkToTau(2) || looseTauData.getSelectedTau().isQuarkToTau(3);
      bHeavyQToTau  = looseTauData.getSelectedTau().isQuarkToTau(4) || looseTauData.getSelectedTau().isQuarkToTau(5);
      bGluonToTau   = looseTauData.getSelectedTau().isGluonToTau();
      bUnknownToTau = looseTauData.getSelectedTau().isUnknownTauDecay();
      bGenuineTau   = looseTauData.getSelectedTau().isGenuineTau();
    }

  // Define a bit to store the source of fake. Note that more than one source might be true (rare but happens)
  // int tauSrcBit =  1*bEleToTau + 2*bMuonToTau + 3*bLightQToTau + 4*bHeavyQToTau + 5*bGluonToTau + 6*bUnknownToTau;// => genuineTau = 0
  int tauSrcBit =  pow(2,0)*bEleToTau + pow(2,1)*bMuonToTau + pow(2,2)*bLightQToTau + pow(2,3)*bHeavyQToTau + pow(2,4)*bGluonToTau + pow(2,5)*bUnknownToTau;// => genuineTau = 0

  // Fake tau in our analysis is "not genuine tau and not e->tau and not mu->tau"
  bool bFakeTau = ( fEvent.isMC() && !bGenuineTau && !(bEleToTau || bMuonToTau) );
  
  // Fill histos ( also sets value for boolean bIsGenuineTau
  //hLepLooseTauMass_AfterTauSelection->Fill(dilepton_invMass);
  hNTau_AfterTauSelection->Fill( looseTauData.getSelectedTaus().size() );  

  // Calculate variables for dilepton system
  double dilepton_invMass = ( muData.getSelectedMuons()[0].p4() + looseTauData.getSelectedTaus()[0].p4() ).M();
  double dilepton_pt      = ( muData.getSelectedMuons()[0].p4() + looseTauData.getSelectedTaus()[0].p4() ).pt();
  double dilepton_eta     = ( muData.getSelectedMuons()[0].p4() + looseTauData.getSelectedTaus()[0].p4() ).eta();
  
  hLepLooseTauMass_BeforeOnZSelection->Fill(dilepton_invMass);
  
  //Apply dR cut for dilepton system?   
  double lep_tau_dR = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),looseTauData.getSelectedTaus()[0].p4());
  if (lep_tau_dR < 0.3) return;
  
  hLepLooseTauMass_AfterDRSelection->Fill(dilepton_invMass);

  // Apply dR cut for dilepton system?
  //  dilepton_dR = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),muData.getSelectedMuons()[1].p4());
  //  if(dilepton_dR < 0.3) return;

  // Fill histos
  hLepLooseTauPt_AfterLeptonSelection->Fill(dilepton_pt);
  hLepLooseTauEta_AfterLeptonSelection->Fill(dilepton_eta);
  hLepLooseTauMass_AfterLeptonSelection->Fill(dilepton_invMass);
  hLeptonN_AfterLeptonSelection ->Fill(muData.getSelectedMuons().size());
  /*
  // For-loop: All Selected leptons
  for(unsigned int i=0; i< muData.getSelectedMuons().size(); i++)
    {
      hLeptonPt_AfterLeptonSelection ->Fill(muData.getSelectedMuons()[i].pt());
      hLeptonEta_AfterLeptonSelection->Fill(muData.getSelectedMuons()[i].eta());
    }
  */

  // Defining the splitting of phase-space for control-plots
  std::vector<float> myFactorisationInfo;
  myFactorisationInfo.push_back( looseTauData.getSelectedTaus()[0].pt() );
  myFactorisationInfo.push_back( looseTauData.getSelectedTaus()[0].eta() );
  fCommonPlots.setFactorisationBinForEvent(myFactorisationInfo);
  
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, looseTauData); // uses: bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();

  // Redefine what the "bGenuineTau" boolean is. N.B: Affects the genuineTau plots filled by control plots!
  if (0) std::cout << "=== Tau Selection: bFakeTau = " << bFakeTau << " (src = " << tauSrcBit << "). NTaus = " << looseTauData.getSelectedTaus().size() << std::endl;
  fCommonPlots.setGenuineTauStatus(!bFakeTau); // CommonPlots.cc (L1292) bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();

  
  // Fill Histos

  hLepLooseTauMass_AfterMetSelection->Fill(dilepton_invMass);
  hNTau_AfterMetSelection->Fill( looseTauData.getSelectedTaus().size() );

  //dilepton_dR = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),muData.getSelectedMuons()[1].p4());   
  // projection MET vector into lepton direction
  // ref: https://root.cern/doc/master/namespaceROOT_1_1Math_1_1VectorUtil.html#ad067df491cad5594499d9a1860a00f4c
  //math::XYZTLorentzVector metOntoMuon; (not compile)
  //TLorentzVector metOntoMuon;
  //metOntoMuon = ROOT::Math::VectorUtil::ProjVector(metData.getMET(),muData.getSelectedMuons()[0].p4());
  //metOntoMuon = ROOT::Math::VectorUtil::ProjVector(looseTauData.getSelectedTaus()[0].p4(),muData.getSelectedMuons()[0].p4()); 
  //double DR = ROOT::Math::VectorUtil::DeltaPhi(metData.getMET(),muData.getSelectedMuons()[0].p4()); // DeltaPhi compile but not DeltaR with metData
  //double costheta = ROOT::Math::VectorUtil::CosTheta(looseTauData.getSelectedTaus()[0].p4(),muData.getSelectedMuons()[0].p4());
  //double costheta = ROOT::Math::VectorUtil::CosTheta(metData.getMET(),muData.getSelectedMuons()[0].p4()); // does not compile as metData.getMET() have only x, y componetn  
  //cout << "MET X and Y component: " << metData.getMET().X() << " & " << metData.getMET().Y() << endl;
  //cout << "muData.getSelectedMuons()[0].px(): " << muData.getSelectedMuons()[0].p4().X() << endl; //compiled
  double dPhi_met_tau = ROOT::Math::VectorUtil::DeltaPhi(metData.getMET(),looseTauData.getSelectedTaus()[0].p4());
  double dPhi_met_lep = ROOT::Math::VectorUtil::DeltaPhi(metData.getMET(),muData.getSelectedMuons()[0].p4());
  double dPhi_tau_lep = ROOT::Math::VectorUtil::DeltaPhi(looseTauData.getSelectedTaus()[0].p4(),muData.getSelectedMuons()[0].p4());
  bool collinerFlag = false;
  if ((dPhi_met_tau+dPhi_met_lep) < dPhi_tau_lep){
    //cout << "colliner MET found" << endl;
    collinerFlag = true;
  }

  math::XYZTLorentzVector newMuon;
  if(collinerFlag){
    newMuon.SetPxPyPzE(muData.getSelectedMuons()[0].p4().X()+metData.getMET().X(),muData.getSelectedMuons()[0].p4().Y()+metData.getMET().Y(),muData.getSelectedMuons()[0].p4().Z(),muData.getSelectedMuons()[0].p4().E());
  }else{
    newMuon.SetPxPyPzE(muData.getSelectedMuons()[0].p4().X(),muData.getSelectedMuons()[0].p4().Y(),muData.getSelectedMuons()[0].p4().Z(),muData.getSelectedMuons()[0].p4().E());
  }
  math::XYZTLorentzVector newTau;
  if(collinerFlag){
    newTau.SetPxPyPzE(looseTauData.getSelectedTaus()[0].p4().X()+metData.getMET().X(),looseTauData.getSelectedTaus()[0].p4().Y()+metData.getMET().Y(),looseTauData.getSelectedTaus()[0].p4().Z(),looseTauData.getSelectedTaus()[0].p4().E());
  }else{
    newTau.SetPxPyPzE(looseTauData.getSelectedTaus()[0].p4().X(),looseTauData.getSelectedTaus()[0].p4().Y(),looseTauData.getSelectedTaus()[0].p4().Z(),looseTauData.getSelectedTaus()[0].p4().E());
  }
  // new invariant mass 
  double new_invMass = (newMuon + newTau).M();
  if (0) cout << "new_invMass: " << new_invMass << endl;
  if(collinerFlag) hLepLooseTauMass_METProjZSelection->Fill(new_invMass);

  // collinear mass
  // CollinearMass( const MuonSelection::Data& muData, const TauSelection::Data& looseTauData, const METSelection::Data& metData);
  double collmass = CollinearMass(muData, looseTauData, metData);
  if (0) cout << "collmass: " << collmass << endl;
  hLepLooseTauCollMass_AfterMetSelection->Fill(collmass);

  // collinear mass
  if (collinerFlag) {
    double collmass_Flag = CollinearMass(muData, looseTauData, metData);
    hLepLooseTauCollMassFlag_AfterMetSelection->Fill(collmass_Flag);
  }

  // 40 < dilepton_invMass < 80
  if (dilepton_invMass < 80 && dilepton_invMass > 40){
    double collmass_mass_40_80 = CollinearMass(muData, looseTauData, metData);
    hLepLooseTauMass40_80_AfterMetSelection->Fill(dilepton_invMass);
    hLepLooseTauCollMass40_80_AfterMetSelection->Fill(collmass_mass_40_80);
  }

  if (collinerFlag) {
    if (dilepton_invMass < 80 && dilepton_invMass > 40){
      double collmass_Flag_mass_40_80 = CollinearMass(muData, looseTauData, metData);
      hLepLooseTauCollMassFlag_mass_40_80_AfterMetSelection->Fill(collmass_Flag_mass_40_80);
    }
  }
  
  hLepLooseTauPtVsLepLooseTauMass_AfterMetSelection->Fill(dilepton_pt,dilepton_invMass);
  hMet_AfterMetSelection->Fill(metData.getMET().R());
  hMetPhi_AfterMetSelection->Fill(metData.getMET().phi());

  // cout << "METValue " << metData.getMET().R() << endl;
  // cout << "METPhi " << metData.getMET().phi() << endl;
  
  // corrected MET value and MET PHI
  // to get runnumber: event.eventID().run() valid for !fEvent.isMC()
  // Function arguments (double uncormet,double uncormet_phi, int runnb, int year, bool isMC, int npv)
  // int runnb = fEvent.eventID().run();
  std::pair<double,double> abcd = METXYCorr_Met_MetPhi(metData.getMET().R(), metData.getMET().phi(), fEvent.eventID().run(), 2016, fEvent.isMC(), nVertices);

  hCorrMet_AfterMetSelection->Fill(abcd.first);
  hCorrMetPhi_AfterMetSelection->Fill(abcd.second);
  //cout << "corrected METValue " << abcd.first << endl;
  //cout << "corrected METPhi " << abcd.second << endl;

  hLepLooseTauPt_AfterMetSelection->Fill(dilepton_pt);
  hDPhiTauMet_AfterMetSelection->Fill(dPhi_met_tau);
  hDPhiTauLep_AfterMetSelection->Fill(dPhi_tau_lep);
  hDPhiMetLep_AfterMetSelection->Fill(dPhi_met_lep);

  // check of if the tau is real or not
  if (0) cout << "Tau flag: " << looseTauData.getSelectedTaus()[0].isGenuineTau() << endl;
  // fill diobjectmass and collmass for GenuineTau
  if ( fEvent.isMC() ){
    if (looseTauData.getSelectedTaus()[0].isGenuineTau()){
      hLepLooseTauMassGenuTau_AfterMetSelection->Fill(dilepton_invMass);
      double collmass_GenuTau = CollinearMass(muData, looseTauData, metData);
      hLepLooseTauCollMassGenuTau_AfterMetSelection->Fill(collmass_GenuTau);
    }
  }
  

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2                                                                                                                       
  // DM = 5*(Nc-1)+Np  [Nc = # charged hadrons, Np = # of pi0s]. DN = 0-11
  // unsigned int tau_dm = tauData.getSelectedTaus()[i].decayMode();

  // cut on dilepton pt
  if(dilepton_pt > 40){
    double collmass_dlpt = CollinearMass(muData, looseTauData, metData);
    hLepLooseTauCollMassDiobjPt_AfterMetSelection->Fill(collmass_dlpt);
    // decay mode 
    unsigned int tau_dm = looseTauData.getSelectedTaus()[0].decayMode();
    double tau_eta      = looseTauData.getSelectedTaus()[0].eta();
    bool tau_inBarrel   = ( fabs(tau_eta) < 1.5);

    hLepLooseTauPtVsLepLooseTauCollMass_AfterMetSelection->Fill(dilepton_pt,collmass_dlpt);
    if (collmass_dlpt > 65 && collmass_dlpt < 115){
      hLepLooseTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection->Fill(collmass_dlpt);
      // make mass for SS and OS
      if (muData.getSelectedMuons()[0].charge() == looseTauData.getSelectedTaus()[0].charge()){
	hLepLooseTauCollMassDiobjPtCollmassBinSS_AfterMetSelection->Fill(collmass_dlpt);
      }else{
	hLepLooseTauCollMassDiobjPtCollmassBinOS_AfterMetSelection->Fill(collmass_dlpt);
      }
      // decay mode (DM) = 0  
      if (tau_dm == 0){
	if(tau_inBarrel) hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_barrel->Fill(collmass_dlpt);
	else hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm0_endcap->Fill(collmass_dlpt);
	// make mass for SS and OS
	if (muData.getSelectedMuons()[0].charge() == looseTauData.getSelectedTaus()[0].charge()){
	  if(tau_inBarrel)  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm0_barrel->Fill(collmass_dlpt);
	  else hLepLooseTauCollMassDiobjPtCollmassBinSS_dm0_endcap->Fill(collmass_dlpt);
	}else{
	  if(tau_inBarrel) hLepLooseTauCollMassDiobjPtCollmassBinOS_dm0_barrel->Fill(collmass_dlpt);
	  else hLepLooseTauCollMassDiobjPtCollmassBinOS_dm0_endcap->Fill(collmass_dlpt);
	}
      } //dm = 0
      // decay mode (DM) = 1
      if (tau_dm == 1){
	if(tau_inBarrel) hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_barrel->Fill(collmass_dlpt);
	else hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm1_endcap->Fill(collmass_dlpt);
	// make mass for SS and OS
	if (muData.getSelectedMuons()[0].charge() == looseTauData.getSelectedTaus()[0].charge()){
	  if(tau_inBarrel)  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm1_barrel->Fill(collmass_dlpt);
	  else hLepLooseTauCollMassDiobjPtCollmassBinSS_dm1_endcap->Fill(collmass_dlpt);
	}else{
	  if(tau_inBarrel) hLepLooseTauCollMassDiobjPtCollmassBinOS_dm1_barrel->Fill(collmass_dlpt);
	  else hLepLooseTauCollMassDiobjPtCollmassBinOS_dm1_endcap->Fill(collmass_dlpt);
	}
      } //dm = 1
      // decay mode (DM) = 10
      if (tau_dm == 10){
	if(tau_inBarrel) hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_barrel->Fill(collmass_dlpt);
	else hLepLooseTauCollMassDiobjPtCollmassBinIncl_dm10_endcap->Fill(collmass_dlpt);
	// make mass for SS and OS
	if (muData.getSelectedMuons()[0].charge() == looseTauData.getSelectedTaus()[0].charge()){
	  if(tau_inBarrel)  hLepLooseTauCollMassDiobjPtCollmassBinSS_dm10_barrel->Fill(collmass_dlpt);
	  else hLepLooseTauCollMassDiobjPtCollmassBinSS_dm10_endcap->Fill(collmass_dlpt);
	}else{
	  if(tau_inBarrel) hLepLooseTauCollMassDiobjPtCollmassBinOS_dm10_barrel->Fill(collmass_dlpt);
	  else hLepLooseTauCollMassDiobjPtCollmassBinOS_dm10_endcap->Fill(collmass_dlpt);
	}
      } //dm = 10
    } // mass range
  } // dilepton pt

  // save the tau pt as well here for FR application
  
  // do the same mass plot with "medium" Tau
  if (!tauData.hasIdentifiedTaus() ) return;
  if(tauData.getSelectedTaus().size() != 1) return; // not the same as line above! 
  
  // Calculate variables for dilepton system
  double lepMedTau_invMass = ( muData.getSelectedMuons()[0].p4() + tauData.getSelectedTaus()[0].p4() ).M();
  double lepMedTau_pt      = ( muData.getSelectedMuons()[0].p4() + tauData.getSelectedTaus()[0].p4() ).pt();
  hLepMedTauMass_AfterMetSelection->Fill(lepMedTau_invMass);
  hLepMedTauPt_AfterMetSelection->Fill(lepMedTau_pt);
  // collmass
  double collmass_lepMedTau = CollinearMass(muData, tauData, metData);
  hLepMedTauCollMass_AfterMetSelection->Fill(collmass_lepMedTau);
  
  // cut on dilepton pt
  if(lepMedTau_pt > 40){
    double collmass_lepMedTau_pt = CollinearMass(muData, tauData, metData);
    hLepMedTauCollMassDiobjPt_AfterMetSelection->Fill(collmass_lepMedTau_pt);
    // decay mode 
    unsigned int tight_tau_dm = tauData.getSelectedTaus()[0].decayMode();
    double tight_tau_eta      = tauData.getSelectedTaus()[0].eta();
    bool tight_tau_inBarrel   = ( fabs(tight_tau_eta) < 1.5);

    hLepMedTauPtVsLepMedTauCollMass_AfterMetSelection->Fill(lepMedTau_pt,collmass_lepMedTau_pt);
    if (collmass_lepMedTau_pt > 65 && collmass_lepMedTau_pt < 115){
      hLepMedTauCollMassDiobjPtCollmassBinIncl_AfterMetSelection->Fill(collmass_lepMedTau_pt);
      // make mass for SS and OS
      if (muData.getSelectedMuons()[0].charge() == tauData.getSelectedTaus()[0].charge()){
	hLepMedTauCollMassDiobjPtCollmassBinSS_AfterMetSelection->Fill(collmass_lepMedTau_pt);
      }else{
	hLepMedTauCollMassDiobjPtCollmassBinOS_AfterMetSelection->Fill(collmass_lepMedTau_pt);
      }
      // decay mode (DM) = 0  
      if (tight_tau_dm == 0){
	if(tight_tau_inBarrel) hLepMedTauCollMassDiobjPtCollmassBinIncl_dm0_barrel->Fill(collmass_lepMedTau_pt);
	else hLepMedTauCollMassDiobjPtCollmassBinIncl_dm0_endcap->Fill(collmass_lepMedTau_pt);
	// make mass for SS and OS
	if (muData.getSelectedMuons()[0].charge() == tauData.getSelectedTaus()[0].charge()){
	  if(tight_tau_inBarrel)  hLepMedTauCollMassDiobjPtCollmassBinSS_dm0_barrel->Fill(collmass_lepMedTau_pt);
	  else hLepMedTauCollMassDiobjPtCollmassBinSS_dm0_endcap->Fill(collmass_lepMedTau_pt);
	}else{
	  if(tight_tau_inBarrel) hLepMedTauCollMassDiobjPtCollmassBinOS_dm0_barrel->Fill(collmass_lepMedTau_pt);
	  else hLepMedTauCollMassDiobjPtCollmassBinOS_dm0_endcap->Fill(collmass_lepMedTau_pt);
	}
      } //dm = 0
      // decay mode (DM) = 1
      if (tight_tau_dm == 1){
	if(tight_tau_inBarrel) hLepMedTauCollMassDiobjPtCollmassBinIncl_dm1_barrel->Fill(collmass_lepMedTau_pt);
	else hLepMedTauCollMassDiobjPtCollmassBinIncl_dm1_endcap->Fill(collmass_lepMedTau_pt);
	// make mass for SS and OS
	if (muData.getSelectedMuons()[0].charge() == tauData.getSelectedTaus()[0].charge()){
	  if(tight_tau_inBarrel)  hLepMedTauCollMassDiobjPtCollmassBinSS_dm1_barrel->Fill(collmass_lepMedTau_pt);
	  else hLepMedTauCollMassDiobjPtCollmassBinSS_dm1_endcap->Fill(collmass_lepMedTau_pt);
	}else{
	  if(tight_tau_inBarrel) hLepMedTauCollMassDiobjPtCollmassBinOS_dm1_barrel->Fill(collmass_lepMedTau_pt);
	  else hLepMedTauCollMassDiobjPtCollmassBinOS_dm1_endcap->Fill(collmass_lepMedTau_pt);
	}
      } //dm = 1
      // decay mode (DM) = 10
      if (tight_tau_dm == 10){
	if(tight_tau_inBarrel) hLepMedTauCollMassDiobjPtCollmassBinIncl_dm10_barrel->Fill(collmass_lepMedTau_pt);
	else hLepMedTauCollMassDiobjPtCollmassBinIncl_dm10_endcap->Fill(collmass_lepMedTau_pt);
	// make mass for SS and OS
	if (muData.getSelectedMuons()[0].charge() == tauData.getSelectedTaus()[0].charge()){
	  if(tight_tau_inBarrel)  hLepMedTauCollMassDiobjPtCollmassBinSS_dm10_barrel->Fill(collmass_lepMedTau_pt);
	  else hLepMedTauCollMassDiobjPtCollmassBinSS_dm10_endcap->Fill(collmass_lepMedTau_pt);
	}else{
	  if(tight_tau_inBarrel) hLepMedTauCollMassDiobjPtCollmassBinOS_dm10_barrel->Fill(collmass_lepMedTau_pt);
	  else hLepMedTauCollMassDiobjPtCollmassBinOS_dm10_endcap->Fill(collmass_lepMedTau_pt);
	}
      } //dm = 1

    } // mass range
  } // dilepton pt
  
  // ekhon etotai gkole

  //================================================================================================
  // Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, looseTauData.getSelectedTau());
  if (!jetData.passedSelection()) return;  
    
  // Fill histos
  fCommonPlots.fillControlPlotsAtJetSelection(fEvent, jetData);
  hLepLooseTauMass_AfterJetSelection->Fill(dilepton_invMass);
  hNTau_AfterJetSelection->Fill( looseTauData.getSelectedTaus().size() );

 
  //================================================================================================  
  // BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;
  
  // Fill histos
  fCommonPlots.fillControlPlotsAtBtagging(fEvent, bjetData);

  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAfterBtagSF(fEvent, jetData, bjetData);
  hLepLooseTauMass_AfterBJetSelection->Fill(dilepton_invMass);
  hNTau_AfterBJetSelection->Fill( looseTauData.getSelectedTaus().size() );


  //================================================================================================
  // All Selections
  //================================================================================================
  if (0) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();
  
  // Fill histos
  //fCommonPlots.fillControlPlotsAfterAllSelections(fEvent); // it was some issues as "fillControlPlotsAfterAllSelections" isnot there (fix me? gkole)
  fCommonPlots.fillControlPlotsAfterPreselections(fEvent, jetData, bjetData, metData);
  
  if (fEvent.isMC() )
    {
      hTauSrc_AfterAllSelections->Fill(tauSrcBit);
      if (tau_dm == 0)  hTauSrcDM0_AfterAllSelections->Fill(tauSrcBit);
      if (tau_dm == 1)  hTauSrcDM1_AfterAllSelections->Fill(tauSrcBit);
      if (tau_dm == 10) hTauSrcDM10_AfterAllSelections->Fill(tauSrcBit);
    }
  hNTau_AfterAllSelections->Fill( looseTauData.getSelectedTaus().size() );
  hNJet_AfterAllSelections ->Fill( jetData.getSelectedJets().size() );
  hNBjet_AfterAllSelections->Fill( bjetData.getSelectedBJets().size() );
  hLepLooseTauPt_AfterAllSelections->Fill(dilepton_pt);
  hLepLooseTauEta_AfterAllSelections->Fill(dilepton_eta);
  hLepLooseTauMass_AfterAllSelections->Fill(dilepton_invMass);
  hMet_AfterAllSelections->Fill(metData.getMET().R());

  // Tau stuff here
  if(looseTauData.getSelectedTaus().size() != 1) return; // why is this necessary? (remove it and get crash)

  // "Tight" Tau
  if (tauData.hasIdentifiedTaus()) 
    {

      if (fEvent.isMC()) 
	{
	  // Apply "tight" tau ID scale factor (SF)
	  fEventWeight.multiplyWeight(tauData.getTauIDSF());
	  cTauSFCounter.increment(); 
	  
	  // Apply "tight" fake tau SF
	  fEventWeight.multiplyWeight(tauData.getTauMisIDSF());
	  cFakeTauSFCounter.increment();
	}

      // Do rest of event selection
      doTightTaus(fEvent, tauData);
      if (fEvent.isMC())
	{ 
	  // Undo tau ID SF!
	  fEventWeight.multiplyWeight(1.0/tauData.getTauIDSF());

	  // Undo tau mis-ID SF!
	  fEventWeight.multiplyWeight(1.0/tauData.getTauMisIDSF());
	}
    }

  // "Loose" Tau
  if (looseTauData.hasIdentifiedTaus()) 
    {
      
      if (fEvent.isMC())
	{
	  // Apply "loose" tau ID SF
	  fEventWeight.multiplyWeight(looseTauData.getTauIDSF());
	  // Apply "loose" tau mis-ID SF
	  fEventWeight.multiplyWeight(looseTauData.getTauMisIDSF());
	}

      // Do rest of event selection
      doLooseTaus(fEvent, looseTauData);
    }

  fEventSaver.save();

  return;
}

double TauFakeRate_mt_MET::CollinearMass( const MuonSelection::Data& muData, const TauSelection::Data& looseTauData, const METSelection::Data& metData) {
  //  TLorentzVector k1 = muData.getSelectedMuons()[0].p4(); // p1->p4();
  TLorentzVector k1(muData.getSelectedMuons()[0].p4().X(),muData.getSelectedMuons()[0].p4().Y(),muData.getSelectedMuons()[0].p4().Z(),muData.getSelectedMuons()[0].p4().E());
  //TLorentzVector k2 = looseTauData.getSelectedTaus()[0].p4(); //p2->p4();
  TLorentzVector k2(looseTauData.getSelectedTaus()[0].p4().X(),looseTauData.getSelectedTaus()[0].p4().Y(),looseTauData.getSelectedTaus()[0].p4().Z(),looseTauData.getSelectedTaus()[0].p4().E());

  TMatrixD K(2, 2);
  K(0, 0) = k1.X();      K(0, 1) = k2.X();
  K(1, 0) = k1.Y();      K(1, 1) = k2.Y();

  if(K.Determinant()==0)
    return -1234.;

  TMatrixD M(2, 1);
  M(0, 0) = metData.getMET().X(); //met->mpx();
  M(1, 0) = metData.getMET().Y(); //met->mpy();

  TMatrixD Kinv = K.Invert();

  TMatrixD X(2, 1);
  X = Kinv * M;

  double X1 = X(0, 0);        double X2 = X(1, 0);
  double x1 = 1. / (1. + X1); double x2 = 1. / (1. + X2);

  TLorentzVector part1 = k1 * (1 / x1);
  TLorentzVector part2 = k2 * (1 / x2);

  double m = (part1 + part2).M();
  return m;

}

void TauFakeRate_mt_MET::doLooseTaus(const Event& event, const TauSelection::Data& tauData) {


  // For-loop: All Selected taus
  for(unsigned int i=0; i<tauData.getSelectedTaus().size(); i++)
    {      
      
      // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
      // DM = 5*(Nc-1)+Np  [Nc = # chargde hadrons, Np = # of pi0s]. DN = 0-11 
      unsigned int tau_dm = tauData.getSelectedTaus()[i].decayMode();
      double tau_pt       = tauData.getSelectedTaus()[i].pt();
      double tau_eta      = tauData.getSelectedTaus()[i].eta();
      bool tau_inBarrel   = ( fabs(tau_eta) < 1.5);
      bool bIsMC          = fEvent.isMC();
      bool bGenuineTau    = false;
      bool bEleToTau      = false;
      bool bMuonToTau     = false;
      bool bFakeTau       = false;
      if (bIsMC)
	{
	  bGenuineTau = tauData.getSelectedTaus()[i].isGenuineTau(); // tauData.isGenuineTau();
	  bEleToTau   = tauData.getSelectedTaus()[i].isElectronToTau();
	  bMuonToTau  = tauData.getSelectedTaus()[i].isMuonToTau();
	  bFakeTau    = (bIsMC && !bGenuineTau && !(bEleToTau || bMuonToTau) );
	}

      // Only if MC  and selected tau is genuine (not fake)
      // if (event.isMC() && tauData.getSelectedTaus()[i].isGenuineTau()) 
      if (bIsMC && !bFakeTau)
	{
	  if (tau_dm == 0)  
	    {
	      hTauPt_den_g_dm0->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_den_g_dm0_barrel->Fill( tau_pt );
	      else hTauPt_den_g_dm0_endcap->Fill( tau_pt );
	    }
	  
	  if (tau_dm==1)
	    {
	      hTauPt_den_g_dm1->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_den_g_dm1_barrel->Fill( tau_pt );
	      else hTauPt_den_g_dm1_endcap->Fill( tau_pt );
	    }
	  
	  if (tau_dm==10) 
	    {
	      hTauPt_den_g_dm10->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_den_g_dm10_barrel->Fill( tau_pt );
	      else hTauPt_den_g_dm10_endcap->Fill( tau_pt );
	    }
	}// (!bFakeTau)
    
      // 1-prong decays; decay mode (DM) = 0
      if (tau_dm==0) 
	{
	  hTauPt_den_dm0 ->Fill( tau_pt  );
	  hTauEta_den_dm0->Fill( tau_eta );

	  if (tau_inBarrel)
	    {
	      hTauPt_den_dm0_barrel ->Fill( tau_pt  );
	      hTauEta_den_dm0_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_den_dm0_endcap ->Fill( tau_pt  );
	      hTauEta_den_dm0_endcap->Fill( tau_eta );
	    }
	}
      
      // 2-prong decays; decay mode (DM) = 1
      if (tau_dm==1) 
	{
	  hTauPt_den_dm1 ->Fill( tau_pt  );
	  hTauEta_den_dm1->Fill( tau_eta );
	  
	  if (tau_inBarrel)
	    {
	      hTauPt_den_dm1_barrel ->Fill( tau_pt  );
	      hTauEta_den_dm1_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_den_dm1_endcap ->Fill( tau_pt  );
	      hTauEta_den_dm1_endcap->Fill( tau_eta );
	    }
	}
      
      // 3-prong decays; decay mode (DM) = 10
      if (tau_dm==10) 
	{
	  hTauPt_den_dm10 ->Fill( tau_pt  );
	  hTauEta_den_dm10->Fill( tau_eta );
	  
	  if (tau_inBarrel)
	    {
	      hTauPt_den_dm10_barrel ->Fill( tau_pt  );
	      hTauEta_den_dm10_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_den_dm10_endcap ->Fill( tau_pt  );
	      hTauEta_den_dm10_endcap->Fill( tau_eta );
	    }
	}
    }
  
  return;
}


void TauFakeRate_mt_MET::doTightTaus(const Event& event, const TauSelection::Data& tightTauData) {

  // For-loop: All selected taus
  for(unsigned int i=0; i<tightTauData.getSelectedTaus().size(); i++)
    {
      
      // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
      // DM = 5*(Nc-1)+Np  [Nc = # charged hadrons, Np = # of pi0s]. DN = 0-11 
      unsigned int tau_dm = tightTauData.getSelectedTaus()[i].decayMode();
      double tau_pt       = tightTauData.getSelectedTaus()[i].pt();
      double tau_eta      = tightTauData.getSelectedTaus()[i].eta();
      bool tau_inBarrel   = ( fabs(tau_eta) < 1.5);
      bool bIsMC          = fEvent.isMC();
      bool bGenuineTau    = false;
      bool bEleToTau      = false;
      bool bMuonToTau     = false;
      bool bFakeTau       = false;
      if (bIsMC)
	{
	  bGenuineTau = tightTauData.getSelectedTaus()[i].isGenuineTau();
	  bEleToTau   = tightTauData.getSelectedTaus()[i].isElectronToTau();
	  bMuonToTau  = tightTauData.getSelectedTaus()[i].isMuonToTau();
	  bFakeTau    = (!bGenuineTau && !(bEleToTau || bMuonToTau) );
	}

      // Only if MC  and selected tau is genuine. (if not genuine tau, reject the events)
      //if (event.isMC() && tightTauData.getSelectedTaus()[i].isGenuineTau()) 
      if (bIsMC && !bFakeTau)
	{
	  if (tau_dm==0) 
	    {
	      hTauPt_num_g_dm0->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_num_g_dm0_barrel->Fill( tau_pt );
	      else hTauPt_num_g_dm0_endcap->Fill( tau_pt );
	    }

	  if (tau_dm==1)
	    {
	      hTauPt_num_g_dm1->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_num_g_dm1_barrel->Fill( tau_pt );
	      else hTauPt_num_g_dm1_endcap->Fill( tau_pt );	   
	    }

	  if (tau_dm==10) 
	    {
	      hTauPt_num_g_dm10->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_num_g_dm10_barrel->Fill( tau_pt );
	      else hTauPt_num_g_dm10_endcap->Fill( tau_pt );
	    }
	}
      
      // 1-prong decays; decay mode (DM) = 0
      if (tau_dm==0) 
	{
	  hTauPt_num_dm0 ->Fill( tau_pt  );
	  hTauEta_num_dm0->Fill( tau_eta );

	  if (tau_inBarrel) 
	    {
	      hTauPt_num_dm0_barrel ->Fill( tau_pt  );
	      hTauEta_num_dm0_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_num_dm0_endcap ->Fill( tau_pt  );
	      hTauEta_num_dm0_endcap->Fill( tau_eta );
	    }
	}
      
      // 2-prong decays; decay mode (DM) = 1
      if (tau_dm==1) 
	{
	  hTauPt_num_dm1 ->Fill( tau_pt  );
	  hTauEta_num_dm1->Fill( tau_eta );
	  if (tau_inBarrel) 
	    {
	      hTauPt_num_dm1_barrel ->Fill( tau_pt  );
	      hTauEta_num_dm1_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_num_dm1_endcap ->Fill( tau_pt  );
	      hTauEta_num_dm1_endcap->Fill( tau_eta );
	    }
	}
      
      // 3-prong decays; decay mode (DM) = 10
      if (tau_dm==10) 
	{
	  hTauPt_num_dm10 ->Fill( tau_pt  );
	  hTauEta_num_dm10->Fill( tau_eta );
	  if (tau_inBarrel) 
	    {
	      hTauPt_num_dm10_barrel ->Fill( tau_pt  );
	      hTauEta_num_dm10_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_num_dm10_endcap ->Fill( tau_pt  );
	      hTauEta_num_dm10_endcap->Fill( tau_eta );
	    }
	}
    }

  return;
}
//**************************************************************************** 
enum TheRunEra{y2016B,y2016C,y2016D,y2016E,y2016F,y2016G,y2016H,y2017B,y2017C,y2017D,y2017E,y2017F,y2018A,y2018B,y2018C,y2018D,y2016MC,y2017MC,y2018MC};
std::pair<double,double> TauFakeRate_mt_MET::METXYCorr_Met_MetPhi(double uncormet,double uncormet_phi, int runnb, int year, bool isMC, int npv){
  std::pair<double,double>  TheXYCorr_Met_MetPhi(uncormet,uncormet_phi);
  if(npv>100) npv=100;
  int runera =-1;
  bool usemetv2 =false;
  if(isMC && year == 2016) runera = y2016MC;
  else if(isMC && year == 2017) {runera = y2017MC; usemetv2 =true;}
  else if(isMC && year == 2018) runera = y2018MC;
  
  else if(!isMC && runnb >=272007 &&runnb<=275376  ) runera = y2016B;
  else if(!isMC && runnb >=275657 &&runnb<=276283  ) runera = y2016C;
  else if(!isMC && runnb >=276315 &&runnb<=276811  ) runera = y2016D;
  else if(!isMC && runnb >=276831 &&runnb<=277420  ) runera = y2016E;
  else if(!isMC && runnb >=277772 &&runnb<=278808  ) runera = y2016F;
  else if(!isMC && runnb >=278820 &&runnb<=280385  ) runera = y2016G;
  else if(!isMC && runnb >=280919 &&runnb<=284044  ) runera = y2016H;
  
  else if(!isMC && runnb >=297020 &&runnb<=299329 ){ runera = y2017B; usemetv2 =true;}
  else if(!isMC && runnb >=299337 &&runnb<=302029 ){ runera = y2017C; usemetv2 =true;}
  else if(!isMC && runnb >=302030 &&runnb<=303434 ){ runera = y2017D; usemetv2 =true;}
  else if(!isMC && runnb >=303435 &&runnb<=304826 ){ runera = y2017E; usemetv2 =true;}
  else if(!isMC && runnb >=304911 &&runnb<=306462 ){ runera = y2017F; usemetv2 =true;}
  
  else if(!isMC && runnb >=315252 &&runnb<=316995 ) runera = y2018A;
  else if(!isMC && runnb >=316998 &&runnb<=319312 ) runera = y2018B;
  else if(!isMC && runnb >=319313 &&runnb<=320393 ) runera = y2018C;
  else if(!isMC && runnb >=320394 &&runnb<=325273 ) runera = y2018D;

  else {
    //Couldn't find data/MC era => no correction applied
    return TheXYCorr_Met_MetPhi;
  }
  //cout << "isMC: " << isMC << endl;
  //cout << "runera: " << runera << endl;
  double METxcorr(0.),METycorr(0.);
  
  if(!usemetv2){//Current recommendation for 2016 and 2018
    //if(runera==y2016B) cout << "in data2016B in !usemetv2 " << endl;
    if(runera==y2016B) METxcorr = -(-0.0478335*npv -0.108032);
    if(runera==y2016B) METycorr = -(0.125148*npv +0.355672);
    if(runera==y2016C) METxcorr = -(-0.0916985*npv +0.393247);
    if(runera==y2016C) METycorr = -(0.151445*npv +0.114491);
    if(runera==y2016D) METxcorr = -(-0.0581169*npv +0.567316);
    if(runera==y2016D) METycorr = -(0.147549*npv +0.403088);
    if(runera==y2016E) METxcorr = -(-0.065622*npv +0.536856);
    if(runera==y2016E) METycorr = -(0.188532*npv +0.495346);
    if(runera==y2016F) METxcorr = -(-0.0313322*npv +0.39866);
    if(runera==y2016F) METycorr = -(0.16081*npv +0.960177);
    if(runera==y2016G) METxcorr = -(0.040803*npv -0.290384);
    if(runera==y2016G) METycorr = -(0.0961935*npv +0.666096);
    if(runera==y2016H) METxcorr = -(0.0330868*npv -0.209534);
    if(runera==y2016H) METycorr = -(0.141513*npv +0.816732);
    if(runera==y2017B) METxcorr = -(-0.259456*npv +1.95372);
    if(runera==y2017B) METycorr = -(0.353928*npv -2.46685);
    if(runera==y2017C) METxcorr = -(-0.232763*npv +1.08318);
    if(runera==y2017C) METycorr = -(0.257719*npv -1.1745);
    if(runera==y2017D) METxcorr = -(-0.238067*npv +1.80541);
    if(runera==y2017D) METycorr = -(0.235989*npv -1.44354);
    if(runera==y2017E) METxcorr = -(-0.212352*npv +1.851);
    if(runera==y2017E) METycorr = -(0.157759*npv -0.478139);
    if(runera==y2017F) METxcorr = -(-0.232733*npv +2.24134);
    if(runera==y2017F) METycorr = -(0.213341*npv +0.684588);
    if(runera==y2018A) METxcorr = -(0.362865*npv -1.94505);
    if(runera==y2018A) METycorr = -(0.0709085*npv -0.307365);
    if(runera==y2018B) METxcorr = -(0.492083*npv -2.93552);
    if(runera==y2018B) METycorr = -(0.17874*npv -0.786844);
    if(runera==y2018C) METxcorr = -(0.521349*npv -1.44544);
    if(runera==y2018C) METycorr = -(0.118956*npv -1.96434);
    if(runera==y2018D) METxcorr = -(0.531151*npv -1.37568);
    if(runera==y2018D) METycorr = -(0.0884639*npv -1.57089);
    //if(runera==y2016MC) cout << "in 2016MC in !usemetv2 " << endl;
    if(runera==y2016MC) METxcorr = -(-0.195191*npv -0.170948);
    if(runera==y2016MC) METycorr = -(-0.0311891*npv +0.787627);
    if(runera==y2017MC) METxcorr = -(-0.217714*npv +0.493361);
    if(runera==y2017MC) METycorr = -(0.177058*npv -0.336648);
    if(runera==y2018MC) METxcorr = -(0.296713*npv -0.141506);
    if(runera==y2018MC) METycorr = -(0.115685*npv +0.0128193);
  }
  else {//these are the corrections for v2 MET recipe (currently recommended for 2017)
    //if(runera==y2016B) cout << "in data2016B in usemetv2 " << endl;
    if(runera==y2016B) METxcorr = -(-0.0374977*npv +0.00488262);
    if(runera==y2016B) METycorr = -(0.107373*npv +-0.00732239);
    if(runera==y2016C) METxcorr = -(-0.0832562*npv +0.550742);
    if(runera==y2016C) METycorr = -(0.142469*npv +-0.153718);
    if(runera==y2016D) METxcorr = -(-0.0400931*npv +0.753734);
    if(runera==y2016D) METycorr = -(0.127154*npv +0.0175228);
    if(runera==y2016E) METxcorr = -(-0.0409231*npv +0.755128);
    if(runera==y2016E) METycorr = -(0.168407*npv +0.126755);
    if(runera==y2016F) METxcorr = -(-0.0161259*npv +0.516919);
    if(runera==y2016F) METycorr = -(0.141176*npv +0.544062);
    if(runera==y2016G) METxcorr = -(0.0583851*npv +-0.0987447);
    if(runera==y2016G) METycorr = -(0.0641427*npv +0.319112);
    if(runera==y2016H) METxcorr = -(0.0706267*npv +-0.13118);
    if(runera==y2016H) METycorr = -(0.127481*npv +0.370786);
    if(runera==y2017B) METxcorr = -(-0.19563*npv +1.51859);
    if(runera==y2017B) METycorr = -(0.306987*npv +-1.84713);
    if(runera==y2017C) METxcorr = -(-0.161661*npv +0.589933);
    if(runera==y2017C) METycorr = -(0.233569*npv +-0.995546);
    if(runera==y2017D) METxcorr = -(-0.180911*npv +1.23553);
    if(runera==y2017D) METycorr = -(0.240155*npv +-1.27449);
    if(runera==y2017E) METxcorr = -(-0.149494*npv +0.901305);
    if(runera==y2017E) METycorr = -(0.178212*npv +-0.535537);
    if(runera==y2017F) METxcorr = -(-0.165154*npv +1.02018);
    if(runera==y2017F) METycorr = -(0.253794*npv +0.75776);
    if(runera==y2018A) METxcorr = -(0.362642*npv +-1.55094);
    if(runera==y2018A) METycorr = -(0.0737842*npv +-0.677209);
    if(runera==y2018B) METxcorr = -(0.485614*npv +-2.45706);
    if(runera==y2018B) METycorr = -(0.181619*npv +-1.00636);
    if(runera==y2018C) METxcorr = -(0.503638*npv +-1.01281);
    if(runera==y2018C) METycorr = -(0.147811*npv +-1.48941);
    if(runera==y2018D) METxcorr = -(0.520265*npv +-1.20322);
    if(runera==y2018D) METycorr = -(0.143919*npv +-0.979328);
    //if(runera==y2016MC) cout << "in 2016MC in usemetv2 " << endl; 
    if(runera==y2016MC) METxcorr = -(-0.159469*npv +-0.407022);
    if(runera==y2016MC) METycorr = -(-0.0405812*npv +0.570415);
    if(runera==y2017MC) METxcorr = -(-0.182569*npv +0.276542);
    if(runera==y2017MC) METycorr = -(0.155652*npv +-0.417633);
    if(runera==y2018MC) METxcorr = -(0.299448*npv +-0.13866);
    if(runera==y2018MC) METycorr = -(0.118785*npv +0.0889588);
  }
  
  double CorrectedMET_x = uncormet *cos( uncormet_phi)+METxcorr;
  double CorrectedMET_y = uncormet *sin( uncormet_phi)+METycorr;
  
  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;
  
  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;
  return TheXYCorr_Met_MetPhi;

  
}
