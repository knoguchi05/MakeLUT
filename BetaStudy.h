//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 28 16:43:39 2021 by ROOT version 6.24/00
// from TTree ntuple/AODVariables
// found on file: /gpfs/fs7001/knoguchi/NSW/Singlemuon/Ntuple/runSingleMu1000k_0727_Merge_Ntuple.root
//////////////////////////////////////////////////////////

#ifndef BetaStudy_h
#define BetaStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class BetaStudy {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventNumber;
   Int_t           runNumber;
   vector<string>  *trigname;
   vector<bool>    *isPassedTrig;
   vector<bool>    *isPassedL1_evt;
   vector<bool>    *isPassedSA_evt;
   vector<bool>    *isPassedCB_evt;
   vector<bool>    *isPassedSAIO_evt;
   vector<bool>    *isPassedCBIO_evt;
   vector<bool>    *isPassedEF_evt;
   vector<vector<bool> > *isPassedL1;
   vector<vector<bool> > *isPassedSA;
   vector<vector<bool> > *isPassedCB;
   vector<vector<bool> > *isPassedSAIO;
   vector<vector<bool> > *isPassedCBIO;
   vector<vector<bool> > *isPassedEF;
   Int_t           n_trig;
   vector<double>  *muon_pt;
   vector<double>  *muon_eta;
   vector<double>  *muon_phi;
   vector<double>  *muon_extEta;
   vector<double>  *muon_extPhi;
   vector<double>  *muon_charge;
   Double_t        segment_x[10];
   Double_t        segment_y[10];
   Double_t        segment_z[10];
   Double_t        segment_px[10];
   Double_t        segment_py[10];
   Double_t        segment_pz[10];
   Double_t        segment_chiSquared[10];
   Double_t        segment_numberDoF[10];
   Double_t        segment_sector[10];
   Double_t        segment_chamberIndex[10];
   Double_t        segment_etaIndex[10];
   Double_t        segment_nPrecisionHits[10];
   Double_t        segment_nPhiLayers[10];
   Double_t        segment_nTrigEtaLayers[10];
   vector<double>  *L1_eta;
   vector<double>  *L1_phi;
   vector<double>  *L1_thrValue;
   vector<int>     *L1_roiNum;
   vector<int>     *L1_roiSector;
   vector<int>     *L1_thrNumber;
   vector<bool>    *L1_isMoreCandInRoI;
   vector<int>     *SA_roiNumber;
   vector<int>     *SA_roiSector;
   vector<double>  *SA_pt;
   vector<double>  *SA_eta;
   vector<double>  *SA_phi;
   vector<double>  *SA_etaMS;
   vector<double>  *SA_phiMS;
   vector<double>  *SA_tgcPt;
   vector<double>  *SA_ptBarrelRadius;
   vector<double>  *SA_ptBarrelSagitta;
   vector<double>  *SA_ptEndcapAlpha;
   vector<double>  *SA_ptEndcapBeta;
   vector<double>  *SA_ptEndcapRadius;
   vector<double>  *SA_ptCSC;
   vector<double>  *SA_EndcapAlpha;
   vector<double>  *SA_EndcapBeta;
   vector<double>  *SA_EndcapRadius;
   vector<double>  *SA_etaBin;
   vector<double>  *SA_phiBin;
   vector<int>     *SA_sAddress;
   vector<float>   *SA_roiEta;
   vector<float>   *SA_roiPhi;
   vector<double>  *SA_superPointR_BI;
   vector<double>  *SA_superPointR_BM;
   vector<double>  *SA_superPointR_BO;
   vector<double>  *SA_superPointR_EI;
   vector<double>  *SA_superPointR_EM;
   vector<double>  *SA_superPointR_EO;
   vector<double>  *SA_superPointR_EE;
   vector<double>  *SA_superPointR_CSC;
   vector<double>  *SA_superPointR_BEE;
   vector<double>  *SA_superPointR_BME;
   vector<double>  *SA_superPointZ_BI;
   vector<double>  *SA_superPointZ_BM;
   vector<double>  *SA_superPointZ_BO;
   vector<double>  *SA_superPointZ_EI;
   vector<double>  *SA_superPointZ_EM;
   vector<double>  *SA_superPointZ_EO;
   vector<double>  *SA_superPointZ_EE;
   vector<double>  *SA_superPointZ_CSC;
   vector<double>  *SA_superPointZ_BEE;
   vector<double>  *SA_superPointZ_BME;
   vector<double>  *SA_superPointSlope_BI;
   vector<double>  *SA_superPointSlope_BM;
   vector<double>  *SA_superPointSlope_BO;
   vector<double>  *SA_superPointSlope_EI;
   vector<double>  *SA_superPointSlope_EM;
   vector<double>  *SA_superPointSlope_EO;
   vector<double>  *SA_superPointSlope_EE;
   vector<double>  *SA_superPointSlope_CSC;
   vector<double>  *SA_superPointSlope_BEE;
   vector<double>  *SA_superPointSlope_BME;
   vector<double>  *SA_superPointIntercept_BI;
   vector<double>  *SA_superPointIntercept_BM;
   vector<double>  *SA_superPointIntercept_BO;
   vector<double>  *SA_superPointIntercept_EI;
   vector<double>  *SA_superPointIntercept_EM;
   vector<double>  *SA_superPointIntercept_EO;
   vector<double>  *SA_superPointIntercept_EE;
   vector<double>  *SA_superPointIntercept_CSC;
   vector<double>  *SA_superPointIntercept_BEE;
   vector<double>  *SA_superPointIntercept_BME;
   vector<double>  *SA_superPointChi2_BI;
   vector<double>  *SA_superPointChi2_BM;
   vector<double>  *SA_superPointChi2_BO;
   vector<double>  *SA_superPointChi2_EI;
   vector<double>  *SA_superPointChi2_EM;
   vector<double>  *SA_superPointChi2_EO;
   vector<double>  *SA_superPointChi2_EE;
   vector<double>  *SA_superPointChi2_CSC;
   vector<double>  *SA_superPointChi2_BEE;
   vector<double>  *SA_superPointChi2_BME;
   vector<vector<float> > *SA_rpcHitX;
   vector<vector<float> > *SA_rpcHitY;
   vector<vector<float> > *SA_rpcHitZ;
   vector<vector<float> > *SA_rpcHitR;
   vector<vector<float> > *SA_rpcHitEta;
   vector<vector<float> > *SA_rpcHitPhi;
   vector<vector<unsigned int> > *SA_rpcHitMeasPhi;
   vector<vector<unsigned int> > *SA_rpcHitLayer;
   vector<vector<string> > *SA_rpcHitStationName;
   vector<vector<float> > *SA_tgcHitZ;
   vector<vector<float> > *SA_tgcHitR;
   vector<vector<float> > *SA_tgcHitEta;
   vector<vector<float> > *SA_tgcHitPhi;
   vector<vector<float> > *SA_tgcHitWidth;
   vector<vector<int> > *SA_tgcHitStationNum;
   vector<vector<bool> > *SA_tgcHitIsStrip;
   vector<vector<int> > *SA_tgcHitBCTag;
   vector<vector<bool> > *SA_tgcHitInRoad;
   vector<vector<int> > *SA_mdtHitIsOutlier;
   vector<vector<int> > *SA_mdtHitChamber;
   vector<vector<float> > *SA_mdtHitR;
   vector<vector<float> > *SA_mdtHitZ;
   vector<vector<float> > *SA_mdtHitPhi;
   vector<vector<float> > *SA_mdtHitResidual;
   vector<vector<float> > *SA_roadAw;
   vector<vector<float> > *SA_roadBw;
   vector<vector<float> > *SA_zMin;
   vector<vector<float> > *SA_zMax;
   vector<vector<float> > *SA_rMin;
   vector<vector<float> > *SA_rMax;
   vector<vector<float> > *SA_etaMin;
   vector<vector<float> > *SA_etaMax;
   vector<vector<float> > *SA_stgcClusterZ;
   vector<vector<float> > *SA_stgcClusterR;
   vector<vector<float> > *SA_stgcClusterEta;
   vector<vector<float> > *SA_stgcClusterPhi;
   vector<vector<float> > *SA_stgcClusterResidualR;
   vector<vector<float> > *SA_stgcClusterResidualPhi;
   vector<vector<int> > *SA_stgcClusterStationEta;
   vector<vector<int> > *SA_stgcClusterStationPhi;
   vector<vector<int> > *SA_stgcClusterStationName;
   vector<vector<int> > *SA_stgcClusterType;
   vector<vector<int> > *SA_stgcClusterIsOutlier;
   vector<vector<unsigned int> > *SA_stgcClusterLayer;
   vector<vector<float> > *SA_mmClusterZ;
   vector<vector<float> > *SA_mmClusterR;
   vector<vector<float> > *SA_mmClusterEta;
   vector<vector<float> > *SA_mmClusterPhi;
   vector<vector<float> > *SA_mmClusterResidualR;
   vector<vector<float> > *SA_mmClusterResidualPhi;
   vector<vector<int> > *SA_mmClusterStationEta;
   vector<vector<int> > *SA_mmClusterStationPhi;
   vector<vector<int> > *SA_mmClusterStationName;
   vector<vector<int> > *SA_mmClusterIsOutlier;
   vector<vector<unsigned int> > *SA_mmClusterLayer;
   vector<double>  *CB_pt;
   vector<double>  *CB_eta;
   vector<double>  *CB_phi;
   vector<double>  *CB_idpt;
   vector<double>  *CB_ideta;
   vector<double>  *CB_idphi;
   vector<int>     *CB_roiNumber;
   vector<int>     *CB_roiSector;
   vector<double>  *CBIO_pt;
   vector<double>  *CBIO_eta;
   vector<double>  *CBIO_phi;
   vector<double>  *CBIO_idpt;
   vector<double>  *CBIO_ideta;
   vector<double>  *CBIO_idphi;
   vector<int>     *CBIO_roiNumber;
   vector<int>     *CBIO_roiSector;
   vector<double>  *EF_pt;
   vector<double>  *EF_eta;
   vector<double>  *EF_phi;

   // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_trigname;   //!
   TBranch        *b_isPassedTrig;   //!
   TBranch        *b_isPassedL1_evt;   //!
   TBranch        *b_isPassedSA_evt;   //!
   TBranch        *b_isPassedCB_evt;   //!
   TBranch        *b_isPassedSAIO_evt;   //!
   TBranch        *b_isPassedCBIO_evt;   //!
   TBranch        *b_isPassedEF_evt;   //!
   TBranch        *b_isPassedL1;   //!
   TBranch        *b_isPassedSA;   //!
   TBranch        *b_isPassedCB;   //!
   TBranch        *b_isPassedSAIO;   //!
   TBranch        *b_isPassedCBIO;   //!
   TBranch        *b_isPassedEF;   //!
   TBranch        *b_n_trig;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_extEta;   //!
   TBranch        *b_muon_extPhi;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_segment_x;   //!
   TBranch        *b_segment_y;   //!
   TBranch        *b_segment_z;   //!
   TBranch        *b_segment_px;   //!
   TBranch        *b_segment_py;   //!
   TBranch        *b_segment_pz;   //!
   TBranch        *b_segment_chiSquared;   //!
   TBranch        *b_segment_numberDoF;   //!
   TBranch        *b_segment_sector;   //!
   TBranch        *b_segment_chamberIndex;   //!
   TBranch        *b_segment_etaIndex;   //!
   TBranch        *b_segment_nPrecisionHits;   //!
   TBranch        *b_segment_nPhiLayers;   //!
   TBranch        *b_segment_nTrigEtaLayers;   //!
   TBranch        *b_L1_eta;   //!
   TBranch        *b_L1_phi;   //!
   TBranch        *b_L1_thrValue;   //!
   TBranch        *b_L1_roiNum;   //!
   TBranch        *b_L1_roiSector;   //!
   TBranch        *b_L1_thrNumber;   //!
   TBranch        *b_L1_isMoreCandInRoI;   //!
   TBranch        *b_SA_roiNumber;   //!
   TBranch        *b_SA_roiSector;   //!
   TBranch        *b_SA_pt;   //!
   TBranch        *b_SA_eta;   //!
   TBranch        *b_SA_phi;   //!
   TBranch        *b_SA_etaMS;   //!
   TBranch        *b_SA_phiMS;   //!
   TBranch        *b_SA_tgcPt;   //!
   TBranch        *b_SA_ptBarrelRadius;   //!
   TBranch        *b_SA_ptBarrelSagitta;   //!
   TBranch        *b_SA_ptEndcapAlpha;   //!
   TBranch        *b_SA_ptEndcapBeta;   //!
   TBranch        *b_SA_ptEndcapRadius;   //!
   TBranch        *b_SA_ptCSC;   //!
   TBranch        *b_SA_EndcapAlpha;   //!
   TBranch        *b_SA_EndcapBeta;   //!
   TBranch        *b_SA_EndcapRadius;   //!
   TBranch        *b_SA_etaBin;   //!
   TBranch        *b_SA_phiBin;   //!
   TBranch        *b_SA_sAddress;   //!
   TBranch        *b_SA_roiEta;   //!
   TBranch        *b_SA_roiPhi;   //!
   TBranch        *b_SA_superPointR_BI;   //!
   TBranch        *b_SA_superPointR_BM;   //!
   TBranch        *b_SA_superPointR_BO;   //!
   TBranch        *b_SA_superPointR_EI;   //!
   TBranch        *b_SA_superPointR_EM;   //!
   TBranch        *b_SA_superPointR_EO;   //!
   TBranch        *b_SA_superPointR_EE;   //!
   TBranch        *b_SA_superPointR_CSC;   //!
   TBranch        *b_SA_superPointR_BEE;   //!
   TBranch        *b_SA_superPointR_BME;   //!
   TBranch        *b_SA_superPointZ_BI;   //!
   TBranch        *b_SA_superPointZ_BM;   //!
   TBranch        *b_SA_superPointZ_BO;   //!
   TBranch        *b_SA_superPointZ_EI;   //!
   TBranch        *b_SA_superPointZ_EM;   //!
   TBranch        *b_SA_superPointZ_EO;   //!
   TBranch        *b_SA_superPointZ_EE;   //!
   TBranch        *b_SA_superPointZ_CSC;   //!
   TBranch        *b_SA_superPointZ_BEE;   //!
   TBranch        *b_SA_superPointZ_BME;   //!
   TBranch        *b_SA_superPointSlope_BI;   //!
   TBranch        *b_SA_superPointSlope_BM;   //!
   TBranch        *b_SA_superPointSlope_BO;   //!
   TBranch        *b_SA_superPointSlope_EI;   //!
   TBranch        *b_SA_superPointSlope_EM;   //!
   TBranch        *b_SA_superPointSlope_EO;   //!
   TBranch        *b_SA_superPointSlope_EE;   //!
   TBranch        *b_SA_superPointSlope_CSC;   //!
   TBranch        *b_SA_superPointSlope_BEE;   //!
   TBranch        *b_SA_superPointSlope_BME;   //!
   TBranch        *b_SA_superPointIntercept_BI;   //!
   TBranch        *b_SA_superPointIntercept_BM;   //!
   TBranch        *b_SA_superPointIntercept_BO;   //!
   TBranch        *b_SA_superPointIntercept_EI;   //!
   TBranch        *b_SA_superPointIntercept_EM;   //!
   TBranch        *b_SA_superPointIntercept_EO;   //!
   TBranch        *b_SA_superPointIntercept_EE;   //!
   TBranch        *b_SA_superPointIntercept_CSC;   //!
   TBranch        *b_SA_superPointIntercept_BEE;   //!
   TBranch        *b_SA_superPointIntercept_BME;   //!
   TBranch        *b_SA_superPointChi2_BI;   //!
   TBranch        *b_SA_superPointChi2_BM;   //!
   TBranch        *b_SA_superPointChi2_BO;   //!
   TBranch        *b_SA_superPointChi2_EI;   //!
   TBranch        *b_SA_superPointChi2_EM;   //!
   TBranch        *b_SA_superPointChi2_EO;   //!
   TBranch        *b_SA_superPointChi2_EE;   //!
   TBranch        *b_SA_superPointChi2_CSC;   //!
   TBranch        *b_SA_superPointChi2_BEE;   //!
   TBranch        *b_SA_superPointChi2_BME;   //!
   TBranch        *b_SA_rpcHitX;   //!
   TBranch        *b_SA_rpcHitY;   //!
   TBranch        *b_SA_rpcHitZ;   //!
   TBranch        *b_SA_rpcHitR;   //!
   TBranch        *b_SA_rpcHitEta;   //!
   TBranch        *b_SA_rpcHitPhi;   //!
   TBranch        *b_SA_rpcHitMeasPhi;   //!
   TBranch        *b_SA_rpcHitLayer;   //!
   TBranch        *b_SA_rpcHitStationName;   //!
   TBranch        *b_SA_tgcHitZ;   //!
   TBranch        *b_SA_tgcHitR;   //!
   TBranch        *b_SA_tgcHitEta;   //!
   TBranch        *b_SA_tgcHitPhi;   //!
   TBranch        *b_SA_tgcHitWidth;   //!
   TBranch        *b_SA_tgcHitStationNum;   //!
   TBranch        *b_SA_tgcHitIsStrip;   //!
   TBranch        *b_SA_tgcHitBCTag;   //!
   TBranch        *b_SA_tgcHitInRoad;   //!
   TBranch        *b_SA_mdtHitIsOutlier;   //!
   TBranch        *b_SA_mdtHitChamber;   //!
   TBranch        *b_SA_mdtHitR;   //!
   TBranch        *b_SA_mdtHitZ;   //!
   TBranch        *b_SA_mdtHitPhi;   //!
   TBranch        *b_SA_mdtHitResidual;   //!
   TBranch        *b_SA_roadAw;   //!
   TBranch        *b_SA_roadBw;   //!
   TBranch        *b_SA_zMin;   //!
   TBranch        *b_SA_zMax;   //!
   TBranch        *b_SA_rMin;   //!
   TBranch        *b_SA_rMax;   //!
   TBranch        *b_SA_etaMin;   //!
   TBranch        *b_SA_etaMax;   //!
   TBranch        *b_SA_stgcClusterZ;   //!
   TBranch        *b_SA_stgcClusterR;   //!
   TBranch        *b_SA_stgcClusterEta;   //!
   TBranch        *b_SA_stgcClusterPhi;   //!
   TBranch        *b_SA_stgcClusterResidualR;   //!
   TBranch        *b_SA_stgcClusterResidualPhi;   //!
   TBranch        *b_SA_stgcClusterStationEta;   //!
   TBranch        *b_SA_stgcClusterStationPhi;   //!
   TBranch        *b_SA_stgcClusterStationName;   //!
   TBranch        *b_SA_stgcClusterType;   //!
   TBranch        *b_SA_stgcClusterIsOutlier;   //!
   TBranch        *b_SA_stgcClusterLayer;   //!
   TBranch        *b_SA_mmClusterZ;   //!
   TBranch        *b_SA_mmClusterR;   //!
   TBranch        *b_SA_mmClusterEta;   //!
   TBranch        *b_SA_mmClusterPhi;   //!
   TBranch        *b_SA_mmClusterResidualR;   //!
   TBranch        *b_SA_mmClusterResidualPhi;   //!
   TBranch        *b_SA_mmClusterStationEta;   //!
   TBranch        *b_SA_mmClusterStationPhi;   //!
   TBranch        *b_SA_mmClusterStationName;   //!
   TBranch        *b_SA_mmClusterIsOutlier;   //!
   TBranch        *b_SA_mmClusterLayer;   //!
   TBranch        *b_CB_pt;   //!
   TBranch        *b_CB_eta;   //!
   TBranch        *b_CB_phi;   //!
   TBranch        *b_CB_idpt;   //!
   TBranch        *b_CB_ideta;   //!
   TBranch        *b_CB_idphi;   //!
   TBranch        *b_CB_roiNumber;   //!
   TBranch        *b_CB_roiSector;   //!
   TBranch        *b_CBIO_pt;   //!
   TBranch        *b_CBIO_eta;   //!
   TBranch        *b_CBIO_phi;   //!
   TBranch        *b_CBIO_idpt;   //!
   TBranch        *b_CBIO_ideta;   //!
   TBranch        *b_CBIO_idphi;   //!
   TBranch        *b_CBIO_roiNumber;   //!
   TBranch        *b_CBIO_roiSector;   //!
   TBranch        *b_EF_pt;   //!
   TBranch        *b_EF_eta;   //!
   TBranch        *b_EF_phi;   //!

   BetaStudy(TTree *tree=0);
   virtual ~BetaStudy();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BetaStudy_cxx
BetaStudy::BetaStudy(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/gpfs/fs7001/knoguchi/NSW/gridSinglemuon10M/Ntuple/gridSingleMuon7M0917_merge.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/gpfs/fs7001/knoguchi/NSW/gridSinglemuon10M/Ntuple/gridSingleMuon7M0917_merge.root");
      }
      f->GetObject("ntuple",tree);

   }
   Init(tree);
}

BetaStudy::~BetaStudy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BetaStudy::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BetaStudy::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BetaStudy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trigname = 0;
   isPassedTrig = 0;
   isPassedL1_evt = 0;
   isPassedSA_evt = 0;
   isPassedCB_evt = 0;
   isPassedSAIO_evt = 0;
   isPassedCBIO_evt = 0;
   isPassedEF_evt = 0;
   isPassedL1 = 0;
   isPassedSA = 0;
   isPassedCB = 0;
   isPassedSAIO = 0;
   isPassedCBIO = 0;
   isPassedEF = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_extEta = 0;
   muon_extPhi = 0;
   muon_charge = 0;
   L1_eta = 0;
   L1_phi = 0;
   L1_thrValue = 0;
   L1_roiNum = 0;
   L1_roiSector = 0;
   L1_thrNumber = 0;
   L1_isMoreCandInRoI = 0;
   SA_roiNumber = 0;
   SA_roiSector = 0;
   SA_pt = 0;
   SA_eta = 0;
   SA_phi = 0;
   SA_etaMS = 0;
   SA_phiMS = 0;
   SA_tgcPt = 0;
   SA_ptBarrelRadius = 0;
   SA_ptBarrelSagitta = 0;
   SA_ptEndcapAlpha = 0;
   SA_ptEndcapBeta = 0;
   SA_ptEndcapRadius = 0;
   SA_ptCSC = 0;
   SA_EndcapAlpha = 0;
   SA_EndcapBeta = 0;
   SA_EndcapRadius = 0;
   SA_etaBin = 0;
   SA_phiBin = 0;
   SA_sAddress = 0;
   SA_roiEta = 0;
   SA_roiPhi = 0;
   SA_superPointR_BI = 0;
   SA_superPointR_BM = 0;
   SA_superPointR_BO = 0;
   SA_superPointR_EI = 0;
   SA_superPointR_EM = 0;
   SA_superPointR_EO = 0;
   SA_superPointR_EE = 0;
   SA_superPointR_CSC = 0;
   SA_superPointR_BEE = 0;
   SA_superPointR_BME = 0;
   SA_superPointZ_BI = 0;
   SA_superPointZ_BM = 0;
   SA_superPointZ_BO = 0;
   SA_superPointZ_EI = 0;
   SA_superPointZ_EM = 0;
   SA_superPointZ_EO = 0;
   SA_superPointZ_EE = 0;
   SA_superPointZ_CSC = 0;
   SA_superPointZ_BEE = 0;
   SA_superPointZ_BME = 0;
   SA_superPointSlope_BI = 0;
   SA_superPointSlope_BM = 0;
   SA_superPointSlope_BO = 0;
   SA_superPointSlope_EI = 0;
   SA_superPointSlope_EM = 0;
   SA_superPointSlope_EO = 0;
   SA_superPointSlope_EE = 0;
   SA_superPointSlope_CSC = 0;
   SA_superPointSlope_BEE = 0;
   SA_superPointSlope_BME = 0;
   SA_superPointIntercept_BI = 0;
   SA_superPointIntercept_BM = 0;
   SA_superPointIntercept_BO = 0;
   SA_superPointIntercept_EI = 0;
   SA_superPointIntercept_EM = 0;
   SA_superPointIntercept_EO = 0;
   SA_superPointIntercept_EE = 0;
   SA_superPointIntercept_CSC = 0;
   SA_superPointIntercept_BEE = 0;
   SA_superPointIntercept_BME = 0;
   SA_superPointChi2_BI = 0;
   SA_superPointChi2_BM = 0;
   SA_superPointChi2_BO = 0;
   SA_superPointChi2_EI = 0;
   SA_superPointChi2_EM = 0;
   SA_superPointChi2_EO = 0;
   SA_superPointChi2_EE = 0;
   SA_superPointChi2_CSC = 0;
   SA_superPointChi2_BEE = 0;
   SA_superPointChi2_BME = 0;
   SA_rpcHitX = 0;
   SA_rpcHitY = 0;
   SA_rpcHitZ = 0;
   SA_rpcHitR = 0;
   SA_rpcHitEta = 0;
   SA_rpcHitPhi = 0;
   SA_rpcHitMeasPhi = 0;
   SA_rpcHitLayer = 0;
   SA_rpcHitStationName = 0;
   SA_tgcHitZ = 0;
   SA_tgcHitR = 0;
   SA_tgcHitEta = 0;
   SA_tgcHitPhi = 0;
   SA_tgcHitWidth = 0;
   SA_tgcHitStationNum = 0;
   SA_tgcHitIsStrip = 0;
   SA_tgcHitBCTag = 0;
   SA_tgcHitInRoad = 0;
   SA_mdtHitIsOutlier = 0;
   SA_mdtHitChamber = 0;
   SA_mdtHitR = 0;
   SA_mdtHitZ = 0;
   SA_mdtHitPhi = 0;
   SA_mdtHitResidual = 0;
   SA_roadAw = 0;
   SA_roadBw = 0;
   SA_zMin = 0;
   SA_zMax = 0;
   SA_rMin = 0;
   SA_rMax = 0;
   SA_etaMin = 0;
   SA_etaMax = 0;
   SA_stgcClusterZ = 0;
   SA_stgcClusterR = 0;
   SA_stgcClusterEta = 0;
   SA_stgcClusterPhi = 0;
   SA_stgcClusterResidualR = 0;
   SA_stgcClusterResidualPhi = 0;
   SA_stgcClusterStationEta = 0;
   SA_stgcClusterStationPhi = 0;
   SA_stgcClusterStationName = 0;
   SA_stgcClusterType = 0;
   SA_stgcClusterIsOutlier = 0;
   SA_stgcClusterLayer = 0;
   SA_mmClusterZ = 0;
   SA_mmClusterR = 0;
   SA_mmClusterEta = 0;
   SA_mmClusterPhi = 0;
   SA_mmClusterResidualR = 0;
   SA_mmClusterResidualPhi = 0;
   SA_mmClusterStationEta = 0;
   SA_mmClusterStationPhi = 0;
   SA_mmClusterStationName = 0;
   SA_mmClusterIsOutlier = 0;
   SA_mmClusterLayer = 0;
   CB_pt = 0;
   CB_eta = 0;
   CB_phi = 0;
   CB_idpt = 0;
   CB_ideta = 0;
   CB_idphi = 0;
   CB_roiNumber = 0;
   CB_roiSector = 0;
   CBIO_pt = 0;
   CBIO_eta = 0;
   CBIO_phi = 0;
   CBIO_idpt = 0;
   CBIO_ideta = 0;
   CBIO_idphi = 0;
   CBIO_roiNumber = 0;
   CBIO_roiSector = 0;
   EF_pt = 0;
   EF_eta = 0;
   EF_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("trigname", &trigname, &b_trigname);
   fChain->SetBranchAddress("isPassedTrig", &isPassedTrig, &b_isPassedTrig);
   fChain->SetBranchAddress("isPassedL1_evt", &isPassedL1_evt, &b_isPassedL1_evt);
   fChain->SetBranchAddress("isPassedSA_evt", &isPassedSA_evt, &b_isPassedSA_evt);
   fChain->SetBranchAddress("isPassedCB_evt", &isPassedCB_evt, &b_isPassedCB_evt);
   fChain->SetBranchAddress("isPassedSAIO_evt", &isPassedSAIO_evt, &b_isPassedSAIO_evt);
   fChain->SetBranchAddress("isPassedCBIO_evt", &isPassedCBIO_evt, &b_isPassedCBIO_evt);
   fChain->SetBranchAddress("isPassedEF_evt", &isPassedEF_evt, &b_isPassedEF_evt);
   fChain->SetBranchAddress("isPassedL1", &isPassedL1, &b_isPassedL1);
   fChain->SetBranchAddress("isPassedSA", &isPassedSA, &b_isPassedSA);
   fChain->SetBranchAddress("isPassedCB", &isPassedCB, &b_isPassedCB);
   fChain->SetBranchAddress("isPassedSAIO", &isPassedSAIO, &b_isPassedSAIO);
   fChain->SetBranchAddress("isPassedCBIO", &isPassedCBIO, &b_isPassedCBIO);
   fChain->SetBranchAddress("isPassedEF", &isPassedEF, &b_isPassedEF);
   fChain->SetBranchAddress("n_trig", &n_trig, &b_n_trig);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_extEta", &muon_extEta, &b_muon_extEta);
   fChain->SetBranchAddress("muon_extPhi", &muon_extPhi, &b_muon_extPhi);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("segment_x", segment_x, &b_segment_x);
   fChain->SetBranchAddress("segment_y", segment_y, &b_segment_y);
   fChain->SetBranchAddress("segment_z", segment_z, &b_segment_z);
   fChain->SetBranchAddress("segment_px", segment_px, &b_segment_px);
   fChain->SetBranchAddress("segment_py", segment_py, &b_segment_py);
   fChain->SetBranchAddress("segment_pz", segment_pz, &b_segment_pz);
   fChain->SetBranchAddress("segment_chiSquared", segment_chiSquared, &b_segment_chiSquared);
   fChain->SetBranchAddress("segment_numberDoF", segment_numberDoF, &b_segment_numberDoF);
   fChain->SetBranchAddress("segment_sector", segment_sector, &b_segment_sector);
   fChain->SetBranchAddress("segment_chamberIndex", segment_chamberIndex, &b_segment_chamberIndex);
   fChain->SetBranchAddress("segment_etaIndex", segment_etaIndex, &b_segment_etaIndex);
   fChain->SetBranchAddress("segment_nPrecisionHits", segment_nPrecisionHits, &b_segment_nPrecisionHits);
   fChain->SetBranchAddress("segment_nPhiLayers", segment_nPhiLayers, &b_segment_nPhiLayers);
   fChain->SetBranchAddress("segment_nTrigEtaLayers", segment_nTrigEtaLayers, &b_segment_nTrigEtaLayers);
   fChain->SetBranchAddress("L1_eta", &L1_eta, &b_L1_eta);
   fChain->SetBranchAddress("L1_phi", &L1_phi, &b_L1_phi);
   fChain->SetBranchAddress("L1_thrValue", &L1_thrValue, &b_L1_thrValue);
   fChain->SetBranchAddress("L1_roiNum", &L1_roiNum, &b_L1_roiNum);
   fChain->SetBranchAddress("L1_roiSector", &L1_roiSector, &b_L1_roiSector);
   fChain->SetBranchAddress("L1_thrNumber", &L1_thrNumber, &b_L1_thrNumber);
   fChain->SetBranchAddress("L1_isMoreCandInRoI", &L1_isMoreCandInRoI, &b_L1_isMoreCandInRoI);
   fChain->SetBranchAddress("SA_roiNumber", &SA_roiNumber, &b_SA_roiNumber);
   fChain->SetBranchAddress("SA_roiSector", &SA_roiSector, &b_SA_roiSector);
   fChain->SetBranchAddress("SA_pt", &SA_pt, &b_SA_pt);
   fChain->SetBranchAddress("SA_eta", &SA_eta, &b_SA_eta);
   fChain->SetBranchAddress("SA_phi", &SA_phi, &b_SA_phi);
   fChain->SetBranchAddress("SA_etaMS", &SA_etaMS, &b_SA_etaMS);
   fChain->SetBranchAddress("SA_phiMS", &SA_phiMS, &b_SA_phiMS);
   fChain->SetBranchAddress("SA_tgcPt", &SA_tgcPt, &b_SA_tgcPt);
   fChain->SetBranchAddress("SA_ptBarrelRadius", &SA_ptBarrelRadius, &b_SA_ptBarrelRadius);
   fChain->SetBranchAddress("SA_ptBarrelSagitta", &SA_ptBarrelSagitta, &b_SA_ptBarrelSagitta);
   fChain->SetBranchAddress("SA_ptEndcapAlpha", &SA_ptEndcapAlpha, &b_SA_ptEndcapAlpha);
   fChain->SetBranchAddress("SA_ptEndcapBeta", &SA_ptEndcapBeta, &b_SA_ptEndcapBeta);
   fChain->SetBranchAddress("SA_ptEndcapRadius", &SA_ptEndcapRadius, &b_SA_ptEndcapRadius);
   fChain->SetBranchAddress("SA_ptCSC", &SA_ptCSC, &b_SA_ptCSC);
   fChain->SetBranchAddress("SA_EndcapAlpha", &SA_EndcapAlpha, &b_SA_EndcapAlpha);
   fChain->SetBranchAddress("SA_EndcapBeta", &SA_EndcapBeta, &b_SA_EndcapBeta);
   fChain->SetBranchAddress("SA_EndcapRadius", &SA_EndcapRadius, &b_SA_EndcapRadius);
   fChain->SetBranchAddress("SA_etaBin", &SA_etaBin, &b_SA_etaBin);
   fChain->SetBranchAddress("SA_phiBin", &SA_phiBin, &b_SA_phiBin);
   fChain->SetBranchAddress("SA_sAddress", &SA_sAddress, &b_SA_sAddress);
   fChain->SetBranchAddress("SA_roiEta", &SA_roiEta, &b_SA_roiEta);
   fChain->SetBranchAddress("SA_roiPhi", &SA_roiPhi, &b_SA_roiPhi);
   fChain->SetBranchAddress("SA_superPointR_BI", &SA_superPointR_BI, &b_SA_superPointR_BI);
   fChain->SetBranchAddress("SA_superPointR_BM", &SA_superPointR_BM, &b_SA_superPointR_BM);
   fChain->SetBranchAddress("SA_superPointR_BO", &SA_superPointR_BO, &b_SA_superPointR_BO);
   fChain->SetBranchAddress("SA_superPointR_EI", &SA_superPointR_EI, &b_SA_superPointR_EI);
   fChain->SetBranchAddress("SA_superPointR_EM", &SA_superPointR_EM, &b_SA_superPointR_EM);
   fChain->SetBranchAddress("SA_superPointR_EO", &SA_superPointR_EO, &b_SA_superPointR_EO);
   fChain->SetBranchAddress("SA_superPointR_EE", &SA_superPointR_EE, &b_SA_superPointR_EE);
   fChain->SetBranchAddress("SA_superPointR_CSC", &SA_superPointR_CSC, &b_SA_superPointR_CSC);
   fChain->SetBranchAddress("SA_superPointR_BEE", &SA_superPointR_BEE, &b_SA_superPointR_BEE);
   fChain->SetBranchAddress("SA_superPointR_BME", &SA_superPointR_BME, &b_SA_superPointR_BME);
   fChain->SetBranchAddress("SA_superPointZ_BI", &SA_superPointZ_BI, &b_SA_superPointZ_BI);
   fChain->SetBranchAddress("SA_superPointZ_BM", &SA_superPointZ_BM, &b_SA_superPointZ_BM);
   fChain->SetBranchAddress("SA_superPointZ_BO", &SA_superPointZ_BO, &b_SA_superPointZ_BO);
   fChain->SetBranchAddress("SA_superPointZ_EI", &SA_superPointZ_EI, &b_SA_superPointZ_EI);
   fChain->SetBranchAddress("SA_superPointZ_EM", &SA_superPointZ_EM, &b_SA_superPointZ_EM);
   fChain->SetBranchAddress("SA_superPointZ_EO", &SA_superPointZ_EO, &b_SA_superPointZ_EO);
   fChain->SetBranchAddress("SA_superPointZ_EE", &SA_superPointZ_EE, &b_SA_superPointZ_EE);
   fChain->SetBranchAddress("SA_superPointZ_CSC", &SA_superPointZ_CSC, &b_SA_superPointZ_CSC);
   fChain->SetBranchAddress("SA_superPointZ_BEE", &SA_superPointZ_BEE, &b_SA_superPointZ_BEE);
   fChain->SetBranchAddress("SA_superPointZ_BME", &SA_superPointZ_BME, &b_SA_superPointZ_BME);
   fChain->SetBranchAddress("SA_superPointSlope_BI", &SA_superPointSlope_BI, &b_SA_superPointSlope_BI);
   fChain->SetBranchAddress("SA_superPointSlope_BM", &SA_superPointSlope_BM, &b_SA_superPointSlope_BM);
   fChain->SetBranchAddress("SA_superPointSlope_BO", &SA_superPointSlope_BO, &b_SA_superPointSlope_BO);
   fChain->SetBranchAddress("SA_superPointSlope_EI", &SA_superPointSlope_EI, &b_SA_superPointSlope_EI);
   fChain->SetBranchAddress("SA_superPointSlope_EM", &SA_superPointSlope_EM, &b_SA_superPointSlope_EM);
   fChain->SetBranchAddress("SA_superPointSlope_EO", &SA_superPointSlope_EO, &b_SA_superPointSlope_EO);
   fChain->SetBranchAddress("SA_superPointSlope_EE", &SA_superPointSlope_EE, &b_SA_superPointSlope_EE);
   fChain->SetBranchAddress("SA_superPointSlope_CSC", &SA_superPointSlope_CSC, &b_SA_superPointSlope_CSC);
   fChain->SetBranchAddress("SA_superPointSlope_BEE", &SA_superPointSlope_BEE, &b_SA_superPointSlope_BEE);
   fChain->SetBranchAddress("SA_superPointSlope_BME", &SA_superPointSlope_BME, &b_SA_superPointSlope_BME);
   fChain->SetBranchAddress("SA_superPointIntercept_BI", &SA_superPointIntercept_BI, &b_SA_superPointIntercept_BI);
   fChain->SetBranchAddress("SA_superPointIntercept_BM", &SA_superPointIntercept_BM, &b_SA_superPointIntercept_BM);
   fChain->SetBranchAddress("SA_superPointIntercept_BO", &SA_superPointIntercept_BO, &b_SA_superPointIntercept_BO);
   fChain->SetBranchAddress("SA_superPointIntercept_EI", &SA_superPointIntercept_EI, &b_SA_superPointIntercept_EI);
   fChain->SetBranchAddress("SA_superPointIntercept_EM", &SA_superPointIntercept_EM, &b_SA_superPointIntercept_EM);
   fChain->SetBranchAddress("SA_superPointIntercept_EO", &SA_superPointIntercept_EO, &b_SA_superPointIntercept_EO);
   fChain->SetBranchAddress("SA_superPointIntercept_EE", &SA_superPointIntercept_EE, &b_SA_superPointIntercept_EE);
   fChain->SetBranchAddress("SA_superPointIntercept_CSC", &SA_superPointIntercept_CSC, &b_SA_superPointIntercept_CSC);
   fChain->SetBranchAddress("SA_superPointIntercept_BEE", &SA_superPointIntercept_BEE, &b_SA_superPointIntercept_BEE);
   fChain->SetBranchAddress("SA_superPointIntercept_BME", &SA_superPointIntercept_BME, &b_SA_superPointIntercept_BME);
   fChain->SetBranchAddress("SA_superPointChi2_BI", &SA_superPointChi2_BI, &b_SA_superPointChi2_BI);
   fChain->SetBranchAddress("SA_superPointChi2_BM", &SA_superPointChi2_BM, &b_SA_superPointChi2_BM);
   fChain->SetBranchAddress("SA_superPointChi2_BO", &SA_superPointChi2_BO, &b_SA_superPointChi2_BO);
   fChain->SetBranchAddress("SA_superPointChi2_EI", &SA_superPointChi2_EI, &b_SA_superPointChi2_EI);
   fChain->SetBranchAddress("SA_superPointChi2_EM", &SA_superPointChi2_EM, &b_SA_superPointChi2_EM);
   fChain->SetBranchAddress("SA_superPointChi2_EO", &SA_superPointChi2_EO, &b_SA_superPointChi2_EO);
   fChain->SetBranchAddress("SA_superPointChi2_EE", &SA_superPointChi2_EE, &b_SA_superPointChi2_EE);
   fChain->SetBranchAddress("SA_superPointChi2_CSC", &SA_superPointChi2_CSC, &b_SA_superPointChi2_CSC);
   fChain->SetBranchAddress("SA_superPointChi2_BEE", &SA_superPointChi2_BEE, &b_SA_superPointChi2_BEE);
   fChain->SetBranchAddress("SA_superPointChi2_BME", &SA_superPointChi2_BME, &b_SA_superPointChi2_BME);
   fChain->SetBranchAddress("SA_rpcHitX", &SA_rpcHitX, &b_SA_rpcHitX);
   fChain->SetBranchAddress("SA_rpcHitY", &SA_rpcHitY, &b_SA_rpcHitY);
   fChain->SetBranchAddress("SA_rpcHitZ", &SA_rpcHitZ, &b_SA_rpcHitZ);
   fChain->SetBranchAddress("SA_rpcHitR", &SA_rpcHitR, &b_SA_rpcHitR);
   fChain->SetBranchAddress("SA_rpcHitEta", &SA_rpcHitEta, &b_SA_rpcHitEta);
   fChain->SetBranchAddress("SA_rpcHitPhi", &SA_rpcHitPhi, &b_SA_rpcHitPhi);
   fChain->SetBranchAddress("SA_rpcHitMeasPhi", &SA_rpcHitMeasPhi, &b_SA_rpcHitMeasPhi);
   fChain->SetBranchAddress("SA_rpcHitLayer", &SA_rpcHitLayer, &b_SA_rpcHitLayer);
   fChain->SetBranchAddress("SA_rpcHitStationName", &SA_rpcHitStationName, &b_SA_rpcHitStationName);
   fChain->SetBranchAddress("SA_tgcHitZ", &SA_tgcHitZ, &b_SA_tgcHitZ);
   fChain->SetBranchAddress("SA_tgcHitR", &SA_tgcHitR, &b_SA_tgcHitR);
   fChain->SetBranchAddress("SA_tgcHitEta", &SA_tgcHitEta, &b_SA_tgcHitEta);
   fChain->SetBranchAddress("SA_tgcHitPhi", &SA_tgcHitPhi, &b_SA_tgcHitPhi);
   fChain->SetBranchAddress("SA_tgcHitWidth", &SA_tgcHitWidth, &b_SA_tgcHitWidth);
   fChain->SetBranchAddress("SA_tgcHitStationNum", &SA_tgcHitStationNum, &b_SA_tgcHitStationNum);
   fChain->SetBranchAddress("SA_tgcHitIsStrip", &SA_tgcHitIsStrip, &b_SA_tgcHitIsStrip);
   fChain->SetBranchAddress("SA_tgcHitBCTag", &SA_tgcHitBCTag, &b_SA_tgcHitBCTag);
   fChain->SetBranchAddress("SA_tgcHitInRoad", &SA_tgcHitInRoad, &b_SA_tgcHitInRoad);
   fChain->SetBranchAddress("SA_mdtHitIsOutlier", &SA_mdtHitIsOutlier, &b_SA_mdtHitIsOutlier);
   fChain->SetBranchAddress("SA_mdtHitChamber", &SA_mdtHitChamber, &b_SA_mdtHitChamber);
   fChain->SetBranchAddress("SA_mdtHitR", &SA_mdtHitR, &b_SA_mdtHitR);
   fChain->SetBranchAddress("SA_mdtHitZ", &SA_mdtHitZ, &b_SA_mdtHitZ);
   fChain->SetBranchAddress("SA_mdtHitPhi", &SA_mdtHitPhi, &b_SA_mdtHitPhi);
   fChain->SetBranchAddress("SA_mdtHitResidual", &SA_mdtHitResidual, &b_SA_mdtHitResidual);
   fChain->SetBranchAddress("SA_roadAw", &SA_roadAw, &b_SA_roadAw);
   fChain->SetBranchAddress("SA_roadBw", &SA_roadBw, &b_SA_roadBw);
   fChain->SetBranchAddress("SA_zMin", &SA_zMin, &b_SA_zMin);
   fChain->SetBranchAddress("SA_zMax", &SA_zMax, &b_SA_zMax);
   fChain->SetBranchAddress("SA_rMin", &SA_rMin, &b_SA_rMin);
   fChain->SetBranchAddress("SA_rMax", &SA_rMax, &b_SA_rMax);
   fChain->SetBranchAddress("SA_etaMin", &SA_etaMin, &b_SA_etaMin);
   fChain->SetBranchAddress("SA_etaMax", &SA_etaMax, &b_SA_etaMax);
   fChain->SetBranchAddress("SA_stgcClusterZ", &SA_stgcClusterZ, &b_SA_stgcClusterZ);
   fChain->SetBranchAddress("SA_stgcClusterR", &SA_stgcClusterR, &b_SA_stgcClusterR);
   fChain->SetBranchAddress("SA_stgcClusterEta", &SA_stgcClusterEta, &b_SA_stgcClusterEta);
   fChain->SetBranchAddress("SA_stgcClusterPhi", &SA_stgcClusterPhi, &b_SA_stgcClusterPhi);
   fChain->SetBranchAddress("SA_stgcClusterResidualR", &SA_stgcClusterResidualR, &b_SA_stgcClusterResidualR);
   fChain->SetBranchAddress("SA_stgcClusterResidualPhi", &SA_stgcClusterResidualPhi, &b_SA_stgcClusterResidualPhi);
   fChain->SetBranchAddress("SA_stgcClusterStationEta", &SA_stgcClusterStationEta, &b_SA_stgcClusterStationEta);
   fChain->SetBranchAddress("SA_stgcClusterStationPhi", &SA_stgcClusterStationPhi, &b_SA_stgcClusterStationPhi);
   fChain->SetBranchAddress("SA_stgcClusterStationName", &SA_stgcClusterStationName, &b_SA_stgcClusterStationName);
   fChain->SetBranchAddress("SA_stgcClusterType", &SA_stgcClusterType, &b_SA_stgcClusterType);
   fChain->SetBranchAddress("SA_stgcClusterIsOutlier", &SA_stgcClusterIsOutlier, &b_SA_stgcClusterIsOutlier);
   fChain->SetBranchAddress("SA_stgcClusterLayer", &SA_stgcClusterLayer, &b_SA_stgcClusterLayer);
   fChain->SetBranchAddress("SA_mmClusterZ", &SA_mmClusterZ, &b_SA_mmClusterZ);
   fChain->SetBranchAddress("SA_mmClusterR", &SA_mmClusterR, &b_SA_mmClusterR);
   fChain->SetBranchAddress("SA_mmClusterEta", &SA_mmClusterEta, &b_SA_mmClusterEta);
   fChain->SetBranchAddress("SA_mmClusterPhi", &SA_mmClusterPhi, &b_SA_mmClusterPhi);
   fChain->SetBranchAddress("SA_mmClusterResidualR", &SA_mmClusterResidualR, &b_SA_mmClusterResidualR);
   fChain->SetBranchAddress("SA_mmClusterResidualPhi", &SA_mmClusterResidualPhi, &b_SA_mmClusterResidualPhi);
   fChain->SetBranchAddress("SA_mmClusterStationEta", &SA_mmClusterStationEta, &b_SA_mmClusterStationEta);
   fChain->SetBranchAddress("SA_mmClusterStationPhi", &SA_mmClusterStationPhi, &b_SA_mmClusterStationPhi);
   fChain->SetBranchAddress("SA_mmClusterStationName", &SA_mmClusterStationName, &b_SA_mmClusterStationName);
   fChain->SetBranchAddress("SA_mmClusterIsOutlier", &SA_mmClusterIsOutlier, &b_SA_mmClusterIsOutlier);
   fChain->SetBranchAddress("SA_mmClusterLayer", &SA_mmClusterLayer, &b_SA_mmClusterLayer);
   fChain->SetBranchAddress("CB_pt", &CB_pt, &b_CB_pt);
   fChain->SetBranchAddress("CB_eta", &CB_eta, &b_CB_eta);
   fChain->SetBranchAddress("CB_phi", &CB_phi, &b_CB_phi);
   fChain->SetBranchAddress("CB_idpt", &CB_idpt, &b_CB_idpt);
   fChain->SetBranchAddress("CB_ideta", &CB_ideta, &b_CB_ideta);
   fChain->SetBranchAddress("CB_idphi", &CB_idphi, &b_CB_idphi);
   fChain->SetBranchAddress("CB_roiNumber", &CB_roiNumber, &b_CB_roiNumber);
   fChain->SetBranchAddress("CB_roiSector", &CB_roiSector, &b_CB_roiSector);
   fChain->SetBranchAddress("CBIO_pt", &CBIO_pt, &b_CBIO_pt);
   fChain->SetBranchAddress("CBIO_eta", &CBIO_eta, &b_CBIO_eta);
   fChain->SetBranchAddress("CBIO_phi", &CBIO_phi, &b_CBIO_phi);
   fChain->SetBranchAddress("CBIO_idpt", &CBIO_idpt, &b_CBIO_idpt);
   fChain->SetBranchAddress("CBIO_ideta", &CBIO_ideta, &b_CBIO_ideta);
   fChain->SetBranchAddress("CBIO_idphi", &CBIO_idphi, &b_CBIO_idphi);
   fChain->SetBranchAddress("CBIO_roiNumber", &CBIO_roiNumber, &b_CBIO_roiNumber);
   fChain->SetBranchAddress("CBIO_roiSector", &CBIO_roiSector, &b_CBIO_roiSector);
   fChain->SetBranchAddress("EF_pt", &EF_pt, &b_EF_pt);
   fChain->SetBranchAddress("EF_eta", &EF_eta, &b_EF_eta);
   fChain->SetBranchAddress("EF_phi", &EF_phi, &b_EF_phi);
   Notify();
}

Bool_t BetaStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BetaStudy::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BetaStudy::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BetaStudy_cxx
