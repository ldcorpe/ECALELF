#include <iostream>
#include <iomanip>
#include <cmath>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TMath.h"
#include <string>
#include <cstring>
#include <sstream>
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THStack.h"
#include "TPaveStats.h"


int selectionPlot() 
{

	//Open Map and both data sources.
	TFile *file_Map = TFile::Open("tmp/Map_d1-louieTest.root");
	TFile *file_53X = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2/190645-193621/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2-190645-193621.root");
	TFile *file_70X  = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-22Jan-v1/190645-193621/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012A-22Jan-v1-190645-193621.root");

	//Define respective TTrees
	TTree *chain7 =  (TTree*)file_70X->Get("selected");
	TTree *chain5 =  (TTree*)file_53X->Get("selected"); 
	TTree *tree_Map =  (TTree*)file_Map->Get("Map"); 

	//------VARIABLES-------------//	
	//////////VARIABLES FOR 70X////////////////
	Double_t selection7;
	Double_t total_weight7=1;
	ULong64_t eventNumber7;	
	// for the energy calculation
	Float_t corrEle_7[2]={1,1};	
	Float_t smearEle_7[2]={1,1};	
	// for the weight calculation
	Float_t weight7=1;//pu weight
	Float_t r9weight7[2]={1,1};	
	Float_t ptweight7[2]={1,1};	
	Float_t mcGenWeight7=1;	
	Float_t smearerCat7[2];	
	// Stored variables
	Float_t energyEle7[2];	
	Float_t etaEle7[2];	
	Float_t etaSCEle7[2];
	Float_t phiSCEle7[2];
	Float_t phiEle7[2];	
	Float_t invMass7;	
	Float_t invMass_raw7;	
	Float_t invMass_corr7;	
	Float_t energySCEle7[2];
	Float_t energySCEle_raw7[2];
	Float_t energySCEle_corr7[2];
	Float_t e5x5SCEle7[2];
	Float_t pModeGsfEle7[2];
	Float_t pAtVtxGsfEle7[2];
	Float_t PtEle7[2];	
	Float_t seedEnergySCEle7[2];	
	Int_t chargeEle7[2];	
	Float_t R9Ele7[2];	
	Float_t seedXSCEle7[2];
	Float_t seedYSCEle7[2];
	//Useful variabled
	Int_t eleID7[2];	
	Char_t HLTfire7;	
	Int_t recoFlagsEle7[2]; 
	//////////VARIABLES FOR 53X////////////////
	Double_t selection5;
	Double_t total_weight5=1;
	ULong64_t eventNumber5;	
	// for the energy calculation
	Float_t corrEle_5[2]={1,1};	
	Float_t smearEle_5[2]={1,1};	
	// for the weight calculation
	Float_t weight5=1;//pu weight
	Float_t r9weight5[2]={1,1};	
	Float_t ptweight5[2]={1,1};	
	Float_t mcGenWeight5=1;	
	Float_t smearerCat5[2];	
	// Stored variables
	Float_t energyEle5[2];	
	Float_t etaEle5[2];	
	Float_t etaSCEle5[2];
	Float_t phiEle5[2];	
	Float_t phiSCEle5[2];	
	Float_t invMass5;	
	Float_t invMass_raw5;	
	Float_t invMass_corr5;	
	Float_t energySCEle5[2];
	Float_t energySCEle_raw5[2];
	Float_t energySCEle_corr5[2];
	Float_t e5x5SCEle5[2];
	Float_t pModeGsfEle5[2];
	Float_t pAtVtxGsfEle5[2];
	Float_t PtEle5[2];	
	Float_t seedEnergySCEle5[2];	
	Int_t chargeEle5[2];	
	Float_t R9Ele5[2];	
	Float_t seedXSCEle5[2];
	Float_t seedYSCEle5[2];
	Float_t fbremEle5[2];
	Float_t fbremEle7[2];
	//Useful variabled
	Int_t eleID5[2];	
	Char_t HLTfire5;	
	Int_t recoFlagsEle5[2]; 
	//////////VARIABLES FOR MAP////////////////
	Int_t eleIndex[2];
	Int_t entryNumber;


	//-----------BRANCHES---------//
	//For 70X//
	TBranch *b_eventNumber7;	
	TBranch *b_weight7;	
	TBranch *b_energyEle7;	
	TBranch *b_etaEle7;	
	TBranch *b_etaSCEle7;
	TBranch *b_phiSCEle7;
	TBranch *b_phiEle7;	
	TBranch *b_invMass7;	
	TBranch *b_invMass_raw7;	
	TBranch *b_invMass_corr7;	
	TBranch *b_energySCEle7;
	TBranch *b_energySCEle_raw7;
	TBranch *b_energySCEle_corr7;
	TBranch *b_e5x5SCEle7;
	TBranch *b_pModeGsfEle7;
	TBranch *b_pAtVtxGsfEle7;
	TBranch *b_PtEle7;	
	TBranch *b_seedEnergySCEle7;	
	TBranch *b_chargeEle7;	
	TBranch *b_R9Ele7;	
	TBranch *b_seedXSCEle7;
	TBranch *b_seedYSCEle7;
	TBranch *b_eleID7;	
	TBranch *b_HLTfire7;	
	TBranch *b_recoFlagsEle7;
	TBranch *b_fbremEle5;
	TBranch *b_fbremEle7;
	//For 53X//
	TBranch *b_eventNumber5;	
	TBranch *b_weight5;	
	TBranch *b_energyEle5;	
	TBranch *b_e5x5SCEle5;	
	TBranch *b_etaEle5;	
	TBranch *b_etaSCEle5;
	TBranch *b_phiSCEle5;
	TBranch *b_phiEle5;	
	TBranch *b_invMass5;	
	TBranch *b_invMass_raw5;	
	TBranch *b_invMass_corr5;	
	TBranch *b_energySCEle5;
	TBranch *b_energySCEle_raw5;
	TBranch *b_energySCEle_corr5;
	TBranch *b_pModeGsfEle5;
	TBranch *b_pAtVtxGsfEle5;
	TBranch *b_PtEle5;	
	TBranch *b_seedEnergySCEle5;	
	TBranch *b_chargeEle5;	
	TBranch *b_R9Ele5;	
	TBranch *b_seedXSCEle5;
	TBranch *b_seedYSCEle5;
	TBranch *b_eleID5;	
	TBranch *b_HLTfire5;	
	TBranch *b_recoFlagsEle5; 
	//For MAP//
	TBranch *b_eleIndex;
	TBranch *b_entryNumber;

	//----------- ADDRESSES----------/
	//For53X//
	chain5->SetBranchAddress("eventNumber", &eventNumber5,&b_eventNumber5);
	chain5->SetBranchAddress("etaEle", etaEle5,&b_etaEle5);
	chain5->SetBranchAddress("fbremEle", fbremEle5,&b_fbremEle5);
	chain7->SetBranchAddress("fbremEle", fbremEle7,&b_fbremEle7);
	chain5->SetBranchAddress("etaSCEle", etaSCEle5,&b_etaSCEle5);
	chain5->SetBranchAddress("phiSCEle", phiSCEle5,&b_phiSCEle5);
	chain5->SetBranchAddress("phiEle", phiEle5,&b_phiEle5);
	chain5->SetBranchAddress("invMass_SC", &invMass5, &b_invMass5);
	chain5->SetBranchAddress("invMass_rawSC", &invMass_raw5, &b_invMass_raw5);
	chain5->SetBranchAddress("invMass_SC_corr", &invMass_corr5, &b_invMass_corr5);
	chain5->SetBranchAddress("pModeGsfEle",pModeGsfEle5,&b_pModeGsfEle5);
	chain5->SetBranchAddress("pAtVtxGsfEle",pAtVtxGsfEle5,&b_pAtVtxGsfEle5);
	//	chain5->SetBranchAddress(invMass_var5.c_str(),&invMass5,&b_invMass5);	
	chain5->SetBranchAddress("energySCEle",energySCEle5,&b_energySCEle5);
	chain5->SetBranchAddress("rawEnergySCEle",energySCEle_raw5,&b_energySCEle_raw5);
	chain5->SetBranchAddress("energySCEle_regrCorrSemiParV5_ele",energySCEle_corr5,&b_energySCEle_corr5);
	chain5->SetBranchAddress("e5x5SCEle",e5x5SCEle5,&b_e5x5SCEle5);
	chain5->SetBranchAddress("PtEle",PtEle5,&b_PtEle5);	
	chain5->SetBranchAddress("seedEnergySCEle",seedEnergySCEle5,&b_seedEnergySCEle5);
	chain5->SetBranchAddress("chargeEle",chargeEle5,&b_chargeEle5);	
	chain5->SetBranchAddress("R9Ele",R9Ele5,&b_R9Ele5);	
	chain5->SetBranchAddress("seedXSCEle",seedXSCEle5,&b_seedXSCEle5);
	chain5->SetBranchAddress("seedYSCEle",seedYSCEle5,&b_seedYSCEle5);
	chain5->SetBranchAddress("eleID",eleID5,&b_eleID5);
	chain5->SetBranchAddress("HLTfire",&HLTfire5,&b_HLTfire5);
	chain5->SetBranchAddress("recoFlagsEle",recoFlagsEle5,&b_recoFlagsEle5);
	TBranch *b_scaleEle5;
	if(chain5->GetBranch("scaleEle")!=NULL){
		chain5->SetBranchAddress("scaleEle", corrEle_5,&b_scaleEle5);
	}
	TBranch *b_smearEle5;
	if(chain5->GetBranch("smearEle")!=NULL){
		chain5->SetBranchAddress("smearEle", smearEle_5,&b_smearEle5);
	}
	TBranch *b_puWeight5;
	if(chain5->GetBranch("puWeight")!=NULL){
		chain5->SetBranchAddress("puWeight", &weight5,&b_puWeight5);
	}
	TBranch *b_r9Weight5;
	if(chain5->GetBranch("r9Weight")!=NULL){
		chain5->SetBranchAddress("r9Weight", r9weight5,&b_r9Weight5);
	}
	TBranch *b_ptWeight5;
	if(chain5->GetBranch("ptWeight")!=NULL){
		chain5->SetBranchAddress("ptWeight", ptweight5,&b_ptWeight5);
	}
	TBranch *b_mcGenWeight5;
	if(chain5->GetBranch("mcGenWeight")!=NULL){
		chain5->SetBranchAddress("mcGenWeight", &mcGenWeight5,&b_mcGenWeight5);
	}
	TBranch *b_smearerCat5;
	if(chain5->GetBranch("smearerCat")!=NULL){
		chain5->SetBranchAddress("smearerCat", smearerCat5,&b_smearerCat5);
	}

	//For73X//
	chain7->SetBranchAddress("eventNumber", &eventNumber7,&b_eventNumber7);
	chain7->SetBranchAddress("etaEle", etaEle7,&b_etaEle7);
	chain7->SetBranchAddress("etaSCEle", etaSCEle7,&b_etaSCEle7);
	chain7->SetBranchAddress("phiSCEle", phiSCEle7,&b_phiSCEle7);
	chain7->SetBranchAddress("phiEle", phiEle7,&b_phiEle7);
	chain7->SetBranchAddress("pModeGsfEle",pModeGsfEle7,&b_pModeGsfEle7);
	chain7->SetBranchAddress("invMass_SC", &invMass7, &b_invMass7);
	chain7->SetBranchAddress("invMass_rawSC", &invMass_raw7, &b_invMass_raw7);
	chain7->SetBranchAddress("invMass_SC_corr", &invMass_corr7, &b_invMass_corr7);
	chain7->SetBranchAddress("pAtVtxGsfEle",pAtVtxGsfEle7,&b_pAtVtxGsfEle7);
	//	chain7->SetBranchAddress(invMass_var7.c_str(),&invMass7,&b_invMass7);	
	chain7->SetBranchAddress("energySCEle",energySCEle7,&b_energySCEle7);
	chain7->SetBranchAddress("rawEnergySCEle",energySCEle_raw7,&b_energySCEle_raw7);
	chain7->SetBranchAddress("energySCEle_corr",energySCEle_corr7,&b_energySCEle_corr7);
	chain7->SetBranchAddress("e5x5SCEle",e5x5SCEle7,&b_e5x5SCEle7);
	chain7->SetBranchAddress("PtEle",PtEle7,&b_PtEle7);	
	chain7->SetBranchAddress("seedEnergySCEle",seedEnergySCEle7,&b_seedEnergySCEle7);
	chain7->SetBranchAddress("chargeEle",chargeEle7,&b_chargeEle7);	
	chain7->SetBranchAddress("R9Ele",R9Ele7,&b_R9Ele7);	
	chain7->SetBranchAddress("seedXSCEle",seedXSCEle7,&b_seedXSCEle7);
	chain7->SetBranchAddress("seedYSCEle",seedYSCEle7,&b_seedYSCEle7);
	chain7->SetBranchAddress("eleID",eleID7,&b_eleID7);
	chain7->SetBranchAddress("HLTfire",&HLTfire7,&b_HLTfire7);
	chain7->SetBranchAddress("recoFlagsEle",recoFlagsEle7,&b_recoFlagsEle7);
	TBranch *b_scaleEle7;
	if(chain7->GetBranch("scaleEle")!=NULL){
		chain7->SetBranchAddress("scaleEle", corrEle_7,&b_scaleEle7);
	}
	TBranch *b_smearEle7;
	if(chain7->GetBranch("smearEle")!=NULL){
		chain7->SetBranchAddress("smearEle", smearEle_7,&b_smearEle7);
	}
	TBranch *b_puWeight7;
	if(chain7->GetBranch("puWeight")!=NULL){
		chain7->SetBranchAddress("puWeight", &weight7,&b_puWeight7);
	}
	TBranch *b_r9Weight7;
	if(chain7->GetBranch("r9Weight")!=NULL){
		chain7->SetBranchAddress("r9Weight", r9weight7,&b_r9Weight7);
	}
	TBranch *b_ptWeight7;
	if(chain7->GetBranch("ptWeight")!=NULL){
		chain7->SetBranchAddress("ptWeight", ptweight7,&b_ptWeight7);
	}
	TBranch *b_mcGenWeight7;
	if(chain7->GetBranch("mcGenWeight")!=NULL){
		chain7->SetBranchAddress("mcGenWeight", &mcGenWeight7,&b_mcGenWeight7);
	}
	TBranch *b_smearerCat7;
	if(chain7->GetBranch("smearerCat")!=NULL){
		chain7->SetBranchAddress("smearerCat", smearerCat7,&b_smearerCat7);
	}
	//MAP//
	tree_Map->SetBranchAddress("eleIndex",&eleIndex);
	tree_Map->SetBranchAddress("entryNumber2",&entryNumber);



	//Define Histograms

	//Simple quantitiy comaprisons
	TH1F *E_hist_EE = new TH1F("energySCEle_EE","energySCEle_EE (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *E3x3_hist_EE = new TH1F("energySCEle3x3_EE","energySCEle3x3_EE (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *rawE_hist_EE = new TH1F("energySCEle_raw_EE","energySCEle_raw_EE (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *corrE_hist_EE = new TH1F("energySCEle_corr_EE","energySCEle_corr_EE (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *E_hist_EB = new TH1F("energySCEle_EB","energySCEle_EB (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *E3x3_hist_EB = new TH1F("energySCEle3x3_EB","energySCEle3x3_EB (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *rawE_hist_EB = new TH1F("energySCEle_raw_EB","energySCEle_raw_EB (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *corrE_hist_EB = new TH1F("energySCEle_corr_EB","energySCEle_corr_EB (_7-_5)/_5",1000,-0.2,0.2);
	TH1F *E5x5_hist = new TH1F("e5x5SCEle","e5x5SCEle (_7-_5)/_5",10000,-0.2,0.2);
	TH1F *R9_hist_5 = new TH1F("R9Ele5","R9Ele5",1000,0,1.2);
	TH1F *R9_hist_7 = new TH1F("R9Ele7","R9Ele7",1000,0,1.2);
	TH1F *R9_hist = new TH1F("R9Ele","R9Ele (_7-_5)",1000,-0.4,0.4);
	TH1F *etaSC_hist = new TH1F("etaSCEle","etaSCEle (_7-_5)",1000,-0.01,0.01);
	TH1F *phiSC_hist = new TH1F("phiSCEle","phiSCEle (_7-_5)",1000,-0.01,0.01);
	TH1F *eta_hist = new TH1F("etaEle","etaEle (_7-_5)",1000,-0.01,0.01);
	TH1F *phi_hist = new TH1F("phiEle","phiEle (_7-_5)",1000,-0.01,0.01);
	TH1F *pModeGsf_hist = new TH1F("pModeGSf","pmode (_7-_5)/_5",1000,-0.02,0.02);
	TH1F *EbyP_hist = new TH1F("EbyP","E/p (_7-_5)",1000,-0.3,0.3);
	TH1F *EbyP_hist_5 = new TH1F("EbyP5","E/p (5)",1000,0.7,1.3);
	TH1F *EbyP_hist_7 = new TH1F("EbyP7","E/p (7)",1000,0.7,1.3);
	TProfile *fbrem = new TProfile("fbrem","fbrem _7/_5",500,-2.6,2.6);
	
	//Z->ee distributions
	TH1F *EEinvMassGold7 = new TH1F("EEinvMassGold7","EEinvMassGold",100,40,140);
	TH1F *EEinvMassBad7 = new TH1F( "EinvMassBad7" ,"EEinvMassBad" ,100,40,140);
	TH1F *EBinvMassGold7 = new TH1F("BinvMassGold7","EBinvMassGold",100,40,140);
	TH1F *EBinvMassBad7 = new TH1F( "BinvMassBad7 ","EBinvMassBad ",100,40,140);
	TH1F *EEinvMassGold5 = new TH1F("EEinvMassGold5","EEinvMassGold",100,40,140);
	TH1F *EEinvMassBad5 = new TH1F( "EinvMassBad5" ,"EEinvMassBad" ,100,40,140);
	TH1F *EBinvMassGold5 = new TH1F("BinvMassGold5","EBinvMassGold",100,40,140);
	TH1F *EBinvMassBad5 = new TH1F( "BinvMassBad5 ","EBinvMassBad ",100,40,140);

	TH1F *EEinvMassGold_raw7 = new TH1F("EEinvMassGold_raw7","EEinvMassGold_raw",100,40,140);
	TH1F *EEinvMassBad_raw7 = new TH1F( "EinvMassBad_raw7" ,"EEinvMassBad_raw" ,100,40,140);
	TH1F *EBinvMassGold_raw7 = new TH1F("BinvMassGold_raw7","EBinvMassGold_raw",100,40,140);
	TH1F *EBinvMassBad_raw7 = new TH1F( "BinvMassBad_raw7 ","EBinvMassBad_raw ",100,40,140);
	TH1F *EEinvMassGold_raw5 = new TH1F("EEinvMassGold_raw5","EEinvMassGold_raw",100,40,140);
	TH1F *EEinvMassBad_raw5 = new TH1F( "EinvMassBad_raw5" ,"EEinvMassBad_raw" ,100,40,140);
	TH1F *EBinvMassGold_raw5 = new TH1F("BinvMassGold_raw5","EBinvMassGold_raw",100,40,140);
	TH1F *EBinvMassBad_raw5 = new TH1F( "BinvMassBad_raw5 ","EBinvMassBad_raw ",100,40,140);

	TH1F *EEinvMassGold_corr7 = new TH1F("EEinvMassGold_corr7","EEinvMassGold_corr",100,40,140);
	TH1F *EEinvMassBad_corr7 = new TH1F( "EinvMassBad_corr7" ,"EEinvMassBad_corr" ,100,40,140);
	TH1F *EBinvMassGold_corr7 = new TH1F("BinvMassGold_corr7","EBinvMassGold_corr",100,40,140);
	TH1F *EBinvMassBad_corr7 = new TH1F( "BinvMassBad_corr7 ","EBinvMassBad_corr ",100,40,140);
	TH1F *EEinvMassGold_corr5 = new TH1F("EEinvMassGold_corr5","EEinvMassGold_corr",100,40,140);
	TH1F *EEinvMassBad_corr5 = new TH1F( "EinvMassBad_corr5" ,"EEinvMassBad_corr" ,100,40,140);
	TH1F *EBinvMassGold_corr5 = new TH1F("BinvMassGold_corr5","EBinvMassGold_corr",100,40,140);
	TH1F *EBinvMassBad_corr5 = new TH1F( "BinvMassBad_corr5 ","EBinvMassBad_corr ",100,40,140);

//iEta and iPhi Plots
	TProfile2D *r9 = new TProfile2D("mean_deltar9","mean deltar9",360,1,360,170,-85,85,"s");
	TProfile2D *e_sc = new TProfile2D("mean_change_e_sc","mean change e_sc",360,1,360,170,-85,85,"s");
	TProfile2D *e_raw = new TProfile2D("mean_change_e_raw","mean change e_raw",360,1,360,170,-85,85,"s");
	TProfile2D *r9_alternative = new TProfile2D("mean_deltar9_alt","mean deltar9",360,-3.14,3.14,170,-1.5,1.5,"s");
	TProfile2D *e_sc_alternative = new TProfile2D("mean_change_e_sc_alt","mean change e_sc",360,-3.14,3.14,170,-1.5,1.5,"s");


	//Sumw2 settings for standard histos.
	E_hist_EB->Sumw2(); 
	E3x3_hist_EB->Sumw2(); 
	corrE_hist_EB->Sumw2();
	rawE_hist_EB->Sumw2();
	E_hist_EE->Sumw2(); 
	E3x3_hist_EE->Sumw2(); 
	corrE_hist_EE->Sumw2();
	rawE_hist_EE->Sumw2();
	E5x5_hist->Sumw2(); 
	R9_hist->Sumw2(); 
	etaSC_hist->Sumw2(); 
	phiSC_hist->Sumw2(); 
	eta_hist->Sumw2(); 
	phi_hist->Sumw2(); 
	pModeGsf_hist ->Sumw2();

EEinvMassGold7 ->Sumw2();
EEinvMassBad7 ->Sumw2();
EBinvMassGold7->Sumw2();
EBinvMassBad7->Sumw2();
EEinvMassGold5->Sumw2();
EEinvMassBad5 ->Sumw2();
EBinvMassGold5->Sumw2();
EBinvMassBad5 ->Sumw2();
                   
EEinvMassGold_raw7->Sumw2();  
EEinvMassBad_raw7 ->Sumw2();
EBinvMassGold_raw7 ->Sumw2();
EBinvMassBad_raw7 ->Sumw2();
EEinvMassGold_raw5 ->Sumw2();
EEinvMassBad_raw5 ->Sumw2();
EBinvMassGold_raw5 ->Sumw2();
EBinvMassBad_raw5 ->Sumw2();
                   
EEinvMassGold_corr7->Sumw2();
EEinvMassBad_corr7 ->Sumw2();
EBinvMassGold_corr7->Sumw2();
EBinvMassBad_corr7 ->Sumw2();
EEinvMassGold_corr5->Sumw2();
EEinvMassBad_corr5 ->Sumw2();
EBinvMassGold_corr5->Sumw2();
EBinvMassBad_corr5 ->Sumw2();

	
	//Cosmetics	
	EEinvMassGold5->SetLineColor(kRed);
	EEinvMassBad5 ->SetLineColor(kRed);
	EBinvMassGold5->SetLineColor(kRed);
	EBinvMassBad5 ->SetLineColor(kRed);
	EEinvMassGold_corr5->SetLineColor(kRed);
	EEinvMassBad_corr5 ->SetLineColor(kRed);
	EBinvMassGold_corr5->SetLineColor(kRed);
	EBinvMassBad_corr5 ->SetLineColor(kRed);
	EEinvMassGold_raw5->SetLineColor(kRed);
	EEinvMassBad_raw5 ->SetLineColor(kRed);
	EBinvMassGold_raw5->SetLineColor(kRed);
	EBinvMassBad_raw5 ->SetLineColor(kRed);
	R9_hist_5->SetLineColor(kRed);
	EbyP_hist_5->SetLineColor(kRed);
	
	
	//Loops over first TTree	
	for(Int_t loop = 0; loop < 

			chain7->GetEntries()
			//20000	
			; loop++)
	
	{
		// Progress tracker
		if (loop%10000==0) 
		{ 
			cout << loop << " / " << chain7->GetEntries() << endl;
		}

		//Get entry in firs TTree and corresponding Map entry	
		tree_Map->GetEntry(loop);
		chain7->GetEntry(loop);

		//Skip Entries where there is no match
		if (entryNumber==-1) continue;

		//Get corresponding entry in TTree 2  using information from Map
		chain5->GetEntry(entryNumber);

		//set selection value for events in 5 and 7
		selection5=((eleID5[0] & 2)==2)*((eleID5[1] & 2)==2)*(HLTfire5==1)*(recoFlagsEle5[0] > 1)*(recoFlagsEle5[1] > 1)*(PtEle5[0]>20)*(PtEle5[1]>20);
		total_weight5=1;
		total_weight5*=weight5*r9weight5[0]*r9weight5[1]*ptweight5[0]*ptweight5[1];//*mcGenWeight;
		selection5*=total_weight5;

		selection7=((eleID7[0] & 2)==2)*((eleID7[1] & 2)==2)*(HLTfire7==1)*(recoFlagsEle7[0] > 1)*(recoFlagsEle7[1] > 1)*(PtEle7[0]>20)*(PtEle7[1]>20);
		total_weight7=1;
		total_weight7*=weight7*r9weight7[0]*r9weight7[1]*ptweight7[0]*ptweight7[1];//*mcGenWeight;
		selection7*=total_weight7;

		//Define iEta and iPhi
		Float_t iPhi5[2] = { floor(phiSCEle5[0]*57.30+0.5)+181, floor(phiSCEle5[1]*57.30+0.5)+181 };
		Float_t iEta5[2] = { ((etaSCEle5[0]>0)- (etaSCEle5[0]<0))*(floor(fabs(etaSCEle5[0]*57.47)+0.5)+1), ((etaSCEle5[1]>0) - (etaSCEle5[1]<0))*(floor(fabs(etaSCEle5[1]*57.46)+0.5)+1 )};

		//Fill Z->ee histograms
		if((fabs(etaSCEle5[0])>1.566 ) && (fabs(etaSCEle5[1])>1.566))
		{

			if (R9Ele5[0]>0.94 && R9Ele5[1]>0.94)
			{
				EEinvMassGold5->Fill(invMass5,selection5);
				EEinvMassGold_raw5->Fill(invMass_raw5,selection5);
				EEinvMassGold_corr5->Fill(invMass_corr5,selection5);
			}	

			if (R9Ele5[0]<0.94 && R9Ele5[1]<0.94)
			{
				EEinvMassBad5->Fill(invMass5,selection5);
				EEinvMassBad_raw5->Fill(invMass_raw5,selection5);
				EEinvMassBad_corr5->Fill(invMass_corr5,selection5);
			}
		}
		if((fabs(etaSCEle5[0])<1.444) && (fabs(etaSCEle5[1])<1.4442)) 
		{

			if (R9Ele5[0]>0.94 && R9Ele5[1]>0.94)
			{
				EBinvMassGold5->Fill(invMass5,selection5);
				EBinvMassGold_raw5->Fill(invMass_raw5,selection5);
				EBinvMassGold_corr5->Fill(invMass_corr5,selection5);
			}	
			if (R9Ele5[0]<0.94 && R9Ele5[1]<0.94)
			{
				EBinvMassBad5->Fill(invMass5,selection5);
				EBinvMassBad_raw5->Fill(invMass_raw5,selection5);
				EBinvMassBad_corr5->Fill(invMass_corr5,selection5);
			}
		}

		if((fabs(etaSCEle7[0])>1.566 ) && (fabs(etaSCEle7[1])>1.566))
		{

			if (R9Ele7[0]>0.94 && R9Ele7[1]>0.94)
			{
				EEinvMassGold7->Fill(invMass7,selection5);
				EEinvMassGold_raw7->Fill(invMass_raw7,selection5);
				EEinvMassGold_corr7->Fill(invMass_corr7,selection5);
			}	
			if (R9Ele7[0]<0.94 && R9Ele7[1]<0.94)
			{
				EEinvMassBad7->Fill(invMass7,selection5);
				EEinvMassBad_raw7->Fill(invMass_raw7,selection5);
				EEinvMassBad_corr7->Fill(invMass_corr7,selection5);
			}
		}
		if((fabs(etaSCEle7[0])<1.444) && (fabs(etaSCEle7[1])<1.4442)) 
		{

			if (R9Ele7[0]>0.94 && R9Ele7[1]>0.94)
			{
				EBinvMassGold7->Fill(invMass7,selection5);
				EBinvMassGold_raw7->Fill(invMass_raw7,selection5);
				EBinvMassGold_corr7->Fill(invMass_corr7,selection5);
			}	
			if (R9Ele7[0]<0.94 && R9Ele7[1]<0.94)
			{
				EBinvMassBad7->Fill(invMass7,selection5);
				EBinvMassBad_raw7->Fill(invMass_raw7,selection5);
				EBinvMassBad_corr7->Fill(invMass_corr7,selection5);
			}
		}
		
		//Fill iEta and iPhi maps
		if (eleIndex[0]==1)
		{

			if(fabs(etaSCEle5[0])<1.4442)
			{
				r9->Fill(iPhi5[0],iEta5[0],R9Ele7[0]-R9Ele5[0]);
				e_sc->Fill(iPhi5[0],iEta5[0],(energySCEle7[0]-energySCEle5[0])/energySCEle5[0]);
				e_raw->Fill(iPhi5[0],iEta5[0],(energySCEle_raw7[0]-energySCEle_raw5[0])/energySCEle_raw5[0]);
				r9_alternative->Fill(phiSCEle5[0],etaSCEle5[0],R9Ele7[0]-R9Ele5[0]);
				e_sc_alternative->Fill(phiSCEle5[0],etaSCEle5[0],(energySCEle7[0]-energySCEle5[0])/energySCEle5[0]);
			}

			if(fabs(etaSCEle5[1])<1.4442)
			{
				r9->Fill(iPhi5[1],iEta5[1],R9Ele7[1]-R9Ele5[1]);
				e_sc->Fill(iPhi5[1],iEta5[1],(energySCEle7[1]-energySCEle5[1])/energySCEle5[1]);
				e_raw->Fill(iPhi5[1],iEta5[1],(energySCEle_raw7[1]-energySCEle_raw5[1])/energySCEle_raw5[1]);
				r9_alternative->Fill(phiSCEle5[1],etaSCEle5[1],R9Ele7[1]-R9Ele5[1]);
				e_sc_alternative->Fill(phiSCEle5[1],etaSCEle5[1],(energySCEle7[1]-energySCEle5[1])/energySCEle5[1]);
			}
		}


		if (eleIndex[0]==2)
		{
			if(fabs(etaSCEle5[0])<1.4442)
			{
				r9->Fill(iPhi5[0],iEta5[0],R9Ele7[1]-R9Ele5[0]);
				e_sc->Fill(iPhi5[0],iEta5[0],(energySCEle7[1]-energySCEle5[0])/energySCEle5[0]);
				e_raw->Fill(iPhi5[0],iEta5[0],(energySCEle_raw7[1]-energySCEle_raw5[0])/energySCEle_raw5[0]);
				r9_alternative->Fill(phiSCEle5[0],etaSCEle5[0],R9Ele7[1]-R9Ele5[0]);
				e_sc_alternative->Fill(phiSCEle5[0],etaSCEle5[0],(energySCEle7[1]-energySCEle5[0])/energySCEle5[0]);
			}

			if(fabs(etaSCEle5[1])<1.4442)
			{
				r9->Fill(iPhi5[1],iEta5[1],R9Ele7[0]-R9Ele5[1]);
				e_sc->Fill(iPhi5[1],iEta5[1],(energySCEle7[0]-energySCEle5[1])/energySCEle5[1]);
				e_raw->Fill(iPhi5[1],iEta5[1],(energySCEle_raw7[0]-energySCEle_raw5[1])/energySCEle_raw5[1]);
				r9_alternative->Fill(phiSCEle5[0],etaSCEle5[0],R9Ele7[1]-R9Ele5[0]);
				e_sc_alternative->Fill(phiSCEle5[0],etaSCEle5[0],(energySCEle7[1]-energySCEle5[0])/energySCEle5[0]);
			}
		}


		//Fill basic histograms
		if (eleIndex[0]==1)
		{
			if(etaSCEle5[0]>1.5)
			{
			E3x3_hist_EE->Fill((energySCEle7[0]*R9Ele7[0]-energySCEle5[0]*R9Ele5[0])/energySCEle5[0]*R9Ele5[0],selection5);
			E_hist_EE->Fill((energySCEle7[0]-energySCEle5[0])/energySCEle5[0],selection5);
			rawE_hist_EE->Fill((energySCEle_raw7[0]-energySCEle_raw5[0])/energySCEle_raw5[0],selection5);
			corrE_hist_EE->Fill((energySCEle_corr7[0]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
			}

			if(etaSCEle5[1]>1.5)
			{
			E3x3_hist_EE->Fill((energySCEle7[1]*R9Ele7[1]-energySCEle5[1]*R9Ele5[1])/energySCEle5[1]*R9Ele5[1],selection5);
			E_hist_EE->Fill((energySCEle7[1]-energySCEle5[1])/energySCEle5[1],selection5);
			rawE_hist_EE->Fill((energySCEle_raw7[1]-energySCEle_raw5[1])/energySCEle_raw5[1],selection5);
			corrE_hist_EE->Fill((energySCEle_corr7[1]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
			}

			if(etaSCEle5[0]<1.5)
			{
			E3x3_hist_EB->Fill((energySCEle7[0]*R9Ele7[0]-energySCEle5[0]*R9Ele5[0])/energySCEle5[0]*R9Ele5[0],selection5);
			E_hist_EB->Fill((energySCEle7[0]-energySCEle5[0])/energySCEle5[0],selection5);
			rawE_hist_EB->Fill((energySCEle_raw7[0]-energySCEle_raw5[0])/energySCEle_raw5[0],selection5);
			corrE_hist_EB->Fill((energySCEle_corr7[0]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
			}

			if(etaSCEle5[1]<1.5)
			{
			E3x3_hist_EB->Fill((energySCEle7[1]*R9Ele7[1]-energySCEle5[1]*R9Ele5[1])/energySCEle5[1]*R9Ele5[1],selection5);
			E_hist_EB->Fill((energySCEle7[1]-energySCEle5[1])/energySCEle5[1],selection5);
			rawE_hist_EB->Fill((energySCEle_raw7[1]-energySCEle_raw5[1])/energySCEle_raw5[1],selection5);
			corrE_hist_EB->Fill((energySCEle_corr7[1]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
			}


			E5x5_hist->Fill((e5x5SCEle7[0]-e5x5SCEle5[0])/e5x5SCEle5[0],selection5);
			E5x5_hist->Fill((e5x5SCEle7[1]-e5x5SCEle5[1])/e5x5SCEle5[1],selection5);
			R9_hist_7->Fill((R9Ele7[0]),selection5);
			R9_hist_7->Fill((R9Ele7[1]),selection5);
			R9_hist_5->Fill((R9Ele5[0]),selection5);
			R9_hist_5->Fill((R9Ele5[1]),selection5);
			R9_hist->Fill((R9Ele7[0]-R9Ele5[0]),selection5);
			R9_hist->Fill((R9Ele7[1]-R9Ele5[1]),selection5);
			etaSC_hist->Fill((etaSCEle7[0]-etaSCEle5[0]),selection5);
			etaSC_hist->Fill((etaSCEle7[1]-etaSCEle5[1]),selection5);
			phiSC_hist->Fill((phiSCEle7[0]-phiSCEle5[0]),selection5);
			phiSC_hist->Fill((phiSCEle7[1]-phiSCEle5[1]),selection5);
			eta_hist->Fill((etaEle7[0]-etaEle5[0]),selection5);
			eta_hist->Fill((etaEle7[1]-etaEle5[1]),selection5);
			phi_hist->Fill((phiEle7[0]-phiEle5[0]),selection5);
			phi_hist->Fill((phiEle7[1]-phiEle5[1]),selection5);
			pModeGsf_hist->Fill((pModeGsfEle7[0]-pModeGsfEle5[0])/pModeGsfEle5[0],selection5);
			pModeGsf_hist->Fill((pModeGsfEle7[1]-pModeGsfEle5[1])/pModeGsfEle5[1],selection5);
			EbyP_hist->Fill((energySCEle5[0]/(PtEle5[0]*std::cosh(etaEle5[0])))-(energySCEle7[0]/(PtEle7[0]*std::cosh(etaEle7[0]))),selection5);
			EbyP_hist->Fill((energySCEle5[1]/(PtEle5[1]*std::cosh(etaEle5[1])))-(energySCEle7[1]/(PtEle7[1]*std::cosh(etaEle7[1]))),selection5);
			EbyP_hist_5->Fill((energySCEle5[0]/(PtEle5[0]*std::cosh(etaEle5[0]))),selection5);
			EbyP_hist_7->Fill((energySCEle7[0]/(PtEle7[0]*std::cosh(etaEle7[0]))),selection5);
			EbyP_hist_5->Fill((energySCEle5[1]/(PtEle5[1]*std::cosh(etaEle5[1]))),selection5);
			EbyP_hist_7->Fill((energySCEle7[1]/(PtEle7[1]*std::cosh(etaEle7[1]))),selection5);
			fbrem->Fill(etaSCEle5[0],fbremEle7[0]/fbremEle5[0],selection5);
			fbrem->Fill(etaSCEle5[1],fbremEle7[1]/fbremEle5[1],selection5);
		}

		if (eleIndex[0]==2)
		{

			if(etaSCEle5[0]>1.5)
			{
			E3x3_hist_EE->Fill((energySCEle7[1]*R9Ele7[1]-energySCEle5[0]*R9Ele5[0])/energySCEle5[0]*R9Ele5[0],selection5);
			E_hist_EE->Fill((energySCEle7[1]-energySCEle5[0])/energySCEle5[0],selection5);
			rawE_hist_EE->Fill((energySCEle_raw7[1]-energySCEle_raw5[0])/energySCEle_raw5[0],selection5);
			corrE_hist_EE->Fill((energySCEle_corr7[1]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
			}

			if(etaSCEle5[1]>1.5)
			{
			E3x3_hist_EE->Fill((energySCEle7[0]*R9Ele7[0]-energySCEle5[1]*R9Ele5[1])/energySCEle5[1]*R9Ele5[1],selection5);
			E_hist_EE->Fill((energySCEle7[0]-energySCEle5[1])/energySCEle5[1],selection5);
			rawE_hist_EE->Fill((energySCEle_raw7[0]-energySCEle_raw5[1])/energySCEle_raw5[1],selection5);
			corrE_hist_EE->Fill((energySCEle_corr7[0]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
			}

			if(etaSCEle5[0]<1.5)
			{
			E3x3_hist_EB->Fill((energySCEle7[1]*R9Ele7[1]-energySCEle5[0]*R9Ele5[0])/energySCEle5[0]*R9Ele5[0],selection5);
			E_hist_EB->Fill((energySCEle7[1]-energySCEle5[0])/energySCEle5[0],selection5);
			rawE_hist_EB->Fill((energySCEle_raw7[1]-energySCEle_raw5[0])/energySCEle_raw5[0],selection5);
			corrE_hist_EB->Fill((energySCEle_corr7[1]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
			}

			if(etaSCEle5[1]<1.5)
			{
			E3x3_hist_EB->Fill((energySCEle7[0]*R9Ele7[0]-energySCEle5[1]*R9Ele5[1])/energySCEle5[1]*R9Ele5[1],selection5);
			E_hist_EB->Fill((energySCEle7[0]-energySCEle5[1])/energySCEle5[1],selection5);
			rawE_hist_EB->Fill((energySCEle_raw7[0]-energySCEle_raw5[1])/energySCEle_raw5[1],selection5);
			corrE_hist_EB->Fill((energySCEle_corr7[0]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
			}

			R9_hist_7->Fill((R9Ele7[0]),selection5);
			R9_hist_7->Fill((R9Ele7[1]),selection5);
			R9_hist_5->Fill((R9Ele5[0]),selection5);
			R9_hist_5->Fill((R9Ele5[1]),selection5);
			R9_hist->Fill((R9Ele7[1]-R9Ele5[0]),selection5);
			R9_hist->Fill((R9Ele7[0]-R9Ele5[1]),selection5);
			etaSC_hist->Fill((etaSCEle7[1]-etaSCEle5[0]),selection5);
			etaSC_hist->Fill((etaSCEle7[0]-etaSCEle5[1]),selection5);
			phiSC_hist->Fill((phiSCEle7[1]-phiSCEle5[0]),selection5);
			phiSC_hist->Fill((phiSCEle7[0]-phiSCEle5[1]),selection5);
			eta_hist->Fill((etaEle7[1]-etaEle5[0]),selection5);
			eta_hist->Fill((etaEle7[0]-etaEle5[1]),selection5);
			phi_hist->Fill((phiEle7[1]-phiEle5[0]),selection5);
			phi_hist->Fill((phiEle7[0]-phiEle5[1]),selection5);
			pModeGsf_hist->Fill((pModeGsfEle7[1]-pModeGsfEle5[0])/pModeGsfEle5[0],selection5);
			pModeGsf_hist->Fill((pModeGsfEle7[0]-pModeGsfEle5[1])/pModeGsfEle5[1],selection5);
			EbyP_hist_5->Fill((energySCEle5[0]/(PtEle5[0]*std::cosh(etaEle5[0]))),selection5);
			EbyP_hist_7->Fill((energySCEle7[0]/(PtEle7[0]*std::cosh(etaEle7[0]))),selection5);
			EbyP_hist_5->Fill((energySCEle5[1]/(PtEle5[1]*std::cosh(etaEle5[1]))),selection5);
			EbyP_hist_7->Fill((energySCEle7[1]/(PtEle7[1]*std::cosh(etaEle7[1]))),selection5);
			EbyP_hist->Fill((energySCEle5[0]/(PtEle5[0]*std::cosh(etaEle5[0])))-(energySCEle7[1]/(PtEle7[1]*std::cosh(etaEle7[1]))),selection5);
			EbyP_hist->Fill((energySCEle5[1]/(PtEle5[1]*std::cosh(etaEle5[1])))-(energySCEle7[0]/(PtEle7[0]*std::cosh(etaEle7[0]))),selection5);
			fbrem->Fill(etaSCEle5[0],fbremEle7[1]/fbremEle5[0],selection5);
			fbrem->Fill(etaSCEle5[1],fbremEle7[0]/fbremEle5[1],selection5);
		}
	//end of loop	
	}
	//Draw basic Histos
	TCanvas *t = new TCanvas("t","t",600,600);
	corrE_hist_EE->Draw();
	t->SaveAs("selectionPlots/corrEnergySC_EE.pdf");
	rawE_hist_EE->Draw();
	t->SaveAs("selectionPlots/rawEnergySC_EE.pdf");
	E_hist_EE->Draw();
	t->SaveAs("selectionPlots/energySC_EE.pdf");
	E3x3_hist_EE->Draw();
	t->SaveAs("selectionPlots/e3x3SC_EE.pdf");
	corrE_hist_EB->Draw();
	t->SaveAs("selectionPlots/corrEnergySC_EB.pdf");
	rawE_hist_EB->Draw();
	t->SaveAs("selectionPlots/rawEnergySC_EB.pdf");
	E_hist_EB->Draw();
	t->SaveAs("selectionPlots/energySC_EB.pdf");
	E3x3_hist_EB->Draw();
	t->SaveAs("selectionPlots/e3x3SC_EB.pdf");
	etaSC_hist->Draw();
	t->SaveAs("selectionPlots/etaSC.pdf");
	phiSC_hist->Draw();
	t->SaveAs("selectionPlots/phiSC.pdf");
	eta_hist->Draw();
	t->SaveAs("selectionPlots/eta.pdf");
	phi_hist->Draw();
	t->SaveAs("selectionPlots/phi.pdf");
	pModeGsf_hist->Draw();
	t->SaveAs("selectionPlots/pModeGsf.pdf");
	R9_hist->Draw();
	t->SaveAs("selectionPlots/R9.pdf");
	R9_hist_5->Draw();
	R9_hist_7->Draw("same");
	t->SaveAs("selectionPlots/R9_dist.pdf");
	EbyP_hist->Draw();
	t->SaveAs("selectionPlots/EbyP.pdf");
	EbyP_hist_5->Draw();
	EbyP_hist_7->Draw("same");
	t->SaveAs("selectionPlots/EbyP_dist.pdf");
//	fbrem->SetMimimum(0.8);
//	fbrem->SetMaximum(1.2);
	fbrem->Draw();
	t->SaveAs("selectionPlots/fbrem2.pdf");
	
	//Draw iEta and iPhi plots
gStyle->SetOptStat(0);
	TCanvas *t2 = new TCanvas("t2","t2",600,600);
	t2->Divide(1,2);
	t2->cd(1);
	r9->SetMinimum(-0.05);
	r9->SetMaximum(0.05);
	r9->Draw("colz");
	t2->cd(2);
	TH2D *r9_error = r9->ProjectionXY("rms","c=e");
	r9_error->SetMinimum(-0.03);
	r9_error->SetMaximum(0.03);
	r9_error->Draw("colz");
	t2->SaveAs("selectionPlots/iEtaPhiR9.pdf");
	t2->cd(1);
	e_sc->SetMinimum(-0.05);
	e_sc->SetMaximum(0.05);
	e_sc->Draw("colz");
	t2->cd(2);
	TH2D *e_sc_error = e_sc->ProjectionXY("esc_rms","c=e");
	e_sc_error->SetMinimum(-0.03);
	e_sc_error->SetMaximum(0.03);
	e_sc_error->Draw("colz");
	t2->SaveAs("selectionPlots/iEtaPhiEsc.pdf");
	t2->cd(1);
	e_raw->SetMinimum(-0.05);
	e_raw->SetMaximum(0.05);
	e_raw->Draw("colz");
	t2->cd(2);
	TH2D *e_raw_error = e_raw->ProjectionXY("e_raw_rms","c=e");
	e_raw_error->SetMinimum(-0.03);
	e_raw_error->SetMaximum(0.03);
	e_raw_error->Draw("colz");
	t2->SaveAs("selectionPlots/iEtaPhiE_raw.pdf");
	t2->cd(1);
	r9_alternative->SetMinimum(-0.05);
	r9_alternative->SetMaximum(0.05);
	r9_alternative->Draw("colz");
	t2->cd(2);
	TH2D *r9_error_alternative = r9_alternative->ProjectionXY("rms_alternative","c=e");
	r9_error_alternative->SetMinimum(-0.03);
	r9_error_alternative->SetMaximum(0.03);
	r9_error_alternative->Draw("colz");
	t2->SaveAs("selectionPlots/iEtaPhiR9_alternative.pdf");
	t2->cd(1);
	e_sc_alternative->SetMinimum(-0.05);
	e_sc_alternative->SetMaximum(0.05);
	e_sc_alternative->Draw("colz");
	t2->cd(2);
	TH2D *e_sc_error_alternative = e_sc_alternative->ProjectionXY("esc_rms_alternative","c=e");
	e_sc_error_alternative->SetMinimum(-0.03);
	e_sc_error_alternative->SetMaximum(0.03);
	e_sc_error_alternative->Draw("colz");
	t2->SaveAs("selectionPlots/iEtaPhiEsc_alternative.pdf");
gStyle->SetOptStat(1);


	//Draw Z->ee istos
	TCanvas *t3 = new TCanvas("t3","t3",600,600);
	t3->Divide(2,2);
	THStack stEEGold ("stEEGold","InvMass of High R9 e's in EE");
	THStack stEBGold ("stEBGold","InvMass of High R9 e's in EB");
	THStack stEEBad ("stEEBad","InvMass of Low R9 e's in EE");
	THStack stEBBad ("stEBBad","InvMass of Low R9 e's in EB");
	stEEGold.Add(EEinvMassGold7);
	stEEGold.Add(EEinvMassGold5);
	stEBGold.Add(EBinvMassGold7);
	stEBGold.Add(EBinvMassGold5);
	stEEBad.Add( EEinvMassBad7);
	stEEBad.Add (EEinvMassBad5);
	stEBBad.Add( EBinvMassBad7);
	stEBBad.Add( EBinvMassBad5);
	THStack stEEGold_raw ("stEEGold_raw","InvMass of High R9 e's in EE");
	THStack stEBGold_raw ("stEBGold_raw","InvMass of High R9 e's in EB");
	THStack stEEBad_raw ("stEEBad_raw","InvMass of Low R9 e's in EE");
	THStack stEBBad_raw ("stEBBad_raw","InvMass of Low R9 e's in EB");
	stEEGold_raw.Add(EEinvMassGold_raw7);
	stEEGold_raw.Add(EEinvMassGold_raw5);
	stEBGold_raw.Add(EBinvMassGold_raw7);
	stEBGold_raw.Add(EBinvMassGold_raw5);
	stEEBad_raw.Add(EEinvMassBad_raw7);
	stEEBad_raw.Add(EEinvMassBad_raw5);
	stEBBad_raw.Add(EBinvMassBad_raw7);
	stEBBad_raw.Add(EBinvMassBad_raw5);
	THStack stEEGold_corr ("stEEGold_corr","InvMass of High R9 e's in EE");
	THStack stEBGold_corr ("stEBGold_corr","InvMass of High R9 e's in EB");
	THStack stEEBad_corr ("stEEBad_corr","InvMass of Low R9 e's in EE");
	THStack stEBBad_corr ("stEBBad_corr","InvMass of Low R9 e's in EB");
	stEEGold_corr.Add(EEinvMassGold_corr7);
	stEEGold_corr.Add(EEinvMassGold_corr5);
	stEBGold_corr.Add(EBinvMassGold_corr7);
	stEBGold_corr.Add(EBinvMassGold_corr5);
	stEEBad_corr.Add(EEinvMassBad_corr7);
	stEEBad_corr.Add(EEinvMassBad_corr5);
	stEBBad_corr.Add(EBinvMassBad_corr7);
	stEBBad_corr.Add(EBinvMassBad_corr5);
	t3->cd(1);
	stEEGold.Draw("nostack");
	t3->cd(2);
	stEBGold.Draw("nostack");
	t3->cd(3);
	stEEBad.Draw("nostack");
	t3->cd(4);
	stEBBad.Draw("nostack");
	t3->SaveAs("selectionPlots/Zee.pdf");

		t3->cd(1);
			stEEGold_raw.Draw("nostack");
			t3->cd(2);
			stEBGold_raw.Draw("nostack");
			t3->cd(3);
			stEEBad_raw.Draw("nostack");
			t3->cd(4);
			stEBBad_raw.Draw("nostack");
			t3->SaveAs("selectionPlots/Zee_raw.pdf");

			t3->cd(1);
			stEEGold_corr.Draw("nostack");
			t3->cd(2);
			stEBGold_corr.Draw("nostack");
			t3->cd(3);
			stEEBad_corr.Draw("nostack");
			t3->cd(4);
			stEBBad_corr.Draw("nostack");
			t3->SaveAs("selectionPlots/Zee_corr.pdf");
	 
std::cout 

<<EEinvMassGold7->GetName()<< " : RMS= " <<EEinvMassGold7->GetRMS() << std::endl
<<EEinvMassGold5->GetName()<< " : RMS= " <<EEinvMassGold5->GetRMS()<< std::endl
<<EBinvMassGold7->GetName()<< " : RMS= " <<EBinvMassGold7->GetRMS()<< std::endl
<<EBinvMassGold5->GetName()<< " : RMS= " <<EBinvMassGold5->GetRMS()<< std::endl
<<EEinvMassBad7 ->GetName()<< " : RMS= " <<EEinvMassBad7 ->GetRMS()<< std::endl
<<EEinvMassBad5->GetName() << " : RMS= " <<EEinvMassBad5->GetRMS() << std::endl
<<EBinvMassBad7->GetName() << " : RMS= " <<EBinvMassBad7->GetRMS() << std::endl
<<EBinvMassBad5->GetName() << " : RMS= " <<EBinvMassBad5->GetRMS() << std::endl

<<EEinvMassGold_corr7->GetName()<< " : RMS= " <<EEinvMassGold_corr7->GetRMS() << std::endl
<<EEinvMassGold_corr5->GetName()<< " : RMS= " <<EEinvMassGold_corr5->GetRMS()<< std::endl
<<EBinvMassGold_corr7->GetName()<< " : RMS= " <<EBinvMassGold_corr7->GetRMS()<< std::endl
<<EBinvMassGold_corr5->GetName()<< " : RMS= " <<EBinvMassGold_corr5->GetRMS()<< std::endl
<<EEinvMassBad_corr7 ->GetName()<< " : RMS= " <<EEinvMassBad_corr7 ->GetRMS()<< std::endl
<<EEinvMassBad_corr5->GetName() << " : RMS= " <<EEinvMassBad_corr5->GetRMS() << std::endl
<<EBinvMassBad_corr7->GetName() << " : RMS= " <<EBinvMassBad_corr7->GetRMS() << std::endl
<<EBinvMassBad_corr5->GetName() << " : RMS= " <<EBinvMassBad_corr5->GetRMS() << std::endl


<<EEinvMassGold_raw7->GetName()<< " : RMS= " <<EEinvMassGold_raw7->GetRMS() << std::endl
<<EEinvMassGold_raw5->GetName()<< " : RMS= " <<EEinvMassGold_raw5->GetRMS()<< std::endl
<<EBinvMassGold_raw7->GetName()<< " : RMS= " <<EBinvMassGold_raw7->GetRMS()<< std::endl
<<EBinvMassGold_raw5->GetName()<< " : RMS= " <<EBinvMassGold_raw5->GetRMS()<< std::endl
<<EEinvMassBad_raw7 ->GetName()<< " : RMS= " <<EEinvMassBad_raw7 ->GetRMS()<< std::endl
<<EEinvMassBad_raw5->GetName() << " : RMS= " <<EEinvMassBad_raw5->GetRMS() << std::endl
<<EBinvMassBad_raw7->GetName() << " : RMS= " <<EBinvMassBad_raw7->GetRMS() << std::endl
<<EBinvMassBad_raw5->GetName() << " : RMS= " <<EBinvMassBad_raw5->GetRMS() << std::endl;

	//end function
	return 0;
}
