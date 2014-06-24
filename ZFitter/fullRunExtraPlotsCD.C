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
#include "TPaveText.h"
#include "TFitResultPtr.h"


// EEFECTIVE SIGMA FUNCYION DEFINITION //
const char* gausSigma_cstr;
double gausSigma_val;
Double_t gausSigma(TH1* h)
{

	Double_t ave = h->GetMean();
	TF1 *g1= new TF1("g1","gaus",ave-2,ave+2);
	h->Fit("g1","NQ");

	cout << ave << " "<< g1->GetParameter(0)<< " "<< g1->GetParameter(1)<< " "<< g1->GetParameter(2)<< " "<<endl;

	Double_t mu= g1->GetParameter(2);


	std::ostringstream ostr1gaus;

	ostr1gaus << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) <<  mu ;
	cout << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) <<  mu << endl;
	std::string gausSigma_str = ostr1gaus.str();
	gausSigma_cstr = gausSigma_str.c_str();
	gausSigma_val = mu;

	return mu;
}

const char* effSigma_cstr;
double effSigma_val;
Double_t effSigma(TH1 * hist)
{

	TAxis *xaxis = hist->GetXaxis();
	Int_t nb = xaxis->GetNbins();
	//	cout << nb << endl;
	if(nb < 10) {
		cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
		return 0.;
	}

	Double_t bwid = xaxis->GetBinWidth(1);
	if(bwid == 0) {
		cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
		return 0.;
	}
	Double_t xmin = xaxis->GetXmin();
	Double_t ave = hist->GetMean();
	Double_t rms = hist->GetRMS();

	Double_t total=0.;
	for(Int_t i=0; i<nb+2; i++) {
		total+=hist->GetBinContent(i);
	}
	// if(total < 100.) {
	// cout << "effsigma: Too few entries " << total << endl;
	// return 0.;
	// }
	Int_t ierr=0;
	Int_t ismin=999;

	Double_t rlim=0.683*total;
	Int_t nrms=rms/(bwid); // Set scan size to +/- rms
	if(nrms > nb/10) nrms=nb/10; // Could be tuned...

	Double_t widmin=9999999.;
	for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
		Int_t ibm=(ave-xmin)/bwid+1+iscan;
		Double_t x=(ibm-0.5)*bwid+xmin;
		Double_t xj=x;
		Double_t xk=x;
		Int_t jbm=ibm;
		Int_t kbm=ibm;
		Double_t bin=hist->GetBinContent(ibm);
		total=bin;
		for(Int_t j=1;j<nb;j++){
			if(jbm < nb) {
				jbm++;
				xj+=bwid;
				bin=hist->GetBinContent(jbm);
				total+=bin;
				if(total > rlim) break;
			}
			else ierr=1;
			if(kbm > 0) {
				kbm--;
				xk-=bwid;
				bin=hist->GetBinContent(kbm);
				total+=bin;
				if(total > rlim) break;
			}
			else ierr=1;
		}
		Double_t dxf=(total-rlim)*bwid/bin;
		Double_t wid=(xj-xk+bwid-dxf)*0.5;
		if(wid < widmin) {
			widmin=wid;
			ismin=iscan;
		}
	}
	if(ismin == nrms || ismin == -nrms) ierr=3;
	if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
	//cout <<widmin << endl;
	std::ostringstream ostr1;

	ostr1 << "#sigma_{eff} = " <<  std::fixed <<  std::setprecision(3) <<  widmin ;
	cout << "#sigma_{eff} = " <<  std::fixed <<  std::setprecision(3) <<  widmin << endl;
	std::string effSigma_str = ostr1.str();
	effSigma_cstr = effSigma_str.c_str();
	//	cout << effSigma_str << endl <<effSigma_cstr << endl;
	effSigma_val=widmin;
	return widmin;

}


int fullRunExtraPlotsCD() 
{
	//Define Histograms
	//Simple quantitiy comaprisons
	TH1F *corrE_hist_EE = new TH1F("energySCEle_corr_EE","energySCEle_corr_EE fractional difference",1000,-0.2,0.2);
	TH1F *corrE_hist_EE_7 = new TH1F("energySCEle_corr_EE_7","energySCEle_corr_EE 7",1000,0,400);
	TH1F *corrE_hist_EE_5 = new TH1F("energySCEle_corr_EE_5","energySCEle_corr_EE 5",1000,0,400);
	TH1F *corrE_hist_EB = new TH1F("energySCEle_corr_EB","energySCEle_corr_EB fractional difference",1000,-0.2,0.2);
	TH1F *corrE_hist_EB_7 = new TH1F("energySCEle_corr_EB_7","energySCEle_corr_EB 7",1000,0,400);
	TH1F *corrE_hist_EB_5 = new TH1F("energySCEle_corr_EB_5","energySCEle_corr_EB 5",1000,0,400);

	//Z->ee distributions

	TH1F *EEinvMassGold_corr7 = new TH1F("EEinvMassGold_corr7","EEinvMassGold_corr",100,40,140);
	TH1F *EEinvMassBad_corr7 = new TH1F( "EinvMassBad_corr7" ,"EEinvMassBad_corr" ,100,40,140);
	TH1F *EBinvMassGold_corr7 = new TH1F("BinvMassGold_corr7","EBinvMassGold_corr",100,40,140);
	TH1F *EBinvMassBad_corr7 = new TH1F( "BinvMassBad_corr7 ","EBinvMassBad_corr ",100,40,140);
	TH1F *EEinvMassGold_corr5 = new TH1F("EEinvMassGold_corr5","EEinvMassGold_corr",100,40,140);
	TH1F *EEinvMassBad_corr5 = new TH1F( "EinvMassBad_corr5" ,"EEinvMassBad_corr" ,100,40,140);
	TH1F *EBinvMassGold_corr5 = new TH1F("BinvMassGold_corr5","EBinvMassGold_corr",100,40,140);
	TH1F *EBinvMassBad_corr5 = new TH1F( "BinvMassBad_corr5 ","EBinvMassBad_corr ",100,40,140);
	TH1F *EE_corr5 = new TH1F("E_corr5","EE_corr",1000,40,140);
	TH1F *EE_corr7 = new TH1F( "E_corr7 ","EE_corr ",1000,40,140);

	//nPV plots
	int histCounter;
	TH1F **Histos_eff_EELoR9;
	TH1F **Histos_eff_EEHiR9;
	TH1F **Histos_eff_EBLoR9;
	TH1F **Histos_eff_EBHiR9;
	Histos_eff_EELoR9 = new TH1F* [80];
	Histos_eff_EEHiR9 = new TH1F* [80];
	Histos_eff_EBLoR9 = new TH1F* [80];
	Histos_eff_EBHiR9 = new TH1F* [80];
	TProfile *clusterSize = new TProfile("size","size",40,0,40);
	TProfile *clusterSize_5 = new TProfile("size_5","size 5",40,0,40);
	TProfile *clusterSize_7 = new TProfile("size_7","size 7",40,0,40);
	TH1F **Histos_eff_pT40;
	TH1F **Histos_eff_pT50;
	TH1F **Histos_eff_pT60;
	Histos_eff_pT40 = new TH1F* [80];
	Histos_eff_pT50 = new TH1F* [80];
	Histos_eff_pT60 = new TH1F* [80];
	TProfile *ratio_profile_pT40 = new TProfile("ratio_pT40","ratio40",40,0,40);
	TProfile *ratio_profile_pT50 = new TProfile("ratio_pT50","ratio50",40,0,40);
	TProfile *ratio_profile_pT60 = new TProfile("ratio_pT50plus","ratio50",40,0,40);
	gROOT->SetBatch();


	//DEFINE
	for (histCounter=0; histCounter<80; histCounter++)
	{

		std::ostringstream name_EELoR9;
		std::ostringstream name_EEHiR9;
		std::ostringstream name_EBLoR9;
		std::ostringstream name_EBHiR9;
		std::ostringstream name_pT40;
		std::ostringstream name_pT50;
		std::ostringstream name_pT60;
		name_EELoR9 << "PV_hist_"<< histCounter << "_EELoR9";
		name_EEHiR9 << "PV_hist_"<< histCounter << "_EEHiR9";
		name_EBLoR9 << "PV_hist_"<< histCounter << "_EBLoR9";
		name_EBHiR9 << "PV_hist_"<< histCounter << "_EBHiR9";
		name_pT40<< "PV_hist_"<< histCounter << "_pT40";
		name_pT50<< "PV_hist_"<< histCounter << "_pT50";
		name_pT60<< "PV_hist_"<< histCounter << "_pT60";

		std::ostringstream title_EELoR9;
		std::ostringstream title_EEHiR9;
		std::ostringstream title_EBLoR9;
		std::ostringstream title_EBHiR9;
		std::ostringstream title_pT40;
		std::ostringstream title_pT50;
		std::ostringstream title_pT60;
		name_EELoR9 << "PV_hist_"<< histCounter << "_EELoR9";
		title_EELoR9 << "PV Histogram " << histCounter << "_EELoR9";
		title_EEHiR9 << "PV Histogram " << histCounter << "_EEHiR9";
		title_EBLoR9 << "PV Histogram " << histCounter << "_EBLoR9";
		title_EBHiR9 << "PV Histogram " << histCounter << "_EBHiR9";
		title_pT40<< "PV Histogram " << histCounter << "_pT40";
		title_pT50<< "PV Histogram " << histCounter << "_pT50";
		title_pT60<< "PV Histogram " << histCounter << "_pT60";

		Histos_eff_EELoR9[histCounter] = new TH1F(title_EELoR9.str().c_str(), name_EELoR9.str().c_str(),100,40,140);
		Histos_eff_EEHiR9[histCounter] = new TH1F(title_EEHiR9.str().c_str(), name_EEHiR9.str().c_str(),100,40,140);
		Histos_eff_EBLoR9[histCounter] = new TH1F(title_EBLoR9.str().c_str(), name_EBLoR9.str().c_str(),100,40,140);
		Histos_eff_EBHiR9[histCounter] = new TH1F(title_EBHiR9.str().c_str(), name_EBHiR9.str().c_str(),100,40,140);
		Histos_eff_pT40[histCounter] = new TH1F(title_pT40.str().c_str(), name_pT40.str().c_str(),100,40,140);
		Histos_eff_pT50[histCounter] = new TH1F(title_pT50.str().c_str(), name_pT50.str().c_str(),100,40,140);
		Histos_eff_pT60[histCounter] = new TH1F(title_pT60.str().c_str(), name_pT60.str().c_str(),100,40,140);

		Histos_eff_EELoR9[histCounter]->Sumw2();
		Histos_eff_EEHiR9[histCounter]->Sumw2();
		Histos_eff_EBLoR9[histCounter]->Sumw2();
		Histos_eff_EBHiR9[histCounter]->Sumw2();
		Histos_eff_pT40[histCounter]->Sumw2();
		Histos_eff_pT50[histCounter]->Sumw2();
		Histos_eff_pT60[histCounter]->Sumw2();

	}


	//Sumw2 settings for standard histos. WEIGHTS
	corrE_hist_EB->Sumw2();

	corrE_hist_EE->Sumw2();
	EEinvMassGold_corr7->Sumw2();
	EEinvMassBad_corr7 ->Sumw2();
	EBinvMassGold_corr7->Sumw2();
	EBinvMassBad_corr7 ->Sumw2();
	EEinvMassGold_corr5->Sumw2();
	EEinvMassBad_corr5 ->Sumw2();
	EBinvMassGold_corr5->Sumw2();
	EBinvMassBad_corr5 ->Sumw2();

	corrE_hist_EB_5 ->Sumw2();
	corrE_hist_EB_7 ->Sumw2();
	corrE_hist_EE_5 ->Sumw2();
	corrE_hist_EE_7 ->Sumw2();

	//Cosmetics	
	corrE_hist_EB_5->SetLineColor(kRed);
	corrE_hist_EE_5->SetLineColor(kRed);
	EEinvMassGold_corr5->SetLineColor(kRed);
	EEinvMassBad_corr5 ->SetLineColor(kRed);
	EBinvMassGold_corr5->SetLineColor(kRed);
	EBinvMassBad_corr5 ->SetLineColor(kRed);
	clusterSize_5->SetLineColor(kRed);

	TFile **Files_53X;
	TFile **Files_70X;
	TFile **Files_Map;
	Files_53X = new TFile*[4];
	Files_70X = new TFile*[4];
	Files_Map = new TFile*[4];



	//Open Map and both data sources.
	Files_Map[0]  = TFile::Open("tmp/Map_d1-louieTest2.root");
	Files_70X[0] = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2/190645-193621/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2-190645-193621.root");
	Files_53X[0]  = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-22Jan-v1/190645-193621/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012A-22Jan-v1-190645-193621.root");
	Files_Map[1]= TFile::Open("tmp/Map_d2-louieTest2.root");
	Files_70X[1]= TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012B-15Apr-v1/193834-196531/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012B-15Apr-v1-193834-196531.root");
	Files_53X[1] = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012B-22Jan-v1/193834-196531/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012B-22Jan-v1-193834-196531.root");
	Files_Map[2]= TFile::Open("tmp/Map_d3-louieTest2.root");
	Files_70X[2]= TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012C-15Apr-v1/198022-203742/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012C-15Apr-v1-198022-203742.root");
	Files_53X[2] = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012C-22Jan-v1/198022-203742/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012C-22Jan-v1-198022-203742.root");
	Files_Map[3]= TFile::Open("tmp/Map_d4-louieTest2.root");
	Files_70X[3]= TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012D-15Apr-v1/203777-208686/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012D-15Apr-v1-203777-208686.root");
	Files_53X[3] = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012D-22Jan-v1/203777-208686/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012D-22Jan-v1-203777-208686.root");

	//LOOPS
	for(Int_t overLoop=2;overLoop<4;overLoop++)
	{


		//Define respective TTrees


		TTree *chain7 =  (TTree*)Files_70X[overLoop]->Get("selected");
		TTree *chain5 =  (TTree*)Files_53X[overLoop]->Get("selected"); 
		TTree *tree_Map =  (TTree*)Files_Map[overLoop]->Get("Map"); 


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
		Int_t nPV7=0;
		Int_t nHitsSCEle7[2]={0,0};
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
		Int_t nPV5=0;
		Int_t nHitsSCEle5[2]={0,0};
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
		TBranch *b_nPV7;
		TBranch *b_nHitsSCEle7;
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
		TBranch *b_nPV5;
		TBranch *b_nHitsSCEle5;
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
		chain5->SetBranchAddress("invMass_SC_regrCorrSemiParV5_ele", &invMass_corr5, &b_invMass_corr5);
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
		chain5->SetBranchAddress("nPV",&nPV5,&b_nPV5);
		chain5->SetBranchAddress("nHitsSCEle",nHitsSCEle5,&b_nPV5);
		chain7->SetBranchAddress("nPV",&nPV7,&b_nPV7);
		chain5->SetBranchAddress("nHitsSCEle",nHitsSCEle5,&b_nHitsSCEle5);
		chain7->SetBranchAddress("nHitsSCEle",nHitsSCEle7,&b_nHitsSCEle7);

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





		//Loops over first TTree	
		for(Int_t loop =0; loop < 

				chain5->GetEntries()
				//	10000	
				; loop++)
			//SIZE
		{
			// Progress tracker
			if (loop%10000==0) 
			{ 
				cout << loop << " / " << chain5->GetEntries() << endl;
	//			cout << PtEle5[0] << " " << PtEle5[1] << " "<< PtEle7[0]<< " " << PtEle7[1]<< endl;
	//			cout << etaSCEle5[0] << " " << etaSCEle5[1] << " "<< etaSCEle7[0]<< " " << etaSCEle7[1]<< endl;
	//			cout << phiSCEle5[0] << " " << phiSCEle5[1] << " "<< phiSCEle7[0]<< " " << phiSCEle7[1]<< endl
			}

			//Get entry in firs TTree and corresponding Map entry	
			tree_Map->GetEntry(loop);
			chain5->GetEntry(loop);

			//Skip Entries where there is no match
			if (entryNumber==-1) continue;
			if (loop ==537558 || loop== 1706308) continue;
			if (loop ==1663201 || loop==1743949|| loop==1807467) continue;
			if (loop ==5380 || loop==9368 || loop ==214927 || loop==238203 || loop==1807467 || loop==1922998 || loop==2224825 || loop==2225083 || loop==2278786 || loop==2360602 || loop==2523867 || loop==2549355 || loop==2565052 || loop==2713956 || loop==3075338) continue;
			if (loop ==111743 || loop==118334 || loop==206385 || loop==252009 || loop==300763 || loop==329568|| loop==330032 || loop ==340427 || loop==428721 || loop==471861 || loop==474983 || loop==490438 || loop==497148 || loop==499118 || loop==543006 || loop==547608 || loop==623201 || loop==630981 || loop==631337 || loop==632624 || loop==780105 || loop==821481 || loop==910510 || loop==916414 || loop==953734 || loop==1007707 || loop==1030291 || loop==1030292 || loop ==1034067 || loop==1076924 || loop==1084153 || loop==1085899 || loop==1160607 || loop==1185241 || loop==1291573 || loop==1294912 || loop==1298424||loop==1302748 || loop==1315395 || loop==1342609 || loop==1343752 || loop==1354931 || loop==1387467 || loop==1388841 || loop==1478640 || loop==1479175|| loop==1488490 || loop==1488601 || loop==1581724 || loop==1614967 || loop==1639044|| loop==1654105 || loop==1657316 || loop==1670642 || loop==1688603 || loop==1699287 || loop==1708662 || loop==1710986 || loop==1715556 || loop==1750385 || loop==1773280 || loop==1774122 || loop==1775257 || loop==1798099 || loop==1803465 || loop==1803588 || loop==1833543 || loop==1842084 || loop==1874487 || loop==1894180 || loop==1906029 || loop==1911685 || loop==1944677 || loop==2048863 || loop==2064783 || loop==2065200 || loop==2072618 || loop==2092166 || loop==2093154 || loop==2094338 || loop==2250628 || loop==2251854 || loop==2339125|| loop==2355749 || loop==2358017 || loop==2387750 || loop==2403797 || loop==2411403 || loop==2479599 || loop==2500319 || loop==2532993 || loop==2539760 || loop==2557339 || loop==2606066 || loop==2609357 || loop==2667359 || loop==2667494 || loop==2686057 || loop==2695787 || loop==2773919 || loop==2892501 || loop==3033892 ||loop==3100742 || loop==3141489) continue;
			//if (PtEle5[0] ==0 || PtEle5[1] ==0 || PtEle7[0]==0 || PtEle7[1]==0) continue;

			//Get corresponding entry in TTree 2  using information from Map
			chain7->GetEntry(entryNumber);

			//	cout << nHitsSCEle5[0] << " " << nHitsSCEle7[0] << endl;
			//		cout << nHitsSCEle5[1] << " " << nHitsSCEle7[1] << endl;

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
			//Double_t iPhi5[2] = { floor(phiSCEle5[0]*57.30+0.5)+181, floor(phiSCEle5[1]*57.30+0.5)+181 };
			//Double_t iEta5[2] = { ((etaSCEle5[0]>0)- (etaSCEle5[0]<0))*(floor(fabs(etaSCEle5[0]*57.47)+0.5)+1), ((etaSCEle5[1]>0) - (etaSCEle5[1]<0))*(floor(fabs(etaSCEle5[1]*57.46)+0.5)+1 )};

			//Fill Z->ee histograms
			if((fabs(etaSCEle5[0])>1.566 ) && (fabs(etaSCEle5[1])>1.566))
			{

				if (R9Ele5[0]>0.94 && R9Ele5[1]>0.94)
				{
					EEinvMassGold_corr5->Fill(invMass_corr5,selection5);
					Histos_eff_EEHiR9[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>30 && PtEle5[0]<40 && PtEle5[1]>30 && PtEle5[0]<40)	Histos_eff_pT40[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>40 && PtEle5[0]<50 && PtEle5[1]>40 && PtEle5[0]<50)	Histos_eff_pT50[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>50 && PtEle5[0]>50)	Histos_eff_pT60[nPV5]->Fill(invMass_corr5,selection5);
				}	

				if (R9Ele5[0]<0.94 && R9Ele5[1]<0.94)
				{
					EEinvMassBad_corr5->Fill(invMass_corr5,selection5);
					Histos_eff_EELoR9[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>30 && PtEle5[0]<40 && PtEle5[1]>30 && PtEle5[0]<40)	Histos_eff_pT40[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>40 && PtEle5[0]<50 && PtEle5[1]>40 && PtEle5[0]<50)	Histos_eff_pT50[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>50 && PtEle5[0]>50)	Histos_eff_pT60[nPV5]->Fill(invMass_corr5,selection5);
				}
			}
			if((fabs(etaSCEle5[0])<1.444) && (fabs(etaSCEle5[1])<1.4442)) 
			{

				if (R9Ele5[0]>0.94 && R9Ele5[1]>0.94)
				{
					EBinvMassGold_corr5->Fill(invMass_corr5,selection5);
					Histos_eff_EBHiR9[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>30 && PtEle5[0]<40 && PtEle5[1]>30 && PtEle5[0]<40)	Histos_eff_pT40[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>40 && PtEle5[0]<50 && PtEle5[1]>40 && PtEle5[0]<50)	Histos_eff_pT50[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>50 && PtEle5[0]>50)	Histos_eff_pT60[nPV5]->Fill(invMass_corr5,selection5);
				}	
				if (R9Ele5[0]<0.94 && R9Ele5[1]<0.94)
				{
					EBinvMassBad_corr5->Fill(invMass_corr5,selection5);
					Histos_eff_EBLoR9[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>30 && PtEle5[0]<40 && PtEle5[1]>30 && PtEle5[0]<40)	Histos_eff_pT40[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>40 && PtEle5[0]<50 && PtEle5[1]>40 && PtEle5[0]<50)	Histos_eff_pT50[nPV5]->Fill(invMass_corr5,selection5);
					if(PtEle5[1]>50 && PtEle5[0]>50)	Histos_eff_pT60[nPV5]->Fill(invMass_corr5,selection5);
				}
			}

			if((fabs(etaSCEle7[0])>1.566 ) && (fabs(etaSCEle7[1])>1.566))
			{

				if (R9Ele7[0]>0.94 && R9Ele7[1]>0.94)
				{
					EEinvMassGold_corr7->Fill(invMass_corr7,selection5);
					Histos_eff_EEHiR9[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>30 && PtEle7[0]<40 && PtEle7[1]>30 && PtEle7[0]<40)	Histos_eff_pT40[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>40 && PtEle7[0]<50 && PtEle7[1]>40 && PtEle7[0]<50)	Histos_eff_pT50[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>50 && PtEle7[0]>50)	Histos_eff_pT60[nPV7+40]->Fill(invMass_corr7,selection5);
				}	
				if (R9Ele7[0]<0.94 && R9Ele7[1]<0.94)
				{
					EEinvMassBad_corr7->Fill(invMass_corr7,selection5);
					Histos_eff_EELoR9[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>30 && PtEle7[0]<40 && PtEle7[1]>30 && PtEle7[0]<40)	Histos_eff_pT40[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>40 && PtEle7[0]<50 && PtEle7[1]>40 && PtEle7[0]<50)	Histos_eff_pT50[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>50 && PtEle7[0]>50)	Histos_eff_pT60[nPV7+40]->Fill(invMass_corr7,selection5);
				}
			}
			if((fabs(etaSCEle7[0])<1.444) && (fabs(etaSCEle7[1])<1.4442)) 
			{

				if (R9Ele7[0]>0.94 && R9Ele7[1]>0.94)
				{
					EBinvMassGold_corr7->Fill(invMass_corr7,selection5);
					Histos_eff_EBHiR9[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>30 && PtEle7[0]<40 && PtEle7[1]>30 && PtEle7[0]<40)	Histos_eff_pT40[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>40 && PtEle7[0]<50 && PtEle7[1]>40 && PtEle7[0]<50)	Histos_eff_pT50[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>50 && PtEle7[0]>50)	Histos_eff_pT60[nPV7+40]->Fill(invMass_corr7,selection5);
				}	
				if (R9Ele7[0]<0.94 && R9Ele7[1]<0.94)
				{
			//	cout << invMass_corr7 << "  "<< selection5<< endl;
				EBinvMassBad_corr7->Fill(invMass_corr7,selection5);
					Histos_eff_EBLoR9[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>30 && PtEle7[0]<40 && PtEle7[1]>30 && PtEle7[0]<40)	Histos_eff_pT40[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>40 && PtEle7[0]<50 && PtEle7[1]>40 && PtEle7[0]<50)	Histos_eff_pT50[nPV7+40]->Fill(invMass_corr7,selection5);
					if(PtEle7[1]>50 && PtEle7[0]>50)	Histos_eff_pT60[nPV7+40]->Fill(invMass_corr7,selection5);
				}
			}

if(fabs(etaSCEle5[0])>1.5 || fabs(etaSCEle5[1]>1.5))
{EE_corr5->Fill(invMass_corr5,selection5);}

if(fabs(etaSCEle7[0])>1.5 || fabs(etaSCEle7[1]>1.5))
{EE_corr7->Fill(invMass_corr7,selection5);}




			//Fill basic histograms
			if (eleIndex[0]==1)
			{
				if(etaSCEle5[0]>1.5)
				{
					corrE_hist_EE->Fill((energySCEle_corr7[0]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
					corrE_hist_EE_5 ->Fill((energySCEle_corr5[0]),selection5);
					corrE_hist_EE_7 ->Fill((energySCEle_corr7[0]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[0]-nHitsSCEle5[0],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[0],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[0],selection5);
				}

				if(etaSCEle5[1]>1.5)
				{
					corrE_hist_EE->Fill((energySCEle_corr7[1]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
					corrE_hist_EE_5 ->Fill((energySCEle_corr5[1]),selection5);
					corrE_hist_EE_7 ->Fill((energySCEle_corr7[1]),selection5);
					clusterSize->Fill(nHitsSCEle7[1]-nHitsSCEle5[1],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[1],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[1],selection5);
				}

				if(etaSCEle5[0]<1.5)
				{
					corrE_hist_EB->Fill((energySCEle_corr7[0]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
					corrE_hist_EB_5 ->Fill((energySCEle_corr5[0]),selection5);
					corrE_hist_EB_7 ->Fill((energySCEle_corr7[0]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[0]-nHitsSCEle5[0],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[1],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[1],selection5);
				}

				if(etaSCEle5[1]<1.5)
				{
					corrE_hist_EB->Fill((energySCEle_corr7[1]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
					corrE_hist_EB_5 ->Fill((energySCEle_corr5[1]),selection5);
					corrE_hist_EB_7 ->Fill((energySCEle_corr7[1]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[1]-nHitsSCEle5[1],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[1],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[1],selection5);
				}


			}

			if (eleIndex[0]==2)
			{

				if(etaSCEle5[0]>1.5)
				{
					corrE_hist_EE->Fill((energySCEle_corr7[1]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
					corrE_hist_EE_5 ->Fill((energySCEle_corr5[0]),selection5);
					corrE_hist_EE_7 ->Fill((energySCEle_corr7[1]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[1]-nHitsSCEle5[0],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[0],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[1],selection5);
				}

				if(etaSCEle5[1]>1.5)
				{
					corrE_hist_EE->Fill((energySCEle_corr7[0]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
					corrE_hist_EE_5 ->Fill((energySCEle_corr5[1]),selection5);
					corrE_hist_EE_7 ->Fill((energySCEle_corr7[0]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[0]-nHitsSCEle5[1],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[1],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[0],selection5);
				}

				if(etaSCEle5[0]<1.5)
				{
					corrE_hist_EB->Fill((energySCEle_corr7[1]-energySCEle_corr5[0])/energySCEle_corr5[0],selection5);
					corrE_hist_EB_5 ->Fill((energySCEle_corr5[0]),selection5);
					corrE_hist_EB_7 ->Fill((energySCEle_corr7[1]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[1]-nHitsSCEle5[0],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[0],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[1],selection5);
				}

				if(etaSCEle5[1]<1.5)
				{
					corrE_hist_EB->Fill((energySCEle_corr7[0]-energySCEle_corr5[1])/energySCEle_corr5[1],selection5);
					corrE_hist_EB_5 ->Fill((energySCEle_corr5[1]),selection5);
					corrE_hist_EB_7 ->Fill((energySCEle_corr7[0]),selection5);
					clusterSize->Fill(nPV5,nHitsSCEle7[0]-nHitsSCEle5[1],selection5);
					clusterSize_5->Fill(nPV5,nHitsSCEle5[1],selection5);
					clusterSize_7->Fill(nPV7,nHitsSCEle7[0],selection5);
				}

			}
			//end of loop	
		}

		//end of overLoop
	}
	//Draw basic Histos_eff

	gStyle->SetOptStat(0);
	//DRAW

	TMultiGraph *mg_gau_EELoR9 = new TMultiGraph();
	TMultiGraph *mg_gau_EEHiR9 = new TMultiGraph();
	TMultiGraph *mg_gau_EBLoR9 = new TMultiGraph();
	TMultiGraph *mg_gau_EBHiR9 = new TMultiGraph();
	TGraph *graph_gau_EELoR9 = new TGraph(40);
	TGraph *graph_gau_EEHiR9 = new TGraph(40);
	TGraph *graph_gau_EBLoR9 = new TGraph(40);
	TGraph *graph_gau_EBHiR9 = new TGraph(40);
	TGraph *graph_gau_EELoR9_7 = new TGraph(40);
	TGraph *graph_gau_EEHiR9_7 = new TGraph(40);
	TGraph *graph_gau_EBLoR9_7 = new TGraph(40);
	TGraph *graph_gau_EBHiR9_7 = new TGraph(40);
	mg_gau_EELoR9->SetTitle("EELoR9");
	mg_gau_EEHiR9->SetTitle("EEHiR9");
	mg_gau_EBLoR9->SetTitle("EBLoR9");
	mg_gau_EBHiR9->SetTitle("EBHiR9");

	for (histCounter =0 ; histCounter < 40; histCounter++)
	{
		//cout << Histos_gau[histCounter]->GetEntries() << " ";

		gausSigma(Histos_eff_EELoR9[histCounter]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EELoR9->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EEHiR9[histCounter]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EEHiR9->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EBLoR9[histCounter]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EBLoR9->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EBHiR9[histCounter]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EBHiR9->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EELoR9[histCounter+40]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EELoR9_7->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EEHiR9[histCounter+40]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EEHiR9_7->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EBLoR9[histCounter+40]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EBLoR9_7->SetPoint(histCounter,histCounter,gausSigma_val);}

		gausSigma(Histos_eff_EBHiR9[histCounter+40]);
		if (fabs(gausSigma_val) < 10000)
		{	graph_gau_EBHiR9_7->SetPoint(histCounter,histCounter,gausSigma_val);}
	}


	for (histCounter =0; histCounter <40 ; histCounter++)
	{
Float_t buffer;

effSigma(Histos_eff_pT40[histCounter]);
buffer= effSigma_val;
effSigma(Histos_eff_pT40[histCounter+40]);
if(fabs(effSigma_val) <10000&& buffer!=0  && fabs(buffer)<10000)
{
ratio_profile_pT40->Fill(histCounter,effSigma_val/buffer);
}
else ratio_profile_pT40->Fill(histCounter,0);

effSigma(Histos_eff_pT50[histCounter]);
buffer= effSigma_val;
effSigma(Histos_eff_pT50[histCounter+40]);
if(fabs(effSigma_val) <10000&& buffer!=0  && fabs(buffer)<10000)
{
ratio_profile_pT50->Fill(histCounter,effSigma_val/buffer);
}
else ratio_profile_pT50->Fill(histCounter,0);

effSigma(Histos_eff_pT60[histCounter]);
buffer= effSigma_val;
effSigma(Histos_eff_pT60[histCounter+40]);
if(fabs(effSigma_val) <10000&& buffer!=0  && fabs(buffer)<10000)
{
ratio_profile_pT60->Fill(histCounter,effSigma_val/buffer);
}
else ratio_profile_pT60->Fill(histCounter,0);

	}

	TCanvas *t1 = new TCanvas("t1","t1",600,600);

	graph_gau_EELoR9->SetMarkerColor(kBlue);
	graph_gau_EEHiR9->SetMarkerColor(kBlue);
	graph_gau_EBLoR9->SetMarkerColor(kBlue);
	graph_gau_EBHiR9->SetMarkerColor(kBlue);
	graph_gau_EELoR9->SetLineColor(kBlue);
	graph_gau_EEHiR9->SetLineColor(kBlue);
	graph_gau_EBLoR9->SetLineColor(kBlue);
	graph_gau_EBHiR9->SetLineColor(kBlue);
	graph_gau_EELoR9->SetMarkerStyle(7);
	graph_gau_EEHiR9->SetMarkerStyle(7);
	graph_gau_EBLoR9->SetMarkerStyle(7);
	graph_gau_EBHiR9->SetMarkerStyle(7);
	graph_gau_EELoR9_7->SetMarkerColor(kRed);
	graph_gau_EEHiR9_7->SetMarkerColor(kRed);
	graph_gau_EBLoR9_7->SetMarkerColor(kRed);
	graph_gau_EBHiR9_7->SetMarkerColor(kRed);
	graph_gau_EELoR9_7->SetLineColor(kRed);
	graph_gau_EEHiR9_7->SetLineColor(kRed);
	graph_gau_EBLoR9_7->SetLineColor(kRed);
	graph_gau_EBHiR9_7->SetLineColor(kRed);
	graph_gau_EELoR9_7->SetMarkerStyle(7);
	graph_gau_EEHiR9_7->SetMarkerStyle(7);
	graph_gau_EBLoR9_7->SetMarkerStyle(7);
	graph_gau_EBHiR9_7->SetMarkerStyle(7);
	mg_gau_EELoR9->Add(graph_gau_EELoR9_7);
	mg_gau_EEHiR9->Add(graph_gau_EEHiR9_7);
	mg_gau_EBLoR9->Add(graph_gau_EBLoR9_7);
	mg_gau_EBHiR9->Add(graph_gau_EBHiR9_7);
	mg_gau_EELoR9->Add(graph_gau_EELoR9);
	mg_gau_EEHiR9->Add(graph_gau_EEHiR9);
	mg_gau_EBLoR9->Add(graph_gau_EBLoR9);
	mg_gau_EBHiR9->Add(graph_gau_EBHiR9);
	mg_gau_EELoR9->SetMaximum(15);
	mg_gau_EEHiR9->SetMaximum(15);
	mg_gau_EBLoR9->SetMaximum(15);
	mg_gau_EBHiR9->SetMaximum(15);
	mg_gau_EELoR9->SetMinimum(0);
	mg_gau_EEHiR9->SetMinimum(0);
	mg_gau_EBLoR9->SetMinimum(0);
	mg_gau_EBHiR9->SetMinimum(0);


	TMultiGraph *mg_eff_EELoR9 = new TMultiGraph();
	TMultiGraph *mg_eff_EEHiR9 = new TMultiGraph();
	TMultiGraph *mg_eff_EBLoR9 = new TMultiGraph();
	TMultiGraph *mg_eff_EBHiR9 = new TMultiGraph();
	TGraph *graph_eff_EELoR9 = new TGraph(40);
	TGraph *graph_eff_EEHiR9 = new TGraph(40);
	TGraph *graph_eff_EBLoR9 = new TGraph(40);
	TGraph *graph_eff_EBHiR9 = new TGraph(40);
	TGraph *graph_eff_EELoR9_7 = new TGraph(40);
	TGraph *graph_eff_EEHiR9_7 = new TGraph(40);
	TGraph *graph_eff_EBLoR9_7 = new TGraph(40);
	TGraph *graph_eff_EBHiR9_7 = new TGraph(40);
	mg_eff_EELoR9->SetTitle("EELoR9");
	mg_eff_EEHiR9->SetTitle("EEHiR9");
	mg_eff_EBLoR9->SetTitle("EBLoR9");
	mg_eff_EBHiR9->SetTitle("EBHiR9");
	mg_eff_EELoR9->SetMaximum(15);
	mg_eff_EEHiR9->SetMaximum(15);
	mg_eff_EBLoR9->SetMaximum(15);
	mg_eff_EBHiR9->SetMaximum(15);
	mg_eff_EELoR9->SetMinimum(0);
	mg_eff_EEHiR9->SetMinimum(0);
	mg_eff_EBLoR9->SetMinimum(0);
	mg_eff_EBHiR9->SetMinimum(0);

	for (histCounter =0 ; histCounter < 40; histCounter++)
	{
		//cout << Histos_eff[histCounter]->GetEntries() << " ";

		effSigma(Histos_eff_EELoR9[histCounter]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EELoR9->SetPoint(histCounter,histCounter,effSigma_val);
//		graph_eff_EELoR9->SetPointError(histCounter,histCounter,effSigma_val);
		else	graph_eff_EELoR9->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EEHiR9[histCounter]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EEHiR9->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EEHiR9->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EBLoR9[histCounter]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EBLoR9->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EBLoR9->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EBHiR9[histCounter]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EBHiR9->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EBHiR9->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EELoR9[histCounter+40]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EELoR9_7->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EELoR9_7->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EEHiR9[histCounter+40]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EEHiR9_7->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EEHiR9_7->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EBLoR9[histCounter+40]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EBLoR9_7->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EBLoR9_7->SetPoint(histCounter,histCounter,0);

		effSigma(Histos_eff_EBHiR9[histCounter+40]);
		if (fabs(effSigma_val) < 10000)
			graph_eff_EBHiR9_7->SetPoint(histCounter,histCounter,effSigma_val);
		else	graph_eff_EBHiR9_7->SetPoint(histCounter,histCounter,0);
	}

	TCanvas *t2 = new TCanvas("t2","t2",600,600);

	graph_eff_EELoR9->SetMarkerColor(kBlue);
	graph_eff_EEHiR9->SetMarkerColor(kBlue);
	graph_eff_EBLoR9->SetMarkerColor(kBlue);
	graph_eff_EBHiR9->SetMarkerColor(kBlue);
	graph_eff_EELoR9->SetLineColor(kBlue);
	graph_eff_EEHiR9->SetLineColor(kBlue);
	graph_eff_EBLoR9->SetLineColor(kBlue);
	graph_eff_EBHiR9->SetLineColor(kBlue);
	graph_eff_EELoR9->SetMarkerStyle(7);
	graph_eff_EEHiR9->SetMarkerStyle(7);
	graph_eff_EBLoR9->SetMarkerStyle(7);
	graph_eff_EBHiR9->SetMarkerStyle(7);
	graph_eff_EELoR9_7->SetMarkerColor(kRed);
	graph_eff_EEHiR9_7->SetMarkerColor(kRed);
	graph_eff_EBLoR9_7->SetMarkerColor(kRed);
	graph_eff_EBHiR9_7->SetMarkerColor(kRed);
	graph_eff_EELoR9_7->SetLineColor(kRed);
	graph_eff_EEHiR9_7->SetLineColor(kRed);
	graph_eff_EBLoR9_7->SetLineColor(kRed);
	graph_eff_EBHiR9_7->SetLineColor(kRed);
	graph_eff_EELoR9_7->SetMarkerStyle(7);
	graph_eff_EEHiR9_7->SetMarkerStyle(7);
	graph_eff_EBLoR9_7->SetMarkerStyle(7);
	graph_eff_EBHiR9_7->SetMarkerStyle(7);
	mg_eff_EELoR9->Add(graph_eff_EELoR9_7);
	mg_eff_EEHiR9->Add(graph_eff_EEHiR9_7);
	mg_eff_EBLoR9->Add(graph_eff_EBLoR9_7);
	mg_eff_EBHiR9->Add(graph_eff_EBHiR9_7);
	mg_eff_EELoR9->Add(graph_eff_EELoR9);
	mg_eff_EEHiR9->Add(graph_eff_EEHiR9);
	mg_eff_EBLoR9->Add(graph_eff_EBLoR9);
	mg_eff_EBHiR9->Add(graph_eff_EBHiR9);

	t1->Divide(2,2);
	t1->cd(1);
	mg_eff_EELoR9->Draw("APL");
	t1->cd(2);
	mg_eff_EEHiR9->Draw("APL");
	t1->cd(3);
	mg_eff_EBLoR9->Draw("APL");
	t1->cd(4);
	mg_eff_EBHiR9->Draw("APL");

	t1->SaveAs("fullRunExtraPlotsCD/effSigma_vs_nPV.pdf");

	t2->Divide(2,2);
	t2->cd(1);
	mg_gau_EELoR9->Draw("APL");
	t2->cd(2);
	mg_gau_EEHiR9->Draw("APL");
	t2->cd(3);
	mg_gau_EBLoR9->Draw("APL");
	t2->cd(4);
	mg_gau_EBHiR9->Draw("APL");

	t2->SaveAs("fullRunExtraPlotsCD/gausSigma_vs_nPV.pdf");

	//Draw Z->ee istos
	TCanvas *t3 = new TCanvas("t3","t3",600,600);


	EEinvMassGold_corr7->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EEinvMassGold_corr5->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EBinvMassGold_corr7->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EBinvMassGold_corr5->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EEinvMassBad_corr7->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EEinvMassBad_corr5->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EBinvMassBad_corr7->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	EBinvMassBad_corr5->GetXaxis()->SetTitle("Invariant Mass (GeV)");

	TPaveText *pt0a =new TPaveText(0.15,0.6,0.35,0.85,"blNDC");
	pt0a->SetFillColor(0);
	pt0a->SetTextAlign(11);
	pt0a->SetBorderSize(0);
	pt0a->SetTextColor(kRed);
	TPaveText *pt0b =new TPaveText(0.65,0.6,0.85,0.85,"blNDC");
	pt0b->SetFillColor(0);
	pt0b->SetTextAlign(11);
	pt0b->SetBorderSize(0);
	pt0b->SetTextColor(kBlue);
	EEinvMassGold_corr7->Draw();
	EEinvMassGold_corr5->Draw("same");

	pt0a->AddText("53X");

	effSigma(EEinvMassGold_corr5);
	gausSigma(EEinvMassGold_corr5);

	std::ostringstream Label1;
	std::ostringstream Label2;
	Label1 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label2 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0a->AddText(Label1.str().c_str());
	pt0a->AddText(Label2.str().c_str());
	pt0a->Draw();

	pt0b->AddText("70X");
	effSigma(EEinvMassGold_corr7);
	gausSigma(EEinvMassGold_corr7);

	std::ostringstream Label3;
	std::ostringstream Label4;
	Label3 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label4 <<"#sigma_{eff} = "<<effSigma_val ; 

	pt0b->AddText(Label3.str().c_str());
	pt0b->AddText(Label4.str().c_str());
	pt0b->Draw();
	t3->SaveAs("fullRunExtraPlotsCD/Zee_corr_1.pdf");
	Label1.str("");
	Label1.clear();
	Label2.str("");
	Label2.clear();
	Label3.str("");
	Label3.clear();
	Label4.str("");
	Label4.clear();
	pt0a->Clear();
	pt0b->Clear();

	EBinvMassGold_corr7->Draw();
	EBinvMassGold_corr5->Draw("same");
	effSigma(EBinvMassGold_corr5);

	gausSigma(EBinvMassGold_corr5);
	pt0a->AddText("53X");

	Label1 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label2 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0a->AddText(Label1.str().c_str());
	pt0a->AddText(Label2.str().c_str());

	pt0a->Draw();
	effSigma(EBinvMassGold_corr7);
	gausSigma(EBinvMassGold_corr7);
	pt0b->AddText("70X");
	Label3 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label4 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0b->AddText(Label3.str().c_str());
	pt0b->AddText(Label4.str().c_str());
	pt0b->Draw();
	t3->SaveAs("fullRunExtraPlotsCD/Zee_corr_2.pdf");
	Label1.str("");
	Label1.clear();
	Label2.str("");
	Label2.clear();
	Label3.str("");
	Label3.clear();
	Label4.str("");
	Label4.clear();
	pt0a->Clear();
	pt0b->Clear();


	EEinvMassBad_corr7->Draw();
	EEinvMassBad_corr5->Draw("same");
	effSigma(EEinvMassBad_corr5);
	gausSigma(EEinvMassBad_corr5);

	pt0a->AddText("53X");
	Label1 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label2 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0a->AddText(Label1.str().c_str());
	pt0a->AddText(Label2.str().c_str());
	pt0a->Draw();

	effSigma(EEinvMassBad_corr7);
	gausSigma(EEinvMassBad_corr7);

	pt0b->AddText("70X");
	Label3 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label4 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0b->AddText(Label3.str().c_str());
	pt0b->AddText(Label4.str().c_str());
	pt0b->Draw();
	t3->SaveAs("fullRunExtraPlotsCD/Zee_corr_3.pdf");

	Label1.str("");
	Label1.clear();
	Label2.str("");
	Label2.clear();
	Label3.str("");
	Label3.clear();
	Label4.str("");
	Label4.clear();
	pt0a->Clear();
	pt0b->Clear();

	EBinvMassBad_corr7->Draw();
	EBinvMassBad_corr5->Draw("same");
	effSigma(EBinvMassBad_corr5);
	gausSigma(EBinvMassBad_corr5);

	pt0a->AddText("53X");
	Label1 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label2 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0a->AddText(Label1.str().c_str());
	pt0a->AddText(Label2.str().c_str());
	pt0a->Draw();


	effSigma(EBinvMassBad_corr7);
	gausSigma(EBinvMassBad_corr7);

	pt0b->AddText("70X");
	Label3 << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) << gausSigma_val <<std::endl;
	Label4 <<"#sigma_{eff} = "<<effSigma_val ; 
	pt0b->AddText(Label3.str().c_str());
	pt0b->AddText(Label4.str().c_str());
	pt0b->Draw();
	t3->SaveAs("fullRunExtraPlotsCD/Zee_corr_4.pdf");

	Label1.str("");
	Label1.clear();
	Label2.str("");
	Label2.clear();
	Label3.str("");
	Label3.clear();
	Label4.str("");
	Label4.clear();
	pt0a->Clear();
	pt0b->Clear();

	TCanvas *t4 = new TCanvas("t4","t4",600,600);
	clusterSize_5->Draw();
	clusterSize_7->Draw("Same");
	t4->SaveAs("fullRunExtraPlotsCD/clusterSize_dist.pdf");

	TCanvas *t5 = new TCanvas("t5","t5",600,600);
	clusterSize->Draw();
	t5->SaveAs("fullRunExtraPlotsCD/clusterSize_profile.pdf");

	TCanvas *t6 = new TCanvas("t6","t6",600,600);
	t6->Divide(1,3);
	t6->cd(1);
	ratio_profile_pT40->Draw();
	t6->cd(2);
	ratio_profile_pT50->Draw();
	t6->cd(3);
	ratio_profile_pT60->Draw();
	t6->SaveAs("fullRunExtraPlotsCD/ratio_profile.pdf");

TCanvas *t7= new TCanvas("t7","t7",600,600);
EE_corr5->SetLineColor(kRed);
EE_corr5->Draw();
EE_corr7->Draw("Same");

gausSigma(EE_corr5);
gausSigma(EE_corr7);
t7->SaveAs("fullRunExtraPlotsCD/t7.pdf");
	//end function
	return 0;
}
