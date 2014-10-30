#include <iostream>
#include <fstream>
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


// EFFECTIVE SIGMA and Gaussian Sigma  FUNCTION DEFINITIONs//

const char* gausSigma_cstr;
double gausSigma_val;
Double_t gausSigma(TH1* h)
{

	Double_t ave = h->GetMean();
	TF1 *g1= new TF1("g1","gaus",ave*0.9,ave*1.1);
	h->Fit("g1","N");

	std::cout << ave << " "<< g1->GetParameter(0)<< " "<< g1->GetParameter(1)<< " "<< g1->GetParameter(2)<< " "<<std::endl;

	Double_t mu= g1->GetParameter(2);


	std::ostringstream ostr1gaus;

	ostr1gaus << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) <<  mu ;
	std::cout << "#sigma_{Gaus} = " <<  std::fixed <<  std::setprecision(3) <<  mu << std::endl;
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
	//	std::cout << nb << std::endl;
	if(nb < 10) {
		std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
		return 0.;
	}

	Double_t bwid = xaxis->GetBinWidth(1);
	if(bwid == 0) {
		std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
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
	// std::cout << "effsigma: Too few entries " << total << std::endl;
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
	if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;
	//std::cout <<widmin << std::endl;
	std::ostringstream ostr1;

	ostr1 << "#sigma_{eff} = " <<  std::fixed <<  std::setprecision(3) <<  widmin ;
	//	std::cout << "#sigma_{eff} = " <<  std::fixed <<  std::setprecision(3) <<  widmin << std::endl;
	std::string effSigma_str = ostr1.str();
	effSigma_cstr = effSigma_str.c_str();
	//	std::cout << effSigma_str << std::endl <<effSigma_cstr << std::endl;
	effSigma_val=widmin;
	return widmin;

}


int fullRunPlots_OctA() 
{
	//Define Histograms
	//Simple quantitiy comaprisons
	TH1F *Delta_eN = new TH1F("Delta_eN","sanity - event",1000,-10000000,10000000); //sanity check that event and run numbers match
	TH1F *Delta_rN = new TH1F("Delta_rN","sanity -run ",1000,-10000000,10000000);

	TH1F *R9_hist_5 = new TH1F("R9Ele5","R9Ele5",1000,0,1.2);
	TH1F *R9_hist_7 = new TH1F("R9Ele7","R9Ele7",1000,0,1.2);
	TH1F *R9_hist = new TH1F("R9Ele","R9Ele difference",1000,-0.4,0.4);
	TH1F *R9_frac_hist = new TH1F("R9Ele","R9Ele frac",1000,-0.4,0.4);
	TH1F *R9_hist_5_un      = new TH1F("R9Ele5_un","R9Ele5_unmatched",1000,0,1.2);
	R9_hist_5       ->Sumw2(); 
	R9_hist_5_un    ->Sumw2(); 
	R9_hist_7       ->Sumw2(); 
	R9_hist->Sumw2(); 
	R9_frac_hist->Sumw2(); 

	TH1F *etaSC_hist = new TH1F("etaSCEle","etaSCEle difference",100,-0.001,0.001);
	TH1F *phiSC_hist = new TH1F("phiSCEle","phiSCEle difference",100,-0.001,0.001);
	TH1F *eta_hist = new TH1F("etaEle","etaEle difference",100,-0.001,0.001);
	TH1F *phi_hist = new TH1F("phiEle","phiEle difference",100,-0.001,0.001);
	TH1F *pModeGsf_hist = new TH1F("pModeGSf","pmode fractional difference",1000,-0.02,0.02);
	etaSC_hist->Sumw2(); 
	phiSC_hist->Sumw2(); 
	eta_hist->Sumw2(); 
	phi_hist->Sumw2(); 
	pModeGsf_hist ->Sumw2();

	TProfile *fbrem = new TProfile("fbrem","f_{brem} ratio vs #eta",500,-2.5,2.5);
	TProfile *corrE_profile = new TProfile("corrE_profile","corrE ratio vs #eta",500,-2.5,2.5);
	TProfile *R9_profile = new TProfile("R9_profile","R9 ratio vs #eta",500,-2.5,2.5);
	TH1F *PtCalcEle_hist5            = new TH1F("PtCalcEle5","PtCalc of electrons in 53X",1000,0,500);
	TH1F *PtCalcEle_hist5_un         = new TH1F("PtCalcEle5_un","PtCalc of unmatched electrons in 53X",1000,0,500);
	TH1F *PtCalcEle_hist7            = new TH1F("PtCalcEle7","PtCalc of electrons in 70X",1000,0,500);
	TH1F *DeltaPtCalcEle        = new TH1F("PtCalcDelta","PtCalc_5 - PtCalc-7 ",1000,-50,50);
	TH1F *FracPtCalcEle         = new TH1F("PtCalcChange","(PtCalc_5 - PtCalc-7 )/PtCalc_5" ,1000,-1,1);
	PtCalcEle_hist5          ->Sumw2(); 
	PtCalcEle_hist5_un       ->Sumw2(); 
	PtCalcEle_hist7          ->Sumw2(); 
	DeltaPtCalcEle      ->Sumw2(); 
	FracPtCalcEle       ->Sumw2(); 

	TH1F *PtEle_hist5            = new TH1F("PtEle5","Pt of electrons in 53X",1000,0,500);
	TH1F *PtEle_hist5_un         = new TH1F("PtEle5_un","Pt of unmatched electrons in 53X",1000,0,500);
	TH1F *PtEle_hist7            = new TH1F("PtEle7","Pt of electrons in 70X",1000,0,500);
	TH1F *DeltaPtEle        = new TH1F("PtDelta","Pt_5 - Pt-7 ",1000,-50,50);
	TH1F *FracPtEle         = new TH1F("PtChange","(Pt_5 - Pt-7 )/Pt_5" ,1000,-1,1);
	PtEle_hist5          ->Sumw2(); 
	PtEle_hist5_un       ->Sumw2(); 
	PtEle_hist7          ->Sumw2(); 
	DeltaPtEle      ->Sumw2(); 
	FracPtEle       ->Sumw2(); 

	TH1F *FBremEle5         = new TH1F("FBremEle5","FBrem of electrons in 53X",1000,0,1);
	TH1F *FBremEle5_un      = new TH1F("FBremEle5_un","FBrem of unmatched electrons in 53X",1000,0,1);
	TH1F *FBremEle7         =  new TH1F("FBremEle7","FBrem of electrons in 70X",1000,0,1);
	TH1F *DeltaFBremEle     = new TH1F("FBremDelta","FBrem_5 - FBrem-7 ",1000,-0.1,0.1);
	TH1F *FracFBremEle      = new TH1F("FBremChange","(FBrem_5 - FBrem-7 )/FBrem_5" ,1000,-0.1,0.1);
	FBremEle5       ->Sumw2(); 
	FBremEle5_un    ->Sumw2(); 
	FBremEle7       ->Sumw2(); 
	DeltaFBremEle   ->Sumw2(); 
	FracFBremEle    ->Sumw2(); 

	TH1F *PhiTrkSCEle5      = new TH1F("PhiTrkSCEle5","PhiTrkSC of electrons in 53X",200,-0.1,0.1);
	TH1F *PhiTrkSCEle5_un   = new TH1F("PhiTrkSCEle5_un","PhiTrkSC of unmatched electrons in 53X",200,-0.1,0.1);
	TH1F *PhiTrkSCEle7      = new TH1F("PhiTrkSCEle7","PhiTrkSC of electrons in 70X",200,-0.1,0.1);
	TH1F *DeltaPhiTrkSCEle  = new TH1F("PhiTrkSCDelta","PhiTrkSC_5 - PhiTrkSC-7 ",1000,-0.02,0.02);
	PhiTrkSCEle5    ->Sumw2(); 
	PhiTrkSCEle5_un ->Sumw2(); 
	PhiTrkSCEle7    ->Sumw2(); 
	DeltaPhiTrkSCEle->Sumw2(); 

	TH1F *EtaTrkSCEle5      = new TH1F("EtaTrkSCEle5","EtaTrkSC of electrons in 53X",200,-0.1,0.1);
	TH1F *EtaTrkSCEle5_un   = new TH1F("EtaTrkSCEle5_un","EtaTrkSC of unmatched electrons in 53X",200,-0.1,0.1);
	TH1F *EtaTrkSCEle7      = new TH1F("EtaTrkSCEle7","EtaTrkSC of electrons in 70X",200,-0.1,0.1);
	TH1F *DeltaEtaTrkSCEle  = new TH1F("EtaTrkSCDelta","EtaTrkSC_5 - EtaTrkSC-7 ",1000,-0.02,0.02);
	EtaTrkSCEle5    ->Sumw2(); 
	EtaTrkSCEle5_un ->Sumw2(); 
	EtaTrkSCEle7    ->Sumw2(); 
	DeltaEtaTrkSCEle->Sumw2(); 

	TH1F *MeeEle5           = new TH1F("MeeEle5","Mee of electrons in 53X",1000,0,200);
	TH1F *MeeEle5_un        = new TH1F("MeeEle5_un","Mee of unmatched electrons in 53X",1000,0,200);
	TH1F *MeeEle7           = new TH1F("MeeEle7","Mee of electrons in 70X",1000,0,200);
	TH1F *DeltaMeeEle       = new TH1F("MeeDelta","Mee_5 - Mee-7 ",5000,-10,10);
	TH1F *FracMeeEle       = new TH1F("MeeFrac","Mee_5 - Mee-7 ",5000,-1,1);
	MeeEle5         ->Sumw2(); 
	MeeEle5_un      ->Sumw2(); 
	MeeEle7         ->Sumw2(); 
	DeltaMeeEle     ->Sumw2(); 
	FracMeeEle     ->Sumw2(); 

	//	TH1F *R9_hist_5         = new TH1F("R9Ele5","R9Ele5",1000,0,1.2);
	TH1F *R9_anomaly_5      = new TH1F("eta_around_r9_0.5","eta aorund r9 0.5",1000,-5,5);
	TH1F *R9_anomaly_7     = new TH1F("eta_around_r9_0.5_70x","eta aorund r9 0.5 for 70x",1000,-5,5);
	R9_anomaly_5->Sumw2();
	R9_anomaly_7->Sumw2();
	//	TH1F *R9_hist_7         = new TH1F("R9Ele7","R9Ele7",1000,0,1.2);


	TH2F *eta_v_phi5_un     = new TH2F("evp5_un","eta v phi unmatched in 53x",100,-3,3,100,-7,7);
	eta_v_phi5_un   ->Sumw2(); 
	TH2F *x_v_y5_un_EEp     = new TH2F("xvy5_un_EEp","X v Y unmatched in 53x EEp",600,-300,300,410,-10,400);
	TH2F *x_v_y5_un_EEm     = new TH2F("xvy5_un_EEm","X v Y unmatched in 53x EEm",600,-300,300,410,-10,400);
	TH2F *x_v_y5_un_EB    = new TH2F("xvy5_un_EB","X v Y unmatched in 53x EB",600,-300,300,410,-10,400);

	//Cosmetics	
	R9_anomaly_5->SetLineColor(kRed);
	PtEle_hist5         ->SetLineColor(kRed);
	PtCalcEle_hist5         ->SetLineColor(kRed);
	FBremEle5      ->SetLineColor(kRed);
	PhiTrkSCEle5   ->SetLineColor(kRed);
	EtaTrkSCEle5   ->SetLineColor(kRed);
	MeeEle5        ->SetLineColor(kRed);
	R9_hist_5      ->SetLineColor(kRed);
	R9_hist_5->SetLineColor(kRed);

	//Declare list of files
	TFile **Files_53X;
	TFile **Files_70X;
	TFile **Files_Map;
	Files_53X = new TFile*[4];
	Files_70X = new TFile*[4];
	Files_Map = new TFile*[4];



	//Open Map and both data sources.
	Files_Map[0]  = TFile::Open("tmp/Map_d1-louieTest3.root"); //map provided my matching code
	//Files_70X[0] = TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2/190645-193621/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2-190645-193621.root");
	Files_70X[0] = TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2/190645-193621/clusteringValidationAND_06Oct14_70x_v1/DoubleElectron-ZHLTSkimPath-RUN2012A-15Apr-v2-190645-193621.root");
	//Files_53X[0]  = TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-22Jan-v1/190645-193621/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012A-22Jan-v1-190645-193621.root");
	Files_53X[0]  = TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkimPath-RUN2012A-22Jan-v1/190645-193621/clusteringValidationAND_06Oct14_53x_v3/DoubleElectron-ZSkimPath-RUN2012A-22Jan-v1-190645-193621.root");
	Files_Map[1]= TFile::Open("tmp/Map_d2-louieTest3.root");
	Files_70X[1]= TFile::Open("~lcorpe/eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012B-15Apr-v1/193834-196531/clusteringValidationAND_06Oct14_70x_v4/DoubleElectron-ZHLTSkimPath-RUN2012B-15Apr-v1-193834-196531.root");
	Files_53X[1] = TFile::Open("~lcorpe/eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkimPath-RUN2012B-22Jan-v1/193834-196531/clusteringValidationAND_06Oct14_53x_v4/DoubleElectron-ZSkimPath-RUN2012B-22Jan-v1-193834-196531.root");
	//Files_70X[1]= TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012B-15Apr-v1/193834-196531/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012B-15Apr-v1-193834-196531.root");
	//Files_53X[1] = TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkimPath-RUN2012B-22Jan-v1/193834-196531/clusteringValidationAND_06Oct14_53x_v4/DoubleElectron-ZSkimPath-RUN2012B-22Jan-v1-193834-196531.root");
	//Files_Map[2]= TFile::Open("tmp/Map_d3-louieTest3.root");
	//Files_70X[2]= TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012C-15Apr-v1/198022-203742/190456-208686-22Jan_v1/DoubleElectron-ZHLTSkimPath-RUN2012C-15Apr-v1-198022-203742.root");
	Files_Map[2]= TFile::Open("tmp/Map_d3-louieTest4.root");
	Files_70X[2]= TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012C-15Apr-v1/198022-203742/clusteringValidationAND_06Oct14_70x_v5/DoubleElectron-ZHLTSkimPath-RUN2012C-15Apr-v1-198022-203742_v2.root");
	Files_53X[2] = TFile::Open(" root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkimPath-RUN2012C-22Jan-v1/198022-203742/clusteringValidationAND_06Oct14_53x_v4/DoubleElectron-ZSkimPath-RUN2012C-22Jan-v1-198022-203742.root");
	Files_Map[3]= TFile::Open("tmp/Map_d4-louieTest3.root");
	Files_70X[3] = TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZHLTSkimPath-RUN2012D-15Apr-v1/203777-208686/clusteringValidationAND_06Oct14_70x_v4/DoubleElectron-ZHLTSkimPath-RUN2012D-15Apr-v1-203777-208686.root");
	Files_53X[3]= TFile::Open("root://eoscms///eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkimPath-RUN2012D-22Jan-v1/203777-208686/clusteringValidationAND_06Oct14_53x_v4/DoubleElectron-ZSkimPath-RUN2012D-22Jan-v1-203777-208686.root");
	//LOOPS
	for(Int_t overLoop=1;overLoop<0;overLoop++) //in theory can run over all ABCD at once by running overLoop 0 to <4. In practice, easier to cherry-pick runs by changing this lines.
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
		Int_t runNumber7;	
		Int_t lumiNumber7;	
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
		Int_t runNumber5;	
		Int_t lumiNumber5;	
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
		TBranch *b_runNumber7;	
		TBranch *b_lumiNumber7;	
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
		TBranch *b_runNumber5;	
		TBranch *b_lumiNumber5;	
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
		chain5->SetBranchAddress("runNumber", &runNumber5,&b_runNumber5);
		chain5->SetBranchAddress("lumiBlock", &lumiNumber5,&b_lumiNumber5);
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
		chain7->SetBranchAddress("runNumber", &runNumber7,&b_runNumber7);
		chain7->SetBranchAddress("lumiBlock", &lumiNumber7,&b_lumiNumber7);
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

		int debug =0;
		int selectionDebug =0;
		int totFillCounter =0;
		int oneMatchCounter =0;
		int twoMatchCounter =0;
		int noMatchCounter =0;
		int DeN_error =0;

		//Loops over first TTree	
		for(Int_t loop = 0; loop < 

				chain5->GetEntries()
				//	10000	
				; loop++)
		{
			//Get entry in first TTree and corresponding Map entry	
			tree_Map->GetEntry(loop);
			chain5->GetEntry(loop);

			// Progress tracker
			if (loop%10000==0) 
			{ //prints out info as a sort of progress tracker. 
				std::cout << loop << " / " << chain5->GetEntries() << std::endl;
				std:: cout << "missing debug " << debug << std::endl;
				std:: cout << "selection  debug " << selectionDebug << std::endl;
				std::cout << " total entries: " << chain5->GetEntries() << std::endl;
				std::cout << "total fill " << totFillCounter << std::endl;
				std::cout << " two match " << twoMatchCounter << " / " << totFillCounter << " = " << (float) 100* twoMatchCounter/totFillCounter << " \%" << std::endl;
				std::cout << " one match " << oneMatchCounter << " / " << totFillCounter << " = " << (float) 100* oneMatchCounter/totFillCounter << " \%" << std::endl;
				std::cout << " no match " << noMatchCounter << " / " << totFillCounter << " = " << (float) 100* noMatchCounter/totFillCounter << " \%" << std::endl;
			}

			//Get corresponding entry in TTree 2  using information from Map
			chain7->GetEntry(entryNumber);

			//set selection value for events in 5 
			selection5=((eleID5[0] & 2)==2)*((eleID5[1] & 2)==2)*(HLTfire5==1)*(recoFlagsEle5[0] > 1)*(recoFlagsEle5[1] > 1)*(PtEle5[0]>20)*(PtEle5[1]>20);
			total_weight5=1;
			total_weight5*=weight5*r9weight5[0]*r9weight5[1]*ptweight5[0]*ptweight5[1];//*mcGenWeight;
			selection5*=total_weight5;
			float PtCalcEle5[2];
			float PtCalcEle7[2];

			//Calcuated version of Pt (or more precisely, Et)
			PtCalcEle5[0]=( energySCEle_corr5[0]/ std::cosh(etaEle5[0]));
			PtCalcEle5[1]=( energySCEle_corr5[1]/ std::cosh(etaEle5[1]));
			PtCalcEle7[0]=( energySCEle_corr7[0]/ std::cosh(etaEle7[0]));
			PtCalcEle7[1]=( energySCEle_corr7[1]/ std::cosh(etaEle7[1]));

			if (selection5 ==0) {selectionDebug++;continue;} //reject event which do not pass selection on 5
			totFillCounter++; 


			if (entryNumber ==-1 ) //if no match in map found
			{
				debug++;
				noMatchCounter++;
				if (noMatchCounter%1000==0)
				{
					//std::cout << "[DEBUG] run"<<  runNumber5 << ", lumi " << lumiNumber5 << std::endl;
				}
			}
			else
			{

				if((eventNumber5 - eventNumber7) != 0) //sanity check that the event number match!
				{
					std::cout << "[ERROR] Delta_eN = "<< eventNumber5 - eventNumber7<< ", NOT FILLING" << std::endl;
					DeN_error++;//Delta eventNumber error
					continue;
				}
				if (eleIndex[0] ==-1 || eleIndex[1] ==-1) {oneMatchCounter++;} // count cases with one matching electron
				if (eleIndex[0] != -1 && eleIndex[1] !=-1) {twoMatchCounter++;} //count cases with two matchign electrons


				if(eleIndex[0] != -1 && eleIndex[1] != -1) //if both electrons have match...
				{

					if(R9Ele5[0] > 0.45 && R9Ele5[0] <0.55)
					{R9_anomaly_5->Fill(etaEle5[0]);}
					if(R9Ele5[1] > 0.45 && R9Ele5[1] <0.55)
					{R9_anomaly_5->Fill(etaEle5[1]);}
					if(R9Ele7[0] > 0.45 && R9Ele7[0] <0.55)
					{R9_anomaly_7->Fill(etaEle7[0]);}
					if(R9Ele7[1] > 0.45 && R9Ele7[1] <0.55)
					{R9_anomaly_7->Fill(etaEle7[1]);}
				
				//Fill basic histograms
					if (eleIndex[0]==1)  //correct order
					{

						Delta_eN->Fill(eventNumber5 - eventNumber7);
						Delta_rN->Fill(runNumber5 - runNumber7);
						PtEle_hist5->Fill(PtEle5[0]); 
						PtEle_hist5->Fill(PtEle5[1]); 
						PtEle_hist7->Fill(PtEle7[0]);        
						PtEle_hist7->Fill(PtEle7[1]);        
						DeltaPtEle   ->Fill(PtEle5[0] - PtEle7[0]); 
						DeltaPtEle   ->Fill(PtEle5[1] - PtEle7[1]); 
						FracPtEle    ->Fill( (PtEle5[0] - PtEle7[0])/PtEle5[0]);
						FracPtEle    ->Fill( (PtEle5[1] - PtEle7[1])/PtEle5[1]);
						PtCalcEle_hist5->Fill(PtCalcEle5[0]); 
						PtCalcEle_hist5->Fill(PtCalcEle5[1]); 
						PtCalcEle_hist7->Fill(PtCalcEle7[0]);        
						PtCalcEle_hist7->Fill(PtCalcEle7[1]);        
						DeltaPtCalcEle   ->Fill(PtCalcEle5[0] - PtCalcEle7[0]); 
						DeltaPtCalcEle   ->Fill(PtCalcEle5[1] - PtCalcEle7[1]); 
						FracPtCalcEle    ->Fill( (PtCalcEle5[0] - PtCalcEle7[0])/PtCalcEle5[0]);
						FracPtCalcEle    ->Fill( (PtCalcEle5[1] - PtCalcEle7[1])/PtCalcEle5[1]);
						FBremEle5->Fill(fbremEle5[0]); 
						FBremEle5->Fill(fbremEle5[1]); 
						FBremEle7->Fill(fbremEle7[0]);        
						FBremEle7->Fill(fbremEle7[1]);        
						DeltaFBremEle   ->Fill(fbremEle5[0] - fbremEle7[0]); 
						DeltaFBremEle   ->Fill(fbremEle5[1] - fbremEle7[1]); 
						FracFBremEle    ->Fill( (fbremEle5[0] - fbremEle7[0])/fbremEle5[0]);
						FracFBremEle    ->Fill( (fbremEle5[1] - fbremEle7[1])/fbremEle5[1]);
						PhiTrkSCEle5 ->Fill(phiEle5[0]-phiSCEle5[0]); 
						PhiTrkSCEle5 ->Fill(phiEle5[1]-phiSCEle5[1]); 
						PhiTrkSCEle7  ->Fill(phiEle7[0]-phiSCEle7[0]); 
						PhiTrkSCEle7  ->Fill(phiEle7[1]-phiSCEle7[1]); 
						DeltaPhiTrkSCEle ->Fill((phiEle5[0]-phiSCEle5[0])-(phiEle7[0]-phiSCEle7[0]));
						DeltaPhiTrkSCEle ->Fill((phiEle5[1]-phiSCEle5[1])-(phiEle7[1]-phiSCEle7[1]));
						EtaTrkSCEle5 ->Fill(etaEle5[0]-etaSCEle5[0]); 
						EtaTrkSCEle5 ->Fill(etaEle5[1]-etaSCEle5[1]); 
						EtaTrkSCEle7  ->Fill(etaEle7[0]-etaSCEle7[0]); 
						EtaTrkSCEle7  ->Fill(etaEle7[1]-etaSCEle7[1]); 
						DeltaEtaTrkSCEle ->Fill((etaEle5[0]-etaSCEle5[0])-(etaEle7[0]-etaSCEle7[0]));
						DeltaEtaTrkSCEle ->Fill((etaEle5[1]-etaSCEle5[1])-(etaEle7[1]-etaSCEle7[1]));
						MeeEle5 ->     Fill(invMass_corr5); 
						MeeEle7 -> Fill(invMass_corr7);      
						DeltaMeeEle   ->Fill(invMass_corr7 - invMass_corr5);
						FracMeeEle   ->Fill((invMass_corr7 - invMass_corr5)/invMass_corr5);


						R9_hist_7->Fill((R9Ele7[0]),selection5);
						R9_hist_7->Fill((R9Ele7[1]),selection5);
						R9_hist_5->Fill((R9Ele5[0]),selection5);
						R9_hist_5->Fill((R9Ele5[1]),selection5);
						R9_hist->Fill((R9Ele7[0]-R9Ele5[0]),selection5);
						R9_hist->Fill((R9Ele7[1]-R9Ele5[1]),selection5);
						R9_frac_hist->Fill((R9Ele7[0]-R9Ele5[0])/R9Ele5[0],selection5);
						R9_frac_hist->Fill((R9Ele7[1]-R9Ele5[1])/R9Ele5[0],selection5);
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
						fbrem->Fill(etaSCEle5[0],fbremEle7[0]/fbremEle5[0],selection5);
						fbrem->Fill(etaSCEle5[1],fbremEle7[1]/fbremEle5[1],selection5);
						corrE_profile->Fill(etaSCEle5[1],energySCEle_corr7[1]/energySCEle_corr5[1],selection5);
						corrE_profile->Fill(etaSCEle5[0],energySCEle_corr7[0]/energySCEle_corr5[0],selection5);
						R9_profile->Fill(etaSCEle5[0],(R9Ele5[0]/R9Ele7[0]),selection5);
						R9_profile->Fill(etaSCEle5[1],(R9Ele5[1]/R9Ele7[1]),selection5);
					}

					if (eleIndex[0]==2) //inverted order
					{


						Delta_eN->Fill(eventNumber5 - eventNumber7);
						Delta_rN->Fill(runNumber5 - runNumber7);
						PtEle_hist5->Fill(PtEle5[0]); 
						PtEle_hist5->Fill(PtEle5[1]); 
						PtEle_hist7->Fill(PtEle7[0]);        
						PtEle_hist7->Fill(PtEle7[1]);        
						DeltaPtEle   ->Fill(PtEle5[1] - PtEle7[0]); 
						DeltaPtEle   ->Fill(PtEle5[0] - PtEle7[1]); 
						FracPtEle    ->Fill( (PtEle5[1] - PtEle7[0])/PtEle5[1]);
						FracPtEle    ->Fill( (PtEle5[0] - PtEle7[1])/PtEle5[0]);
						PtCalcEle_hist5->Fill(PtCalcEle5[0]); 
						PtCalcEle_hist5->Fill(PtCalcEle5[1]); 
						PtCalcEle_hist7->Fill(PtCalcEle7[0]);        
						PtCalcEle_hist7->Fill(PtCalcEle7[1]);        
						DeltaPtCalcEle   ->Fill(PtCalcEle5[1] - PtCalcEle7[0]); 
						DeltaPtCalcEle   ->Fill(PtCalcEle5[0] - PtCalcEle7[1]); 
						FracPtCalcEle    ->Fill( (PtCalcEle5[1] - PtCalcEle7[0])/PtCalcEle5[1]);
						FracPtCalcEle    ->Fill( (PtCalcEle5[0] - PtCalcEle7[1])/PtCalcEle5[0]);
						FBremEle5->Fill(fbremEle5[0]); 
						FBremEle5->Fill(fbremEle5[1]); 
						FBremEle7->Fill(fbremEle7[0]);        
						FBremEle7->Fill(fbremEle7[1]);        
						DeltaFBremEle   ->Fill(fbremEle5[1] - fbremEle7[0]); 
						DeltaFBremEle   ->Fill(fbremEle5[0] - fbremEle7[1]); 
						FracFBremEle    ->Fill( (fbremEle5[0] - fbremEle7[1])/fbremEle5[0]);
						FracFBremEle    ->Fill( (fbremEle5[1] - fbremEle7[0])/fbremEle5[1]);
						PhiTrkSCEle5 ->Fill(phiEle5[0]-phiSCEle5[0]); 
						PhiTrkSCEle5 ->Fill(phiEle5[1]-phiSCEle5[1]); 
						PhiTrkSCEle7  ->Fill(phiEle7[1]-phiSCEle7[1]); 
						PhiTrkSCEle7  ->Fill(phiEle7[0]-phiSCEle7[0]); 
						DeltaPhiTrkSCEle ->Fill((phiEle5[1]-phiSCEle5[1])-(phiEle7[0]-phiSCEle7[0]));
						DeltaPhiTrkSCEle ->Fill((phiEle5[0]-phiSCEle5[0])-(phiEle7[1]-phiSCEle7[1]));
						EtaTrkSCEle5 ->Fill(etaEle5[0]-etaSCEle5[0]); 
						EtaTrkSCEle5 ->Fill(etaEle5[1]-etaSCEle5[1]); 
						EtaTrkSCEle7  ->Fill(etaEle7[0]-etaSCEle7[0]); 
						EtaTrkSCEle7  ->Fill(etaEle7[1]-etaSCEle7[1]); 
						DeltaEtaTrkSCEle ->Fill((etaEle5[1]-etaSCEle5[1])-(etaEle7[0]-etaSCEle7[0]));
						DeltaEtaTrkSCEle ->Fill((etaEle5[0]-etaSCEle5[0])-(etaEle7[1]-etaSCEle7[1]));
						MeeEle5 ->     Fill(invMass_corr5); 
						MeeEle7 -> Fill(invMass_corr7);      
						DeltaMeeEle   ->Fill(invMass_corr7 - invMass_corr5);
						FracMeeEle   ->Fill((invMass_corr7 - invMass_corr5)/invMass_corr5);
						R9_hist_7->Fill((R9Ele7[0]),selection5);
						R9_hist_7->Fill((R9Ele7[1]),selection5);
						R9_hist_5->Fill((R9Ele5[0]),selection5);
						R9_hist_5->Fill((R9Ele5[1]),selection5);
						R9_hist->Fill((R9Ele7[1]-R9Ele5[0]),selection5);
						R9_hist->Fill((R9Ele7[0]-R9Ele5[1]),selection5);
						R9_frac_hist->Fill((R9Ele7[1]-R9Ele5[0])/R9Ele5[0],selection5);
						R9_frac_hist->Fill((R9Ele7[0]-R9Ele5[1])/R9Ele5[1],selection5);
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
						fbrem->Fill(etaSCEle5[0],fbremEle7[1]/fbremEle5[0],selection5);
						fbrem->Fill(etaSCEle5[1],fbremEle7[0]/fbremEle5[1],selection5);
						corrE_profile->Fill(etaSCEle5[1],energySCEle_corr7[0]/energySCEle_corr5[1],selection5);
						corrE_profile->Fill(etaSCEle5[0],energySCEle_corr7[1]/energySCEle_corr5[0],selection5);
						R9_profile->Fill(etaSCEle5[0],(R9Ele5[0]/R9Ele7[1]),selection5);
						R9_profile->Fill(etaSCEle5[1],(R9Ele5[1]/R9Ele7[0]),selection5);
					}
				}
				if(eleIndex[0] == -1 && eleIndex[1] != -1) //one missing electron (0th electron missing)
				{


					PtEle_hist5_un->Fill(PtEle5[0]);     
					PtCalcEle_hist5_un->Fill(PtCalcEle5[0]);     
					FBremEle5_un->  Fill(fbremEle5[0]); 

					PhiTrkSCEle5_un ->Fill((phiEle5[0]-phiSCEle5[0]));
					EtaTrkSCEle5_un ->Fill((etaEle5[0]-etaSCEle5[0]));
					eta_v_phi5_un ->Fill(etaSCEle5[0], phiSCEle5[0]);
					if(etaSCEle5[0] < -1.556){
						x_v_y5_un_EEm ->Fill(seedXSCEle5[0], seedYSCEle5[0]);	
					}
					if(etaSCEle5[0]>1.556){
						x_v_y5_un_EEp ->Fill(seedXSCEle5[0], seedYSCEle5[0]);	
					}
					if(fabs(etaSCEle5[0])<1.4442){
						x_v_y5_un_EB ->Fill(seedXSCEle5[0], seedYSCEle5[0]);	
					}
					MeeEle5_un    ->Fill(invMass_corr5);
					R9_hist_5_un  ->Fill(R9Ele5[0]);

					if (eleIndex[1] ==2) //correct order
					{
						if(R9Ele5[1] > 0.45 && R9Ele5[1] <0.55)
						{R9_anomaly_5->Fill(etaEle5[1]);}
						if(R9Ele7[1] > 0.45 && R9Ele7[1] <0.55)
						{R9_anomaly_7->Fill(etaEle7[1]);}
						PtEle_hist5->Fill(PtEle5[1]); 
						PtEle_hist7->Fill(PtEle7[1]);        
						PtCalcEle_hist5->Fill(PtCalcEle5[1]); 
						PtCalcEle_hist7->Fill(PtCalcEle7[1]);        
						DeltaPtEle ->Fill(PtEle5[1] - PtEle7[1]); 
						FracPtEle    ->Fill( (PtEle5[1] - PtEle7[1])/PtEle5[1]);
						DeltaPtCalcEle ->Fill(PtCalcEle5[1] - PtCalcEle7[1]); 
						FracPtCalcEle    ->Fill( (PtCalcEle5[1] - PtCalcEle7[1])/PtCalcEle5[1]);
						FBremEle5->Fill(fbremEle5[1]); 
						FBremEle7->Fill(fbremEle7[1]);        
						DeltaFBremEle   ->Fill(fbremEle5[1] - fbremEle7[1]); 
						FracFBremEle    ->Fill( (fbremEle5[1] - fbremEle7[1])/fbremEle5[1]);
						PhiTrkSCEle5 ->Fill((phiEle5[1]-phiSCEle5[1])); 
						PhiTrkSCEle7  ->Fill((phiEle7[1]-phiSCEle7[1])); 
						DeltaPhiTrkSCEle ->Fill((phiEle5[1]-phiSCEle5[1])-(phiEle7[1]-phiSCEle7[1]));
						EtaTrkSCEle5 ->Fill((etaEle5[1]-etaSCEle5[1])); 
						EtaTrkSCEle7  ->Fill((etaEle7[1]-etaSCEle7[1])); 
						DeltaEtaTrkSCEle ->Fill((etaEle5[1]-etaSCEle5[1])-(etaEle7[1]-etaSCEle7[1]));
					}
					if (eleIndex[1] ==1) //reversed order
					{
						if(R9Ele5[1] > 0.45 && R9Ele5[1] <0.55)
						{R9_anomaly_5->Fill(etaEle5[1]);}
						if(R9Ele7[0] > 0.45 && R9Ele7[0] <0.55)
						{R9_anomaly_7->Fill(etaEle7[0]);}
						PtEle_hist5->Fill(PtEle5[1]); 
						PtEle_hist7->Fill(PtEle7[0]);        
						DeltaPtEle   ->Fill(PtEle5[1] - PtEle7[0]); 
						FracPtEle    ->Fill( (PtEle5[1] - PtEle7[0])/PtEle5[1]);
						PtCalcEle_hist5->Fill(PtCalcEle5[1]); 
						PtCalcEle_hist7->Fill(PtCalcEle7[0]);        
						DeltaPtCalcEle   ->Fill(PtCalcEle5[1] - PtCalcEle7[0]); 
						FracPtCalcEle    ->Fill( (PtCalcEle5[1] - PtCalcEle7[0])/PtCalcEle5[1]);
						FBremEle5->Fill(fbremEle5[1]); 
						FBremEle7->Fill(fbremEle7[0]);        
						DeltaFBremEle   ->Fill(fbremEle5[1] - fbremEle7[0]); 
						FracFBremEle    ->Fill( (fbremEle5[1] - fbremEle7[0])/fbremEle5[1]);
						PhiTrkSCEle5 ->Fill((phiEle5[1]-phiSCEle5[1])); 
						PhiTrkSCEle7  ->Fill((phiEle7[0]-phiSCEle7[0])); 
						DeltaPhiTrkSCEle ->Fill((phiEle5[1]-phiSCEle5[1])-(phiEle7[0]-phiSCEle7[0]));
						EtaTrkSCEle5 ->Fill((etaEle5[1]-etaSCEle5[1])); 
						EtaTrkSCEle7  ->Fill((etaEle7[0]-etaSCEle7[0])); 
						DeltaEtaTrkSCEle ->Fill((etaEle5[1]-etaSCEle5[1])-(etaEle7[0]-etaSCEle7[0]));
					}

				}


				if(eleIndex[1] == -1 && eleIndex[0] != -1) //electron 1 missing
				{


					PtEle_hist5_un->Fill(PtEle5[1]);     
					PtCalcEle_hist5_un->Fill(PtCalcEle5[1]);     
					FBremEle5_un->  Fill(fbremEle5[1]); 

					PhiTrkSCEle5_un ->Fill((phiEle5[1]-phiSCEle5[1]));
					EtaTrkSCEle5_un ->Fill((etaEle5[1]-etaSCEle5[1]));
					eta_v_phi5_un ->Fill(etaSCEle5[1],phiSCEle5[1]);	
					if(etaSCEle5[1]>1.556){
						x_v_y5_un_EEp ->Fill(seedXSCEle5[1],seedYSCEle5[1]);	
					}
					if(etaSCEle5[1] < -1.556){
						x_v_y5_un_EEm ->Fill(seedXSCEle5[1],seedYSCEle5[1]);	
					}
					if(fabs(etaSCEle5[1])<1.4442){
						x_v_y5_un_EB ->Fill(seedXSCEle5[1],seedYSCEle5[1]);	
					}
					MeeEle5_un    ->Fill(invMass_corr5);
					R9_hist_5_un  ->Fill(R9Ele5[1]);

					if (eleIndex[0] ==1) //regular order
					{
						if(R9Ele5[0] > 0.45 && R9Ele5[0] <0.55)
						{R9_anomaly_5->Fill(etaEle5[0]);}
						if(R9Ele7[0] > 0.45 && R9Ele7[0] <0.55)
						{R9_anomaly_7->Fill(etaEle7[0]);}
						PtEle_hist5->Fill(PtEle5[0]); 
						PtEle_hist7->Fill(PtEle7[0]);        
						DeltaPtEle   ->Fill(PtEle5[0] - PtEle7[0]); 
						FracPtEle    ->Fill( (PtEle5[0] - PtEle7[0])/PtEle5[0]);
						PtCalcEle_hist5->Fill(PtCalcEle5[0]); 
						PtCalcEle_hist7->Fill(PtCalcEle7[0]);        
						DeltaPtCalcEle   ->Fill(PtCalcEle5[0] - PtCalcEle7[0]); 
						FracPtCalcEle    ->Fill( (PtCalcEle5[0] - PtCalcEle7[0])/PtCalcEle5[0]);
						FBremEle5->Fill(fbremEle5[0]); 
						FBremEle7->Fill(fbremEle7[0]);        
						DeltaFBremEle   ->Fill(fbremEle5[0] - fbremEle7[0]); 
						FracFBremEle    ->Fill( (fbremEle5[0] - fbremEle7[0])/fbremEle5[0]);
						PhiTrkSCEle5 ->Fill((phiEle5[0]-phiSCEle5[0])); 
						PhiTrkSCEle7  ->Fill((phiEle7[0]-phiSCEle7[0])); 
						DeltaPhiTrkSCEle ->Fill((phiEle5[0]-phiSCEle5[0])-(phiEle7[0]-phiSCEle7[0]));
						EtaTrkSCEle5 ->Fill((etaEle5[0]-etaSCEle5[0])); 
						EtaTrkSCEle7  ->Fill((etaEle7[0]-etaSCEle7[0])); 
						DeltaEtaTrkSCEle ->Fill((etaEle5[0]-etaSCEle5[0])-(etaEle7[0]-etaSCEle7[0]));
					}
					if (eleIndex[0] ==2) //reverse order
					{
						if(R9Ele5[0] > 0.45 && R9Ele5[0] <0.55)
						{R9_anomaly_5->Fill(etaEle5[0]);}
						if(R9Ele7[1] > 0.45 && R9Ele7[1] <0.55)
						{R9_anomaly_7->Fill(etaEle7[1]);}
						PtEle_hist5->Fill(PtEle5[0]); 
						PtEle_hist7->Fill(PtEle7[1]);        
						DeltaPtEle   ->Fill(PtEle5[0] - PtEle7[1]); 
						FracPtEle    ->Fill( (PtEle5[0] - PtEle7[1])/PtEle5[0]);
						PtCalcEle_hist5->Fill(PtCalcEle5[0]); 
						PtCalcEle_hist7->Fill(PtCalcEle7[1]);        
						DeltaPtCalcEle   ->Fill(PtCalcEle5[0] - PtCalcEle7[1]); 
						FracPtCalcEle    ->Fill( (PtCalcEle5[0] - PtCalcEle7[1])/PtCalcEle5[0]);
						FBremEle5->Fill(fbremEle5[0]); 
						FBremEle7->Fill(fbremEle7[1]);        
						DeltaFBremEle   ->Fill(fbremEle5[0] - fbremEle7[1]); 
						FracFBremEle    ->Fill( (fbremEle5[0] - fbremEle7[1])/fbremEle5[0]);
						PhiTrkSCEle5 ->Fill((phiEle5[0]-phiSCEle5[0])); 
						PhiTrkSCEle7  ->Fill((phiEle7[1]-phiSCEle7[1])); 
						DeltaPhiTrkSCEle ->Fill((phiEle5[0]-phiSCEle5[0])-(phiEle7[1]-phiSCEle7[1]));
						EtaTrkSCEle5 ->Fill((etaEle5[0]-etaSCEle5[0])); 
						EtaTrkSCEle7  ->Fill((etaEle7[1]-etaSCEle7[1])); 
						DeltaEtaTrkSCEle ->Fill((etaEle5[0]-etaSCEle5[0])-(etaEle7[1]-etaSCEle7[1]));
					}

				}

				if(eleIndex[1] == -1 && eleIndex[0] == -1) //both electrons missing
				{
					
					// plot properties of _un ( =unmatched) electrons in 53X
					PtEle_hist5_un->Fill(PtEle5[0]);     
					PtCalcEle_hist5_un->Fill(PtCalcEle5[0]);     
					FBremEle5_un->  Fill(fbremEle5[0]); 

					PhiTrkSCEle5_un ->Fill((phiEle5[0]-phiSCEle5[0]));
					EtaTrkSCEle5_un ->Fill((etaEle5[0]-etaSCEle5[0]));
					eta_v_phi5_un ->Fill(etaSCEle5[0],phiSCEle5[0]);	
					//	MeeEle5_un    ->Fill(invMass_corr5);
					R9_hist_5_un  ->Fill(R9Ele5[0]);

					PtEle_hist5_un->Fill(PtEle5[1]);     
					PtCalcEle_hist5_un->Fill(PtCalcEle5[1]);     
					FBremEle5_un->  Fill(fbremEle5[1]); 

					PhiTrkSCEle5_un ->Fill((phiEle5[1]-phiSCEle5[1]));
					EtaTrkSCEle5_un ->Fill((etaEle5[1]-etaSCEle5[1]));
					eta_v_phi5_un ->Fill(etaSCEle5[1],phiSCEle5[1]);	
					if(etaSCEle5[0]>1.556){
						x_v_y5_un_EEp ->Fill(seedXSCEle5[0],seedYSCEle5[0]);	
					}

					if(etaSCEle5[1]>1.556){
						x_v_y5_un_EEp ->Fill(seedXSCEle5[1],seedYSCEle5[1]);
					}
					if(etaSCEle5[0] < -1.556){
						x_v_y5_un_EEm ->Fill(seedXSCEle5[0],seedYSCEle5[0]);
					}
					if(etaSCEle5[1] < -1.556){
						x_v_y5_un_EEm ->Fill(seedXSCEle5[1],seedYSCEle5[1]);
					}
					if(fabs(etaSCEle5[0])<1.4442){
						x_v_y5_un_EB ->Fill(seedXSCEle5[0],seedYSCEle5[0]);
					}
					if(fabs(etaSCEle5[1])<1.4442){
						x_v_y5_un_EB ->Fill(seedXSCEle5[1],seedYSCEle5[1]);
					}
					MeeEle5_un    ->Fill(invMass_corr5);
					R9_hist_5_un  ->Fill(R9Ele5[1]);
				}

			}//end of loop	
		}
		std:: cout << "missing debug " << debug << std::endl; //summary of errors and matches
		//std:: cout << "lumi debug " << lumiDebug << std::endl;
		std:: cout << "selection  debug " << selectionDebug << std::endl;
		std::cout << " total entries: " << chain5->GetEntries() << std::endl;
		std::cout << "missing % : " << ((float) debug/chain5->GetEntries()) *100 << " \% " << std::endl;

		std::cout << std::endl;
		std::cout << "total fill " << totFillCounter << std::endl;
		std::cout << " two match " << twoMatchCounter << " / " << totFillCounter << " = " << (float) 100* twoMatchCounter/totFillCounter << " \%" << std::endl;
		std::cout << " one match " << oneMatchCounter << " / " << totFillCounter << " = " << (float) 100* oneMatchCounter/totFillCounter << " \%" << std::endl;
		std::cout << " no match " << noMatchCounter << " / " << totFillCounter << " = " << (float) 100* noMatchCounter/totFillCounter << " \%" << std::endl;
		std::cout << " Delta_eN errors " << DeN_error << " / " << totFillCounter << " = " << (float) 100* DeN_error/totFillCounter << " \%" << std::endl;

		//end of overLoop
	}
	//Draw basic Histos

	gStyle->SetOptStat(0);


	TFile f("fullRunPlots_OctA/histosB.root","recreate");	
	R9_anomaly_5->Write();
	R9_anomaly_7->Write();
	PtEle_hist5         ->Write(); 
	PtEle_hist5_un       ->Write(); 
	PtEle_hist7          ->Write(); 
	DeltaPtEle      ->Write(); 
	FracPtEle       ->Write(); 
	PtCalcEle_hist5         ->Write(); 
	PtCalcEle_hist5_un       ->Write(); 
	PtCalcEle_hist7          ->Write(); 
	DeltaPtCalcEle      ->Write(); 
	FracPtCalcEle       ->Write(); 

	FBremEle5       ->Write(); 
	FBremEle5_un    ->Write(); 
	FBremEle7       ->Write(); 
	DeltaFBremEle   ->Write(); 
	FracFBremEle    ->Write(); 

	PhiTrkSCEle5    ->Write(); 
	PhiTrkSCEle5_un ->Write(); 
	PhiTrkSCEle7    ->Write(); 
	DeltaPhiTrkSCEle->Write(); 

	EtaTrkSCEle5    ->Write(); 
	EtaTrkSCEle5_un ->Write(); 
	EtaTrkSCEle7    ->Write(); 
	DeltaEtaTrkSCEle->Write(); 

	MeeEle5         ->Write(); 
	MeeEle5_un      ->Write(); 
	MeeEle7         ->Write(); 
	DeltaMeeEle     ->Write(); 
	FracMeeEle     ->Write(); 

	R9_hist_5       ->Write(); 
	R9_hist_5_un    ->Write(); 
	R9_hist_7       ->Write(); 
	R9_hist     ->Write(); 
	R9_frac_hist     ->Write(); 
	Delta_eN->Write();
	Delta_rN->Write();

	eta_v_phi5_un   ->Write(); 
	x_v_y5_un_EEp   ->Write(); 
	x_v_y5_un_EEm   ->Write(); 
	x_v_y5_un_EB   ->Write(); 
	TCanvas *t = new TCanvas("t","t",600,600);

	TPaveText *pt2 =new TPaveText(0.65,0.7,0.85,0.85,"blNDC");
	pt2->SetFillColor(0);
	pt2->SetTextAlign(11);
	pt2->SetBorderSize(0);


	R9_profile->SetMinimum(0.8);
	R9_profile->SetMaximum(1.2);
	R9_profile->Draw();
	R9_profile->Write();
	R9_profile->GetXaxis()->SetTitle("#eta");
	R9_profile->GetYaxis()->SetTitle("R9_{7}/R9_{5}");
	t->SaveAs("fullRunPlots_OctA/R9_profile.pdf");

	etaSC_hist->Draw("c");
	etaSC_hist->Write();
	pt2->Clear();
	effSigma(etaSC_hist);
	pt2->AddText(effSigma_cstr);
	pt2->Draw();
	etaSC_hist->GetXaxis()->SetTitle("#Delta #eta_{SC}");
	t->SaveAs("fullRunPlots_OctA/etaSC.pdf");

	phiSC_hist->Draw("C");
	phiSC_hist->Write();
	pt2->Clear();
	effSigma(phiSC_hist);
	pt2->AddText(effSigma_cstr);
	pt2->Draw();
	phiSC_hist->GetXaxis()->SetTitle("#Delta #phi_{SC}");
	t->SaveAs("fullRunPlots_OctA/phiSC.pdf");

	eta_hist->Draw("C");
	eta_hist->Write();
	pt2->Clear();
	effSigma(eta_hist);
	pt2->AddText(effSigma_cstr);
	pt2->Draw();
	eta_hist->GetXaxis()->SetTitle("#Delta #eta");
	t->SaveAs("fullRunPlots_OctA/eta.pdf");

	phi_hist->Draw("C");
	phi_hist->Write();
	pt2->Clear();
	effSigma(phi_hist);
	pt2->AddText(effSigma_cstr);
	pt2->Draw();
	phi_hist->GetXaxis()->SetTitle("#Delta #phi");
	t->SaveAs("fullRunPlots_OctA/phi.pdf");

	pModeGsf_hist->Draw();
	pModeGsf_hist->Write();
	pt2->Clear();
	effSigma(pModeGsf_hist);
	pt2->AddText(effSigma_cstr);
	pt2->Draw();

	pModeGsf_hist->GetXaxis()->SetTitle("#Delta pModeGsf");
	t->SaveAs("fullRunPlots_OctA/pModeGsf.pdf");
	//	fbrem->SetMimimum(0.8);
	//	fbrem->SetMaximum(1.2);
	fbrem->Draw();
	fbrem->Write();
	fbrem->GetXaxis()->SetTitle("#eta");
	fbrem->GetYaxis()->SetTitle("f_{brem7}/f_{brem5}");
	t->SaveAs("fullRunPlots_OctA/fbrem2.pdf");
	return 0;
}
