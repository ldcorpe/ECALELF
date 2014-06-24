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
#include "TGaxis.h"


// EEFECTIVE SIGMA FUNCYION DEFINITION //

const char* gausSigma_cstr;
double gausSigma_val;
Double_t gausSigma(TH1* h)
{
using namespace std;
	Double_t ave = h->GetMean();
	TF1 *g1= new TF1("g1","gaus",ave*0.9,ave*1.1);
	h->Fit("g1","N");

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
//	cout << "#sigma_{eff} = " <<  std::fixed <<  std::setprecision(3) <<  widmin << endl;
	std::string effSigma_str = ostr1.str();
	effSigma_cstr = effSigma_str.c_str();
//	cout << effSigma_str << endl <<effSigma_cstr << endl;
	effSigma_val=widmin;
	return widmin;

}


int debug_matching() 
{
//gROOT->SetBatch();
TGaxis::SetMaxDigits(2);
	//DEFINE Histograms
	//Simple quantitiy comaprisons


TH1F *eta1 = new TH1F("eta1","eta1",170,-5,5);
TH1F *phi1 = new TH1F("phi1","phi1",360,-6.3,6.3);
TH1F *pt1 = new TH1F("pt1","pt1",100,0,300);
TH1F *eta2 = new TH1F("eta2","eta2",170,-5,5);
TH1F *phi2 = new TH1F("phi2","phi2",360,-6.3,6.3);
TH1F *pt2 = new TH1F("pt2","pt2",100,0,300);
TH1F *evNo1 = new TH1F("evNo1","evNo1",1000,-1000,3000000000);
TH1F *evNo2 = new TH1F("evNo2","evNo2",1000,-1000,3000000000);
TH1F *runNo1 = new TH1F("runNo1","runNo1",1000,150000,200000);
TH1F *runNo2 = new TH1F("runNo2","runNo2",1000,150000,200000);


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
	for(Int_t overLoop=0;overLoop<1;overLoop++)
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
		for(Int_t loop = 0; loop < 

				chain5->GetEntries()
			//	10000	
				; loop++)
			//SIZE
		{
			// Progress tracker
			if (loop%10000==0) 
			{ 
				cout << loop << " / " << chain5->GetEntries() << endl;
			}

chain5->GetEntry(loop);
eta1->Fill(etaEle5[0]);
eta1->Fill(etaEle5[1]);
phi1->Fill(phiEle5[0]);
phi1->Fill(phiEle5[1]);
pt1->Fill(PtEle5[0]);
pt1->Fill(PtEle5[1]);
evNo1->Fill(eventNumber5);
runNo1->Fill(runNumber5);

			}
			//end of loop	
		for(Int_t loop = 0; loop < 

				chain7->GetEntries()
			//	10000	
				; loop++)
			//SIZE
		{
			// Progress tracker
			if (loop%10000==0) 
			{ 
				cout << loop << " / " << chain7->GetEntries() << endl;
			}

chain7->GetEntry(loop);


eta2->Fill(etaEle7[0]);
eta2->Fill(etaEle7[1]);
phi2->Fill(phiEle7[0]);
phi2->Fill(phiEle7[1]);
pt2->Fill(PtEle7[0]);
pt2->Fill(PtEle7[1]);
evNo2->Fill(eventNumber7);
runNo2->Fill(runNumber7);


			}
		}

		//end of overLoop
	
	//Draw basic Histos
TCanvas *t1 = new TCanvas("t1","",600,600);
eta1->SetLineColor(kRed);
eta2->Draw();
eta1->Draw("same");

TCanvas *t2 = new TCanvas("t2","",600,600);
phi1->SetLineColor(kRed);
phi2->Draw();
phi1->Draw("same");

TCanvas *t3 = new TCanvas("t3","",600,600);
pt1->SetLineColor(kRed);
pt2->Draw();
pt1->Draw("same");


TCanvas *t4 = new TCanvas("t4","",600,600);
evNo1->SetLineColor(kRed);
evNo2->Draw();
evNo1->Draw("same");

TCanvas *t5 = new TCanvas("t5","",600,600);
runNo1->SetLineColor(kRed);
runNo2->Draw();
runNo1->Draw("same");
	//end function
	return 0;
}
