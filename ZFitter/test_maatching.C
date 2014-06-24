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
#include "src/event_class.cc"
#include "../../DataFormats/Math/interface/deltaPhi.h"


int test_maatching() 
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


	Int_t entryNumber2;
	Int_t eleIndex[2];
	Int_t counter=0;
	Int_t counter_reg=0;
	Int_t counter_inv=0;
	Int_t counter_flop=0;


	TFile f("test.root","recreate");
	f.cd();
	TTree *newtree = new TTree("test","test");
	newtree->Branch("entryNumber2",&entryNumber2,"entryNumber2/I");
	newtree->Branch("eleIndex",&eleIndex,"eleIndex[2]/I");

	for(Int_t overLoop=0;overLoop<1;overLoop++)
	{


		//Define respective TTrees
		TTree *secondChain =  (TTree*)Files_70X[overLoop]->Get("selected");
		TTree *originalChain =  (TTree*)Files_53X[overLoop]->Get("selected"); 
		TTree *tree_Map =  (TTree*)Files_Map[overLoop]->Get("Map"); 

		// Declare variables to be extracted from the trees
		Int_t runNumber1=-1;
		Int_t runNumber2=-1;
		ULong64_t eventNumber1=-1;
		ULong64_t eventNumber2=-1;
		Float_t         phiSCEle1[2];
		Float_t         etaSCEle1[2];
		Float_t         phiSCEle2[2];
		Float_t         etaSCEle2[2];
		Float_t dEtaSwitch[2] = {10000, 10000};
		Float_t dPhiSwitch[2] = {10000, 10000};
		Float_t dPhi[2] = {10000, 10000};
		Float_t dEta[2] = {10000, 10000};
		ULong64_t eventNo2=0;
		Float_t PtEle1[2];
		Float_t PtEle2[2];

		// Declare limit values of delta eta and phi to call a match
		Float_t deltaEtaLim=0.3;
		Float_t deltaPhiLim=0.3;

		// Declare vairables to fill in newtree

		// Set Status of all branches to 0, and then only set those used to 1
		originalChain->SetBranchStatus("*",0);
		secondChain->SetBranchStatus("*",0);

		originalChain->SetBranchStatus("runNumber", 1);
		originalChain->SetBranchStatus("eventNumber", 1);
		originalChain->SetBranchStatus("etaSCEle", 1);
		originalChain->SetBranchStatus("phiSCEle", 1);
		originalChain->SetBranchStatus("PtEle", 1);
		secondChain->SetBranchStatus("runNumber", 1);
		secondChain->SetBranchStatus("eventNumber", 1);
		secondChain->SetBranchStatus("etaSCEle", 1);
		secondChain->SetBranchStatus("phiSCEle", 1);
		secondChain->SetBranchStatus("PtEle", 1);

		// Set Branch Addresses
		originalChain->SetBranchAddress("runNumber", &runNumber1);
		secondChain->SetBranchAddress("runNumber", &runNumber2);
		originalChain->SetBranchAddress("eventNumber", &eventNumber1);
		secondChain->SetBranchAddress("eventNumber", &eventNo2);
		originalChain->SetBranchAddress("etaSCEle", &etaSCEle1);
		secondChain->SetBranchAddress("etaSCEle", &etaSCEle2);
		originalChain->SetBranchAddress("phiSCEle", &phiSCEle1);
		secondChain->SetBranchAddress("phiSCEle", &phiSCEle2);
		originalChain->SetBranchAddress("PtEle", &PtEle1);
		secondChain->SetBranchAddress("PtEle", &PtEle2);

		std::cout << "branch addresses set" << std::endl;

		// Declare branches for new tree

		// Debug Histograms
		TH2F *debug_spatial = new TH2F("spatial","spatial",50,-2.5,2.5,32,-3.2,3.2);
		TH1F *debug_pT = new TH1F("pT","pT",100,-100,100);

		// declare vector used to store positions of possible matches in second tree
		std::multimap<event_class, Long64_t> myMap;

		// Fill this vector with all possible psoitions
		for (Long64_t j=0 ; j<secondChain->GetEntries() ; j++)
		{
			secondChain->GetEntry(j);
			//	cout << runNumber2 << " " << eventNo2 << endl;
			event_class event(runNumber2,eventNo2);
			myMap.insert(std::pair<event_class,Long64_t>(event,j));
		}

		std::cout << "second tree entries vector done" << std::endl;

		// begin main "Heavy" loop over original chain entries
		for( Long64_t heavyLoop =0 ; heavyLoop < 

				originalChain->GetEntries()

				; heavyLoop++)
		{
			originalChain->GetEntry(heavyLoop);
	//		cout << eventNumber1 << endl;
			// Intialise variables for each step in the heavy loop
			entryNumber2=-1;// will remain -1 unless match is found
			eleIndex[0]=-1;
			eleIndex[1]=-1;
			// progress tracker spits out to screen (useful for large number of events)
			if(heavyLoop % 1 ==0)
			{
				std::cout << std::endl << heavyLoop <<" / "<< originalChain->GetEntries() << std::endl;
				std::cout << heavyLoop <<" / "<< myMap.size()<< std::endl;
			}

			// Begin subloop over second tree events. Loops over remaining available positions in second tree (using iterator)

			event_class event1(runNumber1,eventNumber1);
			std::pair <std::multimap<event_class,Long64_t>::iterator, std::multimap<event_class,Long64_t>::iterator> range;
			range = myMap.equal_range(event1);

			//		std::cout << heavyLoop<< " - run: " << runNumber1 << " event: "<< eventNumber1 << std::endl; 
			//		std::cout << "heavyLoop: " <<  myMap.count(event1) <<" potential matches" <<std::endl;

			//for (std::multimap<event_class, Long64_t>::iterator subLoop=range.first; subLoop !=range.second ; subLoop++)
			for (Int_t subLoop=0; subLoop <secondChain->GetEntries() ; subLoop++)
			{
				secondChain->GetEntry(subLoop);//->second);
				//		cout << PtEle2[0] <<  " " << PtEle2[1]<< endl;
				if(PtEle2[0]==0 || PtEle2[1]==0) { counter++ ;continue;}


				// Calculate dEta and dPhi, but also switching around electron 1 and 2 in each entry.
				dPhi[0]=fabs(deltaPhi(phiSCEle1[0],phiSCEle2[0]));
				dEta[0]=fabs((etaSCEle1[0]-etaSCEle2[0]));
				dEta[1]=fabs((etaSCEle1[1]-etaSCEle2[1]));
				dPhi[1]=fabs(deltaPhi(phiSCEle1[1],phiSCEle2[1]));


				dPhiSwitch[0]=fabs(deltaPhi(phiSCEle1[0],phiSCEle2[1])); // switched versions
				dPhiSwitch[1]=fabs(deltaPhi(phiSCEle1[1],phiSCEle2[0]));
				dEtaSwitch[0]=fabs(etaSCEle1[0]-etaSCEle2[1]);
				dEtaSwitch[1]=fabs(etaSCEle1[1]-etaSCEle2[0]);


				// check if regular or inverted match & set correct electron index.
				if( dEta[0] < deltaEtaLim && dPhi[0] < deltaPhiLim && dEta[1] < deltaEtaLim && dPhi[1] < deltaPhiLim) 
				{
					eleIndex[0]=1;
					eleIndex[1]=2;
					entryNumber2 = subLoop; //->second;
					counter_reg++;
	//				myMap.erase(subLoop);
					break;

					std::cout << heavyLoop << " | " << subLoop
					<<std::endl;
					//->second << std::endl;
				}


				if( dEtaSwitch[0] < deltaEtaLim && dPhiSwitch[0] < deltaPhiLim && dEtaSwitch[1] < deltaEtaLim && dPhiSwitch[1] < deltaPhiLim) 
				{
					eleIndex[0]=2;
					eleIndex[1]=1;

					entryNumber2 = subLoop; //->second;
					counter_inv++;
				//	std::cout << "Inverted "<<heavyLoop << endl; 
	//				myMap.erase(subLoop);
					break;
				}

				// Match has been found. Now record the entry number from second tree and delete relevant entry in vector
				//		std::cout << "MATCH ENTRY : " << heavyLoop << " "<< entryNumber2 << std:: endl;
				// if a match is found, break the loop

				//	entryNumber2 = subLoop->second;
				//	break;
				counter_flop++;
			}

			// Fill the new tree with relevant branch values.
			newtree->Fill();
		}
		cout << "0 pT: " << counter << endl;
		cout << " reg: " << counter_reg << endl;
		cout << " inv: " << counter_inv << endl;
		cout << "flop: " << counter_flop << endl;
	}

	newtree->Write();
	f.Close();
	//end function
	return 0;
}
