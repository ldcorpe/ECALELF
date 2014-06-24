
{
	TFile *file_Map = TFile::Open("tmp/Map_d1-test.root");
	TFile *file_70X = TFile::Open("../ALCARAW_RECO/ntuple_numEvent1000.root");
	TFile *file_53X = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-22Jan-v1/190645-193621/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012A-22Jan-v1-190645-193621.root");

	TTree *tree_70X =  (TTree*)file_70X->Get("selected"); 
	TTree *tree_53X =  (TTree*)file_53X->Get("selected"); 
	TTree *tree_Map =  (TTree*)file_Map->Get("Map"); 

	Float_t energySCEle_7[2];
	Float_t energySCEle_5[2];
	Float_t R9Ele_7[2];
	Float_t R9Ele_5[2];
	Float_t etaSCEle_7[2];
	Float_t etaSCEle_5[2];
	Float_t phiSCEle_7[2];
	Float_t phiSCEle_5[2];
	Float_t etaEle_7[2];
	Float_t etaEle_5[2];
	Float_t phiEle_7[2];
	Float_t phiEle_5[2];
	Float_t pModeGsfEle_7[2];
	Float_t pModeGsfEle_5[2];
	Int_t eleIndex[2];
	Int_t entryNumber;



	tree_70X->SetBranchAddress("energySCEle",&energySCEle_7);
	tree_53X->SetBranchAddress("energySCEle",&energySCEle_5);
	tree_70X->SetBranchAddress("R9Ele",&R9Ele_7);
	tree_53X->SetBranchAddress("R9Ele",&R9Ele_5);
	tree_70X->SetBranchAddress("etaSCEle",&etaSCEle_7);
	tree_53X->SetBranchAddress("etaSCEle",&etaSCEle_5);
	tree_70X->SetBranchAddress("phiSCEle",&phiSCEle_7);
	tree_53X->SetBranchAddress("phiSCEle",&phiSCEle_5);
	tree_70X->SetBranchAddress("etaEle",&etaEle_7);
	tree_53X->SetBranchAddress("etaEle",&etaEle_5);
	tree_70X->SetBranchAddress("phiEle",&phiEle_7);
	tree_53X->SetBranchAddress("phiEle",&phiEle_5);
	tree_70X->SetBranchAddress("pModeGsfEle",&pModeGsfEle_7);
	tree_53X->SetBranchAddress("pModeGsfEle",&pModeGsfEle_5);
	tree_Map->SetBranchAddress("eleIndex",&eleIndex);
	tree_Map->SetBranchAddress("entryNumber2",&entryNumber);

	TH1F *E_hist = new TH1F("energySCEle","energySCEle (x-y)/x",1000,-0.5,0.5);
	TH1F *R9_hist = new TH1F("R9Ele","R9Ele (x-y)/x",1000,-0.5,0.5);
	TH1F *etaSC_hist = new TH1F("etaSCEle","etaSCEle (x-y)/x",1000,-0.5,0.5);
	TH1F *phiSC_hist = new TH1F("phiSCEle","phiSCEle (x-y)/x",1000,-0.5,0.5);
	TH1F *eta_hist = new TH1F("etaEle","etaEle (x-y)/x",1000,-0.5,0.5);
	TH1F *phi_hist = new TH1F("phiEle","phiEle (x-y)/x",1000,-0.5,0.5);
	TH1F *pModeGsf_hist = new TH1F("pModeGSf","pmode (x-y)/x",1000,-0.5,0.5);

	for(Int_t loop = 0; loop < tree_70X->GetEntries(); loop++)
	{

		tree_Map->GetEntry(loop);
		tree_70X->GetEntry(loop);

		if (entryNumber==-1) continue;

		tree_53X->GetEntry(entryNumber);

		if (eleIndex[0]==1)
		{
			E_hist->Fill((energySCEle_7[0]-energySCEle_5[0])/energySCEle_5[0]);
			E_hist->Fill((energySCEle_7[1]-energySCEle_5[1])/energySCEle_5[1]);
			R9_hist->Fill((R9Ele_7[0]-R9Ele_5[0])/R9Ele_5[0]);
			R9_hist->Fill((R9Ele_7[1]-R9Ele_5[1])/R9Ele_5[1]);
			etaSC_hist->Fill((etaSCEle_7[0]-etaSCEle_5[0])/etaSCEle_5[0]);
			etaSC_hist->Fill((etaSCEle_7[1]-etaSCEle_5[1])/etaSCEle_5[1]);
			phiSC_hist->Fill((phiSCEle_7[0]-phiSCEle_5[0])/phiSCEle_5[0]);
			phiSC_hist->Fill((phiSCEle_7[1]-phiSCEle_5[1])/phiSCEle_5[1]);
			eta_hist->Fill((etaEle_7[0]-etaEle_5[0])/etaEle_5[0]);
			eta_hist->Fill((etaEle_7[1]-etaEle_5[1])/etaEle_5[1]);
			phi_hist->Fill((phiEle_7[0]-phiEle_5[0])/phiEle_5[0]);
			phi_hist->Fill((phiEle_7[1]-phiEle_5[1])/phiEle_5[1]);
			pModeGsf_hist->Fill((pModeGsfEle_7[0]-pModeGsfEle_5[0])/pModeGsfEle_5[0]);
			pModeGsf_hist->Fill((pModeGsfEle_7[1]-pModeGsfEle_5[1])/pModeGsfEle_5[1]);
		}

		if (eleIndex[0]==2)
		{
			E_hist->Fill((energySCEle_7[1]-energySCEle_5[0])/energySCEle_5[0]);
			E_hist->Fill((energySCEle_7[0]-energySCEle_5[1])/energySCEle_5[1]);
			R9_hist->Fill((R9Ele_7[1]-R9Ele_5[0])/R9Ele_5[0]);
			R9_hist->Fill((R9Ele_7[0]-R9Ele_5[1])/R9Ele_5[1]);
			etaSC_hist->Fill((etaSCEle_7[1]-etaSCEle_5[0])/etaSCEle_5[0]);
			etaSC_hist->Fill((etaSCEle_7[0]-etaSCEle_5[1])/etaSCEle_5[1]);
			phiSC_hist->Fill((phiSCEle_7[1]-phiSCEle_5[0])/phiSCEle_5[0]);
			phiSC_hist->Fill((phiSCEle_7[0]-phiSCEle_5[1])/phiSCEle_5[1]);
			eta_hist->Fill((etaEle_7[1]-etaEle_5[0])/etaEle_5[0]);
			eta_hist->Fill((etaEle_7[0]-etaEle_5[1])/etaEle_5[1]);
			phi_hist->Fill((phiEle_7[1]-phiEle_5[0])/phiEle_5[0]);
			phi_hist->Fill((phiEle_7[0]-phiEle_5[1])/phiEle_5[1]);
			pModeGsf_hist->Fill((pModeGsfEle_7[1]-pModeGsfEle_5[0])/pModeGsfEle_5[0]);
			pModeGsf_hist->Fill((pModeGsfEle_7[0]-pModeGsfEle_5[1])/pModeGsfEle_5[1]);
		}
	}

	TCanvas *t = new TCanvas("t","t",600,600);
	E_hist->Draw();
	t->SaveAs("plots/energySC.pdf");
	etaSC_hist->Draw();
	t->SaveAs("plots/etaSC.pdf");
	phiSC_hist->Draw();
	t->SaveAs("plots/phiSC.pdf");
	eta_hist->Draw();
	t->SaveAs("plots/eta.pdf");
	phi_hist->Draw();
	t->SaveAs("plots/phi.pdf");
	pModeGsf_hist->Draw();
	t->SaveAs("plots/pModeGsf.pdf");
	R9_hist->Draw();
	t->SaveAs("plots/R9.pdf");
}
