	{
	TFile *file_Map = TFile::Open("tmp/Map_d1-louieTest.root");
	TFile *file_70X = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-15Apr2014-v2/190645-193621/lumi/DoubleElectron-ZSkim-RUN2012A-15Apr2014-v2-190645-193621.root");
	TFile *file_53X = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-22Jan-v1/190645-193621/190456-208686-22Jan_v1/GainSwitch_v3/DoubleElectron-ZSkim-RUN2012A-22Jan-v1-190645-193621.root");

	TTree *chain7 =  (TTree*)file_70X->Get("selected"); 
	TTree *chain5 =  (TTree*)file_53X->Get("selected"); 
	TTree *tree_Map =  (TTree*)file_Map->Get("Map"); 
    gStyle->SetOptStat(0);
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
	Float_t energy5x5Ele7[2];	
	Float_t etaEle7[2];	
	Float_t etaSCEle7[2];
	Float_t phiSCEle7[2];
	Float_t phiEle7[2];	
	Float_t invMass7;	
	Float_t energySCEle7[2];
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
	Float_t energy5x5Ele5[2];	
	Float_t etaEle5[2];	
	Float_t etaSCEle5[2];
	Float_t phiEle5[2];	
	Float_t phiSCEle5[2];	
	Float_t invMass5;	
	Float_t energySCEle5[2];
	Float_t pModeGsfEle5[2];
	Float_t pAtVtxGsfEle5[2];
	Float_t PtEle5[2];	
	Float_t seedEnergySCEle5[2];	
	Int_t chargeEle5[2];	
	Float_t R9Ele5[2];	
	Float_t seedXSCEle5[2];
	Float_t seedYSCEle5[2];
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
	TBranch *b_energy5x5Ele7;	
	TBranch *b_etaEle7;	
	TBranch *b_etaSCEle7;
	TBranch *b_phiSCEle7;
	TBranch *b_phiEle7;	
	TBranch *b_invMass7;	
	TBranch *b_energySCEle7;
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
//For 53X//
	TBranch *b_eventNumber5;	
	TBranch *b_weight5;	
	TBranch *b_energy5x5Ele5;	
	TBranch *b_etaEle5;	
	TBranch *b_etaSCEle5;
	TBranch *b_phiSCEle5;
	TBranch *b_phiEle5;	
	TBranch *b_invMass5;	
	TBranch *b_energySCEle5;
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
	chain5->SetBranchAddress("etaSCEle", etaSCEle5,&b_etaSCEle5);
	chain5->SetBranchAddress("phiSCEle", phiSCEle5,&b_phiSCEle5);
	chain5->SetBranchAddress("phiEle", phiEle5,&b_phiEle5);
	chain5->SetBranchAddress("pModeGsfEle",pModeGsfEle5,&b_pModeGsfEle5);
	chain5->SetBranchAddress("pAtVtxGsfEle",pAtVtxGsfEle5,&b_pAtVtxGsfEle5);
//	chain5->SetBranchAddress(invMass_var5.c_str(),&invMass5,&b_invMass5);	
	chain5->SetBranchAddress("energySCEle",energySCEle5,&b_energySCEle5);
	chain5->SetBranchAddress("e5x5SCEle",energy5x5Ele5,&b_energy5x5Ele5);
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
	chain7->SetBranchAddress("pAtVtxGsfEle",pAtVtxGsfEle7,&b_pAtVtxGsfEle7);
//	chain7->SetBranchAddress(invMass_var7.c_str(),&invMass7,&b_invMass7);	
	chain7->SetBranchAddress("energySCEle",energySCEle7,&b_energySCEle7);
	chain7->SetBranchAddress("e5x5SCEle",energy5x5Ele7,&b_energy5x5Ele7);
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



/*
	TH1F *E_hist = new TH1F("energySCEle","energySCEle (_7-_5)/_7",10000,-0.4,0.4);
	TH1F *R9_hist = new TH1F("R9Ele","R9Ele (_7 - _5)",10000,-0.4,0.4);
	TH1F *etaSC_hist = new TH1F("etaSCEle","etaSCEle (_7 - _5)",10000,-0.05,0.05);
	TH1F *phiSC_hist = new TH1F("phiSCEle","phiSCEle (_7 - _5)",10000,-0.05,0.05);
	TH1F *eta_hist = new TH1F("etaEle","etaEle (_7-_5)",10000,-0.1,0.1);
	TH1F *phi_hist = new TH1F("phiEle","phiEle (_7-_5)",10000,-0.1,0.1);
	TH1F *pModeGsf_hist = new TH1F("pModeGSf","pmode (_7-_5)/_7",10000,-0.06,0.06);

	E_hist->Sumw2(); 
	R9_hist->Sumw2(); 
	etaSC_hist->Sumw2(); 
	phiSC_hist->Sumw2(); 
	eta_hist->Sumw2(); 
	phi_hist->Sumw2(); 
	pModeGsf_hist ->Sumw2(); 
 */

	//TGraph2D *energySCMap = new TGraph2D("energyMap","; #eta; #phi; energy ",1000,-0.5,0.5,1000,0,1);

	TProfile2D *r9 = new TProfile2D("mean_deltar9","mean deltar9",360,1,360,170,-85,85,"s");
	TProfile2D *e_sc = new TProfile2D("mean_change_e_sc","mean change e_sc",360,1,360,170,-85,85,"s");


			TH1F *hist_ieta = new TH1F("ieta","ieta",171,-85,85);
			TH1F *hist_iphi = new TH1F("iphi","iphi",361,0,360);


	for(Int_t loop = 0; loop < 
	
	chain7->GetEntries(); loop++)
//	10000; loop++)
	{

		if (loop%10000==0) { 
			cout << loop << " / " << chain7->GetEntries() << endl;

		}


		tree_Map->GetEntry(loop);
		chain7->GetEntry(loop);

		if (entryNumber==-1) continue;

		chain5->GetEntry(entryNumber);

	/*	if (entryNumber%10000==0){
			cout << etaEle7[0] << " | "<<etaEle5[0] <<endl;
			cout << etaEle7[1] << " | "<<etaEle5[1] <<endl;}*/

			selection5=((eleID5[0] & 2)==2)*((eleID5[1] & 2)==2)*(HLTfire5==1)*(recoFlagsEle5[0] > 1)*(recoFlagsEle5[1] > 1)*(PtEle5[0]>20)*(PtEle5[1]>20);
			total_weight5=1;
			total_weight5*=weight5*r9weight5[0]*r9weight5[1]*ptweight5[0]*ptweight5[1];//*mcGenWeight;
			//weight is the puWeight: for MC is not 1
			//mcGenWeight is -1 for data => 1 for MC
			selection5*=total_weight5;


			Float_t iPhi5[2] = { floor(phiSCEle5[0]*57.30+0.5)+181, floor(phiSCEle5[1]*57.30+0.5)+181 };
			Float_t iEta5[2] = { ((etaSCEle5[0]>0)- (etaSCEle5[0]<0))*(floor(fabs(etaSCEle5[0]*57.47)+0.5)+1), ((etaSCEle5[1]>0) - (etaSCEle5[1]<0))*(floor(fabs(etaSCEle5[1]*57.46)+0.5)+1 )};
			
		//	cout <<"phi "<< iPhi5[0] << " " << iPhi5[1] << endl;
	//		cout << "eta "<<iEta5[0] << " " << iEta5[1] << endl;


			if (eleIndex[0]==1)
			{

				if(fabs(etaSCEle5[0])<1.4442)
				{
					r9->Fill(iPhi5[0],iEta5[0],R9Ele7[0]-R9Ele5[0]);
					e_sc->Fill(iPhi5[0],iEta5[0],(energySCEle7[0]-energySCEle5[0])/energySCEle5[0]);
					hist_ieta->Fill(iEta5[0]);
					hist_iphi->Fill(iPhi5[0]);
				}

				if(fabs(etaSCEle5[1])<1.4442)
				{
					r9->Fill(iPhi5[1],iEta5[1],R9Ele7[1]-R9Ele5[1]);
					e_sc->Fill(iPhi5[1],iEta5[1],(energySCEle7[1]-energySCEle5[1])/energySCEle5[1]);
					hist_ieta->Fill(iEta5[1]);
					hist_iphi->Fill(iPhi5[1]);
				}
			}


			if (eleIndex[0]==2)
			{
				if(fabs(etaSCEle5[0])<1.4442)
				{
					r9->Fill(iPhi5[0],iEta5[0],R9Ele7[1]-R9Ele5[0]);
					e_sc->Fill(iPhi5[0],iEta5[0],(energySCEle7[1]-energySCEle5[0])/energySCEle5[0]);
					hist_ieta->Fill(iEta5[0]);
					hist_iphi->Fill(iPhi5[0]);
				}

				if(fabs(etaSCEle5[1])<1.4442)
				{
					r9->Fill(iPhi5[1],iEta5[1],R9Ele7[0]-R9Ele5[1]);
					e_sc->Fill(iPhi5[1],iEta5[1],(energySCEle7[0]-energySCEle5[1])/energySCEle5[1]);
					hist_ieta->Fill(iEta5[1]);
					hist_iphi->Fill(iPhi5[1]);
				}
			}
	}
TCanvas *t = new TCanvas("t","t",600,600);
	t->Divide(1,2);
	t->cd(1);
	r9->SetMinimum(-0.05);
	r9->SetMaximum(0.05);
	r9->Draw("COLZ");
	t->cd(2);
	TH2D *r9_error = r9->ProjectionXY("RMS","C=E");
	r9_error->SetMinimum(-0.03);
	r9_error->SetMaximum(0.03);
	r9_error->Draw("COLZ");
	t->SaveAs("selectionPlots/iEtaPhiR9.pdf");
	
	t->cd(1);
	e_sc->SetMinimum(-0.05);
	e_sc->SetMaximum(0.05);
	e_sc->Draw("COLZ");
	t->cd(2);
	TH2D *e_sc_error = e_sc->ProjectionXY("RMS","C=E");
	e_sc_error->SetMinimum(-0.03);
	e_sc_error->SetMaximum(0.03);
	e_sc_error->Draw("COLZ");
	t->SaveAs("selectionPlots/iEtaPhiEsc.pdf");
	return 0;
/*hist_ieta->Draw();
t->cd(2);
hist_iphi->Draw();*/


	}
