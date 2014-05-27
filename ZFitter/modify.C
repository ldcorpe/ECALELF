{

	TFile *_file0 = TFile::Open("ntuple_numEvent100_2.root");
	_file0->cd();
	Int_t runNumber2;
	ULong64_t eventNumber2;
	Float_t         phiSCEle2[2];
	Float_t         etaSCEle2[2];

	Int_t runNumber3;
	ULong64_t eventNumber3;
	Float_t         phiSCEle3[2];
	Float_t         etaSCEle3[2];

	selected->SetBranchAddress("runNumber", &runNumber2);
	selected->SetBranchAddress("eventNumber", &eventNumber2);
	selected->SetBranchAddress("etaSCEle", &etaSCEle2);
	selected->SetBranchAddress("phiSCEle", &phiSCEle2);


	TTree *newtree = new TTree("selected", "ntuple_modified");
	newtree->Branch("runNumber", &runNumber3,"entryNumber2/I");
	newtree->Branch("eventNumber", &eventNumber3,"eventNumber2/I");
	newtree->Branch("etaSCEle", &etaSCEle3,"etaSCEle[2]/F");
	newtree->Branch("phiSCEle", &phiSCEle3,"phiSCEle[2]/F");


	TRandom *R = new TRandom(0);
	for( Long64_t heavyLoop = 0 ; heavyLoop< 111; heavyLoop	++)
	{

	//cout << R << endl;
		selected->GetEntry(heavyLoop);
		if(R->Uniform() >0.9)
		{
			cout << "Switch electrons " << heavyLoop << endl;

		etaSCEle3[0]= etaSCEle2[1];
		etaSCEle3[1]= etaSCEle2[0];
		phiSCEle3[0]= phiSCEle2[1];
		phiSCEle3[1]= phiSCEle2[0];
		}
else
{
		etaSCEle3[0]= etaSCEle2[0];
		etaSCEle3[1]= etaSCEle2[1];
		phiSCEle3[0]= phiSCEle2[0];
		phiSCEle3[1]= phiSCEle2[1];
}


		if(R->Uniform() >0.9)
		{
			cout << "modify energy " <<heavyLoop << endl;

		etaSCEle3[0]= etaSCEle3[0]*(0.5+R->Uniform());
		etaSCEle3[1]= etaSCEle3[1]*(0.5+R->Uniform());
		phiSCEle3[0]= phiSCEle3[0]*(0.5+R->Uniform());
		phiSCEle3[1]= phiSCEle3[1]*(0.5*R->Uniform());
		}
else
{
		etaSCEle3[0]= etaSCEle3[0];
		etaSCEle3[1]= etaSCEle3[1];
		phiSCEle3[0]= phiSCEle3[0];
		phiSCEle3[1]= phiSCEle3[1];
}

		runNumber3 =runNumber2;
		eventNumber3 = eventNumber2;
		newtree->Fill();
	}

	TFile f ("ntuple_modif.root", "recreate");
	f.cd();
	newtree->Write();
	f.Close();

}
