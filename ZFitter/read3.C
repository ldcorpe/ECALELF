{

TFile *_file0 = TFile::Open("tmp/Map_d3-louieTest2.root");
	Int_t entryNumber2;
	Int_t runNumber2;
	Int_t eleIndex[2];
	Float_t counter1=0;
	Int_t counter2=0;
	Map->SetBranchAddress("entryNumber2", &entryNumber2);
	Map->SetBranchAddress("eleIndex", &eleIndex);

//	Map->SetBranchAddress("secondRunNumber", &runNumber2);


	
	for( Long64_t heavyLoop =0 ; heavyLoop < Map->GetEntries() ; heavyLoop++)

	{

	counter1++;
Map->GetEntry(heavyLoop);
if(entryNumber2 ==-1)
{counter2++;
//cout << heavyLoop << " | "<< entryNumber2 << " | eleIndex = {" << eleIndex[0] << ","<< eleIndex[1] << "}" << endl;
}
}
cout << "TOTAL ENTRIES: " << counter1 <<std::endl;
cout << "Matched: " << counter1-counter2 << std::endl;
cout<< "Errors: " << counter2<< " - " << 100*(counter2/counter1) << " percent" << endl;


	}
