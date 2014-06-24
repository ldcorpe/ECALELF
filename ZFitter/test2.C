{

TFile *file_70X = TFile::Open("root://eoscms//eos/cms/store/group/alca_ecalcalib/ecalelf/ntuples/8TeV/ALCARECO/DoubleElectron-ZSkim-RUN2012A-15Apr2014-v2/190645-193621/lumi/DoubleElectron-ZSkim-RUN2012A-15Apr2014-v2-190645-193621.root");



TTree *chain7 =  (TTree*)file_70X->Get("selected"); 


float invMass;

TBranch *b_invmass;


chain7->SetBranchAddress("invMass_SC",&invMass);

for (int i=0; i<100 ; i++)

{
chain7->GetEntry(i);

cout <<invMass<< endl;


}

}


