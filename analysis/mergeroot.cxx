void mergeroot(){

  TList *list = new TList;

  TFile file1("/home/nisiyama/workspace/geant4/msm/result/run0035.root");
  TTree* tree1=(TTree*)file1.Get("sd_data");
  list->Add(tree1);

  TFile file2("/home/nisiyama/workspace/geant4/msm/result/run0036.root");
  TTree* tree2=(TTree*)file2.Get("sd_data");
  list->Add(tree2);

  TFile file3("/home/nisiyama/workspace/geant4/msm/result/run0037.root");
  TTree* tree3=(TTree*)file3.Get("sd_data");
  list->Add(tree3);

  TFile fout ("/home/nisiyama/workspace/geant4/msm/result/run0035-37.root", "recreate");
  TTree *newtree = TTree::MergeTrees(list);
  newtree->SetName("sd_data");
  newtree->Write();

}
