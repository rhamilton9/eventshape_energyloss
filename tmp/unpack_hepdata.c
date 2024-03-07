// A root macro for gathering plots from HEPdata root files
// Made by R. Hamilton on 2/16/2024

void unpack_hepdata(const char* filename) {
  TFile *fin = new TFile(filename, "read");
  if (!fin) {
    std::cout << "File not found!" << std::endl;
    return;
  }
  
  TFile* outfile = new TFile("out.root", "recreate");
  int iTab = 1; int count;
  TDirectoryFile* dir = static_cast<TDirectoryFile*>(fin->Get(Form("Table %i", iTab)));
  while (dir && iTab <= 100) {
    TList* keyList = dir->GetListOfKeys();
    TListIter iter(keyList);
    count = 0;
    while (TObject *name = iter()) {
      TObject *obj = dir->Get(name->GetName());
      
      outfile->cd();
      obj->Write(Form("%s_%i",name->GetName(), iTab));
      ++count;
    }
    std::cout << "Written " << count << " TObjects at iTab = " << iTab << std::endl;
    dir = static_cast<TDirectoryFile*>(fin->Get(Form("Table %i", ++iTab)));
  }
  
  outfile->Close();
  return;
}
