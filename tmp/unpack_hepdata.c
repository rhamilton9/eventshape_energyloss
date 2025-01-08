// A root macro for gathering plots from HEPdata root files
// Made by R. Hamilton on 2/16/2024

// To use succinctly in root, type in command line in the local directory
//
// root -l -q 'unpack_hepdata("[FILENAME]", {[LIST OF INDICES]})'
//
// e.g., to unpack only tables 3, 7:
//
// root -l -q 'unpack_hepdata("hepdata.root", {3, 7})'
//
// The listings on https://www.hepdata.net/ can be used to distinguish which 
// tables are relevant to the analysis you are doing.

void unpack_hepdata(const char* filename, std::vector<int> toKeep = {}) {
  TFile *fin = new TFile(filename, "read");
  if (!fin) {
    std::cout << "File not found!" << std::endl;
    return;
  }
  
  // Default case: run all tables
  if (toKeep.size() == 0)
    for (int i = 1; i < 100; ++i) toKeep.push_back(i);
  
  // Set up file and key pointers
  TFile* outfile = new TFile("out.root", "recreate");
  int count;
  std::vector<int>::iterator tabIndex = toKeep.begin();
  TDirectoryFile* dir = fin->Get<TDirectoryFile>(Form("Table %i", *tabIndex));
  while (dir && tabIndex != toKeep.end()) {
    // TListIterator to loop over the current directory
    TList* keyList = dir->GetListOfKeys();
    TListIter iter(keyList);
    count = 0;
    while (TObject *named = iter()) {
      TObject *obj = dir->Get(named->GetName());
      
      outfile->cd();
      obj->Write(Form("%s_%i",named->GetName(), *tabIndex));
      ++count;
    }
    std::cout << "Written " << count << " TObjects from table #" << *tabIndex << std::endl;
    dir = fin->Get<TDirectoryFile>(Form("Table %i", *(++tabIndex)));
  }
  
  outfile->Close();
  return;
}
