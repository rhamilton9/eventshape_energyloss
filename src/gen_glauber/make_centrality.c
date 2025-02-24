
#include "grace_utils/glauber_tools.h"
#include "glauber_config.h"
#include "grace_utils/color.h"

Color::TextModifier inv_text(Color::TEXT_BLACK);
Color::TextModifier inv_back(Color::BACK_WHITE);
Color::TextModifier def_text(Color::TEXT_DEFAULT);
Color::TextModifier def_back(Color::BACK_DEFAULT);

void make_centrality() {
  char system[15];
  snprintf(system, 15, "%s%s", speciesA, speciesB);
  
  char filename[100];
  snprintf(filename, 100, "out/glauber_%s%.2fTeV_nogrid.root",system,sqrt_s);
  TFile* infile = new TFile(filename, "update");

  Int_t centrality;
  Float_t b_at_cent;
  TTree* cent_tree = new TTree("cent_tree", "cent_tree");
  cent_tree->Branch("cent",  &centrality, "cent/I");
  cent_tree->Branch("bcent", &b_at_cent,  "b/F");

  centrality = 0;
  b_at_cent = 0.0;
  cent_tree->Fill();
  if (use_tabulated_centrality) {
    for (Int_t iCent = 5; iCent <= 100; iCent+=5) {
      centrality = iCent;
      // Reference tabulated centralities in glauber_tools.h
      b_at_cent = getImpactParameterFromCentralityClass(iCent, sqrt_s, speciesA, speciesB);
      cent_tree->Fill();

      std::cout << Form("Centrality %i%%,  \tb = %.5f", iCent, b_at_cent) << std::endl;
    }
  } else {
    // Construct using TTreeIndex (a modular tree sorter)
    TTree* glauber_tree = infile->Get<TTree>("lemon");
    const int nevent = glauber_tree->GetBranch("b")->GetEntries();
    
    TTreeIndex* index = new TTreeIndex(glauber_tree, "0", "10000000*b");
    Long64_t* index_list = index->GetIndex();
    TBranch* bbranch = glauber_tree->GetBranch("b");
    Float_t b;
    glauber_tree->SetBranchAddress("b", &b);
    
    for (int iCent = 1; iCent < 100; ++iCent) {
      glauber_tree->GetEntry(index_list[iCent * (nevent / 100) - 1]);
      centrality = iCent;
      b_at_cent = b;
      cent_tree->Fill();
      if (iCent % 5 == 0) std::cout << inv_back << inv_text;
      std::cout << Form("Centrality %i%%, index %lli", iCent, index_list[iCent * (nevent / 100) - 1]) << "\t";
      std::cout << Form("Impact parameter b = %.5f", b); 
      std::cout << def_text << def_back << std::endl;
    }
    centrality = 100;
    b_at_cent = 20;
    cent_tree->Fill();
  }
  cent_tree->Write("cent_tree", TObject::kOverwrite);
  return;
}


