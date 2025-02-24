
#include "grace_utils/glauber_tools.h"
#include "glauber_config.h"

// Mergemode: which trees to merge.
// 	0: merge finished trees from calc_area.c
// 	1: merge event trees (no grids) from runGlauber_local.c
void mergetree(int mergemode = 0) {
  const Int_t ncent = sizeof(centrality_list)/sizeof(Int_t[2]);
  
  char outfilename[50];
  if (mergemode == 1) {
    snprintf(outfilename, 50, "out/glauber_%s%s%.2fTeV_nogrid.root",speciesA,speciesB,sqrt_s);
  } else { //default
    if ((int)(100 * sqrt_s)%10==0) snprintf(outfilename, 50, "compiled/%s%s%.1fTeV_fullglauber.root",speciesA,speciesB,sqrt_s);
    else                           snprintf(outfilename, 50, "compiled/%s%s%.2fTeV_fullglauber.root",speciesA,speciesB,sqrt_s);
  }TFile* outfile = new TFile(outfilename, "recreate");
  
  // Files for local copying
  Int_t       lemonncoll;
  Int_t       lemonnpart;
  Int_t       lemonnparta;
  Int_t       lemonnpartb;
  Float_t     lemonb;

  Double_t    lemonpsi_arit;
  Double_t    lemonpsi_geom;
  Double_t    lemonareaA;
  Double_t    lemonareaO;
  Double_t    lemonareaW_arit;
  Double_t    lemonareaW_geom;
  
  Float_t     lemoneccgaus[10];
  Float_t     lemonpartpsi[10];
  Int_t       lemonnproj;
  Float_t     lemonxproj[250];
  Float_t     lemonyproj[250];
  bool        lemonwproj[250];
  Int_t       lemonntarg;
  Float_t     lemonxtarg[250];
  Float_t     lemonytarg[250];
  bool        lemonwtarg[250];
  
  TTree* outtree;
  if (mergemode == 0) {
    outtree = new TTree("glauberTree","glauberTree");
    outtree->Branch("npartA",	  &lemonnparta,		"npartA/I");
    outtree->Branch("npartB",	  &lemonnpartb,		"npartB/I");
    outtree->Branch("ncoll",	  &lemonncoll,		"ncoll/I");
    outtree->Branch("b",	  &lemonb,		"b/F");
    outtree->Branch("psi2_arit",  &lemonpsi_arit,	"psi2_arit/D");
    outtree->Branch("psi2_geom",  &lemonpsi_geom,	"psi2_geom/D");
    outtree->Branch("areaA",	  &lemonareaA,		"areaA/D");
    outtree->Branch("areaO",	  &lemonareaO,		"areaO/D");
    outtree->Branch("areaW_arit", &lemonareaW_arit,	"areaW_arit/D");
    outtree->Branch("areaW_geom", &lemonareaW_geom,	"areaW_geom/D");
  } else if (mergemode == 1) {
    outtree = new TTree("lemon","lemon");
    outtree->Branch("npart",	&lemonnpart,	"npart/I");
    outtree->Branch("nparta",	&lemonnparta,	"nparta/I");
    outtree->Branch("npartb",	&lemonnpartb,	"npartb/I");
    outtree->Branch("ncoll",	&lemonncoll,	"ncoll/I");
    outtree->Branch("b",	&lemonb,	"b/F");
    outtree->Branch("eccgaus",	lemoneccgaus,	"eccgaus[10]/F");
    outtree->Branch("partpsi",	lemonpartpsi,	"partpsi[10]/F");
    outtree->Branch("nproj",	&lemonnproj,	"nproj/I");
    outtree->Branch("xproj",	lemonxproj,	"xproj[250]/F");
    outtree->Branch("yproj",	lemonyproj,	"yproj[250]/F");
    outtree->Branch("wproj",	lemonwproj,	"wproj[250]/O");
    outtree->Branch("ntarg",	&lemonntarg,	"ntarg/I");
    outtree->Branch("xtarg",	lemonxtarg,	"xtarg[250]/F");
    outtree->Branch("ytarg",	lemonytarg,	"ytarg[250]/F");
    outtree->Branch("wtarg",	lemonwtarg,	"wtarg[250]/O");
  }
  // Hists for filling
  TH2D* compED_arit[ncent];
  TH2D* compED_geom[ncent];
  Float_t bmax[ncent];
  Int_t tally_cent[ncent];
  if (mergemode == 0) for (int iCent = 0; iCent < ncent; ++iCent) {
    // Construct averaged histograms if merging output from calc_area.c
    tally_cent[iCent] = 0;
    bmax[iCent] = getImpactParameterFromCentralityClass(centrality_list[iCent][1], sqrt_s, speciesA, speciesB);
    compED_arit[iCent] = new TH2D(Form("avgEA_cent%i-%i_total",centrality_list[iCent][0],centrality_list[iCent][1]), ";x;y;E [A.U.]",
                           ogrid_resolution, -grid_edge, grid_edge,
                           ogrid_resolution, -grid_edge, grid_edge);
    compED_geom[iCent] = new TH2D(Form("avgEG_cent%i-%i_total",centrality_list[iCent][0],centrality_list[iCent][1]), ";x;y;E [A.U.]",
                           ogrid_resolution, -grid_edge, grid_edge,
                           ogrid_resolution, -grid_edge, grid_edge);
  }

  for (Int_t iThread = 1; iThread <= 100; ++iThread) {
    char infilename[100];
    char treename[20];
    TFile* tree_merge_file;
    if (mergemode == 1) {
      snprintf(infilename, 100, "out/glauber_%s%s%.2fTeV_nogrid_thread%i.root",speciesA, speciesB, sqrt_s, iThread);
      snprintf(treename, 20, "lemon");
      tree_merge_file = new TFile(infilename);

      // compile this tree 
      if (!tree_merge_file->Get(treename)) return;
      TTreeReader* reader = new TTreeReader(treename, tree_merge_file);
      TTreeReaderValue<Int_t> npart(*reader, "npart");
      TTreeReaderValue<Int_t> nparta(*reader, "nparta");
      TTreeReaderValue<Int_t> npartb(*reader, "npartb");
      TTreeReaderValue<Int_t> ncoll(*reader, "ncoll");
      TTreeReaderValue<Float_t> b(*reader, "b");
      TTreeReaderArray<Float_t> eccgaus(*reader, "eccgaus");
      TTreeReaderArray<Float_t> partpsi(*reader, "partpsi");
      TTreeReaderValue<Int_t> nproj(*reader, "nproj");
      TTreeReaderArray<Float_t> xproj(*reader, "xproj");
      TTreeReaderArray<Float_t> yproj(*reader, "yproj");
      TTreeReaderArray<bool> wproj(*reader, "wproj");
      TTreeReaderValue<Int_t> ntarg(*reader, "ntarg");
      TTreeReaderArray<Float_t> xtarg(*reader, "xtarg");
      TTreeReaderArray<Float_t> ytarg(*reader, "ytarg");
      TTreeReaderArray<bool> wtarg(*reader, "wtarg");

      // Move contents of reader to output tree
      while (reader->Next()) {
        lemonncoll = *ncoll;
        lemonnpart = *npart;
        lemonnparta = *nparta;
        lemonnpartb = *npartb;
        lemonnproj = *nproj;
        lemonntarg = *ntarg;
        lemonb = *b;
        for (int iEcc = 0; iEcc < 10; ++iEcc) {
	  lemoneccgaus[iEcc] = eccgaus[iEcc];
	  lemonpartpsi[iEcc] = partpsi[iEcc];
        } for (int iNucleon = 0; iNucleon < 250; ++iNucleon) {
	  lemonxproj[iNucleon] = xproj[iNucleon];
	  lemonyproj[iNucleon] = yproj[iNucleon];
	  lemonwproj[iNucleon] = wproj[iNucleon];
	  lemonxtarg[iNucleon] = xproj[iNucleon];
	  lemonytarg[iNucleon] = ytarg[iNucleon];
	  lemonwtarg[iNucleon] = wtarg[iNucleon];
        }
        
      	outtree->Fill();
      }
      std::cout << "Finished writing events from thread " << iThread << std::endl;
      tree_merge_file->Close();
      continue;
    }

    // Only arrives here if working with calc_area.c output 
    snprintf(infilename, 100, "out/areadat_%s%s%.2f_thread%i.root",speciesA, speciesB, sqrt_s, iThread);
    snprintf(treename, 20, "areadat_thread%i",iThread);
    tree_merge_file = new TFile(infilename);
    
    // Set up relevant variables to read thread output trees
    if (!tree_merge_file->Get(treename)) return;
    TTreeReader* reader = new TTreeReader(treename, tree_merge_file);
    TTreeReaderValue<Int_t> nparta(*reader, "npartA");
    TTreeReaderValue<Int_t> npartb(*reader, "npartB");
    TTreeReaderValue<Int_t> ncoll(*reader, "ncoll");
    TTreeReaderValue<Float_t> b(*reader, "b");
    TTreeReaderValue<Double_t> psi2_arit(*reader, "psi2_arit");
    TTreeReaderValue<Double_t> psi2_geom(*reader, "psi2_geom");
    TTreeReaderValue<Double_t> areaA(*reader, "areaA");
    TTreeReaderValue<Double_t> areaO(*reader, "areaO");
    TTreeReaderValue<Double_t> areaW_arit(*reader, "areaW_arit");
    TTreeReaderValue<Double_t> areaW_geom(*reader, "areaW_geom");
    
    Int_t count = 0;
    Int_t tally_cent_thistree[ncent];
    for (Int_t iCent = 0; iCent < ncent; ++iCent) tally_cent_thistree[iCent] = 0;

    while (reader->Next()) {
      lemonncoll = *ncoll;
      lemonnparta = *nparta;
      lemonnpartb = *npartb;
      lemonb = *b;
      lemonpsi_arit = *psi2_arit;
      lemonpsi_geom = *psi2_geom;
      lemonareaA = *areaA;
      lemonareaO = *areaO;
      lemonareaW_arit = *areaW_arit;
      lemonareaW_geom = *areaW_geom; 
      outtree->Fill();
      
      // Determine centrality from impact parameter
      Int_t iCent = 0;
      for (Int_t testCent = 0; testCent < ncent; ++testCent) {
        if (*b < bmax[testCent]) break;
        ++iCent;
      }++tally_cent_thistree[iCent];

      ++count;
    }
    
    // Merge histograms
    for (Int_t iCent = 0; iCent < ncent; ++iCent) {
      tally_cent[iCent] += tally_cent_thistree[iCent];
      compED_arit[iCent]->Add(tree_merge_file->Get<TH2D>(Form("avgEA_cent%i-%i_thread%i",
				      centrality_list[iCent][0],centrality_list[iCent][1],iThread)), tally_cent_thistree[iCent]);
      compED_geom[iCent]->Add(tree_merge_file->Get<TH2D>(Form("avgEG_cent%i-%i_thread%i",
				      centrality_list[iCent][0],centrality_list[iCent][1],iThread)), tally_cent_thistree[iCent]);
    } 

    std::cout << "Finished writing events from thread " << iThread << std::endl;
    tree_merge_file->Close();
  }
  
  std::cout << "Finished writing all events." << std::endl;
  outfile->cd();
  outtree->Write(outtree->GetName(), TObject::kOverwrite);
  if (mergemode == 0) for (int iCent = 0; iCent < ncent; ++iCent) {
    compED_arit[iCent]->Scale(1./tally_cent[iCent]);
    compED_geom[iCent]->Scale(1./tally_cent[iCent]);

    compED_arit[iCent]->Write();
    compED_geom[iCent]->Write();
  }

  // Also write the calculated centrality tree after generating
  if (mergemode == 0) {
    TFile* centrality_file = new TFile(Form("out/glauber_%s%s%.2fTeV_nogrid.root",speciesA, speciesB, sqrt_s));
    TTree* cent_tree = centrality_file->Get<TTree>("cent_tree");
    
    outfile->cd();
    cent_tree->CloneTree()->Write();
  }outfile->Close();
  return;
}


