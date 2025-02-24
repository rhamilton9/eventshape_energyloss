/// Root macro for area computation from Glauber output tree.
/// Produces a number of estimates for the transverse area for each event.
/// Written with parallelization via Grace cluster in mind.
/// Expects input in the format of an IP-Glasma TTree.
/// Such a tree can be generated via runAndOutputLemonTree in the Glauber code.
/// A more straightforward interface with this method is implemented in runGlauber_local.c
///
/// Note that this macro is highly optimized for fast execution in C++
/// Calculation time depends heavily on the grid resolution: O(n^2).
/// For this reason, using as coarse a grid as allowable is highly recommended.
/// Written by R. Hamilton on 2024-08-19

#include "glauber_config.h"
#include "grace_utils/glauber_tools.h"
#include "grace_utils/hist_tools.h"
#include <chrono>

using namespace std;

bool debug = false;

void calc_area(Int_t this_thread = 1, 
               Int_t total_threads = 1,
	       bool already_divided = false) {
  if (this_thread < 1) return;
  if (this_thread > total_threads) this_thread = total_threads;
  
  const Double_t sigNN = getNucleonCrossSection(sqrt_s);
  const Double_t interaction_radius = TMath::Sqrt(sigNN/(TMath::Pi()*10));
  const Double_t interaction_threshold = 2.5*interaction_radius;
  const Double_t d2 = sigNN/(TMath::Pi()*10); // in fm^2
  const Double_t r2 = d2/4.;
  
  const Int_t ncent = sizeof(centrality_list)/sizeof(Int_t[2]);
  char species[15];
  snprintf(species, 15, "%s%s", speciesA, speciesB);

  // Functions for computing energy density
  TF2* radius = new TF2("radius", "sqrt(x*x+y*y)", 0, 10, 0, 10);
  TF2* radius_square = new TF2("radius_square", "x*x+y*y", 0, 10, 0, 10);
  TF2* nucleon_smearing2D = new TF2("nucleon_profile", "[0]*exp([1]*(x*x+y*y))", 0, 10*sigs, 0, 10*sigs);
  nucleon_smearing2D->FixParameter(0, 1./(2*TMath::Pi()*sigs*sigs));
  nucleon_smearing2D->FixParameter(1, -0.5/(sigs*sigs));
  
  // Create output file and tree+histograms in that directory
  TFile* outfile = new TFile(Form("out/areadat_%s%.2f_thread%i.root",species,sqrt_s,this_thread), "recreate");
  TH2D* reco_A  = new TH2D("reco_A", ";x;y;T_{A}",
                           ogrid_resolution, -grid_edge, grid_edge,
                           ogrid_resolution, -grid_edge, grid_edge);
  TH2D* reco_B  = new TH2D("reco_B", ";x;y;T_{B}",
                           ogrid_resolution, -grid_edge, grid_edge,
                           ogrid_resolution, -grid_edge, grid_edge);
  TH2D* reco_E = new TH2D("reco_E_arit", ";x;y;E [A.U.]",
                          ogrid_resolution, -grid_edge, grid_edge,
                          ogrid_resolution, -grid_edge, grid_edge);
  TH2D* reco_G = new TH2D("reco_E_geom", ";x;y;E [A.U.]",
                          ogrid_resolution, -grid_edge, grid_edge,
                          ogrid_resolution, -grid_edge, grid_edge);
  TH2D* avg_E[ncent];
  TH2D* avg_G[ncent];
  for (Int_t iCent = 0; iCent < ncent; ++iCent) {
    avg_E[iCent] = new TH2D(Form("avgEA_cent%i-%i_thread%i",centrality_list[iCent][0],centrality_list[iCent][1],this_thread), ";x;y;E [A.U.]",
                           ogrid_resolution, -grid_edge, grid_edge,
                           ogrid_resolution, -grid_edge, grid_edge);
    avg_G[iCent] = new TH2D(Form("avgEG_cent%i-%i_thread%i",centrality_list[iCent][0],centrality_list[iCent][1],this_thread), ";x;y;E [A.U.]",
                           ogrid_resolution, -grid_edge, grid_edge,
                           ogrid_resolution, -grid_edge, grid_edge);
  }
  
  Float_t b;
  Int_t npartA;
  Int_t npartB;
  Int_t ncoll;
  Double_t areaA;
  Double_t areaO;
  Double_t areaW_arit;
  Double_t areaW_geom;
  Double_t ecc2_arit;
  Double_t ecc2_geom;
  Double_t psi2_arit;
  Double_t psi2_geom;
  TTree* areadat_tree = new TTree(Form("areadat_thread%i",this_thread),Form("areadat_thread%i",this_thread));
  areadat_tree->Branch("b",           &b,           "b/F");
  areadat_tree->Branch("npartA",      &npartA,      "npartA/I");
  areadat_tree->Branch("npartB",      &npartB,      "npartB/I");
  areadat_tree->Branch("ncoll",       &ncoll,       "ncoll/I");
  areadat_tree->Branch("areaA",       &areaA,       "AreaA/D");
  areadat_tree->Branch("areaO",       &areaO,       "AreaO/D");
  areadat_tree->Branch("areaW_arit",  &areaW_arit,  "AreaW_arit/D");
  areadat_tree->Branch("areaW_geom",  &areaW_geom,  "AreaW_geom/D");
  areadat_tree->Branch("psi2_arit",   &psi2_arit,   "psi2_arit/D");
  areadat_tree->Branch("psi2_geom",   &psi2_geom,   "psi2_geom/D");
  
  // Access file with centrality reference
  char centfile_name[100];
  snprintf(centfile_name, 100, "out/glauber_%s%.2fTeV_nogrid.root",species,sqrt_s); 
  TFile* centrality_file = new TFile(centfile_name);
  
  // Use the established centrality from make_centrality
  TTree* centrality_tree = centrality_file->Get<TTree>("cent_tree");
  TTreeIndex* cent_index = new TTreeIndex(centrality_tree, "cent", "0");
  Long64_t* cent_sortedlist = cent_index->GetIndex();
  
  Float_t b_cent;
  centrality_tree->SetBranchAddress("bcent", &b_cent);
  Float_t bmax[ncent];
  Int_t tally_cent[ncent];
  for (int iCent = 0; iCent < ncent; ++iCent) {
    centrality_tree->GetEntry(cent_index->GetEntryNumberWithIndex(centrality_list[iCent][1], 0));
    bmax[iCent] = b_cent;
    tally_cent[iCent] = 0;
    std::cout << Form("cent %i%%, b = %.5f",centrality_list[iCent][1], b_cent) << std::endl;
  }
  
  
  // Access file with glauberdata
  TFile* glauber_file;
  if (already_divided) {
    char glauberfile_name[100];
    snprintf(glauberfile_name, 100, "out/glauber_%s%.2fTeV_nogrid_thread%i.root", species, sqrt_s, this_thread);
    glauber_file = new TFile(glauberfile_name);
  } else glauber_file = centrality_file; 
  
  // Get pointer to Data tree
  TTree* glauber_tree = glauber_file->Get<TTree>("lemon");
  const int nevent = glauber_tree->GetBranch("b")->GetEntries();
  
  // Assign reader to access data
  TTreeReader* reader = new TTreeReader("lemon", glauber_file);
  TTreeReaderValue<Float_t> b_base(*reader, "b");
  TTreeReaderValue<Int_t> nA(*reader, "nproj");
  TTreeReaderValue<Int_t> nB(*reader, "ntarg");
  TTreeReaderValue<Int_t> nparta(*reader, "nparta");
  TTreeReaderValue<Int_t> npartb(*reader, "npartb");
  TTreeReaderValue<Int_t> ncoll_base(*reader, "ncoll");
  TTreeReaderArray<Float_t> xA(*reader, "xproj");
  TTreeReaderArray<Float_t> yA(*reader, "yproj");
  TTreeReaderArray<bool> woundedA(*reader, "wproj");
  TTreeReaderArray<Float_t> xB(*reader, "xtarg");
  TTreeReaderArray<Float_t> yB(*reader, "ytarg");
  TTreeReaderArray<bool> woundedB(*reader, "wtarg");
  TTreeReaderArray<Float_t> partpsi(*reader, "partpsi");
  const Int_t max_nucleons = 250;
  
  Int_t total_ecctime = 0;
  
  Int_t count = 0;
  Double_t bincent[ogrid_resolution];
  for (Int_t iBin = 1; iBin <= ogrid_resolution; ++iBin)
    bincent[iBin-1] = reco_A->GetXaxis()->GetBinCenter(iBin);
    
  auto starttime = std::chrono::high_resolution_clock::now();
  while(reader->Next()) {
    if (!already_divided && count < (this_thread - 1)*nevent/total_threads) {++count; continue;}
    if (!already_divided && count >= (this_thread   )*nevent/total_threads) break;
    if (debug) std::cout << "===================================================Event " << count << std::endl;
    auto localtime = std::chrono::high_resolution_clock::now();
    
    // Assign values to new tree that aren't modified
    b = *b_base;
    npartA = *nparta;
    npartB = *npartb;
    ncoll = *ncoll_base;
    psi2_arit = partpsi[2];
    
    // Determine centrality from impact parameter
    Int_t iCent = 0;
    for (Int_t testCent = 0; testCent < ncent; ++testCent) {
      if (b < bmax[testCent]) break;
      ++iCent;
    }++tally_cent[iCent];
    
    // Check that centrality is correct:
    if (debug) {
      if (iCent > 0 )
        std::cout << "b = " << b << ", cent = " << iCent << Form(" corresponding to b-range %.2f-%.2f", bmax[iCent-1],bmax[iCent]) << std::endl;
      else
        std::cout << "b = " << b << ", cent = " << iCent << Form(" corresponding to b-range 0.00-%.2f",bmax[iCent]) << std::endl;
    }
      
    // Construct arrays needed for calculation in each event
    Double_t nucleonsA_x[max_nucleons];
    Double_t nucleonsA_y[max_nucleons];
    Double_t nucleonsB_x[max_nucleons];
    Double_t nucleonsB_y[max_nucleons];
    bool  woundA[max_nucleons];
    bool  woundB[max_nucleons];
    bool  interA[max_nucleons];
    bool  interB[max_nucleons];
    
    // Calculate CM beforehand to align participants for Arit ED
    std::vector<Double_t> cm_frompoint = {0, 0};
    Int_t count_totalwounded = 0;
    for (Int_t iNucleon = 0; iNucleon < *nA; ++iNucleon) {
      if (!woundedA[iNucleon]) continue;
      cm_frompoint.at(0) += xA[iNucleon];
      cm_frompoint.at(1) += yA[iNucleon];
      ++count_totalwounded;
    } for (Int_t iNucleon = 0; iNucleon < *nB; ++iNucleon) {
      if (!woundedB[iNucleon]) continue;
      cm_frompoint.at(0) += xB[iNucleon];
      cm_frompoint.at(1) += yB[iNucleon];
      ++count_totalwounded;
    }
    cm_frompoint[0] /=(Double_t)count_totalwounded;
    cm_frompoint[1] /=(Double_t)count_totalwounded;
    if (debug) {
      std::cout << Form("Arit CM from nucl: (%.5f, %.5f)", cm_frompoint[0],cm_frompoint[1]) << std::endl;
      std::cout << "Part Plane: " << partpsi[2] << std::endl;
    }
    
    // Orient nucleons to zero CM and Psi2 for Arit ED
    for (Int_t i = 0; i < *nA; ++i) {
      nucleonsA_x[i] =  (xA[i] - cm_frompoint[0])*TMath::Cos(partpsi[2]) + (yA[i] - cm_frompoint[1])*TMath::Sin(partpsi[2]);
      nucleonsA_y[i] = -(xA[i] - cm_frompoint[0])*TMath::Sin(partpsi[2]) + (yA[i] - cm_frompoint[1])*TMath::Cos(partpsi[2]);
      nucleonsB_x[i] =  (xB[i] - cm_frompoint[0])*TMath::Cos(partpsi[2]) + (yB[i] - cm_frompoint[1])*TMath::Sin(partpsi[2]);
      nucleonsB_y[i] = -(xB[i] - cm_frompoint[0])*TMath::Sin(partpsi[2]) + (yB[i] - cm_frompoint[1])*TMath::Cos(partpsi[2]);
      woundA[i] = woundedA[i];
      woundB[i] = woundedB[i];
    }
    
    
    // Assign interacting nucleons to avoid performing calculations for long-range spectators
    for (Int_t i = 0; i < *nA; ++i) {
      interA[i] = false;
      if (woundA[i]) interA[i] = true;
      for (Int_t j = 0; j < *nA; ++j) {
        if (radius->Eval(nucleonsA_x[i] - nucleonsB_x[j], nucleonsA_y[i] - nucleonsB_y[j]) <= interaction_threshold) {
          interA[i] = true;
          break;
        }
      }
      
      interB[i] = false;
      if (woundB[i]) interB[i] = true;
      for (Int_t j = 0; j < *nB; ++j) {
        if (radius->Eval(nucleonsB_x[i] - nucleonsA_x[j], nucleonsB_y[i] - nucleonsA_y[j]) <= interaction_threshold) {
          interB[i] = true;
          break;
        }
      }
    }
    
    auto accesstime = std::chrono::high_resolution_clock::now();
    Int_t accessdur = std::chrono::duration_cast<std::chrono::milliseconds>(accesstime - localtime).count();
    if (debug) std::cout << Form("Accessing data from tree: %02d:%02d.%03d", accessdur/60000, accessdur/1000 - 60*(accessdur/60000), accessdur%1000) << std::endl;
    
    
    Double_t energy_threshold = nucleon_smearing2D->Eval(0, TMath::Sqrt(r2));
    TAxis* axis = reco_A->GetXaxis();
    TH2S* flagA = new TH2S("flag_area_A",";x [fm];y [fm]",
                           ogrid_resolution,-grid_edge,grid_edge,
                           ogrid_resolution,-grid_edge,grid_edge);
    TH2S* flagB = new TH2S("flag_area_B",";x [fm];y [fm]",
                           ogrid_resolution,-grid_edge,grid_edge,
                           ogrid_resolution,-grid_edge,grid_edge);
    for (Int_t iNucleon = 0; iNucleon < *nA; ++iNucleon) {  //-----------------Primary loop (A): Building grid
      // Skip any non-interacting nucleons
      if (!interA[iNucleon]) continue;
      Int_t binrange_x[2] = {
        TMath::Max(axis->FindBin(nucleonsA_x[iNucleon] - interaction_threshold), 1),
        TMath::Min(axis->FindBin(nucleonsA_x[iNucleon] + interaction_threshold), ogrid_resolution)
      };
      Int_t binrange_y[2] = {
        TMath::Max(axis->FindBin(nucleonsA_y[iNucleon] - interaction_threshold), 1),
        TMath::Min(axis->FindBin(nucleonsA_y[iNucleon] + interaction_threshold), ogrid_resolution)
      };
      
      if (!woundA[iNucleon]) {
        // Only geometric edensity is affected by non-participant interacting nucleons
        for (Int_t iX = binrange_x[0]; iX <= binrange_x[1]; ++iX) {
          for (Int_t iY = binrange_y[0]; iY <= binrange_y[1]; ++iY) {
            
            // Nuclear Thickness only for geometric calculation
            reco_A->SetBinContent(iX, iY, reco_A->GetBinContent(iX, iY)
                                  + nucleon_smearing2D->Eval(bincent[iX-1] - nucleonsA_x[iNucleon], bincent[iY-1] - nucleonsA_y[iNucleon]));
          }
        }// End of grid loop
      } else {
        // Arithmetic calculations use only participant nucleons
        for (Int_t iX = binrange_x[0]; iX <= binrange_x[1]; ++iX) {
          for (Int_t iY = binrange_y[0]; iY <= binrange_y[1]; ++iY) {
            Double_t hold = nucleon_smearing2D->Eval(bincent[iX-1] - nucleonsA_x[iNucleon], bincent[iY-1] - nucleonsA_y[iNucleon]);
            
            // Nuclear Thickness, arithmetic energy density and area
            reco_A->SetBinContent(iX, iY, reco_A->GetBinContent(iX, iY) + hold);
            reco_E->SetBinContent(iX, iY, reco_E->GetBinContent(iX, iY) + hold);
            
            if (hold >= energy_threshold) flagA->SetBinContent(iX, iY, 1);
          }
        }// End of grid loop
      }
    } for (Int_t iNucleon = 0; iNucleon < *nB; ++iNucleon) { //-----------------Primary loop (B): Building grid
      // Skip any non-interacting nucleons
      if (!interB[iNucleon]) continue;
      Int_t binrange_x[2] = {
        TMath::Max(axis->FindBin(nucleonsB_x[iNucleon] - interaction_threshold), 1),
        TMath::Min(axis->FindBin(nucleonsB_x[iNucleon] + interaction_threshold), ogrid_resolution)
      };
      Int_t binrange_y[2] = {
        TMath::Max(axis->FindBin(nucleonsB_y[iNucleon] - interaction_threshold), 1),
        TMath::Min(axis->FindBin(nucleonsB_y[iNucleon] + interaction_threshold), ogrid_resolution)
      };
      
      if (!woundB[iNucleon]) {
        // Only geometric edensity is affected by non-participant interacting nucleons
        for (Int_t iX = binrange_x[0]; iX <= binrange_x[1]; ++iX) {
          for (Int_t iY = binrange_y[0]; iY <= binrange_y[1]; ++iY) {
            // Nuclear Thickness only for geometric calculation
            reco_B->SetBinContent(iX, iY, reco_B->GetBinContent(iX, iY)
                                  + nucleon_smearing2D->Eval(bincent[iX-1] - nucleonsB_x[iNucleon], bincent[iY-1] - nucleonsB_y[iNucleon]));
          }
        }// End of grid loop
      } else {
        // Arithmetic calculations use only participant nucleons
        for (Int_t iX = binrange_x[0]; iX <= binrange_x[1]; ++iX) {
          for (Int_t iY = binrange_y[0]; iY <= binrange_y[1]; ++iY) {
            Double_t hold = nucleon_smearing2D->Eval(bincent[iX-1] - nucleonsB_x[iNucleon], bincent[iY-1] - nucleonsB_y[iNucleon]);
            
            // Nuclear Thickness, arithmetic energy density and area
            reco_B->SetBinContent(iX, iY, reco_B->GetBinContent(iX, iY) + hold);
            reco_E->SetBinContent(iX, iY, reco_E->GetBinContent(iX, iY) + hold);
            
            if (hold >= energy_threshold) flagB->SetBinContent(iX, iY, 1);
          }
        }// End of grid loop
      }
    }// Finished with grid construction
    
    // Tally areas from grid threshold flag hists
    areaA = 0;
    areaO = 0;
    for (Int_t iX = 1; iX <= ogrid_resolution; ++iX)
      for (Int_t iY = 1; iY <= ogrid_resolution; ++iY) {
        if (flagA->GetBinContent(iX, iY) && flagB->GetBinContent(iX, iY)) ++areaA;
        if (flagA->GetBinContent(iX, iY) || flagB->GetBinContent(iX, iY)) ++areaO;
      }
    Double_t da_ogrid = reco_A->GetXaxis()->GetBinWidth(1)*reco_A->GetYaxis()->GetBinWidth(1);
    areaA *= da_ogrid;
    areaO *= da_ogrid;
    
    if (debug) {
      std::cout << "AreaO: " << areaO << std::endl;
      std::cout << "AreaA: " << areaA << std::endl;
    }
    
    auto looptime = std::chrono::high_resolution_clock::now();
    Int_t loopdur = std::chrono::duration_cast<std::chrono::milliseconds>(looptime - accesstime).count();
    if (debug) std::cout << Form("Time for calculation of A,B,aritED: %02d:%02d.%03d", loopdur/60000, loopdur/1000 - 60*(loopdur/60000), loopdur%1000) << std::endl;
    
    
    // Compute geometric energy density and epsilon/psi, re-align
    reco_G = translateHist_simple(geometricEnergyDensity(reco_A, reco_B), 0, 0, true);
    Double_t avgrad[2] = {0,0};
    Double_t avgsin[2] = {0,0};
    Double_t avgcos[2] = {0,0};
    Double_t sum_w1[2] = {0,0};
    Double_t avg_x1[2] = {0,0};
    Double_t avg_y1[2] = {0,0};
    Double_t avg_x2[2] = {0,0};
    Double_t avg_y2[2] = {0,0};
    Double_t avg_xy[2] = {0,0};
    for (Int_t iX = 1; iX <= ogrid_resolution; ++iX)
      for (Int_t iY = 1; iY <= ogrid_resolution; ++iY) {
        // Calculate averages for epsilon2, psi2
	Double_t cr = radius_square->Eval(bincent[iX], bincent[iY]);
        Double_t ct = 2*TMath::ATan2(bincent[iY], bincent[iX]);
        Double_t ce = reco_E->GetBinContent(iX, iY);
        Double_t cg = reco_G->GetBinContent(iX, iY);
        avgrad[0] += ce*cr;
        avgsin[0] += ce*cr*TMath::Sin(ct);
        avgcos[0] += ce*cr*TMath::Cos(ct);
        avgrad[1] += cg*cr;
        avgsin[1] += cg*cr*TMath::Sin(ct);
        avgcos[1] += cg*cr*TMath::Cos(ct);

	// Calculate means for variance areas
	sum_w1[0] += ce;
	avg_x1[0] += ce*bincent[iX];
	avg_y1[0] += ce*bincent[iY];
	avg_x2[0] += ce*bincent[iX]*bincent[iX];
	avg_y2[0] += ce*bincent[iY]*bincent[iY];
	avg_xy[0] += ce*bincent[iX]*bincent[iY];
	sum_w1[1] += cg;
	avg_x1[1] += cg*bincent[iX];
	avg_y1[1] += cg*bincent[iY];
	avg_x2[1] += cg*bincent[iX]*bincent[iX];
	avg_y2[1] += cg*bincent[iY]*bincent[iY];
	avg_xy[1] += cg*bincent[iX]*bincent[iY];
      }
    ecc2_arit = radius->Eval(avgsin[0], avgcos[0])/avgrad[0];
    ecc2_geom = radius->Eval(avgsin[1], avgcos[1])/avgrad[1];
    psi2_geom = 0.5*(TMath::ATan2(avgsin[1], avgcos[1]) + TMath::Pi());
    
    // Calculate variance areas (areaW)
    Double_t sigmaX2[2];
    Double_t sigmaY2[2];
    Double_t sigmaXY[2];
    for (Int_t iScale = 0; iScale < 2; ++iScale) {
      avg_x1[iScale] /= sum_w1[iScale];
      avg_y1[iScale] /= sum_w1[iScale];
      avg_x2[iScale] /= sum_w1[iScale];
      avg_y2[iScale] /= sum_w1[iScale];
      avg_xy[iScale] /= sum_w1[iScale];
      
      sigmaX2[iScale] = avg_x2[iScale] - TMath::Sq(avg_x1[iScale]);
      sigmaY2[iScale] = avg_y2[iScale] - TMath::Sq(avg_y1[iScale]);
      sigmaXY[iScale] = avg_xy[iScale] - avg_x1[iScale]*avg_y1[iScale];
    }
    areaW_arit = TMath::Sqrt(sigmaX2[0]*sigmaY2[0] - sigmaXY[0]*sigmaXY[0]); 
    areaW_geom = TMath::Sqrt(sigmaX2[1]*sigmaY2[1] - sigmaXY[1]*sigmaXY[1]); 
    
    reco_G = rotateHist2D_simple(reco_G, -psi2_geom, true);
    psi2_geom += psi2_arit; // restore the original rotation performed to nucleon dist to recover psi2 relative to acual reaction plane.
//    psi2_arit = getParticipantPlaneAngle(reco_E);
//    Double_t geompsi_reco = getParticipantPlaneAngle(reco_G);
//    std::cout << Form("Total geom psi2: %.5f", geompsi_reco) << std::endl;
    if (debug) std::cout << Form("Ecc: a=%.5f, g=%.5f", ecc2_arit, ecc2_geom) << std::endl;
    
    auto geomtime = std::chrono::high_resolution_clock::now();
    Int_t geomdur = std::chrono::duration_cast<std::chrono::milliseconds>(geomtime - looptime).count();
    if (debug) std::cout << Form("Time for calculation of eccentricities: %02d:%02d.%03d", geomdur/60000, geomdur/1000 - 60*(geomdur/60000), geomdur%1000) << std::endl;
    total_ecctime += geomdur;
    
    
    
    // Add to event averaged energy densities
    avg_E[iCent]->Add(avg_E[iCent], reco_E);
    avg_G[iCent]->Add(avg_G[iCent], reco_G);
    
    // Clear hists and prepare for next event
    delete flagA;
    delete flagB;
    reco_A->Reset();
    reco_B->Reset();
    reco_E->Reset();
    reco_G->Reset();
    
    Int_t ctime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - starttime).count();
    if (debug) std::cout << Form("Current runtime: %02d:%02d.%03d", ctime/60000, ctime/1000 - 60*(ctime/60000), ctime%1000) << std::endl;
    ++count;
    if (count % 100 == 0) std::cout << "Finished with event " << count << "." << std::endl;
    areadat_tree->Fill();
  }// End of glauber base tree access loop
  glauber_file->Close();
  
  // Store in output tree
  outfile->cd();
  areadat_tree->Write();
  for (int iCent = 0; iCent < ncent; ++iCent) {
    if (tally_cent[iCent] > 0) {
      avg_E[iCent]->Scale(1./tally_cent[iCent]);
      avg_G[iCent]->Scale(1./tally_cent[iCent]);
    }
    avg_E[iCent]->Write();
    avg_G[iCent]->Write();
  }outfile->Close();
  
  total_ecctime /= count;
  std::cout << Form("Average eccentricitiy calc time: %02d:%02d.%03d", total_ecctime/60000, total_ecctime/1000 - 60*(total_ecctime/60000), total_ecctime%1000) << std::endl;
  
  auto endtime = std::chrono::high_resolution_clock::now();
  Int_t elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
  std::cout << Form("Total elapsed time: %02d:%02d.%03d", elapsed_ms/60000, elapsed_ms/1000 - (elapsed_ms/60000)*60, elapsed_ms%1000) << std::endl;
  Int_t avgtime = elapsed_ms/count;
  std::cout << Form("Avg. per event: %02d:%02d.%03d", avgtime/60000, avgtime/1000 - (avgtime/60000)*60, avgtime%1000) << std::endl;
  return;
}

