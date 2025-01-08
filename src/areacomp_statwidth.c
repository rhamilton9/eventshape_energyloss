// A short method for comparing the statwith from EBE against the pheno average.
//

#include "../config.h"
#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"
#include "../utils/glauber_tools.h"

int global_event = 0;

bool killswitch = false;

const int centrality_list_local[9][2] = {
  { 0,  5}, { 5, 10}, {10, 20}, {20, 30}, {30, 40},
  {40, 50}, {50, 60}, {60, 70}, {70, 80}
};

const char phenodesc[2][3][45] = {
  {"T_{A} + T_{B}, no align", "T_{A} + T_{B}, CM align", "T_{A} + T_{B}, CM+#Psi_{2} align"},
  {"#sqrt{T_{A}T_{B}}, no align", "#sqrt{T_{A}T_{B}}, CM align", "#sqrt{T_{A}T_{B}}, CM+#Psi_{2} align"}
};



void areacomp_statwidth() {
  std::cout << "Warning in <src/areacomp_statwidth>: this calculation is slow, due to many nested loops for each TH2D." << std::endl;
  gStyle->SetOptStat(0);
  
  const int ncent = sizeof(centrality_list_local) / sizeof(int[2]);
  std::vector<float> blist;
  std::vector<float> centrality_binedge;
  for (int i = 0; i < ncent; ++i) {
    if (blist.size() == 0 || centrality_binedge.back() != centrality_list_local[i][0]) {
      centrality_binedge.push_back(centrality_list_local[i][0]);
      blist.push_back(getImpactParameterFromCentralityClass(centrality_list_local[i][0], sqrt_s, speciesA, speciesB));
    }
    centrality_binedge.push_back(centrality_list_local[i][1]);
    blist.push_back(getImpactParameterFromCentralityClass(centrality_list_local[i][1], sqrt_s, speciesA, speciesB));
  }
  
  TFile* outfile = new TFile("../ebe_area_hist.root", "recreate");
  
  // Setup hists
  const char scaling_str[2][5] = {"arit", "geom"};
  const char alignstr[3][8] = {"noalign", "cmalign", "aligned"};
  TH1D* areaEBE_hist[6];
  TH1D* areaPhn_hist[6]; 
  for (int i = 0; i < 6; ++i) {
    areaEBE_hist[i] = new TH1D(Form("area_hist_EBE_%s_%s", scaling_str[i/3],alignstr[i%3]),
                               ";Centrality [%];Area [fm^{2}]", ncent, centrality_binedge.data());
    areaPhn_hist[i] = new TH1D(Form("area_hist_Pheno_%s_%s", scaling_str[i/3],alignstr[i%3]),
                               ";Centrality [%];Area [fm^{2}]", ncent, centrality_binedge.data());
    
    areaEBE_hist[i]->SetLineColor(kBlack);
    areaEBE_hist[i]->SetMarkerColor(kBlack);
    areaEBE_hist[i]->SetMarkerStyle(34);
    areaEBE_hist[i]->GetYaxis()->SetRangeUser(2, 175);
    
    areaPhn_hist[i]->SetLineColor(kCyan+2);
    areaPhn_hist[i]->SetMarkerColor(kCyan+2);
    areaPhn_hist[i]->SetMarkerStyle(43);
  }
  
  // Gather data
  for (int iCent = 0; iCent < ncent; ++iCent) {
    
    // Get relevant file for EBE calculation
    float b_bot = getImpactParameterFromCentralityClass(centrality_list_local[iCent][0], sqrt_s, speciesA, speciesB);
    float b_top = getImpactParameterFromCentralityClass(centrality_list_local[iCent][1], sqrt_s, speciesA, speciesB);
    TFile* cEBE_file = new TFile(Form("../data.nosync/%s%s_%.2fTeV/glauber/glauber_withgrid_b%.2f-%.2f.root",
                               speciesA, speciesB, sqrt_s, b_bot, b_top));
    
    // Read hists from file, compute variances for each
    TH2D* chist_e;
    int nevent_cfile = cEBE_file->Get<TTree>("lemon")->GetBranch("npart")->GetEntries();
    double csum[6] = {0, 0, 0, 0, 0, 0};
    for (int iEvent = 1; iEvent < nevent_cfile; ++iEvent) {
      chist_e = cEBE_file->Get<TH2D>(Form("inited_event%i", iEvent));
      
      if (killswitch) return;
      
//      if (iEvent > 100) break;
      
//      // Arithmetic scaling, no alignment, reaction plane
//      areaEBE_hist[0]->Fill(centrality_list_local[iCent][0], varianceArea(chist_e));
//      
//      // Arithmetic scaling, CM alignment, reaction plane
//      chist_e = translateHist_simple(chist_e, 0, 0, true);
//      areaEBE_hist[1]->Fill(centrality_list_local[iCent][0], varianceArea(chist_e));
//      
//      // Arithmetic scaling, CM alignment, rotate to participant plane
//      chist_e = rotateHist2D_simple(chist_e, TMath::Pi()/4 - getParticipantPlaneAngle(chist_e), true);
//      areaEBE_hist[2]->Fill(centrality_list_local[iCent][0], varianceArea(chist_e));
//      
//      // Adjust to Geometric scaling
//      chist_e = geometricEnergyDensity(cEBE_file->Get<TH2D>(Form("initA_event%i", iEvent)),
//                                       cEBE_file->Get<TH2D>(Form("initB_event%i", iEvent)));
//      
//      // Geometric scaling, no alignment, reaction plane
//      areaEBE_hist[3]->Fill(centrality_list_local[iCent][0], varianceArea(chist_e));
//      
//      // Geometric scaling, CM alignment, reaction plane
//      chist_e = translateHist_simple(chist_e, 0, 0, true);
//      areaEBE_hist[4]->Fill(centrality_list_local[iCent][0], varianceArea(chist_e));
//      
//      // Geometric scaling, CM alignment, rotate to participant plane
//      chist_e = rotateHist2D_simple(chist_e, TMath::Pi()/4 - getParticipantPlaneAngle(chist_e), true);
//      areaEBE_hist[5]->Fill(centrality_list_local[iCent][0], varianceArea(chist_e));
      
      global_event = iEvent;
      
      // Arithmetic scaling, no alignment, reaction plane
      csum[0] += varianceArea(chist_e);
      
      // Arithmetic scaling, CM alignment, reaction plane
      chist_e = translateHist_simple(chist_e, 0, 0, true);
      csum[1] += varianceArea(chist_e);
      
      // Arithmetic scaling, CM alignment, rotate to participant plane
      chist_e = rotateHist2D_simple(chist_e, TMath::Pi()/4 - getParticipantPlaneAngle(chist_e), true);
      csum[2] += varianceArea(chist_e);
      
      // Adjust to Geometric scaling
      chist_e = geometricEnergyDensity(cEBE_file->Get<TH2D>(Form("initA_event%i", iEvent)),
                                       cEBE_file->Get<TH2D>(Form("initB_event%i", iEvent)));
      
      // Geometric scaling, no alignment, reaction plane
      csum[3] += varianceArea(chist_e);
      
      // Geometric scaling, CM alignment, reaction plane
      chist_e = translateHist_simple(chist_e, 0, 0, true);
      csum[4] += varianceArea(chist_e);
      
      // Geometric scaling, CM alignment, rotate to participant plane
      chist_e = rotateHist2D_simple(chist_e, TMath::Pi()/4 - getParticipantPlaneAngle(chist_e), true);
      csum[5] += varianceArea(chist_e);
    }
    
    // Rescale bin contents to obtain average area in the current centrality
    for (int iPheno = 0; iPheno < 6; ++iPheno) {
      areaEBE_hist[iPheno]->SetBinContent(iCent+1, csum[iPheno]/nevent_cfile);
//      areaEBE_hist[iPheno]->SetBinContent(iCent+1, areaEBE_hist[iPheno]->GetBinContent(iCent+1)/nevent_cfile);
      areaEBE_hist[iPheno]->SetBinError(iCent+1, 1e-10);
    }
    
    // Get relevant file for Pheno data and fill hist
    TFile* phenoFile = new TFile(Form("../data.nosync/%s%s_%.2fTeV/glauber/glauber_compiled_b%.2f-%.2f.root",
                               speciesA, speciesB, sqrt_s, b_bot, b_top));
    for (int iPheno = 0; iPheno < 6; ++iPheno) {
      areaPhn_hist[iPheno]->SetBinContent(iCent+1, varianceArea(phenoFile->Get<TH2D>(Form("b%.2f-%.2f_eventshape_%s_%s",
                                                                                          b_bot,b_top,alignstr[iPheno%3],scaling_str[iPheno/3]))));
      areaPhn_hist[iPheno]->SetBinError(iCent+1, 1e-10);
    }
    
    std::cout << "Finished with calculations in centrality class " << centrality_list_local[iCent][0] << "-" << centrality_list_local[iCent][1] << "%" << std::endl;
    cEBE_file->Close();
    delete cEBE_file;
    
    phenoFile->Close();
    delete phenoFile;
  }// End of loop over centrality
  
  // Canvas with area calculation information
  TCanvas* canvas_area = new TCanvas();
  canvas_area->SetCanvasSize(1000, 600);
  TPad* pads[2][3];
  double dx = 0.1*0.3/3;
  pads[0][0] = buildPad("pad0_arit", 0, .5, .95/3+2*dx, 1,           0.13, 0, 0, 0.11, false);
  pads[0][1] = buildPad("pad1_arit", .95/3+2*dx, .5, 1.9/3+dx, 1,    0, 0, 0, 0.11, false);
  pads[0][2] = buildPad("pad2_arit", 1.9/3+dx, .5, .95, 1,           0, 0, 0, 0.11, false);
  pads[1][0] = buildPad("pad0_geom", 0, 0, .95/3+2*dx, .5,             0.13, 0, 0.11, 0, false);
  pads[1][1] = buildPad("pad1_geom", .95/3+2*dx, 0, 1.9/3+dx, .5,      0, 0, 0.11, 0, false);
  pads[1][2] = buildPad("pad2_geom", 1.9/3+dx, 0, .95, .5,             0, 0, 0.11, 0, false);
  
  TLegend *leg = new TLegend(0.1, 0.7, 0.83, 0.85);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.043);
  leg->AddEntry(areaEBE_hist[0], "(Co-)Var. EBE Direct (#times4#pi)", "lp");
  leg->AddEntry(areaPhn_hist[0], "(Co-)Var. Phenom. (#times4#pi)", "lp");
  
  outfile->cd();
  for (int iPheno = 0; iPheno < 6; ++iPheno) {
    pads[iPheno/3][iPheno%3]->cd();
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    
    areaEBE_hist[iPheno]->GetXaxis()->SetLabelSize(0.045);
    areaEBE_hist[iPheno]->GetXaxis()->SetTitleSize(0.05);
    areaEBE_hist[iPheno]->GetXaxis()->SetTitleOffset(0.95);
    areaEBE_hist[iPheno]->GetYaxis()->SetLabelSize(0.045);
    areaEBE_hist[iPheno]->GetYaxis()->SetTitleSize(0.05);
    areaEBE_hist[iPheno]->GetYaxis()->SetTitleOffset(1.3);
    
    areaEBE_hist[iPheno]->Draw("hist p e0");
    areaPhn_hist[iPheno]->Draw("hist p e0 same");
    
    areaEBE_hist[iPheno]->Write();
    areaPhn_hist[iPheno]->Write();
    
    if (iPheno/3 == 0) drawText(phenodesc[iPheno/3][iPheno%3], gPad->GetLeftMargin()+0.05, 0.81, false, kBlack, 0.05);
    if (iPheno/3 == 1) drawText(phenodesc[iPheno/3][iPheno%3], gPad->GetLeftMargin()+0.05, 0.89, false, kBlack, 0.05);
    
    if (iPheno == 5) leg->Draw();
    if (iPheno == 3) {
      drawText(Form("#it{TGlauberMC} #bf{%s%s}",speciesA,speciesB), 0.95, 0.91, true, kBlack, 0.05);
      drawText(Form("#sqrt{s_{NN}} = %.2f TeV", sqrt_s), 0.95, 0.84, true, kBlack, 0.05);
    }
//    if (iPheno == 1) {
//      drawText("Events per Centrality", 0.85, 0.72, true, kBlack, 0.05);
//      drawText(Form("Incl. & Excl.: %0.2f", avg_events_per_cent), 0.96, 0.65, true, kBlack, 0.05);
//      drawText(Form("Phenomenological:  %0.2f", 5000.00), 0.95, 0.6, true, kBlack, 0.05);
//    }
  }canvas_area->SaveAs(Form("../plots/%s%s_%.2fTeV/glauber/glaubercomp_statwidth.pdf", speciesA, speciesB, sqrt_s));
  outfile->Close();
  return;
}
