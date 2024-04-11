// Testing macro for all utils
// Runs through all key utils, ensures all are working properly.

#include "root_draw_tools.h"
#include "hist_tools.h"


TH2D* trn_A;

double shift_unf = 0.5;
double shift_nnf = 2;
double shift_ref = 0.1;

void test_cdf();
void test_ks();
void test_translate();
void test_translate_simple();
void test_rotate_simple();

void util_diagnostics() {
  test_ks();
  return;
}

// Test for KS test
void test_ks() {
  TFile *fin = new TFile("test_data/test.root");
  TH1D* testhist_unf = fin->Get<TH1D>("test");
  TH1D* testhist_gaus = fin->Get<TH1D>("test_gaus_2");
  TH1D* testhist_nnf = fin->Get<TH1D>("nonunif_bins");
  TH1D* testhist_sml = fin->Get<TH1D>("test_small");
  TH1D* testhist_ref = fin->Get<TH1D>("test_ref");
  
  
  TFile *infile_ppReference = new TFile("test_data/pp_reference.root");
  TH1F* sigref_pp = static_cast<TH1F*>(infile_ppReference->Get("Hist1D_y2_1"));
  
  TH1* hist_totest_1 = testhist_unf;
//  TH1* hist_totest_2 = static_cast<TH1D*>(testhist_unf->Clone());
  TH1* hist_totest_2 = sigref_pp;
  
//  KS_statistic(testhist_sml, translateHist(testhist_sml, -0.1));
  cout << KS_statistic_new(hist_totest_1, hist_totest_2, 2, 0, 7, true, "test_plots/ks_new.pdf") << endl;
//  cout << KS_statistic(hist_totest_1, hist_totest_2, -4, -4, false, true, "test_plots/ks_old.pdf") << endl;
  
//  TH1* hist1,
//  TH1* hist2,
//  double horizShiftOnHist1 = 0,
//  double comparison_threshold = INT_MIN,
//  bool suppress_CDF_below_threshold = false,
//  bool doPlot = false,
//  char *saveName = (char*)"ks",
//  int iteration = -1
  
  
  
  TH1D* hist1 = new TH1D("test1",";x;y",10, -5, 5);
  TH1D* hist2 = new TH1D("test2",";x;y",10, -5.5, 4.5);
  hist1->FillRandom("gaus", 10000);
  hist2->FillRandom("gaus", 10000);
//  KS_statistic(hist1, translateHist(hist2, -5));
//  KS_statistic(hist1, hist2);
}

// Test for CDF plotting
void test_cdf() {
  TFile *fin = new TFile("test_data/test.root");
  TH1D* testhist_unf = fin->Get<TH1D>("test");
  TH1D* testhist_nnf = fin->Get<TH1D>("nonunif_bins");
  
  TFile *infile_ppReference = new TFile("test_data/pp_reference.root");
  TH1F* sigref_pp = static_cast<TH1F*>(infile_ppReference->Get("Hist1D_y2_1"));
  
  TCanvas *canvas = new TCanvas();
  canvas->SetWindowSize(1000, 500);
  canvas->Divide(2, 1);
  
  double thresh_min = 4;
  double thresh_max = 4.8;
  
  canvas->cd(1);
  testhist_unf->Draw("hist");
  
  canvas->cd(2);
  TMultiGraph* graph = new TMultiGraph();
  TGraph* cgraph;
  cgraph = drawCDF_new(testhist_unf, 0, thresh_min, thresh_max, false, true);
  cgraph->SetLineColor(kBlue);
  graph->Add(cgraph);
  cgraph = drawCDF(testhist_unf, 0, false, true);
  cgraph->SetLineColor(kRed);
  graph->Add(cgraph);
  graph->Draw("al");
  TLine *rangeLine = new TLine();
  setStyleLine(rangeLine, "thin dashed gray");
  rangeLine->DrawLine(thresh_min, 0, thresh_min, 1);
  rangeLine->DrawLine(thresh_max, 0, thresh_max, 1);
  canvas->SaveAs("test_plots/CDF_gaus.pdf");
  
  canvas->cd(1);
  gPad->SetLogy();
  gPad->SetLogx();
  sigref_pp->Draw("hist");
  
  canvas->cd(2);
  graph = new TMultiGraph();
  gPad->SetLogx();
  cgraph = drawCDF_new(sigref_pp, 0, thresh_min, thresh_max, false, true);
  cgraph->SetLineColor(kBlue);
  graph->Add(cgraph);
  cgraph = drawCDF(sigref_pp, 0, false, true);
  cgraph->SetLineColor(kRed);
  graph->Add(cgraph);
  graph->Draw("al");
  rangeLine->DrawLine(thresh_min, 0, thresh_min, 1);
  rangeLine->DrawLine(thresh_max, 0, thresh_max, 1);
  canvas->SaveAs("test_plots/CDF_ppReference.pdf");
}

// Testing for more accurate translate macro
void test_translate() {
  TFile *fin = new TFile("test_data/test.root");
  TFile *infile_ppReference = new TFile("test_data/pp_reference.root");
  TH1D* testhist_unf = static_cast<TH1D*>(fin->Get("test"));
  TH1D* testhist_nnf = static_cast<TH1D*>(fin->Get("nonunif_bins"));
  TH1F* sigref_pp = static_cast<TH1F*>(infile_ppReference->Get("Hist1D_y2_1"));
  
  TH1D* transhist_unf = static_cast<TH1D*>(translateHist(testhist_unf, shift_unf));
  TH1D* transhist_nnf = static_cast<TH1D*>(translateHist(testhist_nnf, shift_nnf));
  TH1F* transhist_ref = static_cast<TH1F*>(translateHist(sigref_pp, shift_ref));
  
  transhist_unf->SetLineColor(kRed);
  transhist_nnf->SetLineColor(kRed);
  transhist_ref->SetLineColor(kRed);
  
  TH1D* transback_unf = static_cast<TH1D*>(translateHist(transhist_unf, -shift_unf));
  TH1D* transback_nnf = static_cast<TH1D*>(translateHist(transhist_nnf, -shift_nnf));
  TH1F* transback_ref = static_cast<TH1F*>(translateHist(transhist_ref, -shift_ref));
  
  transback_unf->SetLineColor(kGray);
  transback_unf->SetLineStyle(8);
  //  transback_unf->SetLineWidth(2);
  transback_nnf->SetLineColor(kGray);
  transback_nnf->SetLineStyle(8);
  //  transback_nnf->SetLineWidth(2);
  transback_ref->SetLineColor(kGray);
  transback_ref->SetLineStyle(8);
  //  transback_ref->SetLineWidth(2);
  
  gStyle->SetOptStat(0);
  const int nbins[3] = {testhist_unf->GetXaxis()->GetNbins(), testhist_nnf->GetXaxis()->GetNbins(), sigref_pp->GetXaxis()->GetNbins()};
  double bins_unf[nbins[0]+1];
  double bins_nnf[nbins[1]+1];
  double bins_ref[nbins[2]+1];
  for (int i = 0; i <= nbins[0]; ++i)
    bins_unf[i] = testhist_unf->GetBinLowEdge(i+1)+shift_unf;
  for (int i = 0; i <= nbins[1]; ++i)
    bins_nnf[i] = testhist_nnf->GetBinLowEdge(i+1)+shift_nnf;
  for (int i = 0; i <= nbins[2]; ++i)
    bins_ref[i] = sigref_pp->GetBinLowEdge(i+1)+shift_ref;
  
  TH1D* pure_unf = new TH1D("","", nbins[0], bins_unf);
  for (int i = 1; i <= nbins[0]; ++i) pure_unf->SetBinContent(i, testhist_unf->GetBinContent(i));
  TH1D* pure_nnf = new TH1D("","", nbins[1], bins_nnf);
  for (int i = 1; i <= nbins[1]; ++i) pure_nnf->SetBinContent(i, testhist_nnf->GetBinContent(i));
  TH1D* pure_ref = new TH1D("","", nbins[2], bins_ref);
  for (int i = 1; i <= nbins[2]; ++i) pure_ref->SetBinContent(i, sigref_pp->GetBinContent(i));
  
  pure_unf->SetLineColor(kAzure+4);
  pure_unf->SetLineWidth(2);
  pure_unf->SetLineStyle(3);
  pure_nnf->SetLineColor(kAzure+4);
  pure_nnf->SetLineWidth(2);
  pure_nnf->SetLineStyle(3);
  pure_ref->SetLineColor(kAzure+4);
  pure_ref->SetLineWidth(2);
  pure_ref->SetLineStyle(3);
  
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->AddEntry(testhist_unf, "Original Hist", "l");
  leg->AddEntry(pure_unf, Form("Translation %.2f (No Rebin)", shift_unf), "l");
  leg->AddEntry(transhist_unf, "After Rebinning", "l");
  leg->AddEntry(transback_unf, Form("Restoring shift -%.2f", shift_unf), "l");
  
  TCanvas *c = new TCanvas();
  c->SetWindowSize(1500, 500);
  c->Divide(3, 1);
  c->cd(1);
  testhist_unf->GetYaxis()->SetRangeUser(0, 1.2*getMaxFromHists(transhist_unf, testhist_unf));
  testhist_unf->Draw("hist");
  pure_unf->Draw("hist same");
  transhist_unf->Draw("hist same");
  transback_unf->Draw("hist same");
  leg->Draw();
  
  leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->AddEntry(testhist_nnf, "Original Hist", "l");
  leg->AddEntry(pure_nnf, Form("Translation %.2f (No Rebin)", shift_nnf), "l");
  leg->AddEntry(transhist_nnf, "After Rebinning", "l");
  leg->AddEntry(transback_nnf, Form("Restoring shift -%.2f", shift_nnf), "l");
  
  c->cd(2);
  testhist_nnf->GetYaxis()->SetRangeUser(0, 1.2*getMaxFromHists(transhist_nnf, testhist_nnf));
  testhist_nnf->Draw("hist");
  pure_nnf->Draw("hist same");
  transhist_nnf->Draw("hist same");
  transback_nnf->Draw("hist same");
  leg->Draw();
  
  leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->AddEntry(sigref_pp, "Original Hist", "l");
  leg->AddEntry(pure_ref, Form("Translation %.2f (No Rebin)", shift_ref), "l");
  leg->AddEntry(transhist_ref, "After Rebinning", "l");
  leg->AddEntry(transback_ref, Form("Restoring shift -%.2f", shift_ref), "l");
  
  c->cd(3);
  gPad->SetLogy();
  gPad->SetLogx();
  sigref_pp->SetLineColor(kBlack);
  sigref_pp->Draw("hist");
  pure_ref->Draw("hist same");
  transhist_ref->Draw("hist same");
  transback_ref->Draw("hist same");
  leg->Draw();
  
  c->SaveAs("test_plots/test_translate.pdf");
  return;
}

// Testing for simple translate macro
void test_translate_simple() {
  TFile *fin_glauber = new TFile("test_data/glauber_withgrid_b0.00-1.00.root");
  TH2D* test_A = static_cast<TH2D*>(fin_glauber->Get("initA_event1"));
  TH2D* test_B = static_cast<TH2D*>(fin_glauber->Get("initB_event0"));
  
  TCanvas* c0 = new TCanvas();
  test_A->Draw("colz");
  std::vector<double> cm = getCM(test_A, 21, true);
  TMarker *cm_marker = new TMarker();
  cm_marker->SetMarkerColor(kRed+2);
  cm_marker->SetMarkerStyle(5);
  cm_marker->DrawMarker(cm.at(0), cm.at(1));

  TCanvas* c1 = new TCanvas();
  TH2D* CM_A = translateHist_simple(test_A);
  CM_A->Draw("colz");
  cm = getCM(CM_A, 21, true);
  cm_marker->DrawMarker(cm.at(0), cm.at(1));

  TCanvas* c2 = new TCanvas();
  trn_A = translateHist_simple(test_A, 3, -2);
  trn_A->Draw("colz");
  cm = getCM(trn_A, 21, true);
  cm_marker->DrawMarker(cm.at(0), cm.at(1));
  return;
}

// Testing for simple rotate macro
void test_rotate_simple() {
  TFile *fin_glauber = new TFile("test_data/glauber_withgrid_b0.00-1.00.root");
  TH2D* test_A = static_cast<TH2D*>(fin_glauber->Get("initA_event1"));
  TH2D* test_B = static_cast<TH2D*>(fin_glauber->Get("initB_event0"));
  
  TCanvas* c0 = new TCanvas();
  test_A->Draw("colz");
  std::vector<double> cm = getCM(test_A, 21, true);
  TMarker *cm_marker = new TMarker();
  cm_marker->SetMarkerColor(kRed+2);
  cm_marker->SetMarkerStyle(5);
  cm_marker->DrawMarker(cm.at(0), cm.at(1));
  
  TCanvas* c1 = new TCanvas();
  TH2D* CM_A = rotateHist2D_simple(translateHist_simple(test_A), 1.24);
  CM_A->Draw("colz");
  cm = getCM(CM_A, 21, true);
  cm_marker->DrawMarker(cm.at(0), cm.at(1));
  
  TCanvas* c2 = new TCanvas();
  trn_A = rotateHist2D_simple(translateHist_simple(test_A, 3, -2), 4.81);
  trn_A->Draw("colz");
  cm = getCM(trn_A, 21, true);
  cm_marker->DrawMarker(cm.at(0), cm.at(1));
  return;
}
