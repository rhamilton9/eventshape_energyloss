// Testing macro for all utils
// Runs through all key utils, ensures all are working properly.

#include "root_draw_tools.h"
#include "hist_tools.h"

void testDrawTools() {
  
}

void testHistTools() {
  
  
}


TH2D* trn_A;

double shift_unf = 0.5;
double shift_nnf = 2;
double shift_ref = 0.1;

void util_diagnostics() {
  // Testing for more accurate translate macro
  TFile *fin = new TFile("test.root");
  TFile *infile_ppReference = new TFile("pp_reference.root");
  TH1D* testhist_unf = static_cast<TH1D*>(fin->Get("test"));
  TH1D* testhist_nnf = static_cast<TH1D*>(fin->Get("nonunif_bins"));
  TH1F* sigref_pp = static_cast<TH1F*>(infile_ppReference->Get("Hist1D_y2_1"));
  
  TH1D* transhist_unf = static_cast<TH1D*>(translateHist(testhist_unf, 11, shift_unf));
  TH1D* transhist_nnf = static_cast<TH1D*>(translateHist(testhist_nnf, 11, shift_nnf));
  TH1F* transhist_ref = static_cast<TH1F*>(translateHist(sigref_pp, 12, shift_ref));
  
  transhist_unf->SetLineColor(kRed);
  transhist_nnf->SetLineColor(kRed);
  transhist_ref->SetLineColor(kRed);
  
  TH1D* transback_unf = static_cast<TH1D*>(translateHist(transhist_unf, 11, -shift_unf));
  TH1D* transback_nnf = static_cast<TH1D*>(translateHist(transhist_nnf, 11, -shift_nnf));
  TH1F* transback_ref = static_cast<TH1F*>(translateHist(transhist_ref, 12, -shift_ref));
  
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
  
  c->SaveAs("test_translate.pdf");
  
  
//  TFile *fin_glauber = new TFile("glauber_withgrid_b0.00-1.00.root");
//  TH2D* test_A = static_cast<TH2D*>(fin_glauber->Get("initA_event1"));
//  TH2D* test_B = static_cast<TH2D*>(fin_glauber->Get("initB_event0"));
//  
//  // Testing for simple translate macro
//  TCanvas* c0 = new TCanvas();
//  test_A->Draw("colz");
//  std::vector<double> cm = getCM(test_A, 21, true);
//  TMarker *cm_marker = new TMarker();
//  cm_marker->SetMarkerColor(kRed+2);
//  cm_marker->SetMarkerStyle(5);
//  cm_marker->DrawMarker(cm.at(0), cm.at(1));
//  
//  TCanvas* c1 = new TCanvas();
//  TH2D* CM_A = translateHist_simple(test_A);
//  CM_A->Draw("colz");
//  cm = getCM(CM_A, 21, true);
//  cm_marker->DrawMarker(cm.at(0), cm.at(1));
//  
//  TCanvas* c2 = new TCanvas();
//  trn_A = translateHist_simple(test_A, 3, -2);
//  trn_A->Draw("colz");
//  cm = getCM(trn_A, 21, true);
//  cm_marker->DrawMarker(cm.at(0), cm.at(1));
  
  //Testing for simple rotate macro
//  TCanvas* c0 = new TCanvas();
//  test_A->Draw("colz");
//  std::vector<double> cm = getCM(test_A, 21, true);
//  TMarker *cm_marker = new TMarker();
//  cm_marker->SetMarkerColor(kRed+2);
//  cm_marker->SetMarkerStyle(5);
//  cm_marker->DrawMarker(cm.at(0), cm.at(1));
//  
//  TCanvas* c1 = new TCanvas();
//  TH2D* CM_A = rotateHist2D_simple(translateHist_simple(test_A), 1.24);
//  CM_A->Draw("colz");
//  cm = getCM(CM_A, 21, true);
//  cm_marker->DrawMarker(cm.at(0), cm.at(1));
//  
//  TCanvas* c2 = new TCanvas();
//  trn_A = rotateHist2D_simple(translateHist_simple(test_A, 3, -2), 4.81);
//  trn_A->Draw("colz");
//  cm = getCM(trn_A, 21, true);
//  cm_marker->DrawMarker(cm.at(0), cm.at(1));
}
