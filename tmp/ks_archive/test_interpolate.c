// A Method for testing some interpolation code
// Each code contains references where applicable
//
// Created by Ryan Hamilton on 4/18/24.

#include "hist_tools.h"
#include "interpolation_tools.h"

void test_interpolate() {
  TFile* fin = new TFile("test_data/test.root");
//  TH1D* hist = fin->Get<TH1D>("test");
  TH1D* hist = new TH1D("testhist_gaus_largestat", ";x;y", 10, -5, 5);
  hist->FillRandom("gaus", 1000);
  hist->Scale(1./hist->Integral("width"));
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kLake);
  
  TCanvas* canvas = new TCanvas();
  canvas->SetWindowSize(1000, 600);
  canvas->Divide(2, 2);
  
  canvas->cd(1);
  hist->Draw("hist");
  
  canvas->cd(2);
  TMultiGraph *all = new TMultiGraph();
  std::vector<double> bins;
  std::vector<double> vals;
  int nbin = hist->GetXaxis()->GetNbins();
  for (int i = 1; i <= nbin; ++i) {
    bins.push_back(hist->GetBinCenter(i));
    vals.push_back(hist->GetBinContent(i));
  }TGraph* smoothedPDF = new TGraph(nbin, bins.data(), vals.data());
  all->Add(smoothedPDF, "c");
  
  bins.clear();
  vals.clear();
  for (int i = 0; i <= 100; ++i) {
    bins.push_back(-5+i*0.1);
    vals.push_back(TMath::Exp(- bins.at(i) * bins.at(i) / 2) / TMath::Sqrt(2 * TMath::Pi()) );
  }TGraph* truePDF = new TGraph(101, bins.data(), vals.data());
  all->Add(truePDF, "c");
  
  all->Draw("a plc");
  
  canvas->cd(3);
  drawCDF_new(hist);
  
  canvas->cd(4);
  TMultiGraph* cdf_all = new TMultiGraph();
  
  // Smooth artificially (only in plot)
  //  cdf_all->Add(drawCDF_new(hist, 0, INT_MIN, INT_MAX, false, true), "c");
  
  // Smooth manually
  std::vector<std::vector<double>> cdf = CDF_FromHist_new(hist, INT_MIN, INT_MAX);
  bins.clear();
  vals.clear();
  int nentries = cdf.size();
  for (int i = 0; i < nentries; ++i) {
    bins.push_back(cdf.at(i).at(0));
    vals.push_back(cdf.at(i).at(1));
  }
  TGraph* interp_cdf = new TGraph(nentries, bins.data(), vals.data());
//  interp_cdf->Draw("a*");
  interp_cdf = Smooth(interp_cdf);
  cdf_all->Add(interp_cdf);
  
  // Reference truth-level CDF
  bins.clear();
  vals.clear();
  for (int i = 0; i <= 100; ++i) {
    bins.push_back(-5+i*0.1);
    vals.push_back(0.5*(1+TMath::Erf(bins.at(i) / TMath::Sqrt(2))) );
  }TGraph *cdf_true = new TGraph(101, bins.data(), vals.data());
  cdf_all->Add(cdf_true, "c");
  cdf_all->Draw("a plc");
  
  canvas->SaveAs("test_plots/proof_of_concept.pdf");
}
