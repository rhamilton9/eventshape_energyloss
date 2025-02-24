T// A Method for testing a possibly novel histogram fitting procedure
// Each code contains references where applicable
//
// Created by Ryan Hamilton on 4/18/24.

#include "root_draw_tools.h"
#include "hist_tools.h"

void test_fit() {
  TFile* fin = new TFile("test_data/test.root");
//  TH1D* hist = fin->Get<TH1D>("test");
  const double lowbound = -10;
  const double highbound = 10;
  
  TH1D* hist = new TH1D("testhist_gaus_largestat", ";x;y", 10, lowbound, highbound);
  
  double mu_truth = 0.1;
  double sig_truth = 1.8;
  int n_counts = 1500;
  
  TF1* parentDist = new TF1("truth_dist", "(1./( [1] * sqrt(2*pi) )) * exp( -0.5*( (x - [0]) / [1] )**2 )", lowbound, highbound);
  parentDist->SetParameters(mu_truth, sig_truth);
  hist->FillRandom("truth_dist", n_counts);
  hist->Scale(1./hist->Integral("width"));
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kLake);
  
  TCanvas* canvas = new TCanvas();
  canvas->SetWindowSize(1000, 600);
  canvas->Divide(2, 2);
  
  canvas->cd(1);
  hist->Draw("hist");
  
  canvas->cd(2);
  //Fit to gaus
  TF1* fitDist = new TF1("fitted_dist", "(1./( [1]*sqrt(2*pi) )) * exp( -0.5*( (x - [0]) / [1] )**2 )", lowbound, highbound);
  fitDist->SetParLimits(0, -3, 10);
  fitDist->SetParLimits(1, 0.2, 3);
  TFitResultPtr pdf_fit = hist->Fit("fitted_dist","SLI");
  setStyleLine(parentDist, "thin dashed gray");
  parentDist->Draw("same");
  
  std::cout << "Abs. Error in mean :: " << TMath::Abs(pdf_fit->Parameter(0) - mu_truth) << std::endl;
  std::cout << "Abs. Error in stdev :: " << TMath::Abs(pdf_fit->Parameter(1) - sig_truth) << std::endl;
  
  canvas->cd(3);
  TGraph* cdf_graph = drawCDF_new(hist);
  
  canvas->cd(4);
  TGraph* to_fit_cdf = static_cast<TGraph*>(cdf_graph->Clone());
  to_fit_cdf->Draw("a* pmc");
  // Fit ERF on CDF, compare fit params
  TF1* cdfDist = new TF1("cdf_dist", "0.5*(1+erf( (x - [0])/([1]*sqrt(2)) ))", lowbound, highbound);
  cdfDist->SetParLimits(0, -3, 10);
  cdfDist->SetParLimits(1, 0.2, 3);
  TFitResultPtr cdf_fit = to_fit_cdf->Fit(cdfDist, "SP", "", lowbound, highbound);
  double store_mu = cdfDist->GetParameter(0);
  double store_sig = cdfDist->GetParameter(1);
  cdfDist->SetParameters(mu_truth, sig_truth);
  setStyleLine(cdfDist, "thin dashed gray");
  cdfDist->Draw("same");
  
  std::cout << "Abs. Error in mean :: " << TMath::Abs(store_mu - mu_truth) << std::endl;
  std::cout << "Abs. Error in stdev :: " << TMath::Abs(store_sig - sig_truth) << std::endl;
  
  canvas->SaveAs("test_plots/fitting_test.pdf");
}
