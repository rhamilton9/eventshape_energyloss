// A simple macro to reproduce the plots in Ref.:
// http://dx.doi.org/10.1016/j.physletb.2013.01.051
// to make sure the data was reconstructed correctly.

#include "../../utils/root_draw_tools.h"

void ALICE_datacheck() {
  TFile* infile_PbSpectra = new TFile("ALICE_PbPb2.76TeV_unpacked.root");
  TFile* infile_ppReference = new TFile("../_ppreference/pp_reference.root");
  
  TH1F* sigref_pp = static_cast<TH1F*>(infile_ppReference->Get("Hist1D_y2_1"));
  sigref_pp->SetLineStyle(7);
  sigref_pp->SetLineColor(kRed);
  sigref_pp->SetLineWidth(2);
  
  TH1F* leadspectra[15];
  for (int i = 0; i < 15; ++i)
    leadspectra[i] = static_cast<TH1F*>(infile_PbSpectra->Get(Form("Hist1D_y1_%i",i+1)));
  
  char spectrum_centrality[15][10] = {
    "0-5%", "5-10%", "10-20%", "20-30%", "30-40%",
    "40-50%", "50-60%", "60-70%", "70-80%", "0-10%",
    "0-20%", "20-40%", "40-60%", "40-80%", "60-80%",
  };
  
  float T_AA[9] = {26.4, 20.6, 14.4, 8.7, 5.0, 2.68, 1.32, 0.59, 0.24};
  Int_t marker[9] = {20, 21, 22, 23, 24, 25, 26, 27, 32};
  
  TCanvas *canvas = new TCanvas();
  canvas->SetWindowSize(650, 1500);
  canvas->SetMargin(0.2, 0.05, 0.1, 0.05);
  canvas->SetLogy();
  canvas->SetLogx();
  gStyle->SetOptStat(0);
  
  TH1F* cref;
  TLegend* leg = new TLegend(0.25, 0.15, 0.7, 0.35);
  leg->SetLineWidth(0);
  for (int i = 8; i >= 0; --i) {
    
    leadspectra[i]->SetMarkerStyle(marker[i]);
    leadspectra[i]->Scale(TMath::Power(10, 16-2*i));
    leadspectra[i]->GetYaxis()->SetRangeUser(1e-16, 1e22);
    leadspectra[i]->SetTitle("");
    leadspectra[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    leadspectra[i]->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) (2^{2} N_{ch}) / (d#eta dp_{T}) [GeV/c]^{-2}");
    leadspectra[i]->Draw("hist p same");
    leg->AddEntry(leadspectra[i], Form("%s (x10^{%i})", spectrum_centrality[i], 16-2*i), "p");
    
    
    cref = static_cast<TH1F*>(sigref_pp->Clone());
    cref->Scale(T_AA[i]*TMath::Power(10, 16-2*i));
    //    cref->GetYaxis()->SetRangeUser(1e-16, 1e22);
    //    cref->GetXaxis()->SetRangeUser(leadspectra[8]->GetXaxis()->GetBinLowEdge(1),
    //                                   leadspectra[8]->GetXaxis()->GetBinLowEdge(leadspectra[8]->GetXaxis()->GetNbins()+1));
    //    cref->SetTitle("");
    //    cref->GetXaxis()->SetTitle("p_{T} [GeV]");
    //    cref->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) (2^{2} N_{ch}) / (d#eta dp_{T}) [GeV/c]^{-2}");
    cref->Draw("hist c same");
    
    
  }
  drawText("#bf{ALICE, Pb-Pb, #sqrt{S_{NN}} = 2.76 TeV}", 0.55, 0.915, false, kBlack, 0.026);
  drawText("Charged particles, |#eta| < 0.8", 0.55, 0.885, false, kBlack, 0.026);
  
  leg->Draw();
  canvas->SaveAs("ALICE_2.76TeV_spectracheck.pdf");
}
