// Testing macro for all utils
// Runs through all key utils, ensures all are working properly.

void make_glauber_plots() {
  TFile *fin = new TFile("glauber_withgrid_b4.91-6.94.root");
  
  TH2D* edensity_default = fin->Get<TH2D>("inited_event2");
  TH2D* thickness_A = fin->Get<TH2D>("initA_event2");
  TH2D* thickness_B = fin->Get<TH2D>("initB_event2");
  
  TCanvas *canvas = new TCanvas();
  canvas->SetWindowSize(500, 500);
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kSunset);
  
  thickness_A->Draw("col");
  canvas->SaveAs("tA.pdf");
  
  thickness_B->Draw("col");
  canvas->SaveAs("tB.pdf");
  
  edensity_default->Draw("col");
  canvas->SaveAs("arith.pdf");
  
  TH2D* edensity_geom = static_cast<TH2D*>(edensity_default->Clone());
  edensity_geom->Reset();
  
  double val;
  for(int i = 1; i <= edensity_geom->GetXaxis()->GetNbins(); ++i) {
    for(int j = 1; j <= edensity_geom->GetYaxis()->GetNbins(); ++j) {
      val = thickness_A->GetBinContent(i, j) * thickness_B->GetBinContent(i, j);
      if (val <= 0) continue;
      edensity_geom->SetBinContent(i, j, TMath::Sqrt(val) );
    }
  }edensity_geom->Draw("col");
  canvas->SaveAs("geometric.pdf");
}
