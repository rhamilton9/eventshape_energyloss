// ROOT macro. Generates energyloss estimate by shifting reference pT spectrum and taking a ratio with original spectrum
// This estimate for energy loss is obtained by fitting this ratio against a known reference R_AA
// Assumes relevant data spectra have been correctly formatted and are accessible in ../data.nosync

#include "../config.h"
#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"

// TODO Implement KS test instead of Chi-square for better results without hist rebinning.

void energyloss_pTspectra() {
  const int n_dpt = (int) max_dpt / dpT_resolution;
  const int n_raa = sizeof(centrality_list) / sizeof(int[2]);
  char species[5];
  snprintf(species, 5, "%s%s", speciesA, speciesB);
  
  TFile *infile_raa = new TFile(Form("../data.nosync/%s_%.2fTeV/spectra/%s_%s%.2fTeV_unpacked.root",
                                species, sqrt_s, experiment, species, sqrt_s));
  TFile *infile_ref = new TFile(Form("../data.nosync/pp_reference_%.2fTeV.root", sqrt_s));
  
  TH1F *pp_reference = static_cast<TH1F*>(infile_ref->Get("Hist1D_y2_1"));
  TH1F *opt_raa;
  TH1F *opt_pt;
  TH1F* lead_raa[n_raa];
  TH1F* lead_totalunc[n_raa];
  for (int i = 0; i < n_raa; ++i) {
    lead_raa[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_%i",i+16)));
    lead_totalunc[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_e1_%i",i+16)));
    lead_totalunc[i]->Add(lead_totalunc[i],
                          static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_e2_%i",i+16))));
    lead_totalunc[i]->Add(lead_totalunc[i],
                          static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_e3_%i",i+16))));
  }
  
  // Open output file or create it if it doesn't exist.
  // Create TTree for energyloss data
  double energyloss;
  double botcentbin_local;
  double topcentbin_local;
  TFile* datfile = new TFile(Form("../data.nosync/%s_%.2fTeV/eventshape_energyloss.root", species, sqrt_s), "update");
  TTree* energyloss_tree = new TTree("tempname", "energyloss_tree");
  energyloss_tree->Branch("centralitybin_low", &botcentbin_local);
  energyloss_tree->Branch("centralitybin_high", &topcentbin_local);
  energyloss_tree->Branch("energyloss", &energyloss);
  
  // Canvas setup
  TCanvas *canvas = new TCanvas();
  canvas->SetWindowSize(1500, 500);
  canvas->Divide(3,1);
  
  const double dpt = -dpT_resolution;
  const int nbin = pp_reference->GetXaxis()->GetNbins();
  
  gStyle->SetOptStat(0);
  
  double cval;
  double cmin_dpt;
  for (int iCent = 0; iCent < n_raa; ++iCent) {
    double cmin_chi2 = 1e50;
    TH1F* chist = static_cast<TH1F*>(pp_reference->Clone());
    TH1F* craa = static_cast<TH1F*>(chist->Clone());
    TH1F* chi2_hist = new TH1F("",";#Delta p_{T};#chi_{2}",n_dpt, 0, n_dpt*(-dpt));
    for (int i = 0; i < n_dpt; ++i) {
      // shift by dpt, take ratio.
      chist = static_cast<TH1F*>(translateHist(pp_reference, 12, dpt*(i+1)));
      craa->Reset();
      for (int ibin = 1; ibin <= nbin; ++ibin) {
        cval = chist->GetBinContent(ibin);
        if (cval == 0) {craa->SetBinContent(ibin, 0); continue;}
        craa->SetBinContent(ibin, cval / pp_reference->GetBinContent(ibin));
      }
      // Compute chi^2 against raa
      double chi2 = 0;
      double cpt;
      for (int ibin = 1; ibin <= lead_raa[iCent]->GetXaxis()->GetNbins(); ++ibin) {
        cpt = lead_raa[iCent]->GetXaxis()->GetBinCenter(ibin);
        if (lead_raa[iCent]->GetXaxis()->GetBinLowEdge(ibin+1) > 5)
          chi2 += TMath::Power(craa->GetBinContent(craa->FindBin(cpt)) - lead_raa[iCent]->GetBinContent(ibin), 2)
                              / lead_totalunc[iCent]->GetBinContent(ibin);
      }chi2_hist->Fill(i*(-dpt), chi2);
      
      if (chi2 < cmin_chi2) {
//        cout << "new best chi2 = " << chi2 << " found at dpt = " << -(i+1)*dpt << endl;
        cmin_chi2 = chi2;
        cmin_dpt = -(i+1)*dpt;
        opt_pt = static_cast<TH1F*>(chist->Clone());
        opt_raa = static_cast<TH1F*>(craa->Clone());
      }
    }
    canvas->cd(1);
//    gPad->SetLogx();
    gPad->SetLogy();
    pp_reference->Draw("hist");
    opt_pt->Draw("hist same");
    
    canvas->cd(2);
    if (dpt < 0.05) chi2_hist->Draw("hist l");
    else chi2_hist->Draw("hist");
    
    canvas->cd(3);
    //  gPad->SetLogx();
    //  gPad->SetLogy();
    lead_raa[iCent]->SetLineColor(kRed);
    lead_raa[iCent]->Draw("hist");
    opt_raa->Draw("hist same");
    
    canvas->cd();
    TLatex* tex = drawText(Form("Cent %i-%i%%, #Delta p_{T} = %.4f",
                                centrality_list[iCent][0], centrality_list[iCent][1], cmin_dpt), 0.4, 0.92);
    canvas->SaveAs(Form("../plots/%s_%.2fTeV/%s/fitting/dptFit_cent%i-%i%%.pdf",
                        species, sqrt_s, experiment, centrality_list[iCent][0], centrality_list[iCent][1]));
    tex->Clear();
    
    // Write to tree
    energyloss = cmin_dpt;
    botcentbin_local = centrality_list[iCent][0];
    topcentbin_local = centrality_list[iCent][1];
    energyloss_tree->Fill();
  }energyloss_tree->Write("energyloss_tree", TObject::kOverwrite);
  return;
}
