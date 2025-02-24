// ROOT macro. Generates pseudo-v2 plot by shifting jet pT spectrum and taking a difference with original spectrum
// Attempts to find mean energy loss by fitting the pseudo-v2 with a reference v2
// Assumes relevant file with input cuts has already been generated
// See ../eventgen/PYTHIA/eventgen_jets.cc for details on output data format

#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"
#include "../config.h"

// TODO lead pT store in TTree....

float pT_shift_best[2] = {16.29, 10.49}; // in GeV, the output of energyloss_pTspectra for the analogous R_AA
float c0_glauber[2][4] = {
  {3.8085, 5.2972, 5.5805, 5.1072},
  {2.5365, 3.4115, 3.0204, 2.9363}
};
float c2_glauber[2][4] = {
  {0.4152, 0.2258, 0.3098, 0.3147},
  {0.5281, 0.7487, 0.8464, 0.7532}
};
int plot_color[4] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1};

// Color and aesthetic settings for root plots
Int_t palette_full = kCherry;
Int_t palette_ch = kLake;
Int_t markerStyle_full = 20;
Int_t markerStyle_ch = 33;
Float_t markerSize_full = 0.75;
Float_t markerSize_ch = 1.0;

void energyloss_anisotropy() {
  char system[10];
  snprintf(system, 10, "%s%s", speciesA, speciesB);
  char data_directory[100];
  snprintf(data_directory, 100, "../data.nosync/%s_%.2fTeV", system, sqrt_s);
  char outplot_directory[100];
  snprintf(outplot_directory, 100, "../plots/%s_%.2fTeV/%s", system, sqrt_s, experiment);
  
  char chargestring[2][10] = {"full","char"};
  
  
  
  
  // TCanvas Setup
  gStyle->SetOptStat(0);
  TCanvas* canvas_check = new TCanvas();
  canvas_check->SetCanvasSize(350,500);
  
  //---------------------------------------------------------- Find data v2 spectra
  
  TFile* data_reference_file = new TFile(Form("%s/%s_%.2fTeV_data.root", data_directory, experiment, sqrt_s));
  
  // Which spectra to get from the data?
  const float data_lowestpT = 25;
  const int ndata_centbin = 2;
  char v2_centbin[ndata_centbin][10] = {"0-5","30-50"};
  TH1F* hist_v2_datreference[ndata_centbin];
  TH1F* hist_shifting[ndata_centbin];
  TH1F* ratio_shifting[ndata_centbin];
  TH1F* hist_baseline[ndata_centbin];
  int color_centbin[2] = {kBlue+1, kRed+1};
  for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent) {
    hist_v2_datreference[i_cent] = data_reference_file->Get<TH1F>(Form("%s_%.2fTeV_v2_cent%s",experiment, sqrt_s, v2_centbin[i_cent]));
    hist_v2_datreference[i_cent]->SetLineColor(color_centbin[i_cent]);
    hist_v2_datreference[i_cent]->SetMarkerColor(color_centbin[i_cent]);
    hist_v2_datreference[i_cent]->SetMarkerStyle(20);
    
    // To later shift spectra into this histogram with the same binning as reference
    hist_shifting[i_cent] = static_cast<TH1F*>(hist_v2_datreference[i_cent]->Clone());
    hist_shifting[i_cent]->Reset();
    ratio_shifting[i_cent] = static_cast<TH1F*>(hist_v2_datreference[i_cent]->Clone());
    ratio_shifting[i_cent]->Reset();
    
    // For taking ratio against pp with same binning
    hist_baseline[i_cent] = static_cast<TH1F*>(hist_v2_datreference[i_cent]->Clone());
    hist_baseline[i_cent]->SetLineColor(color_centbin[i_cent]+2);
    hist_baseline[i_cent]->Reset();
  }
  
  
  // Plot to check we have the right data
  hist_v2_datreference[0]->GetYaxis()->SetRangeUser(0, 0.2);
  hist_v2_datreference[0]->GetXaxis()->SetRangeUser(20, 120);
  hist_v2_datreference[0]->Draw("hist");
  hist_v2_datreference[1]->Draw("hist same");
  canvas_check->SaveAs(Form("%s/checks/v2_inputspectra_check.pdf",outplot_directory));
  
  
  //---------------------------------------------------------- Extract pythia simulated jets with final cuts
  
  TFile* pythia_file = new TFile(Form("%s/pythpp_%.2fTeV_r%.1fjet.root", data_directory, sqrt_s, jet_radius));
  
  // Get the output spectrum from the compiled file to make sure the tree-weighted one looks the same
  TH1D* pythia_spectra_uncut_reference = pythia_file->Get<TH1D>(Form("hist_pT_%sjet",chargestring[flag_use_charged_jets]));
  
  TH1D* pythia_spectra_uncut_check = static_cast<TH1D*>(pythia_spectra_uncut_reference->Clone());
  
  pythia_spectra_uncut_check->Reset();
  
  
  // Reader of tree variables
  TTreeReader* reader = new TTreeReader("pythia_tree", pythia_file);
  TTreeReaderValue<float>                 jetweight(*reader, "binweight");    // Weight of this jet due to pTHard binning
  TTreeReaderValue<int>                   mult(*reader, "ntotal");            // Total event track multiplicity
  TTreeReaderValue<int>                   njet(*reader, "njet");              // Number of full jets in an event
  TTreeReaderValue<std::vector<int>>      jet_mult(*reader, "jet_n");         // Array of full jet constituent multiplicties
  TTreeReaderValue<std::vector<float>>    jet_pT(*reader, "jet_pT");          // Array of full jet p_T
  TTreeReaderValue<std::vector<float>>    jet_y(*reader, "jet_y");            // Array of full jet rapidity
  TTreeReaderValue<std::vector<float>>    jet_eta(*reader, "jet_eta");        // Array of full jet pseudorapidity
  TTreeReaderValue<std::vector<float>>    jet_phi(*reader, "jet_phi");        // Array of full jet azimuth
  TTreeReaderValue<std::vector<float>>    jet_area(*reader, "jet_area");      // Array of full jet area
  TTreeReaderValue<int>                   cmult(*reader, "ncharge");          // Event charged track multiplicity
  TTreeReaderValue<int>                   ncjet(*reader, "nchjet");           // Number of charged jets in an event
  TTreeReaderValue<std::vector<int>>      cjet_mult(*reader, "chjet_n");      // Array of charged jet constituent multiplicties
  TTreeReaderValue<std::vector<float>>    cjet_pT(*reader, "chjet_pT");       // Array of charged jet p_T
  TTreeReaderValue<std::vector<float>>    cjet_y(*reader, "chjet_y");         // Array of charged jet rapidity
  TTreeReaderValue<std::vector<float>>    cjet_eta(*reader, "chjet_eta");     // Array of charged jet pseudorapidity
  TTreeReaderValue<std::vector<float>>    cjet_phi(*reader, "chjet_phi");     // Array of charged jet azimuth
  TTreeReaderValue<std::vector<float>>    cjet_area(*reader, "chjet_area");   // Array of charged jet area
  
  
  // Print some information about the final cuts to be used (should match the R_AA dataset of choice)
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  std::cout << "Beginning energy loss fitting with the following final event cuts:" << std::endl;
  std::cout << Form("PYTHIA pp sqrt_s_NN = %.2f TeV, pT_hadron >= %.2f GeV, |eta_hadron| < %.2f",
                    sqrt_s, pTmin_hadron, max_eta) << std::endl;
  std::cout << Form("FastJet3 R = %.1f %s Jets, |eta_jet| < %.2f, jet_A >= %.2f*pi*R^2",
                    jet_radius,algo_string,max_eta - jet_radius,min_area/(3.141592*jet_radius*jet_radius)) << std::endl;
  std::cout << Form("pT_jet >= %.2f, pT_charged_jet >= %.2f, pT_leading >= %.2f",
                    jetcut_minpT_full,jetcut_minpT_chrg, pTmin_jetcore) << std::endl;
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  
  // Gather jet data and perform the cuts
  std::vector<float> pT_unbinned;
  std::vector<float> weight_unbinned;
  if (!flag_use_charged_jets) { // Use full jets
    while (reader->Next()) for (int i = 0; i < jet_pT->size(); ++i) {
      pythia_spectra_uncut_check->Fill(jet_pT->at(i), *jetweight);
      
      if (jet_eta->at(i) > max_eta - jet_radius) continue;
      if (jet_area->at(i) < min_area) continue;
      if (jet_pT->at(i) < jetcut_minpT_full || jet_pT->at(i) < data_lowestpT) continue;
      // Will have jet core cut here...
      
      // Fill each ref hist separately since the binning of each dataset may be different...
      for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent)
        hist_baseline[i_cent]->Fill(jet_pT->at(i), *jetweight);
      pT_unbinned.push_back(jet_pT->at(i));
      weight_unbinned.push_back(*jetweight);
    }
  } else { // Use charged jets
    while (reader->Next()) {
      for (int i = 0; i < cjet_pT->size(); ++i) {
        pythia_spectra_uncut_check->Fill(cjet_pT->at(i), *jetweight);
        
        if (cjet_eta->at(i) > max_eta - jet_radius) continue;
        if (cjet_area->at(i) < min_area) continue;
        if (cjet_pT->at(i) < jetcut_minpT_chrg) continue;
        // Will have jet core cut here...
        
        // Fill each ref hist separately since the binning of each dataset may be different...
        for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent)
          hist_baseline[i_cent]->Fill(cjet_pT->at(i), *jetweight);
        pT_unbinned.push_back(cjet_pT->at(i));
        weight_unbinned.push_back(*jetweight);
      }
    }
  }
  
  //Spectra are already normalized, except by bin width
  pythia_spectra_uncut_check->Scale(1, "width");
  for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent)
    hist_baseline[i_cent]->Scale(1, "width");
  
  
  std::cout << "Cuts completed. Printing check against reference spectra..." << std::endl;
  gPad->SetLogy();
  pythia_spectra_uncut_reference->Draw("hist");
  pythia_spectra_uncut_check->SetLineColorAlpha(kRed, 0.5);
  pythia_spectra_uncut_check->SetLineStyle(7);
  pythia_spectra_uncut_check->SetLineWidth(2);
  pythia_spectra_uncut_check->Draw("hist same");
  for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent)
    hist_baseline[i_cent]->Draw("hist same");
  canvas_check->SaveAs(Form("%s/checks/check_spectra.pdf",outplot_directory));
  
  
  std::cout << "pT_unbinned size = " << pT_unbinned.size() << ", weight_unbinned.size() = " << weight_unbinned.size() << std::endl;
  
  //---------------------------------------------------------- Perform subtraction
  
  const int ntotal_jets = pT_unbinned.size();
  
  // Reproduce the best fit...
  canvas_check->Clear();
  std::vector<std::vector<TPad*>> pads = divideFlush(gPad, 1, 2);
  
  gPad->SetLogy(0);
  for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent) {
    TH1F* v2_flanking_hist[2];
    TH1F* v2_result_hist[4];
    
    // Perform the shift using unbinned jet data
    hist_shifting[i_cent]->Reset();
    for (int i_jet = 0; i_jet < ntotal_jets; ++i_jet)
      hist_shifting[i_cent]->Fill(pT_unbinned.at(i_jet) - pT_shift_best[i_cent], weight_unbinned.at(i_jet));
    hist_shifting[i_cent]->Scale(1, "width");
    
    pads.at(0).at(0)->cd();
    hist_v2_datreference[i_cent]->GetYaxis()->SetRangeUser(0,0.2);
    hist_v2_datreference[i_cent]->Draw("hist p e0");
    for (int i_algo = 0; i_algo < 4; ++i_algo) {
      pads.at(1).at(0)->cd();
      gPad->SetLogy();
      hist_shifting[i_cent]->Draw("hist");
      // Compute left and right flanking shift hists
      for (int i_flank = 0; i_flank <= 1; ++i_flank) {
        v2_flanking_hist[i_flank] = static_cast<TH1F*>(hist_v2_datreference[i_cent]->Clone());
        v2_flanking_hist[i_flank]->Reset();
        
        
        float flank_delta_pT = - pT_shift_best[i_cent] + (2*i_flank-1) * c2_glauber[i_cent][i_algo]/(0.5*c0_glauber[i_cent][i_algo]);
        for (int i_jet = 0; i_jet < ntotal_jets; ++i_jet)
          v2_flanking_hist[i_flank]->Fill(pT_unbinned.at(i_jet) + flank_delta_pT, weight_unbinned.at(i_jet));
        v2_flanking_hist[i_flank]->Scale(1, "width");
        
        pads.at(1).at(0)->cd();
        v2_flanking_hist[i_flank]->SetLineColor(kGray+i_flank);
        v2_flanking_hist[i_flank]->Draw("hist same");
      }
      
      
      // Compute difference for v2
      v2_result_hist[i_algo] = static_cast<TH1F*>(hist_v2_datreference[i_cent]->Clone());
      v2_result_hist[i_algo]->SetLineColor(plot_color[i_algo]);
      v2_result_hist[i_algo]->Reset();
      
      for (int i_bin = 1; i_bin <= v2_result_hist[i_algo]->GetXaxis()->GetNbins(); ++i_bin) {
        float v2val = ( ( v2_flanking_hist[1]->GetBinContent(i_bin) - v2_flanking_hist[0]->GetBinContent(i_bin) )
                       / (hist_shifting[i_cent]->GetBinContent(i_bin)) );
        std::cout << "i_cent = " <<i_cent << ", i_algo = " << i_algo << " v2val " << v2val << std::endl;
        v2_result_hist[i_algo]->SetBinContent(i_bin,v2val);
      }
      
      pads.at(0).at(0)->cd();
      v2_result_hist[i_algo]->Draw("hist same");
    }
    
    // Save the plot
    canvas_check->SaveAs(Form("%s/v2_results_cent%s.pdf",outplot_directory,v2_centbin[i_cent]));
  }
  
//  float delta_pT = 0.01; // Energy loss resolution
//  float max_delta_pT = 18;
//  const int ntotal_scan = (int)(max_delta_pT/delta_pT);
//  TH1F* chi2_hist[ndata_centbin];
//  float fitmetric_best[ndata_centbin];
//  float pT_shift_best[ndata_centbin];
//  
//  for (int i_cent = 0; i_cent < ndata_centbin; ++i_cent) {
//    fitmetric_best[i_cent] = 1e20; // large value initially
//    
//    std::cout << "Beginning shift finder for centbin" << v2_centbin[i_cent] << "%" << std::endl;
//    
//    chi2_hist[i_cent] = new TH1F(Form("Fit_hist_cent%s",v2_centbin[i_cent]),
//                                      ";#Delta p_{T};MSE", ntotal_scan, delta_pT, max_delta_pT);
//    
//    
//    // To start, a simple scan over the delta pT space
//    // (should replace with a faster method later...)
//    for (int i_delta = 1; i_delta < ntotal_scan; ++i_delta) {
//      
//      // Perform the shift using unbinned jet data
//      hist_shifting[i_cent]->Reset();
//      for (int i_jet = 0; i_jet < ntotal_jets; ++i_jet)
//        hist_shifting[i_cent]->Fill(pT_unbinned.at(i_jet) - i_delta*delta_pT, weight_unbinned.at(i_jet));
//      hist_shifting[i_cent]->Scale(1, "width");
//      
//      // Compute the R_AA against the baseline spectrum
//      ratio_shifting[i_cent]->Reset();
//      for (int i_bin = 1; i_bin <= hist_shifting[i_cent]->GetXaxis()->GetNbins(); ++i_bin) {
//        if (hist_baseline[i_cent]->GetBinContent(i_bin) == 0) continue;
//        ratio_shifting[i_cent]->SetBinContent(i_bin, hist_shifting[i_cent]->GetBinContent(i_bin) / hist_baseline[i_cent]->GetBinContent(i_bin));
//      }
//      
//      
//      // Compare against data RAA using MSE fit metric
//      // Currently separated in case any error analysis from data stat/syst comes into play...
//      float MSE = 0;
//      for (int i_bin = 1; i_bin <= hist_shifting[i_cent]->GetXaxis()->GetNbins(); ++i_bin) {
//        if (ratio_shifting[i_cent]->GetBinContent(i_bin) == 0) continue;
//        if (hist_v2_datreference[i_cent]->GetBinContent(i_bin) == 0) continue;
//        
//        float residual = TMath::Log(hist_v2_datreference[i_cent]->GetBinContent(i_bin)
//                                    / ratio_shifting[i_cent]->GetBinContent(i_bin));
//        MSE += residual*residual;
//      }
//      
//      chi2_hist[i_cent]->Fill(i_delta*delta_pT, MSE);
//      if (MSE < fitmetric_best[i_cent]) {
//        fitmetric_best[i_cent] = MSE;
//        pT_shift_best[i_cent] = i_delta*delta_pT;
//      }
//    }// End of pT Shift loop
//    
//    // Found best fit! Say some things about it...
//    
//    std::cout << "best fit delta_pT = " << pT_shift_best[i_cent] << ", found with MSE = " << fitmetric_best[i_cent] << std::endl;
//    
//    // Draw fit metric over the parameter space...
//    gPad->SetLogy();
//    chi2_hist[i_cent]->Draw("hist");
//    canvas_check->SaveAs(Form("%s/delta_pT_fitmetric_cent%s.pdf",outplot_directory,v2_centbin[i_cent]));
//    
//    // Reproduce the best fit...
//    
//    // Perform the shift using unbinned jet data
//    hist_shifting[i_cent]->Reset();
//    for (int i_jet = 0; i_jet < ntotal_jets; ++i_jet)
//      hist_shifting[i_cent]->Fill(pT_unbinned.at(i_jet) - pT_shift_best[i_cent], weight_unbinned.at(i_jet));
//    hist_shifting[i_cent]->Scale(1, "width");
//    
//    // Compute the R_AA against the baseline spectrum
//    ratio_shifting[i_cent]->Reset();
//    for (int i_bin = 1; i_bin <= hist_shifting[i_cent]->GetXaxis()->GetNbins(); ++i_bin) {
//      if (hist_baseline[i_cent]->GetBinContent(i_bin) == 0) continue;
//      ratio_shifting[i_cent]->SetBinContent(i_bin, hist_shifting[i_cent]->GetBinContent(i_bin) / hist_baseline[i_cent]->GetBinContent(i_bin));
//    }
//    
//    gPad->SetLogy(0);
//    hist_v2_datreference[i_cent]->GetYaxis()->SetRangeUser(0, 0.8);
//    hist_v2_datreference[i_cent]->Draw("hist");
//    ratio_shifting[i_cent]->SetLineColor(kViolet);
//    ratio_shifting[i_cent]->Draw("hist same");
//    canvas_check->SaveAs(Form("%s/delta_pT_fitresult_cent%s.pdf",outplot_directory, v2_centbin[i_cent]));
//    
//  }// End of centbin loop
  
  
  
  return;
}
