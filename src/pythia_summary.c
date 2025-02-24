// ROOT macro. Makes summary sheets to display data from event generation.
// Assumes relevant file with input cuts has already been generated
// See PYTHIA/eventgen_jets.cc for details on output data format

#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"
#include "gen_pythia/config_local.h"
#include "gen_pythia/ptbin.h"

// ---NOTES AND TODO---
// R 0.4 standard and use anti-kt, 10Gev min is what Andrew uses typically in data,
// could use 5GeV min in pythia (for charged fine, for full do 3/2 scale ->7.5 GeV), but this could cause gluon jet contamination
// TODO Maybe move key settings to config file, and make bash script for full gen process??
// Plots to add - neutral energy fraction!


// v2 for high pt tracks
// write the part about areas (appendix)

// Color and aesthetic settings for root plots
Int_t palette_full = kCherry;
Int_t palette_ch = kLake;
Int_t markerStyle_full = 20;
Int_t markerStyle_ch = 33;
Float_t markerSize_full = 0.75;
Float_t markerSize_ch = 1.0;


TPad* pad_main[2];
TPad* pads_left[5];



void pythia_summary(int this_thread = 0) {
  // Attempt to open file with specified parameters
  float collision_energy_TeV = (float) sqrt_s / 1000.;
  char gen_info_string[100];
  snprintf(gen_info_string, 100, "%.2fTeV_r%.1fjet", collision_energy_TeV, jetRadius);
  TFile *fin = new TFile(Form("gen_pythia/out_data/pythpp_%s.root", gen_info_string),"read");
  if (fin->IsZombie()) {
    std::cout << "File with given specifications does not exist." << std::endl;
    std::cout << "Check specfied parameters, or run PYTHIA generator and try again." << std::endl;
    return;
  }fin->ls();
  
  // which thread to run? (should be all threads, or one at a time?)
  const int total_threads = getNbins(sqrt_s);
  if (this_thread >= total_threads) {
    std::cerr << "Error in <eventgen_jets.cc>: requested pTHatBin " << this_thread << " out of bin range " << total_threads << std::endl;
    return;
  }
  pTBinSettings settings = getBinSettings(sqrt_s, this_thread);
  const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
  const float cross_section = settings.cross_section;
  const int nevent = settings.nevent/1000;
  
  
  // Read histogram settings and generator data from file
  TVectorD settings_pT = *(TVectorD*) fin->Get("settings_pT");
  int nbins_pT_trk = (int) settings_pT[0];
  double pTmax_jet = settings_pT[2];
  double pTmax_track = settings_pT[4];
  TVectorD settings_phi_y_eta = *(TVectorD*) fin->Get("settings_phi_y_eta");
  int nbins_phi = (int) settings_phi_y_eta[0];
  int nbins_rap = (int) settings_phi_y_eta[1];
  double max_y = settings_phi_y_eta[2];
  double max_eta = settings_phi_y_eta[3];
  TVectorD settings_PYTHIA  = *(TVectorD*) fin->Get("settings_PYTHIA");
  TVectorD gendata_Ntotal = *(TVectorD*) fin->Get("gendata_Ntotal");
  int nEvent = (int) gendata_Ntotal[0];
  int grandtotal_n[6];
  for (int i = 1; i < 7; ++i) grandtotal_n[i-1] = gendata_Ntotal[i];
  
  // Set up tree readers for extracting jet information
  char tree_name[50];
  snprintf(tree_name, 50, "pythia_tree_pTHard%.f-%.f", pt_binedge[0], pt_binedge[1]);
  TTreeReader *reader = new TTreeReader(tree_name, fin);
  if (!reader) {std::cout << "Tree not found!" << std::endl; return;}
  TTreeReaderValue<Int_t>                     ntotal(*reader, "ntotal");        // Total event track multiplicity
  TTreeReaderValue<Int_t>                     ncharge(*reader, "ncharge");      // Event charged track multiplicity
  TTreeReaderValue<Float_t>                   deltaE(*reader, "deltaE");        // Deviation from energy conservation
  TTreeReaderValue<Int_t>                     njet(*reader, "njet");            // Number of full jets in an event
  TTreeReaderValue<Int_t>                     ncjet(*reader, "ncjet");          // Number of charged jets in an event
  TTreeReaderValue<std::vector<Int_t>>        jet_n(*reader, "jet_n");          // Array of full jet constituent multiplicties
  TTreeReaderValue<std::vector<Int_t>>        cjet_n(*reader, "cjet_n");        // Array of charged jet constituent multiplicities
  TTreeReaderValue<std::vector<Float_t>>      jet_pT(*reader, "jet_pT");        // Array of full jet p_T
  TTreeReaderValue<std::vector<Float_t>>      cjet_pT(*reader, "cjet_pT");      // Array of charged jet p_T
  TTreeReaderValue<std::vector<Float_t>>      jet_y(*reader, "jet_y");          // Array of full jet rapidity
  TTreeReaderValue<std::vector<Float_t>>      cjet_y(*reader, "cjet_y");        // Array of charged jet rapidity
  TTreeReaderValue<std::vector<Float_t>>      jet_eta(*reader, "jet_eta");      // Array of full jet pseudorapidity
  TTreeReaderValue<std::vector<Float_t>>      cjet_eta(*reader, "cjet_eta");    // Array of charged jet pseudorapidity
  TTreeReaderValue<std::vector<Float_t>>      jet_phi(*reader, "jet_phi");      // Array of full jet azimuth
  TTreeReaderValue<std::vector<Float_t>>      cjet_phi(*reader, "cjet_phi");    // Array of charged jet azimuth
  TTreeReaderValue<std::vector<Float_t>>      jet_area(*reader, "jet_area");    // Array of full jet areas
  TTreeReaderValue<std::vector<Float_t>>      cjet_area(*reader, "cjet_area");  // Array of charged jet areas
  
  // Find histograms from file
  // Stored in arrays for later convenience. Index listings:
  //    0-1: jets, full and charged
  //    2-3: constituent level particles, full and charged jet constituents
  //    4-5: track level particles, full and charged tracks
  TH1D *hist_pT_array[6];
  TH2D *hist_y_phi_array[6];
  TH2D *hist_eta_phi_array[6];
  const char specifications[6][7] = {"jet", "cjet", "const", "cconst", "track", "ctrack"};
  const char descriptions[6][15] = {"jet", "Jet", "cst", "Constituent", "trk", "Track"};
  for (int i = 2; i < 6; ++i) {
    hist_pT_array[i] = (TH1D*) fin->Get(Form("hist_pT_%s", specifications[i]));
    hist_y_phi_array[i] = (TH2D*) fin->Get(Form("hist_y_phi_%s", specifications[i]));
    hist_eta_phi_array[i] = (TH2D*) fin->Get(Form("hist_eta_phi_%s", specifications[i]));
  }
  
  // Make similar histograms to the above for jets using tree data
  hist_pT_array[0] = new TH1D("hist_pT_jet",
                               ";#it{p}_{T} [GeV];#frac{1}{N^{event}} #frac{dN^{(ch) jets}}{d#it{p}_{T}}",
                               nbins_pT_jet,pt_binedge[0],pTmax_jet);
  hist_pT_array[1] = new TH1D("hist_pT_cjet",
                               ";#it{p}_{T} [GeV];#frac{1}{N^{event}} #frac{dN^{(ch) jets}}{d#it{p}_{T}}",
                               nbins_pT_jet,pt_binedge[0],pTmax_jet);
  hist_y_phi_array[0] = new TH2D("hist_y_phi_jet",
                               ";y;#phi;#frac{1}{N^{event}} #frac{d^{2}N^{jets}}{dyd#phi}",
                               nbins_rap,-max_y,max_y,nbins_phi,-TMath::Pi(),TMath::Pi());
  hist_y_phi_array[1] = new TH2D("hist_y_phi_cjet",
                               ";y;#phi;#frac{1}{N^{event}} #frac{d^{2}N^{ch jets}}{dyd#phi}",
                               nbins_rap,-max_y,max_y,nbins_phi,-TMath::Pi(),TMath::Pi());
  hist_eta_phi_array[0] = new TH2D("hist_eta_phi_jet",
                               ";y;#phi;#frac{1}{N^{event}} #frac{d^{2}N^{jets}}{d#etad#phi}",
                               nbins_rap,-max_eta,max_eta,nbins_phi,-TMath::Pi(),TMath::Pi());
  hist_eta_phi_array[1] = new TH2D("hist_eta_phi_cjet",
                               ";y;#phi;#frac{1}{N^{event}} #frac{d^{2}N^{ch jets}}{d#etad#phi}",
                               nbins_rap,-max_eta,max_eta,nbins_phi,-TMath::Pi(),TMath::Pi());
  // Histograms for y, eta, phi projections
  TH1D *proj_y_array[6];
  TH1D *proj_eta_array[6];
  TH1D *proj_phi_fromy_array[6];
  TH1D *proj_phi_frometa_array[6];
  for (int i = 0; i < 6; ++i) {
    proj_y_array[i] = new TH1D("","",nbins_rap, -max_y, max_y);
    proj_eta_array[i] = new TH1D("","",nbins_rap, -max_eta, max_eta);
    proj_phi_fromy_array[i] = new TH1D("","",nbins_phi, -TMath::Pi(), TMath::Pi());
    proj_phi_frometa_array[i] = new TH1D("","",nbins_phi, -TMath::Pi(), TMath::Pi());
  }
  // Histograms for general event information
  TH1D *hist_n_array[6];
  hist_n_array[0] = new TH1D("hist_ntotal_jets",
                               ";N^{(ch) jets};#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) jets}}",6,0,6);
  hist_n_array[1] = new TH1D("hist_ncharge_jets",
                               ";N^{(ch) jets};#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) jets}}",6,0,6);
  hist_n_array[2] = new TH1D("hist_ntotal_const",
                               ";N^{(ch) constituents};#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) constituents}}",100,0,100);
  hist_n_array[3] = new TH1D("hist_ncharge_const",
                               ";N^{(ch) constituents};#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) constituents}}",100,0,100);
  hist_n_array[4] = new TH1D("hist_ntotal_tracks",
                               ";N^{(ch) tracks};#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) tracks}}",100,0,600);
  hist_n_array[5] = new TH1D("hist_ncharge_tracks",
                               ";N^{(ch) tracks};#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) tracks}}",100,0,600);
  TH1D *hist_area_jet = new TH1D("hist_area_jet",
                                 ";Jet Area;#frac{1}{N^{event}} #frac{dN^{event}}{dA^{(ch) jet}}"
                                 ,50,0.5*TMath::Pi()*jetRadius*jetRadius,1.5*TMath::Pi()*jetRadius*jetRadius);
  TH1D *hist_area_cjet = new TH1D("hist_area_cjet",
                                 ";Jet Area;#frac{1}{N^{(ch) jet}} #frac{dN^{(ch) jet}}{dA^{(ch) %s}}"
                                  ,50,0.5*TMath::Pi()*jetRadius*jetRadius,1.5*TMath::Pi()*jetRadius*jetRadius);
  float deltaE_xbins[11] = {1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-5};
  TH1D *hist_deltaE = new TH1D("hist_deltaE",
                               ";#DeltaE [GeV];#frac{1}{N^{event}} #frac{dN^{event}}{d#DeltaE}",10,deltaE_xbins);
  
  // Fill jet and general event histograms from tree data
  while (reader->Next()) {
    // DeltaE hist
    hist_deltaE->Fill(*deltaE);
    
    // multiplicity hists
    hist_n_array[0]->Fill(*njet);
    hist_n_array[1]->Fill(*ncjet);
    hist_n_array[4]->Fill(*ntotal);
    hist_n_array[5]->Fill(*ncharge);
    int total_constituents = 0;
    int total_ch_constituents = 0;
    for (int i = 0; i < *njet; ++i)  total_constituents += (*jet_n).at(i);
    for (int i = 0; i < *ncjet; ++i) total_ch_constituents += (*cjet_n).at(i);
    if (total_constituents > 0) hist_n_array[2]->Fill(total_constituents);
    if (total_ch_constituents > 0) hist_n_array[3]->Fill(total_ch_constituents);
      
    // Jet pT hists
    for (int i = 0; i < (*jet_pT).size(); ++i) hist_pT_array[0]->Fill((*jet_pT).at(i));
    for (int i = 0; i < (*cjet_pT).size(); ++i) hist_pT_array[1]->Fill((*cjet_pT).at(i));
    
    // Jet area hists
    for (float f:(*jet_area)) if (f > 0) hist_area_jet->Fill(f);
    for (float f:(*cjet_area)) if (f > 0) hist_area_cjet->Fill(f);
    
    // Jet y/eta-phi hists
    for (int i = 0; i < *njet; ++i) {
      hist_y_phi_array[0]->Fill((*jet_y).at(i), (*jet_phi).at(i));
      hist_eta_phi_array[0]->Fill((*jet_eta).at(i), (*jet_phi).at(i));
    }for (int i = 0; i < *ncjet; ++i) {
      hist_y_phi_array[1]->Fill((*cjet_y).at(i), (*cjet_phi).at(i));
      hist_eta_phi_array[1]->Fill((*cjet_eta).at(i), (*cjet_phi).at(i));
    }
  }
  
  // Renormalize histograms
  hist_deltaE->Scale(1./nEvent, "width");
  hist_area_jet->Scale(1./grandtotal_n[0], "width");
  hist_area_cjet->Scale(1./grandtotal_n[1], "width");
  for (int i = 0; i < 6; ++i) {
    if (i == 2) { // Full Constituents
      hist_n_array[i]->Scale(1./grandtotal_n[0], "width");
    }else if (i == 3) {
      hist_n_array[i]->Scale(1./grandtotal_n[1], "width");
    }else {
      hist_n_array[i]->Scale(1./nEvent, "width");
    }hist_pT_array[i]->Scale(1./nEvent, "width");
    hist_y_phi_array[i]->Scale(1./(nEvent *
                                hist_y_phi_array[i]->GetXaxis()->GetBinWidth(1) *
                                hist_y_phi_array[i]->GetYaxis()->GetBinWidth(1)));
    hist_eta_phi_array[i]->Scale(1./(nEvent *
                                hist_eta_phi_array[i]->GetXaxis()->GetBinWidth(1) *
                                hist_eta_phi_array[i]->GetYaxis()->GetBinWidth(1)));
  }
  
  // Histograms for y, eta, phi projections
  // Better to fill jet proj hists here, rather than in the above loop;
  // this allows an extra check that projections are normalized correctly.
  for (int i = 0; i < 6; ++i) {
    // Loop over relevant 2D hists and fill projections
    double c_phi, c_eta, c_y;
    for (int iPhi = 1; iPhi <= nbins_phi; ++iPhi) {
      c_phi = proj_phi_fromy_array[0]->GetBinCenter(iPhi);
      for (int iRap = 1; iRap <= nbins_rap; ++iRap) {
        c_y = proj_y_array[0]->GetBinCenter(iRap);
        c_eta = proj_eta_array[0]->GetBinCenter(iRap);
        
        proj_y_array[i]->Fill(c_y, hist_y_phi_array[i]->GetBinContent(
                              hist_y_phi_array[i]->GetBin(iRap, iPhi)));
        proj_eta_array[i]->Fill(c_eta, hist_eta_phi_array[i]->GetBinContent(
                              hist_eta_phi_array[i]->GetBin(iRap, iPhi)));
        proj_phi_fromy_array[i]->Fill(c_phi, hist_y_phi_array[i]->GetBinContent(
                              hist_y_phi_array[i]->GetBin(iRap, iPhi)));
        proj_phi_frometa_array[i]->Fill(c_phi, hist_eta_phi_array[i]->GetBinContent(
                              hist_y_phi_array[i]->GetBin(iRap, iPhi)));
      }
    }
  }
  
  // ROOT canvas setup -- one for individual plots and one for summary sheets
  gStyle->SetOptStat(0);
  TCanvas *composite = new TCanvas();
  composite->SetCanvasSize(800,500);
  composite->SetMargin(0, 0, 0, 0);
//  TPad* pad_main[2];
  pad_main[0] = buildPad("pads_main0", 0, 0, 0.385, 0.9, 0, 0, 0, 0);
  pad_main[1] = buildPad("pads_main1", 0.385, 0, 1, 0.9, 0, 0, 0, 0);
//  TPad* pads_left[5]; 
  pad_main[0]->cd();
  pads_left[0] = buildPad("pads_pT", 0, 2./3., 1, 1, 0.05, 0.05, 0.11, 0);
  pads_left[1] = buildPad("pads_misc1", 0, 1./3., 0.5, 2./3., 0.12, 0, 0.1, 0.02);
  pads_left[2] = buildPad("pads_misc2", 0.5, 1./3., 1, 2./3., 0.07, 0.1, 0.1, 0.02);
  pads_left[3] = buildPad("pads_misc3", 0, 0, 0.5, 1./3., 0.12, 0, 0.1, 0.02);
  pads_left[4] = buildPad("pads_misc4", 0.5, 0, 1, 1./3., 0.07, 0.1, 0.1, 0.02);
  TPad* pads_right[10]; pad_main[1]->cd();
  pads_right[0] = buildPad("pads_rap_f", 0, 5.5/9, 0.35, 1,       0, 0, 0, 0, false);
  pads_right[1] = buildPad("pads_eta_f", 0.35, 5.5/9, 0.7, 1,     0, 0, 0, 0, false);
  pads_right[2] = buildPad("pads_rap_c", 0, 2./9, 0.35, 5.5/9,    0, 0, 0, 0, false);
  pads_right[3] = buildPad("pads_eta_c", 0.35, 2./9, 0.7, 5.5/9,  0, 0, 0, 0, false);
  pads_right[4] = buildPad("pads_zscale", 0.9, 2./9*(0.95), 1, 1,        0, 0, 0, 0, false);
  pads_right[5] = buildPad("pads_phi_f", 0.7, 5.5/9, 0.9, 1,      0.01, 0.15, 0, 0);
  pads_right[6] = buildPad("pads_phi_c", 0.7, 2./9*(0.9), 0.9, 5.5/9,   0.01, 0.15, 0.058, 0);
  pads_right[7] = buildPad("pads_yproj", 0, 0, 0.35, 2./9,        0, 0, 0.15, 0.01);
  pads_right[8] = buildPad("pads_eproj", 0.35, 0, 0.7, 2./9,      0, 0, 0.15, 0.01);
  pads_right[9] = buildPad("pads_sinfo", 0.7, 0, 1, 2./9*(0.9),         0, 0, 0, 0);
  TCanvas *indiv = new TCanvas();
  indiv->SetCanvasSize(700,500);
  
  // Color and legend setup
  gStyle->SetPalette(kCherry);
  Int_t color_bg = TColor::GetPalette().At(0);
  gStyle->SetPalette(palette_full);
  Int_t color_full = TColor::GetPalette().At(TColor::GetPalette().GetSize()*1/4);
  gStyle->SetPalette(palette_ch);
  Int_t color_ch = TColor::GetPalette().At(TColor::GetPalette().GetSize()*1/4);
  TLegend *leg = new TLegend(0.75, 0.8, 0.95, 0.95);
  leg->SetLineWidth(0);
  
  //==========================================================================Format and print individual plots
  indiv->cd();
  double rangeY;
  for (int i = 0; i < 6; ++i) {
    if (i % 2 == 0) {
      indiv->SetMargin(0.12, 0.04, 0.08, 0.04);
      indiv->SetFrameFillColor(kWhite);
      // Draw multiplicity hist
      hist_n_array[i]->SetMarkerStyle(markerStyle_full);
      hist_n_array[i]->SetMarkerColor(color_full);
      hist_n_array[i]->SetLineColor(color_full);
      hist_n_array[i]->SetFillColorAlpha(color_full, 0.3);
      hist_n_array[i+1]->SetMarkerStyle(markerStyle_ch);
      hist_n_array[i+1]->SetMarkerColor(color_ch);
      hist_n_array[i+1]->SetLineColor(color_ch);
      hist_n_array[i+1]->SetFillColorAlpha(color_ch, 0.3);
      hist_n_array[i]->GetYaxis()->SetRangeUser(0, 1.1*getMaxFromHists(hist_n_array[i], hist_n_array[i+1]));
      hist_n_array[i]->Draw("hist");
      hist_n_array[i+1]->Draw("hist same");
      leg->Clear();
      if (i/2 != 1) {
        leg->AddEntry(hist_n_array[i], Form("Full %ss", specifications[i]), "f");
        leg->AddEntry(hist_n_array[i+1], Form("Charged %ss", specifications[i]), "f");
      } else {
        leg->AddEntry(hist_n_array[i], "Full constituents", "f");
        leg->AddEntry(hist_n_array[i+1], "Charged constituents", "f");
      }
      leg->Draw();
      indiv->SaveAs(Form("gen_pythia/out_plots/lastgen/multiplicity_%s.pdf", specifications[i]));
      
      // Draw pT hist
      gPad->SetLogy();
      hist_pT_array[i]->SetMarkerStyle(markerStyle_full);
      hist_pT_array[i]->SetMarkerColor(color_full);
      hist_pT_array[i]->SetLineColor(color_full);
      hist_pT_array[i]->SetFillColorAlpha(color_full, 0.3);
      hist_pT_array[i+1]->SetMarkerStyle(markerStyle_ch);
      hist_pT_array[i+1]->SetMarkerColor(color_ch);
      hist_pT_array[i+1]->SetLineColor(color_ch);
      hist_pT_array[i+1]->SetFillColorAlpha(color_ch, 0.3);
      hist_pT_array[i]->Draw("hist");
      hist_pT_array[i+1]->Draw("hist same");
      leg->Draw();
      indiv->SaveAs(Form("gen_pythia/out_plots/lastgen/pT_%s.pdf", specifications[i]));
      gStyle->SetPalette(palette_full);
    } else gStyle->SetPalette(palette_ch);
    // Draw rapidity-phi hist
    gPad->SetLogy(0);
    indiv->SetMargin(0.1, 0.19, 0.1, 0.04);
    indiv->SetFrameFillColor(color_bg);
    hist_y_phi_array[i]->GetXaxis()->SetAxisColor(kWhite);
    hist_y_phi_array[i]->GetYaxis()->SetAxisColor(kWhite);
    hist_y_phi_array[i]->GetZaxis()->SetTitleOffset(1. + 0.2*((7-i)/2));
    hist_y_phi_array[i]->Draw("colz");
    indiv->SaveAs(Form("gen_pythia/out_plots/lastgen/rapidity-phi_%s.pdf", specifications[i]));
    hist_eta_phi_array[i]->GetXaxis()->SetAxisColor(kWhite);
    hist_eta_phi_array[i]->GetYaxis()->SetAxisColor(kWhite);
    hist_eta_phi_array[i]->GetZaxis()->SetTitleOffset(1. + 0.2*((7-i)/2));
    hist_eta_phi_array[i]->Draw("colz");
    indiv->SaveAs(Form("gen_pythia/out_plots/lastgen/eta-phi_%s.pdf", specifications[i]));
  }
  gPad->SetLogy();
  indiv->SetMargin(0.12, 0.04, 0.08, 0.04);
  indiv->SetFrameFillColor(kWhite);
  hist_area_jet->SetMarkerColor(color_full);
  hist_area_jet->SetLineColor(color_full);
  hist_area_jet->SetMarkerSize(markerSize_full);
  hist_area_jet->SetMarkerStyle(markerStyle_full);
  hist_area_jet->Draw("p e0");
  hist_area_cjet->SetMarkerColor(color_ch);
  hist_area_cjet->SetLineColor(color_ch);
  hist_area_cjet->SetMarkerSize(markerSize_ch);
  hist_area_cjet->SetMarkerStyle(markerStyle_ch);
  hist_area_cjet->Draw("p e0 same");
  leg->Clear();
  leg->AddEntry(hist_area_jet, "Full jet area", "lep");
  leg->AddEntry(hist_area_cjet, "Charged jet area", "lep");
  leg->Draw();
  indiv->SaveAs("gen_pythia/out_plots/lastgen/jet-area.pdf");
  
  gPad->SetLogx();
  hist_deltaE->SetMarkerColor(kRed+2);
  hist_deltaE->SetLineColor(kRed+2);
  hist_deltaE->SetMarkerStyle(markerStyle_full);
  hist_deltaE->Draw("p e0");
  indiv->SaveAs("gen_pythia/out_plots/lastgen/deltaE.pdf");
  
  //==========================================================================Generate composite plot/summary sheet
  TString plotfile = Form("gen_pythia/out_plots/%s_summarysheet.pdf", gen_info_string);
  composite->Print(plotfile + "[");
  leg = new TLegend(0.83, 0.7, 0.945, 0.835);
  leg->SetLineWidth(0);
  leg->SetTextAlign(32);
  leg->SetFillStyle(4000);
  TLegend* tinyLeg = new TLegend(0.68, 0.7, 0.895, 0.835);
  tinyLeg->SetTextAlign(32);
  tinyLeg->SetLineWidth(0);
  tinyLeg->SetFillStyle(4000);
  composite->cd();
  drawText(Form("#it{PYTHIA} #bf{pp}, #sqrt{s} = %.f GeV, N^{event} = %i", sqrt_s, nEvent), 1-0.05/3, 0.93, true, kBlack, 0.035);
  for (int i = 0; i < 3; ++i) {
    // Text headers using key generator information
    composite->cd();
    TLatex* headerTex = drawText(Form("Event Generation Summary Sheet (%ss)", descriptions[2*i+1]), 0.05/3, 0.93, false, kBlack, 0.043);
    // Legend setup for general hists
    leg->Clear();
    leg->AddEntry(hist_n_array[2*i], Form("Full %s", descriptions[2*i]), "f");
    leg->AddEntry(hist_n_array[2*i+1], Form("Ch. %s", descriptions[2*i]), "f");
    tinyLeg->Clear();
    tinyLeg->AddEntry(hist_n_array[2*i], Form("Full %s", descriptions[2*i]), "f");
    tinyLeg->AddEntry(hist_n_array[2*i+1], Form("Ch. %s", descriptions[2*i]), "f");
    
    // Draw pT hist
    pads_left[0]->cd();
    pads_left[0]->SetLogy();
    hist_pT_array[2*i]->GetXaxis()->SetTitleSize(0.05);
    hist_pT_array[2*i]->SetTitle(Form(";(Charged) %s #it{p}_{T} [GeV];",descriptions[2*i+1]));
    hist_pT_array[2*i]->Draw("hist");
    hist_pT_array[2*i+1]->Draw("hist same");
    drawText(Form("#frac{1}{N^{event}} #frac{dN^{(ch) %s}}{d#it{p}_{T}}", descriptions[2*i]), 0.94, 0.9, true, kBlack, 0.05);
    leg->Draw();
    
    // Draw N hist
    // Only constituent hist is different -- plotted per jet rather than per event
    pads_left[2]->cd();
    if (i==1) hist_n_array[2*i]->SetTitle(Form(";(Charged) %s Multiplicity, per-jet (N^{(ch) %s});",
                                                      descriptions[2*i+1], descriptions[2*i]));
    else hist_n_array[2*i]->SetTitle(Form(";(Charged) %s, Multiplicity, per-event (N^{(ch) %s});",
                                                      descriptions[2*i+1], descriptions[2*i]));
    hist_n_array[2*i]->Draw("hist");
    hist_n_array[2*i+1]->Draw("hist same");
    if (i==1) drawText(Form("#frac{1}{N^{(ch) jet}} #frac{dN^{(ch) jet)}}{dN^{(ch) %s}}", descriptions[2*i]), 0.89, 0.89, true, kBlack, 0.05);
    else drawText(Form("#frac{1}{N^{event}} #frac{dN^{event}}{dN^{(ch) %s}}", descriptions[2*i]), 0.89, 0.89, true, kBlack, 0.05);
    tinyLeg->Draw();
    
    // Draw extra left hist
    // Some placeholders currently, need to think more on what plots to add
    if (i == 0) { // jets
      pads_left[3]->cd();
      pads_left[3]->SetLogy();
      hist_area_jet->SetTitle(Form(";Reconstructed Jet Area in y-#phi Plane (%s);", algo_string_short));
      hist_area_jet->SetMarkerSize(markerSize_full / 2);
      hist_area_cjet->SetMarkerSize(markerSize_ch / 2);
      hist_area_jet->Draw("p e0");
      hist_area_cjet->Draw("p e0 same");
      TLegend* tinyLeg2 = new TLegend(0.78, 0.7, 0.995, 0.835);
      tinyLeg2->SetTextAlign(32);
      tinyLeg2->SetLineWidth(0);
      tinyLeg2->SetFillStyle(4000);
      tinyLeg2->AddEntry(hist_area_jet, Form("Full %s", descriptions[2*i]), "lep");
      tinyLeg2->AddEntry(hist_area_cjet, Form("Ch. %s", descriptions[2*i]), "lep");
      tinyLeg2->Draw();
      TLine* centerline = new TLine();
      centerline->SetLineStyle(7);
      centerline->SetLineColor(kGray+1);
      centerline->DrawLine(TMath::Pi()*jetRadius*jetRadius,0,
                           TMath::Pi()*jetRadius*jetRadius,10e3);
      drawText("#piR_{jet}^{2}", 0.53, 0.88, false, kGray+1, 0.04)->SetTextAngle(90);
      drawText(Form("#frac{1}{N^{(ch) jet}} #frac{dN^{(ch) jet}}{dA^{(ch) %s}}", descriptions[2*i]), 0.99, 0.89, true, kBlack, 0.05);
      
      pads_left[4]->Clear();
    } else if (i == 1) { // constituents
      pads_left[3]->Clear();
      
      pads_left[4]->Clear();
    } else { // tracks
      pads_left[3]->Clear();
      
      pads_left[4]->cd();
      pads_left[4]->SetLogy();
      pads_left[4]->SetLogx();
      hist_deltaE->GetYaxis()->SetTitle("");
      hist_deltaE->GetXaxis()->SetTitleOffset(1.1);
      hist_deltaE->GetYaxis()->SetNdivisions(4, 3);
      hist_deltaE->GetXaxis()->SetTitle("Deviation from Energy Conservation (#DeltaE) [GeV]");
      hist_deltaE->GetXaxis()->SetNdivisions(5, 3);
      hist_deltaE->SetMarkerSize(markerSize_full/2);
      hist_deltaE->Draw("p e0");
      drawText("#frac{1}{N^{event}} #frac{dN^{event}}{d#DeltaE}", 0.89, 0.89, true, kBlack, 0.05);
    }
    
    // Draw (Pseudo)Rapidity-Azimuth 2D hist
    gStyle->SetPalette(kGreyYellow);
    double max_full = getMaxFromHists(hist_y_phi_array[2*i], hist_eta_phi_array[2*i]);
    double max_ch = getMaxFromHists(hist_y_phi_array[2*i+1], hist_eta_phi_array[2*i+1]);
    double max_all = max_full;
    if (max_full < max_ch) max_all = max_ch;
    
    pads_right[0]->cd();
    pads_right[0]->SetTicks(1, 1);
//    pads_right[0]->SetTickY(1);
    pads_right[0]->SetFrameFillColor(color_bg);
    hist_y_phi_array[2*i]->GetZaxis()->SetRangeUser(0,max_all);
    
    hist_y_phi_array[2*i]->Draw("col");
    pads_right[1]->cd();
    pads_right[1]->SetTicks(1, 1);
//    pads_right[1]->SetTickY(1);
    pads_right[1]->SetFrameFillColor(color_bg);
    hist_eta_phi_array[2*i]->GetZaxis()->SetRangeUser(0,max_all);
    hist_eta_phi_array[2*i]->Draw("col");
    
    pads_right[2]->cd();
    pads_right[2]->SetTicks(1, 1);
//    pads_right[2]->SetTickY(1);
    pads_right[2]->SetFrameFillColor(color_bg);
    hist_y_phi_array[2*i+1]->GetZaxis()->SetRangeUser(0,max_all);
    hist_y_phi_array[2*i+1]->Draw("col");
    pads_right[3]->cd();
    pads_right[3]->SetTicks(1, 1);
//    pads_right[3]->SetTickY(1);
    pads_right[3]->SetFrameFillColor(color_bg);
    hist_eta_phi_array[2*i+1]->GetZaxis()->SetRangeUser(0,max_all);
    hist_eta_phi_array[2*i+1]->Draw("col");
    
    pads_right[4]->Clear();
    pads_right[4]->cd();
    hist_y_phi_array[2*i]->GetZaxis()->SetTitle("");
    hist_y_phi_array[2*i]->GetZaxis()->SetLabelSize(0.15);
    hist_y_phi_array[2*i]->GetZaxis()->SetLabelOffset(0.05);
    hist_y_phi_array[2*i]->GetZaxis()->SetTitleSize(0.15);
    hist_y_phi_array[2*i]->GetZaxis()->SetTitleOffset(1.5);
    TPaletteAxis* zscale = new TPaletteAxis(0.1, 0.015, 0.4, 1, hist_y_phi_array[2*i]);
    zscale->Draw();
    drawText(Form("#frac{1}{N^{event}} #frac{d^{2}N^{(ch) %s}}{dy(#eta)d#phi}",
                  descriptions[2*i]), 0.8, 0.99, true, kBlack, 0.15)->SetTextAngle(90);
    
    // Draw rapidity projection hist
    max_full = getMaxFromHists(proj_y_array[2*i], proj_eta_array[2*i]);
    max_ch = getMaxFromHists(proj_y_array[2*i+1], proj_y_array[2*i+1]);
    max_all = max_full;
    if (max_full < max_ch) max_all = max_ch;
    
    pads_right[7]->cd();
    pads_right[7]->SetTicks(1, 1);
    proj_y_array[2*i]->SetLineColor(color_full);
    proj_y_array[2*i]->SetFillColorAlpha(color_full, 0.3);
    proj_y_array[2*i]->GetXaxis()->SetLabelSize(0.06);
    proj_y_array[2*i]->GetXaxis()->SetTitle("Rapidity (y)");
    proj_y_array[2*i]->GetXaxis()->SetTitleSize(0.06);
    proj_y_array[2*i]->GetYaxis()->SetRangeUser(0, 1.2*max_all);
    proj_y_array[2*i]->Draw("hist");
    proj_y_array[2*i+1]->SetLineColor(color_ch);
    proj_y_array[2*i+1]->SetFillColorAlpha(color_ch, 0.3);
    proj_y_array[2*i+1]->Draw("hist same");
    
    pads_right[8]->cd();
    pads_right[8]->SetTicks(1, 1);
    proj_eta_array[2*i]->SetLineColor(color_full);
    proj_eta_array[2*i]->SetFillColorAlpha(color_full, 0.3);
    proj_eta_array[2*i]->GetXaxis()->SetLabelSize(0.06);
    proj_eta_array[2*i]->GetXaxis()->SetTitle("Pseudorapidity (#eta)");
    proj_eta_array[2*i]->GetXaxis()->SetTitleSize(0.06);
    proj_eta_array[2*i]->GetYaxis()->SetRangeUser(0, 1.2*max_all);
    proj_eta_array[2*i]->Draw("hist");
    proj_eta_array[2*i+1]->SetLineColor(color_ch);
    proj_eta_array[2*i+1]->SetFillColorAlpha(color_ch, 0.3);
    proj_eta_array[2*i+1]->Draw("hist same");
    
    // Draw azimuth projection hist
    max_full = getMaxFromHists(proj_phi_fromy_array[2*i], proj_phi_fromy_array[2*i+1]);
    max_ch = getMaxFromHists(proj_phi_frometa_array[2*i], proj_phi_frometa_array[2*i+1]);
    max_all = max_full;
    if (max_full < max_ch) max_all = max_ch;
    
    pads_right[5]->cd();
    pads_right[5]->SetTicks(2, 2);
    TMultiGraph* proj_phi_full_rotated = new TMultiGraph();
    proj_phi_fromy_array[2*i]->SetLineColor(color_full);
    proj_phi_fromy_array[2*i]->SetFillColorAlpha(color_full, 0.3);
    proj_phi_frometa_array[2*i]->SetLineColor(color_full);
    proj_phi_frometa_array[2*i]->SetFillColorAlpha(color_full, 0.3);
    proj_phi_full_rotated->Add(rotateHist(proj_phi_fromy_array[2*i], 1, 0, false));
    proj_phi_full_rotated->Add(rotateHist(proj_phi_frometa_array[2*i], 1, 0, false));
    proj_phi_full_rotated->Draw();
    proj_phi_full_rotated->GetXaxis()->SetLimits(0, max_all*1.35);
    proj_phi_full_rotated->GetXaxis()->SetNdivisions(4);
    proj_phi_full_rotated->GetYaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
    proj_phi_full_rotated->GetYaxis()->SetLabelOffset(
        3.5*proj_phi_full_rotated->GetYaxis()->GetLabelOffset());
    proj_phi_full_rotated->GetYaxis()->SetLabelSize(0.07);
    proj_phi_full_rotated->Draw("afl rx");
    drawText("Azimuth (#phi) [rad]", 0.99, 0.99, true, kBlack, 0.07)->SetTextAngle(90);
    
    pads_right[6]->cd();
    pads_right[6]->SetTicks(2, 2);
    TMultiGraph* proj_phi_ch_rotated = new TMultiGraph();
    proj_phi_fromy_array[2*i+1]->SetLineColor(color_ch);
    proj_phi_fromy_array[2*i+1]->SetFillColorAlpha(color_ch, 0.3);
    proj_phi_frometa_array[2*i+1]->SetLineColor(color_ch);
    proj_phi_frometa_array[2*i+1]->SetFillColorAlpha(color_ch, 0.3);
    proj_phi_ch_rotated->Add(rotateHist(proj_phi_fromy_array[2*i+1], 1, 0, false));
    proj_phi_ch_rotated->Add(rotateHist(proj_phi_frometa_array[2*i+1], 1, 0, false));
    proj_phi_ch_rotated->Draw();
    proj_phi_ch_rotated->GetXaxis()->SetLimits(0, max_all*1.35);
    proj_phi_ch_rotated->GetXaxis()->SetLabelOffset(0);
    proj_phi_ch_rotated->GetXaxis()->SetNdivisions(4);
    proj_phi_ch_rotated->GetXaxis()->SetLabelSize(0.07);
    proj_phi_ch_rotated->GetYaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
    proj_phi_ch_rotated->GetYaxis()->SetLabelOffset(
        3.5*proj_phi_ch_rotated->GetYaxis()->GetLabelOffset());
    proj_phi_ch_rotated->GetYaxis()->SetLabelSize(0.07);
    proj_phi_ch_rotated->Draw("afl rx");
    drawText("Azimuth (#phi) [rad]", 0.99, 0.99, true, kBlack, 0.07)->SetTextAngle(90);
    
    pads_right[9]->cd();
    pads_right[9]->Clear();
    TBox *box = new TBox(0.025, 0.025, 0.975, 0.975);
    box->SetFillColorAlpha(kGray+1, 0.8);
    box->Draw();
    drawText("#bf{Generator Settings}", 0.05, 0.9, false, kBlack, 0.075);
    drawText("#it{PYTHIA} tune", 0.05, 0.8, false, kBlack, 0.07);
    if ((int)settings_PYTHIA[0] == 6) drawText("4Cx", 0.48, 0.8, true, kBlack, 0.08);
    else drawText("DNE", 0.48, 0.9, true, kBlack, 0.08); // add other tune settings if used
    drawText("#it{PYTHIA} seed", 0.05, 0.7, false, kBlack, 0.07);
    if ((int)settings_PYTHIA[1] == -1) drawText("default", 0.48, 0.7, true, kBlack, 0.07);
    else drawText(Form("%i", (int)settings_PYTHIA[1]), 0.48, 0.7, true, kBlack, 0.08);
    drawText("#it{PYTHIA} #hat{#it{p}}_{T}^{min}", 0.05, 0.6, false, kBlack, 0.07);
    drawText(Form("%.2f", settings_PYTHIA[2]), 0.48, 0.6, true, kBlack, 0.08);
    drawText("#it{FastJet} algo", 0.05, 0.5, false, kBlack, 0.07);
    drawText(Form("%s", algo_string), 0.48, 0.5, true, kBlack, 0.07);
    drawText("#it{p}_{T}^{min, jet}", 0.05, 0.4, false, kBlack, 0.07);
    drawText(Form("%.2f", pt_binedge[0]), 0.48, 0.4, true, kBlack, 0.08);
    drawText("Jet Radius", 0.05, 0.3, false, kBlack, 0.08);
    drawText(Form("%.1f", jetRadius), 0.48, 0.3, true, kBlack, 0.08);
    drawText(Form("#bf{%s-Level Info}", descriptions[2*i+1]), 0.5, 0.9, false, kBlack, 0.075);
    drawText(Form("N^{%s}_{total}", descriptions[2*i]), 0.52, 0.8, false, kBlack, 0.08);
    drawText(Form("%i", grandtotal_n[2*i]), 0.95, 0.8, true, kBlack, 0.08);
    drawText(Form("N^{ch %s}_{total}", descriptions[2*i]), 0.52, 0.7, false, kBlack, 0.08);
    drawText(Form("%i", grandtotal_n[2*i+1]), 0.95, 0.7, true, kBlack, 0.08);
    drawText(Form("%s #it{p}_{T} resolution", descriptions[2*i]), 0.52, 0.6, false, kBlack, 0.07);
    if (i == 0) drawText(Form("%i", nbins_pT_jet), 0.95, 0.6, true, kBlack, 0.08);
    else drawText(Form("%i", nbins_pT_trk), 0.95, 0.6, true, kBlack, 0.08);
    drawText(Form("%s #phi resolution", descriptions[2*i]), 0.52, 0.5, false, kBlack, 0.07);
    drawText(Form("%i", nbins_phi), 0.95, 0.5, true, kBlack, 0.08);
    drawText(Form("%s y,#eta resolution", descriptions[2*i]), 0.52, 0.4, false, kBlack, 0.07);
    drawText(Form("%i", nbins_rap), 0.95, 0.4, true, kBlack, 0.08);
    
    // Write sheet for level i to file
    composite->Print(plotfile);
    composite->SaveAs(Form("gen_pythia/out_plots/lastgen/summarysheet_%ss.pdf", descriptions[2*i+1]));
    headerTex->Clear();
  }composite->Print(plotfile + "]");
  
  
  
  
  // Draw and write to file
//
//    gStyle->SetOptStat(0);
//    gPad->SetLogy();
//
//    // pT Spectra
//    pT_N->SetMarkerStyle(55);
//    pT_N->SetMarkerColor(kRed+1);
//    pT_N->Draw("P E0");
//    pT_NCharge->SetMarkerStyle(56);
//    pT_NCharge->SetMarkerColor(kBlue);
//    pT_NCharge->Draw("P E0 same");
//    canvas->Print(plotfile);
//
//    gPad->SetLogy(0);
//
//    // rapidity/pseudorapidity distributions
//    dndy->SetLineColor(kViolet+7);
//    dndy->Draw("hist");
//    dndy->SetTitle(";y (#eta);#frac{dN}{dy(#eta)}");
//    dndeta->SetLineColor(kOrange+4);
//  //  dndeta->SetLineStyle(7);
//    dndeta->Draw("hist same");
//    canvas->Print(plotfile);
//
//    // Multiplicity Event Distributions
//    nChargeHist->SetLineColor(kBlue-3);
//    nChargeHist->GetYaxis()->SetRangeUser(0, 0.097);
//  //  nChargeHist->SetMarkerStyle(54);
//    nChargeHist->Draw("hist");
//    nTotal->SetLineColor(kOrange-3);
//  //  nTotal->SetMarkerStyle(53);
//    nTotal->Draw("hist same");
//    canvas->Print(plotfile);
//
//    // Energy Deviation Histogram
//    gPad->SetLogx();
//    gPad->SetLogy();
//    DeltaE->SetLineColor(kPink);
//    DeltaE->Draw("hist");
//    canvas->Print(plotfile);
  
//
  
  delete fin;
  return;
}
