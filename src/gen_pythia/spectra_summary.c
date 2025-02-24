// ROOT macro. Should be soon incorperated into pythia_summary
// Just getting an early sense of the spectra shape

#include "../../utils/root_draw_tools.h"
#include "../../utils/hist_tools.h"
#include "config_local.h"
#include "ptbin.h"

// Color and aesthetic settings for root plots
Int_t palette_full = kCherry;
Int_t palette_ch = kLake;
Int_t markerStyle_full = 20;
Int_t markerStyle_ch = 33;
Float_t markerSize_full = 0.75;
Float_t markerSize_ch = 1.0;

// bool local_flag_xsec_scale = false;

TPad* pad_main[2];
TPad* pads_left[5];

void spectra_summary() {
  const int total_threads = getNbins(sqrt_s);
  float collision_energy_TeV = (float) sqrt_s / 1000.;
  char gen_info_string[100];
  snprintf(gen_info_string, 100, "%.2fTeV_r%.1fjet", collision_energy_TeV, jet_radius);
  
  TH1D* hist_fulljet[total_threads];
  TH1D* hist_charjet[total_threads];
  TH1D* hist_alljet[2];
  int max_pt;
  if (sqrt_s < 900) max_pt = 60;
  else              max_pt = 360;
  
  hist_alljet[0] = new TH1D("hist_pT_fulljet",";p_{T};1/N^{event} dN^{jet}/dp_{T}",20, 0, max_pt);
  hist_alljet[1] = new TH1D("hist_pT_charjet",";p_{T};1/N^{event} dN^{charged jet}/dp_{T}",20, 0, max_pt);
  
  TFile* write_file = new TFile(Form("out_data/pythpp_%s.root", gen_info_string),"recreate");
  
  // Merged TTree for storing jet information in each event
  int                     pt_hard_bin;        // Which bin this event belongs to
  int                     mult_merge;         // Total event track multiplicity
  int                     cmult_merge;        // Event charged track multiplicity
  float                   delta_E_merge;      // Deviation from energy conservation
  int                     ntot_jet_merge;     // Number of full jets in an event
  std::vector<int>        ntot_cst_merge;     // Array of full jet constituent multiplicties
  std::vector<float>      jet_pT_merge;       // Array of full jet p_T
  std::vector<float>      jet_y_merge;        // Array of full jet rapidity
  std::vector<float>      jet_eta_merge;      // Array of full jet pseudorapidity
  std::vector<float>      jet_phi_merge;      // Array of full jet azimuth
  std::vector<float>      jet_area_merge;     // Array of full jet area
  int                     ntot_cjet_merge;    // Number of charged jets in an event
  std::vector<int>        ntot_ccst_merge;    // Array of charged jet constituent multiplicities
  std::vector<float>      chjet_pT_merge;     // Array of charged jet p_T
  std::vector<float>      chjet_y_merge;      // Array of charged jet rapidity
  std::vector<float>      chjet_eta_merge;    // Array of charged jet pseudorapidity
  std::vector<float>      chjet_phi_merge;    // Array of charged jet azimuth
  std::vector<float>      chjet_area_merge;   // Array of charged jet area
  
  TTree* pythia_tree_merged = new TTree("pythia_tree","pythia_tree");
  pythia_tree_merged->Branch("pt_bin",     &pt_hard_bin);
  pythia_tree_merged->Branch("ntotal",     &mult_merge);
  pythia_tree_merged->Branch("ncharge",    &cmult_merge);
  pythia_tree_merged->Branch("delta_E",    &delta_E_merge);
  pythia_tree_merged->Branch("njet",       &ntot_jet_merge);
  pythia_tree_merged->Branch("nchjet",     &ntot_cjet_merge);
  pythia_tree_merged->Branch("jet_n",      &ntot_cst_merge);
  pythia_tree_merged->Branch("chjet_n",    &ntot_ccst_merge);
  pythia_tree_merged->Branch("jet_pT",     &jet_pT_merge);
  pythia_tree_merged->Branch("chjet_pT",   &chjet_pT_merge);
  pythia_tree_merged->Branch("jet_y",      &jet_y_merge);
  pythia_tree_merged->Branch("chjet_y",    &chjet_y_merge);
  pythia_tree_merged->Branch("jet_eta",    &jet_eta_merge);
  pythia_tree_merged->Branch("chjet_eta",  &chjet_eta_merge);
  pythia_tree_merged->Branch("jet_phi",    &jet_phi_merge);
  pythia_tree_merged->Branch("chjet_phi",  &chjet_phi_merge);
  pythia_tree_merged->Branch("jet_area",   &jet_area_merge);
  pythia_tree_merged->Branch("chjet_area", &chjet_area_merge);
  
  
  // Print some information about this dataset
  std::cout << Form("PYTHIA pp sqrt_s_NN = %.2f TeV, pT_hadron >= %.2f GeV, |eta_hadron| < %.2f",
                    collision_energy_TeV, pTmin_hadron, max_eta) << std::endl;
  std::cout << Form("FastJet3 R = %.1f %s Jets, |eta_jet| < %.2f, jet_A >= %.2f*pi*R^2",
                    jet_radius,algo_string,max_eta - jet_radius,min_area/(3.141592*jet_radius*jet_radius)) << std::endl;
  std::cout << Form("pT_jet >= %.2f, pT_charged_jet >= %.2f, pT_leading >= %.2f",
                    jetcut_minpT_full,jetcut_minpT_chrg, pTmin_jetcore) << std::endl;
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  // TCanvas Setup
  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas();
  canvas->SetCanvasSize(350,500);
  std::vector<std::vector<TPad*>> pads = divideFlush(canvas, 1, 2);
  TLegend* leg_pTHard = new TLegend(0.7, 0.35, 0.85, 0.95);
  leg_pTHard->SetLineWidth(0);
  
  
  // Colors from root spectra
  gStyle->SetPalette(kBird);
  Int_t color_sample_fulljet[total_threads];
  for (Int_t i = 0; i < total_threads; ++i) color_sample_fulljet[i] = TColor::GetPalette().At(TColor::GetPalette().GetSize() * i/(1.2*total_threads));
  gStyle->SetPalette(kBird);
  Int_t color_sample_charjet[total_threads];
  for (Int_t i = 0; i < total_threads; ++i) color_sample_charjet[i] = TColor::GetPalette().At(TColor::GetPalette().GetSize() * i/(1.2*total_threads));
  
  for (int this_thread = total_threads-1; this_thread >= 0; --this_thread) {
    pt_hard_bin = this_thread;
    pTBinSettings settings = getBinSettings(sqrt_s, this_thread, local_flag_xsec_scale);
    const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
    const float cross_section = settings.cross_section;
    const int nevent = (int) (settings.nevent * nevent_scale_factor); // scaled event number for testing
    
    // File for this thread
    TFile* current_file = new TFile(Form("out_data/pythpp_%s_thread%i.root", gen_info_string, this_thread), "read");
    
    // Find data tree
    char tree_name[50];
    snprintf(tree_name, 50, "pythia_tree_pTHard%.f-%.f", pt_binedge[0], pt_binedge[1]);
    TTreeReader *reader = new TTreeReader(tree_name, current_file);
    if (!reader) {std::cout << "Tree not found!" << std::endl; return;}
    TTreeReaderValue<Int_t>                     ntotal(*reader, "ntotal");          // Total event track multiplicity
    TTreeReaderValue<Int_t>                     ncharge(*reader, "ncharge");        // Event charged track multiplicity
    TTreeReaderValue<Float_t>                   deltaE(*reader, "delta_E");          // Deviation from energy conservation
    TTreeReaderValue<Int_t>                     njet(*reader, "njet");              // Number of full jets in an event
    TTreeReaderValue<Int_t>                     nchjet(*reader, "nchjet");          // Number of charged jets in an event
    TTreeReaderValue<std::vector<Int_t>>        jet_n(*reader, "jet_n");            // Array of full jet constituent multiplicties
    TTreeReaderValue<std::vector<Int_t>>        chjet_n(*reader, "chjet_n");        // Array of charged jet constituent multiplicities
    TTreeReaderValue<std::vector<Float_t>>      jet_pT(*reader, "jet_pT");          // Array of full jet p_T
    TTreeReaderValue<std::vector<Float_t>>      chjet_pT(*reader, "chjet_pT");      // Array of charged jet p_T
    TTreeReaderValue<std::vector<Float_t>>      jet_y(*reader, "jet_y");            // Array of full jet rapidity
    TTreeReaderValue<std::vector<Float_t>>      chjet_y(*reader, "chjet_y");        // Array of charged jet rapidity
    TTreeReaderValue<std::vector<Float_t>>      jet_eta(*reader, "jet_eta");        // Array of full jet pseudorapidity
    TTreeReaderValue<std::vector<Float_t>>      chjet_eta(*reader, "chjet_eta");    // Array of charged jet pseudorapidity
    TTreeReaderValue<std::vector<Float_t>>      jet_phi(*reader, "jet_phi");        // Array of full jet azimuth
    TTreeReaderValue<std::vector<Float_t>>      chjet_phi(*reader, "chjet_phi");    // Array of charged jet azimuth
    TTreeReaderValue<std::vector<Float_t>>      jet_area(*reader, "jet_area");      // Array of full jet areas
    TTreeReaderValue<std::vector<Float_t>>      chjet_area(*reader, "chjet_area");  // Array of charged jet areas
    
    
    
    // Set up histograms
    hist_fulljet[this_thread] = new TH1D(Form("hist_pT_fulljet_thread%i",this_thread),
                                         ";p_{T};1/N^{event} dN^{jet}/dp_{T}",
                                         20, 0, max_pt);
    hist_fulljet[this_thread]->SetLineColor(color_sample_fulljet[this_thread]);
    hist_fulljet[this_thread]->SetLineWidth(2);
    hist_fulljet[this_thread]->SetFillColorAlpha(color_sample_fulljet[this_thread], 0.2);
    hist_charjet[this_thread] = new TH1D(Form("hist_pT_charjet_thread%i",this_thread),
                                         ";p_{T};1/N^{event} dN^{charged jet}/dp_{T}",
                                         20, 0, max_pt);
    hist_charjet[this_thread]->SetLineColor(color_sample_charjet[this_thread]);
    hist_charjet[this_thread]->SetLineWidth(2);
    hist_charjet[this_thread]->SetFillColorAlpha(color_sample_charjet[this_thread], 0.2);
    
    // Gather data
    int ntot_jet = 0;
    int ntot_cjet = 0;
    while (reader->Next()) {
      ntot_jet += *njet;
      ntot_cjet += *nchjet;
      
      // Assign data from current branh into the merged tree
      mult_merge =        *ntotal;
      cmult_merge =       *ncharge;
      delta_E_merge =     *deltaE;
      ntot_jet_merge =    *njet;
      ntot_cst_merge =    *jet_n;
      jet_pT_merge =      *jet_pT;
      jet_y_merge =       *jet_y;
      jet_eta_merge =     *jet_eta;
      jet_phi_merge =     *jet_phi;
      jet_area_merge =    *jet_area;
      ntot_cjet_merge =   *nchjet;
      ntot_ccst_merge =   *chjet_n;
      chjet_pT_merge =    *chjet_pT;
      chjet_y_merge =     *chjet_y;
      chjet_eta_merge =   *chjet_eta;
      chjet_phi_merge =   *chjet_phi;
      pythia_tree_merged->Fill();
      
      // Fill histograms
      for (int i_jet = 0; i_jet < *njet; ++i_jet) {
        hist_fulljet[this_thread]->Fill(jet_pT->at(i_jet));
      } for (int i_jet = 0; i_jet < *nchjet; ++i_jet) {
        hist_charjet[this_thread]->Fill(chjet_pT->at(i_jet));
      }
    }
    
    std::cout << Form("pTHat bin [%.f,%.f]    \tn_event = %i\tn_jet = %i\tn_charged_jet = %i", pt_binedge[0],pt_binedge[1],nevent,ntot_jet,ntot_cjet) << std::endl;
    
    // Normalize by n_event
    float max_range;
    if (sqrt_s < 900) max_range = 0.15;
    else max_range = 0.025;
    hist_fulljet[this_thread]->Scale(1./((float)nevent), "width");
    hist_fulljet[this_thread]->GetYaxis()->SetRangeUser(0, max_range);
    hist_charjet[this_thread]->Scale(1./((float)nevent), "width");
    hist_charjet[this_thread]->GetYaxis()->SetRangeUser(0, max_range);
    
    // Add to net spectrum
    hist_alljet[0]->Add(hist_alljet[0], hist_fulljet[this_thread], 1, cross_section);
    hist_alljet[1]->Add(hist_alljet[1], hist_charjet[this_thread], 1, cross_section);
    
    // Plot spectra
    pads.at(0).at(0)->cd();
    if (this_thread == total_threads-1) {
      gPad->SetTicks(1,1);
      hist_fulljet[this_thread]->Draw("hist");
    } else
      hist_fulljet[this_thread]->Draw("hist same");
    
    
    pads.at(1).at(0)->cd();
    if (this_thread ==  total_threads-1) {
      gPad->SetTicks(1,1);
      hist_charjet[this_thread]->Draw("hist");
    } else
      hist_charjet[this_thread]->Draw("hist same");
    
    leg_pTHard->AddEntry(hist_fulljet[this_thread], Form("[%.f,%.f]",pt_binedge[0],pt_binedge[1]), "f");
    
    //    current_file->Close();
    //    delete current_file;
    //    std::cout << "Finished with thread #" << this_thread << std::endl;
  }
  pads.at(0).at(0)->cd();
  drawText(Form("Full Jet p_{T} #geq %.1f",jetcut_minpT_full), 0.15, 0.82, false);
  drawText(Form("#it{PYTHIA} #bf{pp} #sqrt{s_{NN}} = %.2f", (float) sqrt_s / 1000.), 0.9, 0.82, true);
  drawText(Form("#it{FastJet} %s R = %.1f", algo_string, jet_radius), 0.9, 0.76, true);
  drawText(Form("#left|#eta_{hadron}#right| #leq %.1f", max_eta), 0.9, 0.70, true);
  drawText(Form("p_{T}^{leading} #geq %.1f", pTmin_jetcore), 0.9, 0.62, true);
  
  pads.at(1).at(0)->cd();
  drawText(Form("Charged Jet p_{T} #geq %.1f",jetcut_minpT_chrg), 0.15, 0.92, false);
  leg_pTHard->SetHeader("p_{T}^{Hard} Bin\n");
  leg_pTHard->Draw();
  canvas->SaveAs(Form("out_plots/spectra_normalized_%s.pdf",gen_info_string));
  
  
  
  
  // Draw the sigma scaled version
  for (int this_thread = 0; this_thread < total_threads; ++this_thread) {
    pTBinSettings settings = getBinSettings(sqrt_s, this_thread, local_flag_xsec_scale);
    const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
    const float cross_section = settings.cross_section;
    const int nevent = (int) (settings.nevent * nevent_scale_factor); // scaled event number for testing
    float collision_energy_TeV = (float) sqrt_s / 1000.;
    
    hist_fulljet[this_thread]->SetTitle(";p_{T};#sigma/N^{event} dN^{jet}/dp_{T}");
    hist_charjet[this_thread]->SetTitle(";p_{T};#sigma/N^{event} dN^{charged jet}/dp_{T}");
    
    // scale by cross section
    hist_fulljet[this_thread]->Scale(cross_section);
    hist_fulljet[this_thread]->GetYaxis()->SetRangeUser(1e-9, 1.5e-2);
    hist_charjet[this_thread]->Scale(cross_section);
    hist_charjet[this_thread]->GetYaxis()->SetRangeUser(1e-9, 1.5e-2);
    
    // Plot spectra
    pads.at(0).at(0)->cd();
    if (this_thread == 0) {
      gPad->SetLogy();
      hist_fulljet[this_thread]->Draw("hist");
    } else
      hist_fulljet[this_thread]->Draw("hist same");
    
    
    pads.at(1).at(0)->cd();
    if (this_thread == 0) {
      gPad->SetLogy();
      hist_charjet[this_thread]->Draw("hist");
    } else
      hist_charjet[this_thread]->Draw("hist same");
  }
  
  
  pads.at(0).at(0)->cd();
  hist_alljet[0]->SetLineColor(kBlack);
  hist_alljet[0]->Draw("hist same");
  drawText(Form("Full Jet p_{T} #geq %.1f",jetcut_minpT_full), 0.15, 0.82, false);
  drawText(Form("#it{PYTHIA} #bf{pp} #sqrt{s_{NN}} = %.2f", (float) sqrt_s / 1000.), 0.9, 0.82, true);
  drawText(Form("#it{FastJet} %s R = %.1f", algo_string, jet_radius), 0.9, 0.76, true);
  drawText(Form("#left|#eta_{hadron}#right| #leq %.1f", max_eta), 0.9, 0.70, true);
  drawText(Form("p_{T}^{leading} #geq %.1f", pTmin_jetcore), 0.9, 0.62, true);
  
  pads.at(1).at(0)->cd();
  hist_alljet[1]->SetLineColor(kBlack);
  hist_alljet[1]->Draw("hist same");
  leg_pTHard->AddEntry(hist_alljet[1], "Net Spectrum", "l");
  drawText(Form("Charged Jet p_{T} #geq %.1f",jetcut_minpT_chrg), 0.15, 0.92, false);
  leg_pTHard->SetHeader("p_{T}^{Hard} Bin\n");
  leg_pTHard->Draw();
  canvas->SaveAs(Form("out_plots/spectra_sigmascaled_%s.pdf",gen_info_string));
  
  
  // Save the completed spectra and full tree to a composite output file
  write_file->cd();
  pythia_tree_merged->Print();
  pythia_tree_merged->Write(pythia_tree_merged->GetName(), TObject::kOverwrite);
  
  hist_alljet[0]->Write(hist_alljet[0]->GetName(), TObject::kOverwrite);
  for (int i_thread = 0; i_thread < total_threads; ++i_thread)
    hist_fulljet[i_thread]->Write(hist_fulljet[i_thread]->GetName(), TObject::kOverwrite);
  hist_alljet[1]->Write(hist_alljet[1]->GetName(), TObject::kOverwrite);
  for (int i_thread = 0; i_thread < total_threads; ++i_thread)
    hist_charjet[i_thread]->Write(hist_charjet[i_thread]->GetName(), TObject::kOverwrite);
  
  
  delete write_file;
  return;
}
