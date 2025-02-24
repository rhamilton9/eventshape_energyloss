// ROOT macro. Generates energyloss estimate by shifting reference pT spectrum and taking a ratio with original spectrum
// This estimate for energy loss is obtained by fitting this ratio against a known reference R_AA
// Assumes relevant data spectra have been correctly formatted and are accessible in ../data.nosync

/// ROOT macro for generating an energy loss estimate by comparing pT spectra in PYTHIA monte carlo against data spectra.
/// Expects a spectrum that has been generated with proper weighting to bins of the PYTHIA parameter pTHat.
/// A data spectrum should also be accessible, possibly from HEPdata ( https://www.hepdata.net ).
///
/// For info on proper generation of input spectra, refer to the code for PYTHIA generation in gen-pythia/readme.txt

#include "../config.h"
#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"

void energyloss_pTspectra() {
  const int n_dpt = (int) max_dpt / dpT_resolution;
  const int n_raa = sizeof(centrality_list) / sizeof(int[2]);
//  const int n_raa = 1;
  char species[5];
  snprintf(species, 5, "%s%s", speciesA, speciesB);
  const char string_truefalse[2][10] = {"true","false"};
  
  TFile* infile_raa = new TFile(Form("../data.nosync/%s_%.2fTeV/spectra/%s_%s%.2fTeV_unpacked_%i.root",
                                species, sqrt_s, experiment, species, sqrt_s, dataset));
  TFile* infile_ref = new TFile(Form("../data.nosync/pp_reference_%i.root", dataset));
  
  // Set up data according to inputs
  
  // 2018 dataset.
  //    - 2.76 TeV: 1, 7
  //    - 5.02 TeV: 2, 8
  int pt_table;
  int raa_table;
  // For 2013 reference:
  //    - Hist1D_y2_1 is 2.76 TeV
  // For 2018 reference:
  //    - Hist1D_y1_4 is 5.02 TeV
  //    - Hist1D_y2_4 is 2.76 TeV
  TH1F *pp_reference;
  if (dataset == 2013) {
    pp_reference = infile_ref->Get<TH1F>("Hist1D_y2_1");
    pp_reference->SetTitle(";p_{T} [GeV];1/2#pip_{T} d^{2}#sigma/d#etadp_{T} [mb/GeV^{2}]");
    std::cout << "debug" << std::endl;
  } else if (dataset == 2018) {
    if (sqrt_s == 2.76) {
      pp_reference = infile_ref->Get<TH1F>("Hist1D_y2_4");
      pt_table = 1;
      raa_table = 7;
    } else if (sqrt_s == 5.02) {
      pp_reference = infile_ref->Get<TH1F>("Hist1D_y1_4");
      pt_table = 2;
      raa_table = 8;
    } else {
      std::cout << "Error in <src/energyloss_pTspectra>: Requested CM energy sqrt_s = " <<
        sqrt_s << " is not available in 2018 dataset." << std::endl;
      return;
    }pp_reference->SetTitle(";p_{T} [GeV];1/N^{evt} d^{2}N/D#eta/dp_{T} [C/GeV]");
  } else {
    std::cout<<Form("Error in <src/energyloss_pTspectra>: No data available for %i.",dataset)<<std::endl;
    return;
  }pp_reference->SetLineColor(kBlack);
  
  
  
  
  TH1F* data_pt[n_raa];
  TH1F* data_raa[n_raa];
  TH1F* data_totalunc[n_raa];
  for (int i = 0; i < n_raa; ++i) {
    if (dataset == 2013) {
      data_pt[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_%i",i+1)));
      data_raa[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_%i",i+index_raa_in_file)));
    } else if (dataset == 2018) {
      data_pt[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y%i_%i",i+1, pt_table)));
      data_raa[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y%i_%i",i+1, raa_table)));
    }

    
    data_raa[i]->SetTitle(";p_{T} [GeV];R_{AA}");
    
    if (data_errors_sumbyquadrature) {
      //loop, compute by quadrature
    } else {
      if (dataset == 2013) {
        data_totalunc[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_e1_%i",i+index_raa_in_file)));
        data_totalunc[i]->Add(data_totalunc[i],
                              static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_e2_%i",i+index_raa_in_file))));
        data_totalunc[i]->Add(data_totalunc[i],
                              static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y1_e3_%i",i+index_raa_in_file))));
      } else if (dataset == 2018) {
        data_totalunc[i] = static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y%i_e1_%i",i+1, raa_table)));
        data_totalunc[i]->Add(data_totalunc[i],
                              static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y%i_e2_%i",i+1, raa_table))));
        data_totalunc[i]->Add(data_totalunc[i],
                              static_cast<TH1F*>(infile_raa->Get(Form("Hist1D_y%i_e3_%i",i+1, raa_table))));
      }
    }
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
  // TODO branch other energyloss estimates
  
  // Canvas setup
  TCanvas *canvas = new TCanvas();
//  canvas->SetCanvasSize(1000, 1000);
  canvas->SetWindowSize(1200, 1600);
  
  TPad* mainpad = buildPad("mainpad", 0.07, 0, 1, 0.9, 0, 0, 0, 0);
  mainpad->Divide(4,3);
  gStyle->SetOptStat(0);
  
  
  
  // Energyloss computation
  const double dpt = -dpT_resolution;
  const int nbin = pp_reference->GetXaxis()->GetNbins();
  double cval;
  double cmin_dpt[3];
  TH1F *opt_raa_hist[3];
  TH1F *opt_pt_hist[3];
 
  int ntest = 0;
  for (int iCent = ntest; iCent < ntest+1; ++iCent) {
    // Setup hists for optimizer storage
    TH1F* chi2_hist = new TH1F("",";#Delta p_{T};#chi^{2}",n_dpt, 0, n_dpt*(-dpt));
    TH1F* ks_r_hist = new TH1F("",";#Delta p_{T};KS",      n_dpt, 0, n_dpt*(-dpt));
    TH1F* ks_n_hist = new TH1F("",";#Delta p_{T};KS",      n_dpt, 0, n_dpt*(-dpt));
    
    double cmin_chi2 = 1e50;
    double cmin_ks_r = 1;
    double cmin_ks_n = 1;
    
    //--------------------------------------------------------------Perform fitting
    
    TH1F* placeholder_hist;
    TH1F* placeholder_ratio = static_cast<TH1F*>(pp_reference->Clone());
    for (int ipt = 1; ipt <= n_dpt; ++ipt) {
      placeholder_ratio->Reset();
      
      
      
      //------------------------------------------------------------Rebinning/Ratio
      
      // shift placeholder by dpt, rebin.
      // Take ratio against reference to compute phenom R_AA
      placeholder_hist = static_cast<TH1F*>(translateHist(pp_reference, dpt*ipt, 0, 0, 12));
      for (int ibin = 1; ibin <= nbin; ++ibin) {
        cval = placeholder_hist->GetBinContent(ibin);
        if (cval == 0) {placeholder_ratio->SetBinContent(ibin, 0); continue;}
        placeholder_ratio->SetBinContent(ibin, cval / pp_reference->GetBinContent(ibin));
      }
      
      //------------------------------------------------------------Chi^2 calculation
      
      // Compute chi^2 against data R_AA
      double cpt;
      double chi2 = 0;
      for (int ibin = 1; ibin <= data_raa[iCent]->GetXaxis()->GetNbins(); ++ibin) {
        cpt = data_raa[iCent]->GetXaxis()->GetBinCenter(ibin);
        if (cpt < minpt_comparison_threshold) continue;
        chi2 += TMath::Power(placeholder_ratio->GetBinContent(placeholder_ratio->FindBin(cpt))
                             - data_raa[iCent]->GetBinContent(ibin), 2)
                             / data_totalunc[iCent]->GetBinContent(ibin);
      }chi2_hist->SetBinContent(ipt, chi2);
      
      if (chi2 < cmin_chi2) {
        cmin_chi2 = chi2;
        cmin_dpt[0] = -dpt*ipt;
        opt_pt_hist[0] = static_cast<TH1F*>(placeholder_hist->Clone());
        opt_raa_hist[0] = static_cast<TH1F*>(placeholder_ratio->Clone());
      }
      
      //------------------------------------------------------------KS calculation (rebinned)

      
      //new, with pT
      double ks_r = KS_statistic_new(placeholder_hist, data_pt[iCent], 0, minpt_comparison_threshold, INT_MAX,
                                     false, Form("../tmp/tmpplot/cent_%i-%i%%/rebin/ks_r_Cent%i-%i%%_i%i",
                                                                 centrality_list[iCent][0],centrality_list[iCent][1],
                                                                 centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      // new, with raa
//      double ks_r = KS_statistic_new(placeholder_ratio, data_raa[iCent], 0, minpt_comparison_threshold, INT_MAX,
//                                 (ipt < 20 || ipt%20 == 0), Form("../tmp/tmpplot/cent_%i-%i%%/rebin/ks_r_Cent%i-%i%%_i%i",
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1],
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      //old, with pt, rebin
//      double ks_r = KS_statistic(placeholder_hist, data_pt[iCent], 0, minpt_comparison_threshold, true,
//                                 false, Form("../tmp/tmpplot/cent_%i-%i%%/rebin/ks_r_Cent%i-%i%%_i%i",
//                                             centrality_list[iCent][0],centrality_list[iCent][1],
//                                             centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      //old, with raa, rebin
//      double ks_r = KS_statistic(placeholder_ratio, data_raa[iCent], 0, minpt_comparison_threshold, true,
//                                 false, Form("../tmp/tmpplot/cent_%i-%i%%/rebin/ks_r_Cent%i-%i%%_i%i",
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1],
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      
//      TH1* hist1,
//      TH1* hist2,
//      double horizontal_shift = 0,
//      double thresh_min = INT_MIN,
//      double thresh_max = INT_MAX,
//      bool doPlot = false,
//      const char *saveName = (char*)"ks",
//      int iteration = -1
      
      ks_r_hist->SetBinContent(ipt, ks_r);
      if (ks_r < cmin_ks_r) {
        cmin_ks_r = ks_r;
        cmin_dpt[1] = -dpt*ipt;
        opt_pt_hist[1] = static_cast<TH1F*>(placeholder_hist->Clone());
        opt_raa_hist[1] = static_cast<TH1F*>(placeholder_ratio->Clone());
      }
      
      //------------------------------------------------------------KS calculation (no rebin)
      
      //old, with pt, up
//      double ks_n = KS_statistic(data_pt[iCent], pp_reference, -dpt*ipt, minpt_comparison_threshold-dpt*ipt, true,
//                                 false , Form("../tmp/tmpplot/cent_%i-%i%%/norebin/ks_n_Cent%i-%i%%_i%i",
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1],
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      //old, with pt, down
//      double ks_n = KS_statistic(pp_reference, data_pt[iCent], dpt*ipt, minpt_comparison_threshold, true,
//                                 false , Form("../tmp/tmpplot/cent_%i-%i%%/norebin/ks_n_Cent%i-%i%%_i%i",
//                                              centrality_list[iCent][0],centrality_list[iCent][1],
//                                              centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      //new, with pt, upshift
      double ks_n = KS_statistic_new(data_pt[iCent], pp_reference, -dpt*ipt, minpt_comparison_threshold-dpt*ipt, INT_MAX,
                                     false, Form("../tmp/tmpplot/cent_%i-%i%%/norebin/ks_n_Cent%i-%i%%_i%i",
                                                                 centrality_list[iCent][0],centrality_list[iCent][1],
                                                                 centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      //new, with pt, downshift
//      double ks_n = KS_statistic_new(pp_reference, data_pt[iCent], dpt*ipt, minpt_comparison_threshold, INT_MAX,
//                                 (ipt < 20 || ipt%20 == 0), Form("../tmp/tmpplot/cent_%i-%i%%/norebin/ks_n_Cent%i-%i%%_i%i",
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1],
//                                                                 centrality_list[iCent][0],centrality_list[iCent][1], ipt), ipt);
      ks_n_hist->SetBinContent(ipt, ks_n);
      if (ks_n < cmin_ks_n) {
        cmin_ks_n = ks_n;
        cmin_dpt[2] = -dpt*ipt;
        opt_pt_hist[2] = static_cast<TH1F*>(placeholder_hist->Clone());
        opt_raa_hist[2] = static_cast<TH1F*>(placeholder_ratio->Clone());
      }
      
      
    }// End pT loop
    
    //--------------------------------------------------------------Plotting results
    
    TLine* thresholdLine = new TLine();
    thresholdLine->SetLineColor(kGray+1);
    thresholdLine->SetLineStyle(9);
    
    TLine* optLine = new TLine();
    optLine->SetLineColor(kOrange+7);
    optLine->SetLineStyle(3);
    
    canvas->cd();
    TLatex* tex_title = drawText(Form("Energyloss Fitting Summary for %s #bf{%s}, #sqrt{s} = %.2fTeV",
                                      experiment, species, sqrt_s), 0.02, 0.95);
    TLatex* tex_subtitle;
    tex_subtitle = drawText(Form("Centrality bin %i-%i%%, data_errors_sumbyquadrature = %s",
                                 centrality_list[iCent][0], centrality_list[iCent][1], string_truefalse[!data_errors_sumbyquadrature]),
                            0.02, 0.91, false, kBlack, 0.03);
    drawText("#chi^{2} Fitting", 0.065, 0.7, false, kBlack, 0.03)->SetTextAngle(90);
    drawText("KS Fit, Hist Rebin", 0.065, 0.35, false, kBlack, 0.03)->SetTextAngle(90);
    drawText("KS Fit, No Rebin", 0.065, 0.05, false, kBlack, 0.03)->SetTextAngle(90);
    
    mainpad->cd(1);
//    gPad->SetLogx();
    gPad->SetLogy();
    pp_reference->Draw("hist");
    TH1F* datCompPlot[3];
    datCompPlot[0] = static_cast<TH1F*>(data_pt[iCent]->Clone());
    datCompPlot[0]->SetLineColor(kViolet+5);
    datCompPlot[0]->Draw("hist same");
    datCompPlot[0]->Scale(definiteIntegral(opt_pt_hist[0], minpt_comparison_threshold)/
                          definiteIntegral(datCompPlot[0], minpt_comparison_threshold));
    opt_pt_hist[0]->SetLineColor(kRed+2);
    opt_pt_hist[0]->Draw("hist same");
    thresholdLine->DrawLine(minpt_comparison_threshold, 0, minpt_comparison_threshold,
                            pp_reference->GetMaximum()*1.1);
    TLegend *baseLeg_top = new TLegend(0.5, 0.68, 0.88, 0.89);
    baseLeg_top->SetLineWidth(0);
    baseLeg_top->AddEntry(pp_reference, "#bf{pp} Reference p_{T} Spectrum ", "l");
    baseLeg_top->AddEntry(datCompPlot[0], Form("Rescaled #bf{%s} p_{T} Data", species), "l");
    baseLeg_top->AddEntry(opt_pt_hist[0], Form("Fit Result #Deltap_{T} = %.2f", cmin_dpt[0]), "l");
    baseLeg_top->AddEntry(thresholdLine, "p_{T} Comparison Threshold", "l");
    baseLeg_top->Draw();
    
    mainpad->cd(2);
    if (dpt < 0.05) chi2_hist->Draw("hist l");
    else            chi2_hist->Draw("hist");
    optLine->DrawLine(cmin_dpt[0], 0, cmin_dpt[0], chi2_hist->GetMaximum()*1.05);
    if (cmin_dpt[0] > max_dpt/2) drawText(Form("#Deltap_{T} = %.2f", cmin_dpt[0]), cmin_dpt[0]-0.1,
                                          chi2_hist->GetMaximum()*0.95, true, kBlack, 0.04, 42, false);
    else                         drawText(Form("#Deltap_{T} = %.2f", cmin_dpt[0]), cmin_dpt[0]+0.1,
                                          chi2_hist->GetMaximum()*0.95, false, kBlack, 0.04, 42, false);
    
    mainpad->cd(3);
    //  gPad->SetLogx();
      gPad->SetLogy();
    data_raa[iCent]->SetLineColor(kBlack);
    data_raa[iCent]->Draw("hist");
    data_raa[iCent]->GetYaxis()->SetRangeUser(0.05, 1);
    opt_raa_hist[0]->SetLineColor(kRed+2);
    opt_raa_hist[0]->Draw("hist same");
    thresholdLine->DrawLine(minpt_comparison_threshold, data_raa[iCent]->GetMinimum(),
                            minpt_comparison_threshold, data_raa[iCent]->GetMaximum()*1.05);
    TLegend *baseLeg_low = new TLegend(0.4, 0.15, 0.88, 0.35);
    baseLeg_low->SetLineWidth(0);
    baseLeg_low->SetFillColorAlpha(kWhite, 0.85);
    baseLeg_low->AddEntry(data_raa[iCent], Form("#bf{%s} R_{AA} Data", species), "l");
    baseLeg_low->AddEntry(opt_raa_hist[0], Form("Ratio #frac{#bf{pp} Reference}{Fit Result #Deltap_{T} = %.2f}", cmin_dpt[0]), "l");
    baseLeg_low->Draw();
    
    mainpad->cd(4);
    // Construct double ratio
    TAxis* refaxis = pp_reference->GetXaxis();
    int nbins_doubleratio = refaxis->GetNbins() - refaxis->FindBin(minpt_comparison_threshold+cmin_dpt[0])+1;
    double binedge_ratio_1[nbins_doubleratio+1];
    for (int irbin = 0; irbin <= nbins_doubleratio; ++irbin)
      binedge_ratio_1[irbin] = refaxis->GetBinLowEdge(refaxis->GetNbins() - nbins_doubleratio + irbin + 1);
    TH1D* doubleratio = new TH1D("",";p_{T};Double Ratio (Data / Shifted)",nbins_doubleratio, binedge_ratio_1);
    double raa, praa;
    for (int irbin = 1; irbin <= nbins_doubleratio; ++irbin) {
      raa =  data_raa[iCent]->GetBinContent(data_raa[iCent]->FindBin(doubleratio->GetBinCenter(irbin)));
      praa = opt_raa_hist[0]->GetBinContent(opt_raa_hist[0]->FindBin(doubleratio->GetBinCenter(irbin)));
      std::cout << "bin " << doubleratio->GetBinCenter(irbin) << ", \tRAA = " << raa << ", \tpRAA = " << praa << std::endl;
      if (raa == 0 || praa == 0) continue;
      doubleratio->SetBinContent(irbin, raa/praa);
    }setStyleLine(doubleratio, "violet thin");
    doubleratio->Draw("hist");
    drawUnityLine(doubleratio->GetXaxis());
    
    mainpad->cd(5);
    gPad->SetLogy();
    pp_reference->Draw("hist");
    datCompPlot[1] = static_cast<TH1F*>(data_pt[iCent]->Clone());
    datCompPlot[1]->SetLineColor(kViolet+5);
    datCompPlot[1]->Draw("hist same");
    datCompPlot[1]->Scale(definiteIntegral(opt_pt_hist[1], minpt_comparison_threshold)/
                          definiteIntegral(datCompPlot[1], minpt_comparison_threshold));
    datCompPlot[1]->Draw("hist same");
    opt_pt_hist[1]->SetLineColor(kRed+2);
    opt_pt_hist[1]->Draw("hist same");
    thresholdLine->DrawLine(minpt_comparison_threshold, 0, minpt_comparison_threshold,
                            pp_reference->GetMaximum()*1.1);
    baseLeg_top = new TLegend(0.5, 0.68, 0.88, 0.89);
    baseLeg_top->SetLineWidth(0);
    baseLeg_top->AddEntry(pp_reference, "#bf{pp} Reference p_{T} Spectrum ", "l");
    baseLeg_top->AddEntry(datCompPlot[1], Form("Rescaled #bf{%s} p_{T} Data", species), "l");
    baseLeg_top->AddEntry(opt_pt_hist[1], Form("Fit Result #Deltap_{T} = %.2f", cmin_dpt[1]), "l");
    baseLeg_top->AddEntry(thresholdLine, "p_{T} Comparison Threshold", "l");
    baseLeg_top->Draw();
//    TMultiGraph* bothCDF_rebin = new TMultiGraph();
//    TGraph* cdf_ref = drawCDF(data_raa[iCent], 0, true);
//    TGraph* cdf_rebin = drawCDF(opt_raa_hist[1], 0, true);
//    cdf_ref->SetLineColor(kBlack);
//    cdf_rebin->SetLineColor(kRed+2);
//    bothCDF_rebin->Add(cdf_ref);
//    bothCDF_rebin->Add(cdf_rebin);
//    bothCDF_rebin->Draw("al");
    
    mainpad->cd(6);
    if (dpt < 0.05) ks_r_hist->Draw("hist l");
    else            ks_r_hist->Draw("hist");
    optLine->DrawLine(cmin_dpt[1], 0, cmin_dpt[1], ks_r_hist->GetMaximum()*1.05);
    if (cmin_dpt[1] > max_dpt/2) drawText(Form("#Deltap_{T} = %.2f", cmin_dpt[1]), cmin_dpt[1]-0.1,
                                          ks_r_hist->GetMaximum()*0.95, true, kBlack, 0.04, 42, false);
    else                         drawText(Form("#Deltap_{T} = %.2f", cmin_dpt[1]), cmin_dpt[1]+0.1,
                                          ks_r_hist->GetMaximum()*0.95, false, kBlack, 0.04, 42, false);
    
    mainpad->cd(7);
    gPad->SetLogy();
    data_raa[iCent]->Draw("hist");
    opt_raa_hist[1]->SetLineColor(kRed+2);
    opt_raa_hist[1]->Draw("hist same");
    thresholdLine->DrawLine(minpt_comparison_threshold, data_raa[iCent]->GetMinimum(),
                            minpt_comparison_threshold, data_raa[iCent]->GetMaximum()*1.05);
    baseLeg_low = new TLegend(0.4, 0.15, 0.88, 0.35);
    baseLeg_low->SetLineWidth(0);
    baseLeg_low->SetFillColorAlpha(kWhite, 0.85);
    baseLeg_low->AddEntry(data_raa[iCent], Form("#bf{%s} R_{AA} Data", species), "l");
    baseLeg_low->AddEntry(opt_raa_hist[1], Form("Ratio #frac{#bf{pp} Reference}{Fit Result #Deltap_{T} = %.2f}", cmin_dpt[1]), "l");
    baseLeg_low->Draw();
    
    mainpad->cd(8);
    nbins_doubleratio = refaxis->GetNbins() - refaxis->FindBin(minpt_comparison_threshold+cmin_dpt[1])+1;
    double binedge_ratio_2[nbins_doubleratio+1];
    for (int irbin = 0; irbin <= nbins_doubleratio; ++irbin)
      binedge_ratio_2[irbin] = refaxis->GetBinLowEdge(refaxis->GetNbins() - nbins_doubleratio + irbin + 1);
    doubleratio = new TH1D("",";p_{T};Double Ratio (Data / Shifted)",nbins_doubleratio, binedge_ratio_2);
    for (int irbin = 1; irbin <= nbins_doubleratio; ++irbin) {
      raa =  data_raa[iCent]->GetBinContent(data_raa[iCent]->FindBin(doubleratio->GetBinCenter(irbin)));
      praa = opt_raa_hist[1]->GetBinContent(opt_raa_hist[1]->FindBin(doubleratio->GetBinCenter(irbin)));
      std::cout << "bin " << doubleratio->GetBinCenter(irbin) << ", \tRAA = " << raa << ", \tpRAA = " << praa << std::endl;
      if (raa == 0 || praa == 0) continue;
      doubleratio->SetBinContent(irbin, raa/praa);
    }setStyleLine(doubleratio, "violet thin");
    doubleratio->GetYaxis()->SetRangeUser(0, 2.1);
    doubleratio->Draw("hist");
    drawUnityLine(doubleratio->GetXaxis());
    
    mainpad->cd(9);
    gPad->SetLogy();
    pp_reference->Draw("hist");
    datCompPlot[2] = static_cast<TH1F*>(data_pt[iCent]->Clone());
    datCompPlot[2]->SetLineColor(kViolet+5);
    datCompPlot[2]->Draw("hist same");
    datCompPlot[2]->Scale(definiteIntegral(opt_pt_hist[2], minpt_comparison_threshold)/
                          definiteIntegral(datCompPlot[2], minpt_comparison_threshold));
    datCompPlot[2]->Draw("hist same");
    opt_pt_hist[2]->SetLineColor(kRed+2);
    opt_pt_hist[2]->Draw("hist same");
    thresholdLine->DrawLine(minpt_comparison_threshold, 0, minpt_comparison_threshold,
                            pp_reference->GetMaximum()*1.1);
    baseLeg_top = new TLegend(0.5, 0.68, 0.88, 0.89);
    baseLeg_top->SetLineWidth(0);
    baseLeg_top->AddEntry(pp_reference, "#bf{pp} Reference p_{T} Spectrum ", "l");
    baseLeg_top->AddEntry(datCompPlot[2], Form("Rescaled #bf{%s} p_{T} Data", species), "l");
    baseLeg_top->AddEntry(opt_pt_hist[2], Form("Fit Result #Deltap_{T} = %.2f", cmin_dpt[2]), "l");
    baseLeg_top->AddEntry(thresholdLine, "p_{T} Comparison Threshold", "l");
    baseLeg_top->Draw();
//    TMultiGraph* bothCDF_norebin = new TMultiGraph();
//    TGraph* cdf_ref_pT = drawCDF(data_pt[iCent], 0, true);
//    TGraph* cdf_rebin_pT = drawCDF(pp_reference, cmin_dpt[2], true);
//    cdf_ref_pT->SetLineColor(kBlack);
//    cdf_rebin_pT->SetLineColor(kRed+2);
//    bothCDF_norebin->Add(cdf_ref);
//    bothCDF_norebin->Add(cdf_rebin);
//    bothCDF_norebin->Draw("al");
    
    mainpad->cd(10);
    if (dpt < 0.05) ks_n_hist->Draw("hist l");
    else            ks_n_hist->Draw("hist");
    optLine->DrawLine(cmin_dpt[2], 0, cmin_dpt[2], ks_n_hist->GetMaximum()*1.05);
    if (cmin_dpt[2] > max_dpt/2) drawText(Form("#Deltap_{T} = %.2f", cmin_dpt[2]), cmin_dpt[2]-0.1,
                                          ks_n_hist->GetMaximum()*0.95, true, kBlack, 0.04, 42, false);
    else                         drawText(Form("#Deltap_{T} = %.2f", cmin_dpt[2]), cmin_dpt[2]+0.1,
                                          ks_n_hist->GetMaximum()*0.95, false, kBlack, 0.04, 42, false);
    
    mainpad->cd(11);
    gPad->SetLogy();
    data_raa[iCent]->Draw("hist");
    opt_raa_hist[2]->SetLineColor(kRed+2);
    opt_raa_hist[2]->Draw("hist same");
    thresholdLine->DrawLine(minpt_comparison_threshold, data_raa[iCent]->GetMinimum(),
                            minpt_comparison_threshold, data_raa[iCent]->GetMaximum()*1.05);
    baseLeg_low = new TLegend(0.4, 0.15, 0.88, 0.35);
    baseLeg_low->SetLineWidth(0);
    baseLeg_low->SetFillColorAlpha(kWhite, 0.85);
    baseLeg_low->AddEntry(data_raa[iCent], Form("#bf{%s} R_{AA} Data", species), "l");
    baseLeg_low->AddEntry(opt_raa_hist[2], Form("Ratio #frac{#bf{pp} Reference}{Fit Result #Deltap_{T} = %.2f}", cmin_dpt[2]), "l");
    baseLeg_low->Draw();
    
    mainpad->cd(12);
    nbins_doubleratio = refaxis->GetNbins() - refaxis->FindBin(minpt_comparison_threshold+cmin_dpt[2])+1;
    double binedge_ratio_3[nbins_doubleratio+1];
    for (int irbin = 0; irbin <= nbins_doubleratio; ++irbin)
      binedge_ratio_3[irbin] = refaxis->GetBinLowEdge(refaxis->GetNbins() - nbins_doubleratio + irbin + 1);
    doubleratio = new TH1D("",";p_{T};Double Ratio (Data / Shifted)",nbins_doubleratio, binedge_ratio_3);
    for (int irbin = 1; irbin <= nbins_doubleratio; ++irbin) {
      raa =  data_raa[iCent]->GetBinContent(data_raa[iCent]->FindBin(doubleratio->GetBinCenter(irbin)));
      praa = opt_raa_hist[2]->GetBinContent(opt_raa_hist[2]->FindBin(doubleratio->GetBinCenter(irbin)));
      std::cout << "bin " << doubleratio->GetBinCenter(irbin) << ", \tRAA = " << raa << ", \tpRAA = " << praa << std::endl;
      if (raa == 0 || praa == 0) continue;
      doubleratio->SetBinContent(irbin, raa/praa);
    }setStyleLine(doubleratio, "violet thin");
    doubleratio->GetYaxis()->SetRangeUser(0, 2.1);
    doubleratio->Draw("hist");
    drawUnityLine(doubleratio->GetXaxis());
    
    //--------------------------------------------------------------Export canvas, write data to tree
    
    canvas->SaveAs(Form("../plots/%s_%.2fTeV/%s/fitting/dptFit_cent%i-%i%%.pdf",
                        species, sqrt_s, experiment, centrality_list[iCent][0], centrality_list[iCent][1]));
    tex_title->Clear();
    tex_subtitle->Clear();
    
    energyloss = cmin_dpt[0];
    botcentbin_local = centrality_list[iCent][0];
    topcentbin_local = centrality_list[iCent][1];
    energyloss_tree->Fill();
  }energyloss_tree->Write("energyloss_tree", TObject::kOverwrite);
  return;
}
