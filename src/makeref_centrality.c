// Root macro. Produces distributions and data used to determine centrality class.
// Resulting distributions and data are used for reference by other local macros.

#include "../utils/root_draw_tools.h"

// todo: this file
// Polish plots
// Make global reference file for key config (to do at night)
// generate 1e6 events to have high stats on true ref
// implement interfacing with other files to enable access to centrality informaton across files

// Generator Info
const int maxNucleons = 420;
const char *collisionType = (char*) "PbPb";

// Bin counts
const int nBin = maxNucleons/2;
const int nbins_b = 400;

// Centrality cut
const double centA = 0.6;
const double centB = 0.7;

// Aesthetic Controls
Int_t palette = kBlueRedYellow;

void makeref_centrality() {
  // Read TTree from file
  TFile *fin = new TFile("out_data.nosync/glauber_nogrid_b0.00-20.00.root");
  TTreeReader* reader = new TTreeReader("lemon", fin);
  TTreeReaderValue<Int_t> npart(*reader, "npart");
  TTreeReaderValue<Int_t> ncoll(*reader, "ncoll");
  TTreeReaderValue<Float_t> b(*reader, "b");
  const int nEvent = ( (TTree*) fin->Get("lemon"))->GetBranch("npart")->GetEntries();
  
  // Histograms for event information, to be filled from tree
  TH1D* nPartHist = new TH1D("npart_hist",";N_{part};#frac{1}{N_{event}} #frac{dP}{dN_{part}}",nBin, 2, maxNucleons+2);
  TH1D* nCollHist = new TH1D("ncoll_hist",";N_{coll};#frac{1}{N_{event}} #frac{dP}{dN_{coll}}",nBin, 0, maxNucleons*6);
  TH2D* b_npart2  = new TH2D("b_npart_hist",";b [fm];N_{part};#frac{d^{2}P}{dbdN_{part}}", nbins_b, 0, 20,nBin, 2, maxNucleons+2);
  TH2D* b_ncoll2  = new TH2D("b_ncoll_hist",";b [fm];N_{coll};#frac{d^{2}P}{dbdN_{coll}}", nbins_b, 0, 20,nBin, 2, maxNucleons*6);
  
  // Define arrays for data storage
//  Float_t npart_data[nEvent];
//  Float_t ncoll_data[nEvent];
//  Float_t b_data[nEvent];
//  std::cout << nEvent << std::endl;
  Float_t npart_mean[nbins_b];
  Float_t ncoll_mean[nbins_b];
  Float_t bmean_vals[nbins_b];
  
  // Fill histograms and data arrays from tree
  int iReader = 0;
  
  while (reader->Next()) {
    
    nPartHist->Fill(*npart);
    nCollHist->Fill(*ncoll);
    b_npart2->Fill(*b, *npart);
    b_ncoll2->Fill(*b, *ncoll);
//    npart_data[iReader] = *npart;
//    ncoll_data[iReader] = *ncoll;
//    b_data[iReader] = *b;
    ++iReader;
    
  }reader->Restart();
  
  
  
  // Normalize histograms to unity
  nPartHist->Scale(1./nPartHist->Integral());
  nCollHist->Scale(1./nCollHist->Integral());
  
  //----------------------------------------------------------------------------Centrality line determination
  
  // Find N_part and N_coll bins that correspond to centA, centB
  int binA_ncoll, binB_ncoll, binA_npart, binB_npart;
  for (int iBin = nBin; iBin > 0; --iBin) {
    if (nCollHist->Integral(iBin, nBin) >= centA) {binA_ncoll = iBin; break;}
  }for (int iBin = binA_ncoll; iBin > 0; --iBin) {
    if (nCollHist->Integral(iBin, nBin) >= centB) {binB_ncoll = iBin; break;}
  }for (int iBin = nBin; iBin > 0; --iBin) {
    if (nPartHist->Integral(iBin, nBin) >= centA) {binA_npart = iBin; break;}
  }for (int iBin = binA_npart; iBin > 0; --iBin) {
    if (nPartHist->Integral(iBin, nBin) >= centB) {binB_npart = iBin; break;}
  }
  
  // Store and print values for given centrality cut
  int nPart_centA = (int)nPartHist->GetBinLowEdge(binA_npart);
  int nPart_centB = (int)nPartHist->GetBinLowEdge(binB_npart);
  int nColl_centA = (int)nCollHist->GetBinLowEdge(binA_ncoll);
  int nColl_centB = (int)nCollHist->GetBinLowEdge(binB_ncoll);
  cout << "Computed bins/N for Centrality Measure:" << endl;
  cout << Form("Centrality %i%% (N_coll) :: bin %i, \tN_coll %i; \t%i%%: bin %i, \tN_coll %i",
               (int)(100*centA), binA_ncoll, nColl_centA, (int)(100*centB), binB_ncoll, nColl_centB) << endl;
  cout << Form("Centrality %i%% (N_part) :: bin %i, \tN_part %i; \t%i%%: bin %i, \tN_part %i",
               (int)(100*centA), binA_npart, nPart_centA, (int)(100*centB), binB_npart, nPart_centB) << endl;
  
  // Compute standard centrality bin lines:
  const int nCent = 7;
  double centList[nCent] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
  int centBins_ncoll[nCent];
  centBins_ncoll[0] = nBin;
  int centBins_npart[nCent];
  centBins_npart[0] = nBin;
  for (int iCent = 0; iCent < nCent; ++iCent) {
    for (int iBin = centBins_ncoll[iCent]; iBin > 0; --iBin) {
      if (nCollHist->Integral(iBin, nBin) >= centList[iCent]) {
        centBins_ncoll[iCent] = iBin;
        if (iCent != nCent-1) centBins_ncoll[iCent+1] = iBin;
        break;
      }
    }
    for (int iBin = centBins_npart[iCent]; iBin > 0; --iBin) {
      if (nPartHist->Integral(iBin, nBin) >= centList[iCent]) {
        centBins_npart[iCent] = iBin;
        if (iCent != nCent-1) centBins_npart[iCent+1] = iBin;
        break;
      }
    }
  }
  
  //----------------------------------------------------------------------------Computing b<->npart and b<->ncoll conversion function
  
  // Compute mean N_part, N_coll as a function of b
  double b_centA_fromNpart, b_centB_fromNpart;
  double b_centA_fromNcoll, b_centB_fromNcoll;
  int index_lowb[4] = {0, 0, 0, 0};
  double b_lowedge = 0;
  double b_highedge = 20./nbins_b;
  double nContained, csum_ncoll, csum_npart;
  for (int ib = 0; ib < nbins_b; ++ib) {
    nContained = 0;
    csum_ncoll = 0;
    csum_npart = 0;
    bmean_vals[ib] = (b_lowedge+b_highedge)/2;
    // Loop over events and find all events with b in current range.
    while (reader->Next()) {
      if (*b >= b_lowedge && *b < b_highedge) {
        csum_ncoll += *ncoll;
        csum_npart += *npart;
        ++nContained;
      }
    }reader->Restart();
    
    // Compute and store mean if well-defined.
    if (nContained != 0) {
      ncoll_mean[ib] = csum_ncoll / nContained;
      npart_mean[ib] = csum_npart / nContained;
    } else {
      std::cout << "No data contained at b = " << bmean_vals[ib] << std::endl;
    }
    
    // Check for corresponding b in centrality cut edges
    if (TMath::Abs(ncoll_mean[ib] - nColl_centA) <= TMath::Abs(ncoll_mean[index_lowb[0]] - nColl_centA)) {
      b_centA_fromNcoll = b_lowedge;
      index_lowb[0] = ib;
    }else if (TMath::Abs(ncoll_mean[ib] - nColl_centB) <= TMath::Abs(ncoll_mean[index_lowb[1]] - nColl_centB)) {
      b_centB_fromNcoll = b_lowedge;
      index_lowb[1] = ib;
    }if (TMath::Abs(npart_mean[ib] - nPart_centA) <= TMath::Abs(npart_mean[index_lowb[2]] - nPart_centA)) {
      b_centA_fromNpart = b_lowedge;
      index_lowb[2] = ib;
    }else if (TMath::Abs(npart_mean[ib] - nPart_centB) <= TMath::Abs(npart_mean[index_lowb[3]] - nPart_centB)) {
      b_centB_fromNpart = b_lowedge;
      index_lowb[3] = ib;
    }
    
    // Increment b edges for next loop iteration
    b_lowedge += 20./nbins_b;
    b_highedge += 20./nbins_b;
  }cout << "here" << endl;
  
  cout << Form("Centrality %i%% (N_coll) :: b %.2f", (int)(100*centA), b_centA_fromNpart) << endl;
  cout << Form("Centrality %i%% (N_coll) :: b %.2f", (int)(100*centB), b_centB_fromNpart) << endl;
  
  // Make histograms from b-selections to compare against centrality cut (b-cut means smearing in centrality)
  TH1D *b_smear_npart = (TH1D*) nPartHist->Clone();
  b_smear_npart->Reset();
  TH1D *b_smear_ncoll = (TH1D*) nCollHist->Clone();
  b_smear_ncoll->Reset();
  double curb;
  while (reader->Next()) {
    curb = *b;
    if (curb >= b_centA_fromNpart && curb <= b_centB_fromNpart) b_smear_npart->Fill(*npart);
    if (curb >= b_centA_fromNcoll && curb <= b_centB_fromNcoll) b_smear_ncoll->Fill(*ncoll);
  }b_smear_npart->Scale(1./nEvent);
  b_smear_ncoll->Scale(1./nEvent);
  
  //----------------------------------------------------------------------------Plotting
  
  // Canvas and Pad Setup
  TCanvas *canvas = new TCanvas();
  canvas->SetMargin(0.08, 0.04, 0.1, 0.06);
  canvas->SetWindowSize(1500,1000);
  TCanvas *compiled = new TCanvas();
  compiled->SetMargin(0, 0, 0, 0);
  compiled->SetWindowSize(1250, 750);
  TPad* pads[4];
  pads[0] = buildPad("pad_0", 0, 0.5, 0.625, 1, 0.08, 0.04, 0.1, 0.06, false);
  pads[1] = buildPad("pad_1", 0, 0, 0.625, 0.5, 0.08, 0.04, 0.1, 0.06, false);
  pads[2] = buildPad("pad_2", 0.625, 0.5, 1, 1, 0.12, 0.17, 0.1, 0.06, false);
  pads[3] = buildPad("pad_3", 0.625, 0, 1, 0.5, 0.12, 0.17, 0.1, 0.06, false);
  
  // Color and style setup
  gStyle->SetOptStat(0);
  gStyle->SetPalette(palette);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(3);
  Int_t color_npart = TColor::GetPalette().At(TColor::GetPalette().GetSize() * 0.48);
  Int_t color_ncoll = TColor::GetPalette().At(TColor::GetPalette().GetSize() * 0.23);
  
  // Draw N_part histogram
  nPartHist->GetYaxis()->SetTitleOffset(1.0);
  nPartHist->GetXaxis()->SetTitleOffset(1.1);
  nPartHist->SetLineColor(color_npart);
  
  b_smear_npart->SetLineColor(kGray+2);
  b_smear_npart->SetFillStyle(3345);
  b_smear_npart->SetFillColorAlpha(kGray+2, 0.5);
  
  TH1D *centHist_npart = (TH1D*) nPartHist->Clone();
  centHist_npart->Reset();
  for (int iBin = binB_npart; iBin < binA_npart; ++iBin) centHist_npart->Fill(nPartHist->GetBinCenter(iBin), nPartHist->GetBinContent(iBin));
  centHist_npart->SetFillColorAlpha(color_npart, 0.4);
  centHist_npart->SetLineWidth(0);
  
  TLine *centLine = new TLine();
  centLine->SetLineWidth(2);
  centLine->SetLineStyle(2);
  centLine->SetLineColor(kGray+2);
  double x, y;
  
  TLegend *leg = new TLegend(0.78, 0.85, 0.93, 0.93);
  leg->SetLineWidth(0);
  leg->SetFillColor(kWhite);
  
  canvas->cd();
  for (int i = 0; i < 2; ++i) {
    gPad->SetLogy();
    nPartHist->Draw("hist");
    centHist_npart->Draw("hist same");
    b_smear_npart->Draw("hist same");
    nPartHist->Draw("hist same");
    
    for (int iCent = 0; iCent < nCent; ++iCent) {
      x = nPartHist->GetBinLowEdge(centBins_npart[iCent]);
      y = nPartHist->GetBinContent(centBins_npart[iCent]);
      centLine->DrawLine(x, 0, x, 1.5*y);
      drawText(Form("%i%%", (int)(100*centList[iCent])), x, 1.5*y, false, kBlack, 0.04, 42, false);
    }
    
    drawText(Form("%i Glauber Events", nEvent), 0.95, 0.89, true, kBlack, 0.05);
    drawText(Form("%s, #sigma_{NN} = 76.6 mb", collisionType), 0.95, 0.83, true, kBlack, 0.05);
    drawText(Form("%i-%i%% Centrality", (int)(100*centA), (int)(100*centB)), 0.95, 0.77, true, kBlack, 0.05);
    
    if (i == 0) {
      compiled->cd();
      pads[0]->Draw();
      pads[0]->cd();
    }
  }canvas->SaveAs("outplot_centrality/npart_dist.pdf");
  
  // Draw N_coll histogram
  nCollHist->GetYaxis()->SetTitleOffset(1.0);
  nCollHist->GetXaxis()->SetTitleOffset(1.1);
  nCollHist->SetLineColor(color_ncoll);
  
  b_smear_ncoll->SetLineColor(kGray+2);
  b_smear_ncoll->SetFillStyle(3345);
  b_smear_ncoll->SetFillColorAlpha(kGray+2, 0.5);
  
  TH1D *centHist_ncoll = (TH1D*) nCollHist->Clone();
  centHist_ncoll->Reset();
  for (int iBin = binB_ncoll; iBin < binA_ncoll; ++iBin) centHist_ncoll->Fill(nCollHist->GetBinCenter(iBin), nCollHist->GetBinContent(iBin));
  centHist_ncoll->SetFillColorAlpha(color_ncoll, 0.4);
  centHist_ncoll->SetLineWidth(0);
  
  centLine = new TLine();
  centLine->SetLineWidth(2);
  centLine->SetLineStyle(2);
  centLine->SetLineColor(kGray+2);
  
  canvas->cd();
  for (int i = 0; i < 2; ++i) {
    gPad->SetLogy();
    nCollHist->Draw("hist");
    centHist_ncoll->Draw("hist same");
    b_smear_ncoll->Draw("hist same");
    nCollHist->Draw("hist same");
    
    for (int iCent = 0; iCent < nCent; ++iCent) {
      x = nCollHist->GetBinLowEdge(centBins_ncoll[iCent]);
      y = nCollHist->GetBinContent(centBins_ncoll[iCent]);
      centLine->DrawLine(x, 0, x, 1.5*y);
      drawText(Form("%i%%", (int)(100*centList[iCent])), x, 1.5*y, false, kBlack, 0.04, 42, false);
    }
    
    drawText(Form("%i Glauber Events", nEvent), 0.95, 0.89, true, kBlack, 0.05);
    drawText(Form("%s, #sigma_{NN} = 76.6 mb", collisionType), 0.95, 0.83, true, kBlack, 0.05);
    drawText(Form("%i-%i%% Centrality", (int)(100*centA), (int)(100*centB)), 0.95, 0.77, true, kBlack, 0.05);
    
    if (i == 0) {
      compiled->cd();
      pads[1]->Draw();
      pads[1]->cd();
    }
  }canvas->SaveAs("outplot_centrality/ncoll_dist.pdf");
  
  // Make mean b-npart and b-ncoll correlation graphs
  TGraph* nPartMeanPlot = new TGraph(nbins_b, bmean_vals, npart_mean);
  nPartMeanPlot->SetName("npart_mean");
  nPartMeanPlot->SetLineWidth(2);
  nPartMeanPlot->SetLineColor(kRed-9);
  TGraph* nCollMeanPlot = new TGraph(nbins_b, bmean_vals, ncoll_mean);
  nCollMeanPlot->SetName("ncoll_mean");
  nCollMeanPlot->SetLineWidth(2);
  nCollMeanPlot->SetLineColor(kBlue-9);
  
  canvas->cd();
  gPad->SetLogy(0);
  canvas->SetWindowSize(1000,1000);
  canvas->SetMargin(0.12, 0.06, 0.1, 0.06);
  // Draw N_part(b), N_coll(b) density scatter plot if requested. Warning: large file when data N_event is large.
//  if (doScatterPlot) {
//    // Scatter plots of data points
//    TGraph* b_ncoll_scatter = new TGraph(nEvent, b_data, ncoll_data);
//    b_ncoll_scatter->SetName("b_ncoll_scatter");
//    b_ncoll_scatter->SetMarkerStyle(1);
//    b_ncoll_scatter->SetMarkerColorAlpha(color_ncoll, 0.1);
//    TGraph* b_npart_scatter = new TGraph(nEvent, b_data, npart_data);
//    b_npart_scatter->SetName("b_npart_scatter");
//    b_npart_scatter->SetMarkerStyle(1);
//    b_npart_scatter->SetMarkerColorAlpha(color_npart, 0.1);
//    
//    // Combine all relevant graphs and plot
//    TMultiGraph *full = new TMultiGraph();
//    full->SetTitle(";b [fm]; #color[4]{N_{coll}}, #color[2]{N_{part}}");
//    full->Add(b_ncoll_scatter, "p");
//    full->Add(b_npart_scatter, "p");
//    full->Add(nCollMeanPlot, "l");
//    full->Add(nPartMeanPlot, "l");
//    full->GetXaxis()->SetRangeUser(0, 20);
//    full->GetYaxis()->SetRangeUser(0, 5.5*maxNucleons);
//    full->GetYaxis()->SetTitleOffset(1.4);
//    full->Draw("a");
//    
//    leg->AddEntry("ncoll_mean", "<N_{coll}>(b)", "l");
//    leg->AddEntry("npart_mean", "<N_{part}>(b)", "l");
//    leg->Draw();
//    canvas->SaveAs("outplot_centrality/b_npart_scatter.pdf");
//  }
  
  // Code for TH2 instead of scatter plot
  canvas->SetMargin(0.12, 0.17, 0.1, 0.06);
  canvas->SetFrameFillColor(TColor::GetPalette().At(0));
  centLine = new TLine();
  centLine->SetLineColor(kWhite);
  centLine->SetLineWidth(1);
  centLine->SetLineStyle(9);
  
  b_npart2->GetZaxis()->SetTitleOffset(1.4);
  b_npart2->GetXaxis()->SetAxisColor(kWhite);
  b_npart2->GetYaxis()->SetAxisColor(kWhite);
  nPartMeanPlot->SetLineColor(kWhite);
  leg = new TLegend(0.7, 0.85, 0.82, 0.93);
  leg->AddEntry("npart_mean", "#color[0]{<N_{part}>(b)}", "l");
  leg->SetFillColor(TColor::GetPalette().At(0));
  for (int i = 0; i < 2; ++i) {
    b_npart2->Draw("colz");
    if (i != 0) nPartMeanPlot->SetLineWidth(1);
    nPartMeanPlot->Draw("l");
    leg->Draw();
    
    centLine->DrawLine(0, nPart_centA, b_centA_fromNpart, nPart_centA);
    centLine->DrawLine(b_centA_fromNpart, nPart_centA, b_centA_fromNpart, 0);
    centLine->DrawLine(0, nPart_centB, b_centB_fromNpart, nPart_centB);
    centLine->DrawLine(b_centB_fromNpart, nPart_centB, b_centB_fromNpart, 0);
    drawText(Form("#color[0]{N_{part, %i%%} = %i}", (int)(100*centA), nPart_centA),
             b_centA_fromNpart*0.2, nPart_centA+5, false, kBlack, 0.03, 42, false);
    drawText(Form("#color[0]{N_{part, %i%%} = %i}", (int)(100*centB), nPart_centB),
             b_centA_fromNpart*0.2, nPart_centB+5, false, kBlack, 0.03, 42, false);
    drawText(Form("#color[0]{%.2f}", b_centA_fromNpart),
             b_centA_fromNpart-0.5, nPart_centB*0.8, false, kBlack, 0.03, 42, false)->SetTextAngle(270);
    drawText(Form("#color[0]{%.2f}", b_centB_fromNpart),
             b_centB_fromNpart-0.5, nPart_centB*0.8, false, kBlack, 0.03, 42, false)->SetTextAngle(270);
    
    if (i == 0) {
      compiled->cd();
      pads[2]->SetFrameFillColor(TColor::GetPalette().At(0));
      pads[2]->cd();
    }
  }canvas->SaveAs("outplot_centrality/b_npart_dist2d.pdf");
  
  b_ncoll2->GetZaxis()->SetTitleOffset(1.4);
  b_ncoll2->GetXaxis()->SetAxisColor(kWhite);
  b_ncoll2->GetYaxis()->SetAxisColor(kWhite);
  nCollMeanPlot->SetLineColor(kWhite);
  leg = new TLegend(0.7, 0.85, 0.82, 0.93);
  leg->AddEntry("ncoll_mean", "#color[0]{<N_{coll}>(b)}", "l");
  leg->SetFillColor(TColor::GetPalette().At(0));
  canvas->cd();
  for (int i = 0; i < 2; ++i) {
    b_ncoll2->Draw("colz");
    nCollMeanPlot->Draw("l");
    if (i != 0) nCollMeanPlot->SetLineWidth(1);
    leg->Draw();
    
    centLine->DrawLine(0, nColl_centA, b_centA_fromNcoll, nColl_centA);
    centLine->DrawLine(b_centA_fromNcoll, nColl_centA, b_centA_fromNcoll, 0);
    centLine->DrawLine(0, nColl_centB, b_centB_fromNcoll, nColl_centB);
    centLine->DrawLine(b_centB_fromNcoll, nColl_centB, b_centB_fromNcoll, 0);
    drawText(Form("#color[0]{N_{coll, %i%%} = %i}", (int)(100*centA), nColl_centA),
             b_centA_fromNcoll*0.2, nColl_centA+30, false, kBlack, 0.03, 42, false);
    drawText(Form("#color[0]{N_{coll, %i%%} = %i}", (int)(100*centB), nColl_centB),
             b_centA_fromNcoll*0.2, nColl_centB+30, false, kBlack, 0.03, 42, false);
    drawText(Form("#color[0]{%.2f}", b_centA_fromNcoll),
             b_centA_fromNcoll-0.4, nColl_centB*0.8, false, kBlack, 0.02, 42, false)->SetTextAngle(270);
    drawText(Form("#color[0]{%.2f}", b_centB_fromNcoll),
             b_centB_fromNcoll-0.4, nColl_centB*0.8, false, kBlack, 0.02, 42, false)->SetTextAngle(270);
    
    if (i == 0) {
      compiled->cd();
      pads[3]->SetFrameFillColor(TColor::GetPalette().At(0));
      pads[3]->cd();
    }
  }canvas->SaveAs("outplot_centrality/b_ncoll_dist2d.pdf");
  compiled->SaveAs("outplot_centrality/centrality_compiled.pdf");
  
  // May add mean ncoll/mean npart/mean b across centrality class histogram, if time allows
  return;
}
