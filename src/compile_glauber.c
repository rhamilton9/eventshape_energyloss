// Root macro for compiling Monte Carlo output from TGlauberMC
// Principally averages energy density to produce 2D histogram of average event.
// Expected to run over a determined centrality range

#include "../config.h"
#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"
#include "../utils/glauber_tools.h"

void compile_glauber() {
  const int maxNucleons = getMaxNucleons((char*) speciesA, (char*) speciesB);
  const int nbin_1d = maxNucleons/2;
  char species[5];
  snprintf(species, 5, "%s%s", speciesA, speciesB);
  const double bmin = getImpactParameterFromCentralityClass(centralitybin_low,  sqrt_s, (char*) speciesA, (char*) speciesB);
  const double bmax = getImpactParameterFromCentralityClass(centralitybin_high, sqrt_s, (char*) speciesA, (char*) speciesB);
  
  // Attempt to open glauber data file with correct b-range
  char refstring[20];
  snprintf(refstring, 20, "b%.2f-%.2f",bmin,bmax);
  const char* infilename = Form("../data.nosync/%s_%.2fTeV/glauber/glauber_withgrid_%s.root", species, sqrt_s, refstring);
  TFile *fin = new TFile(infilename);
  if (fin->IsZombie()) {
    std::cout << "File with given specifications does not exist." << std::endl;
    std::cout << "Check specfied parameters, or run TGlauberMC generator and try again." << std::endl;
    return;
  }TFile *fout = new TFile(Form("../data.nosync/%s_%.2fTeV/glauber/glauber_compiled_%s.root", species, sqrt_s, refstring), "update");
  
  // Get key data from file
  const int nEvent = ( (TTree*) fin->Get("lemon"))->GetBranch("npart")->GetEntries();
  TH2D* fullsum_noalign[2];
  TH2D* fullsum_cmalign[2];
  TH2D* fullsum_aligned[2];
  fullsum_aligned[0] = (TH2D*) fin->Get("inited_event0")->Clone();
  fullsum_aligned[0]->GetXaxis()->SetAxisColor(kWhite);
  fullsum_aligned[0]->GetYaxis()->SetAxisColor(kWhite);
  fullsum_aligned[0]->GetXaxis()->SetLabelSize(0.04);
  fullsum_aligned[0]->GetYaxis()->SetLabelSize(0.04);
  fullsum_aligned[0]->GetXaxis()->SetTitle("x [fm]");
  fullsum_aligned[0]->GetYaxis()->SetTitle("y [fm]");
  fullsum_aligned[0]->GetXaxis()->SetTitleSize(0.05);
  fullsum_aligned[0]->GetYaxis()->SetTitleSize(0.05);
  fullsum_aligned[0]->GetXaxis()->SetTitleOffset(0.8);
  fullsum_aligned[0]->GetYaxis()->SetTitleOffset(0.8);
  fullsum_aligned[1] = (TH2D*) fullsum_aligned[0]->Clone();
  
  // Tree reader to get event-level data
  TTreeReader* reader = new TTreeReader("lemon", fin);
  TTreeReaderValue<Int_t> npart(*reader, "npart");
  TTreeReaderValue<Int_t> ncoll(*reader, "ncoll");
  TTreeReaderValue<Float_t> b(*reader, "b");
  
  // Histograms for event information
  TH1D* hist_npart = new TH1D(Form("%s_npart_hist",refstring),
                             ";N_{part};#frac{1}{N_{event}} #frac{dN_{event}}{dN_{part}}",nbin_1d, 2, maxNucleons+2);
  TH1D* hist_ncoll = new TH1D(Form("%s_ncoll_hist", refstring),
                             ";N_{coll};#frac{1}{N_{event}} #frac{dN_{event}}{dN_{coll}}",nbin_1d, 0, maxNucleons*6);
  TH1D* hist_b     = new TH1D(Form("%s_b_hist", refstring),
                             ";b [fm];#frac{1}{N_{event}} #frac{dN_{event}}{dNb}",nbin_1d, 0, 20);
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kGreyYellow);
  Int_t color_bg = TColor::GetPalette().At(0);
  TCanvas *canvas = new TCanvas();
  canvas->SetCanvasSize(1000, 600);
  canvas->SetMargin(0, 0, 0, 0);
  TPad* pads[2][3];
  double dx = 0.1*0.3/3;
  pads[0][0] = buildPad("pad0_arit", 0, .475, .9/3+2*dx, .95,          0.1, 0, 0, 0.1, false);
  pads[0][1] = buildPad("pad1_arit", .9/3+2*dx, .475, 1.8/3+dx, .95,   0, 0, 0, 0.1, false);
  pads[0][2] = buildPad("pad2_arit", 1.8/3+dx, .475, .9, .95,        0, 0, 0, 0.1, false);
  pads[1][0] = buildPad("pad0_geom", 0, 0, .9/3+2*dx, .475,           0.1, 0, 0.1, 0, false);
  pads[1][1] = buildPad("pad1_geom", .9/3+2*dx, 0, 1.8/3+dx, .475,    0, 0, 0.1, 0, false);
  pads[1][2] = buildPad("pad2_geom", 1.8/3+dx, 0, .9, .475,         0, 0, 0.1, 0, false);
  pads[0][0]->SetFrameFillColor(color_bg);
  pads[0][1]->SetFrameFillColor(color_bg);
  pads[0][2]->SetFrameFillColor(color_bg);
  pads[1][0]->SetFrameFillColor(color_bg);
  pads[1][1]->SetFrameFillColor(color_bg);
  pads[1][2]->SetFrameFillColor(color_bg);
  
  TString ebecheck_file = Form("../plots/%s_%.2fTeV/glauber/centbin_%i-%i/ebecheck_%s.pdf",
                               species, sqrt_s, centralitybin_low, centralitybin_high, refstring);
  canvas->Print(ebecheck_file + "[");
  
  // clean histograms
  const char scaling_str[2][5] = {"arit", "geom"};
  for (int iScaling = 0; iScaling < 2; ++iScaling) {
    fullsum_aligned[iScaling]->Reset();
    fullsum_noalign[iScaling] = (TH2D*) fullsum_aligned[iScaling]->Clone();
    fullsum_cmalign[iScaling] = (TH2D*) fullsum_aligned[iScaling]->Clone();
    
    fullsum_noalign[iScaling]->SetName(Form("%s_eventshape_noalign_%s",refstring, scaling_str[iScaling]));
    fullsum_cmalign[iScaling]->SetName(Form("%s_eventshape_cmalign_%s",refstring, scaling_str[iScaling]));
    fullsum_aligned[iScaling]->SetName(Form("%s_eventshape_aligned_%s",refstring, scaling_str[iScaling]));
  }
  
  TH2D* chist;
  double xtext[3] = {0.15, 0.05, 0.05};
  double ytext[2] = {0.82, 0.92};
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    // Loop over MC events to compile average event shape
    bool doDraw = iEvent %(nEvent/10) == 0 || iEvent == nEvent - 1;
  
    // Fill 1D histograms from tree data
    reader->Next();
    hist_npart->Fill(*npart);
    hist_ncoll->Fill(*ncoll);
    hist_b->Fill(*b);
    
    // Loop over both scaling modes
    for (int iScaling = 0; iScaling < 2; ++iScaling) {
      if (!iScaling) chist = (TH2D*) fin->Get(Form("inited_event%i", iEvent));
      else           chist = geometricEnergyDensity((TH2D*) fin->Get(Form("initA_event%i", iEvent)),
                                                    (TH2D*) fin->Get(Form("initB_event%i", iEvent)) );
      
      // Draw current event histogram to event-by-event check plotfile
      std::vector<double> cm;
      double planeAngle = getParticipantPlaneAngle(chist)-TMath::Pi()/4;
      if (doDraw) {
        pads[iScaling][0]->cd();
        pads[iScaling][0]->SetTicks(1, 1);
        chist->GetXaxis()->SetAxisColor(kWhite);
        chist->GetYaxis()->SetAxisColor(kWhite);
        chist->GetXaxis()->SetTitle("x [fm]");
        chist->GetYaxis()->SetTitle("y [fm]");
        chist->GetXaxis()->SetLabelSize(0.04);
        chist->GetYaxis()->SetLabelSize(0.04);
        chist->GetXaxis()->SetTitleSize(0.05);
        chist->GetYaxis()->SetTitleSize(0.05);
        chist->GetXaxis()->SetTitleOffset(0.8);
        chist->GetYaxis()->SetTitleOffset(0.8);
        chist->Draw("col");
        cm = getCM(chist);
        drawText(Form("Event %i: Unmodified", iEvent), xtext[0], ytext[iScaling], false, kWhite);
        drawText(Form("CM :: (%.2f, %.2f)", cm.at(0), cm.at(1)), xtext[0], ytext[iScaling]-0.04, false, kWhite);
        drawText(Form("#Psi_{2} :: %.4f", planeAngle), xtext[0], ytext[iScaling]-0.08, false, kWhite);
        if (!iScaling) {
          drawText("Scaling: T_{A}+T_{B}", 0.95, ytext[iScaling], true, kWhite);
          drawText("(Arithmetic)", 0.95, ytext[iScaling]-0.05, true, kWhite);
        } else {
          drawText("Scaling: #sqrt{T_{A}T_{B}}", 0.95, ytext[iScaling], true, kWhite);
          drawText("(Geometric)", 0.95, ytext[iScaling]-0.05, true, kWhite);
        }
      }
      
      // Align current hist and add to average
      fullsum_noalign[iScaling]->Add(fullsum_noalign[iScaling], chist, 1, 1./nEvent);
      chist = translateHist_simple(chist, 0, 0, true);
      fullsum_cmalign[iScaling]->Add(fullsum_cmalign[iScaling], chist, 1, 1./nEvent);
      chist = rotateHist2D_simple(chist, -planeAngle, true);
      fullsum_aligned[iScaling]->Add(fullsum_aligned[iScaling], chist, 1, 1./nEvent);
      
      // Do something with lostdata_translate, lostdata_rotate?
      
      
      // Draw current aligned histogram to event-by-event check plotfile
      // Draws only fully aligned hist, as this allows checks of both cm and psi_2 alignment.
      if (doDraw) {
        pads[iScaling][1]->cd();
        pads[iScaling][1]->SetTicks(1, 1);
        chist->Draw("col");
        cm = getCM(chist);
        planeAngle = getParticipantPlaneAngle(chist)-TMath::Pi()/4;
        drawText(Form("Event %i: Aligned", iEvent), xtext[1], ytext[iScaling], false, kWhite);
        drawText(Form("CM :: (%.2f, %.2f)", cm.at(0), cm.at(1)), xtext[1], ytext[iScaling]-0.04, false, kWhite);
        drawText(Form("#Psi_{2} :: %.4f", planeAngle), xtext[1], ytext[iScaling]-0.08, false, kWhite);
        if (!iScaling) {
          drawText("Scaling: T_{A}+T_{B}", 0.95, ytext[iScaling], true, kWhite);
          drawText("(Arithmetic)", 0.95, ytext[iScaling]-0.05, true, kWhite);
        } else {
          drawText("Scaling: #sqrt{T_{A}T_{B}}", 0.95, ytext[iScaling], true, kWhite);
          drawText("(Geometric)", 0.95, ytext[iScaling]-0.05, true, kWhite);
        }
        
        pads[iScaling][2]->cd();
        pads[iScaling][2]->SetTicks(1, 1);
        fullsum_aligned[iScaling]->Draw("col");
        cm = getCM(fullsum_aligned[iScaling]);
        planeAngle = getParticipantPlaneAngle(fullsum_aligned[iScaling])-TMath::Pi()/4;
        drawText(Form("Full Average to Event %i", iEvent), xtext[2], ytext[iScaling], false, kWhite);
        drawText(Form("CM :: (%.2f, %.2f)", cm.at(0), cm.at(1)), xtext[2], ytext[iScaling]-0.04, false, kWhite);
        drawText(Form("#Psi_{2} :: %.4f", planeAngle), xtext[2], ytext[iScaling]-0.08, false, kWhite);
        if (!iScaling) {
          drawText("Scaling: T_{A}+T_{B}", 0.95, ytext[iScaling], true, kWhite);
          drawText("(Arithmetic)", 0.95, ytext[iScaling]-0.05, true, kWhite);
        } else {
          drawText("Scaling: #sqrt{T_{A}T_{B}}", 0.95, ytext[iScaling], true, kWhite);
          drawText("(Geometric)", 0.95, ytext[iScaling]-0.05, true, kWhite);
        }
      }
    }if (doDraw) {
      canvas->cd();
      TLatex* headerTex = drawText(Form("#it{TGlauberMC} #bf{%s} event %i", species, iEvent), 0.03, 0.93, false, kBlack, 0.05);
      if (iEvent == 0) {
        drawText(Form("Events Generated with #it{b} #in  [%.2f, %.2f]", bmin, bmax), 0.9, 0.96, true, kBlack, 0.04);
        drawText(Form("Corresponding Centrality Bin: %i-%i%%", centralitybin_low, centralitybin_high), 0.9, 0.92, true, kBlack, 0.04);
        chist->Scale(1./chist->GetMaximum());
        chist->GetZaxis()->SetTitle("Fraction of Local Max Energy Density [a.u.]");
        chist->GetZaxis()->SetNdivisions(6);
        chist->GetZaxis()->SetLabelSize(0.03);
        chist->GetZaxis()->SetTitleOffset(0.85);
        chist->GetZaxis()->SetTitleSize(0.03);
        TPaletteAxis* zscale = new TPaletteAxis(0.91, 0.1*0.475, 0.94, 0.95-0.1*0.475, chist);
        zscale->Draw();

      }canvas->Print(ebecheck_file);
      headerTex->Clear();
    }
  }canvas->Print(ebecheck_file + "]");
  
  // Draw three compiled histograms for each energy scaling
  gStyle->SetPalette(kBlueRedYellow);
  for (int iScaling = 0; iScaling < 2; ++iScaling) {
    std::vector<double> cm;
    double planeAngle;
    pads[iScaling][0]->cd();
    fullsum_noalign[iScaling]->Draw("col");
    cm = getCM(fullsum_noalign[iScaling]);
    planeAngle = getParticipantPlaneAngle(fullsum_noalign[iScaling]) - TMath::Pi()/4;
    drawText("Event Shape (Unmodified)", xtext[0], ytext[iScaling], false, kWhite);
    drawText(Form("CM :: (%.2f, %.2f)", cm.at(0), cm.at(1)), xtext[0], ytext[iScaling]-0.04, false, kWhite);
    drawText(Form("#Psi_{2} :: %.4f", planeAngle), xtext[0], ytext[iScaling]-0.08, false, kWhite);
    if (!iScaling) {
      drawText("Scaling: T_{A}+T_{B}", 0.95, ytext[iScaling], true, kWhite);
      drawText("(Arithmetic)", 0.95, ytext[iScaling]-0.05, true, kWhite);
    } else {
      drawText("Scaling: #sqrt{T_{A}T_{B}}", 0.95, ytext[iScaling], true, kWhite);
      drawText("(Geometric)", 0.95, ytext[iScaling]-0.05, true, kWhite);
    }
    
    pads[iScaling][1]->cd();
    fullsum_cmalign[iScaling]->Draw("col");
    cm = getCM(fullsum_cmalign[iScaling]);
    planeAngle = getParticipantPlaneAngle(fullsum_cmalign[iScaling]) - TMath::Pi()/4;
    drawText("Event Shape (CM Aligned)", xtext[1], ytext[iScaling], false, kWhite);
    drawText(Form("CM :: (%.2f, %.2f)", cm.at(0), cm.at(1)), xtext[1], ytext[iScaling]-0.04, false, kWhite);
    drawText(Form("#Psi_{2} :: %.4f", planeAngle), xtext[1], ytext[iScaling]-0.08, false, kWhite);
    if (!iScaling) {
      drawText("Scaling: T_{A}+T_{B}", 0.95, ytext[iScaling], true, kWhite);
      drawText("(Arithmetic)", 0.95, ytext[iScaling]-0.05, true, kWhite);
    } else {
      drawText("Scaling: #sqrt{T_{A}T_{B}}", 0.95, ytext[iScaling], true, kWhite);
      drawText("(Geometric)", 0.95, ytext[iScaling]-0.05, true, kWhite);
    }
    
    pads[iScaling][2]->cd();
    fullsum_aligned[iScaling]->Draw("col");
    cm = getCM(fullsum_aligned[iScaling]);
    planeAngle = getParticipantPlaneAngle(fullsum_aligned[iScaling]) - TMath::Pi()/4;
    drawText("Event Shape (CM and #Psi_{2} Aligned)", xtext[2], ytext[iScaling], false, kWhite);
    drawText(Form("CM :: (%.2f, %.2f)", cm.at(0), cm.at(1)), xtext[2], ytext[iScaling]-0.04, false, kWhite);
    drawText(Form("#Psi_{2} :: %.4f", planeAngle), xtext[2], ytext[iScaling]-0.08, false, kWhite);
    if (!iScaling) {
      drawText("Scaling: T_{A}+T_{B}", 0.95, ytext[iScaling], true, kWhite);
      drawText("(Arithmetic)", 0.95, ytext[iScaling]-0.05, true, kWhite);
    } else {
      drawText("Scaling: #sqrt{T_{A}T_{B}}", 0.95, ytext[iScaling], true, kWhite);
      drawText("(Geometric)", 0.95, ytext[iScaling]-0.05, true, kWhite);
    }
    
    // Draw main title
    canvas->cd();
    drawText(Form("#it{TGlauberMC} #bf{%s} N^{event} = %i", species, nEvent), 0.03, 0.93, false, kBlack, 0.05);
    
    // Write eventshape histograms to file
    fullsum_noalign[iScaling]->Write(fullsum_noalign[iScaling]->GetName(), TObject::kOverwrite);
    fullsum_cmalign[iScaling]->Write(fullsum_cmalign[iScaling]->GetName(), TObject::kOverwrite);
    fullsum_aligned[iScaling]->Write(fullsum_aligned[iScaling]->GetName(), TObject::kOverwrite);
    
  }canvas->SaveAs(Form("../plots/%s_%.2fTeV/glauber/centbin_%i-%i/eventshape_final_%s.pdf",
                       species, sqrt_s, centralitybin_low, centralitybin_high, refstring));
  
  // Rescale event-level histograms; store these and multiplicity to file and close.
  hist_npart->Scale(1./nEvent, "width");
  hist_ncoll->Scale(1./nEvent, "width");
  hist_b->Scale(1./nEvent, "width");
  hist_npart->Write(hist_npart->GetName(), TObject::kOverwrite);
  hist_ncoll->Write(hist_ncoll->GetName(), TObject::kOverwrite);
  hist_b->Write(hist_b->GetName(), TObject::kOverwrite);
  
  TVectorD gendata_Ntotal(1);
  gendata_Ntotal[0] = nEvent;
  gendata_Ntotal.Write("gendata_Ntotal", TObject::kOverwrite);
  
  fout->Close();
  delete fout;
  return;
}
