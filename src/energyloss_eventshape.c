// A simple root macro for compiling output data from eventshape_summary.c
// and energyloss_pTspectra.c into a single plot.

#include "../config.h"
#include "../utils/root_draw_tools.h"
#include "../utils/hist_tools.h"

void energyloss_eventshape() {
  char species[5];
  snprintf(species, 5, "%s%s", speciesA, speciesB);
  
  // Read data from file.
  char infilename[100];
  snprintf(infilename, 100, "../data.nosync/%s_%.2fTeV/eventshape_energyloss.root", species, sqrt_s);
  TFile* datfile = new TFile(infilename);
  if (datfile->IsZombie()) {
    std::cout << "Error in energyloss_eventshape.c:" << std::endl;
    std::cout << Form("Reference file (%s) does not exist!", infilename) << std::endl;
    std::cout << "Check specfied parameters, or run energyloss_pTspectra.c and eventshape_summary.c to generate this file." << std::endl;
    return;
  }
  
  // Setup tree readers
  TTreeReader* eventshape_reader = new TTreeReader("eventshape_tree", datfile);
  TTreeReaderValue<Double_t> shape_cent_low(*eventshape_reader, "centralitybin_low");
  TTreeReaderValue<Double_t> shape_cent_high(*eventshape_reader, "centralitybin_high");
  TTreeReaderValue<Double_t> eventshape(*eventshape_reader, "eventshape_area");
  
  TTreeReader* energyloss_reader = new TTreeReader("energyloss_tree", datfile);
  TTreeReaderValue<Double_t> energy_cent_low(*energyloss_reader, "centralitybin_low");
  TTreeReaderValue<Double_t> energy_cent_high(*energyloss_reader, "centralitybin_high");
  TTreeReaderValue<Double_t> energyloss(*energyloss_reader, "energyloss");
  
  // Scan trees and gather data
  std::vector<double> eventshape_data;
  std::vector<double> energyloss_data;
  std::vector<double> meancent_reference;
  int count_matches = 0;
  while (eventshape_reader->Next()) {
    while (energyloss_reader->Next()) {
      if (*shape_cent_low != *energy_cent_low) continue;
      if (*shape_cent_high != *energy_cent_high) continue;
      eventshape_data.push_back(*eventshape);
      energyloss_data.push_back(*energyloss);
      meancent_reference.push_back((*shape_cent_low + *shape_cent_high)/2.);
      ++count_matches;
    }energyloss_reader->Restart();
  }
  
  // Sort data into arrays for plotting
  const int ntotal_matches = count_matches;
  double data_eventshape[ntotal_matches];
  double data_energyloss[ntotal_matches];
  for (int i = 0; i < ntotal_matches; ++i) {
    double cmax = 0;
    int grab_index;
    for (int j = 0; j < ntotal_matches; ++j)
      if (meancent_reference.at(j) > cmax) {
        cmax = meancent_reference.at(j);
        grab_index = j;
      }
    data_eventshape[i] = eventshape_data.at(grab_index);
    data_energyloss[i] = energyloss_data.at(grab_index);
    meancent_reference.at(grab_index) = 0;
  }
  
  // Plot and print canvas
  TGraph* result = new TGraph(ntotal_matches, data_eventshape, data_energyloss);
  result->SetTitle(";Averaged Event Area [fm^{2}];Energy Loss [GeV]");
  result->SetMarkerStyle(24);
  result->SetMarkerColor(kBlack);
  result->SetLineColor(kBlack);
  
  
  TCanvas *canvas = new TCanvas();
  result->Draw("alp");
  canvas->SaveAs(Form("../plots/%s_%.2fTeV/eventshape_energyloss.pdf", species, sqrt_s));
}
