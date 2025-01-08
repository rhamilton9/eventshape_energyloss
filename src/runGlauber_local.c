// Root macro for rapid and convenient Glauber Monte Carlo event generation
// Interfaces with TGlauberMC, principally written by C. Loizides.
// See https://arxiv.org/abs/1710.07098 for latest documentation.
//
// Stores generated events continuously, and can be terminated at any time without loss.
// Randomizer end state is stored after each generation so that generation can be resumed.
// Generating events to a non-empty file will add to existing events without overwriting.
//
// Handling of the randomizer storage is done in TGlauberMC-3.2/runglauber_v3.2.C.
// This is due to modifications to this code, and is not a function of the original MC.
//
// Written by R. Hamilton on 2023-08-16

#include "../config.h"
#include "../utils/glauber_tools.h"
#include "TVectorT.h"
#include "TGlauberMC-3.2/runglauber_v3.2.C"


void runGlauber_local(int nEvent, 
                      bool doGrid = false,
                      double b_low = -1,
                      double b_high = -1) {
  // Setup generation parameters from ../config.h
  double sigNN = getNucleonCrossSection(sqrt_s);
  char species[5];
  snprintf(species, 5, "%s%s", speciesA, speciesB);
  
  std::cout << "Running at energy sqrt{s} = " << sqrt_s << " TeV" << std::endl;
  // If b_low and b_high are -1 (i.e. default flag) use config-set centrality
  if (b_low == -1 && b_high == -1) {
    if (useTabulatedCentrality) {
      b_low  = getImpactParameterFromCentralityClass(centralitybin_low, sqrt_s, (char*) speciesA, (char*) speciesB);
      b_high = getImpactParameterFromCentralityClass(centralitybin_high, sqrt_s, (char*) speciesA, (char*) speciesB);
      
      // TODO -- Implement i/o with local centrality reference
    } else { b_low = 0; b_high = 20; }
  }
  
  // Access or create relevant files, and check for existing data.
  const char desc[2][10] = {"withgrid","nogrid"};
  char filename[150];
  snprintf(filename, 150, "../data.nosync/%s_%.2fTeV/glauber/glauber_%s_b%.2f-%.2f.root", species, sqrt_s, desc[!doGrid],b_low,b_high);
  // use snprintf(). Some errors occur with pointer handling at 8th iteration when using Form().
  
  // Attempt to open file with the same parameters.
  // If such a file exists, attempt to recover randomizer state and continue generation.
  int infile_nEvent = 0;
  TFile *test_infile = new TFile(filename,"read");
  if (!test_infile->IsZombie()) {
    // Sample file detected! Attempt to get multiplicity data...
    if (!test_infile->Get("lemon")) {
      std::cout << "File with given specifications already exists, but has no lemon tree." << std::endl;
      std::cout << "Check output file at " << filename << std::endl;
      return;
    }
    
    // Multiplicity data found!
    infile_nEvent = ( (TTree*) test_infile->Get("lemon"))->GetBranch("npart")->GetEntries();
    std::cout << Form("File with given specifications already exists with multiplicity :: %i", infile_nEvent) << std::endl;
    std::cout << "Attemping to get randomizer state from this file and continue generation...\n" << std::endl;
    
    // Attempt to find randomizer state from previous run...
    if (!test_infile->Get("randomizer_endstate")) {
      std::cout << "Randomizer state not found. To ensure reproducabilty of results, generation has been terminated." << std::endl;
      return;
    }std::cout << "Success! Proceeding to Glauber Monte-Carlo...\n" << std::endl;
    
    // NOTE:
    // Randomizer handled in runglauber_v3.2.C
    // Sample implementation :: gRandom = (TRandom*) test_infile->Get("randomizer_endstate")->Clone();
  }delete test_infile;
  
  // Run MC Glauber with reference inputs
  // Execute iteratively to store events as generation occurs.
  int recur_step = 1000;
  if (doGrid) recur_step = 100;
  gSystem->Load("libMathMore");
  for (int iRecur = infile_nEvent+recur_step; iRecur < nEvent+infile_nEvent; iRecur += recur_step) {
    // Reference: Method definition and inputs
    //  void runAndOutputLemonTree(const Int_t n,
    //           const Double_t sigs  = 0.4,
    //           const char *sysA     = "p",
    //           const char *sysB     = "Pb",
    //           const Double_t signn = 67.6,
    //           const Double_t mind  = 0.4,
    //           const Double_t bmin  = 0.0,
    //           const Double_t bmax  = 20.0,
    //           const Bool_t   ogrid = 0,
    //           const char *fname    = 0);
    runAndOutputLemonTree(iRecur, 0.4, speciesA, speciesB, sigNN, 0.4, b_low, b_high, doGrid, filename);
    std::cout << Form("%i Events completed on this generation.", iRecur-infile_nEvent) << std::endl;
    std::cout << "Writing to file..." << std::endl;
  }
  // run final set, separated from loop in case input nEvent not a multiple of recur_step.
  runAndOutputLemonTree(nEvent+infile_nEvent, 0.4, speciesA, speciesB, sigNN, 0.4, b_low, b_high, doGrid, filename);
  std::cout << Form("%i Events completed on this generation.", nEvent) << std::endl;
  std::cout << "All events completed." << std::endl;

  return;
}
