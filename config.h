// Configuration parameters inherited by other code in the package.

//------------------------------------------------ Kinematic Settings

const char experiment[10] = "STAR";       // Experiment string
const char speciesA[3] = "Au";            // Projectile Nucleus (Smaller nucleus if applicable)
const char speciesB[3] = "Au";            // Target Nucleus (Larger nucleus if applicable)
const double sqrt_s = 0.20;               // sqrt{s} CM energy, in TeV
const int dataset = 2013;                 // Two currently implemented datasets: 2013 (2.76 TeV), 2018 (2.76, 5.02 TeV)

//------------------------------------------------ Glauber Settings

// bool for retrieving centrality from tabulated values.
// If set to false, a local centrality reference will be used instead.
// Setting up this local reference must be done manually with runGlauber_local.
// This is necessary for species/energies not included in Ref. arXiv:1710.07098v3
// See ../glauber_tools.h for implemented centrality parameters.
const bool useTabulatedCentrality = true;


// Modify these two settings to control centrality throughout the package!
int centralitybin_low = 70;          // Lower edge of centrality bin on current run
int centralitybin_high = 80;         // Upper edge of centrality bin on current run

// Number of bins in each axis of nuclear thickness histogram.
// Heavily affects computation speed at all levels of Glauber modeling.
// Any value 0-1000 can be chosen.
// Recommended settings:
//  - For initial estimates: 100
//  - For full scale local analysis: 250
//  - For final calculations: 1000 (too slow to be performed locally, must be done on the cluster.)
const int ogrid_resolution = 100;

//------------------------------------------------ Energy Loss Settings

// Algorithm for energy loss extraction
// Currently implemented methods:
//    - CHI2        (chi-square fitting with hist rebinning)
//    - KS-REBIN    (Kalgomorov-Smirnov fitting with rebinning)
//    - KS          (Kalgomorov-Smirnov fitting without rebinning)
const char energyloss_extraction_method[10] = "CHI2";

const bool data_errors_sumbyquadrature = false;
//const bool use_new_ks = false;
//const bool shift_ref_down = false;
//const bool use_pt_spectra = false;

// List of centrality bins corresponding to histograms unpacked from HEP data
// Must be set manually, since this differs between data sets.

// For 2018 dataset
const int centrality_list[9][2] = {
  { 0,  5}, { 5, 10}, {10, 20}, {20, 30}, {30, 40},
  {40, 50}, {50, 60}, {60, 70}, {70, 80}
};

// For 2013 reference
//const int centrality_list[15][2] = {
//  { 0,  5}, { 5, 10}, {10, 20}, {20, 30}, {30, 40},
//  {40, 50}, {50, 60}, {60, 70}, {70, 80}, { 0, 10},
//  { 0, 20}, {20, 40}, {40, 60}, {40, 80}, {60, 80}
//};

const int index_raa_in_file = 16;                   // The first plot index Hist1D_y1_N with an R_AA plot
const double minpt_comparison_threshold = 8;        // Lower bound pT for chi2 comparisons, in GeV
const double dpT_resolution = 0.005;                // Step size for pT shift in energyloss calculation
const double max_dpt = 5;                           // Max pT shift to scan, in GeV

//------------------------------------------------ Event Shape Settings

// Key:
//    0 - No CM align or Psi2 align
//    1 - CM align only
//    2 - CM and Psi2 align (recommended)
const int ebe_alignmode = 2;            // Switch for alignment setting on area calculation

const bool isGeometric = true;          // Toggle for Geometric vs. Arithmetic Energy Density scaling
// Algorithm for Glauber edge extraction
// Currently implemented methods:
//    - radial_mean           (Average radius of energy density distrubution)
//    - gradnorm_mean         (Average radius of energy gradient, i.e. pressure distribution)
//    - gradnorm_max          (Maximal pressure surface (radial contour of largest local descent)
//    - z_contour             (Contour at constant z value, tunable. Half-max is default.)
const char eventshape_edgefind_method[20] = "gradnorm_max";
const int nsample_maxregion = 5;        // Number of samples with which to take the average for the max. For large nsample, max reproduces the average.
const double contour_z = 0.5;           // Fraction of max at which to take z-contour. Half max is possibly best motivated.
const bool use_statwidth_area = true;   // Uses statistical definition of area based on covariances rather than from extracted edge

