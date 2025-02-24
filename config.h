// Configuration parameters inherited by other code in the package.

//------------------------------------------------ Kinematic Settings

const char experiment[10] = "ALICE";       // Experiment string
const char speciesA[3] = "Pb";            // Projectile Nucleus (Smaller nucleus if applicable)
const char speciesB[3] = "Pb";            // Target Nucleus (Larger nucleus if applicable)
const char profile_spec[8] = "";      // Specification of which deformed profile to use (if any); see TGlauberMC reference paper for more info.
const double sqrt_s = 2.76;               // sqrt{s} CM energy, in TeV

const double pTmin_hadron = .150;
// Minimum accepted hadron energy in GeV (due to detector efficiency)
// Note that this depends on the detector and should be:
//    - 150 MeV at ALICE
//    - 200 MeV at STAR


// (Pseudo)Rapidity and Azimuth cuts, should be the detector acceptance
const double max_y = 1;
const double max_eta = 0.7;


// This flag controls the cylinder where jet reconstruction happens.
//    true = rapidity-phi plane, false = pseudorapidity-phi plane.
//    It's more physical to reconstruct in y-phi because y is the lorentz angle
//    However, y is difficult to reconstruct in experiemnt, so it's better to
//    reconstruct in eta/phi when making direct comparisons against experiment.
const bool cylinder_uses_rapidity = false;

//------------------------------------------------ Track Spectra settings


const int dataset = 2013;                 // Two currently implemented datasets: 2013 (2.76 TeV), 2018 (2.76, 5.02 TeV)




//------------------------------------------------ Jet Spectra settings

const bool flag_use_charged_jets = true;

// FastJet3 (jet finder) settings
const double jet_radius = 0.2;                    // Jet Radius -- parameter used by jet clustering algorithm
constexpr char algo_string[20] = "anti-k_{T}";    // String for jet clustering algorithm
constexpr char algo_string_short[5] = "a-kT";     // Shortened version
// Other available algorithms, algo_string, and algo_string_short:
//        kt_algorithm            k_{T}         kT
//        cambridge_algorithm     cambridge     cam
//        antikt_algorithm        anti-k_{T}    a-kT
//        genkt_algorithm         gen-k_{T}[n]  g-kTn     (note n is specified by user)
//        Misc. cone algo       [not implemented or considered in the scope of this macro]

// Jet area cut
const double min_area = .6*3.14159*jet_radius*jet_radius;

// Jet/jet constituents and pT cuts for jet clustering.
const double pTmin_jetcore = 0;         // Minumum pT of the leading hadron in the jet
const double jetcut_minpT_full = 5;    // in GeV
const double jetcut_minpT_chrg = 5;     // in GeV
// Note on general values for these parameters:
//    Generally 10 GeV is a good cut for jets to remove many softer gluon jets
//    Extrapolating to lower pT could be good to use a lower cut (5 GeV)
//    The jetcore cut ensures the jet has at least one piece of hard radiation (not a clump of soft particles)
//    This can be relaxed according to the desired analysis


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

// Number of samples (radial bins) with which to take the average for the max.
// For large nsample, the max algo will converge to the average.
// This is a parameter to slightly smooth the max calculation.
// As small a value as reasonable should be used to assure max != avg
// This parameter should be adjusted depending on the grid resolution.
const int nsample_maxregion = 5;

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
const double contour_z = 0.5;           // Fraction of max at which to take z-contour. Half max is possibly best motivated.
const bool use_statwidth_area = true;   // Uses statistical definition of area based on covariances rather than from extracted edge

