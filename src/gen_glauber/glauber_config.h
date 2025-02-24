/// Config parameters used by all code in local glauber processing
/// Calculation time depends primarily on the grid resolution: O(n^2).
/// For this reason, using as coarse a grid as allowable is highly recommended.

// Kinematic settings
const char speciesA[7] = "Xe2arw";               // Collision species for projectile nucleus (A)
const char speciesB[7] = "Xe2arw";               // Collision species for target nucleus (B)
const Double_t sqrt_s = 5.44;                    // Nucleon-nucleon CM collision energy in TeV

// List of centrality bins to use
const Int_t centrality_list[11][2] = {
  {0, 5}, {5, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 50},
  {50, 60}, {60, 70}, {70, 80}, {80, 90}, {90, 100}
};

// An imporant setting:
//  true: code will use pre-published centrality impact parameter binnings from arXiv:1710.07098v3
//   	- recommended for small numbers of events (less than 1M)
//   	- only works for centrality bins in 5% intervals, and for published systems
//   	- can possibly introduce some systematic error due to sig figs in published b bins
// false: code will evenly divide generated events into centrality classes based on b
// 	- asserts centrality divisions by slicing the generated impact parameter distribution
// 	- recommended for large numbers of events (1M or larger)
// 	- works for any system 
// 	- works for any centrality binning (still assumes integral centrality but this is not necessary in principle)
// 	- avoids error in reported significant digits
const bool use_tabulated_centrality = false;


// Grid settings
const Int_t ogrid_resolution = 200;         // Number of points to sample in each dimension
const Double_t grid_edge = 7.5;              // max grid range from nucleon origins in fm

// Settings for cross sections/interaction threshold.
const Double_t interaction_sigmas = 2.5;     // Interaction cross section to include in geometric grid calculation (less is faster but more error, 2.5 recommended as a sweet spot)
const Double_t sigs = 0.4;                   // Nucleon mass profile stdev/cross section
