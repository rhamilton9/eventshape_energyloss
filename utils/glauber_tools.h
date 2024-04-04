// A collection of small macros for Glauber Monte Carlo Modeling.
// Macros are organized roughly by function.


//========================================================================== Table Retrieval

// Return tabulated nucleon cross section, which depends on collision energy and species.
// Cross section fit reproduced from Ref. arXiv:1710.07098v3:
// https://doi.org/10.48550/arXiv.1710.07098
double getNucleonCrossSection(double sqrts_local) {
  if (sqrts_local == 0.20) return 41.6;
  if (sqrts_local == 2.76) return 61.8;
  if (sqrts_local == 5.02) return 67.6;
  return 29.8 + 0.038 * TMath::Power(TMath::Log(sqrts_local), 2.43);
}

int nucleusLookup(char* species) {
  TString species_t = species;
  // Test for most common types first
  if      (species_t == "Au")   return 197;
  else if (species_t == "Pb")   return 208;
  else if (species_t == "p")    return 1001;
  else if (species_t == "d")    return 1002;
  else if (species_t == "O")    return 16;
  else if (species_t == "Xe")   return 129;
  else if (species_t == "Ar")   return 18040;
  else if (species_t == "Ca")   return 20040;
  else if (species_t == "Cu")   return 63;
  else if (species_t == "U")    return 238;
  else if (species_t == "Ni")   return 58;
  else if (species_t == "t")    return 1003;
  else if (species_t == "H3")   return 1003;
  else if (species_t == "He3")  return 2003;
  else if (species_t == "He4")  return 2004;
  else if (species_t == "C")    return 12;
  else if (species_t == "Al")   return 27;
  else if (species_t == "Si")   return 28;
  else if (species_t == "S")    return 32;
  else if (species_t == "W")    return 186;
  std::cout << "Error in glauber_tools::nucleusLookup: Nucleus not found." << std::endl;
  return 0;
}

// Return the total number of nucleons for a given collision system.
int getMaxNucleons(char* speciesA, char* speciesB) {
  return nucleusLookup(speciesA) % 1000 + nucleusLookup(speciesB) % 1000;
}

// Returns the impact parameter corresponding to a given centrality index.
// Values reproduced from Ref. arXiv:1710.07098v3:
// https://doi.org/10.48550/arXiv.1710.07098
// A little clunky with array handling due to pointer/array troubles. Look into later if time allows.
double getImpactParameterFromCentralityClass(int centrality, double sqrts_local, char* speciesA, char* speciesB) {
  if (centrality % 5 != 0)
    std::cout << "Warning in glauber_tools::getImpactParameterFromCentralityClass: " <<
    "Impact parameter at energy sqrt{s_NN} = " << sqrts_local <<
    " corresponding to centrality " << centrality <<
    " not found in lookup for system " << speciesA << speciesB <<
    ".\nAttempting to replace with nearest known centrality data..." <<std::endl;
  TString speciesA_t = speciesA;
  TString speciesB_t = speciesB;
  bool lookupFlag = false;
  double centralityData[21];
  if (speciesA_t == "Pb" && speciesB_t == "Pb"
      && (TMath::Abs(sqrts_local - 2.76) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 3.47, 4.91, 6.01, 6.94,
      7.76, 8.50, 9.18, 9.81, 10.4,
      11.0, 11.5, 12.0, 12.5, 13.0,
      13.4, 13.9, 14.4, 14.9, 15.6, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Pb" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 5.02) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21]  = {
      0.00, 3.49, 4.93, 6.04, 6.98,
      7.80, 8.55, 9.23, 9.87, 10.5,
      11.0, 11.6, 12.1, 12.6, 13.1,
      13.5, 14.0, 14.4, 14.9, 15.6, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Pb" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 5.5) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 3.50, 4.95, 6.05, 6.99,
      7.81, 8.56, 9.25, 9.89, 10.5,
      11.0, 11.6, 12.1, 12.6, 13.1,
      13.5, 14.0, 14.4, 15.0, 15.6, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Pb" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 10.6) <= 0.5)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 3.51, 4.97, 6.08, 7.02,
      7.85, 8.60, 9.29, 9.93, 10.5,
      11.1, 11.6, 12.2, 12.7, 13.1,
      13.6, 14.1, 14.5, 15.0, 15.7, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Pb" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 39) <= 2)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 3.55, 5.02, 6.15, 7.11,
      7.94, 8.70, 9.40, 10.0, 10.7,
      11.2, 11.8, 12.3, 12.8, 13.3,
      13.8, 14.2, 14.7, 15.2, 15.9, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "p" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 5.02) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 1.82, 2.58, 3.16, 3.65,
      4.08, 4.47, 4.83, 5.16, 5.47,
      5.77, 6.05, 6.32, 6.58, 6.84,
      7.10, 7.36, 7.65, 7.99, 8.49, 14.7 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "p" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 8.16) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 1.84, 2.60, 3.19, 3.68,
      4.11, 4.51, 4.87, 5.20, 5.52,
      5.82, 6.10, 6.38, 6.64, 6.90,
      7.15, 7.42, 7.71, 8.05, 8.55, 14.8 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "p" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 8.8) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 1.84, 2.60, 3.19, 3.68,
      4.12, 4.51, 4.87, 5.21, 5.53,
      5.83, 6.11, 6.38, 6.65, 6.91,
      7.16, 7.43, 7.72, 8.06, 8.56, 14.4 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "p" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 17) <= 0.5)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 1.86, 2.63, 3.23, 3.73,
      4.17, 4.56, 4.93, 5.27, 5.59,
      5.89, 6.18, 6.46, 6.73, 6.98,
      7.24, 7.51, 7.80, 8.14, 8.64, 14.9 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "p" && speciesB_t == "Pb"
             && (TMath::Abs(sqrts_local - 63) <= 1)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 1.91, 2.70, 3.31, 3.81,
      4.26, 4.67, 5.04, 5.39, 5.72,
      6.03, 6.32, 6.60, 6.88, 7.14,
      7.40, 7.67, 7.96, 8.31, 8.80, 14.9 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Xe" && speciesB_t == "Xe"
             && (TMath::Abs(sqrts_local - 5.44) <= 0.05)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 3.01, 4.26, 5.22, 6.02,
      6.73, 7.38, 7.97, 8.52, 9.04,
      9.53, 9.99, 10.4, 10.9, 11.3,
      11.7, 12.1, 12.5, 13.1, 13.8, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Au" && speciesB_t == "Au"
             && (TMath::Abs(sqrts_local - 0.2) <= 0.01)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 3.31, 4.68, 5.73, 6.61,
      7.39, 8.10, 8.75, 9.35, 9.92,
      10.5, 11.0, 11.5, 11.9, 12.4,
      12.8, 13.2, 13.7, 14.2, 14.9, 20.0 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  } else if (speciesA_t == "Cu" && speciesB_t == "Cu"
             && (TMath::Abs(sqrts_local - 0.2) <= 0.01)) {
    lookupFlag = true;
    double tempArray[21] = {
      0.00, 2.34, 3.31, 4.06, 4.68,
      5.24, 5.73, 6.19, 6.62, 7.02,
      7.40, 7.77, 8.11, 8.45, 8.78,
      9.11, 9.47, 9.86, 10.3, 11.0, 19.1 };
    for (int i = 0; i < 21; ++i) centralityData[i] = tempArray[i];
  }
  
  if (lookupFlag) return centralityData[centrality/5];
  std::cout << "Error in glauber_tools::getImpactParameterFromCentralityClass: " <<
    "System " << speciesA << speciesB << " not found in lookup at sqrt{s} = " <<
    sqrt_s << "." << std::endl;
  return 0;
}

//========================================================================== Energy Density Analysis

// Helper method to get 2nd order participant plane angle.
// Takes transverse energy density histogram as input.
double getParticipantPlaneAngle(TH2D* edensity) {
  int nbins_x = edensity->GetXaxis()->GetNbins();
  int nbins_y = edensity->GetYaxis()->GetNbins();
  
  double intUnweighted = 0;
  double intWeightedX  = 0;
  double intWeightedY  = 0;
  double intWeightedX2 = 0;
  double intWeightedY2 = 0;
  double intWeightedXY = 0;
  double cx, cy, ceD;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    cx = edensity->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= nbins_y; ++iy) {
      cy = edensity->GetYaxis()->GetBinCenter(iy);
      ceD = edensity->GetBinContent(ix, iy);
      
      // Add to counters
      intUnweighted += ceD;
      intWeightedX  += ceD*cx;
      intWeightedY  += ceD*cy;
      intWeightedX2 += ceD*cx*cx;
      intWeightedY2 += ceD*cy*cy;
      intWeightedXY += ceD*cx*cy;
    }
  }// end of integration loop

  double avgX  = intWeightedX  / intUnweighted;
  double avgX2 = intWeightedX2 / intUnweighted;
  double avgY  = intWeightedY  / intUnweighted;
  double avgY2 = intWeightedY2 / intUnweighted;
  double avgXY = intWeightedXY / intUnweighted;

  double sigmaX = avgX2 - TMath::Power(avgX, 2);
  double sigmaY = avgY2 - TMath::Power(avgY, 2);
  double sigmaXY = avgXY - avgX*avgY;

  return 0.5*TMath::ATan2(sigmaY - sigmaX, 2*sigmaXY );
}


// Method to produce energy density with Geometric scaling E = sqrt{T_A T_B}
// This scaling is perhaps better motivated than the TGlauberMC default E = T_A + T_B
// Key Ref: Bayesian extraction by Bernhard Et. Al. DOI: 10.1103/PhysRevC.94.024907
// See https://link.aps.org/accepted/10.1103/PhysRevC.94.024907
TH2D* geometricEnergyDensity(TH2D* thicknessA, TH2D* thicknessB) {
  TH2D* edensity = (TH2D*) thicknessA->Clone(); edensity->Reset();
  int nbins_x = edensity->GetXaxis()->GetNbins();
  int nbins_y = edensity->GetYaxis()->GetNbins();
  for (int ix = 1; ix <= nbins_x; ++ix) {
    for (int iy = 1; iy <= nbins_y; ++iy) {
      edensity->SetBinContent(ix, iy, TMath::Sqrt(thicknessA->GetBinContent(ix, iy) *
                                                  thicknessB->GetBinContent(ix, iy)));
    }
  }return edensity;
}

