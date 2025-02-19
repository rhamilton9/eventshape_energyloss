// A collection of small macros for histogram analysis
// Macros are organized roughly by function

// Variables for storing data that may be of interest in the daughter macro

#include "interpolation_tools.h"

double lostdata_translate;
double lostdata_rotate;

//========================================================================== Search/Retrieval/Optimization

// Returns an integer which signals the type of TH1 inheriting class.
// Used to determine qualities of the histogram when not specified.
int getHistType(TH1* hist) {
  // Test most common types first
  if (TH1D* hist_TH1D = dynamic_cast<TH1D*>(hist))  return 11;   // TH1D
  if (TH1F* hist_TH1F = dynamic_cast<TH1F*>(hist))  return 12;   // TH1F
  
  if (TH2D* hist_TH2D = dynamic_cast<TH2D*>(hist))  return 21;   // TH2D
  if (TH2F* hist_TH2F = dynamic_cast<TH2F*>(hist))  return 22;   // TH2F
  
  if (TH3D* hist_TH3D = dynamic_cast<TH3D*>(hist))  return 31;   // TH3D
  if (TH3F* hist_TH3F = dynamic_cast<TH3F*>(hist))  return 32;   // TH3F
  
  // Test for more obscure types
  if (TH1I* hist_TH1I = dynamic_cast<TH1I*>(hist))  return 13;   // TH1I
  if (TH1C* hist_TH1C = dynamic_cast<TH1C*>(hist))  return 14;   // TH1C
  if (TH1S* hist_TH1S = dynamic_cast<TH1S*>(hist))  return 15;   // TH1S
  if (TH1K* hist_TH1K = dynamic_cast<TH1K*>(hist))  return 16;   // TH1K
  
  if (TH2I* hist_TH2I = dynamic_cast<TH2I*>(hist))  return 23;   // TH2I
  if (TH2C* hist_TH2C = dynamic_cast<TH2C*>(hist))  return 24;   // TH2C
  if (TH2S* hist_TH2S = dynamic_cast<TH2S*>(hist))  return 25;   // TH2S
  
  if (TH3I* hist_TH3I = dynamic_cast<TH3I*>(hist))  return 33;   // TH3I
  if (TH3C* hist_TH3C = dynamic_cast<TH3C*>(hist))  return 34;   // TH3C
  if (TH3S* hist_TH3S = dynamic_cast<TH3S*>(hist))  return 35;   // TH3S
  
  // Not recognized. (could be TProfile, TH2Poly, TGLTH3Composition)
  return 0;
}

// Returns the largest value from two input histograms
double getMaxFromHists(TH1* hist1, TH1* hist2) {
    double max1 = hist1->GetMaximum();
    double max2 = hist2->GetMaximum();
    if (max2 > max1) return max2;
    return max1;
}

// Returns the usual calculus definite integral of a hist within a region
// Treats bin contents as uniformly distributed within the bin
double definiteIntegral(TH1* hist, double lowbound = INT_MIN, double topbound = INT_MAX) {
  TAxis* axis = hist->GetXaxis();
  int nbin = axis->GetNbins();
  
  // Edge cases & signed integral considerations
  int int_sign = 1;
  if (axis->FindBin(lowbound) == axis->FindBin(topbound)) return hist->GetBinContent(axis->FindBin(lowbound)) * (topbound - lowbound);
  if (lowbound > topbound) {double temp = topbound; topbound = lowbound; lowbound = temp; int_sign = -1;}
  if (lowbound == INT_MIN || lowbound < axis->GetBinLowEdge(1))       lowbound = axis->GetBinLowEdge(1);
  if (topbound == INT_MAX || topbound > axis->GetBinLowEdge(nbin+1))  topbound = axis->GetBinLowEdge(nbin+1);
  
  // Compute integral
  int startbin = axis->FindBin(lowbound);
  double integral = hist->GetBinContent(startbin)*(axis->GetBinLowEdge(startbin+1) - lowbound);
  for (int ibin = startbin+1; ibin <= nbin; ++ibin) {
    if (axis->GetBinLowEdge(ibin+1) >= topbound) {
      integral += hist->GetBinContent(ibin) * (topbound - axis->GetBinLowEdge(ibin));
      break;
    }integral += hist->GetBinContent(ibin)*axis->GetBinWidth(ibin);
  }return integral*int_sign;
}

// Returns the coordinates for the center-of-mass of an input hist.
// Output is a std::vector<double> of coordinates.
// Length of output varies with histogram dimension.
// Flagging hasUniformBinning increases computation speed.
std::vector<double> getCM(TH1* hist,
                          int histflag = -1,
                          bool hasUniformBinning = false) {
  // Variable setup
  std::vector<double> cm;
  TAxis* ref_axis[3];
  double binVol = 1;
  
  if (histflag == -1) histflag = getHistType(hist);
  int hist_dim = histflag/10;
  switch (histflag) {
    case 31: case 32: case 33: 
      ref_axis[2] = hist->GetZaxis();
      if (hasUniformBinning) binVol *= ref_axis[2]->GetBinWidth(1);
    case 21: case 22: case 23: 
      ref_axis[1] = hist->GetYaxis();
      if (hasUniformBinning) binVol *= ref_axis[1]->GetBinWidth(1);
    case 11: case 12: case 13:
      ref_axis[0] = hist->GetXaxis();
      if (hasUniformBinning) binVol *= ref_axis[0]->GetBinWidth(1);
      break;
    default:
      cout << "Error in getCM: Input Histogram type not recognized." << endl;
      cout << "Input must be of type THnD, THnF or THnI for n = 1,2,3." << endl;
      return cm;
  }
  
  // Compute CM
  double sum;
  double temp;
  double weighted_sum[3] = {0,0,0};
  for (int iBinX = 1; iBinX <= ref_axis[0]->GetNbins(); ++iBinX) {
    if (hist_dim > 1) {
      for (int iBinY = 1; iBinY <= ref_axis[1]->GetNbins(); ++iBinY) {
        if (hist_dim > 2) {
          for (int iBinZ = 1; iBinZ <= ref_axis[2]->GetNbins(); ++iBinZ) {
            // TH3
            if (hasUniformBinning)
              temp = (hist->GetBinContent(iBinX, iBinY, iBinZ)
                      * binVol);
            else
              temp = (hist->GetBinContent(iBinX, iBinY, iBinZ)
                      * ref_axis[0]->GetBinWidth(iBinX)
                      * ref_axis[1]->GetBinWidth(iBinY)
                      * ref_axis[2]->GetBinWidth(iBinZ));
            sum += temp;
            weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
            weighted_sum[1] += temp * ref_axis[1]->GetBinCenter(iBinY);
            weighted_sum[2] += temp * ref_axis[2]->GetBinCenter(iBinZ);
          }
        } else {
          // TH2
          if (hasUniformBinning)
            temp = (hist->GetBinContent(iBinX, iBinY)
                    * binVol);
          else
            temp = (hist->GetBinContent(iBinX, iBinY)
                    * ref_axis[0]->GetBinWidth(iBinX)
                    * ref_axis[1]->GetBinWidth(iBinY));
          sum += temp;
          weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
          weighted_sum[1] += temp * ref_axis[1]->GetBinCenter(iBinY);
        }
      }
    } else {
      // TH1
      if (hasUniformBinning)
        temp = (hist->GetBinContent(iBinX)
                * binVol);
      else
        temp = (hist->GetBinContent(iBinX)
                * ref_axis[0]->GetBinWidth(iBinX));
      sum += temp;
      weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
    }
  }
  
  // Push to CM vector and return
  for (int iDim = 0; iDim < hist_dim; ++iDim) {
    cm.push_back(weighted_sum[iDim] / sum);
  }return cm;
}

//========================================================================== Histogram Modification

// Translates a histogram's contents by the specified amounts.
// Preserves the original binning of the parent histogram.
// If no amount is specified, the method defaults to the CM.
// Produces warning of any data lost to translation over edge.
//
// STILL UNDER DEVELOPMENT. Working as intended for 1D hist
// Use placeholder macro translateHist_simple() for 2D or 3D hist
// Tests indicate that this method is unable to recover the parent hist.
// Further consideration for this implementation (translate, rebin) should be taken.
TH1* translateHist(TH1* hist,
                   double delta_x = 0,
                   double delta_y = 0,
                   double delta_z = 0,
                   int histflag = -1) {
  // Variable setup
  TAxis* ref_axis[3];
  if (histflag == -1) histflag = getHistType(hist);
  int hist_dim = histflag/10;
  int nbins[3];
  switch (histflag) {
    case 31: case 32: 
      ref_axis[2] = hist->GetZaxis();
      nbins[2] = ref_axis[2]->GetNbins();
    case 21: case 22:
      ref_axis[1] = hist->GetYaxis();
      nbins[1] = ref_axis[1]->GetNbins();
    case 11: case 12:
      ref_axis[0] = hist->GetXaxis();
      nbins[0] = ref_axis[0]->GetNbins();
      break;
    default:
      cout << "Error in translateHist: Input Histogram type not recognized." << endl;
      cout << "Input must be of type THnD or THnF for n = 1,2,3." << endl;
      return new TH1D();
  }
  
  TH1D* hist_1D;
  TH1F* hist_1F;
  TH2D* hist_2D;
  TH2F* hist_2F;
  TH3D* hist_3D;
  TH3F* hist_3F;
  double totaldat;
  switch (histflag) {
    case 11: hist_1D = 
      static_cast<TH1D*>(hist->Clone()); hist_1D->Reset();
//      totaldat = hist_1D->Integrate("width"); 
      break;
    case 12:
      hist_1F = static_cast<TH1F*>(hist->Clone()); hist_1F->Reset();
//      totaldat = hist_1F->Integrate("width"); 
      break;
    case 21:
      hist_2D = static_cast<TH2D*>(hist->Clone()); hist_2D->Reset();
//      totaldat = hist_2D->Integrate("width"); 
      break;
    case 22:
      hist_2F = static_cast<TH2F*>(hist->Clone()); hist_2F->Reset();
//      totaldat = hist_2F->Integrate("width"); 
      break;
    case 31:
      hist_3D = static_cast<TH3D*>(hist->Clone()); hist_3D->Reset();
//      totaldat = hist_3D->Integrate("width"); 
      break;
    case 32:
      hist_3F = static_cast<TH3F*>(hist->Clone()); hist_3F->Reset();
//      totaldat = hist_3F->Integrate("width"); 
      break;
  }
  
  // Set shifts to CM if default input given
  if (delta_x == 0 && delta_y == 0 && delta_z == 0) {
    std::vector<double> cm = getCM(hist);
    delta_x = -cm.at(0);
    if (hist_dim > 1) delta_y = -cm.at(1);
    if (hist_dim > 2) delta_z = -cm.at(2);
  }
  
  // Translate
  int fillbin;
  double lostdata = 0;
  double cbin_data;
  double c_lowedge[3];
  double c_topedge[3];
  for (int iBinX = 1; iBinX <= nbins[0]; ++iBinX) {
    if (hist_dim > 1) {
      for (int iBinY = 1; iBinY <= nbins[1]; ++iBinY) {
        if (hist_dim > 2) {
          for (int iBinZ = 1; iBinZ <= nbins[2]; ++iBinZ) {
            // TH3
            
            // TODO
            
          } // End z loop
        } else {
          // TH2
          
          // TODO
          
        }
      } // End y loop
    } else {
      // TH1
      cbin_data = hist->GetBinContent(iBinX);
      
      c_lowedge[0] = ref_axis[0]->GetBinLowEdge(iBinX) + delta_x;
      c_topedge[0] = ref_axis[0]->GetBinLowEdge(iBinX+1) + delta_x;
      while (c_lowedge[0] < c_topedge[0]) {
        // Test edge cases
        if (c_lowedge[0] < ref_axis[0]->GetBinLowEdge(1)) { // Underflow
          // Check for full bin underflow
//          cout << "found undeflow." << endl;
          if (c_topedge[0] <= ref_axis[0]->GetBinLowEdge(1)) {
//            cout << "found total underflow" << endl;
            lostdata += cbin_data * (c_topedge[0] - c_lowedge[0]);
            break; }
          // Otherwise, only partial underflow
          lostdata += cbin_data * (ref_axis[0]->GetBinLowEdge(1) - c_lowedge[0]);
          c_lowedge[0] = ref_axis[0]->GetBinLowEdge(1);
        } else if (c_topedge[0] > ref_axis[0]->GetBinLowEdge(nbins[0]+1)) { // Overflow
          // Check for full bin overflow
          if (c_lowedge[0] >= ref_axis[0]->GetBinLowEdge(nbins[0]+1)) {
            lostdata += cbin_data * (c_topedge[0] - c_lowedge[0]);
            break; }
          // Otherwise, only partial underflow.
          lostdata += cbin_data * (c_topedge[0] - ref_axis[0]->GetBinLowEdge(nbins[0]+1));
          c_topedge[0] = ref_axis[0]->GetBinLowEdge(nbins[0]+1);
        }
        
        // Fill the current bin
        fillbin = ref_axis[0]->FindBin(c_lowedge[0] + 1e-10); // Extra bit to prevent roundoff error
        if (c_topedge[0] >= ref_axis[0]->GetBinLowEdge(fillbin+1)) { // only partial or exact filling
          switch (histflag) {
            case 11: hist_1D->SetBinContent(fillbin,
                                            hist_1D->GetBinContent(fillbin) +
                                            cbin_data * (ref_axis[0]->GetBinLowEdge(fillbin+1) - c_lowedge[0])
                                            / ref_axis[0]->GetBinWidth(fillbin)); break;
            case 12: hist_1F->SetBinContent(fillbin,
                                            hist_1F->GetBinContent(fillbin) +
                                            cbin_data * (ref_axis[0]->GetBinLowEdge(fillbin+1) - c_lowedge[0])
                                            / ref_axis[0]->GetBinWidth(fillbin)); break;
          }c_lowedge[0] = ref_axis[0]->GetBinLowEdge(fillbin+1);
//          cout << Form("Trapped in loop on iBinX = %i, c_topedge = %.2f, c_lowedge = %.2f", iBinX, c_topedge[0], c_lowedge[0]) << endl;
        } else { // Add all bin content to this bin
          switch (histflag) {
            case 11: hist_1D->SetBinContent(fillbin,
                                            hist_1D->GetBinContent(fillbin) +
                                            cbin_data * (c_topedge[0] - c_lowedge[0])
                                            / ref_axis[0]->GetBinWidth(fillbin)); break;
            case 12: hist_1F->SetBinContent(fillbin,
                                            hist_1F->GetBinContent(fillbin) +
                                            cbin_data * (c_topedge[0] - c_lowedge[0])
                                            / ref_axis[0]->GetBinWidth(fillbin)); break;
          }c_lowedge[0] = c_topedge[0];
//          cout << "Ran the break condition" << endl;
        }
      } 
//      cout << "broke out." << endl;
      // End while (distributing bin contents to translated bins)
    }
  } // End x loop
  
  // Report on lost data
//  if (lostdata > 0)
//    cout << Form("Warning in translateHist: %.2f%% of data lost to over/underflow.",
//                 100*lostdata / hist->Integral("width")) << endl;
  lostdata_translate = lostdata;
  
  // Return translated histogram
  switch (histflag) {
    case 12: return hist_1F;
    case 21: return hist_2D;
    case 22: return hist_1F;
    case 31: return hist_3D;
    case 32: return hist_3F;
    default: return hist_1D;
  }
}

// Rotates a histogram counterclockwise through the specified angle
//
// STILL UNDER DEVELOPMENT. Use placeholder macro rotateHist2D_simple()
TH2D* rotateHist2D(TH1* hist, double angle) {
  return new TH2D();
}



//========================================================================== Statistical Tools

/// Computes the Empirical Cumulative Distribution Function for a 1D hist.
/// -----------------------*IMPORTANT*-----------------------
/// The results will vary enormously with the choice of interpolation mode:
///     1) Linear interpolation:
///       Assumes the contents within a bin are uniformly distributed. This is generally a bad assumption
///     2) Cubic Spline Interpolation
///       Not yet implemented
///     3) Root-default curve interpolation algorithm
///       See Ref. https://academic.oup.com/comjnl/article/13/4/392/540402
/// ---------------------------------------------------------
/// Input must be either TH1D or TH1F type; other types will throw an error.
///
/// The output is a vector of ordered pairs (axis value, CDF) each of which
/// is individually formed as a vector.
std::vector<std::vector<double>> CDF_FromHist_new(TH1* hist,
                                                  double thresh_min = INT_MIN,
                                                  double thresh_max = INT_MAX,
                                                  int interpolation_mode = 1) {
  // Verify that input hist is of compatible type
  int histflag = getHistType(hist);
  TAxis* axis = hist->GetXaxis();
  int nbin = axis->GetNbins();
  std::vector<std::vector<double>> cdf;
  switch (histflag) {
    case 11: case 12: break;
    default:
      std::cout<<"Error in <utils/hist_tools::CDF_FromHist>: Input hist type is not TH1F or TH1D."<<std::endl;
      return cdf;
  }
  
  // Test for bad thresholds, throw errors if necessary,
  if (thresh_min > thresh_max) {
    std::cout<<"Warning in <utils/hist_tools::CDF_FromHist>: Axis range is not well defined (min > max)."<<std::endl;
    std::cout<<"Attempting to switch inputs..."<<std::endl;
    double temp = thresh_min;
    thresh_min = thresh_max;
    thresh_max = temp;
  } if (thresh_min > axis->GetBinLowEdge(nbin+1)) {
    std::cout<<"Error in <utils/hist_tools::CDF_FromHist>: Minimum Threshold is above max histogram range."<<std::endl;
    return cdf;
  } if (thresh_max < axis->GetBinLowEdge(1)) {
    std::cout<<"Error in <utils/hist_tools::CDF_FromHist>: Maximum Threshold is below min histogram range."<<std::endl;
    return cdf;
  }
  
  // Get bin count and integral for counting and rescaling
  double norm = definiteIntegral(hist, thresh_min, thresh_max);
  if (thresh_min == INT_MIN || thresh_min < axis->GetBinLowEdge(1))      thresh_min = axis->GetBinLowEdge(1);
  if (thresh_max == INT_MAX || thresh_max > axis->GetBinLowEdge(nbin+1)) thresh_max = axis->GetBinLowEdge(nbin+1);
  
  // Get first CDF bins and test edge cases
  int startbin = axis->FindBin(thresh_min);
  std::vector<double> cpair = {thresh_min, 0};
  cdf.push_back(cpair);
  if (startbin == axis->FindBin(thresh_max)) {
    cpair = {thresh_max, 1};
    cdf.push_back(cpair);
    return cdf;
  }cpair.at(0) = axis->GetBinLowEdge(startbin+1);
  cpair.at(1) = (cpair.at(0) - thresh_min) * hist->GetBinContent(startbin) / norm;
  cdf.push_back(cpair);
  
  // Compute CDF over remaining bin range
  double cdf_counter = cpair.at(1);
  for (int ibin = startbin+1; ibin <= nbin; ++ibin) {
    if (axis->GetBinLowEdge(ibin+1) >= thresh_max) { // Threshold end case
      cdf_counter += hist->GetBinContent(ibin) * (thresh_max - axis->GetBinLowEdge(ibin)) / norm;
      cpair = {thresh_max, cdf_counter};
      cdf.push_back(cpair);
      break;
    }
    
    cdf_counter += hist->GetBinContent(ibin) * axis->GetBinWidth(ibin) / norm;
    cpair = {axis->GetBinLowEdge(ibin+1), cdf_counter};
    cdf.push_back(cpair);
  }
  
  if (interpolation_mode == 1) return cdf;
  if (cdf.size() <= 2) return cdf;
  
  // Perform an interpolation on the CDF if specified.
  std::vector<double> bin_axis;
  std::vector<double> cdf_axis;
  int npoint = cdf.size();
  for (int i = 0; i < npoint; ++i) {
    bin_axis.push_back(cdf.at(i).at(0));
    cdf_axis.push_back(cdf.at(i).at(1));
  }TGraph* interpolated = Smooth(new TGraph(npoint, bin_axis.data(), cdf_axis.data()), true);
  Double_t* interp_binaxis = interpolated->GetX();
  Double_t* interp_cdfaxis = interpolated->GetY();
  npoint = interpolated->GetN();
  cdf.clear();
  for (int i = 0; i < npoint; ++i) {
    if (interp_cdfaxis[i] > 1 || interp_cdfaxis[i] < 0) continue;
    cpair = {interp_binaxis[i], interp_cdfaxis[i]};
//    std::cout << "bin ax : " << interp_binaxis[i] << "\tcdf : " << interp_cdfaxis[i] << std::endl;
    cdf.push_back(cpair);
  }return cdf;
}

// Draws and returns a graphical representation of the CDF
// Result is drawn on the current gPad unless suppressed.
TGraph* drawCDF_new(TH1* hist,
                    double horiz_shift = 0,
                    double thresh_min = INT_MIN,
                    double thresh_max = INT_MAX,
                    bool drawLogy = false,
                    bool suppressDraw = false) {
  std::vector<std::vector<double>> cdf = CDF_FromHist_new(hist, thresh_min, thresh_max);
  int nentry = cdf.size();
  
  // Set up value arrays
  double bin_axis[nentry];
  double cdf_axis[nentry];
  
  // Loop over CDF and assign axis vals to CDF
  for (int ientry = 0; ientry < nentry; ++ientry) {
    bin_axis[ientry] = cdf.at(ientry).at(0);
    if (drawLogy) cdf_axis[ientry] = 1-cdf.at(ientry).at(1);
    else          cdf_axis[ientry] = cdf.at(ientry).at(1);
  }
  
  // Construct TGraph and plot if specified
  TGraph* cdf_graph = new TGraph(nentry, bin_axis, cdf_axis);
  cdf_graph->SetTitle(Form(";%s;eCDF",hist->GetXaxis()->GetTitle()));
  if (!suppressDraw) {
    if (thresh_min < bin_axis[0])        thresh_min = bin_axis[0];
    if (thresh_max > bin_axis[nentry-1]) thresh_max = bin_axis[nentry-1];
    cdf_graph->Draw("al");
    TLine *rangeLine = new TLine();
    setStyleLine(rangeLine, "thin dashed gray");
    rangeLine->DrawLine(thresh_min, cdf_axis[0], thresh_min, cdf_axis[nentry-1]);
    rangeLine->DrawLine(thresh_max, cdf_axis[0], thresh_max, cdf_axis[nentry-1]);
  }return cdf_graph;
}

// Computes the Kolgomorov-Smirnov (KS) statistic for 2 histograms
//
// The variable horizShiftOnHist1 allows for the CDF of hist1 to be
// shifted horizontally by some amount before computing the KS.
// Only horizontal axis values above comparison_threshold are considered.
double KS_statistic_new(TH1* hist1,
                        TH1* hist2,
                        double horizontal_shift = 0,
                        double thresh_min = INT_MIN,
                        double thresh_max = INT_MAX,
                        bool doPlot = false,
                        const char *saveName = (char*)"ks",
                        int iteration = -1) {
  // Handle axis shift if specified
  if (horizontal_shift != 0) {
    TH1* hist = dynamic_cast<TH1*>(hist1->Clone());
    const int nbins_axis = hist1->GetXaxis()->GetNbins();
    double newbins[nbins_axis+1];
    for (int ibin = 0; ibin <= nbins_axis; ++ibin)
      newbins[ibin] = hist1->GetXaxis()->GetBinLowEdge(ibin+1) + horizontal_shift;
    hist->SetBins(nbins_axis, newbins);
    hist1 = hist;
  }
  
  // Gather CDF data
  TAxis* axis1 = hist1->GetXaxis();
  TAxis* axis2 = hist2->GetXaxis();
  if (thresh_min == INT_MIN) {
    double temp = axis1->GetBinLowEdge(1);
    if (temp < axis2->GetBinLowEdge(1)) temp = axis2->GetBinLowEdge(1);
    thresh_min = temp;
  } if (thresh_max == INT_MAX) {
    double temp = axis1->GetBinLowEdge(axis1->GetNbins()+1);
    if (temp > axis2->GetBinLowEdge(axis2->GetNbins()+1))
      temp = axis2->GetBinLowEdge(axis2->GetNbins()+1);
    thresh_max = temp;
  }
  
  std::vector<std::vector<double>> cdf1 = CDF_FromHist_new(hist1, thresh_min, thresh_max);
  std::vector<std::vector<double>> cdf2 = CDF_FromHist_new(hist2, thresh_min, thresh_max);
  double n1 = cdf1.size();
  double n2 = cdf2.size();
  
  //Test edge cases
  if (cdf1.at(0).at(0) > cdf2.at(n2-1).at(0) ||
      cdf2.at(0).at(0) > cdf1.at(n1-1).at(0)) {
    std::cout<<"Warning in <utils/hist_tools::KS_statistic>: Input CDFs have disjoint domains. No meaningful comparison can be made."<<std::endl;
    return 1;
  }
  
  
  // Compute KS
  std::vector<std::vector<double>> ks;
  std::vector<double> opt_ks = {0, 0};
  std::vector<double> cpair_ks(2);
  int i1 = 0; int i2 = 0;
  std::vector<double> cpair1;
  std::vector<double> cpair2;
  std::vector<double> cpairT;
  while (i1 < n1 || i2 < n2) {
    if (i1 < n1) cpair1 = cdf1.at(i1);
    if (i2 < n2) cpair2 = cdf2.at(i2);
//    std::cout << "loop i1 : " << i1 << ", i2 : " << i2 << std::endl;
    
    if (cpair1.at(0) == cpair2.at(0)) {
//      std::cout << "are equal" << std::endl;
      // Same value, no need to interpolate
      cpair_ks.at(0) = cpair1.at(0);
      cpair_ks.at(1) = TMath::Abs(cpair1.at(1) - cpair2.at(1));
      ++i1; ++i2;
    } else if (cpair1.at(0) < cpair2.at(0)) {
//      std::cout << "1 < 2, reldif = " << cpair1.at(0) - cpair2.at(0) << std::endl;
      if (i2 == 0) {
        cpair_ks.at(0) = cpair1.at(0);
        cpair_ks.at(1) = cpair1.at(1);
        ++i1;
      } else if (i1 >= n1) {
        cpair_ks.at(0) = cpair2.at(0);
        cpair_ks.at(1) = 1 - cpair2.at(1);
//        std::cout << "end of cdf 1, looping over cdf2." << std::endl;
        ++i2;
      } else { // interpolate on CDF2
        cpair_ks.at(0) = cpair1.at(0);
//        if (i2 == 0) std::cout << "FLAG!" << std::endl;
        cpairT = cdf2.at(i2-1);
        cpair_ks.at(1) = TMath::Abs( ((cpair2.at(1) - cpairT.at(1)) / (cpair2.at(0) - cpairT.at(0)))
                                    *(cpair1.at(0) - cpairT.at(0)) + cpairT.at(1) - cpair1.at(1));
        ++i1;
//        std::cout << "interp on cdf2 with ks = " << cpair_ks.at(1) << std::endl;
      }
    } else if (cpair1.at(0) > cpair2.at(0)) {
//      std::cout << "1 > 2, reldif = " << cpair1.at(0) - cpair2.at(0) << std::endl;
      if (i1 == 0) {
        cpair_ks.at(0) = cpair2.at(0);
        cpair_ks.at(1) = cpair2.at(1);
        ++i2;
      } else if (i2 >= n2) {
        cpair_ks.at(0) = cpair1.at(0);
        cpair_ks.at(1) = 1 - cpair1.at(1);
//        std::cout << "end of cdf 2, looping over cdf1." << std::endl;
        ++i1;
      } else { // interpolate on CDF1
        cpair_ks.at(0) = cpair2.at(0);
//        if (i1 == 0) std::cout << "FLAG!" << std::endl;
        cpairT = cdf1.at(i1-1);
        cpair_ks.at(1) = TMath::Abs( ((cpair1.at(1) - cpairT.at(1)) / (cpair1.at(0) - cpairT.at(0)))
                                    *(cpair2.at(0) - cpairT.at(0)) + cpairT.at(1) - cpair2.at(1));
        ++i2;
//        std::cout << "interp on cdf1 with ks = " << cpair_ks.at(1) << std::endl;
      }
    }
//    std::cout << opt_ks.at(1) << std::endl;
    
    // ks push back, compare with opt
    if (cpair_ks.at(1) > opt_ks.at(1)) opt_ks = cpair_ks;
    ks.push_back(cpair_ks);
  } // End of KS loop
  
  // Plot if desired
  if (doPlot) {
    TCanvas* c = new TCanvas();
    c->SetWindowSize(500, 500);
    c->SetCanvasSize(1000,1000);
    c->Divide(2, 2);
    
    // Maybe worth trying to use gDirectory to catch the canvas/not keep deleting/remaking canvases.
    
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    TGraph* cdf1_graph = drawCDF_new(hist1, horizontal_shift, thresh_min, thresh_max, false, true);
    TGraph* cdf2_graph = drawCDF_new(hist2, 0, thresh_min, thresh_max, false, true);
    TMultiGraph* both_cdf = new TMultiGraph();
    both_cdf->SetTitle(Form(";%s;CDF", hist1->GetXaxis()->GetTitle()));
    setStyleLine(cdf1_graph, "thin red");
    cdf1_graph->SetMarkerColor(kRed+2);
    setStyleLine(cdf2_graph, "thin blue");
    cdf2_graph->SetMarkerColor(kBlue+2);
    both_cdf->Add(cdf1_graph);
    both_cdf->Add(cdf2_graph);
    double nonlog_rangemin = 0; //0.999 for non renorm
    both_cdf->GetYaxis()->SetRangeUser(nonlog_rangemin, 1);
    both_cdf->Draw("al*");
    TLine* thresholdLine = new TLine();
    setStyleLine(thresholdLine, "gray, dashed");
    thresholdLine->DrawLine(thresh_min, nonlog_rangemin, thresh_min, 1);
    thresholdLine->DrawLine(thresh_max, nonlog_rangemin, thresh_max, 1);
    TLine* ksLine = new TLine();
    setStyleLine(ksLine, "dotted orange thick");
    ksLine->DrawLine(opt_ks.at(0), nonlog_rangemin, opt_ks.at(0), 1);
    
//    TLegend *leg = new TLegend(0.5, 0.15, 0.94, 0.4);
//    leg->SetLineWidth(0);
//    leg->AddEntry(cdf2_graph, "Data p_{T} Spectrum CDF", "lp");
//    leg->AddEntry(cdf1_graph, "Shifted p_{T} #bf{pp} Reference CDF", "lp");
//    leg->AddEntry(ksLine, "KS line", "l");
//    leg->AddEntry(thresholdLine, "p_{T} Comparison Threshold", "l");
//    leg->Draw();
    
    // Make arrays from ks vector for graphing
    double axisval_array[ks.size()];
    double ksval_array[ks.size()];
    for (int iks = 0; iks < ks.size(); ++iks) {
      axisval_array[iks] = ks.at(iks).at(0);
      ksval_array[iks] = ks.at(iks).at(1);
    }
    
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    TGraph* ksgraph = new TGraph(ks.size(), axisval_array, ksval_array);
    ksgraph->SetTitle(Form(";%s;Abs(#Delta_{CDF})", hist1->GetXaxis()->GetTitle()));
    ksgraph->GetYaxis()->SetRangeUser(0, opt_ks.at(1)*1.2);
    ksgraph->Draw("al*");
    thresholdLine->DrawLine(thresh_min, 0, thresh_min, opt_ks.at(1)*1.2);
    thresholdLine->DrawLine(thresh_max, 0, thresh_max, opt_ks.at(1)*1.2);
    ksLine->DrawLine(opt_ks.at(0), 0, opt_ks.at(0), opt_ks.at(1)*1.2);
    drawText(Form("KS = max(#Delta_{CDF}): %.4f", opt_ks.at(1)), 0.9, 0.8, true);
    drawText(Form("KS found at %s = %.2f", hist1->GetXaxis()->GetTitle(), opt_ks.at(0)), 0.9, 0.75, true);
    
    if (iteration >= 0) drawText(Form("#Deltap_{T} KS Fitting Algorithm Iteration %i", iteration), 0.95, 0.93, true);
    
    c->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLogy();
    TGraph* cdf1_graph_log = drawCDF_new(hist1, horizontal_shift, thresh_min, thresh_max, true, true);
    TGraph* cdf2_graph_log = drawCDF_new(hist2, 0, thresh_min, thresh_max, true, true);
    TMultiGraph* both_cdf_log = new TMultiGraph();
    both_cdf_log->SetTitle(Form(";%s;1-CDF", hist1->GetXaxis()->GetTitle()));
    setStyleLine(cdf1_graph_log, "thin red");
    cdf1_graph_log->SetMarkerColor(kRed+2);
    setStyleLine(cdf2_graph_log, "thin blue");
    cdf2_graph_log->SetMarkerColor(kBlue+2);
    both_cdf_log->Add(cdf1_graph_log);
    both_cdf_log->Add(cdf2_graph_log);
    both_cdf_log->GetYaxis()->SetRangeUser(1e-11, 1);
    both_cdf_log->Draw("al*");
    thresholdLine->DrawLine(thresh_min, 0, thresh_min, 1);
    thresholdLine->DrawLine(thresh_max, 0, thresh_max, 1);
    ksLine->DrawLine(opt_ks.at(0), 0, opt_ks.at(0), 1);
    
    c->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLogy();
    TGraph* ksgraph_log = static_cast<TGraph*>(ksgraph->Clone());
    ksgraph_log->GetYaxis()->SetRangeUser(1e-11, 1);
    ksgraph_log->Draw("al*");
    thresholdLine->DrawLine(thresh_min, 0, thresh_min, 1);
    thresholdLine->DrawLine(thresh_max, 0, thresh_max, 1);
    ksLine->DrawLine(opt_ks.at(0), 0, opt_ks.at(0), 1);
    c->SaveAs(Form("%s.pdf",saveName));
//    c->SaveAs(Form("../tmp/tmpplot/%s.pdf",saveName));
    
    delete c;
  }
  
  // return ks
  return opt_ks.at(1);
}

//========================================================================== Placeholder Methods

// Simple translate method for TH2D, assumes uniform binning.
// Only translates histogram within the resolution of bins,
// i.e. moving by the nearest bin division to input
TH2D* translateHist_simple(TH2D* hist, 
                           double delta_x = 0,
                           double delta_y = 0,
                           bool suppress_warnings = false) {
  if (delta_x == 0 && delta_y == 0) {
    std::vector<double> cm = getCM(hist, 21, true);
    delta_x = -cm.at(0);
    delta_y = -cm.at(1);
    if (!suppress_warnings) cout << "Defaulting to CM at x = " << cm.at(0) << ", y = " << cm.at(1) << endl;
  }
  
  double width_x = hist->GetXaxis()->GetBinWidth(1);
  double width_y = hist->GetYaxis()->GetBinWidth(1);
  int nbins_x = hist->GetXaxis()->GetNbins();
  int nbins_y = hist->GetYaxis()->GetNbins();
  double histbounds[2][2] = {
    {hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinLowEdge(nbins_x+1)},
    {hist->GetYaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(nbins_y+1)} };
  if (TMath::Abs(delta_x) <= width_x/2 && TMath::Abs(delta_y) <= width_y/2) {
    if (!suppress_warnings) cout << "Shift too small to affect hist. Returned original hist." << endl;
    return hist;
  }
  
  TH2D* out_hist = static_cast<TH2D*>(hist->Clone());
  out_hist->Reset();
  
  double cx, cy, cfill;
  double lostdata = 0;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    for (int iy = 1; iy <= nbins_y; ++iy) {
      cfill = hist->GetBinContent(ix, iy);
      cx = hist->GetXaxis()->GetBinCenter(ix) + delta_x;
      cy = hist->GetYaxis()->GetBinCenter(iy) + delta_y;
      
      if (cx <= histbounds[0][0] || cx >= histbounds[0][1] ||
          cy <= histbounds[1][0] || cy >= histbounds[1][1]) {
        lostdata += cfill*width_x * width_y;
        continue;
      }out_hist->Fill(cx, cy, cfill);
    }
  }
  
  // Report on lost data, store in public variable
  if (!suppress_warnings && lostdata > 0)
    cout << Form("Warning in translateHist_simple: %.5f%% of data lost to translation over edge.",
                 100.*lostdata / hist->Integral("width")) << endl;
  lostdata_translate = lostdata;
  
  return out_hist;
}


TH2D* rotateHist2D_simple(TH2D* hist, double theta, bool suppress_warnings = false) {
  double width_x = hist->GetXaxis()->GetBinWidth(1);
  double width_y = hist->GetYaxis()->GetBinWidth(1);
  int nbins_x = hist->GetXaxis()->GetNbins();
  int nbins_y = hist->GetYaxis()->GetNbins();
  double histbounds[2][2] = {
    {hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinLowEdge(nbins_x+1)},
    {hist->GetYaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(nbins_y+1)} };
  if (TMath::Abs(theta) <= TMath::Abs(TMath::ATan(width_y / (width_x * nbins_x)))) {
    if (!suppress_warnings) cout << "Warning in rotatehist2D_simple: Angle too small to affect hist; returned original hist." << endl;
    return hist;
  }
  
  TH2D* out_hist = static_cast<TH2D*>(hist->Clone());
  out_hist->Reset();
  
  double cx, cy, cr, ctheta, cfill;
  double lostdata = 0;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    for (int iy = 1; iy <= nbins_y; ++iy) {
      // Tally lost data by manually translating each item
      cfill = hist->GetBinContent(ix, iy);
      cx = hist->GetXaxis()->GetBinCenter(ix);
      cy = hist->GetYaxis()->GetBinCenter(iy);
      cr = TMath::Sqrt(cx*cx + cy*cy);
      ctheta = TMath::ATan2(cy, cx);
      cx = cr*TMath::Cos(ctheta + theta);
      cy = cr*TMath::Sin(ctheta + theta);
      
      if (cx <= histbounds[0][0] || cx >= histbounds[0][1] ||
          cy <= histbounds[1][0] || cy >= histbounds[1][1]) {
        lostdata += cfill * width_x * width_y;
      }
      
      // Fill hist by forward rotating to get bin content
      // This is necessary to avoid 2 bins filling the same rotated bin
      cx = cr*TMath::Cos(ctheta - theta);
      cy = cr*TMath::Sin(ctheta - theta);
      if (cx <= histbounds[0][0] || cx >= histbounds[0][1] ||
          cy <= histbounds[1][0] || cy >= histbounds[1][1]) {
        continue;
      }
      cfill = hist->GetBinContent(hist->GetXaxis()->FindBin(cx),
                                  hist->GetYaxis()->FindBin(cy));
      out_hist->SetBinContent(ix, iy, cfill);
    }
  }
  
  // Report on lost data, store in public variable
  if (!suppress_warnings && lostdata > 0)
    cout << Form("Warning in rotateHist2D_simple: %.5f%% of data lost to rotation over edge.",
                 100.*lostdata / hist->Integral("width")) << endl;
  lostdata_rotate = lostdata;
  
  return out_hist;
}
