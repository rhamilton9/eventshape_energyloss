// A collection of small macros for histogram analysis
// Macros are organized roughly by function

// Variables for storing data that may be of interest in the daughter macro

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

// Returns gradient-norm map of input 2D histogram.
// Gradient is computed using finite difference approximation.
//
// CONSIDERATIONS :: Could be implemented more elegantly, or include
//    higher order corrections to finite difference upon request.
TH2D* gradientNorm(TH2D* hist) {
  TH2D* out_hist = (TH2D*) hist->Clone();
  out_hist->Reset();
  
  const int nbins_x = hist->GetXaxis()->GetNbins();
  const int nbins_y = hist->GetYaxis()->GetNbins();
  
  double ux, uy;
  double cx, cy;
  double lx, ly;
  double hx, hy;
  double dx, dy;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    if (ix == 1) {
      lx = ix;
      ux = ix + 1;
      cx = ix;
      hx = hist->GetXaxis()->GetBinWidth(ix);
    } else if (ix == nbins_x) {
      lx = ix-1;
      cx = ix;
      ux = ix;
      hx = hist->GetXaxis()->GetBinWidth(ix);
    } else {
      lx = ix - 1;
      cx = ix;
      ux = ix + 1;
      hx = hist->GetXaxis()->GetBinWidth(ix-1)
            +hist->GetXaxis()->GetBinWidth(ix);
    }
    
    for (int iy = 1; iy <= nbins_y; ++iy) {
      if (iy == 1) {
        ly = iy;
        uy = iy + 1;
        cy = iy;
        hy = hist->GetYaxis()->GetBinWidth(iy);
      } else if (iy == nbins_y) {
        ly = iy-1;
        cy = iy;
        uy = iy;
        hy = hist->GetYaxis()->GetBinWidth(iy);
      } else {
        ly = iy - 1;
        cy = iy;
        uy = iy + 1;
        hy = hist->GetYaxis()->GetBinWidth(iy-1)
              +hist->GetYaxis()->GetBinWidth(iy);
      }
      
      // compute norm gradient with finite difference
      out_hist->Fill(hist->GetXaxis()->GetBinCenter(ix),
                     hist->GetXaxis()->GetBinCenter(iy),
                     TMath::Sqrt(TMath::Power((hist->GetBinContent(ux, cy) - hist->GetBinContent(lx, cy)) / hx, 2) +
                                 TMath::Power((hist->GetBinContent(cx, uy) - hist->GetBinContent(cx, ly)) / hy, 2)));
    }
  }return out_hist;
}

//========================================================================== Statistical Tools

// Computes the Empirical Cumulative Distribution Function for a 1D hist.
// -----------------------*IMPORTANT*-----------------------
// Assumes the probability distribution within a single bin is uniform.
// ---------------------------------------------------------
// Input must be either TH1D or TH1F type; other types will throw an error.
// The eCDF is stored as a vector that can be used for other analysis.
std::vector<double> CDF_FromHist(TH1* hist) {
  // Verify that input hist is of compatible type
  int histflag = getHistType(hist);
  std::vector<double> cdf;
  switch (histflag) {
    case 11: case 12: break;
    default:
      std::cout<<"Error in CDF_FromHist :: Input hist type is not TH1F or TH1D."<<std::endl;
      return cdf;
  }
  
  // Take integral for later rescaling
  double norm = hist->Integral("width");
  
  // Extract CDF assuming locally uniform prior
  int nbin = hist->GetXaxis()->GetNbins();
  cdf.push_back(0);
  for (int ibin = 1; ibin <= nbin; ++ibin) {
    cdf.push_back(cdf.at(ibin - 1) + (hist->GetBinContent(ibin) * hist->GetXaxis()->GetBinWidth(ibin))/norm);
  }cdf.push_back(1);
  return cdf;
}

// Draws and returns a graphical representation of the CDF
// Result is drawn on the current gPad unless suppressed.
TGraph* drawCDF(TH1* hist,
                double horiz_shift = 0,
                bool drawLogy = false,
                bool suppressDraw = false) {
  std::vector<double> cdf = CDF_FromHist(hist);
  int nbin = hist->GetXaxis()->GetNbins();
  
  // Establish bin threshold + test for edge cases
  double bin_axis[nbin+1];
  double cdf_axis[nbin+1];
  
  // Loop over CDF and assign axis vals to CDF
  for (int ibin = 0; ibin <= nbin; ++ibin) {
    bin_axis[ibin] = hist->GetXaxis()->GetBinLowEdge(ibin+1) + horiz_shift;
    if (drawLogy) cdf_axis[ibin] = 1-cdf.at(ibin);
    else          cdf_axis[ibin] = cdf.at(ibin);
  }
  
  // Construct TGraph and plot if specified
  TGraph* cdf_graph = new TGraph(nbin+1, bin_axis, cdf_axis);
  cdf_graph->SetTitle(Form(";%s;eCDF",hist->GetXaxis()->GetTitle()));
  if (!suppressDraw) cdf_graph->Draw("al");
  return cdf_graph;
}

// Computes the Kalgomorov-Smirnov (KS) statistic for 2 histograms
//
// The variable horizShiftOnHist1 allows for the CDF of hist1 to be
// shifted horizontally by some amount before computing the KS.
// Only horizontal axis values above comparison_threshold are considered.
double KS_statistic(TH1* hist1,
                    TH1* hist2,
                    double horizShiftOnHist1 = 0,
                    double comparison_threshold = INT_MIN,
                    bool suppress_CDF_below_threshold = false,
                    bool doPlot = false,
                    char *saveName = (char*)"ks",
                    int iteration = -1) {
  TAxis* axis1 = hist1->GetXaxis();
  TAxis* axis2 = hist2->GetXaxis();
  int nbin_1 = axis1->GetNbins();
  int nbin_2 = axis2->GetNbins();
  if (horizShiftOnHist1 != 0) {
    //Make new axis.
    double bins[nbin_1 + 1];
    for (int i = 0; i <= nbin_1; ++i)
      bins[i] = axis1->GetBinLowEdge(i+1) + horizShiftOnHist1;
    axis1 = new TAxis(nbin_1, bins);
  }
  
  // If requested, renormalize the computed CDF and recalculate KS above the threshold
  // To do this, create a new histogram with only bin information above the threshold.
  if (suppress_CDF_below_threshold) {
    int startbin_1 = axis1->FindBin(comparison_threshold);
    int startbin_2 = axis2->FindBin(comparison_threshold);
    std::vector<double> binedge1_vector;
    std::vector<double> binedge2_vector;
    binedge1_vector.push_back(comparison_threshold);
    binedge2_vector.push_back(comparison_threshold);
    
    for (int ibin = startbin_1; ibin <= nbin_1; ++ibin) {
      if (axis1->GetBinLowEdge(ibin+1) > axis2->GetBinLowEdge(nbin_2+1)) {
//        binedge1_vector.push_back(axis2->GetBinLowEdge(nbin_2+1));
        break;
      }binedge1_vector.push_back(axis1->GetBinLowEdge(ibin+1));
    }for (int ibin = startbin_2; ibin <= nbin_2; ++ibin) {
      if (axis2->GetBinLowEdge(ibin+1) > axis1->GetBinLowEdge(nbin_1+1)) {
//        binedge2_vector.push_back(axis1->GetBinLowEdge(nbin_1+1));
        break;
      }binedge2_vector.push_back(axis2->GetBinLowEdge(ibin+1));
    }
    
    double binedge1_array[binedge1_vector.size()];
    double binedge2_array[binedge2_vector.size()];
    for (int i = 0; i < binedge1_vector.size(); ++i)
      binedge1_array[i] = binedge1_vector.at(i);
    for (int i = 0; i < binedge2_vector.size(); ++i)
      binedge2_array[i] = binedge2_vector.at(i);
    
//    for (int i = 0; i < binedge2_vector.size(); ++i) cout << binedge2_array[i] << endl;
    
    TH1D* temphist1 = new TH1D("", hist1->GetTitle(), binedge1_vector.size()-1, binedge1_array);
    TH1D* temphist2 = new TH1D("", hist2->GetTitle(), binedge2_vector.size()-1, binedge2_array);
    
    for (int i = 1; i <= binedge1_vector.size()-1; ++i)
      temphist1->SetBinContent(i, hist1->GetBinContent(i+startbin_1-1));
    for (int i = 1; i <= binedge1_vector.size()-2; ++i)
      temphist2->SetBinContent(i, hist2->GetBinContent(i+startbin_2-1));
    
    hist1 = temphist1;
    hist2 = temphist2;
    axis1 = hist1->GetXaxis();
    axis2 = hist2->GetXaxis();
    nbin_1 = axis1->GetNbins();
    nbin_2 = axis2->GetNbins();
  }
  
  
  
  std::vector<double> cdf1 = CDF_FromHist(hist1);
  std::vector<double> cdf2 = CDF_FromHist(hist2);
  cdf1.push_back(1); cdf2.push_back(1);
  
  //debug
//  std::cout << "CDF1 ::" << std::endl;
//  for (double d:cdf1) std::cout<<d<<endl;
//  std::cout << "CDF2 ::" << std::endl;
//  for (double d:cdf2) std::cout<<d<<endl;
//  std::cout << "end CDF" << std::endl;
  
  int start_bin1 = axis1->FindBin(comparison_threshold);
  int start_bin2 = axis2->FindBin(comparison_threshold);
  int i1 = 0;
  int i2 = 0;
  double ks = 0;
  double axis_ks;
  double ks_local, interp;
  std::vector<double> ksvals;
  std::vector<double> axisvals;
  bool checkCompFlag = false;
  double cdf_at_threshold[2] = {0, 0};
  if ((comparison_threshold < axis1->GetBinLowEdge(1) &&
       comparison_threshold < axis2->GetBinLowEdge(1) ) ||
      (comparison_threshold > axis1->GetBinLowEdge(nbin_1 + 1) &&
       comparison_threshold > axis2->GetBinLowEdge(nbin_2 + 1) )) checkCompFlag = true;
  do {
//    cout << i1 << ": " << axis1->GetBinLowEdge(i1+1) << " " << i2 << ": " << axis2->GetBinLowEdge(i2+1) << endl;
    // Check if the comparison threshold is between the current values
    if (!checkCompFlag &&
        axis1->GetBinLowEdge(i1+1) > comparison_threshold &&
        axis2->GetBinLowEdge(i2+1) > comparison_threshold &&
        axis1->GetBinLowEdge(i1) < comparison_threshold &&
        axis2->GetBinLowEdge(i2) < comparison_threshold) {
      cdf_at_threshold[0] = ( ((cdf1.at(i1) - cdf1.at(i1-1)) /
                               (axis1->GetBinLowEdge(i1+1) - axis1->GetBinLowEdge(i1)))
                             *(comparison_threshold - axis1->GetBinLowEdge(i1)) + cdf1.at(i1-1) );
      cdf_at_threshold[1] = ( ((cdf2.at(i2) - cdf2.at(i2-1)) /
                               (axis2->GetBinLowEdge(i2+1) - axis2->GetBinLowEdge(i2)))
                             *(comparison_threshold - axis2->GetBinLowEdge(i2)) + cdf2.at(i2-1) );
      ks_local = TMath::Abs(cdf_at_threshold[1] - cdf_at_threshold[0]);
      axisvals.push_back(comparison_threshold);
      checkCompFlag = true;
//      cout << "added comp with ks = " << ks_local << std::endl;
    } else if (TMath::Abs(axis1->GetBinLowEdge(i1+1) - axis2->GetBinLowEdge(i2+1)) < 1e-10) {
      // Begin regular cases
      
      // Both have same edge; no need to interpolate
      ks_local = TMath::Abs(cdf1.at(i1) - cdf2.at(i2));
      axisvals.push_back(axis1->GetBinLowEdge(i1+1));
//      std::cout << "same!" << std::endl;
//      if (i1 < nbin_1) ++i1;
//      if (i2 < nbin_2) ++i2;
      ++i1;
      ++i2;
    } else if (axis1->GetBinLowEdge(i1+1) < axis2->GetBinLowEdge(i2+1)) {
      // Compute at current hist1 bin edge
      if (i2 == 0) { // Test edge cases (no overlap)
        ks_local = cdf1.at(i1);
        axisvals.push_back(axis1->GetBinLowEdge(i1+1));
        ++i1;
      } else if (i1 > nbin_1) { // not actually using hist1, hist1 is out of entries.
        ks_local = 1 - cdf2.at(i2);
        axisvals.push_back(axis2->GetBinLowEdge(i2+1));
        ++i2;
      } else {
        // Interpolate on CDF2 to match val at CDF1 (this is where uniform assumption appears)
        interp = ( ((cdf2.at(i2) - cdf2.at(i2-1)) /
                     (axis2->GetBinLowEdge(i2+1) - axis2->GetBinLowEdge(i2)))
                  *(axis1->GetBinLowEdge(i1+1) - axis2->GetBinLowEdge(i2)) + cdf2.at(i2-1) );
        ks_local = TMath::Abs(interp - cdf1.at(i1));
        axisvals.push_back(axis1->GetBinLowEdge(i1+1));
        ++i1;
      }
    } else { 
      // Compute at current hist2 bin edge
      if (i1 == 0) { // Test edge cases (no overlap)
        ks_local = cdf2.at(i2);
        axisvals.push_back(axis2->GetBinLowEdge(i2+1));
        ++i2;
      } else if (i2 > nbin_2) { // not actually using hist2, hist2 is out of entries.
        ks_local = 1 - cdf1.at(i1);
        axisvals.push_back(axis1->GetBinLowEdge(i1+1));
        ++i1;
      } else {
        // Interpolate on CDF1 to match val at CDF2 (this is where uniform assumption appears)
        interp = ( ((cdf1.at(i1) - cdf1.at(i1-1)) /
                    (axis1->GetBinLowEdge(i1+1) - axis1->GetBinLowEdge(i1)))
                  *(axis2->GetBinLowEdge(i2+1) - axis1->GetBinLowEdge(i1)) + cdf1.at(i1-1) );
        ks_local = TMath::Abs(interp - cdf2.at(i2));
        axisvals.push_back(axis2->GetBinLowEdge(i2+1));
        ++i2;
      }
    }ksvals.push_back(ks_local);
    
    //debug
//    std::cout<< Form("%.2f,\t\tks: ",axisvals.back()) << ks_local << std::endl;
//    std::cout << Form("axis1: %.2f, axis2: %.2f", axis1->GetBinLowEdge(i1+1), axis2->GetBinLowEdge(i2+1)) << std::endl;
    
    if (axisvals.back() < comparison_threshold) continue;
    
    if (ks_local > ks) {
      ks = ks_local;
      axis_ks = axisvals.back();
    }
  } while (i1 + i2 <= nbin_1 + nbin_2);
  
//  std::cout << "max ks " << ks << std::endl;
  
  ksvals.push_back(0);
  if (axis2->GetBinLowEdge(nbin_2+1) > axis1->GetBinLowEdge(nbin_1+1))
    axisvals.push_back(axis2->GetBinLowEdge(nbin_2+1));
  else
    axisvals.push_back(axis1->GetBinLowEdge(nbin_1+1));
  
  //debug
//  cout << ksvals.size() << endl;
//  cout << axisvals.size() << endl;
  
  // Make plot to check macro is working
  // mainly for debugging purposes, macro would work without it.
  double axisval_array[axisvals.size()];
  for (int i = 0; i < axisvals.size(); ++i) {
    axisval_array[i] = axisvals.at(i);
  }
  double ksval_array[ksvals.size()];
  for (int i = 0; i < ksvals.size(); ++i) {
    ksval_array[i] = ksvals.at(i);
  }
  
  if (doPlot) {
    TCanvas* c = new TCanvas();
    c->SetWindowSize(500, 500);
    c->SetCanvasSize(1000,1000);
    c->Divide(2, 2);
    
    // Maybe worth trying to use gDirectory to catch the canvas/not keep deleting/remaking canvases.
    
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    TGraph* cdf1_graph = drawCDF(hist1, horizShiftOnHist1, false, true);
    TGraph* cdf2_graph = drawCDF(hist2, 0, false, true);
    TMultiGraph* both_cdf = new TMultiGraph();
    both_cdf->SetTitle(Form(";%s;CDF", hist1->GetXaxis()->GetTitle()));
    cdf1_graph->SetLineColor(kRed+2);
    cdf1_graph->SetMarkerColor(kRed+2);
    cdf1_graph->SetMarkerStyle(24);
    cdf2_graph->SetLineColor(kBlue+2);
    cdf2_graph->SetMarkerColor(kBlue+2);
    cdf2_graph->SetMarkerStyle(25);
    both_cdf->Add(cdf1_graph);
    both_cdf->Add(cdf2_graph);
    double nonlog_rangemin = 0; //0.999 for non renorm
    both_cdf->GetYaxis()->SetRangeUser(nonlog_rangemin, 1);
    both_cdf->Draw("al*");
    TLine* thresholdLine = new TLine();
//    thresholdLine->SetLineColor(kGray+1);
//    thresholdLine->SetLineWidth(2);
    TLine* ksLine = new TLine();
    ksLine->SetLineColor(kOrange+10);
    ksLine->SetLineWidth(2);
    ksLine->SetLineStyle(3);
    thresholdLine->DrawLine(comparison_threshold, nonlog_rangemin, comparison_threshold, 1);
    ksLine->DrawLine(axis_ks, nonlog_rangemin, axis_ks, 1);
    
    TLegend *leg = new TLegend(0.5, 0.15, 0.94, 0.4);
    leg->SetLineWidth(0);
    leg->AddEntry(cdf2_graph, "Data p_{T} Spectrum CDF", "lp");
    leg->AddEntry(cdf1_graph, "Shifted p_{T} #bf{pp} Reference CDF", "lp");
    leg->AddEntry(ksLine, "KS line", "l");
    leg->AddEntry(thresholdLine, "p_{T} Comparison Threshold", "l");
    leg->Draw();
    
    
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    TGraph* ksgraph = new TGraph(axisvals.size(), axisval_array, ksval_array);
    ksgraph->SetTitle(Form(";%s;Abs(#Delta_{CDF})", hist1->GetXaxis()->GetTitle()));
    ksgraph->GetYaxis()->SetRangeUser(0, ks*1.2);
    ksgraph->Draw("al*");
    thresholdLine->DrawLine(comparison_threshold, 0, comparison_threshold, ks*1.2);
    ksLine->DrawLine(axis_ks, 0, axis_ks, ks*1.2);
    drawText(Form("KS = max(#Delta_{CDF}): %.4f", ks), 0.9, 0.8, true);
    drawText(Form("KS found at %s = %.2f", hist1->GetXaxis()->GetTitle(), axis_ks), 0.9, 0.75, true);
    
    if (iteration >= 0) drawText(Form("#Deltap_{T} KS Fitting Algorithm Iteration %i", iteration), 0.95, 0.93, true);
    
    c->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLogy();
    TGraph* cdf1_graph_log = drawCDF(hist1, horizShiftOnHist1, true, true);
    TGraph* cdf2_graph_log = drawCDF(hist2, 0, true, true);
    TMultiGraph* both_cdf_log = new TMultiGraph();
    both_cdf_log->SetTitle(Form(";%s;1-CDF", hist1->GetXaxis()->GetTitle()));
    cdf1_graph_log->SetLineColor(kRed+2);
    cdf1_graph_log->SetMarkerColor(kRed+2);
    cdf1_graph_log->SetMarkerStyle(24);
    cdf2_graph_log->SetLineColor(kBlue+2);
    cdf2_graph_log->SetMarkerColor(kBlue+2);
    cdf2_graph_log->SetMarkerStyle(25);
    both_cdf_log->Add(cdf1_graph_log);
    both_cdf_log->Add(cdf2_graph_log);
    both_cdf_log->GetYaxis()->SetRangeUser(1e-11, 1);
    both_cdf_log->Draw("al*");
    thresholdLine->DrawLine(comparison_threshold, 0, comparison_threshold, 1);
    ksLine->DrawLine(axis_ks, 0, axis_ks, 1);
    
    c->cd(4);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLogy();
    TGraph* ksgraph_log = static_cast<TGraph*>(ksgraph->Clone());
    ksgraph_log->GetYaxis()->SetRangeUser(1e-11, 1);
    ksgraph_log->Draw("al*");
    thresholdLine->DrawLine(comparison_threshold, 0, comparison_threshold, 1);
    ksLine->DrawLine(axis_ks, 0, axis_ks, 1);
    c->SaveAs(Form("../tmp/tmpplot/%s.pdf",saveName));
    
    delete c;
  }
  return ks;
}

//========================================================================== KS Rewrite

std::vector<std::vector<double>> CDF_FromHist_new(TH1* hist, double thresh_min = INT_MIN, double thresh_max = INT_MAX) {
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
  } if (thresh_min > thresh_max) {
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
  }return cdf;
}

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

double KS_statistic_new(TH1* hist1,
                        TH1* hist2,
                        double horizShiftOnHist1 = 0,
                        double comparison_threshold = INT_MIN,
                        bool suppress_CDF_below_threshold = false,
                        bool doPlot = false,
                        char *saveName = (char*)"ks",
                        int iteration = -1) {
  std::vector<std::vector<double>> cdf1 = CDF_FromHist_new(hist1);
  std::vector<std::vector<double>> cdf2 = CDF_FromHist_new(hist2);
  double n1 = cdf1.size();
  double n2 = cdf2.size();
  
  // Handle comparison_threshold and CDF suppression
  
  // Compute KS
  std::vector<std::vector<double>> ks;
  std::vector<double> cpair_ks(2);
  int i1 = 0;
  int i2 = 0;
  double caxis_1, caxis_2, cval1, cval2;
  while (i1 < n1 && i2 < n2) {
    caxis_1 = cdf1.at(i1).at(0);
    caxis_2 = cdf2.at(i2).at(0);
  }
  
  // Plot if desired
  
  // return ks
  return 0;
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
