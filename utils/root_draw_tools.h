// A collection of small macros to make certain tasks simpler with root.
// Macros are organized roughly by function

#include "TPad.h"
#include "TH1.h"
#include "TGraph.h"
#include "TLatex.h"

//========================================================================== Construction of TObjects

//Draws text as specified by the input parameters
//Returns pointer for further modification if needed
//  e.g. drawText(...)->SetTextAngle(90);
TLatex* drawText(const char *text,
              float xp,
              float yp,
              bool isRightAlign = 0,
              int textColor=kBlack,
              double textSize=0.04,
              int textFont = 42,
              bool isNDC = true){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  if(isNDC) tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
  return tex;
}

//Constructs, draws and returns a TPad with the given specifications
TPad* buildPad(const char *padName,
               float x1, float y1,
               float x2, float y2,
               float leftMargin = 0.1,
               float rightMargin = 0.1,
               float bottomMargin = 0.1,
               float topMargin = 0.1,
               bool isTransparent = true) {
  TPad *pad = new TPad(padName, "", x1, y1, x2, y2);
  pad->SetLeftMargin(leftMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetTopMargin(topMargin);
  pad->SetBottomMargin(bottomMargin);
  if (isTransparent) {pad->SetFillStyle(4000); pad->SetFrameFillStyle(4000);}
  pad->Draw();
  return pad;
}

//========================================================================== Style Settings

void setStyleLine(TAttLine* line, const char* style) {
  TString str(style);
  str.ToLower();
  // Line width
  if      (str.Contains("thick"))     line->SetLineWidth(3);
  else if (str.Contains("thin"))      line->SetLineWidth(1);
  else                                line->SetLineWidth(2);
  
  // Line style
  if      (str.Contains("dashed"))    line->SetLineStyle(7);
  else if (str.Contains("longdash"))  line->SetLineStyle(9);
  else if (str.Contains("shortdash")) line->SetLineStyle(2);
  else if (str.Contains("dotted"))    line->SetLineStyle(3);
  else if (str.Contains("-."))        line->SetLineStyle(4);
  else if (str.Contains("--."))       line->SetLineStyle(5);
  else if (str.Contains("---."))      line->SetLineStyle(9);
  else if (str.Contains("...-"))      line->SetLineStyle(6);
  else if (str.Contains("..-"))       line->SetLineStyle(8);
  else                                line->SetLineStyle(1);
  
  // Line color
  if (str.Contains("gray")
      || str.Contains("grey")) line->SetLineColor(kGray+1);
  else if (str.Contains("red")) line->SetLineColor(kRed+2);
  else if (str.Contains("green")) line->SetLineColor(kGreen+2);
  else if (str.Contains("blue")) line->SetLineColor(kBlue+2);
  else if (str.Contains("yellow")) line->SetLineColor(kYellow+2);
  else if (str.Contains("orange")) line->SetLineColor(kOrange+10);
  else if (str.Contains("green")) line->SetLineColor(kGreen+3);
  else if (str.Contains("sky")) line->SetLineColor(kAzure+5);
  else if (str.Contains("violet")) line->SetLineColor(kViolet+7);
}

//========================================================================== Misc. Drawing Tools

// Draw a dashed gray line through unity on the input pad.
void drawUnityLine(TAxis* reference) {
  TLine* unityline = new TLine();
  setStyleLine(unityline, "gray, dashed, thin");
  unityline->DrawLine(reference->GetXmin(), 1, reference->GetXmax(), 1);
  return;
}

// Flip a histogram to its side or upside-down by drawing a TGraph
// Inputs give
//      0, 0 - sideways with origin to the left
//      1, 0 - sideways with origin to the right
//      0, 1 - upside-down (reflected about horizontal)
//      1, 1 - rotated 180 degrees (reflected about both axes)
// Input nodraw can be given to supress drawing, if e.g. graph is needed
// as input for TMultiGraph (many rotated graphs are drawn)
// If drawing is supressed, settings for flipping axes (rx, ry) must be
// manually written into the TGraph options when plotting elsewhere.
TGraph* rotateHist(TH1* hist, bool horiz, bool vert, bool dodraw = true) {
  const int nbin = hist->GetXaxis()->GetNbins();
  double x_edge[2*(nbin+1)+1];
  double y_edge[2*(nbin+1)+1];
  
  // Sweep over data to gather hist edges
  y_edge[0] = hist->GetYaxis()->GetBinLowEdge(1);
  for (int i = 0; i < nbin; ++i) {
    x_edge[2*i] = hist->GetBinLowEdge(i+1);
    x_edge[2*i+1] = x_edge[2*i];
    
    y_edge[2*i+1] = hist->GetBinContent(i+1);
    y_edge[2*(i+1)] = y_edge[2*i+1];
  }
  x_edge[2*nbin] = x_edge[2*nbin-1] + hist->GetBinWidth(nbin);
  x_edge[2*nbin+1] = x_edge[2*nbin];
  x_edge[2*nbin+2] = x_edge[0];
  
  y_edge[2*nbin+1] = y_edge[0];
  y_edge[2*nbin+2] = x_edge[0];
  
  // Construct graph with swaps as necessary to reflect input
  TGraph *reflGraph;
  if (!vert) {
    reflGraph = new TGraph(2*(nbin+1)+1, y_edge, x_edge);
    reflGraph->SetTitle(Form("%s;%s;%s",
                             hist->GetTitle(),
                             hist->GetYaxis()->GetTitle(),
                             hist->GetXaxis()->GetTitle()));
  } else {
    reflGraph = new TGraph(2*(nbin+1)+1, x_edge, y_edge);
    reflGraph->SetTitle(Form("%s;%s;%s",
                             hist->GetTitle(),
                             hist->GetXaxis()->GetTitle(),
                             hist->GetYaxis()->GetTitle()));
  }
  
  // Get Color and other settings from original hist
  reflGraph->SetLineWidth(hist->GetLineWidth());
  reflGraph->SetLineColor(hist->GetLineColor());
  reflGraph->SetLineStyle(hist->GetLineStyle());
  reflGraph->SetFillColor(hist->GetFillColor());
  reflGraph->SetFillStyle(hist->GetFillStyle());
  
  // Draw if not supressed
  if (dodraw) {
    switch(2*horiz+3*vert) {
      case 0: // (0, 0)
        reflGraph->Draw("afl"); break;
      case 2: // (1, 0)
        reflGraph->Draw("afl rx"); break;
      case 3: // (0, 1)
        reflGraph->Draw("afl ry"); break;
      case 5: // (1, 1)
        reflGraph->Draw("afl rx ry"); break;
    }
  }return reflGraph;
}
