// A collection of small macros to make certain tasks simpler with root.
// Macros are organized roughly by function

#include "TPad.h"
#include "TH1.h"
#include "TGraph.h"
#include "TLatex.h"

//========================================================================== Root-Default Curve Interpolation Macro

// Reproduced from TGraphPainter.cxx source file:
// https://root.cern/doc/master/TGraphPainter_8cxx_source.html
// Line 4688:
// https://root.cern/doc/master/classTGraphPainter.html#abfe5af15f1920e7b11dd3e15480784ed

// Below is the ROOT local documentation:
////////////////////////////////////////////////////////////////////////////////
/// Smooth a curve given by N points.
///
/// The original code is from an underlaying routine for Draw based on the
/// CERN GD3 routine TVIPTE:
///
/// Author - Marlow etc.   Modified by - P. Ward     Date -  3.10.1973
///
/// This method draws a smooth tangentially continuous curve through
/// the sequence of data points P(I) I=1,N where P(I)=(X(I),Y(I)).
/// The curve is approximated by a polygonal arc of short vectors.
/// The data points can represent open curves, P(1) != P(N) or closed
/// curves P(2) == P(N). If a tangential discontinuity at P(I) is
/// required, then set P(I)=P(I+1). Loops are also allowed.
///
/// Reference Marlow and Powell, Harwell report No.R.7092.1972
/// MCCONALOGUE, Computer Journal VOL.13, NO4, NOV1970P p392 6
///
///  Original publication: "A quasi-intrinsic scheme for passing a smooth curve through a discrete set of points,"
///  https://academic.oup.com/comjnl/article/13/4/392/540402
///
/// -  npoints   : Number of data points.
/// -  x         : Abscissa
/// -  y         : Ordinate
TGraph* Smooth(TGraph *theGraph, bool suppress_info = false) {
  const Int_t npoints = theGraph->GetN();
  Double_t *x = theGraph->GetX();
  Double_t *y = theGraph->GetY();
  Int_t drawtype = 1000;
  
  Int_t i, k, kp, km, npointsMax, banksize, n2, npt;
  Int_t maxiterations, finished;
  Int_t jtype, ktype, closed;
  Double_t sxmin, sxmax, symin, symax;
  Double_t delta;
  Double_t xorg, yorg;
  Double_t ratio_signs, xratio, yratio;
  Int_t flgic, flgis;
  Int_t iw, loptx;
  Double_t p1, p2, p3, p4, p5, p6;
  Double_t w1, w2, w3;
  Double_t a, b, c, r, s=0.0, t, z;
  Double_t co, so, ct, st, ctu, stu, xnt;
  Double_t dx1, dy1, dx2, dy2, dk1, dk2;
  Double_t xo, yo, dx, dy, xt, yt;
  Double_t xa, xb, ya, yb;
  Double_t u1, u2, u3, tj;
  Double_t cc, err;
  Double_t sb, sth;
  Double_t wsign, tsquare, tcube;
  c = t = co = so = ct = st = ctu = stu = dx1 = dy1 = dx2 = dy2 = 0;
  xt = yt = xa = xb = ya = yb = u1 = u2 = u3 = tj = sb = 0;
  
  npointsMax  = npoints*10;
  n2          = npointsMax-2;
  banksize    = n2;
  
  std::vector<Double_t> qlx(npointsMax);
  std::vector<Double_t> qly(npointsMax);
  if (qlx.empty() || qly.empty()) {
    Error("Smooth", "not enough space in memory");
    return new TGraph();
  }
  
  //  Decode the type of curve (draw type).
  
  loptx = kFALSE;
  jtype  = (drawtype%1000)-10;
  if (jtype > 0) { ktype = jtype; loptx = kTRUE; }
  else             ktype = drawtype%1000;
  
  Double_t ruxmin = gPad->GetUxmin();
  Double_t ruymin = gPad->GetUymin();
  if (ktype == 3) {
    xorg = ruxmin;
    yorg = ruymin;
  } else {
    xorg = TMath::Max((Double_t)0,ruxmin);
    yorg = TMath::Min(TMath::Max((Double_t)0,ruymin),gPad->GetUymax());
  }
  
  // delta is the accuracy required in constructing the curve.
  // If it is zero then the routine calculates a value otherwise
  // it uses this value. (default is 0.0)
  
  delta         = 0.00055;
  maxiterations = 20;
  
  //       Scale data to the range 0-ratio_signs in X, 0-1 in Y
  //       where ratio_signs is the ratio between the number of changes
  //       of sign in Y divided by the number of changes of sign in X
  
  sxmin = x[0];
  sxmax = x[0];
  symin = y[0];
  symax = y[0];
  Double_t six   = 1;
  Double_t siy   = 1;
  for (i=1;i<npoints;i++) {
    if (i > 1) {
      if ((x[i]-x[i-1])*(x[i-1]-x[i-2]) < 0) six++;
      if ((y[i]-y[i-1])*(y[i-1]-y[i-2]) < 0) siy++;
    }
    if (x[i] < sxmin) sxmin = x[i];
    if (x[i] > sxmax) sxmax = x[i];
    if (y[i] < symin) symin = y[i];
    if (y[i] > symax) symax = y[i];
  }
  closed = 0;
  Double_t dx1n   = TMath::Abs(x[npoints-1]-x[0]);
  Double_t dy1n   = TMath::Abs(y[npoints-1]-y[0]);
  if (dx1n < 0.01*(sxmax-sxmin) && dy1n < 0.01*(symax-symin))  closed = 1;
  if (sxmin == sxmax) {
    xratio = 1;
  } else {
    if (six > 1) ratio_signs = siy/six;
    else         ratio_signs = 20;
    xratio = ratio_signs/(sxmax-sxmin);
  }
  if (symin == symax) yratio = 1;
  else                yratio = 1/(symax-symin);
  
  qlx[0] = x[0];
  qly[0] = y[0];
  for (i=0;i<npoints;i++) {
    x[i] = (x[i]-sxmin)*xratio;
    y[i] = (y[i]-symin)*yratio;
  }
  
  //          "finished" is minus one if we must draw a straight line from P(k-1)
  //          to P(k). "finished" is one if the last call to PaintPolyLine has < n2
  //          points. "finished" is zero otherwise. npt counts the X and Y
  //          coordinates in work . When npt=n2 a call to IPL is made.
  
  finished = 0;
  npt      = 1;
  k        = 1;
  
  //           Convert coordinates back to original system
  
  //           Separate the set of data points into arcs P(k-1),P(k).
  //           Calculate the direction cosines. first consider whether
  //           there is a continuous tangent at the endpoints.
  
  if (!closed) {
    if (x[0] != x[npoints-1] || y[0] != y[npoints-1]) goto L40;
    if (x[npoints-2] == x[npoints-1] && y[npoints-2] == y[npoints-1]) goto L40;
    if (x[0] == x[1] && y[0] == y[1]) goto L40;
  }
  flgic = kFALSE;
  flgis = kTRUE;
  
  //           flgic is true if the curve is open and false if it is closed.
  //           flgis is true in the main loop, but is false if there is
  //           a deviation from the main loop.
  
  km = npoints - 1;
  
  //           Calculate direction cosines at P(1) using P(N-1),P(1),P(2).
  
  goto L100;
L40:
  flgic = kTRUE;
  flgis = kFALSE;
  
  //           Skip excessive consecutive equal points.
  
L50:
  if (k >= npoints) {
    finished = 1;  // Prepare to clear out remaining short vectors before returning
    if (npt > 1) goto L310;
    goto L390;
  }
  k++;
  if (x[k-1] == x[k-2] && y[k-1] == y[k-2])  goto L50;
L60:
  km = k-1;
  if (k > npoints) {
    finished = 1;  // Prepare to clear out remaining short vectors before returning
    if (npt > 1) goto L310;
    goto L390;
  }
  if (k < npoints) goto L90;
  if (!flgic) { kp = 2; goto L130;}
  
L80:
  if (flgis) goto L150;
  
  //           Draw a straight line from P(k-1) to P(k).
  
  finished = -1;
  goto L170;
  
  //           Test whether P(k) is a cusp.
  
L90:
  if (x[k-1] == x[k] && y[k-1] == y[k]) goto L80;
L100:
  kp = k+1;
  goto L130;
  
  //           Branch if the next section of the curve begins at a cusp.
  
L110:
  if (!flgis) goto L50;
  
  //           Carry forward the direction cosines from the previous arc.
  
L120:
  co = ct;
  so = st;
  k++;
  goto L60;
  
  //           Calculate the direction cosines at P(k).  If k=1 then
  //           N-1 is used for k-1. If k=N then 2 is used for k+1.
  //           direction cosines at P(k) obtained from P(k-1),P(k),P(k+1).
  
L130:
  dx1 = x[k-1]  - x[km-1];
  dy1 = y[k-1]  - y[km-1];
  dk1 = dx1*dx1 + dy1*dy1;
  dx2 = x[kp-1] - x[k-1];
  dy2 = y[kp-1] - y[k-1];
  dk2 = dx2*dx2 + dy2*dy2;
  ctu = dx1*dk2 + dx2*dk1;
  stu = dy1*dk2 + dy2*dk1;
  xnt = ctu*ctu + stu*stu;
  
  //           If both ctu and stu are zero,then default.This can
  //           occur when P(k)=P(k+1). I.E. A loop.
  
  if (xnt < 1.E-25) {
    ctu = dy1;
    stu =-dx1;
    xnt = dk1;
  }
  //           Normalise direction cosines.
  
  ct = ctu/TMath::Sqrt(xnt);
  st = stu/TMath::Sqrt(xnt);
  if (flgis) goto L160;
  
  //           Direction cosines at P(k-1) obtained from P(k-1),P(k),P(k+1).
  
  w3    = 2*(dx1*dy2-dx2*dy1);
  co    = ctu+w3*dy1;
  so    = stu-w3*dx1;
  xnt   = 1/TMath::Sqrt(co*co+so*so);
  co    = co*xnt;
  so    = so*xnt;
  flgis = kTRUE;
  goto L170;
  
  //           Direction cosines at P(k) obtained from P(k-2),P(k-1),P(k).
  
L150:
  w3    = 2*(dx1*dy2-dx2*dy1);
  ct    = ctu-w3*dy2;
  st    = stu+w3*dx2;
  xnt   = 1/TMath::Sqrt(ct*ct+st*st);
  ct    = ct*xnt;
  st    = st*xnt;
  flgis = kFALSE;
  goto L170;
L160:
  if (k <= 1) goto L120;
  
  //           For the arc between P(k-1) and P(k) with direction cosines co,
  //           so and ct,st respectively, calculate the coefficients of the
  //           parametric cubic represented by X(t) and Y(t) where
  //           X(t)=xa*t**3 + xb*t**2 + co*t + xo
  //           Y(t)=ya*t**3 + yb*t**2 + so*t + yo
  
L170:
  xo = x[k-2];
  yo = y[k-2];
  dx = x[k-1] - xo;
  dy = y[k-1] - yo;
  
  //           Initialise the values of X(TI),Y(TI) in xt and yt respectively.
  
  xt = xo;
  yt = yo;
  if (finished < 0) {  // Draw a straight line between (xo,yo) and (xt,yt)
    xt += dx;
    yt += dy;
    goto L300;
  }
  c  = dx*dx+dy*dy;
  a  = co+ct;
  b  = so+st;
  r  = dx*a+dy*b;
  t  = c*6/(TMath::Sqrt(r*r+2*(7-co*ct-so*st)*c)+r);
  tsquare = t*t;
  tcube   = t*tsquare;
  xa = (a*t-2*dx)/tcube;
  xb = (3*dx-(co+a)*t)/tsquare;
  ya = (b*t-2*dy)/tcube;
  yb = (3*dy-(so+b)*t)/tsquare;
  
  //           If the curve is close to a straight line then use a straight
  //           line between (xo,yo) and (xt,yt).
  
  if (.75*TMath::Max(TMath::Abs(dx*so-dy*co),TMath::Abs(dx*st-dy*ct)) <= delta) {
    finished = -1;
    xt += dx;
    yt += dy;
    goto L300;
  }
  
  //           Calculate a set of values 0 == t(0).LTCT(1) <  ...  < t(M)=TC
  //           such that polygonal arc joining X(t(J)),Y(t(J)) (J=0,1,..M)
  //           is within the required accuracy of the curve
  
  tj = 0;
  u1 = ya*xb-yb*xa;
  u2 = yb*co-xb*so;
  u3 = so*xa-ya*co;
  
  //           Given t(J), calculate t(J+1). The values of X(t(J)),
  //           Y(t(J)) t(J) are contained in xt,yt and tj respectively.
  
L180:
  s  = t - tj;
  iw = -2;
  
  //           Define iw here later.
  
  p1 = (2*u1)*tj-u3;
  p2 = (u1*tj-u3)*3*tj+u2;
  p3 = 3*tj*ya+yb;
  p4 = (p3+yb)*tj+so;
  p5 = 3*tj*xa+xb;
  p6 = (p5+xb)*tj+co;
  
  //           Test D(tj,THETA). A is set to (Y(tj+s)-Y(tj))/s.b is
  //           set to (X(tj+s)-X(tj))/s.
  
  cc  = 0.8209285;
  err = 0.1209835;
L190:
  iw -= 2;
L200:
  a   = (s*ya+p3)*s+p4;
  b   = (s*xa+p5)*s+p6;
  
  //           Set z to PSI(D/delta)-cc.
  
  w1 = -s*(s*u1+p1);
  w2 = s*s*u1-p2;
  w3 = 1.5*w1+w2;
  
  //           Set the estimate of (THETA-tj)/s.Then set the numerator
  //           of the expression (EQUATION 4.4)/s. Then set the square
  //           of D(tj,tj+s)/delta. Then replace z by PSI(D/delta)-cc.
  
  if (w3 > 0) wsign = TMath::Abs(w1);
  else        wsign = -TMath::Abs(w1);
  sth = 0.5+wsign/(3.4*TMath::Abs(w1)+5.2*TMath::Abs(w3));
  z   = s*sth*(s-s*sth)*(w1*sth+w1+w2);
  z   = z*z/((a*a+b*b)*(delta*delta));
  z   = (z+2.642937)*z/((.3715652*z+3.063444)*z+.2441889)-cc;
  
  //           Branch if z has been calculated
  
  if (iw > 0) goto L250;
  if (z > err) goto L240;
  goto L220;
L210:
  iw -= 2;
L220:
  if (iw+2 == 0) goto L190;
  if (iw+2 >  0) goto L290;
  
  //           Last part of arc.
  
L230:
  xt = x[k-1];
  yt = y[k-1];
  s  = 0;
  goto L300;
  
  //           z(s). find a value of s where 0 <= s <= sb such that
  //           TMath::Abs(z(s)) < err
  
L240:
  kp = 0;
  c  = z;
  sb = s;
L250:
  theGraph->Zero(kp,0,sb,err,s,z,maxiterations);
  if (kp == 2) goto L210;
  if (kp > 2) {
    Error("Smooth", "Attempt to plot outside plot limits");
    goto L230;
  }
  if (iw > 0) goto L200;
  
  //           Set z=z(s) for s=0.
  
  if (iw < 0) {
    z  = -cc;
    iw = 0;
    goto L250;
  }
  
  //           Set z=z(s) for s=sb.
  
  z  = c;
  iw = 1;
  goto L250;
  
  //           Update tj,xt and yt.
  
L290:
  xt = xt + s*b;
  yt = yt + s*a;
  tj = s  + tj;
  
  //           Convert coordinates to original system
  
L300:
  qlx[npt] = sxmin + xt/xratio;
  qly[npt] = symin + yt/yratio;
  npt++;
  
  //           If a fill area must be drawn and if the banks LX and
  //           LY are too small they are enlarged in order to draw
  //           the filled area in one go.
  
  if (npt < banksize)  goto L320;
  if (drawtype >= 1000 || ktype > 1) {
    Int_t newsize = banksize + n2;
    std::vector<Double_t> qtemp(banksize);
    for (i=0;i<banksize;i++) qtemp[i] = qlx[i];
    qlx.resize(newsize);
    for (i=0;i<banksize;i++) qlx[i]   = qtemp[i];
    for (i=0;i<banksize;i++) qtemp[i] = qly[i];
    qly.resize(newsize);
    for (i=0;i<banksize;i++) qly[i] = qtemp[i];
    banksize = newsize;
    
    goto L320;
  }
  
  //           Store the graph
  
L310:
  if (drawtype >= 1000) {
    gPad->PaintFillArea(npt,qlx.data(),qly.data(), "B");
  } else {
    if (ktype > 1) {
      if (!loptx) {
        qlx[npt]   = qlx[npt-1];
        qlx[npt+1] = qlx[0];
        qly[npt]   = yorg;
        qly[npt+1] = yorg;
      } else {
        qlx[npt]   = xorg;
        qlx[npt+1] = xorg;
        qly[npt]   = qly[npt-1];
        qly[npt+1] = qly[0];
      }
    }
  }theGraph = new TGraph(npt,qlx.data(),qly.data());
  
  npt = 1;
  qlx[0] = sxmin + xt/xratio;
  qly[0] = symin + yt/yratio;
L320:
  if (finished > 0) goto L390;
  if (finished < 0) { finished = 0; goto L110;}
  if (s > 0) goto L180;
  goto L110;
  
  //           Convert coordinates back to original system
  
L390:
  for (i=0;i<npoints;i++) {
    x[i] = sxmin + x[i]/xratio;
    y[i] = symin + y[i]/yratio;
  }
  if (!suppress_info) std::cout <<
    Form("Info in <utils/interpolation_tools::Smooth>: Interpolated curve provided with %i new points.", theGraph->GetN() - npoints) << std::endl;
  return theGraph;
}
