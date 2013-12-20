#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

/******************************************************************************/
/*                                                                            */
/*  Identical aligned parallel plates                                         */
/*                                           +---------------------+          */
/*  Hamilton and Morgan (1952), A-1         /                     /           */
/*  Howell catalog, C-11                   /                   b /            */
/*  Modest (2013), Configuration 38       /         a           /             */
/*                                       +---------------------+              */
/*  Arguments:                           |                                    */
/*    X - a/c                            |                                    */
/*    Y - b/c                          c |   +---------------------+          */
/*  where:                               |  /                     /           */
/*    a - length of plates               | /                   b /            */
/*    b - width of plates                |/         a           /             */
/*    c - separation distance            +---------------------+              */
/*                                                                            */
/******************************************************************************/
double identicalAlignedParallelRectangles(double X, double Y)
{
  double X2,Y2,sqrt1pX2,sqrt1pY2;
  X2 = X*X;
  Y2 = Y*Y;
  sqrt1pX2 = sqrt(1.0+X2);
  sqrt1pY2 = sqrt(1.0+Y2);
  return 2.0*( 0.5*log((1.0+X2)*(1.0+Y2)/(1+X2+Y2)) 
	       + X*sqrt1pY2*atan(X/sqrt1pY2)
	       + Y*sqrt1pX2*atan(Y/sqrt1pX2)
	       - X*atan(X)
	       - Y*atan(Y) )/(M_PI*X*Y);
}

/******************************************************************************/
/*                                                                            */
/*  Rectangles with a common edge at                                          */
/*  90 degree angle                             +                             */
/*                                             /|                             */
/*  Hamilton and Morgan (1952), ?             / |                             */
/*  Howell catalog, C-14                   l /  |                             */
/*  Modest (2013), Configuration 39         /   |                             */
/*                                         +    |         w                   */
/*  Arguments:                             | A2 +---------------------+       */
/*    H = h/l                            h |   /                     /        */
/*    W = w/l                              |  /                     /         */
/*  where:                                 | /         A1          / l        */
/*    l - length of plates                 |/                     /           */
/*    h,w - width of plates                +---------------------+            */
/*                                                                            */
/******************************************************************************/
double perpendicularCommonEdgeRectangles(double H, double W)
{
	double H2 = H*H;
	double W2 = W*W;
	double H2pW2 = H2+W2;
	double H2pW2p1 = H2pW2 + 1;
	double H2p1 = H2+1;
	double W2p1 = W2+1;
	double rootH2pW2 = sqrt(H2pW2);
	return (W*atan(1.0/W)+H*atan(1.0/H)-rootH2pW2*atan(1.0/rootH2pW2)
	        + 0.25*log((W2p1*H2p1/H2pW2p1)
			           *pow(W2*H2pW2p1/(W2p1*H2pW2),W2)
					   *pow(H2*H2pW2p1/(H2p1*H2pW2),H2)))/(M_PI*W);
}

/******************************************************************************/
/*                                                                            */
/*  Rectangles with a common edge at                                          */
/*  90 degree angle                             +                             */
/*                                             /|                             */
/*  Hamilton and Morgan (1952), ?             / |                             */
/*  Howell catalog, C-14                   l /  |                             */
/*  Modest (2013), Configuration 39         /   |                             */
/*                                         +    |         w                   */
/*  Arguments:                             | A2 +---------------------+       */
/*    l - length of plates               h |   /                     /        */
/*    h,w - width of plates                |  /                     /         */
/*                                         | /         A1          / l        */
/*                                         |/                     /           */
/*                                         +---------------------+            */
/*                                                                            */
/******************************************************************************/
double perpendicularCommonEdgeRectangles(double w, double h, double l)
{
	double H, W, H2, W2, H2pW2, H2pW2p1, H2p1, W2p1, rootH2pW2;
	if(w*l*h <= 0.0)
		return 0.0;
	H = h/l;
	W = w/l;
	H2 = H*H;
	W2 = W*W;
	H2pW2 = H2+W2;
	H2pW2p1 = H2pW2 + 1;
	H2p1 = H2+1;
	W2p1 = W2+1;
	rootH2pW2 = sqrt(H2pW2);
	return l*l*(W*atan(1.0/W)+H*atan(1.0/H)-rootH2pW2*atan(1.0/rootH2pW2)
	            + 0.25*log((W2p1*H2p1/H2pW2p1)
			               *pow(W2*H2pW2p1/(W2p1*H2pW2),W2)
					       *pow(H2*H2pW2p1/(H2p1*H2pW2),H2)))/M_PI;
}

/******************************************************************************/
/*                                                                            */
/*  Differential sphere to aligned plane rectangle                            */
/*                                                                            */
/*  (Sphere is on perpendicular ray from one corner of the rectangle, which   */
/*   is a 90 degree angle in Howell's formula)                                */
/*                                                                            */
/*  Hamilton and Morgan (1952), P-5                                           */
/*  Howell catalog, B-106                                                     */
/*                                                                            */
/*  Arguments:                                                                */
/*    a - height of plate                                                     */
/*    b - width of plate                                                      */
/*    c - separation distance                                                 */
/*                                                                            */
/*         /|                                                                 */
/*        / |                                                                 */
/*      b/  |_                                                                */
/*      /   | |                                                               */
/*     /    +-----------------------+                                         */
/*    |    /                     /_/                                          */
/*    |   /                       /      The sphere is placed at "o" and      */
/*   a|  /                       /       the target plane is "c" away and     */
/*    | /_                      /        "a" by "b" in size.                  */
/*    |/ /       c             /                                              */
/*    +-----------------------o                                               */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
double differentialSphereToAlignedPlaneRectangle(double a, double b, double c)
{
  double X = a/c;
  double Y = b/c;
  return 0.25*atan(X*Y/sqrt(1.0+X*X+Y*Y))/M_PI;
}

/******************************************************************************/
/*                                                                            */
/*  Differential planar element to aligned plane rectangle                    */
/*                                                                            */
/*  (Diffential element is aligned with one corner of the rectangle)          */
/*                                                                            */
/*  Hamilton and Morgan (1952), P-1                                           */
/*  Howell catalog, B-3                                                       */
/*                                                                            */
/*  Arguments:                                                                */
/*    a - length of plate                                                     */
/*    b - width of plate                                                      */
/*    c - separation distance                                                 */
/*                                                                            */
/******************************************************************************/
double differentialPlaneToParallelPlaneRectangle(double a, double b, double c)
{
  double A = a/c;
  double B = b/c;
  double rA,rB;
  rA = 1.0/sqrt(1.0+A*A);
  rB = 1.0/sqrt(1.0+B*B);
  return 0.5*(A*atan(B*rA)*rA + B*atan(A*rB)*rB)/M_PI;
}

/******************************************************************************/
/*                                                                            */
/*  Differential planar element to aligned finite cylinder                    */
/*                                                                            */
/*  (Diffential element is aligned with one end of the cylinder)              */
/*                                                                            */
/*  Hamilton and Morgan (1952), P-8                                           */
/*  Howell catalog, B-31                                                      */
/*                                                                            */
/*  Arguments:                                                                */
/*    r - radius of cylinder                                                  */
/*    l - length of cylinder                                                  */
/*    h - distance between plane and cylinder axis                            */
/*                                                                            */
/******************************************************************************/
double differentialPlaneToFiniteCylinder(double r, double l, double h)
{
  double L = l/r;
  double H = h/r;
  double Hm1 = H-1.0;
  double Hp1 = H+1.0;
  double X = Hp1*Hp1 + L*L;
  double Y = Hm1*Hm1 + L*L;
  return (atan(L/sqrt(Hp1*Hm1))/H
      + L*((X-2.0*H)*atan(sqrt(X*Hm1/(Y*Hp1)))/(H*sqrt(X*Y))
           - atan(sqrt(Hm1/Hp1))/H))/M_PI;
}
