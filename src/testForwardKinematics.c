//#include <math.h>
//#include <stdio.h>
//
//static double L, R;
//static double Ax, Ay, Bx, By, Cx, Cy, L2;
//
//#define SQ3    (sqrt(3))
//
//#define SIN_60 (SQ3/2)
//#define COS_60 (.5)
//
//static double sq(double x) { return x*x; }
//static double sqr(double x) { return x*x; }
//
//static void set_geometry(double r_, double l_)
//{
//    if(L == l_ && R == r_) return;
//
//    L = l_;
//    R = r_;
//
//    L2 = sq(L);
//
//    Ax = 0.0;
//    Ay = R;
//
//    Bx = -SIN_60 * R;
//    By = -COS_60 * R;
//
//    Cx = SIN_60 * R;
//    Cy = -COS_60 * R;
//}
//
//int kinematics_forward(double ta, double tb, double tc, double *x, double *y, double *z)
//{
//	double den = (By-Ay)*Cx-(Cy-Ay)*Bx;
//
//	double w1 = Ay*Ay + ta*ta; // n.b. assumption that Ax is 0 all through here
//	double w2 = Bx*Bx + By*By + tb*tb;
//	double w3 = Cx*Cx + Cy*Cy + tc*tc;
//
//	double a1 = (tb-ta)*(Cy-Ay)-(tc-ta)*(By-Ay);
//	double b1 = -((w2-w1)*(Cy-Ay)-(w3-w1)*(By-Ay))/2.0;
//
//	double a2 = -(tb-ta)*Cx+(tc-ta)*Bx;
//	double b2 = ((w2-w1)*Cx - (w3-w1)*Bx)/2.0;
//
//	// a*z^2 + b*z + c = 0
//	double a = a1*a1 + a2*a2 + den*den;
//	double b = 2*(a1*b1 + a2*(b2-Ay*den) - ta*den*den);
//	double c = (b2-Ay*den)*(b2-Ay*den) + b1*b1 + den*den*(ta*ta - L*L);
//
//	double discr = b*b - 4.0*a*c;
//	if (discr < 0) return -1; // non-existing point
//	*z = -0.5*(b+sqrt(discr))/a;
//	*x = (a1* *z + b1)/den;
//	*y = (a2* *z + b2)/den;
//
////	double dr = L;
////	double xa = -0.8660254 * dr, ya = -0.5 * dr, xc = 0;
////	double oa = ta, ob = tb, oc = tc;
////	double r = R;
//	// this should be correct formula for calculating z: from ta tb and tc
////	z=-(6*xa*ya*sqrt((-324*xa^2-81*tb^2+162*ta*tb-81*ta^2)*ya^4+((-72*xa^2-18*tb^2+36*ta*tb-18*ta^2)*xc^2-72*xa^4+
////	(-72*tc^2+(72*tb+72*ta)*tc-54*tb^2+36*ta*tb-54*ta^2+144*r^2)*xa^2+(-18*tb^2+36*ta*tb-18*ta^2)*tc^2+(18*tb^3-18*ta*tb^2-18*ta^2*tb+18*ta^3)*tc
////	-9*tb^4+18*ta*tb^3+(36*r^2-18*ta^2)*tb^2+(18*ta^3-72*r^2*ta)*tb-9*ta^4+36*r^2*ta^2)*ya^2+(-4*xa^2-tb^2+2*ta*tb-ta^2)*xc^4+(8*xa^4+
////	(-8*tc^2+(8*tb+8*ta)*tc-2*tb^2-4*ta*tb-2*ta^2)*xa^2+(-2*tb^2+4*ta*tb-2*ta^2)*tc^2+(2*tb^3-2*ta*tb^2-2*ta^2*tb+2*ta^3)*tc-tb^4+2*ta*tb^3+
////	(4*r^2-2*ta^2)*tb^2+(2*ta^3-8*r^2*ta)*tb-ta^4+4*r^2*ta^2)*xc^2+(((16*ta-16*tb)*tc+8*tb^2-8*ta^2)*xa^3+
////	((-4*tb^3+12*ta*tb^2+(16*r^2-12*ta^2)*tb+4*ta^3-16*r^2*ta)*tc+2*tb^4-4*ta*tb^3-8*r^2*tb^2+4*ta^3*tb-2*ta^4+8*r^2*ta^2)*xa)*xc-4*xa^6+
////	(-8*tc^2+(8*tb+8*ta)*tc-5*tb^2+2*ta*tb-5*ta^2)*xa^4+(-4*tc^4+(8*tb+8*ta)*tc^3+(-6*tb^2-12*ta*tb-6*ta^2+16*r^2)*tc^2+
////	(2*tb^3+6*ta*tb^2+(6*ta^2-16*r^2)*tb+2*ta^3-16*r^2*ta)*tc-tb^4+2*ta*tb^3+(4*r^2-6*ta^2)*tb^2+(2*ta^3+8*r^2*ta)*tb-ta^4+4*r^2*ta^2)*xa^2+
////	(-tb^2+2*ta*tb-ta^2)*tc^4+(2*tb^3-2*ta*tb^2-2*ta^2*tb+2*ta^3)*tc^3+(-tb^4-2*ta*tb^3+6*ta^2*tb^2-2*ta^3*tb-ta^4)*tc^2+
////	(2*ta*tb^4-2*ta^2*tb^3-2*ta^3*tb^2+2*ta^4*tb)*tc-ta^2*tb^4+2*ta^3*tb^3-ta^4*tb^2)+
////	((18*ta-18*tb)*xa*xc+(-36*tc-18*tb-18*ta)*xa^2-9*tb^3+9*ta*tb^2+9*ta^2*tb-9*ta^3)*ya^2+(2*ta-2*tb)*xa*xc^3+
////	((-4*tc+2*tb+2*ta)*xa^2-tb^3+ta*tb^2+ta^2*tb-ta^3)*xc^2+((2*tb-2*ta)*xa^3+((2*ta-2*tb)*tc^2+(2*ta^2-2*tb^2)*tc+2*tb^3-2*ta^3)*xa)*xc+
////	(4*tc-2*tb-2*ta)*xa^4+(-4*tc^3+(2*tb+2*ta)*tc^2+(2*tb^2+2*ta^2)*tc-tb^3-ta*tb^2-ta^2*tb-ta^3)*xa^2)/((72*xa^2+18*tb^2-36*ta*tb+18*ta^2)*ya^2+
////	(2*tb^2-4*ta*tb+2*ta^2)*xc^2+((8*tb-8*ta)*tc-4*tb^2+4*ta^2)*xa*xc+(8*tc^2+(-8*tb-8*ta)*tc+2*tb^2+4*ta*tb+2*ta^2)*xa^2),z=(6*xa*ya*sqrt(
////	(-324*xa^2-81*tb^2+162*ta*tb-81*ta^2)*ya^4+((-72*xa^2-18*tb^2+36*ta*tb-18*ta^2)*xc^2-72*xa^4+
////	(-72*tc^2+(72*tb+72*ta)*tc-54*tb^2+36*ta*tb-54*ta^2+144*r^2)*xa^2+(-18*tb^2+36*ta*tb-18*ta^2)*tc^2+(18*tb^3-18*ta*tb^2-18*ta^2*tb+18*ta^3)*tc
////	-9*tb^4+18*ta*tb^3+(36*r^2-18*ta^2)*tb^2+(18*ta^3-72*r^2*ta)*tb-9*ta^4+36*r^2*ta^2)*ya^2+(-4*xa^2-tb^2+2*ta*tb-ta^2)*xc^4+(8*xa^4+
////	(-8*tc^2+(8*tb+8*ta)*tc-2*tb^2-4*ta*tb-2*ta^2)*xa^2+(-2*tb^2+4*ta*tb-2*ta^2)*tc^2+(2*tb^3-2*ta*tb^2-2*ta^2*tb+2*ta^3)*tc-tb^4+2*ta*tb^3+
////	(4*r^2-2*ta^2)*tb^2+(2*ta^3-8*r^2*ta)*tb-ta^4+4*r^2*ta^2)*xc^2+(((16*ta-16*tb)*tc+8*tb^2-8*ta^2)*xa^3+
////	((-4*tb^3+12*ta*tb^2+(16*r^2-12*ta^2)*tb+4*ta^3-16*r^2*ta)*tc+2*tb^4-4*ta*tb^3-8*r^2*tb^2+4*ta^3*tb-2*ta^4+8*r^2*ta^2)*xa)*xc-4*xa^6+
////	(-8*tc^2+(8*tb+8*ta)*tc-5*tb^2+2*ta*tb-5*ta^2)*xa^4+(-4*tc^4+(8*tb+8*ta)*tc^3+(-6*tb^2-12*ta*tb-6*ta^2+16*r^2)*tc^2+
////	(2*tb^3+6*ta*tb^2+(6*ta^2-16*r^2)*tb+2*ta^3-16*r^2*ta)*tc-tb^4+2*ta*tb^3+(4*r^2-6*ta^2)*tb^2+(2*ta^3+8*r^2*ta)*tb-ta^4+4*r^2*ta^2)*xa^2+
////	(-tb^2+2*ta*tb-ta^2)*tc^4+(2*tb^3-2*ta*tb^2-2*ta^2*tb+2*ta^3)*tc^3+(-tb^4-2*ta*tb^3+6*ta^2*tb^2-2*ta^3*tb-ta^4)*tc^2+
////	(2*ta*tb^4-2*ta^2*tb^3-2*ta^3*tb^2+2*ta^4*tb)*tc-ta^2*tb^4+2*ta^3*tb^3-ta^4*tb^2)+
////	((18*tb-18*ta)*xa*xc+(36*tc+18*tb+18*ta)*xa^2+9*tb^3-9*ta*tb^2-9*ta^2*tb+9*ta^3)*ya^2+(2*tb-2*ta)*xa*xc^3+
////	((4*tc-2*tb-2*ta)*xa^2+tb^3-ta*tb^2-ta^2*tb+ta^3)*xc^2+((2*ta-2*tb)*xa^3+((2*tb-2*ta)*tc^2+(2*tb^2-2*ta^2)*tc-2*tb^3+2*ta^3)*xa)*xc+
////	(-4*tc+2*tb+2*ta)*xa^4+(4*tc^3+(-2*tb-2*ta)*tc^2+(-2*tb^2-2*ta^2)*tc+tb^3+ta*tb^2+ta^2*tb+ta^3)*xa^2)/((72*xa^2+18*tb^2-36*ta*tb+18*ta^2)*ya^2
////	+(2*tb^2-4*ta*tb+2*ta^2)*xc^2+((8*tb-8*ta)*tc-4*tb^2+4*ta^2)*xa*xc+(8*tc^2+(-8*tb-8*ta)*tc+2*tb^2+4*ta*tb+2*ta^2)*xa^2)
//
//	return 0;
//}
//
//int main(int argc, char **argv) {
//
//	double delta_r = 105;
//	double rod = 216;
//	set_geometry(delta_r, rod);
//	printf("delta_r: %f   rod: %f   center z:%f\n", delta_r, rod, -sqrt(sqr(rod) - sqr(delta_r)));
//
//	double q1 = 272.3, q2 = 272.3, q3 = 272.3;
////	q1 = q2 = q3 = 0;
//	double x, y, z;
//	kinematics_forward(q1, q2, q3, &x, &y, &z);
//
//	printf("x:%f   y:%f   z:%f", x, y, z);
//}
//
