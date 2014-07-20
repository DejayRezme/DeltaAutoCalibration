
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "lbfgs.h"

const Real PI = 3.14159265359;

// probe values from original maxima worksheet
Real rod = 150.590;
Real xa = -207.178;
Real ya = 62.015;
Real xc = -41.687;
Real oa = -9.49, ob = 9.49, oc = 9.49;

const Real m[7][3] = {
		{227.34, 227.34, 227.21},
		{258.65, 170.95, 171.55},
		{233.28, 233.28, 127.56},
		{170.94, 258.64, 169.76},
		{127.96, 233.20, 232.52},
		{170.82, 170.82, 258.62},
		{233.19, 127.95, 233.80}
};

//// Dejay's probe values from Kossel Mini, I think they are wrong?
//Real rod = 216;
//Real delta_r = 105;
//Real xa = -0.866 * 105;
//Real ya = -0.5 * 105;
//Real xc = 0;
//Real oa = 0, ob = 0, oc = 0;
//
//const Real m[7][3] = {
//		{-272.3500,	-272.3500,	-272.3500},
//		{-246.3134,	-316.5323,	-316.5323},
//		{-267.0739,	-349.5644,	-267.0749},
//		{-317.4313,	-317.4300,	-247.2136},
//		{-350.5674,	-268.0746,	-268.0746},
//		{-317.5313,	-247.3136,	-317.5300},
//		{-267.2239,	-267.2239,	-349.7144}
//};


static const int sampleCount = 7;

inline Real max(Real a, Real b) {
	return a > b ? a : b;
}

inline Real sqr(Real a) {
	return a*a;
}

static const int varCount = 7;

inline static Real eval(Real r, Real xa, Real ya, Real xc, Real oa, Real ob, Real oc, Real *g) {

	int i;
	Real value = 0;
	for (i = 0; i < varCount; i++)
		g[i] = 0;

	for (i = 0; i < sampleCount; i++) {
		const Real ta = m[i][0];
		const Real tb = m[i][1];
		const Real tc = m[i][2];

		Real tb22obtbta = (sqr(tb)+2*ob*tb-sqr(ta)-2*oa*ta+sqr(ob)-sqr(oa));
		Real yaxctb22ob = (6*sqr(ya)+2*sqr(xc)+(tb22obtbta*xc)/xa-2*sqr(xa)+2*sqr(tc)+4*oc*tc-sqr(tb)-2*ob*tb-sqr(ta)-2*oa*ta+2*sqr(oc)-sqr(ob)-sqr(oa));
		Real yayaxctb22 = (2*ya-yaxctb22ob/(12*ya));
		Real yayaxctb22sq = sqr(yayaxctb22);
		Real xctbobtb = (-xc-tb22obtbta/(4*xa));
		Real xctbobtb2tc = sqr(xctbobtb)+sqr(tc+oc)-sqr(r);

		Real error = sqr(yayaxctb22sq+xctbobtb2tc);
		Real diff_r = -4*r*(yayaxctb22sq+xctbobtb2tc);
		Real diff_xa = 2*((tb22obtbta*xctbobtb)/(2*sqr(xa))-((-(tb22obtbta*xc)/sqr(xa)-4*xa)*yayaxctb22)/(6*ya))*(yayaxctb22sq+xctbobtb2tc);
		Real diff_ya = 4*(yaxctb22ob/(12*sqr(ya))+1)*yayaxctb22*(yayaxctb22sq+xctbobtb2tc);
		Real diff_xc = 2*(-((4*xc+tb22obtbta/xa)*yayaxctb22)/(6*ya)-2*xctbobtb)*(yayaxctb22sq+xctbobtb2tc);
		Real diff_oa = 2*(-((((-2*ta-2*oa)*xc)/xa-2*ta-2*oa)*yayaxctb22)/(6*ya)-((-2*ta-2*oa)*xctbobtb)/(2*xa))*(yayaxctb22sq+xctbobtb2tc);
		Real diff_ob = 2*(-((((2*tb+2*ob)*xc)/xa-2*tb-2*ob)*yayaxctb22)/(6*ya)-((2*tb+2*ob)*xctbobtb)/(2*xa))*(yayaxctb22sq+xctbobtb2tc);
		Real diff_oc = 2*(2*(tc+oc)-((4*tc+4*oc)*yayaxctb22)/(6*ya))*(yayaxctb22sq+xctbobtb2tc);

		value += error;
		g[0] += diff_r;
		g[1] += diff_xa;
		g[2] += diff_ya;
		g[3] += diff_xc;
		g[4] += diff_oa;
		g[5] += diff_ob;
		g[6] += diff_oc;

//		Real error = ((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2)^2;
//		Real diff_r = -4*r*((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		Real diff_xa = 2*(((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))/(2*xa^2)-((-((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa^2-4*xa)*
//		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya))*(
//		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		Real diff_ya = 4*((6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya^2)+1)*
//		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))*(
//		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		Real diff_xc = 2*(-((4*xc+(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/xa)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya)
//		-2*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))*(
//		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		Real diff_oa = 2*(-((((-2*ta-2*oa)*xc)/xa-2*ta-2*oa)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya)-
//		((-2*ta-2*oa)*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))/(2*xa))*((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		Real diff_ob = 2*(-((((2*tb+2*ob)*xc)/xa-2*tb-2*ob)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya)
//		-((2*tb+2*ob)*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))/(2*xa))*((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		Real diff_oc = 2*(2*(tc+oc)-((4*tc+4*oc)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya))*(
//		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
//		(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);

	}

	return value;
}

static Real evaluate(void *instance, const Real *x, Real *g, const int n, const Real step) {

	Real fx = eval(x[0], x[1], x[2], x[3], x[4], x[5], x[6], g);
	return fx;
}

static int progress(void *instance, const Real *x, const Real *g, const Real fx, const Real xnorm, const Real gnorm, const Real step, int n, int k, int ls) {
	printf("Iteration %d:\n", k);
	printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f, x[4] = %f, x[5] = %f, x[6] = %f\n", fx, x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
	printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	printf("\n");
	return 0;
}



int calibrate() {

	int ret = 0;
	Real fx;
	Real *x = lbfgs_malloc(varCount);
	lbfgs_parameter_t param;

	if (x == NULL) {
		printf("ERROR: Failed to allocate a memory block for variables.\n");
		return 1;
	}

	/* Initialize the variables. */
	x[0] = rod;
	x[1] = xa;
	x[2] = ya;
	x[3] = xc;
	x[4] = oa;
	x[5] = ob;
	x[6] = oc;

	for (int i = 0; i < varCount; i++) x[i] *= (0.2 + ((Real) rand()) / RAND_MAX * 1.8);

//	printf("\n==============================================================\n");
//	printf("Starting L-BFGS optimization with initial values:\n");
//	printf("  x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f, x[4] = %f, x[5] = %f, x[6] = %f\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

	/* Initialize the parameters for the L-BFGS optimization. */
	lbfgs_parameter_init(&param);
//	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
//	param.epsilon = 0.01;

	/*
	 Start the L-BFGS optimization; this will invoke the callback functions
	 evaluate() and progress() when necessary.
	 */
	ret = lbfgs(varCount, x, &fx, evaluate, NULL, NULL, &param);

	/* Report the result. */
	printf("L-BFGS optimization terminated with status = %s\n", lbfgs_errorString(ret));
	printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f, x[4] = %f, x[5] = %f, x[6] = %f\n", fx, x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
	printf("  rod: %f   delta_r: %f   oa: %f   ob: %f   oc: %f\n\n", x[0], sqrt(sqr(x[1]) + sqr(x[2])), x[4], x[5], x[6]);

	lbfgs_free(x);
	return 0;
}

int main(int argc, char *argv[]) {

	clock_t startClock = clock();

	srand((unsigned)time(NULL));
	for (int i = 0; i < 1000; i++)
		calibrate();

	printf("Used %0.3f seconds of CPU time. \n", (double) (clock() - startClock) / CLOCKS_PER_SEC);

	return 0;
}
