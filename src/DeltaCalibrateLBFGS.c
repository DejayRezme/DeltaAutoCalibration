
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "lbfgs.h"

const Real PI = 3.14159265359;

typedef struct {
	Real delta_radius; // optional
	Real rod, rod2, rod3; // diagonal rod length, optionally different length for rods.
	Real xa; // x position of column a, ideally should be sin(240) * delta_radius
	Real ya; // y position of column a, ideally should be cos(240) * delta_radius
	Real xc; // x position of column c, ideally should be sin(0) * delta_radius = 0
	Real oa, ob, oc; // home offsets / limit switches
	int probeCount;
	Real probePoints[100][3]; // a list of probed points on the bed in column or motor coordinates
} DeltaParams;

// https://github.com/hercek/Marlin/blob/Marlin_v1/calibration.wxm
DeltaParams p = { // probe values from original maxima worksheet
	123.983, // delta radius
	250.590, 250.590, 250.590, // diagonal rod length
	-107.178, // xa
	62.015, // ya
	-41.687, // xc
	-9.49, -9.49, -9.49, // home offset a, b, c
	7,
	{ // probed points on bed Z = 0 in column coordinates
		{227.34, 227.34, 227.21},
		{258.65, 170.95, 171.55},
		{233.28, 233.28, 127.56},
		{170.94, 258.64, 169.76},
		{127.96, 233.20, 232.52},
		{170.82, 170.82, 258.62},
		{233.19, 127.95, 233.80}
	}
};

//DeltaParams p = { // Dejay's probe values from Kossel Mini, I think they are wrong?
//	105, // delta radius
//	216, 216, 216, // diagonal rod length
//	-0.866 * 105, // xa
//	-0.5 * 105, // ya
//	0, // xc
//	0, 0, 0, // home offset a, b, c
//	7,
//	{ // probed points on bed Z = 0 in column coordinates
//		{-272.3500,	-272.3500,	-272.3500},
//		{-246.3134,	-316.5323,	-316.5323},
//		{-267.0739,	-349.5644,	-267.0749},
//		{-317.4313,	-317.4300,	-247.2136},
//		{-350.5674,	-268.0746,	-268.0746},
//		{-317.5313,	-247.3136,	-317.5300},
//		{-267.2239,	-267.2239,	-349.7144}
//	}
//};

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

	for (i = 0; i < p.probeCount; i++) {
		const Real ta = p.probePoints[i][0];
		const Real tb = p.probePoints[i][1];
		const Real tc = p.probePoints[i][2];

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

static void printDeltaParams(DeltaParams *p) {

	printf("  rod: %.3f  delta_radius: %.3f  xa:%.2f  ya:%.2f  xc:%.2f  oa: %.2f  ob: %.2f  oc: %.2f  probeCount: %i\n  probePoints:", p->rod, p->delta_radius, p->xa, p->ya, p->xc, p->oa, p->ob, p->oc, p->probeCount);
	for (int i = 0; i < p->probeCount; i++)
		printf("(%.2f,%.2f,%.2f) ", p->probePoints[i][0], p->probePoints[i][1], p->probePoints[i][2]);
	printf("\n");
}

static Real random() {
	return (Real) rand() / RAND_MAX;
}

static Real random2(Real min, Real max) {
	return min + random() * (max - min);
}

static int calibrate() {

	int ret = 0;
	Real fx;
	Real *x = lbfgs_malloc(varCount);
	lbfgs_parameter_t param;

	if (x == NULL) {
		printf("ERROR: Failed to allocate a memory block for variables.\n");
		return 1;
	}

	/* Initialize the variables. */
	x[0] = p.rod;
	x[1] = p.xa;
	x[2] = p.ya;
	x[3] = p.xc;
	x[4] = p.oa;
	x[5] = p.ob;
	x[6] = p.oc;

	// randomize starting values to exhaust search space
	for (int i = 0; i < varCount; i++)
		x[i] *= random2(0.1, 2.0);

//	printf("\n==============================================================\n");
	printf("\n");
//	printf("\nStarting L-BFGS optimization with initial values:\n");
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
	printf("LBFGS error: %f   status: %s\n", fx, lbfgs_errorString(ret));
	if (ret == LBFGS_SUCCESS) {
		//	printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f, x[4] = %f, x[5] = %f, x[6] = %f\n  ", fx, x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
		//	printf("  rod: %f   delta_r: %f   oa: %f   ob: %f   oc: %f\n\n", x[0], sqrt(sqr(x[1]) + sqr(x[2])), x[4], x[5], x[6]);
		DeltaParams r = p;
		r.rod = x[0];
		r.delta_radius = sqrt(sqr(x[1]) + sqr(x[2]));
		r.xa = x[1];
		r.ya = x[2];
		r.xc = x[3];
		r.oa = x[4];
		r.ob = x[5];
		r.oc = x[6];
		printDeltaParams(&r);
	}

	lbfgs_free(x);
	return 0;
}

Real asqrt(Real a) {
//	if (a > 0)
		return sqrt(a);
//	else
//		return -sqrt(-a);
}

static DeltaParams* generateDummyDelta(DeltaParams* p) {

	Real psca = 1; // probe x/y scatter
	Real perr = 0.1; // probe z error
	Real terr = 1; // tower error
	Real x, y, z = 0;

	p->delta_radius = random2(50, 500);
	p->rod = p->rod2 = p->rod3 = p->delta_radius * random2(2.01, 2.3);
	p->xa = sin(240 * PI / 180) * p->delta_radius + random2(-terr, terr);
	p->ya = cos(240 * PI / 180) * p->delta_radius + random2(-terr, terr);
	p->xc = 0 + random2(-terr, terr);
	p->oa = -random2(100, 500);
	p->ob = p->oa + random2(-1, 1);
	p->oc = p->oa + random2(-1, 2);
	p->probeCount = 7 + (int) random2(6, 21);

	for (int i = 0; i < p->probeCount; i++) {

		if (i == 0) { // center point
			x = 0;
			y = 0;
		} else if (i < 7) {
			x = sin((i - 1) * 60 * PI / 180) * p->delta_radius + random(-psca, psca);
			y = cos((i - 1) * 60 * PI / 180) * p->delta_radius + random(-psca, psca);
		} else if (i < 13) {
			x = sin((i - 1 + 0.5) * 60 * PI / 180) * p->delta_radius * 0.5 + random(-psca, psca);
			y = cos((i - 1 + 0.5) * 60 * PI / 180) * p->delta_radius * 0.5 + random(-psca, psca);
		} else{
			do { // make sure to generate a x/y inside delta radius
				x = random2(-p->delta_radius, p->delta_radius);
				y = random2(-p->delta_radius, p->delta_radius);
			} while (x*x + y*y > p->delta_radius);
		}
//		printf("%f, %f\n", x, y);


		z = random2(-perr, perr); // generate a bit of probe error or bed level error for the Z value

		p->probePoints[i][0] = asqrt((z - sqr(p->ya) + 2 * y * p->ya - sqr(y) - sqr(p->xa) + 2 * x * p->xa - sqr(x) + sqr(p->rod))) - p->oa;
		p->probePoints[i][1] = asqrt((z - sqr(p->ya) + 2 * y * p->ya - sqr(y) - sqr(p->xa) - 2 * x * p->xa - sqr(x) + sqr(p->rod2))) - p->ob;
		p->probePoints[i][2] = asqrt((z - 4 * sqr(p->ya) - 4 * y * p->ya - sqr(y) - sqr(p->xc) + 2 * x * p->xc - sqr(x) + sqr(p->rod3))) - p->oc;
	}

	return p;
}

int main(int argc, char *argv[]) {

	// initialize random seed
	srand((unsigned)time(NULL));

	// generate randomized delta
	generateDummyDelta(&p);

	printf("Initial delta parameters:\n");
	printDeltaParams(&p);

	clock_t startClock = clock();
	for (int i = 0; i < 100; i++)
		calibrate();

	printf("\nUsed %0.3f seconds of CPU time. \n", (double) (clock() - startClock) / CLOCKS_PER_SEC);

	return 0;
}
