
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "lbfgs.h"

const Real PI = 3.14159265359;

typedef struct {
	Real calibError;
	int calibStatus;
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
	0, 0, // calibError and status
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
//	0, 0, // calibError and status
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

inline Real solvePoly(Real a, Real b, Real c) {
	return (sqrt(b * b - 4 * a * c) + b) / (2 * a);
}

void inverseKinetmatic(Real *ta, Real *tb, Real *tc, const Real x, const Real y, const Real rod1, const Real rod2, const Real rod3, const Real xa, const Real ya, const Real xc, const Real xsa, const Real ysa, const Real xsb, const Real ysb, const Real xsc, const Real ysc, const Real oa, const Real ob, const Real oc) {

	// solve the polynomial for ta, tb, tc
	//	ta^2 + ta*(-2*oa*ysa^2+(-2*ya+2*y)*ysa-2*oa*xsa^2+(-2*xa+2*x)*xsa+2*oa) - oa^2*ysa^2+ya^2-2*y*ya+y^2-oa^2*xsa^2+xa^2-2*x*xa+x^2-rod1^2+oa^2 = 0
	//	tb^2 + tb*(-2*ob*ysb^2+(-2*ya+2*y)*ysb-2*ob*xsb^2+(2*xa+2*x)*xsb+2*ob) - ob^2*ysb^2+ya^2-2*y*ya+y^2-ob^2*xsb^2+xa^2+2*x*xa+x^2-rod2^2+ob^2 = 0
	//	tc^2 + tc*(-2*oc*ysc^2+(4*ya+2*y)*ysc-2*oc*xsc^2+(-2*xc+2*x)*xsc+2*oc) - oc^2*ysc^2+4*ya^2+4*y*ya+y^2-oc^2*xsc^2+xc^2-2*x*xc+x^2-rod3^2+oc^2 = 0

	Real xb = -xa;
	Real yb = ya;
	Real yc = -2*ya;
	*ta = solvePoly(1, -2*oa*sqr(ysa)+(2*y-2*ya)*ysa-2*oa*sqr(xsa)+(2*x-2*xa)*xsa+2*oa, -sqr(oa)*sqr(ysa)+sqr(ya)-2*y*ya+sqr(y)-sqr(oa)*sqr(xsa)+sqr(xa)-2*x*xa+sqr(x)-sqr(rod1)+sqr(oa));
	*tb = solvePoly(1, -2*ob*sqr(ysb)+(2*y-2*yb)*ysb-2*ob*sqr(xsb)+(2*x-2*xb)*xsb+2*ob, -sqr(ob)*sqr(ysb)+sqr(yb)-2*y*yb+sqr(y)-sqr(ob)*sqr(xsb)+sqr(xb)-2*x*xb+sqr(x)-sqr(rod1)+sqr(ob));
	*tc = solvePoly(1, -2*oc*sqr(ysc)+(2*y-2*yc)*ysc-2*oc*sqr(xsc)+(2*x-2*xc)*xsc+2*oc, -sqr(oc)*sqr(ysc)+sqr(yc)-2*y*yc+sqr(y)-sqr(oc)*sqr(xsc)+sqr(xc)-2*x*xc+sqr(x)-sqr(rod1)+sqr(oc));

//	*ta=sqrt(sqr(oa)*pow(ysa,4)+(2*oa*ya-2*oa*y)*pow(ysa,3)+
//	(sqr(ya)-2*y*ya+sqr(y)+2*sqr(oa)*sqr(xsa)+(2*oa*xa-2*oa*x)*xsa-sqr(oa))*sqr(ysa)+
//	((2*oa*sqr(xsa)+(2*xa-2*x)*xsa-2*oa)*ya+(-2*oa*sqr(xsa)+(2*x-2*xa)*xsa+2*oa)*y)*ysa-sqr(ya)+2*y*ya-sqr(y)+sqr(oa)*
//	pow(xsa,4)+(2*oa*xa-2*oa*x)*pow(xsa,3)+(sqr(xa)-2*x*xa+sqr(x)-sqr(oa))*sqr(xsa)+(2*oa*x-2*oa*xa)*xsa-sqr(xa)+2*x*xa-
//	sqr(x)+sqr(rod1)) + oa*sqr(ysa)+(ya-y)*ysa+oa*sqr(xsa)+(xa-x)*xsa-oa;

}

//Real eval(const Real rod1, const Real rod2, const Real rod3, const Real xa, const Real ya, const Real xc, const Real xsa, const Real ysa, const Real xsb, const Real ysb, const Real xsc, const Real ysc, const Real oa, const Real ob, const Real oc, Real *g, const int probeCount, const Real probePoints[][3]) {
void forwardKinematic(Real *x, Real *y, const Real ta, const Real tb, const Real tc, const Real rod1, const Real rod2, const Real rod3, const Real xa, const Real ya, const Real xc, const Real xsa, const Real ysa, const Real xsb, const Real ysb, const Real xsc, const Real ysc, const Real oa, const Real ob, const Real oc) {

	Real xb = -xa;
	Real yb = ya;
//	Real yc = -2*ya;

	Real r2a = -(1-sqr(xsa)-sqr(ysa)) * sqr(ta+oa) + sqr(rod1);
	Real r2b = -(1-sqr(xsa)-sqr(ysa)) * sqr(tb+ob) + sqr(rod2);
//	Real r2c = -(1-sqr(xsa)-sqr(ysa)) * sqr(tc+oc) + sqr(rod3);

	Real xaa=xa - ta * xsa;
	Real xbb=xb - tb * xsb;
//	Real xcc=xc - tc * xsc;
	Real yaa=ya - ta * ysa;
	Real ybb=yb - tb * ysb;
//	Real ycc=yc - tc * ysc;

	Real xaabb = xaa - xbb;
	Real yaabb = yaa - ybb;
	Real x2aabb = sqr(xaabb);
	Real y2aabb = sqr(yaabb);

	Real ra = sqrt(r2a);
	Real rb = sqrt(r2b);
	Real D = -y2aabb * (x2aabb+y2aabb - sqr(ra - rb)) * (x2aabb + y2aabb - sqr(ra + rb));
	Real sqrD = sqrt(D);

	*x = xaa - xaabb * (r2a - r2b + x2aabb + y2aabb) / 2 * (x2aabb + y2aabb) + sqrD / (x2aabb + y2aabb);
	*y = yaa - yaabb * (r2a - r2b + x2aabb + y2aabb) / 2 * (x2aabb + y2aabb) + (-xaabb / yaabb) * sqrD / (x2aabb + y2aabb);
}

inline static Real eval2(Real *g, const int probeCount, const Real probePoints[][3], const Real ta, const Real tb, const Real tc, const Real rod1, const Real rod2, const Real rod3, const Real xa, const Real ya, const Real xc, const Real xsa, const Real ysa, const Real xsb, const Real ysb, const Real xsc, const Real ysc, const Real oa, const Real ob, const Real oc) {

	int i;
	Real error = 0;
	for (i = 0; i < 13; i++)
		g[i] = 0;

	Real x, y;
	forwardKinematic(&x, &y, ta, tb, tc, rod1, rod2, rod3, xa, ya, xc, xsa, ysa, xsb, ysb, xsc, ysc, oa, ob, oc);

	for (i = 0; i < probeCount; i++) {
		const Real ta = probePoints[i][0];
		const Real tb = probePoints[i][1];
		const Real tc = probePoints[i][2];

		Real D = sqr(tc*ysc+2*ya+y)+sqr(tb*ysb-ya+y)+sqr(ta*ysa-ya+y)+sqr(tc+oc)*(-sqr(ysa)-sqr(xsa)+1)+sqr(tb+ob)*(-sqr(ysa)-sqr(xsa)+1)+sqr(ta+oa)*(-sqr(ysa)-sqr(xsa)+1)+sqr(tc*xsc-xc+x)+sqr(tb*xsb+xa+x)+sqr(ta*xsa-xa+x)-sqr(rod3)-sqr(rod2)-sqr(rod1);

		// error function
		error += sqr(D);
		//gradient rod1
		g[0]  += -4*rod1*(D);
		//gradient rod2
		g[1]  += -4*rod2*(D);
		//gradient rod3
		g[2]  += -4*rod3*(D);

		//gradient xa
		g[3]  += 2*(2*(tb*xsb+xa+x)-2*(ta*xsa-xa+x))*(D);
		//gradient ya
		g[4]  += 2*(4*(tc*ysc+2*ya+y)-2*(tb*ysb-ya+y)-2*(ta*ysa-ya+y))*(D);
		//gradient xc

		g[5]  += -4*(tc*xsc-xc+x)*(D);
		//gradient xsa
		g[6]  += 2*(2*ta*(ta*xsa-xa+x)-2*sqr(tc+oc)*xsa-2*sqr(tb+ob)*xsa-2*sqr(ta+oa)*xsa)*(D);
		//gradient ysa
		g[7]  += 2*(2*ta*(ta*ysa-ya+y)-2*sqr(tc+oc)*ysa-2*sqr(tb+ob)*ysa-2*sqr(ta+oa)*ysa)*(D);
		//gradient xsb
		g[8]  += 4*tb*(tb*xsb+xa+x)*(D);
		//gradient ysb
		g[9]  += 4*tb*(tb*ysb-ya+y)*(D);
		//gradient xsc
		g[10]  += 4*tc*(tc*xsc-xc+x)*(D);
		//gradient ysc
		g[11]  += 4*tc*(tc*ysc+2*ya+y)*(D);

		//gradient oa
		g[12]  += 4*(ta+oa)*(-sqr(ysa)-sqr(xsa)+1)*(D);
		//gradient ob
		g[13] += 4*(tb+ob)*(-sqr(ysa)-sqr(xsa)+1)*(D);
		//gradient oc
		g[14] += 4*(tc+oc)*(-sqr(ysa)-sqr(xsa)+1)*(D);
	}

	return error;
}

inline static Real eval(Real r, Real xa, Real ya, Real xc, Real oa, Real ob, Real oc, Real *g) {

	int i;
	Real value = 0;
	for (i = 0; i < varCount; i++)
		g[i] = 0;

	for (i = 0; i < p.probeCount; i++) {
		const Real ta = p.probePoints[i][0];
		const Real tb = p.probePoints[i][1];
		const Real tc = p.probePoints[i][2];

//		Real x = -(sqr(tb)+2*ob*tb-sqr(ta)-2*oa*ta+sqr(ob)-sqr(oa))/(4*xa);
//		Real y = -(6*ya^2+2*xc^2-4*x*xc-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya);

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

	printf("  rod: %.3f  delta_radius: %.3f  xa:%.3f  ya:%.3f  xc:%.3f  oa: %.3f  ob: %.3f  oc: %.3f  probeCount: %i\n", p->rod, p->delta_radius, p->xa, p->ya, p->xc, p->oa, p->ob, p->oc, p->probeCount);
	if (p->probeCount > 0) {
		printf("probePoints: ");
		for (int i = 0; i < p->probeCount; i++)
			printf("(%.2f,%.2f,%.2f) ", p->probePoints[i][0], p->probePoints[i][1], p->probePoints[i][2]);
		printf("\n");
	}
}

static Real random() {
	return (Real) rand() / RAND_MAX;
}

static Real random2(Real min, Real max) {
	return min + random() * (max - min);
}

static int shmequal(Real a, Real b) {
	return abs(a - b) < 0.000001;
}

#define REPEATS 1000
DeltaParams *results[REPEATS];
int resultsFailed = 0;
int resultsConverged = 0;

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
//	printf("\n");
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
//	printf("LBFGS error: %f   status: %s\n", fx, lbfgs_errorString(ret));
	if (ret == LBFGS_SUCCESS) {
		resultsConverged++;
		//	printf("  fx = %f, x[0] = %f, x[1] = %f, x[2] = %f, x[3] = %f, x[4] = %f, x[5] = %f, x[6] = %f\n  ", fx, x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
		//	printf("  rod: %f   delta_r: %f   oa: %f   ob: %f   oc: %f\n\n", x[0], sqrt(sqr(x[1]) + sqr(x[2])), x[4], x[5], x[6]);
		DeltaParams *r = malloc(sizeof(DeltaParams));
		r->calibError = fx;
		r->calibStatus = ret;
		r->rod = x[0];
		r->delta_radius = sqrt(sqr(x[1]) + sqr(x[2]));
		r->xa = x[1];
		r->ya = x[2];
		r->xc = x[3];
		r->oa = x[4];
		r->ob = x[5];
		r->oc = x[6];
//		r->probeCount =

		for (int i = 0; i < REPEATS; i++) {
			if (results[i] == NULL) {
				results[i] = r;
				break;
			} else if (shmequal(results[i]->calibError, r->calibError) && shmequal(results[i]->rod, r->rod)) {
				break;
			} else if (results[i]->calibError > r->calibError) {
				for (int o = REPEATS - 1; o > i; o--) // make space and bubble up
					results[o] = results[o-1];
				results[i] = r;
				break;
			}
		}
	} else {
		resultsFailed++;
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
	Real oerr = 1; // endstop offset error
	Real x, y, z = 0;

	p->delta_radius = random2(50, 500);
	p->rod = p->rod2 = p->rod3 = p->delta_radius * random2(2.01, 2.3);
	p->xa = sin(240 * PI / 180) * p->delta_radius + random2(-terr, terr);
	p->ya = cos(240 * PI / 180) * p->delta_radius + random2(-terr, terr);
	p->xc = 0 + random2(-terr, terr);
	p->oa = random2(100, 500);
	p->ob = p->oa + random2(-oerr, oerr);
	p->oc = p->oa + random2(-oerr, oerr);
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


//		z = random2(-perr, perr); // generate a bit of probe error or bed level error for the Z value
//		z += x * 0.01; // add some tilt

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
	for (int i = 0; i < REPEATS; i++)
		calibrate();

	printf("\nCalibration finished - used %0.3f seconds of CPU time. %i results converged and %i failed \n\n", (double) (clock() - startClock) / CLOCKS_PER_SEC, resultsConverged, resultsFailed);

	for (int i = 0; i < REPEATS; i++)
		if (results[i] != NULL) {
			printf("Solution %i error: %f (%e)\n", i, results[i]->calibError, results[i]->calibError);
			printDeltaParams(results[i]);
		}

	return 0;
}
