//#include <stdio.h>
//#include <math.h>
//#include "lbfgs.h"
//
//const char* lbfgsErrorString(int error) {
//	if (error == LBFGS_SUCCESS) return "LBFGS_SUCCESS";
//	if (error == LBFGS_CONVERGENCE) return "LBFGS_CONVERGENCE";
//	if (error == LBFGS_STOP) return "LBFGS_STOP";
//	if (error == LBFGS_ALREADY_MINIMIZED) return "LBFGS_ALREADY_MINIMIZED - The initial variables already minimize the objective function.";
//	if (error == LBFGSERR_UNKNOWNERROR) return "LBFGSERR_UNKNOWNERROR";
//	if (error == LBFGSERR_LOGICERROR) return "LBFGSERR_LOGICERROR";
//	if (error == LBFGSERR_OUTOFMEMORY) return "LBFGSERR_OUTOFMEMORY";
//	if (error == LBFGSERR_CANCELED) return "LBFGSERR_CANCELED";
//	if (error == LBFGSERR_INVALID_N) return "LBFGSERR_INVALID_N - Invalid number of variables specified.";
//	if (error == LBFGSERR_INVALID_N_SSE) return "LBFGSERR_INVALID_N_SSE - Invalid number of variables (for SSE) specified.";
//	if (error == LBFGSERR_INVALID_X_SSE) return "LBFGSERR_INVALID_X_SSE - The array x must be aligned to 16 (for SSE).";
//	if (error == LBFGSERR_INVALID_EPSILON) return "LBFGSERR_INVALID_EPSILON - Invalid parameter lbfgs_parameter_t::epsilon specified.";
//	if (error == LBFGSERR_INVALID_TESTPERIOD) return "LBFGSERR_INVALID_TESTPERIOD - Invalid parameter lbfgs_parameter_t::past specified.";
//	if (error == LBFGSERR_INVALID_DELTA) return "LBFGSERR_INVALID_DELTA - Invalid parameter lbfgs_parameter_t::delta specified.";
//	if (error == LBFGSERR_INVALID_LINESEARCH) return "LBFGSERR_INVALID_LINESEARCH - Invalid parameter lbfgs_parameter_t::linesearch specified.";
//	if (error == LBFGSERR_INVALID_MINSTEP) return "LBFGSERR_INVALID_MINSTEP - Invalid parameter lbfgs_parameter_t::max_step specified.";
//	if (error == LBFGSERR_INVALID_MAXSTEP) return "LBFGSERR_INVALID_MAXSTEP - Invalid parameter lbfgs_parameter_t::max_step specified.";
//	if (error == LBFGSERR_INVALID_FTOL) return "LBFGSERR_INVALID_FTOL - Invalid parameter lbfgs_parameter_t::ftol specified.";
//	if (error == LBFGSERR_INVALID_WOLFE) return "LBFGSERR_INVALID_WOLFE - Invalid parameter lbfgs_parameter_t::wolfe specified.";
//	if (error == LBFGSERR_INVALID_GTOL) return "LBFGSERR_INVALID_GTOL - Invalid parameter lbfgs_parameter_t::gtol specified.";
//	if (error == LBFGSERR_INVALID_XTOL) return "LBFGSERR_INVALID_XTOL - Invalid parameter lbfgs_parameter_t::xtol specified.";
//	if (error == LBFGSERR_INVALID_MAXLINESEARCH) return "LBFGSERR_INVALID_MAXLINESEARCH - Invalid parameter lbfgs_parameter_t::max_linesearch specified.";
//	if (error == LBFGSERR_INVALID_ORTHANTWISE) return "LBFGSERR_INVALID_ORTHANTWISE - Invalid parameter lbfgs_parameter_t::orthantwise_c specified.";
//	if (error == LBFGSERR_INVALID_ORTHANTWISE_START) return "LBFGSERR_INVALID_ORTHANTWISE_START - Invalid parameter lbfgs_parameter_t::orthantwise_start specified.";
//	if (error == LBFGSERR_INVALID_ORTHANTWISE_END) return "LBFGSERR_INVALID_ORTHANTWISE_END - Invalid parameter lbfgs_parameter_t::orthantwise_end specified.";
//	if (error == LBFGSERR_OUTOFINTERVAL) return "LBFGSERR_OUTOFINTERVAL - The line-search step went out of the interval of uncertainty.";
//	if (error == LBFGSERR_INCORRECT_TMINMAX) return "LBFGSERR_INCORRECT_TMINMAX - A logic error occurred; alternatively, the interval of uncertainty became too small.";
//	if (error == LBFGSERR_ROUNDING_ERROR) return "LBFGSERR_ROUNDING_ERROR - A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.";
//	if (error == LBFGSERR_MINIMUMSTEP) return "LBFGSERR_MINIMUMSTEP - The line-search step became smaller than lbfgs_parameter_t::min_step.";
//	if (error == LBFGSERR_MAXIMUMSTEP) return "LBFGSERR_MAXIMUMSTEP - The line-search step became larger than lbfgs_parameter_t::max_step.";
//	if (error == LBFGSERR_MAXIMUMLINESEARCH) return "LBFGSERR_MAXIMUMLINESEARCH - The line-search routine reaches the maximum number of evaluations.";
//	if (error == LBFGSERR_MAXIMUMITERATION) return "LBFGSERR_MAXIMUMITERATION - The algorithm routine reaches the maximum number of iterations.";
//	if (error == LBFGSERR_WIDTHTOOSMALL) return "LBFGSERR_WIDTHTOOSMALL - Relative width of the interval of uncertainty is at most lbfgs_parameter_t::xtol.";
//	if (error == LBFGSERR_INVALIDPARAMETERS) return "LBFGSERR_INVALIDPARAMETERS - A logic error (negative line-search step) occurred.";
//	if (error == LBFGSERR_INCREASEGRADIENT) return "LBFGSERR_INCREASEGRADIENT - The current search direction increases the objective function value.";
//	return "UNKNOWN!";
//}
//
//const Real PI = 3.14159265359;
//
//
////Real rod = 216;
////Real idelta_r = 105;
////Real xa = sin(240 * PI / 180) * idelta_r;
////Real ya = cos(240 * PI / 180) * idelta_r;
////Real xc = 0;
////Real oa = 0, ob = 0, oc = 0;
////
////const Real m[7][3] = {
////		{-272.3500,	-272.3500,	-272.3500},
////		{-246.3134,	-316.5323,	-316.5323},
////		{-267.0739,	-349.5644,	-267.0749},
////		{-317.4313,	-317.4300,	-247.2136},
////		{-350.5674,	-268.0746,	-268.0746},
////		{-317.5313,	-247.3136,	-317.5300},
////		{-267.2239,	-267.2239,	-349.7144}
////};
//
//
//
//Real rod = 250.590;
//Real xa = -107.178;
//Real ya = -62.015;
//Real xc = -1.687;
//Real oa = 9.49, ob = 9.49, oc = 9.49;
//
//const Real m[7][3] = {
//		{227.34, 227.34, 227.21},
//		{258.65, 170.95, 171.55},
//		{233.28, 233.28, 127.56},
//		{170.94, 258.64, 169.76},
//		{127.96, 233.20, 232.52},
//		{170.82, 170.82, 258.62},
//		{233.19, 127.95, 233.80}
//};
//
//static const int sampleCount = 7;
//
//
//class objective_function {
//protected:
//	Real *m_x;
//
//public:
//
//	objective_function() :
//			m_x(NULL) {
//	}
//
//	virtual ~objective_function() {
//		if (m_x != NULL) {
//			lbfgs_free(m_x);
//			m_x = NULL;
//		}
//	}
//
//	int run() {
//		Real fx;
//		Real *m_x = lbfgs_malloc(varCount);
//
////		(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2=0;
//
//		if (m_x == NULL) {
//			printf("ERROR: Failed to allocate a memory block for variables.\n");
//			return 1;
//		}
//
//		/* Initialize the variables. */
//		m_x[0] = rod;
//		m_x[1] = xa;
//		m_x[2] = ya;
//		m_x[3] = xc;
//		m_x[4] = oa;
//		m_x[5] = ob;
//		m_x[6] = oc;
//
//
//		static lbfgs_parameter_t params = {
//		    6, 1e-5, 0, 1e-5,
//		    0, LBFGS_LINESEARCH_DEFAULT, 40,
//		    1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
//		    0.0, 0, -1,
//		};
//
//		/*
//		 Start the L-BFGS optimization; this will invoke the callback functions
//		 evaluate() and progress() when necessary.
//		 */
//		int ret = lbfgs(varCount, m_x, &fx, _evaluate, _progress, this, &params);
//
//		/* Report the result. */
//		printf("L-BFGS optimization terminated with status = %s\n", lbfgsErrorString(ret));
//		printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);
//
//		Real delta_r = sqrt(sqr(m_x[1]) + sqr(m_x[2]));
//		printf("==============================================================\n");
//		printf("rod: %f   delta_r: %f   oa: %f   ob: %f   oc: %f\n", m_x[0], delta_r, m_x[4], m_x[5], m_x[6]);
//
//		return ret;
//	}
//
//protected:
//	static Real _evaluate(void *instance, const Real *x, Real *g, const int n,
//			const Real step) {
//		return reinterpret_cast<objective_function*>(instance)->evaluate(x, g,
//				n, step);
//	}
//
//	inline static Real max(Real a, Real b) {
//		return a > b ? a : b;
//	}
//
//	inline static Real sqr(Real a) {
//		return a*a;
//	}
//
//	static const int varCount = 7;
//
//	Real evaluate(const Real *x, Real *g, const int n, const Real step) {
//
////		Real fx = 0.0;
////
////		for (int i = 0; i < n; i += 2) {
////			Real t1 = 1.0 - x[i];
////			Real t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
////			g[i + 1] = 20.0 * t2;
////			g[i] = -2.0 * (x[i] * g[i + 1] + t1);
////			fx += t1 * t1 + t2 * t2;
////		}
////		return fx;
//
//		Real fx = eval(x[0], x[1], x[2], x[3], x[4], x[5], x[6], g);
////		Real d[varCount];
////		for (int i = 0; i < varCount; i++)
////			d[i] = x[i];
////		const Real difstep = max(0.0000000001, step);
////		for (int i = 0; i < varCount; i++) {
////
////			d[i] += difstep;
////			g[i] = (eval(d[0], d[1], d[2], d[3], d[4], d[5], d[6], g) - fx) / difstep;
////			d[i] = x[i];
////		}
////		printf("fx: %f   step: %f   x[0] = %f    x[1] = %f          g[0] = %f    g[1] = %f\n", fx, step, x[0], x[1], g[0], g[1]);
//
//		return fx;
//
////		for (int i = 0; i < varCount; i++)
////			g[i] = 0;
////		g[0] = 2*(x[0] - 2);
////		g[1] = 2*(x[1] - 3);
////		return sqr(x[0] - 2) + sqr(x[1] - 3);
//	}
//
//	inline static Real eval(Real r, Real xa, Real ya, Real xc, Real oa, Real ob, Real oc, Real *g) {
//
//		Real value = 0;
//		for (int i = 0; i < varCount; i++)
//			g[i] = 0;
//
//		for (int i = 0; i < sampleCount; i++) {
//			const Real ta = m[i][0];
//			const Real tb = m[i][1];
//			const Real tc = m[i][2];
////			printf("m[%i][] = %f, %f, %f\n", i, m[i][0], m[i][1], m[i][2]);
////			Real error = sqr(sqr(2 * ya - (6 * sqr(ya) + 2 * sqr(xc) + ((sqr(tb) + 2 * ob * tb - sqr(ta) - 2 * oa * ta + sqr(ob) - sqr(oa)) * xc) / xa - 2 * sqr(xa) + 2 * sqr(tc) + 4 * oc * tc - sqr(tb) - 2 * ob * tb - sqr(ta) - 2 * oa * ta + 2 * sqr(oc) - sqr(ob) - sqr(oa)) / (12 * ya)) + sqr(
////									-xc - (sqr(tb) + 2 * ob * tb - sqr(ta) - 2 * oa * ta + sqr(ob) - sqr(oa)) / (4 * xa)) + sqr(tc + oc) - sqr(r));
//
//			Real tb22obtbta = (sqr(tb)+2*ob*tb-sqr(ta)-2*oa*ta+sqr(ob)-sqr(oa));
//			Real yaxctb22ob = (6*sqr(ya)+2*sqr(xc)+(tb22obtbta*xc)/xa-2*sqr(xa)+2*sqr(tc)+4*oc*tc-sqr(tb)-2*ob*tb-sqr(ta)-2*oa*ta+2*sqr(oc)-sqr(ob)-sqr(oa));
//			Real yayaxctb22 = (2*ya-yaxctb22ob/(12*ya));
//			Real yayaxctb22sq = sqr(yayaxctb22);
//			Real xctbobtb = (-xc-tb22obtbta/(4*xa));
//			Real xctbobtb2tc = sqr(xctbobtb)+sqr(tc+oc)-sqr(r);
//
//			Real error = sqr(yayaxctb22sq+xctbobtb2tc);
//			Real diff_r = -4*r*(yayaxctb22sq+xctbobtb2tc);
//			Real diff_xa = 2*((tb22obtbta*xctbobtb)/(2*sqr(xa))-((-(tb22obtbta*xc)/sqr(xa)-4*xa)*yayaxctb22)/(6*ya))*(yayaxctb22sq+xctbobtb2tc);
//			Real diff_ya = 4*(yaxctb22ob/(12*sqr(ya))+1)*yayaxctb22*(yayaxctb22sq+xctbobtb2tc);
//			Real diff_xc = 2*(-((4*xc+tb22obtbta/xa)*yayaxctb22)/(6*ya)-2*xctbobtb)*(yayaxctb22sq+xctbobtb2tc);
//			Real diff_oa = 2*(-((((-2*ta-2*oa)*xc)/xa-2*ta-2*oa)*yayaxctb22)/(6*ya)-((-2*ta-2*oa)*xctbobtb)/(2*xa))*(yayaxctb22sq+xctbobtb2tc);
//			Real diff_ob = 2*(-((((2*tb+2*ob)*xc)/xa-2*tb-2*ob)*yayaxctb22)/(6*ya)-((2*tb+2*ob)*xctbobtb)/(2*xa))*(yayaxctb22sq+xctbobtb2tc);
//			Real diff_oc = 2*(2*(tc+oc)-((4*tc+4*oc)*yayaxctb22)/(6*ya))*(yayaxctb22sq+xctbobtb2tc);
//
//			value += error;
//			g[0] += diff_r;
//			g[1] += diff_xa;
//			g[2] += diff_ya;
//			g[3] += diff_xc;
//			g[4] += diff_oa;
//			g[5] += diff_ob;
//			g[6] += diff_oc;
//
////			Real error = ((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2)^2;
////			Real diff_r = -4*r*((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
////			Real diff_xa = 2*(((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))/(2*xa^2)-((-((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa^2-4*xa)*
////			(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya))*(
////			(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
////			Real diff_ya = 4*((6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya^2)+1)*
////			(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))*(
////			(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
////			Real diff_xc = 2*(-((4*xc+(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/xa)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya)
////			-2*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))*(
////			(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
////			Real diff_oa = 2*(-((((-2*ta-2*oa)*xc)/xa-2*ta-2*oa)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya)-
////			((-2*ta-2*oa)*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))/(2*xa))*((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
////			Real diff_ob = 2*(-((((2*tb+2*ob)*xc)/xa-2*tb-2*ob)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya)
////			-((2*tb+2*ob)*(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa)))/(2*xa))*((2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
////			Real diff_oc = 2*(2*(tc+oc)-((4*tc+4*oc)*(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya)))/(6*ya))*(
////			(2*ya-(6*ya^2+2*xc^2+((tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)*xc)/xa-2*xa^2+2*tc^2+4*oc*tc-tb^2-2*ob*tb-ta^2-2*oa*ta+2*oc^2-ob^2-oa^2)/(12*ya))^2+
////			(-xc-(tb^2+2*ob*tb-ta^2-2*oa*ta+ob^2-oa^2)/(4*xa))^2+(tc+oc)^2-r^2);
//		}
//
//		return value;
////		return sqr(sqr(r - 2) + sqr(xa - 3));
//	}
//
//	static int _progress(void *instance, const Real *x,
//			const Real *g, const Real fx,
//			const Real xnorm, const Real gnorm,
//			const Real step, int n, int k, int ls) {
//		return reinterpret_cast<objective_function*>(instance)->progress(x, g,
//				fx, xnorm, gnorm, step, n, k, ls);
//	}
//
//	int progress(const Real *x, const Real *g,
//			const Real fx, const Real xnorm,
//			const Real gnorm, const Real step, int n,
//			int k, int ls) {
//		printf("Iteration %d:\n", k);
//		printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
//		printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
//		printf("\n");
//		return 0;
//	}
//};
//
//
//int main(int argc, char **argv) {
//	objective_function obj;
//	return obj.run();
//}
