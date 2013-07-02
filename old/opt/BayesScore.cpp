/*
 * BayesScore.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: tiago
 */

#include "BayesScore.h"

BayesScore::BayesScore() {
	SettingsCont *s = SettingsCont::getInstance();
	verb = s->verbose;
}

BayesScore::~BayesScore() {
}

double BayesScore::run(dvec *x) {
	clock_t init = clock();

	s = SettingsCont::getInstance();

	for(uint i=0; i < x->size(); i++) {
		if(x->at(i) > 1) return 10;
		if(x->at(i) < s->bound) return 10;
	}

	Profile *m = new Profile(x);

	paramCont *params = new paramCont;
	params->g = s->bound;
	params->n = s->coop;
	params->p = s->p;
	params->s0 = s->sigma0;
	params->n0 = s->poiss0;
	params->m = m;

	Interpolate *barx = this->createPositionEstimates(params);
	params->xe = barx;

	int nregions, neval, fail;
	double res, err, prob;

	/*CUHRE ROUTINE*/
	Cuhre(2, 1, &(this->integrand), (void *)params, 1e-6, 1e-6,
		  0, 1e4, 1e7, 13, &nregions, &neval, &fail, &res, &err, &prob);
	if(verb == 2) std::cout << res << " +- " << err << " for " << neval << std::endl;

	delete m;
	delete params;
	delete barx;
	double time = (double)(clock()-init) / ((double)CLOCKS_PER_SEC);
	if(verb == 2) std::cout << "Got score " << res << " in time: " << time << std::endl;

	return res;
}
/* return the bayesian estimate for the position given some data
 * point as a function of the data
 * x_estimate(d) = \int dx P(d|x)P(x)x/(\int dx' P(d|x')P(x'))
 * where P(x) = 1
 */
Interpolate* BayesScore::createPositionEstimates(paramCont *p) {
	double x[100];
	double y[100];
	for(uint i = 0; i < 100; i++) { //iterate each data value
		double dp = (double) i / 99;
		double xd = 0;
		double pd = 0;
		for(uint j = 0; j < 100; j++) { //approximate the integral
			double xp = (double) j / 99;
			double mx = p->m->eval(xp);
			double sm = p->s0*(p->p*mx*mx+mx);
			double pdgx = exp(-(dp - mx)*(dp - mx)/(2*sm))/sqrt(2*PI*sm);
			xd += pdgx*xp;
			pd += pdgx;
		} //xd and pd should be divided by 100 but it cancels out
		x[i] = dp;
		y[i] = xd / pd;
	}

	Interpolate *xestimate = new Interpolate(x, y, 100);
	return xestimate;
}

int BayesScore::integrand(const int *ndim, const double *x, const int *ncomp, double *fval, void *params) {
	paramCont *p = (paramCont *)params;

	if(*ndim == 2) {
		double xp = x[0];
		double dp = x[1];

		double mx = p->m->eval(xp);

		double sm = p->s0*(p->p*mx*mx+mx);
		double pdgx = exp(-(dp - mx)*(dp - mx)/(2*sm))/sqrt(2*PI*sm);
		double xdif = xp - p->xe->value(dp);
		(*fval) = pdgx*xdif*xdif;

		return 0;
	}

	else return -999;
}
