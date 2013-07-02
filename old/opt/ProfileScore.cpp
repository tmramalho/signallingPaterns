/*
 * ProfileScore.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 */

#include "ProfileScore.h"

ProfileScore::ProfileScore() {
	SettingsCont *s = SettingsCont::getInstance();
	verb = s->verbose;
}

ProfileScore::~ProfileScore() {
}

double ProfileScore::run(dvec *x) {
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

	int nregions, neval, fail;
	double res, err, prob;

	/*VEGAS ROUTINE
	Vegas(2, 1, &(this->integrand), (void*)params, 1e-5, 1e-4, 0, time(NULL), 1e5, 1e6,
			1000, 100, 1000, 0, NULL, &neval, &fail, &res, &err, &prob);
	if(verb) std::cout << res << " +- " << err << " for " << neval << std::endl;*/

	/*SUAVE ROUTINE
	Suave(2, 1, &(this->integrand), (void*)params, 1e-5, 1e-4, 0, time(NULL), 1e5, 1e6,
			10000, 1, &nregions, &neval, &fail, &res, &err, &prob);
	if(verb) std::cout << res << " +- " << err << " for " << neval << std::endl;*/

	/*DIVONNE ROUTINE
	int ngiven = 1000, ldxgiven = 2;
	double *xgiven = new double[1000];
	for(int i = 0; i< ngiven/2; i++)
		xgiven[2*i] = xgiven[2*i+1] = 2.0*i/100.0;

	Divonne(2, 1, &(this->integrand), (void*)params, 1e-5, 1e-4, 0, time(NULL), 1e4, 1e7,
			30, 30, 1, 4, 0, 10000, 1e-3, ngiven, ldxgiven, xgiven, 0, 0, &nregions, &neval, &fail, &res, &err, &prob);
	if(verb) std::cout << res << " +- " << err << " for " << neval << std::endl;*/

	/*CUHRE ROUTINE*/
	Cuhre(2, 1, &(this->integrand), (void *)params, 1e-6, 1e-6,
		  0, 1e4, 1e7, 13, &nregions, &neval, &fail, &res, &err, &prob);
	if(verb == 2) std::cout << res << " +- " << err << " for " << neval << std::endl;

	delete m;
	delete params;
	double time = (double)(clock()-init) / ((double)CLOCKS_PER_SEC);
	if(verb == 2) std::cout << "Got score " << res << " in time: " << time << std::endl;

	return res;
}

int ProfileScore::integrand(const int *ndim, const double *x, const int *ncomp, double *fval, void *params) {
	paramCont *p = (paramCont *)params;

	if(*ndim == 2) {
		double xp = x[0];
		double yp = x[1]*0.65 + 0.15;//y from 0.15 to 0.8

		double mx = p->m->eval(xp);
		double my = p->m->eval(yp);

		double th = (xp < yp) ? 1 : 0;
		double h = pow(mx, p->n)/(pow(mx, p->n) + pow(my, p->n));
		double dh = pow(mx, p->n-1.0)*pow(my, p->n)*p->n/
				((pow(mx, p->n) + pow(my, p->n))*(pow(mx, p->n) + pow(my, p->n)));

		double dif = h - th;

		double sm = p->s0*(p->p*mx*mx+mx);
		double poisson = p->n0*h;
		double integrand = dif*dif + poisson + dh*dh*sm;

		(*fval) = integrand*0.65;

		return 0;
	}

	else return -999;
}
