/*
 * BayesScore.h
 *
 *  Created on: Sep 26, 2012
 *      Author: tiago
 */

#ifndef BAYESSCORE_H_
#define BAYESSCORE_H_

#include <ctime>
#include <cmath>

#include "cuba.h"
#include "FitnessFunction.h"
#include "../numeric/Profile.h"
#include "../numeric/Interpolate.h"
#include "../helpers/SettingsCont.h"
#include "../helpers/ParamsCont.h"

class BayesScore: public FitnessFunction {
public:
	BayesScore();
	~BayesScore();
	double run(dvec *x);

private:
	int verb;
	SettingsCont *s;

	Interpolate *createPositionEstimates(paramCont *params);
	static int integrand(const int *ndim, const double *x, const int *ncomp, double *fval, void *params);
};

#endif /* BAYESSCORE_H_ */
