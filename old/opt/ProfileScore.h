/*
 * ProfileScore.h
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 */

#ifndef PROFILESCORE_H_
#define PROFILESCORE_H_
#include <ctime>
#include <cmath>

#include "cuba.h"
#include "cubature.h"
#include "FitnessFunction.h"
#include "../numeric/Profile.h"
#include "../helpers/SettingsCont.h"
#include "../helpers/ParamsCont.h"

class ProfileScore : public FitnessFunction {
public:
	ProfileScore();
	~ProfileScore();
	double run(dvec *x);

private:
	int verb;
	SettingsCont *s;

	static int integrand(const int *ndim, const double *x, const int *ncomp, double *fval, void *params);
};

#endif /* PROFILESCORE_H_ */
