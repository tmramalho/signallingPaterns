/*
 * Profile.h
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 */

#ifndef PROFILE_H_
#define PROFILE_H_
#include <gsl/gsl_bspline.h>
#include "dvec.h"

class Profile {
public:
	Profile(dvec *points);
	~Profile();

	double eval(double x);

private:
	dvec *_points;
	int _dim;
	gsl_vector *_x;
	double *_y;
	gsl_bspline_workspace *bw;
	gsl_vector *B;
};

#endif /* PROFILE_H_ */
