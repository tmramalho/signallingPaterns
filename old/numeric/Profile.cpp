/*
 * Profile.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 */

#include "Profile.h"

Profile::Profile(dvec *points) {
	/* Variable knots
	 * first S/2-k/2 points are xpositions for knots
	 * last S/2+k/2 points are ypositions for control points
	 */
	int k = 6;
	_dim = (int)points->size()/2 + k/2;
	_points = points;
	const size_t nbreak = _dim - k + 2;

	_x = gsl_vector_alloc(nbreak);
	_y = new double[_dim];

	(_x->data)[0] = 0;

	for(uint i = 1; i < nbreak-1; i++)
		(_x->data)[i] = points->at(i-1);

	//bubble sort
	bool check = 1;
	while(check) {
		check = 0;
		for(uint i = 1; i < nbreak-2; i++){
			if((_x->data)[i] > (_x->data)[i+1]) {
				double d = (_x->data)[i];
				(_x->data)[i] = (_x->data)[i+1];
				(_x->data)[i+1] = d;
				check = 1;
			}
		}
	}

	(_x->data)[nbreak - 1] = 1;

	for(int i = 0; i < _dim;  i++)
		_y[i] = points->at(i + nbreak - 2);

	bw = gsl_bspline_alloc(k, nbreak);
	B = gsl_vector_alloc(_dim);
	gsl_bspline_knots(_x, bw);
}

Profile::~Profile() {
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_vector_free(_x);
	delete[] _y;
}

double Profile::eval(double x) {
	double acc = 0;

	/* puts the value of each basis spline in B */
	gsl_bspline_eval(x, B, bw);

	/*computes the value of the spline */
	for(int i=0; i<_dim; i++)
		acc += _y[i] * (B->data)[i];

	return acc;
}
