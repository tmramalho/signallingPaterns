#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <gsl/gsl_spline.h>

class Interpolate
{
public:
	Interpolate(double *x, double *y, int dim);
	virtual ~Interpolate();
	double value(double x);
	double deriv(double x);
	double deriv2(double x);
	double integ(double a, double b);
	
private:
	gsl_spline *_spline;
	gsl_interp_accel *_acc;
	
};

#endif /*INTERPOLATE_H_*/
