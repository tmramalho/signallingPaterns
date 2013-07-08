#include "Interpolate.h"

Interpolate::Interpolate(double *x, double *y, int dim)
{
	_acc = gsl_interp_accel_alloc ();
	_spline = gsl_spline_alloc (gsl_interp_akima, dim);
	gsl_spline_init (_spline, x, y, dim);
}

double Interpolate::value(double x) {
	return gsl_spline_eval (_spline, x, _acc);
}

double Interpolate::deriv(double x) {
	return gsl_spline_eval_deriv (_spline, x, _acc);
}

double Interpolate::deriv2(double x) {
	return gsl_spline_eval_deriv2 (_spline, x, _acc);
}

double Interpolate::integ(double a, double b) {
	return gsl_spline_eval_integ (_spline, a, b, _acc);
}

Interpolate::~Interpolate()
{
	gsl_spline_free (_spline);
	gsl_interp_accel_free (_acc);
}
