#include "TestFunction.h"

TestFunction::TestFunction(int dim, int type)
{
	_dim = dim;
	_type = type;
}

TestFunction::~TestFunction()
{
}

double TestFunction::run(dvec *x) {
	switch (_type) {
		case 0:
			if(_dim > 1) {
				return (1 - x->at(0))*(1 - x->at(0)) + 100 * 
					(x->at(1) - x->at(0)*x->at(0))*(x->at(1) - x->at(0)*x->at(0));
			} else return 0;
		case 1:
			if(_dim > 1) {
				return (x->at(0)*x->at(0) + x->at(1) - 11.0)*(x->at(0)*x->at(0) + x->at(1) - 11.0) +
						(x->at(0) + x->at(1)*x->at(1) - 7.0)*(x->at(0) + x->at(1)*x->at(1) - 7.0);
			} else return 0;
		case 2:
			if(_dim > 1) {
				double acc = 0;
				for(int j = 0; j<_dim; j++) {
					acc += x->at(j) * x->at(j) - 10*cos(2*PI*x->at(j));
				}
				return 10 * _dim + acc;
			} else return 0;
		default:
			return 0;
	}
	return 0;
}
