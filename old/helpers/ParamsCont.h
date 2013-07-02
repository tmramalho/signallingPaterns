/*
 * ParamsCont.h
 *
 *  Created on: Sep 26, 2012
 *      Author: tiago
 */

#ifndef PARAMSCONT_H_
#define PARAMSCONT_H_
#include "../numeric/Profile.h"
#include "../numeric/Interpolate.h"

typedef struct paramCont {
	double g;
	double n;
	double p;
	double s0;
	double n0;
	Profile *m;
	Interpolate *xe;
} paramCont;

#endif /* PARAMSCONT_H_ */
