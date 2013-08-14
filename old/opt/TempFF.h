/*
 *  TempFF.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 8/12/13.
 *
 */

# include "FitnessFunction.h"

#ifndef TEMPFF_H
#define TEMPFF_H

class TempFF : public FitnessFunction
{
public:
	TempFF() {}
	virtual ~TempFF() {}
	
	virtual double run(Manager *x);
	
};

#endif