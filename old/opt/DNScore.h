/*
 *  DNScore.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 8/20/13.
 *
 */

/* Fitness function for creating an alternating pattern of proteins between cells.
 */

#ifndef DNSCORE_H
#define DNSCORE_H

#include "FitnessFunction.h"

class DNScore : public FitnessFunction {
public:
	DNScore() {}
	~DNScore() {}
	
	double run(Manager *x,boost::random::mt19937& generator);
};

#endif