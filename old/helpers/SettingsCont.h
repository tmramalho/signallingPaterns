/*
 * SettingsCont.h
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 *
 *  Singleton design pattern
 */

#ifndef SETTINGSCONT_H_
#define SETTINGSCONT_H_

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <string>

#include "../../IntegrationType.h"
#include "../../old/genetic/MutationType.h"
#include "../../Reactions/ReactionType.h"

#define SIG(x,b) ((1.0-b)/(1.0+exp(-x))+b)

class SettingsCont {

private:
	static SettingsCont *_sc_ref;
	SettingsCont();

public:
	
	/* Evolution Parameters */
	int _num_threads, _na, _ng, _mutations_per_gen;
	double _mutation_probability;
	std::vector<MutationType> _mutation_types;

	/* Integration Parameters */
	double _dt, _stochasticity;
	IntegrationType _mode;
	
	/* Tissue Geometry */
	std::vector< std::vector<int>* > _neighbors;
	
	/* User Interface */
	int _verbose;
	
	/* Data Tracking */
	int _do_track;
	
	/* Public Methods */
	void set_neighbors( std::string neighbor_code );
	static SettingsCont *getInstance();
	//void setParameters(int argc, char *argv[]);
	~SettingsCont();

};

#endif /* SETTINGSCONT_H_ */
