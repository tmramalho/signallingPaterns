/*
 *  HillDegReaction.cpp
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#include "HillRepReaction.h"

HillRepReaction::HillRepReaction() {
	
	_type = HILL_REPRESSION;
	
	_num_part = 1;
	
}

HillRepReaction::HillRepReaction( int i_repressor , int i_repressed ,
								   double kinetic , double K , double cooperativity ) {
	
	_type = HILL_REPRESSION;
	
	_num_part = 1;
	
	_i_repressor = i_repressor;
	_i_repressed = i_repressed;
	_kinetic = kinetic;
	_K = K;
	_cooperativity = cooperativity;
	
}

HillRepReaction::~HillRepReaction() {}

int HillRepReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_repressed;
			break;
		default:
			return NEXIST;
			break;
	}
}

void HillRepReaction::react( dmat& curr_tissue , dmat& dx_dt ,
							 std::vector< std::vector<int>* >& neighbors ,
							 int i_curr_cell ) {
	
	double r = 
	_kinetic
	* _K
	/
	( _K + pow( curr_tissue.at(i_curr_cell , i_repressor ) , _cooperativity ));
	
	dx_dt.at(i_curr_cell,_i_repressed) += r;
}

void HillRepReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_repressor >= first_index) _i_repressor += num_insertion;
	if (_i_repressed >= first_index) _i_repressed += num_insertion;
}
