/*
 *  HillPromReaction.cpp
 *	Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#include "HillPromReaction.h"

HillPromReaction::HillPromReaction() {
	
	_type = HILL_PROMOTION;
	
	_num_part = 1;

}

HillPromReaction::HillPromReaction( int i_promoter , int i_promoted ,
								   double kinetic , double K , double cooperativity ) {
	
	_type = HILL_PROMOTION;
	
	_num_part = 1;
	
	_i_promoter = i_promoter;
	_i_promoted = i_promoted;
	_kinetic = kinetic;
	_K = K;
	_cooperativity = cooperativity;
	
}

HillPromReaction::~HillPromReaction() {}

int HillPromReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_promoted;
			break;
		default:
			return NEXIST;
			break;
	}
}

void HillPromReaction::react( dmat& curr_tissue , dmat& dx_dt ,
							 std::vector< std::vector<int>* >& neighbors ,
							 int i_curr_cell ) {

	double r = 
	_kinetic
	* pow( curr_tissue.at(i_curr_cell , i_promoter ) , _cooperativity )
	/
	( _K + pow( curr_tissue.at(i_curr_cell , i_promoter ) , _cooperativity ));
	 
	dx_dt.at(i_curr_cell,_i_promoted) += r;
}

void HillPromReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_promoter >= first_index) _i_promoter += num_insertion;
	if (_i_promoted >= first_index) _i_promoted += num_insertion;
}


