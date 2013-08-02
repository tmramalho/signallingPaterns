/*
 *  PromReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "PromReaction.h"

PromReaction::PromReaction() {
	
	_type = PROMOTION;
	
	_num_part = 1;
	
}

PromReaction::PromReaction( int i_gene , int i_prot , double kinetic ) {
	
	_type = PROMOTION;
	
	_num_part = 1;
	
	_i_gene = i_gene;
	_i_prot = i_prot;
	_kinetic = kinetic;
	
}

PromReaction::~PromReaction() {}

int PromReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_prot;
			break;
		default:
			return NEXIST;
			break;
	}
}

void PromReaction::react( dmat& curr_tissue , dmat& dx_dt ,
						  std::vector< std::vector<int>* >& neighbors,
						  int i_curr_cell ) {
	
	double r = _kinetic * curr_tissue.at( i_curr_cell , _i_gene );
	dx_dt.at( i_curr_cell , _i_prot ) += r;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of _i_gene, _i_prot, given an insertion of size 
 * num_insertion, beginning at first_index, into our dvecs containing 
 * molecule concentrations in our manager.
 *
 * See description for updateIndices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void PromReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_gene >= first_index) _i_gene += num_insertion;
	if (_i_prot >= first_index) _i_prot += num_insertion;
}
