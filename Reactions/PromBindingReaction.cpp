/*
 *  PromBindingReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "PromBindingReaction.h"

PromBindingReaction::PromBindingReaction() {
	
	_type = PROMOTER_BINDING;
	
	_num_part = 2;
	
}

PromBindingReaction::PromBindingReaction( int i_root_gene , int i_promoted_gene , int i_bound_protein ,
										 double forward_kinetic , double backward_kinetic ) {
	
	_type = PROMOTER_BINDING;
	
	_num_part = 2;
	
	_i_root_gene = i_root_gene;
	_i_promoted_gene = i_promoted_gene;
	_i_bound_protein = i_bound_protein;
	_forward_kinetic = forward_kinetic;
	_backward_kinetic = backward_kinetic;
	
}

PromBindingReaction::~PromBindingReaction() {}

int PromBindingReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_root_gene;
			break;
		case 1:
			return _i_promoted_gene;
		default:
			return NEXIST;
			break;
	}
}

void PromBindingReaction::react( dmat& curr_tissue , dmat& dx_dt ,
								std::vector< std::vector<int>* >& neighbors,
								int i_curr_cell ) {
	
	double r =
	_forward_kinetic
	* curr_tissue.at(i_curr_cell , _i_root_gene)
	* curr_tissue.at(i_curr_cell, _i_bound_protein)
	- 
	_backward_kinetic
	* curr_tissue.at(i_curr_cell,_i_promoted_gene);
	
	dx_dt.at(i_curr_cell,_i_root_gene) -= r;
	dx_dt.at(i_curr_cell, _i_promoted_gene) += r;
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of _i_root_gene, _i_promoted_gene, _i_bound_protein, 
 * given an insertion of size num_insertion, beginning at first_index, into 
 * our dvecs containing molecule concentrations in our manager.
 *
 * See description for updateIndices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUM_INSERTION >= 0. 
 */

void PromBindingReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_root_gene >= first_index) _i_root_gene += num_insertion;
	if (_i_promoted_gene >= first_index) _i_promoted_gene += num_insertion;
	if (_i_bound_protein >= first_index) _i_bound_protein += num_insertion;
}



