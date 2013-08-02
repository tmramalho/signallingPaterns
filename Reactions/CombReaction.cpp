/*
 *  CombReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "CombReaction.h"

CombReaction::CombReaction() {
	
	_type = COMBINATION;
	
	_num_part = 3;
	
}

CombReaction::CombReaction( int i_reac_zero , int i_reac_one , int i_product ,
						   double forward_kinetic , double backward_kinetic ) {
	
	_type = COMBINATION;
	
	_num_part = 3;
	
	_i_reac_zero = i_reac_zero;
	_i_reac_one = i_reac_one;
	_i_product = i_product;
	_forward_kinetic = forward_kinetic;
	_backward_kinetic = backward_kinetic;
	
}

CombReaction::~CombReaction() {}

int CombReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_reac_zero;
			break;
		case 1:
			return _i_reac_one;
			break;
		case 2:
			return _i_product;
			break;
		default:
			return NEXIST;
			break;
	}
}

void CombReaction::react( dmat& curr_tissue , dmat& dx_dt ,
						  std::vector< std::vector<int>* >& neighbors ,
						  int i_curr_cell ) {
	
	
	double flow = 
	_forward_kinetic 
	* curr_tissue.at(i_curr_cell,_i_reac_zero)
	* curr_tissue.at(i_curr_cell,_i_reac_one) 
	-
	_backward_kinetic 
	* curr_tissue.at(i_curr_cell,_i_product);
	
	dx_dt.at(i_curr_cell,_i_reac_zero) -= flow;
	dx_dt.at(i_curr_cell,_i_reac_one) -= flow;
	dx_dt.at(i_curr_cell,_i_product) += flow;
		
}

/* Public Method: updateIndices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of iReacZero, iReacOne, iProduct, given an insertion
 * of size numInsertions, beginning at firstIndex, into our dvecs containing
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

void CombReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_reac_zero >= first_index) _i_reac_zero += num_insertion;
	if (_i_reac_one >= first_index) _i_reac_one += num_insertion;
	if (_i_product >= first_index) _i_product += num_insertion;
}


