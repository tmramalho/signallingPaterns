/*
 *  DegReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "DegReaction.h"

DegReaction::DegReaction() {
	
	_type = DEGRADATION;
	
	_num_part = 1;
	
}

DegReaction::DegReaction( int i_reac , double kinetic ) {
	
	_type = DEGRADATION;
	
	_num_part = 1;
	
	_i_reac = i_reac;
	_kinetic = kinetic;
	
}

DegReaction::DegReaction(DegReaction* newOne) {
	
	_type = DEGRADATION;
	
	_num_part = 1;
	
	_i_reac = newOne->_i_reac;
	_kinetic = newOne->_kinetic;
	
}

DegReaction::~DegReaction() {}

int DegReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_reac;
			break;
		default:
			return NEXIST;
			break;
	}
}

void DegReaction::react( dmat& curr_tissue , dmat& dx_dt ,
						 std::vector< std::vector<int>* >& neighbors,
						 int i_curr_cell ) {
	
	double r = _kinetic * curr_tissue.at(i_curr_cell,_i_reac);
	
	dx_dt.at(i_curr_cell,_i_reac) -= r;
	
}

/* Public Method: updateIndices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of _i_reac, given an insertion of size num_insertion, 
 * beginning at first_index, into our dvecs containing molecule 
 * concentrations in our manager.
 *
 * See description for updateIndices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void DegReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_reac >= first_index) _i_reac += num_insertion;
}

void DegReaction::mutate ( boost::random::mt19937& generator ) {
	
	boost::random::uniform_real_distribution<> mutation_factor(0.0,2.0);
	
	double mut_factor = mutation_factor(generator);
	
	_kinetic *= mut_factor;
	
}

void DegReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "ReactionType: Degredation Reaction" << std::endl;
	//std::cout << line_start << "_num_part: " << _num_part << std::endl;
	std::cout << line_start << "Index of Degrading Protein: " << _i_reac << std::endl;
	std::cout << line_start << "Degredation Rate: " << _kinetic << std::endl;
	
}

