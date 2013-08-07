/*
 *  LatPromReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "LatPromReaction.h"

LatPromReaction::LatPromReaction() {
	
	_type = LATERAL_PROMOTION;
	
	_num_part = 1;
	
}

LatPromReaction::LatPromReaction( int i_local_prot , int i_neighbor_prot ,
								 double kinetic ,
								 double K ) {
	
	_type = LATERAL_PROMOTION;
	
	_num_part = 1;
	
	_i_local_prot = i_local_prot;
	_i_neighbor_prot = i_neighbor_prot;
	_kinetic = kinetic;
	_K = K;
	
}

LatPromReaction::LatPromReaction(LatPromReaction* newOne) {
	
	_type = LATERAL_PROMOTION;
	
	_num_part = 1;
	
	_i_local_prot = newOne->_i_local_prot;
	_i_neighbor_prot = newOne->_i_neighbor_prot;
	_kinetic = newOne->_kinetic;
	_K = newOne->_K;
	
}

LatPromReaction::~LatPromReaction() {}

int LatPromReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_local_prot;
			break;
		default:
			return NEXIST;
			break;
	}
}

void LatPromReaction::react( dmat& curr_tissue , dmat& dx_dt , 
							std::vector< std::vector<int>* >& neighbors,
							int i_curr_cell ) {
	
	/* Compute average concentration of neighboring protein */
	double neighbor_sum = 0.0;
	double avg_neighbor_conc;
	for (int i = 0; i < neighbors.at(i_curr_cell)->size() ; i++) {
		neighbor_sum += curr_tissue.at( neighbors.at(i_curr_cell)->at(i) , _i_neighbor_prot );
	}
	avg_neighbor_conc = neighbor_sum/neighbors.at(i_curr_cell)->size();
	
	/* Compute dx_dt */
	double r = _kinetic * (pow(avg_neighbor_conc,2))/(_K+pow(avg_neighbor_conc,2));
	dx_dt.at( i_curr_cell , _i_local_prot ) += r;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of _i_local_prot, _i_neighbor_prot, given an insertion 
 * of size num_insertion, beginning at first_index, into our dvecs containing 
 * molecule concentrations in our manager.
 *
 * See description for update_indices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void LatPromReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_local_prot >= first_index) _i_local_prot += num_insertion;
	if (_i_neighbor_prot >= first_index) _i_neighbor_prot += num_insertion;
}

void LatPromReaction::mutate ( boost::random::mt19937& generator ) {
	
	boost::random::uniform_int_distribution<> which_kinetic(0,1);
	boost::random::uniform_real_distribution<> mutation_factor(0.0,2.0);
	
	int kinetic_num = which_kinetic(generator);
	double mut_factor = mutation_factor(generator);
	
	switch ( kinetic_num ) {
		case 0:
			_kinetic *= mut_factor;
			break;
		case 1:
			_K *= mut_factor;
			break;
		default:
			break;
	}
}

void LatPromReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "ReactionType: Lateral Promotion Reaction" << std::endl;
	//std::cout << line_start << "_num_part: " << _num_part << std::endl;
	std::cout << line_start << "Index of Protein Promoted by Neighbors: " << _i_local_prot << std::endl;
	std::cout << line_start << "Index of Protein Promoting its Neighbors: " << _i_neighbor_prot << std::endl;
	std::cout << line_start << "Kinetic Rate of Reaction: " << _kinetic << std::endl;
	std::cout << line_start << "K: " << _K << std::endl;
	
}

