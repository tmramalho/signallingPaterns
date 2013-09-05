/*
 *  LatPromReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "LatPromReaction.h"

/* Constructor: LatPromReaction(i_promoted_by_neighbors,i_promoting_neighbors,kinetic,K)
 * -------------------------------------------------------------------------- 
 * Constructs LatPromReaction, reading in data passed to it.
 */
LatPromReaction::LatPromReaction( int i_promoting_neighbors , int i_promoted_by_neighbors ,
								 double kinetic , double K ) {
	
	_type = LATERAL_PROMOTION;
	
	_num_part = 1;
	
	_i_promoted_by_neighbors = i_promoted_by_neighbors;
	_i_promoting_neighbors = i_promoting_neighbors;
	_kinetic = kinetic;
	_K = K;
	
}

LatPromReaction::LatPromReaction( std::ifstream& file ) {
	
	for (int i = 0; i < NUM_LINE; i++) {
		switch (i) {
			case I_PROMOTED_BY_NEIGHBORS_LINE:
				file.ignore(256,':');
				file >> _i_promoted_by_neighbors;
				break;
			case I_PROMOTING_NEIGHBORS_LINE:
				file.ignore(256,':');
				file >> _i_promoting_neighbors;
				break;
			case KINETIC_LINE:
				file.ignore(256,':');
				file >> _kinetic;
				break;
			case K_LINE:
				file.ignore(256,':');
				file >> _K;
				break;
			default:
				break;
		}
	}
	
}

Reaction* LatPromReaction::copy() {
	return new LatPromReaction(_i_promoting_neighbors,_i_promoted_by_neighbors,_kinetic,_K);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the protein in the current cell, which will be
 *		promoted by the precense of a protein in neighboring cells.
 *
 * Returns NEXIST for part_num != 0.
 */
int LatPromReaction::get_i_part( int part_num ) const {
	switch (part_num) {
		case 0:
			return _i_promoted_by_neighbors;
			break;
		default:
			return NEXIST;
			break;
	}
}

int LatPromReaction::get_i_dependent_molecule() const {
	return NEXIST;
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this LatPromReaction, occuring in i_curr_cell.
 *		
 *		(d/dt)[proteinZero]=
 *				kinetic*([proteinOne]^2)/(K+[proteinOne]^2)
 *
 * Where [proteinZero] is the concentration in the current cell and 
 * [proteinOne] is the average concentration in its neighbors.
 *
 */
void LatPromReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {
	
	/* Compute average concentration of neighboring protein */
	
	double neighbor_sum = 0.0;
	double avg_neighbor_conc;
	for (int i = 0; i < _sc_ref->_neighbors.at(i_curr_cell)->size() ; i++) {
		neighbor_sum += curr_tissue.at( _sc_ref->_neighbors.at(i_curr_cell)->at(i) , _i_promoting_neighbors );
	}
	avg_neighbor_conc = neighbor_sum/_sc_ref->_neighbors.at(i_curr_cell)->size();
	
	/* Compute dx_dt */
	double r = _kinetic * (pow(avg_neighbor_conc,5))/(_K+pow(avg_neighbor_conc,5));
	dx_dt.at( i_curr_cell , _i_promoted_by_neighbors ) += r;
	
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this LatPromReaction, occuring in i_curr_cell.
 *		
 *		(d/dt)[proteinZero]=
 *				kinetic*([proteinOne]^2)/(K+[proteinOne]^2)
 *
 * Where [proteinZero] is the concentration in the current cell and 
 * [proteinOne] is the average concentration in its neighbors.
 *
 * Stochastic term is
 *
 *		det_flow * random_var * sqrt(q/dt)
 *
 * where q measures the magnitued of the noise.
 *
 * We divde by sqrt(dt) so that dx_dt * dt is the correct term for the 
 * the numerical stochastic integration of time step dt (see Documentation,
 * section INTEGRATION).
 */
void LatPromReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
							boost::random::normal_distribution<>& dist , 
							boost::random::mt19937& generator ,  double q ) {
	
	/* Compute average concentration of neighboring protein */
	double neighbor_sum = 0.0;
	double avg_neighbor_conc;
	for (int i = 0; i < _sc_ref->_neighbors.at(i_curr_cell)->size() ; i++) {
		neighbor_sum += curr_tissue.at( _sc_ref->_neighbors.at(i_curr_cell)->at(i) , _i_promoting_neighbors );
	}
	avg_neighbor_conc = neighbor_sum/_sc_ref->_neighbors.at(i_curr_cell)->size();
	
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = 1 */
	
	double det_flow = _kinetic * (pow(avg_neighbor_conc,5))/(_K+pow(avg_neighbor_conc,5)); 
	
	double rand = dist(generator);
	double stoc_flow = det_flow * rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at( i_curr_cell , _i_promoted_by_neighbors ) += flow;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void LatPromReaction::update_mol_indices( int first_index , int num_insertion ) {
	_i_promoted_by_neighbors = update_index( first_index , num_insertion , _i_promoted_by_neighbors );
	_i_promoting_neighbors = update_index( first_index , num_insertion , _i_promoting_neighbors );
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates either the kinetic rate or K, randomly chosen with equal 
 * probability. Mutates each through multiplication by a random real chosen
 * uniformly between 0 and 2.
 */
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

/* Public Method: print_info(line_start)
 * -------------------------------------------------------------------------- 
 * Prints the data describing the reaction, with line_start beginning each 
 * line.
 *
 * Format:
 *
 *		<line_start> Reaction Type: Lateral Promotion Reaction
 *		<line_start> Index of Protein Promoted by Neighbors: <i>
 *		<line_start> Index of Protein Promoting its Neighbors: <i>
 *		<line_start> Kinetic Rate of Reaction: <r>
 *		<line_start> K: <r>
 *		
 */
void LatPromReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "Reaction Type: Lateral Promotion Reaction" << "\n";
	
	for (int i = 0; i < NUM_LINE; i++) {
		switch (i) {
			case I_PROMOTED_BY_NEIGHBORS_LINE:
				std::cout << line_start << "Index of Protein Promoted by Neighbors: " << _i_promoted_by_neighbors << "\n";
				break;
			case I_PROMOTING_NEIGHBORS_LINE:
				std::cout << line_start << "Index of Protein Promoting its Neighbors: " << _i_promoting_neighbors << "\n";
				break;
			case KINETIC_LINE:
				std::cout << line_start << "Kinetic Rate of Reaction: " << _kinetic << "\n";
				break;
			case K_LINE:
				std::cout << line_start << "K: " << _K << "\n";
				break;
			default:
				break;
		}
	}	
}

/* Public Method: to_file(file,line_start)
 * -------------------------------------------------------------------------- 
 */
void LatPromReaction::to_file ( std::ofstream& file , std::string line_start ) {
	
	file << line_start << "Reaction Type: Lateral Promotion Reaction" << "\n";
	
	for (int i = 0; i < NUM_LINE; i++) {
		switch (i) {
			case I_PROMOTED_BY_NEIGHBORS_LINE:
				file << line_start << "Index of Protein Promoted by Neighbors: " << _i_promoted_by_neighbors << "\n";
				break;
			case I_PROMOTING_NEIGHBORS_LINE:
				file << line_start << "Index of Protein Promoting its Neighbors: " << _i_promoting_neighbors << "\n";
				break;
			case KINETIC_LINE:
				file << line_start << "Kinetic Rate of Reaction: " << _kinetic << "\n";
				break;
			case K_LINE:
				file << line_start << "K: " << _K << "\n";
				break;
			default:
				break;
		}
	}
	
}

