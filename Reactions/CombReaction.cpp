/*
 *  CombReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "CombReaction.h"

/* Constructor: CombReaction(i_reac_zero,i_reac_one,i_product,
 *							 forward_kinetic,backward_kinetic)
 * -------------------------------------------------------------------------- 
 * Constructs CombReaction, reading in data passed to it.
 */
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

/* Public method: copy()
 * -------------------------------------------------------------------------- 
 */

Reaction *CombReaction::copy() {
	return new CombReaction(_i_reac_zero,_i_reac_one,_i_product,_forward_kinetic,_backward_kinetic);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the first uncomplexed reactant.
 *		Participant 1 is the second uncomplexed reactant.
 *		Participant 2 is the complex of participant 0 and 1.
 *
 * Returns NEXIST for part_num != 0,1, or 2.
 */
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

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this CombReaction, occuring in i_curr_cell.
 *
 *		reacZero + reacOne <--> product
 *		
 *		forwardReaction: reactant zero and reactant one combine to
 *		form product with rate forwardKinetic * [reacZero] * [reacOne].
 *		
 *		backwardReaction: product degrades into
 *		reactant zero and reactant one with rate backwardKinetic *
 *		[product].
 */
void CombReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {
	
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

/* Public method: react(curr_tissue,dx_dt,i_curr_cell,
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this CombReaction, occuring in i_curr_cell.
 *
 *		reacZero + reacOne <--> product
 *		
 *		forwardReaction: reactant zero and reactant one combine to
 *		form product with rate forwardKinetic * [reacZero] * [reacOne].
 *		
 *		backwardReaction: product degrades into
 *		reactant zero and reactant one with rate backwardKinetic *
 *		[product].
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

void CombReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
						 boost::random::normal_distribution<>& dist , 
						 boost::random::mt19937& generator , double q ) {
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = det_flow */
	
	double det_flow = 
	_forward_kinetic 
	* curr_tissue.at(i_curr_cell,_i_reac_zero)
	* curr_tissue.at(i_curr_cell,_i_reac_one) 
	-
	_backward_kinetic 
	* curr_tissue.at(i_curr_cell,_i_product);
	
	double rand = dist(generator);
	double stoc_flow = rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at(i_curr_cell,_i_reac_zero) -= flow;
	dx_dt.at(i_curr_cell,_i_reac_one) -= flow;
	dx_dt.at(i_curr_cell,_i_product) += flow; 
	
}



/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void CombReaction::update_indices( int first_index , int num_insertion ) {
	/* Note: Because NEXIST = -1, non-existant participants will never be 
	 * updated by the insertion procedure as desired.
	 */
	if (_i_reac_zero >= first_index) _i_reac_zero += num_insertion;
	if (_i_reac_one >= first_index) _i_reac_one += num_insertion;
	if (_i_product >= first_index) _i_product += num_insertion;
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates either forward kinetic or backward kinetic rate, chosen randomly.
 * Mutates them through multiplication by random real chosen uniformly
 * between 0 and 2.
 */
void CombReaction::mutate( boost::random::mt19937& generator ) {
	
	boost::random::uniform_int_distribution<> which_kinetic(0,1);
	boost::random::uniform_real_distribution<> mutation_factor(0.0,2.0);
	
	int kinetic_num = which_kinetic(generator);
	double mut_factor = mutation_factor(generator);
	
	switch ( kinetic_num ) {
		case 0:
			_forward_kinetic *= mut_factor;
			break;
		case 1:
			_backward_kinetic *= mut_factor;
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
 *		<line_start> Reaction Type: Combination Reaction
 *		<line_start> Index of First Constituent Protein (in Complex): <i>
 *		<line_start> Index of Second Constituent Protein (in Complex): <i>
 *		<line_start> Index of Produced Complex: <i>
 *		<line_start> Forward Kinetic (binding): <r>
 *		<line_start> Backward Kinetic (separation): <r>
 *		
 */
void CombReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "ReactionType: Combination Reaction" << std::endl;
	std::cout << line_start << "Index of First Constituent Protein: " << _i_reac_zero << std::endl;
	std::cout << line_start << "Index of Second Constituent Protein: " << _i_reac_one << std::endl;
	std::cout << line_start << "Index of Produced Complex: " << _i_product << std::endl;
	std::cout << line_start << "Foward Kinetic (binding): " << _forward_kinetic << std::endl;
	std::cout << line_start << "Backward Kinetic (separation): " << _backward_kinetic << std::endl;

}

/* Public Method: print_info(file,line_start)
 * -------------------------------------------------------------------------- 
 */
void CombReaction::to_file ( std::ofstream& file , std::string line_start ) {
	
	file << line_start << "ReactionType: Combination Reaction\n";
	file << line_start << "Index of First Constituent Protein:\n";
	file << line_start << "Index of Second Constituent Protein:\n";
	file << line_start << "Index of Produced Complex: " << _i_product << "\n";
	file << line_start << "Foward Kinetic (binding): " << _forward_kinetic << "\n";
	file << line_start << "Backward Kinetic (separation): " << _backward_kinetic << "\n";
	
}


