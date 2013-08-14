/*
 *  PromBindingReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "PromBindingReaction.h"

/* Constructor: PromBindingReaction(i_root_gene,i_promoted_gene,i_bound_protein
 *									forward_kinetic,backward_kinetic)
 * -------------------------------------------------------------------------- 
 * Constructs PromBindingReaction, reading in data passed to it.
 */
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

Reaction* PromBindingReaction::copy() {
	return new PromBindingReaction(_i_root_gene,_i_promoted_gene,_i_bound_protein,_forward_kinetic,_backward_kinetic);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the un-promoted gene to which the protein will bind.
 *		Participant 1 is the promoted gene, complexed with the protein.
 *
 * Returns NEXIST for part_num != 0 or 1.
 */
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

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this PromBindingReaction, occuring in i_curr_cell.
 *		
 *		a <--> a:B
 *
 *		(d/dt)[a:B] = forwardKinetic * [a][B] - backwardKinetic * [a:B]
 *		(d/dt)[a] = backwardKinetic * [a:B] - forwardKinetic * [a][B]
 *
 */
void PromBindingReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {
	
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

/* Public method: react(curr_tissue,dx_dt,i_curr_cell
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this PromBindingReaction, occuring in i_curr_cell.
 *		
 *		a <--> a:B
 *
 *		(d/dt)[a:B] = forwardKinetic * [a][B] - backwardKinetic * [a:B]
 *		(d/dt)[a] = backwardKinetic * [a:B] - forwardKinetic * [a][B]
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
void PromBindingReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
						 boost::random::normal_distribution<>& dist , 
						 boost::random::mt19937& generator , double q ) {
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = 1 */
	
	double det_flow = 
	_forward_kinetic
	* curr_tissue.at(i_curr_cell , _i_root_gene)
	* curr_tissue.at(i_curr_cell, _i_bound_protein)
	- 
	_backward_kinetic
	* curr_tissue.at(i_curr_cell,_i_promoted_gene); 
	
	double rand = dist(generator);
	double stoc_flow = rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at(i_curr_cell,_i_root_gene) -= flow;
	dx_dt.at(i_curr_cell, _i_promoted_gene) += flow;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */
void PromBindingReaction::update_indices( int first_index , int num_insertion ) {
	/* Note: Because NEXIST = -1, non-existant participants will never be 
	 * updated by the insertion procedure as desired.
	 */		
	if (_i_root_gene >= first_index) _i_root_gene += num_insertion;
	if (_i_promoted_gene >= first_index) _i_promoted_gene += num_insertion;
	if (_i_bound_protein >= first_index) _i_bound_protein += num_insertion;
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates either forward kinetic or backward kinetic rate, chosen randomly.
 * Mutates them through multiplication by random real chosen uniformly
 * between 0 and 2.
 */
void PromBindingReaction::mutate ( boost::random::mt19937& generator ) {
	
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
 *		<line_start> Reaction Type: Promoter Binding Reaction
 *		<line_start> Index of Gene Without Bound Promoter: <i>
 *		<line_start> Index of Protein Binding to Gene: <i>
 *		<line_start> Index of Gene-Protein Complex: <i>
 *		<line_start> Forward Kinetic: <r>
 *		<line_start> Backward Kinetic: <r>
 *		
 */
void PromBindingReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "Reaction Type: Promoter Binding Reaction" << std::endl;
	std::cout << line_start << "Index of Gene Without Bound Promoter: " << _i_root_gene << std::endl;
	std::cout << line_start << "Index of Protein Binding to Gene: " << _i_bound_protein << std::endl;
	std::cout << line_start << "Index of Gene-Protein Complex: " << _i_promoted_gene << std::endl;
	std::cout << line_start << "Forward Kinetic (binding): " << _forward_kinetic << std::endl;
	std::cout << line_start << "Backward Kinetic (degredation): " << _backward_kinetic << std::endl;
	
}

/* Public Method: to_file(file,line_start)
 * -------------------------------------------------------------------------- 
 */
void PromBindingReaction::to_file ( std::ofstream& file ,  std::string line_start ) {
	
	file << line_start << "Reaction Type: Promoter Binding Reaction\n";
	file << line_start << "Index of Gene Without Bound Promoter: " << _i_root_gene << "\n";
	file << line_start << "Index of Protein Binding to Gene: " << _i_bound_protein << "\n";
	file << line_start << "Index of Gene-Protein Complex: " << _i_promoted_gene << "\n";
	file << line_start << "Forward Kinetic (binding): " << _forward_kinetic << "\n";
	file << line_start << "Backward Kinetic (degredation): " << _backward_kinetic << "\n"
	;
	
}



