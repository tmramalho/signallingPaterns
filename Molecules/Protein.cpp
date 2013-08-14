/*
 *  Protein.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Protein.h"

/* Constructor: Protein(i_self,i_root_zero,i_root_one,init_conc)
 * -------------------------------------------------------------------------- 
 */
Protein::Protein( int i_self , int i_root_zero , int i_root_one , double init_conc ) {
	_i_root_zero = i_root_zero;
	_i_root_one = i_root_one;
	_i_self = i_self;
	_init_conc = init_conc;
	/* _reactions will automatically initiate to empty vector */
}

/* Public Method: get_i_root_zero()
 * -------------------------------------------------------------------------- 
 * If this protein is a complex, it returns the index in the dvecs of the 
 * first member protein in this complex. If it is not a complex, returns
 * NEXIST.
 */
int Protein::get_i_root_zero() {
	return _i_root_zero;
}

/* Public Method: get_i_root_one()
 * -------------------------------------------------------------------------- 
 * If this protein is a complex, it returns the index in the dvecs of the 
 * second member protein in this complex. If it is not a complex, returns
 * NEXIST.
 */
int Protein::get_i_root_one() {
	return _i_root_one;
}

/* Private Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of i_self, i_product, i_bound_promoter, i_root, given an 
 * insertion of size num_insertion of molecules beginning at first_index.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void Protein::update_indices( int first_index , int num_insertion ) {
	if (_i_self >= first_index) _i_self += num_insertion;
	if (_i_root_zero >= first_index) _i_root_zero += num_insertion;
	if (_i_root_one >= first_index) _i_root_one += num_insertion;
}

/* Public Method: print_info(line_start)
 * -------------------------------------------------------------------------- 
 * Prints the data describing the protein, with line_start beginning each 
 * line.
 *
 * Format:
 *
 *		<line_start> Index of Self: <i>
 *		<line_start> Index of First Constituent Protein (in Complex): <i>
 *		<line_start> Index of Second Constituent Protein (in Complex): <i>
 *		
 */

void Protein::print_info( std::string line_start ) {

	std::cout << line_start << "Index of Self: " << _i_self << std::endl;
	
	if ( _i_root_zero != NEXIST ) {
		std::cout << line_start << "Index of First Constituent Protein (in Complex): ";
		std::cout << _i_root_zero << std::endl;
	}
	else {
		std::cout << line_start << "Index of First Constituent Protein (in Complex): "; 
		std::cout << "NOT A COMPLEX" << std::endl;
	}
	
	if ( _i_root_one != NEXIST ) {
		std::cout << line_start << "Index of Second Constituent Protein (in Complex): ";
		std::cout << _i_root_one << std::endl;
	}
	else {
		std::cout << line_start << "Index of Second Constituent Protein (in Complex): "; 
		std::cout << "NOT A COMPLEX" << std::endl;
	}
	
}

/* Public Method: to_file(file,line_start)
 * -------------------------------------------------------------------------- 
 */

void Protein::to_file( std::ofstream& file , std::string line_start ) {
	
	file << line_start << "Index of Self: " << _i_self << "\n";
	
	if ( _i_root_zero != NEXIST ) {
		file << line_start << "Index of First Constituent Protein (in Complex): ";
		file << _i_root_zero << "\n";
	}
	else {
		file << line_start << "Index of First Constituent Protein (in Complex): "; 
		file << "NOT A COMPLEX\n";
	}
	
	if ( _i_root_one != NEXIST ) {
		file << line_start << "Index of Second Constituent Protein (in Complex): ";
		file << _i_root_one << "\n";
	}
	else {
		file << line_start << "Index of Second Constituent Protein (in Complex): "; 
		file << "NOT A COMPLEX\n";
	}
	
}

