/*
 *  Protein.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Protein.h"

/* Constructor: Protein()
 * -------------------------------------------------------------------------- 
 */
Protein::Protein( int i_self , int i_root_zero , int i_root_one , double init_conc ) {
	_i_root_zero = i_root_zero;
	_i_root_one = i_root_one;
	_i_self = i_self;
	_init_conc = init_conc;
	
	/* reactions will automatically initiate to empty vector */
}

Protein::Protein(Protein* newOne) {
	_i_root_zero = newOne->_i_root_zero;
	_i_root_one = newOne->_i_root_one;
}

/* Destructor: ~Protein()
 * -------------------------------------------------------------------------- 
 * No heap allocated memory owned by the protein class.
 */
Protein::~Protein() {}

/* Public Method: getIRootZero()
 * -------------------------------------------------------------------------- 
 * If this protein is a complex, it returns the index in the dvecs of the 
 * first member protein in this complex. If it is not a complex, returns
 * NEXIST.
 */
int Protein::get_i_root_zero() {
	return _i_root_zero;
}

/* Public Method: getIRootOne()
 * -------------------------------------------------------------------------- 
 * If this protein is a complex, it returns the index in the dvecs of the 
 * second member protein in this complex. If it is not a complex, returns
 * NEXIST.
 */
int Protein::get_i_root_one() {
	return _i_root_one;
}

/* Private Method: updateIndices(firstIndex,numInsertions)
 * -------------------------------------------------------------------------- 
 * Updates the index of iSelf, iRootZero, and iRootOne, given an insertion
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

void Protein::update_indices( int first_index , int num_insertion ) {
	if (_i_self >= first_index) _i_self += num_insertion;
	if (_i_root_zero >= first_index) _i_root_zero += num_insertion;
	if (_i_root_one >= first_index) _i_root_one += num_insertion;
}





