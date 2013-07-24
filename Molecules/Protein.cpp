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
Protein::Protein( int iSelf , int iRootZero , int iRootOne , double initConc ) {
	this->iRootZero = iRootZero;
	this->iRootOne = iRootOne;
	this->iSelf = iSelf;
	this->initConc = initConc;
	
	/* reactions will automatically initiate to empty vector */
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
int Protein::getIRootZero() {
	return iRootZero;
}

/* Public Method: getIRootOne()
 * -------------------------------------------------------------------------- 
 * If this protein is a complex, it returns the index in the dvecs of the 
 * second member protein in this complex. If it is not a complex, returns
 * NEXIST.
 */
int Protein::getIRootOne() {
	return iRootOne;
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

void Protein::updateIndices( int firstIndex , int numInsertions ) {
	if (iSelf >= firstIndex) iSelf += numInsertions;
	if (iRootZero >= firstIndex) iRootZero += numInsertions;
	if (iRootOne >= firstIndex) iRootOne += numInsertions;
}





