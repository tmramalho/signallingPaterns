/*
 *  Gene.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Gene.h"

Gene::Gene( int iSelf, int iProduct , int iBoundPromoter , int iRoot , double initConc ) {
	this->iProduct = iProduct;
	this->iBoundPromoter = iBoundPromoter;
	this->iRoot = iRoot;
	this->iSelf = iSelf;
	this->initConc = initConc;
}

/* Destructor: ~Gene()
 * -------------------------------------------------------------------------- 
 * No heap allocated memory owned by the gene class.
 */
Gene::~Gene() {}

/* Public Method: getIProduct()
 * -------------------------------------------------------------------------- 
 * Returns the index of this gene in the dvecs held in the cell.
 */
int Gene::getIProduct() {
	return iProduct;
}

/* Public Method: getIBoundPromoter()
 * -------------------------------------------------------------------------- 
 * Returns the index in the dvecs of the protein that our bound to this 
 * gene as a promoter (or represser). If there is no such protein, returns
 * NEXIST.
 */
int Gene::getIBoundPromoter() {
	return iBoundPromoter;
}

/* Public Method: getIRoot()
 * -------------------------------------------------------------------------- 
 * Returns the index in the dvecs of the gene from which this gene derives,
 * if this gene is complexed with a protein. If it is not, it returns 
 * NEXIST.
 *
 * Example: For the gene-protein complex a:B, this would return the index of
 * a
 */
int Gene::getIRoot() {
	return iRoot;
}

/* Private Method: updateIndices(firstIndex,numInsertions)
 * -------------------------------------------------------------------------- 
 * Updates the index of iSelf, iProduct, iBoundPromoter, iRoot, given an 
 * insertion of size numInsertions, beginning at firstIndex, into our dvecs 
 * containing molecule concentrations in our manager.
 *
 * See description for updateIndices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void Gene::updateIndices( int firstIndex , int numInsertions ) {
	if ( iSelf >= firstIndex ) iSelf += numInsertions;
	if ( iProduct >= firstIndex ) iProduct += numInsertions;
	if ( iBoundPromoter >= firstIndex ) iBoundPromoter += numInsertions;
	if ( iRoot >= firstIndex ) iRoot += numInsertions;
}


