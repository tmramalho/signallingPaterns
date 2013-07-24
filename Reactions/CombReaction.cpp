/*
 *  CombReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "CombReaction.h"

CombReaction::CombReaction() {
	
	type = COMBINATION;
	
	numPart = 3;
	
}

CombReaction::CombReaction( int iReacZero , int iReacOne , int iProduct ,
			 double dxReacZero , double dxReacOne , double dxProduct ,
			 double forwardKinetic , double backwardKinetic ) {
	
	type = COMBINATION;
	
	numPart = 3;
	
	this->iReacZero = iReacZero;
	this->iReacOne = iReacOne;
	this->iProduct = iProduct;
	this->dxReacZero = dxReacZero;
	this->dxReacOne = dxReacOne;
	this->dxProduct = dxProduct;
	this->forwardKinetic = forwardKinetic;
	
}

CombReaction::~CombReaction() {}

int CombReaction::getIPart( int partNum ) {
	switch (partNum) {
		case 0:
			return iReacZero;
			break;
		case 1:
			return iReacOne;
			break;
		case 2:
			return iProduct;
			break;
		default:
			return NEXIST;
			break;
	}
}

double CombReaction::getDx( int partNum ) {
	switch (partNum) {
		case 0:
			return dxReacZero;
			break;
		case 1:
			return dxReacOne;
			break;
		case 2:
			return dxProduct;
			break;
		default:
			break;
	}
}

void CombReaction::react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
								int iCurrCell , IntegrationType mode , double dt ) {
	double forwardFlow;
	double backwardFlow;
	
	switch (mode) {
			
		case RK1_DET_TI:
			/* Find forward and backward flows */
			forwardFlow = forwardKinetic * 
			currTissue.at(iCurrCell)->at(iReacZero) *
			currTissue.at(iCurrCell)->at(iReacOne);
			
			backwardFlow = backwardKinetic *
			currTissue.at(iCurrCell)->at(iProduct);
			
			/* Update our dxdt values */
			dxReacZero = dxReacOne = (backwardFlow - forwardFlow) * dt;
			dxProduct = (forwardFlow - backwardFlow) * dt;
			break;
			
		default:
			break;
	}	
	
}

/* Public Method: updateIndices(firstIndex,numInsertions)
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

void CombReaction::updateIndices( int firstIndex , int numInsertions ) {
	if (iReacZero >= firstIndex) iReacZero += numInsertions;
	if (iReacOne >= firstIndex) iReacOne += numInsertions;
	if (iProduct >= firstIndex) iProduct += numInsertions;
}


