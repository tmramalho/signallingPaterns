/*
 *  DegReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "DegReaction.h"

DegReaction::DegReaction() {
	
	type = DEGRADATION;
	
	numPart = 1;
	
}

DegReaction::DegReaction( int iReac , double dxReac , double kinetic ) {
	
	type = DEGRADATION;
	
	numPart = 1;
	
	this->iReac = iReac;
	this->dxReac = dxReac;
	this->kinetic = kinetic;
	
}

DegReaction::~DegReaction() {}

int DegReaction::getIPart( int partNum ) {
	switch (partNum) {
		case 0:
			return iReac;
			break;
		default:
			return NEXIST;
			break;
	}
}

double DegReaction::getDx( int partNum ) {
	switch (partNum) {
		case 0:
			return dxReac;
			break;
		default:
			break;
	}
}

void DegReaction::react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
								int iCurrCell , IntegrationType mode , double dt ) {
	switch (mode) {
		
		case RK1_DET_TI:
			dxReac = - kinetic * dt * currTissue.at(iCurrCell)->at(iReac);
			break;
		
		default:
			break;
	}
	
	
}

/* Public Method: updateIndices(firstIndex,numInsertions)
 * -------------------------------------------------------------------------- 
 * Updates the index of iReac, given an insertion of size numInsertions, 
 * beginning at firstIndex, into our dvecs containing molecule concentrations 
 * in our manager.
 *
 * See description for updateIndices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void DegReaction::updateIndices( int firstIndex , int numInsertions ) {
	if (iReac >= firstIndex) iReac += numInsertions;
}


