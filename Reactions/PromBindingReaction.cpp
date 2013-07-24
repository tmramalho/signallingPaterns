/*
 *  PromBindingReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "PromBindingReaction.h"

PromBindingReaction::PromBindingReaction() {
	
	type = PROMOTER_BINDING;
	
	numPart = 2;
	
}

PromBindingReaction::PromBindingReaction( int iRootGene , int iPromotedGene , int iBoundProtein ,
					double dxLoneGene , double dxPromotedGene ,
					double forwardKinetic , double backwardKinetic ) {
	
	type = PROMOTER_BINDING;
	
	numPart = 2;
	
	this->iRootGene = iRootGene;
	this->iPromotedGene = iPromotedGene;
	this->iBoundProtein = iBoundProtein;
	this->dxRootGene = dxRootGene;
	this->dxPromotedGene = dxPromotedGene;
	this->forwardKinetic = forwardKinetic;
	this->backwardKinetic = backwardKinetic;
	
}

PromBindingReaction::~PromBindingReaction() {}

int PromBindingReaction::getIPart( int partNum ) {
	switch (partNum) {
		case 0:
			return iRootGene;
			break;
		case 1:
			return iPromotedGene;
		default:
			return NEXIST;
			break;
	}
}

double PromBindingReaction::getDx( int partNum ) {
	switch (partNum) {
		case 0:
			return dxRootGene;
			break;
		case 1:
			return dxPromotedGene;
			break;
		default:
			break;
	}
}

void PromBindingReaction::react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
									int iCurrCell , IntegrationType mode , double dt ) {
	
	double forwardFlow;
	double backwardFlow;
	
	switch (mode) {
			
		case RK1_DET_TI:
			forwardFlow = forwardKinetic * 
			currTissue.at(iCurrCell)->at(iRootGene) *
			currTissue.at(iCurrCell)->at(iBoundProtein);
			
			backwardFlow = backwardKinetic *
			currTissue.at(iCurrCell)->at(iPromotedGene);
			
			dxRootGene = (backwardFlow - forwardFlow) * dt;
			dxPromotedGene = (forwardFlow - backwardFlow) * dt;
			break;
			
		default:
			dxRootGene = 0.0;
			dxPromotedGene = 0.0;
			break;
			
	}
}

/* Public Method: updateIndices(firstIndex,numInsertions)
 * -------------------------------------------------------------------------- 
 * Updates the index of iRootGene, iPromotedGene, iBoundProtein, given an 
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

void PromBindingReaction::updateIndices( int firstIndex , int numInsertions ) {
	if (iRootGene >= firstIndex) iRootGene += numInsertions;
	if (iPromotedGene >= firstIndex) iPromotedGene += numInsertions;
	if (iBoundProtein >= firstIndex) iBoundProtein += numInsertions;
}



