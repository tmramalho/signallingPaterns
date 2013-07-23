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

PromBindingReaction::PromBindingReaction( int iLoneGene , int iPromotedGene , int iProt ,
					double dxLoneGene , double dxPromotedGene ,
					double forwardKinetic , double backwardKinetic ) {
	
	type = PROMOTER_BINDING;
	
	numPart = 2;
	
	this->iLoneGene = iLoneGene;
	this->iPromotedGene = iPromotedGene;
	this->iProt = iProt;
	this->dxLoneGene = dxLoneGene;
	this->dxPromotedGene = dxPromotedGene;
	this->forwardKinetic = forwardKinetic;
	this->backwardKinetic = backwardKinetic;
	
}

PromBindingReaction::~PromBindingReaction() {}

int PromBindingReaction::getIPart( int partNum ) {
	switch (partNum) {
		case 0:
			return iLoneGene;
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
			return dxLoneGene;
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
			currTissue.at(iCurrCell)->at(iLoneGene) *
			currTissue.at(iCurrCell)->at(iProt);
			
			backwardFlow = backwardKinetic *
			currTissue.at(iCurrCell)->at(iPromotedGene);
			
			dxLoneGene = (backwardFlow - forwardFlow) * dt;
			dxPromotedGene = (forwardFlow - backwardFlow) * dt;
			break;
			
		default:
			dxLoneGene = 0.0;
			dxPromotedGene = 0.0;
			break;
			
	}
}



