/*
 *  PromReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "PromReaction.h"

PromReaction::PromReaction() {
	
	type = PROMOTION;
	
	numPart = 1;
	
}

PromReaction::PromReaction( int iGene , int iProt , double dxProt , double kinetic ) {
	
	type = PROMOTION;
	
	numPart = 1;
	
	this->iGene = iGene;
	this->iProt = iProt;
	this->dxProt = dxProt;
	this->kinetic = kinetic;
	
}

PromReaction::~PromReaction() {}

int PromReaction::getIPart( int partNum ) {
	switch (partNum) {
		case 0:
			return iProt;
			break;
		default:
			return NEXIST;
			break;
	}
}

double PromReaction::getDx( int partNum ) {
	switch (partNum) {
		case 0:
			return dxProt;
			break;
		default:
			break;
	}
}

void PromReaction::react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
								 int iCurrCell , IntegrationType mode , double dt ) {
	switch (mode) {
			
		case RK1_DET_TI:
			dxProt = kinetic * currTissue.at(iCurrCell)->at(iGene) * dt;
			break;
			
		default:
			break;
			
	}
	
	
}

