/*
 *  LatPromReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "LatPromReaction.h"

LatPromReaction::LatPromReaction() {
	
	type = LATERAL_PROMOTION;
	
	numPart = 1;
	
}

LatPromReaction::LatPromReaction( int iLocalProt , int iNeighborProt ,
								 double dxLocalProt , double kinetic ,
								 double K ) {
	
	type = LATERAL_PROMOTION;
	
	numPart = 1;
	
	this->iLocalProt = iLocalProt;
	this->iNeighborProt = iNeighborProt;
	this->dxLocalProt = dxLocalProt;
	this->kinetic = kinetic;
	this->K = K;
	
}

LatPromReaction::~LatPromReaction() {}

int LatPromReaction::getIPart( int partNum ) {
	switch (partNum) {
		case 0:
			return iLocalProt;
			break;
		default:
			return NEXIST;
			break;
	}
}

double LatPromReaction::getDx( int partNum ) {
	switch (partNum) {
		case 0:
			return dxLocalProt;
			break;
		default:
			break;
	}
}

void LatPromReaction::react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
									int iCurrCell , IntegrationType mode , double dt ) {
	
	double neighborSum;
	double avgNeighborConc;
	
	switch (mode) {
			
		case RK1_DET_TI:
			
			neighborSum = 0.0;
			for (int i = 0; i < neighbors.at(iCurrCell)->size() ; i++) {
				neighborSum += currTissue.at(neighbors.at(iCurrCell)->at(i))->at(iNeighborProt);
			}
			avgNeighborConc = neighborSum/neighbors.at(iCurrCell)->size();
			
			dxLocalProt = kinetic * (pow(avgNeighborConc,2))/(K+pow(avgNeighborConc,2)) * dt;
			
			if ( abs(avgNeighborConc-.5) < .01 ) {
				int k = 2;
			}
			break;
			
		default:
			dxLocalProt = 0.0;
			break;
			
	}
}



