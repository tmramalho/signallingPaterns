/*
 *  LatPromReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

# include "Reaction.h"

#ifndef LATPROMREACTION_H
#define LATPROMREACTION_H

/* These reactions, the concentration of a protein in a cell
 * increases due to the average concentration of a proteins 
 * in neighboring cells.
 *		
 *		(d/dt)[proteinZero]=
 *				kinetic*([proteinOne]^2)/(K+[proteinOne]^2)
 *
 * Where [proteinZero} is the concentration in the current
 * cell and [proteinOne] is the average concentration in its
 * neighbors.
 *
 */

class LatPromReaction : public Reaction {
	
public:
	LatPromReaction();
	LatPromReaction( int iLocalProt , int iNeighborProt ,
					double dxLocalProt , double kinetic ,
					double K );
	~LatPromReaction();
	
	virtual int getIPart( int partNum );
	virtual double getDx( int partNum );
	
	virtual void react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
					   int iCurrCell , IntegrationType mode , double dt );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int iLocalProt;
	int iNeighborProt;
	
	/* What is the current rate of change of concentrations for the participants
	 * in this reaction */
	double dxLocalProt;
	
	/* Kinetic constants for the reaction. */
	double kinetic;
	double K;
	
};

# endif

