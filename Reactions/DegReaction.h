/*
 *  DegReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */


# include "Reaction.h"

#ifndef DEGREACTION_H
#define DEGREACTION_H

/* These reactions are of the form:
 *
 *		reactant --> NOTHING
 *		
 *		forwardReaction: reactant degrades at rate -kinetic * [reactant]
 *
 */

class DegReaction : public Reaction {
	
public:
	DegReaction();
	DegReaction( int iReac , double dxReac , double kinetic );
	~DegReaction();
	
	virtual int getIPart( int partNum );
	virtual double getDx( int partNum );
	
	virtual void react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
					   int iCurrCell , IntegrationType mode , double dt );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int iReac;
	
	/* What is the current rate of change of concentrations for the participants
	 * in this reaction */
	double dxReac;
	
	/* Kinetic constants for the reaction. */
	double kinetic;
	
};

# endif