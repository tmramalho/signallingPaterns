/*
 *  PromReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */


# include "Reaction.h"

#ifndef PROMREACTION_H
#define PROMREACTION_H

/* These reactions are of the form:
 *
 *		gene --> protein
 *		
 *		reaction: protein produced at rate kinetic*[gene]
 *
 */

class PromReaction : public Reaction {
	
public:
	PromReaction();
	PromReaction( int iGene , int iProt , double dxProt , double kinetic );
	~PromReaction();
	
	virtual int getIPart( int partNum );
	virtual double getDx( int partNum );
	
	virtual void react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
					   int iCurrCell , IntegrationType mode , double dt );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int iGene;
	int iProt;
	
	/* What is the current rate of change of concentrations for the participants
	 * in this reaction */
	double dxProt;
	
	/* Kinetic constants for the reaction. */
	double kinetic;
	
};

# endif