/*
 *  PromBindingReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

# include "Reaction.h"

#ifndef PROMBINDINGREACTION_H
#define PROMBINDINGREACTION_H

/* These reactions, a gene turns into a gene bound to a protein, or
 * a gene bound to a protein turns into the lone gene. The concentration
 * of the protein, however, does not change, because genes are in such
 * low concentration in the cell.
 *
 *		a <--> a:B
 *
 *		(d/dt)[a:B] = forwardKinetic * [a][B] - backwardKinetic * [a:B]
 *		(d/dt)[a] = backwardKinetic * [a:B] - forwardKinetic * [a][B]
 *
 */

class PromBindingReaction : public Reaction {
	
public:
	PromBindingReaction();
	PromBindingReaction( int iRootGene , int iPromotedGene , int iBoundProtein ,
						double dxRootGene , double dxPromotedGene ,
						double forwardKinetic , double backwardKinetic );
	~PromBindingReaction();
	
	virtual int getIPart( int partNum );
	virtual double getDx( int partNum );
	
	virtual void react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
					   int iCurrCell , IntegrationType mode , double dt );
	
	virtual void updateIndices( int firstIndex , int numInsertions );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int iRootGene;
	int iPromotedGene;
	int iBoundProtein;
	
	/* What is the current rate of change of concentrations for the participants
	 * in this reaction */
	double dxRootGene;
	double dxPromotedGene;
	
	/* Kinetic constants for the reaction. */
	double forwardKinetic;
	double backwardKinetic;
	
};

# endif

