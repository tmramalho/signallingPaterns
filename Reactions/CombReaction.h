/*
 *  CombReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */


# include "Reaction.h"

#ifndef COMBREACTION_H
#define COMBREACTION_H

/* These reactions are of the form:
 *
 *		reacZero + reacOne <--> product
 *		
 *		forwardReaction: reactant zero and reactant one combine to
 *		form product with rate forwardKinetic * [reacZero] * [reacOne].
 *		
 *		backwardReaction: product degrades into
 *		reactant zero and reactant one with rate backwardKinetic *
 *		[product].
 *
 */

class CombReaction : public Reaction {

public:
	CombReaction();
	CombReaction( int iReacZero , int iReacOne , int iProduct ,
				 double dxReacZero , double dxReacOne , double dxProduct ,
				 double forwardKinetic , double backwardKinetic );
	~CombReaction();

	virtual int getIPart( int partNum );
	virtual double getDx( int partNum );
	
	virtual void react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
						int iCurrCell , IntegrationType mode , double dt );
	
	virtual void updateIndices( int firstIndex , int numInsertions);
	
private:
	
	/* Where in the ODEManager are the participants located */
	int iReacZero;
	int iReacOne;
	int iProduct;
	
	/* What is the current rate of change of concentrations for the participants
	 * in this reaction */
	double dxReacZero;
	double dxReacOne;
	double dxProduct;
	
	/* Kinetic constants for the reaction. */
	double forwardKinetic;
	double backwardKinetic;
};

# endif