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
	CombReaction( int i_reac_zero , int i_reac_one , int i_product ,
				 double forward_kinetic , double backward_kinetic );
	~CombReaction();

	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , 
					   std::vector< std::vector<int>* >& neighbors ,
					   int i_curr_cell );
	
	virtual void update_indices( int first_index , int num_insertions);
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_reac_zero;
	int _i_reac_one;
	int _i_product;
	
	/* Kinetic constants for the reaction. */
	double _forward_kinetic;
	double _backward_kinetic;
};

# endif