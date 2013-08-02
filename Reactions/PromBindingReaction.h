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
	PromBindingReaction( int i_root_gene , int i_promoted_gene , int i_bound_protein ,
						double forward_kinetic , double backward_kinetic );
	~PromBindingReaction();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt ,
					    std::vector< std::vector<int>* >& neighbors,
					    int i_curr_cell );
	
	virtual void update_indices( int first_index , int num_insertion );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_root_gene;
	int _i_promoted_gene;
	int _i_bound_protein;
	
	/* Kinetic constants for the reaction. */
	double _forward_kinetic;
	double _backward_kinetic;
	
};

# endif

