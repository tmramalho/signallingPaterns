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
	PromReaction( int i_gene , int i_prot , double kinetic );
	~PromReaction();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& currTissue , dmat& dx_dt ,
					    std::vector< std::vector<int>* >& neighbors,
					    int i_curr_cell );
	
	virtual void update_indices( int first_index , int num_insertion );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_gene;
	int _i_prot;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;

};

# endif