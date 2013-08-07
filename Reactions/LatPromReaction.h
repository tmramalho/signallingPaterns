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
 * Where [proteinZero] is the concentration in the current
 * cell and [proteinOne] is the average concentration in its
 * neighbors.
 *
 */

class LatPromReaction : public Reaction {
	
public:
	LatPromReaction();
	LatPromReaction( int i_local_prot , int i_neighbor_prot ,
					 double kinetic ,
					 double K );
	LatPromReaction(LatPromReaction* newOne);
	~LatPromReaction();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , 
					    std::vector< std::vector<int>* >& neighbors,
					    int i_curr_cell );
	
	virtual void update_indices( int first_index , int num_insertion );
	
	virtual void mutate ( boost::random::mt19937& generator );
	
	virtual void print_info ( std::string line_start );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_local_prot;
	int _i_neighbor_prot;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;
	double _K;
	
};

# endif

