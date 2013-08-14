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
	DegReaction( int i_reac , double kinetic );
	~DegReaction() {}
	
	virtual Reaction* copy();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell );
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q);	
	
	virtual void update_indices( int first_index , int num_insertion );
	
	virtual void mutate ( boost::random::mt19937& generator );
	
	virtual void print_info ( std::string line_start );
	
	virtual void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_reac;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;
	
	DegReaction(DegReaction* newOne) {}
	DegReaction& operator=( const DegReaction& rhs ) {return *this;}
	
};

# endif