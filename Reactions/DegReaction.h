/*
 *  DegReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#ifndef DEGREACTION_H
#define DEGREACTION_H

# include "Reaction.h"

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
	DegReaction( std::ifstream& file );
	~DegReaction();
	
	Reaction* copy();
	
	int get_i_part( int part_num ) const;
	int get_i_dependent_molecule() const;
	
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell );
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q);	
	
	void update_mol_indices( int first_index , int num_insertion );
	
	void mutate ( boost::random::mt19937& generator );
	
	void print_info ( std::string line_start );
	
	void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_reac;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;
	
	/* File Format Variables */
	static const int NUM_LINE = 2;
	static const int I_REAC_LINE = 0;
	static const int KINETIC_LINE = 1;
	
	DegReaction(DegReaction* newOne) {}
	DegReaction& operator=( const DegReaction& rhs ) {return *this;}
	
};

# endif