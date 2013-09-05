/*
 *  LatPromReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#ifndef LATPROMREACTION_H
#define LATPROMREACTION_H

# include "Reaction.h"

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
	LatPromReaction( int i_promoting_neighbors , int i_promoted_by_neighbors ,
					 double kinetic , double K );
	LatPromReaction( std::ifstream& file );
	~LatPromReaction() {}
	
	Reaction* copy();
	
	int get_i_part( int part_num ) const;
	int get_i_dependent_molecule() const;
	
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell );
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q );
	
	void update_mol_indices( int first_index , int num_insertion );
	
	void mutate ( boost::random::mt19937& generator );
	
	void print_info ( std::string line_start );
	
	void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_promoted_by_neighbors;
	int _i_promoting_neighbors;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;
	double _K;
	
	/* File Format Variables */
	static const int NUM_LINE = 4;
	static const int I_PROMOTED_BY_NEIGHBORS_LINE = 0;
	static const int I_PROMOTING_NEIGHBORS_LINE = 1;
	static const int KINETIC_LINE = 2;
	static const int K_LINE = 3;
	
	LatPromReaction(LatPromReaction* newOne) {}
	LatPromReaction &operator=( const LatPromReaction& rhs ) {return *this;}
	
};

# endif

