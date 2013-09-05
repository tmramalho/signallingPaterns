/*
 *  CombReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#ifndef COMBREACTION_H
#define COMBREACTION_H

# include "Reaction.h"

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
	CombReaction( int i_reac_zero , int i_reac_one , int i_product ,
				 double forward_kinetic , double backward_kinetic );
	CombReaction( std::ifstream& file );
	~CombReaction() {}
	
	Reaction* copy();

	int get_i_part( int part_num ) const;
	int get_i_dependent_molecule() const;
	
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell );
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , 
					   double q );	
	
	void update_mol_indices( int first_index , int num_insertions);
	
	void mutate( boost::random::mt19937& generator );
	
	void print_info ( std::string line_start );
	
	void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* At what index in the manager are the participants located */
	int _i_reac_zero;
	int _i_reac_one;
	int _i_product;
	
	/* Kinetic constants for the reaction. */
	double _forward_kinetic;
	double _backward_kinetic;
	
	/* File Format Variables */
	static const int NUM_LINE = 5;
	static const int I_REAC_ZERO_LINE = 0;
	static const int I_REAC_ONE_LINE = 1;
	static const int I_PRODUCT_LINE = 2;
	static const int FORWARD_KINETIC_LINE = 3;
	static const int BACKWARD_KINETIC_LINE = 4;
	
	CombReaction(CombReaction* newOne) {}
	CombReaction& operator=( const CombReaction& rhs ) {return *this;}
	
};

# endif