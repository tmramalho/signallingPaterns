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
	CombReaction( int i_reac_zero , int i_reac_one , int i_product ,
				 double forward_kinetic , double backward_kinetic );
	~CombReaction() {}
	
	virtual Reaction* copy();

	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell );
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , 
					   double q );	
	
	virtual void update_indices( int first_index , int num_insertions);
	
	virtual void mutate( boost::random::mt19937& generator );
	
	virtual void print_info ( std::string line_start );
	
	virtual void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* At what index in the manager are the participants located */
	int _i_reac_zero;
	int _i_reac_one;
	int _i_product;
	
	/* Kinetic constants for the reaction. */
	double _forward_kinetic;
	double _backward_kinetic;
	
	CombReaction(CombReaction* newOne) {}
	CombReaction& operator=( const CombReaction& rhs ) {return *this;}
};

# endif