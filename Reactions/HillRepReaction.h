/*
 *  HillDegReaction.h
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#ifndef HILLREPREACTION_H
#define HILLREPREACTION_H

# include "Reaction.h"

/* In these reactions one protein's presence causes the repression
 * of another protein.
 *
 *		A repressed B
 *
 *		(d/dt)[B] = kinetic * K / ( K + [A]^n )
 *
 */

class HillRepReaction : public Reaction {
	
public:
	HillRepReaction( int i_repressor , int i_repressed ,
					 double kinetic , double K , double cooperativity );
	HillRepReaction( std::ifstream& file );
	~HillRepReaction() {}
	
	Reaction* copy();
	
	int get_i_part( int part_num ) const ;
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
	
	int _i_repressor;
	int _i_repressed;
	
	double _kinetic;
	double _K;
	double _cooperativity;
	
	/* File Format Variables */
	static const int NUM_LINE =  5;
	static const int I_REPRESSOR_LINE = 0;
	static const int I_REPRESSED_LINE = 1;
	static const int KINETIC_LINE = 2;
	static const int K_LINE = 3;
	static const int COOPERATIVITY_LINE = 4;	
	
	HillRepReaction(HillRepReaction* newOne) {}
	HillRepReaction& operator=( const HillRepReaction& rhs ) {return *this;}
	
};

#endif