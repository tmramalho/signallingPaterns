/*
 *  HillPromReaction.h
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#ifndef HILLPROMREACTION_H
#define HILLPROMREACTION_H

# include "Reaction.h"

/* In these reactions one protein's presence causes the production
 * of another protein according to a Hill function.
 *
 *		A promotes B
 *
 *		(d/dt)[B] = kinetic * ( [A]^n ) / ( K + [A]^n )
 *
 */

class HillPromReaction : public Reaction {

public:
	HillPromReaction( int i_promoter , int i_promoted ,
					 double kinetic , double K , double cooperativity );
	HillPromReaction( std::ifstream& file );
	~HillPromReaction() {}
	
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
	
	int _i_promoter;
	int _i_promoted;
	
	double _kinetic;
	double _K;
	double _cooperativity;
	
	/* File Format Variables */
	static const int NUM_LINE =  5;
	static const int I_PROMOTER_LINE = 0;
	static const int I_PROMOTED_LINE = 1;
	static const int KINETIC_LINE = 2;
	static const int K_LINE = 3;
	static const int COOPERATIVITY_LINE = 4;

	HillPromReaction(HillPromReaction* newOne);
	HillPromReaction& operator=( const HillPromReaction& rhs) {return *this;}
};

#endif

