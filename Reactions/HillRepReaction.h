/*
 *  HillDegReaction.h
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

# include "Reaction.h"

#ifndef HILLREPREACTION_H
#define HILLREPREACTION_H

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
	~HillRepReaction() {}
	
	virtual Reaction* copy();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell );
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q );
	
	virtual void update_indices( int first_index , int num_insertion );
	
	virtual void mutate ( boost::random::mt19937& generator );
	
	virtual void print_info ( std::string line_start );
	
	virtual void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	int _i_repressor;
	int _i_repressed;
	
	double _kinetic;
	double _K;
	double _cooperativity;
	
	HillRepReaction(HillRepReaction* newOne) {}
	HillRepReaction& operator=( const HillRepReaction& rhs ) {return *this;}
	
};

#endif