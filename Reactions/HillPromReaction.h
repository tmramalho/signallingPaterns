/*
 *  HillPromReaction.h
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

# include "Reaction.h"

#ifndef HILLPROMREACTION_H
#define HILLPROMREACTION_H

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
	~HillPromReaction() {}
	
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
	
	int _i_promoter;
	int _i_promoted;
	
	double _kinetic;
	double _K;
	double _cooperativity;

	HillPromReaction(HillPromReaction* newOne);
	HillPromReaction& operator=( const HillPromReaction& rhs) {return *this;}
};

#endif

