/*
 *  Reaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

/* Reaction class
 * -------------------------------------------------------------------------- 
 * The Reaction class is an abstract class.
 *
 
 * All reactions have a number of participants, which is the number of
 * molecules whose concentration gets changed due to the reaction.
 *
 * CAUTION: This is different than the number of molecules whose 
 * concentrations play a role in the rate of the reaction.
 *
 */

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <vector>

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>


# include "../old/numeric/dvec.h"
# include "../old/numeric/dmat.h"
# include "ReactionType.h"
# include "../IntegrationType.h"

#ifndef REACTION_H
#define REACTION_H

class Reaction {

public:

	~Reaction();
	
	int get_num_part();
	ReactionType get_type();
	virtual int get_i_part(int part_num) = 0;
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , 
					    std::vector< std::vector<int>* >& neighbors ,
					    int i_curr_cell ) = 0;
	
	virtual void update_indices( int first_index , int num_insertion) = 0;
	
	virtual void mutate( boost::random::mt19937& generator ) = 0;
	
	virtual void print_info ( std::string line_start ) = 0;
	
protected:
	
	Reaction();
	
	/* For reactions with less than three participants, all indexes and 
	 * locations for out-of-scope participants will be set to NEXIST.
	 */
	static const int NEXIST = -1;
	
	/* Reaction will calculate rates differently depending on its type */
	ReactionType _type;
	
	int _num_part;
	
private:
	
	/* Copy and assignment operators made private */
	Reaction( Reaction *newOne ) {}
	Reaction& operator=( const Reaction& rhs ) { return *this; }
	
};

#endif