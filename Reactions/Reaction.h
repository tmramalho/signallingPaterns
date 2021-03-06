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
 * molecules whose concentration gets changed due to the reaction, NOT the 
 * number of molecules whos concentration is relevant to calculating reaction
 * rate.
 *
 * Example: For a promotion reaction, where a protein is produced at a rate
 * depending upon the concentation of a gene, the number of molecules is one,
 * because only the proteins concentration changes due to the reaction,
 * despite there being two molecules (gene and protein) whose concentration
 * is important to the reaction.
 *
 */


#ifndef REACTION_H
#define REACTION_H

# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <vector>
# include <algorithm>

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>
# include <boost/random/normal_distribution.hpp>


# include "../old/numeric/dvec.h"
# include "../old/numeric/dmat.h"
# include "ReactionType.h"
# include "../IntegrationType.h"
# include "../old/helpers/SettingsCont.h"
# include "../nexist.cpp"

class Reaction {

public:
	virtual ~Reaction();
	
	bool operator==( const Reaction& rhs ) const;
	bool operator!=( const Reaction& rhs ) const;
	
	int get_num_part() const;
	ReactionType get_type();
	
	virtual Reaction* copy() = 0;
	
	virtual int get_i_part( int part_num ) const = 0;
	virtual int get_i_dependent_molecule() const = 0;
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) = 0;
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q ) = 0;
	
	virtual void update_mol_indices( int first_index , int num_insertion) = 0;
	
	virtual void mutate( boost::random::mt19937& generator ) = 0;
	
	virtual void print_info ( std::string line_start ) = 0;
	
	virtual void to_file ( std::ofstream& file , std::string line_start ) = 0;
	
protected:
	
	Reaction();
	
	SettingsCont* _sc_ref;
	
	/* Reaction will calculate rates differently depending on its type */
	ReactionType _type;
	
	int _num_part;
	
	/* Update Index Helper */
	int update_index ( int first_index , int num_insertion , int index );
	
private:
	Reaction(Reaction *newOne);
	
};

#endif