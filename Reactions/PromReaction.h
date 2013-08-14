/*
 *  PromReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */


# include "Reaction.h"

#ifndef PROMREACTION_H
#define PROMREACTION_H

/* These reactions are of the form:
 *
 *		gene --> protein
 *		
 *		reaction: protein produced at rate kinetic*[gene]
 *
 */

class PromReaction : public Reaction {
	
public:
	PromReaction( int i_gene , int i_prot , double kinetic );
	~PromReaction() {}
	
	virtual Reaction* copy();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& currTissue , dmat& dx_dt , int i_curr_cell );
	virtual void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q );	
	
	virtual void update_indices( int first_index , int num_insertion );
	
	virtual void mutate ( boost::random::mt19937& generator );
	
	virtual void print_info ( std::string line_start );
	
	virtual void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_gene;
	int _i_prot;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;
	
	PromReaction(PromReaction* newOne) {}
	PromReaction& operator=( const PromReaction& rhs ) {return *this;}

};

# endif