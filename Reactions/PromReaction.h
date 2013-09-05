/*
 *  PromReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#ifndef PROMREACTION_H
#define PROMREACTION_H

# include "Reaction.h"

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
	PromReaction( std::ifstream& file );
	~PromReaction();
	
	Reaction* copy();
	
	int get_i_part( int part_num ) const;
	int get_i_dependent_molecule() const;
	
	void react( dmat& currTissue , dmat& dx_dt , int i_curr_cell );
	void react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
					   boost::random::normal_distribution<>& dist , 
					   boost::random::mt19937& generator , double q );	
	
	void update_mol_indices( int first_index , int num_insertion );
	
	void mutate ( boost::random::mt19937& generator );
	
	void print_info ( std::string line_start );
	
	void to_file ( std::ofstream& file , std::string line_start );
	
private:
	
	/* Where in the ODEManager are the participants located */
	int _i_gene;
	int _i_prot;
	
	/* Kinetic constants for the reaction. */
	double _kinetic;
	
	/* File Format Variables */
	static const int NUM_LINE = 3;
	static const int I_GENE_LINE = 0;
	static const int I_PROT_LINE = 1;
	static const int KINETIC_LINE = 2;
	
	PromReaction(PromReaction* newOne) {}
	PromReaction& operator=( const PromReaction& rhs ) {return *this;}

};

# endif