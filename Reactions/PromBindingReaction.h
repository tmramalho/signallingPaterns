/*
 *  PromBindingReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#ifndef PROMBINDINGREACTION_H
#define PROMBINDINGREACTION_H

# include "Reaction.h"

/* These reactions, a gene turns into a gene bound to a protein, or
 * a gene bound to a protein turns into the lone gene. The concentration
 * of the protein, however, does not change, because genes are in such
 * low concentration in the cell.
 *
 *		a <--> a:B
 *
 *		(d/dt)[a:B] = forwardKinetic * [a][B] - backwardKinetic * [a:B]
 *		(d/dt)[a] = backwardKinetic * [a:B] - forwardKinetic * [a][B]
 *
 */

class PromBindingReaction : public Reaction {
	
public:
	PromBindingReaction( int i_root_gene , int i_promoted_gene , int i_bound_protein ,
						double forward_kinetic , double backward_kinetic );
	PromBindingReaction( std::ifstream& file );
	~PromBindingReaction() {}
	
	Reaction* copy();
	
	int get_i_part( int part_num ) const;
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
	
	/* Where in the ODEManager are the participants located */
	int _i_root_gene;
	int _i_bound_protein;
	int _i_promoted_gene;
	
	/* Kinetic constants for the reaction. */
	double _forward_kinetic;
	double _backward_kinetic;
	
	/* File Format Variables */
	static const int NUM_LINE = 5;
	static const int I_ROOT_GENE_LINE = 0;
	static const int I_BOUND_PROTEIN_LINE = 1;
	static const int I_PROMOTED_GENE_LINE = 2;
	static const int FORWARD_KINETIC_LINE = 3;
	static const int BACKWARD_KINETIC_LINE = 4;
	
	PromBindingReaction(PromBindingReaction* newOne) {}
	PromBindingReaction& operator=( const PromBindingReaction& rhs ) {return *this;}
	
};

# endif

