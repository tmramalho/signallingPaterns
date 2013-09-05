/*
 *  Manager.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */


# ifndef _Manager_H
# define _Manager_H

# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <set>

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>
# include <boost/random/normal_distribution.hpp>

# include "old/numeric/dvec.h"
# include "old/numeric/dmat.h"
# include "Reactions/Reaction.h"
# include "Reactions/CombReaction.h"
# include "Reactions/DegReaction.h"
# include "Reactions/PromBindingReaction.h"
# include "Reactions/PromReaction.h"
# include "Reactions/LatPromReaction.h"
# include "Reactions/LatRepReaction.h"
# include "Reactions/PromBindingReaction.h"
# include "Reactions/HillPromReaction.h"
# include "Reactions/HillRepReaction.h"
# include "IntegrationType.h"
# include "ConstructionMethod.h"
# include "ReactionType.h"
# include "Molecules/Molecule.h"
# include "Molecules/Protein.h"
# include "Molecules/Gene.h"
# include "old/helpers/SettingsCont.h"
# include "nexist.cpp"
#include "Operations.h"
#include "ManagerDebugger.h"

using namespace std;

/* Manager class
 * -------------------------------------------------------------------------- 
 * This class holds all the necessary information for a particular genome 
 * and the integration of the corresponding ODE.
 *
 */

class Manager {
	
public:
	
	Manager( ConstructionMethod m );
	Manager( Manager *newOne );
	Manager( std::ifstream& file );
	~Manager();
	
	/* Mutation */
	void mutate( boost::random::mt19937& generator );
	void mutate( boost::random::mt19937& generator , int generation );
	
	/* Functions for running the ODE */
	void initialize();
	void integrate( int num_step , boost::random::mt19937& generator );
	void integrate( int num_step , boost::random::mt19937& generator ,
				   ofstream& file );
	
	/* Information Retrieval */
	void get_curr_state(dmat& curr_state);
	int get_num_cell() const;
	int get_num_gene() const;
	int get_num_prot() const;
	
	/* Display Functions */
	void print_state();
	void print_genome();
	
	/* File Writing */
	void genome_to_file( std::ofstream& file , std::string line_start );
	void state_to_file( std::ofstream& file );
	
	/* Debugging */
	bool has_comb_reac();
	SettingsCont *get_sc_ref();
	
	/* Lineage Data Tracking */
	void add_score(double score);
	void add_integration(int generation);
	
	/* Lineage Data Output */
	void scores_to_file( ofstream& file );
	void evaluation_generations_to_file( ofstream& file );
	void mutation_generations_to_file( ofstream& file );
	
	
private:
	
	/* PRIVATE DATA MEMBERS */
	
	/* Pointer to SettingsCont */
	SettingsCont *_sc_ref;
	
	/* State vectors */
	dmat _i_tissue, _curr_tissue;
	
	/* Gene and protein info */
	std::vector<Gene*> _genes;
	std::vector<Protein*> _proteins;
	
	/* Reaction info */
	std::vector< Reaction* > _reactions;
	
	/* Current Time */
	double _time;
	
	/* Number of molecules */
	int _num_mol;
	
	/* Number of cells */
	int _num_cell;
	
	/* Lineage Info */
	std::vector<int> _evaluation_generations;
	std::vector<double> _evaluation_scores;
	std::vector<int> _mutation_generations;
	
	/* PRIVATE METHODS */
	
	/* Assignment operator made private */
	Manager& operator=( const Manager& rhs );
	
	/* Helpers for running of ODEs */
	void react( dmat& xi , dmat& dx_dt);
	void react( dmat& xi , dmat& dx_dt ,
			   boost::random::normal_distribution<>& dist ,
			   boost::random::mt19937& generator , 
			   double q);
	
	/* Helpers for mutating genome */
	void resize_dmats();
	void update_mol_indices( int first_index , int num_insertion );
	void update_reac_indices( int first_index , int num_deletion );
	
	/* Helpers for adding and removing reactions */
	void add_gene( double prod_rate , double deg_rate );
	void add_PromBindingReac( int i_gene_in_genes , int i_prot_in_prots ,
							 double forward_kinetic , double backward_kinetic , 
							 double new_prod_rate );
	void add_CombReaction( int i_prot_zero_in_prots , int i_prot_one_in_prots , 
						  double forward_kinetic , double backward_kinetic );
	void add_LatPromReac( int i_promoting_neighbors_in_prots , int i_promoted_by_neighbors_in_prots ,
						 double kinetic , double K );
	void add_LatRepReac( int i_repressing_neighbors_in_prots , int i_repressed_by_neighbors_in_prots ,
						double kinetic , double K );
	void add_HillPromReac( int i_promoter_in_prots , int i_promoted_in_prots ,
						  double kinetic , double K , double cooperativity );
	void add_HillRepReac( int i_repressor_in_prots , int i_repressed_in_prots ,
						 double kinetic , double K , double cooperativity );
	void remove_gene( int i_gene );
	void remove_reaction( int i_reac );
	void reac_removal_cascade( int i_reac_to_remove ,
							  std::set<int>& reacs_to_delete ,
							  std::set<int>& genes_to_delete ,
							  std::set<int>& prots_to_delete );
	void gene_removal_cascade( int i_gene_to_delete ,
							 std::set<int>& reacs_to_delete ,
							 std::set<int>& genes_to_delete , 
							 std::set<int>& prots_to_delte );
	void prot_removal_cascade( int i_prot_to_delete , 
							  std::set<int>& reacs_to_delete ,
							  std::set<int>& genes_to_delete ,
							  std::set<int>& prots_to_delete);
	void perform_removal_using_sets(std::set<int>& reacs_to_delete , 
									std::set<int>& genes_to_delete , 
									std::set<int>& prots_to_delete);
	
	/* Helpers for performing specific mutations */
	void degredation_mutation ( boost::random::mt19937& generator );
	void kinetic_mutation ( boost::random::mt19937& generator );
	void prom_binding_mutation ( boost::random::mt19937& generator );
	void post_transcript_mutation ( boost::random::mt19937& generator );
	void add_hill_mutation ( boost::random::mt19937& generator );
	void removal_mutation ( boost::random::mt19937& generator );
	
	/* Helpers for integration */
	void rk4_det_ti ( int num_step );
	void rk4_det_ti ( int num_step , ofstream& file );
	void rk4_stc_ti ( int num_step , boost::random::mt19937& generator );
	void rk4_stc_ti ( int num_step , boost::random::mt19937& generator , ofstream& file );
	
	/* Helpers for the constructor */
	void hakim_delta_notch_construct();
	void collier_delta_notch_construct();
	void mutation_construct();
	void one_protein_construct();
	void two_protein_construct();
	void best_genome_construct();
	
	/* FRIEND CLASSES */
	
	/* Debugger */
	friend class ManagerDebugger;
	
};


# endif