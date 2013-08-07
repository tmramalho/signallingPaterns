/*
 *  Manager.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>

# include "old/numeric/dvec.h"
# include "old/numeric/dmat.h"
# include "Reactions/Reaction.h"
# include "Reactions/CombReaction.h"
# include "Reactions/DegReaction.h"
# include "Reactions/PromBindingReaction.h"
# include "Reactions/PromReaction.h"
# include "Reactions/LatPromReaction.h"
# include "Reactions/PromBindingReaction.h"
# include "Reactions/HillPromReaction.h"
# include "Reactions/HillRepReaction.h"
# include "IntegrationType.h"
# include "ReactionType.h"
# include "Molecules/Molecule.h"
# include "Molecules/Protein.h"
# include "Molecules/Gene.h"

# ifndef _Manager_H
# define _Manager_H

using namespace std;

/* Manager class
 * -------------------------------------------------------------------------- 
 * This class holds all the necessary information for a particular genome 
 * and the integration of the corresponding ODE.
 *
 * It also constantly stores information about which integration mode it
 * it should use to integrate, the current time in the integration, and the
 * time step it should use to integrate.
 *
 * This means that we just call integrate() for the integration, and must
 * change these settings before calling integrate() if we want to manage
 * them.
 *
 */

class Manager {
	
public:
	
	Manager();
	Manager( string method );
	~Manager();
	
	void set_mode( IntegrationType mode );
	IntegrationType get_mode();
	
	void set_dt( double dt );
	double get_dt();
	
	/* Mutation */
	void mutate( boost::random::mt19937& generator );
	
	/* Functions for running the ODE */
	void initialize();
	void integrate( int num_step );
	
	/* Display Functions */
	void print_state();
	void print_genome();
	
private:
	
	static const int NEXIST = -1;
	
	/* Neighbors info */
	/* Our tissue data is held in a vector, with each element a separate
	 * cell. Because our genetic network will include membrane reactions,
	 * we would like each cell to be able to know who its neighbors are.
	 * Each element of this vector holds a vector of the indices of the
	 * corresponding cells.
	 *
	 * For example, if we have a 2x2 tissue, we might index the cells
	 * from left to right, top to bottom, as
	 *
	 *		1	2
	 *
	 *		3	4
	 *
	 * Then, our vector of neighbors would look like this:
	 *
	 *		< [2,3] , [1,4] , [1,4] , [2,3] >
	 *
	 * This is how we store the geometry of our lattice, and it can 
	 * accomodate any geometry we desire.
	 * Any use can then easily use the manager class to study any geometry
	 * he or she wants.
	 *
	 */
	std::vector< std::vector<int>* > _neighbors;
	
	/* State vectors */
	/* The class is implemented so that these always up-to-date. In particular,
	 * no public method should return before the dxdt vector has the correct
	 * rates given the state of currTissue. Then, to get the current rates of
	 * change in our ODE, we only need to ask the ODEManager to look in the
	 * dxdt vector. All calculations should have already been done.
	 */
	dmat _i_tissue, _curr_tissue;
	
	/* Gene and protein info */
	std::vector<Gene*> _genes;
	std::vector<Protein*> _proteins;
	
	/* Time variables */
	double _dt; 
	double _time;
	
	/* Number of molecules */
	int _num_mol;
	
	/* Number of cells */
	int _num_cell;
	
	/* Vector of Reactions */
	std::vector< Reaction* > _reactions;
	
	/* Integration Mode */
	IntegrationType _mode;
	
	/* Copy and assignment operators made private */
	Manager( Manager *newOne );
	Manager& operator=( const Manager& rhs );
	
	/* Helpers for running of ODEs */
	void react( dmat& xi , dmat& dx_dt);
	
	/* Helpers for mutating genome */
	void resize_dmats();
	void update_indices( int first_index , int num_insertion );
	
	/* Adding Reactions Helpers */
	void add_gene( double prod_rate , double deg_rate );
	void add_PromBindingReac( int i_gene_in_genes , int i_prot_in_prots ,
							 double forward_kinetic , double backward_kinetic , 
							 double new_prod_rate );
	//void addPhosphReac();
	//void addPartialDegReac();
	void add_CombReaction( int i_prot_zero_in_prots , int i_prot_one_in_prots , 
						  double forward_kinetic , double backward_kinetic );
	//void addPartialCatDegReac();
	void add_LatPromReac( int i_loc_prot_in_prots , int i_neighbor_prot_in_prots ,
						 double kinetic , double K );
	void add_HillPromReac( int i_promoter_in_prots , int i_promoted_in_prots ,
						  double kinetic , double K , double cooperativity );
	void add_HillRepReac( int i_repressor_in_prots , int i_repressed_in_prots ,
						 double kinetic , double K , double cooperativity );
	
	/* Removing Reactions */
	//void remove_reaction( int i_reac );
	
	/* Mutation Types: Helper Functions */
	void degredation_mutation ( boost::random::mt19937& generator );
	void kinetic_mutation ( boost::random::mt19937& generator );
	void prom_binding_mutation ( boost::random::mt19937& generator );
	void post_transcript_mutation ( boost::random::mt19937& generator );
	
	/* Integration Types Methods */
	void rk4_det_ti ( int num_step );
	void rk4_stc_ti ( int num_step );
	
};


# endif