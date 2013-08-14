/*
 *  Manager.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

#include "Manager.h"

/* Constructor: Manager(method) 
 * --------------------------------------------------------------------------
 * Constructor that allows different construction methods, indicated by the 
 * string method.
 */
Manager::Manager( ConstructionMethod m ) {
	
	_sc_ref = SettingsCont::getInstance();
	_num_cell = _sc_ref->_neighbors.size();
	
	try {
		if (m == HAKIM_DELTA_NOTCH) {
			hakim_delta_notch_construct();
		}
		else if (m == COLLIER_DELTA_NOTCH) {
			collier_delta_notch_construct();
		}
		else if (m == MUTATION) {
			mutation_construct();
		}
		else if (m == ONE_PROTEIN) {
			one_protein_construct();
		}
		else { throw -1; }
	}
	catch ( int i ) {
		std::cout << "Invalid construction method\n";
	}
	
}

/* Copy Constructor: Manager(newOne) 
 * --------------------------------------------------------------------------
 * WARNING: We want our Manager to always be consistent with the settings
 * in SettingsCont. We will copy directly from *newOne. If newOne is not
 * consistent (ie, incorrect _sc_pointer or _num_cell doesn't match the
 * _neighbors vector), then neither will this Manager.
 * The reactions retreive _sc_ref through the getInstance() method.
 */

Manager::Manager( Manager *newOne ) {
	
	_sc_ref = newOne->_sc_ref;
	_time = newOne->_time;
	_num_mol = newOne->_num_mol;
	_num_cell = newOne->_num_cell;
	
	/* We do not give Gene a public copy constructor. We instead copy 
	 * manually
	 */
	
	for ( int i_gene = 0 ; i_gene < newOne->_genes.size() ; i_gene++ ) {
		Gene *gene_ref = newOne->_genes.at(i_gene);
		_genes.push_back( new Gene ( gene_ref->get_i_self() , gene_ref->get_i_product() , 
									gene_ref->get_i_bound_promoter() , gene_ref->get_i_root() , 
									gene_ref->get_init_conc() ) );
	}
	
	for ( int i_prot = 0 ; i_prot < newOne->_proteins.size() ; i_prot++ ) {
		Protein *prot_ref = newOne->_proteins.at(i_prot);
		_proteins.push_back( new Protein ( prot_ref->get_i_self() , prot_ref->get_i_root_zero() ,
										  prot_ref->get_i_root_one() , prot_ref->get_init_conc() ) );
	}
	
	for ( int i_reac = 0 ; i_reac < newOne->_reactions.size() ; i_reac++ ) {
		_reactions.push_back( newOne->_reactions.at(i_reac)->copy() );
	}
	
}

/* Destructor: ~Manager() 
 * -------------------------------------------------------------------------- 
 */
Manager::~Manager() {
	_genes.erase(_genes.begin(),_genes.end());
	_proteins.erase(_proteins.begin(),_proteins.end());
	_reactions.erase(_reactions.begin(),_reactions.end());
}

/* Constructor: mutate(generator) 
 * --------------------------------------------------------------------------
 * Performs a random mutation on the genome.
 *
 * Though eventually we want to be able to change the parameters of the
 * mutation process, currently mutations occur as follows:
 *
 * One of the following 5 types of mutations can occur, with an equal
 * probability of each:
 *		
 *		1. The degredation constant of a protein is modified.
 *		2. A kinetic constant for a non-degredation reaction is modified.
 *		3. A gene and its associated protein are added to the genome.
 *		4. A interaction between a protein and a gene is added, creating
 *		a gene-protein complex with a modified protein production rate.
 *		5. A post-transcriptional mutation is added (ie, dimerization, etc)
 * 
 */
void Manager::mutate( boost::random::mt19937& generator ) {
	
	boost::random::uniform_int_distribution<> which_mutation(0,4);
	int choice = which_mutation(generator);
	
	switch (choice) {
		case 0:
			degredation_mutation(generator);
			break;
		case 1:
			kinetic_mutation(generator);
			break;
		case 2:
			add_gene(1.0,1.0); // Production and degredation set at 1.0
			break;
		case 3:
			prom_binding_mutation(generator);
			break;
		case 4:
			post_transcript_mutation(generator);
			break;
		default:
			break;
	}
	
}


/* Public Method: initialize()
 * -------------------------------------------------------------------------- 
 * Initialize by making the dmat the appropriate size according to _num_mol,
 * and then filling each dmat with the initial concentrations of the proteins 
 * and genes, and setting _time = 0.
 */
void Manager::initialize() {
	resize_dmats();
	_time = 0.0;
	
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		for ( int i_gene = 0 ; i_gene < _genes.size() ; i_gene++ ) {
			Gene* gene_ref = _genes.at(i_gene);
			_i_tissue.at(i_cell,gene_ref->get_i_self()) = gene_ref->get_init_conc();
			_curr_tissue.at(i_cell,gene_ref->get_i_self()) = gene_ref->get_init_conc();
		}
		for ( int i_protein = 0 ; i_protein < _proteins.size() ; i_protein++ ) {
			Protein* protein_ref = _proteins.at(i_protein);
			_i_tissue.at(i_cell,protein_ref->get_i_self()) = protein_ref->get_init_conc();
			_curr_tissue.at(i_cell,protein_ref->get_i_self()) = protein_ref->get_init_conc();
		}
	}
}

/* Public Method: integrate(num_step,generator)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the currTissue vector based on iTissue, 
 * reactions, and integration mode.
 */
void Manager::integrate( int num_step , boost::random::mt19937& generator ) {
	
	switch ( _sc_ref->_mode ) {
		
		case RK4_DET_TI:
			rk4_det_ti(num_step);
			break;
			
		case RK4_STC_TI:
			rk4_stc_ti( num_step , generator );
			break;
			
		default:
			std::cout << "No integration occurred: Invalid mode\n";
			break;
	}
}

/* Public Method: integrate(num_step,generator,file)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the currTissue vector based on iTissue, 
 * reactions, and integration mode. Writes state to file at each time step.
 */
void Manager::integrate( int num_step , boost::random::mt19937& generator ,
						ofstream& file ) {
	
	switch ( _sc_ref->_mode ) {
			
		case RK4_DET_TI:
			rk4_det_ti(num_step,file);
			break;
			
		case RK4_STC_TI:
			rk4_stc_ti( num_step , generator , file );
			break;
			
		default:
			std::cout << "No integration occurred: Invalid mode\n";
			break;
	}
}

/* Public Method: get_curr_state()
 * -------------------------------------------------------------------------- 
 */
dmat Manager::get_curr_state() {
	return _curr_tissue;
}

/* Public Method: print_state()
 * -------------------------------------------------------------------------- 
 * Prints the concentration of each molecule in each cell.
 *
 * Format:
 *
 *		time: <t>
 *		Cell: <i>
 *		conc_1 :::: conc_2 :::: conc_3 :::: conc_4 :::: ..... :::: conc_n
 *
 */
void Manager::print_state() {
	std::cout << "time: " << _time << std::endl;
	for ( int i = 0 ; i < _num_cell ; i++ ) {
		std::cout << "Cell: " << i << std::endl;
		for ( int j = 0 ; j < (_genes.size() + _proteins.size() - 1) ; j++ ) {
			std::cout << _curr_tissue.at(i,j) << " :::: ";
		}
		std::cout << _curr_tissue.at(i,(_genes.size() + _proteins.size() - 1)) << std::endl;
	}
	std::cout << std::endl;
}

/* Public Method: print_genome()
 * -------------------------------------------------------------------------- 
 * Lists the current genes, proteins, and reactions in our genome.
 * 
 * Format:
 *
 *		GENES:
 *		Gene <i>:
 *			<info>
 *
 *		PROTEINS:
 *		Protein <i>:
 *			<info>
 *
 *		REACTIONS:
 *		Reaction <i>:
 *			<info>
 * 
 */
void Manager::print_genome() {
	
	std::cout << std::endl;
	std::cout << "GENES: " << std::endl;
	for ( int i_gene = 0 ; i_gene < _genes.size() ; i_gene++ ) {
		std::cout << "Gene " << i_gene << ": " << std::endl;
		_genes.at(i_gene)->print_info("\t");
	}
	std::cout << std::endl << std::endl;
	
	std::cout << "PROTEINS: " << std::endl;
	for ( int i_protein = 0 ; i_protein < _proteins.size() ; i_protein++ ) {
		std::cout<< "Protein " << i_protein + _genes.size() << ":\n";
		_proteins.at(i_protein)->print_info("\t");
	}
	std::cout << std::endl << std::endl;
	
	std::cout << "REACTIONS: " << std::endl;
	for ( int i_reac = 0 ; i_reac < _reactions.size() ; i_reac++ ) {
		std::cout << "Reaction " << i_reac << ": " << std::endl;
		_reactions.at(i_reac)->print_info("\t");
	}
	std::cout << std::endl << std::endl;
}

/* Public Method: genome_to_file(file,line_start)
 * -------------------------------------------------------------------------- 
 */

void Manager::genome_to_file( ofstream& file , std::string line_start ) {

	file << "GENES:\n";
	for ( int i_gene = 0 ; i_gene < _genes.size() ; i_gene++ ) {
		file << "Gene " << i_gene << ":\n";
		_genes.at(i_gene)->to_file(file,line_start);
	}
	file << "\n\n";
	
	file << "PROTEINS:\n";
	for ( int i_protein = 0 ; i_protein < _proteins.size() ; i_protein++ ) {
		file << "Protein " << i_protein + _genes.size() << ":\n";
		_proteins.at(i_protein)->to_file(file,line_start);
	}
	file << "\n\n";
	
	file << "REACTIONS:\n";
	for ( int i_reac = 0 ; i_reac < _reactions.size() ; i_reac++ ) {
		file << "Reaction " << i_reac << ":\n";
		_reactions.at(i_reac)->to_file(file,line_start);
	}
	
}

/* Public Method: state_to_file(file)
 * -------------------------------------------------------------------------- 
 */

void Manager::state_to_file( std::ofstream& file ) {
	
	file << _time << " ";
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		for ( int i_mol = 0 ; i_mol < _num_mol ; i_mol++ ) {
			file << _curr_tissue.at(i_cell,i_mol) << " ";
		}
	}
	file << "\n";
	
}

/* Private Method: react(xi,dx_dt)
 * -------------------------------------------------------------------------- 
 * Given the current reactions in our genome, it fills the dx_dt vector, 
 * passed to it by reference, with the current rates of change of the 
 * molecules, given the concentrations held in the dmat xi.
 *
 * CAUTION: We assume xi and dx_dt are the same size. Nothing in the method
 * checks for this.
 */

void Manager::react(dmat& xi,dmat& dx_dt) {
	
	/* Because the react method for reactions only adds a rate to the current
	 * rates held in dx_dt, we must first set all elements to zero.
	 */
	dx_dt.zero();
	
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		for ( int i_reac = 0 ; i_reac < _reactions.size() ; i_reac++ ) {
			_reactions.at(i_reac)->react( xi , dx_dt , i_cell );
		}
	}
	
}

/* Private Method: react(xi,dx_dt,dist,generator)
 * -------------------------------------------------------------------------- 
 * Given the current reactions in our genome, it fills the dx_dt vector, 
 * passed to it by reference, with the current rates of change of the 
 * molecules, given the concentrations held in the dmat xi.
 *
 * CAUTION: We assume xi and dx_dt are the same size. Nothing in the method
 * checks for this.
 */

void Manager::react( dmat& xi , dmat& dx_dt , 
					boost::random::normal_distribution<>& dist ,
					boost::random::mt19937& generator ,
					double q) {
	
	/* Because the react method for reactions only adds a rate to the current
	 * rates held in dx_dt, we must first set all elements to zero.
	 */
	dx_dt.zero();
	
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		for ( int i_reac = 0 ; i_reac < _reactions.size() ; i_reac++ ) {
			_reactions.at(i_reac)->react( xi , dx_dt , i_cell ,
										 dist , generator , q );
		}
	}
	
}

/* Private Method: resize_dmats()
 * -------------------------------------------------------------------------- 
 */
void Manager::resize_dmats() {
	
	_i_tissue.resize( _num_cell , _num_mol );
	_curr_tissue.resize( _num_cell , _num_mol );
	
}

/* Private Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Goes through every gene, protein, and reaction and instructs them to 
 * update their internally stored data regarding their own index in the
 * dvecs and well as the index of other compounds they store (like their
 * products, roots, promoters, etc.)
 *
 * first_index is the index of the first member of the insertions. A negative
 * num_insertion refers to deletions.
 *
 * CAUTION: we assume first_index is within the range of our current number
 * of molecules, and num_insertions, when negative, does not delete more
 * substances than there are to delete. Could create run-time errors
 * otherwise.
 *
 * EXAMPLES: 
 *
 *		update_indices(3,2) produces
 *
 *		  < 1 , 10 , 3 , 7 , 4 , 11 > --> < 1 , 10 , 3 ,___,___, 7 , 4 , 11 >
 *
 *		update_indices(3,-2) produces
 *
 *		   < 1 , 10 , 3 , 7 , 4 , 11 > --> < 1 , 10 , 3 , 11 >
 *
 *		On the vector above, update_indices(3,-4) would produce run-time 
 *		errors.
 */
void Manager::update_indices( int first_index , int num_insertion ) {
	for ( int i = 0 ; i < _genes.size() ; i++ ) {
		_genes.at(i)->update_indices(first_index,num_insertion);
	}
	for ( int i = 0 ; i < _proteins.size() ; i++ ) {
		_proteins.at(i)->update_indices(first_index,num_insertion);
	}
	for ( int i = 0 ; i < _reactions.size() ; i++ ) {
		_reactions.at(i)->update_indices(first_index,num_insertion);
	}
}

/* Private Method: add_gene(prod_rate)
 * -------------------------------------------------------------------------- 
 * Adds a new gene to the genome, as well as the corresonding protein, and
 * sets the promotion rate to prod_rate. Does nothing to the dmats that will
 * be used for integration.
 *
 * All existing substances and reactions remain unchanged.
 *
 * The tasks that need to be performed are:
 *
 *		1. All internal data of current genes, proteins, and reactions, is 
 *		updated to accomodate the insertion of a new gene at the end of our
 *		list of genes, and a new protein at the end of our list of proteins.
 *		2. A new gene is added to the gene vector, with appropriate data.
 *		3. A new protein is added to the protein vector, with appropriate 
 *		data. 
 *		4. A new reaction is added describing the promotion of the protein
 *		by the gene, with appropriate data.
 *		5. A new reaction is added describing the degredation of the protein.
 *		6. The promotion reaction is added to the reaction vectors of the
 *		gene and protein.
 *		7. The degredation reaction is added to the reaction vector of the
 *		protein.
 *		8. Updates _num_mol.
 */
void Manager::add_gene( double prod_rate , double deg_rate ) {
	
	
	int i_new_gene = _genes.size();
	int i_new_prot = _num_mol+1; /* We are using the old, not yet updated, value 
								  * of numMol 
								  */
	int i_new_gene_in_genes = i_new_gene;
	int i_new_prot_in_prots = i_new_prot - (_genes.size() + 1); // -1 b/c we haven't yet added the new gene to genes
	
	/* Step 1 */
	update_indices(i_new_gene,1);
	
	/* Steps 2 - 5 */
	_genes.push_back(new Gene(i_new_gene,i_new_prot,NEXIST,NEXIST,1.0));
	_proteins.push_back(new Protein(i_new_prot,NEXIST,NEXIST,0.0));
	_reactions.push_back(new PromReaction(i_new_gene,i_new_prot,prod_rate));
	_reactions.push_back(new DegReaction(i_new_prot,deg_rate));
	
	/* Step 6 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-2); // add new PromReaction
	_genes.at(i_new_gene_in_genes)->add_reaction(curr_reac_ref); // to new gene
	_proteins.at(i_new_prot_in_prots)->add_reaction(curr_reac_ref); // to new product protein
	
	Reaction* success;
	success = _genes.at(i_new_gene_in_genes)->get_reaction(0);
	
	/* Step 7 */
	curr_reac_ref = _reactions.at(_reactions.size()-1); // add new DegReaction
	_proteins.at(i_new_prot_in_prots)->add_reaction(curr_reac_ref); // to now product protein
	
	/* Step 8 */
	_num_mol += 2;
	
	std::cout << "GENE ADDITION OCCURRED" << std::endl;
	
}

/* Private Method: addPromBindingReaction(iGeneInGenes,iProtInProts,
 *										 forwardKinetic,backwardKinetic,
 *										 newProdRate)
 * -------------------------------------------------------------------------- 
 * Makes the protein in position iProt in proteins vector promote the gene 
 * in position iGene in the genes vector. A new gene:prot complex is formed, 
 * which promotes the protein gene produces at newProdRate.
 *
 * CAUTION: It is assumed that iGeneInGenes is the index of a gene and 
 * iProtInProts is the index of a protein. No check is performed by the 
 * function to guarantee this.
 *
 * The tasks that need to be performed are:
 *
 *		1. All internal data of current genes, proteins, and reactions, is 
 *		updated to accomodate the insertion of a new gene at the end of our
 *		list of genes. (The gene:prot complex)
 *		2. A new gene is added to the gene vector, with appropriate data.
 *		2. A new reaction is added describing the binding and unbinding of
 *		gene and protein. (PromBindingReaction)
 *		3. A reaction is added describing the production of proteins by the
 *		gene:prot complex. (PromReaction)
 *		5. The PromBindingReaction is added to the internal data of the gene,
 *		binding protein, and gene-protein complex.
 *		6. The PromReaction is added to the interal data of the gene-protein
 *		complex and the promoted gene.
 *		7. Updates numMol.
 */
void Manager::add_PromBindingReac( int i_gene_in_genes , int i_prot_in_prots ,
								  double forward_kinetic , double backward_kinetic , 
								  double new_prod_rate ) {
	
	int i_promoted_gene = _genes.size();
	
	/* Step 1 */
	update_indices(i_promoted_gene,1);
	
	/* We find these indices after step 1 to assure they are consistent with
	 * the new indexing.
	 */
	int i_bound_protein = _proteins.at(i_prot_in_prots)->get_i_self();
	int i_root_gene = _genes.at(i_gene_in_genes)->get_i_self();
	int i_product_protein = _genes.at(i_gene_in_genes)->get_i_product();
	
	/* Steps 2 */
	_genes.push_back(new Gene(i_promoted_gene,i_product_protein,i_bound_protein,
							  i_root_gene,0.0));
	
	/* Steps 3 - 4 */
	_reactions.push_back(new PromBindingReaction(i_root_gene,i_promoted_gene,i_bound_protein,
												 forward_kinetic,backward_kinetic));
	_reactions.push_back(new PromReaction(i_promoted_gene,i_product_protein,new_prod_rate));
	
	/* Step 5 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-2); // add new PromBindingReaction
	_genes.at(i_gene_in_genes)->add_reaction(curr_reac_ref); // to root gene
	_proteins.at(i_prot_in_prots)->add_reaction(curr_reac_ref); // to bound protein
	_genes.at(_genes.size()-1)->add_reaction(curr_reac_ref); // to new gene:prot complex
	
	/* Step 6 */
	int i_product_prot_in_prots = i_product_protein - _genes.size();
	curr_reac_ref = _reactions.at(_reactions.size()-1); // add new PromReaction
	_proteins.at(i_product_prot_in_prots)->add_reaction(curr_reac_ref); // to product protein
	_genes.at(_genes.size()-1)->add_reaction(curr_reac_ref); // to new gene:prot complex
	
	/* Step 7 */
	_num_mol += 1;
	
}

/* Private Method: add_CombReaction(i_prot_zero_in_prots,_i_prot_one_in_prots,
 *								  forward_kinetic,backward_kinetic)
 * -------------------------------------------------------------------------- 
 * Makes the proteins in positions i_prot_zero_in_prots and 
 * i_prot_one_in_prots (which can be the same) combine to form a new protein 
 * complex. With the forward and backward rates as given.
 *
 * CAUTION: It is assumed that i_prot_zero_in_prots and i_prot_one_in_prots 
 * are the indices of proteins. No check is performed by the function to 
 * guarantee this.
 *
 * The tasks that need to be performed are:
 *
 *		1. All internal data of current genes, proteins, and reactions, is 
 *		updated to accomodate the insertion of a new protein at the end of 
 *		our list of genes. (The prot:prot complex)
 *		2. A new protein is added to the protein vector with the appropriate
 *		data.
 *		3. A new reaction is added describing the binding and unbinding of
 *		the proteins. (CombReaction)
 *		4. The CombReaction is added to the internal list of reactions of the 
 *		individual proteins and the new prot:prot complex.
 *		5. Updates numMol.
 */

void Manager::add_CombReaction( int i_prot_zero_in_prots , int i_prot_one_in_prots , 
							   double forward_kinetic , double backward_kinetic ) {
	
	int i_product = _num_mol; // The yet to be updated value of numMol
	
	/* Step 1 */
	/* No updates are necessary because the protein is added to the end of our
	 * list of molecules.
	 */
	
	/* Step 2 */
	int i_reac_zero = i_prot_zero_in_prots + _genes.size();
	int i_reac_one = i_prot_zero_in_prots + _genes.size();
	_proteins.push_back(new Protein(i_product,i_reac_zero,i_reac_one,0.0));
	
	/* Step 3 */
	_reactions.push_back(new CombReaction(i_reac_zero,i_reac_one,i_product,
										  forward_kinetic,backward_kinetic));
	
	/* Step 4 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-1); // add new CombReaction
	_proteins.at(i_prot_zero_in_prots)->add_reaction(curr_reac_ref); // to first participant
	if (i_prot_zero_in_prots != i_prot_one_in_prots) {
		_proteins.at(i_prot_one_in_prots)->add_reaction(curr_reac_ref); // to second participant
	}
	_proteins.at(_proteins.size()-1)->add_reaction(curr_reac_ref); // to their product.
	
	/* Step 5 */
	_num_mol += 1;
}

/* Private Method: add_LatPromReac(i_loc_prot_in_prots,
 *								  i_neighbor_prot_in_prots,
 *								  double kinetic,double K)
 * -------------------------------------------------------------------------- 
 * Makes the protein in position i_loc_prot_in_prots in the protein vector
 * get promoted by the presence of the protein i_neighbor_prot_in_prots in
 * neighboring cells, with the provided kinetic constants.
 * 
 * CAUTION: It is assumed that i_loc_prot_in_prots and 
 * i_neighbor_prot_in_prots is the index of a protein. No check is performed 
 * by the function to guarantee this.
 *
 * The tasks that need to be performed are:
 *
 *		1. A new reaction is added describing the lateral promotion reaction.
 *		2. This reaction is added to the internal data of the participating
 *		proteins.
 */

void Manager::add_LatPromReac( int i_promoted_by_neighbors_in_prots , 
							  int i_promoting_neighbors_in_prots ,
							  double kinetic , double K ) {
	
	int i_promoting_neighbors = _genes.size() + i_promoting_neighbors_in_prots;
	int i_promoted_by_neighbors = _genes.size() + i_promoted_by_neighbors_in_prots;
	
	/* Step 1 */
	_reactions.push_back(new LatPromReaction(i_promoted_by_neighbors,i_promoting_neighbors,
											 kinetic,K));
	
	/* Step 2 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-1); // add new LatPromReaction
	_proteins.at(i_promoted_by_neighbors_in_prots)->add_reaction(curr_reac_ref); // to local protein
	if ( i_promoted_by_neighbors_in_prots != i_promoting_neighbors_in_prots ) {
		_proteins.at(i_promoting_neighbors_in_prots)->add_reaction(curr_reac_ref); // to neighbor protein
	}
	
}

/* Private Method: add_HillPromReac(i_promoter,i_promoted,
 *								  kinetic,K,cooperativity)
 * -------------------------------------------------------------------------- 
 * Makes the protein in position i_promoted_in_prots in the protein vector
 * get promoted by the presence of the protein i_promoter_in_prots in
 * the same cell, with the provided kinetic constants.
 * 
 * CAUTION: It is assumed that i_promoter_in_prots and 
 * i_promoted_in_prots is the index of a protein. No check is performed 
 * by the function to guarantee this.
 *
 * The tasks that need to be performed are:
 *
 *		1. A new reaction is added describing the hill promotion reaction.
 *		2. This reaction is added to the internal data of the participating
 *		proteins.
 */

void Manager::add_HillPromReac( int i_promoter_in_prots , int i_promoted_in_prots ,
							   double kinetic , double K , double cooperativity ) {
	
	int i_promoter = _genes.size() + i_promoter_in_prots;
	int i_promoted = _genes.size() + i_promoted_in_prots;
	
	/* Step 1 */
	_reactions.push_back( new HillPromReaction(i_promoter,i_promoted,
											   kinetic,K,cooperativity) );
	
	/* Step 2 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-1); // add new HillPromReaction
	_proteins.at(i_promoter_in_prots)->add_reaction(curr_reac_ref); // to promoter protein
	if ( i_promoter_in_prots != i_promoted_in_prots ) {
		_proteins.at(i_promoted_in_prots)->add_reaction(curr_reac_ref); // to promoted protein
	}
	
}

/* Private Method: add_HillRepReac(i_promoter,i_promoted,
 *								  kinetic,K,cooperativity)
 * -------------------------------------------------------------------------- 
 * Makes the protein in position i_repressed_in_prots in the protein vector
 * get promoted by the presence of the protein i_repressor_in_prots in
 * the same cell, with the provided kinetic constants.
 * 
 * CAUTION: It is assumed that i_repressor_in_prots and 
 * i_repressed_in_prots is the index of a protein. No check is performed 
 * by the function to guarantee this.
 *
 * The tasks that need to be performed are:
 *
 *		1. A new reaction is added describing the hill repression reaction.
 *		2. This reaction is added to the internal data of the participating
 *		proteins.
 */

void Manager::add_HillRepReac( int i_repressor_in_prots , int i_repressed_in_prots ,
							   double kinetic , double K , double cooperativity ) {
	
	int i_repressor = _genes.size() + i_repressor_in_prots;
	int i_repressed = _genes.size() + i_repressed_in_prots;
	
	/* Step 1 */
	_reactions.push_back(new HillRepReaction(i_repressor,i_repressed,
											  kinetic,K,cooperativity));
	
	/* Step 2 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-1); // add new HillPromReaction
	_proteins.at(i_repressor_in_prots)->add_reaction(curr_reac_ref); // to promoter protein
	if ( i_repressor_in_prots != i_repressed_in_prots ) {
		_proteins.at(i_repressed_in_prots)->add_reaction(curr_reac_ref); // to promoted protein
	}	
}


/* Private Method: degredation_mutation(generator)
 * -------------------------------------------------------------------------- 
 */

void Manager::degredation_mutation(  boost::random::mt19937& generator ) {
	
	/* Vector to tell us indices of degradation reactions */
	std::vector<int> i_deg;
	
	/* Collect indices of degradation reactions */
	for ( int i = 0 ; i < _reactions.size() ; i++ ) {
		if ( _reactions.at(i)->get_type() == DEGRADATION ) {
			i_deg.push_back(i);
		}
	}
	
	/* Randomly pick degradation to mutate */
	if ( i_deg.size() > 0 ) {
		boost::random::uniform_int_distribution<> which_reaction(0,i_deg.size()-1);
		int i_reac = which_reaction(generator);
		_reactions.at( i_deg.at(i_reac) )->mutate(generator);
		std::cout << "DEGREDATION MUTATION OCCURRED (Reaction " << i_deg.at(i_reac) << ")" << std::endl;
	}
	
	
	
}

/* Private Method: kinetic_mutation(generator)
 * --------------------------------------------------------------------------
 * Mutates a kinetic constant of a non-degradation reaction.
 */

void Manager::kinetic_mutation( boost::random::mt19937& generator ) {
	
	/* Vector to tell us indices of non-degradation reactions */
	std::vector<int> i_ndeg;
	
	/* Collect indices of degradation reactions */
	for ( int i = 0 ; i < _reactions.size() ; i++ ) {
		if ( _reactions.at(i)->get_type() != DEGRADATION ) {
			i_ndeg.push_back(i);
		}
	}
	
	/* Randomly pick degradation to mutate */
	if ( i_ndeg.size() > 0 ) {
		boost::random::uniform_int_distribution<> which_reaction(0,i_ndeg.size()-1);
		int i_reac = which_reaction(generator);
		_reactions.at( i_ndeg.at(i_reac) )->mutate(generator);
		std::cout << "MUTATION OF KINETIC CONSTANT OCCURRED (Reaction " << i_ndeg.at(i_reac) << ")" << std::endl;
	}
	
}

/* Private Method: prom_binding_mutation(generator)
 * -------------------------------------------------------------------------- 
 */

void Manager::prom_binding_mutation( boost::random::mt19937& generator ) {
	
	if ( _genes.size() > 0 && _proteins.size() > 0 ) {
		
		boost::random::uniform_int_distribution<> which_gene(0,_genes.size()-1);
		boost::random::uniform_int_distribution<> which_protein(0,_proteins.size()-1);
		boost::random::uniform_real_distribution<> binding_kinetic(0.0,10.0);
		boost::random::uniform_real_distribution<> prod_rate(0.0,2.0);
	
		/* Add the new reaction */
		int i_gene = which_gene(generator);
		int i_prot = which_protein(generator);
		double forward_kinetic = binding_kinetic(generator);
		double backward_kinetic = binding_kinetic(generator);
		double production_rate = prod_rate(generator);
		add_PromBindingReac( i_gene , i_prot ,
							forward_kinetic , backward_kinetic ,
							production_rate );
	
		std::cout << "PROMOTER BINDING MUTATION OCCURRED" << std::endl;
		
	}
										
}

/* Private Method: post_transcript_mutation(generator)
 * -------------------------------------------------------------------------- 
 */

void Manager::post_transcript_mutation( boost::random::mt19937& generator ) {
	
	if ( _proteins.size() > 0 ) {
		boost::random::uniform_int_distribution<> reaction_type(0,1);
		boost::random::uniform_int_distribution<> which_protein(0,_proteins.size()-1);
		boost::random::uniform_real_distribution<> kinetic(0.0,5.0);
	
		int reac_type = reaction_type(generator);
		int i_prot_zero = which_protein(generator);
		int i_prot_one = which_protein(generator);
		double kinetic_zero = kinetic(generator);
		double kinetic_one = kinetic(generator);
	
		switch (reac_type) {
			case 0:
				add_CombReaction( i_prot_zero , i_prot_one ,
								 kinetic_zero , kinetic_one );
				break;
			case 1:
				add_LatPromReac( i_prot_zero , i_prot_one ,
								kinetic_zero , kinetic_one );
				break;
			default:
				break;
		}
		
		std::cout << "POST TRANSCRIPTION MUTATION OCCURRED" << std::endl;
		
	}
	
}

/* Private Method: rk4_det_ti(num_step)
 * -------------------------------------------------------------------------- 
 */

void Manager::rk4_det_ti( int num_step ) {

	dmat k1 ( _num_cell , _num_mol );
	dmat k2 ( _num_cell , _num_mol );
	dmat k3 ( _num_cell , _num_mol );
	dmat k4 ( _num_cell , _num_mol );
	
	dmat x2 ( _num_cell , _num_mol );
	dmat x3 ( _num_cell , _num_mol );
	dmat x4 ( _num_cell , _num_mol );
	
	dmat dx_dt ( _num_cell , _num_mol );
	
	double dt = _sc_ref->_dt;
	
	for ( int step = 0 ; step < num_step ; step++ ) {
		
		/* ----------------------------------------- */
		react ( _curr_tissue , dx_dt );
		k1.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x2.equals_lin_comb(_curr_tissue , 1.0 ,
						   k1 , .5 );
		react ( x2 , dx_dt );
		k2.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x3.equals_lin_comb(_curr_tissue , 1.0 ,
						   k2 , .5 );
		react ( x3 , dx_dt );
		k3.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x4.equals_lin_comb(_curr_tissue,1.0 ,
						   k3,1.0 );
		react ( x4 , dx_dt );
		k4.equals_lin_comb( dx_dt , dt);
		/* ----------------------------------------- */
		_curr_tissue.plus_equals_lin_comb(k1,.166666666667,
										  k2,.333333333333,
										  k3,.333333333333,
										  k4,.166666666667);
		/* ----------------------------------------- */
		_time += dt;
		//print_state();
	}
}

/* Private Method: rk4_det_ti(num_step,file)
 * -------------------------------------------------------------------------- 
 */

void Manager::rk4_det_ti( int num_step , ofstream& file) {
	
	dmat k1 ( _num_cell , _num_mol );
	dmat k2 ( _num_cell , _num_mol );
	dmat k3 ( _num_cell , _num_mol );
	dmat k4 ( _num_cell , _num_mol );
	
	dmat x2 ( _num_cell , _num_mol );
	dmat x3 ( _num_cell , _num_mol );
	dmat x4 ( _num_cell , _num_mol );
	
	dmat dx_dt ( _num_cell , _num_mol );
	
	double dt = _sc_ref->_dt;
	
	for ( int step = 0 ; step < num_step ; step++ ) {
		
		/* ----------------------------------------- */
		react ( _curr_tissue , dx_dt );
		k1.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x2.equals_lin_comb(_curr_tissue , 1.0 ,
						   k1 , .5 );
		react ( x2 , dx_dt );
		k2.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x3.equals_lin_comb(_curr_tissue , 1.0 ,
						   k2 , .5 );
		react ( x3 , dx_dt );
		k3.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x4.equals_lin_comb(_curr_tissue,1.0 ,
						   k3,1.0 );
		react ( x4 , dx_dt );
		k4.equals_lin_comb( dx_dt , dt);
		/* ----------------------------------------- */
		_curr_tissue.plus_equals_lin_comb(k1,.166666666667,
										  k2,.333333333333,
										  k3,.333333333333,
										  k4,.166666666667);
		/* ----------------------------------------- */
		_time += dt;
		
		state_to_file(file);
		
	}
}

/* Private Method: rk4_stc_ti(num_step,generator)
 * -------------------------------------------------------------------------- 
 */

void Manager::rk4_stc_ti( int num_step , boost::random::mt19937& generator ) {
	
	boost::random::normal_distribution<> norm_dist(0.0,_sc_ref->_stochasticity);
	
	double a21;
	double a31;
	double a32;
	double a41;
	double a42;
	double a43;
	double a51;
	double a52;
	double a53;
	double a54;
	dmat k1 ( _num_cell , _num_mol );
	dmat k2 ( _num_cell , _num_mol );
	dmat k3 ( _num_cell , _num_mol );
	dmat k4 ( _num_cell , _num_mol );
	double q1;
	double q2;
	double q3;
	double q4;
	dmat x2 ( _num_cell , _num_mol );
	dmat x3 ( _num_cell , _num_mol );
	dmat x4 ( _num_cell , _num_mol );
	dmat dx_dt ( _num_cell , _num_mol );
	
	a21 =   2.71644396264860;
	a31 = - 6.95653259006152;
	a32 =   0.78313689457981;
	a41 =   0.0;
	a42 =   0.48257353309214;
	a43 =   0.26171080165848;
	a51 =   0.47012396888046;
	a52 =   0.36597075368373;
	a53 =   0.08906615686702;
	a54 =   0.07483912056879;
	
	q1 =   2.12709852335625;
	q2 =   2.73245878238737;
	q3 =  11.22760917474960;
	q4 =  13.36199560336697;
	
	double dt = _sc_ref->_dt;
	
	for ( int step = 0 ; step < num_step ; step++ ) {
		 
		react(_curr_tissue,dx_dt,norm_dist,generator,q1);
		k1.equals_lin_comb(dx_dt,dt);
		/* ----------------------------------------- */
		x2.equals_lin_comb(_curr_tissue,1.0,
						   k1,a21);
		react(x2,dx_dt,norm_dist,generator,q2);
		k2.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x3.equals_lin_comb(_curr_tissue,1.0,
						   k1,a31,
						   k2,a32);
		react(x3,dx_dt,norm_dist,generator,q3);
		k3.equals_lin_comb(dx_dt,dt);
		/* ----------------------------------------- */
		x4.equals_lin_comb(_curr_tissue,1.0,
						   k1,a41,
						   k2,a42);
		react(x4,dx_dt,norm_dist,generator,q4);
		k4.equals_lin_comb(dx_dt,dt);
		/* ----------------------------------------- */
		_curr_tissue.plus_equals_lin_comb(k1,a51,
										  k2,a52,
										  k3,a53,
										  k4,a54);
		
		_time += dt;
		
		
	}
}

/* Private Method: rk4_stc_ti(num_step,generator,file)
 * -------------------------------------------------------------------------- 
 */

void Manager::rk4_stc_ti( int num_step , boost::random::mt19937& generator , std::ofstream& file) {
	
	boost::random::normal_distribution<> norm_dist(0.0,_sc_ref->_stochasticity);
	
	double a21;
	double a31;
	double a32;
	double a41;
	double a42;
	double a43;
	double a51;
	double a52;
	double a53;
	double a54;
	dmat k1 ( _num_cell , _num_mol );
	dmat k2 ( _num_cell , _num_mol );
	dmat k3 ( _num_cell , _num_mol );
	dmat k4 ( _num_cell , _num_mol );
	double q1;
	double q2;
	double q3;
	double q4;
	dmat x2 ( _num_cell , _num_mol );
	dmat x3 ( _num_cell , _num_mol );
	dmat x4 ( _num_cell , _num_mol );
	dmat dx_dt ( _num_cell , _num_mol );
	
	a21 =   2.71644396264860;
	a31 = - 6.95653259006152;
	a32 =   0.78313689457981;
	a41 =   0.0;
	a42 =   0.48257353309214;
	a43 =   0.26171080165848;
	a51 =   0.47012396888046;
	a52 =   0.36597075368373;
	a53 =   0.08906615686702;
	a54 =   0.07483912056879;
	
	q1 =   2.12709852335625;
	q2 =   2.73245878238737;
	q3 =  11.22760917474960;
	q4 =  13.36199560336697;
	
	double dt = _sc_ref->_dt;
	
	for ( int step = 0 ; step < num_step ; step++ ) {
		
		react(_curr_tissue,dx_dt,norm_dist,generator,q1);
		k1.equals_lin_comb(dx_dt,dt);
		/* ----------------------------------------- */
		x2.equals_lin_comb(_curr_tissue,1.0,
						   k1,a21);
		react(x2,dx_dt,norm_dist,generator,q2);
		k2.equals_lin_comb( dx_dt , dt );
		/* ----------------------------------------- */
		x3.equals_lin_comb(_curr_tissue,1.0,
						   k1,a31,
						   k2,a32);
		react(x3,dx_dt,norm_dist,generator,q3);
		k3.equals_lin_comb(dx_dt,dt);
		/* ----------------------------------------- */
		x4.equals_lin_comb(_curr_tissue,1.0,
						   k1,a41,
						   k2,a42);
		react(x4,dx_dt,norm_dist,generator,q4);
		k4.equals_lin_comb(dx_dt,dt);
		/* ----------------------------------------- */
		_curr_tissue.plus_equals_lin_comb(k1,a51,
										  k2,a52,
										  k3,a53,
										  k4,a54);
		
		_time += dt;
		state_to_file(file);
	}
}

void Manager::hakim_delta_notch_construct() {
	_num_mol = 0; /* Before we add anything to the genome. Otherwise has random
				   * initiation.
				   */
	
	/* Initialize gene and protein vectors */
	/* Reference detailed constructors and array of
	 * elements of genome above to understand 
	 * paramaters.
	 */
	add_gene(1.0,1.0); // Delta
	/* < d D > */
	add_gene(0.0,1.0); // Notch
	/* < d n D N > */
	add_CombReaction( 1 , 1 , 50.0 , 5.0 ); // N + N <--> N:N
	/* < d n D N N:N > */
	add_PromBindingReac(0 , 2 , 10.0 , 1.0 , 0.0); // d + N:N <--> d:(N:N)
	/* < d n d:(N:N) D N N:N > */
	add_LatPromReac(1 , 0 , 1.0 , .01 ); // N promoted by neighboring D
	
	/* Initialize */
	initialize();
	
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		if ( i_cell % 2 == 0 ) {
			_curr_tissue.at(i_cell,3) = 1.0;
			_curr_tissue.at(i_cell,4) = 0.0;
		}
		else { // i_cell % 2 == 1
			_curr_tissue.at(i_cell,3) = 0.0;
			_curr_tissue.at(i_cell,4) = 1.0;
			_curr_tissue.at(i_cell,5) = 0.0;
		}
	}
}

void Manager::collier_delta_notch_construct() {
	_num_mol = 0;
	
	add_gene(0.0,1.0); // Delta
	/* < d D > */
	add_gene(0.0,1.0); // Notch
	/* < d n D N > */
	add_LatPromReac(1,0,1.0,.01); // dN/dt = (D_n)^2 / (.01 + (D_n)^2)
	/* < d n D N > */
	add_HillRepReac(1,0,1.0,.01,2); // dD/dt = 1 / (1 + 100N^2)
	/* < d n D N > */
	
	/* Initialize protein concentrations */
	initialize();
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	boost::random::normal_distribution<> norm_dist(0.0,1.0);
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		if ( i_cell % 2 == 0 ) {
			_curr_tissue.at(i_cell,2) = .5 + 0.0 * norm_dist(generator);
			_curr_tissue.at(i_cell,3) = .5 + 0.0 * norm_dist(generator);
		}
		else {
			_curr_tissue.at(i_cell,2) = .5 + 0.0 * norm_dist(generator);
			_curr_tissue.at(i_cell,3) = .5 + 0.0 * norm_dist(generator);
		}
	}
	
}

void Manager::mutation_construct() {
	_num_mol = 0; /* Before we add anything to the genome. Otherwise has random
				   * initiation.
				   */
	
	/* Build genome through random mutations */
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	
	add_gene(1.0,1.0);
	
	for (int i = 0 ; i < 5 ; i++) {
		//std::cout << std::endl << "Mutation " << i << std::endl;
		mutate(generator);
		//print_genome();
	}	
	
	initialize();
}

void Manager::one_protein_construct() {
	_num_mol = 0;
	
	add_gene(1.0,1.0);
	
	initialize();
}

