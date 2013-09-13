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
 * Constructor that allows different construction methods, indicated by 
 * ConstructionMethod variable m. ConstructionMethod is enumerated in the 
 * file ConstructionMethod.h.
 *
 * See individual functions for documentation.
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
		else if (m == TWO_PROTEIN) {
			two_protein_construct();
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
	
	/* Pointer to SettingsCont */
	_sc_ref = newOne->_sc_ref;
	
	/* Gene and protein info */
	for ( int i_gene = 0 ; i_gene < newOne->_genes.size() ; i_gene++ ) {
		Gene *gene_ref = newOne->_genes.at(i_gene);
		_genes.push_back(new Gene(*gene_ref));
	}
	for ( int i_prot = 0 ; i_prot < newOne->_proteins.size() ; i_prot++ ) {
		Protein *prot_ref = newOne->_proteins.at(i_prot);
		_proteins.push_back( new Protein (*prot_ref) );
	}
	
	/* Reaction info */
	for ( int i_reac = 0 ; i_reac < newOne->_reactions.size() ; i_reac++ ) {
		_reactions.push_back( newOne->_reactions.at(i_reac)->copy() );
	}
	
	/* Current Time */
	/* Time is not copied. It is also a state variable, so its intiatializing 
	 * is handled by intialize(), called by the user.
	 */
	
	/* Number of molecules */
	_num_mol = newOne->_num_mol;
	
	/* Number of cells */
	_num_cell = newOne->_num_cell;
	
	
	/* State matrix */
	/* Only the size of this is set as it is its only characteristic that
	 * does not depend on where we are in the integration. The values it
	 * hold is handled by initialize().
	 */
	_curr_tissue.resize(_num_cell, _num_mol);
	
	/* Lineage Info */
	for ( int i = 0 ; i < newOne->_evaluation_generations.size() ; i++ ) {
		_evaluation_generations.push_back(newOne->_evaluation_generations.at(i));
	}
	for ( int i = 0 ; i < newOne->_evaluation_scores.size() ; i++ ) {
		_evaluation_scores.push_back(newOne->_evaluation_scores.at(i));
	}
	for ( int i = 0 ; i < newOne->_mutation_generations.size() ; i++ ) {
		_mutation_generations.push_back(newOne->_mutation_generations.at(i));
	}
}

/* Constructor: Manager(file)
 * --------------------------------------------------------------------------
 * Initialize a Manager with the genetic network specified in the file.
 * It is assumed that file is open and at the beginning of a file that is in
 * the standardized format, matching the format created by the 
 * genome_to_file() method. This format should maintain consistancy.
 *
 * See Documentation.txt for details on the format of the files, and how to
 * easily change some aspects of the format in a way that will maintain 
 * consistancy between file input and output.
 */

Manager::Manager( std::ifstream& file ) {
	
	_sc_ref = SettingsCont::getInstance();
	_time = 0.0;
	_num_cell = _sc_ref->_neighbors.size();
	
	char temp[256];
	file.getline(temp,256); // GENES:
	while ( strncmp("PROTEINS:",temp,9) != 0 ) {
		if ( strncmp(temp, "Gene",4) == 0 ) {
			_genes.push_back( new Gene(file) );
		}
		file.getline(temp,256);
	}
	
	while ( strncmp(temp,"REACTIONS:",9) != 0 ) {
		if ( strncmp(temp,"Protein",7) == 0 ) {
			_proteins.push_back( new Protein(file) );
		}
		file.getline(temp,256);
	}
	
	while (!file.eof()) {
		if ( strncmp(temp,"Reaction",8) == 0 ) {
			file.ignore(256,':');
			file.getline(temp,256);
			if ( strncmp(temp," Combination Reaction",21) == 0 ) {
				_reactions.push_back( new CombReaction(file) );
			}
			else if ( strncmp(temp," Degradation Reaction",21) == 0 ) {
				_reactions.push_back( new DegReaction(file) );
			}
			else if ( strncmp(temp," Promotion Reaction",19) == 0 ) {
				_reactions.push_back( new PromReaction(file) );
			}
			else if ( strncmp(temp," Lateral Promotion Reaction",27) == 0 ) {
				_reactions.push_back( new LatPromReaction(file) );
			}
			else if ( strncmp(temp," Promoter Binding Reaction",26) == 0 ) {
				_reactions.push_back( new PromBindingReaction(file) );
			}
			else if ( strncmp(temp," Hill Promotion Reaction",24) == 0 ) {
				_reactions.push_back( new HillPromReaction(file) );
			}
			else if ( strncmp(temp," Hill Repression Reaction",25) == 0 ) {
				_reactions.push_back( new HillRepReaction(file) );
			}
			else if (strncmp(temp," Lateral Repression Reaction",28) == 0 ) {
				_reactions.push_back( new LatRepReaction(file) );
			}
		}
		file.getline(temp,256);
	}
	
	_num_mol = _genes.size() + _proteins.size();
	
}

/* Destructor: ~Manager() 
 * -------------------------------------------------------------------------- 
 */
Manager::~Manager() {
	int num_gene = _genes.size();
	int num_prot = _proteins.size();
	int num_reac = _reactions.size();
	for ( int i_gene = 0 ; i_gene < num_gene ; i_gene++ ) {
		delete _genes.at(i_gene);
	}
	for ( int i_prot = 0 ; i_prot < num_prot ; i_prot++ ) {
		delete _proteins.at(i_prot);
	}
	for ( int i_reac = 0 ; i_reac < num_reac ; i_reac++ ) {
		delete _reactions.at(i_reac);
	}
	_num_mol = -100;
}

/* Constructor: mutate(generator) 
 * --------------------------------------------------------------------------
 * Performs a random mutation on the genome, picked randomly at uniform from
 * the possibilities listed in the _mutation_types in SettingsCont.
 *
 * See Documentation.txt or the individual methods for details on the
 * mutation types.
 */
void Manager::mutate( boost::random::mt19937& generator ) {
	int num_mutation = _sc_ref->_mutation_types.size();
	boost::random::uniform_int_distribution<> which_mutation(0,num_mutation-1);
	int choice = which_mutation(generator);
	int type = _sc_ref->_mutation_types.at(choice);
	switch (type) {
		case DEGRADATION_M:
			degredation_mutation(generator);
			break;
		case KINETIC:
			kinetic_mutation(generator);
			break;
		case ADD_GENE:
			add_gene(1.0,1.0);
			break;
		case PROM_BINDING:
			prom_binding_mutation(generator);
			break;
		case POST_TRANSCRIPT:
			post_transcript_mutation(generator);
			break;
		case ADD_HILL:
			add_hill_mutation(generator);
			break;
		case REMOVAL:
			removal_mutation(generator);
		default:
			break;
	}
	
}

/* Constructor: mutate(generator,generation) 
 * --------------------------------------------------------------------------
 * Exact copy of mutate method which just takes a generator as the argument,
 * but useful for lineage data tracking. Adds the mutation generation to the
 * _mutations_generations vector keeping track of in which generations 
 * mutations have occurred in the current lineage.
 */
void Manager::mutate( boost::random::mt19937& generator , int generation ) {
	int num_mutation = _sc_ref->_mutation_types.size();
	boost::random::uniform_int_distribution<> which_mutation(0,num_mutation-1);
	int choice = which_mutation(generator);
	int type = _sc_ref->_mutation_types.at(choice);
	switch (type) {
		case DEGRADATION_M:
			degredation_mutation(generator);
			break;
		case KINETIC:
			kinetic_mutation(generator);
			break;
		case ADD_GENE:
			add_gene(1.0,1.0);
			break;
		case PROM_BINDING:
			prom_binding_mutation(generator);
			break;
		case POST_TRANSCRIPT:
			post_transcript_mutation(generator);
			break;
		case ADD_HILL:
			add_hill_mutation(generator);
			break;
		case REMOVAL:
			removal_mutation(generator);
		default:
			break;
	}
	_mutation_generations.push_back(generation);
}

/* Public Method: initialize()
 * -------------------------------------------------------------------------- 
 * Initialize by making the dmat the appropriate size according to _num_mol,
 * and then filling the dmat with the initial concentrations of the proteins 
 * and genes, and setting _time = 0.
 */
void Manager::initialize() {
	resize_dmats();
	_time = 0.0;
	
	int genes_size = _genes.size();
	
	for ( int i_cell = 0 ; i_cell < _num_cell ; i_cell++ ) {
		for ( int i_gene = 0 ; i_gene < genes_size ; i_gene++ ) {
			Gene* gene_ref = _genes.at(i_gene);
			_curr_tissue.at(i_cell,gene_ref->get_i_self()) = gene_ref->get_init_conc();
		}
		for ( int i_protein = 0 ; i_protein < _proteins.size() ; i_protein++ ) {
			Protein* protein_ref = _proteins.at(i_protein);
			_curr_tissue.at(i_cell,protein_ref->get_i_self()) = protein_ref->get_init_conc();
		}
	}
}

/* Public Method: integrate(num_step,generator)
 * -------------------------------------------------------------------------- 
 * Runs the ODE for num_step more time steps from the current time, according
 * to curr_tissue, reactions, and integration mode.
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
 * Runs the ODE for num_step more time steps from the current time, according
 * to curr_tissue, reactions, and integration mode. Writes state to file at 
 * each time step.
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

/* Public Method: get_curr_state(curr_state)
 * -------------------------------------------------------------------------- 
 * We copy _curr_state into the curr_state parameter passed by reference. 
 * We pass by reference rather than return a dmat to avoid shallow copying 
 * due to copy elision or other unwanted procedures used by the compiler to 
 * return a dmat object without using copy constructor or assignment 
 * operator.
 */
void Manager::get_curr_state(dmat& curr_state) {
	int num_row = _curr_tissue.get_num_row();
	int num_col = _curr_tissue.get_num_col();
	curr_state.resize( num_row , num_col );
	for ( int row = 0 ; row < num_row ; row++ ) {
		for ( int col = 0 ; col < num_col ; col++ ) {
			curr_state.at(row,col) = _curr_tissue.at(row,col);
		}
	}
}

/* Public Method: get_num_cell()
 * -------------------------------------------------------------------------- 
 */
int Manager::get_num_cell() const {
	return _num_cell;
}

/* Public Method: get_num_gene()
 * -------------------------------------------------------------------------- 
 */
int Manager::get_num_gene() const {
	return _genes.size();
}

/* Public Method: get_num_prot()
 * -------------------------------------------------------------------------- 
 */
int Manager::get_num_prot() const {
	return _proteins.size();
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
 * Outputs the genetic network information to a console. It is assumed that 
 * file is open and at the beginning of a file. Output follows the 
 * standardized format, matching the format expected by the Manager(file) 
 * constructor (which constructs a Manager according to the genome in a 
 * file). Any change in this format should be made consistently in both 
 * places.
 *
 * See Documentation.txt for details on the format of the files, and how to
 * easily change some aspects of the format in a way that will maintain 
 * consistancy between file input and output.
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

/* Public Method: genome_to_file(file,line_start)
 * --------------------------------------------------------------------------
 * Outputs the genetic network information to a console. It is assumed that 
 * file is open and at the beginning of a file. Output follows the 
 * standardized format, matching the format expected by the Manager(file) 
 * constructor (which constructs a Manager according to the genome in a 
 * file). Any change in this format should be made consistently in both 
 * places.
 *
 * See Documentation.txt for details on the format of the files, and how to
 * easily change some aspects of the format in a way that will maintain 
 * consistancy between file input and output.
 */

void Manager::genome_to_file( ofstream& file , std::string line_start ) const {

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

SettingsCont *Manager::get_sc_ref() {
	return _sc_ref;
}

void Manager::add_score(double score) {
	_evaluation_scores.push_back(score);
}

void Manager::add_integration(int generation) {
	_evaluation_generations.push_back(generation);
}

void Manager::scores_to_file( ofstream& file ) {
	for ( int i_score = 0 ; i_score < _evaluation_scores.size() ; i_score++ ) {
		file << _evaluation_scores.at(i_score) << "\n";
	}
}

void Manager::evaluation_generations_to_file( ofstream& file ) {
	for ( int i_evaluation = 0 ; i_evaluation < _evaluation_generations.size() ; i_evaluation++ ) {
		file << _evaluation_generations.at(i_evaluation) << "\n";
	}
}

void Manager::mutation_generations_to_file( ofstream& file ) {
	for ( int i_generation ; i_generation < _mutation_generations.size() ; i_generation++ ) {
		file << _mutation_generations.at(i_generation) << "\n";
	}
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
	
	_curr_tissue.resize( _num_cell , _num_mol );
	
}

/* Private Method: update_mol_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Goes through every gene, protein, and reaction, and instructs them to 
 * update their internally stored data regarding their own index (for genes
 * and proteins) as well as the index of other molecules they are related to 
 * (eg. proteins in a complex, or for reactions, the reactants and products).
 *
 * The actual update method is performed by the gene and proteins themselves.
 */
void Manager::update_mol_indices( int first_index , int num_insertion ) {
	for ( int i = 0 ; i < _genes.size() ; i++ ) {
		_genes.at(i)->update_mol_indices(first_index,num_insertion);
	}
	for ( int i = 0 ; i < _proteins.size() ; i++ ) {
		_proteins.at(i)->update_mol_indices(first_index,num_insertion);
	}
	for ( int i = 0 ; i < _reactions.size() ; i++ ) {
		_reactions.at(i)->update_mol_indices(first_index,num_insertion);
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
 *
 * NB Because reactions only added at the end, no indices of existing
 * reactions require updating.
 */
void Manager::add_gene( double prod_rate , double deg_rate ) {
	
	
	int i_new_gene = _genes.size();
	int i_new_prot = _num_mol+1; /* We are using the old, not yet updated, value 
								  * of numMol 
								  */
	int i_new_gene_in_genes = i_new_gene;
	int i_new_prot_in_prots = i_new_prot - (_genes.size() + 1); // -1 b/c we haven't yet added the new gene to genes
	
	/* Step 1 */
	update_mol_indices(i_new_gene,1);
	
	/* Steps 2 - 5 */
	_genes.push_back(new Gene(i_new_gene,i_new_prot,NEXIST,NEXIST,1.0));
	_proteins.push_back(new Protein(i_new_prot,NEXIST,NEXIST,0.5));
	_reactions.push_back(new PromReaction(i_new_gene,i_new_prot,prod_rate));
	_reactions.push_back(new DegReaction(i_new_prot,deg_rate));
	
	/* Step 6 */
	_genes.at(i_new_gene_in_genes)->add_reaction(_reactions.size()-2); // to new gene
	_proteins.at(i_new_prot_in_prots)->add_reaction(_reactions.size()-2); // to new product protein
	/* Step 7 */
	_proteins.at(i_new_prot_in_prots)->add_reaction(_reactions.size()-1); // to now product protein
	
	/* Step 8 */
	_num_mol += 2;
	
	if ( _sc_ref->_verbose == 1 ) {
		std::cout << "GENE ADDITION OCCURRED" << std::endl;
	}
	
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
 *		7. Updates _num_mol.
 *
 * NB Because reactions only added at the end, no indices of existing
 * reactions require updating.
 */
void Manager::add_PromBindingReac( int i_gene_in_genes , int i_prot_in_prots ,
								  double forward_kinetic , double backward_kinetic , 
								  double new_prod_rate ) {
	
	int i_promoted_gene = _genes.size();
	
	/* Step 1 */
	update_mol_indices(i_promoted_gene,1);
	
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
	_genes.at(i_gene_in_genes)->add_reaction(_reactions.size()-2); // to root gene
	_proteins.at(i_prot_in_prots)->add_reaction(_reactions.size()-2); // to bound protein
	_genes.at(_genes.size()-1)->add_reaction(_reactions.size()-2); // to new gene:prot complex
	
	/* Step 6 */
	int i_product_prot_in_prots = i_product_protein - _genes.size();
	_proteins.at(i_product_prot_in_prots)->add_reaction(_reactions.size()-1); // to product protein
	_genes.at(_genes.size()-1)->add_reaction(_reactions.size()-1); // to new gene:prot complex
	
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
 *		5. Updates _num_mol.
 *
 * NB Because reactions only added at the end, no indices of existing
 * reactions require updating.
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
	int i_reac_one = i_prot_one_in_prots + _genes.size();
	_proteins.push_back(new Protein(i_product,i_reac_zero,i_reac_one,0.0));
	
	/* Step 3 */
	_reactions.push_back(new CombReaction(i_reac_zero,i_reac_one,i_product,
										  forward_kinetic,backward_kinetic));
	
	/* Step 4 */
	_proteins.at(i_prot_zero_in_prots)->add_reaction(_reactions.size()-1); // to first participant
	if (i_prot_zero_in_prots != i_prot_one_in_prots) {
		_proteins.at(i_prot_one_in_prots)->add_reaction(_reactions.size()-1); // to second participant
	}
	_proteins.at(_proteins.size()-1)->add_reaction(_reactions.size()-1); // to their product.
	
	/* Step 5 */
	_num_mol += 1;
}

/* Private Method: add_LatPromReac(i_promoting_neighbors_in_prots,
 *								  i_promoted_by_neighbors_in_prots,
 *								  kinetic,K)
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
 *
 * NB Because no genes or proteins are added, and reactions only added at the 
 * end, no indices of existing reactions require updating.
 */

void Manager::add_LatPromReac( int i_promoting_neighbors_in_prots , 
							  int i_promoted_by_neighbors_in_prots ,
							  double kinetic , double K ) {
	
	int i_promoting_neighbors = _genes.size() + i_promoting_neighbors_in_prots;
	int i_promoted_by_neighbors = _genes.size() + i_promoted_by_neighbors_in_prots;
	
	/* Step 1 */
	_reactions.push_back(new LatPromReaction(i_promoting_neighbors,i_promoted_by_neighbors,
											 kinetic,K));
	
	/* Step 2 */
	_proteins.at(i_promoting_neighbors_in_prots)->add_reaction(_reactions.size()-1); // to local protein
	if ( i_promoting_neighbors_in_prots != i_promoted_by_neighbors_in_prots ) {
		_proteins.at(i_promoted_by_neighbors_in_prots)->add_reaction(_reactions.size()-1); // to neighbor protein
	}
	
}

/* Private Method: add_LatRepReac(i_repressing_neighbors_in_prots,
 *								  i_repressed_by_neighbors_in_prots,
 *								  kinetic,K)
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
 *
 * NB Because no genes or proteins are added, and reactions only added at the 
 * end, no indices of existing reactions require updating.
 */

void Manager::add_LatRepReac( int i_repressing_neighbors_in_prots , 
							  int i_repressed_by_neighbors_in_prots ,
							  double kinetic , double K ) {
	
	int i_repressing_neighbors = _genes.size() + i_repressing_neighbors_in_prots;
	int i_repressed_by_neighbors = _genes.size() + i_repressed_by_neighbors_in_prots;
	
	/* Step 1 */
	_reactions.push_back(new LatRepReaction(i_repressing_neighbors,i_repressed_by_neighbors,
											 kinetic,K));
	
	/* Step 2 */
	_proteins.at(i_repressing_neighbors_in_prots)->add_reaction(_reactions.size()-1); // to local protein
	if ( i_repressing_neighbors_in_prots != i_repressed_by_neighbors_in_prots ) {
		_proteins.at(i_repressed_by_neighbors_in_prots)->add_reaction(_reactions.size()-1); // to neighbor protein
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
 *
 * NB Because no genes or proteins are added, and reactions only added at the 
 * end, no indices of existing reactions require updating.
 */

void Manager::add_HillPromReac( int i_promoter_in_prots , int i_promoted_in_prots ,
							   double kinetic , double K , double cooperativity ) {
	
	int i_promoter = _genes.size() + i_promoter_in_prots;
	int i_promoted = _genes.size() + i_promoted_in_prots;
	
	/* Step 1 */
	_reactions.push_back( new HillPromReaction(i_promoter,i_promoted,
											   kinetic,K,cooperativity) );
	
	/* Step 2 */
	_proteins.at(i_promoter_in_prots)->add_reaction(_reactions.size()-1); // to promoter protein
	if ( i_promoter_in_prots != i_promoted_in_prots ) {
		_proteins.at(i_promoted_in_prots)->add_reaction(_reactions.size()-1); // to promoted protein
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
 *
 * NB Because no genes or proteins are added, and reactions only added at the 
 * end, no indices of existing reactions require updating.
 */

void Manager::add_HillRepReac( int i_repressor_in_prots , int i_repressed_in_prots ,
							  double kinetic , double K , double cooperativity ) {
	
	int i_repressor = _genes.size() + i_repressor_in_prots;
	int i_repressed = _genes.size() + i_repressed_in_prots;
	
	/* Step 1 */
	_reactions.push_back(new HillRepReaction(i_repressor,i_repressed,
											 kinetic,K,cooperativity));
	
	/* Step 2 */
	_proteins.at(i_repressor_in_prots)->add_reaction(_reactions.size()-1); // to promoter protein
	if ( i_repressor_in_prots != i_repressed_in_prots ) {
		_proteins.at(i_repressed_in_prots)->add_reaction(_reactions.size()-1); // to promoted protein
	}
}

/* Private Method: remove_gene( int i_gene )
 * -------------------------------------------------------------------------- 
 * Removes the gene indexed by i_gene, as well as all things dependent on it.
 *
 * Finding dependencies is handled by gene_removal_cascade. See 
 * Documentation.txt for more info on reaction removal and the gene removal
 * cascade.
 */

void Manager::remove_gene( int i_gene ) {
	std::set<int> reacs_to_delete;
	std::set<int> genes_to_delete;
	std::set<int> prots_to_delete;
	
	genes_to_delete.insert(i_gene);
	gene_removal_cascade(i_gene,reacs_to_delete,genes_to_delete,prots_to_delete);
	
	perform_removal_using_sets(reacs_to_delete,genes_to_delete,prots_to_delete);
}

/* Private Method: remove_reaction(i_reac)
 * --------------------------------------------------------------------------
 * WARNING: As this method stands right now, the genes and proteins will not
 * maintain accurate _reactions vectors. Also, will cause problems if i_reac
 * not a valid index.
 */

void Manager::remove_reaction( int i_reac ) {
	
	std::set<int> reacs_to_delete;
	std::set<int> genes_to_delete;
	std::set<int> prots_to_delete;
	
	reacs_to_delete.insert(i_reac);
	reac_removal_cascade(i_reac,reacs_to_delete,genes_to_delete,prots_to_delete);
	
	perform_removal_using_sets(reacs_to_delete,genes_to_delete,prots_to_delete);
}

/* Private Method: reac_removal_cascade(i_reac_to_remove,reacs_to_delete,
 *										genes_to_delete,prots_to_delete)
 * --------------------------------------------------------------------------
 * Adds to the reacs_to_delete, genes_to_delete, and prots_to_delete sets the
 * indices of all genes, proteins, and reactions that need to be deleted due
 * to the removal of reaction indexed by i_reac_to_remove.
 *
 * Performs via recursion, using the gene_removal_cascade and
 * prot_removal_cascade methods.
 *
 * We assume i_reac_to_remove has already been added to the reacs_to_delete
 * method when the function is called. Thus, the function only deals with
 * dependencies.
 *
 * First checks if any molecules must be added due to dependencies on this
 * reaction.
 *
 * If so, it, based on whether the molecules is a gene or a protein, performs
 * the following method:
 *
 *		1. Check whether this molecule has already been added as a molecule
 *		to delete. Check by seeing if adding the molecule to the 
 *		corresponding set (gene or protein) changes the size of the set.
 *		2. If it has not yet been added (the set size did change at its
 *		addition), we call the appropriate cascade for that molecule.
 *
 * Note that because the next molecule is added to the corresponding set
 * before we call its cascade, we maintain the convention that the cascade
 * is only called once the substance it is cascading off of has been added
 * to the corresponding to_delete set. This assures that the recursion will
 * eventually return.
 *
 * See Documetnation.txt for the details of the dependencies that govern
 * the cascade.
 */

void Manager::reac_removal_cascade( int i_reac_to_remove , 
								   std::set<int>& reacs_to_delete ,
								   std::set<int>& genes_to_delete , 
								   std::set<int>& prots_to_delete ) {
	
	int i_mol_to_delete = _reactions.at(i_reac_to_remove)->get_i_dependent_molecule();
	if (i_mol_to_delete == NEXIST) return;
	
	if (i_mol_to_delete < _genes.size()) { // if molecule is a gene
		
		int genes_size = genes_to_delete.size();
		genes_to_delete.insert(i_mol_to_delete);
		/* Return if gene was already in genes_to_delete */
		if (genes_to_delete.size() == genes_size) return;
		
		else gene_removal_cascade( i_mol_to_delete , reacs_to_delete , 
								 genes_to_delete , prots_to_delete );
		
	}
	
	else {
		
		int i_prot_to_delete = i_mol_to_delete - _genes.size();
		int prots_size = prots_to_delete.size();
		prots_to_delete.insert(i_prot_to_delete);
		/* Return if prot was already in prots_to_delete */
		if (prots_to_delete.size() == prots_size) return;
		
		else prot_removal_cascade( i_prot_to_delete , reacs_to_delete , 
								 genes_to_delete , prots_to_delete );
		
	}
	
}

/* Private Method: gene_removal_cascade(i_gene_to_remove,reacs_to_delete,
 *										genes_to_delete,prots_to_delete)
 * --------------------------------------------------------------------------
 * Adds to the reacs_to_delete, genes_to_delete, and prots_to_delete sets the
 * indices of all genes, proteins, and reactions that need to be deleted due
 * to the removal of gene indexed by i_gene_to_remove.
 *
 * Performs via recursion, using the reac_removal_cascade and
 * prot_removal_cascade methods.
 *
 * We assume i_gene_to_remove has already been added to the genes_to_delete
 * method when the function is called. Thus, the function only deals with
 * dependencies.
 *
 * See Documetnation.txt for the details of the dependencies that govern
 * the cascade.
 */

void Manager::gene_removal_cascade( int i_gene_to_delete ,
								  std::set<int>& reacs_to_delete ,
								  std::set<int>& genes_to_delete ,
								  std::set<int>& prots_to_delete ) {

	for ( std::set<int>::iterator it = _genes.at(i_gene_to_delete)->reacs_begin() ; 
		 it != _genes.at(i_gene_to_delete)->reacs_end(); 
		 it++ ) {
		
		int reacs_size = reacs_to_delete.size();
		reacs_to_delete.insert(*it);
		/* If *it wasn't already among those to delete, we cascade */
		if (reacs_size != reacs_to_delete.size()) reac_removal_cascade(*it,
																	   reacs_to_delete,
																	   genes_to_delete,
																	   prots_to_delete);
	}
	
	/* If this is a base gene, then deleting it should delete the protein as well. */
	/* Because the protein is not considered dependent on the promotion reaction, 
	 * we need to eliminate the protein separately */
	if (_genes.at(i_gene_to_delete)->get_i_root() == NEXIST && 
		_genes.at(i_gene_to_delete)->get_i_bound_promoter() == NEXIST) {
		int i_prot_to_delete = _genes.at(i_gene_to_delete)->get_i_product()-_genes.size();
		int prots_size = prots_to_delete.size();
		prots_to_delete.insert(i_prot_to_delete);
		if (prots_size != prots_to_delete.size()) prot_removal_cascade(_genes.at(i_gene_to_delete)->get_i_product()-_genes.size(),
																	   reacs_to_delete,
																	   genes_to_delete,
																	   prots_to_delete);
	}
	
}

/* Private Method: prot_removal_cascade(i_prot_to_remove,reacs_to_delete,
 *										genes_to_delete,prots_to_delete)
 * --------------------------------------------------------------------------
 * Adds to the reacs_to_delete, genes_to_delete, and prots_to_delete sets the
 * indices of all genes, proteins, and reactions that need to be deleted due
 * to the removal of the protein indexed by i_prot_to_remove.
 *
 * i_prot_to_remove is the index of the protein in the proteins vector, and
 * not in the entire dvec of the tissue.
 *
 * Performs via recursion, using the reac_removal_cascade and
 * gene_removal_cascade methods.
 *
 * We assume i_prot_to_remove has already been added to the prots_to_delete
 * method when the function is called. Thus, the function only deals with
 * dependencies.
 *
 * See Documetnation.txt for the details of the dependencies that govern
 * the cascade.
 */
void Manager::prot_removal_cascade( int i_prot_to_delete ,
								   std::set<int>& reacs_to_delete ,
								   std::set<int>& genes_to_delete ,
								   std::set<int>& prots_to_delete ) {
	
	
	std::set<int>::iterator it = _proteins.at(i_prot_to_delete)->reacs_begin();
	int test = *it;
	test = test;
	
	for ( std::set<int>::iterator it = _proteins.at(i_prot_to_delete)->reacs_begin() ;
		 it != _proteins.at(i_prot_to_delete)->reacs_end() ;
		 it++ ) {
		
		int reacs_size = reacs_to_delete.size();
		reacs_to_delete.insert(*it);
		/* If *it wasn't already among those to delete, we cascade */
		if (reacs_size != reacs_to_delete.size()) reac_removal_cascade(*it,
																	   reacs_to_delete,
																	   genes_to_delete,
																	   prots_to_delete);
	}
	
}

/* Private Method: perform_removal_using_sets(i_prot_to_remove,reacs_to_delete,
 *											genes_to_delete,prots_to_delete)
 * --------------------------------------------------------------------------
 * Removes the reactions, genes, and proteins from the genome using the sets
 * passed to it. 
 *
 * This method does nothing regarding dependencies. It is assumed that reacs,
 * genes, and prots to delete have already been found according to a removal
 * cascade so that all dependencies are already taken care of. 
 *
 * It removes reactions or molecules one at a time, and in between updates
 * all necessary indices, including those in teh reacs_to_delete,
 * genes_to_delete, and prots_to_delete sets. 
 *
 * Updates are made using the update index method or function class from
 * the Operations namespace. This, referenced from all part of the code,
 * makes it easy to edit the update index method as necessary and maintain
 * consistancy throughout the code. Useful for correcting errors or adding
 * more complicated update procedures.
 */
void Manager::perform_removal_using_sets(std::set<int>& reacs_to_delete,
										  std::set<int>& genes_to_delete,
										  std::set<int>& prots_to_delete) {
	/* Delete reactions */
	while (!reacs_to_delete.empty()) {
		int i_deletion = *reacs_to_delete.begin();
		Operations::UpdateIndex updater(i_deletion,-1);
		for ( int i_gene = 0 ; i_gene < _genes.size() ; i_gene++ ) {
			_genes.at(i_gene)->update_reac_indices(i_deletion,-1);
		}
		for ( int i_prot = 0 ; i_prot < _proteins.size() ; i_prot++ ) {
			_proteins.at(i_prot)->update_reac_indices(i_deletion,-1);
		}
		std::set<int> updates;
		for (std::set<int>::iterator it = reacs_to_delete.begin(); it != reacs_to_delete.end(); it++) {
			updates.insert(updater(*it));
		}
		delete _reactions.at(i_deletion);
		reacs_to_delete = updates;
		reacs_to_delete.erase(NEXIST);
		for ( int i = i_deletion + 1 ; i < _reactions.size() ; i++ ) {
			_reactions.at(i-1) = _reactions.at(i);
		}
		_reactions.pop_back();
		
	}
	
	/* Delete genes */
	while (!genes_to_delete.empty()) {
		int i_deletion = *genes_to_delete.begin();
		update_mol_indices(i_deletion, -1);
		std::set<int> updates;
		Operations::UpdateIndex updater(i_deletion,-1);
		for (std::set<int>::iterator it = genes_to_delete.begin(); it != genes_to_delete.end(); it++) {
			updates.insert(updater(*it));
		}
		_num_mol--;
		delete _genes.at(i_deletion);
		for	( int i = i_deletion +1 ; i < _genes.size() ; i++ ) {
			_genes.at(i-1) = _genes.at(i);
		}
		_genes.pop_back();
		genes_to_delete = updates;
		genes_to_delete.erase(NEXIST);
	}
	
	/* Delete proteins */
	while (!prots_to_delete.empty()) {
		int i_deletion = *prots_to_delete.begin();
		int i_deletion_as_mol = i_deletion + _genes.size();
		update_mol_indices(i_deletion_as_mol, -1);
		std::set<int> updates;
		Operations::UpdateIndex updater(i_deletion,-1);
		for (std::set<int>::iterator it = prots_to_delete.begin(); it != prots_to_delete.end(); it++) {
			updates.insert(updater(*it));
		}
		delete _proteins.at(i_deletion);
		_num_mol--;
		for (int i = i_deletion+1; i < _proteins.size() ; i++) {
			_proteins.at(i-1) = _proteins.at(i);
		}
		_proteins.pop_back();
		prots_to_delete = updates;
		prots_to_delete.erase(NEXIST);
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
		if ( _sc_ref->_verbose == 1 ) {
			std::cout << "DEGREDATION MUTATION OCCURRED (Reaction " << i_deg.at(i_reac) << ")" << std::endl;
		}
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
		if ( _sc_ref->_verbose == 1 ) {
			std::cout << "MUTATION OF KINETIC CONSTANT OCCURRED (Reaction " << i_ndeg.at(i_reac) << ")" << std::endl;
		}
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
		if ( _sc_ref->_verbose == 1 ) {
			std::cout << "PROMOTER BINDING MUTATION OCCURRED" << std::endl;
		}
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
		
		if ( _sc_ref->_verbose == 1 ) {
			std::cout << "POST TRANSCRIPTION MUTATION OCCURRED" << std::endl;
		}
	}
	
}

/* Private Method: add_hill_mutation(generator)
 * -------------------------------------------------------------------------- 
 */
void Manager::add_hill_mutation( boost::random::mt19937& generator ) {
	
	if ( _proteins.size() == 0 ) return;
	
	boost::random::uniform_int_distribution<> reaction_type(0,3);
	boost::random::uniform_int_distribution<> which_protein(0,_proteins.size()-1);
	boost::random::uniform_real_distribution<> rand_real(0,2.0);
	
	for ( int trial_num = 0 ; trial_num < 5 ; trial_num++ ) {
	
		int type = reaction_type(generator);
		double kinetic = .5+rand_real(generator);
		double K = rand_real(generator);
		double cooperativity = 5.0; // We know cooperativities should general be more than one
		int actor_protein = which_protein(generator);
		int receiver_protein = which_protein(generator);
		
		switch (type) {
			case 0:
				add_HillRepReac(actor_protein,receiver_protein,kinetic,K,cooperativity);
				break;
			case 1:
				add_HillPromReac(actor_protein,receiver_protein,kinetic,K,cooperativity);
				break;
			case 2:
				add_LatPromReac(actor_protein, receiver_protein, kinetic, K);
				break;
			case 3:
				add_LatRepReac(actor_protein, receiver_protein, kinetic, K);
			default:
				break;
		}
		
		/* Remove the added reaction if it is a duplicate */
		int temp_trial_num = trial_num;
		trial_num = 5; // We only run loop again if we find no duplicates of this reaction
		int num_reac = _reactions.size();
		for ( int i_reac = 0 ; i_reac < num_reac-1 ; i_reac++ ) {
			if ( *(_reactions.at(i_reac)) == *(_reactions.at(num_reac-1)) ) {
				remove_reaction(num_reac-1);
				/* Set looping parameters so we return */
				trial_num = temp_trial_num;
				i_reac = num_reac-1;
			}
		}
	}
	
	if ( _sc_ref->_verbose == 1 ) {
		std::cout << "HILL MUTATION OCCURRED\n";
	}

}

/* Private Method: remove_reaction(generator)
 * --------------------------------------------------------------------------
 */

void Manager::removal_mutation( boost::random::mt19937& generator ) {
	std::vector<int> i_candidates;
	for ( int i = 0 ; i < _reactions.size() ; i++ ) {
		if ( (_reactions.at(i)->get_type()) != DEGRADATION && (_reactions.at(i)->get_type()) != PROMOTION ) {
			i_candidates.push_back(i);
		}
	}

	if (i_candidates.size() > 0) {
		boost::random::uniform_int_distribution<> which_reaction(0,i_candidates.size()-1);
		remove_reaction(i_candidates.at(which_reaction(generator)));
	}
	if ( _sc_ref->_verbose == 1 ) {
		std::cout << "REMOVAL MUTATION OCCURRED\n";
	}
	
}

/* Private Method: rk4_det_ti(num_step)
 * --------------------------------------------------------------------------
 * At each step, performs the 4th order deterministic runge kutta update. 
 * Calls the Manager::react method which given the curr_tissue reads into
 * the dx_dt dmat passed by reference the appropriate rates of change of 
 * reactants at this time step.
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
 * Equivalent to the ordinary 4th order deterministic method, but reads into
 * the provided file the concentrations of each substance in each cell at
 * each point in time.
 *
 * file is considered to be open and pointing to the beginning of the
 * desired output file. The file is NOT closed before returning.
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
 * ---------------------------------------------------------------------------
 * At each step, performs the 4th order stochastic runge kutta update designed
 * by Jeremy Kasdin at Princeton University. 
 *
 * See http://people.sc.fsu.edu/~jburkardt/m_src/stochastic_rk/stochastic_rk.html
 *
 * At each step it reads it calls the stochastic Manager::react method, which 
 * reads into the dx_dt dmat passed by reference the appropriate rates of
 * chante at that time step given the curr_state dmat. 
 *
 * For a clear specification of what the dx_dt dmat represents in the
 * stochastic case, see Documentation.txt.
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
 * Same as orindary 4th order stochastic runge kutta method, but reads into
 * the provided file the concentrations of each substance in each cell at
 * each point in time.
 *
 * file is considered to be open and pointing to the beginning of the
 * desired output file. The file is NOT closed before returning.
 * step reads into t
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

/* Private Method: hakim_delta_notch_construct()
 * --------------------------------------------------------------------------
 * HAKIM_DELTA_NOTCH constructor helper. 
 * The constructed network attempts to make a Hakim Delta Notch network, 
 * though given further tests, the appropriate altnerating pattern only was
 * found for initial concentrations that already had the desired alternating.
 */

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

/* Private Method: collier_delta_notch_construct()
 * --------------------------------------------------------------------------
 * COLLER_DELTA_NOTCH constructor helper. 
 * Creates the exact genetic network studies in Collier et al.'s paper.
 * (ie. "Pattern Formation by Lateral Inhibition with Feedback: a 
 * Mathematical Model of Delta-Notch Intercellular Signalling.")
 */
void Manager::collier_delta_notch_construct() {
	_num_mol = 0;
	
	add_gene(0.0,1.0); // Delta
	/* < d D > */
	add_gene(0.0,1.0); // Notch
	/* < d n D N > */
	add_LatPromReac(0,1,1.0,.01); // dN/dt = (D_n)^2 / (.01 + (D_n)^2)
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

/* Private Method: mutate_construct()
 * --------------------------------------------------------------------------
 * MUTATION constructor helper. 
 * Creates genome by beginning with a gene, and then performing 5 mutations. 
 */
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

/* Private Method: one_protein_construct()
 * --------------------------------------------------------------------------
 * ONE_PROTEIN constructor helper. 
 * Creates genome consisting of just a gene and a protein.
 */
void Manager::one_protein_construct() {
	_num_mol = 0;
	
	add_gene(1.0,1.0);
	
	initialize();
}

/* Private Method: two_protein_construct()
 * --------------------------------------------------------------------------
 * TWO_PROTEIN constructor helper. 
 * Creates genome consisting of two gene and their related proteins.
 *
 * promotion rate = 0. Degredation rate = 1.
 */
void Manager::two_protein_construct() {
	_num_mol = 0;
	add_gene(0.0,1.0);
	add_gene(0.0,1.0);
	initialize();
}
