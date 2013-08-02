/*
 *  Manager.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

#include "Manager.h"

/* Constructor: Manager()
 * -------------------------------------------------------------------------- 
 * Default implementation of constructor for use during development.
 *
 * Currently attempts to create a Delta-Notch like network.
 */

Manager::Manager() {
	
	_num_mol = 0; /* Before we add anything to the genome. Otherwise has random
				 * initiation.
				 */
	
	_num_cell = 3;
	
	/* Create neighbors vector: Linear tissue three long */
	/* < [2] , [1,3] , [2] > */
	_neighbors.push_back(new vector<int>(1,1));
	vector<int>* addition = new vector<int>;
	addition->push_back(0);
	addition->push_back(2);
	_neighbors.push_back(addition);
	addition = new vector<int>(1,1);
	_neighbors.push_back(addition);
	
	/* Initialize tissue vectors */
	/* The substances in our genome are, in order, 
	 * 
	 *		[ d   n   d:N   d:N:N   D   N]
	 *
	 * Where lower-case refers to a gene, upper-case to a protein.
	 *
	 * Initial concentrations are
	 *
	 *		Cell 1: [ 1.0  1.0  0.0  0.0  0.0  1.0 ]
	 *		Cell 2: [ 1.0  1.0  0.0  0.0  1.0  0.0 ]
	 *		Cell 3: [ 1.0  1.0  0.0  0.0  0.0  1.0 ]
	 *
	 */
	/*
	double nums[5];
	nums[0] = 1;
	nums[1] = 1;
	nums[2] = 0;
	nums[3] = 0;
	nums[4] = 0.0;
	nums[5] = 1.0;
	iTissue.push_back(new dvec(6,nums));
	nums[4] = 1.0;
	nums[5] = 0.0;
	iTissue.push_back(new dvec(6,nums));
	nums[4] = 0.0;
	nums[5] = 1.0;
	iTissue.push_back(new dvec(6,nums));
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
	
	/*
	genes.push_back(new Gene(0,4,-1,-1));
	genes.push_back(new Gene(1,5,-1,-1));
	genes.push_back(new Gene(2,4,5,0));
	genes.push_back(new Gene(3,4,5,2));
	proteins.push_back(new Protein(4,-1,-1));
	proteins.push_back(new Protein(5,-1,-1));
	*/
	
	/* Initialize time variables */
	_dt = .001;
	_time = 0.0;
	
	/* Initialize Integration Mode */
	_mode = RK4_DET_TI;
	
	/* Initialize */
	initialize();
	
	
	_curr_tissue.at(0,3) = 1.0;
	_curr_tissue.at(0,4) = 0.0;
	_curr_tissue.at(1,3) = 0.0;
	_curr_tissue.at(1,4) = 1.0;
	_curr_tissue.at(1,5) = 0.0;
	_curr_tissue.at(2,3) = 1.0;
	_curr_tissue.at(2,3) = 0.0;
	 
	/*
	_curr_tissue.at(0,0) = .942203;
	_curr_tissue.at(0,1) = 1.0;
	_curr_tissue.at(0,2) = .0577972;
	_curr_tissue.at(0,3) = .942328;
	_curr_tissue.at(0,4) = .00112699;
	_curr_tissue.at(0,5) = .00614686;
	
	_curr_tissue.at(1,0) = .00330837;
	_curr_tissue.at(1,1) = 1.0;
	_curr_tissue.at(1,2) = .996692;
	_curr_tissue.at(1,3) = .00335278;
	_curr_tissue.at(1,4) = .2987579;
	_curr_tissue.at(1,5) = 20.1277;
	
	_curr_tissue.at(2,0) = .942203;
	_curr_tissue.at(2,1) = 1.0;
	_curr_tissue.at(2,2) = .0577972;
	_curr_tissue.at(2,3) = .942328;
	_curr_tissue.at(2,4) = .00112699;
	_curr_tissue.at(2,5) = .00614686;
	*/
	
}

/* Destructor: ~Manager() 
 * -------------------------------------------------------------------------- 
 */
Manager::~Manager() {

	_reactions.erase(_reactions.begin(),_reactions.end());
	_genes.erase(_genes.begin(),_genes.end());
	_proteins.erase(_proteins.begin(),_proteins.end());
	_neighbors.erase(_neighbors.begin(),_neighbors.end());
	
}

/* Public Method: set_mode(mode)
 * -------------------------------------------------------------------------- 
 * Sets integration mode of Manager.
 */
void Manager::set_mode( IntegrationType mode ) {
	_mode = mode;
}

/* Public Method: get_mode()
 * -------------------------------------------------------------------------- 
 * Returns current integration mode of Manager.
 */
IntegrationType Manager::get_mode() {
	return _mode;
}

/* Public Method: set_dt(dt)
 * -------------------------------------------------------------------------- 
 * Sets time step size to be used in integration.
 */
void Manager::set_dt( double dt ) {
	_dt = dt;
}

/* Public Method: get_dt()
 * -------------------------------------------------------------------------- 
 * Returns current time step size used in integration.
 */
double Manager::get_dt() {
	return _dt;
}

/* Public Method: add_gene(prod_rate)
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
	
}

/* Public Method: addPromBindingReaction(iGeneInGenes,iProtInProts,
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

/* Public Method: add_CombReaction(i_prot_zero_in_prots,_i_prot_one_in_prots,
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

/* Public Method: add_LatPromReac(i_loc_prot_in_prots,
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

void Manager::add_LatPromReac( int i_loc_prot_in_prots , 
							  int i_neighbor_prot_in_prots ,
							  double kinetic , double K ) {
	
	int i_local_prot = _genes.size() + i_loc_prot_in_prots;
	int i_neighbor_prot = _genes.size() + i_neighbor_prot_in_prots;
	
	/* Step 1 */
	_reactions.push_back(new LatPromReaction(i_local_prot,i_neighbor_prot,
											 kinetic,K));
	
	/* Step 2 */
	Reaction* curr_reac_ref = _reactions.at(_reactions.size()-1); // add new LatPromReaction
	_proteins.at(i_loc_prot_in_prots)->add_reaction(curr_reac_ref); // to local protein
	if ( i_loc_prot_in_prots != i_neighbor_prot_in_prots ) {
		_proteins.at(i_neighbor_prot_in_prots)->add_reaction(curr_reac_ref); // to neighbor protein
	}
	
}

/* Public Method: initialize()
 * -------------------------------------------------------------------------- 
 * Initialize by making the dmat the appropriate size according to _num_mol,
 * and then filling each dmat with the initial concentrations of the proteins 
 * and genes, and setting the time = 0.
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

/* Public Method: integrate(mode,dt,num_step)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the currTissue vector based on iTissue, 
 * reactions, and mode.
 */
void Manager::integrate( int num_step ) {
	
	switch (_mode) {
		
		case RK4_DET_TI:
			rk4_det_ti(num_step);
			break;
			
		case RK4_STC_TI:
			rk4_stc_ti(num_step);
			break;
			
		default:
			break;
	}
}

/* Public Method: print_state()
 * -------------------------------------------------------------------------- 
 * For the user.
 * Currently goes through each cell in the vector and prints out the 
 * protein concentrations in each cell, (assuming there are 6 of them).
 */
void Manager::print_state() {
	std::cout << "time: " << _time << std::endl;
	for ( int i = 0 ; i < _num_cell ; i++ ) {
		std::cout << "Cell: " << i << std::endl;
		std::cout << _curr_tissue.at(i,0) << " :::: " <<
		_curr_tissue.at(i,1) << " :::: " << _curr_tissue.at(i,2) <<
		" :::: " << _curr_tissue.at(i,3) << " :::: " << _curr_tissue.at(i,4) <<
		" :::: " << _curr_tissue.at(i,5) << std::endl;
	}
	std::cout << std::endl;
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
			_reactions.at(i_reac)->react( xi , dx_dt , _neighbors , i_cell );
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
	
	for ( int step = 0 ; step < num_step ; step++ ) {
		
		/* ----------------------------------------- */
		react ( _curr_tissue , dx_dt );
		k1.equals_lin_comb( dx_dt , _dt );
		/* ----------------------------------------- */
		x2.equals_lin_comb(_curr_tissue , 1.0 ,
						   k1 , .5 );
		react ( x2 , dx_dt );
		k2.equals_lin_comb( dx_dt , _dt );
		/* ----------------------------------------- */
		x3.equals_lin_comb(_curr_tissue , 1.0 ,
						   k2 , .5 );
		react ( x3 , dx_dt );
		k3.equals_lin_comb( dx_dt , _dt );
		/* ----------------------------------------- */
		x4.equals_lin_comb(_curr_tissue,1.0 ,
						   k3,1.0 );
		react ( x4 , dx_dt );
		k4.equals_lin_comb( dx_dt , _dt);
		/* ----------------------------------------- */
		_curr_tissue.plus_equals_lin_comb(k1,.166666666667,
										  k2,.333333333333,
										  k3,.333333333333,
										  k4,.166666666667);
		/* ----------------------------------------- */
		_time += _dt;
		print_state();
	}
}

/* Private Method: rk4_stc_ti(num_step)
 * -------------------------------------------------------------------------- 
 */

void Manager::rk4_stc_ti( int num_step ) {
	
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
	double t1;
	double t2;
	double t3;
	double t4;
	double w1;
	double w2;
	double w3;
	double w4;
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
	
	for ( int step = 0 ; step < num_step ; step++ ) {
		
		
		
		_time += _dt;
		print_state();
	}
}

