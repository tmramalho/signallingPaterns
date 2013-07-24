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
	
	numMol = 0; /* Before we add anything to the genome. Otherwise has random
				 * initiation.
				 */
	
	/* Create neighbors vector: Linear tissue three long */
	/* < [2] , [1,3] , [2] > */
	neighbors.push_back(new vector<int>(1,1));
	vector<int>* addition = new vector<int>;
	addition->push_back(0);
	addition->push_back(2);
	neighbors.push_back(addition);
	addition = new vector<int>(1,1);
	neighbors.push_back(addition);
	
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
	
	for ( int i = 0 ; i < 3 ; i++ ) {
		iTissue.push_back(new dvec());
		currTissue.push_back(new dvec());
		dx.push_back(new dvec());
	}
	
	
	/* Initialize gene and protein vectors */
	/* Reference detailed constructors and array of
	 * elements of genome above to understand 
	 * paramaters.
	 */
	addGene(1.0,1.0); // Delta
	/* < d D > */
	addGene(0.0,2.0); // Notch
	/* < d n D N > */
	addPromBindingReac(0 , 1 , 1.0 , 0.1 , 1.0); // d + N --> d:N
	/* < d n d:N D N > */
	addPromBindingReac(2 , 1 , 1.0 , 0.1 , 0.0); // d:N + N --> d:N:N
	/* < d n d:N d:N:N D N > */
	addLatPromReac(1 , 0 , 1.0 , .25 ); // N promoted by neighboring D
	
	/*
	genes.push_back(new Gene(0,4,-1,-1));
	genes.push_back(new Gene(1,5,-1,-1));
	genes.push_back(new Gene(2,4,5,0));
	genes.push_back(new Gene(3,4,5,2));
	proteins.push_back(new Protein(4,-1,-1));
	proteins.push_back(new Protein(5,-1,-1));
	*/
	
	/* Initialize time variables */
	dt = .01;
	
	/* Initialize Integration Mode */
	mode = RK1_DET_TI;
	
	/* Initialize */
	initialize();
	
}

/* Destructor: ~Manager() 
 * -------------------------------------------------------------------------- 
 */
Manager::~Manager() {
	
	for ( int i = 0 ; i < iTissue.size() ; i++ ) {
		delete iTissue.at(i);
		delete currTissue.at(i);
		delete dx.at(i);
	}
	
	for ( int i = 0 ; i < reactions.size() ; i++ ) {
		delete reactions.at(i);
	}
	
	for ( int i = 0 ; i < genes.size() ; i++ ) {
		delete genes.at(i);
	}
	
	for ( int i = 0 ; i < proteins.size() ; i++ ) {
		delete proteins.at(i);
	}
	
	neighbors.erase(neighbors.begin(),neighbors.begin()+2);
	
}

/* Public Method: setMode(mode)
 * -------------------------------------------------------------------------- 
 * Sets integration mode of Manager.
 */
void Manager::setMode( IntegrationType mode ) {
	this->mode = mode;
}

/* Public Method: getMode()
 * -------------------------------------------------------------------------- 
 * Returns current integration mode of Manager.
 */
IntegrationType Manager::getMode() {
	return mode;
}

/* Public Method: setDt(dt)
 * -------------------------------------------------------------------------- 
 * Sets time step size to be used in integration.
 */
void Manager::setDt( double dt ) {
	this->dt = dt;
}

/* Public Method: getDt()
 * -------------------------------------------------------------------------- 
 * Returns current time step size used in integration.
 */
double Manager::getDt() {
	return dt;
}

/* Public Method: addGene(prodRate)
 * -------------------------------------------------------------------------- 
 * Adds a new gene to the genome, as well as the corresonding protein, and
 * sets the promotion rate to prodRate.
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
 *		8. Updates numMol.
 */
void Manager::addGene( double prodRate , double degRate ) {
	
	
	int iNewGene = genes.size();
	int iNewProt = numMol+1; /* We are using the old, not yet updated, value 
							  * of numMol 
							  */
	int iNewGeneInGenes = iNewGene;
	int iNewProtInProts = iNewProt - (genes.size() + 1); // -1 b/c we haven't yet added the new gene to genes
	
	/* Step 1 */
	updateIndices(iNewGene,1);
	
	/* Steps 2 - 5 */
	genes.push_back(new Gene(iNewGene,iNewProt,NEXIST,NEXIST,1.0));
	proteins.push_back(new Protein(iNewProt,NEXIST,NEXIST,0.0));
	reactions.push_back(new PromReaction(iNewGene,iNewProt,0.0,prodRate));
	reactions.push_back(new DegReaction(iNewProt,0.0,degRate));
	
	/* Step 6 */
	Reaction* currReac = reactions.at(reactions.size()-2); // add new PromReaction
	genes.at(iNewGeneInGenes)->addReaction(currReac); // to new gene
	proteins.at(iNewProtInProts)->addReaction(currReac); // to new product protein
	
	/* Step 7 */
	currReac = reactions.at(reactions.size()-1); // add new DegReaction
	proteins.at(iNewProtInProts)->addReaction(currReac); // to now product protein
	
	/* Step 8 */
	numMol += 2;
	
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
void Manager::addPromBindingReac( int iGeneInGenes , int iProtInProts ,
								 double forwardKinetic , double backwardKinetic , 
								 double newProdRate ) {
	
	int iPromotedGene = genes.size();
	
	/* Step 1 */
	updateIndices(iPromotedGene,1);
	
	/* We find these indices after step 1 to assure they are consistent with
	 * the new indexing.
	 */
	int iBoundProtein = proteins.at(iProtInProts)->getISelf();
	int iRootGene = genes.at(iGeneInGenes)->getISelf();
	int iProductProtein = genes.at(iGeneInGenes)->getIProduct();
	
	/* Steps 2 */
	genes.push_back(new Gene(iPromotedGene,iProductProtein,iBoundProtein,
							 iRootGene,0.0));
	
	/* Steps 3 - 4 */
	reactions.push_back(new PromBindingReaction(iRootGene,iPromotedGene,iBoundProtein,
												0.0,0.0,forwardKinetic,backwardKinetic));
	reactions.push_back(new PromReaction(iPromotedGene,iProductProtein,0.0,newProdRate));
	
	/* Step 5 */
	Reaction* currReac = reactions.at(reactions.size()-2); // add new PromBindingReaction
	genes.at(iGeneInGenes)->addReaction(currReac); // to root gene
	proteins.at(iProtInProts)->addReaction(currReac); // to bound protein
	genes.at(genes.size()-1)->addReaction(currReac); // to new gene:prot complex
	
	/* Step 6 */
	int iProductProtInProts = iProductProtein - genes.size();
	currReac = reactions.at(reactions.size()-1); // add new PromReaction
	proteins.at(iProductProtInProts)->addReaction(currReac); // to product protein
	genes.at(genes.size()-1)->addReaction(currReac); // to new gene:prot complex
	
	/* Step 7 */
	numMol += 1;
	
}

/* Public Method: addCombReaction(iProtZeroInProts,iProtOneInProts,
 *								  forwardKinetic,backwardKinetic)
 * -------------------------------------------------------------------------- 
 * Makes the proteins in positions iProtZeroInProts and iProtOneInProts
 * (which can be the same) combine to form a new protein complex. With the
 * forward and backward rates as given.
 *
 * CAUTION: It is assumed that iProtZeroInProts and iProtOneInProts are the 
 * indices of proteins No check is performed by the function to guarantee 
 * this.
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

void Manager::addCombReaction( int iProtZeroInProts , int iProtOneInProts , 
							  double forwardKinetic , double backwardKinetic ) {
	
	int iProduct = numMol; // The yet to be updated value of numMol
	
	/* Step 1 */
	/* No updates are necessary because the protein is added to the end of our
	 * dvces.
	 */
	
	/* Step 2 */
	int iReacZero = iProtZeroInProts + genes.size();
	int iReacOne = iProtZeroInProts + genes.size();
	proteins.push_back(new Protein(iProduct,iReacZero,iReacOne,0.0));
	
	/* Step 3 */
	reactions.push_back(new CombReaction(iReacZero,iReacOne,iProduct,
										 0.0,0.0,0.0,
										 forwardKinetic,backwardKinetic));
	
	/* Step 4 */
	Reaction* currReac = reactions.at(reactions.size()-1); // add new CombReaction
	proteins.at(iProtZeroInProts)->addReaction(currReac); // to first participant
	if (iProtZeroInProts != iProtOneInProts) {
		proteins.at(iProtOneInProts)->addReaction(currReac); // to second participant
	}
	proteins.at(proteins.size()-1)->addReaction(currReac); // to their product.
	
	/* Step 5 */
	numMol += 1;
}

/* Public Method: addLatPromReac(iLocProtInProts,iNeighborProtInProts,
 *								 double kinetic,double K)
 * -------------------------------------------------------------------------- 
 * Makes the protein in position iLocProtInProts in the protein vector
 * get promoted by the presence of the protein iNeighborProtInProts in
 * neighboring cells, with the provided kinetic constants.
 * 
 * CAUTION: It is assumed that iLocProtInProts and iNeighborProtInProts is 
 * the index of a protein. No check is performed by the function to guarantee 
 * this.
 *
 * The tasks that need to be performed are:
 *
 *		1. A new reaction is added describing the lateral promotion reaction.
 *		2. This reaction is added to the internal data of the participating
 *		proteins.
 */

void Manager::addLatPromReac( int iLocProtInProts , int iNeighborProtInProts ,
							 double kinetic , double K ) {
	
	int iLocalProt = genes.size() + iLocProtInProts;
	int iNeighborProt = genes.size() + iNeighborProtInProts;
	
	/* Step 1 */
	reactions.push_back(new LatPromReaction(iLocalProt,iNeighborProt,
											0.0,kinetic,K));
	
	/* Step 2 */
	Reaction* currReac = reactions.at(reactions.size()-1); // add new LatPromReaction
	proteins.at(iLocProtInProts)->addReaction(currReac); // to local protein
	if ( iLocProtInProts != iNeighborProtInProts ) {
		proteins.at(iNeighborProtInProts)->addReaction(currReac); // to neighbor protein
	}
	
}

/* Public Method: initialize()
 * -------------------------------------------------------------------------- 
 * Initialize by filling each dvec with the initial concentrations of the
 * proteins and genes, finding current Dx values, and setting the time = 0.
 */
void Manager::initialize() {
	resizeDVecs();
	time = 0.0;
	for ( int iCell = 0 ; iCell < iTissue.size() ; iCell++ ) {
		for ( int iGene = 0 ; iGene < genes.size() ; iGene++ ) {
			Gene* currGene = genes.at(iGene);
			iTissue.at(iCell)->at(currGene->getISelf()) = currGene->getInitConc();
			currTissue.at(iCell)->at(currGene->getISelf()) = currGene->getInitConc();
		}
		for ( int iProtein = 0 ; iProtein < proteins.size() ; iProtein++ ) {
			Protein* currProtein = proteins.at(iProtein);
			iTissue.at(iCell)->at(currProtein->getISelf()) = currProtein->getInitConc();
			currTissue.at(iCell)->at(currProtein->getISelf()) = currProtein->getInitConc();
		}
	}
	updateDx();
}

/* Public Method: integrate(mode,dt,numSteps)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the currTissue vector based on iTissue, 
 * reactions, and mode.
 */
void Manager::integrate( int numSteps ) {
	for ( int step = 0 ; step < numSteps ; step++ ) {
		
		/* For the User */
		printState();
		
		/* Update currTissue according to current rates */
		/* Our dx vector always should have the appropriate update for the 
		 * next step.
		 */
		for ( unsigned int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
			dvec* currCell = currTissue.at(iCell);
			for ( unsigned int iMol = 0; iMol < currCell->size(); iMol++ ) {
				double update = dx.at(iCell)->at(iMol);
				currCell->at(iMol) += update;
			}
		}
		
		/* Never return a public method without dx vector being up to date */
		updateDx();
		
		/* Progress time */
		time += dt;
	}
}

/* Private Method: updateDx()
 * -------------------------------------------------------------------------- 
 */
void Manager::updateDx() {
	
	/* Reset all rates to zero */
	for ( int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		dx.at(iCell)->zero();
	}
	
	for (unsigned int iCell = 0; iCell < iTissue.size() ; iCell++) {
		/* Update dx vector with new rates, which the Reaction objects calculates
		 * for us.
		 */
		for ( unsigned int iReac = 0 ; iReac < reactions.size() ; iReac++ ) {
		
			/* For this reaction, we go through the participants and update our dxdt
			 * dvecs accordingly.
			 */
			Reaction *currReac = reactions.at(iReac);
			
			/* Have reaction object calculate change */
			currReac->react(currTissue,neighbors,iCell,mode,dt); 
			
			/* Transfer changes to dx vector */
			for ( int partNum = 0 ; partNum < reactions.at(iReac)->getNumPart() ; partNum++ ) {
				double update = currReac->getDx(partNum);
				dx.at(iCell)->at(currReac->getIPart(partNum)) += update;
			}
		}
	}
}

/* Public Method: printState()
 * -------------------------------------------------------------------------- 
 * For the user.
 * Currently goes through each cell in the vector and prints out the 
 * protein concentrations in each cell, (assuming there are 6 of them).
 */
void Manager::printState() {
	std::cout << "time: " << time << std::endl;
	for ( int i = 0 ; i < currTissue.size() ; i++ ) {
		std::cout << "Cell: " << i << std::endl;
		std::cout << currTissue.at(i)->at(0) << " :::: " <<
		currTissue.at(i)->at(1) << " :::: " << currTissue.at(i)->at(2) <<
		" :::: " << currTissue.at(i)->at(3) << " :::: " << currTissue.at(i)->at(4) <<
		" :::: " << currTissue.at(i)->at(5) << std::endl;
	}
	std::cout << std::endl;
}

/* Private Method: addMolecules(num)
 * -------------------------------------------------------------------------- 
 * Replaces dvecs with with new ones incremented in length by num, and
 * updates the numMol variable by incrementing by num.
 *
 * CAUTION: For num < 0, we decrease the number of molecules in the cell.
 * In this case, we assume |num| <= numMol, or else we will have errors.
 */
void Manager::resizeDVecs() {
	
	for ( int iCell = 0 ; iCell < iTissue.size() ; iCell++ ) {
		dvec* currVec = iTissue.at(iCell);
		iTissue.at(iCell) = new dvec(numMol);
		delete currVec;
		currVec = currTissue.at(iCell);
		currTissue.at(iCell) = new dvec(numMol);
		delete currVec;
		currVec = dx.at(iCell);
		dx.at(iCell) = new dvec(numMol);
		delete currVec;
	}
	
}

/* Private Method: updateIndices(firstIndex,numInsertions)
 * -------------------------------------------------------------------------- 
 * Goes through every gene, protein, and reaction and instructs them to 
 * update their internally stored data regarding their own index in the
 * dvecs and well as the index of other compounds they store (like their
 * products, roots, promoters, etc.)
 *
 * firstIndex is the index of the first member of the insertions. A negative
 * numInsertions refers to deletions.
 *
 * CAUTION: we assume firstIndex is within the range of our current number
 * of molecules, and numInsertions, when negative, does not delete more
 * substances than there are to delete. Errors will be produced otherwise.
 *
 * EXAMPLES: 
 *
 *		updateIndices(3,2) produces
 *
 *		  < 1 , 10 , 3 , 7 , 4 , 11 > --> < 1 , 10 , 3 ,___,___, 7 , 4 , 11 >
 *
 *		updateIndices(3,-2) produces
 *
 *		   < 1 , 10 , 3 , 7 , 4 , 11 > --> < 1 , 10 , 3 , 11 >
 *
 *		On the vector above, updateIndices(3,-4) would produce run-time errors.
 */
void Manager::updateIndices( int firstIndex , int numInsertions ) {
	for ( int i = 0 ; i < genes.size() ; i++ ) {
		genes.at(i)->updateIndices(firstIndex,numInsertions);
	}
	for ( int i = 0 ; i < proteins.size() ; i++ ) {
		proteins.at(i)->updateIndices(firstIndex,numInsertions);
	}
	for ( int i = 0 ; i < reactions.size() ; i++ ) {
		reactions.at(i)->updateIndices(firstIndex,numInsertions);
	}
}

