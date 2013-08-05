/*
 *  HillDegReaction.h
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

# include "Reaction.h"

#ifndef HILLREPREACTION_H
#define HILLREPREACTION_H

class HillRepReaction : public Reaction {
	
public:
	HillRepReaction();
	HillRepReaction( int i_repressor , int i_repressed ,
					 double kinetic , double K , double cooperativity );
	~HillRepReaction();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt ,
					   std::vector< std::vector<int>* >& neighbors ,
					   int i_curr_cell );
	
	virtual void update_indices( int first_index , int num_insertion );
	
private:
	
	int _i_repressor;
	int _i_repressed;
	
	double _kinetic;
	double _K;
	double _cooperativity;
	
};

#endif