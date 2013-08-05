/*
 *  HillPromReaction.h
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

# include "Reaction.h"

#ifndef HILLPROMREACTION_H
#define HILLPROMREACTION_H

class HillPromReaction : public Reaction {

public:
	HillPromReaction();
	HillPromReaction( int i_promoter , int i_promoted ,
					 double kinetic , double K , double cooperativity );
	~HillPromReaction();
	
	virtual int get_i_part( int part_num );
	
	virtual void react( dmat& curr_tissue , dmat& dx_dt ,
					   std::vector< std::vector<int>* >& neighbors ,
					   int i_curr_cell );
	
	virtual void update_indices( int first_index , int num_insertion );
	
private:
	
	int _i_promoter;
	int _i_promoted;
	
	double _kinetic;
	double _K;
	double _cooperativity;

};

#endif

