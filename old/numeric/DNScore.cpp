/*
 *  DNScore.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 8/20/13.
 *
 */

#include "DNScore.h"

double DNScore::run(Manager* x,boost::random::mt19937& generator ) {
	x->initialize();
	x->integrate(10000,generator);
	dmat final_state = x->get_curr_state();
	
	double score = 0.0;
	
	int num_cell = x->get_num_cell();
	
	for ( int i = 0 ; i < num_cell-1 ; i++ ) {
		score -= abs( final_state.at(i,2)-final_state.at(i+1,2));
		score += exp(final_state.at(i,2)/2);
	}
	
	score += exp(final_state.at(num_cell-1,2));
	return score;
}