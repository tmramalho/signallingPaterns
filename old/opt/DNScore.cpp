/*
 *  DNScore.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 8/20/13.
 *
 */

#include "DNScore.h"

/* See Documentation.txt for information on the DNScore objective function */

double DNScore::run(Manager* x,boost::random::mt19937& generator ) {
	double score = 0.0;
	x->initialize();
	x->integrate(10000,generator);
	dmat final_state;
	x->get_curr_state(final_state);
	
	int num_cell = x->get_num_cell();
	
	int next_sign = 1;
	if (final_state.at(0,x->get_num_gene()) < final_state.at(1,x->get_num_gene())) {
		next_sign = -1;
	}
	
	double tot1 = 0.0;
	double tot2 = 0.0;
	double avg1, avg2;
	double score_update = 0.0;
	for ( int i = 0 ; i < num_cell-1 ; i++ ) {
		if ( final_state.at(i,x->get_num_gene()) < 0 || final_state.at(i,x->get_num_gene()+1) < 0 ) score += 100;
		else {
			score_update = 5*next_sign*(final_state.at(i+1,x->get_num_gene())-final_state.at(i,2))/std::max(final_state.at(i,x->get_num_gene()), final_state.at(i+1,2));
			score += score_update;
			if ( score_update < 0 ) next_sign *= -1; // concentrations change in expected direction, meaning we want the next one to switch.
			tot1 += final_state.at(i,x->get_num_gene());
			tot2 += final_state.at(i,x->get_num_gene()+1);
		}
	}
	if ( final_state.at(num_cell-1,x->get_num_gene()) < 0 || final_state.at(num_cell-1,x->get_num_gene()+1) < 0 ) score += 100;
	else {
		tot1 += final_state.at(num_cell-1,x->get_num_gene());
		tot2 += final_state.at(num_cell-1,x->get_num_gene()+1);
	}
	
	avg1 = tot1/(x->get_num_cell());
	avg2 = tot2/(x->get_num_cell());
	if ( avg1 < .2 || avg2 < .2 || avg1 > 2 || avg2 > 2) score += 50.0;
	if ( avg1 < .1 || avg2 < .1 || avg1 > 3 || avg2 > 3) score += 50.0;
	score += 1/(8*avg1) + avg1/2;
	return score;
}