/*
 *  dmat.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/25/13.
 *
 */

# ifndef DMAT_H
# define DMAT_H

# include "dvec.h"

class dmat : public dvec {

public:
	dmat();
	dmat( unsigned int num_row , unsigned int num_col );
	~dmat();
	
	double &at(int row, int col);
	
	void resize( unsigned int num_row , unsigned int num_col );
	
	
private:
	
	int _num_row;
	int _num_col;

};

# endif