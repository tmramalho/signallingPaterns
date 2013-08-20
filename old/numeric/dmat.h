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
	dmat(dmat* newOne);
	dmat( unsigned int num_row , unsigned int num_col );
	~dmat();
	
	double &at(int row, int col);
	dmat &operator=(const dmat &other);
	
	int get_num_row() const;
	int get_num_col() const;
	
	void resize( unsigned int num_row , unsigned int num_col );
	
	
private:
	
	int _num_row;
	int _num_col;

};

# endif