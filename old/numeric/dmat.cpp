/*
 *  dmat.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/25/13.
 *
 */

#include "dmat.h"

dmat::dmat() {
	_size = 1;
	_main = new double[1];
	_main[0] = 0;
	_num_row = 1;
	_num_col = 1;
}

dmat::dmat( unsigned int num_row , unsigned int num_col ) {
	_size = num_row*num_col;
	_main = new double[_size];
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = 0;
	}
	_num_row = num_row;
	_num_col = num_col;
}

dmat::~dmat() {}

double &dmat::at( int row , int col ) {
	return dvec::at( row * _num_col + col );
}

void dmat::resize( unsigned int num_row , unsigned int num_col ) {
	delete[] _main;
	_size = num_row*num_col;
	_main = new double[_size];
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = 0;
	}
	_num_row = num_row;
	_num_col = num_col;
}


