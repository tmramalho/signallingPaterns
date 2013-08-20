/*
 *  dmat.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/25/13.
 *
 */

#include "dmat.h"

dmat::dmat() {
	//std::cout << "dmat default constructor\n";
	/* dvec() will initialize to dvec of size 1 */
	_num_row = 1;
	_num_col = 1;
}

dmat::dmat(dmat* newOne) : dvec(newOne) {
	//std::cout << "dmat copy constructor\n";
	_num_row = newOne->_num_row;
	_num_col = newOne->_num_col;
}

dmat::dmat( unsigned int num_row , unsigned int num_col ) : dvec(num_row*num_col) {
	//std::cout << "dmat size constructor ( " << num_row << " , " << num_col << " )\n";
	_num_row = num_row;
	_num_col = num_col;
}

dmat::~dmat() {
	std::cout << "dmat destructor\n";
}

double &dmat::at( int row , int col ) {
	return dvec::at( row * _num_col + col );
}

dmat &dmat::operator=(const dmat &other) {
	/* dvec assignment operator called automatically */ 
	//std::cout << "dmat operator=\n";
	dvec::operator=(static_cast<dvec const&>(other));
	_num_row = other._num_row;
	_num_col = other._num_col;
	return *this;
}

int dmat::get_num_row() const {
	return _num_row;
}

int dmat::get_num_col() const {
	return _num_col;
}

void dmat::resize( unsigned int num_row , unsigned int num_col ) {
	dvec::resize(num_row*num_col);
	_num_row = num_row;
	_num_col = num_col;
}


