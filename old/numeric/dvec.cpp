#include "dvec.h"

dvec::dvec() {
	_main = new double[1];
	_size = 1;
	_main[0] = 0;
}

dvec::dvec(unsigned int size) {
	_main = new double[size];
	_size = size;
	for(unsigned int i = 0 ; i < size ; i++) {
		_main[i] = 0;
	}
	
}

dvec::dvec(unsigned int size, double value) {
	_main = new double[size];
	_size = size;
	for(unsigned int i = 0 ; i < size ; i++) {
		_main[i] = value;
	}
	
}

dvec::dvec(unsigned int size, double *nums) {
	_main = new double[size];
	_size = size;
	for(unsigned int i=0;i<size;i++) {
		_main[i] = nums[i];
	}
	
}

dvec::dvec(unsigned int size, dvec *x) {
	_size = x->size() > size ? size : x->size();
	_main = new double[_size];
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = (*x)[i];
	}
	
}

/*dvec::dvec(dvec &newOne) {
 _size = newOne.size();
 _main = new double[_size];
 //_main = new double[size];
 for(unsigned int i=0;i<_size;i++) {
 _main[i] = newOne[i];
 }
 
 }*/

dvec::dvec(dvec *newOne) {
	_size = newOne->size();
	_main = new double[_size];
	//_main = new double[size];
	for(unsigned int i=0;i<_size;i++) {
		_main[i] = (*newOne)[i];
	}
	
}

dvec::~dvec() {
	delete[] _main;
}

unsigned int dvec::size() const {
	return _size;
}

double *dvec::getData() {
	return _main;
}

void dvec::zero() {
	memset(_main, 0, _size * sizeof(double));
}

void dvec::set(double num) {
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = num;
	}
}

void dvec::set(const dvec &other) {
	// Only do assignment if RHS is a different object from this.
	if (this != &other) {
		delete[] _main;
		unsigned int size = other.size();
		_main = new double[size];
		_size = size;
		for(unsigned int i = 0 ; i < _size ; i++) {
			_main[i] = other.getValue(i);
		}
    }
}

void dvec::copy(dvec &other) {
	unsigned int size = _size < other.size() ? _size : other.size();
	for(unsigned int i = 0 ; i < size ; i++)
		_main[i] = other[i];
}

void dvec::copy(dvec *other) {
	unsigned int size = _size < other->size() ? _size : other->size();
	for(unsigned int i = 0 ; i < size ; i++)
		_main[i] = other->at(i);
}

double &dvec::at(int pos) {
	if(pos > 0 && (unsigned int)pos <= _size) return _main[pos];
	else return _main[0];
}

double &dvec::operator[](int pos) {
	if(pos > 0 && (unsigned int)pos <= _size) return _main[pos];
	else return _main[0];
}

double dvec::getValue(int pos) const {
	if(pos > 0 && (unsigned int)pos <= _size) return _main[pos];
	else return _main[0];
}

dvec &dvec::operator=(const dvec &other) {
	// Only do assignment if RHS is a different object from this.
	if (this != &other) {
		delete[] _main;
		unsigned int size = other.size();
		_main = new double[size];
		_size = size;
		for(unsigned int i = 0 ; i < _size ; i++) {
			_main[i] = other.getValue(i);
		}
    }
	
	return *this;
}

dvec &dvec::operator+=(const dvec &other) {
	if (this != &other) {
		unsigned int size = other.size();
		if(_size <= size) {
			for(unsigned int i = 0 ; i < _size ; i++) {
				_main[i] += other.getValue(i);
			}
		} else {
			for(unsigned int i = 0 ; i < size ; i++) {
				_main[i] += other.getValue(i);
			}
		}
	}
	
	return *this;
}

dvec &dvec::operator-=(const dvec &other) {
	if (this != &other) {
		unsigned int size = other.size();
		if(_size <= size) {
			for(unsigned int i = 0 ; i < _size ; i++) {
				_main[i] -= other.getValue(i);
			}
		} else {
			for(unsigned int i = 0 ; i < size ; i++) {
				_main[i] -= other.getValue(i);
			}
		}
	}
	
	return *this;
}

dvec &dvec::operator*=(double rhs) {
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] *= rhs;
	}
	
	return *this;
}

dvec dvec::operator+(dvec &other) {
	dvec result = *this;
	result += other;
	return result;
}

dvec dvec::operator-(dvec &other) {
	dvec result = *this;
	result -= other;
	return result;
}

dvec dvec::operator*(double constant) {
	dvec result = *this;
	result *= constant;
	return result;
}

double dvec::operator*(dvec &other) {
	double output = 0;
	for(unsigned int i = 0 ; i < _size ; i++) {
		output += _main[i] * other[i];
	}
	return output;
}

std::ostream &operator <<(std::ostream &output, dvec &vector) {
	for(unsigned int i=0;i<vector._size;i++) {
		output << vector._main[i];
		if(i<(vector._size-1)) output << ", ";
	}
	return output;
}
