#include "dvec.h"

dvec::dvec() {
	//std::cout << "dvec default constructor\n";
	_main = new double[1];
	_size = 1;
	_main[0] = 0;
}

dvec::dvec(unsigned int size) {
	//std::cout << "dvec size constructor ( " << size << " )\n";
	_main = new double[size];
	_size = size;
	for(unsigned int i = 0 ; i < size ; i++) {
		_main[i] = 0;
	}
	
}

dvec::dvec(unsigned int size, double value) {
	//std::cout << "dvec size value constructor ( " << size << " )\n";
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
	//std::cout << "dvec copy constructor\n";
	_size = newOne->size();
	_main = new double[_size];
	//_main = new double[size];
	for(unsigned int i=0;i<_size;i++) {
		_main[i] = (*newOne)[i];
	}
	
}

dvec::~dvec() {
	std::cout << "dvec destructor\n"; /*<< _size << " : ";
	for ( int i = 0 ; i < _size ; i++ ) {
		std::cout << _main[i] << " , ";
	}
	std::cout << ")\n";*/
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
	if(pos >= 0 && (unsigned int)pos < _size) return _main[pos];
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
	//std::cout << "dvec operator=\n";
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
/*
double dvec::operator*(dvec &other) {
	double output = 0;
	for(unsigned int i = 0 ; i < _size ; i++) {
		output += _main[i] * other[i];
	}
	return output;
}
*/
std::ostream &operator <<(std::ostream &output, dvec &vector) {
	for(unsigned int i=0;i<vector._size;i++) {
		output << vector._main[i];
		if(i<(vector._size-1)) output << ", ";
	}
	return output;
}

void dvec::resize(int size) {
	if (size == _size) return;
	else {
		_size = size;
		delete[] _main;
		_main = new double[_size];
	}
}

/* CAUTION: MEMORY ERROR IF DVEC NOT ALREADY ALLOCATED THE RIGHT SIZE */
void dvec::equals_lin_comb(const dvec &vec1, double c1) {
	
	if (vec1.size() != _size) {
		std::cout << "Size Error: equals_lin_comb 1\n";
	}
	
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = vec1.getValue(i)*c1;
	}
}

void dvec::equals_lin_comb(const dvec &vec1, double c1,
						   const dvec &vec2, double c2) {
	
	if (vec1.size() != _size) {
		std::cout << "Size Error: equals_lin_comb 2\n";
	}
	
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = 
		vec1.getValue(i)*c1 + 
		vec2.getValue(i)*c2;
	}
}

void dvec::equals_lin_comb(const dvec &vec1, double c1,
						   const dvec &vec2, double c2,
						   const dvec &vec3, double c3) {
	if (vec1.size() != _size) {
		std::cout << "Size Error: equals_lin_comb 3\n";
	}
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = 
		vec1.getValue(i)*c1 + 
		vec2.getValue(i)*c2 +
		vec3.getValue(i)*c3;
	}
}

void dvec::equals_lin_comb(const dvec &vec1, double c1,
						   const dvec &vec2, double c2,
						   const dvec &vec3, double c3,
						   const dvec &vec4, double c4) {
	if (vec1.size() != _size) {
		std::cout << "Size Error: equals_lin_comb 4\n";
	}
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] = 
		vec1.getValue(i)*c1 + 
		vec2.getValue(i)*c2 +
		vec3.getValue(i)*c3 +
		vec4.getValue(i)*c4;
	}
}

void dvec::plus_equals_lin_comb(const dvec &vec1, double c1) {
	if (vec1.size() != _size) {
		std::cout << "Size Error: plus_equals_lin_comb 1\n";
	}
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] += vec1.getValue(i)*c1;
	}
}

void dvec::plus_equals_lin_comb(const dvec &vec1, double c1,
								const dvec &vec2, double c2) {
	if (vec1.size() != _size) {
		std::cout << "Size Error: plus_equals_lin_comb 2\n";
	}
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] += 
		(vec1.getValue(i)*c1 +
		 vec2.getValue(i)*c2);
	}
}

void dvec::plus_equals_lin_comb(const dvec &vec1, double c1,
								const dvec &vec2, double c2,
								const dvec &vec3, double c3) {
	if (vec1.size() != _size) {
		std::cout << "Size Error: plus_equals_lin_comb 3\n";
	}
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] += 
		(vec1.getValue(i)*c1 +
		 vec2.getValue(i)*c2 +
		 vec3.getValue(i)*c3);
	}
}

void dvec::plus_equals_lin_comb(const dvec &vec1, double c1,
								const dvec &vec2, double c2,
								const dvec &vec3, double c3,
								const dvec &vec4, double c4) {
	if (vec1.size() != _size) {
		std::cout << "Size Error: plus_equals_lin_comb 4\n";
	}
	for(unsigned int i = 0 ; i < _size ; i++) {
		_main[i] += 
		(vec1.getValue(i)*c1 +
		 vec2.getValue(i)*c2 +
		 vec3.getValue(i)*c3 +
		 vec4.getValue(i)*c4);
	}
}




