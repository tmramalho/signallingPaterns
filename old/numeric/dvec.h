#ifndef _DVECTOR_H
#define _DVECTOR_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

class dvec {
	
public:
	dvec();
	dvec(unsigned int size);
	dvec(unsigned int size, double value);
	dvec(unsigned int size, double *nums);
	dvec(unsigned int size, dvec *x);
	dvec(dvec *newOne);
	~dvec();
	
	unsigned int size() const;
	void copy(dvec &other);
	void copy(dvec *other);
	double &at(int pos);
	void zero();
	void set(double num);
	void set(const dvec &other);
	
	double *getData();
	
	dvec &operator=(const dvec &other);
	dvec &operator+=(const dvec &rhs);
	dvec &operator-=(const dvec &rhs);
	dvec &operator*=(double rhs);
	dvec operator+(dvec &other);
	dvec operator-(dvec &other);
	dvec operator*(double constant);
	double &operator[] (int pos);
	//double operator* (dvec &other);
	friend std::ostream &operator <<(std::ostream &output, dvec &vector);
	
	void resize( int size );
	
	void equals_lin_comb(const dvec &vec1, double c1);
	void equals_lin_comb(const dvec &vec1, double c1,
						 const dvec &vec2, double c2);
	void equals_lin_comb(const dvec &vec1, double c1,
						 const dvec &vec2, double c2,
						 const dvec &vec3, double c3);
	void equals_lin_comb(const dvec &vec1, double c1,
						 const dvec &vec2, double c2,
						 const dvec &vec3, double c3,
						 const dvec &vec4, double c4);
	
	void plus_equals_lin_comb(const dvec &vec1, double c1);
	void plus_equals_lin_comb(const dvec &vec1, double c1,
							  const dvec &vec2, double c2);
	void plus_equals_lin_comb(const dvec &vec1, double c1,
							  const dvec &vec2, double c2,
							  const dvec &vec3, double c3);
	void plus_equals_lin_comb(const dvec &vec1, double c1,
							  const dvec &vec2, double c2,
							  const dvec &vec3, double c3,
							  const dvec &vec4, double c4);
	
protected:
	double getValue(int pos) const;
	
	unsigned int _size;
	
private:
	double *_main;
};

#endif
