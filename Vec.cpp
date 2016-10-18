#include "Vec.h"

using namespace std;

Vec::Vec(vector<float> data)
{
	if (data.size() != SIZE) {
		//cerr << "Invalid Vector Size:" << data.size() << endl;
		throw exception();
	}
	this->data = data;
}

Vec::Vec()
{
	data = vector<float>(SIZE, 0);
}

Vec & Vec::operator=(const Vec& other)
{
	this->data = other.data;
	return *this;
}

string Vec::toString()
{
	string s;
	for (float f : data) {
		s += to_string(f) + " ";
	}
	return s;
}

Vec Vec::operator+(const Vec & other)
{
	Vec vec;
	for (int i = 0; i <SIZE; i++) {
		vec.data[i] = this->data[i] + this->data[i];
	}
	return vec;
}

Vec Vec::operator/(float other)
{
	Vec vec;
	for (int i = 0; i < SIZE; i++) {
		vec.data[i] /= other;
	}
	return vec;
}

void Vec::operator+=(const Vec & other)
{
	for (int i = 0; i < SIZE; ++i) {
		this->data[i] += other.data[i];
	}
}

void Vec::operator/=(float other)
{
	for (int i = 0; i < SIZE; i++) {
		data[i] /= other;
	}
}

double Vec::dot(const Vec& other) const{
    double r = 0;
    for(int i = 0; i < SIZE; ++i){
        r += this->data[i] * other.data[i];
    }
    return r;
}

double Vec::magnitude() const{
    double r = 0;
    for(float d : data)
        r += d * d;
    return sqrt(r);
}

double Vec::cosDist(const Vec& other) const{
    return 1 - (this->dot(other) / (this->magnitude() * other.magnitude()));
}

float Vec::get(int i) const{
    return this->data[i];
}