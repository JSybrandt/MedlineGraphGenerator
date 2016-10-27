#pragma once

#include <vector>
#include<exception>
#include<string>
#include<iostream>

#include<cmath>

using std::vector;
using std::string;




class Vec {
private:
	vector<float> data;
        const int SIZE = 500;

public:
	Vec(vector<float> data);
	Vec();
	Vec(const Vec&) = default;
	Vec& operator=(const Vec&);
	std::string toString();
	Vec operator+(const Vec& other);
	Vec operator/(float other);
	void operator+=(const Vec& other);
	void operator/=(float other);
        double dot(const Vec& other) const;
        double magnitude() const;
        double cosDist(const Vec& other) const;
        float get(int i) const;
};
