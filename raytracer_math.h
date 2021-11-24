#ifndef __HW1__VECTOR__OPERATIONS__
#define __HW1__VECTOR__OPERATIONS__

#include "parser.h"

#define VECTOR_LENGTH(v) sqrt(v.x*v.x + v.y*v.y + v.z*v.z)
#define clampLigth()

#define MAX(a,b) ((a>b) ? a : b)
#define MIN(a,b) ((a<b) ? a : b)

parser::Vec3f vectorScalerMult(float c, parser::Vec3f v);

parser::Vec3f vectorSum(parser::Vec3f a, parser::Vec3f b);

parser::Vec3i vectorSum(parser::Vec3i a, parser::Vec3i b);

float dotProduct(parser::Vec3f v, parser::Vec3f u);

int dotProduct(parser::Vec3i v, parser::Vec3i u);

parser::Vec3f crossProduct(parser::Vec3f a, parser::Vec3f b);

parser::Vec3i crossProduct(parser::Vec3i a, parser::Vec3i b);

parser::Vec3f elementViseMultiply(const parser::Vec3f& a, const parser::Vec3f& b);

parser::Vec3f normalize(const parser::Vec3f& v);// normalizes a vector

parser::Vec3f clampColor(const parser::Vec3f& color);

//operator overloads
parser::Vec3f operator+(const parser::Vec3f& a, const parser::Vec3f& b);
parser::Vec3f operator+=(parser::Vec3f& a, const parser::Vec3f& b);
parser::Vec3f operator-(const parser::Vec3f& a, const parser::Vec3f& b);
parser::Vec3f operator-=(parser::Vec3f& a, const parser::Vec3f& b);
parser::Vec3f operator*(const float& a, const parser::Vec3f& b);


std::ostream& operator<<(std::ostream& out, const parser::Vec3f& v);

#endif
