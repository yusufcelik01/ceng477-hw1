#ifndef __HW1__VECTOR__OPERATIONS__
#define __HW1__VECTOR__OPERATIONS__

#include "parser.h"

parser::Vec3f vectorScalerMult(float c, parser::Vec3f v)
{
    v.x *= c;
    v.y *= c;
    v.z *= c;

    return v;
}

//vector summation
parser::Vec3f vectorSum(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result;

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    
    return result;
}

parser::Vec3i vectorSum(parser::Vec3i a, parser::Vec3i b)
{
    parser::Vec3i result;

    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    
    return result;
}
//dot product

float dotProduct(parser::Vec3f v, parser::Vec3f u)
{
    float sum = 0;

    sum += v.x * u.x;
    sum += v.y * u.y;
    sum += v.z * u.z;

    return sum;
}

int dotProduct(parser::Vec3i v, parser::Vec3i u)
{
    int sum = 0;

    sum += v.x * u.x;
    sum += v.y * u.y;
    sum += v.z * u.z;

    return sum;
}


//cross product

parser::Vec3f crossProduct(parser::Vec3f a, parser::Vec3f b)
{
    parser::Vec3f result;

    result.x = a.y*b.z - a.z*b.y;
    result.y = a.z*b.x - a.x*b.z;
    result.z = a.x*b.y - a.y*b.x;

    return result;
}
parser::Vec3i crossProduct(parser::Vec3i a, parser::Vec3i b)
{
    parser::Vec3i result;

    result.x = a.y*b.z - a.z*b.y;
    result.y = a.z*b.x - a.x*b.z;
    result.z = a.x*b.y - a.y*b.x;

    return result;
}
#endif
