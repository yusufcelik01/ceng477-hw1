#include "parser.h"
#include "raytracer_math.h"
#include <cmath>

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

// float getDeterminant(float arr[3][])//3x3 matrix
// {
//     float det = 0;
//     det += arr[0][0] * arr[1][1] * arr[2][2];
//     det += arr[0][1] * arr[1][2] * arr[2][0];
//     det += arr[0][2] * arr[1][0] * arr[2][1];

//     det -= arr[0][0] * arr[][] * arr[][];
//     det -= arr[0][0] * arr[][] * arr[][];
//     det -= arr[0][0] * arr[][] * arr[][];
    

// }


///////////////////////////////////////
