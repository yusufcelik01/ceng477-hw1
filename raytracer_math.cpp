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

float getDeterminant(float arr[3][])//3x3 matrix
{
    float det = 0;
    det += arr[0][0] * arr[1][1] * arr[2][2];
    det += arr[0][1] * arr[1][2] * arr[2][0];
    det += arr[0][2] * arr[1][0] * arr[2][1];

    det -= arr[0][0] * arr[][] * arr[][];
    det -= arr[0][0] * arr[][] * arr[][];
    det -= arr[0][0] * arr[][] * arr[][];
    

}


///////////////////////////////////////

parser::Intersection parser::Scene::intersectRaySphere(Ray eye_ray, int sphere_id)
{
    Intersection intersect;

    parser::Vec3f c;
    float radius;
    float temp;

    c = vertex_data[spheres[sphere_id].center_vertex_id -1 ];
    radius = spheres[sphere_id].radius;


    parser::Vec3f e_c; //e-c and d^2 is freq used so assign it to a variable
    float d_sqr;

    e_c = vectorSum(eye_ray.e, vectorScalerMult(-1.0, c));
    d_sqr = dotProduct(eye_ray.d, eye_ray.d);


    temp = dotProduct(eye_ray.d, e_c);
    temp *= temp;//(d.(e-c))^2

    intersect.discriminant = temp;
    intersect.discriminant -= d_sqr * (dotProduct(e_c, e_c) - radius*radius);

    if(intersect.discriminant >= 0)
    {
        intersect.t1 = -dotProduct(eye_ray.d, e_c);
        intersect.t2 = intersect.t1;

        intersect.t1 += sqrt(intersect.discriminant);//-b + sqrt delta
        intersect.t1 /= d_sqr;
        intersect.t2 -= sqrt(intersect.discriminant);//-b - sqrt delta
        intersect.t2 /= d_sqr;

        return intersect;
    }
    else
    {
        intersect.t1 = -1;
        intersect.t2 = -1;

        return intersect;
    }


}
