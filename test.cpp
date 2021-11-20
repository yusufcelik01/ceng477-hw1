#include<iostream>

#include "raytracer_math.h"

int main()
{
    parser::Vec3f f1,f2 ,fr;
    parser::Vec3i i1,i2 ,ir;
    float ff;

    f1.x = 1.5; f1.y = 0; f1.z = 0;
    f2.x = 0  ; f2.y = 2; f2.z = 0;
    i1.x = 1  ; i1.y = 0; i1.z = 0;
    i2.x = 0  ; i2.y = 1; i2.z = 0;

    fr = vectorScalerMult(3.5,f2);
    std::cout << fr.x <<" " <<fr.y <<" " << fr.z  << std::endl;

    fr = vectorSum(f1,f2);
    std::cout << fr.x <<" " <<fr.y <<" " << fr.z  << std::endl;

    fr = crossProduct(f2,f1);
    std::cout <<"crossProduct "<< fr.x <<" " <<fr.y <<" " << fr.z  << std::endl;

    ff = dotProduct(f1,f1);
    std::cout << "dotProduct " << ff << std::endl;
    return 0;
}
