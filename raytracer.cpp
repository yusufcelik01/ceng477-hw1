#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <cmath>
#include "raytracer_math.h"

typedef unsigned char RGB[3];

enum IntersecType {none, sphere, triangle, mash};

int RayIntersecObj(const parser::Scene &scene,parser::Ray ray, IntersecType &i_type);
parser::Intersection intersectRaySphere(const parser::Scene &scene, parser::Ray &eye_ray, int sphere_id);

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    int i,j,k,bar;
    int width, height;
    int cam_id; 
    int num_of_cameras;

    
    parser::Vec3f temp_vec;
    float temp;

    num_of_cameras = scene.cameras.size();

    for(cam_id = 0; cam_id < num_of_cameras; cam_id++)
    {
        parser::Camera cam;
        parser::Vec3f e;//origin
        parser::Vec3f w,v,u;//coordinate system

        parser::Vec3f ray;
        parser::Vec3f q;// top corner coordinate
        parser::Vec3f s;// pixel coordinate
        parser::Vec3f d;// d in the eq r(t)= e+dt
        float pixel_width, pixel_height;
        

        cam = scene.cameras[cam_id];

        width = cam.image_width;
        height = cam.image_height;


        unsigned char* image = new unsigned char [width * height * 3];


        e = cam.position;
        w = vectorScalerMult(-1.0, cam.gaze);
        v = cam.up;
        u = crossProduct(v, w);

        /*
        std::cout << "e = " << e.x << "," << e.y << "," << e.z << std::endl
                  << "w = " << w.x << "," << w.y << "," << w.z << std::endl
                  << "v = " << v.x << "," << v.y << "," << v.z << std::endl
                  << "u = " << u.x << "," << u.y << "," << u.z << std::endl;

        std::cout << "\nNear plane things: " << cam.near_plane.x << " "
                                            << cam.near_plane.y << " "
                                            << cam.near_plane.z << " "
                                            << cam.near_plane.w << std::endl;

        */
        pixel_width = (cam.near_plane.y - cam.near_plane.x) / cam.image_width;
        pixel_height = (cam.near_plane.w - cam.near_plane.z) / cam.image_height;

        
        temp_vec = vectorScalerMult(cam.near_plane.x, u);//l.u

        q = vectorSum(cam.position, temp_vec);//e + l.u
        
        temp_vec = vectorScalerMult(cam.near_plane.w, v);//t.v
        q = vectorSum(q, temp_vec); // += t.v

        temp_vec = vectorScalerMult(cam.near_distance, cam.gaze);
        q = vectorSum(q, temp_vec);
        //q has the position of top left corner

        /*
        std::cout << "pixel width " << pixel_height << std::endl;
        std::cout << "pixel height " << pixel_height << std::endl;

        std::cout << "q = " << q.x << "," << q.y << "," << q.z << std::endl;
        */

        i = 0; //pixels' color value
        float p_height = pixel_height*0.5;
        float p_width = pixel_width*0.5;
        for(int y=0; y < cam.image_height; y++){
            for(int x=0; x < cam.image_width; x++){
                //first compute pixel coordinates;
                temp_vec = vectorScalerMult(x*pixel_width+p_width, u);

                s = vectorSum(q, temp_vec);//s = q + s_u.u

                temp_vec = vectorScalerMult( -1.0 * (y*pixel_height + p_height), v);// - s_v . v

                s = vectorSum(s, temp_vec);
                //s = q + s_u.u - s_v.v

                d = vectorSum(s, vectorScalerMult(-1.0, e));//s-e

                parser::Ray ray;
                ray.e = e;
                ray.d = d;

                //float intersect.t1, intersect.t2;//different solutions of the equation

                //calculate spheres' closest

                IntersecType i_type = none;//whether it intersects or not
                int closest_obj_id = RayIntersecObj(scene,ray,i_type);

                //find the colour of that material

                if(!i_type){
                    // does not intersect get backround colour
                    parser::Vec3i bg;//backgroung
                    bg = scene.background_color;
                    image[i++] = bg.x;
                    image[i++] = bg.y;
                    image[i++] = bg.z;
                }
                else{
                    // if there is a sphere just write green
                    parser::Material material = scene.materials[scene.spheres[closest_obj_id].material_id - 1];
                    image[i++] = 255 * material.diffuse.x;//R
                    image[i++] = 255 * material.diffuse.y;//G
                    image[i++] = 255 * material.diffuse.z;//B
                }
                //write that colour to the array

            }
        }
        //print to ppm

        //write to file
        write_ppm(scene.cameras[cam_id].image_name.c_str(), image, width, height);

    }
    //write_ppm("test.ppm", image, width, height);
}

int RayIntersecObj(const parser::Scene &scene,parser::Ray ray, IntersecType &i_type){
    float t_min;//closest objects parameter
    int closest_obj_id;
    int numOfSpheres = scene.spheres.size();
    for(int k = 0; k < numOfSpheres; k++){
        parser::Intersection intersect = intersectRaySphere(scene,ray, k);

        if(intersect.discriminant >= 0){
            //meaning they intersect
            if(i_type==none){
                i_type = sphere;
                //assign the smallest
                if(intersect.t1 < intersect.t2)
                    t_min = intersect.t1;
                else
                    t_min = intersect.t2;
                closest_obj_id = k;
            }
            else{
                if(intersect.t1 < t_min){
                    t_min = intersect.t1;
                    closest_obj_id = k;
                }
                if(intersect.t2 < t_min);{
                    t_min = intersect.t2;
                    closest_obj_id = k;
                }
            }
        }
    }
    return closest_obj_id;
}

parser::Intersection intersectRaySphere(const parser::Scene &scene, parser::Ray &eye_ray, int sphere_id){
    parser::Intersection intersect;

    parser::Vec3f c;
    float radius;
    float temp;

    c = scene.vertex_data[scene.spheres[sphere_id].center_vertex_id -1 ];
    radius = scene.spheres[sphere_id].radius;


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