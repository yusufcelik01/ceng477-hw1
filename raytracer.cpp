#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <cmath>
#include "raytracer_math.h"
#include <limits>

typedef unsigned char RGB[3];

enum IntersecType {none, sphere, triangle, mesh};

int RayIntersecObj(const parser::Scene &scene,parser::Ray ray, IntersecType &i_type);
parser::Intersection intersectRaySphere(const parser::Scene &scene, parser::Ray &eye_ray, int sphere_id);
parser::Vec3f intersectRayFace(const parser::Scene &scene, parser::Ray &eye_ray, const parser::Face &face);
float Determinant();

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
        float pixel_width, pixel_height, p_height, p_width;
        

        cam = scene.cameras[cam_id];

        width = cam.image_width;
        height = cam.image_height;


        unsigned char* image = new unsigned char [width * height * 3];


        e = cam.position;
        w = vectorScalerMult(-1.0, cam.gaze);
        v = cam.up;
        u = crossProduct(v, w);

        pixel_width = (cam.near_plane.y - cam.near_plane.x) / cam.image_width;
        pixel_height = (cam.near_plane.w - cam.near_plane.z) / cam.image_height;
        
        temp_vec = vectorScalerMult(cam.near_plane.x, u);//l.u

        q = vectorSum(cam.position, temp_vec);//e + l.u
        
        temp_vec = vectorScalerMult(cam.near_plane.w, v);//t.v
        q = vectorSum(q, temp_vec); // += t.v

        temp_vec = vectorScalerMult(cam.near_distance, cam.gaze);
        q = vectorSum(q, temp_vec);
        //q has the position of top left corner

        i = 0; //pixels' color value
        p_height = pixel_height*0.5;
        p_width = pixel_width*0.5;
        for(int y=0; y < cam.image_height; y++){
            for(int x=0; x < cam.image_width; x++){
                //first compute pixel coordinates;
                parser::Ray ray;
                IntersecType i_type = none;//whether it intersects or not
                parser::Material material;
                int closest_obj_id = -1;
                temp_vec = vectorScalerMult(x*pixel_width+p_width, u);

                s = vectorSum(q, temp_vec);//s = q + s_u.u

                temp_vec = vectorScalerMult( -1.0 * (y*pixel_height + p_height), v);// - s_v . v

                s = vectorSum(s, temp_vec);
                //s = q + s_u.u - s_v.v

                d = vectorSum(s, vectorScalerMult(-1.0, e));//s-e

                ray.e = e;
                ray.d = d;

                //float intersect.t1, intersect.t2;//different solutions of the equation

                //calculate spheres' closest

                closest_obj_id = RayIntersecObj(scene,ray,i_type);

                //find the colour of that material
                switch (i_type){
                case none:
                    parser::Vec3i bg;//backgroung
                    bg = scene.background_color;
                    image[i++] = bg.x;
                    image[i++] = bg.y;
                    image[i++] = bg.z;
                    break;
                case sphere:
                    material = scene.materials[scene.spheres[closest_obj_id].material_id - 1];
                    image[i++] = 255 * material.diffuse.x;//R
                    image[i++] = 255 * material.diffuse.y;//G
                    image[i++] = 255 * material.diffuse.z;//B
                    break;
                case triangle:
                    material = scene.materials[scene.triangles[closest_obj_id].material_id - 1];
                    image[i++] = 255 * material.diffuse.x;//R
                    image[i++] = 255 * material.diffuse.y;//G
                    image[i++] = 255 * material.diffuse.z;//B
                    break;
                case mesh:
                    material = scene.materials[scene.meshes[closest_obj_id].material_id - 1];
                    image[i++] = 255 * material.diffuse.x;//R
                    image[i++] = 255 * material.diffuse.y;//G
                    image[i++] = 255 * material.diffuse.z;//B
                    break;
                // default:
                //     break;
                }
            }
        }
        //print to ppm

        //write to file
        write_ppm(scene.cameras[cam_id].image_name.c_str(), image, width, height);

    }
    //write_ppm("test.ppm", image, width, height);
}
// TODO
// add face id in return along as mesh id
int RayIntersecObj(const parser::Scene &scene,parser::Ray ray, IntersecType &i_type){
    float t_min = __FLT_MAX__;//closest objects parameter
    int closest_obj_id = -1;
    int numOfSpheres = scene.spheres.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();
    for(int i = 0; i < numOfSpheres; i++){
        parser::Intersection intersect = intersectRaySphere(scene,ray, i);

        if(intersect.discriminant >= 0){
            //meaning they intersect
            if(intersect.t1 < t_min){
                t_min = intersect.t1;
                closest_obj_id = i;
                i_type = sphere;
            }
            if(intersect.t2 < t_min);{
                t_min = intersect.t2;
                closest_obj_id = i;
                i_type = sphere;
            }
        }
    }
    for(int i = 0; i<numOfTriangles;i++){
        parser::Vec3f b_g_t = intersectRayFace(scene,ray,scene.triangles[i].indices);
        if(b_g_t.z < t_min && b_g_t.x >= 0 && b_g_t.y >= 0){
            if(b_g_t.y <= 1 && b_g_t.x <= (1 - b_g_t.y)){
                t_min = b_g_t.z;
                closest_obj_id = i;
                i_type = triangle;
            }
        }
    }
    for (size_t j = 0; j < numOfMeshes; j++){
        parser::Mesh currmesh = scene.meshes[j];
        int numOfFaces = currmesh.faces.size();
        for(int i = 0; i<numOfFaces;i++){
            parser::Vec3f b_g_t = intersectRayFace(scene,ray,currmesh.faces[i]);
            if(b_g_t.z < t_min && b_g_t.x >= 0 && b_g_t.y >= 0){
                if(b_g_t.y <= 1 && b_g_t.x <= (1 - b_g_t.y)){
                    t_min = b_g_t.z;
                    closest_obj_id = j;
                    i_type = mesh;
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

float Determinant(float matrix[][3]){
    float det = 0;
    det += matrix[0][0] * matrix[1][1] * matrix[2][2];
    det += matrix[0][1] * matrix[1][2] * matrix[2][0];
    det += matrix[0][2] * matrix[1][0] * matrix[2][1];

    det -= matrix[2][0] * matrix[1][1] * matrix[0][2];
    det -= matrix[2][1] * matrix[1][2] * matrix[0][0];
    det -= matrix[2][2] * matrix[1][0] * matrix[0][1];
    return det;
}

parser::Vec3f intersectRayFace(const parser::Scene &scene, parser::Ray &eye_ray, const parser::Face &face){
    parser::Vec3f A,B,C,res;
    float a,b,c,d,e,f,g,h,i,j,k,l;
    float ei_hf,gf_di,dh_eg,ak_jb,jc_al,bl_kc;
    float M = -1;
    // if do not intersec
    res.x = -1; 
    res.y = -1;
    res.z = __FLT_MAX__;

    A = scene.vertex_data[face.v0_id-1];
    B = scene.vertex_data[face.v1_id-1];
    C = scene.vertex_data[face.v2_id-1];

    a = A.x - B.x;
    b = A.y - B.y;
    c = A.z - B.z;
    d = A.x - C.x;
    e = A.y - C.y;
    f = A.z - C.z;
    g = eye_ray.d.x;
    h = eye_ray.d.y;
    i = eye_ray.d.z;
    j = A.x - eye_ray.e.x;
    k = A.y - eye_ray.e.y;
    l = A.z - eye_ray.e.z;

    gf_di = g*f - d*i;
    dh_eg = d*h - e*g;
    ei_hf = e*i - h*f;
    jc_al = j*c - a*l;
    bl_kc = b*l - k*c;
    ak_jb = a*k - j*b;
    // TODO
    // try parallization

    M = a*ei_hf + b*gf_di + c*dh_eg;
    if(M==0.0)
        return res;
    res.x = (j*ei_hf + k*gf_di + l*dh_eg)/M; // beta
    res.y = (i*ak_jb + h*jc_al + g*bl_kc)/M; // gama
    res.z = -(f*ak_jb + e*jc_al + d*bl_kc)/M; // t
    return res;
}
