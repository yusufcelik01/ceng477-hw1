#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <cmath>
#include "raytracer_math.h"
#include <limits>
#include <thread>

typedef unsigned char RGB[3];
static const int Num_Th = 8;

enum IntersecType {none, sphere, triangle, mesh};

typedef struct 
{
    IntersecType obj_type;
    float t;
    int obj_id;
    int face_id;//if it is a mesh we need to find which face it is

}IntersectionData;

struct ARGS
{
    parser::Camera cam;
    parser::Vec3f e;//origin
    parser::Vec3f w,v,u;//coordinate system

    //parser::Vec3f ray;
    parser::Vec3f q;// top corner coordinate
    parser::Vec3f s;// pixel coordinate
    parser::Vec3f d;// d in the eq r(t)= e+dt

    unsigned char* image;
    float pixel_width, pixel_height, p_height, p_width;
    int width, height;
    int t_heigth;
    const parser::Scene *scene;
    std::vector<std::vector<parser::Vec3f>> *normals_of_meshes;
} arg;

struct thArg
{
    ARGS arg;
    int start;
    unsigned char* image;

} tharg[Num_Th];

void foo(int arg);
void RayIntersecObj(const parser::Scene &scene,parser::Ray ray, IntersectionData &closest_obj_data);
parser::Intersection intersectRaySphere(const parser::Scene &scene, parser::Ray &eye_ray, int sphere_id);
parser::Vec3f intersectRayFace(const parser::Scene &scene, parser::Ray &eye_ray, const parser::Face &face);
float Determinant();
parser::Vec3f* getRayColor(const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, parser::Ray ray, int recursion_depth, bool isEyeRay);
parser::Vec3f getNormalOfFace(const parser::Scene& scene, const parser::Face& face);

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file

    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    size_t i,j;
    int cam_id; 
    size_t num_of_cameras;
    size_t num_of_triangles;
    size_t num_of_meshes;

    std::vector<std::vector<parser::Vec3f>> normals_of_meshes;
    
    //precompute normals of all triangles
    
    num_of_triangles = scene.triangles.size();
    num_of_meshes = scene.meshes.size();
    num_of_cameras = scene.cameras.size();

    for(i=0; i<num_of_triangles; i++)
        scene.triangles[i].norm = getNormalOfFace(scene, scene.triangles[i].indices); 

    for(i=0; i<num_of_meshes; i++){   
        size_t num_of_faces = scene.meshes[i].faces.size();
        std::vector<parser::Vec3f> face_normals;
        normals_of_meshes.push_back(face_normals);

        for(j=0; j<num_of_faces; j++)
            normals_of_meshes[i].push_back(getNormalOfFace(scene, scene.meshes[i].faces[j]));
    }

    for(cam_id = 0; cam_id < num_of_cameras; cam_id++){
        arg.cam = scene.cameras[cam_id];
        arg.width = arg.cam.image_width;
        arg.height = arg.cam.image_height;
        arg.e = arg.cam.position;
        arg.w = vectorScalerMult(-1.0, arg.cam.gaze);
        arg.v = arg.cam.up;
        arg.u = crossProduct(arg.v, arg.w);
        arg.pixel_width = (arg.cam.near_plane.y - arg.cam.near_plane.x) / arg.cam.image_width;
        arg.pixel_height = (arg.cam.near_plane.w - arg.cam.near_plane.z) / arg.cam.image_height;
        arg.normals_of_meshes = &normals_of_meshes;
        arg.scene = &scene;
        arg.q = arg.cam.position + arg.cam.near_plane.x * arg.u + arg.cam.near_plane.w * arg.v + arg.cam.near_distance * arg.cam.gaze;
        arg.p_height = arg.pixel_height*0.5;
        arg.p_width = arg.pixel_width*0.5;
        arg.image = new unsigned char [arg.width * arg.height * 3];
        arg.t_heigth = arg.height/(Num_Th);
        int height = arg.height/(Num_Th);

        for(i=0; i<Num_Th;i++){
            tharg[i].arg.cam = scene.cameras[cam_id];
            tharg[i].arg.width = arg.width;
            tharg[i].arg.height = arg.height;
            tharg[i].arg.e = arg.cam.position;
            tharg[i].arg.w = vectorScalerMult(-1.0, arg.cam.gaze);
            tharg[i].arg.v = arg.cam.up;
            tharg[i].arg.u = crossProduct(arg.v, arg.w);
            tharg[i].arg.pixel_height = arg.pixel_height;
            tharg[i].arg.pixel_width = arg.pixel_width;
            tharg[i].arg.normals_of_meshes = &normals_of_meshes;
            tharg[i].arg.scene = &scene;
            tharg[i].arg.q = arg.q;
            tharg[i].arg.p_height = arg.p_height;
            tharg[i].arg.p_width = arg.p_width;
            tharg[i].arg.image = new unsigned char [arg.width * arg.height * 3];
            tharg[i].arg.t_heigth = arg.t_heigth;
            tharg[i].start = i*height;
            tharg[i].image = new unsigned char[(height*arg.width*3)];
        }
        
        std::thread th[Num_Th];
        for (size_t i = 0; i < Num_Th; i++)
            th[i] = std::thread(foo,i);
        for (size_t i = 0; i < Num_Th; i++){
            if(th[i].joinable())
                th[i].join();
        }

        //write to file
        int m = 0;
        for(i=0; i<Num_Th;i++){
            int n = 0;
            for(j=0;j<height;j++){
                for(int k=0;k<arg.width;k++){
                    arg.image[m++]= tharg[i].image[n++];
                    arg.image[m++]= tharg[i].image[n++];
                    arg.image[m++]= tharg[i].image[n++];
                }
            }
        }
        write_ppm(scene.cameras[cam_id].image_name.c_str(), arg.image, arg.width, arg.height);
    }
    //write_ppm("test.ppm", image, width, height);

}
// TODO
// add face id in return along as mesh id
void RayIntersecObj(const parser::Scene &scene,parser::Ray ray, IntersectionData &closest_obj_data){
    int numOfSpheres = scene.spheres.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();

    float t_min = __FLT_MAX__;

    closest_obj_data.t = __FLT_MAX__;//closest objects intersection parameter
    closest_obj_data.obj_id = -1;
    closest_obj_data.obj_type = none;


    for(int i = 0; i < numOfSpheres; i++){
        parser::Intersection intersect = intersectRaySphere(scene,ray, i);

        if(intersect.discriminant >= 0 ){
            //meaning they intersect
            if(intersect.t1 < intersect.t2 && intersect.t1 > 0 && intersect.t1 < t_min)
            {
                t_min = intersect.t1;
                closest_obj_data.obj_id = i;
                closest_obj_data.obj_type = sphere;
                    
            }
            else if(intersect.t2 > 0 && intersect.t2 < t_min)
            {
                t_min = intersect.t2;
                closest_obj_data.obj_id = i;
                closest_obj_data.obj_type = sphere;
            }
        }
    }
    for(int i = 0; i<numOfTriangles;i++){
        parser::Vec3f b_g_t = intersectRayFace(scene,ray,scene.triangles[i].indices);
        if(b_g_t.z > 0 && b_g_t.z < t_min && b_g_t.x >= 0 && b_g_t.y >= 0){
            if(b_g_t.y <= 1 && b_g_t.x <= (1 - b_g_t.y)){
                t_min = b_g_t.z;
                closest_obj_data.obj_id = i;
                closest_obj_data.obj_type = triangle;
            }
        }
    }
    for (size_t j = 0; j < numOfMeshes; j++){
        parser::Mesh currmesh = scene.meshes[j];
        int numOfFaces = currmesh.faces.size();
        for(int i = 0; i<numOfFaces;i++){
            parser::Vec3f b_g_t = intersectRayFace(scene,ray,currmesh.faces[i]);
            if(b_g_t.z > 0 && b_g_t.z < t_min && b_g_t.x >= 0 && b_g_t.y >= 0){
                if(b_g_t.y <= 1 && b_g_t.x <= (1 - b_g_t.y)){
                    t_min = b_g_t.z;
                    closest_obj_data.obj_id = j;
                    closest_obj_data.obj_type = mesh;
                    closest_obj_data.face_id = i;

                }
            }
        }
    }
    closest_obj_data.t = t_min;
    //return closest_obj_data.obj_id;
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

bool isShadow(const parser::Scene &scene, parser::Ray ray, float t=1)//t is the parameter for lights position
{
    int numOfSpheres = scene.spheres.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();


    for(int i = 0; i < numOfSpheres; i++){
        parser::Intersection intersect = intersectRaySphere(scene,ray, i);

        if(intersect.discriminant >= 0 ){
            //meaning they intersect
            if(intersect.t1 < intersect.t2 && intersect.t1 > 0 && intersect.t1 < t)
            {
                return true; 
            }
            else if(intersect.t2 > 0 && intersect.t2 < t)
            {
                return true;
            }
        }
    }
    for(int i = 0; i<numOfTriangles;i++){
        parser::Vec3f b_g_t = intersectRayFace(scene,ray,scene.triangles[i].indices);
        if(b_g_t.z > 0 && b_g_t.z < t && b_g_t.x >= 0 && b_g_t.y >= 0){
            if(b_g_t.y <= 1 && b_g_t.x <= (1 - b_g_t.y)){
                return true;
            }
        }
    }
    for (size_t j = 0; j < numOfMeshes; j++){
        parser::Mesh currmesh = scene.meshes[j];
        int numOfFaces = currmesh.faces.size();
        for(int i = 0; i<numOfFaces;i++){
            parser::Vec3f b_g_t = intersectRayFace(scene,ray,currmesh.faces[i]);
            if(b_g_t.z > 0 && b_g_t.z < t && b_g_t.x >= 0 && b_g_t.y >= 0){
                if(b_g_t.y <= 1 && b_g_t.x <= (1 - b_g_t.y)){
                    return true;
                }
            }
        }
    }
    //return closest_obj_data.obj_id;
    return false;


}

parser::Vec3f* getRayColor(const parser::Scene& scene, const std::vector<std::vector<parser::Vec3f>>& mesh_normals, parser::Ray ray, int recursion_depth, bool isEyeRay)
{
    if(recursion_depth < 0)
    {
        return NULL;
    }
    IntersectionData closest_obj_data;
    parser::Vec3f* pixel_color, *reflection_color;
    parser::Vec3f reflection_point, normal_vector;
    parser::Ray reflecting_ray;

    parser::Vec3f center;
    parser::Material material;
    parser::Ray light_ray;
    parser::Vec3f half_vector;

    float light_distance = 0;

    pixel_color = new parser::Vec3f;
    pixel_color->x = 0;     
    pixel_color->y = 0;
    pixel_color->z = 0;

    RayIntersecObj(scene, ray, closest_obj_data);
    int numberOfLightSources = scene.point_lights.size();
    reflection_point = ray.e + closest_obj_data.t * ray.d;

    switch (closest_obj_data.obj_type){
        case none: 
            delete pixel_color;
            return NULL;
            break;

        case sphere:

            material = scene.materials[scene.spheres[closest_obj_data.obj_id].material_id - 1];
            center = scene.vertex_data[scene.spheres[closest_obj_data.obj_id].center_vertex_id -1];

            
            normal_vector = reflection_point - center;
            normal_vector = normalize(normal_vector);
            
            

            //TODO possible error point epsilon


            
            break;

        case triangle:
            material = scene.materials[scene.triangles[closest_obj_data.obj_id].material_id -1];
            normal_vector = scene.triangles[closest_obj_data.obj_id].norm;

            if(dotProduct(normal_vector, ray.d) > 0)//if the normal vector is inverted
            {
                normal_vector = -normal_vector;//invert it back
            }


            break;
        case mesh:
            material = scene.materials[scene.meshes[closest_obj_data.obj_id].material_id -1];
            normal_vector = mesh_normals[closest_obj_data.obj_id][closest_obj_data.face_id];

            if(dotProduct(normal_vector, ray.d) > 0)//if the normal vector is inverted
            {
                normal_vector = -normal_vector;//invert it back
            }

            break;
    }

    reflection_point +=  scene.shadow_ray_epsilon * normal_vector; //add offset to reflection point by epsilon
    reflecting_ray.e = reflection_point;
    reflecting_ray.d = ray.d - 2* dotProduct(ray.d, normal_vector) *normal_vector;

    for(int i=0; i<numberOfLightSources; i++)
    {
        parser::PointLight point_light = scene.point_lights[i];
        float cosine_theta = 0;

        light_ray.e = reflection_point;
        light_ray.d = point_light.position - reflection_point;

        //check if the point is in shadow
        //if(false)
        if(isShadow(scene, light_ray))
        {
            continue;//if in shadow diffuse and specular components are not calculated
        }
        light_distance = VECTOR_LENGTH(light_ray.d);
        light_ray.d = normalize(light_ray.d);//this must be called AFTER isShaow is called

        cosine_theta = dotProduct(light_ray.d, normal_vector);
        cosine_theta = MAX(0, cosine_theta);//TODO be carefull if the normal is negative

        point_light.intensity = point_light.intensity /(light_distance*light_distance);// I/r^2

        //*pixel_color += clampColor(cosine_theta * elementViseMultiply(material.diffuse, point_light.intensity));
        *pixel_color += cosine_theta * elementViseMultiply(material.diffuse, point_light.intensity);
        //TODO add specular component

        half_vector = normalize(-ray.d) + light_ray.d;
        half_vector = normalize(half_vector);

        cosine_theta = dotProduct(half_vector, normal_vector);
        cosine_theta = MAX(0, cosine_theta);
        cosine_theta = pow(cosine_theta, material.phong_exponent);

        //*pixel_color += clampColor(cosine_theta * elementViseMultiply(point_light.intensity, material.specular));
        *pixel_color += cosine_theta * elementViseMultiply(point_light.intensity, material.specular);
    }

    if(material.is_mirror){
        reflection_color = getRayColor(scene, mesh_normals, reflecting_ray, recursion_depth-1, false);
        if(reflection_color != NULL)
        {
            //*pixel_color += clampColor(elementViseMultiply(material.mirror, *reflection_color));
            *pixel_color += elementViseMultiply(material.mirror, *reflection_color);
        }
    }

    *pixel_color += elementViseMultiply(scene.ambient_light, material.ambient);

    return pixel_color;
}
parser::Vec3f getNormalOfFace(const parser::Scene& scene, const parser::Face& face)
{
    parser::Vec3f A,B,C;//corners of the face
    parser::Vec3f a,b;  //edges of the triangle
    parser::Vec3f normal_vector;

    A = scene.vertex_data[face.v0_id-1];
    B = scene.vertex_data[face.v1_id-1];
    C = scene.vertex_data[face.v2_id-1];

    a = B - C;
    b = A - C;

    normal_vector = crossProduct(a, b);

    normal_vector = normalize(normal_vector);
    return normal_vector;
}

void foo(int thno){
    int i = 0;
    for(int y=0; y < tharg[thno].arg.t_heigth; y++){
        int tmp = tharg[thno].start + y;
        for(int x=0; x < tharg[thno].arg.width; x++){
            //first compute pixel coordinates;
            parser::Ray r;
            //parser::Material material;
            parser::Vec3f* pixel_color;

            parser::Vec3f temp_vec = (x*tharg[thno].arg.pixel_width+tharg[thno].arg.p_width)* tharg[thno].arg.u;

            tharg[thno].arg.s = tharg[thno].arg.q + temp_vec;//s = q + s_u.u

            temp_vec =  (tmp*tharg[thno].arg.pixel_height + tharg[thno].arg.p_height) * tharg[thno].arg.v;// temp_vec =  s_v . v

            tharg[thno].arg.s -= temp_vec;
            //s = q + s_u.u - s_v.v

            tharg[thno].arg.d = tharg[thno].arg.s - tharg[thno].arg.e;//s-e

            r.e = tharg[thno].arg.e;
            r.d = tharg[thno].arg.d;

            //float intersect.t1, intersect.t2;//different solutions of the equation

            pixel_color = getRayColor(*(tharg[thno].arg.scene), *(tharg[thno].arg.normals_of_meshes), r, tharg[thno].arg.scene->max_recursion_depth, true);

            if(pixel_color == NULL)//meaning ray didn't hit any objects
            {
                pixel_color = new parser::Vec3f;
                pixel_color->x = tharg[thno].arg.scene->background_color.x;
                pixel_color->y = tharg[thno].arg.scene->background_color.y;
                pixel_color->z = tharg[thno].arg.scene->background_color.z;
            }
            *pixel_color = clampColor(*pixel_color);
                
            tharg[thno].image[i++] = pixel_color->x;
            tharg[thno].image[i++] = pixel_color->y;
            tharg[thno].image[i++] = pixel_color->z;
        }
    }
}
