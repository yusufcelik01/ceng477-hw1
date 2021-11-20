#include <iostream>
#include "parser.h"
#include "ppm.h"

#include <cmath>
#include "vector_op.h"

typedef unsigned char RGB[3];

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);


    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int i,j,foo,bar;
    int width, height;

    /*
    width = 640;
    height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    i = 0;
        float w, v, u;
    
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }

    write_ppm("test.ppm", image, width, height);
    */

    
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
        for(int y=0; y < cam.image_height; y++)
        {
            for(int x=0; x < cam.image_width; x++)
            {
                //first compute pixel coordinates;
                temp_vec = vectorScalerMult((x+0.5)*pixel_width, u);

                s = vectorSum(q, temp_vec);//s = q + s_u.u

                temp_vec = vectorScalerMult( -1.0 * (y+0.5)*pixel_height, v);// - s_v . v

                s = vectorSum(s, temp_vec);
                //s = q + s_u.u - s_v.v


                d = vectorSum(s, vectorScalerMult(-1.0, e));//s-e

                float t_min;//closest objects parameter
                float t_1, t_2;//different solutions of the equation

                //calculate spheres' closest

                int numOfSpheres = scene.spheres.size();
                int intersects = 0;//whether it intersects or not
                int closest_sphere;//id

                for(foo = 0; foo < numOfSpheres; foo++)
                {
                    parser::Vec3f c;
                    float radius;

                    c = scene.vertex_data[scene.spheres[foo].center_vertex_id -1 ];
                    radius = scene.spheres[foo].radius;

                    //std::cout << "center " << c.x << "," << c.y << "," << c.z  << " radius: " << radius << std::endl;

                    float discriminant;
                    parser::Vec3f e_c; //e-c and d^2 is freq used so assign it to a variable
                    float d_sqr;

                    e_c = vectorSum(e, vectorScalerMult(-1.0, c));
                    d_sqr = dotProduct(d, d);


                    temp = dotProduct(d, e_c);
                    temp *= temp;//(d.(e-c))^2

                    discriminant = temp;
                    discriminant -= d_sqr * (dotProduct(e_c, e_c) - radius*radius);
                    //d^2(e-c)^2-r^2

                    if(discriminant >= 0)//meaning they intersect
                    {
                        t_1 = -dotProduct(d, e_c);
                        t_2 = t_1;

                        t_1 += sqrt(discriminant);//-b + sqrt delta
                        t_1 /= d_sqr;
                        t_2 -= sqrt(discriminant);//-b - sqrt delta
                        t_2 /= d_sqr;

                        if(!intersects)
                        {
                            intersects = 1;
                            //assign the smallest
                            if(t_1 < t_2)
                            {
                                t_min = t_1;
                            }
                            else
                            {
                                t_min = t_2;
                            }
                            closest_sphere = foo;
                        }
                        else
                        {
                            if(t_1 < t_min)
                            {
                                t_min = t_1;
                                closest_sphere = foo;
                            }
                            if(t_2 < t_min);
                            {
                                t_min = t_2;
                                closest_sphere = foo;
                            }

                        }
                        
                    }



                }

                //find the colour of that material

                if(!intersects)// does not intersect get backround colour
                {
                    parser::Vec3i bg;//backgroung
                    bg = scene.background_color;
                    image[i++] = bg.x;
                    image[i++] = bg.y;
                    image[i++] = bg.z;
                }
                else// if there is a sphere just write green
                {

                    parser::Material material = scene.materials[scene.spheres[closest_sphere].material_id - 1];


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
