#include <iostream>
#include "parser.h"
#include "ppm.h"

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

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
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

//checking given structures
//*************************

    //std::cout<< "cameras" ;


//*************************
//checking given structures
    
    int cam_id; 
    int num_of_cameras;

    

    num_of_cameras = scene.cameras.size();

    for(cam_id = 0; cam_id < num_of_cameras; cam_id++)
    {
        parser::Camera cam;
        parser::Vec3f ray;
        parser::Vec3f q;// top corner coordinate
        parser::Vec3f s;// pixel coordinate
        parser::Vec3f w,v,u;//coordinate system


        float pixel_width, pixel_height;

        cam = scene.cameras[cam_id];
        w = vectorScalerMult(-1.0, cam.gaze);
        v = cam.up;
        u = crossProduct(v, w);



        pixel_width = (cam.near_plane.y - cam.near_plane.x) / cam.image_width;
        pixel_height = (cam.near_plane.w - cam.near_plane.z) / cam.image_height;

        
        parser::Vec3f temp_vec;
        temp_vec = vectorScalerMult(cam.near_plane.x, u);//l.u

        q = vectorSum(cam.position, temp_vec);//e + l.u
        
        temp_vec = vectorScalerMult(cam.near_plane.w, v);//t.v
        q = vectorSum(q, temp_vec); // += t.v

        temp_vec = vectorScalerMult(cam.near_distance, cam.gaze);
        q = vectorSum(q, temp_vec);
        //q has the position of top left corner



        for(int y=0; y < cam.image_height; y++)
        {
            for(int x=0; x < cam.image_width; x++)
            {
                //first compute pixel coordinates;
                temp_vec = vectorScalerMult((x+0.5)*pixel_width, u);

                s = vectorSum(q, temp_vec);//s = q + s_u.u

                temp_vec = vectorScalerMult( (y+0.5)*pixel_height, v);

                s = vectorSum(s, temp_vec);
                //s = q + s_u.u + s_v.v



                


                //find ray equation
                //ray = cam.position + ()



            }
        }

        parser::Vec3f e;//origin of ray
        

        //compute near plane's pixels
        //send ray

            //for all objects check if they intersect
            //if so then find the colour


        //compute image

         
        //write to file
        write_ppm(scene.cameras[cam_id].image_name.c_str(), image, width, height);

    }


    

    write_ppm("test.ppm", image, width, height);

}
