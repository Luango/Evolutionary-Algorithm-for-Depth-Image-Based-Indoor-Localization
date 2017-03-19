//
//  Camera.h
//  Project
//
//  Created by Lu Ang on 2/7/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//

#ifndef __Project__Camera__
#define __Project__Camera__
#include <iostream>
#include <stdio.h>
#include "Eigen/Geometry"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

///// WITHOUT QUATERNION
/*
class Camera
{
public:
    double eye_x;
    double eye_y;
    double eye_z;
    double center_x;
    double center_y;
    double center_z;
    double ssd_score; // Each camera has a SSD score compared with input image.
    
    Camera(double x1, double y1, double z1, double x2, double y2, double z2, double score);
    Camera(Camera const &a);
    Camera();
    void operator=(Camera a);
    void copy(Camera a);
    
    
    void show(){
        cout<<"eye.x = "<<eye_x<<endl;
        cout<<"eye.y = "<<eye_y<<endl;
        cout<<"eye.z = "<<eye_z<<endl;
        
        cout<<"center.x = "<<center_x<<endl;
        cout<<"center.y = "<<center_y<<endl;
        cout<<"center.z = "<<center_z<<endl;
    }
};
*/

/////// WITH QUATERNION
class Camera
{
public:
    Quaternion<double> quaternionCamera;
    double eye_x;
    double eye_y;
    double eye_z;
    double center_x;
    double center_y;
    double center_z;
    double up_x;
    double up_y;
    double up_z;
    double ssd_score; // Each camera has a SSD score compared with input image.
    double mutationVec[7];
    
    //Camera(double x1, double y1, double z1, double x2, double y2, double z2,double w2, double score);
    Camera(Camera inputCamera, double posRange, double rotRange);
    Camera(Camera const &a);
    Camera();
    Camera(double x, double y, double z, double qx, double qy, double qz, double qw);
    void operator=(Camera a);
    void copy(Camera a);
    
    
    void show() const
    {
        cout<<"eye.x = "<<eye_x<<endl;
        cout<<"eye.y = "<<eye_y<<endl;
        cout<<"eye.z = "<<eye_z<<endl;
        cout<<"quaternion_x = "<<quaternionCamera.x()<<endl;
        cout<<"quaternion_y = "<<quaternionCamera.y()<<endl;
        cout<<"quaternion_z = "<<quaternionCamera.z()<<endl;
        cout<<"quaternion_w = "<<quaternionCamera.w()<<endl;
        cout<<"up_x = "<<up_x<<endl;
        cout<<"up_y = "<<up_y<<endl;
        cout<<"up_z = "<<up_z<<endl<<endl;
        /*
        cout<<"mutation Vec 1 = "<<mutationVec[0]<<endl;
        cout<<"mutation Vec 2 = "<<mutationVec[1]<<endl;
        cout<<"mutation Vec 3 = "<<mutationVec[2]<<endl;
        cout<<"mutation Vec 4 = "<<mutationVec[3]<<endl;
        cout<<"mutation Vec 5 = "<<mutationVec[4]<<endl;
        cout<<"mutation Vec 6 = "<<mutationVec[5]<<endl;
        cout<<"mutation Vec 7 = "<<mutationVec[6]<<endl;
         */
    }
};
#endif /* defined(__Project__Camera__) */
