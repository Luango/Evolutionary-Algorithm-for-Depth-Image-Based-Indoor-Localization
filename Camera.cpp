//
//  Camera.cpp
//  Project
//
//  Created by Lu Ang on 2/7/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//
#include "Header.h"

using namespace cv;
using namespace std;


// With quaternion
// Generate random camera around input camera
Camera::Camera(Camera inputCamera, double posRange, double rotRange)
{
    RNG rng(time(NULL));
    // Initialize camera in position close to goal position
    
    eye_x = inputCamera.eye_x + rng.gaussian(posRange);
    eye_y = inputCamera.eye_y + rng.gaussian(posRange);
    eye_z = inputCamera.eye_z + rng.gaussian(posRange);
    
    
    // Initialize camera rotation near the right angle (shot (1: 0, -1, -, 24))
    quaternionCamera.x() = inputCamera.quaternionCamera.x() + rng.gaussian(rotRange);
    quaternionCamera.y() = inputCamera.quaternionCamera.y() + rng.gaussian(rotRange);
    quaternionCamera.z() = inputCamera.quaternionCamera.z() + rng.gaussian(rotRange);
    quaternionCamera.w() = inputCamera.quaternionCamera.w() + rng.gaussian(rotRange);
    quaternionCamera.normalize(); // Normalize
    
    // To test different pose how it works
    
    // Rotation Matrix
    MatrixXd m(3,3);
    m = quaternionCamera.toRotationMatrix();
    
    // Center point
    Vector3d a;
    a<<0,0,-1;
    Vector3d b;
    b = m * a;
    
    center_x = eye_x + b(0);
    center_y = eye_y + b(1);
    center_z = eye_z + b(2);
    
    // Up vector
    a<<0,1,0;
    Vector3d c;
    c = m * a;
    
    up_x = c(0);
    up_y = c(1);
    up_z = c(2);
    
    ssd_score = 0;
    
    for (int i=0; i<7; i++) {
        mutationVec[i]=0;
    }
}
Camera::Camera()
{
    RNG rng(time(NULL));
    
    // Initialize camera in position
    eye_x = 0;
    eye_y = 0;
    eye_z = 0;
    
    // Initialize camera rotation
    quaternionCamera.x() = 0;
    quaternionCamera.y() = 1;
    quaternionCamera.z() = 0;
    quaternionCamera.w() = 0;
    quaternionCamera.normalize(); // Normalize
    
    // Rotation Matrix
    MatrixXd m(3,3);
    m = quaternionCamera.toRotationMatrix();
    
    // Center point
    Vector3d a;
    a<<0,0,-1;
    Vector3d b;
    b = m * a;
    
    center_x = eye_x + b(0);
    center_y = eye_y + b(1);
    center_z = eye_z + b(2);
    
    // Up vector
    a<<0,1,0;
    Vector3d c;
    c = m * a;
    
    up_x = c(0);
    up_y = c(1);
    up_z = c(2);
    
    ssd_score = 0;
    
    for (int i=0; i<7; i++) {
        mutationVec[i] = 0;
    }
}

Camera::Camera(double x, double y, double z, double qx, double qy, double qz, double qw){
    eye_x = x;
    eye_y = y;
    eye_z = z;

    // Initialize camera rotation near the right angle (shot (1: 0, -1, -, 24))
    quaternionCamera.x() = qx;
    quaternionCamera.y() = qy;
    quaternionCamera.z() = qz;
    quaternionCamera.w() = qw;
    quaternionCamera.normalize(); // Normalize
    
    // Rotation Matrix
    MatrixXd m(3,3);
    m = quaternionCamera.toRotationMatrix();
    
    // Center point
    Vector3d a;
    a<<0,0,-1;
    Vector3d b;
    b = m * a;
    
    center_x = eye_x + b(0);
    center_y = eye_y + b(1);
    center_z = eye_z + b(2);
    
    // Up vector
    a<<0,1,0;
    Vector3d c;
    c = m * a;
    
    up_x = c(0);
    up_y = c(1);
    up_z = c(2);
    
    ssd_score = 0;
    
    for (int i=0; i<7; i++) {
        mutationVec[i] = 0;
    }
}

Camera::Camera(Camera const &a)
{
    eye_x = a.eye_x;
    eye_y = a.eye_y;
    eye_z = a.eye_z;
    center_x = a.center_x;
    center_y = a.center_y;
    center_z = a.center_z;
    quaternionCamera = a.quaternionCamera;
    up_x = a.up_x;
    up_y = a.up_y;
    up_z = a.up_z;
    
    for (int i=0; i<7; i++) {
        mutationVec[i] = a.mutationVec[i];
    }
    
    ssd_score = a.ssd_score;
}
void Camera::operator=(Camera a)
{
    eye_x = a.eye_x;
    eye_y = a.eye_y;
    eye_z = a.eye_z;
    center_x = a.center_x;
    center_y = a.center_y;
    center_z = a.center_z;
    quaternionCamera = a.quaternionCamera;
    up_x = a.up_x;
    up_y = a.up_y;
    up_z = a.up_z;
    for (int i=0; i<7; i++) {
        mutationVec[i] = a.mutationVec[i];
    }
    ssd_score = a.ssd_score;
}