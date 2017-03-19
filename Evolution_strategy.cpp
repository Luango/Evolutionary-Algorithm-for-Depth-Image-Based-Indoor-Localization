//
//  Evolution_strategy.cpp
//  Project
//
//  Created by Lu Ang on 2/7/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//
#include "Header.h"

using namespace cv;
using namespace std;

// Updating camera pose 6-DOF quaternion (ES)
void evolution_strategy(Camera *camera, RNG *rng, double posStepFactor, double rotStepFactor)
{
    // Random update variable
    double update_eye_x, update_eye_y, update_eye_z;
    double update_center_x, update_center_y, update_center_z;
    double update_up_x, update_up_y, update_up_z;
    
    double position_dotProduct;
    // Update the eye position
    update_eye_x = rng->gaussian(1.0);
    update_eye_y = rng->gaussian(1.0);
    update_eye_z = rng->gaussian(1.0);
    camera->mutationVec[0] = update_eye_x;
    camera->mutationVec[1] = update_eye_y;
    camera->mutationVec[2] = update_eye_z;
    
    // Update camera position
    camera->eye_x += update_eye_x*posStepFactor;
    camera->eye_y += update_eye_y*posStepFactor;
    camera->eye_z += update_eye_z*posStepFactor;
    
    // Update camera center
    camera->center_x += update_eye_x*posStepFactor;
    camera->center_y += update_eye_y*posStepFactor;
    camera->center_z += update_eye_z*posStepFactor;
    
    // Update quaternion rotation
    // Using PPSN 2014
    Quaternion<double> update_quaternion;
    Quaternion<double> w;
    w.x() = rng->gaussian(1.0);
    w.y() = rng->gaussian(1.0);
    w.z() = rng->gaussian(1.0);
    w.w() = rng->gaussian(1.0);
    // Dot product
    double dot_product;
    Quaternion<double> z;
    
    dot_product = w.dot(camera->quaternionCamera);
    
    z.x() = w.x() - dot_product * camera->quaternionCamera.x();
    z.y() = w.y() - dot_product * camera->quaternionCamera.y();
    z.z() = w.z() - dot_product * camera->quaternionCamera.z();
    z.w() = w.w() - dot_product * camera->quaternionCamera.w();
    
    camera->mutationVec[3] = z.x();
    camera->mutationVec[4] = z.y();
    camera->mutationVec[5] = z.z();
    camera->mutationVec[6] = z.w();
    
    double cosFactor;
    double sinFactor;
    cosFactor = cos(sqrt(z.dot(z))*rotStepFactor);
    sinFactor = (sin(sqrt(z.dot(z))*rotStepFactor))/(sqrt(z.dot(z)));
    
    camera->quaternionCamera.x() = cosFactor * camera->quaternionCamera.x() + sinFactor * z.x();
    camera->quaternionCamera.y() = cosFactor * camera->quaternionCamera.y() + sinFactor * z.y();
    camera->quaternionCamera.z() = cosFactor * camera->quaternionCamera.z() + sinFactor * z.z();
    camera->quaternionCamera.w() = cosFactor * camera->quaternionCamera.w() + sinFactor * z.w();
    camera->quaternionCamera.normalize();
    
    // Rotation matrix
    MatrixXd m(3,3);
    m = camera->quaternionCamera.toRotationMatrix();
    
    // Update center point
    Vector3d directionVector;
    directionVector<<0,0,-1;
    Vector3d b;
    b = m * directionVector;
    update_center_x = camera->eye_x + b(0);
    update_center_y = camera->eye_y + b(1);
    update_center_z = camera->eye_z + b(2);
    
    // Update up vector
    directionVector<<0,1,0;
    Vector3d c;
    c = m * directionVector;
    update_up_x = c(0);
    update_up_y = c(1);
    update_up_z = c(2);
    
    // Update camera rotation
    camera->center_x = update_center_x;
    camera->center_y = update_center_y;
    camera->center_z = update_center_z;
    camera->up_x = update_up_x;
    camera->up_y = update_up_y;
    camera->up_z = update_up_z;
}