//
//  Evolution_strategy.h
//  Project
//
//  Created by Lu Ang on 2/7/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//

#ifndef __Project__Evolution_strategy__
#define __Project__Evolution_strategy__

#include <stdio.h>
#include "Camera.h"
using namespace cv;

void evolution_strategy(Camera *camera, RNG *rng, double posStepFactor, double rotStepFactor);
void evolution_strategy(Camera *camera, RNG *rng, double stepFactor, Vector3d &sx, Quaternion<double> &sq, double &alpha);
void evolution_strategy(Camera *camera, RNG *rng, double posStepFactor, double rotStepFactor, Quaternion<double> &w);
void evolution_strategy(Camera *camera, RNG *rng, double stepFactor, Vector3d &sx, Quaternion<double> &sq, double &alpha);
#endif /* defined(__Project__Evolution_strategy__) */
