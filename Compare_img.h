//
//  Compare_img.h
//  Project
//
//  Created by Lu Ang on 30/6/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//

#ifndef __Project__Compare_img__
#define __Project__Compare_img__

#include <opencv2/imgproc/imgproc.hpp>
#include <stdio.h>
using namespace cv;
double SSD(Mat image1, Mat image2);
double SSD(float * img1, float * img2, int w, int h);
double ImageMin(Mat img1, Mat img2);
void diff_img(float * img1, float * img2, int w, int h, float *img3);

#endif /* defined(__Project__Compare_img__) */
