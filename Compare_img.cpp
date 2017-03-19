//
//  Compare_img.cpp
//  Project
//
//  Created by Lu Ang on 30/6/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//

#include "Compare_img.h"
#include "Header.h"
/*
// SSD is the scoring function
double SSD(Mat img1, Mat img2)
{
    double SSD_score=0;
    int row;
    int col;
    row = img1.rows;
    col = img1.cols;
    
    for (int i=0; i<row; i++) {
        for (int j=0; j<col; j++) {
            Scalar intensity1 = img1.at<uchar>(i, j);
            Scalar intensity2 = img2.at<uchar>(i, j);
            // Sum of squared difference
            SSD_score +=  ((intensity1[0]-intensity2[0])*(intensity1[0]-intensity2[0]));
        }
    }
    return SSD_score;
}
 */
double SSD(float * img1, float * img2, int w, int h)
{
    double SSD_score=0.0;
    
    for (int i=0; i<w*h; i++) {
        //img1[i] = 2 * 60 * 0.5 / (60 + 0.5 - (60 - 0.5) * ( 2 * img1[i] - 1));
        //img2[i] = 2 * 60 * 0.5 / (60 + 0.5 - (60 - 0.5) * ( 2 * img2[i] - 1));
        SSD_score += (img1[i]-img2[i])*(img1[i]-img2[i]);
    }
    return SSD_score;
}

void diff_img(float * img1, float * img2, int w, int h, float *img3)
{
    for (int i=0; i<w*h; i++) {
        //float im1 = ((20/img1[i])-20-0.5)/(-2*(20-0.5))+1;
        //float im2 = ((20/img2[i])-20-0.5)/(-2*(20-0.5))+1;
        //img3[i] = sqrt(((im1-im2)*(im1-im2)));
        img3[i] = sqrt(((img1[i]-img2[i])*(img1[i]-img2[i])))/10;
    }
}