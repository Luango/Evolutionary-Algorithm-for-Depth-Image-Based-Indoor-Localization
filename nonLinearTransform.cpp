//
//  Logarithm.cpp
//  Project2
//
//  Created by Lu Ang on 22/8/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//
#include "Header.h"
#include "nonLinearTransform.h"

void nonLinearTransform(float *img[], int w, int h)
{
    
    for (int i=0; i<w*h; i++) {
        //Concatenation
        //(*img)[i] = ( 2 * 0.1 ) / ( 10.0 + 0.1 - (*img)[i] * ( 10.0 - 0.1 ) );
        //(*img)[i] = 2 * 2.0 * 0.5 / (2.0 + 0.5 - (2.0 - 0.5) * ( 2 * (*img)[i] - 1));
        
        //World coordinate transformation
        (*img)[i] = 2 * 20 * 0.1 / (20 + 0.1 + (20 - 0.1) * ( 1 - 2*(*img)[i]));
    }
}
