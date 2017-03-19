//
//  Model_OBJ.h
//  Project
//
//  Created by Lu Ang on 30/6/15.
//  Copyright (c) 2015 Lu Ang. All rights reserved.
//

#ifndef __Project__Model_OBJ__
#define __Project__Model_OBJ__

#include <stdio.h>
/***************************************************************************
 OBJ Loading
 ***************************************************************************/

class Model_OBJ
{
public:
    Model_OBJ();
    float* calculateNormal(float* coord1,float* coord2,float* coord3 );
    int Load(char *filename);       // Loads the model
    void Draw();					// Draws the model on the screen
    void Release();                 // Release the model
    
    float* normals;							// Stores the normals
    float* Faces_Triangles;					// Stores the triangles
    float* vertexBuffer;					// Stores the points which make the object
    long TotalConnectedPoints;				// Stores the total number of connected verteces
    long TotalConnectedTriangles;			// Stores the total number of connected triangles
    
};

#endif /* defined(__Project__Model_OBJ__) */
