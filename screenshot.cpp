//
//  screenshot.cpp
//  Project
//
// Part from http://www.david-amador.com/2012/09/how-to-take-screenshot-in-opengl/ , include glPixelStorei and glReadPixels
//
#include "Header.h"
#include "screenshot.h"
bool save_screenshot(string filename, int w, int h)
{
    // Using OpenCV to store image
    Mat img(h, w, CV_8UC3);
    Mat gray_image;
    glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
    glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
    cvtColor(img, gray_image, CV_BGR2GRAY);
    flip(gray_image, gray_image, 0);
    imwrite(filename, gray_image);
    
    namedWindow( "Display window", WINDOW_AUTOSIZE ); // Create a window for display.
    imshow( "Display window", gray_image );                // Show our image inside it.
    return true;
}