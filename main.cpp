#include "Header.h"
#define KEY_ESCAPE 27
#define KEY_GOOD 103
#define KEY_BAD 98
#define addRotStep 97
#define addPosStep 99

const double PI = 3.14159265395;
int start_time;
RNG rng(time(NULL));

// Process num
int process_num=0;
// Stopping criteria count
int stopCounter = 0;
// Count number
int counter = 0;

// Output file
ofstream ssdValueFile[20];
ofstream sigmaLocFile[20],sigmaRotFile[20],locDifFile[20], rotDifFile[20];
ofstream iterationCountFile;
ofstream ratioFile;
ofstream rotSSDFile, posSSDFile;
ofstream successRateFile;
ofstream rotDisFile,posDisFile;
ofstream cvFile, meanFile, sdFile;
ofstream subProcessNumFile;

// window
typedef struct {
    int width;
    int height;
    char* title;
    
    float field_of_view_angle;
    float z_near;
    float z_far;
} glutWindow;

// Calculate Path
Vector3d sx;
Quaternion<double> sq;
VectorXd stepSizePath(7);
double mutationVector[7]={0};

Model_OBJ obj; // The Model
glutWindow win; // Window

Mat show_input_image; // Using openCV to show image
float * input_image;  // This array save the input image
float * current_image; // This array save the current image that rendered from camera
float * squar_diff_image; // This is squared of difference image

// Cameras
Camera one_camera; // one_camera is the current camera that used for mutate
Camera save_camera; // save_camera is used to save the previous camera, if current camera is not suppose to exist, use save_camera overwrite it.
//Camera inputCamera(0.3,1.5,-5,0,1,0,0); //This is where input camera is initialised

//Classroom 4 initialization
//Camera inputCamera(1.5,2.3,-7,0.3,1,0.2,0.1);
//Camera inputCamera(-0.5,2.8,-3,0.1,-0.9,-0.3,0.0);
//Camera inputCamera(2.0,2.6,-4,0.1,0.9,0.1,0.9);
Camera inputCamera(1,1.9,-5.5,-0.1,-0.7,-0.1,0.2);
// Familia
//Camera inputCamera(0.0,0.3,-1.5,0.1,1,0.1,0.1); //This is where input camera is initialised
Quaternion<double> w; // This quaternion w is used for PPSN 2014 parallel transport purpose
Camera iniCamera[10000];

// Previous SSD
double best_ssd = 20000000000000;
double previous_ssd; // used to save previous ssd value which is the best ssd
double current_ssd;  // current ssd value

// Image sequence
int sequence_number = 0;
string location = "/Users/luang/Desktop/";

// Step size
double positionStepFactor; // position step size
double rotationStepFactor; // rotation step size
double stepFactor;  // one single step size
double stepRatio;
const double tau = 0.3;
double sigma;
const int populationSize = 10;
const int selectSize = 3;

// Better SSD
bool betterSSD = false;

// Check if first time. if it is their first time, initialise the input camera.
bool isFirstTime = true;

// Output file
ofstream ssd_output;
ofstream iteration_num;

// Counting
int counting = 0;
int subCounting = 0;
int subProcessNum = 0;

const double eyeDistanceWrite()
{
    // Calculate and write eye distance value
    double delta_x, delta_y, delta_z;
    delta_x = abs(one_camera.eye_x - inputCamera.eye_x);
    delta_y = abs(one_camera.eye_y - inputCamera.eye_y);
    delta_z = abs(one_camera.eye_z - inputCamera.eye_z);
    double eyeDistanceValue = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
    //cout<<"pos : "<<eyeDistanceValue<<endl;
    
    return eyeDistanceValue;
}
const double rotationDistanceWrite()
{
    // Calculate and write rotation value
    double x,y,z,w;
    double difference2q;
    double angle2q;
    x = inputCamera.quaternionCamera.x();
    y = inputCamera.quaternionCamera.y();
    z = inputCamera.quaternionCamera.z();
    w = inputCamera.quaternionCamera.w();
    difference2q = (one_camera.quaternionCamera.x()*x+one_camera.quaternionCamera.y()*y+one_camera.quaternionCamera.z()*z+one_camera.quaternionCamera.w()*w);
    
    angle2q = 2*acos(fabs(difference2q));
    //cout<<"rot : "<<angle2q<<endl;
    
    return angle2q;
}

// Display 3D model and set depth map.
const void display(float *depthImg[])
{
    int w = win.width;
    int h = win.height;
    
    // Display 3D model
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, win.width, win.height);
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    double ar = w / static_cast< double >( h );
    const float zNear = 0.1;
    const float zFar = 10.0;
    gluPerspective( 60, ar, zNear, zFar );
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // This camera has all 6 DOF (Degree of freedom) 9 parameters
    gluLookAt(one_camera.eye_x, one_camera.eye_y, one_camera.eye_z, one_camera.center_x, one_camera.center_y, one_camera.center_z, one_camera.up_x, one_camera.up_y, one_camera.up_z);

    glPushMatrix();
    obj.Draw();
    glPopMatrix();
    
    //vector< GLfloat > depth( w * h, 0 );
    float *depth = new float [w*h];
    glReadPixels( 0, 0, w, h, GL_DEPTH_COMPONENT, GL_FLOAT, depth );
    *depthImg = depth;
    
    /*
    static GLuint tex = 0;
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_LUMINANCE, w, h, 0, GL_LUMINANCE, GL_FLOAT, &depth[0] );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    
    glOrtho( 0, w, 0, h, -1, 1 );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glEnable( GL_TEXTURE_2D );
    glColor3ub( 255, 255, 255 );
    glScalef( 1, 1, 1 );
    glBegin( GL_QUADS );
    glTexCoord2i( 0, 0 );
    glVertex2i( 0, 0 );
    glTexCoord2i( 1, 0 );
    glVertex2i( w, 0 );
    glTexCoord2i( 1, 1 );
    glVertex2i( w, h);
    glTexCoord2i( 0, 1 );
    glVertex2i( 0, h );
    glEnd();
    glutSwapBuffers();
     */
}

// Render depth image and using OpenCV to store present image
const void render(Mat *gray_image)
{
    int w = win.width;
    int h = win.height;
    // Render depth image and using OpenCV to store present image
    Mat img(h, w, CV_8UC3);
    
    location = "/Users/luang/Desktop/sequence/" + to_string(sequence_number)+".jpg";
    glPixelStoref(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
    glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
    cvtColor(img, *gray_image, CV_BGR2GRAY);
    flip(*gray_image, *gray_image, 0);
    // Save image
    imwrite(location, *gray_image);
    sequence_number++;
}
void ini(Camera *ini){
    //CAMERA classroom
    ini->eye_x = rng.uniform(-1.0, 3.0);
    ini->eye_y = rng.uniform(1.0, 2.8);
    ini->eye_z = rng.uniform(-9.0, -2.0);
    
    ini->quaternionCamera.x() = rng.gaussian(1.0);
    ini->quaternionCamera.y() = rng.gaussian(1.0);
    ini->quaternionCamera.z() = rng.gaussian(1.0);
    ini->quaternionCamera.w() = rng.gaussian(1.0);
    /*
    ini->quaternionCamera.x() = rng.uniform(0.0,1.0);
    ini->quaternionCamera.y() = rng.uniform(0.0,1.0);
    ini->quaternionCamera.z() = rng.uniform(0.0,1.0);
    ini->quaternionCamera.w() = rng.uniform(0.0,1.0);
    ini->quaternionCamera.normalize();
    */
    /*
    //CAMERA church
    //inputCamera(0.0,2.6,0.1,0.5,1,0.2,0.2)
    ini->eye_x = rng.uniform(-0.5, 0.5);
    ini->eye_y = rng.uniform(-0.2, 0.8);
    ini->eye_z = rng.uniform(-1.0, -0.5);
    
    ini->quaternionCamera.x() = rng.gaussian(1.0);
    ini->quaternionCamera.y() = rng.gaussian(1.0);
    ini->quaternionCamera.z() = rng.gaussian(1.0);
    ini->quaternionCamera.w() = rng.gaussian(1.0);
    ini->quaternionCamera.normalize();
    */
    
    // Rotation Matrix
    MatrixXd m(3,3);
    m = ini->quaternionCamera.toRotationMatrix();
    
    // Center point
    Vector3d a;
    a<<0,0,-1;
    Vector3d b;
    b = m * a;
    
    ini->center_x = ini->eye_x + b(0);
    ini->center_y = ini->eye_y + b(1);
    ini->center_z = ini->eye_z + b(2);
    
    // Up vector
    a<<0,1,0;
    Vector3d c;
    c = m * a;
    
    ini->up_x = c(0);
    ini->up_y = c(1);
    ini->up_z = c(2);
    
    ini->ssd_score = 0;
    
    for (int i=0; i<7; i++) {
        ini->mutationVec[i] = 0;
    }
}
const void parameterInitialize()
{
    sigma = 0.2;
    previous_ssd = 2000000000000000;
    stepRatio = 1;
    
    // To initialize camera that used as input camera.
    //one_camera = *new Camera(inputCamera, 0.2, 0.1);
    //one_camera = *new Camera(1.26242,1.42547,-5.05309,0.312132,0.823215,0.219588,0.320322);
    //one_camera = *new Camera(1.3,1.8,-0.5,0,1,0,0);
    //one_camera = *new Camera(0.6,1.6,-8.5,0,1,0,0);
    
    //one_camera = *new Camera(1.0,1.6,-5.5,0.4,1,0.6,0);
    //one_camera = *new Camera(1.1,1.8,-5.5,0.1,1.2,0.6,0);
    //one_camera = *new Camera(1.3,1.8,-0.5,0,1,0,0);
    //one_camera = *new Camera(0.6,1.8,-1,0.3,1,0.5,0);
    //one_camera = *new Camera(1.0,1.6,-5.5,0,1,0,0);
    //one_camera = *new Camera(1.2,1.8,-5.1,0.1,1,0.2,0.2);
    //one_camera = *new Camera(1.9,1.9,-5.7,0.1,1,-0.2,0.2);
    //one_camera = *new Camera(1.9,1.9,-8.7,0.1,1,-0.1,0.1);
    //one_camera = *new Camera(0.6,1.9,-6.7,0.1,1.2,-0.1,0.1);
    //one_camera = *new Camera(0.2,1.9,-5.8,-0.1,1.0,-0.1,0.1);
    //one_camera = *new Camera(2.3,2.2,-3.7,0.1,1,-0.2,0.2);
    //one_camera = *new Camera(2.3,2.5,-1.7,0.1,1.2,0.1,0.1);
    //one_camera = *new Camera(2.3,2.5,-9.7,0.1,1.2,0.1,0.1);
    
    //GOAL (1.0,1.6,-5,0.0,1,0.0,0.0);
    
    // Initialize camera
    ini(&one_camera);
    
    save_camera = one_camera; // update save_camera
}
void processing()
{
    Camera cameraPopulation[populationSize];
    double stepRatioPopulation[populationSize];
    double cv = 0; // Coefficient variant
    
    // Record camera previous camera rotation
    Quaternion<double> x;  // x is quaternion before update
    Quaternion<double> next_x;   // next_x is quaternion after update
    x = one_camera.quaternionCamera;
    
    // First time render and
    if (isFirstTime) {
        start_time = time(NULL);
        one_camera = inputCamera;
        display(&input_image); // Render depth image inside 3D model
        nonLinearTransform(&input_image, win.height, win.width);
        
        render(&show_input_image);
        parameterInitialize();
        
        // No more first time
        isFirstTime = false;
    }
    else {
        ///////////***********************///////////
        ///                                       ///
        ////////    Step ratio population    ////////
        ///                                       ///
        ///////////***********************///////////
        for (int i = 0; i<populationSize; i++) {
            stepRatioPopulation[i] = 0;
        }
        for (int i = 0; i<populationSize; i++) {
            //stepRatioPopulation[i] = stepRatio*exp(tau*rng.gaussian(1));
            stepRatioPopulation[i] = rng.gaussian(1);
        }
        ///////////***********************///////////
        ///                                       ///
        ////////    Initialize population    ////////
        ///                                       ///
        ///////////***********************///////////
        for (int i = 0; i<populationSize; i++) {
            cameraPopulation[i] = one_camera;
        }

        ///////////***********************///////////
        ///                                       ///
        ////////     Update population       ////////
        ///                                       ///
        ///////////***********************///////////
        // Random set camera poses
        for (int i = 0; i<populationSize; i++) {
            /* save the mutation vector */
            evolution_strategy(&cameraPopulation[i], &rng, sigma*sqrt(stepRatio*exp(tau*stepRatioPopulation[i])), sigma/sqrt(stepRatio*exp(tau*stepRatioPopulation[i])));
        }
    
        ///////////***********************///////////
        ///                                       ///
        ////////       Calculate SSD         ////////
        ///                                       ///
        ///////////***********************///////////
        save_camera = one_camera;
        double posSSD = 0;
        double rotSSD = 0;
        double sum = 0;
        for (int i = 0; i<populationSize; i++) {
            one_camera = cameraPopulation[i];
            
            display(&current_image);
            //render(&camera_image);
            nonLinearTransform(&current_image, win.height, win.width);
            cameraPopulation[i].ssd_score = SSD(input_image, current_image, win.width, win.height);
            sum += cameraPopulation[i].ssd_score;
            one_camera = save_camera;
        }
        double mean = sum/populationSize;
        double standardDeviation=0;
        for (int i=0; i<populationSize; i++) {
            standardDeviation += (cameraPopulation[i].ssd_score-mean)*(cameraPopulation[i].ssd_score-mean);
        }
        standardDeviation = sqrt(standardDeviation/(populationSize-1));
        cv = standardDeviation/mean;
        cvFile<<cv<<endl;
        meanFile<<mean<<endl;
        sdFile<<standardDeviation<<endl;
        
        one_camera = save_camera;
        ///////////***********************///////////
        ///                                       ///
        ////////       Find best SSDs        ////////
        ///       Sort the population array       ///
        ///////////***********************///////////
        for (int i = 0; i<selectSize; i++) {
            for (int j = i; j<populationSize; j++) {
                if (cameraPopulation[j].ssd_score<cameraPopulation[i].ssd_score) {
                    Camera replace = cameraPopulation[j];
                    cameraPopulation[j] = cameraPopulation[i];
                    cameraPopulation[i] = replace;
                    
                    double replaceRatio = stepRatioPopulation[j];
                    stepRatioPopulation[j] = stepRatioPopulation[i];
                    stepRatioPopulation[i] = replaceRatio;
                }
            }
        }
        double ratioSum = 0;
        for (int i = 0; i<selectSize; i++) {
            ratioSum += stepRatioPopulation[i];
        }
        //stepRatio = ratioSum/selectSize;
        stepRatio *= exp(tau*ratioSum/selectSize);
        //cout<<"step ratio : "<<stepRatio<<endl;
        ratioFile<<stepRatio<<endl;
        
        for (int i=0; i<7; i++) {
            mutationVector[i] = 0;
        }
        for (int i=0; i<selectSize; i++) {
            for (int j=0; j<7; j++) {
                mutationVector[j] += cameraPopulation[i].mutationVec[j];
            }
        }
        for (int i=0; i<7; i++) {
            mutationVector[i] /= selectSize;
        }
        
        ///////////***********************///////////
        ///                                       ///
        ////////       Update camera         ////////
        ///                                       ///
        ///////////***********************///////////
        /*
        one_camera.eye_x = one_camera.eye_x+mutationVector[0]*sigma*sqrt(stepRatio);
        one_camera.eye_y = one_camera.eye_y+mutationVector[1]*sigma*sqrt(stepRatio);
        one_camera.eye_z = one_camera.eye_z+mutationVector[2]*sigma*sqrt(stepRatio);
        one_camera.center_x = one_camera.center_x+mutationVector[0]*sigma*sqrt(stepRatio);
        one_camera.center_y = one_camera.center_y+mutationVector[1]*sigma*sqrt(stepRatio);
        one_camera.center_z = one_camera.center_z+mutationVector[2]*sigma*sqrt(stepRatio);
         */
        one_camera.eye_x = 0;
        one_camera.eye_y = 0;
        one_camera.eye_z = 0;
        one_camera.center_x = 0;
        one_camera.center_y = 0;
        one_camera.center_z = 0;
        for (int i=0; i<selectSize; i++) {
            one_camera.eye_x += cameraPopulation[i].eye_x;
            one_camera.eye_y += cameraPopulation[i].eye_y;
            one_camera.eye_z += cameraPopulation[i].eye_z;
            one_camera.center_x += cameraPopulation[i].center_x;
            one_camera.center_y += cameraPopulation[i].center_y;
            one_camera.center_z += cameraPopulation[i].center_z;
        }
        one_camera.eye_x /= selectSize;
        one_camera.eye_y /= selectSize;
        one_camera.eye_z /= selectSize;
        one_camera.center_x /= selectSize;
        one_camera.center_y /= selectSize;
        one_camera.center_z /= selectSize;
        
        // Update quaternion rotation
        // Using PPSN 2014
        double cosFactor;
        double sinFactor;
        double dotProduct;
        double rotationMutationVector[4] = {0};
        for (int i=0; i<selectSize; i++) {
            rotationMutationVector[0] += (sigma/sqrt(stepRatio*exp(tau*stepRatioPopulation[i])))*cameraPopulation[i].mutationVec[3];
            rotationMutationVector[1] += (sigma/sqrt(stepRatio*exp(tau*stepRatioPopulation[i])))*cameraPopulation[i].mutationVec[4];
            rotationMutationVector[2] += (sigma/sqrt(stepRatio*exp(tau*stepRatioPopulation[i])))*cameraPopulation[i].mutationVec[5];
            rotationMutationVector[3] += (sigma/sqrt(stepRatio*exp(tau*stepRatioPopulation[i])))*cameraPopulation[i].mutationVec[6];
        }
        for (int i=0; i<4; i++) {
            rotationMutationVector[i] /= selectSize;
        }
        
        dotProduct = rotationMutationVector[0]*rotationMutationVector[0]+rotationMutationVector[1]*rotationMutationVector[1]+rotationMutationVector[2]*rotationMutationVector[2]+rotationMutationVector[3]*rotationMutationVector[3];
        cosFactor = cos(sqrt(dotProduct));
        sinFactor = sin((sqrt(dotProduct)))/(sqrt(dotProduct));
        
        one_camera.quaternionCamera.x() = cosFactor * save_camera.quaternionCamera.x() + sinFactor * rotationMutationVector[0];
        one_camera.quaternionCamera.y() = cosFactor * save_camera.quaternionCamera.y() + sinFactor * rotationMutationVector[1];
        one_camera.quaternionCamera.z() = cosFactor * save_camera.quaternionCamera.z() + sinFactor * rotationMutationVector[2];
        one_camera.quaternionCamera.w() = cosFactor * save_camera.quaternionCamera.w() + sinFactor * rotationMutationVector[3];
        one_camera.quaternionCamera.normalize();
        
        // Rotation matrix
        MatrixXd m(3,3);
        m = one_camera.quaternionCamera.toRotationMatrix();
        
        // Update center point
        Vector3d directionVector;
        directionVector<<0,0,-1;
        Vector3d b;
        b = m * directionVector;
        
        one_camera.center_x = one_camera.eye_x + b(0);
        one_camera.center_y = one_camera.eye_y + b(1);
        one_camera.center_z = one_camera.eye_z + b(2);
        
        // Update up vector
        directionVector<<0,1,0;
        Vector3d c;
        c = m * directionVector;
        one_camera.up_x = c(0);
        one_camera.up_y = c(1);
        one_camera.up_z = c(2);
    }
    ///////////***********************///////////
    ///                                       ///
    ////////       Scoring function      ////////
    ///                                       ///
    ///////////***********************///////////
    // Calculate SSD score
    double current_ssd;
    display(&current_image);
    nonLinearTransform(&current_image, win.height, win.width);
    current_ssd = SSD(input_image, current_image, win.width, win.height);
    //cout<<"current ssd : "<<current_ssd<<endl;
    ///////////***********************///////////
    ///                                       ///
    ////////      Step size control      ////////
    ///                                       ///
    ///////////***********************///////////
    // Calculate the evolution paths
    // s <-- (1-c)*s+sqrt(mu*c*(2-c))*z
    // calculate z
    // Calculate path with 7 dimensions.
    const double c_path = 1.0/sqrt(6.0);
    const double mu = selectSize;
    // use the mutation vector
    stepSizePath(0) = (1.0-c_path)*(stepSizePath(0)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[0]);
    stepSizePath(1) = (1.0-c_path)*(stepSizePath(1)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[1]);
    stepSizePath(2) = (1.0-c_path)*(stepSizePath(2)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[2]);
    
    stepSizePath(3) = (1.0-c_path)*(stepSizePath(3)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[3]);
    stepSizePath(4) = (1.0-c_path)*(stepSizePath(4)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[4]);
    stepSizePath(5) = (1.0-c_path)*(stepSizePath(5)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[5]);
    stepSizePath(6) = (1.0-c_path)*(stepSizePath(6)) + sqrt(c_path*(2.0-c_path)*mu)*(mutationVector[6]);
    
    // PPSN 2014
    // Transplant the path w (update w)
    next_x = one_camera.quaternionCamera;  // next_x is updated quaternion x'
    Quaternion<double> v;   // v in formula (5)
    double v_scale;
    v_scale = x.dot(next_x);    // <x, x'> in formula (4)
    
    Quaternion<double> v_upper;     // upper part in formula (5)
    v_upper.x() = next_x.x() - v_scale * x.x();
    v_upper.y() = next_x.y() - v_scale * x.y();
    v_upper.z() = next_x.z() - v_scale * x.z();
    v_upper.w() = next_x.w() - v_scale * x.w();
    
    double v_bottom;    // bottom part in formula (5)
    v_bottom = sqrt(v_upper.dot(v_upper));
    v.x() = v_upper.x()/v_bottom;
    v.y() = v_upper.y()/v_bottom;
    v.z() = v_upper.z()/v_bottom;
    v.w() = v_upper.w()/v_bottom;
    
    // Update rotation searching path using parallel transport
    double qdotv = stepSizePath(3)*v.x() + stepSizePath(4)*v.y() + stepSizePath(5)*v.z() + stepSizePath(6)*v.w();
    stepSizePath(3) = stepSizePath(3) - qdotv * ((1.0-v_scale)*v.x() + v.dot(next_x)*x.x());
    stepSizePath(4) = stepSizePath(4) - qdotv * ((1.0-v_scale)*v.y() + v.dot(next_x)*x.y());
    stepSizePath(5) = stepSizePath(5) - qdotv * ((1.0-v_scale)*v.z() + v.dot(next_x)*x.z());
    stepSizePath(6) = stepSizePath(6) - qdotv * ((1.0-v_scale)*v.w() + v.dot(next_x)*x.w());
    
    // calculate & update step size
    sigma = sigma * exp(((stepSizePath.dot(stepSizePath)) - 6)/(2*10*6));
    positionStepFactor = sigma*sqrt(stepRatio);
    rotationStepFactor = sigma/sqrt(stepRatio);
    //cout<<positionStepFactor<<endl;
    //cout<<"position size : "<<positionStepFactor<<endl;
    //cout<<"rotation size : "<<rotationStepFactor<<endl;
    if (rotationStepFactor>0.5) {
        rotationStepFactor=0.5;
    }
    sigmaLocFile[process_num]<<positionStepFactor<<endl;
    sigmaRotFile[process_num]<<rotationStepFactor<<endl;
    locDifFile[process_num]<<eyeDistanceWrite()<<endl;
    rotDifFile[process_num]<<rotationDistanceWrite()<<endl;
    ///////////***********************///////////
    ///                                       ///
    ////////      Best SSD update        ////////
    ///                                       ///
    ///////////***********************///////////
    if (current_ssd < previous_ssd) {
        previous_ssd = current_ssd; // Update local best SSD
        stopCounter=0;
    }
    else{
        // Use previous camera for next
        betterSSD = false;
        stopCounter++;
    }
    ssdValueFile[process_num]<<current_ssd<<endl;
    counting++;
    subCounting++;
    
    // stopping criteria
    const int subCountingThreshold = 500;
    const int countingThreshold = 500;
    const double cvThreshold = 0.03;
    //if ((eyeDistanceWrite()<0.05 && rotationDistanceWrite()<2*PI/180) || (counting>countingThreshold) || (cv<cvThreshold&&subCounting>subCountingThreshold)) {
    if ((eyeDistanceWrite()<0.05 && rotationDistanceWrite()<2*PI/180) || (counting>countingThreshold)) {
        cout<<"processing time : "<<time(NULL)-start_time<<endl;
        cout<<"best ssd : "<<previous_ssd<<endl;
        cout<<"counting number : "<<counting<<endl;
        cout<<"process num : "<<process_num<<endl;
        cout<<"subCounting : "<<subCounting<<endl;
        cout<<"sub process num : "<<subProcessNum<<endl;
        if (counting>countingThreshold) {
            successRateFile<<"0"<<endl;
            process_num++;
            iterationCountFile<<counting<<endl;
            posDisFile<<eyeDistanceWrite()<<endl;
            rotDisFile<<rotationDistanceWrite()<<endl;
            subProcessNumFile<<subProcessNum<<endl;
            subProcessNum = 0;
            subCounting = 0;
            counting = 0;
        }
        /*
        // Start again
        else if (cv<cvThreshold&&subCounting>subCountingThreshold){
            subProcessNum++;
            subCounting = 0;
            cout<<"one sub finish"<<endl;
        }
         */
        else if (eyeDistanceWrite()<0.05 && rotationDistanceWrite()<2*PI/180){
            successRateFile<<"1"<<endl;
            process_num++;
            iterationCountFile<<counting<<endl;
            posDisFile<<eyeDistanceWrite()<<endl;
            rotDisFile<<rotationDistanceWrite()<<endl;
            subProcessNumFile<<subProcessNum<<endl;
            subCounting = 0;
            counting = 0;
            subProcessNum = 0;
        }
        isFirstTime = true;
        //ssdValueFile<<previous_ssd<<endl;
    }
    if (process_num>19) {
        exit(0);
    }
    
    ///////////***********************///////////
    ///                                       ///
    ////////     Repeat this process     ////////
    ///                                       ///
    ///////////***********************///////////
    glutPostRedisplay();
}
void keyboard ( unsigned char key, int x, int y )
{
    switch ( key ) {
        case KEY_ESCAPE:
            save_screenshot("/Users/luang/Desktop/Interpolate.jpg", win.width, win.height);
            cout<<"saved"<<endl;
            cout<<"GOOD"<<endl;
            one_camera.show();
            exit ( 0 );
            break;

        default:
            break;
    }
}
int main(int argc, char **argv)
{
    // Set window parameters
    win.width = 160;
    win.height = 120;
    win.title = "Match and Search";
    win.field_of_view_angle = 45;
    win.z_near = 1.0f;
    win.z_far = 500.0f;
    
    input_image = new float[win.width * win.height];
    current_image = new float[win.width * win.height];
    squar_diff_image = new float[win.width * win.height];
    
    ratioFile.open("ratio.txt");
    rotSSDFile.open("rotSSDFile.txt");
    posSSDFile.open("posSSDFile.txt");
    successRateFile.open("successRate.txt");
    iterationCountFile.open("iterationCount.txt");
    posDisFile.open("posDis.txt");
    rotDisFile.open("rotDis.txt");
    cvFile.open("coefficentVariance.txt");
    meanFile.open("mean.txt");
    sdFile.open("standardDeviation.txt");
    subProcessNumFile.open("subPfocessNum.txt");
    
    stringstream sstm;
    for (int i=0;i<20 ;i++)
    {
        sstm.str("");
        sstm << "ssdValue" << i<<".txt";
        ssdValueFile[i].open(sstm.str());
    }
    for (int i=0;i<20 ;i++)
    {
        sstm.str("");
        sstm << "sigmaLoc" << i<<".txt";
        sigmaLocFile[i].open(sstm.str());
    }
    for (int i=0;i<20 ;i++)
    {
        sstm.str("");
        sstm << "sigmaRot" << i<<".txt";
        sigmaRotFile[i].open(sstm.str());
    }
    for (int i=0;i<20 ;i++)
    {
        sstm.str("");
        sstm << "locDif" << i<<".txt";
        locDifFile[i].open(sstm.str());
    }
    for (int i=0;i<20 ;i++)
    {
        sstm.str("");
        sstm << "rotDif" << i<<".txt";
        rotDifFile[i].open(sstm.str());
    }
    
    // Start processing
    glutInit(&argc, argv);
    glutInitWindowSize(win.width,win.height);					// set window size
    glutCreateWindow(win.title);								// create Window
    glutDisplayFunc(processing);								// register Display Function
    glutKeyboardFunc(keyboard);
    obj.Load("/Users/luang/Desktop/Research/ThesisProject/Model/Maya_Classroom.obj");        // load 3D models
    //obj.Load("/Users/luang/Desktop/Research/ThesisProject/Model/familia_blender2.obj");        // load 3D models
    glEnable( GL_DEPTH_TEST );
    glutMainLoop();												// run GLUT mainloop
    return 0;
}