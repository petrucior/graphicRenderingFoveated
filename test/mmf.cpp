#include <stdio.h>
#include <iostream>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/xfeatures2d/nonfree.hpp"
#include "../mmf.h"

using namespace std;
using namespace cv;

int main(int argc, char** argv){
  
  Mat image = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  if ( image.empty() ){
    printf("Can't read one of the images\n");
    return -1;
  }

  MMF mmf;
  
  //Mat imgFoveated = mmf.MMF_CPU( image, 4, 4, Point(40, 40), Point(image.cols, image.rows), Point(0, 0) );
  Mat imgFoveated = mmf.foveated( image, 1000, Point(60, 60), Point(image.cols, image.rows), Point(0, 0), 0 );

  imshow("Foveated image", imgFoveated);
  waitKey(0);
  
  return 0;
}
