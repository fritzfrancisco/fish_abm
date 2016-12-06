#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

void paint(int x, int y, Mat image);

int main(){

  // Attraction heatmap:
  char heatmap_window[] = "Attraction Heatmap";
  Mat heatmap_image_grey = Mat::zeros(410,510,CV_8UC3);
  Mat heatmap_image_color;

  Point center(205,205);
  // Point corner_1(455,5);
  // Point corner_2(475,405);

  for(int x = 0;x<410;x++){
    for(int y = 0; y < 410;y++){
      paint(x,y,heatmap_image_grey);
    }
  }

  applyColorMap(heatmap_image_grey,heatmap_image_color,11);
  circle(heatmap_image_color,center,15,(0,0,0),1,CV_AA,0);
  circle(heatmap_image_color,center,100,(0,0,0),1,CV_AA,0);
  circle(heatmap_image_color,center,200,(0,0,0),1,CV_AA,0);
  // rectangle(heatmap_image_color,corner_1,corner_2,(255,255,255),CV_FILLED,0);

  imshow(heatmap_window,heatmap_image_color);

  waitKey();

  return 0;
}

void paint(int x, int y,Mat image){
  double s = 0;
  double d = sqrt(pow((205-x),2)+pow((205-y),2));
  if(d < 15){
    s = 1 / ( 0.01 * pow((d),2) + 4 );
    image.at<Vec3b>(x,y) = Vec3b(127 +(255*s),127+(255*s),127+(255*s));
  }
  else if(d >= 15 && d < 100){
    s = 1 / (0.01 * pow((d-15),2) + 2);
    image.at<Vec3b>(x,y) =  Vec3b(127 +(255*s),127+(255*s),127+(255*s));
  }
  else if(d >= 100 && d < 200){
    s = (1 / ( 0.01 * pow((d-200),2) + 2 ));
    image.at<Vec3b>(x,y) =  Vec3b(127 +(255*s),127+(255*s),127+(255*s));
  }
}
