#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <random>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

void paint(int x, int y, Mat image);

int main(){

  // Attraction heatmap:
  char heatmap_window[] = "Social heat map";
  Mat heatmap_image_grey = Mat::zeros(460,510,CV_8UC3); // (y,x,CV_8UC3)
  Mat heatmap_image_color;
  Mat invert = Mat::ones(460,510,CV_8UC3);
  invert = Scalar(255,255,255);
  Mat invert_out = Mat::zeros(460,510,CV_8UC3);

  Point center(205,205);
  // Point outer(207-cos(0.25*M_PI)*200,207 - sin(0.25*M_PI)*200);
  Point outer(207-cos(M_PI)*200,207 - sin(M_PI)*200);


  for(int x = 0;x<410;x++){
    for(int y = 0; y < 410;y++){
      paint(y,x,heatmap_image_grey);
    }
  }

  for(int x = 435;x < 455; x++){
    double s = 0.5;
    for(int y = 35 ;y < 436;y++){
        // double d = abs(y-405)/2;
        // s = 1 / ( 0.01 * pow((d-200),2) + 2);
        heatmap_image_grey.at<Vec3b>(y,x) = Vec3b((510*s),(510*s),(510*s));
        s = s - (0.5/400);
    }
  }

  Point triangle[1][3];
  triangle[0][0] = Point(205,205);
  triangle[0][1] = Point(271,405);
  triangle[0][2] = Point(139,405);

  const Point* ppt[1] = { triangle[0] };
  int npt[] = { 3 };


  Scalar w(255,255,255);
  Scalar s(0,0,0);

  // draw shapes:
  fillPoly(heatmap_image_grey,ppt,npt,1,Scalar(0,0,0),CV_AA);
  circle(heatmap_image_grey,center,15,w,1,CV_AA,0);
  circle(heatmap_image_grey,center,100,w,1,CV_AA,0);
  circle(heatmap_image_grey,center,200,w,1,CV_AA,0);
  line(heatmap_image_grey,Point(205,220),Point(205,250),w,1,CV_AA);
  line(heatmap_image_grey,Point(205,305),Point(205,335),w,1,CV_AA);
  line(heatmap_image_grey,Point(205,405),Point(205,435),w,1,CV_AA);
  // fish
  // circle(heatmap_image_grey,Point(40,430),15,w,1,CV_AA,0);
  // circle(heatmap_image_grey,Point(40,430),100,w,1,CV_AA,0);
  circle(heatmap_image_grey,Point(40,425),2,w,1,CV_AA,0);
  arrowedLine(heatmap_image_grey,Point(40,435),Point(40,415),w,1,CV_AA,0,0.2);

  // annotation
  putText(heatmap_image_grey,"0.75",Point(188,265),1,0.9,w,1,CV_AA);
  putText(heatmap_image_grey,"BL",Point(198,280),1,0.9,w,1,CV_AA);
  putText(heatmap_image_grey,"5 BL",Point(193,350),1,0.9,w,1,CV_AA);
  putText(heatmap_image_grey,"10 BL",Point(184,450),1,0.9,w,1,CV_AA);
  putText(heatmap_image_grey,"Body Length",Point(50,430),1,0.9,w,1,CV_AA);
  putText(heatmap_image_grey,"Social Force",Point(400,20),1,0.9,w,1,CV_AA);

  for(double i = 0; i < 6;i++){
    Point pt1(435, 35 + i *80);
    Point pt2(455, 35 + i * 80);
    Point pt3(460, 40 + i * 80);
    line(heatmap_image_grey,pt1,pt2,s,1,8,0);
    char s_value[10];
    sprintf(s_value, "%g", 0.5-i/10);

    putText(heatmap_image_grey,s_value,pt3,1,0.9,w,1,CV_AA,false);
}
  applyColorMap(heatmap_image_grey,heatmap_image_color,11);

  subtract(invert,heatmap_image_grey,invert_out);

  imshow(heatmap_window,invert_out);

  // imwrite( "Social heat map.png", invert_out);

  waitKey();

  return 0;
}

void paint(int x, int y,Mat image){
  double s = 0;
  double d = sqrt(pow((205-x),2)+pow((205-y),2));
  if(d < 15){
    s = 1 / ( 0.1 * pow((d),2) + 2 );
    image.at<Vec3b>(y,x) = Vec3b((510*s),(510*s),(510*s));
  }
  else if(d >= 15 && d < 100){
    s = 1 / (0.004 * pow((d-15),2) + 2);
    image.at<Vec3b>(y,x) = Vec3b((510*s),(510*s),(510*s));
  }
  else if(d >= 100 && d < 200){
    s = (1 / ( 0.004 * pow((d-200),2) + 2 ));
    image.at<Vec3b>(y,x) = Vec3b((510*s),(510*s),(510*s));
  }
}
