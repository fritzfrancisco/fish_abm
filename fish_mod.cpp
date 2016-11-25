#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <math.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"
#include <opencv2/imgproc.hpp>
#include "opencv2/core/core_c.h"
#include "opencv2/highgui/highgui_c.h"

using namespace std;
using namespace cv;

#define w 550

class individuals{

public:

  // int n = 3;
  // int dimensions = 2;
  int n;
  int dim;
  double lag [40]; // length n
  double coords [40][2]; // height n * width dim
  bool state [40]; // state of movement (directed/random)
  double dir [40]; // array of agles between 0 - 360 degrees for each individual
  int colors[40][3]; // individual color in RGB

};


void move(individuals& inds);
void printcoords(individuals inds);
void initialize(individuals& inds,int n, int dim);
void drawfish(individuals& inds, Mat fish_image);
int follow(individuals inds,int id);

int main(){

  srand(time(0));

  individuals fish;
  initialize(fish,40,2);

  // Creating window for imaging
  char fish_window[] = "Fish";
  Mat fish_image = Mat::zeros(w,2*w,CV_8UC3);

  // printcoords(fish);
  // printcoords(fish);

  for(int z = 0; z < 100; z++){
  // Creating black empty window
  // if(z % 4 == 0){
  //   Mat fish_image = Mat::zeros(w, w,CV_8UC3);
  // }
  Mat fish_image = Mat::zeros(w,2*w,CV_8UC3);
  drawfish(fish,fish_image);
  move(fish);

  imshow(fish_window,fish_image);
  waitKey(100);
  }
  return(0);
}


void move(individuals& inds){
  for(int i = 0; i < inds.n; i++){
    if(follow(inds,i) > 20 && inds.lag[i] == 0){

      inds.state[i] = 1;
      // inds.coords[i][0] = inds.coords[i][0] + 2;
    }

    else{

      inds.state[i] = 0;
    }
  }

  for(int i = 0; i < inds.n;i++){

    if(inds.state[i] == 1 && inds.lag[i] == 0){
      inds.coords[i][0] = inds.coords[i][0] + 40; // swim distance
      // inds.coords[i][1] = inds.coords[i][1] + 10;
      // inds.coords[i][0] = inds.coords[i][0] + rand()%15;
      // inds.coords[i][1] = inds.coords[i][1] + rand()%15;
    }

    else{

      if(inds.lag[i] == 0){
        inds.lag[i] = 20;
      }
      else{
      inds.lag[i] = inds.lag[i]-1;
      }

      for(int u = 0; u < inds.dim; u++){
        inds.coords[i][u] = inds.coords[i][u] + rand()%5 - 2; // initialize random walk
      }
    }
  }
}

void printcoords(individuals inds){
  for(int i = 0; i < inds.n; i++){
  }
}

void initialize(individuals& inds, int n ,int dim){
  inds.n = n;
  inds.dim = dim;

  for(int i = 0; i < inds.n; i++){

    inds.lag[i] = 0;
    inds.colors[i][0] = 70 + rand()%155;
    inds.colors[i][1] = 70 + rand()%155;
    inds.colors[i][2] = 70 + rand()%155;
    inds.dir[i] = 0;
    // inds.dir[i] = rand() % 360;

    for(int u = 0; u < inds.dim; u++){
      if(u == 0){
        inds.coords[i][u] = rand()%100;
      }
      else{
        inds.coords[i][u] = w/2 + rand()%100;
        // inds.coords[i][u] = rand()%100;
      }
    }
  }
}

void drawfish(individuals& inds, Mat fish_image){
  for(int i = 0; i < inds.n; i++){
    // Scalar fishcol(100,inds.coords[i][1],inds.coords[i][0]);
    // Scalar fishcol(0,0,255);
    Scalar fishcol(inds.colors[i][0],inds.colors[i][1],inds.colors[i][2]);
    Point center(inds.coords[i][0],inds.coords[i][1]);
    Point tail(inds.coords[i][0] - 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] - 10 * sin(inds.dir[i] * M_PI / 180));
    Point head(inds.coords[i][0] + 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] + 10 * sin(inds.dir[i] * M_PI / 180));
    arrowedLine(fish_image,tail,head,fishcol,1,20,0,0.1);
    circle(fish_image,center,2,fishcol,1,20,0);
  }
}

int follow(individuals inds,int id){
  // count number of inds in front (x-axis) of id individual
  int c = 0;

  double idcoords[2] = {inds.coords[id][0],inds.coords[id][1]}; // includes corrected coordinates for each fish per fish [id]

  for (int i = 0; i < inds.n;i++){

    double dist[40]; // absolute distance between fish [id] and fish [i]
    double angle[40]; // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis
    double dirvector[40][2]; // angle of fish[id] is facing in respect to original x-axis

    inds.coords[i][0] = inds.coords[i][0] - idcoords[0]; // set fish id x to 0
    inds.coords[i][1] = inds.coords[i][1] - idcoords[1]; // set fish id y to 0

    dist[i] = sqrt(pow(inds.coords[i][0],2) + pow(inds.coords[i][1],2));

    dirvector[i][0] = cos(inds.dir[i] * M_PI / 180);
    dirvector[i][1] = sin(inds.dir[i] * M_PI / 180);

    if(i != id){
      angle[i] = acos((dirvector[i][0]*inds.coords[i][0]+dirvector[i][1]*inds.coords[i][1])/dist[i])*180/M_PI;
    }
    else{
      angle[i] = 180;
    }

    if (angle[i] < 90){
      c++;
      }
    }
  return c;
}
