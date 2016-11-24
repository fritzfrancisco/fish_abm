#include <iostream>
#include <stdlib.h>
#include <ctime>

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
  double lag [30]; // length n
  double coords [30][2]; // height n * width dim
  bool state [30]; // state of movement (directed/random)

};


void move(individuals& inds);
void printcoords(individuals inds);
void initialize(individuals& inds,int n, int dim);
void drawfish(individuals& inds, Mat fish_image);
int follow(individuals inds,int id);


// int colors[4][3] = {255,0,0,0,255,0,0,0,255,120,0,120};


int main(){

  srand(time(0));

  // Creating window for imaging
  char fish_window[] = "Fish";
  Mat fish_image = Mat::zeros(w,w,CV_8UC3);


  individuals fish;
  initialize(fish,30,2);

  // printcoords(fish);
  // printcoords(fish);

  for(int z = 0; z < 1000; z++){
  // Creating black empty window
  // if(z % 4 == 0){
  //   Mat fish_image = Mat::zeros(w, w,CV_8UC3);
  // }
  Mat fish_image = Mat::zeros(w,w,CV_8UC3);

  drawfish(fish,fish_image);
  move(fish);

  imshow(fish_window,fish_image);
  waitKey(30);
  }
  return(0);
}

void move(individuals& inds){
  for(int i = 0; i < inds.n; i++){
    if(follow(inds,i) > 0 && inds.lag[i] == 0){

      inds.state[i] = 1;
      // inds.coords[i][0] = inds.coords[i][0] + 2;
    }

    else{

      inds.state[i] = 0;
    }
  }

  for(int i = 0; i < inds.n;i++){

    if(inds.state[i] == 1 && inds.lag[i] == 0){
      inds.coords[i][0] = inds.coords[i][0] + 15;
      inds.coords[i][1] = inds.coords[i][1] + 15;
      // inds.coords[i][0] = inds.coords[i][0] + rand()%15;
      // inds.coords[i][1] = inds.coords[i][1] + rand()%15;
    }

    else{

      if(inds.lag[i] == 0){
        inds.lag[i] = 10;
      }
      else{
      inds.lag[i] = inds.lag[i]-1;
      }

      for(int u = 0; u < inds.dim; u++){
        inds.coords[i][u] = inds.coords[i][u] + rand()%5 - 2;
      }
    }
  }
}

void printcoords(individuals inds){
  for(int i = 0; i < inds.n; i++){
    cout << inds.coords[i][0] << "," << inds.coords[i][1] << "\n";
  }
}

void initialize(individuals& inds, int n ,int dim){
  inds.n = n;
  inds.dim = dim;

  for(int i = 0; i < inds.n; i++){
    inds.lag[i] = 0;
    for(int u = 0; u < inds.dim; u++){
      inds.coords[i][u] = 10;
    }
  }
}

void drawfish(individuals& inds, Mat fish_image){
  for(int i = 0; i < inds.n; i++){
    // Scalar fishcol(100,inds.coords[i][1],inds.coords[i][0]);
    Scalar fishcol(0,0,255);
    Point pos(inds.coords[i][0],inds.coords[i][1]);
    circle(fish_image,pos,2,fishcol,1,20,0);
  }
}

int follow(individuals inds,int id){
  // count number of inds in front (x-axis) of id individual
  int c = 0;
  for (int i = 0; i < inds.n;i++){
    if (inds.coords[i][0] > inds.coords[id][0]){
      c++;
    }
  }
  return c;
}
