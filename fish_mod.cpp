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

#define w 700
#define h 500
#define num 100

class individuals{

public:
  // int n = 3;
  // int dimensions = 2;
  int n;
  int dim;
  int colors[num][3];     // individual color in RGB
  int nfollow;

  double dspeed;          // dspeed+1 = lag-phase / swim speed
  double speed;           // movement speed
  double lag [num];       // length n
  double coords [num][2]; // height n * width dim
  double dir [num];       // array of agles between 0 - 360 degrees for each individual

  bool state [num];       // state of movement (directed/random)
};

int follow(individuals inds,int id, int fov); // fov: field of view as angle

void move(individuals& inds);
void printcoords(individuals inds);
void initialize(individuals& inds,int n, int dim, double speed, double dspeed, int nfollow);
void drawfish(individuals& inds, Mat fish_image);

double avoidance(individuals inds,int id,int fov);
double collision(individuals inds, int id, double newdir);

bool state(individuals inds, int id);

int main(){

  srand(time(0));

  individuals fish;
  initialize(fish,num,2,20,4,50); // initialize: object fish,num individuals, dimensions,speed,dspeed,n individuals to follow

  // Creating window for imaging
  char fish_window[] = "Fish";
  Mat fish_image = Mat::zeros(h,w,CV_8UC3);

  for(int z = 0; z < 1000; z++){
  Mat fish_image = Mat::zeros(h,w,CV_8UC3);

  move(fish);
  drawfish(fish,fish_image);
  imshow(fish_window,fish_image);
  waitKey(30);
  }

  return(0);
}


void move(individuals& inds){
  for(int id = 0; id < inds.n;id++){

    double dir = inds.dir[id];
    inds.state[id] = state(inds, id);

    if(inds.state[id] == 0){  // feeding, initialize random walk
      dir = dir + rand()%90+315;

      if(inds.lag[id] == 0){  // individual lag phase
        inds.lag[id] = rand()%40;
      }
      else{
      inds.lag[id] = inds.lag[id]-1;
      }
    }
    double newdir = avoidance(inds,id,330);
    newdir = collision(inds,id,dir);

    inds.coords[id][0] =  inds.coords[id][0] + (inds.state[id]*inds.dspeed+1) * inds.speed * cos(newdir * M_PI / 180);
    inds.coords[id][1] = inds.coords[id][1] + (inds.state[id]*inds.dspeed+1) * inds.speed * sin(newdir * M_PI/ 180);
    inds.dir[id] = newdir;
  }
}

void printcoords(individuals inds){
  for(int i = 0; i < inds.n; i++){
  }
}

void initialize(individuals& inds, int n ,int dim, double speed, double dspeed, int nfollow){
  inds.n = n;
  inds.dim = dim;
  inds.speed = speed;
  inds.dspeed = dspeed;
  inds.nfollow = nfollow;

  for(int i = 0; i < inds.n; i++){

    inds.lag[i] = 0;
    inds.colors[i][0] = 75 + rand()%155;
    inds.colors[i][1] = 75 + rand()%155;
    inds.colors[i][2] = 75 + rand()%155;
    inds.dir[i] = 0;
    // inds.dir[i] = rand()%361; // all individuals facing x-axis
    // inds.dir[i] = rand() % 20 + 350; // field of vision

    for(int u = 0; u < inds.dim; u++){
      if(u == 0){
        inds.coords[i][u] = w/4 + rand()%100-50;
      }
      else{
        inds.coords[i][u] = h/2 + rand()%100-50;
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

    arrowedLine(fish_image,tail,head,fishcol,1,20,0,0.1); // draw fish
    circle(fish_image,center,2,fishcol,1,20,0);
  }
}

int follow(individuals inds,int id,int fov){
  // count number of inds in front (x-axis) of id individual
  int c = 0;
  // double idcoords[2] = {inds.coords[id][0],inds.coords[id][1]}; // includes corrected coordinates for each fish per fish [id]
  double dirvector[inds.dim];   // angle of fish[id] is facing in respect to original x-axis
  double dist[inds.n];          // absolute distance between fish [id] and fish [i]
  double angle[inds.n];         // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

  for (int i = 0; i < inds.n;i++){

    double univec[inds.n][inds.dim];

    inds.coords[i][0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id x to 0
    inds.coords[i][1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id y to 0

    dist[i] = sqrt(pow(inds.coords[i][0],2) + pow(inds.coords[i][1],2)); // absolute value

    // for(int u = 0; u < inds.dim-1; u++){ // creating unit vector
    // univec[i][u] = inds.coords[i][u]/dist[i];
    // }

    if(i != id){
      angle[i] = acos((dirvector[0]*inds.coords[i][0]+dirvector[1]*inds.coords[i][1])/dist[i])*180/M_PI;
    }

    else{
      angle[i] = 180; // never see yourself
    }

    if (angle[i] < (fov/2) && dist[i] < 200){ // visual field with angle and radius
      c++;
      }
    }
  return c;
}

bool state(individuals inds, int id){
  bool state = 0;
  if(follow(inds,id,180) >= inds.nfollow && inds.lag[id] == 0){
    state = 1;
  }
  return state;
}

double collision(individuals inds, int id, double newdir){
	double bound = 20;     // boundary in which to detect the border
	double leftcoords[2];  // new possible coords after left turn
	double rightcoords[2]; // new possible coords after right turn

	leftcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdir * M_PI / 180); // initialize new coords after turn with original values
	leftcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdir * M_PI / 180);

	rightcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdir * M_PI / 180);
	rightcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdir * M_PI / 180);

	double newdirleft = newdir; // new direction after turning (left/right), initialized with original direction
	double newdirright = newdir;

	// check collision after turn (left or right, 0 for no collision, 1 for collision); starting with original direction..
	bool turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
	bool turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));

	// .. then adjust turn direction and check collision again; repeat until one collision statement (left/right) is 0 (no collision):
	while(turnleft && turnright){
		newdirleft = newdirleft - 1;
		leftcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirleft * M_PI / 180); // new x coord after left turn
		leftcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirleft * M_PI / 180); // new y coord

		newdirright = newdirright + 1;
		rightcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirright * M_PI / 180); // new x coord after right turn
		rightcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirright * M_PI / 180); // new y coord

		turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
		turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));
	}

	// choose new direction based on false collision statement (no collision in this direction):
	if (turnleft == 0){
		newdir = newdirleft;
	}
	else{
		newdir = newdirright;
	}
	return newdir;
}

double avoidance(individuals inds,int id,int fov){

  int c = follow(inds,id,330);  // counter for fish in field of vision

  double newdir;                // angle after avoidance decision
  double repulsion[inds.n];     // state of avoidance (0 = swim away, 1 = swim towards group)
  double dirvector[inds.dim];   // angle of fish[id] is facing in respect to original x-axis
  double dist[inds.n];          // absolute distance between fish [id] and fish [i]
  double angle[inds.n];         // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis
  double coordsum[2] ={0,0};
  bool avoid = 0;               // switch for avoidance behaviour (0 = no avoidance)

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

  for(int i = 0; i < inds.n;i++){

    repulsion[i] = 0;
    inds.coords[i][0] = inds.coords[i][0] - inds.coords[id][0];          // set fish id x to 0
    inds.coords[i][1] = inds.coords[i][1] - inds.coords[id][1];          // set fish id y to 0
    dist[i] = sqrt(pow(inds.coords[i][0],2) + pow(inds.coords[i][1],2)); // absolute value

    if(i != id){
      angle[i] = acos((dirvector[0]*inds.coords[i][0]+dirvector[1]*inds.coords[i][1])/dist[i])*180/M_PI;
    }

    else{
      angle[i] = 180; // never see yourself
    }

    if(angle[i]< fov/2 && dist[i] < 200){
      repulsion[i] = 1 / (0.01*pow(dist[i],2)+1); // weighted repulsion impact of fish in proximity (200) (f(x) = 1/(0.01x^2+1))
    }

    if(angle[i]< fov/2 && dist[i] < 20){
      avoid = 1;
    }

    if(i != id){
    coordsum[0] = coordsum[0] + inds.coords[i][0] * repulsion[i];
    coordsum[1] = coordsum[1] + inds.coords[i][1] * repulsion[i];
    }
  }
  coordsum[0] = coordsum[0]/c + inds.coords[id][0];
  coordsum[1] = coordsum[1]/c + inds.coords[id][1];

  double leftcoords[2];  // new possible coords after left turn
  double rightcoords[2]; // new possible coords after right turn

  leftcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(inds.dir[id] - 90 * M_PI / 180); // new x coord after left turn of 90 degrees
  leftcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(inds.dir[id] - 90 * M_PI / 180); // new y coord

  rightcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(inds.dir[id] - 90 * M_PI / 180); // new x coord after right turn of 90 degrees
  rightcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(inds.dir[id] - 90 * M_PI / 180); // new y coord

  double avoiddistl = sqrt(pow(leftcoords[0] - coordsum[0],2) + pow(leftcoords[1] - coordsum[1],2));
  double avoiddistr = sqrt(pow(rightcoords[0] - coordsum[0],2) + pow(rightcoords[1] - coordsum[1],2));

  if(avoiddistl > avoiddistr){
    newdir = inds.dir[id] - 90;
  }
  else if(avoiddistl == avoiddistr){
    newdir = inds.dir[id] + rand()%181 - 90;
  }
  else{
    newdir = inds.dir[id] + 90;
  }
  return newdir;
}
