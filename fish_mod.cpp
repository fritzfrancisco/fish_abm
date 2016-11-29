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

#define w 1400
#define h 700
#define num 200

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
void initialize(individuals& inds,int n, int dim, double speed, double dspeed, int nfollow);
void drawfish(individuals& inds, Mat fish_image);

double social(individuals inds,int id,int fov);
double collision(individuals inds, int id, double newdir);

bool state(individuals inds, int id);

int main(){

  srand(time(0));

  // Creating window for imaging
  char fish_window[] = "Fish";
  Mat fish_image = Mat::zeros(h,w,CV_8UC3);

  individuals fish;
  initialize(fish,num,2,5,4,5); // initialize: object fish,num individuals, dimensions,speed,dspeed,n individuals to follow

  for(int z = 0; z < 10000; z++){
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
    double socialdir = social(inds,id,330);

    if(inds.state[id] == 0 && socialdir == inds.dir[id]){  // feeding, initialize random walk
      socialdir = socialdir + rand()%90+315;

      if(inds.lag[id] == 0){  // individual lag phase
        inds.lag[id] = rand()%5;
      }
      else{
      inds.lag[id] = inds.lag[id]-1;
      }
    }

    double coldir = collision(inds,id,socialdir);

    inds.coords[id][0] =  inds.coords[id][0] + (inds.state[id]*inds.dspeed+1) * inds.speed * cos(coldir * M_PI / 180);
    inds.coords[id][1] = inds.coords[id][1] + (inds.state[id]*inds.dspeed+1) * inds.speed * sin(coldir * M_PI/ 180);

    inds.dir[id] = coldir;
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
    inds.colors[i][0] = 80 + rand()%155;
    inds.colors[i][1] = 80 + rand()%155;
    inds.colors[i][2] = 80 + rand()%155;
    inds.dir[i] = 0;
    // inds.dir[i] = rand()%361; // all individuals facing x-axis
    // inds.dir[i] = rand() % 20 + 350; // field of vision

    for(int u = 0; u < inds.dim; u++){
      if(u == 0){
        inds.coords[i][u] = w/2 + rand()%201 - 100;
      }
      else{
        inds.coords[i][u] = h/2 + rand()%201 - 100;
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

    circle(fish_image,center,1,fishcol,1,CV_AA,0);            // visualize center
    // circle(fish_image,center,50,fishcol,1,CV_AA,0);           // visualize comfort zone
    // circle(fish_image,center,200,fishcol,1,CV_AA,0);       // visualize visual distance
    arrowedLine(fish_image,tail,head,fishcol,1,CV_AA,0,0.1);  // draw fish

  }
}

int follow(individuals inds,int id,int fov){
  // count number of inds in front (x-axis) of id individual
  int c = 0;
  // double idcoords[2] = {inds.coords[id][0],inds.coords[id][1]}; // includes corrected coordinates for each fish per fish [id]
  double dirvector[inds.dim];   // angle of fish[id] is facing in respect to original x-axis
  double dist;          // absolute distance between fish [id] and fish [i]
  double angle;         // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

  for (int i = 0; i < inds.n;i++){

    inds.coords[i][0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id x to 0
    inds.coords[i][1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id y to 0

    dist = sqrt(pow(inds.coords[i][0],2) + pow(inds.coords[i][1],2)); // absolute value

    if(i != id){
      angle = acos((dirvector[0]*inds.coords[i][0]+dirvector[1]*inds.coords[i][1])/dist)*180/M_PI;
    }

    else{
      angle = 180; // never see yourself
    }

    if (angle < (fov/2) && dist < 200){ // visual field with angle and radius
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
		newdirleft = newdirleft - 10;
		leftcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirleft * M_PI / 180); // new x coord after left turn
		leftcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirleft * M_PI / 180); // new y coord

		newdirright = newdirright + 10;
		rightcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirright * M_PI / 180); // new x coord after right turn
		rightcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirright * M_PI / 180); // new y coord

		turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
		turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));
	}

	// choose new direction based on false collision statement (no collision in this direction):
  if(turnleft == 0 && turnright == 0){
    if(rand()%2 == 0){
      newdir = newdirleft;
    }
    else{
      newdir = newdirright;
    }
  }
  else if (turnleft == 0 && turnright != 0){
		newdir = newdirleft;
	}
	else{
		newdir = newdirright;
	}
	return newdir;
}

double social(individuals inds,int id,int fov){ // function for social behaviour

int c = follow(inds,id,fov);                    // counter for fish in field of vision

  double newdir = inds.dir[id];                 // angle after decision
  double repulsion;                     // state of avoidance (0 = swim away, 1 = swim towards group)
  double dirvector[inds.dim];                   // angle of fish[id] is facing in respect to original x-axis
  double dist;                          // absolute distance between fish [id] and fish [i]
  double angle;                         // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis
  double coordsum[2] ={0,0};
  double nearestvisdist = 200;                      // distance in which to find closest visible neighbor

  int nearestvisid;                             // id of nearest visible neighbor

  bool socialfactor = 2;                        // switch for social behaviour (0 = avoid, 1 = correct angle, 2 = reconnect)

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

// check distances to all other fish and determine social factor:
  for(int i = 0; i < inds.n;i++){

      repulsion = 0;
      inds.coords[i][0] = inds.coords[i][0] - inds.coords[id][0];          // set fish id x to 0
      inds.coords[i][1] = inds.coords[i][1] - inds.coords[id][1];          // set fish id y to 0
      dist = sqrt(pow(inds.coords[i][0],2) + pow(inds.coords[i][1],2)); // absolute value

      if(i != id){
        angle = acos((dirvector[0]*inds.coords[i][0]+dirvector[1]*inds.coords[i][1])/dist)*180/M_PI;
        if(dist < nearestvisdist){
          nearestvisdist = dist;
          nearestvisid = i;
        }
      }

      else{
        angle = 180; // never see yourself
      }

      if(angle < (fov/2) && dist < 200){
        repulsion = 1 / (0.01*pow(dist,2)+1); // weighted repulsion impact of fish in proximity (200) (f(x) = 1/(0.01x^2+1))
      }

      if(angle < (fov/2) && dist < 50){       // if no individual is detected in radius group attraction takes place
        socialfactor = 1;
      }

      else if(angle < (fov/2) && dist < 30){
        socialfactor = 0;
      }

      if(i != id){
        coordsum[0] = coordsum[0] + inds.coords[i][0] * repulsion;
        coordsum[1] = coordsum[1] + inds.coords[i][1] * repulsion;
      }
    }
 // get new direction according to social factor:
 double leftcoords[2];  // new possible coords after left turn
 double rightcoords[2]; // new possible coords after right turn

if(socialfactor == 0){
    if(c != 0){
      coordsum[0] = coordsum[0]/c;
      coordsum[1] = coordsum[1]/c;
    }
// left turn:
    leftcoords[0] = (state(inds, id) * inds.dspeed + 1) * inds.speed * cos((inds.dir[id] - 90) * M_PI / 180); // new x coord after left turn of 90 degrees
    leftcoords[1] = (state(inds, id) * inds.dspeed + 1) * inds.speed * sin((inds.dir[id] - 90) * M_PI / 180); // new y coord
// turn right:
    rightcoords[0] = (state(inds, id) * inds.dspeed + 1) * inds.speed * cos((inds.dir[id] + 90) * M_PI / 180); // new x coord after right turn of 90 degrees
    rightcoords[1] = (state(inds, id) * inds.dspeed + 1) * inds.speed * sin((inds.dir[id] + 90) * M_PI / 180); // new y coord

    double avoiddistl = sqrt(pow(leftcoords[0] - coordsum[0],2) + pow(leftcoords[1] - coordsum[1],2));
    double avoiddistr = sqrt(pow(rightcoords[0] - coordsum[0],2) + pow(rightcoords[1] - coordsum[1],2));

    if(avoiddistl > avoiddistr){
      newdir = inds.dir[id] - 45;
    }
    else if(avoiddistl == avoiddistr){
      newdir = inds.dir[id] + (rand()%2 - 0.5) * 90; // +- 90 degrees
    }
    else{
      newdir = inds.dir[id] + 45;
    }
  }

if(socialfactor == 2){
    double neighbordistl[2];
    neighbordistl[0] = sqrt(pow(inds.coords[nearestvisid][0],2) + pow(inds.coords[nearestvisid][1],2));
    neighbordistl[1] = neighbordistl[0];

    double neighbordistr[2];
    neighbordistr[0] = sqrt(pow(inds.coords[nearestvisid][0],2) + pow(inds.coords[nearestvisid][1],2));
    neighbordistr[1] = neighbordistr[0];

    double newdirleft = newdir;
    double newdirright = newdir;

    bool turnleft = 0;
    bool turnright = 0;

    while(turnleft == 0 && turnright ==0){
      newdirleft = newdirleft - 10;
      leftcoords[0] = (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirleft * M_PI / 180);
      leftcoords[1] = (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirleft * M_PI / 180);

      newdirright = newdirright + 10;
      rightcoords[0] = (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirright * M_PI / 180);
      rightcoords[1] = (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirright * M_PI / 180);

      neighbordistl[1] = sqrt(pow((leftcoords[0] - inds.coords[nearestvisid][0]),2) + pow((leftcoords[1] - inds.coords[nearestvisid][1]),2));
      neighbordistr[1] = sqrt(pow((rightcoords[0] - inds.coords[nearestvisid][0]),2) + pow((rightcoords[1] - inds.coords[nearestvisid][1]),2));

      if(neighbordistl[1] < neighbordistl[0]){
        neighbordistl[0] = neighbordistl[1];
      }else{
        turnleft = 1;
      }

      if(neighbordistr[1] < neighbordistr[0]){
        neighbordistr[0] = neighbordistr[1];
      }else{
        turnright = 1;
      }
    }
    if(turnleft == 1 && turnright == 1){
      if(rand()%2 == 0){
        newdir = newdirleft;
      }
      else{
        newdir = newdirright;
      }
    }
    else if (turnleft == 1 && turnright == 0){
  		newdir = newdirleft;
  	}
  	else{
  		newdir = newdirright;
  	}
  }
  return newdir;
}
