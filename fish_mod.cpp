#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <vector>
// #include <fstream>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"
#include <opencv2/imgproc.hpp>
#include "opencv2/core/core_c.h"
#include "opencv2/highgui/highgui_c.h"

using namespace std;
using namespace cv;

#define w 1000
#define h 700
#define num 80

const int sp_max = 50;
int sp;
const int dsp_max = 15;
int dsp;
const int nf_max = num;
int nf;

class individuals{

public:

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

int follow(individuals inds,int id, int fov,double r); // fov: field of view as angle
int on_trackbar_nfollow(int,void*);
int lag(bool state, int lag);

void move(individuals& inds);
void initialize(individuals& inds,int n, int dim, double speed, double dspeed, int nfollow);
void drawfish(individuals inds, Mat fish_image);

double on_trackbar_dspeed(int,void*);
double on_trackbar_speed(int,void*);
double correctangle(double dir);
double social(individuals inds,int id,int fov);
double collision(individuals inds, int id, double newdir);
double getangle(individuals inds, int id, double* position);
double averageturn(vector<double> vec);

bool state(individuals inds, int id);

int main(){
  // cout << "Created using: OpenCV version 2.4.13" << endl;
  // cout << "Currently using: " << "OpenCV version : " << CV_VERSION << endl;
  srand(time(0));

  // Creating window for imaging
  char fish_window[] = "Fish";
  Mat fish_image = Mat::zeros(h,w,CV_8UC3);

  individuals fish;
  fish.speed = 1;
  fish.dspeed = 5;
  fish.nfollow = 100;
  fish.dim = 2;

  initialize(fish,num,fish.dim,fish.speed,fish.dspeed,fish.nfollow); // initialize: object fish,num individuals, dimensions,speed,dspeed,n individuals to follow

  for(int z = 0; z < 10000; z++){
      Mat fish_image = Mat::zeros(h,w,CV_8UC3);

      move(fish);
      drawfish(fish,fish_image);
      imshow(fish_window,fish_image);

      char speed_trackbar[10];  // create trackbar in "Fish" window for changing speed
      sprintf(speed_trackbar,"%g",fish.speed);
      createTrackbar("Speed","Fish",&sp,sp_max);
      fish.speed = on_trackbar_speed(fish.speed,0);
      //
      // char dspeed_trackbar[10];  // create trackbar in "Fish" window for changing speed
      // sprintf(dspeed_trackbar,"%g",fish.dspeed);
      // createTrackbar("DSpeed","Fish",&dsp,dsp_max);
      // fish.dspeed = on_trackbar_dspeed(fish.dspeed,0);
      //
      // char nfollow_trackbar[10];  // create trackbar in "Fish" window for changing speed
      // sprintf(nfollow_trackbar,"%i",fish.nfollow);
      // createTrackbar("N-Follow","Fish",&nf,nf_max);
      // fish.nfollow = on_trackbar_nfollow(fish.nfollow,0);
      // cout << "speed: "<< fish.speed << " dspeed: "<< fish.dspeed << " N-Follow: "<< fish.nfollow << "\n";

      waitKey(30);
}
  return(0);
}

// int on_trackbar_nfollow(int,void*)
// {
// return nf;
// }
//
// double on_trackbar_dspeed(int,void*)
// {
// return dsp;
// }

double on_trackbar_speed(int,void*)
{
return sp;
}

void move(individuals& inds){
  for(int id = 0; id < inds.n;id++){

    inds.state[id] = state(inds, id);
    double dir = correctangle(inds.dir[id]);
    double socialdir = dir + social(inds,id,330);

    // if(inds.state[id] == 0 && socialdir == dir){  // feeding, initialize random walk
      // socialdir = socialdir + rand()%91-45;
      // socialdir = socialdir + rand()%21-10;
      // socialdir = correctangle(socialdir);
    // }

    socialdir = socialdir + rand()%21-10;
    inds.lag[id] = lag(inds.state[id],inds.lag[id]);

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
    inds.dir[i] = rand()%361; // all individuals facing x-axis
    // inds.dir[i] = rand() % 20 + 350; // field of vision

    for(int u = 0; u < inds.dim; u++){
      if(u == 0){
        inds.coords[i][u] = w/8 + rand()%201 - 100;
      }
      else{
        inds.coords[i][u] = h/2 + rand()%201 - 100;
        // inds.coords[i][u] = rand()%100;
      }
    }
  }
}

void drawfish(individuals inds, Mat fish_image){
  for(int i = 0; i < inds.n; i++){
    // Scalar fishcol(100,inds.coords[i][1],inds.coords[i][0]);
    // Scalar fishcol(0,0,255);
    Scalar fishcol(inds.colors[i][0],inds.colors[i][1],inds.colors[i][2]);

    Point center(inds.coords[i][0],inds.coords[i][1]);
    Point tail(inds.coords[i][0] - 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] - 10 * sin(inds.dir[i] * M_PI / 180));
    Point head(inds.coords[i][0] + 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] + 10 * sin(inds.dir[i] * M_PI / 180));

    // circle(fish_image,center,30,fishcol,1,CV_AA,0);            // visualize center
    // circle(fish_image,center,100,fishcol,1,CV_AA,0);           // visualize comfort zone
    // circle(fish_image,center,200,fishcol,1,CV_AA,0);           // visualize visual distance
    // char angle[10];
    // sprintf(angle,"%g",inds.dir[i]);
    // putText(fish_image,angle,center,1,1,fishcol,1,8);
    arrowedLine(fish_image,tail,head,fishcol,1,CV_AA,0,0.2);   // draw fish
  }
}



int follow(individuals inds,int id,int fov,double r){
  // count number of inds in front (x-axis) of id individual
  int c = 0;

  double dirvector[inds.dim]; // angle of fish[id] is facing in respect to original x-axis
  double dist; // absolute distance between fish [id] and fish [i]
  double angle; // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

  for (int i = 0; i < inds.n;i++){

    double coords[2];
    coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id x to 0
    coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id y to 0

    dist = sqrt(pow(coords[0],2) + pow(coords[1],2)); // absolute value

    if(i != id){
      angle = acos((dirvector[0]*coords[0]+dirvector[1]*coords[1])/dist)*180/M_PI;
    }

    else{
      angle = 180; // never see yourself
    }

    if (angle < (fov/2) && dist < r){ // visual field with angle and radius
      c++;
      }
    }
  return c;
}

bool state(individuals inds, int id){
  bool state = 0;
  if(follow(inds,id,180,200) >= inds.nfollow && inds.lag[id] == 0){
    state = 1;
  }
  return state;
}

int lag(bool state, int lag){
  if(state == 0){
    if(lag == 0){  // individual lag phase
      lag = rand()%10;
    }
    else{
      lag = lag-1;
    }
  }
  return lag;
}

double collision(individuals inds, int id, double newdir){
	double bound = 20; // boundary in which to detect the border
	double leftcoords[2]; // new possible coords after left turn
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

double social(individuals inds, int id, int fov){
	double dirvector[2]; // direction vector (fish id) of lenth 1
  double repulsion; // repulsion value of other fish..
  double newdir = inds.dir[id]; // angle after avoidance behavior
  double avoidir; // direciton in whicht to avoid nearestvisid

  int c_avoid = 0;  // counter to initialize array
  int c_attract = 0;
  int c_align = 0;

  int c = follow(inds, id, fov, 200); // counter for fish in visual radius
  int nearestvisid; // id of nearest seeable neighbor
  int socialfactor = 2; // switch for behavior: 0 = avoidance, 1 = alignment, 2 = attraction

  int n_avoid = follow(inds, id, fov, 30); // number of individals to avoid. [id,dist,angle,angle to neighbor]
  int n_align = follow(inds, id, fov, 100) - n_avoid; // number of individals to align to [id,dist,angle between directions]
  int n_attract = follow(inds, id, fov, 200) - (n_avoid + n_align); // number of individals to be attracted to [id,dist,angle to neighbor, turning angle]

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

  vector<double> neighbors_avoid_id(n_avoid);
  vector<double> neighbors_avoid_dist(n_avoid);
  vector<double> neighbors_avoid_angle(n_avoid);
  vector<double> neighbors_avoid_turn(n_avoid);

  vector<double> neighbors_align_id(n_align);
  vector<double> neighbors_align_dist(n_align);
  vector<double> neighbors_align_angle(n_align);
  vector<double> neighbors_align_turn(n_align);

  vector<double> neighbors_attract_id(n_attract);
  vector<double> neighbors_attract_dist(n_attract);
  vector<double> neighbors_attract_angle(n_attract);
  vector<double> neighbors_attract_turn(n_attract);

	// check distances for all other fish and determine social factor (0 = avoid, 1 = correct angle, 2 = reconnect with others)
	for (int i = 0; i < inds.n; i++){
    if(i != id){
      double coords[2];
      	coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
      	coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0
        double dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id
        double angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / M_PI;

      		if (angle < (fov / 2) && dist < 200){
            if(dist < 30){
              neighbors_avoid_id.at(c_avoid) = i; // first column: id of fish
              neighbors_avoid_dist.at(c_avoid) = dist;
              neighbors_avoid_angle.at(c_avoid) = getangle(inds,id,inds.coords[i]);

              if(neighbors_avoid_angle.at(c_avoid) > 0){
                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) - 180;
              }
              else if(neighbors_avoid_angle.at(c_avoid) < 0){
                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + 180;
              }
              else{
                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + (rand()%2 - 0.5) * 360;
              }
              neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_turn.at(c_avoid) * (1 / ( 0.01 * pow((dist),2) + 4 )); // change in angle with algebraic sign (+,-)
              c_avoid++;
            }
            else if(dist < 100){
              neighbors_align_id.at(c_align) = i; // id of neighbor fish
              neighbors_align_dist.at(c_align) = dist; // distance between neighbor and focal fish

              double dir[2];
              dir[0] = cos(inds.dir[i] * M_PI / 180) + inds.coords[id][0]; // unit vector of neighbor i, starting at focal fish
              dir[1] = sin(inds.dir[i] * M_PI / 180) + inds.coords[id][1];

              neighbors_align_angle.at(c_align) = getangle(inds,id,dir); // angle to neighbor
              neighbors_align_turn.at(c_align) = neighbors_align_angle.at(c_align) * (1 / ( 0.01 * pow((dist-30),2) + 2 )); // turning angle, weighted according to distance of fish i
              c_align++;
            }
            else{
              neighbors_attract_id.at(c_attract) = i; // id of neighbor fish
              neighbors_attract_dist.at(c_attract) = dist; // distance between neighbor and focal fish
              neighbors_attract_angle.at(c_attract) = getangle(inds,id,inds.coords[i]);
              neighbors_attract_turn.at(c_attract) = neighbors_attract_angle.at(c_attract) * (1 / ( 0.01 * pow(( dist-200),2) + 2 )); // turning angle, weighted according to distance fish i
              c_attract++;
              }
            }
      		}
        }
        double turn = 0;

        if(n_avoid > 0){
          turn = averageturn(neighbors_avoid_turn);
        }
        else{
          if(n_align > 0 && n_attract > 0){
            vector<double> comb_al_at(2);
            comb_al_at.at(0) = averageturn(neighbors_align_turn);
            comb_al_at.at(1) = averageturn(neighbors_attract_turn);
            turn = averageturn(comb_al_at);
          }
          else if(n_align > 0){
            turn = averageturn(neighbors_align_turn);
          }
          else if(n_attract > 0){
            turn = averageturn(neighbors_attract_turn);
          }
          else{
            turn = rand()%91 - 45; // random walk if no social interactions present
          }
        }
      return turn;

// 	if (socialfactor == 0){
//     // newdir = newdir - getangle(inds,id,coordsum) * avoidsum;
//     if(getangle(inds,id,inds.coords[nearestvisid]) < 0){
//       avoidir = getangle(inds,id,inds.coords[nearestvisid]) + 180;
//     }
//     else if(getangle(inds,id,inds.coords[nearestvisid]) == 0){
//       if(rand()%2 == 0){
//         avoidir = getangle(inds,id,inds.coords[nearestvisid]) + 180;
//       }
//       else{
//         avoidir = getangle(inds,id,inds.coords[nearestvisid]) - 180;
//       }
//     }
//     else{
//       avoidir = getangle(inds,id,inds.coords[nearestvisid]) - 180;
//     }
//     newdir = newdir + (1 / ( 0.01 * pow((nearestvisdist),2) + 4 ))* avoidir;
// 	}
//   else if(socialfactor == 1){
//     double dir[2];
//     dir[0] = cos(inds.dir[nearestvisid] * M_PI / 180); // unit vector of nearestvisid
//     dir[1] = sin(inds.dir[nearestvisid] * M_PI / 180); // unit vector of nearestvisid
//
//     dir[0] = dir[0] + inds.coords[id][0];
//     dir[1] = dir[1] + inds.coords[id][1];
//
//     newdir = newdir + (1 / ( 0.01 * pow((nearestvisdist-30),2) + 2 )) * getangle(inds,id,dir);
//   }
//
// 	else if (socialfactor == 2 && c > 0){
//       // newdir = newdir + ((nearestvisdist)/400) * getangle(inds, id, inds.coords[nearestvisid]); // angle corrected to swim towards nearest visible individual with visible field of 200 and multiplied by 0.5
//       newdir = newdir + (1 / ( 0.01 * pow((nearestvisdist-200),2) + 2 )) * getangle(inds, id, inds.coords[nearestvisid]);
// 	}
// 	return newdir;
}

double averageturn(vector<double> vec){
  double x_sum = 0;
  double y_sum = 0;
  for(int i = 0; i < vec.size();i++){
      x_sum = x_sum + cos(vec.at(i) * M_PI / 180); // unit vector of neighbor i, starting at focal fish
      y_sum = y_sum + sin(vec.at(i) * M_PI / 180);
  }
  double x = x_sum/vec.size();
  double y = y_sum/vec.size();
  double angle = acos(x / sqrt(pow(x, 2) + pow(y, 2))) * 180 / M_PI; // angle to average turnpoint from x-axis == 0 degrees (0 - 180°)

  if (y < 0){
		angle = 360 - angle;     // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
	}
	else if (y == 0){   // correct 0 angle and 0,0 position
		if (x >= 0){ // reference point in front of focal point.don't change angle
			angle = 0;
		}
		else{ // reference point behind focal point
			angle = 180;
		}
	}
  return angle;
}

double correctangle(double dir){
	if (dir < 0){
		dir = dir + 360;
	}
	else if (dir >= 360){
		dir = dir - 360;
	}
	return dir;
}

double getangle(individuals inds, int id, double* position){
	double point[2];
  point[0] = position[0] - inds.coords[id][0]; // move point in respect to setting fish id to 0,0
	point[1] = position[1] - inds.coords[id][1];
	double angle = acos(point[0] / sqrt(pow(point[0], 2) + pow(point[1], 2))) * 180 / M_PI; // arccos of dot product between id direction and point position

	if (point[1] < 0){
		angle = 360 - angle;     // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
	}
	else if (point[1] == 0){   // correct 0 angle and 0,0 position
		if (point[0] >= 0){
			angle = 0;
		}
		else{
			angle = 180;
		}
	}

	angle = angle - inds.dir[id];          // angle difference between id direction and point position, can be negative
	double newangle = correctangle(angle); // set angle difference between 0 < x =< 360

	if (newangle > 180){
		newangle = newangle - 360; // left is negative, right is positive; angles between -180 < angle =< 180
	}
	return newangle;
}
