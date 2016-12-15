#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <stdlib.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"
#include <opencv2/imgproc.hpp>
#include "opencv2/core/core_c.h"
#include "opencv2/highgui/highgui_c.h"

using namespace std;
using namespace cv;

#define latitude 14
#define longitude 7
#define resolution 4 // must be 100 > resolution > 1
#define res_factor 25 // should be: res_factor * resolution = 100

#define w (latitude * resolution * res_factor)
#define h (longitude * resolution * res_factor)
#define rows (longitude*resolution)
#define cols (latitude*resolution)

#define num 45

const int sp_max = 50;
int sp;
const int dsp_max = 15;
int dsp;
const int nf_max = num;
int nf;

random_device rd;
uniform_int_distribution<int>distrib_choice(0,1);
uniform_int_distribution<int>distrib_error(-10,10);
uniform_int_distribution<int>distrib_randomwalk(-45,45);
uniform_int_distribution<int>distrib_quality(0,10);
uniform_int_distribution<int>distrib_seed(0,10);
uniform_int_distribution<int>distrib_sample(20,100);


class individuals{

  public:

    int n;
    int dim;
    int lag[num];
    int colors[num][3];     // individual color in RGB
    int nfollow;
    int sample_rate [num];
    int speed_factor[num];

    double speed[num];      // movement speed
    double coords [num][2]; // height n * width dim
    double dir [num];       // array of agles between 0 - 360 degrees for each individual

    bool sample[num];

  };

int in_zone(individuals inds,int id, int fov,double r); // fov: field of view as angle
// int on_trackbar_nfollow(int,void*);
int get_quality(individuals inds, int id, int quality[][cols]);
int sum_array(int quality[][cols]);

void move(individuals& inds, int quality[][cols]);
void initialize(individuals& inds,int n, int dim, double speed);
void drawfish(individuals inds, Mat fish_image);
void correct_coords(double* coords);
void create_environment(Mat& environment,int quality[][cols]);
void feed(individuals& inds,int id, int quality[][cols],Mat& environment);
void sample(individuals& inds, int quality[][cols], Mat& environment);

double correctangle(double dir);
double social(individuals inds,int id,int fov);
double collision(individuals inds, int id, double newdir);
double getangle(double* focal_coords, double focal_dir, double* position);
double averageturn(vector<double> vec, bool comb_weight);
double get_speed(individuals inds,int id);

int main(){

  srand(time(0));

  ofstream sum_quality;
  sum_quality.open("sum_quality.csv");

  while(sum_quality){
    for(int n = 0;n < 10;n++){
  // cout << "Created using: OpenCV version 2.4.13" << endl;
  // cout << "Currently using: " << "OpenCV version : " << CV_VERSION << endl;

  // Creating window for imaging
  char fish_window[] = "Fish";
  Mat fish_image = Mat::zeros(h,w,CV_8UC3);

  // Create "food" environment
  char env_window[] = "Environment";
  Mat environment = Mat::zeros(h,w,CV_8UC3);
  Mat fishinenvironment = Mat::zeros(h,w,CV_8UC3);
  // Establish mixed pattern of environment with quality[h]x[w] squares
  int quality[rows][cols];
  create_environment(environment,quality);

  // imshow(env_window,environment);

  individuals fish;

  initialize(fish,num,2,3); // initialize: object fish,num individuals, dimensions,speed

  // sum_quality << sum_array(quality) << ",";
  sum_quality << endl;
  for(int z = 0; z < 1000; z++){
      // plot fish on black background
      // fish_image = Mat::zeros(h,w,CV_8UC3);

      sample(fish,quality,environment);
      move(fish,quality);
      // drawfish(fish,fish_image);
      // addWeighted(fish_image,1,environment,0.1,10,fishinenvironment);
      // imshow(fish_window,fishinenvironment);

      sum_quality << sum_array(quality) << ",";

      // waitKey(30);
      }
    }
    sum_quality.close();
    }
    return(0);
  }

void move(individuals& inds,int quality[][cols]){
  for(int id = 0; id < inds.n;id++){

    if(inds.lag[id] == 0){
      double dir = correctangle(inds.dir[id]);
      double turn = social(inds,id,330) + distrib_error(rd); // maximum error of movement is +- 5°

      dir = dir + turn;

      inds.speed_factor[id] = get_speed(inds,id);

      inds.coords[id][0] =  inds.coords[id][0] + inds.speed[id] * cos(dir * M_PI / 180) * inds.speed_factor[id];
      inds.coords[id][1] = inds.coords[id][1] + inds.speed[id] * sin(dir * M_PI/ 180) * inds.speed_factor[id];
      correct_coords(inds.coords[id]);
      inds.dir[id] = dir;
    }
    else{
      inds.dir[id] = inds.dir[id] + distrib_error(rd);
    }
  }
}

void initialize(individuals& inds, int n ,int dim, double speed){
  inds.n = n;
  inds.dim = dim;

  for(int i = 0; i < inds.n; i++){
    inds.sample_rate[i] = 100;
    inds.speed_factor[i] = 1;
    inds.speed[i] = speed;
    inds.colors[i][0] = 80 + rand()%155;
    inds.colors[i][1] = 80 + rand()%155;
    inds.colors[i][2] = 80 + rand()%155;
    inds.dir[i] = distrib_randomwalk(rd);
    inds.lag[i] = 0; // equals handling time
    // inds.dir[i] = rand() % 20 + 350; // field of vision

    for(int u = 0; u < inds.dim; u++){
      if(u == 0){
        inds.coords[i][u] = w/2 + rand()%201 - 100;
      }
      else{
        inds.coords[i][u] = h/2 + rand()%201 - 100;
      }
    }
  }
}

void drawfish(individuals inds, Mat fish_image){
  for(int i = 0; i < inds.n; i++){
    // Scalar fishcol(5,inds.coords[i][1],inds.coords[i][0]);
    // Scalar fishcol(0,0,255);
    Scalar fishcol(inds.colors[i][0],inds.colors[i][1],inds.colors[i][2]);

    Point center(inds.coords[i][0],inds.coords[i][1]);
    Point tail(inds.coords[i][0] - 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] - 10 * sin(inds.dir[i] * M_PI / 180));
    Point head(inds.coords[i][0] + 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] + 10 * sin(inds.dir[i] * M_PI / 180));

    circle(fish_image,center,2,fishcol,1,CV_AA,0);            // visualize center
    // circle(fish_image,center,100,fishcol,1,CV_AA,0);       // visualize comfort zone
    // circle(fish_image,center,200,fishcol,1,CV_AA,0);       // visualize visual distance
    // char angle[10];
    // sprintf(angle,"%g",inds.dir[i]);
    // putText(fish_image,angle,center,1,1,fishcol,1,8);
    arrowedLine(fish_image,tail,head,fishcol,1,CV_AA,0,0.2);   // draw fish
  }
}

int in_zone(individuals inds,int id,int fov,double r){
  // count number of inds in front (x-axis) of id individual
  int c = 0;

  double dirvector[inds.dim]; // angle of fish[id] is facing in respect to original x-axis
  double dist; // absolute distance between fish [id] and fish [i]
  double angle; // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis

  dirvector[0] = cos(inds.dir[id] * M_PI / 180);
  dirvector[1] = sin(inds.dir[id] * M_PI / 180);

  for (int i = 0; i < inds.n;i++){

    if(i != id){

    double coords[2];
    coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id x to 0
    coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id y to 0

    if(coords[0] > w/2){
      coords[0] = coords[0] - w;
    }
    else if(coords[0] < -w/2){
      coords[0] = coords[0] + w;
    }

    if(coords[1] > h/2){
      coords[1] = coords[1] - h;
    }
    else if(coords[1] < -h/2){
      coords[1] = coords[1] + h;
    }

    dist = sqrt(pow(coords[0],2) + pow(coords[1],2)); // absolute value

    angle = acos((dirvector[0]*coords[0]+dirvector[1]*coords[1])/dist)*180/M_PI; // dot prodcut (direction vector of id, position vector of i)

    if (angle < (fov/2) && dist < r){ // visual field with angle and radius
      c++;
        }
      }
    }
  return c;
}

double social(individuals inds, int id, int fov){
	double dirvector[2]; // direction vector (fish id) of lenth 1
  double repulsion; // repulsion value of other fish..
  double newdir = inds.dir[id]; // angle after avoidance behavior
  double avoidir; // direciton in whicht to avoid nearestvisid

  int c_avoid = 0;  // counter to initialize array
  int c_attract = 0;
  int c_align = 0;

  int c = in_zone(inds, id, fov, 200); // counter for fish in visual radius
  int socialfactor = 2; // switch for behavior: 0 = avoidance, 1 = alignment, 2 = attraction

  int n_avoid = in_zone(inds, id, fov, 15); // number of individals to avoid. [id,dist,angle,angle to neighbor]
  int n_align = in_zone(inds, id, fov, 100) - n_avoid; // number of individals to align to [id,dist,angle between directions]
  int n_attract = in_zone(inds, id, fov, 200) - (n_avoid + n_align); // number of individals to be attracted to [id,dist,angle to neighbor, turning angle]

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

        if(coords[0] > w/2){
          coords[0] = coords[0] - w;
        }
        else if(coords[0] < -w/2){
          coords[0] = coords[0] + w;
        }

        if(coords[1] > h/2){
          coords[1] = coords[1] - h;
        }
        else if(coords[1] < -h/2){
          coords[1] = coords[1] + h;
        }

        double dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id
        double angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / M_PI;

      		if (angle < (fov / 2) && dist < 200){
            if(dist < 15){
              neighbors_avoid_id.at(c_avoid) = i; // first column: id of fish
              neighbors_avoid_dist.at(c_avoid) = dist; // distance between focal fish and neighbor
              neighbors_avoid_angle.at(c_avoid) = getangle(inds.coords[id],inds.dir[id],inds.coords[i]); // angle to neighbor

              if(neighbors_avoid_angle.at(c_avoid) > 0){
                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) - 180;
              }
              else if(neighbors_avoid_angle.at(c_avoid) < 0){
                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + 180;
              }
              else{
                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + (distrib_choice(rd)-0.5) * 360;
              }
              neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_turn.at(c_avoid) * (1 / ( 0.1 * pow((dist),2) + 2 )); // change in angle with algebraic sign (+,-)
              c_avoid++;
            }
            else if(dist < 100){
              neighbors_align_id.at(c_align) = i; // id of neighbor fish
              neighbors_align_dist.at(c_align) = dist; // distance between neighbor and focal fish

              double dir[2];
              dir[0] = cos(inds.dir[i] * M_PI / 180) + inds.coords[id][0]; // unit vector of neighbor i, starting at focal fish
              dir[1] = sin(inds.dir[i] * M_PI / 180) + inds.coords[id][1];

              neighbors_align_angle.at(c_align) = getangle(inds.coords[id],inds.dir[id],dir); // angle to neighbor
              neighbors_align_turn.at(c_align) = neighbors_align_angle.at(c_align) * (1 / ( 0.004 * pow((dist-15),2) + 2 )); // turning angle, weighted according to distance of fish i
              c_align++;
            }
            else{
              neighbors_attract_id.at(c_attract) = i; // id of neighbor fish
              neighbors_attract_dist.at(c_attract) = dist; // distance between neighbor and focal fish
              neighbors_attract_angle.at(c_attract) = getangle(inds.coords[id],inds.dir[id],inds.coords[i]);
              neighbors_attract_turn.at(c_attract) = neighbors_attract_angle.at(c_attract) * (1 / ( 0.004 * pow(( dist-200),2) + 2 )); // turning angle, weighted according to distance fish i
              c_attract++;
              }
            }
      		}
        }
        double turn = 0;

        if(n_avoid > 0){
          turn = averageturn(neighbors_avoid_turn,0);
        }
        else{
          if(n_align > 0 && n_attract > 0){
            vector<double> comb_al_at(2);
            comb_al_at.at(0) = averageturn(neighbors_align_turn,0);
            comb_al_at.at(1) = averageturn(neighbors_attract_turn,0);
            turn = averageturn(comb_al_at,1);
          }
          else if(n_align > 0){
            turn = averageturn(neighbors_align_turn,0);
          }
          else if(n_attract > 0){
            turn = averageturn(neighbors_attract_turn,0);
          }
          else{
            turn = distrib_randomwalk(rd); // random walk if no social interactions present
          }
        }
      return turn;
}

double averageturn(vector<double> vec, bool comb_weight){
  double x_sum = 0;
  double y_sum = 0;
  for(int i = 0; i < vec.size();i++){
    if(comb_weight){
      // weighting difference between attraction and alignment
      x_sum = x_sum + cos(vec.at(i) * M_PI / 180) * (1.2 - i); // unit vector of neighbor i, starting at focal fish
      y_sum = y_sum + sin(vec.at(i) * M_PI / 180) * (1.2 - i);
    }
    else{
      x_sum = x_sum + cos(vec.at(i) * M_PI / 180); // unit vector of neighbor i, starting at focal fish
      y_sum = y_sum + sin(vec.at(i) * M_PI / 180);
    }
  }
  double x = x_sum/vec.size();
  double y = y_sum/vec.size();
  double angle = acos(x / sqrt(pow(x, 2) + pow(y, 2))) * 180 / M_PI; // angle to average turnpoint from x-axis == 0 degrees (0 - 180°)

  if (y < 0){
		angle = 0 - angle;     // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
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

double getangle(double* focal_coords, double focal_dir, double* position){
  // stream cout to file
	double point[2];
  point[0] = position[0] - focal_coords[0]; // move point in respect to setting fish id to 0,0
	point[1] = position[1] - focal_coords[1];

  if(point[0] > w/2){
    point[0] = point[0] - w;
  }
  else if(point[0] < -w/2){
    point[0] = point[0] + w;
  }

  if(point[1] > h/2){
    point[1] = point[1] - h;
  }
  else if(point[1] < -h/2){
    point[1] = point[1] + h;
  }

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

	angle = angle - focal_dir; // angle difference between id direction and point position, can be negative
	double newangle = correctangle(angle); // set angle difference between 0 < x =< 360

	if (newangle > 180){
		newangle = newangle - 360; // left is negative, right is positive; angles between -180 < angle =< 180
	}
  else if(newangle == 180){
    newangle = (distrib_choice(rd) - 0.5) * 360; // random choice between left or right turn to 180°
  }
	return newangle;
}

void correct_coords(double* coords){

  if(coords[0] > w){
    coords[0] = coords[0] - w;
  }
  else if(coords[0] < 0){
    coords[0] = coords[0] + w;
  }

  if(coords[1] > h){
    coords[1] = coords[1] - h;
  }
  else if(coords[1] < 0){
    coords[1] = coords[1] + h;
  }
}

void create_environment(Mat& environment,int quality[][cols]){
  int seed = distrib_seed(rd);
  for(int u=0;u<rows;u++){
     for(int i =0;i<cols;i++){
       quality[u][i] = 0;
     }
   }

  if(seed > 0){
  int seed_coords[seed][2];
  uniform_int_distribution<int>distrib_rows(0,rows);
  uniform_int_distribution<int>distrib_cols(0,cols);

  for(int i = 0; i < seed;i++){
    seed_coords[i][0] = distrib_rows(rd);
    seed_coords[i][1] = distrib_cols(rd);
    quality[seed_coords[i][0]][seed_coords[i][1]] = 10;
    }
  }
  for(int p = 0;p < 20; p++){
      for(int u=0;u<rows;u++){
        for(int i =0;i<cols;i++){
          if(u < rows-1 && i < cols-1){
            if(quality[u+1][i] > p || quality[u+1][i+1] > p || quality[u][i+1]){
              quality[u][i] = distrib_quality(rd);
          }
        }
      }
    }

  for(int u = rows - 1;u > 0;u--){
      for(int i = cols - 1;i > 0;i--){
          if(quality[u-1][i] > p || quality[u-1][i-1] > p || quality[u][i-1]){
          quality[u][i] = distrib_quality(rd);
        }
      }
    }
  }

  // Color squares according to quality
  for(int x = 0;x < w;x++){
    for(int y = 0; y < h;y++){
      int q_color = quality[y/(h/rows)][x/(w/cols)];
          environment.at<Vec3b>(y,x)= Vec3b(q_color*(255/10),q_color*(255/10),q_color*(255/10));
    }
  }
}

int get_quality(individuals inds, int id, int quality[][cols]){
  int x = inds.coords[id][0];
  int y = inds.coords[id][1];
  int q = quality[y/(h/rows)][x/(w/cols)];
  return q;
}

void sample(individuals& inds, int quality[][cols], Mat& environment){
  for(int id = 0; id < inds.n;id++){
    int n_avoid = in_zone(inds, id, 330, 15); // number of individals to avoid. [id,dist,angle,angle to neighbor]
    int n_align = in_zone(inds, id, 330, 100) - n_avoid; // number of individals to align to [id,dist,angle between directions]
    int n_attract = in_zone(inds, id, 330, 200) - (n_avoid + n_align); // number of individals to be attracted to [id,dist,angle to neighbor, turning angle]
    bool alone = (n_avoid == 0 && n_align == 0 && n_attract > 0);
    // avoid being left behind.lag decreases with increasing distance to group
    if(n_avoid == 0 && n_align == 0 && inds.lag[id] != 0 && n_attract > 0){
      inds.lag[id] = inds.lag[id] - 1;
    }

    if(inds.lag[id] == 0 && inds.sample_rate[id] > distrib_sample(rd) && alone == 0){
        inds.lag[id] = 5;
        feed(inds,id,quality,environment);
      }
    else if(inds.lag[id] == 0){
      inds.sample_rate[id] = inds.sample_rate[id] + 1;
    }
    else{
      inds.lag[id] = inds.lag[id] - 1;
      inds.sample_rate[id] = inds.sample_rate[id] + 1;
    }
  }
}

void feed(individuals& inds,int id, int quality[][cols],Mat& environment){
    int x = inds.coords[id][0];
    int y = inds.coords[id][1];
    int q =  quality[y/(h/rows)][x/(w/cols)];
    inds.sample_rate[id] = 10 * q;

    if(q>0){
      quality[y /(h/rows)][x/(w/cols)] = q - 1;
      q = q - 1;
    }

    // determining box in quality matrix
    int y_min = y/(h/rows);
    y_min = y_min * (h/rows);
    int y_max = y_min + (h/rows);

    int x_min = x/(w/cols);
    x_min = x_min * (w/cols);
    int x_max = x_min + (w/cols);

    // draw environment boxes as filled rectangles (255/max. quality)
    // rectangle(environment,Point(x_min,y_min),Point(x_max-1,y_max-1),Scalar(q*(255/10),q*(255/10),q*(255/10)),-1,8,0);
}

double get_speed(individuals inds,int id){
  int count = in_zone(inds,id,180,100);
  double speed_factor = 1 + count/(num/2);
  return speed_factor;
}

int sum_array(int quality[][cols]){
  int sum = 0;
  for(int i = 0; i < cols;i++){
    for(int u = 0;u < rows;u++){
      sum = sum + quality[u][i];
    }
  }
  return sum;
}
