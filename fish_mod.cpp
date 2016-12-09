#include <iostream>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdio.h>
#include <string>
#include <map>
#include <random>
#include <fstream>

#define w 1200
#define h 600
#define PI 3.14159265
#define num 60

#define rows 4
#define cols 8

const int sp_max = 30;
int sp = 2;
const int dsp_max = 20;
int dsp = 6;
const int nf_max = num;
int nf;

using namespace std;
using namespace cv;

random_device rd;
uniform_int_distribution<int> distrib_choice(0, 1);
uniform_int_distribution<int> distrib_error(-10, 10);
uniform_int_distribution<int> distrib_rwalk(-45, 45);

class Individuals {
	public:
		double coords[num][2]; // height n and width dim
		double dir[num]; // fish direction between 0 and 360°
		int lag[num]; // with length n
		bool state[num]; // 0 for random, 1 for directed
		int color[num][3];
		int n;
		int dim;
		double speed; // movement speed
		double dspeed; //dspeed + 1 = swimming speed / feeding speed
		int nfollow;
};

void move(Individuals& inds, int quality[][cols]);
void initialize(Individuals& inds, int n, int dim, double speed, double dspeed, int nfollow);
void drawfish(Individuals inds, Mat image);
int in_zone(Individuals inds, int id, int fov, double r); // fov: field of view
bool state(Individuals inds, int id);
double collision(Individuals inds, int id, double newdir);
double correct_angle(double dir);
double get_angle(double* focal_coords, double focal_dir, double* pos);
double social(Individuals inds, int id, int fov);
double distance(double* coordsa, double* coordsb);
int lag(bool state, int lag);

double on_trackbar_dspeed(int, void*);
double on_trackbar_speed(int, void*);
int on_trackbar_nfollow(int, void*);

double average_turn(vector<double> vec);
void correct_coords(double* coords);

void create_environment(Mat& environment, int quality[][cols]);
int get_quality(Individuals inds, int id, int quality[][cols]);

int main() {
	srand(time(0));

	char fish_window[] = "Fish";
	Mat fish_image = Mat::zeros(h, w, CV_8UC3);
	// Mat fish_mask = Mat::zeros(h, w, CV_8UC3);

	char environment_window[] = "Environment";
	Mat environment = Mat::zeros(h, w, CV_8UC3);

	int quality[rows][cols];
	create_environment(environment, quality);

	Mat fish_environment;

	Individuals fish;
	fish.speed = 1;
	fish.dspeed = 1;
	fish.nfollow = num;

	// initializes object fish of class Individuals with num individuals, dimensions, speed, dspeed, nfollow
	initialize(fish, num, 2, fish.speed, fish.dspeed, fish.nfollow);

	for (int z = 0; z < 2000; z++) {
		fish_image = Mat::zeros(h, w, CV_8UC3);

		drawfish(fish, fish_image);

		// cvtColor(fish_image, fish_mask, CV_BGR2GRAY, 0);
		// threshold(fish_image, fish_mask, 60, 255, THRESH_BINARY);
		//
		// subtract(environment, fish_mask, fish_environment);

		addWeighted(fish_image, 1, environment, 0.1, 0.0, fish_environment);

		imshow(fish_window, fish_environment);
		waitKey( 30 );

		move(fish, quality);

		char speed_trackbar[10];  // create trackbar in "Fish" window for changing speed
	    sprintf(speed_trackbar,"%g",fish.speed);
	    createTrackbar("Speed","Fish",&sp,sp_max);
	    fish.speed = on_trackbar_speed(fish.speed,0);

	}
	return(0);
}

void move(Individuals& inds, int quality[][cols]) {
	for (int id = 0; id < inds.n; id++) {
		double dir = correct_angle(inds.dir[id]);

		inds.state[id] = 1 - get_quality(inds, id, quality); //state(inds, id);
		double small_error = distrib_error(rd);

		// f_small_error << small_error << ",\n";

		double turn = social(inds, id, 330) + small_error;

		dir = dir + turn;
		inds.lag[id] = lag(inds.state[id], inds.lag[id]);

		inds.coords[id][0] = inds.coords[id][0] + (inds.state[id] * inds.dspeed + 1) * inds.speed * cos(dir * PI / 180);
		inds.coords[id][1] = inds.coords[id][1] + (inds.state[id] * inds.dspeed + 1) * inds.speed * sin(dir * PI / 180);
		correct_coords(inds.coords[id]);

		inds.dir[id] = dir;
	}
}

void initialize(Individuals& inds, int n, int dim, double speed, double dspeed, int nfollow) {
	inds.n = n;
	inds.dim = dim;
	inds.speed = speed;
	inds.dspeed = dspeed;
	inds.nfollow = nfollow;
	for (int i = 0; i < inds.n; i++) {
		inds.lag[i] = 0;
		inds.color[i][0] = 120 + rand() % 120;
		inds.color[i][1] = 120 + rand() % 120;
		inds.color[i][2] = 120 + rand() % 120;
		inds.dir[i] = rand() % 360; // direction
		for (int u = 0; u < inds.dim; u++) {
			if (u == 0) {
				inds.coords[i][u] = w / 2 + rand() % 201 - 100; // x coord
			}
			else {
				inds.coords[i][u] = h / 2 + rand() % 201 - 100; // y coord
			}
		}
	}
}

void drawfish(Individuals inds, Mat image) {
	for (int i = 0; i < inds.n; i++) {
		Scalar fishcolor(inds.color[i][0], inds.color[i][1], inds.color[i][2]);
		Point center(inds.coords[i][0], inds.coords[i][1]);
		Point tail(inds.coords[i][0] - 10 * cos(inds.dir[i] * PI / 180), inds.coords[i][1] - 10 * sin(inds.dir[i] * PI / 180));
		Point head(inds.coords[i][0] + 10 * cos(inds.dir[i] * PI / 180), inds.coords[i][1] + 10 * sin(inds.dir[i] * PI / 180));
		circle(image, center, 2, fishcolor, 1, CV_AA, 0);
		// circle(image, center, 30, fishcolor, 1, CV_AA, 0);
		// circle(image, center, 200, fishcolor, 1, CV_AA, 0);
		// circle(image, center, 100, fishcolor, 1, CV_AA, 0);
		arrowedLine(image, tail, head, fishcolor, 1, CV_AA, 0, 0.2);
		// char angle[10];
		// sprintf(angle, "%g", inds.dir[i]);
		// putText(image, angle, head, 1, 1, fishcolor, 1, CV_AA, false );
	}
}

int in_zone(Individuals inds, int id, int fov, double r) {
	int c =  0; // count num of inds in front of id
	double dirvector[2]; // direction vector (fish id) of lenth 1
	dirvector[0] = cos(inds.dir[id] * PI / 180);
	dirvector[1] = sin(inds.dir[id] * PI / 180);
	double dist; // absolut distance between fish i and fish id // does not need to be an array
	double angle; // angle between fish id direction vector and fish i pos vector
	for (int i = 0; i < inds.n; i++) {
		if (i != id) {
			double coords[2];
			coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0 and correct for continous environment
			if (coords[0] > w / 2) {
				coords[0] = coords[0] - w;
			}
			else if (coords[0] < - w / 2) {
				coords[0] = coords[0] + w;
			}
			coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
			if (coords[1] > h / 2) {
				coords[1] = coords[1] - h;
			}
			else if (coords[1] < - h / 2) {
				coords[1] = coords[1] + h;
			}

			dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id:

			// dot product(direction vector of id, position vector of i)
			angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;

			// f_all_dotproduct << angle << ",\n";

			if (angle < (fov / 2) && dist < r) {
				c++;
			}
		}
	}
	return c;
}

double social(Individuals inds, int id, int fov) {
	double dirvector[2]; // direction vector (focal fish) of lenth 1
	dirvector[0] = cos(inds.dir[id] * PI / 180);
	dirvector[1] = sin(inds.dir[id] * PI / 180);

	int n_avoid = in_zone(inds, id, fov, 15); // number of neighbors in respective zones
	int n_align = in_zone(inds, id, fov, 100) - n_avoid;
	int n_attract = in_zone(inds, id, fov, 200) - (n_avoid + n_align);

	vector<double> neighbors_avoid_id(n_avoid); // identies of all neighbors in zone of repulsion: id, dist, angle to neighbor, newdir of focal fish
	vector<double> neighbors_avoid_dist(n_avoid);
	vector<double> neighbors_avoid_angle(n_avoid);
	vector<double> neighbors_avoid_turn(n_avoid);

	vector<double> neighbors_align_id(n_align); // identies of all neighbors in zone of alignment: id, dist, angle between directions, newdir of focal fish
	vector<double> neighbors_align_dist(n_align);
	vector<double> neighbors_align_angle(n_align);
	vector<double> neighbors_align_turn(n_align);

	vector<double> neighbors_attract_id(n_attract); // identies of all neighbors in zone of attraction: id, dist, angle to neighbor, newdir of focal fish
	vector<double> neighbors_attract_dist(n_attract);
	vector<double> neighbors_attract_angle(n_attract);
	vector<double> neighbors_attract_turn(n_attract);

	int c_avoid = 0; // counter for array initialization
	int c_align = 0;
	int c_attract = 0;

	for (int i = 0; i < inds.n; i++) {
		if (i != id) {
			double coords[2]; // coords of comparison fish
			coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set focal fish to x = 0 and correct for continous environment
			if (coords[0] > w / 2) {
				coords[0] = coords[0] - w;
			}
			else if (coords[0] < - w / 2) {
				coords[0] = coords[0] + w;
			}
			coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
			if (coords[1] > h / 2) {
				coords[1] = coords[1] - h;
			}
			else if (coords[1] < - h / 2) {
				coords[1] = coords[1] + h;
			}

			double dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id

			// angle between fish id direction vector and fish i pos vector, dot product(direction vector of id, position vector of i)
			double angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;

			// f_all_dotproduct << angle << ",\n";

			if (dist < 200 && angle < (fov / 2)) {
				if (dist < 15) {
					neighbors_avoid_id.at(c_avoid) = i; // id of neighbor fish
					neighbors_avoid_dist.at(c_avoid) = dist; // distance between neighbor and focal fish
					neighbors_avoid_angle.at(c_avoid) = get_angle(inds.coords[id], inds.dir[id], inds.coords[i]); // angle to neighbor

					if (neighbors_avoid_angle.at(c_avoid) > 0) {
						neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) - 180;
					}
					else if (neighbors_avoid_angle.at(c_avoid) < 0) {
						neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + 180;
					}
					else {
						neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + (distrib_choice(rd) - 0.5) * 360;
					}
					neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_turn.at(c_avoid) * (1 / (0.1 * pow((dist), 2) + 4)); // turning angle, change in direction (+/-), weighted with distance

					// f_all_repulsion << neighbors_avoid_turn.at(c_avoid) << ",\n";
					c_avoid++;
				}
				else if (dist < 100) {
					neighbors_align_id.at(c_align) = i; // id of neighbor fish
					neighbors_align_dist.at(c_align) = dist; // distance between neighbor and focal fish

					double dir[2];
					dir[0] = cos(inds.dir[i] * PI / 180) + inds.coords[id][0]; // unit vector of neighbor i, starting at focal fish
					dir[1] = sin(inds.dir[i] * PI / 180) + inds.coords[id][1];

					neighbors_align_angle.at(c_align) = get_angle(inds.coords[id], inds.dir[id], dir); // angle to neighbor
					neighbors_align_turn.at(c_align) = neighbors_align_angle.at(c_align) * (1 / (0.004 * pow((dist - 15), 2) + 2)); //turning angle, weighted with distance

					// f_all_alignment << neighbors_align_turn.at(c_align) << ",\n";
					c_align++;
				}
				else {
					neighbors_attract_id.at(c_attract) = i; // id of neighbor fish
					neighbors_attract_dist.at(c_attract) = dist; // distance between neighbor and focal fish
					neighbors_attract_angle.at(c_attract) = get_angle(inds.coords[id], inds.dir[id], inds.coords[i]); // angle to neighbor
					neighbors_attract_turn.at(c_attract) = neighbors_attract_angle.at(c_attract) * (1 / (0.004 * pow((dist - 200), 2) + 2)); // weighted turning angle

					// f_all_attraction << neighbors_attract_turn.at(c_attract) << ",\n";
					c_attract++;
				}
			}
		}
	}

	double turn = 0;
	if (n_avoid > 0) {
		turn = average_turn(neighbors_avoid_turn);
		// f_mean_repulsion << turn << ",\n";
	}
	else {
		if (n_align > 0 && n_attract > 0) {
			vector<double> comb_al_at(2);
			comb_al_at.at(0) = average_turn(neighbors_align_turn);
			comb_al_at.at(1) = average_turn(neighbors_attract_turn);
			turn = average_turn(comb_al_at);
			// f_mean_combination << turn << ",\n";
		}
		else if (n_align > 0) {
			turn = average_turn(neighbors_align_turn);
			// f_mean_alignment << turn << ",\n";
		}
		else if (n_attract > 0) {
			turn = average_turn(neighbors_attract_turn);
			// f_mean_attraction << turn << ",\n";
		}
		else {
			turn = distrib_rwalk(rd); // random walk if no social interactions
			// f_random_walk << turn << ",\n";
		}
	}
	return turn;
}

double average_turn(vector<double> vec) {
	double x_sum = 0;
	double y_sum = 0;
	for (int i = 0; i < vec.size(); i++) {
		x_sum = x_sum + cos(vec.at(i) * PI / 180);
		y_sum = y_sum + sin(vec.at(i) * PI / 180);
	}
	double x = x_sum / vec.size();
	double y = y_sum / vec.size();

	double angle = acos(x / sqrt(pow(x, 2) + pow(y, 2))) * 180 / PI; // angle to average turn destination from 0° (0-180°)

	if (y < 0) { // in our coordinates, +y is down, -y is up, switched from > to <
		angle = 0 - angle; // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
	}
	else if (y == 0) { // correct 0 angle and 0,0 position
		if (x >= 0) {
			angle = 0;
		}
		else {
			angle = 180;
		}
	}
	return angle;
}

bool state(Individuals inds, int id) {
	bool state = 0;
	if (in_zone(inds, id, 180, 200) >= inds.nfollow && inds.lag[id] == 0) {
		state = 1; // swimming state
	}
	return state;
}

double correct_angle(double dir) {
	if (dir < 0) {
		dir = dir + 360;
	}
	else if (dir >= 360) {
		dir = dir - 360;
	}
	return dir;
}

double get_angle(double* focal_coords, double focal_dir, double* pos) {
	double point[2];
	point[0] = pos[0] - focal_coords[0]; // move point in respect to setting fish id to 0,0
	if (point[0] > w / 2) {
		point[0] = point[0] - w;
	}
	else if (point[0] < - w / 2) {
		point[0] = point[0] + w;
	}
	point[1] = pos[1] - focal_coords[1];
	if (point[1] > h / 2) {
		point[1] = point[1] - h;
	}
	else if (point[1] < - h / 2) {
		point[1] = point[1] + h;
	}

	double angle = acos(point[0] / sqrt(pow(point[0], 2) + pow(point[1], 2))) * 180 / PI; // arccos of dot product between x-axis and point position
	// f_all_dotproduct << angle << ",\n";

	if (point[1] < 0) { // in our coordinates, +y is down, -y is up
		angle = 360 - angle; // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
	}
	else if (point[1] == 0) { // correct 0 angle and 0,0 position
		if (point[0] >= 0) {
			angle = 0;
		}
		else {
			angle = 180;
		}
	}

	angle = angle - focal_dir; // angle difference between id direction and point position, can be negative
	double newangle = correct_angle(angle); // set angle difference between 0 < x =< 360

	if (newangle > 180) {
		newangle = newangle - 360; // left is negative, right is positive; angles between -180 < angle =< 180
	}
	else if (newangle == 180) {
		newangle = (distrib_choice(rd) - 0.5) * 360; // when turn is 180, random left/right decision
	}
	return newangle;
}

double distance(double* coordsa, double* coordsb) {
	double distance = sqrt(pow(coordsa[0] - coordsb[0], 2) + pow(coordsa[1] - coordsb[1], 2));
	return distance;
}

int on_trackbar_nfollow (int, void*) {
	return nf;
}

double on_trackbar_dspeed (int, void*) {
	return dsp;
}

double on_trackbar_speed (int, void*) {
	return sp;
}

int lag(bool state, int lag) {
	if (state == 0) {
		if (lag == 0) {
			lag = 40;
		}
		else {
			lag = lag - 1;
		}
	}
	return lag;
}

void correct_coords(double* coords) {
	if (coords[0] > w) {
		coords[0] = coords[0] - w;
	}
	else if (coords[0] < 0) {
		coords[0] = coords[0] + w;
	}

	if (coords[1] > h) {
		coords[1] = coords[1] - h;
	}
	else if (coords[1] < 0) {
		coords[1] = coords[1] + h;
	}
}

void create_environment(Mat& environment, int quality[][cols]) {
	for (int i = 0; i < cols; i++) {
		for (int u = 0; u < rows; u++) {
			quality[u][i] = distrib_choice(rd);
		}
	}

	for (int i = 0; i < h; i++) {
		for (int u = 0; u < w; u++) {
			int q = quality[i / (h / rows)][u / (w / cols)];
			environment.at<Vec3b>(i, u) = Vec3b(q * 255, q * 255, q * 255);
		}
	}
}

int get_quality(Individuals inds, int id, int quality[][cols]) {
	int x = inds.coords[id][0];
	int y = inds.coords[id][1];

	int q = quality[y / (h / rows)][x / (w / cols)];
	return q;
}
