#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdio.h>

#define w 1600
#define h 800
#define PI 3.14159265
#define num 5

const int sp_max = 10;
int sp = 2;
const int dsp_max = 20;
int dsp;
const int nf_max = num;
int nf;

using namespace std;
using namespace cv;

class Individuals {
	public:
		double coords[num][2]; // height n and width dim
		double dir[num]; // fish direction between 0 and 360°
		double socialdir[num];
		int lag[num]; // with length n
		bool state[num]; // 0 for random, 1 for directed
		int color[num][3];
		int n;
		int dim;
		double speed; // movement speed
		double dspeed; //dspeed + 1 = swimming speed / feeding speed
		int nfollow;
		int fov;
		double r;
};

void move(Individuals& inds);
void initialize(Individuals& inds, int n, int dim, double speed, double dspeed, int nfollow);
void drawfish(Individuals inds, Mat image);
int follow(Individuals inds, int id, int fov, double r); // fov: field of view; r: radius
bool state(Individuals inds, int id);
double collision(Individuals inds, int id, double newdir);
double correctangle(double dir);
double get_angle(Individuals inds, int id, double* point);
// double social(Individuals inds, int id);
double get_distance(double* coords_a, double* coords_b);
int lag(bool state, int lag);

//new function for social behavior
int get_social_state(Individuals inds, int id);
double get_avoid_angle(Individuals inds, int id);
double get_align_angle(Individuals inds, int id);
double get_attract_angle(Individuals inds, int id);
double socialize(Individuals inds, int id);

double on_trackbar_dspeed(int, void*);
double on_trackbar_speed(int, void*);
int on_trackbar_nfollow(int, void*);

int main() {
	srand(time(0));

	char fish_window[] = "Fish";
	Mat fish_image = Mat::zeros(h, w, CV_8UC3);

	Individuals fish;
	fish.speed = 1;
	fish.dspeed = 1;
	fish.nfollow = num;

	// initializes object fish of class Individuals with num individuals, dimensions, speed, dspeed, nfollow
	initialize(fish, num, 2, fish.speed, fish.dspeed, fish.nfollow);

	for (int z = 0; z < 2000; z++) {

		Mat fish_image = Mat::zeros(h, w, CV_8UC3);
		drawfish(fish, fish_image);

		move(fish);
		imshow(fish_window, fish_image);
		waitKey( 30 );

	char speed_trackbar[10];  // create trackbar in "Fish" window for changing speed
    sprintf(speed_trackbar,"%g",fish.speed);
    createTrackbar("Speed","Fish",&sp,sp_max);
    fish.speed = on_trackbar_speed(fish.speed,0);

    // char dspeed_trackbar[10];  // create trackbar in "Fish" window for changing speed
    // sprintf(dspeed_trackbar,"%g",fish.dspeed);
    // createTrackbar("DSpeed","Fish",&dsp,dsp_max);
    // fish.dspeed = on_trackbar_dspeed(fish.dspeed,0);
	//
    // char nfollow_trackbar[10];  // create trackbar in "Fish" window for changing speed
    // sprintf(nfollow_trackbar,"%i",fish.nfollow);
    // createTrackbar("N-Follow","Fish",&nf,nf_max);
    // fish.nfollow = on_trackbar_nfollow(fish.nfollow,0);

	}
	return(0);
}

void move(Individuals& inds) {
	for (int id = 0; id < num; id++) {
		inds.socialdir[id] = correctangle(socialize(inds, id));

	}

	for (int id = 0; id < inds.n; id++) {
		double dir = correctangle(inds.dir[id] + inds.socialdir[id]);

		inds.state[id] = state(inds, id);

		if (inds.state[id] == 0 && dir == inds.dir[id]) {
			// dir = dir + rand() % 91 - 45; // feeding, new fish direction
			dir = correctangle(dir);
		}
		// dir = dir + rand() % 21 - 10;
		inds.lag[id] = lag(inds.state[id], inds.lag[id]);

		double coldir = collision(inds, id, dir);

		inds.coords[id][0] = inds.coords[id][0] + (inds.state[id] * inds.dspeed + 1) * inds.speed * cos(coldir * PI / 180);
		inds.coords[id][1] = inds.coords[id][1] + (inds.state[id] * inds.dspeed + 1) * inds.speed * sin(coldir * PI / 180);

		inds.dir[id] = coldir;
	}
}

void initialize(Individuals& inds, int n, int dim, double speed, double dspeed, int nfollow) {
	inds.n = n;
	inds.dim = dim;
	inds.speed = speed;
	inds.dspeed = dspeed;
	inds.nfollow = nfollow;
	inds.fov = 330;
	inds.r = 200;
	for (int i = 0; i < inds.n; i++) {
		inds.lag[i] = 0;
		inds.color[i][0] = 120 + rand() % 120;
		inds.color[i][1] = 120 + rand() % 120;
		inds.color[i][2] = 120 + rand() % 120;
		inds.dir[i] = 0; // rand() % 360; // direction
		inds.socialdir[i] = 0;
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
		// circle(image, center, 200 / 6, fishcolor, 1, CV_AA, 0);
		circle(image, center, 200, fishcolor, 1, CV_AA, 0);
		circle(image, center, 100, fishcolor, 1, CV_AA, 0);
		arrowedLine(image, tail, head, fishcolor, 1, CV_AA, 0, 0.2);
		// char idtext[10];
		// sprintf(idtext, "%d", i);
		// putText(image, idtext, head, 1, 1, fishcolor, 1, CV_AA, false );
	}
}

int follow(Individuals inds, int id, int fov, double r) {
	int c =  0; // count num of inds in front of id
	double dirvector[2]; // direction vector (fish id) of lenth 1
	dirvector[0] = cos(inds.dir[id] * PI / 180);
	dirvector[1] = sin(inds.dir[id] * PI / 180);
	double dist; // absolut distance between fish i and fish id // does not need to be an array
	double angle; // angle between fish id direction vector and fish i pos vector
	for (int i = 0; i < inds.n; i++) {
		double coords[2];
		coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0
		coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
		dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id

		if (i != id) {
			// dot product(direction vector of id, position vector of i)
			angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;
		}
		else {
			angle = 180; // never see yourself!
		}

		if (angle < (fov / 2) && dist < r) {
			c++;
		}

	}
	return c;
}

// double social(Individuals inds, int id) {
// 	double dirvector[2]; // direction vector (fish id) of lenth 1
// 	dirvector[0] = cos(inds.dir[id] * PI / 180);
// 	dirvector[1] = sin(inds.dir[id] * PI / 180);
// 	int c = follow(inds, id, inds.fov, 200); // counter for fish in visual radius
// 	int socialfactor = 2; // switch for avoidance behavior, 0 for no avoidance
// 	double newdir = inds.dir[id]; // angle after avoidance behavior
// 	double nearestvisdist = 200; // initializes distance of nearest seeable neighbour
// 	int nearestvisid; // id of nearest seeable neighbor
// 	double avoidsum = 0; // sum of all repulsion values for weighting turning angle
//
// 	// check distances for all other fish and determine social factor (0 = avoid, 1 = correct angle, 2 = reconnect with others)
// 	for (int i = 0; i < inds.n; i++) {
// 		double coords[2]; // coords of comparison fish
// 		double angle; // angle between fish id direction vector and fish i pos vector
// 		double dist; // absolut distance between fish i and fish id
// 		coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0
// 		coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
// 		dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id
//
// 		if (i != id) {
// 			// dot product(direction vector of id, position vector of i)
// 			angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;
// 			if (dist < nearestvisdist) {
// 				nearestvisdist = dist;
// 				nearestvisid = i;
// 			}
// 		}
// 		else {
// 			angle = 180; // never see yourself!
// 		}
//
// 		if (angle < (inds.fov / 2) && dist < 30) {
// 			socialfactor = 0;
// 		}
// 		else if (angle < (inds.fov / 2) && dist < 100) {
// 			socialfactor = 1;
// 		}
// 	}
//
// 	// get a new direction according to social factor
// 	if (socialfactor == 0) {
// 		double avoiddir;
// 		if (get_angle(inds, id, inds.coords[nearestvisid]) < 0) {
// 			avoiddir = get_angle(inds, id, inds.coords[nearestvisid]) + 180;
// 		}
// 		else if (get_angle(inds, id, inds.coords[nearestvisid]) == 0) {
// 			if (rand() % 2 == 0) {
// 				avoiddir = get_angle(inds, id, inds.coords[nearestvisid]) + 180;
// 			}
// 			else {
// 				avoiddir = get_angle(inds, id, inds.coords[nearestvisid]) - 180;
// 			}
// 		}
// 		else {
// 			avoiddir = get_angle(inds, id, inds.coords[nearestvisid]) - 180;
// 		}
//
// 		newdir = newdir + (1 / (0.01 * pow((nearestvisdist), 2) + 4)) * avoiddir;
// 	}
// 	else if (socialfactor == 1) {
// 		double dir[2];
// 		dir[0] = cos(inds.dir[nearestvisid] * PI / 180); // unit vector of nearestvisid
// 		dir[1] = sin(inds.dir[nearestvisid] * PI / 180);
//
// 		dir[0] = dir[0] + inds.coords[id][0]; // shift with coords of id
// 		dir[1] = dir[1] + inds.coords[id][1];
//
// 		newdir = newdir + (1 / (0.01 * pow((nearestvisdist - 30), 2) + 2)) * get_angle(inds, id, dir);
// 	}
// 	else if (socialfactor == 2 && c > 0) {
// 		newdir = newdir + (1 / (0.01 * pow((nearestvisdist - 200), 2) + 2)) * get_angle(inds, id, inds.coords[nearestvisid]);
// 	}
// 	return newdir;
// }

bool state(Individuals inds, int id) {
	bool state = 0;
	if (follow(inds, id, 180, 200) >= inds.nfollow && inds.lag[id] == 0) {
		state = 1; // swimming state
	}
	return state;
}

double collision(Individuals inds, int id, double newdir) {
	double bound = 20; // boundary in which to detect the border

	double leftcoords[2]; // new possible coords after left turn
	double rightcoords[2]; // new possible coords after right turn

	leftcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdir * PI / 180); // initialize new coords after turn with original values
	leftcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdir * PI / 180);

	rightcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdir * PI / 180);
	rightcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdir * PI / 180);

	double newdirleft = newdir; // new direction after turning (left/right), initialized with original direction
	double newdirright = newdir;

	// check collision after turn (left or right, 0 for no collision, 1 for collision); starting with original direction..
	bool turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
	bool turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));

	// .. then adjust turn direction and check collision again; repeat until one collision statement (left/right) is 0 (no collision):
	while(turnleft && turnright) {
		newdirleft = newdirleft - 10;
		leftcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirleft * PI / 180); // new x coord after left turn
		leftcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirleft * PI / 180); // new y coord

		newdirright = newdirright + 10;
		rightcoords[0] = inds.coords[id][0] + (state(inds, id) * inds.dspeed + 1) * inds.speed * cos(newdirright * PI / 180); // new x coord after right turn
		rightcoords[1] = inds.coords[id][1] + (state(inds, id) * inds.dspeed + 1) * inds.speed * sin(newdirright * PI / 180); // new y coord

		turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
		turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));
	}

	// choose new direction based on false collision statement (no collision in this direction):
	if (turnleft == 0 && turnright == 0) {
		if (rand() % 2 == 0) {
			newdir = newdirleft;
		}
		else {
			newdir = newdirright;
		}
	}
	else if (turnleft == 0 && turnright == 1) {
		newdir = newdirleft;
	}
	else {
		newdir = newdirright;
	}
	return newdir;
}

double correctangle(double dir) {
	if (dir < 0) {
		dir = dir + 360;
	}
	else if (dir >= 360) {
		dir = dir - 360;
	}
	return dir;
}

double get_angle(Individuals inds, int id, double* point) {
	double pos[2];
	pos[0] = point[0] - inds.coords[id][0]; // move point in respect to setting fish id to 0,0
	pos[1] = point[1] - inds.coords[id][1];

	double angle = acos(pos[0] / sqrt(pow(pos[0], 2) + pow(pos[1], 2))) * 180 / PI; // arccos of dot product between id direction and point position

	if (pos[1] < 0) {
		angle = 360 - angle; // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
	}
	else if (pos[1] == 0) { // correct 0 angle and 0,0 position
		if (pos[0] >= 0) {
			angle = 0;
		}
		else {
			angle = 180;
		}
	}

	angle = angle - inds.dir[id]; // angle difference between id direction and point position, can be negative
	double newangle = correctangle(angle); // set angle difference between 0 < x =< 360


	if (newangle > 180) {
		newangle = newangle - 360; // left is negative, right is positive; angles between -180 < angle =< 180
	}

	return newangle;
}

double get_distance(double* coords_a, double* coords_b) {
	double distance = sqrt(pow(coords_a[0] - coords_b[0], 2) + pow(coords_a[1] - coords_b[1], 2));
	return distance;
}

int on_trackbar_nfollow(int, void*) {
	return nf;
}

double on_trackbar_dspeed(int, void*) {
	return dsp;
}

double on_trackbar_speed(int, void*) {
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

// social behavior:
// function 1: get social state based on neighbors
// function 2: avoid behavior returns turning angle respecting all neighbors in zone of repulsion
// function 3: alignment also respecting all neighbors -> angle
// function 4: preference for all (distant) neighbors -> angle
// function 5: contains functions 1-4, combines angles if 3 and 4 both return an angle -> returns final angle

// function 1: social state
int get_social_state(Individuals inds, int id) {
	int social_state = 4; // switch for social states, initialized with 4: NO social interactions

	int n_avoid = follow(inds, id, inds.fov, inds.r / 6); // counter for fish in zone of repulsion of focal fish
	int n_align = follow(inds, id, inds.fov, inds.r / 2) - n_avoid; // counter for fish in zone of alignment
	int n_attract = follow(inds, id, inds.fov, inds.r) - n_avoid - n_align; // counter for fish in zone of attraction

	if (n_attract > 0) {
		social_state = 3; // ONLY attraction to distant fish
	}

	if (n_attract > 0 && n_align > 0) {
		social_state = 2; // BOTH attraction and alignment to neighbors
	}

	if (n_attract == 0 && n_align > 0) {
		social_state = 1; // ONLY alignment to neighbors
	}

	if (n_avoid > 0) {
		social_state = 0; // ONLY avoidance of neighbors in zone of repulsion
	}

	return social_state;
}

// function 2: calculate avoid turn angle based on all neighbors
double get_avoid_angle(Individuals inds, int id) {
	double dir_focal[2]; // direction (unit) vector of focal fish
	dir_focal[0] = cos(inds.dir[id] * PI / 180);
	dir_focal[1] = sin(inds.dir[id] * PI / 180);

	int n_avoid = follow(inds, id, inds.fov, inds.r / 6); // counter for fish in zone of repulsion of focal fish

	double pos_neighbors[n_avoid][2]; // position vectors of all fish in zone of repusion
	int u = 0; //counter for for loop array initialization
	for (int i = 0; i < inds.n; i++) {
		if ((i != id) && (get_distance(inds.coords[id], inds.coords[i]) < inds.r / 6) && (get_angle(inds, id, inds.coords[i]) < inds.fov / 2)) {
			pos_neighbors[u][0] = inds.coords[i][0];
			pos_neighbors[u][1] = inds.coords[i][1];
			u++;
		}
	}

	double angle_n[n_avoid]; // angles to all neighbors in zone of repulsion; if turn should be weighted based on distance, here is the possibility
	for (int i = 0; i < n_avoid; i++) {
		double angle = get_angle(inds, id, pos_neighbors[i]);
		if (angle < 0) {
			angle_n[i] = angle + 180;
		}
		else if (angle == 0) {
			if (rand() % 2 == 0) {
				angle_n[i] = angle + 180;
			}
			else {
				angle_n[i] = angle - 180;
			}
		}
		else {
			angle_n[i] = angle - 180;
		}
	}

	double dir_n_turn[n_avoid][2]; // calculate direction vectors of turns for all neighbors and respective angles
	for (int i = 0; i < n_avoid; i++) {
		dir_n_turn[i][0] = cos(angle_n[i] * PI / 180);
		dir_n_turn[i][1] = sin(angle_n[i] * PI / 180);
	}

	double dir_turn[2]; // calculate mean direction vector
	dir_turn[0] = 0;
	dir_turn[1] = 0;
	for (int i = 0; i < n_avoid; i++) {
		dir_turn[0] = dir_turn[0] + dir_n_turn[i][0];
		dir_turn[1] = dir_turn[1] + dir_n_turn[i][1];
	}
	dir_turn[0] = dir_turn[0] / n_avoid;
	dir_turn[1] = dir_turn[1] / n_avoid;

	double angle_turn = acos(dir_turn[0]) * 180 / PI; // same in degree
	if (dir_turn[1] < 0) {
		angle_turn = - angle_turn;
	}

	double turn = inds.dir[id] - angle_turn;
	return turn;
}

// function 3: align to all neighbors in alignment zone
double get_align_angle(Individuals inds, int id) {
	double dir_focal[2]; // direction (unit) vector of focal fish
	dir_focal[0] = cos(inds.dir[id] * PI / 180);
	dir_focal[1] = sin(inds.dir[id] * PI / 180);

	int n_align = follow(inds, id, inds.fov, inds.r / 2) - follow(inds, id, inds.fov, inds.r / 6); // counter for fish in zone of alignment

	double angle_n[n_align]; // direction angles of all fish in zone of alignment
	int u = 0; //counter for for loop array initialization
	for (int i = 0; i < inds.n; i++) {
		if (get_distance(inds.coords[id],inds.coords[i]) >= inds.r / 6 && get_distance(inds.coords[id],inds.coords[i]) < inds.r / 2 && get_angle(inds, id, inds.coords[i]) < inds.fov / 2) {
			angle_n[u] = inds.dir[i];
			u++;
		}
	}

	double dir_n_turn[n_align][2]; // calculate direction vectors of turns for all neighbors and respective angles
	for (int i = 0; i < n_align; i++) {
		dir_n_turn[i][0] = cos(angle_n[i] * PI / 180);
		dir_n_turn[i][1] = sin(angle_n[i] * PI / 180);
	}

	double dir_turn[2]; // calculate mean direction vector
	dir_turn[0] = 0;
	dir_turn[1] = 0;
	for (int i = 0; i < n_align; i++) {
		dir_turn[0] = dir_turn[0] + dir_n_turn[i][0];
		dir_turn[1] = dir_turn[1] + dir_n_turn[i][1];
	}
	dir_turn[0] = dir_turn[0] / n_align;
	dir_turn[1] = dir_turn[1] / n_align;

	double angle_turn = acos(dir_turn[0]) * 180 / PI; // same in degree
	if (dir_turn[1] < 0) {
		angle_turn = - angle_turn;
	}

	double turn = inds.dir[id] + angle_turn;
	return turn;
}

// function 4: get attracted by all neighbors in attraction zone
double get_attract_angle(Individuals inds, int id) {
	double dir_focal[2]; // direction (unit) vector of focal fish
	dir_focal[0] = cos(inds.dir[id] * PI / 180);
	dir_focal[1] = sin(inds.dir[id] * PI / 180);

	int n_attract = follow(inds, id, inds.fov, inds.r) - follow(inds, id, inds.fov, inds.r / 2); // counter for fish in zone of attraction

	double pos_neighbors[n_attract][2]; // position vectors of all fish in zone of repusion
	int u = 0; //counter for for loop array initialization
	for (int i = 0; i < inds.n; i++) {
		if (get_distance(inds.coords[id],inds.coords[i]) >= inds.r / 2 && get_distance(inds.coords[id],inds.coords[i]) < inds.r && get_angle(inds, id, inds.coords[i]) < inds.fov / 2) {
			pos_neighbors[u][0] = inds.coords[i][0];
			pos_neighbors[u][1] = inds.coords[i][1];
			u++;
		}
	}

	double angle_n[n_attract]; // angles to all neighbors in zone of repulsion; if turn should be weighted based on distance, here is the possibility
	for (int i = 0; i < n_attract; i++) {
		double angle = get_angle(inds, id, pos_neighbors[i]);
		if (angle < 0) {
			angle_n[i] = angle - 180;
		}
		else if (angle == 0) {
			if (rand() % 2 == 0) {
				angle_n[i] = angle + 180;
			}
			else {
				angle_n[i] = angle - 180;
			}
		}
		else {
			angle_n[i] = angle + 180;
		}
	}

	double dir_n_turn[n_attract][2]; // calculate direction vectors of turns for all neighbors and respective angles
	for (int i = 0; i < n_attract; i++) {
		dir_n_turn[i][0] = cos(angle_n[i] * PI / 180);
		dir_n_turn[i][1] = sin(angle_n[i] * PI / 180);
	}

	double dir_turn[2]; // calculate mean direction vector
	dir_turn[0] = 0;
	dir_turn[1] = 0;
	for (int i = 0; i < n_attract; i++) {
		dir_turn[0] = dir_turn[0] + dir_n_turn[i][0];
		dir_turn[1] = dir_turn[1] + dir_n_turn[i][1];
	}
	dir_turn[0] = dir_turn[0] / n_attract;
	dir_turn[1] = dir_turn[1] / n_attract;

	double angle_turn = acos(dir_turn[0]) * 180 / PI; // same in degree
	if (dir_turn[1] < 0) {
		angle_turn = - angle_turn;
	}

	double turn = inds.dir[id] - angle_turn;
	return turn;
}

// function 5, container function for social behavior and turn angle calculations
double socialize(Individuals inds, int id) {
	int social_state = get_social_state(inds, id);
	double turn = 0;
	double turn_align;
	double turn_attract;

	switch (social_state) {
		case 0:
			// turn = get_avoid_angle(inds, id);
			// cout << "only avoid angle: " << turn << "\n";
			break;

		case 1:
			turn = get_align_angle(inds, id);
			// cout << "only align angle: " << turn << "\n";
			break;

		case 2:
			turn_align = get_align_angle(inds, id);
			turn_attract = get_attract_angle(inds, id);

			double dir_turn[2]; // direction (unit) vector of alignment turn
			dir_turn[0] = (cos(turn_align * PI / 180) + cos(turn_attract * PI / 180)) / 2;
			dir_turn[1] = (sin(turn_align * PI / 180) + sin(turn_attract * PI / 180)) / 2;

			turn = acos(dir_turn[0]) * 180 / PI; // same in degree
			if (dir_turn[1] < 0) {
				turn = - turn;
			}
			// cout << "combined angle: " << turn << "\n";
			break;

		case 3:
			turn = get_attract_angle(inds, id);
			// cout << "only attract angle: " << turn << "\n";
			break;

		case 4:
			// turn = rand() % 91 - 45; // random walk if no social interactions
			// cout << "random walk mode\n";
			break;
	}
	return turn / 20;
}
