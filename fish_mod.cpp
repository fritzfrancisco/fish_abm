#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#define w 500
#define h 500
#define PI 3.14159265
#define num 4

using namespace std;
using namespace cv;

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

void move(Individuals& inds);
void initialize(Individuals& inds, int n, int dim, double speed, double dspeed, int nfollow);
void drawfish(Individuals inds, Mat image);
int follow(Individuals inds, int id, int fov); // fov: field of view
bool state(Individuals inds, int id);
double collision(Individuals inds, int id, double newdir);
double correctangle(double dir);
double getangle(Individuals inds, int id, double point[2]);
double social(Individuals inds, int id, int fov);

int main(){
	
	srand(time(0));

	char fish_window[] = "Fish";
	Mat fish_image = Mat::zeros(h, w, CV_8UC3);
		
	Individuals fish;
	// initializes object fish of class Individuals with num individuals, dimensions, speed, dspeed, nfollow 
	initialize(fish, num, 2, 3, 4, 1000);
	
	for (int z = 0; z < 666; z++){
		cout << "step " << z << "\n";
		Mat fish_image = Mat::zeros(h, w, CV_8UC3);
		drawfish(fish, fish_image);
		move(fish);
		imshow(fish_window, fish_image);
		waitKey(  );

	}
	return(0);
}

void move(Individuals& inds){
	for (int id = 0; id < inds.n; id++){
		double dir = inds.dir[id];
		
		inds.state[id] = state(inds, id);
		double socialdir = social(inds, id, 330);
 
		if(inds.state[id] == 0 && socialdir == inds.dir[id]){
			socialdir = socialdir + rand() % 91 - 45; // feeding, new fish direction
			if(inds.lag[id] == 0){
				inds.lag[id] = 10;
			}
			else{
				inds.lag[id] = inds.lag[id] - 1;			
			}
		}		
	//	correctangle(dir);	

		double coldir = collision(inds, id, socialdir);		

		inds.coords[id][0] = inds.coords[id][0] + (inds.state[id] * inds.dspeed + 1) * inds.speed * cos(coldir * PI / 180);
		inds.coords[id][1] = inds.coords[id][1] + (inds.state[id] * inds.dspeed + 1) * inds.speed * sin(coldir * PI / 180);
		
		inds.dir[id] = coldir;
	}
}

void initialize(Individuals& inds, int n, int dim, double speed, double dspeed, int nfollow){
	inds.n = n;
	inds.dim = dim;
	inds.speed = speed;
	inds.dspeed = dspeed;
	inds.nfollow = nfollow;
	for (int i = 0; i < inds.n; i++){
		inds.lag[i] = 0;
		inds.color[i][0] = 120 + rand() % 120;
		inds.color[i][1] = 120 + rand() % 120;
		inds.color[i][2] = 120 + rand() % 120;
		inds.dir[i] = 0; // rand() % 360; // direction
		for (int u = 0; u < inds.dim; u++){
			if (u == 0) {
				inds.coords[i][u] = w / 2 + rand() % 201 - 100; // x coord
			}
			else{
				inds.coords[i][u] = h / 2 + rand() % 201 - 100; // y coord
			}
		} 
	}
}

void drawfish(Individuals inds, Mat image){
	for (int i = 0; i < inds.n; i++){
		Scalar fishcolor(inds.color[i][0], inds.color[i][1], inds.color[i][2]);
		Point center(inds.coords[i][0], inds.coords[i][1]);
		Point tail(inds.coords[i][0] - 10 * cos(inds.dir[i] * PI / 180), inds.coords[i][1] - 10 * sin(inds.dir[i] * PI / 180));
		Point head(inds.coords[i][0] + 10 * cos(inds.dir[i] * PI / 180), inds.coords[i][1] + 10 * sin(inds.dir[i] * PI / 180));
		circle(image, center, 2, fishcolor, 1, CV_AA, 0);
		circle(image, center, 30, fishcolor, 1, CV_AA, 0);
		circle(image, center, 200, fishcolor, 1, CV_AA, 0);
		arrowedLine(image, tail, head, fishcolor, 1, CV_AA, 0, 0.2);
	}
}

int follow(Individuals inds, int id, int fov){
	int c =  0; // count num of inds in front of id
	double dirvector[2]; // direction vector (fish id) of lenth 1
	dirvector[0] = cos(inds.dir[id] * PI / 180);
	dirvector[1] = sin(inds.dir[id] * PI / 180);
	double dist; // absolut distance between fish i and fish id // does not need to be an array
	double angle; // angle between fish id direction vector and fish i pos vector
	for (int i = 0; i < inds.n; i++){
		double coords[2];
		coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0
		coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
		dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id
		
		if (i != id){
			// dot product(direction vector of id, position vector of i)
			angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI; 
		}
		else{
			angle = 180; // never see yourself!	
		}

		if (angle < (fov / 2) && dist < 200){
			c++;		
		}

	}
	return c;
}

double social(Individuals inds, int id, int fov){
	double dirvector[2]; // direction vector (fish id) of lenth 1
	dirvector[0] = cos(inds.dir[id] * PI / 180);
	dirvector[1] = sin(inds.dir[id] * PI / 180);
	double angle; // angle between fish id direction vector and fish i pos vector
	double dist; // absolut distance between fish i and fish id
	int c = follow(inds, id, fov); // counter for fish in visual radius
	double repulsion; // repulsion value of other fish..
	double coordsum[2] = {0, 0};
	int socialfactor = 2; // switch for avoidance behavior, 0 for no avoidance
	double newdir = inds.dir[id]; // angle after avoidance behavior
	double nearestvisdist = 200; // initializes distance of nearest seeable neighbour
	int nearestvisid; // id of nearest seeable neighbor
	double avoidsum = 0; // sum of all repulsion values for weighting turning angle

	// check distances for all other fish and determine social factor (0 = avoid, 1 = correct angle, 2 = reconnect with others)
	for (int i = 0; i < inds.n; i++){
		double coords[2]; // coords of comparison fish
		repulsion = 0;
		coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0
		coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
		dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id

		if (i != id){
			// dot product(direction vector of id, position vector of i)
			angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;
			if (dist < nearestvisdist){
				nearestvisdist = dist;
				nearestvisid = i;			
			}
		}
		else{
			angle = 180; // never see yourself!	
		}

		if (angle < (fov / 2) && dist < 200){
			repulsion = 1 / (0.01 * pow(dist, 2) + 1); // f(x) = 1 / (0.001 * x^2 + 1), also see google
		}
		
		if (angle < (fov / 2) && dist < 30){
			socialfactor = 0;
		}
		//else if (angle < (fov / 2) && dist < 30){
		//	socialfactor = 0;
		//}
		
		if (i != id){
			coordsum[0] = coordsum[0] + coords[0] * repulsion;
			coordsum[1] = coordsum[1] + coords[1] * repulsion;
			avoidsum = avoidsum + repulsion;		
		}
	}
	
	// get a new direction according to social factor
	double leftcoords[2]; // new possible coords after left turn
	double rightcoords[2]; // new possible coords after right turn
	
	if (socialfactor == 0){
		//cout << "avoid; ";
		if (c != 0){
			coordsum[0] = coordsum[0] / c;
			coordsum[1] = coordsum[1] / c;
			avoidsum = avoidsum / c;
		}

		leftcoords[0] = (state(inds, id) * inds.dspeed + 1) * inds.speed * cos((inds.dir[id] - 90) * PI / 180);
		leftcoords[1] = (state(inds, id) * inds.dspeed + 1) * inds.speed * sin((inds.dir[id] - 90) * PI / 180);

		rightcoords[0] = (state(inds, id) * inds.dspeed + 1) * inds.speed * cos((inds.dir[id] + 90) * PI / 180);
		rightcoords[1] = (state(inds, id) * inds.dspeed + 1) * inds.speed * sin((inds.dir[id] + 90) * PI / 180);

		double avoiddistl = sqrt(pow(leftcoords[0] - coordsum[0], 2) + pow(leftcoords[1] - coordsum[1], 2));
		double avoiddistr = sqrt(pow(rightcoords[0] - coordsum[0], 2) + pow(rightcoords[1] - coordsum[1], 2));

		if (avoiddistl > avoiddistr){
			newdir = inds.dir[id] - 90 * avoidsum;
		}
		else if (avoiddistl == avoiddistr){
			newdir = inds.dir[id] + (rand() % 2 - 0.5 ) * 180 * avoidsum;	
		}
		else{
			newdir = inds.dir[id] + 90 * avoidsum;
		}
	}

	if (socialfactor == 2 && c > 0){
		cout << "I am: " << id << "; nearest vis id: " << nearestvisid << "\n";
		cout << "nearest vis id coords: " << inds.coords[nearestvisid][0] << ", " << inds.coords[nearestvisid][1] << "\n";
		newdir = newdir + getangle(inds, id, inds.coords[nearestvisid]);
		cout << "nearest vis id coords: " << inds.coords[nearestvisid][0] << ", " << inds.coords[nearestvisid][1] << "\n";
	}
	return newdir;	
}

bool state(Individuals inds, int id){
	bool state = 0;
	if (follow(inds, id, 180) >= inds.nfollow && inds.lag[id] == 0){ 
		state = 1; // swimming state
	}
	return state;
}

double collision(Individuals inds, int id, double newdir){
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
	while(turnleft && turnright){
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
	if (turnleft == 0 && turnright == 0){
		if (rand() % 2 == 0){
			newdir = newdirleft;
		}
		else{
			newdir = newdirright;					
		}
	}
	else if (turnleft == 0 && turnright == 1){
		newdir = newdirleft;
	}
	else{
		newdir = newdirright;		
	}
	return newdir;
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

double getangle(Individuals inds, int id, double point[2]){
	point[0] = point[0] - inds.coords[id][0]; // move point in respect to setting fish id to 0,0
	point[1] = point[1] - inds.coords[id][1];
	double angle = acos(point[0] / sqrt(pow(point[0], 2) + pow(point[1], 2))) * 180 / PI; // arccos of dot product between id direction and point position

	if (point[1] > 0){
		angle = 360 - angle; // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
	}
	else if (point[1] == 0){ // correct 0 angle and 0,0 position
		if (point[0] >= 0){
			angle = 0;		
		}
		else{
			angle = 180;		
		}
	}
	
	angle = angle - inds.dir[id]; // angle difference between id direction and point position, can be negative
	double newangle = correctangle(angle); // set angle difference between 0 < x =< 360


	if (newangle > 180){
		newangle = newangle - 360; // left is negative, right is positive; angles between -180 < angle =< 180
	}

	return newangle;
}
