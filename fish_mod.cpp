#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#define w 700
#define h 700
#define PI 3.14159265
#define speed 2

using namespace std;
using namespace cv;

class Individuals {
	public:
		double coords[20][2]; // height n and width dim
		double dir[20]; // fish direction between 0 and 360°
		int lag[20]; // with length n
		// bool state[20]; // 0 for random, 1 for directed
		int color[20][3];
		int n;
		int dim;
};

void move(Individuals& inds);
void printcoords(Individuals inds);
void initialize(Individuals& inds, int n, int dim);
void drawfish(Individuals inds, Mat image);
int follow(Individuals inds, int id);
bool state(Individuals inds, int id);
double border(Individuals inds, int id, double newdir);
void correctangle(double& dir);

int main(){
	
	srand(time(0));

	char fish_window[] = "Fish";
	Mat fish_image = Mat::zeros(h, w, CV_8UC3);
		
	Individuals fish;
	initialize(fish, 20, 2);
	
	for (int z = 0; z < 900; z++){
		Mat fish_image = Mat::zeros(h, w, CV_8UC3);
		drawfish(fish, fish_image);
		move(fish);
		imshow(fish_window, fish_image);
		waitKey( 30 );

	}
	return(0);
}



void move(Individuals& inds){
	for (int i = 0; i < inds.n; i++){
		double dir = inds.dir[i];
		if(state(inds, i) == 0){
			dir = inds.dir[i] + rand() % 91 - 45; // feeding, new fish direction
			if(inds.lag[i] == 0){
				inds.lag[i] = 10;
			}
			else{
				inds.lag[i] = inds.lag[i] - 1;			
			}
		}		
		correctangle(dir);		
		double newdir = border(inds, i, dir);		

		inds.coords[i][0] = inds.coords[i][0] + (state(inds, i) * 4 + 1) * speed * cos(newdir * PI / 180);
		inds.coords[i][1] = inds.coords[i][1] + (state(inds, i) * 4 + 1) * speed * sin(newdir * PI / 180);
		
		inds.dir[i] = newdir;		
	}
}

void printcoords(Individuals inds){
	cout << "individuals are at coordinates:\n";
	for (int i = 0; i < inds.n; i++){
		cout << inds.coords[i][0] << ", " << inds.coords[i][1] << "\n";
	}
}

void initialize(Individuals& inds, int n, int dim){
	inds.n = n;
	inds.dim = dim;
	for (int i = 0; i < inds.n; i++){
		inds.lag[i] = 10;
		inds.color[i][0] = 120 + rand() % 120;
		inds.color[i][1] = 120 + rand() % 120;
		inds.color[i][2] = 120 + rand() % 120;
		inds.dir[i] = rand() % 360; // direction
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
		circle(image, center, 2, fishcolor, 1, 8, 0);
		arrowedLine(image, tail, head, fishcolor, 1, 8, 0, 0.1);
	}
}

int follow(Individuals inds, int id){
	int c =  0; // count number of inds in front of id
	double idcoords[2] = {inds.coords[id][0], inds.coords[id][1]};
	for (int i = 0; i < inds.n; i++){
		double dist[20]; // absolut distance between fish i and fish id
		double dirvector[20][2]; // direction vector of lenth 1
		double angle[20]; // angle between fish id direction vector and fish i pos vector

		inds.coords[i][0] = inds.coords[i][0] - idcoords[0]; // set fish id to x = 0
		inds.coords[i][1] = inds.coords[i][1] - idcoords[1]; // set fish id to y = 0
		dist[i] = sqrt(pow(inds.coords[i][0], 2) + pow(inds.coords[i][1], 2));
		dirvector[i][0] = cos(inds.dir[i] * PI / 180);
		dirvector[i][1] = sin(inds.dir[i] * PI / 180);
		
		if (i != id){
			angle[i] =  acos((dirvector[i][0] * inds.coords[i][0] + dirvector[i][1] * inds.coords[i][1]) / dist[i]) * 180 / PI;
		}
		else{
			angle[i] = 180; // never see yourself!	
		}

		if (angle[i] < 90 && dist[i] < 200){
			c++;		
		}

	}
	return c;
}

bool state(Individuals inds, int id){
	bool state = 0;
	if (follow(inds, id) > 2 && inds.lag[id] == 0){ 
		state = 1; // swimming state
	}
	return state;
}

double border(Individuals inds, int id, double newdir){
	double bound = 20;
	double leftcoords[2];	
	double rightcoords[2];	
	leftcoords[0] = inds.coords[id][0] + (state(inds, id) * 4 + 1) * speed * cos(newdir * PI / 180); // "new" x coord after left turn
	leftcoords[1] = inds.coords[id][1] + (state(inds, id) * 4 + 1) * speed * sin(newdir * PI / 180); // "new" y coord
	rightcoords[0] = inds.coords[id][0] + (state(inds, id) * 4 + 1) * speed * cos(newdir * PI / 180); // "new" x coord after right turn
	rightcoords[1] = inds.coords[id][1] + (state(inds, id) * 4 + 1) * speed * sin(newdir * PI / 180); // "new" y coord

	double newdirleft = newdir;
	double newdirright = newdir;

	bool turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
	bool turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));

	while(turnleft && turnright){
		newdirleft = newdirleft - 1;	
		leftcoords[0] = inds.coords[id][0] + (state(inds, id) * 4 + 1) * speed * cos(newdirleft * PI / 180); // "new" x coord after left turn
		leftcoords[1] = inds.coords[id][1] + (state(inds, id) * 4 + 1) * speed * sin(newdirleft * PI / 180); // "new" y coord

		newdirright = newdirright + 1;	
		rightcoords[0] = inds.coords[id][0] + (state(inds, id) * 4 + 1) * speed * cos(newdirright * PI / 180); // "new" x coord after right turn
		rightcoords[1] = inds.coords[id][1] + (state(inds, id) * 4 + 1) * speed * sin(newdirright * PI / 180); // "new" y coord

		turnleft = (leftcoords[0] < bound || leftcoords[0] > (w - bound) || leftcoords[1] < bound || leftcoords[1] > (h - bound));
		turnright = (rightcoords[0] < bound || rightcoords[0] > (w - bound) || rightcoords[1] < bound || rightcoords[1] > (h - bound));
	}
	if (turnright){
		newdir = newdirleft;
	}
	else{
		newdir = newdirright;		
	}
	return newdir;
}

void correctangle(double& dir){
	if (dir < 0){
		dir = dir + 360;	
	}
	else if (dir > 360){
		dir = dir - 360;
	}
}
