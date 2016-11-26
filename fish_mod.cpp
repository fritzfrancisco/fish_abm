#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#define w 400
#define PI 3.14159265

using namespace std;
using namespace cv;

class Individuals {
	public:
		double coords[20][2]; // height n and width dim
		double dir[20]; // fish direction between 0 and 360°
		int lag[20]; // with length n
		bool state[20]; // 0 for random, 1 for directed
		int color[20][3];
		int n;
		int dim;
};

void move(Individuals& inds);
void printcoords(Individuals inds);
void initialize(Individuals& inds, int n, int dim);
void drawfish(Individuals inds, Mat image);
int follow(Individuals inds, int id);

int main(){
	
	srand(time(0));

	char fish_window[] = "Fish";
	Mat fish_image = Mat::zeros(w, w * 2, CV_8UC3);
		
	Individuals fish;
	initialize(fish, 20, 2);
	
	for (int z = 0; z < 200; z++){
		
		Mat fish_image = Mat::zeros(w, w * 2, CV_8UC3);
		drawfish(fish, fish_image);
		move(fish);
		imshow(fish_window, fish_image);
		waitKey( 30 );

	}
	return(0);
}



void move(Individuals& inds){
	for (int i = 0; i < inds.n; i++){
		if (follow(inds, i) > 1 && inds.lag[i] == 0){ 
			inds.state[i] = 1; // swimming state
		}
		else{
			inds.state[i] = 0; // feeding state
		}
	}

	for (int i = 0; i < inds.n; i++){
		if (inds.state[i] == 1 && inds.lag[i] == 0){ 
			for (int u = 0; u < inds.dim - 1; u++){ //only x dim
				inds.coords[i][u] = inds.coords[i][u] + 20;
			}
		}
		else{
			if (inds.lag[i] == 0){
					inds.lag[i] = inds.lag[i] + 20;
				}
				else{
					inds.lag[i] = inds.lag[i] - 1;
				}
			for (int u = 0; u < inds.dim; u++){
				inds.coords[i][u] = inds.coords[i][u] + rand() % 5 - 2;
			}
		}
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
		inds.lag[i] = 0;
		inds.color[i][0] = 120 + rand() % 120;
		inds.color[i][1] = 120 + rand() % 120;
		inds.color[i][2] = 120 + rand() % 120;
		inds.dir[i] = 0; // rand() % 360; // direction
		for (int u = 0; u < inds.dim; u++){
			if (u == 0) {
				inds.coords[i][u] = w / 8 + rand() % 100 - 50; // x coord
			}
			else{
				inds.coords[i][u] = w / 2 + rand() % 100 - 50; // y coord
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
		}else{
			angle[i] = 180; // never see yourself!	
		}

		if (angle[i] < 90){
			c++;		
		}

	}
	return c;
}


