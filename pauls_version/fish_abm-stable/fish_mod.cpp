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

#define PI 3.14159265
#define num 60
#define n_cols 96 // must be Structure.longitude * resolution
#define n_rows 54 // must be Structure.latitude * resolution

using namespace std;
using namespace cv;

random_device rd;
uniform_int_distribution<int> distrib_choice(0, 1);
uniform_int_distribution<int> distrib_error(-10, 10);
uniform_int_distribution<int> distrib_rwalk(-45, 45);
uniform_int_distribution<int> distrib_quality(0, 10);
uniform_int_distribution<int>distrib_seed(0,0);
uniform_int_distribution<int>distrib_sample(40, 100);

class Individuals {
public:
double coords[num][2];   // height n and width dim
double dir[num];   // fish direction between 0 and 360°
int color[num][3];
int n;
int dim;
double speed[num];   // movement speed
double speed_factor[num];
int sample_rate[num];
int lag[num];
int food_intake[num];
};

class Structure {
public:
int latitude = 48;
int longitude = 27;
int resolution = 2;   // must be 100 > resolution > 1
int res_factor = 20;   // should be: res_factor * resolution = ~ 100

int w = latitude * resolution * res_factor;
int h = longitude * resolution * res_factor;
int cols = latitude * resolution;
int rows = longitude * resolution;
};

void move(Individuals& inds, Structure struc);
void initialize(Individuals& inds, int n, int dim, double speed, Structure struc);
void correct_coords(double* coords, Structure struc);
void drawfish(Individuals inds, Mat image);
void create_environment(Mat& environment, int quality[][n_cols], Structure struc);
void get_environment(Mat& environment, int quality[][n_cols], Structure struc);
void feed(Individuals& inds, int id, int quality[][n_cols], Mat& environment, Structure struc);
void sample(Individuals& inds, int quality[][n_cols], Mat& environment, Structure struc);
void get_speed(Individuals& inds, int id, Structure struc);


int in_zone(Individuals inds, int id, int fov, double r, Structure struc); // fov: field of view
int get_quality(Individuals inds, int id, int quality[][n_cols], Structure struc);
int sum_array(int quality[][n_cols], Structure struc);

double collision(Individuals inds, int id, double newdir);
double correct_angle(double dir);
double get_angle(double* focal_coords, double focal_dir, double* pos, Structure struc);
double social(Individuals inds, int id, int fov, Structure struc);
double average_turn(vector<double> vec, bool comb_weight);

int main() {
								VideoWriter output_video;
								output_video.open("fish_mod.avi", CV_FOURCC('M','P','4','2'), 30, Size(1920, 1080), true);
								if (!output_video.isOpened()) {
								        cout << "video writer not initialized\n";
								}
								else {
								        cout << "video writer successfully initialized\n";
								}

								ofstream sum_quality;
								sum_quality.open("sum_quality.csv", fstream::in | fstream::out | fstream::app);

								srand(time(0));

								Structure fish_struc;

								cout << "environment structure object initialized\n";
								Individuals fish;
								// initializes object fish of class Individuals with num individuals, dimensions, speed and initializes structure
								initialize(fish, num, 2, 5, fish_struc);

								cout << "fish individuals object initialized\n";

								char fish_window[] = "Fish";
								char environment_window[] = "Environment";
								Mat fish_image = Mat::zeros(fish_struc.h, fish_struc.w, CV_8UC3);
								Mat environment = Mat::zeros(fish_struc.h, fish_struc.w, CV_8UC3);
								Mat fish_environment;

								int quality[n_rows][n_cols];

								cout << "2d quality array..\n";
								get_environment(environment, quality, fish_struc);

								// sum_quality << sum_array(quality, fish_struc) << ",";

								cout << "starting simulation..\n";
								for (int z = 0; z < 800; z++) {
																if (z % 100 == 0) {
																								cout << ".. step: " << z << "\n";
																}

																move(fish, fish_struc);
																sample(fish, quality, environment, fish_struc);

																if(z % 10 == 0) {
																								for(int id = 0; id < fish.n; id++) {
																																if(id < (fish.n-1)) {
																																								sum_quality << fish.food_intake[id] << ",";
																																}
																																else{
																																								sum_quality << fish.food_intake[id] << endl;
																																}
																								}
																}

																fish_image = Mat::zeros(fish_struc.h, fish_struc.w, CV_8UC3);
																drawfish(fish, fish_image);
																addWeighted(fish_image, 1, environment, 0.1, 0.0, fish_environment);
																// imshow(fish_window, fish_environment);

																output_video.write(fish_environment);

																// if(z < 4999){
																// 	sum_quality << sum_array(quality, fish_struc) << ",";
																// }
																// else{
																// 	sum_quality << sum_array(quality, fish_struc) << endl;
																// }

																// waitKey(33);

								}
								sum_quality.close();
								cout << "simulation complete\n";
								return(0);
}

void move(Individuals& inds, Structure struc) {
								for (int id = 0; id < inds.n; id++) {
																if (inds.lag[id] == 0) {
																								double dir = correct_angle(inds.dir[id]);

																								double small_error = distrib_error(rd);
																								double turn = social(inds, id, 330, struc) + small_error;

																								dir = dir + turn;

																								inds.coords[id][0] = inds.coords[id][0] + inds.speed[id] * cos(dir * PI / 180) * inds.speed_factor[id];
																								inds.coords[id][1] = inds.coords[id][1] + inds.speed[id] * sin(dir * PI / 180) * inds.speed_factor[id];
																								get_speed(inds, id, struc);
																								correct_coords(inds.coords[id], struc);

																								inds.dir[id] = dir;
																}
																else {
																								inds.dir[id] = inds.dir[id] + distrib_error(rd);
																}
								}
}

void initialize(Individuals& inds, int n, int dim, double speed, Structure struc) {
								inds.n = n;
								inds.dim = dim;
								for (int i = 0; i < inds.n; i++) {
																inds.color[i][0] = 120 + rand() % 120;
																inds.color[i][1] = 120 + rand() % 120;
																inds.color[i][2] = 120 + rand() % 120;
																inds.dir[i] = distrib_rwalk(rd); // direction
																inds.speed[i] = speed;
																inds.speed_factor[i] = 1;
																inds.sample_rate[i] = 100;
																inds.lag[i] = 0; // 17 equals half a second
																inds.food_intake[i] = 0;
																for (int u = 0; u < inds.dim; u++) {
																								if (u == 0) {
																																inds.coords[i][u] = struc.w / 2 + rand() % 201 - 100; // x coord
																								}
																								else {
																																inds.coords[i][u] = struc.h / 2 + rand() % 201 - 100; // y coord
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

int in_zone(Individuals inds, int id, int fov, double r, Structure struc) {
								int c = 0; // count num of inds in front of id
								double dirvector[2]; // direction vector (fish id) of lenth 1
								dirvector[0] = cos(inds.dir[id] * PI / 180);
								dirvector[1] = sin(inds.dir[id] * PI / 180);
								double dist; // absolut distance between fish i and fish id // does not need to be an array
								double angle; // angle between fish id direction vector and fish i pos vector
								for (int i = 0; i < inds.n; i++) {
																if (i != id) {
																								double coords[2];
																								coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0 and correct for continous environment
																								if (coords[0] > struc.w / 2) {
																																coords[0] = coords[0] - struc.w;
																								}
																								else if (coords[0] < -struc.w / 2) {
																																coords[0] = coords[0] + struc.w;
																								}
																								coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
																								if (coords[1] > struc.h / 2) {
																																coords[1] = coords[1] - struc.h;
																								}
																								else if (coords[1] < -struc.h / 2) {
																																coords[1] = coords[1] + struc.h;
																								}

																								dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id:

																								// dot product(direction vector of id, position vector of i)
																								angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;
																								if (angle < (fov / 2) && dist < r) {
																																c++;
																								}
																}
								}
								return c;
}

double social(Individuals inds, int id, int fov, Structure struc) {
								double dirvector[2]; // direction vector (focal fish) of lenth 1
								dirvector[0] = cos(inds.dir[id] * PI / 180);
								dirvector[1] = sin(inds.dir[id] * PI / 180);

								int n_avoid = in_zone(inds, id, fov, 15, struc); // number of neighbors in respective zones
								int n_align = in_zone(inds, id, fov, 100, struc) - n_avoid;
								int n_attract = in_zone(inds, id, fov, 200, struc) - (n_avoid + n_align);

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
																								if (coords[0] > struc.w / 2) {
																																coords[0] = coords[0] - struc.w;
																								}
																								else if (coords[0] < -struc.w / 2) {
																																coords[0] = coords[0] + struc.w;
																								}
																								coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
																								if (coords[1] > struc.h / 2) {
																																coords[1] = coords[1] - struc.h;
																								}
																								else if (coords[1] < -struc.h / 2) {
																																coords[1] = coords[1] + struc.h;
																								}

																								double dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id

																								// angle between fish id direction vector and fish i pos vector, dot product(direction vector of id, position vector of i)
																								double angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / PI;

																								if (dist < 200 && angle < (fov / 2)) {
																																if (dist < 15) {
																																								neighbors_avoid_id.at(c_avoid) = i; // id of neighbor fish
																																								neighbors_avoid_dist.at(c_avoid) = dist; // distance between neighbor and focal fish
																																								neighbors_avoid_angle.at(c_avoid) = get_angle(inds.coords[id], inds.dir[id], inds.coords[i], struc); // angle to neighbor

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
																																								c_avoid++;
																																}
																																else if (dist < 100) {
																																								neighbors_align_id.at(c_align) = i; // id of neighbor fish
																																								neighbors_align_dist.at(c_align) = dist; // distance between neighbor and focal fish

																																								double dir[2];
																																								dir[0] = cos(inds.dir[i] * PI / 180) + inds.coords[id][0]; // unit vector of neighbor i, starting at focal fish
																																								dir[1] = sin(inds.dir[i] * PI / 180) + inds.coords[id][1];

																																								neighbors_align_angle.at(c_align) = get_angle(inds.coords[id], inds.dir[id], dir, struc); // angle to neighbor
																																								neighbors_align_turn.at(c_align) = neighbors_align_angle.at(c_align) * (1 / (0.004 * pow((dist - 15), 2) + 2)); //turning angle, weighted with distance
																																								c_align++;
																																}
																																else {
																																								neighbors_attract_id.at(c_attract) = i; // id of neighbor fish
																																								neighbors_attract_dist.at(c_attract) = dist; // distance between neighbor and focal fish
																																								neighbors_attract_angle.at(c_attract) = get_angle(inds.coords[id], inds.dir[id], inds.coords[i], struc); // angle to neighbor
																																								neighbors_attract_turn.at(c_attract) = neighbors_attract_angle.at(c_attract) * (1 / (0.004 * pow((dist - 200), 2) + 2)); // weighted turning angle
																																								c_attract++;
																																}
																								}
																}
								}

								double turn = 0;
								if (n_avoid > 0) {
																turn = average_turn(neighbors_avoid_turn, 0);
								}
								else {
																if (n_align > 0 && n_attract > 0) {
																								vector<double> comb_al_at(2);
																								comb_al_at.at(0) = average_turn(neighbors_align_turn, 0);
																								comb_al_at.at(1) = average_turn(neighbors_attract_turn, 0);
																								turn = average_turn(comb_al_at, 1);
																}
																else if (n_align > 0) {
																								turn = average_turn(neighbors_align_turn, 0);
																}
																else if (n_attract > 0) {
																								turn = average_turn(neighbors_attract_turn, 0);
																}
																else {
																								turn = distrib_rwalk(rd); // random walk if no social interactions
																}
								}
								return turn;
}

double average_turn(vector<double> vec, bool comb_weight) {
								double x_sum = 0;
								double y_sum = 0;
								for (int i = 0; i < vec.size(); i++) {
																if (comb_weight) {
																								x_sum = x_sum + cos(vec.at(i) * PI / 180) * (2 - i);
																								y_sum = y_sum + sin(vec.at(i) * PI / 180) * (2 - i);
																}
																else {
																								x_sum = x_sum + cos(vec.at(i) * PI / 180);
																								y_sum = y_sum + sin(vec.at(i) * PI / 180);
																}
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

double correct_angle(double dir) {
								if (dir < 0) {
																dir = dir + 360;
								}
								else if (dir >= 360) {
																dir = dir - 360;
								}
								return dir;
}

double get_angle(double* focal_coords, double focal_dir, double* pos, Structure struc) {
								double point[2];
								point[0] = pos[0] - focal_coords[0]; // move point in respect to setting fish id to 0,0
								if (point[0] > struc.w / 2) {
																point[0] = point[0] - struc.w;
								}
								else if (point[0] < -struc.w / 2) {
																point[0] = point[0] + struc.w;
								}
								point[1] = pos[1] - focal_coords[1];
								if (point[1] > struc.h / 2) {
																point[1] = point[1] - struc.h;
								}
								else if (point[1] < -struc.h / 2) {
																point[1] = point[1] + struc.h;
								}

								double angle = acos(point[0] / sqrt(pow(point[0], 2) + pow(point[1], 2))) * 180 / PI; // arccos of dot product between x-axis and point position
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

void correct_coords(double* coords, Structure struc) {
								if (coords[0] > struc.w) {
																coords[0] = coords[0] - struc.w;
								}
								else if (coords[0] < 0) {
																coords[0] = coords[0] + struc.w;
								}

								if (coords[1] > struc.h) {
																coords[1] = coords[1] - struc.h;
								}
								else if (coords[1] < 0) {
																coords[1] = coords[1] + struc.h;
								}
}

void create_environment(Mat& environment, int quality[][n_cols], Structure struc) {
								int seed = distrib_seed(rd);

								cout << ".. number of seeds: " << seed << "\n";

								for (int u = 0; u < struc.rows; u++) {
																for (int i =0; i < struc.cols; i++) {
																								quality[u][i] = 0;
																}
								}

								cout << ".. array initialized with zeros\n";

								if(seed > 0) {
																int seed_coords[seed][2];
																uniform_int_distribution<int>distrib_rows(0,struc.rows);
																uniform_int_distribution<int>distrib_cols(0,struc.cols);

																for(int i = 0; i < seed; i++) {
																								seed_coords[i][0] = distrib_rows(rd);
																								seed_coords[i][1] = distrib_cols(rd);
																								quality[seed_coords[i][0]][seed_coords[i][1]] = 10;
																}
								}

								cout << ".. seed coordinates initialized\n";

								for(int p = 0; p < 20; p++) {
																for(int u = 0; u < struc.rows; u++) {
																								for(int i = 0; i < struc.cols; i++) {
																																if (u < struc.rows - 1 && i < struc.cols - 1) {
																																								if(quality[u + 1][i] > p || quality[u + 1][i + 1]  > p || quality[u][i + 1] > p) {
																																																quality[u][i] = distrib_quality(rd);
																																								}
																																}
																								}
																}
																for(int u = struc.rows - 1; u > 0; u--) {
																								for(int i = struc.cols - 1; i > 0; i--) {
																																if(quality[u - 1][i] > 5 || quality[u - 1][i - 1]  > 5 || quality[u][i - 1] > 5) {
																																								quality[u][i] = distrib_quality(rd);
																																}
																								}
																}
								}

								for (int y = 0; y < struc.h; y++) {
																for (int x = 0; x < struc.w; x++) {
																								int q = quality[y / (struc.h / struc.rows)][x / (struc.w / struc.cols)];
																								environment.at<Vec3b>(y, x) = Vec3b(q * 25.5, q * 25.5, q * 25.5);
																}
								}

								cout << ".. seeds successfully spread\n";
}

int get_quality(Individuals inds, int id, int quality[][n_cols], Structure struc) {
								int x = inds.coords[id][0];
								int y = inds.coords[id][1];

								int q = quality[y / (struc.h / struc.rows)][x / (struc.w / struc.cols)];
								return q;
}

void feed(Individuals& inds, int id, int quality[][n_cols], Mat& environment, Structure struc) {
								int x = inds.coords[id][0];
								int y = inds.coords[id][1];

								int q = quality[y / (struc.h / struc.rows)][x / (struc.w / struc.cols)];
								inds.sample_rate[id] = 10 * q;

								if (q > 0) {
																inds.speed_factor[id] = inds.speed_factor[id]/(q);

																quality[y / (struc.h / struc.rows)][x / (struc.w / struc.cols)] = q - 1;
																q = q - 1;
								}

        int y_min = y / (struc.h / struc.rows);
        y_min = y_min * (struc.h / struc.rows);
        int y_max = y_min + (struc.h / struc.rows) - 1;

        int x_min = x / (struc.w / struc.cols);
        x_min = x_min * (struc.w / struc.cols);
        int x_max = x_min + (struc.w / struc.cols) - 1;

        Point pt_min(x_min, y_min);
        Point pt_max(x_max, y_max);
        Scalar q_color(q * (255 / 10), q * (255 / 10), q * (255 / 10));
        rectangle(environment, pt_min, pt_max, q_color, -1, 8, 0);
}

void sample(Individuals& inds, int quality[][n_cols], Mat& environment, Structure struc) {
								for (int id = 0; id < inds.n; id++) {
																int n_avoid = in_zone(inds, id, 330, 15, struc); // number of individals to avoid. [id,dist,angle,angle to neighbor]
																int n_align = in_zone(inds, id, 330, 100, struc) - n_avoid; // number of individals to align to [id,dist,angle between directions]
																int n_attract = in_zone(inds, id, 330, 200, struc) - (n_avoid + n_align); // number of individals to be attracted to [id,dist,angle to neighbor, turning angle]
																bool alone = (n_avoid == 0 && n_align == 0);

																if (inds.lag[id] == 0 && inds.sample_rate[id] > distrib_sample(rd) && alone == 0) {
																								inds.lag[id] = 10;
																								feed(inds, id, quality, environment, struc);
																								inds.food_intake[id] = inds.food_intake[id] + 1;
																}
																else if (inds.lag[id] == 0) {
																								inds.sample_rate[id] = inds.sample_rate[id] + 1;
																}
																else {
																								inds.lag[id] = inds.lag[id] - 1;
																								inds.sample_rate[id] = inds.sample_rate[id] + 1;
																}
								}
}

void get_speed(Individuals& inds, int id, Structure struc) {
								int count = in_zone(inds, id, 180, 200, struc);
								inds.speed_factor[id] = 1 + 2 * (double) count/ (double) num;
}

int sum_array(int quality[][n_cols], Structure struc) {
								int sum = 0;
								for (int i = 0; i < struc.rows; i++) {
																for (int u = 0; u < struc.cols; u++) {
																								sum = sum + quality[i][u];
																}
								}
								return sum;
}

void get_environment(Mat& environment, int quality[][n_cols], Structure struc) {
								using namespace std;

								ifstream in("quality.csv");

								string line, field;

								vector< vector<string> > array; // the 2D array
								vector<string> v;         // array of values for one line only

								while ( getline(in,line) ) // get next line in file
								{
								v.clear();
								stringstream ss(line);

								while (getline(ss,field,','))  // break line into comma delimitted fields
								{
																v.push_back(field); // add each field to the 1D array
								}

								array.push_back(v);  // add the 1D array to the 2D array
								}

								// print out what was read in

								for (size_t i=0; i<array.size(); ++i)
								{
																for (size_t j=0; j<array[i].size(); ++j)
																{
																								quality[i][j] = stoi(array[i][j]); // (separate fields by |) -> write into quality
																}
								}

								for (int y = 0; y < struc.h; y++) {
								        for (int x = 0; x < struc.w; x++) {
								                int q = quality[y / (struc.h / struc.rows)][x / (struc.w / struc.cols)];
								                environment.at<Vec3b>(y, x) = Vec3b(q * 25.5, q * 25.5, q * 25.5);
								        }
								}
}
