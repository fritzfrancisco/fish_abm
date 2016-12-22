#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <map>
#include <string>
#include <cmath>
#include <random>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/imgproc.hpp>
#include "opencv2/core/core_c.h"
#include "opencv2/highgui/highgui_c.h"

using namespace std;
using namespace cv;

#define num 60
#define n_cols 250 // must be Structure.latitude * Structure.resolution
#define n_rows 250 // must be Structure.longitude * Structure.resolution

random_device rd;
uniform_int_distribution<int>distrib_choice(0,1);
uniform_int_distribution<int>distrib_error(-10,10);
uniform_int_distribution<int>distrib_randomwalk(-45,45);
uniform_int_distribution<int>distrib_quality(0,10);
uniform_int_distribution<int>distrib_seed(0,0);
uniform_int_distribution<int>distrib_sample(0,35);

class individuals {
public:
int n;
int dim;
int lag[num];
int colors[num][3];         // individual color in RGB
int nfollow;
int sample_rate [num];
int food_intake[num];
double speed[num];          // movement speed
double coords [num][2];     // height n * width dim
double dir [num];           // array of agles between 0 - 360 degrees for each individual
bool sample[num];
double speed_factor[num];
};

class Structure {
public:
int latitude = 100;
int longitude = 100;
int resolution = 1;
int res_factor = 20;

int w = (latitude * resolution * res_factor);
int h = (longitude * resolution * res_factor);
int rows = (longitude * resolution);
int cols = (latitude * resolution);
};

int in_zone(individuals inds,int id, int fov,double r,Structure struc); // fov: field of view as angle
int get_quality(individuals inds, int id, int quality[][n_cols],Structure struc);
int sum_array(int quality[][n_cols],Structure struc);

void move(individuals& inds,Structure struc);
void initialize(individuals& inds,int n, int dim, double speed,Structure struc);
void drawfish(individuals inds, Mat fish_image);
void correct_coords(double* coords,Structure struc);
void create_environment(Mat& environment,int quality[][n_cols], Structure struc);
void feed(individuals& inds,int id, int quality[][n_cols],Mat& environment,Structure struc);
void sample(individuals& inds, int quality[][n_cols], Mat& environment,Structure struc);
void get_environment(Mat& environment, int quality[][n_cols], Structure struc);
void get_speed(individuals& inds,int id,Structure struc);

double correctangle(double dir);
double social(individuals inds,int id,int fov,Structure struc);
double collision(individuals inds, int id, double newdir);
double getangle(double* focal_coords, double focal_dir, double* position,Structure struc);
double averageturn(vector<double> vec, bool comb_weight);

int main(){
        cout << "Created using: OpenCV version 2.4.13" << endl;
        cout << "Currently using: " << "OpenCV version : " << CV_VERSION << "\n";
        cout << "-------------------------------------------" << endl;
        VideoWriter output_video;
        output_video.open("fish_mod.avi", CV_FOURCC('M','P','4','2'), 30, Size(2000, 2000), true);

        if (!output_video.isOpened()) {
                cout << "video writer not initialized\n";
        }
        else {
                cout << "video writer successfully initialized\n";
        }
        // ofstream sum_quality;
        // ofstream ind_quality;
        // sum_quality.open("sum_quality.csv", fstream::in | fstream::out | fstream::app);
        // ind_quality.open("ind_quality.csv", fstream::in | fstream::out | fstream::app);

        srand(time(0));

        Structure fish_struc;
        cout << "environment structure object initialized\n";
        // Creating window for imaging
        char fish_window[] = "Fish";
        // Create matrices
        // Create "food" environment
        Mat fish_image = Mat::zeros(fish_struc.h,fish_struc.w,CV_8UC3);
        Mat environment = Mat::zeros(fish_struc.h,fish_struc.w,CV_8UC3);
        Mat fishinenvironment;
        // imshow(env_window,environment);

        individuals fish;
        initialize(fish,num,2,5,fish_struc); // initialize: object fish,num individuals, dimensions,speed
        cout << "fish individuals object initialized\n";
        // Establish mixed pattern of environment with quality[h]x[w] squares
        int quality[n_rows][n_cols];
        cout << "2d quality array..\n";
        get_environment(environment,quality,fish_struc);

        // sum_quality << sum_array(quality) << ",";
        // sum_quality << "\n";
        cout << "starting simulation..\n";
        for(int z = 0; z < 500; z++) {
                if (z % 100 == 0) {
                        cout << ".. step: " << z << "\n";
                }
                move(fish,fish_struc);
                sample(fish,quality,environment,fish_struc);
                // if(z % 10 == 0) {
                //                 for(int id = 0; id < fish.n; id++) {
                //                                 if(id < (fish.n-1)) {
                //                                                 sum_quality << fish.food_intake[id] << ",";
                //                                 }
                //                                 else{
                //                                                 sum_quality << fish.food_intake[id] << endl;
                //                                 }
                //                 }
                // }
                fish_image = Mat::zeros(fish_struc.h,fish_struc.w,CV_8UC3);
                drawfish(fish,fish_image);
                addWeighted(fish_image,1,environment,0.1,0.0,fishinenvironment);
                // imshow(fish_window,fishinenvironment);

                output_video.write(fishinenvironment);

                // if(z < 1999){
                //  sum_quality << sum_array(quality, fish_struc) << ",";
                // }
                // else{
                //  sum_quality << sum_array(quality, fish_struc) << endl;
                // }
                // waitKey(30);
        }
        // sum_quality.close();
        // ind_quality.close();
        cout << "simulation complete\n";
        return(0);
}

void move(individuals& inds,Structure struc){
        for(int id = 0; id < inds.n; id++) {
                if(inds.lag[id] == 0) {
                        double dir = correctangle(inds.dir[id]);
                        double turn = social(inds,id,330,struc) + distrib_error(rd); // maximum error of movement is +- 5°

                        dir = dir + turn;

                        inds.coords[id][0] = inds.coords[id][0] + inds.speed[id] * cos(dir * M_PI / 180) * inds.speed_factor[id];
                        inds.coords[id][1] = inds.coords[id][1] + inds.speed[id] * sin(dir * M_PI/ 180) * inds.speed_factor[id];
                        get_speed(inds, id, struc);
                        correct_coords(inds.coords[id],struc);
                        inds.dir[id] = dir;
                }
                else{
                        inds.dir[id] = inds.dir[id] + distrib_error(rd);
                }
        }
}

void initialize(individuals& inds, int n,int dim, double speed,Structure struc){
        inds.n = n;
        inds.dim = dim;
        for(int i = 0; i < inds.n; i++) {
                inds.lag[i] = 0; // equals handling time
                inds.speed[i] = speed;
                inds.sample_rate[i] = 0;
                inds.speed_factor[i] = 1;
                inds.food_intake[i] = 0;

                inds.colors[i][0] = 80 + rand()%155;
                inds.colors[i][1] = 80 + rand()%155;
                inds.colors[i][2] = 80 + rand()%155;
                inds.dir[i] = distrib_randomwalk(rd);
                for(int u = 0; u < inds.dim; u++) {
                        if(u == 0) {
                                inds.coords[i][u] = struc.w/2 + rand()%201 - 100;
                        }
                        else{
                                inds.coords[i][u] = struc.h/2 + rand()%201 - 100;
                        }
                }
        }
}

void drawfish(individuals inds, Mat fish_image){
        for(int i = 0; i < inds.n; i++) {
                // Scalar fishcol(5,inds.coords[i][1],inds.coords[i][0]);
                // Scalar fishcol(0,0,255);
                Scalar fishcol(inds.colors[i][0],inds.colors[i][1],inds.colors[i][2]);
                Point center(inds.coords[i][0],inds.coords[i][1]);
                Point tail(inds.coords[i][0] - 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] - 10 * sin(inds.dir[i] * M_PI / 180));
                Point head(inds.coords[i][0] + 10 * cos(inds.dir[i] * M_PI / 180),inds.coords[i][1] + 10 * sin(inds.dir[i] * M_PI / 180));
                circle(fish_image,center,2,fishcol,1,CV_AA,0); // visualize center
                // circle(fish_image,center,100,fishcol,1,CV_AA,0);       // visualize comfort zone
                // circle(fish_image,center,200,fishcol,1,CV_AA,0);       // visualize visual distance
                // char angle[10];
                // sprintf(angle,"%g",inds.dir[i]);
                // putText(fish_image,angle,center,1,1,fishcol,1,8);
                arrowedLine(fish_image,tail,head,fishcol,1,CV_AA,0,0.2); // draw fish
        }
}

int in_zone(individuals inds,int id,int fov,double r,Structure struc){
        // count number of inds in front (x-axis) of id individual
        int c = 0;
        double dirvector[inds.dim]; // angle of fish[id] is facing in respect to original x-axis
        double dist; // absolute distance between fish [id] and fish [i]
        double angle; // angular distance between fish[id] direction vector and fish[i] position vector with fish [id] oriented forward along x-axis
        dirvector[0] = cos(inds.dir[id] * M_PI / 180);
        dirvector[1] = sin(inds.dir[id] * M_PI / 180);
        for (int i = 0; i < inds.n; i++) {
                if(i != id) {
                        double coords[2];
                        coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id x to 0
                        coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id y to 0
                        if(coords[0] > struc.w/2) {
                                coords[0] = coords[0] - struc.w;
                        }
                        else if(coords[0] < -struc.w/2) {
                                coords[0] = coords[0] + struc.w;
                        }
                        if(coords[1] > struc.h/2) {
                                coords[1] = coords[1] - struc.h;
                        }
                        else if(coords[1] < -struc.h/2) {
                                coords[1] = coords[1] + struc.h;
                        }

                        dist = sqrt(pow(coords[0],2) + pow(coords[1],2)); // absolute value
                        angle = acos((dirvector[0]*coords[0]+dirvector[1]*coords[1])/dist)*180/M_PI; // dot prodcut (direction vector of id, position vector of i)
                        if (angle < (fov/2) && dist < r) { // visual field with angle and radius
                                c++;
                        }
                }
        }
        return c;
}

double social(individuals inds, int id, int fov,Structure struc){
        double dirvector[2]; // direction vector (fish id) of lenth 1
        dirvector[0] = cos(inds.dir[id] * M_PI / 180);
        dirvector[1] = sin(inds.dir[id] * M_PI / 180);

        int socialfactor = 2; // switch for behavior: 0 = avoidance, 1 = alignment, 2 = attraction

        int c_avoid = 0; // counter to initialize array
        int c_attract = 0;
        int c_align = 0;

        int n_avoid = in_zone(inds, id, fov, 15,struc); // number of individals to avoid. [id,dist,angle,angle to neighbor]
        int n_align = in_zone(inds, id, fov, 100,struc) - n_avoid; // number of individals to align to [id,dist,angle between directions]
        int n_attract = in_zone(inds, id, fov, 200,struc) - (n_avoid + n_align); // number of individals to be attracted to [id,dist,angle to neighbor, turning angle]

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
        for (int i = 0; i < inds.n; i++) {
                if(i != id) {
                        double coords[2];
                        coords[1] = inds.coords[i][1] - inds.coords[id][1]; // set fish id to y = 0
                        coords[0] = inds.coords[i][0] - inds.coords[id][0]; // set fish id to x = 0
                        if(coords[0] > struc.w/2) {
                                coords[0] = coords[0] - struc.w;
                        }
                        else if(coords[0] < -struc.w/2) {
                                coords[0] = coords[0] + struc.w;
                        }

                        if(coords[1] > struc.h/2) {
                                coords[1] = coords[1] - struc.h;
                        }
                        else if(coords[1] < -struc.h/2) {
                                coords[1] = coords[1] + struc.h;
                        }

                        double dist = sqrt(pow(coords[0], 2) + pow(coords[1], 2)); // calculate absolute distances of respective i to id
                        // angle between fish id direction vector and fish i pos vector, dot product(direction vector of id, position vector of i)
                        double angle =  acos((dirvector[0] * coords[0] + dirvector[1] * coords[1]) / dist) * 180 / M_PI;
                        if (angle < (fov / 2) && dist < 200) {
                                if(dist < 15) {
                                        neighbors_avoid_id.at(c_avoid) = i; // first column: id of fish
                                        neighbors_avoid_dist.at(c_avoid) = dist; // distance between focal fish and neighbor
                                        neighbors_avoid_angle.at(c_avoid) = getangle(inds.coords[id],inds.dir[id],inds.coords[i],struc); // angle to neighbor

                                        if(neighbors_avoid_angle.at(c_avoid) > 0) {
                                                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) - 180;
                                        }
                                        else if(neighbors_avoid_angle.at(c_avoid) < 0) {
                                                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + 180;
                                        }
                                        else{
                                                neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_angle.at(c_avoid) + (distrib_choice(rd)-0.5) * 360;
                                        }
                                        neighbors_avoid_turn.at(c_avoid) = neighbors_avoid_turn.at(c_avoid) * (1 / ( 0.1 * pow((dist),2) + 2 )); // change in angle with algebraic sign (+,-)
                                        c_avoid++;
                                }
                                else if(dist < 100) {
                                        neighbors_align_id.at(c_align) = i; // id of neighbor fish
                                        neighbors_align_dist.at(c_align) = dist; // distance between neighbor and focal fish

                                        double dir[2];
                                        dir[0] = cos(inds.dir[i] * M_PI / 180) + inds.coords[id][0]; // unit vector of neighbor i, starting at focal fish
                                        dir[1] = sin(inds.dir[i] * M_PI / 180) + inds.coords[id][1];

                                        neighbors_align_angle.at(c_align) = getangle(inds.coords[id],inds.dir[id],dir,struc); // angle to neighbor
                                        neighbors_align_turn.at(c_align) = neighbors_align_angle.at(c_align) * (1 / ( 0.004 * pow((dist-15),2) + 2 )); // turning angle, weighted according to distance of fish i
                                        c_align++;
                                }
                                else{
                                        neighbors_attract_id.at(c_attract) = i; // id of neighbor fish
                                        neighbors_attract_dist.at(c_attract) = dist; // distance between neighbor and focal fish
                                        neighbors_attract_angle.at(c_attract) = getangle(inds.coords[id],inds.dir[id],inds.coords[i],struc);
                                        neighbors_attract_turn.at(c_attract) = neighbors_attract_angle.at(c_attract) * (1 / ( 0.004 * pow(( dist-200),2) + 2 )); // turning angle, weighted according to distance fish i
                                        c_attract++;
                                }
                        }
                }
        }
        double turn = 0;
        if(n_avoid > 0) {
                turn = averageturn(neighbors_avoid_turn,0);
        }
        else{
                if(n_align > 0 && n_attract > 0) {
                        vector<double> comb_al_at(2);
                        comb_al_at.at(0) = averageturn(neighbors_align_turn,0);
                        comb_al_at.at(1) = averageturn(neighbors_attract_turn,0);
                        turn = averageturn(comb_al_at,1);
                }
                else if(n_align > 0) {
                        turn = averageturn(neighbors_align_turn,0);
                }
                else if(n_attract > 0) {
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
        for(int i = 0; i < vec.size(); i++) {
                if(comb_weight) {
                        // weighting difference between attraction and alignment
                        x_sum = x_sum + cos(vec.at(i) * M_PI / 180) * (1.7 - i); // unit vector of neighbor i, starting at focal fish
                        y_sum = y_sum + sin(vec.at(i) * M_PI / 180) * (1.7 - i);
                }
                else{
                        x_sum = x_sum + cos(vec.at(i) * M_PI / 180); // unit vector of neighbor i, starting at focal fish
                        y_sum = y_sum + sin(vec.at(i) * M_PI / 180);
                }
        }
        double x = x_sum/vec.size();
        double y = y_sum/vec.size();
        double angle = acos(x / sqrt(pow(x, 2) + pow(y, 2))) * 180 / M_PI; // angle to average turnpoint from x-axis == 0 degrees (0 - 180°)
        if (y < 0) {
                angle = 0 - angle; // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
        }
        else if (y == 0) { // correct 0 angle and 0,0 position
                if (x >= 0) { // reference point in front of focal point.don't change angle
                        angle = 0;
                }
                else{ // reference point behind focal point
                        angle = 180;
                }
        }
        return angle;
}

double correctangle(double dir){
        if (dir < 0) {
                dir = dir + 360;
        }
        else if (dir >= 360) {
                dir = dir - 360;
        }
        return dir;
}

double getangle(double* focal_coords, double focal_dir, double* position,Structure struc){
        // stream cout to file
        double point[2];
        point[0] = position[0] - focal_coords[0]; // move point in respect to setting fish id to 0,0
        point[1] = position[1] - focal_coords[1];

        if(point[0] > struc.w/2) {
                point[0] = point[0] - struc.w;
        }
        else if(point[0] < -struc.w/2) {
                point[0] = point[0] + struc.w;
        }

        if(point[1] > struc.h/2) {
                point[1] = point[1] - struc.h;
        }
        else if(point[1] < -struc.h/2) {
                point[1] = point[1] + struc.h;
        }

        double angle = acos(point[0] / sqrt(pow(point[0], 2) + pow(point[1], 2))) * 180 / M_PI; // arccos of dot product between id direction and point position

        if (point[1] < 0) {
                angle = 360 - angle; // correct angle from 0 =< angle < 180 to 0 =< angle < 360 based on y coordinates
        }
        else if (point[1] == 0) { // correct 0 angle and 0,0 position
                if (point[0] >= 0) {
                        angle = 0;
                }
                else{
                        angle = 180;
                }
        }

        angle = angle - focal_dir; // angle difference between id direction and point position, can be negative
        double newangle = correctangle(angle); // set angle difference between 0 < x =< 360

        if (newangle > 180) {
                newangle = newangle - 360; // left is negative, right is positive; angles between -180 < angle =< 180
        }
        else if(newangle == 180) {
                newangle = (distrib_choice(rd) - 0.5) * 360; // random choice between left or right turn to 180°
        }
        return newangle;
}

void correct_coords(double* coords,Structure struc){

        if(coords[0] > struc.w) {
                coords[0] = coords[0] - struc.w;
        }
        else if(coords[0] < 0) {
                coords[0] = coords[0] + struc.w;
        }
        if(coords[1] > struc.h) {
                coords[1] = coords[1] - struc.h;
        }
        else if(coords[1] < 0) {
                coords[1] = coords[1] + struc.h;
        }
}

void create_environment(Mat& environment,int quality[][n_cols],Structure struc){
        int seed = distrib_seed(rd);
        cout << ".. number of seeds: " << seed << "\n";
        for(int u=0; u<struc.rows; u++) {
                for(int i =0; i<struc.cols; i++) {
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
                for(int u=0; u<struc.rows; u++) {
                        for(int i =0; i<struc.cols; i++) {
                                if(u < struc.rows-1 && i < struc.cols-1) {
                                        if(quality[u+1][i] > p || quality[u+1][i+1] > p || quality[u][i+1]) {
                                                quality[u][i] = distrib_quality(rd);
                                        }
                                }
                        }
                }
                for(int u = struc.rows - 1; u > 0; u--) {
                        for(int i = struc.cols - 1; i > 0; i--) {
                                if(quality[u-1][i] > p || quality[u-1][i-1] > p || quality[u][i-1]) {
                                        quality[u][i] = distrib_quality(rd);
                                }
                        }
                }
        }

        // Color squares according to quality
        for(int x = 0; x < struc.w; x++) {
                for(int y = 0; y < struc.h; y++) {
                        int q_color = quality[y/(struc.h/struc.rows)][x/(struc.w/struc.cols)];
                        // colour depends on maximum quality: 255/maximum quality
                        environment.at<Vec3b>(y,x)= Vec3b(q_color*(255/2),q_color*(255/2),q_color*(255/2));
                }
        }
        cout << ".. seeds successfully spread\n";
}

int get_quality(individuals inds, int id, int quality[][n_cols],Structure struc){
        int x = inds.coords[id][0];
        int y = inds.coords[id][1];
        int q = quality[y/(struc.h/struc.rows)][x/(struc.w/struc.cols)];
        return q;
}

void sample(individuals& inds, int quality[][n_cols], Mat& environment,Structure struc){
        for(int id = 0; id < inds.n; id++) {
                // cout << "in sample:l 537" << "\n";
                int n_avoid = in_zone(inds, id, 330, 15,struc); // number of individals to avoid. [id,dist,angle,angle to neighbor]
                int n_align = in_zone(inds, id, 330, 100,struc) - n_avoid; // number of individals to align to [id,dist,angle between directions]
                int n_attract = in_zone(inds, id, 330, 200,struc) - (n_avoid + n_align); // number of individals to be attracted to [id,dist,angle to neighbor, turning angle]
                bool alone = (n_avoid == 0 && n_align < 5);
                if(inds.lag[id] == 0 && inds.sample_rate[id] > distrib_sample(rd) && alone == 0) {
                        inds.lag[id] = 10;
                        feed(inds,id,quality,environment,struc);
                }
                else if(inds.lag[id] == 0) {
                        inds.sample_rate[id] = inds.sample_rate[id] + 1;
                }
                else{
                        inds.lag[id] = inds.lag[id] - 1;
                }
        }
}

void feed(individuals& inds,int id, int quality[][n_cols],Mat& environment,Structure struc){
        int x = inds.coords[id][0];
        int y = inds.coords[id][1];
        int q =  quality[y/(struc.h/struc.rows)][x/(struc.w/struc.cols)];
        if(q>0) {
                // make speed dependent on food quality sampled
                inds.speed_factor[id] = 0.5;
                quality[y /(struc.h/struc.rows)][x/(struc.w/struc.cols)] = q - 1;
                q = q - 1;
                inds.food_intake[id] = inds.food_intake[id] + 1;
                inds.sample_rate[id] = 30;
        }
        else{
                inds.sample_rate[id] = 0;
        }
        // determining box in quality matrix
        int y_min = y/(struc.h/struc.rows);
        y_min = y_min * (struc.h/struc.rows);
        int y_max = y_min + (struc.h/struc.rows);

        int x_min = x/(struc.w/struc.cols);
        x_min = x_min * (struc.w/struc.cols);
        int x_max = x_min + (struc.w/struc.cols);
        // draw environment boxes as filled rectangles (255/max. quality)
        // color dependent on maximum quality: 255/maximum quality
        rectangle(environment,Point(x_min,y_min),Point(x_max-1,y_max-1),Scalar(q*(255/2),q*(255/2),q*(255/2)),-1,8,0);
}

void get_speed(individuals& inds,int id,Structure struc){
        int count = in_zone(inds,id,180,200,struc);
        inds.speed_factor[id] = 1 + 2 * (double) count / (double) num;
}

int sum_array(int quality[][n_cols],Structure struc){
        int sum = 0;
        for(int i = 0; i < struc.cols; i++) {
                for(int u = 0; u < struc.rows; u++) {
                        sum = sum + quality[u][i];
                }
        }
        return sum;
}

void get_environment(Mat& environment, int quality[][n_cols], Structure struc) {
        ifstream in("quality.csv");

        string line, field;

        vector< vector<string> > array;   // the 2D array
        vector<string> v;                 // array of values for one line only

        while ( getline(in,line) )        // get next line in file
        {
        v.clear();
        stringstream ss(line);

        while (getline(ss,field,','))     // break line into comma delimitted fields
        {
                v.push_back(field);       // add each field to the 1D array
        }
        array.push_back(v);               // add the 1D array to the 2D array
        }

        for (size_t i=0; i<array.size(); ++i)
        {
                for (size_t j=0; j<array[i].size(); ++j)
                {
                        quality[i][j] = stoi(array[i][j]);                         // (separate fields by |) -> write into quality
                }
        }

        for (int y = 0; y < struc.h; y++) {
                for (int x = 0; x < struc.w; x++) {
                        int q = quality[y / (struc.h / struc.rows)][x / (struc.w / struc.cols)];
                        environment.at<Vec3b>(y, x) = Vec3b(q * (255/2), q * (255/2), q * (255/2));
                }
        }
}
