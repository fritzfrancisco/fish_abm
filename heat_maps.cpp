#include <iostream>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "opencv2/core/core.hpp"
#include <opencv2/imgproc.hpp>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"


#include <stdio.h>

using namespace std;
using namespace cv;

void paint(int x, int y, Mat image);

int main() {

    char hm_window[] = "Social heat map";
	Mat hm_image_g = Mat::zeros(460, 510, CV_8UC3);
    Mat hm_image_invert = Mat::zeros(460, 510, CV_8UC3);
    hm_image_invert = Scalar(255, 255, 255);
    Mat hm_image_out;

    for (int x = 0; x < 410; x++) {
        for (int y = 0; y < 410; y++) {
            paint(x, y, hm_image_g);
        }
    }

    Point triangle[1][3];
    triangle[0][0] = Point(205, 205);
    triangle[0][1] = Point(271, 405);
    triangle[0][2] = Point(139, 405);

    const Point* ppt[1] = { triangle[0] };
    int npt[] = { 3 };

    fillPoly(hm_image_g, ppt, npt, 1, Scalar(0, 0, 0 ), CV_AA);

    for (int x = 430; x < 451; x++) {
        double s = 0.5;
        for (int y = 40; y < 441; y++) {
            hm_image_g.at<Vec3b>(y, x) = Vec3b(s * 510, s * 510, s * 510);
            s = s - (0.5 / 400);
        }
    }

    Point center(205, 205);
    Scalar s(0, 0, 0);
    Scalar w(255, 255, 255);


    circle(hm_image_g, center, 15, w, 1, CV_AA, 0);
    circle(hm_image_g, center, 100, w, 1, CV_AA, 0);
    circle(hm_image_g, center, 200, w, 1, CV_AA, 0);

    for (double i = 0; i < 6; i++) {
        Point pt1(430, 40 + i * 80);
        Point pt2(450, 40 + i * 80);
        Point pt3(455, 45 + i * 80);
        line(hm_image_g, pt1, pt2, s, 1, 8, 0);

        char s_value[10];
        sprintf(s_value, "%g", 0.5 - i / 10);
        putText(hm_image_g, s_value, pt3, 1, 1, w, 1, CV_AA, false);
    }

    line(hm_image_g, Point(205, 220), Point(205, 250), w, 1, CV_AA ,0);
    line(hm_image_g, Point(205, 305), Point(205, 335), w, 1, CV_AA ,0);
    line(hm_image_g, Point(205, 405), Point(205, 435), w, 1, CV_AA ,0);

    putText(hm_image_g, "0.75BL", Point(178, 265), 1, 1, w, 1, CV_AA, false);
    putText(hm_image_g, "5BL", Point(192, 350), 1, 1, w, 1, CV_AA, false);
    putText(hm_image_g, "10BL", Point(188, 450), 1, 1, w, 1, CV_AA, false);

    Point fish(40, 425);
    Point tail(40, 435);
    Point head(40, 415);
    circle(hm_image_g, fish, 2, w, 1, CV_AA, 0);
    arrowedLine(hm_image_g, tail, head, w, 1, CV_AA, 0, 0.2);
    putText(hm_image_g, "Body Length", Point(50, 430), 1, 1, w, 1, CV_AA, false);
    putText(hm_image_g, "Social Force", Point(390, 20), 1, 1, w, 1, CV_AA, false);


    subtract(hm_image_invert, hm_image_g, hm_image_out);
    imshow(hm_window, hm_image_out);

    waitKey(  );

    return(0);
}

void paint(int x, int y, Mat image) {
    double d = sqrt(pow((205 - x), 2) + pow((205 - y), 2));
    double s = 0;
    if (d < 15) {
        s = 1 / (0.1 * pow((d), 2) + 2);
        image.at<Vec3b>(y, x) = Vec3b(s * 510, s * 510, s * 510);
    }
    else if (d >= 15 && d < 100) {
        s = 1 / (0.004 * pow((d - 15), 2) + 2);
        image.at<Vec3b>(y, x) = Vec3b(s * 510, s * 510, s * 510);
    }
    else if (d >= 100 && d < 200) {
        s = 1 / (0.004 * pow((d - 200), 2) + 2);
        image.at<Vec3b>(y, x) = Vec3b(s * 510, s * 510, s * 510);
    }
}
