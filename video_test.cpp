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

using namespace std;
using namespace cv;

int main () {

    Mat test = Mat::zeros(1920, 1080, CV_8UC3);

    int fcc = CV_FOURCC('D', 'I', 'V', '3');

    VideoWriter output_video;
    output_video = VideoWriter("test.AVI", -1, 30.0, test.size(), true);
    if (!output_video.isOpened()) {
        cout << "video writer not initialized\n";
    }
    else {
        cout << "video writer successfully initialized\n";
    }



    for (int i = 0; i < 100; i++) {
        output_video << test;
    }

    return(0);
}
