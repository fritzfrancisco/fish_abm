#include <iostream>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <string>
#include <map>
#include <random>
#include <fstream>

using namespace std;

#define n_cols 100 // must be Structure.longitude * resolution
#define n_rows 100 // must be Structure.latitude * resolution

random_device rd;
uniform_int_distribution<int> distrib_quality(0, 2);

class Structure {
	public:
		int latitude = 100;
		int longitude = 100;
		int resolution = 1; // must be 100 > resolution > 1
		int res_factor = 40; // should be: res_factor * resolution = ~ 100

		int w = latitude * resolution * res_factor;
		int h = longitude * resolution * res_factor;
		int cols = latitude * resolution;
		int rows = longitude * resolution;
};

int main() {
  ofstream store_quality;
  store_quality.open("quality.csv", fstream::in | fstream::out | fstream::app);

  int quality[n_cols][n_rows];
  Structure struc;

  int seed = 10; // set number of seeds

  cout << ".. number of seeds: " << seed << "\n";

    for (int u = 0; u < struc.rows; u++) {
    for (int i =0; i < struc.cols; i++) {
      quality[u][i] = 0;
    }
  }

  cout << ".. array initialized with zeros\n";

  if(seed > 0){
    int seed_coords[seed][2];
    uniform_int_distribution<int>distrib_rows(0,struc.rows);
    uniform_int_distribution<int>distrib_cols(0,struc.cols);

    for(int i = 0; i < seed; i++){
    seed_coords[i][0] = distrib_rows(rd);
    seed_coords[i][1] = distrib_cols(rd);
    quality[seed_coords[i][0]][seed_coords[i][1]] = 2;
    }
  }

  cout << ".. seed coordinates initialized\n";

  for(int p = 0; p < 20; p++){
    for(int u = 0; u < struc.rows; u++){
      for(int i = 0; i < struc.cols; i++){
        if (u < struc.rows - 1 && i < struc.cols - 1) {
          if(quality[u + 1][i] > 1 || quality[u + 1][i + 1]  > 1 || quality[u][i + 1] > 1){
            quality[u][i] = distrib_quality(rd);
          }
        }
      }
    }
    for(int u = struc.rows - 1; u > 0; u--){
      for(int i = struc.cols - 1; i > 0; i--){
        if(quality[u - 1][i] > 1 || quality[u - 1][i - 1]  > 1 || quality[u][i - 1] > 1){
          quality[u][i] = distrib_quality(rd);
        }
      }
		}
  }
  cout << ".. seeds successfully spread\n";

  for (int u = 0; u < struc.rows; u++) {
		for (int i = 0; i < struc.cols; i++) {
			if (i == struc.cols - 1) {
				store_quality << quality[u][i] << endl;
			}
			else {
				store_quality << quality[u][i] << ",";
			}
	  }
  }
  store_quality.close();
  cout << "quality array stored to file quality.csv\n";
}
