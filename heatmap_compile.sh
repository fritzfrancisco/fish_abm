g++ heatmap.cpp -o heatmap.out -std=c++11 `pkg-config --cflags --libs opencv`
./heatmap.out
rm heatmap.out
