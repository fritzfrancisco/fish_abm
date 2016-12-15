g++ video_test.cpp -o video_test.out -std=c++11 `pkg-config --cflags --libs opencv`
./video_test.out
rm video_test.out
