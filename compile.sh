g++ fish_mod.cpp -o opencv.out `pkg-config --cflags --libs opencv`
./opencv.out
rm opencv.out
