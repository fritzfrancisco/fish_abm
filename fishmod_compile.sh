g++ fish_mod.cpp -o fish_abm.out -std=c++11 `pkg-config --cflags --libs opencv`
./fish_abm.out
rm fish_abm.out

# g++ heatmap.cpp -o heatmap.out `pkg-config --cflags --libs opencv`
# ./heatmap.out
# rm heatmap.out
