#include "verlet.h"
#include <iostream>

int main(int argc, char *argv[]) {
    const int frames = 30;
    int nb_atoms = 10;
    Positions_t positions(3, nb_atoms);
    positions(2, 1) = 1.0;

    for (int i = 0; i < frames; i++) {
    }
}