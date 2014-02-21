#include <ctime>
#include <iostream>

#include "algorithm.h"
#include "mesh.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

clock_t f;
#define START_TIME() f = clock();
#define END_TIME(s) cout << "Operation " << s << " " << clock() - f << "ms" << endl;

void usage() {
    cout << "./Geo <in_obj> <out_obj> <operation_1> <operation_2> ..." << endl;
    cout << "operations:" << endl;
    cout << "l : remove face that the largest component remains" << endl;
    cout << "c : center by average vertex positions" << endl;
    cout << "r : fix entire mesh orientation" << endl;
    cout << "example: ./Geo fengkan_10000.obj fengkan_20000.obj c r l" << endl;
    exit(-1);
}

int main(int argc, char** argv) {
    if (argc < 4)
        usage();
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    clock_t f;
    Mesh m;
    m.load(in);
    in.close();

    int t = 3;
    while (t < argc) {
        START_TIME()
        if (strcmp(argv[t], "l") == 0)
            remove_face_by_largest_component(m);
        else if (strcmp(argv[t], "c") == 0)
            center_positioning_by_averaging_vertex(m);
        else if (strcmp(argv[t], "r") == 0)
            rotate_mesh(m);
        END_TIME(argv[t])
        ++t;
    }

    m.save(out);
    out.close();
    return 0;
}
