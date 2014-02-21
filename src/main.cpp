#include <ctime>
#include <iostream>

#include "algorithm.h"
#include "mesh.h"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

clock_t f;
#define START_TIME(s) f = clock(); cout << s << endl;
#define END_TIME(s) cout << s << clock() - f << "ms" << endl;

int main(int argc, char** argv) {
    if (argc < 3)
        usage();
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    clock_t f;
    Mesh m;

    START_TIME("Loading...")
    m.load(in);
    END_TIME("Loading...");

    START_TIME("Executing positioning...")
    center_positioning_by_averaging_vertex(m);
    END_TIME("Executing positioning...")

    START_TIME("Executing rotating...")
    rotate_mesh(m);
    END_TIME("Executing rotating...")

    START_TIME("Executing remaining...")
    remove_face_by_largest_component(m);
    END_TIME("Executing remaining...")

    START_TIME("Writing...")
    m.save(out);
    END_TIME("Writing...")
    in.close();
    out.close();
    return 0;
}
