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
    cout << "u : unify face normal orientation" << endl;
    cout << "a <a> <b> <c> <d> : partition by plane ax + by + cz + d = 0" << endl;
    cout << "p <a> <b> <c> <d> : project by plane ax + by + cz + d = 0" << endl;
    cout << "f <a> <b> <c> <d> : fill max border face by plane ax + by + cz + d = 0" << endl;
    cout << "e <a> <b> <c> <d> : remove face by plane ax + by + cz + d = 0" << endl;
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
        char* cur_op = argv[t];
        if (strcmp(argv[t], "l") == 0)
            remove_face_by_largest_component(m);
        else if (strcmp(argv[t], "c") == 0)
            center_positioning_by_averaging_vertex(m);
        else if (strcmp(argv[t], "r") == 0)
            rotate_mesh(m);
        else if (strcmp(argv[t], "u") == 0)
            unify_face_normals(m);
        else if (strcmp(argv[t], "a") == 0) {
            Plane p;
            p.a = atof(argv[++t]);
            p.b = atof(argv[++t]);
            p.c = atof(argv[++t]);
            p.d = atof(argv[++t]);
            partition_by_plane(m, p);
        } else if (strcmp(argv[t], "p") == 0) {
            Plane p;
            p.a = atof(argv[++t]);
            p.b = atof(argv[++t]);
            p.c = atof(argv[++t]);
            p.d = atof(argv[++t]);
            project_by_plane(m, p);
        } else if (strcmp(argv[t], "f") == 0) {
            Plane p;
            p.a = atof(argv[++t]);
            p.b = atof(argv[++t]);
            p.c = atof(argv[++t]);
            p.d = atof(argv[++t]);
            fill_max_border_face_by_plane(m, p);
        } else if (strcmp(argv[t], "e") == 0) {
            Plane p;
            p.a = atof(argv[++t]);
            p.b = atof(argv[++t]);
            p.c = atof(argv[++t]);
            p.d = atof(argv[++t]);
            remove_face_by_plane(m, p);
        } else
            usage();
        END_TIME(cur_op)
        ++t;
    }

    m.save(out);
    out.close();
    return 0;
}
