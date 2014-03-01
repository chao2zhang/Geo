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
    cout << "\tl : remove face that the largest component remains" << endl;
    cout << "\tc : center by average vertex positions" << endl;
    cout << "\tr : fix entire mesh orientation" << endl;
    cout << "\tu : unify face normal orientation" << endl;
    cout << "\ty <y0> : partition by plane y = y0" << endl;
    cout << "\tz <z0> : partition by plane z = z0" << endl;
    cout << "\tya <ratio_y> : ratio_y in (0, 1), partition by plane y = ymin * (1 - ratio) + ymax * ratio" << endl;
    cout << "\tza <ratio_z> : ratio_z in (0, 1), partition by plane z = zmin * (1 - ratio) + zmax * ratio" << endl;
    cout << "\tzh <ratio> : partition by plane z = zmax - ratio * (ymax - ymin)" << endl;
    cout << "\txr <angle> : rotate around x-axis by angle" << endl;
    cout << "\tyr <angle> : rotate around y-axis by angle" << endl;
    cout << "\tzr <angle> : rotate around z-axis by angle" << endl;
    cout << "\tyd : remove face by plane y = ymin" << endl;
    cout << "\tzp <z0> : project by plane z = zmin - z0" << endl;
    cout << "\tzf : fill max border face by plane z = zmin" << endl;
    cout << "\ta : analyze along z-axis" << endl;
    cout << "\ts <file> : save to file" << endl;
    cout << "return value: 0 for success, 1 for failure" << endl;
    cout << "example: ./Geo fengkan_10000.obj c r l s fengkan_20000.obj" << endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
    if (argc <= 3)
        usage();
    ifstream in(argv[1]);
    clock_t f;
    Mesh m;
    m.load(in);
    in.close();

    int t = 2;
    while (t < argc) {
        START_TIME()
        char* cur_op = argv[t];
        if (strcmp(argv[t], "l") == 0)
            remove_face_by_largest_component(m);
        else if (strcmp(argv[t], "c") == 0)
            center_positioning_by_averaging_vertex(m);
        else if (strcmp(argv[t], "r") == 0)
            auto_rotate_mesh(m);
        else if (strcmp(argv[t], "u") == 0)
            unify_face_normals(m);
        else if (strcmp(argv[t], "xr") == 0) {
            float a = atof(argv[++t]);
            rotate_mesh(m, Point3f(1, 0, 0), deg_to_rad(a));
        } else if (strcmp(argv[t], "yr") == 0) {
            float a = atof(argv[++t]);
            rotate_mesh(m, Point3f(0, 1, 0), deg_to_rad(a));
        } else if (strcmp(argv[t], "zr") == 0) {
            float a = atof(argv[++t]);
            rotate_mesh(m, Point3f(0, 0, 1), deg_to_rad(a));
        } else if (strcmp(argv[t], "y") == 0)
            partition_by_plane(m, Plane(0, -1, 0, atof(argv[++t])));
        else if (strcmp(argv[t], "z") == 0)
            partition_by_plane(m, Plane(0, 0, -1, atof(argv[++t])));
        else if (strcmp(argv[t], "ya") == 0) {
            float r = atof(argv[++t]);
            partition_by_plane(m, Plane(0, -1, 0, m.bbox[0].x[1] * (1 - r) + m.bbox[1].x[1] * r));
        } else if (strcmp(argv[t], "za") == 0) {
            float r = atof(argv[++t]);
            partition_by_plane(m, Plane(0, 0, -1, m.bbox[0].x[2] * (1 - r) + m.bbox[1].x[2] * r));
        } else if (strcmp(argv[t], "zh") == 0) {
            float r = atof(argv[++t]);
            partition_by_plane(m, Plane(0, 0, -1, m.bbox[1].x[2] - r * (m.bbox[1].x[1] - m.bbox[0].x[1])));
        } else if (strcmp(argv[t], "yd") == 0)
            remove_face_by_plane(m, Plane(0, -1, 0, m.bbox[0].x[1]));
        else if (strcmp(argv[t], "zp") == 0)
            project_by_plane(m, Plane(0, 0, -1, m.bbox[0].x[2] - atof(argv[++t])));
        else if (strcmp(argv[t], "zf") == 0)
            fill_max_border_face_by_plane(m, Plane(0, 0, -1, m.bbox[0].x[2]));
        else if (strcmp(argv[t], "a") == 0)
            analyze_z(m);
        else if (strcmp(argv[t], "s") == 0) {
            ofstream out(argv[++t]);
            m.save(out);
            out.close();
        } else
            usage();
        END_TIME(cur_op)
        ++t;
    }
    return EXIT_SUCCESS;
}
