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
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    clock_t f;
    Mesh m;

    START_TIME("Loading...")
    m.load(in);
    END_TIME("Loading...");

    START_TIME("Counting spikes...")
    count_spike(m);
    END_TIME("Counting spikes...");

//    START_TIME("Executing laplacian HC smooth...")
//    laplacian_hc_smooth(m, 5);
//    END_TIME("Executing laplacian HC smooth...")

    START_TIME("Executing positioning...")
    center_positioning_by_averaging_vertex(m);
    END_TIME("Executing positioning...")

//    START_TIME("Executing shell...")
//    mesh_offset(m, -0.05);
//    END_TIME("Executing shell...")

    START_TIME("Executing partitioning...")
    partition_by_plane(m, Plane(0.586, 0.980, -0.993, 1));
    END_TIME("Executing partitioning..")

    START_TIME("Executing partitioning...")
    partition_by_plane(m, Plane(0.837, -0.234, 0.263, -0.8));
    END_TIME("Executing partitioning..")

    START_TIME("Executing projecting...")
    project_by_plane(m, Plane(0.586, 0.980, -0.993, 1.1));
    END_TIME("Executing projecting...")

    START_TIME("Executing removing...")
    remove_face_by_plane(m, Plane(0.837, -0.234, 0.263, -0.8));
    END_TIME("Executing removing...")

    START_TIME("Executing filling...")
    fill_max_border_face_by_plane(m, Plane(0.586, 0.980, -0.993, 1.1));
    END_TIME("Executing filling...")

    START_TIME("Executing unify normals...")
    unify_face_normals(m);
    END_TIME("Executing unify normals...")

    START_TIME("Writing...")
    m.save(out);
    END_TIME("Writing...")
    in.close();
    out.close();
    return 0;
}
