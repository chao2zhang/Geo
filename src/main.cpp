#include <ctime>
#include <iostream>

#include "algorithm.h"
#include "object.h"

using namespace std;

clock_t t;
#define START_TIME(s) t = clock(); cout << s << endl;
#define END_TIME(s) cout << s << clock() - t << "ms" << endl;

int main(int argc, char** argv) {
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    clock_t t;
    Object o;

    START_TIME("Loading...")
    o.load(in);
    END_TIME("Loading...");

    count_spikes(o);

    START_TIME("Executing laplacian HC smooth...")
    laplacian_hc_smooth(o,3);
    END_TIME("Executing laplacian HC smooth...")

//    START_TIME("Executing partitioning")
//    partition_by_plane(o,Plane(1, 0, 0, 0));
//    END_TIME()

    START_TIME("Executing positioning...")
    center_positioning(o);
    END_TIME("Executing positioning...")

    START_TIME("Executing shell...")
    shell(o, -0.1);
    END_TIME("Executing shell...")

//    START_TIME("Executing unify normals...")
//    unify_face_normals(o);
//    END_TIME("Executing unify normals...")

    count_spikes(o);

    START_TIME("Writing...")
    o.save(out);
    END_TIME("Writing...")
    in.close();
    out.close();
    return 0;
}
