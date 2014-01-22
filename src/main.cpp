#include <ctime>
#include <iostream>

#include "algorithm.h"
#include "object.h"

using namespace std;

int main(int argc, char** argv) {
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    clock_t t;
    Object o;

    t = clock();
    cout << "Loading..." << endl;
    o.load(in);
    cout << "Loaded " << clock() - t << "ms" << endl;

    count_spikes(o);

    t = clock();
    cout << "Executing laplacian HC smooth" << endl;
    laplacian_hc_smooth(o,3);
    cout << "Executed laplacian HC smooth " << clock() - t << "ms" << endl;

//    t = clock();
//    cout << "Executing partitioning" << endl;
//    partition_by_plane(o,Plane(1, 0, 0, 0));
//    cout << "Executed partitioning " << clock() - t << "ms" << endl;

    t = clock();
    cout << "Executing positioning" << endl;
    center_positioning(o);
    cout << "Executed positioning " << clock() - t << "ms" << endl;

    t = clock();
    cout << "Executing shell" << endl;
    shell(o, -0.1);
    cout << "Executed shell " << clock() - t << "ms" << endl;

    t = clock();
    cout << "Executing unify normals" << endl;
    unify_face_normals(o);
    cout << "Executed unify normals " << clock() - t << "ms" << endl;

    count_spikes(o);

    t = clock();
    cout << "Writing..." << endl;
    o.save(out);
    cout << "Writed " << clock() - t << "ms" << endl;
    in.close();
    out.close();
    return 0;
}
