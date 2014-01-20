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

    t = clock();
    cout << "Executing laplacian HC smooth" << endl;
    laplacian_hc_smooth(o,3);
    partition_by_plane(o,1,0,0,0);
    center_positioning(o);
    cout << "Executed laplacian HC smooth " << clock() - t << "ms" << endl;

    t = clock();
    cout << "Writing..." << endl;
    o.save(out);
    cout << "Writed " << clock() - t << "ms" << endl;
    in.close();
    out.close();
    return 0;
}
