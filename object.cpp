#include "object.h"
#include "algorithm.h"

#include <ctime>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    clock_t t;
    Object o;

    t = clock();
    cout << "Loading..." << endl;
    o.load(in);
    cout << "Loaded " << clock - t << "ms" << endl;

    t = clock();
    cout << "Executing laplacian HC smooth"
    laplacian_hc_smooth(o,3);
    cout << "Executed laplacian HC smooth " << clock - t << "ms" << endl;
    
    t = clock();
    cout << "Writing..."
    o.save(out);
    cout << "Writed " << clock() - t << "ms" << endl;
    in.close();
    out.close();
    return 0;
}