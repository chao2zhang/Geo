#include "algorithm.h"
#include "object.h"

#include <vector>

using namespace std;

void laplacian_hc_smooth(Object& obj, int times=1, float alpha=0, float beta=0.5) {
    vector<Point3f> o = obj.vertex;
    vector<Point3f> p = o;
    for (int i = 0; i < times; i++) {
        vector<Point3f> b;
        vector<Point3f> q = p;
        for (int i = 0; i < o.size(); i++) {
            if (obj.adjacentList[i].size()) {
                Point3f s;
                for (j = 0; j < obj.adjacentList[i].size(); j++)
                    s += q[obj.adjacentList[i][j];
                p[i]  = s / obj.adjacentList[i].size();
            }
            b[i] = p[i] - (alpha * o[i] + (1 - alpha) * q[i]);
        }
        for (int i = 0; i < o.size(); i++) {
            if (obj.adjacentList[i].size()) {
                Point3f s;
                for (j = 0; j < obj.adjacentList[i].size(); j++)
                    s += b[obj.adjacentList[i][j];
                p[i] = p[i] - (beta * b[i] + (1 - beta) * s / obj.adjacentList[i].size());
            }
        }
    }
    obj.vertex = p;
    obj.update();
}