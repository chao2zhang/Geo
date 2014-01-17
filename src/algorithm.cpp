#include "algorithm.h"

#include <vector>

using namespace std;

void laplacian_hc_smooth(Object& obj, int times, float alpha, float beta) {
    vector<Point3f> o(obj.vertex);
    vector<Point3f> p(o);
    for (int i = 0; i < times; i++) {
        vector<Point3f> b(o.size());
        vector<Point3f> q(p);
        for (int i = 0; i < o.size(); i++) {
            if (obj.adjacentList[i].size()) {
                Point3f s(0, 0, 0);
                for (int j = 0; j < obj.adjacentList[i].size(); j++)
                    s += q[obj.adjacentList[i][j]];
                p[i] = s / obj.adjacentList[i].size();
            }
            b[i] = p[i] - (o[i] * alpha + q[i] * (1 - alpha));
        }
        for (int i = 0; i < o.size(); i++) {
            if (obj.adjacentList[i].size()) {
                Point3f s(0, 0, 0);
                for (int j = 0; j < obj.adjacentList[i].size(); j++)
                    s += b[obj.adjacentList[i][j]];
                p[i] = p[i] - (b[i] * beta + s * (1 - beta) / obj.adjacentList[i].size());
            }
        }
    }
    obj.vertex.assign(p.begin(), p.end());
    obj.update();
}
