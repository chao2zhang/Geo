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

void center_positioning(Object& o) {
    Point3f u(o.vertex[0]);
    Point3f d(o.vertex[0]);
    for (int i = 0; i < o.vertex.size(); i++) {
        if (o.vertex[0].x[0] > u.x[0]) u.x[0] = o.vertex[0].x[0];
        if (o.vertex[0].x[0] < d.x[0]) d.x[0] = o.vertex[0].x[0];
        if (o.vertex[0].x[1] > u.x[1]) u.x[1] = o.vertex[0].x[1];
        if (o.vertex[0].x[1] < d.x[1]) d.x[1] = o.vertex[0].x[1];
        if (o.vertex[0].x[2] > u.x[2]) u.x[2] = o.vertex[0].x[2];
        if (o.vertex[0].x[2] < d.x[2]) d.x[2] = o.vertex[0].x[2];
    }
    Point3f mid = (u + d) / 2;
    for (int i = 0; i < o.vertex.size(); i++)
        o.vertex[i] -= mid;
}
