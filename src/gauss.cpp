#include "gauss.h"

#include <cmath>

using std::vector;

bool gauss(const vector<vector<float> >& matrix, const vector<float>& value, vector<float>& ans) {
    int n = value.size();
    if (matrix.size() != n)
        return false;
    for (vector<float>& v : matrix)
        if (v.size() != n)
            return false;
    if (ans.size() != n)
        return false;
    vector<vector<float> > m(matrix);
    for (int i = 0; i < n; ++i) {
        m[i].push_back(value[i]);
    }
    for (int i = 0; i < n; ++i) {
        int maxi = 0;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(m[j][i]) > fabs(m[maxi][i]))
                maxi = j;
        }
        if (fabs(m[maxi][i]) < eps)
            return false;
        m[maxi].swap(m[i]);
        float v = m[i][i];
        for (int j = i; j < n; ++j)
            m[i][j] /= v;
        for (int j = i + 1; j < n; ++j) {
            v = m[j][i];
            for (int k = i; k <= n; ++k)
                m[j][k] -= v * m[i][k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        ans[i] = m[i][n];
        for (int j = i + 1; j < n; ++j) {
            ans[i] -= ans[j] * m[i][j];
        }
    }
    return true;
}
