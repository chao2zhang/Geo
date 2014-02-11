#ifndef GAUSS_H
#define GAUSS_H

#include <vector>
#include "base.h"

bool gauss(const std::vector<std::vector<float> >& matrix,
           const std::vector<float>& value,
           std::vector<float>& ans);

#endif // GAUSS_H
