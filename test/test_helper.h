#ifndef TEST_HELPER_H
#define TEST_HELPER_H

#include "ProjectDefs.h"
#include <vector>

ArrayType<Real, 2> make_mat(int m, int n, const std::vector<Real>& vals);

ArrayType<Real, 1> make_vec(const std::vector<Real>& vals);


#endif
