#include "../../Vendor/cpp/src/linalg.h"
#include "../../GeomFigures.h"

alglib::real_2d_array compute(const alglib::real_2d_array& data, const double& eps, alglib::real_1d_array & centre);
Ellipsoid* KhachiyanAlgo(const vector<Point> & points, const double & eps);