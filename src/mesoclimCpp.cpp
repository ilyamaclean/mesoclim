#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for delineating basins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// Written 19th Feb 2024 by Ilya Maclean
// Identify lowest pixel in dm2 where dun = 0
IntegerVector whichmin(NumericMatrix& dm2, IntegerMatrix& dun) {
    int rows = dm2.nrow(); // Get number of rows of dm2
    int cols = dm2.ncol(); // Get number of columns of dm2
    int rw = -1; // Initialize row index
    int cl = -1; // Initialize column index
    double d = dm2(0, 0); // Initialize d with the first element of dm2
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (dm2(i, j) < d && dun(i, j) < 1) {
                rw = i; // Update row index
                cl = j; // Update column index
                d = dm2(i, j); // Update minimum value
            }
        }
    }
    IntegerVector out(2);
    out[0] = rw;
    out[1] = cl;
    return out;
}
// Identify lowest pixel within same basin in dm2 where dun = 0
IntegerVector whichmin2(NumericMatrix& dm2, IntegerMatrix& b, IntegerMatrix& dun, int bn) {
    int rows = dm2.nrow(); // Get number of rows of dm2
    int cols = dm2.ncol(); // Get number of columns of dm2
    int rw = -1; // Initialize row index
    int cl = -1; // Initialize column index
    double d = dm2(0, 0); // Initialize d with the first element of dm2
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (dm2(i, j) < d && b(i, j) == bn && dun(i, j) == 0) {
                rw = i; // Update row index
                cl = j; // Update column index
                d = dm2(i, j); // Update minimum value
            }
        }
    }
    IntegerVector out(2);
    out[0] = rw;
    out[1] = cl;
    return out;
}
// Select 3 x 3 matrix surrounding target focal cell given by s (numeric)
NumericMatrix sel3n(NumericMatrix& dm2, IntegerVector s) {
    NumericMatrix m3(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m3(i, j) = dm2(s(0) - 1 + i, s(1) - 1 + j);
        }
    }
    return m3;
}
// Select 3 x 3 matrix surrounding target focal cell given by s (integer)
IntegerMatrix sel3i(IntegerMatrix& dm2, IntegerVector s) {
    IntegerMatrix m3(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m3(i, j) = dm2(s(0) - 1 + i, s(1) - 1 + j);
        }
    }
    return m3;
}
// Assign all grid cells of equal height or higher surrounding a focal grid cell to same basin
IntegerMatrix assignhigher(NumericMatrix& m3, IntegerMatrix& b3, int bn) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (m3(i, j) >= m3(1, 1) && b3(i, j) != 0 && m3(i, j) != 9999) {
                b3(i, j) = bn;
            }
        }
    }
    return b3;
}
// slot b3 into bsn based on s
IntegerMatrix slotin(IntegerMatrix& b, IntegerMatrix b3, IntegerVector s) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            b(s[0] + i - 1, s[1] + j - 1) = b3(i, j);
        }
    }
    return b;
}
// Function that does the basin delineation
// 'basinCpp
// @export
// [[Rcpp::export]]
IntegerMatrix basinCpp(NumericMatrix& dm2, IntegerMatrix& bsn, IntegerMatrix& dun) {
    int bn = 1;
    int tsta = 1;
    while (tsta == 1) {
        // Initial iteration
        IntegerVector s = whichmin(dm2, dun); // lowest undone pixel
        if (s[1] > -1) {
            // assign basin and dun
            bsn(s[0], s[1]) = bn;
            dun(s[0], s[1]) = 1;
            // select 3 x 3 matrix around dm2 and bsn
            NumericMatrix m3 = sel3n(dm2, s);
            IntegerMatrix b3 = sel3i(bsn, s);
            // identify which grid cells in 3 x 3 undone and higher and assign to basin bn
            b3 = assignhigher(m3, b3, bn);
            // Slot in b3 into basin
            bsn = slotin(bsn, b3, s);
        }
        else {
            tsta = 0;
        }
        // Subsequent iteration
        int tst = 1;
        while (tst == 1) {
            s = whichmin2(dm2, bsn, dun, bn); // lowest undone pixel
            if (s[1] > -1) {
                // assign basin and dun
                bsn(s[0], s[1]) = bn;
                dun(s[0], s[1]) = 1;
                // select 3 x 3 matrix around dm2 and bsn
                NumericMatrix m3 = sel3n(dm2, s);
                IntegerMatrix b3 = sel3i(bsn, s);
                // identify which grid cells in 3 x 3 undone and higher and assign to basin bn
                b3 = assignhigher(m3, b3, bn);
                // Slot in b3 into basin
                bsn = slotin(bsn, b3, s);
            }
            else {
                tst = 0;
            }
        }
        bn++;
    }
    return bsn;
}
// Function used to renumber basins sequentially
// [[Rcpp::export]]
IntegerVector renumberbasin(IntegerVector& m, IntegerVector u) {
    for (int i = 0; i < u.size(); ++i) {
        for (int j = 0; j < m.size(); ++j) {
            if (m[j] == u[i]) {
                m[j] = i + 1;
            }
        }
    }
    return m;
}
// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~ Function used for calculating coastal exposure ~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// Based on code written by James P Duffy
// ' invls_calc
//'
//' This function returns a numeric matrix identifying
//' matching values from a raster that intersect with
//' a set of points
//' @export
// [[Rcpp::export]]
NumericMatrix invls_calc(NumericMatrix lsm, double resolution,
    double xmin, double ymax, NumericVector s,
    int direction, NumericMatrix slr,
    double slr_xmin, double slr_xmax,
    double slr_ymin, double slr_ymax) {

    double pi = 3.141593;
    int lsm_row = lsm.nrow();
    int lsm_col = lsm.ncol();

    //set up collection matrix
    NumericMatrix lsw(lsm_row, lsm_col);

    //for each cell in the land/sea mask
    for (int yy = 0; yy < lsm_row; ++yy) {
        for (int xx = 0; xx < lsm_col; ++xx) {

            //if it is not a sea pixel...
            if (lsm(yy, xx) != 0) {

                double x = (xx + 1) * resolution + xmin - resolution / 2;
                double y = ymax + resolution / 2 - (yy + 1) * resolution;

                NumericVector xdist = round(s * sin(direction * pi / 180), 0);
                NumericVector ydist = round(s * cos(direction * pi / 180), 0);

                //update existing
                xdist = xdist + x;
                ydist = ydist + y;

                //create matrix of xy coords (points)
                NumericMatrix xy(xdist.size(), 2);
                xy(_, 0) = xdist;
                xy(_, 1) = ydist;

                //collector
                NumericVector lsc(xy.nrow());

                //for every point (in steps away from focal cell)
                for (int point = 0; point < xy.nrow(); point++) {
                    double x2 = xy(point, 0);
                    double y2 = xy(point, 1);

                    //check position of point falls within raster
                    if (x2 > slr_xmin && slr_xmax > x2 && y2 > slr_ymin && slr_ymax > y2) {

                        //work out index of position
                        double rdist = slr_ymax - y2;
                        int row = floor(rdist / resolution);
                        double cdist = x2 - slr_xmin;
                        int col = floor(cdist / resolution);

                        //pull value from slr that matches the point
                        lsc(point) = slr(row, col);
                    }
                    else {
                        //if no point matches add NA to collector
                        lsc(point) = NA_REAL;
                    }
                }

                lsc(0) = 1;
                lsc = na_omit(lsc);

                if (lsc.size() > 0) {

                    lsw(yy, xx) = sum(lsc) / lsc.size();

                }
                else {

                    lsw(yy, xx) = 1;

                }

            }
            else {
                lsw(yy, xx) = NA_REAL;
            }
        }
    }
    return lsw;
}
// ============================================================================================= #
// ~~~~~~~~~~~~~~~~~~~~~~~~~ Functions used for adjusting rainfall ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================================= #
// [[Rcpp::export]]
std::vector<double> rainadjustv(std::vector<double> rain, std::vector<double> rrain, double rfrac, double rtot)
{
    // Calculate expected rianfall days/hours
    if (!std::isnan(rain[0])) {
        double pr = rfrac * rain.size();
        int prd = std::round(pr);
        int ard = 0;
        // Calculate actual rainfall days/hours
        for (size_t i = 0; i < rain.size(); ++i) if (rain[i] > 0) ard++;
        // if actual rain days/hours < predicted rain days/hours give some random no rain days a small amount of rain
        if (ard < prd) {
            // Calculate number of zero days/hours to assign rain to
            int rta = prd - ard;
            // Calculate minimum none-zero in rain
            double nz = 5001;
            for (size_t i = 0; i < rain.size(); ++i) if (rain[i] > 0 && rain[i] < nz) nz = rain[i];
            // Create a vector of indices corresponding to zeros in rain
            std::vector<int> zeroIndices;
            for (int i = 0; i < rain.size(); ++i) {
                if (rain[i] == 0) {
                    zeroIndices.push_back(i);
                }
            }
            // Sort zeroIndices based on corresponding values in rrain
            std::sort(zeroIndices.begin(), zeroIndices.end(), [&rrain](int i, int j) { return rrain[i] < rrain[j]; });
            // Assign the lowest rta indices a value nz
            for (int i = 0; i < rta && i < zeroIndices.size(); ++i) {
                rain[zeroIndices[i]] = nz; // or any desired non-zero value
            }
        }
        // if actual rain days > predicted rain days give lowest rainfall days zero rain
        if (ard > prd) {
            // Find number of days/hours that should be zero
            int nrd = rain.size() - prd;
            // Create a vector of indices ordered by the values in rain
            std::vector<size_t> indices(rain.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&rain](size_t i, size_t j) { return rain[i] < rain[j]; });
            // Assign 0 to the lowest nrd elements
            for (int i = 0; i < nrd && i < rain.size(); ++i) {
                rain[indices[i]] = 0;
            }
        }
        // adjust rainfall to match total
        double rsum = rain[0];
        for (size_t i = 1; i < rain.size(); ++i) rsum = rsum + rain[i];
        double mu = rtot / rsum;
        for (size_t i = 0; i < rain.size(); ++i) rain[i] = rain[i] * mu;
    }
    return rain;
}
std::vector<std::vector<double>> convertoCppmatrix(NumericMatrix mat) {
    std::vector<std::vector<double>> result(mat.nrow(), std::vector<double>(mat.ncol()));
    for (int i = 0; i < mat.nrow(); ++i) {
        for (int j = 0; j < mat.ncol(); ++j) {
            result[i][j] = mat(i, j);
        }
    }
    return result;
}
NumericMatrix convertoRmatrix(std::vector<std::vector<double>>& mat) {
    int nrow = mat.size();
    int ncol = mat[0].size(); // Assuming all inner vectors have the same size
    NumericMatrix mat2(nrow, ncol);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            mat2(i, j) = mat[i][j];
        }
    }
    return mat2;
}
// [[Rcpp::export]]
NumericMatrix rainadjustm(NumericMatrix rainm, std::vector<double> rrain, std::vector<double> rfrac, std::vector<double> rtot)
{
    std::vector<std::vector<double>> rainc = convertoCppmatrix(rainm);
    for (size_t i = 0; i < rainc.size(); ++i) {
        rainc[i] = rainadjustv(rainc[i], rrain, rfrac[i], rtot[i]);
    }
    rainm = convertoRmatrix(rainc);
    return rainm;
}
