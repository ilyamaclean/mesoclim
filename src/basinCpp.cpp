#include <Rcpp.h>
using namespace Rcpp;

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