#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix basinCpp(NumericMatrix dm2, IntegerMatrix bsn, IntegerMatrix dun) {
  int bn = 1;
  int tsta = 1;
  
  int rows = dm2.nrow(); // Get number of rows of dm2
  int cols = dm2.ncol(); // Get number of columns of dm2

  while (tsta == 1) {
    // Initial iteration
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
    if (rw > -1) {
      // assign basin and dun
      bsn(rw, cl) = bn;
      dun(rw, cl) = 1;
      
      // identify which grid cells in 3 x 3 undone and higher and assign to basin bn
      for (int i = rw - 1; i <= rw + 1; ++i) {
        for (int j = cl - 1; j <= cl + 1; ++j) {
          if (i >= 0 && i < rows && j >= 0 && j < cols && dm2(i, j) >= dm2(rw, cl) && bsn(i, j) != 0 && dm2(i, j) != 9999) { 
            bsn(i, j) = bn;
          }
        }
      }
    } else {
      tsta = 0;
    }

    // Subsequent iteration
    int tst = 1;
    while (tst == 1) {
      rw = -1;
      cl = -1;
      d = dm2(0, 0);
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
          if (dm2(i, j) < d && bsn(i, j) == bn && dun(i, j) == 0) {
            rw = i;
            cl = j;
            d = dm2(i, j);
          }
        }
      }
      if (rw > -1) {
        bsn(rw, cl) = bn;
        dun(rw, cl) = 1;
        for (int i = rw - 1; i <= rw + 1; ++i) {
          for (int j = cl - 1; j <= cl + 1; ++j) {
            if (i >= 0 && i < rows && j >= 0 && j < cols && dm2(i, j) >= dm2(rw, cl) && bsn(i, j) != 0 && dm2(i, j) != 9999) { 
              bsn(i, j) = bn;
            }
          }
        }
      } else {
        tst = 0;
      }
    }
    bn++;
  }
  return bsn;
}