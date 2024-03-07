#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;
// ============================================================================================ = #
// ~~~~~~~~~~~~~~~~~~~ Functions used for converting between R and C++ matrices ~~~~~~~~~~~~~~~~~ #
// ============================================================================================== #
// ** Convert R matrix to C++ matrix ** //
std::vector<std::vector<double>> convertoCppmatrix(NumericMatrix mat) {
    std::vector<std::vector<double>> result(mat.nrow(), std::vector<double>(mat.ncol()));
    for (int i = 0; i < mat.nrow(); ++i) {
        for (int j = 0; j < mat.ncol(); ++j) {
            result[i][j] = mat(i, j);
        }
    }
    return result;
}
// ** Convert C++ matrix to R matrix ** //
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
// ** Applies a specified function to an hourly 3D array to return daily data ** //
// 'hourtodayCpp
// @export
// [[Rcpp::export]]
NumericVector hourtodayCpp(NumericVector a, std::string fun) {
    IntegerVector dims = a.attr("dim");
    int dim1 = dims[0];
    int dim2 = dims[1];
    int dim3 = dims[2];
    int ndays = dim3 / 24; // Number of days
    NumericVector daily(dim1 * dim2 * ndays, NA_REAL);
    // Calculate daily mean, max, or min
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            double x = a[j * dim1 + i];
            if (!std::isnan(x)) {
                for (int day = 0; day < ndays; ++day) {
                    double out = 0;
                    if (fun == "mean" || fun == "sum") {
                        for (int hr = 0; hr < 24; ++hr) {
                            double v = a[hr * dim1 * dim2 + day * dim1 * dim2 * 24 + j * dim1 + i];
                            out = out + v;
                        }
                    }
                    else if (fun == "max") {
                        out = a[day * dim1 * dim2 * 24 + j * dim1 + i];
                        for (int hr = 1; hr < 24; ++hr) {
                            double v = a[hr * dim1 * dim2 + day * dim1 * dim2 * 24 + j * dim1 + i];
                            if (v > out) out = v;
                        }
                    }
                    else if (fun == "min") {
                        out = a[day * dim1 * dim2 * 24 + j * dim1 + i];
                        for (int hr = 1; hr < 24; ++hr) {
                            double v = a[hr * dim1 * dim2 + day * dim1 * dim2 * 24 + j * dim1 + i];
                            if (v < out) out = v;
                        }
                    }
                    if (fun == "mean") out = out / 24;
                    daily[day * dim1 * dim2 + j * dim1 + i] = out;
                }
            }
        }
    }
    // Reshape the 1D vector to a 3D array
    daily.attr("dim") = IntegerVector::create(dim1, dim2, ndays);
    return daily;
}
// ============================================================================================ = #
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions for calculating solar variables  ~~~~~~~~~~~~~~~~~ #
// ============================================================================================== #
int juldayCpp(int year, int month, int day)
{
    double dd = day + 0.5;
    int madj = month + (month < 3) * 12;
    int yadj = year + (month < 3) * -1;
    double j = std::trunc(365.25 * (yadj + 4716)) + std::trunc(30.6001 * (madj + 1)) + dd - 1524.5;
    int b = 2 - std::trunc(yadj / 100) + std::trunc(std::trunc(yadj / 100) / 4);
    int jd = static_cast<int>(j + (j > 2299160) * b);
    return jd;
}
// 'juldayvCpp
// @export
// [[Rcpp::export]]
IntegerVector juldayvCpp(IntegerVector year, IntegerVector month, IntegerVector day)
{
    IntegerVector jd = clone(day);
    for (int i = 0; i < jd.size(); ++i) {
        jd[i] = juldayCpp(year[i], month[i], day[i]);
    }
    return jd;
}
// ** Calculates solar time ** //
double soltimeCpp(int jd, double lt, double lon)
{
    double m = 6.24004077 + 0.01720197 * (jd - 2451545);
    double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
    double st = lt + (4 * lon + eot) / 60;
    return st;
}
// ** Calculates day length ** //
double daylengthCpp(int jd, double lat) 
{
    double declin = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double latr = lat * M_PI / 180;
    double hc = -0.01453808 / (cos(latr) * cos(declin)) - tan(latr) * tan(declin);
    double dl = 0;
    if (hc < -1) {
        dl = 24;
    }
    else if (hc < 1) {
        double ha = (acos(hc)) * 180 / M_PI;
        double m = 6.24004077 + 0.01720197 * (jd - 2451545);
        double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
        double sr = (720 - 4 * ha - eot) / 60;
        double ss = (720 + 4 * ha - eot) / 60;
        dl = ss - sr;
    }
    return(dl);
}
// ** Calculates solar zenith in radians ** //
double solzenCpp(double jd, double st, double lat, double lon)
{
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double coh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double z = acos(coh);
    return z;
}
// ** Calculates solar azimuth in radians ** //
double solaziCpp(double jd, double st, double lat, double lon)
{
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double sh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double hh = atan(sh / sqrt(1 - sh * sh));
    double sazi = cos(dec) * sin(tt) / cos(hh);
    double cazi = (sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec)) /
        sqrt(pow(cos(dec) * sin(tt), 2) + pow(sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec), 2));
    double sqt = 1 - sazi * sazi;
    if (sqt < 0) sqt = 0;
    double azi = M_PI + atan(sazi / sqrt(sqt));
    if (cazi < 0) {
        if (sazi < 0) {
            azi = M_PI - azi;
        }
        else {
            azi = 3 * M_PI - azi;
        }
    }
    return azi;
}
// ** Calculates clear sky radiation ** //
double clearskyradCpp(int jd, double lt, double lat, double lon, double tc = 15.0, double rh = 80.0, double pk = 101.3)
{
    double st = soltimeCpp(jd, lt, lon);
    double z = solzenCpp(jd, st, lat, lon);
    double Ic = 0.0;
    if (z <= M_PI / 2) {
        double m = 35 * cos(z) * pow(1224 * cos(z) * cos(z) + 1, -0.5);
        double TrTpg = 1.021 - 0.084 * sqrt(m * 0.00949 * pk + 0.051);
        double xx = log(rh / 100) + ((17.27 * tc) / (237.3 + tc));
        double Td = (237.3 * xx) / (17.27 - xx);
        double u = exp(0.1133 - log(3.78) + 0.0393 * Td);
        double Tw = 1 - 0.077 * pow(u * m, 0.3);
        double Ta = 0.935 * m;
        double od = TrTpg * Tw * Ta;
        Ic = 1352.778 * cos(z) * od;
    }
    return Ic;
}
// ** Calculates diffuse fraction ** //
double difpropCpp(double swrad, int jd, double lt, double lat, double lon)
{
    double d = 1.0;
    if (swrad > 0) {
        double st = soltimeCpp(jd, lt, lon);
        double z = solzenCpp(jd, st, lat, lon);
        if (z < M_PI / 2) {
            double zd = z * 180 / M_PI;
            double k1 = 0.83 - 0.56 * exp(-0.06 * (90 - zd));
            double si = 0.0;
            if (z <= M_PI / 2) si = cos(z);
            double k = 0.0;
            if (si > 0) k = swrad / (1352 * si);
            if (k > k1) k = k1;
            if (k < 0) k = 0;
            double rho = k / k1;
            double sigma3 = 0;
            if (rho > 1.04) {
                sigma3 = 0.12 + 0.65 * (rho - 1.04);
            }
            else {
                sigma3 = 0.021 + 0.397 * rho - 0.231 * pow(rho, 2) - 0.13 * exp(-1 * pow((rho - 0.931) / 0.134, 2) * 0.834);
            }
            double k2 = 0.95 * k1;
            double d1 = 1.0;
            if (zd < 88.6) d1 = 0.07 + 0.046 * zd / (93 - zd);
            double K = 0.5 * (1 + sin(M_PI * (k - 0.22) / (k1 - 0.22) - M_PI / 2));
            double d2 = 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K * K));
            double d3 = (d2 * k2) * (1 - k) / (k * (1 - k2));
            double alpha = pow(1 / cos(z), 0.6);
            double kbmax = pow(0.81, alpha);
            double kmax = (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2));
            double dmax = (d2 * k2) * (1 - kmax) / (kmax * (1 - k2));
            d = 1 - kmax * (1 - dmax) / k;
            if (k <= kmax) d = d3;
            if (k <= k2) d = d2;
            if (k <= 0.22) d = 1;
            double kX = 0.56 - 0.32 * exp(-0.06 * (90 - zd));
            double kL = (k - 0.14) / (kX - 0.14);
            double kR = (k - kX) / 0.71;
            double delta = (k >= 0.14 && k < kX) ? (-3 * pow(kL, 2) * (1 - kL) * pow(sigma3, 1.3)) : 0;
            if (k >= kX && k < (kX + 0.71)) delta = 3 * kR * pow((1 - kR), 2) * pow(sigma3, 0.6);
            if (sigma3 > 0.01) d = d + delta;
        }
    }
    return d;
}
// ** Calculates diffuse fraction on a vector of data ** //
std::vector<double> difpropvCpp(std::vector<double> swrad, std::vector<int> jd, std::vector<double> lt, double lat, double lon)
{
    std::vector<double> d(lt.size());
    for (size_t i = 0; i < d.size(); ++i) {
        d[i] = difpropCpp(swrad[i], jd[i], lt[i], lat, lon);
    }
    return d;
}
// ** Calculates clear sky radiation for an hourly vector ** //
std::vector<double> clearskyradhCpp(std::vector<int>& jd, std::vector<double>& lt, double lat, double lon)
{
    std::vector<double> Ic(lt.size());
    for (size_t i = 0; i < Ic.size(); ++i) {
        Ic[i] = clearskyradCpp(jd[i], lt[i], lat, lon, 15.0, 80.0, 101.3);
    }
    return Ic;
}
// ** Calculates clear sky radiation for a daily vector ** //
std::vector<double> clearskyraddCpp(std::vector<int>& jd, double lat)
{
    std::vector<double> Ic(jd.size());
    for (size_t i = 0; i < Ic.size(); ++i) {
        Ic[i] = 0;
        for (int j = 0; j < 1440; ++j) {
            double lt = j / 60;
            Ic[i] = Ic[i] + clearskyradCpp(jd[i], lt, lat, 0.0, 15.0, 80.0, 101.3);
        }
        Ic[i] = Ic[i] / 1440;
    }
    return Ic;
}
// Calculate clear sky radiation with matrices ** //
// 'clearskyradmCpp
// @export
// [[Rcpp::export]]
NumericMatrix clearskyradmCpp(std::vector<int> jd, std::vector<double> lt, std::vector<double> lat, std::vector<double> lon,
    bool hourly = true) {
    int rows = lat.size();
    int cols = jd.size();
    std::vector<std::vector<double>> Ic(rows, std::vector<double>(cols, std::nan("")));
    for (size_t i = 0; i < lat.size(); ++i) {
        if (!std::isnan(lat[i])) {
            if (hourly) {
                Ic[i] = clearskyradhCpp(jd, lt, lat[i], lon[i]);
            }
            else {
                Ic[i] = clearskyraddCpp(jd, lat[i]);
            }
        }
    }
    NumericMatrix ICm = convertoRmatrix(Ic);
    return ICm;
}
// ** Calculates diffuse fraction with matrices ** //
// 'difpropmCpp
// @export
// [[Rcpp::export]]
NumericMatrix difpropmCpp(NumericMatrix swrad, std::vector<int> jd, std::vector<double> lt, std::vector<double> lat, std::vector<double> lon)
{
    std::vector<std::vector<double>> swradc = convertoCppmatrix(swrad);
    int rows = lat.size();
    int cols = jd.size();
    std::vector<std::vector<double>> d(rows, std::vector<double>(cols, std::nan("")));
    for (size_t i = 0; i < lat.size(); ++i) {
        if (!std::isnan(swradc[i][1])) {
            d[i] = difpropvCpp(swradc[i], jd, lt, lat[i], lon[i]);
        }
    }
    NumericMatrix dr = convertoRmatrix(d);
    return dr;
}
// ** Calculates solar zenith for vector ** //
std::vector<double> solzenvCpp(std::vector<int> jd, std::vector<double> lt, double lat, double lon)
{
    std::vector<double> z(lt.size());
    for (size_t i = 0; i < z.size(); ++i) {
        double st = soltimeCpp(jd[i], lt[i], lon);
        z[i] = solzenCpp(jd[i], st, lat, lon);
        z[i] = z[i] * 180 / M_PI;
    }
    return z;
}
// ** Calculates solar zenith for matrices ** //
// 'solzenmCpp
// @export
// [[Rcpp::export]]
NumericMatrix solzenmCpp(std::vector<int> jd, std::vector<double> lt, std::vector<double> lat, std::vector<double> lon)
{
    int nrow = lat.size();
    int ncol = jd.size();
    std::vector<std::vector<double>> z(nrow, std::vector<double>(ncol, std::nan("")));
    for (size_t i = 0; i < lat.size(); ++i) {
        if (!std::isnan(lat[i])) {
            z[i] = solzenvCpp(jd, lt, lat[i], lon[i]);
        }
    }
    NumericMatrix zd = convertoRmatrix(z);
    return(zd);
}
// ** Calculates 24 hourly temperatures for daily data
std::vector<double> tempintdayCpp(double tmn, double tmnn, double tmx, double dl, double stt, double lat, double lon, double srte = 0.09)
{
    // Calculate predicted night fraction
    double ngtp = 0.04187957 * ((tmx - tmn) * (1 - dl / 24)) + 0.4372056;
    if (ngtp < 0.01) ngtp = 0.01;
    if (ngtp > 0.99) ngtp = 0.99;
    // Calculate sunrise time
    double sr = 12 - 0.5 * dl;
    // Calculate solar time after sunrise
    std::vector<double> thour(24);
    for (int lt = 0; lt < 24; lt++) {
        double st = lt + stt - sr; // solar time after sunrise
        if (st < 0) st = st + 24;
        if (st > 24) st = st - 24;
        if (dl == 24) {
            thour[lt] = (tmx - tmn) * sin((M_PI * st) / 28) + tmn;
            double gr = (tmnn - tmn) * (st / 24);
            thour[lt] = thour[lt] + gr;
        }
        else if (dl == 0) {
            st = lt + stt; // solar time after sunrise
            if (st < 0) st = st + 24;
            if (st > 24) st = st - 24;
            thour[lt] = (tmx - tmn) * sin((M_PI * st) / 24) + tmn;
            double gr = (tmnn - tmn) * (st / 24);
            thour[lt] = thour[lt] + gr;
        }
        else {
            double k = -(24 - dl) / log(srte / ngtp);
            double ph = -0.5 * dl * ((M_PI / (asin(ngtp) - M_PI)) + 1);
            double rho = dl + 2 * ph;
            if (st > dl) {
                thour[lt] = ngtp * exp(-(st - dl) / k);
            }
            else {
                thour[lt] = sin((M_PI * st) / rho);
            }
            // Adjust by tmax and tmin
            thour[lt] = (tmx - tmn) * thour[lt] + tmn;
            // Apply gradient to times after sunset as tmn is next day
            if (lt + stt > dl) {
                double gr = (tmnn - tmn) * (st - dl) / (24 - dl);
                thour[lt] = thour[lt] + gr;    
            }
        }
    }
    // Adjust to ensure tmax and tmin match
    double ptmx = thour[0];
    double ptmn = thour[0];
    for (int lt = 1; lt < 24; lt++) {
        if (thour[lt] > ptmx) ptmx = thour[lt];
        if (thour[lt] < ptmn) ptmn = thour[lt];
    }
    double b = (ptmx - ptmn) / (tmx - tmn);
    double a = b * tmn - ptmn;
    for (int lt = 0; lt < 24; lt++) thour[lt] = (thour[lt] + a) / b;
    return(thour);
}
// ** Interpolates hourly temperature (vector) ** //
// 'hourlytemp
// @export
// [[Rcpp::export]]
std::vector<double> hourlytempv(std::vector<double> tmn, std::vector<double> tmx,
    std::vector<int> year, std::vector<int> month, std::vector<int> day,
    double lat, double lon, double srte = 0.09)
{
    int nrow = tmn.size();
    std::vector<std::vector<double>> thour(nrow, std::vector<double>(24, 0.0));
    for (int i = 0; i < nrow; ++i) {
        // Calculate jd
        int jd = juldayCpp(year[i], month[i], day[i]);
        // Calculate sunrise time
        double dl = daylengthCpp(jd, lat);
        double stt = soltimeCpp(jd, 0, lon);
        double sr = stt + 12 - 0.5 * dl;
        if (sr > 0) { // Sunrise in current day
            double tmnn = tmn[i+1];
            if (i == (nrow - 1)) tmnn = tmn[nrow - 1];
            thour[i] = tempintdayCpp(tmn[i], tmnn, tmx[i], dl, stt, lat, lon, srte);
        }
        else { // sunrise in previous day
            double tmnp = tmn[0];
            if (i > 0) tmnp = tmn[i - 1];
            thour[i] = tempintdayCpp(tmnp, tmn[i], tmx[i], dl, stt, lat, lon, srte);
        }

    }
    std::vector<double> thourv;
    thourv.reserve(nrow * 24);
    for (const auto& row : thour) {
        thourv.insert(thourv.end(), row.begin(), row.end());
    }
    return thourv;
}
// ** Interpolates hourly temperature (matrix) ** //
// 'hourlytemp
// @export
// [[Rcpp::export]]
NumericMatrix hourlytempm(NumericMatrix tmn, NumericMatrix tmx, std::vector<int> year, std::vector<int> month, std::vector<int> day,
    std::vector<int> lat, std::vector<int> lon, double srte = 0.09)
{
    int nrow = lat.size();
    int ncol = year.size() * 24;
    // Convert tmn and tmx to c++ type matrices
    std::vector<std::vector<double>> tmnc = convertoCppmatrix(tmn);
    std::vector<std::vector<double>> tmxc = convertoCppmatrix(tmx);
    // Initalise output
    std::vector<std::vector<double>> thour(nrow, std::vector<double>(ncol, std::nan("")));
    for (int i = 0; i < nrow; ++i) {
        if (!std::isnan(lat[i])) {
            thour[i] = hourlytempv(tmnc[i], tmxc[i], year, month, day, lat[i], lon[i], srte);
        }
    }
    NumericMatrix thourm = convertoRmatrix(thour);
    return(thourm);
}
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
