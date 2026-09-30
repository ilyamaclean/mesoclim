// Coastal exposure: land fraction in the upwind direction.
//
// Ported from the terravars package (I. Maclean, src/wind.cpp @ 68d99cf,
// https://github.com/ilyamaclean/terravars), which generalises the earlier
// invls_calc() (microclima / mesoclim) to multiple resolution levels.
// Keep in step with terravars if the method changes there.
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <cmath>
#include <cstdint>

// jitter_offset
//
// Deterministic pseudo-random value in [-1, 1], derived purely from
// integer bit-mixing of (row, col, point) -- no external state, no RNG
// calls, so it's safe to call from any thread inside the OpenMP parallel
// region below (R's own RNG is not thread-safe and must never be called
// from inside a parallel region, per the pattern used throughout this
// package). Deterministic in the inputs means the same (row, col, point)
// always yields the same offset, so a call with jitter enabled is fully
// reproducible for identical inputs, even though there's no user-facing
// seed to set.
// ------------------------------------------------------------
static inline double jitter_offset(int row, int col, int point) {
    uint32_t h = static_cast<uint32_t>(row)   * 374761393u
                + static_cast<uint32_t>(col)   * 668265263u
                + static_cast<uint32_t>(point) * 2654435761u;
    h = (h ^ (h >> 13)) * 1274126177u;
    h = h ^ (h >> 16);
    return (static_cast<double>(h) / 4294967295.0) * 2.0 - 1.0;
}

// coastal_exposure_cpp
//
// For each land cell in `lsm` (the cropped target-region land/sea
// matrix), samples land/sea values at a set of distances `s` along
// azimuth `direction`, using whichever supplied resolution level is
// finest among those that actually cover each sampled point (levels are
// tried in the order supplied -- the R wrapper sorts them finest-first),
// falling back to coarser levels only when a point falls outside a
// finer level's extent. A point outside every level's extent is skipped.
// The mean of all successfully-sampled 0/1 (sea/land) values is this
// cell's coefficient.
//
// This is a port of the reference implementation `invls_calc()`,
// generalised from a single fixed-resolution land/sea raster to a list
// of levels at different resolutions. The "inverse distance" weighting
// happens entirely via the *spacing* of `s`, built by the R wrapper with
// a tunable power -- there is no explicit distance-based weight applied
// here, matching the original.
//
// Parameters:
//   lsm        - binary (0/1) land/sea matrix for the CROPPED target
//                region only; determines which cells get a value
//                computed (sea cells -> NA) and the output's shape
//   resolution - resolution (metres) of `lsm` -- i.e. the finest level's
//                resolution
//   xmin, ymax - origin of the TARGET region (used to place each output
//                cell's own coordinate; not the same as any level's own
//                xmin/ymax, which describe where to look things up)
//   s          - sample distances (metres) along `direction`, including
//                a leading 0 for the focal cell itself
//   direction  - azimuth (degrees, N=0, clockwise) searched directly --
//                i.e. the direction the wind is blowing FROM, with no
//                +180 adjustment (matching horizon_angle_cpp()'s own
//                "search toward this azimuth" convention)
//   masks      - list of binary (0/1) land/sea matrices, one per
//                resolution level, ORDERED FINEST TO COARSEST
//   mask_reso, mask_xmin, mask_xmax, mask_ymin, mask_ymax
//              - parallel vectors (one entry per `masks` element) giving
//                each level's own resolution and extent
//   jitter_deg - 0 (default): sample exactly along `direction`, as
//                before. > 0: each sample point's azimuth is perturbed
//                independently by a deterministic pseudo-random amount
//                in [-jitter_deg, jitter_deg] (see jitter_offset() above)
//                before its (x2, y2) location is computed, fanning the
//                transect into a narrow cone instead of an exact line.
//
// No Rcpp object allocation happens inside the parallel loop below
// (only plain double/int arithmetic and reads from pre-built plain
// vectors/matrices) -- matches the pattern used elsewhere in this
// package for OpenMP-parallel loops, since R's allocator isn't safe to
// call concurrently from multiple threads.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix coastal_exposure_cpp(Rcpp::NumericMatrix lsm, double resolution,
    double xmin, double ymax, Rcpp::NumericVector s, double direction,
    Rcpp::List masks, Rcpp::NumericVector mask_reso,
    Rcpp::NumericVector mask_xmin, Rcpp::NumericVector mask_xmax,
    Rcpp::NumericVector mask_ymin, Rcpp::NumericVector mask_ymax,
    double jitter_deg = 0.0) {

    const double pi = 3.141593;
    int lsm_row  = lsm.nrow();
    int lsm_col  = lsm.ncol();
    int n_levels = masks.size();
    int n_s      = s.size();

    // Cache each level's matrix (and dimensions) once, up front -- all
    // Rcpp/R-API work happens here, before the parallel region.
    std::vector<Rcpp::NumericMatrix> level_mats(n_levels);
    std::vector<int> level_nrow(n_levels), level_ncol(n_levels);
    for (int lev = 0; lev < n_levels; ++lev) {
        level_mats[lev] = Rcpp::as<Rcpp::NumericMatrix>(masks[lev]);
        level_nrow[lev] = level_mats[lev].nrow();
        level_ncol[lev] = level_mats[lev].ncol();
    }

    Rcpp::NumericMatrix lsw(lsm_row, lsm_col);

    double sin_d = std::sin(direction * pi / 180.0);
    double cos_d = std::cos(direction * pi / 180.0);

    const double* s_ptr         = s.begin();
    const double* mask_reso_ptr = mask_reso.begin();
    const double* mask_xmin_ptr = mask_xmin.begin();
    const double* mask_xmax_ptr = mask_xmax.begin();
    const double* mask_ymin_ptr = mask_ymin.begin();
    const double* mask_ymax_ptr = mask_ymax.begin();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int yy = 0; yy < lsm_row; ++yy) {
        for (int xx = 0; xx < lsm_col; ++xx) {
            if (lsm(yy, xx) != 0) {
                double x = (xx + 1) * resolution + xmin - resolution / 2;
                double y = ymax + resolution / 2 - (yy + 1) * resolution;

                double sum_val = 0.0;
                int    count   = 0;

                for (int point = 0; point < n_s; ++point) {
                    double val;
                    if (point == 0) {
                        val = 1.0;   // the focal cell itself is land
                    } else {
                        double sp = sin_d, cp = cos_d;
                        if (jitter_deg > 0.0) {
                            double off = jitter_offset(yy, xx, point) * jitter_deg;
                            double ang = (direction + off) * pi / 180.0;
                            sp = std::sin(ang);
                            cp = std::cos(ang);
                        }
                        double x2 = x + std::round(s_ptr[point] * sp);
                        double y2 = y + std::round(s_ptr[point] * cp);

                        val = NA_REAL;
                        for (int lev = 0; lev < n_levels; ++lev) {
                            if (x2 > mask_xmin_ptr[lev] && mask_xmax_ptr[lev] > x2 &&
                                y2 > mask_ymin_ptr[lev] && mask_ymax_ptr[lev] > y2) {
                                int row = static_cast<int>(
                                    std::floor((mask_ymax_ptr[lev] - y2) / mask_reso_ptr[lev]));
                                int col = static_cast<int>(
                                    std::floor((x2 - mask_xmin_ptr[lev]) / mask_reso_ptr[lev]));
                                if (row >= 0 && row < level_nrow[lev] &&
                                    col >= 0 && col < level_ncol[lev]) {
                                    val = level_mats[lev](row, col);
                                    break;   // finest covering level wins
                                }
                            }
                        }
                    }
                    if (!std::isnan(val)) {
                        sum_val += val;
                        ++count;
                    }
                }
                lsw(yy, xx) = (count > 0) ? (sum_val / count) : 1.0;
            } else {
                lsw(yy, xx) = NA_REAL;
            }
        }
    }
    return lsw;
}
