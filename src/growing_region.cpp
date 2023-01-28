#include <Rcpp.h>
#include <cmath>
#include <forward_list>
using namespace Rcpp;

//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

inline int _nrow(List r) {
    return(as<int>(r["nrow"]));
};

inline int _ncol(List r) {
    return(as<int>(r["ncol"]));
};

inline int _size(List r) {
    return(as<int>(r["size"]));
}

inline int _px_i(List r, int px) {
    return(px % _nrow(r));
};

inline int _px_j(List r, int px) {
    return(px / _nrow(r));
};

inline int _ij_px(List r, int i, int j) {
    if (i < 0 || i >= _nrow(r)) return(-1);
    if (j < 0 || j >= _ncol(r)) return(-1);
    return(j * _nrow(r) + i);
};

inline NumericMatrix _avg(List r) {
    return(as<NumericMatrix>(r["avg"]));
};

inline IntegerMatrix _reg(List r) {
    return(as<IntegerMatrix>(r["reg"]));
};

inline IntegerMatrix _count(List r) {
    return(as<IntegerMatrix>(r["count"]));
};

// [[Rcpp::export]]
List reg_setup(int nrow, int ncol, NumericMatrix& avg) {
    int size = nrow * ncol;
    if (avg.nrow() != size)
        stop("Invalid matrix dimentions");
    IntegerMatrix reg(size, 1);
    IntegerMatrix count(size, 1);
    for (int i = 0; i < size; i++) {
        reg(i, 0) = i;
        count(i, 0) = 1;
    }
    List r = List::create(
        Named("nrow") = nrow,
        Named("ncol") = ncol,
        Named("size") = size,
        Named("reg") = reg,
        Named("count") = count,
        Named("avg") = avg
    );
    return(r);
};

int _get_reg(List r, int px) {
    int r0 = _reg(r)(px, 0);
    if (r0 < 0 || r0 >= _size(r)) return(-1); // not valid region
    if (r0 == px) return(px);
    return(_get_reg(r, r0));
};

void _set_reg(List r, int px0, int px1) {
    int size = _size(r);
    int r0 = _reg(r)(px0, 0), r1 = _get_reg(r, px1);
    if (r0 < 0 || r0 >= size) return; // not valid region
    if (r1 < 0 || r1 >= size) return; // not valid region
    if (r0 != px0) _set_reg(r, r0, r1);
    // update reg, count, and avg
    _reg(r)(px0, 0) = r1;
    _count(r)(px0, 0) = _count(r)(r1, 0);
    _avg(r).row(px0) = _avg(r).row(r1);
};

inline int _von_neumann(List r, int px, int ni, int nj) {
    int i = _px_i(r, px), j = _px_j(r, px);
    if ((ni != 0 && nj != 0) || (ni == 0 && nj == 0)) return(-1);
    return(_ij_px(r, i + ni, j + nj));
}

int best_neigh(List r, int px, double h) {
    int size = _size(r);
    _set_reg(r, px, px); // update px
    int r0 = _reg(r)(px, 0), px_min = -1;
    if (r0 < 0 || r0 >= size) return(-1); // not valid region
    double d_min = -1.0;
    NumericVector avg0 = _avg(r).row(px), avg1;
    for (int ni = -1; ni <= 1; ni++) {
        for (int nj = -1; nj <= 1; nj++) {
            // von-neumann neighborhood
            int px1 = _von_neumann(r, px, ni, nj);
            if (px1 < 0) continue; // non neighbor or outside raster
            // update neighbor;
            _set_reg(r, px1, px1);
            int r1 = _reg(r)(px1, 0);
            if (r1 < 0 || r1 >= size) continue; // not valid region
            if (r0 == r1) continue; // in the same region
            double d = sum(pow(_avg(r).row(px1) - avg0, 2));
            if (std::isnan(d)) continue; // invalid distance
            // get minimum
            if (d_min < 0 || d_min > d) {
                d_min = d;
                px_min = px1;
            }
        }
    }
    if (d_min > h) return(-1);
    return(px_min);
};

void merge(List r, int px0, int px1) {
    int c0 = _count(r)(px0, 0), c1 = _count(r)(px1, 0), c2 = c0 + c1;
    _avg(r).row(px0) = (_avg(r).row(px0) * c0 + _avg(r).row(px1) * c1) / c2;
    _count(r)(px0, 0) = c2;
    _set_reg(r, px1, px0);
};

// [[Rcpp::export]]
bool merge_best_neigh(List r, int px, double h) {
    int size = _size(r);
    if (px < 0 || px >= size) return(false);
    if (h <= 0) return(false);
    int px1 = best_neigh(r, px, h);
    if (px1 < 0) return(false); // no best neighbor
    merge(r, px, px1);
    return(true);
};

// [[Rcpp::export]]
void segment_region_growing(List r, double a, double h0) {
    int size = _size(r);
    for (int t = 0; t < 100; t++) {
        bool merged = false;
        double h = pow(a, t) * h0;
        for (int px = 0; px < size; px++) {
            merged |= merge_best_neigh(r, px, h);
        }
        if (!merged) break;
    }
};
