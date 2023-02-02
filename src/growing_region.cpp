#include <Rcpp.h>
#include <cmath>
#include <map>
using namespace Rcpp;

// ---- type definitions ----

typedef std::tuple<int, int, double> neigh_dist; // neigh, seed, dist
template<class T>
using neighborhood = std::map<int, T>;
template<class T>
inline int _px(T& neigh) {
    return(std::get<0>(neigh));
};
template<class T>
inline int _seed(T& neigh) {
    return(std::get<1>(neigh));
};
template<class T>
inline int _dist(T& neigh) {
    return(std::get<2>(neigh));
};
template<class T>
inline bool _no_neigh(T& neigh) {
    return(_px(neigh) < 0);
};

// ---- accessors ----

inline int _nrow(List& r) {
    return(as<int>(r["nrow"]));
};

inline int _ncol(List& r) {
    return(as<int>(r["ncol"]));
};

inline int _size(List& r) {
    return(as<int>(r["size"]));
}

inline NumericMatrix _avg(List& r) {
    return(as<NumericMatrix>(r["avg"]));
};

inline IntegerMatrix _reg(List& r) {
    return(as<IntegerMatrix>(r["reg"]));
};

inline IntegerMatrix _count(List& r) {
    return(as<IntegerMatrix>(r["count"]));
};

inline IntegerVector _iter(List& r) {
    return(as<IntegerVector>(r["iter"]));
}

inline neigh_dist _get_dist(neighborhood<neigh_dist>& neigh, int r0) {
    neighborhood<neigh_dist>::iterator n_it = neigh.find(r0);
    if (n_it == neigh.end()) return(neigh_dist(-1, -1, -1.0));
    return(n_it->second);
}

int _get_seed(List& r, int px) {
    int r0 = _reg(r)(px, 0);
    if (r0 < 0) return(-1); // not valid region
    if (r0 == px) return(r0);
    return(_get_seed(r, r0));
};

void _set_reg(List& r, int px, int reg) {
    int r0 = _get_seed(r, px);
    if (r0 < 0) return; // not a valid region
    // update reg, count, and avg
    _reg(r)(r0, 0) = reg;
    _count(r)(r0, 0) = _count(r)(reg, 0);
    _avg(r).row(r0) = _avg(r).row(reg);
};

// ---- utils ----

void _flat_reg(List& r) {
    int size = _size(r);
    for (int px = 0; px < size; px++) {
        int r0 = _get_seed(r, px);
        if (r0 < 0) continue; // invalid region
        _reg(r)(px, 0) = r0;
        _count(r)(px, 0) = _count(r)(r0, 0);
        _avg(r).row(px) = _avg(r).row(r0);
    }
}

inline int _px_i(List& r, int px) {
    return(px / _ncol(r));
};

inline int _px_j(List& r, int px) {
    return(px % _ncol(r));
};

inline int _ij_px(List& r, int i, int j) {
    return(i * _ncol(r) + j);
};

inline int _von_neumann(List& r, int px, int ni, int nj) {
    int i = _px_i(r, px), j = _px_j(r, px);
    if ((ni != 0 && nj != 0) || (ni == 0 && nj == 0)) return(-1);
    if (i + ni < 0 || j + nj < 0 ||
        i + ni >= _nrow(r) || j + nj >= _ncol(r)) return(-1);
    return(_ij_px(r, i + ni, j + nj));
}

// ---- growing region (best neighbor) ----

// find best distance neighbor
neigh_dist _best_neigh(List& r, int px, double h) {
    // starting negative => to be updated
    int size = _size(r);
    double d_min = -1.0;
    int r0 = _get_seed(r, px), r_min = -1, px_min = -1;
    if (r0 < 0) return(neigh_dist(-1, -1, -1.0)); // not a valid region
    NumericVector avg0 = _avg(r).row(r0);
    for (int ni = -1; ni <= 1; ni++) {
        for (int nj = -1; nj <= 1; nj++) {
            // von-neumann neighborhood
            int px1 = _von_neumann(r, px, ni, nj);
            if (px1 < 0) continue; // outside raster
            int r1 = _get_seed(r, px1);
            if (r1 < 0) continue; // not a valid region
            if (r0 == r1) continue; // they are in the same region
            double d = sqrt(sum(pow(_avg(r).row(r1) - avg0, 2)));
            if (std::isnan(d)) continue; // invalid distance
            if (d > h) continue; // above the minimum threshold
            // get minimum distance above or equal zero
            if (r_min < 0 || d_min > d) {
                px_min = px1;
                r_min = r1; // seed
                d_min = d;
            }
        }
    }
    return(neigh_dist(px_min, r_min, d_min));
};

// find best distance neighbors
neighborhood<neigh_dist> _find_best_neighbors(List& r, double h) {
    int size = _size(r);
    neighborhood<neigh_dist> neigh;
    for (int px = 0; px < size; px++) {
        int r0 = _get_seed(r, px);
        if (r0 < 0) continue; // not a valid region
        neigh_dist r1_min = _best_neigh(r, px, h); // reference neighbor
        // no best neighbor
        if (_no_neigh(r1_min)) continue;
        // not the best neighbor
        neigh_dist r0_neigh = _get_dist(neigh, r0);
        if (_dist(r0_neigh) >= 0 &&
            _dist(r0_neigh) <= _dist(r1_min)) continue;
        // add reference neighbor
        neigh[r0] = r1_min;
    }
    return(neigh);
};

bool _neigh_reciprocal(neighborhood<neigh_dist>& neigh, int r0, int r1) {
    neigh_dist r0_neigh = _get_dist(neigh, r0);
    return(!_no_neigh(r0_neigh) && r1 == _seed(r0_neigh));
}

inline void _merge(List& r, int r0, int r1) {
    int c0 = _count(r)(r0, 0), c1 = _count(r)(r1, 0), c2 = c0 + c1;
    _avg(r).row(r0) = (_avg(r).row(r0) * c0 + _avg(r).row(r1) * c1) / c2;
    _count(r)(r0, 0) = c2;
    _set_reg(r, r1, r0);
}

// merge best neighboring regions
bool _merge_best_neighbors(List& r, neighborhood<neigh_dist>& neigh) {
    bool merged = false;
    // check if best neighbor is reciprocal and merge it
    for (auto n_it = neigh.begin(); n_it != neigh.end(); n_it++) {
        // check reciprocal best neighbor
        int r0 = _get_seed(r, n_it->first);
        int r1 = _get_seed(r, _seed(n_it->second));
        if (r0 < 0 || r1 < 0) continue; // invalid neigh
        // not reciprocal neighbor
        if (!_neigh_reciprocal(neigh, r0, r1)) continue;
        // merge them
        _merge(r, r0, r1);
        merged = true;
    }
    return(merged);
};

// [[Rcpp::export]]
List reg_setup(int nrow, int ncol, NumericMatrix x) {
    int size = nrow * ncol;
    if (x.nrow() != size)
        stop("Invalid matrix dimentions");
    IntegerMatrix reg(size, 1);
    IntegerMatrix count(size, 1);
    for (int i = 0; i < size; i++) {
        reg(i, 0) = i;
        count(i, 0) = 1;
    }
    IntegerVector iter = {0};
    List r = List::create(
        Named("nrow") = nrow,
        Named("ncol") = ncol,
        Named("size") = size,
        Named("reg") = reg,
        Named("count") = count,
        Named("avg") = x,
        Named("iter") = iter
    );
    return(r);
};

// [[Rcpp::export]]
void segment_region_growing(List& r, double a, double h0, int dt) {
    if (a <= 0.0 || h0 <= 0.0 || dt < 1) return;
    int t, t0 = _iter(r)(0);
    // update
    for (t = t0; t < t0 + dt; t++) {
        double h = pow(a, t) * h0;
        neighborhood<neigh_dist> neigh_best = _find_best_neighbors(r, h);
        bool merged = _merge_best_neighbors(r, neigh_best);
        if (!merged) break;
    }
    _flat_reg(r);
    _iter(r)(0) = t;
};


// find and merge best neighboring pixels
bool _merge_best_pixels(List& r, double h, int t) {
    int size = _size(r);
    bool merged = false;
    for (int px = 0; px < size; px++) {
        int r0 = _get_seed(r, px);
        if (r0 < 0) continue; // not a valid region
        neigh_dist r1_min = _best_neigh(r, px, h);
        // no best neighbor
        if (_no_neigh(r1_min)) continue;
        // not the best neighbor
        neigh_dist r0_min = _best_neigh(r, _px(r1_min), h);
        // no best neighbor
        if (_no_neigh(r0_min)) continue;
        // not reciprocal
        if (_px(r0_min) != px) continue;
        // merge them
        int r1 = _seed(r1_min);
        _merge(r, r0, r1);
        merged = true;
    }
    return(merged);
};

// [[Rcpp::export]]
void segment_region_growing_pix(List& r, double a, double h0, int dt) {
    if (a <= 0.0 || h0 <= 0.0 || dt < 1) return;
    int t, t0 = _iter(r)(0);
    // update
    for (t = t0; t < t0 + dt; t++) {
        double h = pow(a, t) * h0;
        bool merged = _merge_best_pixels(r, h, t);
        if (!merged) break;
    }
    _flat_reg(r);
    _iter(r)(0) = t;
};

