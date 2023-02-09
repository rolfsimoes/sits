#include <Rcpp.h>
#include <cmath>
#include <map>
#include <vector>
#include <set>
using namespace Rcpp;

// ---- type definitions ----

struct Regions {
    Regions(NumericMatrix data, int nrow, int ncol) :
        nrow_(nrow), ncol_(ncol), size_(nrow * ncol), data_(data) {
        if (size_ != data.nrow()) {
            stop("Invalid matrix dimentions");
        }
        for (int i = 0; i < size_; i++) {
            area_.insert({i, {i}});
            neigh_.insert({i, neighbors(i)});
        }
    }

    std::set<int> neighbors(int index) {
        std::set<int> neigh;
        int dx[] = { 0, 1, 0, -1};
        int dy[] = {-1, 0, 1,  0};
        for (int i = 0; i < 4; i++) {
            int x = index % ncol_ + dx[i];
            int y = index / ncol_ + dy[i];
            if (x >= 0 && y >= 0 && x < ncol_ && y < nrow_) {
                neigh.insert(y * ncol_ + x);
            };
        }
        return(neigh);
    }

    void add_neigh(int r1, int r2) {
        neigh_[r1].insert(r2);
        neigh_[r2].insert(r1);
    }

    void del_neigh(int r1, int r2) {
        neigh_[r1].erase(r2);
        neigh_[r2].erase(r1);
    }

    void merge_regions(int r1, int r2) {
        int n1 = neigh_[r1].size(), n2 = neigh_[r2].size();
        data_.row(r1) = (data_.row(r1) * n1 + data_.row(r2) * n2) / (n1 + n2);
        del_neigh(r1, r2);
        for (int neighbor : neigh_[r2]) {
            neigh_[neighbor].erase(r2); // remove r2 from neighbors
            add_neigh(r1, neighbor);
        }
        neigh_[r2].clear(); // remove neighbors from r2
        area_[r1].insert(area_[r2].begin(), area_[r2].end());
        area_[r2].clear();
    }

    void flat() {
        for (auto r1 : area_) {
            for (int r2 : r1.second) {
                data_.row(r2) = data_.row(r1.first);
            }
        }
    }

    int best_neigh(int r1, double threshold) {
        int r_min = -1;
        double d_min = -1;
        for (int r2 : neigh_[r1]) {
            int d = sqrt(sum(pow(data_.row(r1) - data_.row(r2), 2)));
            if (std::isnan(d)) continue;
            if (d > threshold) continue;
            if (d < d_min || d_min < 0) {
                d_min = d;
                r_min = r2;
            }
        }
        return(r_min);
    }

    std::vector<int> merge_best_neighbors(double threshold) {
        std::vector<int> merged;
        for (auto region : neigh_) {
            if (!region.second.size()) continue;
            int r1 = region.first;
            int r2 = best_neigh(r1, threshold);
            // no best neighbor
            if (r2 < 0) continue;
            // not reciprocal
            int r1_min = best_neigh(r2, threshold);
            if (r1_min != r1) continue;
            // merge them
            merge_regions(r1, r2);
            merged.push_back(r2);
        }
        return(merged);
    };

    int nrow_, ncol_, size_;
    NumericMatrix data_;
    std::map<int, std::set<int>> area_;
    std::map<int, std::set<int>> neigh_;
};

// ---- accessors ----

inline int _nrow(List& r) {
    return(as<int>(r["nrow"]));
}

inline int _ncol(List& r) {
    return(as<int>(r["ncol"]));
}

inline int _size(List& r) {
    return(as<int>(r["size"]));
}

inline NumericMatrix _avg(List& r) {
    return(as<NumericMatrix>(r["avg"]));
}

inline IntegerVector _iter(List& r) {
    return(as<IntegerVector>(r["iter"]));
}


// [[Rcpp::export]]
List reg_setup(int nrow, int ncol, NumericMatrix x) {
    int size = nrow * ncol;
    if (x.nrow() != size)
        stop("Invalid matrix dimentions");
    IntegerVector iter = {0};
    List r = List::create(
        Named("nrow") = nrow,
        Named("ncol") = ncol,
        Named("size") = size,
        Named("avg") = x,
        Named("iter") = iter
    );
    return(r);
}

// [[Rcpp::export]]
void segment_region_growing(List& r, double a, double h0, int dt) {
    if (a <= 0.0 || h0 <= 0.0 || dt < 1) return;
    int t, t0 = _iter(r)(0);
    Regions reg(_avg(r), _nrow(r), _ncol(r));
    // update
    for (t = t0; t < t0 + dt; t++) {
        std::vector<int> merged = reg.merge_best_neighbors(pow(a, t) * h0);
        if (!merged.size()) break;
        for (auto m : merged) {
            reg.neigh_.erase(m);
        }
    }
    reg.flat();
    _avg(r) = reg.data_;
    _iter(r)(0) = t;
}

