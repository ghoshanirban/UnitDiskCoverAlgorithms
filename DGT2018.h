#ifndef DGT2018_H
#define DGT2018_H

#include <chrono>
#include <list>
#include <vector>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/spatial_sort.h>

#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<double>::Point_2 Point;

class DGT2018 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;
    typedef std::pair<int,int> intPair;
    typedef std::unordered_map< intPair, std::list<Point>, boost::hash<intPair> > HashMap;
    const double sqrt2 = std::sqrt(2);

    bool inline static isCoveredBy(const Point &p, const HashMap &H, const HashMap::iterator it) {
        if( it != H.end() )
            for(const Point &q : it->second)
                if(CGAL::squared_distance(p,q) <= 1)
                    return true;
        return false;
    }

    void inline checkIfCoveredElseInsertCenter(HashMap &H, const Point &p, std::list<Point> &C) const {
        int v = floor(p.x()/sqrt2), h = floor(p.y()/sqrt2);
        auto it = H.find(std::make_pair(v,h));

        if(isCoveredBy(p,H,it)) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v,  h+1)))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v-1,h  )))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v+1,h  )))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v,  h-1)))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v-1,h-1)))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v-1,h+1)))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v+1,h-1)))) return;
        if(isCoveredBy(p,H,H.find(std::make_pair(v+1,h+1)))) return;

        if( it == H.end() ) {
            std::list<Point> L;
            L.push_back(p);
            H[std::make_pair(v,h)] = L;
        }
        else
            it->second.emplace_back(p);

        C.emplace_back(p);
    }

public:
    DGT2018(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {

        assert(!P.empty());

        auto start = std::chrono::high_resolution_clock::now();
        HashMap H;

        for(const Point &p : P)
            checkIfCoveredElseInsertCenter(H, p, diskCenters);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // DGT2018_H
