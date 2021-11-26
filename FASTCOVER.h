//
// Created by Ghosh, Anirban on 11/11/21.
//

#ifndef UDC_FASTCOVER_H
#define UDC_FASTCOVER_H

#include <chrono>
#include <list>
#include <vector>
#include <unordered_set>

#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<double>::Point_2 Point;

class FASTCOVER {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;

    typedef std::pair<int,int> intPair;
    typedef std::unordered_set<intPair,boost::hash<intPair>> SetOfCells;
    const double sqrt2 = std::sqrt(2);
    const double additiveFactor = sqrt2/2;

    public:
    FASTCOVER(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {
        assert(!P.empty());

        auto start = std::chrono::high_resolution_clock::now();

        SetOfCells S;
        for(const Point &p : P)
            S.insert(std::make_pair(floor(p.x()/sqrt2),floor(p.y()/sqrt2)));

        for(auto pair : S)
            diskCenters.emplace_back(Point(pair.first*sqrt2+additiveFactor,pair.second*sqrt2+additiveFactor));

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif //UDC_FASTCOVER_H
