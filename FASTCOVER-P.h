//
// Created by Ghosh, Anirban on 11/11/21.
//

#ifndef UDC_FASTCOVER_P_H
#define UDC_FASTCOVER_P_H

#include <chrono>
#include <list>
#include <vector>
#include <unordered_set>

#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<double>::Point_2 Point;

class FASTCOVER_P {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;

    typedef std::pair<int,int> intPair;
    typedef std::unordered_set<intPair,boost::hash<intPair>> SetOfCells;
    const double sqrt2 = std::sqrt(2);
    const double additiveFactor = sqrt2/2;
    const double sqrt2TimesOnePointFiveMinusOne = (sqrt2 * 1.5)-1;
    const double sqrt2TimesZeroPointFivePlusOne = (sqrt2 * 0.5)-1;

public:
    FASTCOVER_P(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {
        assert(!P.empty());

        auto start = std::chrono::high_resolution_clock::now();

        SetOfCells S;
        for(const Point &p : P) {
            int v = floor(p.x()/sqrt2), h = floor(p.y()/sqrt2);
            double verticalTimesSqrtTwo = v*sqrt2, horizontalTimesSqrt2 = h*sqrt2;

            auto it = S.find(std::make_pair(v,h));
            if( it != S.end() )
                continue;

            if( (p.x() >= verticalTimesSqrtTwo + sqrt2TimesOnePointFiveMinusOne) ) {
                it = S.find(std::make_pair(v+1,h));
                if( it != S.end() && (CGAL::squared_distance(p,Point(sqrt2*(v+1)+additiveFactor,horizontalTimesSqrt2+additiveFactor))<=1))
                    continue;
            }

            if( (p.x() <= verticalTimesSqrtTwo - sqrt2TimesZeroPointFivePlusOne) ) {
                it = S.find(std::make_pair(v-1,h));
                if( it != S.end() && (CGAL::squared_distance(p,Point(sqrt2*(v-1)+additiveFactor,horizontalTimesSqrt2+additiveFactor))<=1))
                    continue;
            }

            if( (p.y() <= horizontalTimesSqrt2 + sqrt2TimesOnePointFiveMinusOne) ) {
                it = S.find(std::make_pair(v,h-1));
                if( it != S.end() && (CGAL::squared_distance(p,Point(verticalTimesSqrtTwo+additiveFactor,sqrt2*(h-1)+additiveFactor))<=1))
                    continue;
            }

            if( (p.y() >= horizontalTimesSqrt2 - sqrt2TimesZeroPointFivePlusOne) ) {
                it = S.find(std::make_pair(v,h+1));
                if(it != S.end() && (CGAL::squared_distance(p,Point(verticalTimesSqrtTwo+additiveFactor,sqrt2*(h+1)+additiveFactor))<=1))
                    continue;
            }
            S.insert(std::make_pair(v,h));
        }

        for(auto pair : S)
            diskCenters.emplace_back(Point(pair.first*sqrt2+additiveFactor,pair.second*sqrt2+additiveFactor));

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif //UDC_FASTCOVER_P_H
