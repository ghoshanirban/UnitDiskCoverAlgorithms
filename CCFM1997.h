//
// Created by Ghosh, Anirban on 11/9/21.
//

#ifndef CCFM1997_H
#define CCFM1997_H

#include <chrono>
#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/spatial_sort.h>

#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<double>::Point_2 Point;

class CCFM1997 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;
    long double sqrt3 = std::sqrt(3), sqrt3Over2 = sqrt3/2;
    typedef std::pair<int,int> intPair;
    typedef std::unordered_map< intPair, std::list<Point>, boost::hash<intPair> > HashMap;
    const double sqrt2 = std::sqrt(2);

    inline void insertIntoMap(HashMap &H, const Point &p) const {
        int v = floor(p.x()/sqrt2), h = floor(p.y()/sqrt2);
        auto it = H.find(std::make_pair(v,h));

        if( it != H.end() )
            it->second.emplace_back(p);
        else {
            std::list<Point> L;
            L.push_back(p);
            H[std::make_pair(v, h)] = L;
        }
    }

    inline void removeFromMap(HashMap &H, const Point &p) const {
        int v = floor(p.x()/sqrt2), h = floor(p.y()/sqrt2);
        auto it = H.find(std::make_pair(v,h));

        if( it != H.end() )
            return;
        else {
           auto L = it->second;
           for(auto itForL = L.begin(); itForL != L.end(); ++itForL)
               if( *itForL == p ) {
                   L.erase(itForL);
                   return;
               }
        }
    }

    inline void static nnInsideACell(const std::list<Point> &L, const Point &p, double &minDistance, Point &nn) {
        for(const Point &q : L) {
            double dist = CGAL::squared_distance(p,q);
            if( dist < minDistance ) {
                minDistance = dist;
                nn = q;
            }
        }
    }

    inline Point findNN(HashMap &H, const Point &p) const {
        double minDistance = INFINITY;
        Point nn(INFINITY,INFINITY);
        int v = floor(p.x()/sqrt2), h = floor(p.y()/sqrt2);

        auto it = H.find(std::make_pair(v,h));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v-1,h+1));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v,h+1));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v+1,h+1));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v-1,h));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v+1,h));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v-1,h-1));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v,h-1));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        it = H.find(std::make_pair(v+1,h-1));
        if( it != H.end() && !it->second.empty() )
            nnInsideACell(it->second,p,minDistance,nn);

        return nn;
    }

    inline void newPoint(const Point &p, HashMap &activeDiskCenters, HashMap &inactiveDiskCenters)  {
        insertIntoMap(activeDiskCenters,p);
        insertIntoMap(inactiveDiskCenters,Point(p.x()+sqrt3,       p.y()));
        insertIntoMap(inactiveDiskCenters,Point(p.x()+sqrt3Over2,  p.y()+1.5));
        insertIntoMap(inactiveDiskCenters,Point(p.x()+sqrt3Over2,  p.y()-1.5));
        insertIntoMap(inactiveDiskCenters,Point(p.x()-sqrt3Over2,  p.y()+1.5));
        insertIntoMap(inactiveDiskCenters,Point(p.x()-sqrt3,       p.y()));
        insertIntoMap(inactiveDiskCenters,Point(p.x()-sqrt3Over2,  p.y()-1.5));
    }

public:
    CCFM1997(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {

        assert(!P.empty());

        auto start = std::chrono::high_resolution_clock::now();

        HashMap activeDiskCenters, inactiveDiskCenters;
        newPoint(P[0],activeDiskCenters,inactiveDiskCenters);
        diskCenters.emplace_back(P[0]);

        for(unsigned i = 1; i <  P.size(); i++) {
            Point nearestDiskCenter = findNN(activeDiskCenters,P[i]);

            if(nearestDiskCenter != Point(INFINITY,INFINITY) && squared_distance(P[i],nearestDiskCenter) <=1)
                continue;

            if( inactiveDiskCenters.empty() ) {
                newPoint(P[i],activeDiskCenters,inactiveDiskCenters);
                diskCenters.emplace_back(P[i]);
                continue;
            }

            nearestDiskCenter = findNN(inactiveDiskCenters,P[i]);

            if( nearestDiskCenter == Point(INFINITY,INFINITY) || squared_distance(P[i],nearestDiskCenter) > 1) {
                newPoint(P[i], activeDiskCenters, inactiveDiskCenters);
                diskCenters.emplace_back(P[i]);
            }
            else {
                removeFromMap(inactiveDiskCenters,nearestDiskCenter);
                insertIntoMap(activeDiskCenters,nearestDiskCenter);
                diskCenters.emplace_back(nearestDiskCenter);
            }
        }

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // CCFM1997_H

