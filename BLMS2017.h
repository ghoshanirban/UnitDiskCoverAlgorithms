#ifndef BLMS2017_H
#define BLMS2017_H

#include <chrono>
#include <list>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

class BLMS2017 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;

    struct BSTcomp {
        bool operator () ( const Point &lhs, const Point &rhs ) const {
            return lhs.y() < rhs.y()
                   || ( lhs.y() == rhs.y() && lhs.x() < rhs.x() );
        }
    };

    struct Event {
        Point point;
        Point deletion = this->point;

        bool operator < ( const Event b ) const {
            return point.x() < b.point.x()
                   || ( point.x() == b.point.x() && point.y() < b.point.y() );
        }
    };

    bool isCovered( Point p, Point diskCenter ) {
        return CGAL::squared_distance( p, diskCenter ) <= 4;
    }

public:
    BLMS2017(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {
        assert(P.size() > 0);

        auto start = std::chrono::high_resolution_clock::now();

        const double sqrt3 = std::sqrt(3), sqrt3Over2 = sqrt3/2;
        std::vector<Event> E;
        std::set<Point,BSTcomp> BST;

        for( Point p : P ) {
            E.push_back( { p } );
            E.push_back( { Point( p.x()+2, p.y() ), p } );
        }

        std::sort( E.begin(), E.end() );

        double x = E.at(0).point.x();
        double y = E.at(0).point.y();

        diskCenters.push_back( Point( x, y ) );
        diskCenters.push_back( Point( x + sqrt3, y ) );
        diskCenters.push_back( Point( x + sqrt3Over2, y + 1.5 ) );
        diskCenters.push_back( Point( x + sqrt3Over2, y - 1.5 ) );

        BST.insert( E.at(0).point );
        E.erase( E.begin() );

        for( Event e : E ) {

            if( e.deletion == e.point ) {

                std::set<Point,BSTcomp>::iterator iter = BST.lower_bound( e.point );
                if( iter != BST.end() ) {
                    if( isCovered( e.point, *(iter) ) ) {
                        continue;
                    }
                    if( ++iter != BST.end() ) {
                        if( isCovered( e.point, *(iter) ) )  {
                            continue;
                        }
                    }
                    --iter;
                }
                if( iter != BST.begin() ) {
                    --iter;
                    if( isCovered( e.point, *(iter) ) ) {
                        continue;
                    }
                    if( iter != BST.begin() ) {
                        --iter;
                        if( isCovered( e.point, *(iter) ) ) {
                            continue;
                        }
                    }
                }

                double x = e.point.x();
                double y = e.point.y();

                diskCenters.push_back( Point( x, y ) );
                diskCenters.push_back( Point( x + sqrt3, y ) );
                diskCenters.push_back( Point( x + sqrt3Over2, y + 1.5 ) );
                diskCenters.push_back( Point( x + sqrt3Over2, y - 1.5 ) );

                BST.insert( e.point );

            } else {
                BST.erase( e.deletion );
            }
        }

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // BLMS2017_H


