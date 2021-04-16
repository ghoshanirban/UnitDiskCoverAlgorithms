#ifndef CCFM1997_H
#define CCFM1997_H

#include <chrono>
#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Delaunay_triangulation_2<K>  DelaunayTriangulation;
typedef DelaunayTriangulation::Vertex_handle Vertex_handle;

class CCFM1997 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;

    long double sqrt3 = std::sqrt(3), sqrt3Over2 = sqrt3/2;

    inline void newPoint(const Point &p, DelaunayTriangulation &activeDiskCenters, DelaunayTriangulation &inactiveDiskCenters) {
        activeDiskCenters.insert(p);
        inactiveDiskCenters.insert(Point(p.x()+sqrt3,       p.y()));
        inactiveDiskCenters.insert(Point(p.x()+sqrt3Over2,  p.y()+1.5));
        inactiveDiskCenters.insert(Point(p.x()+sqrt3Over2,  p.y()-1.5));
        inactiveDiskCenters.insert(Point(p.x()-sqrt3Over2,  p.y()+1.5));
        inactiveDiskCenters.insert(Point(p.x()-sqrt3,       p.y()));
        inactiveDiskCenters.insert(Point(p.x()-sqrt3Over2,  p.y()-1.5));
    }

public:
    CCFM1997(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {

        assert(P.size() > 0);
        auto start = std::chrono::high_resolution_clock::now();

        DelaunayTriangulation activeDiskCenters, inactiveDiskCenters;
        newPoint(P[0],activeDiskCenters,inactiveDiskCenters);

        for(unsigned i = 1; i <  P.size(); i++) {

            Vertex_handle handleToTheNearestDiskCenter = activeDiskCenters.nearest_vertex(P[i]);
            if(!(squared_distance(P[i],handleToTheNearestDiskCenter->point()) > 1))
                continue;

            if(inactiveDiskCenters.number_of_vertices() == 0) {
                newPoint(P[i],activeDiskCenters,inactiveDiskCenters);
                continue;
            }

            handleToTheNearestDiskCenter = inactiveDiskCenters.nearest_vertex(P[i]);

            if( squared_distance(P[i],handleToTheNearestDiskCenter->point()) > 1)
                newPoint(P[i],activeDiskCenters,inactiveDiskCenters);
            else {
                inactiveDiskCenters.remove(handleToTheNearestDiskCenter);
                activeDiskCenters.insert(handleToTheNearestDiskCenter->point());
            }
        }

        for (auto it = activeDiskCenters.finite_vertices_begin(); it != activeDiskCenters.finite_vertices_end(); it++)
            diskCenters.push_back(it->point());

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // CCFM1997_H



