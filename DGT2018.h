#ifndef DGT2018_H
#define DGT2018_H

#include <chrono>
#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Delaunay_triangulation_2<K>  DelaunayTriangulation;
typedef DelaunayTriangulation::Vertex_handle Vertex_handle;

class DGT2018 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;

public:
    DGT2018(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {
        auto start = std::chrono::high_resolution_clock::now();

        DelaunayTriangulation T;
        bool isEmptyTriangulation = true;

        for( Point p : P ) {
            if(isEmptyTriangulation) {
                isEmptyTriangulation = false;
                diskCenters.push_back(p);
                T.insert(p);
            } else {
                Vertex_handle handleToTheNearestDiskCenter = T.nearest_vertex(p);
                if(squared_distance(p,handleToTheNearestDiskCenter->point()) > 1) {
                    diskCenters.push_back(p);
                    T.insert(p);
                }
            }
        }
        T.clear();

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // DGT2018_H
