#ifndef FCB2001_H
#define FCB2001_H

#include <chrono>
#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Exact_circular_kernel_2             Circular_k;
typedef CGAL::Point_2<Circular_k>                 PointInCircularKernel;
typedef CGAL::Circle_2<Circular_k>                CircleInCircularKernel;
typedef CGAL::Circular_arc_2<Circular_k>          Circular_arc;
typedef CGAL::CK2_Intersection_traits<Circular_k, CircleInCircularKernel, CircleInCircularKernel>::type Intersection_result;
typedef CGAL::Delaunay_triangulation_2<K>  DelaunayTriangulation;
typedef DelaunayTriangulation::Vertex_handle Vertex_handle;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

class FCB2001 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;
    int l;

    static inline unsigned min(const int a, const int b) {
        return (a < b )? a : b;
    }

    inline void appendPoints(std::list<Point> &dest, std::list<Point> &source) {
        for(Point p: source)
            dest.push_back(p);
    }

    bool isCovered(const std::vector<Point> &P,  const std::list<Point> &selectedCenters) {
        DelaunayTriangulation T(selectedCenters.begin(),selectedCenters.end());

        for(Point p : P) {
            Vertex_handle handleToTheNearestDiskCenter = T.nearest_vertex(p);
            if(squared_distance(p,handleToTheNearestDiskCenter->point()) > 1)//1.0000000001)
                return false;
        }
        return true;
    }

    const double sqrt2 = std::sqrt(2), oneOverSqrt2 = 1/sqrt2, gridSize =  4/(5*sqrt2);

    void findSolInSquare( std::vector<Point> &P, std::list<Point> &optCenters,  const int ell,
                          const double minX, const double maxX, const double minY, const double maxY) {

        std::vector<Point> candidateDiskCenters;
        candidateDiskCenters.reserve(2*P.size()*P.size());

        double xOfLowestLeftmostLatticeDiskCenter = floor(minX/gridSize)*gridSize;
        double yOfLowestLeftmostLatticeDiskCenter = floor(minY/gridSize)*gridSize;

        for( double x = xOfLowestLeftmostLatticeDiskCenter; x <=  maxX; x += gridSize )
            for( double y = yOfLowestLeftmostLatticeDiskCenter; y <=  maxY; y += gridSize )
                candidateDiskCenters.push_back(Point(x,y));

        for(unsigned currentCoverSize = 1; currentCoverSize <= FCB2001::min(2*ell*ell-1,candidateDiskCenters.size()); currentCoverSize++) {

            std::string bitmask(currentCoverSize, '1');
            bitmask.resize(candidateDiskCenters.size(), '0');

            do {
                std::list<Point> selectedCenters;

                for( unsigned pos = 0; pos < bitmask.size(); pos++)
                    if( bitmask[pos] == '1')
                        selectedCenters.push_back(candidateDiskCenters[pos]);

                if( isCovered(P,selectedCenters) ) {
                    appendPoints(optCenters,selectedCenters);
                    return;
                }
            } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
        }


        std::list<Point> tempC;
        xOfLowestLeftmostLatticeDiskCenter = floor(minX/sqrt2)*sqrt2;
        yOfLowestLeftmostLatticeDiskCenter = floor(minY/sqrt2)*sqrt2;
        for( double x = minX; x <= maxX; x += sqrt2 )
            for( double y = minY; y <= maxY; y += sqrt2 )
                tempC.push_back(Point(x + oneOverSqrt2, y + oneOverSqrt2 ));

        for(Point p : tempC)
            optCenters.push_back(p);
    }

public:
    FCB2001(std::vector<Point> &P, std::list<Point> &diskCenters, int ell) : P(P), diskCenters(diskCenters), l(ell) { }

    double execute() {

        assert(l > 0);
        assert(P.size() > 0);

        if(P.size() == 1) {
            diskCenters.push_back(P[0]);
            return 0.0;
        }

        auto start = std::chrono::high_resolution_clock::now();

        long double minX = P[0].x(), maxX = P[0].x();
        for(Point p : P) {
            if( p.x() < minX ) minX = p.x();
            if( p.x() > maxX ) maxX = p.x();
        }

        unsigned numberOfBinsForX = 1 + ceil((maxX - minX) / (2*l));
        unsigned answer = UINT_MAX;

        for(int shiftInX = 0; shiftInX < l; shiftInX++) {

            std::list<Point> tempCForThisShift;
            long double rightOfFirstGroup = minX + (shiftInX*2);
            std::vector<std::list<Point>> bin(numberOfBinsForX);

            for(Point p : P)
                if( p.x() < rightOfFirstGroup)
                    bin[0].push_back(p);
                else {
                    unsigned indexWhereToBeInserted = min(1 + floor((p.x()-rightOfFirstGroup)/(2*l)),bin.size()-1);
                    bin[indexWhereToBeInserted].push_back(p);
                }

            for(unsigned pos = 0; pos < bin.size(); pos++) {

                if( bin[pos].size() == 0)
                    continue;

                std::list<Point> bestSoFarForThisVerticalStrip;
                unsigned answerY = UINT_MAX;
                long double minY = bin[pos].front().y(), maxY = bin[pos].front().y();

                for(Point p : bin[pos]) {
                    if( p.y() < minY )  minY = p.y();
                    if( p.y() > maxY )  maxY = p.y();
                }

                unsigned numberOfSquares =  1 + ceil((maxY - minY) / (2*l));

                for(int shiftInY = 0; shiftInY < l; shiftInY++) {

                    std::list<Point> centersForThisShiftinY;
                    std::vector<std::vector<Point>> squares(numberOfSquares);
                    long double topOfFirstSquare = minY + (shiftInY*2);

                    for( Point p : bin[pos])
                        if( p.y() < topOfFirstSquare)
                            squares[0].push_back(p);
                        else {
                            unsigned indexWhereToBeInserted = min(1 + floor((p.y()-topOfFirstSquare)/(2*l)),squares.size()-1);
                            squares[indexWhereToBeInserted].push_back(p);
                        }
                    for( unsigned i = 0; i < squares.size(); i++) {
                        if( squares[i].empty() )
                            continue;
                        std::list<Point> optCoveringForThisSquare;
                        double minimumX = (pos == 0) ? rightOfFirstGroup-2*l : rightOfFirstGroup+(pos-1)*2*l;
                        double maximumX = minimumX + 2*l;
                        double minimumY = (i == 0) ? topOfFirstSquare-2*l : topOfFirstSquare+(i-1)*2*l;
                        double maximumY = minimumY + 2*l;

                        findSolInSquare(squares[i], optCoveringForThisSquare, l, minimumX, maximumX, minimumY, maximumY);
                        appendPoints(centersForThisShiftinY,optCoveringForThisSquare);
                    }

                    if( centersForThisShiftinY.size() < answerY ) {
                        answerY = centersForThisShiftinY.size();
                        bestSoFarForThisVerticalStrip.clear();
                        appendPoints(bestSoFarForThisVerticalStrip,centersForThisShiftinY);
                    }
                }
                appendPoints(tempCForThisShift,bestSoFarForThisVerticalStrip);
            }

            if(tempCForThisShift.size() < answer) {
                answer = tempCForThisShift.size();
                diskCenters.clear();
                appendPoints(diskCenters,tempCForThisShift);
            }
        }

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};
#endif // FCB2001_H
