#ifndef HM1985_H
#define HM1985_H

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
class HM1985 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;
    int l;

    inline void appendPoints(std::list<Point> &dest, const std::list<Point> &source) {
        for(Point p: source)
            dest.push_back(p);
    }

    bool isCover(const std::vector<Point> &P,  const std::list<Point> &selectedCenters) {

        DelaunayTriangulation T(selectedCenters.begin(),selectedCenters.end());

        for(Point p : P) {
            Vertex_handle handleToTheNearestDiskCenter = T.nearest_vertex(p);
            if(squared_distance(p,handleToTheNearestDiskCenter->point()) > 1.0000000001)
                return false;
        }
        return true;
    }

    void findOPTCoverinSquare(std::vector<Point> &P, std::list<Point> &optCenters,  const int ell) {

        if(P.size() == 1) {
            optCenters.push_back(P[0]);
            return;
        }

        std::vector<Point> candidateDiskCenters;
        candidateDiskCenters.reserve(2*P.size()*P.size());
        std::vector<bool> hasConsidered(P.size());

        for(unsigned i = 0; i < P.size() - 1 ; i++)
            for(unsigned j = i + 1; j < P.size(); j++) {

                PointInCircularKernel pi(P[i].x(),P[i].y()), pj(P[j].x(),P[j].y());
                CircleInCircularKernel c1(pi,1), c2(pj,1);

                std::vector<Intersection_result> intersectionPoints;
                intersection(c1,c2,std::back_inserter(intersectionPoints));

                if(!intersectionPoints.empty())
                    hasConsidered[i] = hasConsidered[j] = true;

                using boostRetVal = std::pair<CGAL::Circular_arc_point_2<CGAL::Filtered_bbox_circular_kernel_2<CGAL::Circular_kernel_2<CGAL::Cartesian<CGAL::Gmpq>,
                                                                         CGAL::Algebraic_kernel_for_circles_2_2<CGAL::Gmpq> > > >, unsigned>;

                for(const auto& element : intersectionPoints) {
                    CGAL::Circular_arc_point_2<Circular_k> p = get<0>( boost::get<boostRetVal>(element) );
                    candidateDiskCenters.push_back(Point(to_double(p.x()), to_double(p.y())));
                }
            }

        for( unsigned i = 0; i < hasConsidered.size(); i++ )
            if( !hasConsidered[i] )
                candidateDiskCenters.push_back(P[i]);

        for(unsigned currentCoverSize = 1; currentCoverSize <= HM1985::min(2*ell*ell,candidateDiskCenters.size()); currentCoverSize++) {

            std::string bitmask(currentCoverSize, '1');
            bitmask.resize(candidateDiskCenters.size(), '0');

            do {
                std::list<Point> selectedCenters;

                for( unsigned pos = 0; pos < bitmask.size(); pos++)
                    if( bitmask[pos] == '1')
                        selectedCenters.push_back(candidateDiskCenters[pos]);

                if( isCover(P,selectedCenters) ) {
                    appendPoints(optCenters,selectedCenters);
                    return;
                }
            } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
        }
    }

public:
    HM1985(std::vector<Point> &P, std::list<Point> &diskCenters, int ell) : P(P), diskCenters(diskCenters), l(ell) { }

    static inline unsigned min(const int a, const int b) {
        return (a < b )? a : b;
    }

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
        unsigned answer = P.size() + 1;

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
                unsigned answerY = bin[pos].size() + 1;

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
                    for( std::vector<Point> pointsInSquare : squares) {
                        if( pointsInSquare.empty() )
                            continue;
                        std::list<Point> optCoveringForThisSquare;
                        findOPTCoverinSquare(pointsInSquare,optCoveringForThisSquare,l);
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

#endif // HM1985_H
