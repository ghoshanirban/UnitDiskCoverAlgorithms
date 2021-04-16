#ifndef G1991_H
#define G1991_H

#include <chrono>
#include <list>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

class G1991 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;
    long double D = std::sqrt(2);

    void fastCLPF( const std::list<Point> &P, std::list<Point> &solutionCover) {
        if(P.empty())
            return;

        std::map<int,std::list<Point>> pointsGroupedByI1;

        for ( Point p : P ) {
            int i1ValueOfp = floor(p.x()/D);
            auto it = pointsGroupedByI1.find(i1ValueOfp);

            if(it == pointsGroupedByI1.end()) {
                std::list<Point> newS;
                newS.push_back(p);
                pointsGroupedByI1[i1ValueOfp] = newS;
            } else
                it->second.push_back(p);
        }

        std::vector<std::list<Point>> S;
        S.reserve(pointsGroupedByI1.size());
        for(auto &it : pointsGroupedByI1) {
            S.push_back(it.second);
        }

        std::vector<std::pair<unsigned,Point>> R;
        R.reserve(2*(P.size()/pointsGroupedByI1.size()));
        std::vector<unsigned> countElementsFromSiPresentinR(S.size(),0);
        long unsigned numberOfSiPresentinR = 0, j = 1;

        for(Point p : S[0]) {
            R.push_back(std::make_pair(0,p));
            countElementsFromSiPresentinR[0]++;
        }
        numberOfSiPresentinR++;

        if( S.size () > 1) {
            for(Point p : S[1]) {
                R.push_back(std::make_pair(1,p));
                countElementsFromSiPresentinR[1]++;
            }
            numberOfSiPresentinR++;
        } else
            j = 0;

        while( !R.empty() ) {
            Point q = R.front().second;
            for(auto aPair : R )
                if(aPair.second.x() < q.x())
                    q = aPair.second;

            std::vector<std::pair<unsigned,Point>> tempR;
            tempR.reserve(R.size());
            long double lowerYofSquareTobePlaced = D*floor(q.y()/D), upperYofSquareTobePlaced = lowerYofSquareTobePlaced + D;
            for(auto aPair : R )
                if( aPair.second.x() - q.x() > D ||
                        aPair.second.y() > upperYofSquareTobePlaced ||
                        aPair.second.y() < lowerYofSquareTobePlaced )
                    tempR.push_back(aPair);
                else {
                    countElementsFromSiPresentinR[aPair.first]--;
                    if(countElementsFromSiPresentinR[aPair.first] == 0)
                        numberOfSiPresentinR--;
                }
            R = tempR;

            solutionCover.push_back( Point( (2*q.x() + D)/2, (D/2)+(D*floor(q.y()/D))  ));

            while(j < S.size()-1 && numberOfSiPresentinR <= 1) {
                j++;
                for(Point p : S[j]) {
                    R.push_back(std::make_pair(j,p));
                    countElementsFromSiPresentinR[j]++;
                }
                if(S[j].size() > 0)
                    numberOfSiPresentinR++;
            }
        }
    }


public:
    G1991(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {
        assert(P.size() > 0);
        auto start = std::chrono::high_resolution_clock::now();

        std::list<Point> solution1, solution2, subProblemR1, subProblemR2;

        for(Point p : P)
            if( ((int)floor(p.y()/D))%2 != 0 )
                subProblemR1.push_back(p);
            else
                subProblemR2.push_back(p);

        fastCLPF(subProblemR1,solution1);
        fastCLPF(subProblemR2,solution2);

        for(Point p : solution1)
            diskCenters.push_back(p);

        for(Point p : solution2)
            diskCenters.push_back(p);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // G1991_H



