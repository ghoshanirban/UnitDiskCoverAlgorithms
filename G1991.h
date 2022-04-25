#ifndef G1991_2_H
#define G1991_2_H

#include <chrono>
#include <list>
#include <vector>
#include <unordered_map>

#include <CGAL/Cartesian.h>

typedef CGAL::Cartesian<double>::Point_2 Point;

class G1991 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;
    long double D = std::sqrt(2);

   void fastCLPF(const std::list<Point> &Q, std::list<Point> &solutionCover) const {
        if(Q.empty())
            return;

        std::unordered_map<int,std::list<Point>> pointsGroupedByI1;

        for ( const Point &p : Q ) {
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
        R.reserve(2*(Q.size() / pointsGroupedByI1.size()));
        std::vector<unsigned> countElementsFromSiPresentinR(S.size(),0);
        long unsigned numberOfSiPresentinR = 0, j = 1;

        for(const Point &p : S[0]) {
            R.emplace_back(std::make_pair(0,p));
            countElementsFromSiPresentinR[0]++;
        }
        numberOfSiPresentinR++;

        if( S.size () > 1) {
            for(const Point &p : S[1]) {
                R.emplace_back(std::make_pair(1,p));
                countElementsFromSiPresentinR[1]++;
            }
            numberOfSiPresentinR++;
        } else
            j = 0;

        while( !R.empty() ) {
            Point q = R.front().second;
            for(const auto &aPair : R )
                if(aPair.second.x() < q.x())
                    q = aPair.second;

            std::vector<std::pair<unsigned,Point>> tempR;
            tempR.reserve(R.size());
            long double lowerYofSquareTobePlaced = D*floor(q.y()/D),
                        upperYofSquareTobePlaced = lowerYofSquareTobePlaced + D;
            for(const auto &aPair : R )
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

            solutionCover.emplace_back( Point( (2*q.x() + D)/2, (D/2)+(D*floor(q.y()/D))  ));

            while(j < S.size()-1 && numberOfSiPresentinR <= 1) {
                j++;
                for(const Point &p : S[j]) {
                    R.emplace_back(std::make_pair(j,p));
                    countElementsFromSiPresentinR[j]++;
                }
                if(!S[j].empty())
                    numberOfSiPresentinR++;
            }
        }
    }


public:
    G1991(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) { }

    double execute() {
        assert(!P.empty());
        auto start = std::chrono::high_resolution_clock::now();

        std::list<Point> subProblemR1, subProblemR2;

        std::unordered_map<int,std::list<Point>> pointsGroupedByI2;

        for ( const Point &p : P ) {
            int i2ValueOfp = floor(p.y()/D);
            auto it = pointsGroupedByI2.find(i2ValueOfp);

            if(it == pointsGroupedByI2.end()) {
                std::list<Point> newS;
                newS.push_back(p);
                pointsGroupedByI2[i2ValueOfp] = newS;
            } else
                it->second.push_back(p);
        }

        for(auto &it : pointsGroupedByI2) {
            fastCLPF( it.second, diskCenters);
        }

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};

#endif // G1991_2_H



