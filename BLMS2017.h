#ifndef BLMS2017_2_H
#define BLMS2017_2_H

#include <chrono>
#include <list>
#include <vector>
#include <set>

#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<double>::Point_2 Point;

class BLMS2017 {
    std::vector<Point> &P;
    std::list<Point> &diskCenters;

    struct BSTcomp {
        bool operator()(const std::pair<Point, unsigned> &lhs, const std::pair<Point, unsigned> &rhs) const {
            return lhs.first.y() < rhs.first.y()
                   || (lhs.first.y() == rhs.first.y() && lhs.first.x() < rhs.first.x());
        }
    };

    struct Event {
        Point point;
        Point deletion = this->point;

        bool operator<(const Event &b) const {
            return point.x() < b.point.x() || (point.x() == b.point.x() && point.y() < b.point.y());
        }
    };

    static inline bool isCovered(const Point &p, const Point &diskCenter) {
        return CGAL::squared_distance(p, diskCenter) <= 1;
    }

    std::vector<std::pair<Point, bool> > potentialCenters;

    inline bool checkIfCoveredByAnExistingDisk( Point &p, const unsigned index) {
        if (isCovered(p, potentialCenters[4 * index].first)) {
            return true;
        }
        if (isCovered(p, potentialCenters[4 * index + 1].first)) {
            potentialCenters[4 * index + 1].second = true;
            return true;
        }
        if (isCovered(p, potentialCenters[4 * index + 2].first)) {
            potentialCenters[4 * index + 2].second = true;
            return true;
        }
        if (isCovered(p, potentialCenters[4 * index + 3].first)) {
            potentialCenters[4 * index + 3].second = true;
            return true;
        }
        return false;
    }


public:
    BLMS2017(std::vector<Point> &P, std::list<Point> &diskCenters) : P(P), diskCenters(diskCenters) {}

    double execute() {
        assert(!P.empty());

        auto start = std::chrono::high_resolution_clock::now();
        const double sqrt3 = std::sqrt(3), sqrt3Over2 = sqrt3 / 2;
        std::vector<Event> E;
        std::set<std::pair<Point, int>, BSTcomp> BST;

        for (const Point &p: P) {
            E.push_back({p});
            E.push_back({Point(p.x() + 2, p.y()), p});
        }

        std::sort(E.begin(), E.end());

        double x = E.at(0).point.x(), y = E.at(0).point.y();
        unsigned count = 0;

        potentialCenters.emplace_back(std::make_pair(Point(x, y), true));
        potentialCenters.emplace_back(std::make_pair(Point(x + sqrt3, y), false));
        potentialCenters.emplace_back(std::make_pair(Point(x + sqrt3Over2, y + 1.5), false));
        potentialCenters.emplace_back(std::make_pair(Point(x + sqrt3Over2, y - 1.5), false));

        BST.insert(std::make_pair(E.at(0).point, count++));
        E.erase(E.begin());

        for (Event &e: E) {
            if (e.deletion == e.point) {
                auto iter = BST.lower_bound(std::make_pair(e.point, 0));
                if (iter != BST.end()) {
                    unsigned index = iter->second;

                    if( checkIfCoveredByAnExistingDisk(e.point, index) )
                        continue;

                    if (++iter != BST.end()) {
                        index = iter->second;

                        if( checkIfCoveredByAnExistingDisk(e.point, index) )
                            continue;
                    }
                    --iter;
                }

                if (iter != BST.begin()) {
                    --iter;
                    unsigned index = iter->second;
                    if( checkIfCoveredByAnExistingDisk(e.point, index) )
                        continue;

                    if (iter != BST.begin()) {
                        --iter;
                        index = iter->second;
                        if( checkIfCoveredByAnExistingDisk(e.point, index) )
                            continue;
                    }
                }
                x = e.point.x(); y = e.point.y();
                potentialCenters.emplace_back(std::make_pair(Point(x, y), true));
                potentialCenters.emplace_back(std::make_pair(Point(x + sqrt3, y), false));
                potentialCenters.emplace_back(std::make_pair(Point(x + sqrt3Over2, y + 1.5), false));
                potentialCenters.emplace_back(std::make_pair(Point(x + sqrt3Over2, y - 1.5), false));
                BST.insert(std::make_pair(e.point, count++));
            }
            else{
                BST.erase(std::make_pair(e.deletion, 0));
            }
        }

        for(const auto &center : potentialCenters)
            if (center.second)
                diskCenters.push_back(center.first);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;
        return duration.count();
    }
};
#endif // BLMS2017_2_H