#include "Algorithms.h"
#include <iostream>

/**********************************************************************/
/* cgal and boost libraries are needed to compile and run the project */
/**********************************************************************/

/*
	We advise using the latest g++ compiler (v9 and above) with the -O3 flag
	for the fastest execution times.
*/

int main() {
    /* It is assumed that an input text file contains 'n' lines exactly, where
       'n' denotes the size of the input point set. Every line in the file should
       contain two doubles: x and y-coordinates of the ith point, separated by a single space.
    */

    std::ifstream input("world.txt"); // this point set is zipped and uploaded to this repository
    if(input.fail()){
        std::cout << "File does not exist!" << std::endl;
        exit(1);
    }

    // Read the input file to the vector 'points'
    std::vector<Point> P;
    long double x,y;
    while(input >> x >> y )
        P.emplace_back(Point(x,y));

    input.close();

    std::list<Point> C; // a list to hold the disk centers
    FASTCOVER ob(P,C);
    std::cout << "Time taken: " << ob.execute() << " seconds, "; // prints the execution time
    std::cout << "Disks placed: " << C.size() << std::endl;  // prints the number disks the algorithm has placed

    std::ofstream output("cover.txt"); // send the disk centers to an output file
    for(const Point &c : C)
        output << c.x() << " " << c.y() << std::endl;

    output.close();

    return EXIT_SUCCESS;
}

