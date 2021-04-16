#include "BLMS2017.h"
#include "DGT2018.h"
#include "HM1985.h"
#include "LL2014.h"
#include "CCFM1997.h"
#include "G1991.h"
#include "FASTCOVER.h"
#include "FCB2001.h"

#include <iostream>
#include <fstream>

/*
For every engineered algorithm, we have created a different class (refer to the header files). 
The constructors for these classes need two parameters: a vector of the input
points and an output list where the disk centers (a pair of doubles) will be sent to
by the algorithm. This interface works for all the engineered algorithms. 
However, for the algorithms HM1985 and FCB2001, an extra integer parameter
'l' is required (warning: these two algorithms may run forever for large point sets). 
*/

/*
	We advise using the latest g++ compiler (v9 and above) with the -O3 flag
	for the fastest execution times. 
*/

int main() {
	/* Please make sure the input file is present and the name is exact
       It is assumed that an input text file contains 'n' lines exactly, where
       'n' denotes the size of the input point set. Every line in the file should
       contain two doubles: x-coordinate of the point i and y-coordinate of
       the point i, separated by a single space.
    */
    
    std::ifstream input("world.txt");   
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
    
    FASTCOVER ob(P,C); // choose an algorithm (for HM1985/FCB2001 an extra third integer parameter is needed)
    
    std::cout << "Time taken: " << ob.execute() << " seconds, "; // prints the execution time
    std::cout << "Disks placed: " << C.size() << std::endl;  // prints the number disks the algorithm has placed
    
    std::ofstream output("cover.txt"); // send the disk centers to an output file
    for(Point c : C)
    	output << c.x() << " " << c.y() << std::endl;
    
    output.close();
       
    return EXIT_SUCCESS;
}

