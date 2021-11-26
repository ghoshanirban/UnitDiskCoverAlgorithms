//
// Created by Ghosh, Anirban on 11/9/21.
//

#ifndef UDC_EXPERIMENT_H
#define UDC_EXPERIMENT_H

#include "Algorithms.h"
#include <CGAL/point_generators_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/random_convex_set_2.h>

#include <random>

typedef CGAL::Cartesian<double>::Point_2 Point;

struct result {
    long double disks;
    double time;
};

using namespace std;

void generatePointsInAnnulus(const unsigned n, const double r2, const double r1, std::vector<Point> &points) {
    assert(r2 > r1);
    assert(n > 1);

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distributionR(r1,r2), distributionT(0,1);

    for(unsigned i = 0; i < n; i++) {
        double t = 2 * M_PI * distributionT(generator);
        double r = distributionR(generator);
        points.emplace_back(Point(r*cos(t),r*sin(t)));
    }
}

void plotGraphAll(const unsigned type, const string &gName, const unsigned area, const unsigned orderOf,
                  const unsigned numberOfTimes, const unsigned N, const bool displayPDF, double r1 = 0.95) {
    unsigned nALG = 9;

    unsigned ALG_G1991     = 0;
    unsigned ALG_CCFM1997  = 1;
    unsigned ALG_LL2014    = 2;
    unsigned ALG_LL2014_P1 = 3;
    unsigned ALG_BLMS2017  = 4;
    unsigned ALG_DGT2018   = 5;
    unsigned ALG_FASTCOVER_NAIVE = 6;
    unsigned ALG_FASTCOVER_H1 = 7;
    unsigned ALG_FASTCOVER = 8;
    //unsigned N = 10;

    result plotData[nALG][N];
    for(unsigned i = 0; i < nALG; i++)
        for(unsigned j = 0; j < N; j++) {
            plotData[i][j].disks = 0;
            plotData[i][j].time = 0.0;
        }

    for(unsigned i = 0; i < N; i++) {
        for(unsigned sample = 0; sample <= numberOfTimes; sample++) {
            cout << "===========================" << endl;
            std::vector<Point> points;

            if(type == 0) {
                typedef CGAL::Creator_uniform_2<double,Point>  Creator;
                CGAL::Random_points_in_square_2<Point,Creator> g(std::sqrt(area)/2);
                std::copy_n( g,(i+1)*orderOf, std::back_inserter(points));
            }
            else if(type == 1) {
                typedef CGAL::Creator_uniform_2<double,Point>  Creator;
                CGAL::Random_points_in_disc_2<Point,Creator> g(std::sqrt(area/M_PI));
                std::copy_n( g,(i+1)*orderOf, std::back_inserter(points));
            }
            else if(type == 2){
                typedef CGAL::Creator_uniform_2<double,Point>  Creator;
                CGAL::Random_points_in_square_2<Point,Creator> g(std::sqrt(area)/2);
                CGAL::random_convex_set_2((i+1)*orderOf,std::back_inserter(points),g);
            }
            else if(type == 3) {
                generatePointsInAnnulus((i+1)*orderOf,1000,r1*1000,points);
            }

            std::vector<Point> P(points.size());

            cout << i+1 << " G-1991  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                G1991 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_G1991][i].time  += time;
                if(sample > 0 ) plotData[ALG_G1991][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " CCFM-1997  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                CCFM1997 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_CCFM1997][i].time  += time;
                if(sample > 0 ) plotData[ALG_CCFM1997][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " BLMS-2017  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                BLMS2017 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_BLMS2017][i].time  += time;
                if(sample > 0 ) plotData[ALG_BLMS2017][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " LL2014-1PASS  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                LL2014_1PASS ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_LL2014_P1][i].time  += time;
                if(sample > 0 ) plotData[ALG_LL2014_P1][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " LL-2014  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                LL2014 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_LL2014][i].time  += time;
                if(sample > 0 ) plotData[ALG_LL2014][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " DGT-2018  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                DGT2018 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_DGT2018][i].time  += time;
                if(sample > 0 ) plotData[ALG_DGT2018][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " FAST-COVER  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                FASTCOVER ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_FASTCOVER_NAIVE][i].time  += time;
                if(sample > 0 ) plotData[ALG_FASTCOVER_NAIVE][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " FAST-COVER-P  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                FASTCOVER_P ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_FASTCOVER_H1][i].time  += time;
                if(sample > 0 ) plotData[ALG_FASTCOVER_H1][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " FAST-COVER-PP  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                FASTCOVER_PP ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_FASTCOVER][i].time  += time;
                if(sample > 0 ) plotData[ALG_FASTCOVER][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

        }
    }

    // taking the average over all runs
    for(unsigned i = 0; i < nALG; i++)
        for(unsigned j = 0; j < N; j++) {
            plotData[i][j].disks = plotData[i][j].disks/numberOfTimes;
            plotData[i][j].time = plotData[i][j].time/numberOfTimes;
        }

    ///////////////////////////////////////////
    string fName = gName + "Time.tex";
    FILE *fp = fopen(fName.c_str(),"w");
    fprintf(fp,"\\documentclass[tikz]{standalone} \n\\usepackage{pgfplots} \n\\usetikzlibrary{plotmarks} \n \n\n\\begin{document}\n");
    fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");
    fprintf(fp,"\\begin{axis}[legend pos=north west,xlabel=$n$ (in millions), xtick={1,2,3,4,5,6,7,8,9,10},ylabel near ticks, ylabel={Average execution time (in seconds)},legend style={nodes={scale=0.6, transform shape}}]");

    for(unsigned i = 0; i < nALG; i++) {
        if(i == ALG_G1991)
            fprintf(fp,"\n\n\\addplot[ color=blue,mark=star] coordinates {\n");
        else if(i==ALG_BLMS2017)
            fprintf(fp,"\n\n\\addplot[ color=black,mark=square*] coordinates {\n");
        else if(i==ALG_LL2014)
            fprintf(fp,"\n\n\\addplot[ color=orange,mark=triangle*] coordinates {\n");
        else if(i==ALG_CCFM1997)
            fprintf(fp,"\n\n\\addplot[ color=red,mark=diamond*] coordinates {\n");
        else if(i==ALG_DGT2018)
            fprintf(fp,"\n\n\\addplot[ color=teal,mark=pentagon*] coordinates {\n");
        else if(i==ALG_FASTCOVER)
            fprintf(fp,"\n\n\\addplot[ color=purple,mark=*] coordinates {\n");
        else if(i==ALG_FASTCOVER_H1)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=oplus*] coordinates {\n");
        else if(i==ALG_FASTCOVER_NAIVE)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=otimes*] coordinates {\n");
        else if(i==ALG_LL2014_P1)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=square*] coordinates {\n");

        for(unsigned j = 0; j < N; j++)
            fprintf(fp,"( %u, %.4f )\n", (j+1), plotData[i][j].time);
        fprintf(fp,"}node [pos=1.15, above left] {};");
    }
    fprintf(fp,"\n\\end{axis} \n legend style={at={(axis cs:0.5,1)},anchor=south west}");
    fprintf(fp,"\n\n\\end{tikzpicture}");
    fprintf(fp,"\n\n\\end{document}");
    fclose(fp);
    if(displayPDF) {
        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        cout << "PDF generation terminated...\n";
        command = "atril " + gName + "Time.pdf &";
        system(command.c_str());
    }
    ////////////////////////////////////////////////////
    fName = gName + "TimeNoLL6P.tex";
    fp = fopen(fName.c_str(),"w");
    fprintf(fp,"\\documentclass[tikz]{standalone} \n\\usepackage{pgfplots} \n\\usetikzlibrary{plotmarks} \n \n\n\\begin{document}\n");
    fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");
    fprintf(fp,"\\begin{axis}[legend pos=north west,xlabel=$n$ (in millions), xtick={1,2,3,4,5,6,7,8,9,10},ylabel near ticks, ylabel={Average execution time (in seconds)},legend style={nodes={scale=0.6, transform shape}}]");

    for(unsigned i = 0; i < nALG; i++) {
        if(i == ALG_G1991)
            fprintf(fp,"\n\n\\addplot[ color=blue,mark=star] coordinates {\n");
        else if(i==ALG_BLMS2017)
            fprintf(fp,"\n\n\\addplot[ color=black,mark=square*] coordinates {\n");
        else if(i==ALG_LL2014)
            continue;
        else if(i==ALG_CCFM1997)
            fprintf(fp,"\n\n\\addplot[ color=red,mark=diamond*] coordinates {\n");
        else if(i==ALG_DGT2018)
            fprintf(fp,"\n\n\\addplot[ color=teal,mark=pentagon*] coordinates {\n");
        else if(i==ALG_FASTCOVER)
            fprintf(fp,"\n\n\\addplot[ color=purple,mark=*] coordinates {\n");
        else if(i==ALG_FASTCOVER_H1)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=oplus*] coordinates {\n");
        else if(i==ALG_FASTCOVER_NAIVE)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=otimes*] coordinates {\n");
        else if(i==ALG_LL2014_P1)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=square*] coordinates {\n");

        for(unsigned j = 0; j < N; j++)
            fprintf(fp,"( %u, %.4f )\n", (j+1), plotData[i][j].time);
        fprintf(fp,"}node [pos=1.15, above left] {};");
    }
    fprintf(fp,"\n\\end{axis} \n legend style={at={(axis cs:0.5,1)},anchor=south west}");
    fprintf(fp,"\n\n\\end{tikzpicture}");
    fprintf(fp,"\n\n\\end{document}");
    fclose(fp);
    if(displayPDF) {
        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        cout << "PDF generation terminated...\n";
        command = "atril " + gName + "Time.pdf &";
        system(command.c_str());
    }
    //////////////
    fName = gName + "Disks.tex";
    fp = fopen(fName.c_str(),"w");
    fprintf(fp,"\\documentclass[tikz]{standalone} \n\\usepackage{pgfplots} \n\\usetikzlibrary{plotmarks} \n \n\n\\begin{document}\n");
    fprintf(fp,"\n\n\n\\begin{tikzpicture}\n\n");
    fprintf(fp,"\\begin{axis}[legend pos=north west,xlabel=$n$ (in millions), xtick={1,2,3,4,5,6,7,8,9,10},ylabel near ticks, ylabel={Disks placed on average},legend style={nodes={scale=0.6, transform shape}}	]");

    for(unsigned i = 0; i < nALG; i++) {
        if(i == ALG_G1991)
            fprintf(fp,"\n\n\\addplot[ color=blue,mark=star] coordinates {\n");
        else if(i==ALG_BLMS2017)
            fprintf(fp,"\n\n\\addplot[ color=black,mark=square*] coordinates {\n");
        else if(i==ALG_LL2014)
            fprintf(fp,"\n\n\\addplot[ color=orange,mark=triangle*] coordinates {\n");
        else if(i==ALG_CCFM1997)
            fprintf(fp,"\n\n\\addplot[ color=red,mark=diamond*] coordinates {\n");
        else if(i==ALG_DGT2018)
            fprintf(fp,"\n\n\\addplot[ color=teal,mark=pentagon*] coordinates {\n");
        else if(i==ALG_FASTCOVER)
            fprintf(fp,"\n\n\\addplot[ color=purple,mark=*] coordinates {\n");
        else if(i==ALG_FASTCOVER_H1)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=oplus*] coordinates {\n");
        else if(i==ALG_FASTCOVER_NAIVE)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=otimes*] coordinates {\n");
        else if(i==ALG_LL2014_P1)
            fprintf(fp,"\n\n\\addplot[ color=black,mark options={fill=yellow},mark=square*] coordinates {\n");

        for(unsigned j = 0; j < N; j++)
            fprintf(fp,"( %u, %.4Lf )\n", (j+1), plotData[i][j].disks);
        fprintf(fp,"}node [pos=1.15, above left] {};");
    }
    fprintf(fp,"\n\\end{axis} \n legend style={at={(axis cs:0.5,1)},anchor=south west}");
    fprintf(fp,"\n\n\\end{tikzpicture}");
    fprintf(fp,"\n\n\\end{document}");
    fclose(fp);

    if(displayPDF) {
        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        cout << "PDF generation terminated...\n";
        command = "atril " + gName + "Disks.pdf &";
        system(command.c_str());
    }

    //////////////
    fName = gName + "Table.tex";
    fp = fopen(fName.c_str(),"w");
    fprintf(fp,"\\documentclass{standalone} \\begin{document} \\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c| } \n \\hline \n ");
    fprintf(fp," $n$  & G-1991 & CCFM-1997 & LL-2014 & LL-2014-1P & BLMS-2017 & DGT-2018 & \\textsc{FastCover} & \\textsc{FastCover\\texttt{+}} & \\textsc{FastCover\\texttt{++}} \\\\ \n \\hline \\hline \n");

    for(unsigned j = 0; j < N; j++) {
        fprintf(fp, "$%u M$ ", (j+1));

        for(unsigned i = 0; i < nALG; i++)
            fprintf(fp," & %.2Lf, %.2f ",plotData[i][j].disks,plotData[i][j].time);

        fprintf(fp,"\\\\ \n \\hline \n");
    }

    fprintf(fp,"\\end{tabular} \n \\end{document}");
    fclose(fp);

    if(displayPDF) {
        cout << "\nOutput PDF generation started...\n";
        string command = "pdflatex " + fName + " > /dev/null";
        system(command.c_str());
        cout << "PDF generation terminated...\n";
        command = "atril " + gName + "Table.pdf &";
        system(command.c_str());
    }
}

void loadPointsFrom(const std::string &file, std::vector<Point> &P) {
    std::ifstream input(file);
    if(input.fail()){
        std::cout << "File does not exist!" << std::endl;
        exit(1);
    }

    long double x,y;
    while(input >> x >> y )
        P.emplace_back(Point(x,y));

    input.close();
}

void plotReal(std::string gName) {
    unsigned nALG = 9;

    unsigned ALG_G1991     = 0;
    unsigned ALG_CCFM1997  = 1;
    unsigned ALG_LL2014    = 2;
    unsigned ALG_LL2014_P1 = 3;
    unsigned ALG_BLMS2017  = 4;
    unsigned ALG_DGT2018   = 5;
    unsigned ALG_FASTCOVER_NAIVE = 6;
    unsigned ALG_FASTCOVER_H1 = 7;
    unsigned ALG_FASTCOVER = 8;
    unsigned N = 11;

    result plotData[nALG][N];
    for(unsigned i = 0; i < nALG; i++)
        for(unsigned j = 0; j < N; j++) {
            plotData[i][j].disks = 0;
            plotData[i][j].time = 0.0;
        }

    unsigned numberOfTimes = 5;
    for(unsigned i = 0; i < N; i++) {
        for(unsigned sample = 0; sample <= numberOfTimes; sample++) {
            cout << "===========================" << endl;
            std::vector<Point> points;

            if(i==0)
                loadPointsFrom("pointsets/birch3.txt",points);
            else if(i==1)
                loadPointsFrom("pointsets/monalisa.txt",points);
            else if(i==2)
                loadPointsFrom("pointsets/usa.txt",points);
            else if(i==3)
                loadPointsFrom("pointsets/KDDCU2D.txt",points);
            else if(i==4)
                loadPointsFrom("pointsets/europe.txt",points);
            else if(i==5)
                loadPointsFrom("pointsets/wildfires.txt",points);
            else if(i==6)
                loadPointsFrom("pointsets/world.txt",points);
            else if(i==7)
                loadPointsFrom("pointsets/china.txt",points);
            else if(i==8)
                loadPointsFrom("pointsets/nyctaxi.txt",points);
            else if(i==9)
                loadPointsFrom("pointsets/uber.txt",points);
            else if(i==10)
                loadPointsFrom("pointsets/hail2015.txt",points);

            std::vector<Point> P(points.size());

            cout << i+1 << " G-1991  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                G1991 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_G1991][i].time  += time;
                if(sample > 0 ) plotData[ALG_G1991][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " CCFM-1997  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                CCFM1997 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_CCFM1997][i].time  += time;
                if(sample > 0 ) plotData[ALG_CCFM1997][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " BLMS-2017  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                BLMS2017 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_BLMS2017][i].time  += time;
                if(sample > 0 ) plotData[ALG_BLMS2017][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " LL2014-1PASS  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                LL2014_1PASS ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_LL2014_P1][i].time  += time;
                if(sample > 0 ) plotData[ALG_LL2014_P1][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " LL-2014  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                LL2014 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_LL2014][i].time  += time;
                if(sample > 0 ) plotData[ALG_LL2014][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " DGT-2018  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                DGT2018 ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_DGT2018][i].time  += time;
                if(sample > 0 ) plotData[ALG_DGT2018][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " FAST-COVER  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                FASTCOVER ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_FASTCOVER_NAIVE][i].time  += time;
                if(sample > 0 ) plotData[ALG_FASTCOVER_NAIVE][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " FAST-COVER+  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                FASTCOVER_P ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_FASTCOVER_H1][i].time  += time;
                if(sample > 0 ) plotData[ALG_FASTCOVER_H1][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

            cout << i+1 << " FAST-COVER++  ";
            {
                P.clear();
                list<Point> centersOfPlacedDisks;
                for(const Point &p : points)
                    P.emplace_back(p);
                FASTCOVER_PP ob(P,centersOfPlacedDisks);
                double time = ob.execute();
                if(sample > 0 ) plotData[ALG_FASTCOVER][i].time  += time;
                if(sample > 0 ) plotData[ALG_FASTCOVER][i].disks += centersOfPlacedDisks.size();
                cout <<  time << ",  " << centersOfPlacedDisks.size() << endl;
            }

        }
    }

    // taking the average over all runs
    for(unsigned i = 0; i < nALG; i++)
        for(unsigned j = 0; j < N; j++) {
            plotData[i][j].disks = plotData[i][j].disks/numberOfTimes;
            plotData[i][j].time = plotData[i][j].time/numberOfTimes;
        }

    ///////////////////////////////////////////
    string fName = gName + "Table.tex";
    FILE *fp = fopen(fName.c_str(),"w");
    fprintf(fp,"\\documentclass{standalone} \\begin{document} \\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c| } \n \\hline \n ");
    fprintf(fp," Dataset  & G-1991 & CCFM-1997 & LL-2014 & LL-2014-1P & BLMS-2017 & DGT-2018 & \\textsc{FastCover} & \\textsc{FastCover\\texttt{+}} & \\textsc{FastCover\\texttt{++}} \\\\ \n \\hline \\hline \n");

    for(unsigned j = 0; j < N; j++) {

        if(j==0)
            fprintf(fp, "\\texttt{birch3}");
        else if(j==1)
            fprintf(fp, "\\texttt{monalisa}");
        else if(j==2)
            fprintf(fp, "\\texttt{usa}");
        else if(j==3)
            fprintf(fp, "\\texttt{KDDCU2D}");
        else if(j==4)
            fprintf(fp, "\\texttt{europe}");
        else if(j==5)
            fprintf(fp, "\\texttt{wildfires}");
        else if(j==6)
            fprintf(fp, "\\texttt{world}");
        else if(j==7)
            fprintf(fp, "\\texttt{china}");
        else if(j==8)
            fprintf(fp, "\\texttt{nyctaxi}");
        else if(j==9)
            fprintf(fp, "\\texttt{uber}");
        else if(j==10)
            fprintf(fp, "\\texttt{hail2015}");

        for(unsigned i = 0; i < nALG; i++)
            fprintf(fp," & %.2Lf, %.2f ",plotData[i][j].disks,plotData[i][j].time);

        fprintf(fp,"\\\\ \n \\hline \n");
    }

    fprintf(fp,"\\end{tabular} \n \\end{document}");
    fclose(fp);
}
#endif //UDC_EXPERIMENT_H
