## Experiments with Unit Disk Cover Algorithms for Covering Massive Pointsets

Given a set of n points in the plane, the Unit Disk Cover (UDC) problem asks to compute
the minimum number of unit disks required to cover the points, along with a placement
of the disks. The problem is NP-hard and several approximation algorithms have been
designed over the last three decades. In this paper, we have engineered and experimentally
compared practical performances of some of these algorithms on massive pointsets. The
goal is to investigate which algorithms run fast and give good approximation in practice.
We present a simple 7-approximation algorithm for UDC that runs in O(n) expected
time and uses O(s) extra space, where s denotes the size of the generated cover. In our
experiments, it turned out to be the speediest of all. We also present two heuristics to
reduce the sizes of covers generated by it without slowing it down by much.
To our knowledge, this is the first work that experimentally compares algorithms for the
UDC problem. Experiments with them using massive pointsets (in the order of millions)
throw light on their practical uses. We share the engineered algorithms via GitHub for
broader uses and future research in the domain of geometric optimization.

Journal: Computational Geometry, Volume 109, February 2023, 101925.

@article{friederich2022experiments,
  title={Experiments with Unit Disk Cover Algorithms for Covering Massive Pointsets},
  author={Friederich, Rachel and Ghosh, Anirban and Graham, Matthew and Hicks, Brian and Shevchenko, Ronald},
  journal={Computational Geometry},
  pages={101925},
  year={2022},
  publisher={Elsevier}
}

## Built with

* [GCC](https://gcc.gnu.org/)
* [C++17](https://en.cppreference.com/w/cpp/17)
* [CGAL](https://www.cgal.org/)
* [GMP](https://gmplib.org/)
* [MPFR](https://www.mpfr.org/)


## Acknowledgement

Research on this project was supported by the University of North Florida Academic Technology Grant and partially by NSF Award CCF-1947887.
