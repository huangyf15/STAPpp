# stappp
[![build_status](https://travis-ci.com/gwy15/STAPpp.svg?token=xVpBdFRd1VWbmgq6LXmh&branch=master)](https://travis-ci.com/gwy15/STAPpp/)

STAP++ is a C++ finite element method (FEM) code whose input/output data files are the same as STAP90. STAP90 is a FEM code in Fortran 90 provided by our textbook (Xiong Zhang, Tianshu Wang. Computational Dynamics, Tsinghau University Press, 2007; Xiong Zhang, Tianshu Wang, Yan Liu. Computational Dynamics (2nd edition), Tsinghua University Press, 2015)

STAP++ is developed for the course "Finite Element Method" delivered by Professor Xiong Zhang (xzhang@tsinghua.edu.cn) in the School of Aerospace Engineering at Tsinghua University. It helps students to understand the basic implementation techniques of the FEM, and servers as a starting point for students to practice programming the FEM.

STAP++ is developed and maintained by the Computational Dynamics Laboratory (http://www.comdyn.cn/), School of Aerospace Engineering, Tsinghua University, China. Your feedbacks are welcome.

The documentation of STAP++ can be found at https://xzhang66.github.io/stappp/index.html.

## before cmake

### dependencies

+ eigen-3
    
    Please download source code [here](http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz), unzip and move the directory to `src/eigen-3`

+ MKL

    Please download setup package [here](https://registrationcenter.intel.com/en/forms/?productid=2558&licensetype=2).

    