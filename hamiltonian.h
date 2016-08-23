//**************************************************************************
//************************* Hamiltonian Class ******************************
//**************************************************************************
// Giacomo Torlai, Apr. 2015
//**************************************************************************
//**************************************************************************

#ifndef _HAMILTONIAN_h
#define _HAMILTONIAN_h

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include "MersenneTwister.h"
#include <string>


using namespace std;
using namespace Eigen;

class Hamiltonian {
    
public:
    
    int L;
    string BC;
    MatrixXd A;
    MatrixXd B;
    
    Hamiltonian(int L_,const string & BC_);
    
    void build_Ising(double h, double J);
    void build_RandomIsing(double h_, double J_, MTRand & random);
};

#endif
