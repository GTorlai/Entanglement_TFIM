//**************************************************************************
//**************************************************************************
//************************* Entanglement Class *****************************
//**************************************************************************
//**************************************************************************
//*******   Calculate the Von Neumann etropy and the entanglement    *******
//*******   spectrum for free fermions in 1d.                        *******
//*******                                                            *******
//*******   Giacomo Torlai, Apr. 2015                                *******
//**************************************************************************
//**************************************************************************

#ifndef _entanglement_h
#define _entanglement_h

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <bitset>
#include "hamiltonian.cpp"


using namespace std;
using namespace Eigen;

//2D Square Lattice classs
class entanglement {

public:
    
    int L;
    
    vector<double> ES; 
    double VN; 
    double SG; 
    
    int bound;
    int n_eig;
    int eig_kept;
    int rod;

    MatrixXd Gamma;
    MatrixXd Pi;
    VectorXd Energy;

    //Contructore
    entanglement(int L_, char n_eig_);

    //Main Functions
    vector<double> get_ES(const Hamiltonian & H);
    void measure(vector<double> & ES); 
    void accumulate(vector<double> & es); 

    void record(ofstream & file,double h);

    //Utilities
    int myPow (int x,int p);
    string build_fileName(string & model,double h);
};

#endif
