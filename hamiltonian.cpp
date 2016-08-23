#include "hamiltonian.h"
#include <vector>

//***********************************************************************
// Hamiltoanian class consctructor
//***********************************************************************

Hamiltonian::Hamiltonian(int L_, const string & BC_) {
    
    L = L_;
    BC.assign(BC_);

    cout << "***********************" <<endl;
    cout << "Lattice Size: " << L << endl;
    cout << "Boundary conditions: " << BC << endl << endl;
    
}

//***********************************************************************
// Build transverse field Ising model
//***********************************************************************

void Hamiltonian::build_Ising(double h, double J) {
    
    A.setZero(L,L);
    B.setZero(L,L);
    cout << "Magnetic field: " << h << endl;
    cout << "Interaction strength: " << J << endl;
    cout << "..building transverse field Ising model" << endl;
    for(int i=0;i<L-2;i++) {
        A(i,i) = -h;
        A(i,i+1) = -0.5*J;
        A(i+1,i) = A(i,i+1);
        B(i,i+1) = -0.5*J;
        B(i+1,i) = -B(i,i+1);
    }
    
    A(L-1,L-1) = -h;
    A(L-2,L-2) = -h;
    A(L-2,L-1) = -0.5*J;
    A(L-1,L-2) = -0.5*J;
    B(L-2,L-1) = -0.5*J;
    B(L-1,L-2) = 0.5*J;
    
    if(BC.compare("Periodic") == 0) {
        A(L,1) = -0.5*J;
        A(1,L) = -0.5*J;
        B(L,1) = 0.5*J;
        B(1,L) = -0.5*J;
    }
    
    cout << "..done" << endl << endl;
}

//***********************************************************************
// Build random-field random bond-Ising model
//***********************************************************************

void Hamiltonian::build_RandomIsing(double h_, double J_, MTRand & random) {
    
    A.setZero(L,L);
    B.setZero(L,L);
    cout << "Magnetic field range: [0," << h_ << "]" << endl;
    cout << "Interaction strength bound: [0," << J_ << "]" << endl;
    cout << "..building random-field random-bond Ising model" << endl;
    vector<double> h;
    vector<double> J;

    for(int i=0;i<L;i++) {
        J.push_back(random.randDblExc(J_));
        h.push_back(random.randDblExc(h_));
    }
    
    //Initialize the matrices A and B with NN interactions only
    for(int i=0;i<L-1;i++) {
        A(i,i) = -h[i];
        A(i,i+1) = -0.5*J[i];
        A(i+1,i) = A(i,i+1);
        B(i,i+1) = -0.5*J[i];
        B(i+1,i) = -B(i,i+1);
    }
    A(L-1,L-1) = -h[L-1];
    A(L-2,L-1) = -0.5*J[L-2];
    A(L-1,L-2) = -0.5*J[L-2];
    B(L-2,L-1) = -0.5*J[L-2];
    B(L-1,L-2) = 0.5*J[L-2];
    
    if(BC.compare("Periodic") == 0) {
        A(L,1) = -0.5*J[L-1];
        A(1,L) = -0.5*J[L-1];
        B(L,1) = 0.5*J[L-1];
        B(1,L) = -0.5*J[L-1];
    }
    
    cout << "..done" << endl << endl;

}




