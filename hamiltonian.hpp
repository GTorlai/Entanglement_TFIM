#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <random>


class Hamiltonian{

  int L_;
  std::string boundary_conditions_;
  std::mt19937 rgen_;

public:
  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  
  Hamiltonian(int L,std::string & boundary_conditions):L_(L),
                boundary_conditions_(boundary_conditions){
    rgen_.seed(13579);
    std::cout << "***********************" <<std::endl;
    std::cout << "Lattice Size: " << L << std::endl;
    std::cout << "Boundary conditions: " << boundary_conditions_ << std::endl << std::endl;
  }

  void Ising(double h){

    double J = 1.0;
    A_.setZero(L_,L_);
    B_.setZero(L_,L_);
    std::cout << "Magnetic field: " << h << std::endl;
    std::cout << "..building transverse field Ising model" << std::endl;
    for(int i=0;i<L_-2;i++) {
      A_(i,i) = -h;
      A_(i,i+1) = -0.5*J;
      A_(i+1,i) = A_(i,i+1);
      B_(i,i+1) = -0.5*J;
      B_(i+1,i) = -B_(i,i+1);
    }
    
    A_(L_-1,L_-1) = -h;
    A_(L_-2,L_-2) = -h;
    A_(L_-2,L_-1) = -0.5*J;
    A_(L_-1,L_-2) = -0.5*J;
    B_(L_-2,L_-1) = -0.5*J;
    B_(L_-1,L_-2) = 0.5*J;
    
    if(boundary_conditions_ == "periodic"){//.compare("periodic") == 0) {
      A_(L_,1) = -0.5*J;
      A_(1,L_) = -0.5*J;
      B_(L_,1) = 0.5*J;
      B_(1,L_) = -0.5*J;
    }


  }


//***********************************************************************
// Build random-field random bond-Ising model
//***********************************************************************

//void Hamiltonian::build_RandomIsing(double h_, double J_, MTRand & random) {
//    
//    A.setZero(L,L);
//    B.setZero(L,L);
//    cout << "Magnetic field range: [0," << h_ << "]" << endl;
//    cout << "Interaction strength bound: [0," << J_ << "]" << endl;
//    cout << "..building random-field random-bond Ising model" << endl;
//    vector<double> h;
//    vector<double> J;
//
//    for(int i=0;i<L;i++) {
//        J.push_back(random.randDblExc(J_));
//        h.push_back(random.randDblExc(h_));
//    }
//    
//    //Initialize the matrices A and B with NN interactions only
//    for(int i=0;i<L-1;i++) {
//        A(i,i) = -h[i];
//        A(i,i+1) = -0.5*J[i];
//        A(i+1,i) = A(i,i+1);
//        B(i,i+1) = -0.5*J[i];
//        B(i+1,i) = -B(i,i+1);
//    }
//    A(L-1,L-1) = -h[L-1];
//    A(L-2,L-1) = -0.5*J[L-2];
//    A(L-1,L-2) = -0.5*J[L-2];
//    B(L-2,L-1) = -0.5*J[L-2];
//    B(L-1,L-2) = 0.5*J[L-2];
//    
//    if(BC.compare("Periodic") == 0) {
//        A(L,1) = -0.5*J[L-1];
//        A(1,L) = -0.5*J[L-1];
//        B(L,1) = 0.5*J[L-1];
//        B(1,L) = -0.5*J[L-1];
//    }
//    
//    cout << "..done" << endl << endl;
//
//}

};

#endif


