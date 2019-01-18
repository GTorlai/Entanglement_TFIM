//**************************************************************************
//**************************************************************************
//***** ******************** Entanglement Class *****************************
//**************************************************************************
//**************************************************************************
//*******   Calculate the Von Neumann etropy and the entanglement    *******
//*******   spectrum for free fermions in 1d.                        *******
//*******                                                            *******
//*******   Giacomo Torlai, Apr. 2015                                *******
//**************************************************************************
//**************************************************************************

#ifndef ENTANGLEMENT_HPP
#define ENTANGLEMENT_HPP

#include <boost/format.hpp>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <bitset>
#include "hamiltonian.hpp"

class Entanglement{
  
  int L_;
  int bound_;
  int n_eig_;
  int eig_keep_;
  
  Eigen::MatrixXd Gamma_;
  Eigen::MatrixXd Pi_;
  Eigen::VectorXd Energy_;
 
  std::vector<double> ES_;
  double VN_;
  double SG_;
  
public:
  Entanglement(int L,int n_eig,int eig_keep):L_(L),n_eig_(n_eig),eig_keep_(eig_keep){
    bound_ = pow(2,n_eig);
    std::cout<<bound_<<std::endl;
    ES_.assign(bound_,0.0);
  }



 
  void ComputeEntanglementSpectrum(const Hamiltonian &H) {
    
    double Norm=0;
    double epsilon=0.000000000000001;
    double dot_prod;
    
    std::vector<double> eigen0;
    eigen0.assign(L_/2,0.0);
    std::vector<double> eigen1;
    eigen1.assign(n_eig_,0.0);
    std::vector<char> bin;
    bin.assign(n_eig_,0);
    std::vector<double> lambda;
    lambda.assign(bound_,0.0);
 
    Gamma_.setZero(L_,L_);
    Pi_.setZero(L_,L_);
    Energy_.setZero(L_);

    Eigen::MatrixXd C_tot(L_,L_);
    Eigen::MatrixXd F_tot(L_,L_);
    Eigen::MatrixXd C_sub(L_/2,L_/2);
    Eigen::MatrixXd F_sub(L_/2,L_/2);
    
    Eigen::MatrixXd temp(L_/2,L_/2);
   
    std::vector<double> spectrum;
     
    Eigen::JacobiSVD<Eigen::MatrixXd> E(H.A_+H.B_, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    Gamma_  =   0.5 *(E.matrixV().transpose()-E.matrixU().transpose());
    Pi_     =   0.5 *(E.matrixV().transpose()+E.matrixU().transpose());
    Energy_ = E.singularValues();
    
    C_tot = 0.25 * (E.matrixV()-E.matrixU())*(E.matrixV().transpose()-E.matrixU().transpose());
    F_tot = 0.25 * (E.matrixV()-E.matrixU())*(E.matrixV().transpose()+E.matrixU().transpose());
    
    C_sub = C_tot.block(0,0,L_/2,L_/2);
    F_sub = F_tot.block(0,0,L_/2,L_/2);
    
    temp=(2*C_sub - Eigen::MatrixXd::Identity(L_/2,L_/2) - 2*F_sub);
    
    Eigen::JacobiSVD<Eigen::MatrixXd> eps(temp);
    
    for (int i=0;i<L_/2;i++) {
        eigen0[i]=eps.singularValues()(i);
        if (eigen0[i]>=1.0)
            eigen0[i]=1.0-epsilon;        
    }

    std::sort(eigen0.begin(),eigen0.begin()+L_/2);

    for (int i=0;i<n_eig_;i++) {
        eigen1[i]=2*atanhl(eigen0[i]);
 
    }
    
    for (int i=0;i<bound_;i++) {
        
        dot_prod=0;
        std::bitset<16> bin(i);
        
        for(int j=0;j<n_eig_;j++) {
            dot_prod += bin[n_eig_-1-j]*eigen1[j];
        }
        
        lambda[i] = exp(-dot_prod);
        Norm += lambda[i];
    }

    for (int i=0;i<bound_;i++) {
        spectrum.push_back(lambda[i]/Norm);
    }
    std::sort(spectrum.begin(),spectrum.begin()+bound_);
    ES_ = spectrum; 
  }

  void Measure() {
    VN_ = 0.0;
    for (int k=0;k<bound_;k++) {
        VN_ += - ES_[k] * (log(ES_[k]));
    }
    SG_ = ES_[bound_-1]-ES_[bound_-2];
    std::cout<<"vN Entropy = " << VN_ <<std::endl;
    std::cout<<"Schmidt gap = " << SG_ << std::endl<<std::endl;
  }


//string entanglement::build_fileName(string & model, double h) {
//
//
//    string name;
//     
//    //name = "data/Ising/L";
//    //name += "data/randomIsing/L";
//    //name += boost::str(boost::format("%d") % L);
//    //name += "/";
//    name += model;
//
//    name += "_L";
//    name += boost::str(boost::format("%d") % L);
//
//    name += "_h";
//    name += boost::str(boost::format("%.4f") % h);
//    
//    //name += ".txt";
//    name += "data.txt";
//    return name;
//}
//
//void entanglement::record(ofstream & file,double h) {
//    
//    //file << h  << "    ";
//    file << SG << "    ";
//    file << VN << "    ";
//    for (int i=0; i<8; i++) {
//
//        file << ES[bound-i-1] << "  ";
//    }
//    file << endl;
//}

};

#endif
