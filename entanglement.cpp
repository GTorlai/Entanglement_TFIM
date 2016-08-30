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
#include "entanglement.h"
#include <boost/format.hpp>

entanglement::entanglement(int L_, char n_eig_) {
    
    L=L_;
    n_eig=n_eig_;
    eig_kept = 8;
    rod = 1;

    bound = myPow(2,n_eig);
    
    ES.assign(bound,0.0);
}


vector<double> entanglement::get_ES(const Hamiltonian & H) {
    
    double Norm=0;
    double epsilon=0.000000000000001;
    double dot_prod;
    
    vector<double> eigen0;
    eigen0.assign(L/2,0.0);
    vector<double> eigen1;
    eigen1.assign(n_eig,0.0);
    vector<char> bin;
    bin.assign(n_eig,0);
    vector<double> lambda;
    lambda.assign(bound,0.0);
 
    Gamma.setZero(L,L);
    Pi.setZero(L,L);
    Energy.setZero(L);

    MatrixXd C_tot(L,L);
    MatrixXd F_tot(L,L);
    MatrixXd C_sub(L/2,L/2);
    MatrixXd F_sub(L/2,L/2);
    
    MatrixXd temp(L/2,L/2);
   
    vector<double> spectrum;
     
    JacobiSVD<MatrixXd> E(H.A+H.B, ComputeFullU | ComputeFullV);
    
    Gamma  =   0.5 *(E.matrixV().transpose()-E.matrixU().transpose());
    Pi     =   0.5 *(E.matrixV().transpose()+E.matrixU().transpose());
    Energy = E.singularValues();
    
    C_tot = 0.25 * (E.matrixV()-E.matrixU())*(E.matrixV().transpose()-E.matrixU().transpose());
    F_tot = 0.25 * (E.matrixV()-E.matrixU())*(E.matrixV().transpose()+E.matrixU().transpose());
    
    C_sub = C_tot.block(0,0,L/2,L/2);
    F_sub = F_tot.block(0,0,L/2,L/2);
    
    temp=(2*C_sub - MatrixXd::Identity(L/2,L/2) - 2*F_sub);
    
    JacobiSVD<MatrixXd> eps(temp);
    
    for (int i=0;i<L/2;i++) {
        eigen0[i]=eps.singularValues()(i);
        if (eigen0[i]>=1.0)
            eigen0[i]=1.0-epsilon;        
    }

    sort(eigen0.begin(),eigen0.begin()+L/2);

    for (int i=0;i<n_eig;i++) {
        eigen1[i]=2*atanhl(eigen0[i]);
 
    }
    
    for (int i=0;i<bound;i++) {
        
        dot_prod=0;
        bitset<8> bin(i);
        
        for(int j=0;j<n_eig;j++) {
            dot_prod += bin[n_eig-1-j]*eigen1[j];
        }
        
        lambda[i] = exp(-dot_prod);
        Norm += lambda[i];
    }

    for (int i=0;i<bound;i++) {
        spectrum.push_back(lambda[i]/Norm);
    }
    sort(spectrum.begin(),spectrum.begin()+bound);
    
    return spectrum;
}

void entanglement::measure(vector<double> & ES) {
    
    VN = 0.0;

    for (int k=0;k<bound;k++) {
        VN += - ES[k] * (log(ES[k]) / log(2));
    }
    
    SG = ES[bound-1]-ES[bound-2];
}


string entanglement::build_fileName(string & model, double h) {


    string name;
     
    //name = "data/Ising/L";
    name += "data/randomIsing/L";
    name += boost::str(boost::format("%d") % L);
    name += "/";
    name += model;

    name += "_L";
    name += boost::str(boost::format("%d") % L);

    name += "_h";
    name += boost::str(boost::format("%.4f") % h);
    
    //name += ".txt";
    name += "data.txt";
    return name;
}

void entanglement::record(ofstream & file,double h) {
    
    //file << h  << "    ";
    file << SG << "    ";
    file << VN << "    ";
    for (int i=0; i<8; i++) {

        file << ES[bound-i-1] << "  ";
    }
    file << endl;
}

int entanglement::myPow (int x, int p) {
    int i = 1;
    for (int j = 1; j <= p; j++)  i *= x;
    return i;
}

