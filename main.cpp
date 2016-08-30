/*
 *  2D_Ising_MC.cpp
 *
 *  Monte Carlo Simulation of a 2D classical Ising model
 *
 *  Created by Giacomo on 10/10/2014.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *bla
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "entanglement.cpp"
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    
    int L;
    double h;
    int sweeps;
   
    //string model = argv[1];

    for (int i=1;i<argc;i++) {
        if (strcmp(argv[i],"--L") == 0) {
            L = atoi(argv[i+1]);
            cout << "L is " << L << endl;
        }
        else if (strcmp(argv[i], "--h") == 0) {
            float hF = atof(argv[i+1]);
            h = double(hF);
            cout << "h is " << h << endl;
        }
        else if (strcmp(argv[i], "--sw") == 0) {
            sweeps = atoi(argv[i+1]);
        }
    }

    
    Hamiltonian H(L,"Open");
    entanglement E(L,8);
    //string model = "Collapse_Ising";
    string model = "randomIsing";
    string fileName = E.build_fileName(model,h);
    
    ofstream fout(fileName);
    
    //Homogeneous
    //H.build_Ising(h,1.0);
    //E.ES = E.get_ES(H);
    //E.measure(E.ES);
    //E.record(fout,h);
    
    //Random
    MTRand random;
    for (int k=0; k<sweeps; k++) {
        H.build_RandomIsing(h,1.0,random);
        E.ES = E.get_ES(H);
        E.measure(E.ES);
        E.record(fout,h); 
    } 
//    
//    Par parameter;
//    
//    MTRand random;
//    
//    switch(alg) {
//            
//        case 1: {
//            model_Hamiltonian Hamilt(L);
//            
//            entanglement_1d ESFF(L,parameter.save_eig,parameter.n_eig);
//            
//            measure experiment(parameter);
//            
//            experiment.create_file(L,sweeps,g_S,m_flag);
//            ofstream file_data(experiment.ss.c_str());
//            
//            for(int k=0;k<sweeps;k++) {
//                
//                switch(m_flag) {
//                    case 0:
//                        Hamilt.TFIM(g);
//                        break;
//                    case 1:
//                        Hamilt.RBRFIM(g,random);
//                        break;
//                }
//                
//                ESFF.entanglement_spectrum(Hamilt);
//                experiment.update(ESFF);
//                
//            }
//            
//            experiment.record(g,sweeps,file_data);
//            file_data.close();
//        }
//            break;
//            
//        case 2: {
//            
//            model_Hamiltonian Hamilt0(L);
//            model_Hamiltonian Hamilt1(L);
//            
//            entanglement_1d ESFF0(L,parameter.save_eig,parameter.n_eig);
//            entanglement_1d ESFF1(L,parameter.save_eig,parameter.n_eig);
//            
//            dynamics quench(L,parameter.n_eig);
//            
//            measure experiment(parameter);
//            
//            experiment.create_file(L,sweeps,t_S,m_flag);
//            ofstream file_data(experiment.ss.c_str());
//            
//            for(int k=0;k<sweeps;k++) {
//                
//                switch(m_flag) {
//                    case 0:
//                        Hamilt0.TFIM(g0);
//                        Hamilt1.TFIM(g1);
//                        break;
//                    case 1:
//                        Hamilt0.RBRFIM(g0,random);
//                        Hamilt1.RBRFIM(g1,random);
//                        break;
//                }
//                
//                ESFF0.entanglement_spectrum(Hamilt0);
//                ESFF1.entanglement_spectrum(Hamilt1);
//                
//                quench.instant_quench(ESFF0,ESFF1,t);
//                
//                experiment.update_dyn(quench);
//
//                
//            }
//            experiment.record(t,sweeps,file_data);
//
//            file_data.close();
//
//        }
//            break;
//    }
    
}
