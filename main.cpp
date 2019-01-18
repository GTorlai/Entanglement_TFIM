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
#include "entanglement.hpp"
#include <cmath>
#include "parameters.hpp"

int main(int argc, char *argv[]) {
  
  Parameters pars;
  pars.ReadParameters(argc,argv);

  Hamiltonian H(pars.L_,pars.boundary_conditions_);
  Entanglement entanglement(pars.L_,pars.n_eig_,pars.eig_keep_);
  std::string model = "ising";
  H.Ising(pars.h_);
 
  entanglement.ComputeGreensFunction(H);
  for (int l=1;l<pars.L_;l++){
    entanglement.ComputeEntanglementSpectrum(l);
    entanglement.Measure();
    std::cout<<"Entanglement entropy at bond "<< l <<" : S = " << std::setprecision(10)<<entanglement.VN_ <<std::endl;
  }
  //entanglement.ComputeEntanglementSpectrum(1);
  //entanglement.Measure();
  //std::cout<<"Entanglement entropy = " << setprecision(10)<<entanglement.VN_ <<std::endl;
  //entanglement.ComputeEntanglementSpectrum(2);
  //entanglement.Measure();
  //std::cout<<"Entanglement entropy = " << setprecision(10)<<entanglement.VN_ <<std::endl;
  //entanglement.ComputeEntanglementSpectrum(3);
  //entanglement.Measure();
  //std::cout<<"Entanglement entropy = " << setprecision(10)<<entanglement.VN_ <<std::endl;
  //entanglement.ComputeEntanglementSpectrum(4);
  //entanglement.Measure();
  //std::cout<<"Entanglement entropy = " << setprecision(10)<<entanglement.VN_ <<std::endl;

    //Hamiltonian H(L,boundary_conditions);
    //entanglement E(L,8);
    ////string model = "Collapse_Ising";
    //string model = "ising";
    //string fileName = E.build_fileName(model,h);
    //
    //ofstream fout(fileName);
    //
    ////Homogeneous
    //H.build_Ising(h,1.0);
    //E.ES = E.get_ES(H);
    //E.measure(E.ES);
    //E.record(fout,h);
    
    //Random
    //MTRand random;
    //for (int k=0; k<sweeps; k++) {
    //    H.build_RandomIsing(h,1.0,random);
    //    E.ES = E.get_ES(H);
    //    E.measure(E.ES);
    //    E.record(fout,h); 
    //} 
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
