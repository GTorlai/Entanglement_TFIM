#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <stdio.h>
#include <stdlib.h>

// Parameter Class
class Parameters{

public:
  
  int L_; 
  double h_;
  int seed_;          
  int eigen_keep_;
  std::string boundary_conditions_; 

  Parameters() {
    L_ = 8;
    eigen_keep_=4;
    seed_ = 16382;
    boundary_conditions_ = "open";
  }
    
  // Read parameters from the command line
  void ReadParameters(int argc,char** argv){
    std::string flag;
    
    flag = "-L";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) L_=atoi(argv[i+1]);
    }
    flag = "-h";
    for(int i=1;i<argc;i++){
      if(flag==argv[i]) h_=double(atof(argv[i+1]));
    }
  }
};
#endif
