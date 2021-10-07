#include "matrix.h"//Need to include my custom matrix class for this.
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
int main(int argc, char *argv[]){
  int M, N;
  double tolerance, ut, ub, ur, ul;
  //Grid points
  M = atoi(argv[1]);//1st input is the number of rows.
  N = atoi(argv[2]);//2nd input is the number of columns.
  //Dirichlet boundary conditions.
  ut = std::stod(argv[3]);//3rd input is the top temperature.
  ul = std::stod(argv[4]);//4th input is the left temperature.
  ur = std::stod(argv[5]);//5th input is the right temperature.
  ub = std::stod(argv[6]);//6th input is the bottom temperature.
  //Convergence criterion
  tolerance = std::stod(argv[7]);//7th input is the tolerance.
  std::cout<<"INPUT PARAMETERS: "<<std::endl;
  std::cout<<"================="<<std::endl;
  std::cout<<"# Rows       = "<<M<<std::endl;
  std::cout<<"# Columns    = "<<N<<std::endl;
  std::cout<<"Top Temp.    = "<<ut<<std::endl;
  std::cout<<"Bottom Temp. = "<<ub<<std::endl;
  std::cout<<"Right Temp.  = "<<ur<<std::endl;
  std::cout<<"Left Temp.   = "<<ul<<std::endl;
  std::cout<<"Tolerance    = "<<tolerance<<"\n"<<std::endl;

  //DOMAIN INITIALIZATION
  double u_avg = (ut*(N-2)+ub*N+ur*(M-1)+ul*(M-1))/(2*(M+N)-4);//Average temperature of  edges.
  matrix U(M,N,u_avg);//Initialize matrix with u_avg as default value of all entries.
  matrix P(M,N);//Declare a second matrix to store solution at "n-1" time step
  //Initialize the top edge
  for(int i=0; i<N;i++){
      U.DATA[i] = ut;
  }//end "i" loop
  //Initialize both the left and right edges
  for(int i=0; i<M;i++){
    U.DATA[i*N] = ul;//Assign values to the left edge.
    U.DATA[(i+1)*N - 1] = ur;//Assign values to the right edge.
  }//end "i" loop
  //Initialize the bottom edge
  for(int i=0; i<N;i++){
    U.DATA[N*(M-1) + i] = ub;
  }//end "i" loop

  //SOLUTION PARAMETERS
  int iter=0; //Iteraction counter.

  //Initialization
  int index;//Variable to loop through matrix indices.
  double diff;//Variable to compute node-wise residuals.
  double max_diff=tolerance*1.01;//store maximum residual.
  int power = 0;//Powers of 2 worth of iterations spanned.
  int next = 1;//Will equal 2^(power) inside the loop.

  std::cout<<"NOTES:"<<std::endl;
  std::cout<<"================="<<std::endl;
  std::cout<<"Maximum Iterations are hardcoded at 1000!"<<std::endl;
  std::cout<<"Temperature inputs are of type double!\n"<<std::endl;
  std::cout<<"Begin Iterations!"<<std::endl;
  std::cout<<"================="<<std::endl;
  while(max_diff>tolerance && iter <= 1000){
    iter += 1;//Update iteration counter.

    if(max_diff > tolerance){
      max_diff = 0;//Reset the max difference variable
    }//end else
    P = U;//Copy contents of updated matrix to old one.

    //Update the U-matrix
    for(int j=1; j<N-1;j++){//Column looper
      for(int k=1; k<M-1;k++){//Row looper
        index = j+k*N;//Linear index of internal matrix element.
        U.DATA[index] = (P.DATA[index+1] + P.DATA[index-1] + P.DATA[index+N] + P.DATA[index-N])/4;//Update the value
        diff = abs(U.DATA[index]-P.DATA[index]);//Compute residual
        if(diff>max_diff){//If this residual greater than current max
          max_diff = diff;//Label new residual as maximum.
        }//end if
      }//end "k" loop
    }//end "j" loop

    if(iter == next){
      power += 1;//Update the number of powers worth of iterations spanned.
      next *= 2;//Compute the next iteration at which to report.
      //Report residuals
      std::cout<<"Iteration #"<<setw(5)<<iter<<" max_diff = "<<setw(15)<<max_diff<<std::endl;
    }//end if

  }//end while loop
  if(tolerance > max_diff){
  std::cout<<"!!!CONVERGENCE CRITERIA MET!!!"<<std::endl;
}// end if
  else{
    std::cout<<"!!!WARNING: MAXIMUM ITERATIONS REACHED!!!"<<std::endl;
  }// end else
  std::cout<<"Iteration #"<<setw(5)<<iter<<" max_diff = "<<setw(15)<<max_diff<<std::endl;

  return 0;
}//end main()
