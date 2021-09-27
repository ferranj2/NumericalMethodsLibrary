#include "matrix.h"//Need to include my custom matrix class for this.
#include <cmath>
#include <iostream>
int main(){

  //SIMULATION PARAMETERS:
  //Grid points
  int M = 5;//Number of rows.
  int N = 5;//Number of columns.
  //Dirichlet boundary conditions.
  double ut = 100.0;//Temperature of the top row.
  double ub = 400.0;//Temperature of the bottom row.
  double ur = 300.0;//Temperature of the right column.
  double ul = 200.0;//Temperature of the left column.

  //DOMAIN INITIALIZATION
  double u_avg = (ut*(N-2)+ub*N+ur*(M-1)+ul*(M-1))/(2*(M+N)-4);//Average temperature of  edges.
  matrix U(M,N,u_avg);//Initialize matrix with u_avg as default value of all entries.
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
  std::cout<<"Initial matrix"<<std::endl;
  U.display();

  //SOLUTION PARAMETERS
  int n = 4;//Number of iterations as a power of 2.
  int max_iter = pow(2,n);//Actual number of iterations.
  int iter=0;//Maximum allowed iterations.
  double tolerance = 0.0001;//Smallest tolerable residual

  //std::cout<<iter<<std::endl;

  int index;//Variable to loop through matrix indices.
  double diff;//Variable to compute node-wise residuals.
  double old_val;//Variable to store results of previous iteration.
  double max_diff=tolerance*1.01;//store maximum residual.
  double min_diff=abs(ut-ub);//store minimum residual.
  int power = 0;//Powers of 2 worth of iterations spanned.
  int next = 1;//Will equal 2^(power) inside the loop.

  while(max_diff>tolerance && iter <= max_iter){
  //for(int i=0;i<max_iter;i++){//Iteration looper (Change to while loop)

    iter += 1;//Update iteration counter.

    //Update the U-matrix
    for(int j=1; j<N-1;j++){//Column looper
      for(int k=1; k<M-1;k++){//Row looper
        index = j+k*N;//Linear index of internal matrix element.
        old_val = U.DATA[index];//Copy current value.
        U.DATA[index] = (U.DATA[index+1] + U.DATA[index-1] + U.DATA[index+N] + U.DATA[index-N])/4;//Update the value
        diff = abs(U.DATA[index]-old_val);//Compute residual
        if(diff>max_diff){//If this residual greater than current max
          max_diff = diff;//Label new residual as maximum.
        }//end if
        else if(diff<min_diff){//If this residual less than current min
          min_diff = diff;//Label new residual as minimum
        }//end else if
      }//end "k" loop
    }//end "j" loop

    if(iter == next){
      power += 1;//Update the number of powers worht of iterations spanned.
      next *= 2;//Compute the next iteration at which to report.
      //Report residuals
      std::cout<<"Iteration #"<<iter+1<<" max_diff = "<<max_diff<<" min_diff = "<<min_diff<<std::endl;
      //U.display();
    }//end if

    max_diff = tolerance*1.01;//Reset the max difference variable
    min_diff = abs(ut-ub);//Reset the min difference variable
  }//end while loop

  std::cout<<"Iteration #"<<iter+1<<" max_diff = "<<max_diff<<" min_diff = "<<min_diff<<std::endl;

  return 0;
}//end main()
