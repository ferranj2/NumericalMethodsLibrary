#include <iostream>//Standard Library for input/output
#include <iomanip>//Standrd Library for I/0 manipulation (need setw() function)
#include <cstdlib>//Allows use of the rand() function
#include <ctime>//Allows querying numbers from the clock.
#include <string>
#include <limits>
using namespace std;
//TO-DO:
//Create a class specifically for solving linear systems.
//Make an LU factorization function?
//Cholesky decomposition?

class matrix{
//CLASS MEMBERS
public:
  int R;//Rows of the matrix
  int C;//Columns of the matrix
  int N;//Number of elements in the matrix
  double *DATA;//Matrix coefficients

//CONSTRUCTORS
public:
  //Default Constructor
  matrix(){
    this->R=1;//Default number of rows is 1.
    this->C=1;//Default number of columns is 1.
    this->N=1;//Default number of elements is 1.
    this->DATA = new double[1];//Default data is a 1x1 0-matrix
  }//end constructor

  //Constructor for simple memory allocation
  matrix(int r, int c){
    this->R=r;
    this->C=c;
    this->N=r*c;
    this->DATA = new double[N];
  }//end constructor

  //Constructor for predefined values (i.e., MATLAB's zeros and ones)
  matrix(int r, int c, double val){
    this->R=r;
    this->C=c;
    this->N=r*c;
    this->DATA=new double[N];
    for(int i=0; i<N;i++){
      DATA[i]=val;//Assign the same value to all matrix entries.
    }//end "i" loop
  }//end constructor

  //Constructor for special matrices
  matrix(int r, int c, string a){
    this->R=r;
    this->C=c;
    this->N=r*c;
    this->DATA=new double[N];

    //Random integer matrix
    if(a=="randi"){
      srand(time(NULL));//Use clock time to generate random seeds.
      for(int i=0; i<N;i++){
        DATA[i] = rand() % 10;//Assign pseudo-random integers (0-9) to vector B
      }//end "i" loop
    }//end if

    //Identity matrix
    else if(a=="eye"){
      int pivot;
      for(int i=0;i<R;i++){
        pivot = i*(1+C);//Index of the pivot element
        *(DATA+pivot) = 1.0;//Set pivot elements equal to 1.
        for(int j =pivot+1;j<pivot+C;j++){//All elements until next pivot...
          *(DATA+j) = 0.0;//Set them to zero
        }//end "j" loop
      }//end "i" loop
    }//end else if

  }//end constructor

public://METHODS
  //Compute the trace of a matrix
  double trace(void){
    if(R==C){
      double sum = 0;//Running sum variable for the trace.
      for(int i=0;i<C;i++){
        sum += DATA[i*(C+1)];//Reference diagonal elements
      }//end "i" loop
      return sum;
    }//end if
    else{
      cout<<"WARNING: Trace undefined for non-square matrix!"<<endl;
      return std::numeric_limits<double>::quiet_NaN();
    }//end else
  }
  //Compute the eigenvalues of a square matrix.
  double eigen(void){
    if(R==C){
      cout<<"WIP"<<endl;
    }//end if
    else{
      cout<<"WARNING: Eigenvalues undefined for non-square matrix!"<<endl;
      return std::numeric_limits<double>::quiet_NaN();
    }//end else
  }
  //Naive forward elimination
  void fel(void){
    int swap = 0;//Counter variable for number of partial pivots.
    for(int i=0; i<(R-1); i++){//Scan all rows expcet the very last
      for(int k=i+1; k<R; k++){
        //Conditional Partial pivoting.
        if(*(DATA+i*C+i) == 0){
          double temp[C];//temporary vector to store row with pivot  = 0;
          for(int l=0; l<R;l++){
            temp[l] = *(DATA+i*C+l); //Copy elements of 0-pivot row to temporary vector
          }//end "l" loop
          for(int l=0; l<R;l++){
            *(DATA+i*C+l) = *(DATA+N-C+l);//Overwrite current row with a copy of the last row.
          }//end "l" loop
          for(int l=0; l<R;l++){//end constructor
            *(DATA+N-C+l) = temp[l];//Overwrite last row with data stored in the temporary vector.
          }//end "l" loop
          swap += 1;
        }//end if statement
        double ratio = *(DATA+k*C+i)/(*(DATA+i*C+i));
        for(int m=1; m<C; m++){//loop for eliminating the remainder of the row.
          *(DATA+k*C+m) = *(DATA+k*C+m) - ratio*(*(DATA+i*C+m));
          *(DATA+k*C+i) = 0.0;
        }//end "m" loop
      }//end "k" statement
    }//end "i" loop
  }
  //Return number of elements in matrix
  int nel(void){
    return N;
  }
  //Return number of rows
  int cols(void){
    return C;
  }
  //Return number of columns
  int rows(void){
    return R;
  }
  //Print matrix to console
  void display(void){
    for(int i=0; i<R; i++){//Scan the rows.
      cout <<"[";//Announce beginning of row by using a square bracket.
      for(int j=0; j<C; j++){//scan the columns.
        cout<< setw(10);//This is a reasonable width for outputting matrix elements.
        cout<< *(DATA+i*C+j) <<' ';//Output values of the matrix without jumping lines.
      }//end "j" for loop
      cout<<"]"<<endl;//Jump to the next line.
    }//end "i" for loop
  }

  //OPERATOR OVERLOADING

  //Matrix addition with a constant
  matrix operator+(double B){
    matrix ApB(this->R,this->C);
    for(int i=0;i<N;i++){
      ApB.DATA[i] = *(DATA+i) + B;
    }//end "i" loop
    return ApB;
  }//end matA + B

  //Matrix multiplication with a constant
  matrix operator*(double B){
    matrix AtB(this->R,this->C);
    for(int i=0;i<N;i++){
      AtB.DATA[i] = (*(DATA+i))*B;
    }//end "i" loop
    return AtB;
  }// matA*B

  //Matrix addition with another matrix
  matrix operator+(const matrix& B){
    if (R == B.R && C == B.C){
      matrix ApB(B.R,B.C);
      for(int i=0; i<N; i++){
        ApB.DATA[i] = *(DATA+i) + B.DATA[i];//Add elements of both matrices together.
      }//end "i" lloop
      return ApB;
    }//end if
    else{
      cout<<"ERROR: matrices are of unequal size"<<endl;
      return matrix();
    }//end else
  }//end matA +_matB

  //Matrix multiplication with another matrix
  matrix operator*(const matrix& B){
    if(C==B.R){
      matrix AtB(R,B.C);//Allocate space for the output matrix
      double sum = 0.0;//Running sum variable for elementwise computation
      for(int i=0;i<R;i++){//For all rows of the left matrix
        for(int j=0;j<B.C;j++){//For all columns of the right matrix
          for(int k=0;k<B.R;k++){//For all elements of the jth column's elements
            sum += (*(DATA+i*C+k))*B.DATA[j+k*C];//
          }//end "k" loop
          AtB.DATA[j+i*C] = sum;//Assign element to the output matrix.
          sum = 0.0;//Reset running sum variable
        }//end "j" loop
      }//end "i" loop
      return AtB;
    }//end if
    else{
      cout<<"ERROR: matrices are of unequal size"<<endl;
      return matrix();
    }//end else
  }//end matA*matB


};//end matrix CLASS

/*
//Test program
int main(){

  matrix A(5,5,"rand");
  A.display();
  matrix B(5,5,"rand");
  int rows = A.rows();
  cout<<"rows of A "<<rows<<endl;
  int cols = A.cols();
  cout<<"columns of A "<<cols<<endl;
  int elements = A.nel();
  cout<<"# elements of A "<<elements<<endl;
  double trac = A.trace();
  cout<<"trace of A "<<trac<<endl;
  A.fel();
  A.display();
  cout<<"matrix B is "<<endl;
  B.display();
  matrix ApB = A+B;
  cout<<"matrix B plus A.fel() "<<endl;
  ApB.display();
matrix A(3,3,"randi");
cout<<"Matrix A:"<<endl;
A.display();
matrix B(3,3,"randi");
cout<<"Matrix B:"<<endl;
B.display();
matrix C = A*B;
cout<<"Matrix C = A*B:"<<endl;
C.display();
cout<<"Matrix D = C + 1:"<<endl;
matrix D = C+1.0;
D.display();
cout<<"Matrix E = D*3:"<<endl;
matrix E = D*3.0;
E.display();
cout<<"Matrix I = Identity:"<<endl;
matrix I(3,3,"eye");
I.display();
  return 0;
}//end main program
*/
