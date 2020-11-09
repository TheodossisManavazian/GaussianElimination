#include <iostream>
#include "Eigen/Core"
#include "Eigen/LU"


using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void backSubstitution(MatrixXd&);
void swapRows(MatrixXd&,int,int);

int main(){

    int rows, columns;
    cout << "Enter rows and columns" << endl;
    cin >> rows;
    cin >> columns;

    MatrixXd m(rows,columns);

    cout << "Enter equations for augmented matrix" << endl;

    double num;
    for(int i = 0; i < rows; i ++){
        for(int j = 0; j < columns; j++){
            cin >> num;
            m(i,j) = num;
        }
    }

    MatrixXd A = m.block(0,0,m.rows(), m.cols()-1);

    // ensures that the first row has a leading 1
    //////////////////////////////////////
    bool flag = false;
    for(int i = 0; i < m.rows(); i ++){
        if (m(i, 0) == 1){
            swapRows(m,0,i);
            flag = true;
        }
    }
    if (!flag){
        double divi = m(0,0);
        for(int col = 0; col < m.cols(); col++){
            m(0,col) /= divi;
        }
    }
    //////////////////////////////////////

    // performs elimination, makes sure the matrix is upperr triangular
    //////////////////////////////////////
    for(int col = 0; col < m.cols()-1; col ++){
        // go through the matrix, make sure that the first leading value that isnt a 0, is a 1.
        if (m(col, col) != 1){
            m.row(col) /= m(col, col);
        }
        // set the temp vector = to the row before the one were trying to evaluate 
        VectorXd temp = m.row(col);
        for(int row = col; row < m.rows()-1; row ++){
            double mult = m(row+1,col) * -1;
            temp *= mult; 
            m.row(row+1) += temp;
            // reset the temp vector so the next row gets evaluated correctly
            temp = m.row(col);
        }   
    }
    //////////////////////////////////////

    backSubstitution(m);


    //////////////////////////////////////
    // compute B vector by plugging in the x values back into their original equations
    VectorXd x = m.col(columns-1);
    VectorXd b(m.rows());
    for(int row = 0; row < A.rows(); row ++){
        for(int col = 0; col < A.cols(); col++){
            A(row, col) *= x(col);
        }
        b(row) = A.row(row).sum();
    }
    //////////////////////////////////////

    cout << "M" << endl << m << endl;
    cout << "B" << endl << b << endl;


    return 0;
}
void backSubstitution(MatrixXd& m){
    // start at end of the matrix and subsitute the b vector into the matrix so we can get 1's along the diangonal and 0s in every other place
    for (int i = m.rows()-1; i >= 0; i--){
        for (int j = i+1; j < m.cols()-1; j++){
            m(i,j) *= m(j,m.cols()-1); 
            m(i,m.cols()-1) -= m(i,j);
            m(i,j) = 0;
        }
        // divide row i by m(i, i) so we can eliminate that variable
        m.row(i) /= m(i,i); 
    }
}

void swapRows(MatrixXd& m, int i, int j){
    // create a temp variable so we can swap 2 rows. pass matrix by reference so oour changes stick outside the function
    VectorXd temp = m.row(i);
    m.row(i) = m.row(j);
    m.row(j) = temp;

}