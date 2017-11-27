#include <iostream>
#include <vector>
#include <windef.h>
#include <sys/time.h>
#include <random>

#define SECONDS_TO_NANO_SECONDS 1000000000
#define MICRO_SECONDS_TO_NANO_SECONDS 1000

typedef std::vector<std::vector<int> > matrix;
int beforeStatus, afterStatus;
struct timeval before, after;



/**
 * convert to nano second
 * @param before the time before the operation
 * @param after the time after the operation
 * @return return the converted time
 */
inline double conversionToNanoSecond(timeval &before, timeval &after){
    return ((after.tv_sec-before.tv_sec)*SECONDS_TO_NANO_SECONDS)
           +((after.tv_usec- before.tv_usec)*MICRO_SECONDS_TO_NANO_SECONDS);
}

matrix naiveMatrixMultiplication(matrix &a, matrix &b,unsigned int numRows, unsigned int numCols){

    matrix c;
    c.resize(numCols,std::vector<int>(numRows,0));
    // go row by row
    for(unsigned int row =0;row < a.size(); row++){
        //sum the multiplication for each b column and a row
        for(unsigned int i = 0; i<b[0].size(); i++){
            for(unsigned int j = 0; j< a[0].size(); j++){
                c.at(row).at(i)+= a.at(row).at(j)*b.at(j).at(i);
            }
        }
    }

    return c;
}

matrix addition(matrix& a, matrix &b){
    matrix c;
    c.resize(a[0].size(),std::vector<int>(a.size(),0));

    for(int i =0; i<a.size();i++){
        for(int j=0;j<a[0].size();j++){
            c[i][j] = a[i][j]+b[i][j];
        }
    }
    return c;
}

matrix getBlock( matrix& a, int rowStart,int rowEnd, int colStart, int colEnd,unsigned int n)
{
    matrix block;
    block.resize(n,std::vector<int>(n,0));
    for(int i = rowStart; i < rowEnd; i++ ) {
        for (int j = colStart; j < colEnd; j++) {
            block[i - rowStart][j - colStart] = a[i][j];
        }
    }
    return block;
}

matrix unifyMatrix(matrix &c11, matrix &c12, matrix &c21, matrix &c22, unsigned int n){
    matrix c;
    c.resize(n,std::vector<int>(n,0));
    int half = n/2;
    for (int row=0;row < n;row++ ){
        for(int col = 0; col < n; col++){
            if (row < half){
                if(col< half){
                    c[row][col]= c11[row][col];
                }else{
                    c[row][col] = c12[row][col-half];
                }
            }else{
                if(col< half){
                    c[row][col] = c21[row-half][col];
                }else{
                    c[row][col] = c22[row-half][col-half];
                }
            }
        }
    }
    return c;
}

matrix recursiveMatrixMultiplication(matrix &a, matrix &b, unsigned int n){

    //base case
    if(n ==1){
        matrix c;
        c.resize(n,std::vector<int>(n,0));
        c[0][0]= a[0][0]*b[0][0];
        return c;
    }
    // get the blocks
    matrix a11 = getBlock(a,0, n/2, 0,n/2,n/2);
    matrix b11 = getBlock(b,0, n/2, 0,n/2,n/2);
    matrix a12 = getBlock(a,0, n/2, n/2,n,n/2);
    matrix b12 = getBlock(b,0, n/2, n/2,n,n/2);
    matrix a21 = getBlock(a,n/2, n, 0,n/2,n/2);
    matrix b21 = getBlock(b,n/2, n, 0,n/2,n/2);
    matrix a22 = getBlock(a,n/2, n, n/2,n,n/2);
    matrix b22 = getBlock(b,n/2, n, n/2,n,n/2);

    matrix part1 = recursiveMatrixMultiplication(a11, b11, n/2 );
    matrix part2 = recursiveMatrixMultiplication(a12,b21,n/2);
    matrix part3 = recursiveMatrixMultiplication(a11, b12, n/2 );
    matrix part4 = recursiveMatrixMultiplication(a12,b22,n/2);
    matrix part5 = recursiveMatrixMultiplication(a21, b11, n/2 );
    matrix part6 = recursiveMatrixMultiplication(a22,b21,n/2);
    matrix part7 = recursiveMatrixMultiplication(a21, b12, n/2 );
    matrix part8 = recursiveMatrixMultiplication(a22,b22,n/2);

    matrix c11 = addition(part1,part2);
    matrix c12 = addition(part3,part4);
    matrix c21 = addition(part5,part6);
    matrix c22 = addition(part7,part8);

    return unifyMatrix(c11, c12, c21, c22, n);

}

matrix blockMatrixMultiplication(matrix &a, matrix &b, int blockSize, int n){

    matrix c;
    c.resize(n,std::vector<int>(n,0));

    int sum = 0 ;
    for(int kk=0;  kk <n; kk += blockSize) {
        for (int jj = 0; jj < n; jj += blockSize) {
            for (int i = 0; i < n; i++) {
                for (int j = jj; j < jj + blockSize; j++) {
                    sum = c[i][j];
                    for (int k = kk; k < kk + blockSize; k++) {
                        sum += a[i][k] * b[k][j];
                    }
                    c[i][j] = sum;
                }
            }
        }
    }

    return c;
}

int main(){

    unsigned int n=10000;
    matrix a,b,c;
    a.resize(n,std::vector<int>(n,0));
    b.resize(n,std::vector<int>(n,0));
    c.resize(n,std::vector<int>(n,0));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 10);

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            a[i][j]= dis(gen);
        }
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            b[i][j]= dis(gen);
        }
    }

    beforeStatus = gettimeofday(&before,NULL);
    c = naiveMatrixMultiplication(a,b,4,4);
    afterStatus= gettimeofday(&after,NULL);
    std::cout<<conversionToNanoSecond(before,after)<<std::endl;

    beforeStatus = gettimeofday(&before,NULL);
    c = recursiveMatrixMultiplication(a,b,4);
    afterStatus= gettimeofday(&after,NULL);
    std::cout<<conversionToNanoSecond(before,after)<<std::endl;

    beforeStatus = gettimeofday(&before,NULL);
    c = blockMatrixMultiplication(a,b,2,4);
    afterStatus= gettimeofday(&after,NULL);
    std::cout<<conversionToNanoSecond(before,after)<<std::endl;

    return 0;

}
