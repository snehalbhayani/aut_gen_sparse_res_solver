// Reduce the monomial basis and the corresponding coefficient matrix
// by finding the nonzero independent entries in the coefficient matrix.
// i - spans the columns
// j - spans the rows
#include "mex.h"
#include "math.h"

bool checkIsMatrix0(double *matrixToTest, int m, int n) {
    // to be used to check if a given matrix has all elements 0
    int i=0,j=0, isMatrix0=true, ofs=0;
    
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            if(*(matrixToTest + ofs) != 0.0) {
                isMatrix0 = false;
                break;
            }
            ofs++;
        }
        if(!isMatrix0 ) {
            break;
        }
    }
    return isMatrix0;
}

void outerProduct(double *M, double *N, double *outerProductResult, int sizeOfEachVector,  int degreeExponentOfHiddenVar, double *zeroCoeffMatrixIndices, double coeffToRemove) {
//     double *outerProductResult = (double*)malloc(sizeOfEachVector*sizeOfEachVector*sizeof(double));
    int resultOfs = 0,i,j, offsetToAddTo = sizeOfEachVector * sizeOfEachVector * degreeExponentOfHiddenVar;
//      mexPrintf("\n .....RowCnt %d ........ offsetToAddTo %d ..... degreeExpOfHiddenVar %d",sizeOfEachVector,offsetToAddTo,degreeExponentOfHiddenVar );
    for(i=0;i<sizeOfEachVector;i++) {
        for(j=0;j<sizeOfEachVector;j++) {
            *(outerProductResult + offsetToAddTo + resultOfs) -=   (*(M + j)) * (*(N + i)) /coeffToRemove ;
//              mexPrintf("\n ..... %f ... %d .... %d",*(outerProductResult + offsetToAddTo + resultOfs), i, j);
//             mexPrintf("\n ..... %f ... %d .... %d",*(M + j), i, j);
            resultOfs++;
        }
    }
    if(checkIsMatrix0(outerProductResult + offsetToAddTo, sizeOfEachVector,sizeOfEachVector)){
        *(zeroCoeffMatrixIndices + degreeExponentOfHiddenVar) = 0;
    } else {
        *(zeroCoeffMatrixIndices + degreeExponentOfHiddenVar) = 1;
    }
}

void normalizeMatrix(int m, int n, double *normalizedMatrix, double *matrixToNormalize) {
    int i,j,coefOfs=0;
    float sum=0;
    for(i=0;i<n;i++) {
        for (j=0;j<m;j++) {
            sum += pow(*(matrixToNormalize + coefOfs),2);
            coefOfs++;
        }
    }
    
    // Getting the magnitude
    sum = pow(sum, 0.5);
    coefOfs = 0;
    for(i=0;i<n;i++) {
        for (j=0;j<m;j++) {
            *(normalizedMatrix + coefOfs) = *(matrixToNormalize + coefOfs) / sum;
            coefOfs ++;
        }
    }
    
}

int findSortedInd(double *ptr, int size, int indx) {
    int sortedIndOfs = 0,i;
    
    for(i=0;i<size;i++) {
        if (*(ptr+i) <= indx) {
            sortedIndOfs++;
        }
    }
    
    return sortedIndOfs;
    
}
void mexFunction(int nlhs , mxArray *plhs[] , int nrhs , const mxArray *prhs[]){
    
    int h,i,j,i1,j1,k,l,redInd,totalNoOfPossiblyReducedCoeffMatrices;
    int monomialPosI, monomialPosJ, newMonPosI, newMonPosJ, coeffMatrixEle;
    int coeffOfs = 0,redCoeffOfs =0, tempCoeffOfs=0;
    int rowCntCoeffMatrix = mxGetM(prhs[0]);
    int colCntCoeffMatrix = mxGetN(prhs[0]);
    int basisRCnt = mxGetM(prhs[1]);
    int basisCCnt = mxGetN(prhs[1]);
    int noOfCoeffsMatrices = colCntCoeffMatrix/rowCntCoeffMatrix, origNoOfCoeffsMatrices = noOfCoeffsMatrices;
    int sizeOfEachCoeffMatrix = rowCntCoeffMatrix * rowCntCoeffMatrix;
    int reducedMatrixRowCnt = rowCntCoeffMatrix - 1;

    // Various pointers we shall be using
    double *coeffMatrix = mxGetPr(prhs[0]);
    double *basis = mxGetPr(prhs[1]);
    // Parameter that denotes, the amount of reduction to be attempted for the monomial basis;
    double *basisRedSize = mxGetPr(prhs[2]);
    // parameter that denotes the maximum numeric tolerance  to be used while selecting an element to be inverted and reduced.
    double *thresholdForInvertibleElement = mxGetPr(prhs[3]);
    double *indicesToSkip = mxGetPr(prhs[4]);
    double *reducedCoeffMatrix, *reducedMonomialBasis, *zeroCoeffMatrixIndices, *rowsToRemove, *colsToRemove;
    double *indxTrackMat;
    double *M, *N, *origCoeffMatrix;
    double coeffToRemove;
    int redSize = *(basisRedSize), rowToRemoveIndx = 0, origCsSize = rowCntCoeffMatrix, coeffOfsOfs;
    
    // We need to keep a track of the original coefficient matrix, to be used 
    // for tracking the matrix elements to be removed in a cumulative fashion.
    origCoeffMatrix = (double*)malloc(rowCntCoeffMatrix*colCntCoeffMatrix*sizeof(double));
    memcpy(origCoeffMatrix, coeffMatrix, rowCntCoeffMatrix*colCntCoeffMatrix*sizeof(double));
    
    plhs[3] = mxCreateDoubleMatrix(redSize, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(redSize, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(origCsSize, 2 * origCsSize, mxREAL);
    rowsToRemove = mxGetPr(plhs[3]);
    colsToRemove = mxGetPr(plhs[4]);
    newMonPosI = -1;
    newMonPosJ = -1;
    
//     indxTrackMat = (int*)malloc(2 * origCsSize*origCsSize*sizeof(int));
    indxTrackMat = mxGetPr(plhs[5]);
    coeffOfs = 0;
    for(i=0;i<origCsSize;i++) {
        for(j=0;j<origCsSize;j++) {
            *(indxTrackMat + coeffOfs) = j;
            *(indxTrackMat + origCsSize*origCsSize+ coeffOfs) = i;
            coeffOfs++;
        }
    }
    
    coeffOfs = 0;
    for(redInd=0;redInd<redSize;redInd++) {
        mexPrintf(" reduction index :: %d \n", redInd);
        
        M = (double*)malloc(reducedMatrixRowCnt*sizeof(double));
        N = (double*)malloc(reducedMatrixRowCnt*sizeof(double));
        
        // variables related to output params
        totalNoOfPossiblyReducedCoeffMatrices = noOfCoeffsMatrices*2 - 1;
        plhs[0] = mxCreateDoubleMatrix(reducedMatrixRowCnt, reducedMatrixRowCnt*totalNoOfPossiblyReducedCoeffMatrices, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(basisRCnt, basisCCnt - 1, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(totalNoOfPossiblyReducedCoeffMatrices, 1, mxREAL);
        reducedCoeffMatrix = mxGetPr(plhs[0]);
        reducedMonomialBasis = mxGetPr(plhs[1]);
        zeroCoeffMatrixIndices = mxGetPr(plhs[2]);
        // We need to iterate through the coefficient matrix and find the first non-zero element that we can find and use it to
        // reduce the size of the coefficient matrix.
        coeffOfs = sizeOfEachCoeffMatrix-1;
        monomialPosI=-1;
        coeffToRemove =1.0;
        // We should not scan all columns. As the columns that are mapped to those monomials which are to be used
        // in estimating values of un-hiden vars, should not be removed during a reduction step.
        
        for(i=rowCntCoeffMatrix-1;i>*(indicesToSkip);i--) {
            for(j=rowCntCoeffMatrix-1;j>=0;j--) {
    
                // Checking for whether the elements at corresponding positions for multiples f higher powers of hidden vars 
                // are 0.
                l = 1;
                for(k=1;k<noOfCoeffsMatrices;k++) {                    
                    if(*(coeffMatrix + coeffOfs + k * sizeOfEachCoeffMatrix) == 0.0) {
                        l = l * 1;
                    } else {
                        l = l * 0;
                    }
                }
//                 mexPrintf(" the value is %f \n", *(coeffMatrix + coeffOfs));
                // && fabs(*(coeffMatrix + coeffOfs)) < *(thresholdForInvertibleElement) && fabs(*(coeffMatrix + coeffOfs)) > 0.00005
                if ( l == 1 && *(coeffMatrix + coeffOfs) != 0.0  ) {
                    
                    coeffToRemove = *(coeffMatrix + coeffOfs);
                    // We need to check if the new indices to remove are greater than or equal to the the previously removed ones.
                    // If yes, we need to add 1 to the apporpriate type of indices ( row or column )
                    monomialPosI = i;
                    monomialPosJ = j;
//                     mexPrintf("Testing position : (%d,%d) \n",monomialPosJ+1,monomialPosI+1);
                    // With an updated tracker matrix for the indices,
                    // we extract the orignal indices that map to the current indices we are trying to remove. This is to be done, as the current indices we are trying to remove
                    // are
                    tempCoeffOfs = 0;
                    for(i1=0;i1<origCsSize;i1++) {
                        for(j1=0;j1<origCsSize;j1++) {
                            if(j1 == monomialPosJ && i1 == monomialPosI) {
                                newMonPosJ =  *(indxTrackMat + tempCoeffOfs) + 1;
                                newMonPosI =  *(indxTrackMat + origCsSize*origCsSize + tempCoeffOfs) + 1;
                                break;
                            }
                            tempCoeffOfs++;
                        }
                    }
                    // With the corresponding indices extracted for the original coefficient matrix,
                    // we now check if the corresponding matrix to be extracted at this stage has elements independent of the hidden variable.
                    // We need to check that the non-zero element that we want to remove
                    // has corresponding elements from the row nad column corresponding to the the previously removed elements
                    // are independent of the hidden variable.
                    for(k=redInd-1;k>=0;k--) {
                        int colInd = *(colsToRemove + k) - 1, rowInd = *(rowsToRemove + k) - 1;
                        for(h=1;h<origNoOfCoeffsMatrices;h++) {
                            if(*(origCoeffMatrix + h * origCsSize*origCsSize + (newMonPosJ  - 1) + colInd * origCsSize ) == 0.0) {
                                l = l * 1;
                            } else {
                                l = l * 0;
                            }
                            if(*(origCoeffMatrix + h * origCsSize*origCsSize + (newMonPosI -1) * origCsSize + rowInd ) == 0.0) {
                                l = l * 1;
                                
                            } else {
                                l = l * 0;
                            }
                            
                        }   
//                         mexPrintf("                   Checking for the old removed position : (%d:%d) ---> %d \n",rowInd+1,colInd+1,l );

                    }

                    if (l==1){
                        // we update the matrix with elements being the corresponding indices in the original matrix.
                        tempCoeffOfs = 0;
                        coeffOfsOfs = 0;                       
                        for(i1=0;i1<origCsSize;i1++) {
                            for(j1=0;j1<origCsSize;j1++) {
                                if(monomialPosJ<=j1)
                                {
                                    *(indxTrackMat + tempCoeffOfs) = *(indxTrackMat + tempCoeffOfs + 1);
                                } else {
                                    *(indxTrackMat + tempCoeffOfs) = *(indxTrackMat + tempCoeffOfs);
                                }
                                
                                if(monomialPosI<=i1)
                                {
//                                      *(indxTrackMat + origCsSize*origCsSize + tempCoeffOfs2 )
                                    *(indxTrackMat + origCsSize*origCsSize + tempCoeffOfs ) = *(indxTrackMat + origCsSize*origCsSize + tempCoeffOfs  + origCsSize);
                                } else {
                                    *(indxTrackMat + origCsSize*origCsSize + tempCoeffOfs ) = *(indxTrackMat + origCsSize*origCsSize + tempCoeffOfs );
                                }
                                tempCoeffOfs++;
                                if(monomialPosI == i1-1 || monomialPosJ == j1-1) {
                                    coeffOfsOfs++;
                                }
                            }
                        }
//                         mexPrintf(" Element posiiton found. Original pos : (%d:%d) and current position : (%d:%d) \n", newMonPosJ,newMonPosI,monomialPosJ+1,monomialPosI+1 );
                        break;
                    } else {
//                         mexPrintf(" Element posiiton tried. Original pos : (%d:%d) and current position : (%d:%d) \n", newMonPosJ,newMonPosI,monomialPosJ+1,monomialPosI+1 );
                        monomialPosI = -1;
                        monomialPosJ = -1;                        
                    }                    
                }
                coeffOfs--;
            }
            if (monomialPosI != -1) {
                break;
            }
        }
        
        // In case we find an element to remove then only we do the next steps. Else skip to exit..
        if (monomialPosI != -1) {
            // While saving the indices we add one to the index selected, as these indices will be used by MATLAB where indices start from 1
            // instead of here, where C compiler considers indices from 0;
            *(colsToRemove + rowToRemoveIndx) = newMonPosI;
            *(rowsToRemove + rowToRemoveIndx) = newMonPosJ;
            rowToRemoveIndx++;
            
//             mexPrintf("\n ............. \n row no: %d >>> col no: %d >>> The element to be removed: %f\n", monomialPosJ, monomialPosI, coeffToRemove);
            // With the position to remove found, we now move to the task of determining the coefficient matrices
            // The first part is the M11 matrix
            
            coeffOfs = 0;
            redCoeffOfs = 0;
            for (k=0;k<noOfCoeffsMatrices;k++){
//             mexPrintf(" \n ... %f, %d \n",coeffToRemove,k);
                coeffOfs = 0;
                redCoeffOfs = 0;
                for(i=0; i<rowCntCoeffMatrix; i++){
                    for(j=0; j<rowCntCoeffMatrix; j++){
                        if(monomialPosI != i && monomialPosJ != j) {
                            *(reducedCoeffMatrix+ k * (reducedMatrixRowCnt * reducedMatrixRowCnt) + redCoeffOfs) = (*(coeffMatrix + k*(sizeOfEachCoeffMatrix) +coeffOfs));
                            redCoeffOfs++;
                        }
                        coeffOfs++;
                    }
                }
            }
            
            for (k=0;k<noOfCoeffsMatrices;k++){
                
                // With the position found for the best element to be used for monomial
                // reduction, we attempt to recreate the reduced coeffmatrix and the reduced
                // monomial basis.
                coeffOfs = 0;
                for(i=0;i<rowCntCoeffMatrix;i++){
                    if(monomialPosJ != i) {
                        *(M+coeffOfs) = *(coeffMatrix + rowCntCoeffMatrix*rowCntCoeffMatrix*k + monomialPosI*rowCntCoeffMatrix + i );
                        coeffOfs++;
                    }
                }
                for (l=0;l<noOfCoeffsMatrices;l++){
                    coeffOfs = 0;
                    for(j=0;j<rowCntCoeffMatrix;j++){
                        if(monomialPosI != j) {
                            *(N+coeffOfs) = *(coeffMatrix + rowCntCoeffMatrix*rowCntCoeffMatrix*l + rowCntCoeffMatrix*j + monomialPosJ );
                            coeffOfs++;
                        }
                    }
//              (double *M, double *N, int sizeOfEachVector, double*outerProductResult, int offsetToAddTo, double divisor)
                    outerProduct(M,N, reducedCoeffMatrix, reducedMatrixRowCnt, k+l, zeroCoeffMatrixIndices, coeffToRemove);
                }
            }
        }
        
        // Finally with the coefficient matrix truncated, we further truncate the monomial basis as well.
        redCoeffOfs = 0;
        coeffOfs = 0;
        for (i=0;i<basisCCnt;i++){
            for (j=0;j<basisRCnt;j++){
                if(i != monomialPosI){
                    *(reducedMonomialBasis + redCoeffOfs) = *(basis+ coeffOfs);
                    redCoeffOfs++;
                }
                coeffOfs++;
            }
        }
        
        free(M);
        free(N);
        
        // Variables arelated to input parameters
        coeffOfs = 0;
        redCoeffOfs =0;
        rowCntCoeffMatrix = reducedMatrixRowCnt;
        basisRCnt = basisRCnt;
        basisCCnt = basisCCnt - 1 ;
        noOfCoeffsMatrices = 0;
        for(i=0;i<totalNoOfPossiblyReducedCoeffMatrices;i++){
            if(*(zeroCoeffMatrixIndices + i)==1) {
                noOfCoeffsMatrices++;
            }
        }
        // In case no reduction was possible we break off from here.
        if (monomialPosI == -1 || monomialPosJ == -1 ) {
            mexPrintf("breaking");
            break;
        } else {
            
        }
        
        sizeOfEachCoeffMatrix = rowCntCoeffMatrix * rowCntCoeffMatrix;
        reducedMatrixRowCnt = rowCntCoeffMatrix - 1;
        basis = reducedMonomialBasis;
        coeffMatrix = reducedCoeffMatrix;
    }
    
    free(origCoeffMatrix);
    mexPrintf(" \n The reduction in monomial bases could be achieved upto %d monomials. \n ", redInd);            
}