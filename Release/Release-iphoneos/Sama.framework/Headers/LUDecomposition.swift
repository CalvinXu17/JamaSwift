//
//  LUDecomposition.swift
//  Sama
//
//  Created by Calvin on 2019/7/13.
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

public class LUDecomposition: Codable {
    
    /* ------------------------
     Class variables
     * ------------------------ */
    
    /** Array for internal storage of decomposition.
     @serial internal array storage.
     */
    private var LU: [[Double]]
    
    /** Row and column dimensions, and pivot sign.
     @serial column dimension.
     @serial row dimension.
     @serial pivot sign.
     */
    private var m, n, pivsign: Int
    
    /** Internal storage of pivot vector.
     @serial pivot vector.
     */
    private var piv: [Int]
    
    /* ------------------------
     Constructor
     * ------------------------ */
    
    /** LU Decomposition
     Structure to access L, U and piv.
     @param  A Rectangular matrix
     */
    
    public init(A: Matrix) {
        
        // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
        
        LU = A.getArrayCopy()
        m = A.getRowDimension()
        n = A.getColumnDimension()
        piv =  [Int].init(repeating: 0, count: self.m)
        for i in 0 ..< m {
            piv[i] = i
        }
        pivsign = 1
        var LUcolj = [Double].init(repeating: 0.0, count: self.m)
        
        // Outer loop.
        
        for j in 0 ..< n {
            
            // Make a copy of the j-th column to localize references.
            
            for i in 0 ..< m {
                LUcolj[i] = LU[i][j]
            }
            
            // Apply previous transformations.
            
            for i in 0 ..< m {
                
                // Most of the time is spent in the following dot product.
                
                let kmax: Int = min(i,j)
                var s: Double = 0.0
                for k in 0 ..< kmax {
                    s += LU[i][k] * LUcolj[k]
                }
                LUcolj[i] -= s
                LU[i][j] = LUcolj[i]
            }
            
            // Find pivot and exchange if necessary.
            
            var p: Int = j
            for i in stride(from: j + 1, to: m, by: 1) {
                if (abs(LUcolj[i]) > abs(LUcolj[p])) {
                    p = i
                }
            }
            if (p != j) {
                for k in 0 ..< n {
                    let t: Double = LU[p][k]
                    LU[p][k] = LU[j][k]
                    LU[j][k] = t
                }
                let k: Int = piv[p]
                piv[p] = piv[j]
                piv[j] = k
                pivsign = -pivsign
            }
            
            // Compute multipliers.
            
            if j <= LU.count - 1 && j <= LU[j].count - 1 && ((j < self.m ? 1 : 0) & (self.LU[j][j] != 0.0 ? 1 : 0)) != 0 {
                for i in j + 1 ..< self.m {
                    if i <= LU.count - 1 && j <= LU[i].count - 1 && j <= LU.count - 1 {
                        LU[i][j] /= LU[j][j]
                    }
                }
            }
        }
    }
    
    public func isNonsingular() -> Bool {
        for j in 0 ..< n {
            if (LU[j][j] == 0) {
                return false
            }
        }
        return true
    }
    
    /** Return lower triangular factor
     @return     L
     */
    
    public func getL() -> Matrix {
        let X = Matrix(paramInt1: m,paramInt2: n)
        for i in 0 ..< m {
            for j in 0 ..< n {
                if (i > j) {
                    X.A[i][j] = LU[i][j]
                } else if (i == j) {
                    X.A[i][j] = 1.0
                } else {
                    X.A[i][j] = 0.0
                }
            }
        }
        return X
    }
    
    /** Return upper triangular factor
     @return     U
     */
    
    public func getU() -> Matrix {
        let X = Matrix(paramInt1: n,paramInt2: n)
        for i in 0 ..< n {
            for j in 0 ..< n {
                if i <= j {
                    X.A[i][j] = LU[i][j]
                } else {
                    X.A[i][j] = 0.0
                }
            }
        }
        return X
    }
    
    /** Return pivot permutation vector
     @return     piv
     */
    
    public func getPivot() -> [Int]{
        var p = [Int].init(repeating: 0, count: self.m)
        for i in 0 ..< self.m {
            p[i] = piv[i]
        }
        return p
    }
    
    /** Return pivot permutation vector as a one-dimensional double array
     @return     (double) piv
     */
    
    public func getDoublePivot() -> [Double] {
        var vals = [Double].init(repeating: 0.0, count: self.m)
        for i in 0 ..< m {
            vals[i] = Double(piv[i])
        }
        return vals
    }
    
    /** Determinant
     @return     det(A)
     @exception  IllegalArgumentException  Matrix must be square
     */
    
    public func det() throws -> Double {
        if (m != n) {
            throw MatrixError.IllegalArgumentException(info: "Matrix must be square.")
        }
        var d: Double = Double(pivsign)
        
        for j in 0 ..< n {
            d *= LU[j][j]
        }
        return d
    }
    
    /** Solve A*X = B
     @param  B   A Matrix with as many rows as A and any number of columns.
     @return     X so that L*U*X = B(piv,:)
     @exception  IllegalArgumentException Matrix row dimensions must agree.
     @exception  RuntimeException  Matrix is singular.
     */
    
    public func solve(B: Matrix) throws -> Matrix {
        if B.getRowDimension() != m {
            throw MatrixError.IllegalArgumentException(info: "Matrix row dimensions must agree.")
        }
        if !self.isNonsingular() {
            throw MatrixError.RuntimeException(info: "Matrix is singular.")
        }
        
        // Copy right hand side with pivoting
        let nx = B.getColumnDimension()
        let Xmat = try B.getMatrix(paramArrayOfInt: piv,int: 0,int: nx-1)
        
        // Solve L*Y = B(piv,:)
        for k in 0 ..< n {
            for i in k + 1 ..< n {
                for j in 0 ..< nx {
                    Xmat.A[i][j] -= Xmat.A[k][j] * LU[i][k]
                }
            }
        }
        // Solve U*X = Y
        for k in stride(from: n - 1, through: 0, by: -1) {
            for j in 0 ..< nx {
                Xmat.A[k][j] /= LU[k][k]
            }
            for i in 0 ..< k {
                for j in 0 ..< nx {
                    Xmat.A[i][j] -= Xmat.A[k][j] * LU[i][k]
                }
            }
        }
        return Xmat
    }
}
