//
//  CholeskyDecomposition.swift
//  Sama
//
//  Created by Calvin on 2019/7/12.
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

public class CholeskyDecomposition: Codable
{
    private var L: [[Double]]
    private var n: Int
    private var isspd: Bool
    
    public init(paramMatrix: Matrix)
    {
        var arrayOfDouble: [[Double]] = paramMatrix.getArray()
        self.n = paramMatrix.getRowDimension()
        self.L = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: self.n), count: self.n)
        self.isspd = (paramMatrix.getColumnDimension() == self.n)
        
        for i in 0 ..< self.n {
            var arrayOfDouble1: [Double] = self.L[i]
            var d1: Double = 0.0
            for j in 0 ..< i {
                var arrayOfDouble2: [Double] = self.L[j]
                var d2: Double = 0.0
                for k in 0 ..< j {
                    d2 += arrayOfDouble2[k] * arrayOfDouble1[k]
                }
                d2 = (arrayOfDouble[i][j] - d2) / self.L[j][j]
                arrayOfDouble1[j] = d2
                L[i][j] = d2
                self.L[i][j] = d2
                d1 += d2 * d2
                self.isspd = self.isspd && ( arrayOfDouble[j][i] == arrayOfDouble[i][j] )
            }
            d1 = arrayOfDouble[i][i] - d1
            self.isspd = self.isspd && ( d1 > 0.0 )
            self.L[i][i] = sqrt(Double(max(d1, 0.0)))
            for j in i + 1 ..< self.n {
                self.L[i][j] = 0.0
            }
        }
    }
    
    public func isSPD() -> Bool
    {
        return self.isspd
    }
    
    public func getL() -> Matrix
    {
        return Matrix(paramArrayOfDouble: self.L, paramInt1: self.n, paramInt2: self.n)
    }
    
    public func solve(paramMatrix: Matrix) throws -> Matrix
    {
        
        if (paramMatrix.getRowDimension() != self.n) {
            throw MatrixError.IllegalArgumentException(info: "Matrix row dimensions must agree.")
        }
        if (!self.isspd) {
            throw MatrixError.RuntimeException(info: "Matrix is not symmetric positive definite.")
        }
        var arrayOfDouble: [[Double]] = paramMatrix.getArrayCopy()
        let i: Int = paramMatrix.getColumnDimension()
        for j in 0 ..< self.n {
            for k in 0 ..< i {
                for m in 0 ..< j {
                    arrayOfDouble[j][k] -= arrayOfDouble[m][k] * self.L[j][m]
                }
                arrayOfDouble[j][k] /= self.L[j][j]
            }
        }
        for j in (0 ... self.n - 1).reversed()  {
            for k in 0 ..< i {
                for m in j + 1 ..< self.n {
                    arrayOfDouble[j][k] -= arrayOfDouble[m][k] * self.L[m][j]
                }
                arrayOfDouble[j][k] /= self.L[j][j]
            }
        }
        return Matrix(paramArrayOfDouble: arrayOfDouble, paramInt1: self.n, paramInt2: i)
    }
}
