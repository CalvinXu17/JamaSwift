//
//  QRDecomposition.swift
//  Sama
//
//  Created by Calvin on 2019/7/12.
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

public class QRDecomposition: Codable
{
    private var QR: [[Double]]
    private var m: Int
    private var n: Int
    private var Rdiag: [Double]
    
    public init(paramMatrix: Matrix)
    {
        self.QR = paramMatrix.getArrayCopy()
        self.m = paramMatrix.getRowDimension()
        self.n = paramMatrix.getColumnDimension()
        self.Rdiag = [Double].init(repeating: 0.0, count: self.n)
        
        for i in stride(from: 0, to: self.n, by: 1)
        {
            var d1: Double = 0.0
            for j in stride(from: i, to: self.m, by: 1) {
                d1 = hypot(d1, self.QR[j][i])
            }
            if (d1 != 0.0)
            {
                if (self.QR[i][i] < 0.0) {
                    d1 = -d1
                }
                for j in stride(from: i, to: self.m, by: 1) {
                    self.QR[j][i] /= d1
                }
                self.QR[i][i] += 1.0
                for j in stride(from: i + 1, to: self.n, by: 1)
                {
                    var d2: Double = 0.0
                    for k in stride(from: i, to: self.m, by: 1) {
                        d2 += self.QR[k][i] * self.QR[k][j]
                    }
                    d2 = -d2 / self.QR[i][i]
                    for k in stride(from: i, to: self.m, by: 1) {
                        self.QR[k][j] += d2 * self.QR[k][i]
                    }
                }
            }
            self.Rdiag[i] = (-d1)
        }
    }
    
    public func isFullRank() -> Bool
    {
        for i in 0 ..< self.n {
            if (self.Rdiag[i] == 0.0) {
                return false
            }
        }
        return true
    }
    
    public func getH() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                if (i >= j) {
                    localMatrix.A[i][j] = self.QR[i][j]
                } else {
                    localMatrix.A[i][j] = 0.0
                }
            }
        }
        return localMatrix
    }
    
    public func getR() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.n, paramInt2: self.n)
        for i in 0 ..< self.n {
            for j in 0 ..< self.n {
                if (i < j) {
                    localMatrix.A[i][j] = self.QR[i][j]
                } else if (i == j) {
                    localMatrix.A[i][j] = self.Rdiag[i]
                } else {
                    localMatrix.A[i][j] = 0.0
                }
            }
        }
        return localMatrix
    }
    
    public func getQ() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for  i in stride(from: self.n - 1, through: 0, by: -1)
        {
            for j in 0 ..< self.m {
                localMatrix.A[j][i] = 0.0
            }
            localMatrix.A[i][i] = 1.0
            for j in i ..< self.n {
                if (self.QR[i][i] != 0.0)
                {
                    var d:Double = 0.0
                    for  k in i ..< self.m {
                        d += self.QR[k][i] * localMatrix.A[k][j]
                    }
                    d = -d / self.QR[i][i]
                    for k in i ..< self.m {
                        localMatrix.A[k][j] += d * self.QR[k][i]
                    }
                }
            }
        }
        return localMatrix
    }
    
    public func solve(paramMatrix: Matrix) throws -> Matrix
    {
        if (paramMatrix.getRowDimension() != self.m) {
            throw MatrixError.IllegalArgumentException(info: "Matrix row dimensions must agree.")
        }
        if (!isFullRank()) {
            throw MatrixError.RuntimeException(info: "Matrix is rank deficient.")
        }
        let i = paramMatrix.getColumnDimension()
        var arrayOfDouble = paramMatrix.getArrayCopy()
        for j in 0 ..< self.n {
            for k in 0 ..< i
            {
                var  d: Double = 0.0
                for i2 in j ..< self.m {
                    d += self.QR[i2][j] * arrayOfDouble[i2][k]
                }
                d = -d / self.QR[j][j]
                for i2 in j ..< self.m {
                    arrayOfDouble[i2][k] += d * self.QR[i2][j]
                }
            }
        }
        for j in stride(from: self.n - 1, through: 0, by: -1)
        {
            for k in 0 ..< i {
                arrayOfDouble[j][k] /= self.Rdiag[j]
            }
            for k in 0 ..< j {
                for i1 in 0  ..< i {
                    arrayOfDouble[k][i1] -= arrayOfDouble[j][i1] * self.QR[k][j]
                }
            }
        }
        return try Matrix(paramArrayOfDouble: arrayOfDouble, paramInt1: self.n, paramInt2: i).getMatrix(paramInt1: 0, paramInt2: self.n - 1, paramInt3: 0, paramInt4: i - 1)
    }
}
