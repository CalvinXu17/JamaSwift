//
//  SingularValueDecomposition.swift
//  Sama
//
//  Created by Calvin on 2019/7/12.
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

public class SingularValueDecomposition : Codable
{
    private var U: [[Double]]
    private var V: [[Double]]
    private var s: [Double]
    private var m: Int
    private var n: Int
    
    public init(paramMatrix: Matrix)
    {
        var arrayOfDouble = paramMatrix.getArrayCopy()
        self.m = paramMatrix.getRowDimension()
        self.n = paramMatrix.getColumnDimension()
        
        let i = min(self.m, self.n)
        self.s = [Double].init(repeating: 0.0, count: min(self.m + 1, self.n))
        self.U = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: i), count: self.m)
        self.V = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: self.n), count: self.n)
        var arrayOfDouble1 = [Double].init(repeating: 0.0, count: self.n)
        var arrayOfDouble2 = [Double].init(repeating: 0.0, count: self.m)
        var j = 1
        var k = 1
        
        let i1 = min(self.m - 1, self.n)
        let i2 = max(0, min(self.n - 2, self.m))
        for  i3 in 0 ..< max(i1, i2)
        {
            if (i3 < i1)
            {
                self.s[i3] = 0.0
                for i4 in i3 ..< self.m {
                    self.s[i3] = hypot(self.s[i3], arrayOfDouble[i4][i3])
                }
                if (self.s[i3] != 0.0)
                {
                    if (arrayOfDouble[i3][i3] < 0.0) {
                        self.s[i3] = (-self.s[i3])
                    }
                    for i4 in i3 ..< self.m {
                        arrayOfDouble[i4][i3] /= self.s[i3]
                    }
                    arrayOfDouble[i3][i3] += 1.0
                }
                self.s[i3] = (-self.s[i3])
            }
            
            for i4 in i3 + 1 ..< self.n
            {
                if (((i3 < i1 ? 1 : 0) & (self.s[i3] != 0.0 ? 1 : 0)) != 0)
                {
                    var d1:Double = 0.0
                    for i7 in i3 ..< self.m {
                        d1 += arrayOfDouble[i7][i3] * arrayOfDouble[i7][i4]
                    }
                    d1 = -d1 / arrayOfDouble[i3][i3]
                    for i7 in i3 ..< self.m {
                        arrayOfDouble[i7][i4] += d1 * arrayOfDouble[i7][i3]
                    }
                }
                arrayOfDouble1[i4] = arrayOfDouble[i3][i4]
            }
            if ((j & (i3 < i1 ? 1 : 0)) != 0) {
                for i4 in i3 ..< self.m {
                    self.U[i4][i3] = arrayOfDouble[i4][i3]
                }
            }
            if (i3 < i2)
            {
                arrayOfDouble1[i3] = 0.0
                for i4 in i3 + 1 ..< self.n {
                    arrayOfDouble1[i3] = hypot(arrayOfDouble1[i3], arrayOfDouble1[i4])
                }
                if (arrayOfDouble1[i3] != 0.0)
                {
                    if (arrayOfDouble1[(i3 + 1)] < 0.0) {
                        arrayOfDouble1[i3] = (-arrayOfDouble1[i3])
                    }
                    for i4 in i3 + 1 ..< self.n {
                        arrayOfDouble1[i4] /= arrayOfDouble1[i3]
                    }
                    arrayOfDouble1[(i3 + 1)] += 1.0
                }
                arrayOfDouble1[i3] = (-arrayOfDouble1[i3])
                if (((i3 + 1 < self.m ? 1 : 0) & (arrayOfDouble1[i3] != 0.0 ? 1 : 0)) != 0)
                {
                    for i4 in i3 + 1 ..< self.m {
                        arrayOfDouble2[i4] = 0.0
                    }
                    for i4 in i3 + 1 ..< self.n {
                        for i5 in i3 + 1 ..< self.m {
                            arrayOfDouble2[i5] += arrayOfDouble1[i4] * arrayOfDouble[i5][i4]
                        }
                    }
                    for i4 in i3 + 1 ..< self.n
                    {
                        let d2 = -arrayOfDouble1[i4] / arrayOfDouble1[(i3 + 1)]
                        for i7 in i3 + 1 ..< self.m {
                            arrayOfDouble[i7][i4] += d2 * arrayOfDouble2[i7]
                        }
                    }
                }
                if (k != 0) {
                    for i4 in i3 + 1 ..< self.n {
                        self.V[i4][i3] = arrayOfDouble1[i4]
                    }
                }
            }
        }
        var i3 = min(self.n, self.m + 1)
        if (i1 < self.n) {
            self.s[i1] = arrayOfDouble[i1][i1]
        }
        if (self.m < i3) {
            self.s[(i3 - 1)] = 0.0
        }
        if (i2 + 1 < i3) {
            arrayOfDouble1[i2] = arrayOfDouble[i2][(i3 - 1)]
        }
        arrayOfDouble1[(i3 - 1)] = 0.0
        if (j != 0)
        {
            for i4 in i1 ..< i
            {
                for i6 in 0 ..< self.m {
                    self.U[i6][i4] = 0.0
                }
                self.U[i4][i4] = 1.0
            }
            for i4 in stride(from: i1 - 1, through: 0, by: -1) {
                if (self.s[i4] != 0.0)
                {
                    for i6 in i4 + 1 ..< i
                    {
                        var d3: Double = 0.0
                        for i8 in i4 ..< self.m {
                            d3 += self.U[i8][i4] * self.U[i8][i6]
                        }
                        d3 = -d3 / self.U[i4][i4]
                        for i8 in i4  ..< self.m {
                            self.U[i8][i6] += d3 * self.U[i8][i4]
                        }
                    }
                    for i6 in i4 ..< self.m {
                        self.U[i6][i4] = (-self.U[i6][i4])
                    }
                    self.U[i4][i4] = (1.0 + self.U[i4][i4])
                    if i4 - 1 >= 0 {
                        for i6 in 0 ..< i4 - 1 {
                            self.U[i6][i4] = 0.0
                        }
                    }
                }
                else
                {
                    for i6 in 0 ..< self.m {
                        self.U[i6][i4] = 0.0
                    }
                    self.U[i4][i4] = 1.0
                }
            }
        }
        if (k != 0) {
            for i4 in stride(from: self.n - 1, through: 0, by: -1)
            {
                if (((i4 < i2 ? 1 : 0) & (arrayOfDouble1[i4] != 0.0 ? 1 : 0)) != 0) {
                    for i6 in i4 + 1 ..< i
                    {
                        var d3: Double = 0.0
                        for i8 in i4 + 1 ..< self.n {
                            d3 += self.V[i8][i4] * self.V[i8][i6]
                        }
                        d3 = -d3 / self.V[(i4 + 1)][i4]
                        for i8 in i4 + 1 ..< self.n {
                            self.V[i8][i6] += d3 * self.V[i8][i4]
                        }
                    }
                }
                for i6 in 0 ..< self.n {
                    self.V[i6][i4] = 0.0
                }
                self.V[i4][i4] = 1.0
            }
        }
        let i4 = i3 - 1
        var i6 = 0
        let d3 = Double(pow(2.0, -52.0))
        let d4 = Double(pow(2.0, -966.0))
        while (i3 > 0)
        {
            var i9 = 0
            for i99 in stride(from: i3 - 2, through: -1, by: -1)
            {
                i9 = i99
                if (i99 == -1) {
                    break
                }
                if (abs(arrayOfDouble1[i99]) <= d4 + d3 * (abs(self.s[i99]) + abs(self.s[(i99 + 1)])))
                {
                    arrayOfDouble1[i99] = 0.0
                    break
                }
            }
            var i10 = 0
            if (i9 == i3 - 2)
            {
                i10 = 4
            }
            else
            {
                var i11 = 0
                for i111 in stride(from: i3 - 1, through: i9, by: -1)
                {
                    i11 = i111
                    if (i111 == i9) {
                        break
                    }
                    let d7:Double  = (i111 != i3 ? abs(arrayOfDouble1[i111]) : 0.0) + (i111 != i9 + 1 ? abs(arrayOfDouble1[(i111 - 1)]) : 0.0)
                    if (abs(self.s[i111]) <= d4 + d3 * d7)
                    {
                        self.s[i111] = 0.0
                        break
                    }
                }
                if (i11 == i9)
                {
                    i10 = 3
                }
                else if (i11 == i3 - 1)
                {
                    i10 = 1
                }
                else
                {
                    i10 = 2
                    i9 = i11
                }
            }
            i9 += 1
            var d5: Double = 0.0
            var d9: Double = 0.0
            var d11: Double = 0.0
            var d13: Double = 0.0
            switch i10
            {
            case 1:
                d5 = arrayOfDouble1[(i3 - 2)]
                arrayOfDouble1[(i3 - 2)] = 0.0
                for i13 in stride(from: i3 - 2, through: i9, by: -1)
                {
                    d9 = hypot(self.s[i13], d5)
                    d11 = self.s[i13] / d9
                    d13 = d5 / d9
                    self.s[i13] = d9
                    if (i13 != i9)
                    {
                        d5 = -d13 * arrayOfDouble1[(i13 - 1)]
                        arrayOfDouble1[(i13 - 1)] = (d11 * arrayOfDouble1[(i13 - 1)])
                    }
                    if (k != 0) {
                        for i15 in 0 ..< self.n
                        {
                            d9 = d11 * self.V[i15][i13] + d13 * self.V[i15][(i3 - 1)]
                            self.V[i15][(i3 - 1)] = (-d13 * self.V[i15][i13] + d11 * self.V[i15][(i3 - 1)])
                            self.V[i15][i13] = d9
                        }
                    }
                }
                break
            case 2:
                d5 = arrayOfDouble1[(i9 - 1)]
                arrayOfDouble1[(i9 - 1)] = 0.0
                for i13 in i9 ..< i3
                {
                    d9 = hypot(self.s[i13], d5)
                    d11 = self.s[i13] / d9
                    d13 = d5 / d9
                    self.s[i13] = d9
                    d5 = -d13 * arrayOfDouble1[i13]
                    arrayOfDouble1[i13] = (d11 * arrayOfDouble1[i13])
                    if (j != 0) {
                        for i15 in 0 ..< self.m
                        {
                            d9 = d11 * self.U[i15][i13] + d13 * self.U[i15][(i9 - 1)]
                            self.U[i15][(i9 - 1)] = (-d13 * self.U[i15][i13] + d11 * self.U[i15][(i9 - 1)])
                            self.U[i15][i13] = d9
                        }
                    }
                }
                break
            case 3:
                d5 = max(max(max(max(abs(self.s[(i3 - 1)]), abs(self.s[(i3 - 2)])), abs(arrayOfDouble1[(i3 - 2)])), abs(self.s[i9])), abs(arrayOfDouble1[i9]))
                
                let d8 = self.s[(i3 - 1)] / d5
                let d10 = self.s[(i3 - 2)] / d5
                let d12 = arrayOfDouble1[(i3 - 2)] / d5
                let d14 = self.s[i9] / d5
                let d15 = arrayOfDouble1[i9] / d5
                let d16 = ((d10 + d8) * (d10 - d8) + d12 * d12) / 2.0
                let d17 = d8 * d12 * (d8 * d12)
                var d18:Double = 0.0
                if (((d16 != 0.0 ? 1 : 0) | (d17 != 0.0 ? 1 : 0)) != 0)
                {
                    d18 = sqrt(d16 * d16 + d17)
                    if (d16 < 0.0) {
                        d18 = -d18
                    }
                    d18 = d17 / (d16 + d18)
                }
                var d19: Double = (d14 + d8) * (d14 - d8) + d18
                var d20: Double = d14 * d15
                for i16 in i9 ..< i3 - 1
                {
                    var d21: Double = hypot(d19, d20)
                    var d22 = d19 / d21
                    var d23 = d20 / d21
                    if (i16 != i9) {
                        arrayOfDouble1[(i16 - 1)] = d21
                    }
                    d19 = d22 * self.s[i16] + d23 * arrayOfDouble1[i16]
                    arrayOfDouble1[i16] = (d22 * arrayOfDouble1[i16] - d23 * self.s[i16])
                    d20 = d23 * self.s[(i16 + 1)]
                    self.s[(i16 + 1)] = (d22 * self.s[(i16 + 1)])
                    if (k != 0) {
                        for i17 in 0 ..< self.n
                        {
                            d21 = d22 * self.V[i17][i16] + d23 * self.V[i17][(i16 + 1)]
                            self.V[i17][(i16 + 1)] = (-d23 * self.V[i17][i16] + d22 * self.V[i17][(i16 + 1)])
                            self.V[i17][i16] = d21
                        }
                    }
                    d21 = hypot(d19, d20)
                    d22 = d19 / d21
                    d23 = d20 / d21
                    self.s[i16] = d21
                    d19 = d22 * arrayOfDouble1[i16] + d23 * self.s[(i16 + 1)]
                    self.s[(i16 + 1)] = (-d23 * arrayOfDouble1[i16] + d22 * self.s[(i16 + 1)])
                    d20 = d23 * arrayOfDouble1[(i16 + 1)]
                    arrayOfDouble1[(i16 + 1)] = (d22 * arrayOfDouble1[(i16 + 1)])
                    if ((j != 0) && (i16 < self.m - 1)) {
                        for i17 in 0 ..< self.m
                        {
                            d21 = d22 * self.U[i17][i16] + d23 * self.U[i17][(i16 + 1)]
                            self.U[i17][(i16 + 1)] = (-d23 * self.U[i17][i16] + d22 * self.U[i17][(i16 + 1)])
                            self.U[i17][i16] = d21
                        }
                    }
                }
                arrayOfDouble1[(i3 - 2)] = d19
                i6 += 1
                
                break
            case 4:
                if (self.s[i9] <= 0.0)
                {
                    self.s[i9] = (self.s[i9] < 0.0 ? -self.s[i9] : 0.0)
                    if (k != 0) {
                        for i12 in 0 ... i4 {
                            self.V[i12][i9] = (-self.V[i12][i9])
                        }
                    }
                }
                while ((i9 < i4) && (self.s[i9] < self.s[(i9 + 1)]))
                {
                    var d6: Double = self.s[i9]
                    self.s[i9] = self.s[(i9 + 1)]
                    self.s[(i9 + 1)] = d6
                    if ((k != 0) && (i9 < self.n - 1)) {
                        for i14 in 0 ..< self.n
                        {
                            d6 = self.V[i14][(i9 + 1)]
                            self.V[i14][(i9 + 1)] = self.V[i14][i9]
                            self.V[i14][i9] = d6
                        }
                    }
                    if ((j != 0) && (i9 < self.m - 1)) {
                        for i14 in  0 ..< self.m
                        {
                            d6 = self.U[i14][(i9 + 1)]
                            self.U[i14][(i9 + 1)] = self.U[i14][i9]
                            self.U[i14][i9] = d6
                        }
                    }
                    i9 += 1
                }
                i6 = 0
                i3 -= 1
            default:
                break
            }
        }
    }
    
    public func getU() -> Matrix
    {
        return Matrix(paramArrayOfDouble: self.U, paramInt1: self.m, paramInt2: min(self.m + 1, self.n))
    }
    
    public func getV() -> Matrix
    {
        return Matrix(paramArrayOfDouble: self.V, paramInt1: self.n, paramInt2: self.n)
    }
    
    public func getSingularValues() -> [Double]
    {
        return self.s
    }
    
    public func getS() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.n, paramInt2: self.n)
        for  i in 0 ..< self.n
        {
            for  j in 0 ..< self.n {
                localMatrix.A[i][j] = 0.0
            }
            localMatrix.A[i][i] = self.s[i]
        }
        return localMatrix
    }
    
    public func norm2() -> Double
    {
        return self.s[0]
    }
    
    public func cond() -> Double
    {
        return self.s[0] / self.s[(min(self.m, self.n) - 1)]
    }
    
    public func rank() -> Int
    {
        let d1: Double = pow(2.0, -52.0)
        let d2: Double = Double(max(self.m, self.n)) * self.s[0] * d1
        var i = 0
        for j in 0 ..< self.s.count{
            if (self.s[j] > d2) {
                i += 1
            }
        }
        return i
    }
}
