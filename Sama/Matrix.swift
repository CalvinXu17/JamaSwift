//
//  Matrix.swift
//  Sama
//
//  Created by Calvin on 2019/7/12.
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

public class Matrix: Codable {
    public var A: [[Double]]
    private var m: Int
    private var n: Int
    
    public init(paramInt1: Int, paramInt2: Int)
    {
        self.m = paramInt1
        self.n = paramInt2
        self.A = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: paramInt2), count: paramInt1)
    }
    
    public init(paramInt1: Int, paramInt2: Int, paramDouble: Double)
    {
        self.m = paramInt1
        self.n = paramInt2
        self.A = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: paramInt2), count: paramInt1)
        for i in 0 ..< paramInt1 {
            for j in 0 ..< paramInt2 {
                self.A[i][j] = paramDouble
            }
        }
    }
    
    public init(paramArrayOfDouble: [[Double]]) throws
    {
        self.m = paramArrayOfDouble.count
        self.n = paramArrayOfDouble[0].count
        for i in stride(from: 0, to: self.m, by: 1) {
            if (paramArrayOfDouble[i].count != self.n) {
                throw MatrixError.IllegalArgumentException(info: "All rows must have the same length.")
            }
        }
        self.A = paramArrayOfDouble
    }
    
    public init(paramArrayOfDouble: [[Double]], paramInt1: Int, paramInt2: Int)
    {
        self.A = paramArrayOfDouble
        self.m = paramInt1
        self.n = paramInt2
    }
    
    public init(paramArrayOfDouble: [Double], paramInt: Int) throws
    {
        self.m = paramInt
        self.n = (paramInt != 0 ? paramArrayOfDouble.count / paramInt : 0)
        if (paramInt * self.n != paramArrayOfDouble.count) {
            throw MatrixError.IllegalArgumentException(info: "Array length must be a multiple of m.")
        }
        self.A =  [[Double]].init(repeating: [Double].init(repeating: 0.0, count: self.n), count: paramInt)
        
        for i in stride(from: 0, to: paramInt, by: 1) {
            for j in stride(from: 0, to: self.n, by: 1) {
                self.A[i][j] = paramArrayOfDouble[(i + j * paramInt)]
            }
        }
    }
    
    public static func constructWithCopy(paramArrayOfDouble: [[Double]]) throws -> Matrix
    {
        let i = paramArrayOfDouble.count
        let j = paramArrayOfDouble[0].count
        let localMatrix = Matrix(paramInt1: i, paramInt2: j)
        for k in stride(from: 0, to: i, by: 1)
        {
            if (paramArrayOfDouble[k].count != j) {
                throw MatrixError.IllegalArgumentException(info: "All rows must have the same length.")
            }
            for i1 in 0 ..< j {
                localMatrix.A[k][i1] = paramArrayOfDouble[k][i1]
            }
        }
        return localMatrix
    }
    
    public func copy() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for  i in 0 ..< self.m {
            for j in 0 ..< self.n {
                localMatrix.A[i][j] = self.A[i][j]
            }
        }
        return localMatrix
    }
    
    public func clone() -> Any
    {
        return copy()
    }
    
    public func getArray() -> [[Double]]
    {
        return self.A
    }
    
    public func getArrayCopy() -> [[Double]]
    {
        var arrayOfDouble: [[Double]] = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: self.n), count: self.m)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                arrayOfDouble[i][j] = self.A[i][j]
            }
        }
        return arrayOfDouble
    }
    
    public func getColumnPackedCopy() -> [Double]
    {
        var arrayOfDouble = [Double].init(repeating: 0.0, count: self.m * self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                arrayOfDouble[(i + j * self.m)] = self.A[i][j]
            }
        }
        return arrayOfDouble
    }
    
    public func getRowPackedCopy() -> [Double]
    {
        var arrayOfDouble = [Double].init(repeating: 0.0, count: self.m * self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                arrayOfDouble[(i * self.n + j)] = self.A[i][j]
            }
        }
        return arrayOfDouble
    }
    
    public func getRowDimension() -> Int
    {
        return self.m
    }
    
    public func getColumnDimension() -> Int
    {
        return self.n
    }
    
    public func get(paramInt1: Int, paramInt2: Int) -> Double
    {
        return self.A[paramInt1][paramInt2]
    }
    
    public func getMatrix(paramInt1: Int, paramInt2: Int, paramInt3: Int, paramInt4: Int) throws -> Matrix
    {
        let localMatrix = Matrix(paramInt1: paramInt2 - paramInt1 + 1, paramInt2: paramInt4 - paramInt3 + 1)
        
        for i in paramInt1...paramInt2 {
            for j in paramInt3...paramInt4 {
                
                if i > A.count - 1 || j > A[i].count - 1 || (i - paramInt1) > localMatrix.A.count - 1 || (j - paramInt3) > localMatrix.A[(i - paramInt1)].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                localMatrix.A[(i - paramInt1)][(j - paramInt3)] = self.A[i][j]
            }
        }
        return localMatrix
    }
    
    public func getMatrix(paramArrayOfInt1: [Int], paramArrayOfInt2: [Int]) throws -> Matrix
    {
        let localMatrix = Matrix(paramInt1: paramArrayOfInt1.count, paramInt2: paramArrayOfInt2.count)
        
        for i in 0 ..< paramArrayOfInt1.count {
            for j in 0 ..< paramArrayOfInt2.count {
                if paramArrayOfInt1[i] > A.count - 1 || paramArrayOfInt2[j] > A[paramArrayOfInt1[i]].count - 1 || i > localMatrix.A.count - 1 || j > localMatrix.A[i].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                localMatrix.A[i][j] = self.A[paramArrayOfInt1[i]][paramArrayOfInt2[j]]
            }
        }
        return localMatrix
    }
    
    public func getMatrix(paramInt1: Int, paramInt2: Int, paramArrayOfInt: [Int]) throws -> Matrix
    {
        let localMatrix = Matrix(paramInt1: paramInt2 - paramInt1 + 1, paramInt2: paramArrayOfInt.count)
        
        for i in paramInt1 ... paramInt2 {
            for j in 0 ..< paramArrayOfInt.count {
                if i > A.count - 1 || paramArrayOfInt[j] > A[i].count - 1 || (i - paramInt1) > localMatrix.A.count - 1 || j > localMatrix.A[(i - paramInt1)].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                localMatrix.A[(i - paramInt1)][j] = self.A[i][paramArrayOfInt[j]]
            }
        }
        
        return localMatrix
    }
    
    public func getMatrix(paramArrayOfInt: [Int], int paramInt1: Int, int paramInt2: Int) throws -> Matrix
    {
        let localMatrix = Matrix(paramInt1: paramArrayOfInt.count, paramInt2: paramInt2 - paramInt1 + 1)
        
        for i in 0 ..< paramArrayOfInt.count {
            for j in paramInt1 ... paramInt2 {
                if paramArrayOfInt[i] > A.count - 1 || j > A[paramArrayOfInt[i]].count || i > localMatrix.A.count - 1 || (j - paramInt1) > localMatrix.A[i].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                localMatrix.A[i][(j - paramInt1)] = self.A[paramArrayOfInt[i]][j]
            }
        }
        return localMatrix
    }
    
    public func set(paramInt1: Int, paramInt2: Int, paramDouble: Double)
    {
        self.A[paramInt1][paramInt2] = paramDouble
    }
    
    public func setMatrix(paramInt1: Int,paramInt2: Int,paramInt3: Int,paramInt4: Int, paramMatrix: Matrix)
    {
        for  i in paramInt1 ... paramInt2 {
            for j in paramInt3 ... paramInt4 {
                
                if i > A.count - 1 || j > A[i].count {
                    return
                }
                
                self.A[i][j] = paramMatrix.get(paramInt1: i - paramInt1, paramInt2: j - paramInt3)
            }
        }
        
    }
    
    public func setMatrix(paramArrayOfInt1: [Int], paramArrayOfInt2: [Int],paramMatrix: Matrix) throws
    {
        
        for  i in 0 ..< paramArrayOfInt1.count {
            for  j in 0 ..< paramArrayOfInt2.count {
                if paramArrayOfInt1[i] > A.count - 1 || paramArrayOfInt2[j] > A[paramArrayOfInt1[i]].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                self.A[paramArrayOfInt1[i]][paramArrayOfInt2[j]] = paramMatrix.get(paramInt1: i, paramInt2: j)
            }
        }
    }
    
    public func setMatrix(paramArrayOfInt: [Int],paramInt1: Int, paramInt2: Int,paramMatrix: Matrix) throws {
        
        
        for  i in 0 ..< paramArrayOfInt.count {
            for  j in paramInt1 ... paramInt2 {
                if paramArrayOfInt[i] > A.count - 1 || j > A[paramArrayOfInt[i]].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                self.A[paramArrayOfInt[i]][j] = paramMatrix.get(paramInt1: i, paramInt2: j - paramInt1)
            }
        }
        
    }
    
    public func setMatrix(paramInt1: Int, paramInt2: Int,paramArrayOfInt: [Int], paramMatrix: Matrix) throws
    {
        
        for i in paramInt1 ... paramInt2 {
            for  j in 0 ..< paramArrayOfInt.count {
                
                if i > A.count - 1 || paramArrayOfInt[j] > A[i].count - 1 {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "Submatrix indices")
                }
                
                self.A[i][paramArrayOfInt[j]] = paramMatrix.get(paramInt1: i - paramInt1, paramInt2: j)
            }
        }
        
    }
    
    public func transpose() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.n, paramInt2: self.m)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                localMatrix.A[j][i] = self.A[i][j]
            }
        }
        return localMatrix
    }
    
    public func norm1() -> Double
    {
        var d1: Double = 0.0
        for i in 0 ..< self.n
        {
            var d2:Double = 0.0
            for  j in 0 ..< self.m {
                d2 += abs(self.A[j][i])
            }
            d1 = max(d1, d2)
        }
        return d1
    }
    
    public func norm2() -> Double
    {
        return SingularValueDecomposition(paramMatrix: self).norm2()
    }
    
    public func normInf() -> Double
    {
        var d1: Double = 0.0
        for i in 0 ..< self.m
        {
            var d2: Double = 0.0
            for j in 0 ..< self.n {
                d2 += abs(self.A[i][j])
            }
            d1 = max(d1, d2)
        }
        return d1
    }
    
    public func normF() -> Double
    {
        var d : Double = 0.0
        for  i in 0 ..< self.m {
            for j in 0 ..< self.n {
                d = hypot(d, self.A[i][j])
            }
        }
        return d
    }
    
    public func uminus() -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                localMatrix.A[i][j] = (-self.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func plus(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m  {
            for j in 0 ..< self.n {
                localMatrix.A[i][j] = (self.A[i][j] + paramMatrix.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func plusEquals(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                self.A[i][j] += paramMatrix.A[i][j]
            }
        }
        return self
    }
    
    public func minus(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                localMatrix.A[i][j] = (self.A[i][j] - paramMatrix.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func minusEquals(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                self.A[i][j] -= paramMatrix.A[i][j]
            }
        }
        return self
    }
    
    public func arrayTimes(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                localMatrix.A[i][j] = (self.A[i][j] * paramMatrix.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func arrayTimesEquals(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                self.A[i][j] *= paramMatrix.A[i][j]
            }
        }
        return self
    }
    
    public func arrayRightDivide(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for  i in 0 ..< self.m {
            for  j in 0 ..< self.n {
                localMatrix.A[i][j] = (self.A[i][j] / paramMatrix.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func arrayRightDivideEquals(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                self.A[i][j] /= paramMatrix.A[i][j]
            }
        }
        return self
    }
    
    public func arrayLeftDivide(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                localMatrix.A[i][j] = (paramMatrix.A[i][j] / self.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func arrayLeftDivideEquals(paramMatrix: Matrix) throws -> Matrix
    {
        try checkMatrixDimensions(paramMatrix: paramMatrix)
        
        for  i in 0 ..< self.m {
            for j in 0 ..< self.n {
                paramMatrix.A[i][j] /= self.A[i][j]
            }
        }
        return self
    }
    
    public func times(paramDouble: Double) -> Matrix
    {
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: self.n)
        for i in 0 ..< self.m {
            for j in  0 ..< self.n {
                localMatrix.A[i][j] = (paramDouble * self.A[i][j])
            }
        }
        return localMatrix
    }
    
    public func timesEquals(paramDouble: Double) -> Matrix
    {
        for i in 0 ..< self.m {
            for j in 0 ..< self.n {
                self.A[i][j] = (paramDouble * self.A[i][j])
            }
        }
        return self
    }
    
    public func times(paramMatrix: Matrix) throws -> Matrix
    {
        if (paramMatrix.m != self.n) {
            throw MatrixError.IllegalArgumentException(info: "Matrix inner dimensions must agree.")
        }
        let localMatrix = Matrix(paramInt1: self.m, paramInt2: paramMatrix.n)
        var arrayOfDouble1 = [Double].init(repeating: 0.0, count: self.n)
        for i in 0 ..< paramMatrix.n
        {
            for j in 0 ..< self.n {
                arrayOfDouble1[j] = paramMatrix.A[j][i]
            }
            for j in 0 ..< self.m
            {
                var arrayOfDouble2 = self.A[j]
                var d = 0.0
                for k in 0 ..< self.n {
                    d += arrayOfDouble2[k] * arrayOfDouble1[k]
                }
                localMatrix.A[j][i] = d
            }
        }
        return localMatrix
    }
    
    public func lu() -> LUDecomposition
    {
        return LUDecomposition(A: self)
    }
    
    public func qr() -> QRDecomposition
    {
        return QRDecomposition(paramMatrix: self)
    }
    
    public func chol() -> CholeskyDecomposition
    {
        return CholeskyDecomposition(paramMatrix: self)
    }
    
    public func svd() -> SingularValueDecomposition
    {
        return SingularValueDecomposition(paramMatrix: self)
    }
    
    public func eig() throws -> EigenvalueDecomposition
    {
        return try EigenvalueDecomposition(Arg: self)
    }
    
    public func solve(paramMatrix: Matrix) throws -> Matrix
    {
        return try self.m == self.n ? LUDecomposition(A: self).solve(B: paramMatrix) : QRDecomposition(paramMatrix: self).solve(paramMatrix: paramMatrix)
    }
    
    public func solveTranspose(paramMatrix: Matrix) throws -> Matrix
    {
        return try transpose().solve(paramMatrix: paramMatrix.transpose())
    }
    
    public func inverse() throws -> Matrix
    {
        return try self.solve(paramMatrix: .identity(paramInt1: self.m, paramInt2: self.m))
    }
    
    public func det() throws -> Double
    {
        return try LUDecomposition(A: self).det()
    }
    
    public  func rank() -> Int
    {
        return SingularValueDecomposition(paramMatrix: self).rank()
    }
    
    public func cond() -> Double
    {
        return SingularValueDecomposition(paramMatrix: self).cond()
    }
    
    public func trace() -> Double
    {
        var d: Double = 0.0
        for i in 0 ..< min(self.m, self.n) {
            d += self.A[i][i]
        }
        return d
    }
    
    public static func random(int paramInt1: Int, int paramInt2: Int) -> Matrix
    {
        let localMatrix = Matrix(paramInt1: paramInt1, paramInt2: paramInt2)
        for i in 0 ..< paramInt1 {
            for j in 0 ..< paramInt2 {
                localMatrix.A[i][j] = Double(1 + arc4random_uniform(10)) / 10.0
            }
        }
        return localMatrix
    }
    
    public static func identity(paramInt1: Int, paramInt2: Int) -> Matrix {
        let localMatrix = Matrix(paramInt1: paramInt1, paramInt2: paramInt2)
        for i in 0 ..< paramInt1 {
            for j in 0 ..< paramInt2 {
                localMatrix.A[i][j] = (i == j ? 1.0 : 0.0)
            }
        }
        return localMatrix
    }
    
    private func checkMatrixDimensions(paramMatrix: Matrix) throws
    {
        if ((paramMatrix.m != self.m) || (paramMatrix.n != self.n)) {
            throw MatrixError.IllegalArgumentException(info: "Matrix dimensions must agree.")
        }
    }
    
    //    public void print(int paramInt1, int paramInt2)
    //{
    //    print(new PrintWriter(System.out, true), paramInt1, paramInt2)
    //    }
    //
    //    public void print(PrintWriter paramPrintWriter, int paramInt1, int paramInt2)
    //{
    //    DecimalFormat localDecimalFormat = new DecimalFormat()
    //    localDecimalFormat.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US))
    //    localDecimalFormat.setMinimumIntegerDigits(1)
    //    localDecimalFormat.setMaximumFractionDigits(paramInt2)
    //    localDecimalFormat.setMinimumFractionDigits(paramInt2)
    //    localDecimalFormat.setGroupingUsed(false)
    //    print(paramPrintWriter, localDecimalFormat, paramInt1 + 2)
    //    }
    //
    //    public void print(NumberFormat paramNumberFormat, int paramInt)
    //{
    //    print(new PrintWriter(System.out, true), paramNumberFormat, paramInt)
    //    }
    //
    //    public func print(PrintWriter paramPrintWriter, NumberFormat paramNumberFormat, int paramInt)
    //    {
    //    paramPrintWriter.println()
    //    for (int i = 0 i < self.m i++)
    //    {
    //    for (int j = 0 j < self.n j++)
    //    {
    //    String str = paramNumberFormat.format(self.A[i][j])
    //    int k = Math.max(1, paramInt - str.length())
    //    for (int i1 = 0 i1 < k i1++) {
    //    paramPrintWriter.print(' ')
    //    }
    //    paramPrintWriter.print(str)
    //    }
    //    paramPrintWriter.println()
    //    }
    //    paramPrintWriter.println()
    //    }
    
    //    public static Matrix read(BufferedReader paramBufferedReader)
    //    throws IOException
    //    {
    //    StreamTokenizer localStreamTokenizer = new StreamTokenizer(paramBufferedReader)
    //
    //    localStreamTokenizer.resetSyntax()
    //    localStreamTokenizer.wordChars(0, 255)
    //    localStreamTokenizer.whitespaceChars(0, 32)
    //    localStreamTokenizer.eolIsSignificant(true)
    //    Vector localVector1 = new Vector()
    //    while (localStreamTokenizer.nextToken() == 10) {}
    //    if (localStreamTokenizer.ttype == -1) {
    //    throw new IOException("Unexpected EOF on matrix read.")
    //    }
    //    do
    //    {
    //    localVector1.addElement(Double.valueOf(localStreamTokenizer.sval))
    //    } while (localStreamTokenizer.nextToken() == -3)
    //    int i = localVector1.size()
    //    double[] arrayOfDouble = new double[i]
    //    for (int j = 0 j < i j++) {
    //    arrayOfDouble[j] = ((Double)localVector1.elementAt(j)).doubleValue()
    //    }
    //    Vector localVector2 = new Vector()
    //    localVector2.addElement(arrayOfDouble)
    //    while (localStreamTokenizer.nextToken() == -3)
    //    {
    //    localVector2.addElement(arrayOfDouble = new double[i])
    //    k = 0
    //    do
    //    {
    //    if (k >= i) {
    //    throw new IOException("Row " + localVector2.size() + " is too long.")
    //    }
    //    arrayOfDouble[(k++)] = Double.valueOf(localStreamTokenizer.sval).doubleValue()
    //    } while (localStreamTokenizer.nextToken() == -3)
    //    if (k < i) {
    //    throw new IOException("Row " + localVector2.size() + " is too short.")
    //    }
    //    }
    //    int k = localVector2.size()
    //    double[][] arrayOfDouble1 = new double[k][]
    //    localVector2.copyInto(arrayOfDouble1)
    //    return new Matrix(arrayOfDouble1)
    //    }
    
}
