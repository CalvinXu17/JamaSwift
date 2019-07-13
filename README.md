# JamaSwift
Jama for Swift
Swift for Linear Algebra(Vector/Matrix).
It is translated from Jama(Java Matrix Package).
Thanks for Jama's source and I translate it to Swift.

# Usage
It's Class and Method are the same as Jama, so you can google Jama for it's usage.

# 欢迎加群905506227交流讨论
Jama Swift是基于Jama(Java Matrix Package)翻译成swift版本的一款可以实现线性代数向量、矩阵运算的库
感谢Jama团队开放的源代码
Jama Swift所有类与方法均与Jama一致，具体用法直接搜索Jama即可

不需要第三方库的支持，均使用Swift原生编写

# Demo
let aa:[[Double]] = [[10.0,20.0,30.0],[12.0,24.0,46.0],[34.0,58.0,12.0]]
do {
    let mat = try Matrix(paramArrayOfDouble: aa)
    print("cond:  \(mat.cond())") 
    print("rank:  \(mat.rank())") // 秩
    print("chol:  \(mat.chol().getL().getArray())")
    print("trace: \(mat.trace())") // 迹
    print("transpose:\(mat.transpose().getArray())") // 转置
    let r = try mat.det() // 行列式的值
    print("det:   \(r)")
    let c = try mat.lu().det()
    let d = try mat.inverse().getArray() // 逆矩阵
    print("inverse: \(d)")
    print("特征值: \(mat.eig().getD().getArray())") // 特征值
    print("特征向量: \(mat.eig().getV().getArray())") // 特征向量
    // mat.times(Matrix) 矩阵乘法
    } catch {
    print(error)
}
