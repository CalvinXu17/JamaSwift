//
//  EigenvalueDecomposition.swift
//  Sama
//
//  Created by Calvin on 2019/7/13. 
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

public class EigenvalueDecomposition: Codable {
    
    /* ------------------------
     Class variables
     * ------------------------ */
    
    /** Row and column dimension (square matrix).
     @serial matrix dimension.
     */
    private var n: Int
    
    /** Symmetry flag.
     @serial internal symmetry flag.
     */
    private var issymmetric: Bool
    
    /** Arrays for internal storage of eigenvalues.
     @serial internal storage of eigenvalues.
     */
    private var d, e: [Double]
    
    /** Array for internal storage of eigenvectors.
     @serial internal storage of eigenvectors.
     */
    private var V: [[Double]]
    
    /** Array for internal storage of nonsymmetric Hessenberg form.
     @serial internal storage of nonsymmetric Hessenberg form.
     */
    private var H: [[Double]] = [[]]
    
    /** Working storage for nonsymmetric algorithm.
     @serial working storage for nonsymmetric algorithm.
     */
    private var ort: [Double] = []
    
    /* ------------------------
     Private Methods
     * ------------------------ */
    
    // Symmetric Householder reduction to tridiagonal form.
    
    private func tred2() {
        
        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        
        for j in 0 ..< n {
            d[j] = V[n - 1][j]
        }
        
        // Householder reduction to tridiagonal form.
        for i in stride(from: n - 1, through: 1, by: -1) {
            
            // Scale to avoid under/overflow.
            
            var scale: Double = 0.0
            var h: Double = 0.0
            for k in stride(from: 0, to: i, by: 1) {
                scale = scale + abs(d[k])
            }
            
            if (scale == 0.0) {
                e[i] = d[i-1]
                for j in stride(from: 0, to: i, by: 1) {
                    d[j] = V[i-1][j]
                    V[i][j] = 0.0
                    V[j][i] = 0.0
                }
            } else {
                
                // Generate Householder vector.
                
                for k in stride(from: 0, to: i, by: 1) {
                    d[k] /= scale
                    h += d[k] * d[k]
                }
                var f = d[i-1]
                var g = sqrt(h)
                if (f > 0) {
                    g = -g
                }
                e[i] = scale * g
                h = h - f * g
                d[i-1] = f - g
                for j in 0 ..< i {
                    e[j] = 0.0
                }
                
                // Apply similarity transformation to remaining columns.
                
                for j in stride(from: 0, to: i, by: 1) {
                    f = d[j]
                    V[j][i] = f
                    g = e[j] + V[j][j] * f
                    for k in stride(from: j + 1, through: i - 1, by: 1) {
                        g += V[k][j] * d[k]
                        e[k] += V[k][j] * f
                    }
                    e[j] = g
                }
                f = 0.0
                for j in stride(from: 0, to: i, by: 1) {
                    e[j] /= h
                    f += e[j] * d[j]
                }
                let hh: Double = f / (h + h)
                for j in 0 ..< i {
                    e[j] -= hh * d[j]
                }
                for j in stride(from: 0, to: i, by: 1) {
                    f = d[j]
                    g = e[j]
                    for k in stride(from: j, through: i - 1, by: 1) {
                        V[k][j] -= (f * e[k] + g * d[k])
                    }
                    d[j] = V[i - 1][j]
                    V[i][j] = 0.0
                }
            }
            d[i] = h
        }
        
        // Accumulate transformations.
        
        for i in 0 ..< n - 1 {
            V[n-1][i] = V[i][i]
            V[i][i] = 1.0
            let h = d[i + 1]
            if (h != 0.0) {
                for k in stride(from: 0, through: i, by: 1) {
                    d[k] = V[k][i+1] / h
                }
                for j in stride(from: 0, through: i, by: 1) {
                    var g: Double = 0.0
                    for k in 0 ... i {
                        g += V[k][i+1] * V[k][j]
                    }
                    for k in stride(from: 0, through: i, by: 1) {
                        V[k][j] -= g * d[k]
                    }
                }
            }
            for k in stride(from: 0, through: i, by: 1) {
                V[k][i + 1] = 0.0
            }
        }
        for j in stride(from: 0, to: n, by: 1) {
            d[j] = V[n - 1][j]
            V[n - 1][j] = 0.0
        }
        V[n - 1][n - 1] = 1.0
        e[0] = 0.0
    }
    
    // Symmetric tridiagonal QL algorithm.
    
    private func tql2() {
        
        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        
        for i in stride(from: 1, to: n, by: 1) {
            e[i-1] = e[i]
        }
        e[n-1] = 0.0
        
        var f: Double = 0.0
        var tst1: Double = 0.0
        let eps: Double = pow(2.0,-52.0)
        for l in 0 ..< n {
            
            // Find small subdiagonal element
            
            tst1 = max(tst1,abs(d[l]) + abs(e[l]))
            var m: Int = l
            while (m < n) {
                if (abs(e[m]) <= eps * tst1) {
                    break
                }
                m += 1
            }
            
            // If m == l, d[l] is an eigenvalue,
            // otherwise, iterate.
            
            if (m > l) {
                var iter: Int = 0
                repeat {
                    iter = iter + 1  // (Could check iteration count here.)
                    
                    // Compute implicit shift
                    
                    var g: Double = d[l]
                    var p: Double = (d[l+1] - g) / (2.0 * e[l])
                    var r: Double = hypot(p,1.0)
                    if (p < 0) {
                        r = -r
                    }
                    d[l] = e[l] / (p + r)
                    d[l+1] = e[l] * (p + r)
                    let dl1 = d[l+1]
                    var h = g - d[l]
                    for i in l + 2 ..< n {
                        d[i] -= h
                    }
                    f = f + h
                    
                    // Implicit QL transformation.
                    
                    p = d[m]
                    var c: Double = 1.0
                    var c2: Double = c
                    var c3: Double = c
                    let el1 = e[l + 1]
                    var s: Double = 0.0
                    var s2: Double = 0.0
                    for i in stride(from: m - 1, through: l, by: -1) {
                        c3 = c2
                        c2 = c
                        s2 = s
                        g = c * e[i]
                        h = c * p
                        r = hypot(p,e[i])
                        e[i+1] = s * r
                        s = e[i] / r
                        c = p / r
                        p = c * d[i] - s * g
                        d[i+1] = h + s * (c * g + s * d[i])
                        
                        // Accumulate transformation.
                        
                        for k in stride(from: 0, to: n, by: 1) {
                            h = V[k][i+1]
                            V[k][i+1] = s * V[k][i] + c * h
                            V[k][i] = c * V[k][i] - s * h
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1
                    e[l] = s * p
                    d[l] = c * p
                    
                    // Check for convergence.
                    
                } while (abs(e[l]) > eps * tst1)
            }
            d[l] = d[l] + f
            e[l] = 0.0
        }
        
        // Sort eigenvalues and corresponding vectors.
        
        for i in stride(from: 0, to: n - 1, by: 1) {
            var k: Int = i
            var p: Double = d[i]
            for j in stride(from: i + 1, to: n, by: 1) {
                if (d[j] < p) {
                    k = j
                    p = d[j]
                }
            }
            if (k != i) {
                d[k] = d[i]
                d[i] = p
                for j in stride(from: 0, to: n, by: 1) {
                    p = V[j][i]
                    V[j][i] = V[j][k]
                    V[j][k] = p
                }
            }
        }
    }
    
    // Nonsymmetric reduction to Hessenberg form.
    
    private func orthes() {
        
        //  This is derived from the Algol procedures orthes and ortran,
        //  by Martin and Wilkinson, Handbook for Auto. Comp.,
        //  Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutines in EISPACK.
        
        let low = 0
        let high = n - 1
        
        for m in stride(from: low + 1, through: high - 1, by: 1) {
            
            // Scale column.
            
            var scale: Double = 0.0
            for i in stride(from: m, through: high, by: 1) {
                scale = scale + abs(H[i][m-1])
            }
            if (scale != 0.0) {
                
                // Compute Householder transformation.
                
                var h: Double = 0.0
                for i in stride(from: high, through: m, by: -1) {
                    ort[i] = H[i][m-1] / scale
                    h += ort[i] * ort[i]
                }
                var g: Double = sqrt(h)
                if (ort[m] > 0) {
                    g = -g
                }
                h = h - ort[m] * g
                ort[m] = ort[m] - g
                
                // Apply Householder similarity transformation
                // H = (I-u*u'/h)*H*(I-u*u')/h)
                
                for j in m ..< n {
                    var f: Double = 0.0
                    for i in stride(from: high, through: m, by: -1) {
                        f += ort[i]*H[i][j]
                    }
                    f = f/h
                    for i in stride(from: m, through: high, by: 1) {
                        H[i][j] -= f * ort[i]
                    }
                }
                
                for i in stride(from: 0, through: high, by: 1) {
                    var f: Double = 0.0
                    for j in stride(from: high, through: m, by: -1) {
                        f += ort[j] * H[i][j]
                    }
                    f = f/h
                    for j in stride(from: m, through: high, by: 1) {
                        H[i][j] -= f * ort[j]
                    }
                }
                ort[m] = scale * ort[m]
                H[m][m-1] = scale * g
            }
        }
        
        // Accumulate transformations (Algol's ortran).
        
        for i in 0 ..< n {
            for  j in 0 ..< n{
                V[i][j] = (i == j ? 1.0 : 0.0)
            }
        }
        
        for  m in stride(from: high - 1, through: low + 1, by: -1) {
            if (H[m][m-1] != 0.0) {
                for i in stride(from: m + 1, through: high, by: 1) {
                    ort[i] = H[i][m-1]
                }
                for j in stride(from: m, through: high, by: 1) {
                    var g: Double = 0.0
                    for i in stride(from: m, through: high, by: 1) {
                        g += ort[i] * V[i][j]
                    }
                    // Double division avoids possible underflow
                    g = (g / ort[m]) / H[m][m-1]
                    for i in stride(from: m, through: high, by: 1) {
                        V[i][j] += g * ort[i]
                    }
                }
            }
        }
    }
    
    
    // Complex scalar division.
    
    private var cdivr: Double = 0.0
    private var cdivi: Double = 0.0
    
    private func cdiv(xr: Double, xi: Double, yr: Double, yi: Double) {
        var r: Double = 0.0
        var d: Double = 0.0
        if (abs(yr) > abs(yi)) {
            r = yi / yr
            d = yr + r * yi
            cdivr = (xr + r * xi) / d
            cdivi = (xi - r * xr) / d
        } else {
            r = yr / yi
            d = yi + r * yr
            cdivr = (r * xr + xi) / d
            cdivi = (r * xi - xr) / d
        }
    }
    
    
    // Nonsymmetric reduction from Hessenberg to real Schur form.
    
    private func hqr2 () {
        
        //  This is derived from the Algol procedure hqr2,
        //  by Martin and Wilkinson, Handbook for Auto. Comp.,
        //  Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        
        // Initialize
        
        let nn = self.n
        var n = nn - 1
        let low = 0
        let high = nn - 1
        let eps: Double = pow(2.0, -52.0)
        var exshift: Double = 0.0
        var p: Double = 0.0
        var q: Double = 0.0
        var r: Double = 0.0
        var s: Double = 0.0
        var z: Double = 0.0
        var t: Double = 0.0
        var w: Double = 0.0
        var x: Double = 0.0
        var y: Double = 0.0
        
        // Store roots isolated by balanc and compute matrix norm
        
        var norm: Double = 0.0
        for i in stride(from: 0, to: nn, by: 1) {
            if (i < low || i > high) {
                d[i] = H[i][i]
                e[i] = 0.0
            }
            for j in max(i - 1,0) ..< nn {
                norm = norm + abs(H[i][j])
            }
        }
        
        // Outer loop over eigenvalue index
        
        var iter = 0
        while (n >= low) {
            
            // Look for single small sub-diagonal element
            
            var l = n
            while (l > low) {
                s = abs(H[l-1][l-1]) + abs(H[l][l])
                if (s == 0.0) {
                    s = norm
                }
                if (abs(H[l][l-1]) < eps * s) {
                    break
                }
                l -= 1
            }
            
            // Check for convergence
            // One root found
            
            if (l == n) {
                H[n][n] = H[n][n] + exshift
                d[n] = H[n][n]
                e[n] = 0.0
                n -= 1
                iter = 0
                
                // Two roots found
                
            } else if (l == n-1) {
                w = H[n][n-1] * H[n-1][n]
                p = (H[n-1][n-1] - H[n][n]) / 2.0
                q = p * p + w
                z = sqrt(abs(q))
                H[n][n] = H[n][n] + exshift
                H[n-1][n-1] = H[n-1][n-1] + exshift
                x = H[n][n]
                
                // Real pair
                
                if (q >= 0) {
                    if (p >= 0) {
                        z = p + z
                    } else {
                        z = p - z
                    }
                    d[n-1] = x + z
                    d[n] = d[n-1]
                    if (z != 0.0) {
                        d[n] = x - w / z
                    }
                    e[n-1] = 0.0
                    e[n] = 0.0
                    x = H[n][n-1]
                    s = abs(x) + abs(z)
                    p = x / s
                    q = z / s
                    r = sqrt(p * p+q * q)
                    p = p / r
                    q = q / r
                    
                    // Row modification
                    
                    for j in stride(from: n - 1, to: nn, by: 1) {
                        z = H[n-1][j]
                        H[n-1][j] = q * z + p * H[n][j]
                        H[n][j] = q * H[n][j] - p * z
                    }
                    
                    // Column modification
                    
                    for i in stride(from: 0, through: n, by: 1) {
                        z = H[i][n-1]
                        H[i][n-1] = q * z + p * H[i][n]
                        H[i][n] = q * H[i][n] - p * z
                    }
                    
                    // Accumulate transformations
                    
                    for i in stride(from: low, through: high, by: 1) {
                        z = V[i][n-1]
                        V[i][n-1] = q * z + p * V[i][n]
                        V[i][n] = q * V[i][n] - p * z
                    }
                    
                    // Complex pair
                    
                } else {
                    d[n-1] = x + p
                    d[n] = x + p
                    e[n-1] = z
                    e[n] = -z
                }
                n = n - 2
                iter = 0
                
                // No convergence yet
                
            } else {
                
                // Form shift
                
                x = H[n][n]
                y = 0.0
                w = 0.0
                if (l < n) {
                    y = H[n-1][n-1]
                    w = H[n][n-1] * H[n-1][n]
                }
                
                // Wilkinson's original ad hoc shift
                
                if (iter == 10) {
                    exshift += x
                    for i in stride(from: low, through: n, by: 1) {
                        H[i][i] -= x
                    }
                    s = abs(H[n][n-1]) + abs(H[n-1][n-2])
                    x = 0.75 * s
                    y = 0.75 * s
                    w = -0.4375 * s * s
                }
                
                // MATLAB's new ad hoc shift
                
                if (iter == 30) {
                    s = (y - x) / 2.0
                    s = s * s + w
                    if (s > 0) {
                        s = sqrt(s)
                        if (y < x) {
                            s = -s
                        }
                        s = x - w / ((y - x) / 2.0 + s)
                        for i in low ... n {
                            H[i][i] -= s
                        }
                        exshift += s
                        x = 0.964
                        y = 0.964
                        w = 0.964
                    }
                }
                
                iter = iter + 1   // (Could check iteration count here.)
                
                // Look for two consecutive small sub-diagonal elements
                
                var m: Int = n - 2
                while (m >= l) {
                    z = H[m][m]
                    r = x - z
                    s = y - z
                    p = (r * s - w) / H[m+1][m] + H[m][m+1]
                    q = H[m+1][m+1] - z - r - s
                    r = H[m+2][m+1]
                    s = abs(p) + abs(q) + abs(r)
                    p = p / s
                    q = q / s
                    r = r / s
                    if (m == l) {
                        break
                    }
                    if (abs(H[m][m-1]) * (abs(q) + abs(r)) < eps * (abs(p) * (abs(H[m-1][m-1]) + abs(z) + abs(H[m+1][m+1])))) {
                        break
                    }
                    m -= 1
                }
                
                for i in m + 2 ... n {
                    H[i][i-2] = 0.0
                    if (i > m+2) {
                        H[i][i-3] = 0.0
                    }
                }
                
                // Double QR step involving rows l:n and columns m:n
                
                
                for k in stride(from: m, through: n - 1, by: 1) {
                    let notlast: Bool = (k != n - 1)
                    if (k != m) {
                        p = H[k][k-1]
                        q = H[k+1][k-1]
                        r = (notlast ? H[k+2][k-1] : 0.0)
                        x = abs(p) + abs(q) + abs(r)
                        if (x == 0.0) {
                            continue
                        }
                        p = p / x
                        q = q / x
                        r = r / x
                    }
                    
                    s = sqrt(p * p + q * q + r * r)
                    if (p < 0) {
                        s = -s
                    }
                    if (s != 0) {
                        if (k != m) {
                            H[k][k-1] = -s * x
                        } else if (l != m) {
                            H[k][k-1] = -H[k][k-1]
                        }
                        p = p + s
                        x = p / s
                        y = q / s
                        z = r / s
                        q = q / p
                        r = r / p
                        
                        // Row modification
                        
                        for j in stride(from: k, to: nn, by: 1) {
                            p = H[k][j] + q * H[k+1][j]
                            if (notlast) {
                                p = p + r * H[k+2][j]
                                H[k+2][j] = H[k+2][j] - p * z
                            }
                            H[k][j] = H[k][j] - p * x
                            H[k+1][j] = H[k+1][j] - p * y
                        }
                        
                        // Column modification
                        for  i in stride(from: 0, through: min(n, k + 3), by: 1) {
                            p = x * H[i][k] + y * H[i][k+1]
                            if (notlast) {
                                p = p + z * H[i][k+2]
                                H[i][k+2] = H[i][k+2] - p * r
                            }
                            H[i][k] = H[i][k] - p
                            H[i][k+1] = H[i][k+1] - p * q
                        }
                        
                        // Accumulate transformations
                        
                        for i in stride(from: low, through: high, by: 1) {
                            p = x * V[i][k] + y * V[i][k+1]
                            if (notlast) {
                                p = p + z * V[i][k+2]
                                V[i][k+2] = V[i][k+2] - p * r
                            }
                            V[i][k] = V[i][k] - p
                            V[i][k+1] = V[i][k+1] - p * q
                        }
                    }  // (s != 0)
                }  // k loop
            }  // check convergence
        }  // while (n >= low)
        
        // Backsubstitute to find vectors of upper triangular form
        
        if (norm == 0.0) {
            return
        }
        
        for n in stride(from: nn - 1, through: 0, by: -1) {
            p = d[n]
            q = e[n]
            
            // Real vector
            
            if (q == 0) {
                var l = n
                H[n][n] = 1.0
                if n - 1 >= 0 {
                    for i in stride(from: n - 1, through: 0, by: -1) {
                        w = H[i][i] - p
                        r = 0.0
                        for j in stride(from: l, through: n, by: 1) {
                            r = r + H[i][j] * H[j][n]
                        }
                        if (e[i] < 0.0) {
                            z = w
                            s = r
                        } else {
                            l = i
                            if (e[i] == 0.0) {
                                if (w != 0.0) {
                                    H[i][n] = -r / w
                                } else {
                                    H[i][n] = -r / (eps * norm)
                                }
                                
                                // Solve real equations
                                
                            } else {
                                x = H[i][i+1]
                                y = H[i+1][i]
                                q = (d[i] - p) * (d[i] - p) + e[i] * e[i]
                                t = (x * s - z * r) / q
                                H[i][n] = t
                                if (abs(x) > abs(z)) {
                                    H[i+1][n] = (-r - w * t) / x
                                } else {
                                    H[i+1][n] = (-s - y * t) / z
                                }
                            }
                            
                            // Overflow control
                            
                            t = abs(H[i][n])
                            if ((eps * t) * t > 1) {
                                for j in stride(from: i, through: n, by: 1) {
                                    H[j][n] = H[j][n] / t
                                }
                            }
                        }
                    }
                }
                
                // Complex vector
                
            } else if (q < 0) {
                var l = n - 1
                
                // Last vector component imaginary so matrix is triangular
                
                if (abs(H[n][n-1]) > abs(H[n-1][n])) {
                    H[n-1][n-1] = q / H[n][n-1]
                    H[n-1][n] = -(H[n][n] - p) / H[n][n-1]
                } else {
                    cdiv(xr: 0.0,xi: -H[n-1][n],yr: H[n-1][n-1]-p,yi: q)
                    H[n-1][n-1] = cdivr
                    H[n-1][n] = cdivi
                }
                H[n][n-1] = 0.0
                H[n][n] = 1.0
                for i in stride(from: n - 2, through: 0, by: -1) {
                    var ra: Double = 0.0
                    var sa: Double = 0.0
                    var vr: Double = 0.0
                    var vi: Double = 0.0
                    for j in stride(from: l, through: n, by: 1) {
                        ra = ra + H[i][j] * H[j][n-1]
                        sa = sa + H[i][j] * H[j][n]
                    }
                    w = H[i][i] - p
                    
                    if (e[i] < 0.0) {
                        z = w
                        r = ra
                        s = sa
                    } else {
                        l = i
                        if (e[i] == 0) {
                            cdiv(xr: -ra,xi: -sa,yr: w,yi: q)
                            H[i][n-1] = cdivr
                            H[i][n] = cdivi
                        } else {
                            
                            // Solve complex equations
                            
                            x = H[i][i+1]
                            y = H[i+1][i]
                            vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q
                            vi = (d[i] - p) * 2.0 * q
                            if (vr == 0.0 && vi == 0.0) {
                                vr = eps * norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z))
                            }
                            cdiv(xr: x * r - z * ra + q * sa, xi: x * s - z * sa - q * ra, yr: vr, yi: vi)
                            H[i][n-1] = cdivr
                            H[i][n] = cdivi
                            if (abs(x) > (abs(z) + abs(q))) {
                                H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x
                                H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x
                            } else {
                                cdiv(xr: -r - y * H[i][n-1], xi: -s - y * H[i][n], yr: z, yi: q)
                                H[i+1][n-1] = cdivr
                                H[i+1][n] = cdivi
                            }
                        }
                        
                        // Overflow control
                        
                        t = max(abs(H[i][n-1]), abs(H[i][n]))
                        if ((eps * t) * t > 1) {
                            for j in stride(from: i, through: n, by: 1) {
                                H[j][n-1] = H[j][n-1] / t
                                H[j][n] = H[j][n] / t
                            }
                        }
                    }
                }
            }
        }
        
        // Vectors of isolated roots
        
        for i in stride(from: 0, to: nn, by: 1) {
            if (i < low || i > high) {
                for j in stride(from: i, to: nn, by: 1) {
                    V[i][j] = H[i][j]
                }
            }
        }
        
        // Back transformation to get eigenvectors of original matrix
        
        for j in stride(from: nn - 1, through: low, by: -1) {
            for i in stride(from: low, through: high, by: 1) {
                z = 0.0
                for k in stride(from: low, through: min(j,high), by: 1) {
                    z = z + V[i][k] * H[k][j]
                }
                V[i][j] = z
            }
        }
    }
    
    
    /* ------------------------
     Constructor
     * ------------------------ */
    
    /** Check for symmetry, then construct the eigenvalue decomposition
     Structure to access D and V.
     @param Arg    Square matrix
     */
    
    public init (Arg: Matrix) throws {
        self.n = Arg.getColumnDimension()
        self.V = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: self.n), count: self.n)
        self.d = [Double].init(repeating: 0.0, count: self.n)
        self.e = [Double].init(repeating: 0.0, count: self.n)
        self.issymmetric = true
        
        
        var jj = 0
        var ii = 0
        
        while jj < n && (issymmetric == true) {
            ii = 0
            while ii < n && (issymmetric == true) {
                if ii < Arg.A.count && jj < Arg.A[ii].count && jj < Arg.A.count && ii < Arg.A[jj].count {
                    self.issymmetric = (Arg.A[ii][jj] == Arg.A[jj][ii])
                } else {
                    throw MatrixError.ArrayIndexOutOfBoundsException(info: "index outof range")
                }
                ii += 1
            }
            jj += 1
        }
        
        //        for j = 0 (j < n) && issymmetric j++) {
        //            for (int i = 0 (i < n) & issymmetric i++) {
        //                issymmetric = (A[i][j] == A[j][i])
        //            }
        //        }
        
        if self.issymmetric == true {
            for i in 0 ..< n {
                for j in 0 ..< n {
                    if i < Arg.A.count && j < Arg.A[i].count {
                        V[i][j] = Arg.A[i][j]
                    } else {
                        throw MatrixError.ArrayIndexOutOfBoundsException(info: "index outof range")
                    }
                }
            }
            
            // Tridiagonalize.
            tred2()
            
            // Diagonalize.
            tql2()
            
        } else {
            self.H = [[Double]].init(repeating: [Double].init(repeating: 0.0, count: self.n), count: self.n)
            self.ort = [Double].init(repeating: 0.0, count: self.n)
            
            for j in 0 ..< n {
                for i in 0 ..< n {
                    H[i][j] = Arg.A[i][j]
                }
            }
            
            // Reduce to Hessenberg form.
            self.orthes()
            
            // Reduce Hessenberg to real Schur form.
            self.hqr2()
        }
    }
    
    /* ------------------------
     Public Methods
     * ------------------------ */
    
    /** Return the eigenvector matrix
     @return     V
     */
    
    public func getV() -> Matrix {
        return Matrix(paramArrayOfDouble: V,paramInt1: n,paramInt2: n)
    }
    
    /** Return the real parts of the eigenvalues
     @return     real(diag(D))
     */
    
    public func getRealEigenvalues() -> [Double] {
        return d
    }
    
    /** Return the imaginary parts of the eigenvalues
     @return     imag(diag(D))
     */
    
    public func getImagEigenvalues() -> [Double] {
        return e
    }
    
    /** Return the block diagonal eigenvalue matrix
     @return     D
     */
    
    public func getD() -> Matrix {
        let X = Matrix(paramInt1: n,paramInt2: n)
        for i in 0 ..< n {
            for j in 0 ..< n {
                X.A[i][j] = 0.0
            }
            X.A[i][i] = d[i]
            if (e[i] > 0) {
                X.A[i][i+1] = e[i]
            } else if (e[i] < 0) {
                X.A[i][i-1] = e[i]
            }
        }
        return X
    }
}
