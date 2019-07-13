//
//  ViewController.swift
//  test
//
//  Created by Calvin on 2019/7/13.
//  Copyright © 2019 Calvin. All rights reserved.
//

import UIKit
import Sama

class ViewController: UIViewController {

    override func viewDidLoad() {
        super.viewDidLoad()
        let aa:[[Double]] = [[10.0,20.0,30.0],[12.0,24.0,46.0],[34.0,58.0,12.0]]
        //let aa:[[Double]] = [[1.0,2.0],[3.0,4.0]]
        print("sss")
        do {
            let mat = try Matrix(paramArrayOfDouble: aa)
            print("cond:  \(mat.cond())")
            print("rank:  \(mat.rank())")
            print("chol:  \(mat.chol().getL().getArray())")
            print("trace: \(mat.trace())")
            print("transpose:\(mat.transpose().getArray())")
            
            let r = try mat.det()
            print("det:   \(r)")
            let c = try mat.lu().det()
            print("lu + det: \(c)")
            let d = try mat.inverse().getArray()
            print("inverse: \(d)")
            print("特征值: \(mat.eig().getD().getArray())")
            print("特征向量: \(mat.eig().getV().getArray())")
            print("hello")
            
        } catch {
            print(error)
        }
        do {
            let a1 = try Matrix(paramArrayOfDouble: [[1,2],[2,3]])
            let a2 = try Matrix(paramArrayOfDouble: [[2,4],[5,7]])
            let a3 = try a1.plus(paramMatrix: a2)  //times(paramMatrix: a2)
            print("times: \(a3.getArray())")
        }catch {
            print(error)
        }
        // Do any additional setup after loading the view.
    }


}

