//
//  MatrixError.swift
//  Sama
//
//  Created by Calvin on 2019/7/12.
//  Copyright © 2019 Calvin. All rights reserved.
//

import Foundation

enum MatrixError: Error {
    case IllegalArgumentException(info: String)
    case RuntimeException(info: String)
    case ArrayIndexOutOfBoundsException(info: String)
}
