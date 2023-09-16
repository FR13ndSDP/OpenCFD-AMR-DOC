# 1. 射流燃烧

> 开启化学反应，单层网格, 160x64x64, 计算100步

## 1.1 RTX 3080
Run Time total        = **79.97765838**
Run Time init         = 0.00442036
Run Time advance      = 79.97317812

## 1.2 Tesla V-100
Run Time total        = **20.14456375**
Run Time init         = 0.002265513
Run Time advance      = 20.14219416

## 1.3 A10
Run Time total        = **86.84865471**
Run Time init         = 0.004822282
Run Time advance      = 86.84375563

## 1.4 A10 1.6x grid
Run Time total        = **139.4827673**
Run Time init         = 0.006830545
Run Time advance      = 139.4758513

## 1.5 A10X4 (mpirun -n 4)
Run Time total        = **30.5222507**
Run Time init         = 0.002297403
Run Time advance      = 30.51993597

## 1.6 A10X4 (mpirun -n 4) 1.6x grid
Run Time total        = **39.89737012**
Run Time init         = 0.002684094
Run Time advance      = 39.89466335

## 1.7 i9-12900k (mpirun -n 8)
Run Time total        = **325.6690344**
Run Time init         = 0.064384584
Run Time advance      = 325.6030875