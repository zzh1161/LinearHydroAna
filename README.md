# Linear Hydrostatic Analysis
## 文档说明
+ `include` 目录下包含所有的头文件
  + `info.hpp` 包含了一些宏定义
  
  + `LinSysSolver.hpp`是抽象基类，以CSR格式存储稀疏矩阵，用于实现刚度矩阵的装配和线性系统的求解
  
    `EigenLibSolver.hpp`和`CHOLMODSolver.hpp`是相应继承类，代表采用Eigen库或CHOLMOD库求解
    
  + `readTetMesh.hpp` 定义了读取四面体网格.msh文件的函数
  
    `readHexMesh.hpp` 定义了读取六面体网格.in文件的函数，其中.in文件由.vtk文件解析而来
  
    `readCondition.hpp` 定义了读取边界条件和其他信息的函数，对输入文件格式有要求
  
  + `shapeInterface.hpp` 定义网格抽象基类
  
    `tet4nodeEle.hpp` 是继承而来的4节点四面体类
  
    `hex8nodeEle.hpp` 是继承而来的8节点六面体类
  
  + `tet4nodeSolver.hpp` 将四面体网格模型的有限元分析过程封装起来
  
    `hex8nodeSolver.hpp` 将六面体网格模型的有限元分析过程封装起来
  
+ `tests` 目录是所有的测试样例

  + `tet_cube` 测试了一个8个顶点8个四面体的简单网格模型，并使用matlab程序验证正确性
  + `tet_bunny` 测试了bunny.msh四面体网格模型
  + `tet_armadillo` 测试了Armadillo219K.msh四面体网格模型，此模型一共有219k个节点，使用Eigen库求解大约需要30分钟，使用CHOLMOD库求解仅需5分钟
  + `hex_smallCube` 测试了一个小立方体模型，此模型由27个小六面体构成，并使用matlab验证了正确性
  + `hex_cube` 测试了更大的立方体模型，此模型由512个六面体、共729个节点构成

## 编译与运行

代码作者本地具有如下环境：

+ Ubuntu 20.04
+ g++(GCC) 11.2.0 (C++20)
+ Python 3.8.10
+ CMake 3.16.3
+ Eigen 3.3.7
+ SuiteSparse 5.7.1

编译：

+ 在hex测试样例中，执行

  ```bash
  mkdir build && cd build
  cmake ..
  make
  ```

  完成编译，执行 `./main` 运行程序

+ 在tet测试样例中，执行 `make` 完成编译，执行 `./main` 运行程序

