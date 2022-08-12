#include <iostream>

#include <python3.8/Python.h>

int main()
{
    Py_Initialize();

    // 执行import语句,把当前路径加入路径中,为了找到readvtk.py
    PyRun_SimpleString("import os,sys");
    PyRun_SimpleString("sys.path.append('./')");
    PyRun_SimpleString("sys.path.append('../')");

    PyRun_SimpleString("import readvtk");
    PyRun_SimpleString("readvtk.transHexFile()");

    Py_Finalize();

    return 0;
}