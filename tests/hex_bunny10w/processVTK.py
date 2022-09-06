import os
import meshio
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(__file__))

def transHexFile():
    mesh = meshio.vtk.read("./input/bunny10w.vtk")
    with open("./input/bunny10w.in","w") as fout:
        fout.write(str(len(mesh.points))+"\n")
        for arr in mesh.points:
            fout.write(str(arr[0])+" "+str(arr[1])+" "+str(arr[2])+"\n")
        fout.write(str(len(mesh.cells[0]))+"\n")
        for ele in mesh.cells[0].data:
            for i in range(0,7):
                fout.write(str(ele[i])+" ")
            fout.write(str(ele[7])+"\n")

def outputHexVTK():
    mesh = meshio.vtk.read("./input/bunny10w.vtk")
    result_uvw = np.array(pd.DataFrame(columns=range(3)))
    with open('./output/result_uvw.txt','r') as fin:
        data = fin.readlines()
        for line in data:
            get_numbers = line.split()
            get_numbers[0] = float(get_numbers[0])
            get_numbers[1] = float(get_numbers[1])
            get_numbers[2] = float(get_numbers[2])
            get_numbers = np.array(get_numbers)
            result_uvw = np.append(result_uvw, [get_numbers], axis=0)
    points = mesh.points + result_uvw
    cells = mesh.cells
    points = np.array(points, dtype=np.float64)
    ano_mesh = meshio.Mesh(points, cells)
    ano_mesh.write('./output/result_bunny10w.vtk')
