import os
import meshio
os.chdir(os.path.dirname(__file__))

def transHexFile():
    mesh = meshio.vtk.read("./input/cube.vtk")
    with open("./input/cube.in","w") as fout:
        fout.write(str(len(mesh.points))+"\n")
        for arr in mesh.points:
            fout.write(str(arr[0])+" "+str(arr[1])+" "+str(arr[2])+"\n")
        fout.write(str(len(mesh.cells[0]))+"\n")
        for ele in mesh.cells[0].data:
            for i in range(0,7):
                fout.write(str(ele[i])+" ")
            fout.write(str(ele[7])+"\n")

def readHexFile(TV, TT):
    mesh = meshio.vtk.read("./input/cube.vtk")
    TV = mesh.points
    TT = mesh.cells[0].data