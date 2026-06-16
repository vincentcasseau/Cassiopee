import Converter.PyTree as C
import Converter.Internal as Internal

import os


def createVTUFile():
    import pyvista as pv
    import numpy as np

    points = np.array([
        [0.,0.,0.],
        [1.,0.,0.],
        [0.,1.,0.]
    ])

    cells = np.array([3,0,1,2])  # 3 points in triangle
    celltypes = np.array([5])    # VTK_TRIANGLE

    grid = pv.UnstructuredGrid(cells, celltypes, points)
    grid.save("mesh.vtu")

def readHLAero():
    folder = "/stck/peter/PUBLICC/HighLiftAeroML_LHC017_AoA_6"


if __name__ == "__main__":
    for eltType in ["BAR", "TRI", "QUAD", "TETRA", "PYRA", "PENTA", "HEXA"]:
        a = C.convertFile2PyTree(f"./vtu/mesh_{eltType}.vtu", api=1)
        C.convertPyTree2File(a, f"./vtu/mesh_{eltType}.cgns")

