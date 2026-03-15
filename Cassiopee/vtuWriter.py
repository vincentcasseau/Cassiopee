import pyvista as pv
import numpy as np

def createVTUFile_1D():
    # 3 lines (BAR_2, VTK_LINE=3)
    points = np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [2., 0., 0.],
        [3., 0., 0.],
    ])
    cells = np.array([
        2, 0, 1,
        2, 1, 2,
        2, 2, 3,
    ])
    celltypes = np.array([3, 3, 3])
    grid = pv.UnstructuredGrid(cells, celltypes, points)
    grid.save("mesh_1d.vtu")
    print("mesh_1d.vtu: 3 lines, 4 points")

def createVTUFile_2D():
    # mix of triangles (5) and quads (9)
    points = np.array([
        [0., 0., 0.],  # 0
        [1., 0., 0.],  # 1
        [2., 0., 0.],  # 2
        [3., 0., 0.],  # 3
        [0., 1., 0.],  # 4
        [1., 1., 0.],  # 5
        [2., 1., 0.],  # 6
        [3., 1., 0.],  # 7
        [0., 2., 0.],  # 8
        [1., 2., 0.],  # 9
        [2., 2., 0.],  # 10
    ])
    cells = np.array([
        3, 0, 1, 4,       # TRI
        3, 1, 5, 4,       # TRI
        4, 1, 2, 6, 5,    # QUAD
        4, 2, 3, 7, 6,    # QUAD
        3, 4, 5, 8,       # TRI
        4, 5, 6, 10, 9,   # QUAD
        3, 5, 9, 8,       # TRI
    ])
    celltypes = np.array([5, 5, 9, 9, 5, 9, 5])
    grid = pv.UnstructuredGrid(cells, celltypes, points)
    grid.save("mesh_2d.vtu")
    print("mesh_2d.vtu: 4 tris + 3 quads, 11 points")

def createVTUFile_3D():
    # mix of tetra (10), pyramid (14), wedge (13), hexa (12)
    points = np.array([
      # bottom layer z=0
      [0., 0., 0.],  # 0
      [1., 0., 0.],  # 1
      [2., 0., 0.],  # 2
      [3., 0., 0.],  # 3
      [0., 1., 0.],  # 4
      [1., 1., 0.],  # 5
      [2., 1., 0.],  # 6
      [3., 1., 0.],  # 7
      # top layer z=1
      [0., 0., 1.],  # 8
      [1., 0., 1.],  # 9
      [2., 0., 1.],  # 10
      [3., 0., 1.],  # 11
      [0., 1., 1.],  # 12
      [1., 1., 1.],  # 13
      [2., 1., 1.],  # 14
      [3., 1., 1.],  # 15
      # top layer z=2 (for second hexa)
      [0., 0., 2.],  # 16
      [1., 0., 2.],  # 17
      [0., 1., 2.],  # 18
      [1., 1., 2.],  # 19
      # apex for pyramid
      [1.5, 0.5, 1.5], # 20
    ])
    cells = np.array([
      # TETRA (4 nodes) x2
      4,  0,  1,  4,  8,
      4,  1,  5,  4,  9,
      # PYRAMID (5 nodes): quad base [10,11,15,14] + apex 20
      5, 10, 11, 15, 14, 20,
      # WEDGE/PENTA (6 nodes): triangles [0,1,4] and [8,9,12]
      6,  0,  1,  4,  8,  9, 12,
      # HEXA (8 nodes) x2
      8,  1,  2,  6,  5,  9, 10, 14, 13,
      8,  8,  9, 13, 12, 16, 17, 19, 18,
    ])
    celltypes = np.array([10, 10, 14, 13, 12, 12])
    grid = pv.UnstructuredGrid(cells, celltypes, points)
    grid.save("mesh_3d.vtu")
    print("mesh_3d.vtu: 2 tetras + 1 pyra + 1 wedge + 2 hexa, 21 points")
    
    grid = pv.read("mesh_3d.vtu")
    print("points:\n", grid.points)
    print("cells:\n", grid.cells)
    print("celltypes:\n", grid.celltypes)

if __name__ == "__main__":
    createVTUFile_1D()
    createVTUFile_2D()
    createVTUFile_3D()
