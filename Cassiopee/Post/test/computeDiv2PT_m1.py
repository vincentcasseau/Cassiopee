# - computeDiv2 (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Connector.PyTree as X
import Connector.Mpi as Xmpi
import Generator.PyTree as G
import Post.PyTree as P
import Post.Mpi as Pmpi
import KCore.test as test
import os
import numpy

N = 10
LOCAL = test.getLocal()

def _createTest(filepath, meshType="STRUCT", matchNormalTo="X", plane=None, api=3):

    def _addData(t, plane="XY"):
        def f(x, y, z): return 3.*x + 2.*y + z
        def g(x, y, z): return -1.5*x - 2.5*y -3.5*z
        def h(x, y, z): return 0.25*x + 0.75*y + 1.25*z
        for ax, func in zip(["X", "Y", "Z"], [f, g, h]):
            if plane is None or ax in plane:
                C._initVars(
                    t,
                    'centers:Velocity' + ax.upper(),
                    func,
                    ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'],
                    isVectorized=True
                )
    
    dim = 3 if plane is None else 2
    if Cmpi.master:
        zones = []
        for rank in range(2):
            xblhc = float(rank)*(N-1) if matchNormalTo == "X" else 0.
            yblhc = float(rank)*(N-1) if matchNormalTo == "Y" else 0.
            zblhc = float(rank)*(N-1) if matchNormalTo == "Z" else 0.
            Nx = N if (dim == 3 or "X" in plane) else 1
            Ny = N if (dim == 3 or "Y" in plane) else 1
            Nz = N if (dim == 3 or "Z" in plane) else 1
            xtrhc = 1. if (dim == 3 or "X" in plane) else 0.
            ytrhc = 1. if (dim == 3 or "Y" in plane) else 0.
            ztrhc = 1. if (dim == 3 or "Z" in plane) else 0.
            print(Cmpi.rank, (xblhc, yblhc, zblhc), (xtrhc, ytrhc, ztrhc), (Nx, Ny, Nz))
            if meshType.upper() == "NGON":
                a = G.cartNGon((xblhc, yblhc, zblhc), (xtrhc, ytrhc, ztrhc), (Nx, Ny, Nz), api=api)
            elif meshType.upper() == "STRUCT":
                a = G.cart((xblhc, yblhc, zblhc), (xtrhc, ytrhc, ztrhc), (Nx, Ny, Nz))
            Cmpi._setProc(a, rank)
            zones.append(a)
        t = C.newPyTree(["Base", *zones])
        _addData(t, plane=plane)
        C.convertPyTree2File(t, filepath)
    Cmpi.barrier()

def runTest(filepath, meshType="STRUCT"):
    if meshType.upper() == "NGON":
        h = Filter.Handle(filepath)
        t = h.loadFromProc()
        zones = Internal.getZones(t)
        Xmpi._connectMatchNGon(zones[0])
        Pmpi._computeDiv2(t, var='centers:Velocity', rmVar=False)
    elif meshType.upper() == "STRUCT":
        h = Filter.Handle(filepath)
        t = h.loadFromProc()
        t = Xmpi.connectMatch(t)
        Pmpi._computeDiv2(t, var='centers:Velocity', rmVar=False)
    Cmpi.convertPyTree2File(t, filepath)
    #if Cmpi.master: os.remove(filepath)
    Cmpi.barrier()
    return t


filepath = os.path.join(LOCAL, "out.cgns")

# --- Without BCDataSets --- #
# 2D STRUCT
meshType = "STRUCT"
"""_createTest(filepath, meshType, matchNormalTo="X", plane="XY")
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 1) # WRONG REF!

_createTest(filepath, meshType, matchNormalTo="Y", plane="YZ")
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 2) # WRONG REF!

_createTest(filepath, meshType, matchNormalTo="Z", plane="XZ")
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 3) # WRONG REF!"""

# 3D STRUCT
_createTest(filepath, meshType, matchNormalTo="X")
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 4)

_createTest(filepath, meshType, matchNormalTo="Y")
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 5)

_createTest(filepath, meshType, matchNormalTo="Z")
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 6)

# 2D NGON
meshType = "NGON"
_createTest(filepath, meshType, matchNormalTo="X", plane="XY", api=1)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 7) # CRASH!
exit()

#_createTest(filepath, meshType, matchNormalTo="Y", plane="YZ", api=3)
#t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 8) # CRASH!

#_createTest(filepath, meshType, matchNormalTo="Z", plane="XZ", api=3)
#t = runTest(filepath, meshType)
#if Cmpi.master: test.testT(t, 9) # CRASH!

# 3D NGON
_createTest(filepath, meshType, matchNormalTo="X", api=1)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 10)

_createTest(filepath, meshType, matchNormalTo="Y", api=3)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 11)

_createTest(filepath, meshType, matchNormalTo="Z", api=3)
t = runTest(filepath, meshType)
if Cmpi.master: test.testT(t, 12)
