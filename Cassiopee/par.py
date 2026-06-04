import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Converter.Filter2 as Filter2
import Connector.PyTree as X
import Connector.Mpi as Xmpi
import Generator.PyTree as G
import Post.PyTree as P
import Post.Mpi as Pmpi

import sys
import numpy


N = 20
meshType = "NGON"

def density(x, y, z):
    return x #numpy.cos(x)*numpy.sin(y)

if Cmpi.rank == 0:
    zones = []
    for rank in range(4):
        xblhc = float(rank/2>=1)*(N-1)
        yblhc = float(rank%2)*(N-1)
        if meshType == "NGON":
            a = G.cartNGon((xblhc, yblhc, 0.), (1., 1., 1.), (N, N, N), api=3)
        else:
            a = G.cart((xblhc, yblhc, 0.), (1., 1., 1.), (N, N, N))
        if Cmpi.size == 4: Cmpi._setProc(a, rank)
        zones.append(a)
    t = C.newPyTree(["Base", *zones])
    C._initVars(
        t,
        'centers:Density',
        density,
        ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'],
        isVectorized=True
    )
    C._initVars(
        t,
        'centers:VelocityX',
        numpy.cos,
        ['centers:CoordinateX'],
        isVectorized=True
    )
    C._initVars(
        t,
        'centers:VelocityY',
        numpy.sin,
        ['centers:CoordinateY'],
        isVectorized=True
    )
    C._initVars(
        t,
        'centers:VelocityZ',
        numpy.cos,
        ['centers:CoordinateZ'],
        isVectorized=True
    )
    C.convertPyTree2File(t, f"parCases/par{meshType}.cgns")

if Cmpi.size == 1:
    t = C.convertFile2PyTree(f"parCases/par{meshType}.cgns")
    t = X.connectMatch(t)
    P._computeGrad2(t, var='centers:Density')
    P._computeDiv2(t, var='centers:Velocity')
    C.convertPyTree2File(t, f"parCases/par{meshType}2_REF.cgns")
else:
    if meshType == "NGON":
        #t = Filter2.loadAsChunks(f"parCases/par{meshType}.cgns")
        t, _ = Filter2.loadAndSplit(f"parCases/par{meshType}.cgns")
        #Xmpi._connectMatchNGon(t)
        Pmpi._computeDiv2(t, var='centers:Velocity', rmVar=False)
        Pmpi._computeGrad2(t, var='centers:Density')
        Cmpi.convertPyTree2File(t, f"par{meshType}2.cgns")
    else:  # STRUCT
        h = Filter.Handle(f"parCases/par{meshType}.cgns")
        t = h.loadAndDistribute()
        t = Xmpi.connectMatch(t)
        Pmpi._computeDiv2(t, var='centers:Velocity', rmVar=False)
        Pmpi._computeGrad2(t, var='centers:Density')
        Cmpi.convertPyTree2File(t, f"parCases/par{meshType}2.cgns")
