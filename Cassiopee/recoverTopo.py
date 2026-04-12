# - convertArray2NGon(pyTree): topologic -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import numpy as np

# Simple box, Elements, with BCs
N = 10

def _addBCsBySubzone(t, nxblocks=1, is3D=True, quadFaces=True):
    Nz = N if is3D else 1
    Lz = 1. if is3D else 0.
    cartFunc = G.cartHexa if quadFaces else G.cartTetra
    subzone = cartFunc((0., 0., 0.), (0., 1., Lz), (1, N, Nz))
    C._addBC2Zone(t, 'inlet', 'BCInflow', subzone=subzone)
    if is3D: subzone = cartFunc((0., 0., 0.), (1., 1., 0.), (nxblocks*(N - 1) + 1, N, 1))
    else: subzone = cartFunc((0., 0., 0.), (1., 0., 0.), (nxblocks*(N - 1) + 1, 1, 1))
    C._addBC2Zone(t, 'wall', 'BCWallInviscid', subzone=subzone)
    if is3D: subzone = cartFunc((0., 0., Nz - 1.), (1., 1., 0.), (nxblocks*(N - 1) + 1, N, 1))
    else: subzone = cartFunc((0., N - 1., 0.), (1., 0., 0.), (nxblocks*(N - 1) + 1, 1, 1))
    C._addBC2Zone(t, 'farfield', 'BCFarfield', subzone=subzone)
    subzone = cartFunc(((N - 1.) * nxblocks, 0., 0.), (0., 1., Lz), (1, N, Nz))
    C._addBC2Zone(t, 'outlet', 'BCOutflow', subzone=subzone)
    return None

def _addData(t):
    C._initVars(t, "Density=1.05")
    C._initVars(t, "centers:Pressure=1.")

def _addBCDataAtNodes(t):
    b = Internal.getBCNodesFromType(t, bndType='BCInflow')[0]
    er = Internal.getNodeFromName(b, Internal.__ELEMENTRANGE__)
    er = Internal.getValue(er)[0]
    npts = N #er[1] - er[0] + 1
    d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                              gridLocation='Vertex', parent=b)
    d = Internal.newBCData('BCDirichlet', parent=d)
    """C._initBCDataSet(t, 'Density', 1.1, bndType='BCInflow')"""
    Internal.newDataArray('Pressure', value=npts*[1.1], parent=d)
    """Internal.newDataArray('MomentumX', value=npts*[2.], parent=d)
    Internal.newDataArray('MomentumY', value=npts*[0.1], parent=d)
    Internal.newDataArray('MomentumZ', value=npts*[0.1], parent=d)
    Internal.newDataArray('EnergyStagnationDensity', value=npts*[2.], parent=d)"""

def _addBCDataAtFaceCenters(t):
    b = Internal.getBCNodesFromType(t, bndType='BCWall*')[0]
    er = Internal.getNodeFromName(b, Internal.__ELEMENTRANGE__)
    er = Internal.getValue(er)[0]
    nfaces = er[1] - er[0] + 1
    d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                              gridLocation='FaceCenter', parent=b)
    d = Internal.newBCData('BCNeumann', parent=d)
    """C._initBCDataSet(t, 'Density', 1.05, bndType='BCWall*')
    C._initBCDataSet(t, 'MomentumX', 1., bndType='BCWall*')
    C._initBCDataSet(t, 'MomentumY', 0., bndType='BCWall*')
    C._initBCDataSet(t, 'MomentumZ', 0., bndType='BCWall*')
    C._initBCDataSet(t, 'EnergyStagnationDensity', 1., bndType='BCWall*')"""
    Internal.newDataArray('Density', value=nfaces*[1.05], parent=d)
    Internal.newDataArray('MomentumX', value=nfaces*[1.], parent=d)
    Internal.newDataArray('MomentumY', value=nfaces*[0.], parent=d)
    Internal.newDataArray('MomentumZ', value=nfaces*[0.], parent=d)
    Internal.newDataArray('EnergyStagnationDensity', value=nfaces*[1.], parent=d)

def _convertStruct2NGon(method="geometric", addDataSets=False, api=1):
    C.clearAllNames()
    t = G.cart((0., 0., 0.), (1., 1., 1.), (N, N, N))
    C._addBC2Zone(t, 'inlet', 'BCInflow', wrange='imin')
    C._addBC2Zone(t, 'wall', 'BCWallInviscid', wrange='imax')
    C._addBC2Zone(t, 'farfield', 'BCFarfield', wrange='jmin')
    C._addBC2Zone(t, 'outlet', 'BCOutflow', wrange='jmax')
    _addData(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    C.convertPyTree2File(t, 'recoCases/out_3dSTRUCT.cgns')
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=api)
    if method == "geometric":
        C.convertPyTree2File(tng, 'recoCases/out_3dSTRUCT_ng_REF.cgns')
    else:
        C.convertPyTree2File(tng, 'recoCases/out_3dSTRUCT_ng.cgns')

def _convertQuad2NGon(method="geometric", addDataSets=False, api=1):
    C.clearAllNames()
    t = G.cartHexa((0., 0., 0.), (1., 1., 0.), (N, N, 1))
    _addBCsBySubzone(t, nxblocks=1, quadFaces=True, is3D=False)
    _addData(t)
    # _addBCDataAtNodes(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    C.convertPyTree2File(t, 'recoCases/out_2dBE.cgns')
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=api)
    if method == "geometric":
        C.convertPyTree2File(tng, 'recoCases/out_2dBE_ng_REF.cgns')
    else:
        C.convertPyTree2File(tng, 'recoCases/out_2dBE_ng.cgns')

def _convertHexa2NGon(method="geometric", addDataSets=False, api=1):
    C.clearAllNames()
    t = G.cartHexa((0., 0., 0.), (1., 1., 1.), (N, N, N))
    _addBCsBySubzone(t, nxblocks=1, quadFaces=True)
    _addData(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    C.convertPyTree2File(t, 'recoCases/out_3dBE.cgns')
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=api)
    if method == "geometric":
        C.convertPyTree2File(tng, 'recoCases/out_3dBE_ng_REF.cgns')
    else:
        C.convertPyTree2File(tng, 'recoCases/out_3dBE_ng.cgns')

def _convertME2NGon(method="geometric", addDataSets=False):
    C.clearAllNames()
    a = G.cartPyra((0., 0., 0.), (1., 1., 1.), (N, N, N))
    b = G.cartHexa((N-1., 0., 0.), (1., 1., 1.), (N, N, N))
    t = C.mergeConnectivity(a, b, boundary=0)
    _addBCsBySubzone(t, nxblocks=2, quadFaces=True)
    _addData(t)
    if addDataSets: _addBCDataAtFaceCenters(t)
    C.convertPyTree2File(t, 'recoCases/out_3dME.cgns')
    tng = C.convertArray2NGon(t, recoverBC=True, method=method, api=3)
    if method == "geometric":
        C.convertPyTree2File(tng, 'recoCases/out_3dME_ng_REF.cgns')
    else:
        C.convertPyTree2File(tng, 'recoCases/out_3dME_ng.cgns')

# First create the ref using the geometric method
# Then check against the ref using the topologic method
# 3D STRUCT
#_convertStruct2NGon(method="geometric", addDataSets=False, api=3)
#_convertStruct2NGon(method="topologic", addDataSets=False, api=3)

# 2D BE
#_convertQuad2NGon(method="geometric", addDataSets=True, api=3)
#_convertQuad2NGon(method="topologic", addDataSets=False, api=3)

# 3D BE
#_convertHexa2NGon(method="geometric", addDataSets=True, api=1)
_convertHexa2NGon(method="topologic", addDataSets=True, api=1)

# 3D ME
#_convertME2NGon(method="geometric", addDataSets=True)
#_convertME2NGon(method="topologic", addDataSets=False)
