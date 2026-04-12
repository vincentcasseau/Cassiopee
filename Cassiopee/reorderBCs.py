# reorder sur BC NGON
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import KCore.Vector as Vector

# test case
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
b = G.cartNGon((0,0,0), (1,1,1), (10,10,1), api=3)
C._addBC2Zone(a, 'wall', 'BCWall', subzone=b)
C.convertPyTree2File(a, 'out.cgns')

# subzone a BC
bc = Internal.getNodeFromType(a, 'BC_t')
PL = Internal.getNodeFromName1(bc, 'PointList')[1]
zp = T.subzone(a, PL, type='faces')

# reorder
zp = T.reorder(zp, (1,))

# get pointers
Internal._adaptNFace2PE(a)
PE = Internal.getNodeFromName(a, 'ParentElements')[1]
NFACE = Internal.getNFaceNode(a)
NGON = Internal.getNGonNode(a)
off1 = Internal.getNodeFromName1(NGON, 'ElementStartOffset')[1]
con1 = Internal.getNodeFromName1(NGON, 'ElementConnectivity')[1]
off2 = Internal.getNodeFromName1(NFACE, 'ElementStartOffset')[1]
con2 = Internal.getNodeFromName1(NFACE, 'ElementConnectivity')[1]
xp = Internal.getNodeFromName2(a, 'CoordinateX')[1]
yp = Internal.getNodeFromName2(a, 'CoordinateY')[1]
zp = Internal.getNodeFromName2(a, 'CoordinateZ')[1]

# face 0 et elt 0
face0 = PL[0,0]
elt0 = PE[face0-1,0]
print(face0, elt0)

# barycenter de la face 0
np = off1[face0]-off1[face0-1]+1
PF = (0., 0., 0.)
for i in range(np):
    off = off1[face0-1]
    ind = con1[off+i]-1
    PF = Vector.add(PF, (xp[ind],yp[ind],zp[ind]))
PF = Vector.mul(1./np, PF)
print(PF)

# normale a la face
nF = (0., 0., 0.)
off = off1[face0-1]
for i in range(np-1):
    ind1 = con1[off+i]-1
    ind2 = con1[off+i+1]-1
    P1 = (xp[ind1],yp[ind1],zp[ind1])
    P2 = (xp[ind2],yp[ind2],zp[ind2])
    v = Vector.cross(Vector.sub(P1,PF), Vector.sub(P2,PF))
    nF = Vector.add(nF, v)    
nF = Vector.mul(1./np, nF)
print(nF)

# barycenter de l'element 0
PE = (0., 0., 0.)
cur = 0
nf = off2[elt0]-off2[elt0-1]+1
off = off2[elt0-1]
for j in range(nf):
    facei = con2[off+j]
    np = off1[facei]-off1[facei-1]+1
    offf = off1[facei-1]
    for i in range(np):
        ind = con1[offf+i]-1
        PE = Vector.add(PE, (xp[ind],yp[ind],zp[ind]))
        cur += 1
PE = Vector.mul(1./cur, PE)
print(PE)

# produit scalaire
prod = Vector.dot(nF, Vector.sub(PE, PF))
if prod < 0: T._reorder(zp, (-1,))

