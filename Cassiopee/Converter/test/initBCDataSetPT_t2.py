# - initBCDataSet (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import KCore.test as test

N = 10

# -- initBCDataSet without having to create a BCDataSet node first.
# A family is also defined.
a = G.cart((0,0,0), (1,1,1), (N,N,N))
C._addBC2Zone(a, 'sym1', 'BCSymmetryPlane', 'imin')
C._addBC2Zone(a, 'sym2', 'BCSymmetryPlane', 'imax')
C._addBC2Zone(a, 'fuselage', 'FamilySpecified:AIRCRAFT', 'jmin')
C._addBC2Zone(a, 'wing', 'FamilySpecified:AIRCRAFT', 'jmax')
C._addBC2Zone(a, 'inlet', 'BCInflow', 'kmin')
C._addBC2Zone(a, 'outlet', 'BCOutflow', 'kmax')
t = C.newPyTree(['Base', a])
C._addFamily2Base(t[2][1], 'AIRCRAFT', bndType='BCWallViscous')

# Add a volume field
C._initVars(t, '{MomentumX}=0.')
# Init all BC nodes
C._initBCDataSet(t, '{MomentumX}=2.')
C._initBCDataSet(t, '{var}=4.*{MomentumX}')
# Add another volume field - does not zero out BCDatasets for that field
C._initVars(t, '{var}=0.')

# Init a single BC node
C._initBCDataSet(t, 'MomentumY', 3., bndName='wall1')
# Edit BC nodes by type
C._initBCDataSet(t, 'MomentumX', 1., bndType='BCWallViscous')
test.testT(t, 1)
