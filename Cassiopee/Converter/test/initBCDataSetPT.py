# - initBCDataSet (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G

N = 10; nfaces = (N-1)*(N-1)

a = G.cart((0,0,0), (1,1,1), (N,N,N))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=nfaces*[1.], parent=d)

# Init all BC nodes
C._initBCDataSet(a, 'MomentumX=2.*numpy.min({Density},0)')
# Init a single BC node only
C._initBCDataSet(a, 'MomentumY=3.', bndName='wall')
C.convertPyTree2File(a, 'out.cgns')
