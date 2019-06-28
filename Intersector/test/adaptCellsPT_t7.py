import Geom.PyTree as D
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Converter.Internal as Internal
import time
import KCore.test as test

a = D.triangle((0,0,0), (0,1,0), (1, 0, 0))
a = C.convertArray2Hexa(a)
a = G.quad2Pyra(a)
a = C.convertArray2NGon(a)
a = G.close(a)

#C.convertPyTree2File(a,'z_nonreg7.plt')

b = G.cart((0.4,0.4,-0.2), (0.01,0.01,-0.01),(5,5,5)) 
b = C.convertArray2NGon(b)
b = G.close(b)

#C.convertPyTree2File(b,'z_nonreg7_source.cgns')

a = C.fillEmptyBCWith(a, 'wall', 'BCWall')

m1 = XOR.adaptCells(a,b, sensor_type=0)
m1 = XOR.closeOctalCells(m1)
#C.convertPyTree2File(m1, 'out7.cgns')
test.testT(m1,1)

m2 = XOR.adaptCells(a,b, sensor_type=2)
m2 = XOR.closeOctalCells(m2)
#C.convertPyTree2File(m2, 'out71.cgns')
test.testT(m2,2)