# - quad2Pyra (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test
test.TOLERANCE=1e-8

a = D.sphere6((0,0,0), 1, N=10)
a = C.convertArray2Hexa(a)
a = T.join(a)
a = G.close(a)
b = G.quad2Pyra(a, hratio=1.)
C.convertPyTree2File(b, 'out2.cgns')
test.testT(b,1)


# a = G.cartHexa((0,0,0), (1,1,0), (5, 5, 1))
# b = G.quad2Pyra(a, hratio=1.)
# G._getVolumeMap(b)
# C.convertPyTree2File(b, 'out.cgns')
# test.testT(b,2)
