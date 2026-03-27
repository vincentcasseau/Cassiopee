# - identifyElements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

# 3D ME: pyra - penta - hexa
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c], None)
#C.convertPyTree2File(a, 'out.cgns')
#C.fillEmptyBCWith(a, 'wall', 'BCWall')
wins = C.getEmptyBC(a)
print(wins)
#n_bcs = C.getBCs(a)
#print(n_bcs)

"""
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 7)"""

