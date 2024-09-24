# - getEmptyBCPT (pyTree) -
# 2D NGon
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test
import testers as Q

# Sur une zone
a = G.cartNGon((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 1))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
#a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
#a = C.addBC2Zone(a, 'wall2', 'BCWall', 'imax')
wins = C.getEmptyBC(a, dim=2)
#test.testO(wins, 1)
a = C.addBC2Zone(a, 'wall1', 'BCWall', faceList=wins)
print(wins)
C.convertPyTree2File(a, 'gg.cgns')
wins = C.getEmptyBC(a, dim=2)
print(wins)
#Q.viewer(a)
#test.testO(wins, 2)

# Sur un arbre
b = G.cartNGon((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
#b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
wins = C.getEmptyBC(t, dim=2)
#test.testO(wins, 2)

# Sur une liste de zones
wins = C.getEmptyBC([a,b], dim=2)
#test.testO(wins, 3)

# Sur un arbre
b = G.cartNGon((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
#b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
t = C.fillEmptyBCWith(t, 'ov', 'BCOverlap', dim=2)
wins = C.getEmptyBC(t, dim=2)
#test.testO(wins, 4)
