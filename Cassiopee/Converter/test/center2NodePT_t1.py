# - center2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x,y): return 2*x+y

def H(x,y):
    if (x+y > 5): return 0
    else: return 1

# -- STRUCT
# center2Node: cree une nouvelle zone
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
C._initVars(a,'centers:F',1.)
C._initVars(a, 'Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'cellN', H, ['CoordinateX', 'CoordinateY'])
b = C.center2Node(a); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 1)

# center2Node: modifie une variable
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX', 'CoordinateY'])
a = C.rmVars(a,['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 2)

# -- BE
# center2Node: cree une nouvelle zone
ni = 30; nj = 40; nk = 2
a = G.cartTetra((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
C._initVars(a, 'centers:F', 1.)
C._initVars(a, 'Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'cellN', H, ['CoordinateX', 'CoordinateY'])
b = C.center2Node(a,'centers:F'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 3)

# center2Node: modifie une variable
ni = 30; nj = 40; nk = 2
a = G.cartTetra((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX', 'CoordinateY'])
a = C.rmVars(a,['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 4)

# NGON, api 1
# center2Node: cree une nouvelle zone
ni = 30; nj = 40; nk = 2
a = G.cartNGon((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk), api=1)
C._initVars(a, 'centers:F', 1.)
C._initVars(a, 'Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'cellN', H, ['CoordinateX', 'CoordinateY'])
b = C.center2Node(a,'centers:F'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 5)

# center2Node: modifie une variable
ni = 30; nj = 40; nk = 2
a = G.cartNGon((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk), api=1)
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX', 'CoordinateY'])
a = C.rmVars(a, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 6)

# NGON, api 3
# center2Node: modifie une variable
ni = 30; nj = 40; nk = 2
a = G.cartNGon((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk), api=3)
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX', 'CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX', 'CoordinateY'])
a = C.rmVars(a, ['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 7)

# -- ME
# center2Node: cree une nouvelle zone
ni = 30; nj = 40; nk = 2
a = G.cartTetra((0,(nj-1)/4.,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
b = G.cartPyra(((ni-1)/3.,(nj-1)/4.,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
c = G.cartPenta((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
d = G.cartHexa(((ni-1)/3.,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
a = C.mergeConnectivity([a, b, c, d], None)
C._initVars(a, 'centers:F', 1.)
C._initVars(a, 'Density', F, ['CoordinateX', 'CoordinateY'])
#C._initVars(a, 'cellN', H, ['CoordinateX', 'CoordinateY'])
b = C.center2Node(a, 'centers:F'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 8)

# center2Node: modifie une variable
ni = 30; nj = 40; nk = 2
a = G.cartTetra((0,(nj-1)/4.,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
b = G.cartPyra(((ni-1)/3.,(nj-1)/4.,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
c = G.cartPenta((0,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
d = G.cartHexa(((ni-1)/3.,0,0), (10./(ni-1), 10./(nj-1),1), (ni,nj,nk))
a = C.mergeConnectivity([a, b, c, d], None)
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX', 'CoordinateY'])
#C._initVars(a, 'centers:cellN', H, ['CoordinateX', 'CoordinateY'])
a = C.rmVars(a,['centers:CoordinateX', 'centers:CoordinateY', 'centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1', 3, b])
test.testT(t, 9)
