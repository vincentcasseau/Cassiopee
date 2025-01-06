# - check if mesh is cartesian -
import Transform.PyTree as T
import KCore.test as test
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Generator.IBM as G_IBM
import random

a = G.cart((0.,0.,0.), (0.1,0.1,1), (40,40,40))
t = C.newPyTree(['CARTESIAN', a])
C._addState(t, 'EquationDimension', 3)

t = T.splitNParts(t, 10)

##Test 1 - Cartesian
cartesian=G_IBM.checkCartesian(t)
isCartesian=0
if cartesian:isCartesian=1
Internal._createUniqueChild(t, 'TMP_Node', 'UserDefinedData_t')
Internal._createUniqueChild(Internal.getNodeFromName1(t, 'TMP_Node'), 'Cartesian', 'DataArray_t', value=isCartesian)
test.testT(t,1)
tsave = Internal.copyTree(t)
print('Is Cartesian::',cartesian,isCartesian)
#C.convertPyTree2File(t,'t_test1.cgns')

## Test 2 - Z direction is non homogenous --> non Cartesian mesh
Internal._rmNode(t, Internal.getNodeFromName(t, 'TMP_Node'))
coord = Internal.getNodeFromName(Internal.getZones(t)[0],'CoordinateZ')[1]
for i in range(6): coord[:,:,i]=coord[:,:,i]*i/10+coord[:,:,i]
cartesian=G_IBM.checkCartesian(t)

isCartesian=0
if cartesian:isCartesian=1
Internal._createUniqueChild(t, 'TMP_Node', 'UserDefinedData_t')
Internal._createUniqueChild(Internal.getNodeFromName1(t, 'TMP_Node'), 'Cartesian', 'DataArray_t', value=isCartesian)
test.testT(t,2)
print('Is Cartesian::',cartesian,isCartesian)
#C.convertPyTree2File(t,'t_test2.cgns')

## Test 3 - Y direction is non homogenous --> non Cartesian mesh
t = Internal.copyTree(tsave)
Internal._rmNode(t, Internal.getNodeFromName(t, 'TMP_Node'))
coord = Internal.getNodeFromName(Internal.getZones(t)[0],'CoordinateY')[1]
for i in range(6): coord[:,i,:]=coord[:,i,:]*i/10+coord[:,i,:]
cartesian=G_IBM.checkCartesian(t)

isCartesian=0
if cartesian:isCartesian=1
Internal._createUniqueChild(t, 'TMP_Node', 'UserDefinedData_t')
Internal._createUniqueChild(Internal.getNodeFromName1(t, 'TMP_Node'), 'Cartesian', 'DataArray_t', value=isCartesian)
test.testT(t,3)
print('Is Cartesian::',cartesian,isCartesian)
#C.convertPyTree2File(t,'t_test3.cgns')

## Test 4 - Z direction is non homogenous --> non Cartesian mesh
t = Internal.copyTree(tsave)
Internal._rmNode(t, Internal.getNodeFromName(t, 'TMP_Node'))
coord = Internal.getNodeFromName(Internal.getZones(t)[0],'CoordinateX')[1]
for i in range(6): coord[i,:,:]=coord[i,:,:]*i/10+coord[i,:,:]
cartesian=G_IBM.checkCartesian(t)

isCartesian=0
if cartesian:isCartesian=1
Internal._createUniqueChild(t, 'TMP_Node', 'UserDefinedData_t')
Internal._createUniqueChild(Internal.getNodeFromName1(t, 'TMP_Node'), 'Cartesian', 'DataArray_t', value=isCartesian)
test.testT(t,4)
print('Is Cartesian::',cartesian,isCartesian)
#C.convertPyTree2File(t,'t_test4.cgns')
