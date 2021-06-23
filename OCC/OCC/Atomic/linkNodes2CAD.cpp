/*    
    Copyright 2013-2021 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "occ.h"
#include <TopoDS_Face.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRep_Tool.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <StdFail_NotDone.hxx>

# include <string>
# include <sstream> 

# include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/BbTree.h"

#include <iostream>
//#include <memory>

#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

#include <chrono>


using namespace std;
using namespace K_FLD;
using namespace NUGA;



E_Int check_is_NGON(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{

  E_Int ni, nj, nk;
  
  E_Int res = K_ARRAY::getFromArray(arr, varString, f1, ni, nj, nk,
                                    cn1, eltType);

     
  bool err = (res !=2);
  err &= (strcmp(eltType, "NGON") != 0);

  if (err)
  {
    std::stringstream o;
    o << "input error : " << eltType << " is an invalid array, must be a NGON array." ;
    PyErr_SetString(PyExc_TypeError, o.str().c_str());
    return 1;
  }

  // Check coordinates.
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    PyErr_SetString(PyExc_TypeError, "input error : can't find coordinates in array.");//fixme  conformUnstr
    return 1;
  }
  
  return 0;
}

std::vector<int> getFacesVector (PyObject* arrf){
  
  std::vector<E_Int> Flist;
  {
    FldArrayI* inds=NULL;
    E_Int res=0;
    if (arrf != Py_None)
      res = K_NUMPY::getFromNumpyArray(arrf, inds, true);

    std::unique_ptr<FldArrayI> pL(inds); // to avoid to call explicit delete at several places in the code.
  
    if ((res == 1) && (inds != NULL)  && (inds->getSize() != 0))
    {
      size_t nb_faces = (size_t)inds->getSize();
      Flist.resize(nb_faces);
      for (size_t i = 0; i < nb_faces; ++i) 
      {
        Flist[i]=(*inds)[i]-1;
        //std::cout << Flist[i] << std::endl;
      }
    }
  }
  return Flist;
}


// ============================================================================

//    linkNodes2CAD

// ============================================================================
PyObject* K_OCC::linkNodes2CAD(PyObject* self, PyObject* args)
{
  PyObject *arr, *arrf;
  PyObject* hook;
  PyObject* dhx; PyObject* dhy; PyObject* dhz; PyObject* ncad; 

  if (!PyArg_ParseTuple(args, "OOOOOOO", &arr, &arrf, &hook, &dhx, &dhy, &dhz, &ncad)) return NULL;

  void** packet = NULL;
  #if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    packet = (void**) PyCObject_AsVoidPtr(hook);
  #else
    packet = (void**) PyCapsule_GetPointer(hook, NULL);
  #endif

  E_Float tol = 1e-13;
  
// Get first component of each list 
  arr = PyList_GetItem(arr, 0);
  dhx = PyList_GetItem(dhx, 0);
  dhy = PyList_GetItem(dhy, 0);
  dhz = PyList_GetItem(dhz, 0);
  ncad = PyList_GetItem(ncad, 0);

// array (mesh)
  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString;  char* eltType;
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;

  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);


// faces vector
  std::vector<E_Int>  Flist = getFacesVector(arrf);
    
    
  E_Int erreur = 0;
  
// array with coordinates
  FldArrayF* fi; E_Int ni, nj, nk;
  FldArrayI* c;
  E_Int ret = K_ARRAY::getFromArray2(arr, varString, fi, ni, nj, nk, c, eltType);
  if (ret != 1 && ret != 2) erreur = 1;

// array hx, hy, hz, ncad
  FldArrayF* fhx; FldArrayF* fhy; FldArrayF* fhz;
  FldArrayF* CADfi;
  ret = K_ARRAY::getFromArray2(ncad, CADfi); if (ret != 1 && ret != 2) erreur = 1;
  ret = K_ARRAY::getFromArray2(dhx, fhx); if (ret != 1 && ret != 2) erreur = 1;
  ret = K_ARRAY::getFromArray2(dhy, fhy); if (ret != 1 && ret != 2) erreur = 1;
  ret = K_ARRAY::getFromArray2(dhz, fhz); if (ret != 1 && ret != 2) erreur = 1;
  
  if (erreur == 1 )
  {
    PyErr_SetString(PyExc_TypeError, "linkNodes2CAD: invalid arrays input.");
    return NULL;
  }
  
  E_Float* px = fi->begin(1); // fix
  E_Float* py = fi->begin(2);
  E_Float* pz = fi->begin(3);
  E_Int npts = fi->getSize();
  
  E_Float* pncad = CADfi->begin(1);
  E_Float* phx = fhx->begin(1);
  E_Float* phy = fhy->begin(1);
  E_Float* phz = fhz->begin(1); 


// array (nPoints of faces)
  std::set<E_Int> fPoints;

// Unique points of the faces: 
  for (E_Int f = 0; f < Flist.size(); ++f)
  {
    E_Int PGi = Flist[f];
    E_Int nb_nodes = ngi.PGs.stride(PGi);
    E_Int* p_nodes = ngi.PGs.get_facets_ptr(PGi);

    for (E_Int n = 0; n < nb_nodes; ++n) fPoints.insert(p_nodes[n]-1);
  }
      

// liste des no des faces sur lesquelles on projete
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  
  E_Int nfaces = surfaces.Extent(); 
  std::vector<E_Int> faces(nfaces);
  for (E_Int i = 0; i < nfaces; i++) faces[i] = i+1;

  TopExp_Explorer expl;

// BoundingBoxesTree
  Vector_t<K_SEARCH::BBox3D*> boxes(nfaces);
  for (E_Int j=0; j < nfaces; j++)
  {
	Bnd_Box box;
    BRepBndLib::Add(surfaces(faces[j]),box);
    gp_Pnt PointMin = box.CornerMin();
    gp_Pnt PointMax = box.CornerMax();
    
    E_Float minB[] = { PointMin.X(), PointMin.Y(), PointMin.Z()};
    E_Float maxB[] = { PointMax.X(), PointMax.Y(), PointMax.Z()};
    
    K_SEARCH::BBox3D* Kbox = new K_SEARCH::BBox3D(minB, maxB); boxes[j] = Kbox; 
    boxes[j]-> enlarge(0.0001);
  }
  K_SEARCH::BbTree3D tree(boxes, EPSILON);
   

// Projection des points sur les faces
  E_Float* ptx = new E_Float [npts];
  E_Float* pty = new E_Float [npts];
  E_Float* ptz = new E_Float [npts];
  
  E_Float* dist = new E_Float [npts];

  for (E_Int i = 0; i < npts; i++) dist[i] = K_CONST::E_MAX_FLOAT;
  
  gp_Pnt Point;
  E_Float dx,dy,dz,d;
  
  E_Int nb_pro = 0; 
  E_Int nb_point_pro = 0;
   
  E_Int nb_point_i1 = 0; 
  E_Int nb_point_i2 = 0;
  E_Int nb_point_i3 = 0;
  
  E_Int nb_pro_i1 = 0;
  E_Int nb_pro_i2 = 0;
  
  chrono::steady_clock sc;
  auto start = sc.now();
  

  for (auto index : fPoints)
  {
    
    E_Int k_size = 1;
    
    bool bbtree = false; 
    
    
    if (phx[index] > 1e9){
        
        nb_point_pro += 1; 
        
        
        // BBtree : andidate faces for point i
        Vector_t<E_Int> captured_boxes;
        if (pncad[index] == -1) 
        {
          gp_Pnt PointA;
          PointA.SetCoord(px[index], py[index], pz[index]);
          E_Float pA[] = { PointA.X(), PointA.Y(), PointA.Z()};

          tree.getOverlappingBoxes(pA, pA, captured_boxes);

          k_size = (captured_boxes.size() != 0) ? captured_boxes.size() : nfaces;
          //bbtree = true;
        }
        
        if(pncad[index] == -1)
        {
          nb_point_i1 += 1; 
          nb_pro_i1 += k_size;
        }
            
        if(pncad[index] != -1){
            nb_point_i2 += 1; 
            nb_pro_i2 += k_size;
        }
        
        nb_pro += k_size;
     
        // Projection on candidate faces
        for (E_Int k = 0; k < k_size; k++)
        {
          E_Int j;
          
          if (pncad[index] != -1 && k_size == 1) 
            j = pncad[index]-1;
          else
            j = (captured_boxes.size() != 0) ? captured_boxes[k] : k;

          if (j < 0 || j >= nfaces)
          {
            std::cout << "BUG TO FIX" << std::endl;
            continue;
          }

          const TopoDS_Face& F = TopoDS::Face(surfaces(faces[j]));
          Handle(Geom_Surface) face = BRep_Tool::Surface(F);
      
          Point.SetCoord(px[index], py[index], pz[index]);

          try
          {  
            GeomAPI_ProjectPointOnSurf o(Point, face, Extrema_ExtAlgo_Tree);
            gp_Pnt Pj = o.NearestPoint();
            
            ptx[index] = Pj.X(); pty[index] = Pj.Y(); ptz[index] = Pj.Z();
          }
          catch (StdFail_NotDone& e) 
          { 
            printf("FAIL for point %g %g %g\n", px[index],py[index],pz[index]); 
            ptx[index] = K_CONST::E_MAX_FLOAT;
            pty[index] = K_CONST::E_MAX_FLOAT;
            ptz[index] = K_CONST::E_MAX_FLOAT;
          }
          
          dx = ptx[index]-px[index];
          dy = pty[index]-py[index];
          dz = ptz[index]-pz[index];
          d = dx*dx+dy*dy+dz*dz;
        
          // Result of the closest face
          if (d < dist[index])
          { 
            phx[index] = (abs(dx) <= tol) ? 0 : dx; 
            phy[index] = (abs(dy) <= tol) ? 0 : dy;
            phz[index] = (abs(dz) <= tol) ? 0 : dz; 
            dist[index] = d; pncad[index]=faces[j]; 
            //printf("Point %i : hx: %g, hy: %g, hz: %g      CADid:%i\n",index+1,  dx,dy,dz,pncad[index]);
          }
        }
    }
    else nb_point_i3 += 1;
  }
  auto end = sc.now(); 
  auto time_span = static_cast<chrono::duration<double>>(end - start);
  cout<<"Proyections took: "<<time_span.count()<<" seconds ! \n";
  
  printf("Nombre de points dans BC : %i \n", fPoints.size());
  printf("Nombre de projections  : %i \n", nb_pro);
  
  printf("  nb_point_i1: %i , nb_pro_i1: %i \n", nb_point_i1, nb_pro_i1);
  printf("  nb_point_i2: %i , nb_pro_i2: %i \n", nb_point_i2, nb_pro_i2);
  printf("  nb_point_i3: %i \n", nb_point_i3);

  delete [] ptx; delete [] pty; delete [] ptz;
  boxes.clear();

  Py_DECREF(Py_None);
  return Py_None;

}
// ============================================================================

//    updateFcadidFromNcadid

// ============================================================================
PyObject* K_OCC::updateFcadidFromNcadid(PyObject* self, PyObject* args)
{
  PyObject *arr, *arrf;
  PyObject* ncad; PyObject* fcad;
  if (!PyArg_ParseTuple(args, "OOOO", &arr, &arrf, &ncad, &fcad)) return NULL;
  
// Get first component of each input list 
  arr = PyList_GetItem(arr, 0);
  ncad = PyList_GetItem(ncad, 0);

// array (mesh)
  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString;  char* eltType;
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;


  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

// faces vector
  std::vector<E_Int>  Flist = getFacesVector(arrf);

// ncadid array
  E_Int res;
  FldArrayF* CADfi;
  res = K_ARRAY::getFromArray2(ncad, CADfi); if (res != 1 && res != 2) return NULL;
  E_Float* pncad = CADfi->begin(1);
  
// fcadid array
  FldArrayI* Fid;
  res = K_NUMPY::getFromNumpyArray(fcad, Fid, true); if (res != 1 && res != 2) return NULL;
  E_Int* pFid = Fid->begin(1);

// Update fcadid
  for (E_Int f = 0; f < Flist.size(); ++f)
  {
    E_Int PGi = Flist[f];
    E_Int nb_nodes = ngi.PGs.stride(PGi);
    E_Int* p_nodes = ngi.PGs.get_facets_ptr(PGi);
    //printf("i: %i, PGI: %i, nb_nodes: %i,  \n", f,PGi, nb_nodes);
    
    E_Int nCADid_0 = -2;
    bool same_ncadid = true;

    for (E_Int n = 0; n < nb_nodes; ++n)
    {
      E_Int Ni = p_nodes[n]-1;
      E_Int nCADid = pncad[Ni];
      if ( n == 0 ) nCADid_0 = nCADid;
      if ( nCADid != nCADid_0 ) same_ncadid = false;
      //printf("              n: %i,  Ni: %2i,  CADid(NI): %i\n",n, Ni, nCADid);
    }
    
    if ( same_ncadid ) pFid[Flist[f]] = nCADid_0;
    //printf("              Face: %2i, id: %f \n", PGi, pFid[f]);
  }
 
  
  Py_DECREF(Py_None);
  return Py_None;

}

// ============================================================================

//    updateNcadidFromFcadid

// ============================================================================
PyObject* K_OCC::updateNcadidFromFcadid(PyObject* self, PyObject* args)
{
  PyObject *arr, *arrf;
  PyObject* ncad; PyObject* fcad;
  if (!PyArg_ParseTuple(args, "OOOO", &arr, &arrf, &ncad, &fcad)) return NULL;

// Get first component of each input list 
  arr = PyList_GetItem(arr, 0);
  ncad = PyList_GetItem(ncad, 0);

// array (mesh)
  K_FLD::FloatArray* f(0);
  K_FLD::IntArray* cn(0);
  char* varString;  char* eltType;
  E_Int err = check_is_NGON(arr, f, cn, varString, eltType);
  if (err) return NULL;


  K_FLD::FloatArray & crd = *f;
  K_FLD::IntArray & cnt = *cn;
  
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  ngon_type ngi(cnt);

// faces vector
  std::vector<E_Int>  Flist = getFacesVector(arrf);
  
  E_Int res;
// Array ncadid
  FldArrayF* arr_ncad;
  res = K_ARRAY::getFromArray2(ncad, arr_ncad); if (res != 1 && res != 2) return NULL;
  E_Float* pncad = arr_ncad->begin(1);
  
// Array fcadid
  FldArrayI* arr_fcad;
  res = K_NUMPY::getFromNumpyArray(fcad, arr_fcad, true); if (res != 1 && res != 2) return NULL;
  E_Int* pfcad = arr_fcad->begin(1);
  

// Update ncadid
  std::set<E_Int> fPoints;

  for (E_Int f = 0; f < Flist.size(); ++f)
  {
    E_Int PGi = Flist[f];
    E_Int nb_nodes = ngi.PGs.stride(PGi);
    E_Int* p_nodes = ngi.PGs.get_facets_ptr(PGi);
    
    
    
    //printf("i: %3i, PGI: %i, nb_nodes: %i, fcadid: %i  \n", f,PGi, nb_nodes, pfcad[PGi]);
    
    if (pfcad[PGi] != -1) {
        for (E_Int n = 0; n < nb_nodes; ++n)
        {
          E_Int Ni = p_nodes[n]-1;
          
          if (fPoints.count(Ni) == 0 && pncad[Ni] == -1 ){
            fPoints.insert(Ni);
            E_Int nCADid = pncad[Ni];
            if (pfcad[PGi] != -1) pncad[Ni] = pfcad[PGi];
            //printf("              n: %i,  Ni: %3i,  CADid(NI) before: %2i, CADid(NI) now: %2.0f\n",n, Ni, nCADid, pncad[Ni]);
          } //else printf("              n: %i,  Ni %3i is already updated \n",n, Ni);
        }
    }
  }
 
  Py_DECREF(Py_None);
  return Py_None;

}
