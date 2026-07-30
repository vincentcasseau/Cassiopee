/*
    Copyright 2013-2026 ONERA.

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
# include "post.h"

using namespace K_FLD;

//=============================================================================
/* Compute the divergence of a set of vector fields given in cell centers
   The divergence is given on cell centers. */
//=============================================================================
PyObject* K_POST::computeDiv2(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayc;
  PyObject* volc; PyObject* cellNc;
  PyObject* indices; PyObject* fieldX; PyObject* fieldY; PyObject* fieldZ;
  if (!PYPARSETUPLE_(args, OOOO_ OOOO_, &array, &arrayc, &volc, &cellNc,
                      &indices, &fieldX, &fieldY, &fieldZ)) return NULL;

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn,
                                     eltType);

  E_Int dim = 0;
                                     
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: invalid array.");
    return NULL;
  }
  else if (res == 1)
  {
    if (ni > 1) dim += 1;
    if (nj > 1) dim += 1;
    if (nk > 1) dim += 1;
    std::cout << "ni" << ni << std::endl;
    std::cout << "nj" << nj << std::endl;
    std::cout << "nk" << nk << std::endl;
  }
  else if (res == 2 && K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    dim = cn->getDim();
  }
  else
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: only for STRUCT and NGon arrays.");
    return NULL;
  }

  if (dim < 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: not valid for 1D arrays.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;

  // Check arrayc
  char* varStringc; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int nic, njc, nkc;
  E_Int resc = K_ARRAY::getFromArray3(arrayc, varStringc, fc, nic, njc, nkc, cnc,
                                      eltTypec);

  // Extract cellN if given
  E_Float* cellNp = NULL;
  E_Int nelts = fc->getSize();
  if (cellNc != Py_None) K_NUMPY::getFromNumpyArray(cellNc, cellNp, nelts);

  // Number of vector fields whose divergence to compute
  E_Int nfld = fc->getNfld();  // number of scalar fields
  E_Int nfld2 = 1; // TODO nfld/dim;  // number of vector fields
  std::vector<char*> vars;
  K_ARRAY::extractVars(varStringc, vars);
  // if (nfld % dim != 0)
  // {
  //   PyErr_SetString(PyExc_TypeError,
  //                   "computeDiv2: not all components were found for each vector field.");
  //   RELEASESHAREDB(res, array, f, cn);
  //   RELEASESHAREDB(resc, arrayc, fc, cnc);
  //   if (cellNc != Py_None) Py_DECREF(cellNc);
  //   return NULL;
  // }

  // Check xyz-plane based only on the first two of the given vector fields
  std::vector<char*> varStrings;
  E_Int ixyz = -1; // in 2D, = 0 XY-plane, = 1 XZ-plane, = 2 YZ-plane

  for (E_Int i = 0; i < nfld2; i++)
  {
    // for (E_Int m = 0; m < dim-1; m++)
    // {
    //   std::cout << m << "; " << vars[dim*i+m] << ", " << vars[dim*i+m+1] << ", " << strlen(vars[dim*i+m])-1 << std::endl;
    //   if (K_STRING::cmp(vars[dim*i+m], strlen(vars[dim*i+m])-1, vars[dim*i+m+1]) != 0)
    //   {
    //     RELEASESHAREDB(res, array, f, cn);
    //     RELEASESHAREDB(resc, arrayc, fc, cnc);
    //     PyErr_SetString(PyExc_TypeError,
    //                     "computeDiv2: invalid names for vector component fields.");
    //     return NULL;
    //   }
    // }

    if (dim == 3)
    {
      char* sv0 = vars[dim*i]; char s0 = sv0[strlen(sv0)-1];
      char* sv1 = vars[dim*i+1]; char s1 = sv1[strlen(sv1)-1];
      char* sv2 = vars[dim*i+2]; char s2 = sv2[strlen(sv2)-1];
      if (s0 != 'X' || s1 != 'Y' || s2 != 'Z')
      {
        if (cellNc != Py_None) Py_DECREF(cellNc);
        RELEASESHAREDB(res, array, f, cn);
        RELEASESHAREDB(resc, arrayc, fc, cnc);
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv2: error with the order of given scalar fields.");
        return NULL;
      }
    }
    else // dim == 2
    {
      // char* sv0 = vars[dim*i]; char s0 = sv0[strlen(sv0)-1];
      // char* sv1 = vars[dim*i+1]; char s1 = sv1[strlen(sv1)-1];
      // if (ixyz == -1)
      // {
      //   if (s0 == 'X' && s1 == 'Y') ixyz = 0;
      //   else if (s0 == 'X' && s1 == 'Z') ixyz = 1;
      //   else if (s0 == 'Y' && s1 == 'Z') ixyz = 2;
      //   else 
      //   {
      //     PyErr_SetString(PyExc_TypeError,
      //                     "computeDiv2: error with the order of given scalar fields.");
      //     if (cellNc != Py_None) Py_DECREF(cellNc);
      //     RELEASESHAREDB(res, array, f, cn);
      //     RELEASESHAREDB(resc, arrayc, fc, cnc);
      //     return NULL;
      //   }
      // }
      // else
      // {
      //   if (
      //     (s0 == 'X' && s1 == 'Y' && ixyz != 0)
      //     || (s0 == 'X' && s1 == 'Z' && ixyz != 1)
      //     || (s0 == 'Y' && s1 == 'Z' && ixyz != 2)
      //   )
      //   {
      //     PyErr_SetString(PyExc_TypeError,
      //       "computeDiv2: inconsistant scalar field components.");
      //     if (cellNc != Py_None) Py_DECREF(cellNc);
      //     RELEASESHAREDB(res, array, f, cn);
      //     RELEASESHAREDB(resc, arrayc, fc, cnc);
      //     return NULL;
      //   }
      // }
    }

    char* local;
    computeDivVarsString(vars[i*dim], local);
    varStrings.push_back(local);
  }
  
  E_Int size = 0;
  for (E_Int i = 0; i < nfld2; i++) size += strlen(varStrings[i]) + 1;
  char* varStringOut = new char [size];
  char* pt = varStringOut;
  for (E_Int i = 0; i < nfld2; i++)
  {
    char* v = varStrings[i];
    for (size_t j = 0; j < strlen(v); j++) { *pt = v[j]; pt++; }
    *pt = ','; pt++;
    delete [] varStrings[i];
  }
  pt--; *pt = '\0';
  for (size_t i = 0; i < vars.size(); i++) delete [] vars[i];

  PyObject* tpl = NULL;
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  if (res == 1)
  {
    E_Int nicnjc = nic*njc;
    E_Int ninjc = ni*njc;
    E_Int nicnj = nic*nj;
    E_Int nbIntI = ninjc*nkc;
    E_Int nbIntJ = nicnj*nkc;
    E_Int nbIntK = nicnjc*nk;
    E_Int nbIntIJ = nbIntI+nbIntJ;
    E_Int nbIntTot = nbIntIJ+nbIntK;
    
    // Build
    FldArrayF faceField(nbIntTot, nfld); faceField.setAllValuesAtNull();
    FldArrayI voisins(nbIntTot, 2); voisins.setAllValuesAt(-1);
    E_Int* cellG = voisins.begin(1); E_Int* cellD = voisins.begin(2);
    
    
    if (dim == 2)
    {
      tpl = computeDiv2Struct2D(
        ni, nj, nk, nic, njc, nkc, ixyz, varStringOut, cellNp,
        xt, yt, zt, *fc, faceField, cellG, cellD,
        indices, fieldX, fieldY, fieldZ
      );
    }
    else
    {
      tpl = computeDiv2Struct3D(
        ni, nj, nk, nic, njc, nkc, varStringOut, cellNp,
        xt, yt, zt, *fc, faceField, cellG, cellD,
        indices, fieldX, fieldY, fieldZ
      );
    }
  }
  else
  {
    // Compute volume of each element
    FldArrayF vol(nelts);
    E_Float* volp = vol.begin();
    if (volc == Py_None) K_METRIC::compVolNGon(xt, yt, zt, *cn, volp);
    else
    {
      FldArrayF* vols = NULL; 
      K_NUMPY::getFromNumpyArray(volc, vols);
      volp = vols->begin();
      RELEASESHAREDN(volc, vols);
    }

    tpl = computeDiv2NGon(cn, varStringOut, volp, cellNp, 
                          xt, yt, zt, *fc, indices, fieldX, fieldY, fieldZ);

    if (volc != Py_None) Py_DECREF(volc);
  }

  delete [] varStringOut;
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(resc, arrayc, fc, cnc);
  if (cellNc != Py_None) Py_DECREF(cellNc);
  return tpl;
}

PyObject* K_POST::computeDiv2NGon(
  FldArrayI* cn, const char* varStringOut, E_Float* volp, E_Float* cellNp,
  E_Float* xt, E_Float* yt, E_Float* zt, 
  FldArrayF& fc, PyObject* indices,
  PyObject* fieldX, PyObject* fieldY, PyObject* fieldZ
)
{
  E_Int dim = cn->getDim();
  E_Int nfaces = cn->getNFaces();
  E_Int nelts = cn->getNElts();
  E_Int nfld = fc.getNfld();  // number of scalar fields
  E_Int nfld2 = 1; //nfld/3;  // number of vector fields
  E_Int api = fc.getApi();
  std::cout << "nfld = " << nfld << std::endl;
  std::cout << "nfld2 = " << nfld2 << std::endl;

  // Compute FE connectivity
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  // Compute field on element faces
  FldArrayF faceField(nfaces, nfld);
  if (cellNp == NULL)
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = faceField.begin(n);
      E_Float* s = fc.begin(n);
      #pragma omp parallel
      {
        E_Int i1, i2;
        #pragma omp for
        for (E_Int i = 0; i < nfaces; i++)
        {
          i1 = cFE1[i] - 1; i2 = cFE2[i] - 1;
          if (i2 != -1) fp[i] = 0.5*(s[i1] + s[i2]);
          else fp[i] = s[i1];
        }
      }
    }
  }
  else // cellN
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = faceField.begin(n);
      E_Float* s = fc.begin(n);
      #pragma omp parallel
      {
        E_Int i1, i2;
        #pragma omp for
        for (E_Int i = 0; i < nfaces; i++)
        {
          i1 = cFE1[i] - 1; i2 = cFE2[i] - 1;
          if (i2 != -1)
          {
            if (cellNp[i1] == 0) fp[i] = s[i2];
            else fp[i] = 0.5*(s[i1] + s[i2]);
          }
          else fp[i] = s[i1];
        }
      }
    }
  }

  // Replace DataSet
  std::cout << "H" << std::endl;
  FldArrayI* inds = NULL; FldArrayF* bfieldX = NULL;
  FldArrayF* bfieldY = NULL; FldArrayF* bfieldZ = NULL;
  if (indices != Py_None && fieldX != Py_None && fieldY != Py_None && fieldZ != Py_None)
  {
    K_NUMPY::getFromNumpyArray(indices, inds);
    K_NUMPY::getFromNumpyArray(fieldX, bfieldX);
    K_NUMPY::getFromNumpyArray(fieldY, bfieldY);
    K_NUMPY::getFromNumpyArray(fieldZ, bfieldZ);
    std::cout << "H1" << std::endl;

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pind = inds->begin();

    #pragma omp parallel
    {
      E_Int ind;
      for (E_Int n = 0; n < nfld2; n++)
      {
        std::cout << "H2" << std::endl;
        E_Float* bfx = bfieldX->begin(n+1);
        E_Float* bfy = bfieldY->begin(n+1);
        E_Float* bfz = bfieldZ->begin(n+1);
        std::cout << "H3" << std::endl;
        E_Float* fpx = faceField.begin(3*n+1);
        E_Float* fpy = faceField.begin(3*n+2);
        E_Float* fpz = faceField.begin(3*n+3);
        std::cout << "H4" << std::endl;

        #pragma omp for
        for (E_Int i = 0; i < ninterfaces; i++)
        {
          ind = pind[i]-1;
          fpx[ind] = bfx[i]; fpy[ind] = bfy[i]; fpz[ind] = bfz[i];
        }
        std::cout << "H5" << std::endl;
      }
    }

    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(fieldX, bfieldX);
    RELEASESHAREDN(fieldY, bfieldY);
    RELEASESHAREDN(fieldZ, bfieldZ);
  }
  std::cout << "K" << std::endl;
  
  // Build unstructured NGON array from existing connectivity
  PyObject* tpl = K_ARRAY::buildArray3(nfld2, varStringOut, nelts,
                                       *cn, "NGON", true, api, true);
  FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2);
  f2->setAllValuesAtNull();

  FldArrayF surf(nfaces, 4);
  E_Float* sxp = surf.begin(1);
  E_Float* syp = surf.begin(2);
  E_Float* szp = surf.begin(3);
  E_Float* snp = surf.begin(4);
  K_METRIC::compNGonFacesSurf(xt, yt, zt, *cn, sxp, syp, szp, snp, &cFE);

  // divergence
  std::cout << "L" << std::endl;
  E_Float ffx, ffy, ffz;
  E_Int i1, i2;
  for (E_Int n = 0; n < nfld2; n++)
  {
    E_Float* f2p = f2->begin(n+1);
    E_Float* fpx = faceField.begin(3*n+1);
    E_Float* fpy = faceField.begin(3*n+2);
    E_Float* fpz = faceField.begin(3*n+3);
    for (E_Int i = 0; i < nfaces; i++)
    {
      i1 = cFE1[i] - 1; i2 = cFE2[i] - 1;
      ffx = fpx[i]; ffy = fpy[i]; ffz = fpz[i];
      if (i1 != -1) f2p[i1] += ffx*sxp[i] + ffy*syp[i] + ffz*szp[i];
      if (i2 != -1) f2p[i2] -= ffx*sxp[i] + ffy*syp[i] + ffz*szp[i];
    }
  }

  // free mem
  surf.malloc(0); faceField.malloc(0);
  std::cout << "M" << nfld2 << std::endl;

  #pragma omp parallel
  {
    for (E_Int n = 1; n <= nfld2; n++)
    {
      E_Float* f2p = f2->begin(n);
      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        f2p[i] /= K_FUNC::E_max(volp[i], K_CONST::E_MIN_VOL);
      }
    }
  }
  std::cout << "N" << std::endl;

  return tpl;
}

//=============================================================================
PyObject* K_POST::computeDiv2Struct3D(
    E_Int ni, E_Int nj, E_Int nk, E_Int nic, E_Int njc, E_Int nkc,
    const char* varStringOut, E_Float* cellNp,
    E_Float* xt, E_Float* yt, E_Float* zt,
    FldArrayF& fc, FldArrayF& faceField, E_Int* cellG, E_Int* cellD,
    PyObject* indices, PyObject* fieldX, PyObject* fieldY, PyObject* fieldZ
)
{
  E_Int nicnjc = nic*njc;
  E_Int ninjc = ni*njc;
  E_Int nicnj = nic*nj;
  
  E_Int nbIntI = ninjc*nkc;
  E_Int nbIntJ = nicnj*nkc;
  E_Int nbIntK = nicnjc*nk;
  E_Int nbIntIJ = nbIntI + nbIntJ;
  E_Int nbIntTot = nbIntIJ + nbIntK;
  
  E_Int api = fc.getApi();
  E_Int nfld = fc.getNfld(); // number of scalar fields
  E_Int nfld2 = nfld/3; // number of vector fields
  E_Int ncells = nicnjc*nkc;

  if (cellNp == NULL)
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fcn = fc.begin(n);
      E_Float* fintp = faceField.begin(n);
      
      #pragma omp parallel
      {
        E_Int indf, indcellg, indcelld;
        
        // faces en i
        #pragma omp for nowait collapse(3)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 1; i < nic; i++)
        {
          indf = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          indcelld = indcellg + 1;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait collapse(2)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
        {
          E_Int i = 0;
          indf = i + j*ni + k*ninjc;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          i = nic;
          indf = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en j
        #pragma omp for nowait collapse(3)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 1; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          indf = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          indcelld = indcellg + nic;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for nowait collapse(2)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int i = 0; i < nic; i++)
        {
          E_Int j = 0;
          indf = i + j*nic + k*nicnj + nbIntI;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indf] = -1;
          cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          j = njc;
          indf = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en k
        #pragma omp for nowait collapse(3)
        for (E_Int k = 1; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          indf = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          indcelld = indcellg + nicnjc;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en k
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          E_Int k = 0;
          indf = i + j*nic + k*nicnjc + nbIntIJ;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld];

          k = nkc;
          indf = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg];
        }
      }
    }
  }
  else // cellN
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fcn = fc.begin(n);
      E_Float* fintp = faceField.begin(n);
      
      #pragma omp parallel
      {
        E_Int indf, indcellg, indcelld;
        
        // faces en i
        #pragma omp for nowait collapse(3)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 1; i < nic; i++)
        {
          indf = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          indcelld = indcellg + 1;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indf] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indf] = fcn[indcellg];
          else fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait collapse(2)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
        {
          E_Int i = 0;
          indf = i + j*ni + k*ninjc;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          i = nic;
          indf = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en j
        #pragma omp for nowait collapse(3)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int j = 1; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          indf = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          indcelld = indcellg + nic;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indf] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indf] = fcn[indcellg];
          else fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for nowait collapse(2)
        for (E_Int k = 0; k < nkc; k++)
        for (E_Int i = 0; i < nic; i++)
        {
          E_Int j = 0;
          indf = i + j*nic + k*nicnj + nbIntI;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indf] = -1;
          cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          j = njc;
          indf = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en k
        #pragma omp for nowait collapse(3)
        for (E_Int k = 1; k < nkc; k++)
        for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          indf = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          indcelld = indcellg + nicnjc;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indf] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indf] = fcn[indcellg];
          else fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en k
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          E_Int k = 0;
          indf = i + j*nic + k*nicnjc + nbIntIJ;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld];

          k = nkc;
          indf = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg];
        }
      }
    }
  }

  // Replace DataSet
  if (indices != Py_None && fieldX != Py_None && fieldY != Py_None && fieldZ != Py_None)
  {
    FldArrayI* inds = NULL; FldArrayF* bfieldX = NULL;
    FldArrayF* bfieldY = NULL; FldArrayF* bfieldZ = NULL;
    K_NUMPY::getFromNumpyArray(indices, inds);
    K_NUMPY::getFromNumpyArray(fieldX, bfieldX);
    K_NUMPY::getFromNumpyArray(fieldY, bfieldY);
    K_NUMPY::getFromNumpyArray(fieldZ, bfieldZ);

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pindint = inds->begin();

    E_Float* pf[3];
    pf[0] = bfieldX->begin();
    pf[1] = bfieldY->begin();
    pf[2] = bfieldZ->begin();

    #pragma omp parallel
    {
      E_Int indf;
      for (E_Int n = 0; n < nfld; n++)
      {
        E_Float* fintp = faceField.begin(n+1);
        #pragma omp for
        for (E_Int noint = 0; noint < ninterfaces; noint++)
        {
          indf = pindint[noint];
          fintp[indf] = pf[n][noint];
        }
      }
    }
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(fieldX, bfieldX);
    RELEASESHAREDN(fieldY, bfieldY);
    RELEASESHAREDN(fieldZ, bfieldZ);
  }

  // Build empty array
  PyObject* tpl = K_ARRAY::buildArray3(nfld2, varStringOut, nic, njc, nkc, api);
  FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2);
  f2->setAllValuesAtNull();

  FldArrayF surf(nbIntTot,3);
  FldArrayF centerInt(nbIntTot,3);
  E_Float* sxp = surf.begin(1);
  E_Float* syp = surf.begin(2);
  E_Float* szp = surf.begin(3);
  FldArrayF surfnorm(nbIntTot);
  E_Float* snp = surfnorm.begin();
  FldArrayF vol(ncells); E_Float* volp = vol.begin();

  K_METRIC::compMetricStruct(
    ni, nj, nk, nbIntI, nbIntJ, nbIntK,
    xt, yt, zt,
    volp, sxp, syp, szp, snp,
    centerInt.begin(1), centerInt.begin(2), centerInt.begin(3)
  );
  centerInt.malloc(0); surfnorm.malloc(0);

  // divergence
  E_Int indcellg, indcelld;
  E_Float ffx, ffy, ffz;
  for (E_Int n = 0; n < nfld2; n++)
  {
    E_Float* f2p = f2->begin(n+1);
    E_Float* fpx = faceField.begin(3*n+1);
    E_Float* fpy = faceField.begin(3*n+2);
    E_Float* fpz = faceField.begin(3*n+3);
    for (E_Int i = 0; i < nbIntTot; i++)
    {
      indcellg = cellG[i]; indcelld = cellD[i];
      ffx = fpx[i]; ffy = fpy[i]; ffz = fpz[i];
      if (indcellg != -1) f2p[indcellg] += ffx*sxp[i] + ffy*syp[i] + ffz*szp[i];
      if (indcelld != -1) f2p[indcelld] -= ffx*sxp[i] + ffy*syp[i] + ffz*szp[i];
    }
  }

  // free mem
  surf.malloc(0); faceField.malloc(0);

  #pragma omp parallel
  {
    for (E_Int n = 1; n <= nfld2; n++)
    {
      E_Float* f2p = f2->begin(n);
      #pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        f2p[i] /= K_FUNC::E_max(volp[i], K_CONST::E_MIN_VOL);
      }
    }
  }

  RELEASESHAREDS(tpl, f2);
  return tpl;
}
//=============================================================================
PyObject* K_POST::computeDiv2Struct2D(
    E_Int ni, E_Int nj, E_Int nk, E_Int nic, E_Int njc, E_Int nkc,
    E_Int ixyz, const char* varStringOut, E_Float* cellNp, 
    E_Float* xt, E_Float* yt, E_Float* zt,
    FldArrayF& fc, FldArrayF& faceField, E_Int* cellG, E_Int* cellD,
    PyObject* indices, PyObject* fieldX, PyObject* fieldY, PyObject* fieldZ
)
{
  E_Int n1, n2, n1c, n2c;

  if (ixyz == 0)  // XY-plane
  {
    n1 = ni; n2 = nj; n1c = nic; n2c = njc;
  }
  else if (ixyz == 1)  // XZ-plane
  {
    n1 = ni; n2 = nk; n1c = nic; n2c = nkc;
  }
  else  // YZ-plane
  {
    n1 = nj; n2 = nk; n1c = njc; n2c = nkc;
  }

  E_Int n1n2c = n1*n2c;
  E_Int nbInt1 = n1n2c;
  E_Int nbInt2 = n1c*n2;
  E_Int nbIntTot = nbInt1 + nbInt2;
  E_Int ncells = n1c*n2c;

  std::cout << "ni = " << ni << ", nj = " << nj << ", nk = " << nk << std::endl;
  std::cout << "n1 = " << n1 << ", n2 = " << n2 << std::endl;
  std::cout << "nbInt1 = " << nbInt1 << ", nbInt2 = " << nbInt2 << ", nbIntTot = " << nbIntTot << std::endl;

  E_Int api = fc.getApi();
  E_Int nfld = fc.getNfld();  // number of scalar fields
  E_Int nfld2 = 1; // TODO nfld/2;  // number of vector fields

  FldArrayF sint(nbIntTot,3); sint.setAllValuesAtNull();
  E_Float* sxint = sint.begin(1);
  E_Float* syint = sint.begin(2);
  E_Float* szint = sint.begin(3);

  // Compute length of faces
  #pragma omp parallel
  {
    E_Int indf, indm, indp;

    #pragma omp for nowait collapse(2)
    for (E_Int j = 0; j < n2c; j++)
    for (E_Int i = 0; i < n1; i++)
    {
      indm = i + j*ni; indp = indm + ni;
      d13x = xt[indp] - xt[indm];
      d13y = yt[indp] - yt[indm];
      d13z = 1;
      d24x = xt[indm] - xt[indp];
      d24y = yt[indm] - yt[indp];
      d24z = 1;

      sxint[indm] = 0.5*(d13y*d24z - d13z*d24y);
      syint[indm] = 0.5*(d13z*d24x - d13x*d24z);
      szint[indm] = 0.5*(d13x*d24y - d13y*d24x);
    }

    #pragma omp for collapse(2)
    for (E_Int j = 0; j < n2; j++)
    for (E_Int i = 0; i < n1c; i++)
    {
      indm = i + j*ni; indp = indm + 1;
      d13x = xt[indp] - xt[indm];
      d13y = yt[indp] - yt[indm];
      d13z = 1;
      d24x = xt[indp] - xt[indm];
      d24y = yt[indp] - yt[indm];
      d24z = -1;

      indf = i + j*nic + ninjc;
      sxint[indf] = 0.5*(d13y*d24z - d13z*d24y);
      syint[indf] = 0.5*(d13z*d24x - d13x*d24z);
      szint[indf] = 0.5*(d13x*d24y - d13y*d24x);
    }
  }

  if (cellNp == NULL)
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fcn = fc.begin(n);
      E_Float* fintp = faceField.begin(n);

      #pragma omp parallel
      {
        E_Int indf, indcellg, indcelld;

        // faces en i internes
        #pragma omp for nowait collapse(2)
        for (E_Int j = 0; j < n2c; j++)
        for (E_Int i = 1; i < n1c; i++)
        {
          indf = i + j*n1;
          indcellg = (i - 1) + j*n1c; indcelld = indcellg + 1;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }

        // bords des faces en i
        #pragma omp for nowait
        for (E_Int j = 0; j < n2c; j++)
        {
          // faces i = 0
          indf = j*n1;
          indcelld = j*n1c;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          // faces i = n1
          indf = (n1 - 1) + j*n1;
          indcellg = (n1c - 1) + j*n1c;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }

        // faces en j internes
        #pragma omp for nowait collapse(2)
        for (E_Int j = 1; j < n2c; j++)
        for (E_Int i = 0; i < n1c; i++)
        {
          indf = i + j*n1c + nbInt1;
          indcellg = i + (j - 1)*n1c; indcelld = indcellg + n1c;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          fintp[indf] = 0.5*(fcn[indcellg]+fcn[indcelld]);
        }

        // bords des faces en j
        #pragma omp for
        for (E_Int i = 0; i < n1c; i++)
        {
          // faces j = 0
          indf = i + nbInt1;
          indcelld = i;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          // faces j = jmax
          indf = i + n2c*n1c + nbInt1;
          indcellg = i + (n2c - 1)*n1c;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }
      }
    }
  }
  else // cellN
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fcn = fc.begin(n);
      E_Float* fintp = faceField.begin(n);
        
      #pragma omp parallel
      {
        E_Int indf, indcellg, indcelld;
        
        // faces en i internes
        #pragma omp for nowait collapse(2)
        for (E_Int j = 0; j < n2c; j++)
        for (E_Int i = 1; i < n1c; i++)
        {
          indf = i + j*n1;
          indcellg = (i - 1) + j*n1c; indcelld = indcellg + 1;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indf] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indf] = fcn[indcellg];
          else fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait
        for (E_Int j = 0; j < n2c; j++)
        {
          // faces i = 0
          indf = j*n1;
          indcelld = j*n1c;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          // faces i = n1
          indf = (n1 - 1) + j*n1;
          indcellg = (n1c - 1) + j*n1c;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }

        // faces en j internes
        #pragma omp for nowait collapse(2)
        for (E_Int j = 1; j < n2c; j++)
        for (E_Int i = 0; i < n1c; i++)
        {
          indf = i + j*n1c + nbInt1;
          indcellg = i + (j - 1)*n1c; indcelld = indcellg + n1c;
          cellG[indf] = indcellg; cellD[indf] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indf] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indf] = fcn[indcellg];
          else fintp[indf] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for
        for (E_Int i = 0; i < n1c; i++)
        {
          // faces j = 0
          indf = i + nbInt1;
          indcelld = i;
          cellG[indf] = -1; cellD[indf] = indcelld;
          fintp[indf] = fcn[indcelld]; // Extrapolation de l interieur

          // faces j = jmax
          indf = i + n2c*n1c + nbInt1;
          indcellg = i + (n2c - 1)*n1c;
          cellG[indf] = indcellg; cellD[indf] = -1;
          fintp[indf] = fcn[indcellg]; // Extrapolation de l interieur
        }
      }
    }
  }

  // Replace DataSet
  if (indices != Py_None && (fieldX != Py_None || fieldY != Py_None || fieldZ != Py_None))
  {
    FldArrayI* inds = NULL; FldArrayF* bfieldX = NULL;
    FldArrayF* bfieldY = NULL; FldArrayF* bfieldZ = NULL;
    K_NUMPY::getFromNumpyArray(indices, inds);
    if (ixyz != 2) K_NUMPY::getFromNumpyArray(fieldX, bfieldX);
    if (ixyz != 1) K_NUMPY::getFromNumpyArray(fieldY, bfieldY);
    if (ixyz != 0) K_NUMPY::getFromNumpyArray(fieldZ, bfieldZ);

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pindint = inds->begin();

    E_Float* pf[2];
    if (ixyz == 0) { pf[0] = bfieldX->begin(); pf[1] = bfieldY->begin(); }
    if (ixyz == 1) { pf[0] = bfieldX->begin(); pf[1] = bfieldZ->begin(); }
    if (ixyz == 2) { pf[0] = bfieldY->begin(); pf[1] = bfieldZ->begin(); }

    #pragma omp parallel
    {
      E_Int indf;
      for (E_Int n = 0; n < nfld; n++)
      {
        E_Float* fintp = faceField.begin(n+1);
        #pragma omp for
        for (E_Int noint = 0; noint < ninterfaces; noint++)
        {
          indf = pindint[noint];
          fintp[indf] = pf[n][noint];
        }
      }
    }
    RELEASESHAREDN(indices, inds);
    if (ixyz != 2) RELEASESHAREDN(fieldX, bfieldX);
    if (ixyz != 1) RELEASESHAREDN(fieldY, bfieldY);
    if (ixyz != 0) RELEASESHAREDN(fieldZ, bfieldZ);
  }

  // Build empty array
  PyObject* tpl = K_ARRAY::buildArray3(nfld2, varStringOut, nic, njc, nkc, api);
  FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2);
  f2->setAllValuesAtNull();

  // divergence
  E_Int indcellg, indcelld;
  E_Float ffi, ffj;
  for (E_Int n = 0; n < nfld2; n++)
  {
    E_Float* f2p = f2->begin(n+1);
    E_Float* fpi = faceField.begin(2*n+1);
    E_Float* fpj = faceField.begin(2*n+2);
    for (E_Int i = 0; i < nbIntTot; i++)
    {
      indcellg = cellG[i]; indcelld = cellD[i];
      ffi = fpi[i]; ffj = fpj[i];
      std::cout << "indcellg = " << indcellg << ", indcelld = " << indcelld << ", ffi = " << ffi << ", ffj = " << ffj << std::endl;
      if (indcellg != -1) f2p[indcellg] += ffi*s1int[i] + ffj*s2int[i];
      if (indcelld != -1) f2p[indcelld] -= ffi*s1int[i] + ffj*s2int[i];
      if (indcellg != -1) std::cout << "f2p[indcellg] = " << f2p[indcellg] << std::endl;
      if (indcelld != -1) std::cout << "f2p[indcelld] = " << f2p[indcelld] << std::endl;
    }
  }
  // free mem
  sint.malloc(0); faceField.malloc(0);
  
  #pragma omp parallel
  {
    E_Float voli;
    for (E_Int n = 1; n <= nfld2; n++)
    {
      E_Float* f2p = f2->begin(n);
      #pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        K_METRIC::compVolOfStructCell2D(n1, n2, i, -1, xt, yt, zt, voli);
        voli = 1./K_FUNC::E_max(voli, K_CONST::E_MIN_VOL);
        f2p[i] *= voli;
      }
    }
  }
  
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

//=============================================================================
/* From the initial chain of variables: (x,y,z,var1,var2,...)
   Create the chain (divvar1, divvar2, ....)
   This routine allocates varStringOut */
//=============================================================================
void K_POST::computeDivVarsString(char* varString, char*& varStringOut)
{
  std::vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int c = -1;
  E_Int varsSize = vars.size();
  E_Int sizeVarStringOut = 0;
  for (E_Int v = 0; v < varsSize; v++)
  {
    E_Int vsize = strlen(vars[v]);
    sizeVarStringOut += vsize+4;
  }
  varStringOut = new char [sizeVarStringOut];

  for (E_Int v = 0; v < varsSize; v++)
  {
    char*& var0 = vars[v];
    if (strcmp(var0, "x") != 0 &&
        strcmp(var0, "y") != 0 &&
        strcmp(var0, "z") != 0)
    {
      if (c == -1)
      {
        strcpy(varStringOut, "div");
        c = 1;
      }
      else strcat(varStringOut, ",div");
      char* nme = new char [strlen(var0)+1];
      strcpy(nme, var0);
      nme[strlen(var0)-1] = '\0';
      strcat(varStringOut, nme);
      delete [] nme;
    }
  }
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];
}
