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

# include "transform.h"

//=============================================================================
// splitSharpEdges:
// Decoupe un array TRI, QUAD, BAR en partie lisses (angle entre elements
// inferieur a alphaRef)
// La connectivite doit etre propre.
//=============================================================================
PyObject* K_TRANSFORM::splitSharpEdges(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float alphaRef;
  E_Float dirVect[3];
  dirVect[0] = 0.; dirVect[1] = 0.; dirVect[2] = 1.;

  if (!PYPARSETUPLE_(args, O_ R_, &array, &alphaRef)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a structured array.");
    return NULL;
  }
  else if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: unknown type of array.");
    return NULL;
  }

  if (K_STRING::cmp(eltType, "NODE") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a NODE array.");
    return NULL;
  }
  else if (K_STRING::cmp(eltType, "MIXED") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a MIXED array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  PyObject* tpl;
  if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    tpl = splitSharpEdgesNGon(f, cn, varString,
                              posx, posy, posz, alphaRef, dirVect);
  }
  else
  {
    tpl = splitSharpEdgesBasics(f, cn, eltType, varString,
                                posx, posy, posz, alphaRef, dirVect);
  }

  RELEASESHAREDU(array, f, cn);
  return tpl;
}

//=============================================================================
// Split sharp edges pour un array a elements basiques (BAR,TRI,QUAD)
//=============================================================================
PyObject* K_TRANSFORM::splitSharpEdgesBasics(
  FldArrayF* f, FldArrayI* cn,
  char* eltType, char* varString, E_Int posx, E_Int posy, E_Int posz,
  E_Float alphaRef, E_Float* dirVect)
{
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  E_Int api = f->getApi();
  E_Int npts = f->getSize();

  E_Int nc = cn->getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  E_Int dim = K_CONNECT::getDimME(eltTypes);
  if (dim == 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a 3D ME array.");
    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
    return NULL;
  }

  // Compute total number of elements
  E_Int ntotElts = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }

  std::vector<std::vector<E_Int> > cEEN(ntotElts);
  K_CONNECT::connectEV2EENbrs(eltType, npts, *cn, cEEN);

  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ntotElts, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ntotElts * sizeof(E_Int));
  E_Int mbv, ie, elt, curr, p;
  std::vector<FldArrayI*> components;

  E_Float alpha;
  E_Int ind;
  E_Float pts1[4][3]; E_Float pts2[4][3];
  
  std::vector<E_Int> type;
  if (dim == 1) type.push_back(0);
  else
  {
    for (size_t ic = 0; ic < eltTypes.size(); ic++)
    {
      if (K_STRING::cmp(eltTypes[ic], "TRI") == 0) type.push_back(1);
      else type.push_back(2);  // QUAD
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  K_FLD::FldArrayI& cm = *(cn->getConnect(0));
  E_Int nvpe = cm.getNfld();

  mbv = 0;
  while (nev < ntotElts)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant
    FldArrayI* cn2 = new FldArrayI(ntotElts, nvpe);

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    curr = 0;

    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      for (E_Int i = 0; i < nvpe; i++)
      {
        ind = cm(elt,i+1)-1;
        pts1[i][0] = x[ind]; pts1[i][1] = y[ind]; pts1[i][2] = z[ind];
        (*cn2)(curr,i+1) = ind+1;
      }
      curr++;

      for (size_t iv = 0; iv < cEEN[elt].size(); iv++)
      {
        ie = cEEN[elt][iv];
        if (isVisited[ie] == 0)
        {
          // Calcul de alpha
          for (E_Int i = 0; i < nvpe; i++)
          {
            ind = (*cn)(ie,i+1)-1;
            pts2[i][0] = x[ind]; pts2[i][1] = y[ind]; pts2[i][2] = z[ind];
          }
          switch (type[0])
          {
            case 0:
              alpha = K_COMPGEOM::getAlphaAngleBetweenBars(
              pts1[0], pts1[1],
              pts2[0], pts2[1], dirVect);
              break;

            case 1:
              alpha = K_COMPGEOM::getAlphaAngleBetweenTriangles(
                pts1[0], pts1[1], pts1[2],
                pts2[0], pts2[1], pts2[2]); break;

            case 2:
              alpha = K_COMPGEOM::getAlphaAngleBetweenQuads(
                pts1[0], pts1[1], pts1[2], pts1[3],
                pts2[0], pts2[1], pts2[2], pts2[3]); break;

            default: alpha = 180.; break;
          }

          if (alpha == -1000. || K_FUNC::E_abs(alpha-180.) <= alphaRef)
          {
            mustBeVisited[mbv] = ie;
            mbv++; nev++;
            isVisited[ie] = 1;
          }
        }
      }
    }
    cn2->reAllocMat(curr, nvpe);
    components.push_back(cn2);
  }
  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  E_Int size = components.size();

  for (E_Int i = 0; i < size; i++)
  {
    FldArrayF* f0 = new FldArrayF(*f);
    FldArrayF& fp = *f0;
    FldArrayI& cnp = *components[i];
    K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, fp, cnp);
    tpl = K_ARRAY::buildArray3(fp, varString, cnp, eltType, api);
    delete &fp; delete &cnp;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
    // tpl = K_CONNECT::V_cleanConnectivity(varString, *f, *components[i], eltType, 1.e-10);
    // PyList_Append(l, tpl); Py_DECREF(tpl);
    // delete components[i];
  }
  return l;
}

//=============================================================================
// Split sharp edges pour un array NGON 1D ou 2D
//==============================================================================
PyObject* K_TRANSFORM::splitSharpEdgesNGon(
  FldArrayF* f, FldArrayI* cn, char* varString,
  E_Int posx, E_Int posy, E_Int posz, E_Float alphaRef, E_Float* dirVect)
{
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  E_Int dim = cn->getDim();
  if (dim == 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a 3D NGON array.");
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int* ngon = cn->getNGon(); E_Int* nface = cn->getNFace();
  E_Int* indPG = cn->getIndPG(); E_Int* indPH = cn->getIndPH();
  E_Int nelts = cn->getNElts(); E_Int nfaces = cn->getNFaces();
  E_Int ngonType = cn->getNGonType();
  E_Int shift = 1; if (ngonType == 3) shift = 0;

  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  // Commence par calculer alpha
  FldArrayF alphat(nfaces);
  E_Float* alphap = alphat.begin();

  #pragma omp parallel
  {
    E_Int e1, e2, ind, nv, nf, indf, nvert;
    E_Float alpha;
    std::vector<E_Int> indices; indices.reserve(32);
    std::vector<E_Float*> pts1; pts1.reserve(32);
    std::vector<E_Float*> pts2; pts2.reserve(32);

    #pragma omp for
    for (E_Int i = 0; i < nfaces; i++)
    {
      alpha = -1000.;
      e1 = cFE1[i]-1; e2 = cFE2[i]-1;
      if (e1 != -1 && e2 != -1)
      {
        // Get vertex indices of the element e1 (unique, sorted)
        indices.clear();
        E_Int* elt = cn->getElt(e1, nf, nface, indPH);
        for (E_Int f = 0; f < nf; f++)
        {
          indf = elt[f]-1;
          E_Int* face = cn->getFace(indf, nv, ngon, indPG);
          for (E_Int k = 0; k < nv; k++) indices.push_back(face[k]);
        }
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
        nvert = indices.size();

        for (E_Int k = 0; k < nvert; k++)
        {
          E_Float* pt = new E_Float[3];
          ind = indices[k]-1;
          pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
          pts1.push_back(pt);
        }

        // Get vertex indices of the element e2 (unique, sorted)
        indices.clear();
        elt = cn->getElt(e2, nf, nface, indPH);
        for (E_Int f = 0; f < nf; f++)
        {
          indf = elt[f]-1;
          E_Int* face = cn->getFace(indf, nv, ngon, indPG);
          for (E_Int k = 0; k < nv; k++) indices.push_back(face[k]);
        }
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
        nvert = indices.size();

        for (E_Int k = 0; k < nvert; k++)
        {
          E_Float* pt = new E_Float[3];
          ind = indices[k]-1;
          pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
          pts2.push_back(pt);
        }

        alpha = 180.;
        if (dim == 1)
        {
          if (pts1.size() == 2 && pts2.size() == 2)
          {
            alpha = K_COMPGEOM::getAlphaAngleBetweenBars(
              pts1[0], pts1[1],
              pts2[0], pts2[1], dirVect
            );
          }
        }
        else if (dim == 2)
        {
          alpha = K_COMPGEOM::getAlphaAngleBetweenPolygons(pts1, pts2);
        }

        for (size_t k = 0; k < pts1.size(); k++) delete [] pts1[k];
        for (size_t k = 0; k < pts2.size(); k++) delete [] pts2[k];
        pts1.clear(); pts2.clear();
      }

      alphap[i] = alpha;
    }
  }

  // Assign a component id to each element
  E_Int* mustBeVisited = (E_Int*)malloc(nelts * sizeof(E_Int));
  std::vector<E_Int> componentId(nelts, -1);
  E_Int ncomponents = 0;
  E_Int nev = 0;  // number of elements visited already
  E_Int nextSeed = 0;  // first element that has not been visited yet
  E_Int mbv = 0;  // short for must be visited
  E_Int e1, e2, nei, nf, nv, seed, ie, indf;
  E_Float alpha;

  while (nev < nelts)
  {
    // Find first unvisited element
    for (seed = nextSeed; componentId[seed] != -1; seed++);
    nextSeed = seed + 1;

    mustBeVisited[mbv++] = seed;
    componentId[seed] = ncomponents;
    nev++;

    while (mbv > 0)
    {
      ie = mustBeVisited[--mbv];
      E_Int* elt = cn->getElt(ie, nf, nface, indPH);

      for (E_Int j = 0; j < nf; j++)
      {
        indf = elt[j] - 1;
        e1 = cFE1[indf] - 1;
        e2 = cFE2[indf] - 1;
        nei = (e1 == ie) ? e2 : e1;

        if (nei == -1 || componentId[nei] != -1) continue;

        alpha = alphap[indf];
        if (alpha == -1000. || K_FUNC::E_abs(alpha-180.) <= alphaRef)
        {
          mustBeVisited[mbv++] = nei;
          componentId[nei] = ncomponents;
          nev++;
        }
      }
    }

    ncomponents++;
  }
  std::cout << "dim = " << dim << std::endl;

  // Determine output sizes (pc = per component), incl. duplicates
  std::vector<E_Int> neltspc(ncomponents, 0), nfacespc(ncomponents, 0);
  std::vector<E_Int> sizeFNpc(ncomponents, 0), sizeEFpc(ncomponents, 0);

  E_Int cid;
  E_Int sizeFNIncr = (dim == 1) ? 1 : 2;

  for (E_Int i = 0; i < nelts; i++)
  {
    cid = componentId[i];
    cn->getElt(i, nf, nface, indPH);

    neltspc[cid]++;
    nfacespc[cid] += nf;
    sizeEFpc[cid] += nf + shift;
    for (E_Int j = 0; j < nf; j++) sizeFNpc[cid] += sizeFNIncr + shift;
  }

  // Build new arrays
  std::cout << "ncomponents = " << ncomponents << std::endl;
  PyObject* l = PyList_New(0);
  PyObject* tpl; PyObject* tplClean;
  FldArrayF* f2; FldArrayI* cn2;

  for (E_Int cid = 0; cid < ncomponents; cid++)
  {
    std::cout << "neltspc[cid] = " << neltspc[cid] << std::endl;
    std::cout << "nfacespc[cid] = " << nfacespc[cid] << std::endl;
    std::cout << "sizeFNpc[cid] = " << sizeFNpc[cid] << std::endl;
    std::cout << "sizeEFpc[cid] = " << sizeEFpc[cid] << std::endl;
    tpl = K_ARRAY::buildArray3(*f, varString, neltspc[cid], nfacespc[cid],
                               "NGON", sizeFNpc[cid], sizeEFpc[cid],
                               ngonType, false, api);

    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (ngonType == 2 || ngonType == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }

    E_Int c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    if (ngonType == 2 || ngonType == 3) { indPG2[0] = 0; indPH2[0] = 0; }
    for (E_Int i = 0; i < nelts; i++)
    {
      if (componentId[i] != cid) continue;  // TODO build and use an inverse map instead
      E_Int* elt = cn->getElt(i, nf, nface, indPH);

      nface2[c2] = nf;
      if (c4 < neltspc[cid]-1 && (ngonType == 2 || ngonType == 3))
      {
        indPH2[c4+1] = indPH2[c4] + nf + shift;
      }
      for (E_Int j = 0; j < nf; j++)
      {
        indf = elt[j];  // 1-based face index

        // Fill NFACE
        nface2[c2+j+shift] = indf;

        // Fill NGON
        E_Int* face = cn->getFace(indf-1, nv, ngon, indPG);
        ngon2[c1] = nv;
        if (c3 < nfacespc[cid]-1 && (ngonType == 2 || ngonType == 3))
        {
          indPG2[c3+1] = indPG2[c3] + nv + shift;
        }
        for (E_Int k = 0; k < nv; k++)
        {
          ngon2[c1+k+shift] = face[k];
        }
        c1 += nv+shift; c3++;
      }
      c2 += nf+shift; c4++;
    }

    std::cout << "ngon2 = " << std::endl;
    for (E_Int j = 0; j < sizeFNpc[cid]; j++) std::cout << ngon2[j] << ", ";
    std::cout << std::endl;
    std::cout << "nface2 = " << std::endl;
    for (E_Int j = 0; j < sizeEFpc[cid]; j++) std::cout << nface2[j] << ", ";
    std::cout << std::endl;
    if (ngonType == 2 || ngonType == 3)
    {
      std::cout << "indPG2 = " << std::endl;
      for (E_Int j = 0; j < nfacespc[cid]; j++) std::cout << indPG2[j] << ", ";
      std::cout << std::endl;
      std::cout << "indPH2 = " << std::endl;
      for (E_Int j = 0; j < neltspc[cid]; j++) std::cout << indPH2[j] << ", ";
      std::cout << std::endl;
    }

    tplClean = K_CONNECT::V_cleanConnectivity(
        varString, *f2, *cn2, "NGON", 1.e-12,
        false, true,  // rmOrphans
        true, false,  // rmDuplicatedFaces
        false, false
    );

    std::cout << "Done... " << cid << std::endl;

    RELEASESHAREDU(tpl, f2, cn2); Py_DECREF(tpl);
    PyList_Append(l, tplClean); Py_DECREF(tplClean);
  }
  std::cout << "Exiting sharp edges..." << std::endl;
  return l;
}

//=============================================================================
// Split sharp edges pour un array NGON de dim=2 ou dim=1
// Input : liste d'index et une sortie liste d'index splittee
// Chaque index est associe a un element
//==============================================================================
PyObject* K_TRANSFORM::splitSharpEdgesList(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayI;
  E_Float alphaRef;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &array, &arrayI, &alphaRef))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: can not be used on a structured array.");
    return NULL;
  }
  if (res == 2 && strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: only for NGON array.");
    return NULL;
  }

  FldArrayI* indexI;
  res = K_NUMPY::getFromNumpyArray(arrayI, indexI);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: index numpy is invalid.");
    return NULL;
  }

  E_Float dirVect[3];
  dirVect[0] = 0.; dirVect[1] = 0.; dirVect[2] = 1.;
  E_Int* index = indexI->begin();

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  E_Int nelts = cn->getNElts();
  E_Int nfaces = cn->getNFaces();
  E_Int dim = cn->getDim();
  E_Int* ngon = cn->getNGon(); E_Int* nface = cn->getNFace();
  E_Int* indPG = cn->getIndPG(); E_Int* indPH = cn->getIndPH();

  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);

  E_Int se = cn->getSizeNFace();
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(nelts, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(nelts * sizeof(E_Int));
  E_Int mbv, p, ie, elt, necurr, lt;
  std::vector<FldArrayI*> components; // liste des index par bloc

  // Commence par calculer alpha
  FldArrayF alphat(nfaces);
  E_Float* alphap = alphat.begin();
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  #pragma omp parallel
  {
    E_Int e1, e2;
    std::vector<E_Int> indices;
    std::vector<E_Float*> pts1; std::vector<E_Float*> pts2;
    E_Float alpha;
    E_Int ind, nvert;

    #pragma omp for
    for (E_Int i = 0; i < nfaces; i++)
    {
      e1 = cFE1[i]-1; e2 = cFE2[i]-1;
      if (e1 == -1) alpha = -1000.;
      else if (e2 == -1) alpha = -1000.;
      else
      {
        K_CONNECT::getVertexIndices(*cn, ngon, nface, indPG, indPH, e1, indices);
        nvert = indices.size();
        pts1.reserve(nvert);
        for (E_Int k = 0; k < nvert; k++)
        {
          E_Float* pt = new E_Float[3];
          ind = indices[k]-1;
          pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
          pts1.push_back(pt);
        }

        K_CONNECT::getVertexIndices(*cn, ngon, nface, indPG, indPH, e2, indices);
        nvert = indices.size();
        pts2.reserve(nvert);
        for (E_Int k = 0; k < nvert; k++)
        {
          E_Float* pt = new E_Float[3];
          ind = indices[k]-1;
          pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
          pts2.push_back(pt);
        }

        alpha = 180.;
        if (dim == 2)
        {
          alpha = K_COMPGEOM::getAlphaAngleBetweenPolygons(pts1, pts2);
        }
        else if (dim == 1)
        {
          if (pts1.size() == 2 && pts2.size() == 2)
          {
            alpha = K_COMPGEOM::getAlphaAngleBetweenBars(
                      pts1[0], pts1[1],
                      pts2[0], pts2[1], dirVect);
          }
        }

        nvert = pts1.size();
        for (E_Int k = 0; k < nvert; k++) delete [] pts1[k];
        pts1.clear();
        nvert = pts2.size();
        for (E_Int k = 0; k < nvert; k++) delete [] pts2[k];
        pts2.clear();
      }
      alphap[i] = alpha;
    }
  }

  E_Int e1, e2, nf;
  E_Float alpha;

  // split
  mbv = 0;
  while (nev < nelts)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant (morceau de liste)
    FldArrayI* cn2 = new FldArrayI(se+1);
    E_Int* pc = cn2->begin(); // current pointer

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    necurr = 0;

    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];

      // copie index de l'element
      E_Int* elem = cn->getElt(elt, lt, nface, indPH);
      (*pc) = index[elt]; pc += 1;
      necurr++;

      for (E_Int iv = 0; iv < lt; iv++)
      {
        nf = elem[iv]-1;
        e1 = cFE1[nf]-1; e2 = cFE2[nf]-1;
        if (e1 == elt) ie = e2;
        else ie = e1;

        if (ie != -1 && isVisited[ie] == 0)
        {
          alpha = alphap[nf];
          if (alpha == -1000. || K_FUNC::E_abs(alpha-180.) <= alphaRef)
          {
            mustBeVisited[mbv] = ie;
            mbv++; nev++;
            isVisited[ie] = 1;
          }
        }
      }
    }

    cn2->reAlloc(necurr);
    components.push_back(cn2);
  }

  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  E_Int size = components.size();

  for (E_Int i = 0; i < size; i++)
  {
    FldArrayI* c = components[i];
    E_Int nd = c->getSize();
    E_Int* pc = c->begin();
    tpl = K_NUMPY::buildNumpyArray(nd, 1, 1, 1);
    E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
    for (E_Int k = 0; k < nd; k++) pos[k] = pc[k];
    delete components[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDN(arrayI, indexI);
  return l;
}
