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

// selectCellCenters
# include "stdio.h"
# include "post.h"
# include "Nuga/include/ngon_t.hxx"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Selectionne les cellules d'un array dont les centres sont tagues */
// ============================================================================
PyObject* K_POST::selectCellCenters(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* taga;
  PyObject* PE; E_Int cleanConnectivity;
  if (!PYPARSETUPLE_(args, OOO_ I_, &array, &taga, &PE, &cleanConnectivity)) return NULL;

  // Extract array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res, ni, nj, nk;
  res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, 
                               eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCellCenters: array is invalid.");
    return NULL;
  }

  // Extract tag
  char* varStringa; char* eltTypea;
  FldArrayF* tag; FldArrayI* cnpa;
  E_Int nia, nja, nka;
  E_Int resa = K_ARRAY::getFromArray3(taga, varStringa, 
                                      tag, nia, nja, nka, cnpa, 
                                      eltTypea);

  if (resa != 1 && resa != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCellCenters: tag is invalid.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }

  if (res != resa)
  {  
    PyErr_SetString(PyExc_TypeError,
                    "selectCellCenters: tag and array must represent the same grid.");
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resa, taga, tag, cnpa);
    return NULL;
  }

  if (tag->getNfld() != 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(resa, taga, tag, cnpa);
    PyErr_SetString(PyExc_TypeError,
                    "selectCellCenters: tag must have only one variable.");
    return NULL;
  }

  // Check coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;

  const E_Float eps = 1.e-10;
  const E_Float oneEps = 1.0 - eps;
  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  E_Float* tagp = tag->begin();
  PyObject* l = PyList_New(0);
  PyObject* tpl = NULL;
  FldArrayF* f2; FldArrayI* cn2;
  char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH]; eltType2[0] = '\0';

  if (res == 1) // Create connectivity (BAR, QUAD or HEXA)
  {
    E_Int dim0 = 3;
    if (ni == 1 || nj == 1 || nk == 1) dim0 = 2;
    else if (nj == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nj == 1) dim0 = 1;
    
    if (dim0 == 3) strcpy(eltType2, "HEXA");
    else if (dim0 == 2) strcpy(eltType2, "QUAD");
    else strcpy(eltType2, "BAR");

    E_Int ni1 = E_max(1, E_Int(ni)-1);
    E_Int nj1 = E_max(1, E_Int(nj)-1);
    E_Int nk1 = E_max(1, E_Int(nk)-1);
    E_Int ninj = ni*nj;
    E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees

    if (tag->getSize() != ncells)
    {
      RELEASESHAREDB(res, array, f, cn);
      RELEASESHAREDB(resa, taga, tag, cnpa);
      PyErr_SetString(PyExc_TypeError,
                      "selectCellCenters: dimensions of tag are not valid.");
      return NULL;
    }

    // Create a cell indirection mask: 1 (tagged), 0 (not tagged)
    std::vector<E_Int> eindir(ncells);
    #pragma omp parallel for
    for (E_Int i = 0; i < ncells; i++)
    {
      if (tagp[i] >= oneEps) eindir[i] = 1;
      else eindir[i] = 0;
    }

    // Transform the cell indirection mask of zeros and ones into an element map
    // from old to new connectivities, and get the number of new elements, nelts2
    // Careful: eindir starts at 1
    E_Int nelts2 = K_CONNECT::prefixSum(eindir);

    // Build new BE connectivity
    tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts2, eltType2, false, api);
    
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    FldArrayI& cm2 = *(cn2->getConnect(0));
    
    #pragma omp parallel
    {
      E_Int c, d, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
      
      if (dim0 == 1)
      {
        if (nj1 == 1 && nk1 == 1)
        {
          #pragma omp for
          for (E_Int i = 0; i < ni1; i++)
          {
            c = eindir[i]-1;
            if (c == 0) continue;
            // starts at 1
            ind1 = i + 1;           // (i,1,1)
            ind2 = ind1 + 1;        // (i+1,1,1)
            cm2(c,1) = ind1; cm2(c,2) = ind2;
          }
        }
        else if (ni1 == 1 && nk1 == 1)
        {
          #pragma omp for
          for (E_Int j = 0; j < nj1; j++)
          {
            c = eindir[j]-1;
            if (c == 0) continue;
            ind1 = j*ni + 1;        // (1,j,1)
            ind2 = ind1 + ni;       // (1,j+1,1)
            cm2(j,1) = ind1; cm2(j,2) = ind2;        
          }
        }
        else
        {
          #pragma omp for
          for (E_Int k = 0; k < nk1; k++)
          {
            c = eindir[k]-1;
            if (c == 0) continue;
            ind1 = 1 + k*ninj;      // (1,1,k)
            ind2 = ind1 + ninj;     // (1,1,k+1)
            cm2(k,1) = ind1; cm2(k,2) = ind2;
          }
        }
      }
      else if (dim0 == 2)
      {
        if (nk1 == 1)
        {
          #pragma omp for collapse(2)
          for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            d = j*ni1 + i;
            c = eindir[d]-1;
            if (c == 0) continue;
            // starts at 1
            ind1 = i + j*ni + 1;    // (i,j,1)
            ind2 = ind1 + 1;        // (i+1,j,1)
            ind3 = ind2 + ni;       // (i+1,j+1,1)
            ind4 = ind3 - 1;        // (i,j+1,1)
            cm2(i,1) = ind1; cm2(c,2) = ind2;
            cm2(c,3) = ind3; cm2(c,4) = ind4;
          }
        }
        else if (nj1 == 1)
        {
          #pragma omp for collapse(2)
          for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            d = k*ni1 + i;
            c = eindir[d]-1;
            if (c == 0) continue;
            ind1 = i + k*ninj + 1;  // (i,1,k)
            ind2 = ind1 + ninj;     // (i,1,k+1)
            ind3 = ind2 + 1;        // (i+1,1,k+1)
            ind4 = ind3 - 1;        // (i,1,k+1)
            cm2(c,1) = ind1; cm2(c,2) = ind2;
            cm2(c,3) = ind3; cm2(c,4) = ind4;
          }
        }
        else // i1 = 1 
        {
          #pragma omp for collapse(2)
          for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj1; j++)
          {
            d = k*nj1 + j;
            c = eindir[d]-1;
            if (c == 0) continue;
            ind1 = 1 + j*ni + k*ninj;  // (1,j,k)
            ind2 = ind1 + ni;          // (1,j+1,k)
            ind3 = ind2 + ninj;        // (1,j+1,k+1)
            ind4 = ind3 - ni;          // (1,j,k+1)
            cm2(c,1) = ind1; cm2(c,2) = ind2;
            cm2(c,3) = ind3; cm2(c,4) = ind4;
          }
        }
      }
      else 
      { 
        #pragma omp for collapse(3)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          d = k*nj1*ni1 + j*ni1 + i;
          c = eindir[d]-1;
          if (c == 0) continue;
          ind1 = 1 + i + j*ni + k*ninj; // A(  i,  j,k)
          ind2 = ind1 + 1;              // B(i+1,  j,k)
          ind3 = ind2 + ni;             // C(i+1,j+1,k)
          ind4 = ind3 - 1;              // D(  i,j+1,k)
          ind5 = ind1 + ninj;           // E(  i,  j,k+1)
          ind6 = ind2 + ninj;           // F(i+1,  j,k+1)
          ind7 = ind3 + ninj;           // G(i+1,j+1,k+1)
          ind8 = ind4 + ninj;           // H(  i,j+1,k+1) 
          cm2(c,1) = ind1; cm2(c,2) = ind2;
          cm2(c,3) = ind3; cm2(c,4) = ind4;
          cm2(c,5) = ind5; cm2(c,6) = ind6;
          cm2(c,7) = ind7; cm2(c,8) = ind8;
        }
      }
    }
  }
  else if (K_STRING::cmp(eltType, "NGON") != 0) // BE/ME
  {
    E_Int nc = cn->getNConnect();
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);

    std::vector<E_Int> nepc(nc);
    std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;  // cumulative number of elements per conn.
    // Compute total number of elements across all connectivities
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      E_Int nelts = cm.getSize();
      nepc[ic] = nelts;
      cumnepc[ic+1] = cumnepc[ic] + nelts;
    }
    E_Int ntotElts = cumnepc[nc];

    // Create a cell indirection mask: 1 (tagged), 0(not tagged)
    std::vector<E_Int> eindir(ntotElts);
    #pragma omp parallel for
    for (E_Int i = 0; i < ntotElts; i++)
    {
      if (tagp[i] >= oneEps) eindir[i] = 1;
      else eindir[i] = 0;
    }

    // Transform the cell indirection mask of zeros and ones into an element map
    // from old to new connectivities, and get the number of new elements per
    // connectivity, tmp_nepc2
    std::vector<E_Int> tmp_nepc2 = K_CONNECT::prefixSum(eindir, nepc);

    // Build new eltType from connectivities that have at least one element
    E_Int nc2 = 0;
    std::map<E_Int, E_Int> outConnId;  // map ic to ic2
    for (E_Int ic = 0; ic < nc; ic++)
    {
      if (tmp_nepc2[ic] > 0)
      {
        outConnId[ic] = nc2; nc2++;
        if (eltType2[0] == '\0') strcpy(eltType2, eltTypes[ic]);
        else
        {
          strcat(eltType2, ",");
          strcat(eltType2, eltTypes[ic]);
        }
      }
    }
    if (nc2 > 1) api = 3;
    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

    // Compress the number of elements per connectivity of the output ME, ie,
    // drop connectivities containing no tagged elements
    std::vector<E_Int> nepc2(nc2);
    nc2 = 0;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      if (tmp_nepc2[ic] > 0) { nepc2[nc2] = tmp_nepc2[ic]; nc2++; }
    }

    // Build new ME connectivity
    tpl = K_ARRAY::buildArray3(nfld, varString, npts, nepc2, eltType2, false, api);
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    #pragma omp parallel
    {
      E_Int c, d, nelts, nvpe, offset = 0;

      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn->getConnect(ic));
        FldArrayI& cm2 = *(cn2->getConnect(ic));
        nelts = cm.getSize();
        nvpe = cm.getNfld();

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          d = i + offset;
          c = eindir[d]-1;
          if (c == 0) continue;
          for (E_Int j = 1; j <= nvpe; j++) cm2(c,j) = cm(i,j);
        }

        offset += nelts;
      }
    }
  }
  else  // NGON
  {
    std::cout << "B" << std::endl;
    E_Int* cnpp = cn->begin();
    E_Int sizeFN = cn->getSizeNGon();            // taille de l'ancienne connectivite Face/Noeuds
    E_Int nbElements = cn->getNElts(); // nombre d'elts de l'ancienne connectivite
    E_Int sizeEF = cn->getSizeNFace();     // taille de l'ancienne connectivite Elmt/Faces
    FldArrayI cn2 = FldArrayI(sizeEF); // nouvelle connectivite Elmt/Faces
    E_Int* cn2p = cn2.begin();         // pointeur sur la nouvelle connectivite
    E_Int size2 = 0;                   // compteur pour la nouvelle connectivite
    E_Int next = 0;                      // nbre d'elts selectionnes
    E_Int nbFaces = cn->getNFaces();
    E_Int ii = 0;

    E_Int* ngon = cn->getNGon(); E_Int* nface = cn->getNFace();
    E_Int* indPG = cn->getIndPG(); E_Int* indPH = cn->getIndPH();
 
    E_Int newNumFace = 0;  
    FldArrayI new_pg_ids;
    FldArrayI keep_pg;
    FldArrayI new_ph_ids;
    
    if (PE != Py_None)
    {
      new_pg_ids.malloc(nbFaces);    // Tableau d'indirection des faces (pour maj PE)
      keep_pg.malloc(nbFaces);       // Flag de conservation des faces 
      new_ph_ids.malloc(nbElements); // Tableau d'indirection des elmts (pour maj PE)

      new_pg_ids = -1;
      new_ph_ids = -1;
      keep_pg    = -1; 
      
      // Boucle sur le nombre d'elements
      for (E_Int i = 0; i < nbElements; i++)
      {
        E_Int* elt = cn->getElt(i, nbFaces, nface, indPH);
        if (tagp[i] >= oneEps)
        {
          cn2p[0] = nbFaces; size2 +=1;
          for (E_Int n = 0; n < nbFaces; n++)
          {
            cn2p[n+1] = elt[n];
            keep_pg[elt[n]-1] = 1;
          }
          size2 += nbFaces; cn2p += nbFaces+1; next++;
      
          // Selection elt 
          new_ph_ids[i] = ii;
          ii++;
        }
      }
  
      E_Int nn = 0; 
      for (E_Int n = 0; n<new_pg_ids.getSize(); n++)
      {
        if (keep_pg[n]>0){ new_pg_ids[n] = nn; nn++; newNumFace++;}
      }
  
      cn2.reAlloc(size2);
    }
    else // PE == Py_None - pas de creation de tab d'indirection 
    {
      // Boucle sur le nombre d elements
      for (E_Int i = 0; i < nbElements; i++)
      {
        E_Int* elt = cn->getElt(i, nbFaces, nface, indPH);
        if (tagp[i] >= oneEps)
        {
          cn2p[0] = nbFaces; size2 +=1;
          for (E_Int n = 0; n < nbFaces; n++)
          {
            cn2p[n+1] = elt[n];
          }
          size2 += nbFaces; cn2p += nbFaces+1; next++;
        }
      }
      cn2.reAlloc(size2);
    }

    // Cree la nouvelle connectivite complete
    E_Int coutsize = sizeFN+4+size2;
    FldArrayI* cout = new FldArrayI(coutsize);
    E_Int* coutp = cout->begin();
    cn2p = cn2.begin();
    for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = cnpp[i];
    coutp += sizeFN+2;
    coutp[0] = next;
    coutp[1] = size2; coutp += 2;
    for (E_Int i = 0; i < size2; i++) coutp[i] = cn2p[i];

    RELEASESHAREDB(resa, taga, tag, cnpa);

    if (PE != Py_None)
    {
      // Check numpy (parentElement)
      FldArrayI* cFE;
      E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE);
      
      if (res == 0)
      {
        RELEASESHAREDN(PE, cFE);
        PyErr_SetString(PyExc_TypeError, "selectCellsCenter: PE numpy is invalid.");
        return NULL;
      }
      
      ngon_t<K_FLD::FldArrayI> ng(*cout); // construction d'un ngon_t à partir d'un FldArrayI
      
      FldArrayI* cFEp_new = new FldArrayI(newNumFace,2);
      FldArrayI& cFE_new  = *cFEp_new ; 

      E_Int* cFEl = cFE_new.begin(1);
      E_Int* cFEr = cFE_new.begin(2);

      E_Int* cFEl_old = cFE->begin(1);
      E_Int* cFEr_old = cFE->begin(2);
      
      E_Int old_ph_1, old_ph_2;

      nbFaces = cnpp[0];

      for (E_Int pgi = 0; pgi < nbFaces; pgi++)
      {
        if (new_pg_ids[pgi]>=0)
        {
          old_ph_1 = cFEl_old[pgi]-1;
          old_ph_2 = cFEr_old[pgi]-1;

          if (old_ph_1 >= 0) // l'elmt gauche existe 
          {
            cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_1]+1;
 
            if (old_ph_2 >= 0) // l'elmt droit existe
            {
              cFEr[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
            }
            else
              cFEr[new_pg_ids[pgi]] = 0;
          } 
          else // l'elmt gauche a disparu - switch droite/gauche
          {
            cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
            cFEr[new_pg_ids[pgi]] = 0;
            // reverse
            E_Int s = ng.PGs.stride(pgi);
            E_Int* p = ng.PGs.get_facets_ptr(pgi);
            std::reverse(p, p + s);
          }
        }
      } // boucle pgi      
      
      // export ngon
      ng.export_to_array(*cout);

      // objet Python de sortie
      PyObject* pyPE = K_NUMPY::buildNumpyArray(cFE_new, 1);

      PyList_Append(l, pyPE);
      
      delete cFEp_new;
      RELEASESHAREDN(PE, cFE);
    }
  }

  // Remove orphans
  if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
  {
    PyObject* tplc = K_CONNECT::V_cleanConnectivity(
      varString, *f2, *cn2, eltType2, eps,
      false, true, false, false, false, false, false);

    RELEASESHAREDU(tpl, f2, cn2);
    Py_DECREF(tpl);
    PyList_Append(l, tplc); Py_DECREF(tplc);
  }
  else
  {
    RELEASESHAREDU(tpl, f2, cn2);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  delete [] eltType2;
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(resa, taga, tag, cnpa);
  return l;
}

//=============================================================================
/* Selectionne les cellules d'un array dont les centres sont tagues */
// ============================================================================
PyObject* K_POST::selectCellCentersBoth(PyObject* self, PyObject* args) 
{
  PyObject* arrayNodes; PyObject* arrayCenters; PyObject* taga;
  PyObject* PE; E_Int cleanConnectivity;
  if (!PYPARSETUPLE_(args, OOOO_ I_, &arrayNodes, &arrayCenters, &taga, &PE, &cleanConnectivity)) return NULL;

  // Extract arrayNodes 
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cnp;
  E_Int res, ni, nj, nk;
  res = K_ARRAY::getFromArray3(arrayNodes, varString, f, ni, nj, nk, cnp, 
                               eltType);


  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: array is invalid.");
    return NULL;
  }

 // Extract arrayCenters 
  char* varStringC; char* eltTypeC;
  FldArrayF* fC; FldArrayI* cnpC;
  E_Int resC, niC, njC, nkC;
  resC = K_ARRAY::getFromArray3(arrayCenters, varStringC, fC, niC, njC, nkC, cnpC, 
                                eltTypeC);
  
  if (resC != 1 && resC != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: arrayCenters is invalid.");
    return NULL;
  }

  
  // Extract tag
  char* varStringa; char* eltTypea;
  FldArrayF* tag; FldArrayI* cnpa;
  E_Int nia, nja, nka;
  E_Int resa = K_ARRAY::getFromArray3(taga, varStringa, 
                                      tag, nia, nja, nka, cnpa, 
                                      eltTypea);

  if (resa != 1 && resa != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag is invalid.");
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    return NULL;
  }

  if (res != resa)
  {  
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag and array must represent the same grid.");
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(resa, taga, tag, cnpa);
    return NULL;
  }

  if (tag->getNfld() != 1)
  {
    RELEASESHAREDB(res, arrayNodes, f, cnp);
    RELEASESHAREDB(resa, taga, tag, cnpa);
    PyErr_SetString(PyExc_TypeError,
                    "selectCells2: tag must have only one variable.");
    return NULL;
  }

  const E_Float oneEps = 1. - 1.e-10;
  E_Float* tagp = tag->begin();
  // no check of coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  E_Int nfld = f->getNfld();
  if (res == 1) // Create connectivity cnp (HEXA or QUAD)
  {
    E_Int dim0 = 3;
    if (ni == 1 || nj == 1 || nk == 1) dim0 = 2;
    if (nj == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nk == 1) dim0 = 1;
    else if (ni == 1 && nj == 1) dim0 = 1;
    eltType = new char [128];
    if (dim0 == 3) strcpy(eltType, "HEXA");
    else if (dim0 == 2) strcpy(eltType, "QUAD");
    else strcpy(eltType, "BAR");
    E_Int ni1 = E_max(1, E_Int(ni)-1);
    E_Int nj1 = E_max(1, E_Int(nj)-1);
    E_Int nk1 = E_max(1, E_Int(nk)-1);
    E_Int ninj = ni*nj;
    E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
    if (tag->getSize() != ncells)
    {
      RELEASESHAREDB(res, arrayNodes, f, cnp);
      RELEASESHAREDB(resa, taga, tag, cnpa);
      PyErr_SetString(PyExc_TypeError,
                      "selectCells2: dimensions of tag are not valid.");
      return NULL;
    }
    cnp = new FldArrayI();
    FldArrayI& cn = *cnp;
    E_Int nelts; // nb d'elements non structures
    
    E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
    E_Int c = 0;
    
    if (dim0 == 1)
    {
      nelts = ncells;
      cn.malloc(nelts, 2);
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);

      if (nj1 == 1 && nk1 == 1)
      {
        for (E_Int i = 0; i < ni1; i++)
        {
          // starts from 1
          ind1 = i + 1; //(i,1,1)
          ind2 = ind1 + 1;  //(i+1,1,1)
          cn1[c] = ind1; cn2[c] = ind2;  
          c++;
        }
      }
      else if (ni1 == 1 && nk1 == 1)
      {
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = j*ni + 1;  //(1,j,1)
          ind2 = ind1 + ni; //(1,j+1,1)
          cn1[c] = ind1; cn2[c] = ind2;
          c++;         
        }
      }
      else
      {
        for (E_Int k = 0; k < nk1; k++)
        {
          ind1 = 1 + k*ninj; //(1,1,k)
          ind2 = ind1 + ninj;   //(1,1,k+1)
          cn1[c] = ind1; cn2[c] = ind2;
          c++;
        }
      }
    }
    else if (dim0 == 2)
    {
      nelts = ncells;
      cn.malloc(nelts, 4);
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      E_Int* cn4 = cn.begin(4);
      if (nk1 == 1)
      {
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            //starts from 1
            ind1 = i + j*ni + 1; //(i,j,1)
            ind2 = ind1 + 1;  //(i+1,j,1)
            ind3 = ind2 + ni; //(i+1,j+1,1)
            ind4 = ind3 - 1;  //(i,j+1,1)
            cn1[c] = ind1; cn2[c] = ind2;
            cn3[c] = ind3; cn4[c] = ind4;
            c++;
          }
      }
      else if (nj1 == 1)
      {
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + k*ninj + 1;  //(i,1,k)
            ind2 = ind1 + ninj; //(i,1,k+1)
            ind3 = ind2 + 1;    //(i+1,1,k+1)
            ind4 = ind3 - 1;    //(i,1,k+1)
            cn1[c] = ind1; cn2[c] = ind2;
            cn3[c] = ind3; cn4[c] = ind4;
            c++;         
          }
      }
      else // i1 = 1 
      {
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj1; j++)
          {
            ind1 = 1 + j*ni + k*ninj; //(1,j,k)
            ind2 = ind1 + ni;   //(1,j+1,k)
            ind3 = ind2 + ninj; //(1,j+1,k+1)
            ind4 = ind3 - ni;   //(1,j,k+1)
            cn1[c] = ind1; cn2[c] = ind2;
            cn3[c] = ind3; cn4[c] = ind4;
            c++;
          }
      }// i1 = 1
    }//dim 2
    else 
    { 
      nelts = ncells;
      cn.malloc(nelts,8);
      E_Int* cn1 = cn.begin(1);
      E_Int* cn2 = cn.begin(2);
      E_Int* cn3 = cn.begin(3);
      E_Int* cn4 = cn.begin(4);
      E_Int* cn5 = cn.begin(5);
      E_Int* cn6 = cn.begin(6);
      E_Int* cn7 = cn.begin(7);
      E_Int* cn8 = cn.begin(8);
      
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = 1 + i + j*ni + k*ninj; //A(  i,  j,k)
            ind2 = ind1 + 1;              //B(i+1,  j,k)
            ind3 = ind2 + ni;             //C(i+1,j+1,k)
            ind4 = ind3 - 1;              //D(  i,j+1,k)
            ind5 = ind1 + ninj;           //E(  i,  j,k+1)
            ind6 = ind2 + ninj;           //F(i+1,  j,k+1)
            ind7 = ind3 + ninj;           //G(i+1,j+1,k+1)
            ind8 = ind4 + ninj;           //H(  i,j+1,k+1) 
            
            cn1[c] = ind1; cn2[c] = ind2;
            cn3[c] = ind3; cn4[c] = ind4;
            cn5[c] = ind5; cn6[c] = ind6;
            cn7[c] = ind7; cn8[c] = ind8;
            c++;
          }
    } //dim = 3
  }

  // Selection
  
  PyObject* l = PyList_New(0); 
  PyObject* tpl;
  PyObject* tplc;
  E_Int api = fC->getApi();

  if (strcmp(eltType, "NGON") != 0) // tous les elements sauf NGON
  {
    FldArrayF* fout  = new FldArrayF(*f);
    FldArrayF* foutC = new FldArrayF();
    FldArrayF& fcenter = *foutC; 
    FldArrayI* acn = new FldArrayI();
    FldArrayI& cn = *acn;
    E_Int nt = cnp->getNfld();
    E_Int ne = cnp->getSize();

    E_Int nthreads = __NUMTHREADS__;
    E_Int net = ne/nthreads+1;
    E_Int** ptr   = new E_Int* [nthreads];
    E_Float** ptrF = new E_Float* [nthreads];
    E_Int* nes = new E_Int [nthreads];
    E_Int* prev = new E_Int [nthreads];
    E_Int nfldC = fC->getNfld();
    
    for (E_Int i = 0; i < nthreads; i++)
    {
      ptr[i]  = new E_Int   [net*nt];
      ptrF[i] = new E_Float [nfldC*net*nt];
    }

#pragma omp parallel default(shared)
    {
      E_Int ithread = __CURRENT_THREAD__;
      nes[ithread] = 0;  
      E_Int cprev  = 0;
      E_Int cprev2 = 0;
      E_Int* cnt        = ptr[ithread];
      E_Float* ftcenter = ptrF[ithread];
      E_Int ii = 0 ; 
      
#pragma omp for
      for (E_Int i = 0; i < ne; i++)
      {
        if (tagp[i] >= oneEps)
        {
          for (E_Int n = 1; n <= nt; n++)
          {
            cnt[cprev+(n-1)]      = (*cnp)(i,n);
          }

          // boucle sur les champs
          for (E_Int k = 1; k <= nfldC; k++)
          {
            E_Float* fcpl = fC->begin(k);
            ftcenter[cprev2+(k-1)]  = fcpl[i];
          }
          ii++;
          cprev += nt; cprev2 +=nfldC; nes[ithread]++;
        }
      }
    }

    RELEASESHAREDB(resa, taga, tag, cnpa);
    if (res == 1) delete cnp;

    // total
    E_Int nntot = 0;
    for (E_Int i = 0; i < nthreads; i++) { prev[i] = nntot; nntot += nes[i]; }

    // Compact 
    cn.malloc(nntot, nt);

    fcenter.malloc(nntot, nfldC);

#pragma omp parallel default(shared)
    {
      E_Int ithread     = __CURRENT_THREAD__;
      E_Int* cnt        = ptr[ithread];
      E_Float* ftcenter = ptrF[ithread];
      E_Int p = prev[ithread];
      for (E_Int i = 0; i < nes[ithread]; i++)
      {
        for (E_Int n = 1; n <= nt; n++)
        {
          cn(i+p, n) = cnt[i*nt+(n-1)];
        }

        for (E_Int k = 1; k <= nfldC; k++)
        {
          fcenter(i+p,k) = ftcenter[i*nfldC+(k-1)];	    
        }
      }
    }

    delete [] nes; delete [] prev;
    for (E_Int i = 0; i < nthreads; i++){ delete [] ptr[i]; delete [] ptrF[i]; }
    delete [] ptr; delete [] ptrF;
    
    if (nntot == 0) fout->reAllocMat(0, nfld);
    else 
    {
      if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
        K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, 
                                     *fout, *acn);
    }
    
    tpl = K_ARRAY::buildArray3(*fout, varString, *acn, eltType, api);
    tplc = K_ARRAY::buildArray3(*foutC, varStringC, *acn, eltType, api);
    delete acn; delete fout; delete foutC;
    if (res == 1) delete[] eltType;
  }
  else // elements NGON
  {
    FldArrayF* fout     = new FldArrayF(*f);
    FldArrayF* foutC    = new FldArrayF(*fC); // pour recuperer la solution aux centres 
    FldArrayF& fcenter  = *foutC; 
    FldArrayF& fcenter0 = *fC; 
    
    E_Int nfldC = fC->getNfld(); // nombre de champs en centres 
    E_Int* cnpp = cnp->begin();
    E_Int sizeFN = cnpp[1];            // taille de l'ancienne connectivite Face/Noeuds
    E_Int nbElements = cnpp[sizeFN+2]; // nombre d'elts de l'ancienne connectivite
    E_Int sizeEF = cnpp[sizeFN+3];     // taille de l'ancienne connectivite Elmt/Faces
    E_Int* cnEFp = cnpp+4+sizeFN;    // pointeur sur l'ancienne connectivite Elmt/Faces
    FldArrayI cn2 = FldArrayI(sizeEF); // nouvelle connectivite Elmt/Faces
    E_Int* cn2p = cn2.begin();         // pointeur sur la nouvelle connectivite
    E_Int size2 = 0;                   // compteur pour la nouvelle connectivite
    E_Int next=0;                      // nbre d'elts selectionnes
    E_Int nbFaces = cnpp[0];       // nombre de faces de l'ancienne connectivite;
    E_Int ii = 0 ;

    E_Int newNumFace = 0;  
    FldArrayI new_pg_ids;
    FldArrayI keep_pg;
    FldArrayI new_ph_ids;

    if (PE != Py_None)
    {
      new_pg_ids.malloc(nbFaces);    // Tableau d'indirection des faces (pour maj PE)
      keep_pg.malloc(nbFaces);       // Flag de conservation des faces 
      new_ph_ids.malloc(nbElements); // Tableau d'indirection des elmts (pour maj PE)

      new_pg_ids = -1;
      new_ph_ids = -1;
      keep_pg    = -1; 
      
      // Boucle sur le nombre d'elements
      for (E_Int i = 0; i < nbElements; i++)
      {
        nbFaces = cnEFp[0];
        if (tagp[i] >= oneEps)
        {
          cn2p[0] = nbFaces; size2 +=1;
          for (E_Int n = 1; n <= nbFaces; n++)
          {
            cn2p[n]             = cnEFp[n];
            keep_pg[cnEFp[n]-1] = +1;
          }
          size2 += nbFaces; cn2p += nbFaces+1; next++;

          // Selection champs en centres 
          for (E_Int k = 1; k <= nfldC; k++)
          { 
            fcenter(ii,k) = fcenter0(i,k);
          }
          new_ph_ids[i] = ii;
          ii++;
        }    
        cnEFp += nbFaces+1;
      }

      E_Int nn = 0; 
      for (E_Int n = 0; n<new_pg_ids.getSize(); n++)
      {
        if (keep_pg[n]>0){ new_pg_ids[n] = nn; nn++; newNumFace++;}
      }
    } 
    else // PE == Py_None - pas de creation de tab d'indirection 
    {
        // Boucle sur le nombre d elements
        for (E_Int i = 0; i < nbElements; i++)
        {
          nbFaces = cnEFp[0];
          if (tagp[i] >= oneEps)
          {
            cn2p[0] = nbFaces; size2 +=1;
            for (E_Int n = 1; n <= nbFaces; n++)
            {
              cn2p[n] = cnEFp[n];
            }
            size2 += nbFaces; cn2p += nbFaces+1; next++;

            // Selection champs en centres 
            for (E_Int k = 1; k <= nfldC; k++)
            { 
              fcenter(ii,k) = fcenter0(i,k);
            }
            ii++;
          }    
          cnEFp += nbFaces+1;
        }
    } // (E == Py_None) 

    cn2.reAlloc(size2);
    
    foutC->reAllocMat(ii,nfldC);
    
    // Cree la nouvelle connectivite complete
    E_Int coutsize = sizeFN+4+size2;
    FldArrayI* cout = new FldArrayI(coutsize);
    E_Int* coutp = cout->begin();
    cn2p = cn2.begin();
    for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = cnpp[i];
    coutp += sizeFN+2;
    coutp[0] = next;
    coutp[1] = size2; coutp += 2;
    for (E_Int i = 0; i < size2; i++) coutp[i] = cn2p[i];

    RELEASESHAREDB(resa, taga, tag, cnpa);
    
    if (PE != Py_None)
    {
      // Check numpy (parentElement)
      FldArrayI* cFE;
      E_Int res = K_NUMPY::getFromNumpyArray(PE, cFE);
      
      if (res == 0)
      {
        RELEASESHAREDN(PE, cFE);
        PyErr_SetString(PyExc_TypeError, "selectCellsBoth: PE numpy is invalid.");
        return NULL;
      }
      
      ngon_t<K_FLD::FldArrayI> ng(*cout); // construction d'un ngon_t à partir d'un FldArrayI
      
      FldArrayI* cFEp_new = new FldArrayI(newNumFace,2);
      FldArrayI& cFE_new  = *cFEp_new ; 

      E_Int* cFEl = cFE_new.begin(1);
      E_Int* cFEr = cFE_new.begin(2);

      E_Int* cFEl_old = cFE->begin(1);
      E_Int* cFEr_old = cFE->begin(2);
      
      E_Int old_ph_1, old_ph_2;

      nbFaces = cnpp[0];

      for (E_Int pgi = 0; pgi < nbFaces; pgi++)
      {
        if (new_pg_ids[pgi]>=0)
        {
          old_ph_1 = cFEl_old[pgi]-1;
          old_ph_2 = cFEr_old[pgi]-1;

          if (old_ph_1 >= 0) // l'elmt gauche existe 
          {
            cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_1]+1;
 
            if (old_ph_2 >= 0) // l'elmt droit existe
            {
              cFEr[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
            }
            else
              cFEr[new_pg_ids[pgi]] = 0;
          } 
          else // l'elmt gauche a disparu - switch droite/gauche
          {
            cFEl[new_pg_ids[pgi]] = new_ph_ids[old_ph_2]+1;
            cFEr[new_pg_ids[pgi]] = 0;
            // reverse
            E_Int s = ng.PGs.stride(pgi);
            E_Int* p = ng.PGs.get_facets_ptr(pgi);
            std::reverse(p, p + s);
          }
        }
      } // boucle pgi      
      
      // export ngon
      ng.export_to_array(*cout);

      // objet Python de sortie
      PyObject* pyPE = K_NUMPY::buildNumpyArray(cFE_new, 1);

      PyList_Append(l,pyPE);
      
      delete cFEp_new;
      RELEASESHAREDN(PE, cFE);
    }
    
    // close
    if (cleanConnectivity == 1 && posx > 0 && posy > 0 && posz > 0)
    {
      K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10, *fout, *cout);
      
    }

    cout->setNGonType(1);
    tpl  = K_ARRAY::buildArray3(*fout,  varString, *cout, eltType, api);
    tplc = K_ARRAY::buildArray3(*foutC, varStringC, *cout, eltType, api);
    delete cout; delete fout; delete foutC; 
  }

  PyList_Append(l, tpl); Py_DECREF(tpl);
  PyList_Append(l, tplc); Py_DECREF(tplc);
  
  RELEASESHAREDB(res, arrayNodes, f, cnp);
  RELEASESHAREDB(resC, arrayCenters, fC, cnpC);
  
  return l;
}
//=============================================================================

