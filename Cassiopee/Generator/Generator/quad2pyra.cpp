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

# include "generator.h"

//=============================================================================
// Create a pyramid for each input quad
//=============================================================================
PyObject* K_GENERATOR::quad2Pyra(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float hratio = 0.5;

  if (!PYPARSETUPLE_(args, O_ R_, &array, &hratio)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "quad2Pyra: input array must be QUAD, not STRUCT.");
    RELEASESHAREDS(array, f);
    return NULL;
  }
  else if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "quad2Pyra: unknown type of array.");
    return NULL;
  }
  else if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "quad2Pyra: input array must be QUAD, not NGON.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }

  // Check that only QUADs are present in the input array
  E_Int nc = cn->getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    if (K_STRING::cmp(eltTypes[ic], "QUAD") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                     "quad2Pyra: input array must be QUAD only.");
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "quad2Pyra: coords must be present in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float *xp = f->begin(posx), *yp = f->begin(posy), *zp = f->begin(posz);

  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();

  // Compute total number of elements across all connectivities, ntotElts
  std::vector<E_Int> nepc(nc);
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;  // cumulative number of elts per conn.
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic] = nelts;
    cumnepc[ic+1] = cumnepc[ic] + nelts;
  }
  E_Int ntotElts = cumnepc[nc];

  // Build new connectivity
  E_Int npts2 = npts + ntotElts;  // one new vertex for each QUAD
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2,
                                       nepc, "PYRA", false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  // Coordinates of the tip of the new PYRA
  std::vector<E_Float> x5(ntotElts), y5(ntotElts), z5(ntotElts);
  // Unit normals of each QUAD (pointing outside the PYRA)
  std::vector<E_Float> nx(ntotElts), ny(ntotElts), nz(ntotElts);

  #pragma omp parallel
  {
    E_Int nelts, ioffset, ind1, ind2, ind3, ind4;
    E_Float l14x, l14y, l14z, l12x, l12y, l12z, mag;
    E_Float l13x, l13y, l13z, diagLength, height;

    // Connectivities and coords of the tip of the PYRA
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      FldArrayI& cm2 = *(cn2->getConnect(ic));
      nelts = cm.getSize();
      #pragma omp for schedule(static)
      for (E_Int i = 0; i < nelts; i++)
      {
        ioffset = cumnepc[ic] + i;
        ind1 = cm(i, 1); ind2 = cm(i, 2); ind3 = cm(i, 3); ind4 = cm(i, 4);
        cm2(i, 1) = ind1; cm2(i, 2) = ind2;
        cm2(i, 3) = ind3; cm2(i, 4) = ind4;
        cm2(i, 5) = npts + ioffset + 1;
        ind1 -= 1; ind2 -= 1; ind3 -= 1; ind4 -= 1;

        // Barycenter of the element
        x5[ioffset] = K_CONST::ONE_FOURTH * (xp[ind1] + xp[ind2] + xp[ind3] + xp[ind4]);
        y5[ioffset] = K_CONST::ONE_FOURTH * (yp[ind1] + yp[ind2] + yp[ind3] + yp[ind4]);
        z5[ioffset] = K_CONST::ONE_FOURTH * (zp[ind1] + zp[ind2] + zp[ind3] + zp[ind4]);

        // Compute 'base' unit normal pointing outside the PYRA
        l12x = xp[ind2] - xp[ind1]; l14x = xp[ind4] - xp[ind1];
        l12y = yp[ind2] - yp[ind1]; l14y = yp[ind4] - yp[ind1];
        l12z = zp[ind2] - zp[ind1]; l14z = zp[ind4] - zp[ind1];
        K_MATH::cross(
          l12x, l12y, l12z,
          l14x, l14y, l14z,
          nx[ioffset], ny[ioffset], nz[ioffset]
        );
        mag = std::sqrt(
            nx[ioffset] * nx[ioffset]
          + ny[ioffset] * ny[ioffset]
          + nz[ioffset] * nz[ioffset]
        );
        if (!K_FUNC::fEqualZero(mag))
        {
          nx[ioffset] /= mag; ny[ioffset] /= mag; nz[ioffset] /= mag;
        }

        // Compute the diagonal length of the 'base'
        // l13x = xp[ind3] - xp[ind1];
        // l13y = yp[ind3] - yp[ind1];
        // l13z = zp[ind3] - zp[ind1];
        // diagLength = std::sqrt(l13x * l13x + l13y * l13y + l13z * l13z);
        l13x = xp[ind4] - xp[ind2];
        l13y = yp[ind4] - yp[ind2];
        l13z = zp[ind4] - zp[ind2];
        diagLength = std::sqrt(l13x * l13x + l13y * l13y + l13z * l13z);
        // diagLength *= K_CONST::ONE_HALF;

        // Move barycenter along the normal direction
        height = K_CONST::ONE_HALF * diagLength * hratio;
        x5[ioffset] += height * nx[ioffset];
        y5[ioffset] += height * ny[ioffset];
        z5[ioffset] += height * nz[ioffset];
      }
    }

    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];

      if (n == posx)
      {
        #pragma omp for nowait
        for (E_Int i = npts; i < npts2; i++) f2p[i] = x5[i-npts];
      }
      else if (n == posy)
      {
        #pragma omp for nowait
        for (E_Int i = npts; i < npts2; i++) f2p[i] = y5[i-npts];
      }
      else if (n == posz)
      {
        #pragma omp for nowait
        for (E_Int i = npts; i < npts2; i++) f2p[i] = z5[i-npts];
      }
      else
      {
        #pragma omp for nowait
        for (E_Int i = npts; i < npts2; i++) f2p[i] = 0.;
      }
    }
  }

  RELEASESHAREDU(tpl, f2, cn2);
  RELEASESHAREDU(array, f, cn);
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return tpl;
}
