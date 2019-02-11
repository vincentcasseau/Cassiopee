/*    
    Copyright 2013-2019 Onera.

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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//=============================================================================
/*
  Display une zone en solid ou en material.
*/
//=============================================================================
void DataDL::renderGPUUSolidHOZone( UnstructZone *zonep, int zone, int zonet )
{
    // Style
    float color1[ 3 ];
    float color2[ 3 ];

    // Colormap
    float r, g, b;
    void ( *getrgb )( Data * data, double, float *, float *, float * );
    getrgb = _plugins.colorMap->next->f;

    E_Float nz = 1. / _numberOfUnstructZones;
#include "solidStyles.h"

    // Ecrasement si renderTag
    if ( zonep->colorR > -0.5 ) {
        color1[ 0 ] = zonep->colorR;
        color1[ 1 ] = zonep->colorG;
        color1[ 2 ] = zonep->colorB;
    }

#include "selection.h"

    bool is1D = ( ( zonep->eltType == 1 ) | ( zonep->eltType == 10 && zonep->nelts1D > 0 ) );
    if ( is1D == true && ptrState->mode == RENDER )
        glLineWidth( 1. + 5 * zonep->shaderParam1 );
    else if ( is1D == true )
        glLineWidth( 3. );
    else
        glLineWidth( 1. );

    // scale
    E_Float s = MAX( zonep->xmax - zonep->xmin, zonep->ymax - zonep->ymin );
    s = MAX( s, zonep->zmax - zonep->zmin );
    s = 100. / ( s + 1.e-12 );

    // Only for textured rendering, we use vect display =======================
    if ( ptrState->mode == RENDER && zonep->material == 14 && zonep->nfield >= 3 )  // Textured rendering
    {
#ifdef __SHADERS__
        triggerShader( *zonep, zonep->material, s, color1 );
#endif
        int nofield1 = 0;
        int nofield2 = 1;
        int nofield3 = 2;
        int ff;
        double offb = 0.;
        int ret1, ret2, ret3, ret4, i, n1, n2, n3, n4;
#undef PLOT
#include "displayUVectSolidZone.h"
        glLineWidth( 1. );
        return;
    }
    // END Textured rendering ============================================

#ifdef __SHADERS__
    bool must_define_outer_and_inner = true;
    // Activation du shader de tesselation :
    int ishader = 0;
    if ( ( zonep->eltType == UnstructZone::TRI ) and ( zonep->eltSize == 6 ) )
        ishader = 1;  // OK, element de type Tri_6
    // CONTINUER DE MEME POUR LES AUTRES TYPES DE HO
    if ( not this->_shaders.has_tesselation() ) {
        this->_shaders.set_tesselation( ishader );
        must_define_outer_and_inner = true;
    }
    if ( ptrState->mode == RENDER ) {
        if ( zonep->selected == 1 && zonep->active == 1 )
            triggerShader( *zonep, zonep->material, s, color2 );
        else
            triggerShader( *zonep, zonep->material, s, color1 );
    } else {
        if ( zonep->selected == 1 && zonep->active == 1 )
            triggerShader( *zonep, 0, s, color2 );
        else
            triggerShader( *zonep, 0, s, color1 );
    }
    // Pour eviter de tracer le lo order sans faire expres :-)
    if ( must_define_outer_and_inner == true ) {
        unsigned short idShader = this->_shaders.currentShader();
        int t_inner = this->ptrState->inner_tesselation;
        int t_outer = this->ptrState->outer_tesselation;
        this->_shaders[ idShader ]->setUniform( "uInner", (float)t_inner );
        this->_shaders[ idShader ]->setUniform( "uOuter", (float)t_outer );
    }

#endif

    ZoneImplDL *zImpl = static_cast<ZoneImplDL *>( zonep->ptr_impl );
    glCallList( zImpl->_DLsolid );
    glLineWidth( 1. );
}