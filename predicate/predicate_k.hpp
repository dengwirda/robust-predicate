
    /*
    --------------------------------------------------------
     * PREDICATE-k: robust geometric predicates in E^k.
    --------------------------------------------------------
     *
     * Compute "robust" geometric predicates using filtered
     * floating-point + multi-precision expansions.
     *
     * The sign-correctness of each predicate is guaranteed
     * --- using exact arithmetic where necessary to
     * eliminate floating-point round-off. See Shewchuk for
     * additional detail
     *
     * J. R. Shewchuk (1997), Adaptive Precision Floating-
     * Point Arithmetic & Fast Robust Geometric Predicates
     * Discrete & Computational Geometry, 18, pp. 305-363.
     *
     * A translational version of BFS semi-static filtering
     * is employed, adapted from, e.g.
     *
     * C. Burnikel, S. Funke, and M. Seel (2001), Exact 
     * geometric computation using cascading.
     * IJCGA (Special issue) 11 (3), pp. 245–266.
     *
     * O. Devillers and S. Pion (2002), Efficient Exact 
     * Geometric Predicates for Delaunay Triangulations.
     * RR-4351, INRIA. inria-00072237
     *
    --------------------------------------------------------
     *
     * This program may be freely redistributed under the
     * condition that the copyright notices (including this
     * entire header) are not removed, and no compensation
     * is received through use of the software.  Private,
     * research, and institutional use is free.  You may
     * distribute modified versions of this code UNDER THE
     * CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE
     * TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE
     * ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE
     * MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR
     * NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution
     * of this code as part of a commercial system is
     * permissible ONLY BY DIRECT ARRANGEMENT WITH THE
     * AUTHOR.  (If you are not directly supplying this
     * code to a customer, and you are instead telling them
     * how they can obtain it for free, then you are not
     * required to make any arrangement with me.)
     *
     * Disclaimer:  Neither I nor: Columbia University, The
     * Massachusetts Institute of Technology, The
     * University of Sydney, nor The National Aeronautics
     * and Space Administration warrant this code in any
     * way whatsoever.  This code is provided "as-is" to be
     * used at your own risk.
     *
    --------------------------------------------------------
     *
     * Last updated: 09 April, 2020
     *
     * Copyright 2020--
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __PREDICATE_K__
#   define __PREDICATE_K__

    namespace geompred {

#   define REAL_TYPE mp_float::real_type
#   define INDX_TYPE mp_float::indx_type

    namespace mp=mp_float;

#   include "orient_k.hpp"
#   include "bisect_k.hpp"
//  include "linear_k.hpp"
#   include "inball_k.hpp"

    __inline_call REAL_TYPE orient2d (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc
        )
    {
    /*------------ orient2d predicate, "filtered" version */
        REAL_TYPE _FT, _rr;

        _rr = orient2d_f(               // "float" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;

        _rr = orient2d_e(               // "exact" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
    }

    __inline_call REAL_TYPE orient3d (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd
        )
    {
    /*------------ orient3d predicate, "filtered" version */
        REAL_TYPE _FT, _rr;

        _rr = orient3d_f(               // "float" kernel
            _pa, _pb, _pc, _pd, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;

        _rr = orient3d_e(               // "exact" kernel
            _pa, _pb, _pc, _pd, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
    }

    __inline_call REAL_TYPE bisect2d (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc        
        )
    {
    /*------------ bisect2d predicate, "filtered" version */
        REAL_TYPE _FT, _rr;

        _rr = bisect2d_f(               // "float" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;
        
        _rr = bisect2d_e(               // "exact" kernel
            _pa, _pb, _pc, _FT
            ) ;
        
        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
    }

    __inline_call REAL_TYPE bisect2w (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc        
        )
    {
    /*------------ bisect2w predicate, "filtered" version */
        if (_pa [ 2] == _pb [ 2] )
        {                   // equal weights, do bisect2d
        return bisect2d(_pa, _pb, _pc) ;
        }
        else
        {
        REAL_TYPE _FT, _rr; // given weights, full kernel

        _rr = bisect2w_f(               // "float" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;
        
        _rr = bisect2w_e(               // "exact" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
        }
    }

    __inline_call REAL_TYPE bisect3d (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc    
        )
    {
    /*------------ bisect3d predicate, "filtered" version */
        REAL_TYPE _FT, _rr;

        _rr = bisect3d_f(               // "float" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;
        
        _rr = bisect3d_e(               // "exact" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
    }

    __inline_call REAL_TYPE bisect3w (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc
        )
    {
    /*------------ bisect2w predicate, "filtered" version */
        if (_pa [ 3] == _pb [ 3] )
        {                   // equal weights, do bisect3d
        return bisect3d(_pa, _pb, _pc) ;
        }
        else
        {
        REAL_TYPE _FT, _rr; // given weights, full kernel

        _rr = bisect3w_f(               // "float" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;
        
        _rr = bisect3w_e(               // "exact" kernel
            _pa, _pb, _pc, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
        }
    }

    __inline_call REAL_TYPE inball2d (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd
        )
    {
    /*------------ inball2d predicate, "filtered" version */
        REAL_TYPE _FT, _rr;

        _rr = inball2d_f(               // "float" kernel
            _pa, _pb, _pc, _pd, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;

        _rr = inball2d_e(               // "exact" kernel
            _pa, _pb, _pc, _pd, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
    }

    __inline_call REAL_TYPE inball2w (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd
        )
    {
    /*------------ inball2w predicate, "filtered" version */
        if (_pa [ 2] == _pb [ 2] &&
            _pb [ 2] == _pc [ 2] &&
            _pc [ 2] == _pd [ 2] )
        {
        return inball2d (   // equal weights, do inball2d
            _pa, _pb, _pc, _pd
            ) ;
        }
        else
        {
        REAL_TYPE _FT, _rr; // given weights, full kernel

        _rr = inball2w_f(               // "float" kernel
            _pa, _pb, _pc, _pd, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;

        _rr = inball2w_e(               // "exact" kernel
            _pa, _pb, _pc, _pd, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
        }
    }

    __inline_call REAL_TYPE inball3d (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
      __const_ptr(REAL_TYPE) _pe
        )
    {
    /*------------ inball3d predicate, "filtered" version */
        REAL_TYPE _FT, _rr;

        _rr = inball3d_f(               // "float" kernel
            _pa, _pb, _pc, _pd, _pe, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;

        _rr = inball3d_e(               // "exact" kernel
            _pa, _pb, _pc, _pd, _pe, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
    }

    __inline_call REAL_TYPE inball3w (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
      __const_ptr(REAL_TYPE) _pd ,
      __const_ptr(REAL_TYPE) _pe
        )
    {
    /*------------ inball3w predicate, "filtered" version */
        if (_pa [ 3] == _pb [ 3] &&
            _pb [ 3] == _pc [ 3] &&
            _pc [ 3] == _pd [ 3] &&
            _pd [ 3] == _pe [ 3] )
        {
        return inball3d (   // equal weights, do inball3d
            _pa, _pb, _pc, _pd, _pe
            ) ;
        }
        else
        {
        REAL_TYPE _FT, _rr; // given weights, full kernel

        _rr = inball3w_f(               // "float" kernal
            _pa, _pb, _pc, _pd, _pe, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            if (std::isnormal(_rr))
            return _rr ;

        _rr = inball3w_e(               // "exact" kernel
            _pa, _pb, _pc, _pd, _pe, _FT
            ) ;

        if (_rr > _FT || _rr < -_FT)
            return _rr ;

        return (REAL_TYPE) +0.0E+00;
        }
    }

#   undef REAL_TYPE
#   undef INDX_TYPE


    }

#   endif//__PREDICATE_K__



