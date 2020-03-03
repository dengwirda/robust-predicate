
    /*
    --------------------------------------------------------
     * MPFLOAT: multi-precision floating-point arithmetic.
    --------------------------------------------------------
     *
     * These are the low-level multi-precision kernels ---
     * computing elementary operations on "expansions" of
     * floating-point numbers such that rounding error is
     * eliminated. See Shewchuk for more detail:
     *
     * J. R. Shewchuk (1997): Adaptive Precision Floating-
     * Point Arithmetic & Fast Robust Geometric Predicates
     * Discrete & Computational Geometry, 18, pp. 305-363.
     *
     * This header is adapted from Shewchuk's original C89
     * source (predicates.c).
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
     * Last updated: 01 March, 2020
     *
     * Copyright 2020--
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

#   pragma once

#   ifndef __MP_BASIC__
#   define __MP_BASIC__

    namespace mp_float {

#   define REAL_TYPE mp_float::real_type
#   define INDX_TYPE mp_float::indx_type

    /*
    --------------------------------------------------------
     * multi-precision initialisation, a'la shewchuk
    --------------------------------------------------------
     */

    REAL_TYPE _splitter;
    REAL_TYPE _epsilon ;

    __normal_call void exactinit (
        )
    {
    /*-------------- find machine eps, etc, a'la shewchuk */
        INDX_TYPE _alternate = +1 ;
        REAL_TYPE _lastcheck ;
        REAL_TYPE _halve = +0.5;
        REAL_TYPE _check = +1.0;

    /*-------------- find eps: bisect until 1. + eps ~ 1. */

        _epsilon = _splitter = +1.00 ;

        do {
            _lastcheck = _check;
            _epsilon  *= _halve;

            if (_alternate)
            _splitter *= +2.00 ;

            _alternate = !_alternate ;

            _check = 1.00 + _epsilon ;
        }
        while (_check != +1.00 &&
               _check != _lastcheck) ;

        _splitter += 1.00 ;
    }

    /*
    --------------------------------------------------------
     * multi-precision "add" routines, a'la shewchuk
    --------------------------------------------------------
     */

    __inline_call void one_one_add_fast (
        REAL_TYPE  _aa, REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _bvirt;
        _x1 = _aa + _bb;
        _bvirt = _x1 - _aa;
        _x0 = _bb - _bvirt;
    }

    __inline_call void one_one_add_full (
        REAL_TYPE  _aa, REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _bvirt, _avirt;
        _x1 = _aa + _bb;
        _bvirt = _x1 - _aa;
        _avirt = _x1 - _bvirt;

        REAL_TYPE _bround, _around;
        _bround = _bb - _bvirt;
        _around = _aa - _avirt;
        _x0 = _around + _bround;
    }

    __inline_call void two_one_add_full (
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _bb,
        REAL_TYPE &_x2, REAL_TYPE &_x1,
        REAL_TYPE &_x0
        )
    {
        REAL_TYPE _tt;
        one_one_add_full(_a0, _bb, _tt, _x0
            ) ;
        one_one_add_full(_a1, _tt, _x2, _x1
            ) ;
    }

    __inline_call void two_one_add_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t0, _t1 ;
        one_one_add_full(_a1, _bb, _t1, _t0
            ) ;

        _t0 = _t0 + _a0 ;

        one_one_add_fast(_t1, _t0, _x1, _x0
            ) ;
    }

    __inline_call void two_two_add_full (
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _b1, REAL_TYPE  _b0,
        REAL_TYPE &_x3, REAL_TYPE &_x2,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t1, _t0 ;
        two_one_add_full(_a1, _a0, _b0, _t1,
            _t0, _x0
            ) ;
        two_one_add_full(_t1, _t0, _b1, _x3,
            _x2, _x1
            ) ;
    }

    __inline_call void two_two_add_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _b1, REAL_TYPE  _b0,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t1, _t0 ;
        REAL_TYPE _s1, _s0 ;
        REAL_TYPE _w1, _w0 ;
        one_one_add_full(_a1, _b1, _s1, _s0
            ) ;
        one_one_add_full(_a0, _b0, _t1, _t0
            ) ;

        _s0 = _s0 + _t1 ;

        one_one_add_fast(_s1, _s0, _w1, _w0
            ) ;

        _w0 = _w0 + _t0 ;

        one_one_add_fast(_w1, _w0, _x1, _x0
            ) ;
    }

    /*
    --------------------------------------------------------
     * multi-precision "sub" routines, a'la shewchuk
    --------------------------------------------------------
     */

    __inline_call void one_one_sub_fast (
        REAL_TYPE  _aa, REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _bvirt;
        _x1 = _aa - _bb;
        _bvirt = _aa - _x1;
        _x0 = _bvirt - _bb;
    }

    __inline_call void one_one_sub_full (
        REAL_TYPE  _aa, REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _bvirt, _avirt;
        _x1 = _aa - _bb;
        _bvirt = _aa - _x1;
        _avirt = _x1 + _bvirt;

        REAL_TYPE _bround, _around;
        _bround = _bvirt - _bb;
        _around = _aa - _avirt;
        _x0 = _around + _bround;
    }

    __inline_call void two_one_sub_full (
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _bb,
        REAL_TYPE &_x2, REAL_TYPE &_x1,
        REAL_TYPE &_x0
        )
    {
        REAL_TYPE _tt;
        one_one_sub_full(_a0, _bb, _tt, _x0
            ) ;
        one_one_add_full(_a1, _tt, _x2, _x1
            ) ;
    }

    __inline_call void two_one_sub_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t0, _t1 ;
        one_one_sub_full(_a1, _bb, _t1, _t0
            ) ;

        _t0 = _t0 + _a0 ;

        one_one_add_fast(_t1, _t0, _x1, _x0
            ) ;
    }

    __inline_call void two_two_sub_full (
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _b1, REAL_TYPE  _b0,
        REAL_TYPE &_x3, REAL_TYPE &_x2,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t1, _t0 ;
        two_one_sub_full(_a1, _a0, _b0, _t1,
            _t0, _x0
            ) ;
        two_one_sub_full(_t1, _t0, _b1, _x3,
            _x2, _x1
            ) ;
    }

    __inline_call void two_two_sub_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _b1, REAL_TYPE  _b0,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _s0, _s1 ;
        REAL_TYPE _t0, _t1 ;
        REAL_TYPE _w0, _w1 ;
        one_one_sub_full(_a1, _b1, _s1, _s0
            ) ;
        one_one_sub_full(_a0, _b0, _t1, _t0
            ) ;

        _s0 = _s0 + _t1 ;

        one_one_add_fast(_s1, _s0, _w1, _w0
            ) ;

        _w0 = _w0 + _t0 ;

        one_one_add_fast(_w1, _w0, _x1, _x0
            ) ;
    }

    /*
    --------------------------------------------------------
     * multi-precision "mul" routines, a'la shewchuk
    --------------------------------------------------------
     */

    __inline_call void one_split (
        REAL_TYPE  _aa,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _cc, _ab ;
        _cc = _aa * _splitter;
        _ab = _cc - _aa;
        _x1 = _cc - _ab;
        _x0 = _aa - _x1;
    }

    __inline_call void one_one_mul_full (
        REAL_TYPE  _aa, REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _ah, _al, _bh, _bl;
        _x1 = _aa * _bb;
        one_split (_aa, _ah, _al);
        one_split (_bb, _bh, _bl);

        REAL_TYPE _err1, _err2, _err3;
        _err1 = _x1 - (_ah * _bh);
        _err2 = _err1 - (_al * _bh);
        _err3 = _err2 - (_ah * _bl);
        _x0 = (_al * _bl) - _err3;
    }

    __inline_call void one_one_mul_full (
        REAL_TYPE  _aa,
        REAL_TYPE  _bb, REAL_TYPE  _bh,
        REAL_TYPE  _bl,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _ah, _al;
        _x1 = _aa * _bb;
        one_split (_aa, _ah, _al);

        REAL_TYPE _err1, _err2, _err3;
        _err1 = _x1 - (_ah * _bh);
        _err2 = _err1 - (_al * _bh);
        _err3 = _err2 - (_ah * _bl);
        _x0 = (_al * _bl) - _err3;
    }

    __inline_call void one_one_mul_full (
        REAL_TYPE  _aa, REAL_TYPE  _ah,
        REAL_TYPE  _al,
        REAL_TYPE  _bb, REAL_TYPE  _bh,
        REAL_TYPE  _bl,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        _x1 = _aa * _bb;

        REAL_TYPE _err1, _err2, _err3;
        _err1 = _x1 - (_ah * _bh);
        _err2 = _err1 - (_al * _bh);
        _err3 = _err2 - (_ah * _bl);
        _x0 = (_al * _bl) - _err3;
    }

    __inline_call void one_square_full (
        REAL_TYPE  _aa,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _ah, _al;
        _x1 = _aa * _aa;
        one_split (_aa, _ah, _al);

        REAL_TYPE _err1, _err3;
        _err1 = _x1 - (_ah * _ah);
        _err3 = _err1 - ((_ah + _ah) * _al);
        _x0 = (_al * _al) - _err3;
    }

    __inline_call void one_square_full (
        REAL_TYPE  _aa, REAL_TYPE  _ah,
        REAL_TYPE  _al,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        _x1 = _aa * _aa;

        REAL_TYPE _err1, _err3;
        _err1 = _x1 - (_ah * _ah);
        _err3 = _err1 - ((_ah + _ah) * _al);
        _x0 = (_al * _al) - _err3;
    }

    __inline_call void two_one_mul_full (
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _bb,
        REAL_TYPE &_x3, REAL_TYPE &_x2,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _bh, _bl;
        REAL_TYPE _t0, _t1, _t2, _t3 ;
        one_split(_bb, _bh, _bl) ;

        one_one_mul_full(_a0, _bb, _bh, _bl,
            _t2, _x0
            ) ;
        one_one_mul_full(_a1, _bb, _bh, _bl,
            _t1, _t0
            ) ;

        one_one_add_full(_t2, _t0, _t3, _x1
            ) ;
        one_one_add_fast(_t1, _t3, _x3, _x2
            ) ;
    }

    __inline_call void two_one_mul_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _bb,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t0, _t1, _ss ;
        one_one_mul_full(_a1, _bb, _t1, _t0
            ) ;

        _ss = _a0 * _bb ;
        _t0 = _t0 + _ss ;

        one_one_add_fast(_t1, _t0, _x1, _x0
            ) ;
    }

    __inline_call void two_two_mul_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _b1, REAL_TYPE  _b0,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t0, _t1, _ss ;
        one_one_mul_full(_a1, _b1, _t1, _t0
            ) ;

        _ss = _a1 * _b0 +
              _a0 * _b1 ;
        _t0 = _t0 + _ss ;

        one_one_add_fast(_t1, _t0, _x1, _x0
            ) ;
    }

    __inline_call void two_two_div_fast (   // dd_flt
        REAL_TYPE  _a1, REAL_TYPE  _a0,
        REAL_TYPE  _b1, REAL_TYPE  _b0,
        REAL_TYPE &_x1, REAL_TYPE &_x0
        )
    {
        REAL_TYPE _t0, _t1, _ee;
        _t1 = _a1 / _b1 ;

        REAL_TYPE _r0, _r1 ;
        REAL_TYPE _w0, _w1 ;
        two_one_mul_fast(_b1, _b0, _t1,
            _r1, _r0                    // rr = bb * t1
            ) ;
        two_two_sub_fast(_a1, _a0, _r1, _r0,
            _w1, _w0                    // ww = aa - rr
            ) ;

        _t0 = _w1 / _b1 ;

        REAL_TYPE _u0, _u1 ;
        two_one_mul_fast(_b1, _b0, _t0,
            _r1, _r0                    // rr = bb * t0
            ) ;
        two_two_sub_fast(_w1, _w0, _r1, _r0,
            _u1, _u0                    // uu = ww - rr
            ) ;

        _ee = _u1 / _b1 ;

        REAL_TYPE _q0, _q1 ;
        one_one_add_fast(_t1, _t0, _q1, _q0
            ) ;
        two_one_add_fast(_q1, _q0, _ee,
            _t1, _t0
            ) ;

        _x1 = _t1 ; _x0 = _t0 ;
    }

#   undef REAL_TYPE
#   undef INDX_TYPE


    }

#   endif//__MP_BASIC__



