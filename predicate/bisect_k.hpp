
    /*
    --------------------------------------------------------
     * PREDICATE-k: robust geometric predicates in E^k.
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
     * Last updated: 07 April, 2020
     *
     * Copyright 2020--
     * Darren Engwirda
     * de2363@columbia.edu
     * https://github.com/dengwirda/
     *
    --------------------------------------------------------
     */

    // from predicate_k.hpp...


    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using 
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the unweighted "bisect" predicate in E^2.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect2d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect2d predicate, "exact" version */
        mp::expansion< 2 > _ac_xx_, _ac_yy_, 
                           _bc_xx_, _bc_yy_;
        mp::expansion< 16> _ac_sq_, _bc_sq_;
        mp::expansion< 32> _absum_;
        
        _FT = (REAL_TYPE) +0.0E+00;

    /*----------------------------------- compute: d(p,q) */
        _ac_xx_.from_sub(_pc[0], _pa[0]);
        _ac_yy_.from_sub(_pc[1], _pa[1]);

        _bc_xx_.from_sub(_pc[0], _pb[0]);
        _bc_yy_.from_sub(_pc[1], _pb[1]);
    
        mp::expansion_dot(_ac_xx_, _ac_xx_,
                          _ac_yy_, _ac_yy_,
                          _ac_sq_) ;

        mp::expansion_dot(_bc_xx_, _bc_xx_,
                          _bc_yy_, _bc_yy_,
                          _bc_sq_) ;

    /*----------------------------------- d(a,c) - d(b,c) */
        mp::expansion_sub(
            _ac_sq_, _bc_sq_, _absum_) ;

    /*----------------------------------- return signed d */
        return _absum_[_absum_._xlen - 1] ;
    }

    __normal_call REAL_TYPE bisect2d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect2d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  5. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy ;
        REAL_TYPE _bcx, _bcy ;
        REAL_TYPE _acsqr, _bcsqr ;

        REAL_TYPE _ACSQR, _BCSQR ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;

        _acsqr = _acx * _acx +
                 _acy * _acy ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy ;

        _ACSQR = std::abs(_acsqr);
        _BCSQR = std::abs(_bcsqr);

        _FT  = _ACSQR + _BCSQR ;          // roundoff tol
        _FT *= _ER ;

        return _acsqr - _bcsqr ;          // d_ab - d_bc
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using 
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the weighted "bisect" predicate in E^2.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect2w_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect2w predicate, "exact" version */
        mp::expansion< 2 > _ac_xx_, _ac_yy_, 
                           _bc_xx_, _bc_yy_;
        mp::expansion< 16> _ac_sq_, _bc_sq_;
        mp::expansion< 17> _a_sum_, _b_sum_;
        mp::expansion< 34> _absum_;
        
        _FT = (REAL_TYPE) +0.0E+00;

        mp::expansion< 1 > _pa_ww_(_pa[ 2]);
        mp::expansion< 1 > _pb_ww_(_pb[ 2]);

    /*----------------------------------- compute: d(p,q) */
        _ac_xx_.from_sub(_pc[0], _pa[0]);
        _ac_yy_.from_sub(_pc[1], _pa[1]);

        _bc_xx_.from_sub(_pc[0], _pb[0]);
        _bc_yy_.from_sub(_pc[1], _pb[1]);
    
        mp::expansion_dot(_ac_xx_, _ac_xx_,
                          _ac_yy_, _ac_yy_,
                          _ac_sq_) ;

        mp::expansion_dot(_bc_xx_, _bc_xx_,
                          _bc_yy_, _bc_yy_,
                          _bc_sq_) ;

        mp::expansion_sub(
            _ac_sq_, _pa_ww_, _a_sum_) ;

        mp::expansion_sub(
            _bc_sq_, _pb_ww_, _b_sum_) ;

    /*----------------------------------- d(a,c) - d(b,c) */
        mp::expansion_sub(
            _a_sum_, _b_sum_, _absum_) ;

    /*----------------------------------- return signed d */
        return _absum_[_absum_._xlen - 1] ;
    }

    __normal_call REAL_TYPE bisect2w_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect2w predicate, "float" version */
        REAL_TYPE static const _ER =
        +  6. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy ;
        REAL_TYPE _bcx, _bcy ;
        REAL_TYPE _acsqr, _bcsqr ;
        REAL_TYPE _a_sum, _b_sum ;

        REAL_TYPE _A_SUM, _B_SUM ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;

        _acsqr = _acx * _acx +
                 _acy * _acy ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy ;

        _a_sum = _acsqr - _pa[2] ;
        _b_sum = _bcsqr - _pb[2] ;

        _A_SUM = std::abs(_acsqr)
               + std::abs(_pa[2]);
        _B_SUM = std::abs(_bcsqr)
               + std::abs(_pb[2]);

        _FT  = _A_SUM + _B_SUM ;          // roundoff tol
        _FT *= _ER ;

        return _a_sum - _b_sum ;          // d_ab - d_bc
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using 
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the unweighted "bisect" predicate in E^3.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect3d_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect3d predicate, "exact" version */
        mp::expansion< 2 > _ac_xx_, _ac_yy_, 
                           _ac_zz_,
                           _bc_xx_, _bc_yy_,
                           _bc_zz_;
        mp::expansion< 24> _ac_sq_, _bc_sq_;
        mp::expansion< 48> _absum_;
        
        _FT = (REAL_TYPE) +0.0E+00;

    /*----------------------------------- compute: d(p,q) */
        _ac_xx_.from_sub(_pc[0], _pa[0]);
        _ac_yy_.from_sub(_pc[1], _pa[1]);
        _ac_zz_.from_sub(_pc[2], _pa[2]);

        _bc_xx_.from_sub(_pc[0], _pb[0]);
        _bc_yy_.from_sub(_pc[1], _pb[1]);
        _bc_zz_.from_sub(_pc[2], _pb[2]);
    
        mp::expansion_dot(_ac_xx_, _ac_xx_,
                          _ac_yy_, _ac_yy_,
                          _ac_zz_, _ac_zz_,
                          _ac_sq_) ;

        mp::expansion_dot(_bc_xx_, _bc_xx_,
                          _bc_yy_, _bc_yy_,
                          _bc_zz_, _bc_zz_,
                          _bc_sq_) ;

    /*----------------------------------- d(a,c) - d(b,c) */
        mp::expansion_sub(
            _ac_sq_, _bc_sq_, _absum_) ;

    /*----------------------------------- return signed d */
        return _absum_[_absum_._xlen - 1] ;
    }

    __normal_call REAL_TYPE bisect3d_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect3d predicate, "float" version */
        REAL_TYPE static const _ER =
        +  6. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy, _acz ;
        REAL_TYPE _bcx, _bcy, _bcz ;
        REAL_TYPE _acsqr, _bcsqr ;

        REAL_TYPE _ACSQR, _BCSQR ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;
        _acz = _pa [2] - _pc [2] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;
        _bcz = _pb [2] - _pc [2] ;

        _acsqr = _acx * _acx +
                 _acy * _acy +
                 _acz * _acz ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy +
                 _bcz * _bcz ;

        _ACSQR = std::abs(_acsqr);
        _BCSQR = std::abs(_bcsqr);

        _FT  = _ACSQR + _BCSQR ;          // roundoff tol
        _FT *= _ER ;

        return _acsqr - _bcsqr ;          // d_ab - d_bc
    }

    /*
    --------------------------------------------------------
     *
     * Compute an exact orientation wrt. bisector using 
     * multi-precision expansions, a'la shewchuk
     *
     *   |c-a|**2 - wa = |c-b|**2 - wb
     *
     * This is the weighted "bisect" predicate in E^3.
     *
    --------------------------------------------------------
     */

    __normal_call REAL_TYPE bisect3w_e (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect3w predicate, "exact" version */
        mp::expansion< 2 > _ac_xx_, _ac_yy_, 
                           _ac_zz_,
                           _bc_xx_, _bc_yy_,
                           _bc_zz_;
        mp::expansion< 24> _ac_sq_, _bc_sq_;
        mp::expansion< 25> _a_sum_, _b_sum_;
        mp::expansion< 50> _absum_;
        
        _FT = (REAL_TYPE) +0.0E+00;

        mp::expansion< 1 > _pa_ww_(_pa[ 3]);
        mp::expansion< 1 > _pb_ww_(_pb[ 3]);

    /*----------------------------------- compute: d(p,q) */
        _ac_xx_.from_sub(_pc[0], _pa[0]);
        _ac_yy_.from_sub(_pc[1], _pa[1]);
        _ac_zz_.from_sub(_pc[2], _pa[2]);

        _bc_xx_.from_sub(_pc[0], _pb[0]);
        _bc_yy_.from_sub(_pc[1], _pb[1]);
        _bc_zz_.from_sub(_pc[2], _pb[2]);
    
        mp::expansion_dot(_ac_xx_, _ac_xx_,
                          _ac_yy_, _ac_yy_,
                          _ac_zz_, _ac_zz_,
                          _ac_sq_) ;

        mp::expansion_dot(_bc_xx_, _bc_xx_,
                          _bc_yy_, _bc_yy_,
                          _bc_zz_, _bc_zz_,
                          _bc_sq_) ;

        mp::expansion_sub(
            _ac_sq_, _pa_ww_, _a_sum_) ;

        mp::expansion_sub(
            _bc_sq_, _pb_ww_, _b_sum_) ;

    /*----------------------------------- d(a,c) - d(b,c) */
        mp::expansion_sub(
            _a_sum_, _b_sum_, _absum_) ;

    /*----------------------------------- return signed d */
        return _absum_[_absum_._xlen - 1] ;
    }

    __normal_call REAL_TYPE bisect3w_f (
      __const_ptr(REAL_TYPE) _pa ,
      __const_ptr(REAL_TYPE) _pb ,
      __const_ptr(REAL_TYPE) _pc ,
        REAL_TYPE &_FT
        )
    {
    /*--------------- bisect3w predicate, "float" version */
        REAL_TYPE static const _ER =
        +  7. * std::pow(mp::_epsilon, 1) ;

        REAL_TYPE _acx, _acy, _acz ;
        REAL_TYPE _bcx, _bcy, _bcz ;
        REAL_TYPE _acsqr, _bcsqr ;
        REAL_TYPE _a_sum, _b_sum ;

        REAL_TYPE _A_SUM, _B_SUM ;

        _acx = _pa [0] - _pc [0] ;        // coord. diff.
        _acy = _pa [1] - _pc [1] ;
        _acz = _pa [2] - _pc [2] ;

        _bcx = _pb [0] - _pc [0] ;
        _bcy = _pb [1] - _pc [1] ;
        _bcz = _pb [2] - _pc [2] ;

        _acsqr = _acx * _acx +
                 _acy * _acy +
                 _acz * _acz ;
        _bcsqr = _bcx * _bcx +
                 _bcy * _bcy +
                 _bcz * _bcz ;

        _a_sum = _acsqr - _pa[3] ;
        _b_sum = _bcsqr - _pb[3] ;

        _A_SUM = std::abs(_acsqr)
               + std::abs(_pa[3]);
        _B_SUM = std::abs(_bcsqr)
               + std::abs(_pb[3]);

        _FT  = _A_SUM + _B_SUM ;          // roundoff tol
        _FT *= _ER ;

        return _a_sum - _b_sum ;          // d_ab - d_bc
    }



