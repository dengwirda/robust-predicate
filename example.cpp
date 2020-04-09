
// g++ -std=c++17 -pedantic -Wall -O3 -flto -DNDEBUG 
// example.cpp -oexample

#include <iostream>

#include "geompred.hpp"

int main () {

    // Initialise the internal library state. Call at the 
    // start of any program.

    mp_float::exactinit() ;

/*-------------------------------- test prediactes in E^2 */

    double _pa[3] = {           // (d+1) coord. is weight
        +0.0, +0.0, +0.0
        } ;
    double _pb[3] = {
        +1.0, +0.0, +0.1
        } ;
    double _pc[3] = {
        +1.0, +1.0, +0.2
        } ;

    double _qq[3] = {
        +0.5, +0.5, +0.0
        } ;

    double _rr ;

    // Test the orienation of the point QQ wrt. the line
    // PA, PB in E^2. 

    _rr = geompred::orient2d (
        _pa, _pb, _qq
        ) ;

    std::cout << std::showpos;

    std::cout << "orient2d: " << _rr;
    std::cout << std::endl;

    // Test the orienation of the point QQ wrt. the half-
    // space of PA, PB. 

    // This is the unweighted case in E^2.

    _rr = geompred::bisect2d (
        _pa, _pb, _qq
        ) ;

    std::cout << "bisect2d: " << _rr;
    std::cout << std::endl;

    // Test the orienation of the point QQ wrt. the half-
    // space of PA, PB. 

    // This is the "weighted" case in E^2.

    _rr = geompred::bisect2w (
        _pa, _pb, _qq
        ) ;

    std::cout << "bisect2w: " << _rr;
    std::cout << std::endl;

    // Test whether the point QQ is contained within the
    // circumscribing ball associated with the unweighted 
    // simplex PA, PB, PC in E^2.
 
    // This is the unweighted case, so only the geometric
    // coordinates PP[0..1] are used.

    _rr = geompred::inball2d (
        _pa, _pb, _pc, _qq
        ) ;

    std::cout << "inball2d: " << _rr;
    std::cout << std::endl;

    // Test whether the point QQ is contained within the
    // "orthogonal" ball associated with the weighted 
    // simplex PA, PB, PC in E^2.

    // This is the "weighted" case, so the full (x, y, w)
    // coordinates PP[0..2] are used.

    _rr = geompred::inball2w (
        _pa, _pb, _pc, _qq
        ) ;

    std::cout << "inball2w: " << _rr;
    std::cout << std::endl;

/*-------------------------------- test prediactes in E^3 */

    double _PA[4] = {           // (d+1) coord. is weight
        +0.0, +0.0, +0.0, +0.0
        } ;
    double _PB[4] = {
        +1.0, +0.0, +0.0, +0.1
        } ;
    double _PC[4] = {
        +1.0, +1.0, +0.0, +0.2
        } ;
    double _PD[4] = {
        +1.0, +1.0, +1.0, +0.3
        } ;

    double _QQ[4] = {
        +0.5, +0.5, +0.5, +0.0
        } ;

    // Test the orienation of the point QQ wrt. the plane
    // PA, PB, PC in E^3. 

    _rr = geompred::orient3d (
        _PA, _PB, _PC, _QQ
        ) ;

    std::cout << "orient3d: " << _rr;
    std::cout << std::endl;

    // Test the orienation of the point QQ wrt. the half-
    // space of PA, PB. 

    // This is the unweighted case in E^3.

    _rr = geompred::bisect3d (
        _PA, _PB, _QQ
        ) ;

    std::cout << "bisect3d: " << _rr;
    std::cout << std::endl;

    // Test the orienation of the point QQ wrt. the half-
    // space of PA, PB. 

    // This is the "weighted" case in E^3.

    _rr = geompred::bisect3w (
        _PA, _PB, _QQ
        ) ;

    std::cout << "bisect3w: " << _rr;
    std::cout << std::endl;

    // Test whether the point QQ is contained within the
    // circumscribing ball associated with the unweighted 
    // simplex PA, PB, PC, PD in E^3.
 
    // This is the unweighted case, so only the geometric
    // coordinates PP[0..2] are used.

    _rr = geompred::inball3d (
        _PA, _PB, _PC, _PD, _QQ
        ) ;

    std::cout << "inball3d: " << _rr;
    std::cout << std::endl;

    // Test whether the point QQ is contained within the
    // "orthogonal" ball associated with the weighted 
    // simplex PA, PB, PC, PD in E^3.

    // This is the "weighted" case, so the full (x,y,z,w)
    // coordinates PP[0..3] are used.

    _rr = geompred::inball3w (
        _PA, _PB, _PC, _PD, _QQ
        ) ;

    std::cout << "inball3w: " << _rr;
    std::cout << std::endl;

    return 0 ;
}



