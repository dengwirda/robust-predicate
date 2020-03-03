
// g++ -std=c++17 -pedantic -Wall -03 
// -flto -DNDEBUG example.cpp -oexample

#include "geompred.hpp"

int main () {

    double _pa[3] = {
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

    _rr = geompred::orient2d (
        _pa, _pb, _qq
        ) ;

    _rr = geompred::inball2d (
        _pa, _pb, _pc, _qq
        ) ;

    _rr = geompred::inball2w (
        _pa, _pb, _pc, _qq
        ) ;

    return 0 ;
}



