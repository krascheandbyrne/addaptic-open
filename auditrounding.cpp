/////////////////////////////////////////////////////////////////
//
//	auditrounding.cpp
//
//	The driver half of auditrounding.py: reads "key argument"
//	lines on stdin and prints "key|argument|result|code", using
//	the same one-letter keys AddapticSweep dispatches on.
//
//	Nothing in the build needs this. See auditrounding.py.
//
//	Copyright (c) 2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include "LZNumber.h"
#include "LZMath.h"
static long MeanOfOneAnd( char k, CLZNumber* a, CLZNumber* r ) {
    CLZNumber one; one.Assign( 1L ); CLZNumber* list[] = { &one, a };
    return (k == 'u') ? lz_amean(list,2,r) : (k == 'v') ? lz_gmean(list,2,r) :
           (k == 'q') ? lz_agm(&one,a,r) : lz_hmean(list,2,r); }
// The same partners as AddapticSweep's WithPartner...
static long WithPartner( char k, CLZNumber* a, CLZNumber* r ) {
    CLZNumber p;
    switch ( k ) {
    case 'P': p.Assign(3L);    return lz_exp_n(&p,a,r);   case 'i': p.Assign(3L);    return lz_log_n(&p,a,r);
    case 'p': p.Assign("2.5"); return lz_power(a,&p,r);   case 'r': p.Assign("2.5"); return lz_root_n(a,&p,r);
    case 'Q': p.Assign("-1.5"); return lz_arctangent2(a,&p,r);
    default:  p.Assign(200L);  return lz_raise_n(a,&p,r); } }
static long Invoke( char k, CLZNumber* a, CLZNumber* r )
{
    switch ( k ) {
    case 'R': return lz_root_2(a,r);          case '3': return lz_root_3(a,r);
    case 'L': return lz_log_e(a,r);           case 'E': return lz_exp_e(a,r);
    case 'V': return lz_reciprocal(a,r);      case 's': return lz_sine(a,r);
    case 'c': return lz_cosine(a,r);          case 'S': return lz_hypersine(a,r);
    case 'H': return lz_hypercosine(a,r);     case 'T': return lz_hypertangent(a,r);
    case 'j': return lz_archypersine(a,r);    case 'h': return lz_archypercosine(a,r);
    case 'a': return lz_arctangent(a,r);      case 'A': return lz_arcsine(a,r);
    case 'C': return lz_arccosine(a,r);       case 'g': return lz_archypertangent(a,r);
    case 'o': return lz_cotangent(a,r);       case 'x': return lz_secant(a,r);
    case 'n': return lz_cosecant(a,r);        case 'O': return lz_hypercotangent(a,r);
    case 'X': return lz_hypersecant(a,r);     case 'N': return lz_hypercosecant(a,r);
    case 'k': return lz_arccotangent(a,r);    case 'y': return lz_arcsecant(a,r);
    case 'm': return lz_arccosecant(a,r);     case 'K': return lz_archypercotangent(a,r);
    case 'Y': return lz_archypersecant(a,r);  case 'M': return lz_archypercosecant(a,r);
    case 'G': return lz_gamma(a,r);           case 'F': return lz_log_gamma(a,r);
    case 'u': case 'v': case 'w': case 'q': return MeanOfOneAnd(k,a,r);
    case 'B': return lz_exp_b(a,r);           case 'l': return lz_log_b(a,r);
    case 'e': return lz_raise_e(a,r);         case 't': return lz_tangent(a,r);
    case 'P': case 'i': case 'p': case 'r': case 'z': case 'Q': return WithPartner(k,a,r);
    default:  return kFailed; }
}
int main(){
    char line[600];
    lzu_setup();
    while ( fgets(line,sizeof(line),stdin) ) {
        // The key and its argument, two words -- not by sscanf, which MSVC deprecates...
        std::istringstream words(line);
        std::string key, arg;
        if ( !(words >> key >> arg) ) continue;
        CLZNumber x; x.Assign(arg.c_str());
        CLZNumber y;
        long c = Invoke(key[0],&x,&y);
        // The text is the caller's to free, and was not -- a string a case, until the sanitizer job ran this...
        char* sz = y.GetStringRep();
        printf("%s|%s|%s|%ld\n", key.c_str(), arg.c_str(), sz ? sz : "(null)", c);
        ::free(sz);
    }
    lzu_takedown(); return 0;
}
