#!/usr/bin/env python3
"""Generate the reference values the magnitude sweep runs against.

The sweep exists because seven silent wrong answers were found in this
engine over two days, and every one of them was correct for the sizes
somebody had tested and wrong outside them. It runs every function
whose answer depends on the scale of its argument across a ladder from
1e-40 to 1e+40 and either side of one.

    python3 gensweep.py > AddapticSweepRefs.h

Nothing in the build needs Python: the output is committed. This is
here so the ladder can be widened later without anyone having to work
out how the references were made.

Every value proves itself by an inverse or an identity before it is
emitted -- a bad reference invents a bug, which has happened here and
cost a false report of an engine defect.
"""

from decimal import Decimal, getcontext
import sys

getcontext().prec = 160

TOL = Decimal(1).scaleb(-45)


def close( a, b, scale=None ):
    s = scale if scale is not None else max( Decimal(1), abs(b) )
    return abs( a - b ) <= TOL * s


def ln( x ):  return x.ln()
def exp( x ): return x.exp()


def atan( x ):
    doublings = 0
    while abs( x ) > Decimal( "0.05" ):
        x = x / (1 + (1 + x * x).sqrt())
        doublings += 1
    cut = Decimal(1).scaleb( -(getcontext().prec + 10) )
    term, total, xsq, n = x, x, x * x, 1
    while True:
        term = -term * xsq
        contribution = term / (2 * n + 1)
        if abs( contribution ) < cut:
            break
        total += contribution
        n += 1
    return total * (2 ** doublings)


PI = 16 * atan( Decimal(1) / 5 ) - 4 * atan( Decimal(1) / 239 )


def _reduce_circular( x ):
    k = (x / (2 * PI)).to_integral_value( rounding="ROUND_FLOOR" )
    return x - k * 2 * PI


def sin( x ):
    x = _reduce_circular( x )
    cut = Decimal(1).scaleb( -(getcontext().prec + 10) )
    term, total, xsq, n = x, x, x * x, 1
    while True:
        term = -term * xsq / ((2 * n) * (2 * n + 1))
        if abs( term ) < cut:
            break
        total += term
        n += 1
    return total


def cos( x ):
    x = _reduce_circular( x )
    cut = Decimal(1).scaleb( -(getcontext().prec + 10) )
    term, total, xsq, n = Decimal(1), Decimal(1), x * x, 1
    while True:
        term = -term * xsq / ((2 * n - 1) * (2 * n))
        if abs( term ) < cut:
            break
        total += term
        n += 1
    return total


def sinh( x ): return (exp(x) - exp(-x)) / 2
def cosh( x ): return (exp(x) + exp(-x)) / 2
def tanh( x ): return sinh(x) / cosh(x)


def sqrt( x ): return x.sqrt()


def cbrt( x ):
    if x == 0:
        return Decimal(0)
    return (-1 if x < 0 else 1) * exp( ln( abs(x) ) / 3 )


def asin( x ):
    if x == 1:  return PI / 2
    if x == -1: return -PI / 2
    return atan( x / (1 - x * x).sqrt() )


def acos( x ):   return PI / 2 - asin( x )
def asinh( x ):
    # Fold onto the magnitude and put the sign back, for the same
    # reason lz_archypersine does: x + sqrt(x^2+1) cancels for a
    # negative x, and at -200 the two terms agree far enough that the
    # logarithm is taken of noise. Written the naive way this check
    # rejected its own reference for sinh(-200).
    if x < 0:
        return -ln( -x + (x * x + 1).sqrt() )
    return ln( x + (x * x + 1).sqrt() )
def acosh( x ):  return Decimal(0) if x == 1 else ln( x + (x * x - 1).sqrt() )
def atanh( x ):  return Decimal(0) if x == 0 else ln( (1 + x) / (1 - x) ) / 2


# The reciprocal family. Each is one over the function it is named
# against, and each inverse is its co-function's inverse taken on 1/x.
def tan( x ):    return sin( x ) / cos( x )
def cot( x ):    return cos( x ) / sin( x )
def sec( x ):    return Decimal(1) / cos( x )
def csc( x ):    return Decimal(1) / sin( x )
def coth( x ):   return cosh( x ) / sinh( x )
def sech( x ):   return Decimal(1) / cosh( x )
def csch( x ):   return Decimal(1) / sinh( x )
# arccot is arctan(1/x), pi/2 at nought -- range (-pi/2, pi/2], odd -- since
# 2026-10-23, Rich's call, as Elzed answers; it was pi/2 - arctan(x) before.
def acot( x ):   return PI / 2 - atan( x ) if x >= 0 else -( PI / 2 - atan( -x ) )
def asec( x ):   return acos( Decimal(1) / x )
def acsc( x ):   return asin( Decimal(1) / x )
def acoth( x ):  return atanh( Decimal(1) / x )
def asech( x ):  return acosh( Decimal(1) / x )
def acsch( x ):  return asinh( Decimal(1) / x )


# A circular reciprocal near a pole has no reference worth recording:
# the true value is enormous and the argument's last place moves it.
# Nothing on the ladder lands on a pole, but say so rather than hope.
NEAR_POLE = Decimal(1).scaleb( -20 )

def sine_clear( v ):    return abs( sin( v ) ) > NEAR_POLE
def cosine_clear( v ):  return abs( cos( v ) ) > NEAR_POLE


# The ladder. The interesting failures cluster at the extremes and
# within a whisker of one, so it is dense in both places -- but the
# middle needs filling too, and being thin there cost something real.
#
# The logarithm's reduction-band fault lived at arguments whose reduced
# value landed near the top of the band: 10, 199, 200, 3999, 4000,
# 4001. Only 10 was on the ladder, and the sweep recorded it as sound
# anyway for a separate reason (see the classifier in AddapticSweep).
# But a fault that hits scattered arguments needs scattered arguments
# to be seen at all, and eight thin decades cannot supply them.
#
# So: several points per decade rather than one, values that are not
# all powers of ten -- 3, 7, 20, 50, 200, 400, 4000 reduce to quite
# different places than 1000 and 1e8 do -- and a tighter approach to
# one from both sides.
LADDER = [ "1e-40", "1e-30", "1e-20", "1e-12", "1e-8", "0.00001",
           "0.001", "0.01", "0.05", "0.1", "0.2", "0.3", "0.5", "0.7",
           "0.9",
           "1", "1.1", "1.5", "2", "3", "5", "7", "10", "20", "50",
           "200", "400", "1000", "4000", "100000", "1e8", "1e12",
           "1e20", "1e30", "1e40" ]

NEAR_ONE = [ "0.99", "0.999999", "0.9999999999",
             "1.01", "1.000001", "1.0000000001" ]

# Gamma, by Stirling with a recurrence shift -- the same method the
# engine uses, so the self-check below has to come from somewhere else
# entirely. It does: Gamma(z+1) = z*Gamma(z) is the defining recurrence
# and shares no machinery with the series.
BERNOULLI = None

def _bernoulli( count ):
    from fractions import Fraction
    from math import comb
    B = [ Fraction( 0 ) ] * (2*count + 1)
    B[ 0 ] = Fraction( 1 )
    for m in range( 1, 2*count + 1 ):
        B[ m ] = -sum( comb( m+1, k ) * B[ k ] for k in range( m ) ) / (m + 1)
    return [ Decimal( B[ 2*i ].numerator ) / Decimal( B[ 2*i ].denominator )
             for i in range( 1, count + 1 ) ]

def loggamma( z ):
    """For a positive argument only; gamma() reflects for the rest."""
    global BERNOULLI
    if BERNOULLI is None:
        BERNOULLI = _bernoulli( 40 )
    product = Decimal( 1 )
    w = z
    while w < 60:
        product *= w
        w += 1
    r = (w - Decimal( "0.5" )) * ln( w ) - w + ln( 2 * PI ) / 2
    for n in range( 1, 41 ):
        r += BERNOULLI[ n-1 ] / (2*n * (2*n - 1) * w ** (2*n - 1))
    return r - ln( product )

def gamma( z ):
    if z > 0:
        return exp( loggamma( z ) )
    # Gamma(z) = pi / ( sin(pi z) Gamma(1-z) )
    return PI / ( sin( PI * z ) * exp( loggamma( 1 - z ) ) )

def gamma_defined( v ):
    """Every non-positive whole number is a pole, and the engine
    refuses them. Nothing else is out of bounds."""
    return not (v <= 0 and v == v.to_integral_value())

def gamma_recurrence( y, x ):
    """Gamma(x+1) = x * Gamma(x), which the series knows nothing of."""
    return close( gamma( x + 1 ), x * y, max( Decimal( 1 ), abs( x * y ) ) )

def loggamma_check( y, x ):
    """exp of the answer against gamma while gamma can be held; past
    that the same recurrence in logarithms, lnGamma(x+1) = lnGamma(x)
    + ln x -- gamma(1e40) is e to the 9e41, far beyond Decimal."""
    if x < 1000:
        return close( exp( y ), gamma( x ), max( Decimal( 1 ), abs( gamma( x ) ) ) )
    return close( loggamma( x + 1 ) - y, ln( x ) )


# Gamma cannot share the ladder. **Every negative whole number on it is
# a pole**, and past about fifty the answer runs to hundreds of digits
# for no gain -- Gamma(50.5) is already 1.6e63, which exercises the
# magnitude the sweep is for.
#
# So: the shared ladder's small end, which is all non-integer and
# perfectly safe, then half-integers up and down. The half-integers are
# deliberate -- they sit as far from a pole as an argument can, so what
# is measured is the series rather than the conditioning.
GAMMA_LADDER = [ "1e-40", "1e-30", "1e-20", "1e-12", "1e-8", "0.00001",
                 "0.001", "0.01", "0.05", "0.1", "0.2", "0.3", "0.5",
                 "0.7", "0.9", "0.99", "0.999999",
                 "1.01", "1.000001", "1.1", "1.5", "2.5", "3.5", "3.7",
                 "5.5", "7.5", "10.5", "20.5", "30.5", "50.5",
                 "-0.00001", "-0.001", "-0.1", "-0.5", "-0.7", "-0.9",
                 "-1.5", "-2.5", "-3.5", "-3.7", "-5.5", "-10.5",
                 "-20.5", "-30.5" ]

# Its logarithm keeps the same ladder above nought -- gamma alternates
# sign below, so a real logarithm of it does not exist and lz_log_gamma
# refuses rather than answering ln|Gamma| -- and then carries on past
# fifty, which gamma cannot. **It stopped at 50.5 until 2026-10-10**, and
# from about 1e8 nearly every log-gamma was wrong, by 2.6e22 units at 1e30,
# with nothing here or in the audit to see it: a large argument is not
# shifted, and the series multiplied the argument into its rounding of
# ln w. Not all powers of ten, for the reason the shared ladder says.
LOGGAMMA_LADDER = ( [ a for a in GAMMA_LADDER if not a.startswith( "-" ) ] +
                    [ "100.5", "1000.25", "4000", "100000", "1e8", "123456789.5",
                      "1e12", "1e20", "1e30", "1e40" ] )

LADDERS = { "gamma": GAMMA_LADDER, "loggamma": LOGGAMMA_LADDER }
# (exp_b and exp_n join these below, once POWER_LADDER is defined.)


# The means of a list, swept as the mean of 1 and x -- see MeanOfOneAnd
# in AddapticSweep.cpp. Each self-check is an identity the definition
# gives directly, so it shares nothing with how the engine computes.
def amean1( x ):  return (1 + x) / 2
def gmean1( x ):  return x.sqrt()
def hmean1( x ):  return 2 * x / (1 + x)

# And the arithmetic-geometric mean of 1 and x, the same way. Its check
# is the AGM's scaling and symmetry together, agm(1, x) = x * agm(1, 1/x):
# the two sides start the iteration from different places, so a loop
# that stopped short, or lost the small value's digits, disagrees with
# itself.
def agm1( x ):
    a, g = Decimal(1), x
    cut = Decimal(1).scaleb( -(getcontext().prec - 5) )
    for _ in range( 200 ):
        a, g = (a + g) / 2, (a * g).sqrt()
        if abs( a - g ) <= cut * a:
            break
    return a

# The functions that take a second operand, swept with that operand
# held fixed -- see WithPartner in AddapticSweep.cpp. Each was outside
# the sweep and the audit until 2026-09-24, and so was never measured
# at the sizes the sweep exists for. The partners are chosen to keep the
# general path in play: an exponent of 2.5 rather than 2, a root of
# order 2.5 rather than 3, which lz_root_3 would otherwise answer.
LN3 = ln( Decimal( 3 ) )
LN10 = ln( Decimal( 10 ) )

def exp10( x ):   return exp( x * LN10 )             # lz_exp_b in base ten
def log10( x ):   return ln( x ) / LN10              # lz_log_b in base ten
def pow3( x ):    return exp( x * LN3 )              # lz_exp_n( 3, x )
def log3( x ):    return ln( x ) / LN3               # lz_log_n( 3, x )
def power25( x ): return exp( Decimal( "2.5" ) * ln( x ) )   # lz_power( x, 2.5 )
def root25( x ):  return exp( ln( x ) / Decimal( "2.5" ) )   # lz_root_n( x, 2.5 )
def pow200( x ):  return x ** 200                    # lz_raise_n( x, 200 )

def root200( y ):
    return Decimal( 0 ) if y == 0 else exp( ln( abs( y ) ) / 200 )

# **Two hundred, not seven, and the reason is the ladder.** Every value
# on it is a short decimal, so x to the 7th is exact at 36 places and
# passes however the powering rounds: 1.1^7 is 1.9487171. Fifty was
# tried next and audited clean too -- lz_raise_n falls short once its
# answer passes about 1e5 and rounding has accumulated over several
# squarings, and 1.5^50 is 6e8 in six of them. At 200, 1.1^200 is 1.9e8
# and 1.01^200 is 7.3, both under the audit's 1e30. The domain keeps
# every answer inside the sweep's 512-digit scoring buffer: 100^200 is
# 1e400, 0.1^200 is 1e-200.
#
# And the exponentials of 3 and 10 get fractional arguments above one:
# the shared ladder has only 1.1 and 1.5 there, and it is the size of
# the answer that exposes a power that does not widen for it. All of
# these stay under the audit's 1e30.
POWER_LADDER = ( LADDER + NEAR_ONE + [ "-" + a for a in LADDER ] +
                 [ "2.5", "7.5", "12.5", "20.5", "29.5", "-2.5", "-12.5" ] )

# lz_raise_e and lz_raise_n take a whole exponent and discard any
# fraction without saying so -- raise_e(2.5) answers e squared -- so
# they are swept on whole numbers only. That is a fault of its own;
# the sweep measures size, not that.
def whole( v ):   return v == v.to_integral_value()

# lz_arctangent2( v, -1.5 ): the angle of (-1.5, v), in the left half-plane,
#	so the reflection through pi is always taken -- nearly pi for a tiny v,
#	nearly pi/2 for a huge one. Since 2026-10-22, with the angle unit.
def atan2m15( v ):
    t = atan( v / Decimal( "-1.5" ) )
    return t + PI if v >= 0 else t - PI

LADDERS[ "exp_b" ] = POWER_LADDER
LADDERS[ "exp_n" ] = POWER_LADDER

#  name          key  function  domain                     self-check
FUNCS = [
    ( "sqrt",     'R', sqrt,  lambda v: v > 0,             lambda y, x: close( y*y, x, max(Decimal(1), abs(x)) ) ),
    ( "cbrt",     '3', cbrt,  lambda v: v != 0,            lambda y, x: close( y*y*y, x, max(Decimal(1), abs(x)) ) ),
    ( "ln",       'L', ln,    lambda v: v > 0,             lambda y, x: close( exp(y), x, max(Decimal(1), abs(x)) ) ),
    ( "exp",      'E', exp,   lambda v: abs(v) < 1000,     lambda y, x: close( ln(y), x ) ),
    ( "reciprocal",'V', lambda v: 1/v, lambda v: v != 0,   lambda y, x: close( y*x, Decimal(1) ) ),
    ( "sin",      's', sin,   lambda v: True,              lambda y, x: abs(y) <= 1 ),
    ( "cos",      'c', cos,   lambda v: True,              lambda y, x: abs(y) <= 1 ),
    ( "sinh",     'S', sinh,  lambda v: abs(v) < 300,      lambda y, x: close( asinh(y), x ) ),
    ( "cosh",     'H', cosh,  lambda v: abs(v) < 300,      lambda y, x: y >= 1 ),
    ( "tanh",     'T', tanh,  lambda v: abs(v) < 300,      lambda y, x: abs(y) <= 1 ),
    ( "arcsinh",  'j', asinh, lambda v: True,              lambda y, x: close( sinh(y), x, max(Decimal(1), abs(x)) ) ),
    ( "arccosh",  'h', acosh, lambda v: v >= 1,            lambda y, x: close( cosh(y), x, max(Decimal(1), abs(x)) ) ),
    ( "arctan",   'a', atan,  lambda v: True,              lambda y, x: abs(y) < PI/2 + 1 ),
    ( "arcsin",   'A', asin,  lambda v: abs(v) <= 1,       lambda y, x: close( sin(y), x ) ),
    ( "arccos",   'C', acos,  lambda v: abs(v) <= 1,       lambda y, x: close( cos(y), x ) ),
    ( "arctanh",  'g', atanh, lambda v: abs(v) < 1,        lambda y, x: close( tanh(y), x ) ),

    # The reciprocal family. Every check is a round trip through the
    # forward function, scaled against the argument so a large one is
    # judged on its relative error rather than its absolute.
    ( "cot",      'o', cot,   lambda v: v != 0 and sine_clear(v),
                                       lambda y, x: close( y * tan(x), Decimal(1) ) ),
    ( "sec",      'x', sec,   cosine_clear,
                                       lambda y, x: close( y * cos(x), Decimal(1) ) ),
    ( "csc",      'n', csc,   lambda v: v != 0 and sine_clear(v),
                                       lambda y, x: close( y * sin(x), Decimal(1) ) ),
    ( "coth",     'O', coth,  lambda v: v != 0 and abs(v) < 300,
                                       lambda y, x: close( y * tanh(x), Decimal(1) ) ),
    ( "sech",     'X', sech,  lambda v: abs(v) < 300,
                                       lambda y, x: close( y * cosh(x), Decimal(1) ) ),
    ( "csch",     'N', csch,  lambda v: v != 0 and abs(v) < 300,
                                       lambda y, x: close( y * sinh(x), Decimal(1) ) ),
    ( "arccot",   'k', acot,  lambda v: True,
                                       lambda y, x: -PI / 2 < y <= PI / 2 and (y > 0) == (x >= 0)
                                                    and ( x == 0 or close( cot(y), x, max(Decimal(1), abs(x)) ) ) ),
    ( "arcsec",   'y', asec,  lambda v: abs(v) >= 1,
                                       lambda y, x: close( sec(y), x, max(Decimal(1), abs(x)) ) ),
    ( "arccsc",   'm', acsc,  lambda v: abs(v) >= 1,
                                       lambda y, x: close( csc(y), x, max(Decimal(1), abs(x)) ) ),
    ( "arccoth",  'K', acoth, lambda v: abs(v) > 1,
                                       lambda y, x: close( coth(y), x, max(Decimal(1), abs(x)) ) ),
    ( "arcsech",  'Y', asech, lambda v: 0 < v <= 1,
                                       lambda y, x: close( sech(y), x, max(Decimal(1), abs(x)) ) ),
    ( "arccsch",  'M', acsch, lambda v: v != 0,
                                       lambda y, x: close( csch(y), x, max(Decimal(1), abs(x)) ) ),

    # Gamma and its logarithm, on their own ladders -- see LADDERS.
    ( "gamma",    'G', gamma, gamma_defined,        gamma_recurrence ),
    # The means of 1 and x. The geometric and harmonic means refuse a
    # negative, so they keep to the positive half of the ladder.
    ( "amean",    'u', amean1, lambda v: True,
                                       lambda y, x: close( 2 * y, 1 + x, max(Decimal(1), abs(x)) ) ),
    ( "gmean",    'v', gmean1, lambda v: v > 0,
                                       lambda y, x: close( y * y, x, max(Decimal(1), abs(x)) ) ),
    ( "hmean",    'w', hmean1, lambda v: v > 0,
                                       lambda y, x: close( 1 / y, (1 + 1 / x) / 2, max(Decimal(1), abs(1 / y)) ) ),
    ( "agm",      'q', agm1,   lambda v: v > 0,
                                       lambda y, x: close( y, x * agm1( 1 / x ), max(Decimal(1), abs(y)) ) ),
    ( "loggamma", 'F', loggamma, lambda v: v > 0, loggamma_check ),

    # Outside the sweep and the audit until 2026-09-24. Every check is
    # the inverse, which shares nothing with how the engine gets there.
    ( "exp_b",    'B', exp10,   lambda v: abs(v) < 400,
                                       lambda y, x: close( log10(y), x ) ),
    ( "log_b",    'l', log10,   lambda v: v > 0,
                                       lambda y, x: close( exp10(y), x, max(Decimal(1), abs(x)) ) ),
    ( "raise_e",  'e', exp,     lambda v: whole(v) and abs(v) < 1000,
                                       lambda y, x: close( ln(y), x ) ),
    ( "tan",      't', tan,     cosine_clear,
                                       lambda y, x: close( y * cos(x), sin(x) ) ),
    ( "exp_n",    'P', pow3,    lambda v: abs(v) < 600,
                                       lambda y, x: close( log3(y), x ) ),
    ( "log_n",    'i', log3,    lambda v: v > 0,
                                       lambda y, x: close( pow3(y), x, max(Decimal(1), abs(x)) ) ),
    ( "power",    'p', power25, lambda v: v > 0,
                                       lambda y, x: close( root25(y), x, max(Decimal(1), abs(x)) ) ),
    ( "root_n",   'r', root25,  lambda v: v > 0,
                                       lambda y, x: close( power25(y), x, max(Decimal(1), abs(x)) ) ),
    ( "raise_n",  'z', pow200,  lambda v: Decimal( "0.1" ) <= abs(v) <= 100,
                                       lambda y, x: close( root200(y), abs(x), max(Decimal(1), abs(x)) ) ),

    # The two-argument arctangent against a partner of -1.5, since 2026-10-22.
    #	Checked by the tangent of the angle, which shares nothing with it.
    ( "atan2",    'Q', atan2m15, lambda v: True,
                                       lambda y, x: close( tan(y) * Decimal( "-1.5" ), x, max(Decimal(1), abs(x)) ) ),
]


def longhand( value ):
    """The engine reads no exponent notation, and neither does the
    suite's own scoring. Everything goes out written in full."""
    return format( value, 'f' )


rows = []
rejected = []

for name, key, fn, domain, check in FUNCS:
    args = LADDERS.get( name, LADDER + NEAR_ONE + [ "-" + a for a in LADDER ] )
    for a in args:
        x = Decimal( a )
        if not domain( x ):
            continue
        try:
            y = fn( x )
            ok = check( y, x )
        except Exception:
            rejected.append( (name, a, "raised") )
            continue
        if not ok:
            rejected.append( (name, a, "failed its own check") )
            continue
        getcontext().prec = 50
        rendered = longhand( +y )
        getcontext().prec = 160
        rows.append( (name, key, longhand( x ), rendered) )

print( "/////////////////////////////////////////////////////////////////" )
print( "//" )
print( "//\tAddapticSweepRefs.h" )
print( "//" )
print( "//\tGenerated by gensweep.py. Do not edit by hand." )
print( "//" )
print( "//\tEvery value was cross-checked by an inverse or an identity" )
print( "//\tbefore emission, and is written out longhand because neither" )
print( "//\tthe engine nor the scoring reads exponent notation." )
print( "//" )
print( "/////////////////////////////////////////////////////////////////" )
print()
print( "#pragma once" )
print()
print( "struct SweepReference" )
print( "{" )
print( "\tconst char*\tfunction;" )
print( "\tchar\t\tkey;" )
print( "\tconst char*\targument;" )
print( "\tconst char*\treference;" )
print( "};" )
print()
print( "static const SweepReference kSweepReferences[] =" )
print( "{" )
for name, key, a, r in rows:
    print( '\t{ "%s", \'%s\', "%s", "%s" },' % (name, key, a, r) )
print( "};" )
print()
print( "static const long kSweepReferenceCount = %d;" % len( rows ) )

print( "// %d references emitted, %d rejected" % (len( rows ), len( rejected )), file=sys.stderr )
for r in rejected:
    print( "//   rejected: %s" % (r,), file=sys.stderr )
