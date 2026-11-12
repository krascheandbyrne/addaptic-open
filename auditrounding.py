#!/usr/bin/env python3
"""Audit every sweep case against the correctly rounded answer.

    c++ -std=c++17 -O2 -I. auditrounding.cpp LZMath.cpp LZNumber.cpp -o /tmp/auditdump
    python3 auditrounding.py /tmp/auditdump

AddapticSweep asks a weaker question than this one. It counts agreeing
digits and compares that count against what the engine managed last
time, which makes it a regression gate. This asks the absolute
question instead: is the answer the correctly rounded value at the
working precision, and if not, by how much and why.

That distinction found something once already. The sweep recorded
ln(10) as sound at 34 digits for as long as the logarithm's
reduction-band fault existed, because 34 was better than the flat
threshold it was compared against. Nothing in a regression gate can
report that today's answer is already wrong.

Nothing in the build needs this. It is here rather than in a scratch
directory because the last such harness was nearly lost that way.

Four outcomes, and the last is the only one worth chasing:

  correctly rounded   the answer at 36 places is the true value
                      correctly rounded to 36 places
  within a few ulp    off in the last place or two, which is the
                      series' own accumulated rounding
  series baseline     agreeing to at least kBaselineDigits significant
                      digits, but not correctly rounded and not
                      explained by either of the two above. Good
                      answers that are short of what their size
                      affords.

                      **Read what lands here.** It is meant to stay
                      small enough to be read case by case, the way
                      the sweep's stderr list is, and it is listed at
                      the end of the run for that reason.

                      What is in it today is four gamma cases: its
                      answers stop at 52 significant digits, the
                      Bernoulli table's fifty-digit width stopping them
                      before kGammaMaxWidening does.

                      It held 51 more for one day, 2026-09-24: power,
                      log_n, exp_b, exp_n, raise_n and root_n, added to
                      the sweep that morning, did not widen for the
                      size of their answer, so a large one carried 36
                      to 39 significant digits and fell short of its
                      36th place by as much as 4e29 units. This bucket
                      took them because it asks only for 35
                      significant digits. They widen now.

                      It held five exponential cases until the
                      magnitude widening went into lz_raise_e, and
                      those are worth remembering: they were the exact
                      five that kRaiseEGuardDigits at 32 used to buy at
                      1.7 times the cost of every hyperbolic. Widening
                      by the answer's own size lifted them for
                      nothing, and 32 now buys nothing at all.

                      The bucket used to describe itself as the
                      trigonometric series settling at 35 digits, with
                      sin(1) as the example. That stopped being true
                      when the circular functions got their guard
                      digits, and nothing noticed for some time,
                      because nothing listed the contents.
  inherited           a reciprocal no better than the small thing it
                      divides by. A cotangent of 1e-12 divides by a
                      sine of about 1e-12, and a sine held to w places
                      carries about w - 12 significant digits, so the
                      quotient can carry no more.

                      Corroborated where it can be: the answer is
                      scored against the exact quotient of the engine's
                      OWN numerator and denominator, and if that agrees
                      far better than the answer does with the truth,
                      the division itself is blameless.

                      That corroboration only works while this script
                      can reproduce the intermediate. Once lz_cotangent
                      began guarding across its own division it stopped
                      being able to: the sine this script asks for is
                      the public one at 36 places, while the cotangent
                      is dividing a wider one it made privately. So the
                      ceiling below is the internal width, PLACES plus
                      the guard, and the corroboration is used when it
                      agrees and not required when it cannot run. The
                      engine improving is not the same as the audit
                      breaking, and the first time this fired it was
                      the audit that was out of date.
"""

import collections
import contextlib
import io
import re
import subprocess
import sys

from decimal import Decimal as D, getcontext, ROUND_HALF_EVEN

with contextlib.redirect_stdout( io.StringIO() ), \
     contextlib.redirect_stderr( io.StringIO() ):
    import gensweep as G

getcontext().prec = 260

PLACES = 36
QUANTUM = D( 1 ).scaleb( -PLACES )

# What the series actually deliver, measured rather than assumed:
# sin(1) agrees with the truth to 35 significant digits, and every
# function standing on the trigonometric or exponential series lands
# on the same number. Raise this if the series ever get guard digits.
BASELINE_DIGITS = 35

# The width the circular functions work at internally, from
# kTrigGuardDigits in LZMath.h. A reciprocal of a small argument is
# limited by the divisor computed at this width, not at PLACES.
TRIG_GUARD = 4

# The same key-to-function map the sweep dispatches on.
FUNCTIONS = {
    'R': G.sqrt,  '3': G.cbrt,  'L': G.ln,    'E': G.exp,
    'V': lambda v: D( 1 ) / v,
    's': G.sin,   'c': G.cos,   'S': G.sinh,  'H': G.cosh,  'T': G.tanh,
    'j': G.asinh, 'h': G.acosh, 'a': G.atan,  'A': G.asin,
    'C': G.acos,  'g': G.atanh,
    'o': G.cot,   'x': G.sec,   'n': G.csc,
    'O': G.coth,  'X': G.sech,  'N': G.csch,
    'k': G.acot,  'y': G.asec,  'm': G.acsc,
    'K': G.acoth, 'Y': G.asech, 'M': G.acsch,
    'G': G.gamma, 'F': G.loggamma,
    'u': G.amean1, 'v': G.gmean1, 'w': G.hmean1, 'q': G.agm1,
    'B': G.exp10,  'l': G.log10,  'e': G.exp,    't': G.tan,
    'P': G.pow3,   'i': G.log3,   'p': G.power25, 'r': G.root25, 'z': G.pow200,
    'Q': G.atan2m15,
}

# The reciprocals, and what each divides: (key, numerator, denominator).
# A numerator of None means one.
RECIPROCALS = {
    'o': ( 'c', 's' ),  'n': ( None, 's' ),  'x': ( None, 'c' ),
    'O': ( 'H', 'S' ),  'N': ( None, 'S' ),  'X': ( None, 'H' ),
}


def read_cases():
    text = open( "AddapticSweepRefs.h" ).read()
    return [ (m.group( 1 ), m.group( 2 ), m.group( 3 ))
             for m in re.finditer(
                 r'\{ "([^"]+)", \'(.)\', "([^"]*)", "([^"]*)" \}', text ) ]


def read_kinds():
    text = open( "AddapticSweepExpect.h" ).read()
    return { (m.group( 1 ), m.group( 2 )): int( m.group( 4 ) )
             for m in re.finditer(
                 r'\{ "([^"]+)", "([^"]*)", (-?\d+), (\d+) \}', text ) }


def run( binary, wanted ):
    """wanted is a list of (key, argument). Returns {(key, arg): (value, code)}."""
    stdin = "".join( "%s %s\n" % pair for pair in wanted )
    result = subprocess.run( [ binary ], input=stdin,
                             capture_output=True, text=True )

    # A driver that failed fails the audit, whatever it printed first -- a
    # crash, or a sanitizer's report. LeakSanitizer reports at exit, after
    # every answer, so read as answers alone a leak passed unseen; and
    # until 2026-10-17 the exit status was not read at all.
    if result.returncode != 0:
        sys.stderr.write( result.stderr[ -6000: ] )
        sys.exit( "the audit's driver exited %d" % result.returncode )

    out = result.stdout
    answers = {}
    for line in out.splitlines():
        parts = line.split( "|" )
        if len( parts ) == 4:
            answers[ (parts[ 0 ], parts[ 1 ]) ] = ( parts[ 2 ], int( parts[ 3 ] ) )
    return answers


def leading_zeros( text ):
    """Zeros between the radix mark and the first significant digit."""
    body = text.lstrip( '-' )
    if '.' not in body:
        return 0
    whole, fraction = body.split( '.', 1 )
    if whole.lstrip( '0' ):
        return 0
    count = 0
    for ch in fraction:
        if ch != '0':
            break
        count += 1
    return count


def agreement( a, b ):
    if a == b:
        return 99
    return int( -(abs( a - b ) / abs( b )).log10() )


def main():
    if len( sys.argv ) < 2:
        sys.exit( __doc__ )

    binary = sys.argv[ 1 ]
    cases = read_cases()
    kinds = read_kinds()

    wanted = [ (key, arg) for _, key, arg in cases ]

    # The reciprocals want their own numerator and denominator too, so
    # the inherited check has something to compare against.
    for _, key, arg in cases:
        if key in RECIPROCALS:
            num, den = RECIPROCALS[ key ]
            if num:
                wanted.append( (num, arg) )
            wanted.append( (den, arg) )

    answers = run( binary, wanted )

    exact = 0
    near = collections.Counter()
    baseline = []
    inherited = []
    unexplained = []
    skipped = 0

    for name, key, arg in cases:
        # Poles and the 2*pi table are accounted for elsewhere.
        if kinds.get( (name, arg) ) in ( 2, 3 ):
            skipped += 1
            continue

        rendered, code = answers[ (key, arg) ]
        if code != 0:
            skipped += 1
            continue

        true = FUNCTIONS[ key ]( D( arg ) )

        # Past this the answer has more integer digits than 36 places
        # of fraction, and "correctly rounded to 36 places" stops being
        # a reasonable thing to ask of a series.
        if abs( true ) > D( 10 ) ** 30:
            skipped += 1
            continue

        got = D( rendered )
        want = true.quantize( QUANTUM, rounding=ROUND_HALF_EVEN )
        have = got.quantize( QUANTUM, rounding=ROUND_HALF_EVEN )
        ulps = int( (have - want).scaleb( PLACES ) )

        if ulps == 0:
            exact += 1
            continue

        if abs( ulps ) <= 2:
            near[ ulps ] += 1
            continue

        # Further out than a couple of places. If it is a reciprocal,
        # ask whether the loss is its divisor's rather than its own.
        if key in RECIPROCALS:
            num_key, den_key = RECIPROCALS[ key ]
            denominator = D( answers[ (den_key, arg) ][ 0 ] )
            numerator = ( D( answers[ (num_key, arg) ][ 0 ] )
                          if num_key else D( 1 ) )

            if denominator != 0:
                true_den = FUNCTIONS[ den_key ]( D( arg ) )
                against_truth = agreement( got, true )
                against_own = agreement( got, numerator / denominator )
                divisor_itself = agreement( denominator, true_den )

                # What the argument's own size allows, at the width
                # the function works to internally.
                ceiling = PLACES + TRIG_GUARD - leading_zeros( arg )

                # Blameless when the answer sits AT that ceiling --
                # not merely under it. A genuine fault would land well
                # below and still be reported, which is what keeps this
                # bucket from becoming a place to hide things.
                at_ceiling = ( ceiling - 2 ) <= against_truth <= ( ceiling + 1 )

                # Corroboration where it can still be had: see the note
                # at the top on why it cannot always run.
                corroborated = ( against_own > against_truth + 3 or
                                 abs( divisor_itself - against_truth ) <= 1 )

                # The direct evidence, which needs no model of the
                # internal width at all: the answer agrees with the
                # truth exactly as well as its divisor does, and
                # dividing the engine's own numerator by its own
                # denominator reproduces the answer. The loss is the
                # divisor's and the division is blameless.
                #
                # This has to stand on its own rather than behind the
                # ceiling, because the ceiling is pessimistic wherever
                # a function works wider than TRIG_GUARD -- the
                # hyperbolic reciprocals stand on the exponential and
                # its sixteen places of guard, so csch at 1e-20 holds
                # 40 digits where the model predicts 21 and fell
                # through to the baseline bucket for it.
                #
                # It cannot hide a bad divisor: every divisor here is
                # a sweep case in its own right and is scored on its
                # own row.
                matches_divisor = ( abs( divisor_itself - against_truth ) <= 1 and
                                    against_own > against_truth + 3 )

                if ( matches_divisor or at_ceiling or
                     ( against_truth <= ceiling and corroborated ) ):
                    inherited.append( (name, arg, against_truth) )
                    continue

        digits = agreement( got, true )

        if digits >= BASELINE_DIGITS:
            baseline.append( (name, arg, digits) )
            continue

        unexplained.append( (name, arg, ulps, digits) )

    total = len( cases ) - skipped
    print( "%d cases examined (%d poles, tables and huge answers set aside)\n"
           % (total, skipped) )
    print( "  correctly rounded at %d places : %d" % (PLACES, exact) )
    print( "  within two places of it        : %d" % sum( near.values() ) )
    print( "  at the series baseline (%d dig) : %d"
           % (BASELINE_DIGITS, len( baseline )) )
    print( "  inherited from a small divisor : %d" % len( inherited ) )
    print( "  unexplained                    : %d" % len( unexplained ) )

    if near:
        print( "\nlast-place differences:" )
        for u in sorted( near ):
            print( "   %+d ulp : %d" % (u, near[ u ]) )

    if inherited:
        print( "\ninherited, worst first "
               "(the answer is as good as its divisor and no better):" )
        for name, arg, digits in sorted( inherited, key=lambda r: r[ 2 ] )[ :12 ]:
            print( "   %-9s %-32s %d significant digits" % (name, arg, digits) )

    # Listed, not merely counted. A bucket nobody reads is where a
    # fault goes to be mistaken for a known limit -- the same argument
    # the sweep's stderr review list is built on. Ten cube roots sat
    # here losing a digit per power of the base for as long as this
    # only printed the number.
    if baseline:
        print( "\nat the series baseline -- good, but short of what "
               "their size affords:" )
        for name, arg, digits in sorted( baseline, key=lambda r: (r[ 0 ], r[ 1 ]) ):
            print( "   %-9s %-34s %d significant digits" % (name, arg, digits) )

    if unexplained:
        print( "\nUNEXPLAINED -- these are the ones to look at:" )
        for name, arg, ulps, digits in sorted(
                unexplained, key=lambda r: r[ 3 ] )[ :40 ]:
            print( "   %-9s %-32s %d digits (%+d ulp)"
                   % (name, arg, digits, ulps) )

    return 1 if unexplained else 0


if __name__ == "__main__":
    sys.exit( main() )
