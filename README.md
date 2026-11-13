# Addaptic

**An arbitrary-precision math library for C++.**

Copyright (c) 2000-2026 by Krasche & Byrne LLC. Released under the MIT License;
see [LICENSE](LICENSE).

Addaptic does exact arithmetic on numbers of unlimited length. Numbers are
stored as arrays of digits rather than as machine floats, so a value keeps
every digit it is given — no binary representation error, no 53-bit ceiling —
and every operation up through the transcendental functions is built from
digit-by-digit primitives. **Every answer is correctly rounded** to the working
precision, in any base from 2 to 16.

```
first operand  >> 99999999999999999999
operator       >> *
second operand >> 99999999999999999999

99999999999999999999 * 99999999999999999999 = 9999999999999999999800000000000000000001
```

## The number: `CLZNumber`

A number is an array of `long`, one digit of its base per element, stored
least-significant-first, plus the state that gives those digits meaning:

| Member | Meaning |
|---|---|
| `m_lSign` | `kPositive`, `kZero`, or `kNegative` — sign is held separately from the digits |
| `m_lRadix` | how many of the digits are fractional (where the point sits, counted from the right) |
| `m_lBase` | any base from **2 to 16**, with live conversion via `SetBase()` |
| `m_lPrecision` | per-number working precision, default 36 places |
| `m_rmMode`, `m_rdDigit` | rounding policy (below) |
| `m_lFlags` | zero, infinity, quiet NaN, signalling NaN, divide-by-zero |
| `m_bConstant` | write-protects the number as an operation target |

Storage grows on demand — `AlignStorage`, `PadLeft`/`PadRight`,
`ShiftLeft`/`ShiftRight`, `TrimLeading`/`TrimTrailing`/`TrimAll`.

`SetBase()` re-expresses the value in the new base. Integers convert
exactly. Fractions cannot always — a tenth has no terminating form in
base two — so they are carried to the working precision scaled by the
ratio of the logs, and rounded half-even at the cut. A value that does
not terminate survives a round trip to that precision rather than to
the last digit. A long conversion is done by multiplication rather than a
digit at a time, so it costs about what a product of that length does.

**Twelve rounding modes** are available: down, half-up, half-even, half-down,
half-odd, ceiling, floor, up, truncate, alternate, random, and no-op — applied
at one of three points (zero, last digit, or the precision limit). Half-even
is the default, as in the decimal standards: of the modes in common use it is
the only one that does not push a long series of results steadily in one
direction.

Numbers can be constructed or assigned from a string (`const char*`), a `long`,
a `double`, or another `CLZNumber`.

## The math: `LZMath`

Three internal primitives do the real work on digit arrays — `lzi_add_op`,
`lzi_subtract_op`, and `lzi_divide_op` (long division, several digits to a
limb once the operands are long enough to pay for it). Multiplication does the
same, packing digits into limbs past the smallest products, in every base;
in base ten it splits them Karatsuba's way from 856 digits and takes a number
theoretic transform from 130,870. Everything else is bootstrapped from those:

- **Arithmetic** — `lz_add`, `lz_subtract`, `lz_multiply` (digit by digit for
  the smallest products, packed limbs from there, Karatsuba on those limbs for
  long ones, a transform for the largest), `lz_divide` (long division on
  packed limbs, or for a long divisor a Newton reciprocal and a product),
  `lz_divide_integer`, `lz_modulo`, `lz_remainder`,
  `lz_remainder_near`, `lz_abs`, `lz_negate`, `lz_plus`, `lz_reciprocal`
- **Comparison** — `lz_compare`, `lz_max`, `lz_min`, `lz_same_quantum`,
  each answering `kUnordered` when either side is a NaN or an infinity
- **Powers and roots** — `lz_square`, `lz_cube`, `lz_power`, `lz_raise_n`
  and `lz_raise_e` (by repeated squaring, so the cost follows the number
  of digits in the exponent rather than its size — and `lz_power`, `lz_exp_n`
  and `lz_exp_b` square too for a whole exponent, taking no logarithm),
  `lz_root_2` (Newton–Raphson, and at wide widths the
  arithmetic-geometric mean's own square root), `lz_root_3` (Newton, at wide
  widths on the reciprocal cube root), `lz_root_n` (a power when one over the
  order is whole, otherwise e to the logarithm over the order). A root of
  whole order exactly on a half is decided exactly, never by a coin toss.
- **Exponentials** — `lz_exp_e` (Taylor series against precomputed
  reciprocal-factorial constants, and at wide widths binary splitting a
  chunk of the fraction at a time), `lz_exp_b`, `lz_exp_n`
- **Logarithms** — `lz_log_e` (a series after argument reduction, sized from
  the argument it is left with; at wide widths Newton on the exponential, or
  the arithmetic-geometric mean once a width is in use), `lz_log_b`, `lz_log_n`
- **Trigonometry** — `lz_sine`, `lz_cosine`, `lz_tangent`, `lz_secant`,
  `lz_cosecant`, `lz_cotangent` (Taylor series after the argument is folded
  into [0, π/4], and at wide widths sine and cosine together by binary
  splitting, the tangent and cotangent from one such run); `lz_arcsine`,
  `lz_arccosine`, `lz_arctangent`, `lz_arccotangent`, `lz_arcsecant`,
  `lz_arccosecant`, all through the arcsine (a 23-coefficient series with the
  angle halved first, and at wide widths Newton on the sine and cosine) —
  `lz_arccotangent` is arctan(1/x), π/2 at nought, so its range is (−π/2,
  π/2] and arccot(−1) is −π/4, as Elzed answers; and
  `lz_arctangent2`, the angle of the point (x, y) — y first, as C's `atan2`
  takes them — in (−π, π], 0 at the origin. **All fourteen read or write
  their angle in radians, degrees, grads or millicles**, as a host sets
  `geAngleUnit` for its thread; outside radians an angle is reduced by
  whole turns exactly, so sin(30°) is exactly 0.5 and tan(90°) a pole.
  `lz_sine_cosine` gives the sine and the cosine of one angle together, each
  exactly what its own function answers, usually for about half the cost of
  the two
- **Hyperbolics** — `lz_hypersine`, `lz_hypercosine`, `lz_hypertangent`,
  `lz_hypercotangent`, `lz_hypersecant`, `lz_hypercosecant`, each from a single
  e^x; and their inverses, `lz_archypersine`, `lz_archypercosine`,
  `lz_archypertangent`, `lz_archypercotangent`, `lz_archypersecant`,
  `lz_archypercosecant`, from the logarithm and the square root, with a series
  where those would lose digits. `lz_hypersine_hypercosine` gives sinh and
  cosh together from one exponential, each exactly what its own function
  answers
- **Gamma** — `lz_gamma` and `lz_log_gamma`: a whole number exactly, as a
  factorial; anything else above nought from Stirling's series after a shift;
  below nought by reflection; the poles refused
- **Combinatorics** — `lz_factorial`, `lz_combinations`, `lz_permutations`,
  all on the whole numbers and exact; choosing or arranging more than there
  are is nought
- **Means** — `lz_amean`, `lz_gmean`, `lz_hmean` (arithmetic, geometric and
  harmonic means of a list of any length), and `lz_agm` (the
  arithmetic-geometric mean of two values that are not negative)
- **Decimal-standard operations** — `lz_quantize`, `lz_normalize`,
  `lz_round`, `lz_round_integral`, `lz_round_last_digit`, `lz_integer`,
  `lz_fraction`
- **Complex numbers, stages 1 to 6** — a `CLZComplex` holds two
  `CLZNumber`s, the real and imaginary parts in one base; `lzc_add`,
  `lzc_subtract`, `lzc_negate`, `lzc_conjugate`, `lzc_multiply`,
  `lzc_square`, `lzc_divide` and `lzc_abs`, each part correctly rounded from
  exact intermediates — the four plain operations in each target part's own
  rounding mode, the square and the modulus as functions do — and `lzc_arg`,
  the argument, correctly rounded in the thread's angle unit; and `lzc_polar`,
  a value from its modulus and angle, each part correctly rounded and exact
  where the mathematics is — 1 at 180° is exactly −1 + 0i; and `lzc_log_e`,
  `lzc_exp_e` and `lzc_root_2`, the principal values, each part correctly
  rounded — √(3 + 4i) is exactly 2 + i and ln(−1) is πi; and the twelve
  circular and hyperbolic functions, `lzc_sine` to `lzc_hypercosecant`, the
  circular six over the whole angle in the thread's unit; and `lzc_power`,
  `lzc_raise_n`, `lzc_root_n`, `lzc_root_3`, `lzc_exp_b`, `lzc_log_n` and
  `lzc_log_b`, exponents, orders and bases complex — (2 + i)³ is exactly 2 +
  11i and the cube root of −8 is −2; and the eight inverses that take square
  roots, `lzc_arcsine` to `lzc_archypercosecant`, the circular four over the
  whole angle in the thread's unit — in degrees arcsin(1 + 1.5i) has a real
  part of exactly 30; and the four that take logarithms, `lzc_arctangent`,
  `lzc_arccotangent`, `lzc_archypertangent` and `lzc_archypercotangent`,
  arccot being arctan(1/z) — in degrees arctan(0.4 + 0.2i) has a real part of
  exactly 22.5.

Every function takes an optional target: pass one and the result lands there,
omit it and the left operand is updated in place. Constant-flagged numbers are
never written to. The three means of a list are the exception — they take an
array of values and a count, have no left operand to update, and so require the
target:

```cpp
CLZNumber* values[] = { &a, &b, &c };
lz_gmean( values, 3, &result );
```

`lz_exp_b` raises the operand's own base and answers in it, so `b` squared
comes back as `100` whatever `b` is. It works in that base throughout, with
ln b worked out there, as every function does with the constants it needs.

Two operands need not share a base. The arithmetic operators, the division and
remainder family, and `lz_compare`, `lz_max` and `lz_min` bring the pair into
the left operand's base first, so `0x15 + 0b1010` is `0x1F` and `0b1010 + 0x15`
is `0b11111` — the same thirty-one, written the way the left-hand side was. The
right operand is converted on a copy and never modified.

**Every answer comes back in the left operand's base** — the operand's, for a
function of one, and the first value's, for a mean of a list — whatever base the
target was in before, and in place or not. So the sine of a hexadecimal value is
hexadecimal, and a precision of 36 means 36 places of that base: about eleven
decimal digits in base two, forty-three in base sixteen. `lz_quantize` and
`lz_same_quantum` read no digits, only how many places each operand has after
its point, so for them a place is counted in each operand's own base.

Call `lzu_setup()` before use and `lzu_takedown()` when finished — they build
and release the shared table of high-precision constants.

An operation that has no answer to give says so twice: it zeroes the result and
sets a flag on it — `kfInfinity`, `kfQNaN`, `kfDiv0` — and returns a code.
`kFailed` for an argument outside the function's domain, and `kDivideByZero`
and `kUndefinedResult` for division. So `ln(0)`, `ln(-5)`, `sqrt(-9)`, `1/0`
and `tan(pi/2)` all report rather than returning a plausible-looking number.

Those flags then carry through the arithmetic and the functions alike, so a
failure early in a calculation is still visible at the end of it: `(1/0) * 2 + 2`
comes back flagged infinite rather than as a plain zero, and a NaN stays a NaN
through every operation that touches it.

Infinities here are stored as zero plus a flag and so carry no sign, which
decides the awkward cases: two infinities added or subtracted give a NaN,
because the cancelling case cannot be told from the reinforcing one. Their
product is infinite, and anything finite divided by an infinity is zero. The
same question settles the functions — `sinh`, `cosh` and the cube root are
infinite whichever sign went in, while `exp`, `tanh`, `sqrt`, `ln` and the
factorial would answer differently for each, so they give a NaN instead of
guessing.

An answer too large to hold is refused rather than worked out. Each function
sizes its answer before it starts, and one whose order of magnitude would pass
100,000 places — the whole places of a large answer, or the noughts after the
point of a small quotient, which keeps its significant digits — zeroes the
result and returns `kTooLarge`, flagged `kfInfinity` if it was too large. So
`exp(1E24)`, `factorial(1E6)` and `2^(1E9)` are refused at once, where they
used to run until memory gave out. An answer that rounds to nought at the
working precision is nought at once — `exp(-1E24)` — and `tanh` and `coth` far
from nought are ±1 without working out e^x. `lzu_set_answer_limit()` changes
the limit for the whole process; set it before starting threads.

### Threads

Any number of threads may compute at once. Each keeps its own working
precision, which starts at 36 places in every new thread whatever another has
set, so set `glWorkingPrecision` in the thread that computes. Each also keeps
its own cache of constants converted into other bases, released when the
thread ends.

A number may be read by several threads at once, as an operand, but not
written by one while another reads it — give each thread its own results.
Call `lzu_setup()` once before starting threads and `lzu_takedown()` after
they have finished. The `lzi_set_*` functions, which only the tests use,
change the whole process: call them before a second thread starts.

## The bench

`Addaptic` is a console harness for exercising the engine. It prompts for a
first operand, an operator, and — when the operator needs one — a second
operand. Entering `r` as the first operand reuses the previous result, so
calculations can be chained. A blank first operand exits.

Each answer is printed with the digits-array length, digits used, and radix
position, which is usually the interesting part.

**Binary operators**

| Key | Operation | Key | Operation |
|---|---|---|---|
| `+` | add | `?` | compare |
| `-` | subtract | `>` | maximum |
| `*` | multiply | `<` | minimum |
| `/` | divide | `q` | quantize |
| `^` | raise to power | `=` | same quantum |
| `%` | modulo | `:` | trim (`1` leading, `2` trailing, else all) |
| `R` | round | `0` | round integral |
| `9` | round last digit | | |

The three rounding operators take a mode letter as the second operand:
`D` down, `u` half-up, `e` half-even, `d` half-down, `o` half-odd,
`C` ceiling, `F` floor, `U` up, `A` alternate, `R` random.

**Unary operators**

| Key | Operation | Key | Operation |
|---|---|---|---|
| `\|` | absolute value | `n` | negate |
| `p` | plus (no-op copy) | `~` | normalize |
| `L` | natural log | `l` | log in current base |
| `M` | order of magnitude | `B` | set base — the base goes in the *first* operand (`2`, `8`, `10`, `12`, `16`) |
| `[` | shift left | `]` | shift right |
| `(` | pad left | `)` | pad right |
| `F` | is fraction? | `I` | is integer? |
| `N` | is negative? | `Z` | is zero? |

| Key | Circular | Key | Inverse | Key | Hyperbolic | Key | Inverse |
|---|---|---|---|---|---|---|---|
| `s` | sine | `x` | arcsine | `S` | sinh | `j` | arcsinh |
| `c` | cosine | `y` | arccosine | `C` | cosh | `h` | arccosh |
| `t` | tangent | `X` | arctangent | `T` | tanh | `g` | arctanh |
| `o` | cotangent | `v` | arccotangent | `O` | coth | `W` | arccoth |
| `e` | secant | `V` | arcsecant | `E` | sech | `u` | arcsech |
| `k` | cosecant | `w` | arccosecant | `K` | csch | `U` | arccsch |

The bench reaches only part of the library; `lz_factorial`, `lz_combinations`,
`lz_permutations`, `lz_gamma`, `lz_log_gamma`, the roots, `lz_reciprocal`,
`lz_agm` and the means of a list are callable from code but have no key bound
to them.

## Building

```bash
cmake -G Ninja -S . -B build
ninja -C build
./build/Addaptic
```

Every target builds with `-Wall -Wextra`, or `/W4 /permissive-` under MSVC, and
nothing — engine, tests or bench — draws a warning from clang, GCC or MSVC. It
builds and passes every gate on macOS with clang, on Linux with GCC, and with
Visual Studio on x64 and ARM64, every gate's output identical on all four.
Configuring with `-DADDAPTIC_ASAN=ON` adds AddressSanitizer to each — under
MSVC too — and makes the test run abort on the first memory error: worth doing
after any change to the digit arrays, since the accuracy suite measures digits
and a stray read or write costs none. Configuring with `-DADDAPTIC_TSAN=ON`
adds ThreadSanitizer instead, for the races the threads gate cannot see.

C++17, no external dependencies. Three translation units:

| File | Contents |
|---|---|
| `Addaptic.cpp` | the console bench |
| `LZMath.{h,cpp}` | the math API and its high-precision constants |
| `LZNumber.{h,cpp}` | `CLZNumber` — representation, storage, rounding, base conversion |
| `LZAddaptic.h` | build configuration |

The test programs beside them — `AddapticTests.cpp`, `AddapticSweep.cpp`,
`AddapticThreads.cpp` and the audit's `auditrounding.cpp` — are described
under Tests.

## Accuracy

**Every answer is correctly rounded at the working precision**: the value it
would have if worked out exactly and rounded once, in the number's own base and
rounding mode. The audit checks that absolutely, against Python's `decimal`, on
every case of the magnitude sweep — 2,328 answers from 1e-40 to 1e40 across the
functions, all correctly rounded at 36 places, with the 47 poles, tables and
answers too large to hold set aside.

Three things make it so:

- **Each function works past the answer** by guard digits sized for it, and is
  built to be right to about a unit there.
- **When the final rounding lands within a hair of a half**, the function is
  worked out again wider until the answer is clear of it. That reaches a near
  tie as many places past the answer's last digit as the answer itself has —
  36 more, for a 36-place answer. One closer still is rounded as the widest
  attempt stands.
- **A value exactly on a half is decided exactly** wherever one can arise — the
  roots, a power to a fraction, a logarithm that is a fraction — with whole
  numbers rather than another, wider estimate.

A precision counts places after the point, so a small answer — `arcsin(.05)`
among them — has fewer significant digits than one starting at the first
place. A quotient is the exception: below a tenth — a base-th, in another
base — it keeps the precision as significant digits, the noughts after the
point added.

Every constant table in `LZMath.h` carries 50 significant digits, and past
them `e`, π and the logarithm of the base are worked out at the width asked
for and kept.

## Tests

```bash
ctest --test-dir build --output-on-failure    # the suite, the sweep, the threads gate
python3 auditrounding.py build/auditrounding  # the audit, from the repository root
```

Four programs, each a gate:

- **`AddapticTests`**, 2,805 checks, exercises the API directly, so it reaches
  functions the bench has no key for: exact arithmetic, division and its
  rounding in every mode, every function against references, correct rounding
  at near ties and exact ones, every base, the answer limit, and identities such
  as `sin² + cos² = 1` that need no reference at all. Accuracy is scored by
  aligned decimal subtraction — the ratio of the reference to the difference —
  rather than by matching leading digits, so values straddling a power of ten
  are judged on how far apart they are. References are embedded, correctly
  rounded to 50 digits, so the suite needs nothing beyond a compiler.
- **`AddapticSweep`** runs every function across magnitudes, 1e-40 to 1e40 and
  toward each singular point, 2,375 cases, against the digits each achieved
  when last recorded in `AddapticSweepExpect.h` — a regression gate for
  answers at sizes the suite does not reach. `AddapticSweep --calibrate`
  regenerates the expectations.
- **`auditrounding.py`** asks the absolute question: is each sweep answer the
  correctly rounded value? It drives `auditrounding` and works out every
  reference itself in Python's `decimal`, so it needs no other library.
- **`AddapticThreads`** runs every function eight threads at a time, two to
  each of four bases and widths, and requires each answer to match the one a
  single thread gives — with each thread's own operands, with operands two
  threads share, and with division taking the reciprocal path. It also checks
  that a new thread starts at 36 places. Most races leave every answer right on
  most runs, so run it under `-DADDAPTIC_TSAN=ON` too after any change to state
  the engine keeps.

**CI runs them on every push that touches code** — `.github/workflows/ci.yml`:
GCC on Linux and MSVC on Windows run all four and require the sweep's
calibration to reproduce `AddapticSweepExpect.h` exactly; clang on Linux runs
the suite, the sweep, the threads gate and the audit under AddressSanitizer,
UBSan and LeakSanitizer; and any warning fails the build. MSVC's
AddressSanitizer runs the same weekly, and by hand. The checks themselves are
`ci/checks.sh`, which runs against any build directory locally too.

## The write-up

[`docs/every-digit-it-is-given.html`](docs/every-digit-it-is-given.html) is a
survey of arbitrary-precision arithmetic — the design axes an engine chooses
between, how those choices play out in C++, and where this engine falls among
the alternatives. It is the context this README does not give: what the
competition is, and what Addaptic is trading for what.
