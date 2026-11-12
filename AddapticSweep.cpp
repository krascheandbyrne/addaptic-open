/////////////////////////////////////////////////////////////////
//
//	AddapticSweep.cpp
//
//	The magnitude sweep: every function whose answer depends on
//	the scale of its argument, run from 1e-40 to 1e+40 and either
//	side of one.
//
//	Copyright (c) 2006-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

/*
	Why this exists.

	Seven silent wrong answers were found in this engine over two days,
	and every one of them had the same shape: correct for the sizes
	somebody had tested, wrong outside them, returning a plausible
	number with no flag and no error code. The square root of 1e40 came
	back as 5.4e20. The hyperbolics were right only for whole numbers
	below 2*pi. The logarithm hung below 1e-37. The accuracy suite had
	402 tests and not one above 1e30 or below 1e-7.

	This is the gate that shape of fault has to get past now. It is not
	an accuracy suite -- AddapticTests is that, and it is the finer
	instrument. This one only asks whether each function still answers
	as well as it did at every size, which is the question nobody was
	asking.

	Four outcomes, and they must not be conflated:

	  sound            the answer carries the full working width
	  floor            short of that, but no shorter than the size of
	                   the answer or of the argument allows. A value of
	                   about 1e-k held to 36 places carries about 36 - k
	                   significant digits whatever the algorithm does,
	                   and a reciprocal is no better than the small
	                   thing it divides by
	  constant width   short of everything the representation explains.
	                   **Empty as of the kTwoPi widening.** It held
	                   thirty cases, the circular functions at 1e20 and
	                   above running out of 2*pi; the constant now
	                   carries 120 places and they are all sound. The
	                   bucket stays because it is the one that finds
	                   this kind of fault
	  pole             a refusal, and the right answer. **Empty since
	                   2026-09-25.** It held coth and csch at +-1e-40,
	                   a division by a sinh rounded to nought; an
	                   answer under 0.1 that is a quotient now keeps 36
	                   significant digits, and the two work wider by
	                   twice their answer's size, so they answer 1e40,
	                   correctly rounded. Their only pole is nought,
	                   pinned in AddapticTests with the answers at
	                   +-1e-45

	Only the first is worth improving. Recording the others is what
	stops somebody "fixing" them.

	Run it with --calibrate to print a fresh expectations table after a
	deliberate change; paste that over kSweepExpectations below.

	**Read what --calibrate writes to stderr before you trust the
	table.** Calibration records what the engine does, not what is
	right, so it makes a regression gate and not a correctness one: it
	can tell you tomorrow's answer got worse, never that today's is
	already wrong. That blindness had a cost -- ln(10) sat in this
	table as "sound" at 34 digits, three below its neighbours, for as
	long as the logarithm's reduction-band fault existed. The stderr
	list names every case the representation does not account for, and
	is meant to stay short enough to get read. It is **empty** now: the
	thirty it used to hold were the circular functions running out of
	2*pi, and widening kTwoPi to 120 places answered all thirty. Had
	the list existed earlier, ln at .1, 10, 200 and 4000 would have
	been sitting in it, and they are the four arguments the
	reduction-band fault hit. An empty list is the expected state --
	anything appearing in it is new.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "LZNumber.h"
#include "LZMath.h"

#include "AddapticSweepRefs.h"
#include "AddapticSweepExpect.h"


/////////////////////////////////////////////////////////////////
//	Scoring -- the same measure AddapticTests uses, so a number
//	here means what a number there means.
/////////////////////////////////////////////////////////////////

static void Normalize( const char* szValue, int& iSign, int& iExponent, char* szDigits, int iMaxDigits )
{
	const char*	p = szValue;
	int			iIntegerDigits = 0;
	int			iCount = 0;
	int			iLeadingZeros = 0;

	iSign = 1;

	if ( *p == '-' )		{ iSign = -1; p++; }
	else if ( *p == '+' )	{ p++; }

	for ( const char* q = p; *q && (*q != '.'); ++q )
		iIntegerDigits++;

	for ( const char* q = p; *q; ++q )
	{
		if ( *q == '.' )
			continue;

		if ( iCount < (iMaxDigits - 1) )
			szDigits[ iCount++ ] = *q;
	}

	szDigits[ iCount ] = '\0';
	iExponent = iIntegerDigits;

	while ( szDigits[ iLeadingZeros ] == '0' )
		iLeadingZeros++;

	if ( szDigits[ iLeadingZeros ] == '\0' )
	{
		iSign = 0;
		iExponent = 0;
		szDigits[ 0 ] = '\0';
		return;
	}

	if ( iLeadingZeros )
	{
		::memmove( szDigits, szDigits + iLeadingZeros, (iCount - iLeadingZeros) + 1 );
		iExponent -= iLeadingZeros;
	}
}

const int kExactAgreement = 9999;

static int AgreementDigits( const char* szGot, const char* szReference )
{
	static char	szGotDigits[ 512 ];
	static char	szRefDigits[ 512 ];

	int	iGotSign, iGotExp;
	int	iRefSign, iRefExp;

	Normalize( szGot, iGotSign, iGotExp, szGotDigits, 512 );
	Normalize( szReference, iRefSign, iRefExp, szRefDigits, 512 );

	if ( (iGotSign == 0) && (iRefSign == 0) )
		return kExactAgreement;

	if ( (iGotSign == 0) || (iRefSign == 0) || (iGotSign != iRefSign) )
		return 0;

	const int	kWidth = 700;
	int			iAlignExp = (iGotExp > iRefExp) ? iGotExp : iRefExp;

	static char	szGotAligned[ kWidth + 1 ];
	static char	szRefAligned[ kWidth + 1 ];

	::memset( szGotAligned, '0', kWidth );
	::memset( szRefAligned, '0', kWidth );
	szGotAligned[ kWidth ] = '\0';
	szRefAligned[ kWidth ] = '\0';

	for ( int i = 0; szGotDigits[ i ]; i++ )
	{
		int iPos = (iAlignExp - iGotExp) + i;
		if ( (iPos >= 0) && (iPos < kWidth) ) szGotAligned[ iPos ] = szGotDigits[ i ];
	}

	for ( int i = 0; szRefDigits[ i ]; i++ )
	{
		int iPos = (iAlignExp - iRefExp) + i;
		if ( (iPos >= 0) && (iPos < kWidth) ) szRefAligned[ iPos ] = szRefDigits[ i ];
	}

	if ( ::memcmp( szGotAligned, szRefAligned, kWidth ) == 0 )
		return kExactAgreement;

	const char*	pHigh = szGotAligned;
	const char*	pLow = szRefAligned;

	if ( ::memcmp( szGotAligned, szRefAligned, kWidth ) < 0 )
	{
		pHigh = szRefAligned;
		pLow = szGotAligned;
	}

	static char	szDiff[ kWidth + 1 ];
	int			iBorrow = 0;

	szDiff[ kWidth ] = '\0';

	for ( int i = kWidth - 1; i >= 0; i-- )
	{
		int d = (pHigh[ i ] - '0') - (pLow[ i ] - '0') - iBorrow;

		if ( d < 0 )	{ d += 10; iBorrow = 1; }
		else			{ iBorrow = 0; }

		szDiff[ i ] = static_cast<char>('0' + d);
	}

	for ( int i = 0; i < kWidth; i++ )
		if ( szDiff[ i ] != '0' )
			return iRefExp - (iAlignExp - i);

	return kExactAgreement;
}


/////////////////////////////////////////////////////////////////
//	Dispatch
/////////////////////////////////////////////////////////////////

/*
	The means of a list, swept as the mean of 1 and x. The sweep has one
	argument per case, and pairing it with 1 keeps every magnitude on the
	ladder in play: at 1e40 the three answers are about 5e39, 1e20 and 2,
	and at 1e-40 about 0.5, 1e-20 and 2e-40 -- the sizes the geometric
	mean's exponent scaling exists for.

	The arithmetic-geometric mean goes the same way, as agm(1, x): about
	1.7e38 at 1e40 and 0.017 at 1e-40, where its first product used to
	underflow.
*/

static long MeanOfOneAnd( char cKey, CLZNumber* pArgument, CLZNumber* pResult )
{
	CLZNumber	one;

	one.Assign( 1L );

	CLZNumber*	list[] = { &one, pArgument };

	switch ( cKey )
	{
	case 'u':	return lz_amean( list, 2, pResult );
	case 'v':	return lz_gmean( list, 2, pResult );
	case 'q':	return lz_agm( &one, pArgument, pResult );
	default:	return lz_hmean( list, 2, pResult );
	}
}

/*
	The functions of two operands, swept with the second held fixed, so
	that the one argument a case carries drives the size of the answer:
	3 to the x, x to the 2.5, the log of x to base 3, the 2.5th root of
	x, and x to the 200th. The partners keep the general path in play -- a
	root of order 3 would go to lz_root_3 -- and gensweep.py computes
	each reference with the same partner.

	None of these, nor lz_exp_b, lz_log_b, lz_raise_e and lz_tangent,
	was in the sweep or the audit until 2026-09-24, and those that form
	a product at relative precision had never been asked for 36 places
	of a large answer.
*/

static long WithPartner( char cKey, CLZNumber* pArgument, CLZNumber* pResult )
{
	CLZNumber	partner;

	switch ( cKey )
	{
	case 'P':	partner.Assign( 3L );		return lz_exp_n( &partner, pArgument, pResult );
	case 'i':	partner.Assign( 3L );		return lz_log_n( &partner, pArgument, pResult );
	case 'p':	partner.Assign( "2.5" );	return lz_power( pArgument, &partner, pResult );
	case 'r':	partner.Assign( "2.5" );	return lz_root_n( pArgument, &partner, pResult );
	case 'Q':	partner.Assign( "-1.5" );	return lz_arctangent2( pArgument, &partner, pResult );
	default:	partner.Assign( 200L );	return lz_raise_n( pArgument, &partner, pResult );
	}
}

static long Invoke( char cKey, CLZNumber* pArgument, CLZNumber* pResult )
{
	switch ( cKey )
	{
	case 'R':	return lz_root_2( pArgument, pResult );
	case '3':	return lz_root_3( pArgument, pResult );
	case 'L':	return lz_log_e( pArgument, pResult );
	case 'E':	return lz_exp_e( pArgument, pResult );
	case 'V':	return lz_reciprocal( pArgument, pResult );
	case 's':	return lz_sine( pArgument, pResult );
	case 'c':	return lz_cosine( pArgument, pResult );
	case 'S':	return lz_hypersine( pArgument, pResult );
	case 'H':	return lz_hypercosine( pArgument, pResult );
	case 'T':	return lz_hypertangent( pArgument, pResult );
	case 'j':	return lz_archypersine( pArgument, pResult );
	case 'h':	return lz_archypercosine( pArgument, pResult );
	case 'a':	return lz_arctangent( pArgument, pResult );
	case 'A':	return lz_arcsine( pArgument, pResult );
	case 'C':	return lz_arccosine( pArgument, pResult );
	case 'g':	return lz_archypertangent( pArgument, pResult );

	// The reciprocal family...
	case 'o':	return lz_cotangent( pArgument, pResult );
	case 'x':	return lz_secant( pArgument, pResult );
	case 'n':	return lz_cosecant( pArgument, pResult );
	case 'O':	return lz_hypercotangent( pArgument, pResult );
	case 'X':	return lz_hypersecant( pArgument, pResult );
	case 'N':	return lz_hypercosecant( pArgument, pResult );
	case 'k':	return lz_arccotangent( pArgument, pResult );
	case 'y':	return lz_arcsecant( pArgument, pResult );
	case 'm':	return lz_arccosecant( pArgument, pResult );
	case 'K':	return lz_archypercotangent( pArgument, pResult );
	case 'Y':	return lz_archypersecant( pArgument, pResult );
	case 'M':	return lz_archypercosecant( pArgument, pResult );

	// Gamma and its logarithm, which run on their own ladder -- every
	//	negative whole number is a pole and the shared one is full of
	//	them. See GAMMA_LADDER in gensweep.py...
	case 'G':	return lz_gamma( pArgument, pResult );
	case 'F':	return lz_log_gamma( pArgument, pResult );

	// The three means of a list, and the AGM, as the mean of 1 and x...
	case 'u':
	case 'v':
	case 'w':
	case 'q':	return MeanOfOneAnd( cKey, pArgument, pResult );

	// Outside the sweep until 2026-09-24...
	case 'B':	return lz_exp_b( pArgument, pResult );
	case 'l':	return lz_log_b( pArgument, pResult );
	case 'e':	return lz_raise_e( pArgument, pResult );
	case 't':	return lz_tangent( pArgument, pResult );
	case 'P':
	case 'i':
	case 'p':
	case 'r':
	case 'z':
	case 'Q':	return WithPartner( cKey, pArgument, pResult );

	default:	return kFailed;
	}
}


/////////////////////////////////////////////////////////////////
//	Expectations
/////////////////////////////////////////////////////////////////

static const SweepExpectation* FindExpectation( const char* szFunction, const char* szArgument )
{
	for ( long i = 0; i < kSweepExpectationCount; i++ )
	{
		if ( (::strcmp( kSweepExpectations[ i ].function, szFunction ) == 0) &&
			 (::strcmp( kSweepExpectations[ i ].argument, szArgument ) == 0) )
			return &(kSweepExpectations[ i ]);
	}

	return nullptr;
}

/*
	How many zeros a value carries between its radix mark and its first
	significant digit -- which is how many of the 36 places it holds are
	spent before any information arrives.
*/

static int LeadingZeros( const char* szValue )
{
	const char*	pPoint = ::strchr( szValue, '.' );

	if ( pPoint == nullptr )
		return 0;

	int	iZeros = 0;

	for ( const char* q = pPoint + 1; *q; ++q )
	{
		if ( *q != '0' )
			break;

		iZeros++;
	}

	// An integer part that is not itself nought means no zeros are
	//	being spent at all: 1.0004 starts at its first place...
	for ( const char* q = szValue; q < pPoint; ++q )
		if ( (*q >= '1') && (*q <= '9') )
			return 0;

	return iZeros;
}

static const char* KindName( int iKind )
{
	switch ( iKind )
	{
	case kSweepFloor:			return "floor";
	case kSweepConstantWidth:	return "constant width";
	case kSweepPole:			return "pole";
	default:					return "sound";
	}
}


/////////////////////////////////////////////////////////////////
//	The Main Line
/////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
	bool	bCalibrate = (argc > 1) && (::strcmp( argv[ 1 ], "--calibrate" ) == 0);

	lzu_setup();

	long	lPassed = 0;
	long	lFailed = 0;
	long	lUnknown = 0;
	long	lReview = 0;

	if ( bCalibrate )
	{
		printf( "// Generated by AddapticSweep --calibrate. What the engine\n" );
		printf( "//\tachieves today, so that a drop anywhere is a regression.\n\n" );
		printf( "#pragma once\n\n" );
		printf( "const int kSweepSound = 0;\n" );
		printf( "const int kSweepFloor = 1;\n" );
		printf( "const int kSweepConstantWidth = 2;\n" );
		printf( "const int kSweepPole = 3;\n\n" );
		printf( "struct SweepExpectation\n{\n" );
		printf( "\tconst char*\tfunction;\n\tconst char*\targument;\n" );
		printf( "\tint\t\t\tleast;\n\tint\t\t\tkind;\n};\n\n" );
		printf( "static const SweepExpectation kSweepExpectations[] =\n{\n" );

		fprintf( stderr, "// Cases the representation does not account for.\n" );
		fprintf( stderr, "// Read these before trusting the table: calibration\n" );
		fprintf( stderr, "// records what the engine does, not what is right.\n" );
	}
	else
	{
		printf( "\nMagnitude sweep -- %ld cases\n\n", kSweepReferenceCount );
	}

	for ( long i = 0; i < kSweepReferenceCount; i++ )
	{
		const SweepReference&	theCase = kSweepReferences[ i ];

		CLZNumber	a( theCase.argument );
		CLZNumber	r;

		long	lCode = Invoke( theCase.key, &a, &r );

		int		iDigits = 0;
		bool	bRefused = (lCode != kNoError);

		if ( !bRefused )
		{
			r.TrimTrailing();

			char*	szGot = r.GetStringRep();

			iDigits = AgreementDigits( (szGot != nullptr) ? szGot : "", theCase.reference );

			if ( iDigits == kExactAgreement )
				iDigits = 50;

			if ( iDigits < 0 )
				iDigits = 0;

			::free( szGot );
		}

		if ( bCalibrate )
		{
			/*
				Classify. A case that scores well is sound. One that
				does not is either the floor -- the true answer being
				too small for the precision to hold, which shows up as
				no agreement at all on an argument at the bottom of the
				ladder -- or the width of a constant, which is what
				limits the circular functions on a huge argument. Both
				of those are properties of the representation, not
				faults, and both are recorded so nobody tries to fix
				them.
			*/

			int	iKind = kSweepSound;

			/*
				A refusal can be the right answer. A reciprocal of a
				function that vanishes with its argument -- the
				cotangent, cosecant and their hyperbolic pair -- meets
				a divisor of nought once the argument is small enough
				that the sine underneath rounds away: sin(1e-40) is
				nought at 36 places, correctly, and so cot(1e-40) is a
				division by zero and says so. Recorded as the expected
				outcome, so that answering a plausible number there
				later is the regression it would be.
			*/

			if ( bRefused )
			{
				iKind = kSweepPole;
				iDigits = 0;
			}
			else
			{
				/*
					Three ways to score badly without being wrong.

					The floor: the answer itself is too small for the
					precision. A value of about 1e-k held to 36 places
					carries about 36 - k significant digits whatever
					the algorithm does.

					The inherited floor: the answer is large but is a
					reciprocal of something small, and that something
					was itself held to 36 places. cot(1e-8) divides by
					a sine of about 1e-8 known to 36 places, so to 28
					significant digits, and the quotient can be no
					better. Measured: the answer agrees with the
					quotient of the engine's own sine and cosine to 45
					digits, so the division adds nothing -- the loss is
					entirely the sine's. The argument's size is what
					predicts it, not the answer's.

					The width of a constant: what limits the circular
					functions on a huge argument, the reduction modulo
					2*pi needing more digits of 2*pi than the table
					holds. That one **was** a fault of a sort, and was
					fixed by widening kTwoPi -- but only after the
					sweep had spent long enough recording it to show
					the falloff was exactly one digit per digit of the
					argument, which is what identified the constant.

					The other two are not faults, and each is recorded
					so nobody sets out to fix them. Telling them apart
					matters: widening kTwoPi did nothing at all for
					cot(1e-8), which is the one that most looks like it
					wants it.

					This test used to run only when a case scored under
					30 digits, and everything above that was called
					sound without being looked at. That is how ln(10)
					came to be recorded as sound at 34 digits while its
					neighbours managed 36 and 37 -- a three-digit
					shortfall from a real fault in the logarithm's
					reduction band, sitting in the expectations table
					as the expected outcome for as long as it existed.
					A score is measured against what its own answer can
					hold now, whatever that comes to.

					The slack is one digit, not three. Three was enough
					to swallow ln(10)'s shortfall even once the score
					was looked at, which is the second half of the same
					mistake.
				*/

				const int	kSweepSlack = 1;
				const int	kFullWidth = 36;

				int	iAnswerZeros = LeadingZeros( theCase.reference );
				int	iArgumentZeros = LeadingZeros( theCase.argument );

				int	iAffordable = kFullWidth - iAnswerZeros;
				int	iInherited = kFullWidth - iArgumentZeros;

				if ( iDigits >= (kFullWidth - kSweepSlack) )
				{
					// Carrying the full working width. Nothing to
					//	explain, whatever the sizes involved...
					iKind = kSweepSound;
				}
				else if ( (iDigits >= (iAffordable - kSweepSlack)) ||
						  (iDigits >= (iInherited - kSweepSlack)) )
				{
					// Short of full width, but no shorter than the
					//	answer's own size or the argument's allows...
					iKind = kSweepFloor;
				}
				else
				{
					// Short of everything the representation explains.
					//	The interesting bucket, and the one the review
					//	list below exists for...
					iKind = kSweepConstantWidth;
				}
			}

			/*
				Anything not sound goes to stderr as well, where it can
				be read without disturbing the header on stdout.

				Calibration records what the engine achieves, which
				makes it a regression gate and not a correctness one:
				it cannot tell you that today's answer is already
				wrong, only that tomorrow's got worse. The one defence
				against baking in a fault is for somebody to read the
				cases that scored under what they could have, and
				decide. This list is that, and it is short -- a dozen
				or so out of fifteen hundred -- precisely so it gets
				read.
			*/

			if ( iKind == kSweepConstantWidth )
			{
				lReview++;

				fprintf( stderr, "//   %-12s %-44s %3d digits\n",
						 theCase.function, theCase.argument, iDigits );
			}

			printf( "\t{ \"%s\", \"%s\", %d, %d },\n",
					theCase.function, theCase.argument, iDigits, iKind );
			continue;
		}

		const SweepExpectation*	pExpected = FindExpectation( theCase.function, theCase.argument );

		if ( pExpected == nullptr )
		{
			lUnknown++;
			printf( "  UNKNOWN  %-11s %-44s no expectation recorded\n",
					theCase.function, theCase.argument );
			continue;
		}

		if ( bRefused )
		{
			// A refusal is right where the divisor underneath has
			//	rounded to nought, and wrong everywhere else...
			if ( pExpected->kind == kSweepPole )
			{
				lPassed++;
			}
			else
			{
				lFailed++;
				printf( "  FAIL     %-11s %-44s refused with code %ld\n",
						theCase.function, theCase.argument, lCode );
			}
		}
		else if ( pExpected->kind == kSweepPole )
		{
			// ...and an answer where a refusal was expected is the
			//	regression the recorded pole exists to catch...
			lFailed++;
			printf( "  FAIL     %-11s %-44s answered where the pole was expected\n",
					theCase.function, theCase.argument );
		}
		else if ( iDigits < pExpected->least )
		{
			lFailed++;
			printf( "  FAIL     %-11s %-44s %d digits, had %d (%s)\n",
					theCase.function, theCase.argument, iDigits,
					pExpected->least, KindName( pExpected->kind ) );
		}
		else
		{
			lPassed++;
		}
	}

	if ( bCalibrate )
	{
		printf( "};\n\n" );
		printf( "static const long kSweepExpectationCount =\n" );
		printf( "\t\tsizeof(kSweepExpectations) / sizeof(kSweepExpectations[ 0 ]);\n" );

		fprintf( stderr, "// %ld of %ld cases want reading; the rest are sound.\n",
				 lReview, kSweepReferenceCount );

		lzu_takedown();
		return 0;
	}

	printf( "\n-------------------------------------------------\n" );
	printf( "  held %ld, dropped %ld, unrecorded %ld\n\n", lPassed, lFailed, lUnknown );

	lzu_takedown();

	return ((lFailed == 0) && (lUnknown == 0)) ? 0 : 1;
}
