/////////////////////////////////////////////////////////////////
//
//	LZNumber.cpp
//
//	Part of Addaptic
//
//	Copyright (c) 2006-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <new>
#include <vector>

#include "LZNumber.h"
#include "LZMath.h"


/////////////////////////////////////////////////////////////////
//  Initialize static members...
/////////////////////////////////////////////////////////////////

thread_local bool CLZNumber::ms_bRoundAlternate = false;
thread_local long CLZNumber::ms_lInternalWork = 0;
thread_local long CLZNumber::ms_lRecordRoundings = 0;
thread_local double CLZNumber::ms_dLastTail = -1.0;
thread_local long CLZNumber::ms_lLastPlaces = 0;
thread_local long CLZNumber::ms_lConversionsByProducts = 0;

// How far from the half a recorded tail stops reading -- the widest doubt
//	window with room to spare: a tail that far from the half is in no
//	window, and its first digit or two say so. How many digits it reads
//	at most is kRecordedTailDigits, in LZNumber.h...
const double	kRecordedTailReach = 0.01;


/////////////////////////////////////////////////////////////////
//  Globals...
/////////////////////////////////////////////////////////////////

thread_local long glWorkingPrecision = kInitialPrecision;


/////////////////////////////////////////////////////////////////
//  Reading and writing digits. The two directions have to
//  agree, so they live together.
/////////////////////////////////////////////////////////////////

/*
	What a digit is worth in this base, or -1 if the character is not
	a digit at all.

	**The only place that decides this on the way in.** Initialize's
	decode loop used to spell the same mapping out over sixteen cases
	of its own, base twelve's D and E special-cased in both -- the two
	had to be changed together and nothing enforced it. The loop now
	handles the sign and the radix mark and asks this for everything
	else, so there is one answer rather than two that agree.

	The one that catches people: **base twelve spells ten and eleven D
	and E**, not A and B, and the decoder has said so for as long as
	there have been base tests. A and B are taken as ten and eleven
	there too, by the ordinary hex path, so both spellings go in and D
	and E come out. C is refused in base twelve, being worth twelve.
*/

static long lzi_digit_value( char cDigit, long lBase )
{
	if ( (cDigit >= '0') && (cDigit <= '9') )
		return cDigit - '0';

	if ( lBase == 12 )
	{
		if ( (cDigit == 'D') || (cDigit == 'd') )	return 10;
		if ( (cDigit == 'E') || (cDigit == 'e') )	return 11;
	}

	if ( (cDigit >= 'A') && (cDigit <= 'F') )	return (cDigit - 'A') + 10;
	if ( (cDigit >= 'a') && (cDigit <= 'f') )	return (cDigit - 'a') + 10;

	return -1;
}


/*
	Is this text a number, in the base the value is held in?

	Nothing asked this before, and the decoder below is a switch with
	no default case, so anything it did not recognise was skipped in
	silence: "12x34" became 1234. A sign was taken wherever it
	appeared rather than only at the front, so the minus in "1.5E-8"
	made the whole value negative. A second radix mark was accepted
	and simply moved the first. And A through F were read as ten
	through fifteen whatever the base, so "6.02e23" put a digit of
	fourteen into a base-ten number.

	Every one of those returned kNoError with a plausible-looking
	value attached, which is the shape of fault this engine has spent
	its life turning up. They are refused now.

	Exponent notation is refused rather than understood. Reading it is
	a feature and wants designing alongside the rendering side, which
	cannot write it either; what matters here is that "6.02e23" stops
	being silently accepted as something else entirely.

	The rule: an optional leading sign, at most one radix mark, at
	least one digit, every digit valid for the base, and nothing else
	anywhere.
*/

/*
	Does this character mark an exponent, in this base?

	**This is Elzed's convention and Addaptic follows it**, the two
	being destined to merge. Elzed extends scientific notation twice
	over to survive multiple bases:

	  - X marks an exponent as well as E. X is not a digit in any base
	    this engine allows, so it always works; E is a digit in twelve,
	    fifteen and sixteen, and only marks an exponent where it is
	    not one. So 0x1.0X5 is sixteen to the fifth and 1.0E5 is ten to
	    the fifth, and there is a way to write either in any base.

	  - The base raised is the number's own. 0x1.0X5 is 1.0 sixteen
	    times sixteen to the fifth, not ten to the fifth.

	Both were nearly invented differently here before anyone checked
	the sibling's documentation, which would have left the two spelling
	the same idea two ways at integration. Whether E is free is asked
	of the digit mapping rather than answered with a list of bases, so
	the two cannot drift.
*/

static bool lzi_is_exponent_mark( char cMark, long lBase )
{
	if ( (cMark == 'X') || (cMark == 'x') )
		return true;

	if ( (cMark != 'E') && (cMark != 'e') )
		return false;

	long	l_value = lzi_digit_value( 'E', lBase );

	return (l_value < 0) || (l_value >= lBase);
}




/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////

CLZNumber::CLZNumber()
{
    m_bInitialized = false;
	m_bConstant = false;

	Initialize();
}

CLZNumber::CLZNumber( const char *cpTextNumber )
{
    m_bInitialized = false;
	m_bConstant = false;

	Initialize( cpTextNumber );
}

CLZNumber::CLZNumber( long lLongNumber )
{
    char	cpTextNumber[ 64 ];

	m_bInitialized = false;
	m_bConstant = false;

	snprintf( cpTextNumber, 64, "%ld", lLongNumber );
	Initialize( cpTextNumber );
}

CLZNumber::CLZNumber( double dFloatNumber )
{
    char	cpTextNumber[ 64 ];

	m_bInitialized = false;
	m_bConstant = false;

	snprintf( cpTextNumber, 64, "%f", dFloatNumber );
	Initialize( cpTextNumber );
}

CLZNumber::CLZNumber( const CLZNumber& sourceNumber )
{
    m_bInitialized = true;
	m_bConstant = sourceNumber.m_bConstant;

	m_lFlags = sourceNumber.m_lFlags;
	m_lSign = sourceNumber.m_lSign;
	m_lBase = sourceNumber.m_lBase;
	m_lMaxDigit = sourceNumber.m_lMaxDigit;
	m_lHalfDigit = sourceNumber.m_lHalfDigit;
	m_lAdjust = sourceNumber.m_lAdjust;
	m_lPrecision = sourceNumber.m_lPrecision;
	m_rmMode = sourceNumber.m_rmMode;
	m_rdDigit = sourceNumber.m_rdDigit;

	m_lDigitsArrayLength = sourceNumber.m_lDigitsArrayLength;
	m_lDigitsLength = sourceNumber.m_lDigitsLength;
	m_lRadix = sourceNumber.m_lRadix;
	m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof( long )) );

	::memcpy( m_lpDigits, sourceNumber.m_lpDigits, m_lDigitsArrayLength * sizeof(long) );
}

CLZNumber::~CLZNumber()
{
	::free( m_lpDigits );
	m_lpDigits = nullptr;
}


long CLZNumber::AlignStorage( unsigned long lDigitsLength, unsigned long lExponent, bool bPreserve )
{
	long*	lp_newStorage = nullptr;
	long	l_digitsAdjust = 0;
	long	l_exponentAdjust = 0;

	// Determine increase in whole digits...
	if ( static_cast<long>(lDigitsLength - lExponent) > (m_lDigitsLength - m_lRadix) )
		l_digitsAdjust = static_cast<long>(lDigitsLength - lExponent) - (m_lDigitsLength - m_lRadix);

	// Determine increase in exponent size...
	if ( static_cast<long>(lExponent) > m_lRadix )
		l_exponentAdjust = lExponent - m_lRadix;

	if ( l_digitsAdjust || l_exponentAdjust )
	{
		// Calculate new size...
		lDigitsLength = m_lDigitsLength + l_digitsAdjust + l_exponentAdjust;

		if ( bPreserve )
		{
			// Re-allocate as required...
			if ( static_cast<long>(lDigitsLength) > m_lDigitsArrayLength )
			{
				// Create new storage & copy in digits...
				m_lDigitsArrayLength = lDigitsLength;
				lp_newStorage = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof(long)) );
				::memcpy( lp_newStorage + l_exponentAdjust, m_lpDigits, m_lDigitsLength * sizeof(long) );

				// Free the old array...
				::free( m_lpDigits );

				// Save the new info...
				m_lpDigits = lp_newStorage;
			}

			else
			{
				::memmove( m_lpDigits + l_exponentAdjust, m_lpDigits, m_lDigitsLength * sizeof(long) );
				::memset( m_lpDigits, 0, l_exponentAdjust * sizeof(long) );
			}
		}

		else
		{
			// Re-allocate as required...
			if ( static_cast<long>(lDigitsLength) > m_lDigitsArrayLength )
			{
				m_lDigitsArrayLength = lDigitsLength;

				// Free the old array...
				::free( m_lpDigits );

				// Create new storage...
				m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof(long)) );
			}

			else
				::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );
		}

		// Save the new info...
		m_lDigitsLength = lDigitsLength;
		m_lRadix += l_exponentAdjust;
		m_lAdjust = 0;
	}

	return kNoError;
}


long CLZNumber::Assign( const char *cpTextNumber )
{
	Initialize( cpTextNumber );

	/*
		Initialize cannot report -- the constructors call it too -- so
		it leaves a not-a-number behind on text it could not read, and
		this turns that into a code for the caller. Before this,
		Assign returned kNoError for anything at all, including
		"12x34" and "6.02e23".
	*/

	if ( FlagIsOn( kfQNaN ) )
		return kFailed;

	return kNoError;
}


/*
	The text these two build is decimal, so it has to be read as
	decimal and converted -- not read in whatever base the value
	happens to be in.

	Reading it in the value's own base put a 2 into a base-two number.
	Before the text was validated that stored a digit the base does
	not have and the arithmetic went quietly wrong; after it, the
	assignment failed and left a not-a-number. Either way lz_log_e of
	a base-two value was lost, because it builds its own constants
	with Assign( 2L ) -- it hung before the validation and answered
	nought after it.

	The conversion only runs when the base is not ten, so the ordinary
	call pays a comparison.
*/

long CLZNumber::Assign( long lLongNumber )
{
	char	cpTextNumber[ 64 ];

	snprintf( cpTextNumber, 64, "%ld", lLongNumber );

	long	l_wantedBase = m_bInitialized ? m_lBase : 10;

	m_lBase = 10;
	m_lMaxDigit = 9;
	m_lHalfDigit = 5;

	Initialize( cpTextNumber );

	if ( l_wantedBase != 10 )
		SetBase( l_wantedBase, true );

	// A number this wide cannot fail to read, but the flag is the
	//	carrier and the caller is owed the same answer either way...
	if ( FlagIsOn( kfQNaN ) )
		return kFailed;

	return kNoError;
}


/*
	**Every whole digit is kept, and an infinity is one** -- since
	2026-10-15, item 13 of that day's review. The text was "%f" into 64
	bytes, which cut a double past 63 characters short with no error and
	no flag -- 1e308 assigned about 1e62 -- and wrote an infinity as
	"inf", which reads as malformed and came back NaN. The text is now
	as long as "%f" needs, and an infinity or a NaN is flagged as one; an
	infinity here carries no sign, as none in this engine does. "%f"
	still writes six places after the point, a known limit:
	nothing in the engine assigns a double, so that bites only a host.
*/

long CLZNumber::Assign( double dFloatNumber )
{
	long	l_wantedBase = m_bInitialized ? m_lBase : 10;

	m_lBase = 10;
	m_lMaxDigit = 9;
	m_lHalfDigit = 5;

	if ( std::isnan( dFloatNumber ) || std::isinf( dFloatNumber ) )
	{
		Initialize( "0" );

		if ( l_wantedBase != 10 )
			SetBase( l_wantedBase, true );

		SetFlag( std::isnan( dFloatNumber ) ? kfQNaN : kfInfinity );

		return std::isnan( dFloatNumber ) ? kFailed : kNoError;
	}

	int					i_length = snprintf( nullptr, 0, "%f", dFloatNumber );
	std::vector<char>	v_text( static_cast<size_t>( (i_length > 0) ? i_length : 0 ) + 1, '\0' );

	snprintf( v_text.data(), v_text.size(), "%f", dFloatNumber );

	Initialize( v_text.data() );

	if ( l_wantedBase != 10 )
		SetBase( l_wantedBase, true );

	if ( FlagIsOn( kfQNaN ) )
		return kFailed;

	return kNoError;
}


long CLZNumber::Assign( CLZNumber* p_sourceNumber )
{
	if ( p_sourceNumber != this )
	{
		m_bInitialized = p_sourceNumber->m_bInitialized;

		m_lFlags = p_sourceNumber->m_lFlags;
		m_lSign = p_sourceNumber->m_lSign;
		m_lBase = p_sourceNumber->m_lBase;
		m_lMaxDigit = p_sourceNumber->m_lMaxDigit;
		m_lHalfDigit = p_sourceNumber->m_lHalfDigit;
		m_lAdjust = p_sourceNumber->m_lAdjust;
		m_lPrecision = p_sourceNumber->m_lPrecision;
		m_rmMode = p_sourceNumber->m_rmMode;
		m_rdDigit = p_sourceNumber->m_rdDigit;

		long	l_oldLength = m_lDigitsLength;

		m_lDigitsLength = p_sourceNumber->m_lDigitsLength;
		m_lRadix = p_sourceNumber->m_lRadix;

		/*
			**The digits are copied, not the array** -- since 2026-10-15,
			performance item 8 of that day's review. It cleared the whole
			array and copied the source's whole array, slack included:
			9-12% of a wide logarithm. Everything past a number's length is
			nought -- a carry writes into the digit past it on that
			understanding -- so copying the source's digits and clearing
			this one's past them leaves the same array...
		*/
		if ( m_lDigitsArrayLength < p_sourceNumber->m_lDigitsArrayLength )
		{
			m_lDigitsArrayLength = p_sourceNumber->m_lDigitsArrayLength;

			::free( m_lpDigits );

			m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof( long )) );
		}
		else if ( l_oldLength > m_lDigitsLength )
		{
			long	l_clear = ((l_oldLength < m_lDigitsArrayLength) ? l_oldLength : m_lDigitsArrayLength) - m_lDigitsLength;

			::memset( m_lpDigits + m_lDigitsLength, 0, l_clear * sizeof(long) );
		}

		// Then copy source to target...
		::memcpy( m_lpDigits, p_sourceNumber->m_lpDigits, m_lDigitsLength * sizeof(long) );

		CheckForZero();
	}

	return kNoError;
}


/*
	Re-express the digits, which are held in lOldBase, in m_lBase.
	SetBase has already moved m_lBase and its digit limits on, so the
	base being left behind has to be passed in.

	The integer part goes by short division: divide the value by the
	new base doing the arithmetic in the old, and each remainder is
	the next digit of the answer, least significant first, which is
	the order the array wants anyway. That is exact.

	The fraction goes the other way, by repeated multiplication:
	multiply by the new base in the old, and the whole part carried
	out is the next digit, this time most significant first. That is
	not exact, and cannot be -- a fraction that terminates in one base
	need not terminate in another, a tenth being the standard example
	against base two. So it needs a budget and a rounding rule:

	  - The budget is the working precision scaled by the ratio of the
	    logs, because precision here counts places after the radix and
	    a place is worth less in a smaller base. Thirty-six decimal
	    places is about a hundred and twenty binary ones. m_lPrecision
	    moves with it, or the first rounding would undo the work.

	  - The rounding is half-even, as everywhere else in the engine.
	    Truncating would pull every converted value toward zero, which
	    is the bias the default rounding mode was changed to avoid.

	Neither part needs anything from LZMath -- below a floor. Both are
	quadratic, and past it the conversion goes by products instead; see
	BaseConvertByProducts below.
*/

/*
	Log of each base, to base two, scaled by a thousand, so the budget
	can be worked out without reaching for floating point.

	**Every base SetBase accepts has one.** Until 2026-09-24 only 2, 8,
	10, 12 and 16 did, and a zero here makes BaseConvertValue return
	without converting -- after SetBase has already moved m_lBase on.
	So 21.5 put into base seven kept the digits 2, 1 and 5 under a new
	label, and every mixed-base operation and every function reaching
	a base-ten constant answered a wrong value in the other ten bases,
	without a code or a flag. The conversion itself never depended on
	the base; only this table did.
*/

static long BaseLogScaled( long lBase )
{
	switch ( lBase )
	{
	case 2:		return 1000;
	case 3:		return 1585;
	case 4:		return 2000;
	case 5:		return 2322;
	case 6:		return 2585;
	case 7:		return 2807;
	case 8:		return 3000;
	case 9:		return 3170;
	case 10:	return 3322;
	case 11:	return 3459;
	case 12:	return 3585;
	case 13:	return 3700;
	case 14:	return 3807;
	case 15:	return 3907;
	case 16:	return 4000;
	default:	return 0;
	}
}

/*
	How many places of lToBase are worth lPlaces places of lFromBase,
	rounded up: the ratio of the logs, which is how a conversion sizes
	its budget and moves the precision with it. Thirty-six decimal
	places are 120 binary ones and 30 hexadecimal.

	**Anything the engine sizes in places means a different number of
	digits in each base**, so a count that was measured in decimal --
	a guard, a cap, a threshold -- has to go through this before it is
	applied to a number in another base. lz_gamma's widening cap is
	the case that showed: thirty places were thirty digits in base ten
	and nine in base two.
*/

long CLZNumber::EquivalentPlaces( long lPlaces, long lFromBase, long lToBase )
{
	long	l_fromLog = BaseLogScaled( lFromBase );
	long	l_toLog = BaseLogScaled( lToBase );

	if ( (l_fromLog == 0) || (l_toLog == 0) || (lFromBase == lToBase) )
		return lPlaces;

	return ((lPlaces * l_fromLog) + l_toLog - 1) / l_toLog;
}


/*
	A fraction's digits, least significant first, against one half:
	-1 below, 0 equal, 1 above. In an even base a half is the digit
	base / 2 and noughts after it. In an odd base it is the middle
	digit repeated forever, so a finite fraction is never equal to it
	and the first digit that is not the middle one decides; a run of
	nothing but middle digits falls short.
*/

static long FractionAgainstHalf( const long* lpDigits, long lLength, long lBase )
{
	long	l_half = lBase / 2;

	if ( (lBase % 2) == 0 )
	{
		if ( lLength <= 0 )
			return -1;

		if ( lpDigits[ lLength - 1 ] != l_half )
			return (lpDigits[ lLength - 1 ] > l_half) ? 1 : -1;

		for ( long i = 0; i < (lLength - 1); i++ )
		{
			if ( lpDigits[ i ] != 0 )
				return 1;
		}

		return 0;
	}

	for ( long i = lLength - 1; i >= 0; i-- )
	{
		if ( lpDigits[ i ] != l_half )
			return (lpDigits[ i ] > l_half) ? 1 : -1;
	}

	return -1;
}


/*
	**Past a floor, a conversion is worked by products** -- since
	2026-10-16, performance item 8(g) of the review of 2026-10-15. Both
	loops below are quadratic: SetBase from ten to two took 1,362 ms at
	10,000 places, half a whole part and the rest fraction, and 122 ms at
	3,000. By products it takes 2.8 ms and 0.62.

	**The value is a whole number over a power of its base**, W / b^m,
	m the places after the point, so what the loops write -- the value at
	a budget of N places of the new base B, rounded half-even -- is the
	whole number nearest W B^N / b^m. That is worked in the old base,
	where dividing by b^m is a shift:

		W B^N = Z b^m + R,	0 <= R < b^m

	Z is the answer before rounding, and R against half of b^m is exactly
	what the loops' guard digit and the leftover past it decide between
	them, in either parity of either base: 2R above b^m rounds up, below
	it down, and equal -- which only an even old base can write -- is the
	tie, settled by the kept digit's parity as the loops settle it. Then Z
	is written in B by halves: the high half and the low converted alone,
	put together with one product by a power of b worked in B, the powers
	made by squaring as the halves want them. Every product goes through
	lz_multiply's own dispatch, transform and Karatsuba included, so the
	whole is some log n products in place of n^2 digit steps.

	**What it writes is the loops' answer to the digit**: the same trailing
	noughts trimmed, a whole part of nought kept as the one digit "0" when
	the value had a whole part and as none when it had not, the same
	precision afterwards. Proved by a census of 78,200 conversions against
	the loops, every pair of bases, 4,368 exact ties, carries into the
	whole part, and every value against an exact reference.

	**A fraction that ends in the new base is worked at the places it ends
	within** rather than the budget -- a hex fraction of k places within 4k
	binary ones -- since the loops stop there too, and B^N past it is all
	cost: 100 hex places at a budget of 109,000 binary ones took 3.8 ms
	by products against 0.14 by the loops before this.

	**The floor is the loops' work against the digits the products write**,
	measured on 16,800 shapes over all 210 pairs of bases: by products
	where the loops would take 5,000 digit steps or more and at least four
	per digit written, times its bit length. Below either, the loops are
	quicker: a value of a few dozen digits pays for the powers and the
	vectors, and a short fraction at a wide budget -- three places at
	4,096 -- costs the loops 3N steps where the products write all N
	digits. Checked on 3,000 fresh shapes: every one the rule sends by
	products is quicker, the worst at 0.70 of the loops' time.
*/

static long gConversionFloor = CLZNumber::kConversionFloorMeasured;

constexpr int64_t	kConversionWorkFloor = 5000;
constexpr int64_t	kConversionWorkPerDigit = 4;

long CLZNumber::SetConversionFloor( long lWork )
{
	long	l_previous = gConversionFloor;

	gConversionFloor = lWork;

	return l_previous;
}

long CLZNumber::ConversionsByProducts()
{
	return ms_lConversionsByProducts;
}

// A digit array's length less its leading noughts...
static long DigitsLength( const long* lpDigits, long lLength )
{
	while ( (lLength > 0) && (lpDigits[ lLength - 1 ] == 0) )
		lLength--;

	return lLength;
}

// Two digit arrays' product, by lz_multiply's dispatch...
static std::vector<long> DigitsProduct( const std::vector<long>& vLeft, const std::vector<long>& vRight, long lBase )
{
	std::vector<long>	v_product;

	if ( vLeft.empty() || vRight.empty() )
		return v_product;

	v_product.assign( vLeft.size() + vRight.size(), 0 );

	lzi_multiply_digit_arrays( vLeft.data(), static_cast<long>( vLeft.size() ),
							   vRight.data(), static_cast<long>( vRight.size() ), v_product.data(), lBase );

	v_product.resize( DigitsLength( v_product.data(), static_cast<long>( v_product.size() ) ) );

	return v_product;
}

// A digit array times a small whole number...
static void DigitsTimesSmall( std::vector<long>& vDigits, long lFactor, long lBase )
{
	int64_t	l_carry = 0;

	for ( long& l_digit : vDigits )
	{
		int64_t	l_current = (static_cast<int64_t>( l_digit ) * lFactor) + l_carry;

		l_digit = static_cast<long>( l_current % lBase );
		l_carry = l_current / lBase;
	}

	while ( l_carry > 0 )
	{
		vDigits.push_back( static_cast<long>( l_carry % lBase ) );
		l_carry /= lBase;
	}
}

// lpDigits, written in lFrom, written in lTo by Horner's rule, as many
//	digits a step as keep the step's multiplier under 2^31...
static std::vector<long> DigitsInBaseByDigits( const long* lpDigits, long lLength, long lFrom, long lTo )
{
	int64_t	l_most = 1;
	long	l_perStep = 0;

	while ( (l_most * lFrom) < (static_cast<int64_t>( 1 ) << 31) )
	{
		l_most *= lFrom;
		l_perStep++;
	}

	std::vector<long>	v_out;

	v_out.reserve( static_cast<size_t>( CLZNumber::EquivalentPlaces( lLength, lFrom, lTo ) ) + 2 );

	long	l_next = lLength;

	while ( l_next > 0 )
	{
		long	l_take = (l_next < l_perStep) ? l_next : l_perStep;
		int64_t	l_multiplier = 1;
		int64_t	l_carry = 0;

		for ( long k = 0; k < l_take; k++ )
		{
			l_multiplier *= lFrom;
			l_carry = (l_carry * lFrom) + lpDigits[ --l_next ];
		}

		for ( long& l_digit : v_out )
		{
			int64_t	l_current = (static_cast<int64_t>( l_digit ) * l_multiplier) + l_carry;

			l_digit = static_cast<long>( l_current % lTo );
			l_carry = l_current / lTo;
		}

		while ( l_carry > 0 )
		{
			v_out.push_back( static_cast<long>( l_carry % lTo ) );
			l_carry /= lTo;
		}
	}

	return v_out;
}

// lFrom to the power of lChunk times 2^level, written in lTo, made as wanted...
struct SConversionPowers
{
	long							lFrom;
	long							lTo;
	long							lChunk;
	std::vector<std::vector<long>>	vPowers;
};

// lpDigits, written in lFrom, written in lTo: the high half and the low
//	converted alone, and put together by one product...
static std::vector<long> DigitsInBase( const long* lpDigits, long lLength, SConversionPowers& powers )
{
	lLength = DigitsLength( lpDigits, lLength );

	if ( lLength <= powers.lChunk )
		return DigitsInBaseByDigits( lpDigits, lLength, powers.lFrom, powers.lTo );

	size_t	l_level = 0;
	long	l_split = powers.lChunk;

	while ( (2 * l_split) < lLength )
	{
		l_split *= 2;
		l_level++;
	}

	while ( powers.vPowers.size() <= l_level )
		powers.vPowers.push_back( DigitsProduct( powers.vPowers.back(), powers.vPowers.back(), powers.lTo ) );

	std::vector<long>	v_low = DigitsInBase( lpDigits, l_split, powers );
	std::vector<long>	v_out = DigitsProduct( DigitsInBase( lpDigits + l_split, lLength - l_split, powers ),
											   powers.vPowers[ l_level ], powers.lTo );

	if ( v_out.size() < v_low.size() )
		v_out.resize( v_low.size(), 0 );

	long	l_carry = 0;

	for ( size_t i = 0; (i < v_low.size()) || (l_carry != 0); i++ )
	{
		if ( i == v_out.size() )
			v_out.push_back( 0 );

		long	l_current = v_out[ i ] + ((i < v_low.size()) ? v_low[ i ] : 0) + l_carry;

		l_carry = (l_current >= powers.lTo) ? 1 : 0;
		v_out[ i ] = l_current - (l_carry * powers.lTo);
	}

	v_out.resize( DigitsLength( v_out.data(), static_cast<long>( v_out.size() ) ) );

	return v_out;
}

// lValue to the power lExponent, written in lBase...
static std::vector<long> DigitsPower( long lValue, long lExponent, long lBase )
{
	std::vector<long>	v_power( 1, 1 );
	int					l_bit = 62;

	while ( (l_bit >= 0) && (((static_cast<int64_t>( lExponent ) >> l_bit) & 1) == 0) )
		l_bit--;

	for ( ; l_bit >= 0; l_bit-- )
	{
		v_power = DigitsProduct( v_power, v_power, lBase );

		if ( ((static_cast<int64_t>( lExponent ) >> l_bit) & 1) != 0 )
			DigitsTimesSmall( v_power, lValue, lBase );
	}

	return v_power;
}

/*
	Places of lToBase within which a fraction of lPlaces places of lFromBase
	ends, when every prime of lFromBase divides lToBase -- a binary fraction
	of k places ends within k decimal ones, a hex one within 4k binary --
	and -1 when some fraction of lFromBase never ends there.
*/
static long PlacesWithinWhichItEnds( long lPlaces, long lFromBase, long lToBase )
{
	long	l_most = 0;

	for ( long l_prime = 2; l_prime <= lFromBase; l_prime++ )
	{
		bool	b_prime = (lFromBase % l_prime) == 0;

		for ( long l_factor = 2; b_prime && ((l_factor * l_factor) <= l_prime); l_factor++ )
			b_prime = (l_prime % l_factor) != 0;

		if ( !b_prime )
			continue;

		long	l_from = 0;
		long	l_to = 0;

		for ( long l_left = lFromBase; (l_left % l_prime) == 0; l_left /= l_prime )
			l_from++;

		for ( long l_left = lToBase; (l_left % l_prime) == 0; l_left /= l_prime )
			l_to++;

		if ( l_to == 0 )
			return -1;

		long	l_places = ((lPlaces * l_from) + l_to - 1) / l_to;

		if ( l_places > l_most )
			l_most = l_places;
	}

	return l_most;
}

/*
	The value times m_lBase^budget, whole, rounded to a whole number by the
	same rule as the loops below, and written out as they would write it.
*/
bool CLZNumber::BaseConvertByProducts( long lOldBase, long lFractionBudget, long lIntegerLength )
{
	try
	{
		long				l_length = DigitsLength( m_lpDigits, m_lDigitsLength );
		long				l_budget = (m_lRadix > 0) ? lFractionBudget : 0;
		long				l_ends = PlacesWithinWhichItEnds( m_lRadix, lOldBase, m_lBase );

		// A fraction that ends sooner is exact there, and the budget past
		//	it would only add noughts...
		if ( (l_ends >= 0) && (l_ends < l_budget) )
			l_budget = l_ends;

		long				l_against = -1;
		std::vector<long>	v_scaled;

		if ( m_lRadix > 0 )
		{
			std::vector<long>	v_whole( m_lpDigits, m_lpDigits + l_length );
			std::vector<long>	v_product = DigitsProduct( v_whole, DigitsPower( m_lBase, l_budget, lOldBase ), lOldBase );

			if ( static_cast<long>( v_product.size() ) < m_lRadix )
				v_product.resize( m_lRadix, 0 );

			l_against = FractionAgainstHalf( v_product.data(), m_lRadix, lOldBase );
			v_scaled.assign( v_product.begin() + m_lRadix, v_product.end() );
		}
		else
			v_scaled.assign( m_lpDigits, m_lpDigits + l_length );

		SConversionPowers	powers{ lOldBase, m_lBase, 64, {} };
		std::vector<long>	v_one( powers.lChunk + 1, 0 );

		v_one[ powers.lChunk ] = 1;
		powers.vPowers.push_back( DigitsInBaseByDigits( v_one.data(), powers.lChunk + 1, lOldBase, m_lBase ) );

		std::vector<long>	v_new = DigitsInBase( v_scaled.data(), static_cast<long>( v_scaled.size() ), powers );

		if ( (l_against > 0) || ((l_against == 0) && !v_new.empty() && ((v_new[ 0 ] % 2) != 0)) )
		{
			size_t	i = 0;

			while ( (i < v_new.size()) && (v_new[ i ] == (m_lBase - 1)) )
				v_new[ i++ ] = 0;

			if ( i == v_new.size() )
				v_new.push_back( 1 );
			else
				v_new[ i ]++;
		}

		long	l_newDigits = static_cast<long>( v_new.size() );
		long	l_fracLength = 0;

		for ( long i = 0; (i < l_budget) && (i < l_newDigits); i++ )
		{
			if ( v_new[ i ] != 0 )
			{
				l_fracLength = l_budget - i;
				break;
			}
		}

		long	l_intLength = (l_newDigits > l_budget) ? (l_newDigits - l_budget) : ((lIntegerLength > 0) ? 1 : 0);
		long	l_newLength = l_fracLength + l_intLength;

		if ( l_newLength > m_lDigitsArrayLength )
			SetDigitsArraySize( l_newLength + 2, kNoPreserve );

		::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );

		for ( long i = 0; i < l_newLength; i++ )
		{
			long	l_from = (l_budget - l_fracLength) + i;

			m_lpDigits[ i ] = (l_from < l_newDigits) ? v_new[ l_from ] : 0;
		}

		m_lDigitsLength = l_newLength;
		m_lRadix = l_fracLength;
	}
	catch ( const std::bad_alloc& )
	{
		return false;
	}

	ms_lConversionsByProducts++;

	return true;
}

// The floor: the digit steps the loops would take, against the digits the
//	products would write -- see above...
static bool ConversionByProductsWanted( long lOldBase, long lNewBase, long lIntegerLength, long lFractionLength, long lFractionBudget )
{
	long	l_ends = PlacesWithinWhichItEnds( lFractionLength, lOldBase, lNewBase );
	long	l_steps = ((l_ends >= 0) && (l_ends < (lFractionBudget + 1))) ? l_ends : (lFractionBudget + 1);
	int64_t	l_work = (static_cast<int64_t>( lIntegerLength ) * CLZNumber::EquivalentPlaces( lIntegerLength, lOldBase, lNewBase )) +
					 (static_cast<int64_t>( l_steps ) * lFractionLength);

	if ( gConversionFloor != CLZNumber::kConversionFloorMeasured )
		return l_work >= gConversionFloor;

	// The products' cost goes with the digits they write, in the new base...
	long	l_budget = ((l_ends >= 0) && (l_ends < lFractionBudget)) ? l_ends : lFractionBudget;
	int64_t	l_written = CLZNumber::EquivalentPlaces( lIntegerLength, lOldBase, lNewBase ) + ((lFractionLength > 0) ? l_budget : 0);
	int64_t	l_bits = 0;

	for ( int64_t l_left = l_written + 1; l_left > 0; l_left >>= 1 )
		l_bits++;

	return (l_work >= kConversionWorkFloor) && (l_work >= (kConversionWorkPerDigit * l_written * l_bits));
}


void CLZNumber::BaseConvertValue( long lOldBase )
{
	long	l_oldLog = BaseLogScaled( lOldBase );
	long	l_newLog = BaseLogScaled( m_lBase );

	if ( (lOldBase == m_lBase) || (l_oldLog == 0) || (l_newLog == 0) )
		return;

	CheckForZero();

	if ( IsZero() )
		return;

	long	l_integerLength = m_lDigitsLength - m_lRadix;
	long	l_fractionLength = m_lRadix;

	// Dividing by as little as two lengthens the integer part by at
	//	most the ratio of the logs, which is four across these bases...
	long	l_intLimit = (l_integerLength * 4) + 8;

	// The fraction gets the working precision, scaled the same way...
	long	l_fractionBudget = 0;

	if ( l_fractionLength > 0 )
		l_fractionBudget = EquivalentPlaces( m_lPrecision, lOldBase, m_lBase );

	if ( (l_integerLength >= 0) && ConversionByProductsWanted( lOldBase, m_lBase, l_integerLength, l_fractionLength, l_fractionBudget ) &&
		 BaseConvertByProducts( lOldBase, l_fractionBudget, l_integerLength ) )
	{
		if ( l_fractionLength > 0 )
			m_lPrecision = EquivalentPlaces( m_lPrecision, lOldBase, m_lBase );

		CheckForZero();

		return;
	}

	long*	lp_work = static_cast<long*>(::calloc( m_lDigitsLength + 1, sizeof(long) ));
	long*	lp_int = static_cast<long*>(::calloc( l_intLimit, sizeof(long) ));
	long*	lp_frac = static_cast<long*>(::calloc( l_fractionBudget + 2, sizeof(long) ));

	if ( (lp_work == nullptr) || (lp_int == nullptr) || (lp_frac == nullptr) )
	{
		::free( lp_work );
		::free( lp_int );
		::free( lp_frac );
		return;
	}

	// The integer part, by short division...
	long	l_intLength = 0;
	long	l_workLength = l_integerLength;

	::memcpy( lp_work, m_lpDigits + m_lRadix, l_integerLength * sizeof(long) );

	while ( l_workLength > 0 )
	{
		long	l_remainder = 0;

		// Most significant digit down to least...
		for ( long i = l_workLength - 1; i >= 0; i-- )
		{
			long	l_current = (l_remainder * lOldBase) + lp_work[ i ];

			lp_work[ i ] = l_current / m_lBase;
			l_remainder = l_current % m_lBase;
		}

		lp_int[ l_intLength++ ] = l_remainder;

		// Drop any leading zeros the division just produced...
		while ( (l_workLength > 0) && (lp_work[ l_workLength - 1 ] == 0) )
			l_workLength--;
	}

	// The fraction, by repeated multiplication. One digit past the
	//	budget is taken, to round on...
	long	l_fracLength = 0;
	bool	b_remainderLeft = false;

	if ( l_fractionLength > 0 )
	{
		::memset( lp_work, 0, (m_lDigitsLength + 1) * sizeof(long) );
		::memcpy( lp_work, m_lpDigits, l_fractionLength * sizeof(long) );

		while ( l_fracLength <= l_fractionBudget )
		{
			long	l_carry = 0;

			// Least significant digit up to most...
			for ( long i = 0; i < l_fractionLength; i++ )
			{
				long	l_current = (lp_work[ i ] * m_lBase) + l_carry;

				lp_work[ i ] = l_current % lOldBase;
				l_carry = l_current / lOldBase;
			}

			lp_frac[ l_fracLength++ ] = l_carry;

			// Nothing left to carry out means the fraction ended...
			bool b_allZero = true;

			for ( long i = 0; i < l_fractionLength; i++ )
			{
				if ( lp_work[ i ] != 0 )
				{
					b_allZero = false;
					break;
				}
			}

			if ( b_allZero )
				break;
		}

		for ( long i = 0; i < l_fractionLength; i++ )
		{
			if ( lp_work[ i ] != 0 )
			{
				b_remainderLeft = true;
				break;
			}
		}
	}

	/*
		Round half to even on the guard digit, if one was taken. A
		guard exactly half the base rounds up only when something is
		left beyond it, or when the digit it would land on is odd.

		**In an odd base the guard cannot be exactly half**: its middle
		digit, three in base seven, is under half a place. So there it
		is the leftover beyond the guard that decides -- above half
		rounds up, below rounds down, and exactly half, which a
		leftover can be (a half in base ten is .5), is the tie. The
		leftover is still in the old base, so it is compared with a
		half written in that one.
	*/

	if ( l_fracLength > l_fractionBudget )
	{
		long	l_guard = lp_frac[ --l_fracLength ];
		bool	b_roundUp = false;

		if ( l_guard > m_lHalfDigit )
			b_roundUp = true;

		else if ( l_guard == m_lHalfDigit )
		{
			long	l_beyond = b_remainderLeft ? 1 : 0;

			if ( (m_lBase % 2) != 0 )
				l_beyond = FractionAgainstHalf( lp_work, l_fractionLength, lOldBase );

			if ( l_beyond > 0 )
				b_roundUp = true;
			else if ( l_beyond < 0 )
				b_roundUp = false;
			else if ( l_fracLength > 0 )
				b_roundUp = ((lp_frac[ l_fracLength - 1 ] % 2) != 0);
			else
				b_roundUp = ((lp_int[ 0 ] % 2) != 0);
		}

		// Carry the increment back up through the fraction, and into
		//	the integer part if it runs off the top...
		if ( b_roundUp )
		{
			long	l_index = l_fracLength - 1;

			while ( l_index >= 0 )
			{
				if ( ++lp_frac[ l_index ] <= (m_lBase - 1) )
					break;

				lp_frac[ l_index ] = 0;
				l_index--;
			}

			if ( l_index < 0 )
			{
				long	l_carryIndex = 0;

				while ( l_carryIndex < l_intLength )
				{
					if ( ++lp_int[ l_carryIndex ] <= (m_lBase - 1) )
						break;

					lp_int[ l_carryIndex ] = 0;
					l_carryIndex++;
				}

				if ( l_carryIndex == l_intLength )
					lp_int[ l_intLength++ ] = 1;
			}
		}
	}

	// Trailing zeros in the new fraction carry no information...
	while ( (l_fracLength > 0) && (lp_frac[ l_fracLength - 1 ] == 0) )
		l_fracLength--;

	// Assemble: fraction below the radix, least significant first,
	//	then the integer part above it...
	long	l_newLength = l_fracLength + l_intLength;

	if ( l_newLength > m_lDigitsArrayLength )
		SetDigitsArraySize( l_newLength + 2, kNoPreserve );

	::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );

	for ( long i = 0; i < l_fracLength; i++ )
		m_lpDigits[ i ] = lp_frac[ l_fracLength - 1 - i ];

	::memcpy( m_lpDigits + l_fracLength, lp_int, l_intLength * sizeof(long) );

	m_lDigitsLength = l_newLength;
	m_lRadix = l_fracLength;

	// Precision counts places after the radix, and a place in this
	//	base is not worth what it was in the last one...
	if ( l_fractionLength > 0 )
		m_lPrecision = EquivalentPlaces( m_lPrecision, lOldBase, m_lBase );

	::free( lp_work );
	::free( lp_int );
	::free( lp_frac );

	CheckForZero();
}


void CLZNumber::CheckForZero()
{
	// Decided before anything is written, and then only what changes is
	//	written. Every function that asks whether an operand is nought
	//	asks here, an operand may be shared between threads, and a number
	//	that is already right is then only read...
	bool	b_zero = true;

	for ( long i = m_lDigitsLength-1; i >= 0; i-- )
	{
		if ( m_lpDigits[ i ] != 0 )
		{
			b_zero = false;
			break;
		}
	}

	if ( b_zero != FlagIsOn( kfZero ) )
	{
		if ( b_zero )
			SetFlag( kfZero );
		else
			ClearFlag( kfZero );
	}

	// Set sign to zero...
	if ( b_zero && (m_lSign != kZero) )
		m_lSign = kZero;
}


long CLZNumber::Clear( unsigned long lDigitsLength )
{
	m_bInitialized = false;

	m_lFlags = 0;
	m_lSign = kPositive;
	m_lBase = 10;
	m_lMaxDigit = 9;
	m_lHalfDigit = 5;
	m_lAdjust = 0;
	m_lPrecision = glWorkingPrecision;
	m_rmMode = kDefaultRoundingMode;
	m_rdDigit = RD_digit_precision;
	m_lDigitsLength = lDigitsLength;
	m_lRadix = 0;

	// Resize and/or clear the digits array...
	if ( m_lDigitsArrayLength < static_cast<long>(lDigitsLength) )
	{
		m_lDigitsArrayLength = lDigitsLength;

		::free( m_lpDigits );
		m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof(long)) );
	}
	else
		::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );

	m_bInitialized = true;

	return kNoError;
}


long CLZNumber::ClearFlag( unsigned long lFlag )
{
	m_lFlags &= ( kAllFlagsOn ^ lFlag );

	return kNoError;
}


// Signed: two negatives are ordered the other way round from their
//	magnitudes. It returned the magnitudes' order for any matching signs
//	until 2026-10-15, so -3 compared below -5...
long CLZNumber::Compare( CLZNumber* p_compareNumber )
{
	if ( m_lSign == p_compareNumber->m_lSign )
		return ( (m_lSign == kNegative) ? -CompareABS( p_compareNumber ) : CompareABS( p_compareNumber ) );

	if ( m_lSign == kZero )
		return ( -(p_compareNumber->m_lSign) );

	return m_lSign;
}


/*
	**A value in another base is brought into this one before the
	digits are compared** -- the same move CLZMatchedBase makes for
	lz_compare, and for the same reason: a digit array carries no base,
	so 0x15 and 21 are one number and two different strings of digits.

	Until 2026-09-24 these two compared the raw digits whatever the
	bases. lz_compare had been fixed and these had not, and the engine
	uses them itself against constants it makes in base ten. The one
	that showed was gamma's shift, which runs its argument up until it
	reaches 24: in base two 4.1 is 100.0001..., "100" beat "24", and
	the series was summed at 4 rather than 24 -- five digits where 36
	binary places hold eleven. The exponential's table reach, the
	logarithm's reduction and lz_raise_e's widening probe made the same
	comparison and got away with it.

	The copy is made only when the bases differ, so the ordinary call
	pays one comparison. The argument is never converted in place: it
	is not ours to modify, and callers pass constants. CompareRaw is
	left as it is -- it is the digit comparison, and division uses it
	on two numbers it has already put in one base.
*/

/*
	Kept out of line, and it has to be told. Written inside CompareABS,
	or left for the compiler to inline, the copy it makes gave the
	ordinary same-base call a 224-byte frame and a stack-protector
	check it never needed: 3.1 ns to 3.8 on the comparison and 1.5 per
	cent on lz_add, measured, and none of it with the function kept
	apart.
*/

#if defined( _MSC_VER )
	#define LZN_OUT_OF_LINE		__declspec( noinline )
#else
	#define LZN_OUT_OF_LINE		__attribute__(( noinline ))
#endif

LZN_OUT_OF_LINE static long CompareInThisBase( CLZNumber* p_this, CLZNumber* p_compareNumber )
{
	CLZNumber	matched;

	matched.Assign( p_compareNumber );
	matched.SetBaseForWork( p_this->m_lBase );

	return p_this->CompareABS( &matched );
}

long CLZNumber::CompareABS( CLZNumber* p_compareNumber )
{
	if ( p_compareNumber->m_lBase != m_lBase )
		return CompareInThisBase( this, p_compareNumber );

	long	l_compareResult = kEqualTo;

	/*
		Magnitude is decided by where a value's leading significant
		digit sits, not by how many digits precede the radix mark.
		Counting raw whole digits makes the answer depend on how the
		value was written: "0.5" carries a whole digit that happens to
		be zero, and would otherwise be judged larger than ".84".
	*/

	long	l_thisTop = m_lDigitsLength;
	long	l_passedTop = p_compareNumber->m_lDigitsLength;

	while ( (l_thisTop > 0) && (m_lpDigits[ l_thisTop - 1 ] == 0) )
		l_thisTop--;

	while ( (l_passedTop > 0) && (p_compareNumber->m_lpDigits[ l_passedTop - 1 ] == 0) )
		l_passedTop--;

	long	l_thisDigits = l_thisTop - m_lRadix;
	long	l_passedDigits = l_passedTop - p_compareNumber->m_lRadix;

	// A value with no significant digits is zero, however it is signed...
	bool	b_thisIsZero = (m_lSign == kZero) || (l_thisTop == 0);
	bool	b_passedIsZero = (p_compareNumber->m_lSign == kZero) || (l_passedTop == 0);

	// Start by checking for zeroes...
	if ( b_thisIsZero && b_passedIsZero )
		l_compareResult = kEqualTo;
	else if ( b_thisIsZero )
		l_compareResult = kLessThan;
	else if ( b_passedIsZero )
		l_compareResult = kGreaterThan;

	// Next compare the leading significant digit's magnitude...
	else if ( l_thisDigits > l_passedDigits )
		l_compareResult = kGreaterThan;
	else if ( l_thisDigits < l_passedDigits )
		l_compareResult = kLessThan;

	// If we get here the two values are of the same magnitude, so we
	//   have to compare the digit arrays element by element. Both are
	//   walked down from their leading significant digit, which keeps
	//   them aligned on the same decimal places...
	else
	{
		l_compareResult = kEqualTo;

		for ( long l_offset = 0; ; l_offset++ )
		{
			long	l_thisIndex = (l_thisTop - 1) - l_offset;
			long	l_passedIndex = (l_passedTop - 1) - l_offset;

			if ( (l_thisIndex < 0) && (l_passedIndex < 0) )
				break;

			// A value that has run out of digits contributes zeroes...
			long	l_thisDigit = (l_thisIndex < 0) ? 0 : m_lpDigits[ l_thisIndex ];
			long	l_passedDigit = (l_passedIndex < 0) ? 0 : p_compareNumber->m_lpDigits[ l_passedIndex ];

			if ( l_thisDigit > l_passedDigit )
			{
				l_compareResult = kGreaterThan;
				break;
			}

			if ( l_thisDigit < l_passedDigit )
			{
				l_compareResult = kLessThan;
				break;
			}
		}
	}

	return l_compareResult;
}


long CLZNumber::CompareRaw( CLZNumber* p_compareNumber )
{
	long	l_compareResult = kEqualTo;

	long	l_thisDigits = m_lDigitsLength;
	long	l_passedDigits = p_compareNumber->m_lDigitsLength;
	long	l_compareDigits = 0;

	// Check digit lengths first...
	if ( l_thisDigits > l_passedDigits )
		l_compareResult = kGreaterThan;
	else if ( l_thisDigits < l_passedDigits )
		l_compareResult = kLessThan;

	// If we get here, there are equal numbers of whole digits, so we
	//   have to compare the digit arrays element by element...
	else
	{
		l_thisDigits = m_lDigitsLength;
		l_passedDigits = p_compareNumber->m_lDigitsLength;

		if ( l_thisDigits > l_passedDigits )
		{
			l_compareDigits = l_passedDigits;
			l_compareResult = kGreaterThan;
		}
		else if ( l_thisDigits == l_passedDigits )
		{
			l_compareDigits = l_passedDigits;
			l_compareResult = kEqualTo;
		}
		else if ( l_thisDigits < l_passedDigits )
		{
			l_compareDigits = l_thisDigits;
			l_compareResult = kLessThan;
		}

		while ( l_compareDigits-- > 0 )
		{
			l_thisDigits--;
			l_passedDigits--;

			if ( m_lpDigits[ l_thisDigits ] == p_compareNumber->m_lpDigits[ l_passedDigits ] )
				continue;

			if ( m_lpDigits[ l_thisDigits ] > p_compareNumber->m_lpDigits[ l_passedDigits ] )
			{
				l_compareResult = kGreaterThan;
				break;
			}

			if ( m_lpDigits[ l_thisDigits ] < p_compareNumber->m_lpDigits[ l_passedDigits ] )
			{
				l_compareResult = kLessThan;
				break;
			}
		}
	}

	return l_compareResult;
}


bool CLZNumber::FlagIsOff( unsigned long lFlag )
{
    return ( (m_lFlags & lFlag) == kAllFlagsOff?true:false );
}


bool CLZNumber::FlagIsOn( unsigned long lFlag )
{
    return ( (m_lFlags & lFlag) == kAllFlagsOff?false:true );
}


long CLZNumber::GetOrderOfMagnitude()
{
	long	l_leadingZeroes = 0;
	long	l_testIndex = m_lDigitsLength - 1;

	// Check for leading zeroes...
	while ( l_testIndex >= 0 )
	{
		if ( m_lpDigits[ l_testIndex ] > 0 )
			break;

		l_leadingZeroes++;
		l_testIndex--;
	}

    long l_orderOfMagnitude = m_lDigitsLength - m_lRadix - l_leadingZeroes;

	return (l_orderOfMagnitude < 0 ? 0:l_orderOfMagnitude);
}


/*
	The character a digit is written with: the only place that decides
	this on the way out, as lzi_digit_value is the only one on the way
	in. The two have to agree, or what is written does not read back as
	itself; the round-trip tests are what hold them to it.

	GetStringRep carried a sixteen-case copy of this until the two were
	collapsed. Base twelve spells ten and eleven D and E.
*/

static char lzi_digit_char( long lValue, long lBase )
{
	if ( (lValue >= 0) && (lValue <= 9) )
		return static_cast<char>( '0' + lValue );

	if ( lBase == 12 )
	{
		if ( lValue == 10 )	return 'D';
		if ( lValue == 11 )	return 'E';
	}

	if ( (lValue >= 10) && (lValue <= 15) )
		return static_cast<char>( 'A' + (lValue - 10) );

	return '?';
}


/*
	Which mark introduces an exponent when writing in this base: E
	where the base leaves it free, X where E is a digit. The reader
	takes either, so this only has to pick the conventional one.
*/

static char lzi_exponent_char( long lBase )
{
	return lzi_is_exponent_mark( 'E', lBase ) ? 'E' : 'X';
}


/*
	The index of the most significant digit that is not a nought, or -1
	where there is none. Both exponent renderers measure the exponent
	from it, and the display renderer counts significant digits down
	from it.
*/

static long lzi_high_digit( const CLZNumber* pNumber )
{
	for ( long i = pNumber->m_lDigitsLength - 1; i >= 0; i-- )
	{
		if ( pNumber->m_lpDigits[ i ] != 0 )
			return i;
	}

	return -1;
}


/*
	The digit at an index, with noughts beyond either end of the array.

	A display asks for a fixed number of digits and the value may not
	have that many: SCI 9 of a value carrying four figures wants five
	noughts that are nowhere in the array. Reading through this saves
	the writer from having to know where the digits stop.

	It is also what pads a short value out to the width, so nothing
	else has to. GetDisplayRep called PadRight for that until a
	mutation removing the call changed no test's answer: every nought
	it added was one this was already supplying, and the padding was
	an allocation and a memmove for nothing.
*/

static long lzi_digit_at( const CLZNumber* pNumber, long lIndex )
{
	// The bound is the digits the value has, not the room the array
	//	holds. Widening it to m_lDigitsArrayLength changes no answer
	//	today, because the slack above the digits is always zero -- but
	//	that is a property of how the array is kept, and the value
	//	genuinely stops at m_lDigitsLength...
	if ( (lIndex < 0) || (lIndex >= pNumber->m_lDigitsLength) )
		return 0;

	return pNumber->m_lpDigits[ lIndex ];
}


/*
	lBefore digits, a radix mark, then lAfter digits, reading down the
	array from lTop. Returns where it left off.

	The mark is written only when something follows it, so FIX 0 reads
	123 rather than 123 with a full stop hanging off the end.
*/

static long lzi_write_digits( const CLZNumber* pNumber, char* cpOut, long lOut,
							  long lTop, long lBefore, long lAfter )
{
	long	l_index = lTop;

	for ( long i = 0; i < lBefore; i++ )
		cpOut[ lOut++ ] = lzi_digit_char( lzi_digit_at( pNumber, l_index-- ), pNumber->m_lBase );

	if ( lAfter > 0 )
	{
		cpOut[ lOut++ ] = '.';

		for ( long i = 0; i < lAfter; i++ )
			cpOut[ lOut++ ] = lzi_digit_char( lzi_digit_at( pNumber, l_index-- ), pNumber->m_lBase );
	}

	return lOut;
}


/*
	Write a display of the value, lDigits wide -- see GetStringRep in
	the header for what separates this from writing the value itself.

	FIX counts places after the radix mark, which is what the engine's
	own precision counts, so Round does the work directly. SCI and ENG
	count significant digits, and for a value of any size that means
	rounding somewhere above the mark, which Round has no way to be
	asked for: its precision is places after the mark and cannot go
	negative. Moving the mark to just behind the leading digit turns
	the one question into the other. Nothing about the value moves --
	only where this copy claims the point sits -- and the exponent is
	worked back out afterwards from where the leading digit ended up,
	which also catches the carry that makes 9.99 into 1.0E1.
*/

char* CLZNumber::GetDisplayRep( RenderStyle rsStyle, long lDigits )
{
	if ( lDigits < 0 )
		lDigits = 0;

	// The rounding is done to a copy, and the copy inherits the
	//	constant flag that Round refuses to work through...
	CLZNumber	shown( *this );

	shown.SetConstant( false );

	long	l_high = lzi_high_digit( &shown );
	long	l_exponent = 0;
	long	l_top = 0;
	long	l_before = 1;
	long	l_after = lDigits;

	if ( rsStyle == RS_positional )
	{
		shown.Round( lDigits );

		l_high = lzi_high_digit( &shown );

		/*
			Writing starts at the leading digit, or at the mark where
			the value has nothing above it -- that is the nought in
			front of the point, and it is written rather than
			suppressed. A rounded-away value reads 0.000 and a half
			reads 0.50.
		*/

		l_top = (l_high > shown.m_lRadix) ? l_high : shown.m_lRadix;
		l_before = (l_top - shown.m_lRadix) + 1;
	}

	// Nought has no leading digit to measure an exponent from, and
	//	reads 0.00E0 in both of the exponent styles...
	else if ( l_high >= 0 )
	{
		long	l_scale = l_high - shown.m_lRadix;

		shown.m_lRadix = l_high;
		shown.Round( lDigits );

		l_high = lzi_high_digit( &shown );
		l_exponent = (l_high - shown.m_lRadix) + l_scale;
		l_top = l_high;

		if ( rsStyle == RS_engineering )
		{
			/*
				Down to the multiple of three at or below the exponent,
				taking the remainder toward negative infinity rather
				than toward zero the way C does it, or -4 would group
				as -3 and put four digits before the mark.

				The digits that move in front of the mark come out of
				the count behind it: the display shows lDigits+1
				figures whichever side of the mark they fall on, which
				is what makes ENG and SCI agree on how much of the
				value is being shown.
			*/

			long	l_remainder = l_exponent % 3;

			if ( l_remainder < 0 )
				l_remainder += 3;

			l_before = l_remainder + 1;
			l_exponent -= l_remainder;
			l_after = (lDigits > l_remainder) ? (lDigits - l_remainder) : 0;
		}
	}

	// Sign, the digits either side of the mark, the exponent mark and
	//	its own sign and digits, and the terminator...
	char*	cp_result = static_cast<char*>( ::calloc( shown.m_lDigitsLength + lDigits + 32, sizeof(char) ) );

	if ( cp_result == nullptr )
		return nullptr;

	long	l_out = 0;

	// A value that rounded away to nothing has had its sign cleared by
	//	CheckForZero, so a small negative shows as 0.000 and not -0.000...
	if ( shown.m_lSign == kNegative )
		cp_result[ l_out++ ] = '-';

	l_out = lzi_write_digits( &shown, cp_result, l_out, l_top, l_before, l_after );

	if ( rsStyle != RS_positional )
	{
		cp_result[ l_out++ ] = lzi_exponent_char( shown.m_lBase );

		// The exponent counts positions, so its own digits are decimal
		//	whatever base the value is in...
		::snprintf( cp_result + l_out, 20, "%ld", l_exponent );
	}

	return cp_result;
}


/*
	Write the value with an exponent -- see RenderStyle.

	Scientific puts one digit before the radix mark. Engineering moves
	the exponent down to a multiple of three and lets one, two or
	three digits stand before the mark, which is what makes the
	exponent read as kilo, mega, milli and so on.

	Every digit the value holds is written. Rounding belongs to the
	caller: see the note at RenderStyle.
*/

char* CLZNumber::GetStringRep( RenderStyle rsStyle )
{
	if ( rsStyle == RS_positional )
		return GetStringRep();

	// The exponent is measured from the most significant digit that is
	//	not a nought, not from the top of the array: above it there may
	//	be leading noughts the positional form suppresses...
	long	l_high = lzi_high_digit( this );

	// Sign, the digits, a radix mark, the exponent mark and its own
	//	sign and digits, and the terminator...
	char*	cp_result = static_cast<char*>( ::calloc( m_lDigitsLength + 24, sizeof(char) ) );

	if ( cp_result == nullptr )
		return nullptr;

	// Nought has no exponent worth writing...
	if ( l_high < 0 )
	{
		cp_result[ 0 ] = '0';
		return cp_result;
	}

	long	l_exponent = l_high - m_lRadix;
	long	l_before = 1;

	if ( rsStyle == RS_engineering )
	{
		/*
			Down to the multiple of three at or below the exponent.
			The remainder has to be taken toward negative infinity,
			not toward zero the way C does it, or -4 would group as -3
			and put four digits before the mark.
		*/

		long	l_remainder = l_exponent % 3;

		if ( l_remainder < 0 )
			l_remainder += 3;

		l_before = l_remainder + 1;
		l_exponent -= l_remainder;
	}

	long	l_out = 0;

	if ( m_lSign == kNegative )
		cp_result[ l_out++ ] = '-';

	long	l_written = 0;

	for ( long i = l_high; i >= 0; i-- )
	{
		if ( l_written == l_before )
			cp_result[ l_out++ ] = '.';

		cp_result[ l_out++ ] = lzi_digit_char( m_lpDigits[ i ], m_lBase );
		l_written++;
	}

	// An engineering group can want more digits before the mark than
	//	the value has left to give...
	while ( l_written < l_before )
	{
		cp_result[ l_out++ ] = '0';
		l_written++;
	}

	cp_result[ l_out++ ] = lzi_exponent_char( m_lBase );

	// The exponent counts positions, so its own digits are decimal
	//	whatever base the value is in...
	::snprintf( cp_result + l_out, 20, "%ld", l_exponent );

	return cp_result;
}


char* CLZNumber::GetStringRep()
{
	long	l_digitIndex = 0;
	long	l_stringIndex = 0;
	/*
		Room for the digits, a sign, a radix mark, the nought written
		in front of the mark, and the terminator.

		Three was right while that nought did not exist, and is still
		right whenever the value carries a suppressed digit above the
		mark to pay for it. **A fraction the engine computes usually
		does not**: m_lRadix equals m_lDigitsLength for 1/2, for
		sin(1), for ln(2), for anything read from ".5". Those overrun
		by exactly one byte at three, which ASan reports and a plain
		build does not.
	*/

    char*	cp_result = static_cast<char*>(::calloc( m_lDigitsLength + 4, sizeof(char) ));

	// Handle allocation failure...
	if ( cp_result == nullptr )
		return nullptr;

    if (m_lDigitsLength == 0)
    {
        cp_result[ l_digitIndex++ ] = '0';
    }

    else
    {
        // Apply sign...
        if (m_lSign == kNegative)
            cp_result[l_digitIndex++] = '-';

        /*
        	A leading zero is one with no digit yet written before it,
        	and that has to be judged from where the digits start
        	rather than from the start of the string. The test below
        	used to ask whether anything had been written at all,
        	which a minus sign satisfies -- so a negative value kept
        	every leading zero it had and sixteen less a hundred read
        	-084.
        */

        long	l_digitsBegin = l_digitIndex;

        // Write out string...
        for (long i = 0; i < m_lDigitsLength; i++)
        {
            // Append radix mark...
            if (((m_lDigitsLength - i) == m_lRadix) && (m_lRadix > 0))
            {
                /*
                	A value with nothing above the mark still shows the
                	nought in front of it. Every digit up here was a
                	leading zero and was suppressed, which is right for
                	the ones further up and wrong for this one: a half
                	is 0.5, not .5. The leading-zero test cannot make
                	the distinction on its own -- it has no idea where
                	the mark is going to fall -- so the mark says so on
                	its way past.
                */

                if (l_digitIndex == l_digitsBegin)
                    cp_result[l_digitIndex++] = '0';

                cp_result[l_digitIndex++] = '.';
            }

            /*
            	Append the digit.

            	The sixteen cases this replaces spelled out the same
            	mapping lzi_digit_char holds, base twelve's D and E
            	included, and the two had to be kept in step by hand.
            	The only part that was ever this loop's own business is
            	the leading nought below.
            */

            long	l_value = m_lpDigits[m_lDigitsLength - (++l_stringIndex)];

            /*
            	A leading nought is one with no digit written before
            	it, and that is judged from where the digits start
            	rather than from the start of the string -- a minus
            	sign is not a digit. Getting that wrong is how sixteen
            	less a hundred came to read -084.

            	A value one digit wide is all nought and keeps it.
            */

            if ( (l_value != 0) || (l_digitIndex > l_digitsBegin) || (m_lDigitsLength == 1) )
                cp_result[l_digitIndex++] = lzi_digit_char( l_value, m_lBase );
        }

        /*
        	A value that is all zeros and more than one digit wide had
        	every one of them suppressed as a leading zero and came
        	back as the empty string: Assign("00") rendered as nothing
        	at all. Nought is written with a nought.

        	Only reachable through text -- the arithmetic trims -- which
        	is why it survived: nothing wrote "00" until exponent
        	notation began expanding 0X10 into eleven of them.
        */

        if ( l_digitIndex == l_digitsBegin )
            cp_result[ l_digitIndex++ ] = '0';
    }

	return cp_result;
}


/*
	The base a two-character prefix names, or 0 for anything that is
	not one. **Elzed's convention**, which Addaptic follows:

	    0b   binary        0b10101
	    0o   octal         0o25
	    0d   decimal       0d21, and the default with no prefix at all
	    0t   duodecimal    0t19     t for twelve; its extra digits are
	                                D for dec and E for el
	    0x   hexadecimal   0x15     either x or h is taken
	    0h   hexadecimal   0h15

	The leading character is always a nought and the identifier always
	lower case, which is what keeps the prefix apart from everything
	else that could start a number. It matters here: this engine reads
	a lower-case x as an exponent mark, so without the position rule
	"0x15" would be nought times ten to the fifteenth rather than
	twenty-one. **A prefix at the front wins; an x anywhere else is an
	exponent.** "0X15" is not a prefix -- the identifier is lower case
	-- and does mean nought times ten to the fifteenth.

	Addaptic differs from Elzed on one point deliberately. Elzed reads
	unprefixed text as base ten always and says the default input base
	cannot be changed; here a CLZNumber carries a base of its own, and
	unprefixed text is read in it, which is what SetBase has always
	meant. A prefix overrides that and takes the value with it: after
	Assign("0x15") the number is a base-sixteen one.
*/

static long lzi_prefix_base( char cIdentifier )
{
	switch ( cIdentifier )
	{
		case 'b':	return 2;
		case 'o':	return 8;
		case 'd':	return 10;
		case 't':	return 12;
		case 'x':	return 16;
		case 'h':	return 16;
		default:	return 0;
	}
}


/*
	The exponent counting the number's own base makes it exactly a move
	of the radix mark, which is why expanding the text below is exact
	and cheap.

	Its own digits are read in base ten whatever the value's base: it
	is a count of positions, not a quantity in the number.
*/

const long kMaxExponentDigits = 7;
const long kMaxExpansion = 100000;


/*
	Is this text a number, in the base the value is held in?

	Nothing asked this before, and the decoder below is a switch with
	no default case, so anything it did not recognise was skipped in
	silence: "12x34" became 1234. A sign was taken wherever it
	appeared rather than only at the front, so the minus in "1.5E-8"
	made the whole value negative. A second radix mark was accepted
	and simply moved the first. And A through F were read as ten
	through fifteen whatever the base, so "6.02e23" put a digit of
	fourteen into a base-ten number.

	Every one of those returned kNoError with a plausible-looking
	value attached, which is the shape of fault this engine has spent
	its life turning up.

	The rule: an optional leading sign, at most one radix mark, at
	least one digit, every digit valid for the base, an optional
	exponent introduced by X or by E where the base leaves E free, and
	nothing else anywhere.

	The pieces are handed back so the caller can expand an exponent
	without parsing the text a second time and risking a different
	answer.
*/

struct SNumberText
{
	const char*	cpMantissa;			// first character after any sign
	long		lIntegerDigits;
	long		lFractionDigits;
	bool		bNegative;
	bool		bHasExponent;
	long		lExponent;			// already signed
	bool		bHasPrefix;
	long		lBase;				// from the prefix, or the default
};

static bool lzi_text_is_a_number( const char *cpTextNumber, long lDefaultBase,
								  SNumberText* pParts = nullptr )
{
	if ( cpTextNumber == nullptr )
		return false;

	SNumberText	theParts = { nullptr, 0, 0, false, false, 0, false, lDefaultBase };

	const char*	p = cpTextNumber;

	if ( (*p == '+') || (*p == '-') )
	{
		theParts.bNegative = (*p == '-');
		p++;
	}

	// A base prefix, if there is one, and only at the front -- see
	//	lzi_prefix_base for why the position matters...
	if ( *p == '0' )
	{
		long	l_prefixBase = lzi_prefix_base( p[ 1 ] );

		if ( l_prefixBase > 0 )
		{
			theParts.bHasPrefix = true;
			theParts.lBase = l_prefixBase;
			p += 2;
		}
	}

	long	lBase = theParts.lBase;

	theParts.cpMantissa = p;

	bool	b_seenRadix = false;

	for ( ; *p != '\0'; p++ )
	{
		if ( *p == '.' )
		{
			if ( b_seenRadix )
				return false;

			b_seenRadix = true;
			continue;
		}

		if ( lzi_is_exponent_mark( *p, lBase ) )
			break;

		long	l_value = lzi_digit_value( *p, lBase );

		if ( (l_value < 0) || (l_value >= lBase) )
			return false;

		if ( b_seenRadix )	theParts.lFractionDigits++;
		else				theParts.lIntegerDigits++;
	}

	if ( (theParts.lIntegerDigits + theParts.lFractionDigits) == 0 )
		return false;

	// Anything left is the exponent, and it has its own small grammar...
	if ( *p != '\0' )
	{
		theParts.bHasExponent = true;
		p++;

		bool	b_exponentNegative = false;

		if ( (*p == '+') || (*p == '-') )
		{
			b_exponentNegative = (*p == '-');
			p++;
		}

		long	l_digits = 0;
		long	l_value = 0;

		for ( ; *p != '\0'; p++ )
		{
			if ( (*p < '0') || (*p > '9') )
				return false;

			// Bounded before it can overflow, and well before the
			//	expansion below could be built at all...
			if ( ++l_digits > kMaxExponentDigits )
				return false;

			l_value = (l_value * 10) + (*p - '0');
		}

		if ( l_digits == 0 )
			return false;

		theParts.lExponent = b_exponentNegative ? -l_value : l_value;
	}

	if ( pParts != nullptr )
		*pParts = theParts;

	return true;
}


/*
	Rewrite exponent notation as the plain positional text it stands
	for, so the decoder below needs to know nothing about it.

	The whole value is written out: this representation has no exponent
	field, only a radix position, so 1e300 really is three hundred and
	one digits and there is no cheaper way to hold it. That is what
	kMaxExpansion bounds -- a typo of 1e99999999 should be refused
	rather than allocated.

	Returns a string the caller frees, or nullptr if it will not fit.
*/

static char* lzi_expand_exponent( const SNumberText& theParts )
{
	long	l_mantissaDigits = theParts.lIntegerDigits + theParts.lFractionDigits;

	// Where the radix mark lands once the exponent has moved it...
	long	l_pointAt = theParts.lIntegerDigits + theParts.lExponent;

	long	l_leadingZeros = (l_pointAt < 0) ? -l_pointAt : 0;
	long	l_trailingZeros = (l_pointAt > l_mantissaDigits)
							  ? (l_pointAt - l_mantissaDigits) : 0;

	long	l_total = l_leadingZeros + l_mantissaDigits + l_trailingZeros;

	if ( l_total > kMaxExpansion )
		return nullptr;

	// Sign, a possible leading nought, the radix mark, the digits and
	//	the terminator...
	char*	cp_plain = static_cast<char*>( ::calloc( l_total + 4, sizeof(char) ) );

	if ( cp_plain == nullptr )
		return nullptr;

	long	l_out = 0;

	if ( theParts.bNegative )
		cp_plain[ l_out++ ] = '-';

	if ( l_pointAt <= 0 )
		cp_plain[ l_out++ ] = '.';

	for ( long i = 0; i < l_leadingZeros; i++ )
		cp_plain[ l_out++ ] = '0';

	// The mantissa's own digits, skipping the radix mark it came with...
	long	l_written = 0;

	for ( const char* q = theParts.cpMantissa; l_written < l_mantissaDigits; q++ )
	{
		if ( *q == '.' )
			continue;

		if ( (l_written == l_pointAt) && (l_pointAt > 0) )
			cp_plain[ l_out++ ] = '.';

		cp_plain[ l_out++ ] = *q;
		l_written++;
	}

	for ( long i = 0; i < l_trailingZeros; i++ )
		cp_plain[ l_out++ ] = '0';

	cp_plain[ l_out ] = '\0';

	return cp_plain;
}


void CLZNumber::Initialize( const char *cpTextNumber, unsigned long lDigitsLength )
{
	bool	b_wasInitialized = m_bInitialized;
	long	l_digitIndex = 0;
	bool	b_appliedSign = false;
	bool	b_appliedRadix = false;
	bool	b_allZeroDigits = true;
	bool	b_malformed = false;

	m_bInitialized = false;

	m_lFlags = kAllFlagsOff;
	m_lSign = kPositive;
	m_lAdjust = 0;
	m_lPrecision = glWorkingPrecision;
	m_rmMode = kDefaultRoundingMode;
	m_rdDigit = RD_digit_precision;

	m_lDigitsLength = lDigitsLength;
	m_lRadix = 0;

	if ( !b_wasInitialized )
	{
		m_lBase = 10;
		m_lMaxDigit = 9;
		m_lHalfDigit = 5;

		m_lDigitsArrayLength = kInitialDigitsLength;
		m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof(long)) );
	}

	/*
		Text that is not a number becomes nought carrying a
		not-a-number flag, rather than whatever the decoder below
		would have made of it. Taken here, after the base is settled,
		because what counts as a digit depends on it -- see
		lzi_text_is_a_number.
	*/

	char*	cp_expanded = nullptr;

	if ( cpTextNumber != nullptr )
	{
		SNumberText	theParts;

		if ( !lzi_text_is_a_number( cpTextNumber, m_lBase, &theParts ) )
		{
			b_malformed = true;
			cpTextNumber = nullptr;
		}

		else
		{
			/*
				A prefix takes the value's base with it, so this has to
				land before the decoder runs -- the switch below reads
				D and E differently in base twelve.
			*/

			if ( theParts.lBase != m_lBase )
			{
				m_lBase = theParts.lBase;
				m_lMaxDigit = m_lBase - 1;
				m_lHalfDigit = m_lBase / 2;
			}

			/*
				Rewritten as plain positional text whenever it is not
				already: the decoder understands neither prefixes nor
				exponents and does not need to.
			*/

			if ( theParts.bHasPrefix || theParts.bHasExponent )
			{
				cp_expanded = lzi_expand_exponent( theParts );

				if ( cp_expanded == nullptr )
				{
					b_malformed = true;
					cpTextNumber = nullptr;
				}
				else
				{
					cpTextNumber = cp_expanded;
				}
			}
		}
	}

	// Set up the digits array...
	if ( cpTextNumber != nullptr )
	{
		// A text of two billion characters is past anything the engine is
		//	asked to read, so its length fits a long even where a long is
		//	32 bits, as on Windows...
		m_lDigitsLength = static_cast<long>( strlen( cpTextNumber ) );

		// Reallocate, if required...
		if ( m_lDigitsLength >= m_lDigitsArrayLength )
		{
			m_lDigitsArrayLength = m_lDigitsLength + 2;

			::free( m_lpDigits );
			m_lpDigits = nullptr;

			m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof(long)) );
		}

		else
		{
			::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );
		}

		// Decode the input string...
		for ( long i=0; i<m_lDigitsLength; i++ )
		{
			switch ( cpTextNumber[ m_lDigitsLength - i - 1 ] )
			{
				case '+':
					m_lSign = kPositive;
					b_appliedSign = true;
					break;

				case '-':
					m_lSign = kNegative;
					b_appliedSign = true;
					break;

				case '.':
					m_lRadix = i;	// Update exponent...
					b_appliedRadix = true;
					break;

				/*
					Everything that is not the sign or the mark is a
					digit, and what it is worth is lzi_digit_value's
					question alone.

					This was sixteen cases spelling the mapping out a
					second time, with base twelve's D and E special-cased
					in both places. The comment on lzi_digit_value asked
					for them to be collapsed into a call; this is that.

					A character that is no digit at all is skipped, which
					is what the switch did by having no default case.
					Nothing reaches here to be skipped -- the text was
					checked by lzi_text_is_a_number before the decode
					began -- but the two agreeing on what a digit is
					worth is precisely the point, so it costs nothing to
					say so.
				*/

				default:
				{
					long	l_value = lzi_digit_value( cpTextNumber[ m_lDigitsLength - i - 1 ], m_lBase );

					if ( l_value >= 0 )
					{
						m_lpDigits[ l_digitIndex++ ] = l_value;

						if ( l_value != 0 )
							b_allZeroDigits = false;
					}

					break;
				}
			}
		}

		// Remove radix mark from digit count...
		if ( b_appliedRadix )
			m_lDigitsLength--;

		// Remove sign from digit count...
		if ( b_appliedSign )
			m_lDigitsLength--;
	}

	else
	{
		m_lDigitsLength = 1;

		if ( m_lpDigits != nullptr )
			::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );

		else
		{
			m_lDigitsArrayLength = kInitialDigitsLength;
			m_lpDigits = static_cast<long*>(::calloc( m_lDigitsArrayLength, sizeof(long)) );
		}
	}

	if ( b_allZeroDigits )
	{
		m_lFlags = kfZero;
		m_lSign = kZero;
	}

	// After the zero test, which assigns m_lFlags outright rather than
	//	adding to it...
	if ( b_malformed )
		SetFlag( kfQNaN );

	::free( cp_expanded );

	m_bInitialized = true;
}


bool CLZNumber::IsConstant()
{
	return m_bConstant;
}


bool CLZNumber::IsFraction()
{
	bool b_isFraction = true;

	// The integer digits are m_lRadix .. m_lDigitsLength-1; scan them
	// downward, so start on the last valid digit, not one past it...
	long l_compareDigits = m_lDigitsLength - 1;

	if ( IsZero() )
		b_isFraction = false;

	else while ( l_compareDigits >= m_lRadix )
	{
		if ( m_lpDigits[ l_compareDigits-- ] == 0 )
			continue;

		b_isFraction = false;
		break;
	}

	return b_isFraction;
}


bool CLZNumber::IsInfinity()
{
	return ( (m_lFlags & kfInfinity) == kAllFlagsOff?false:true );
}


bool CLZNumber::IsInteger()
{
	bool b_isInteger = true;
	long l_compareDigits = 0;

	while ( l_compareDigits < m_lRadix )
	{
		if ( m_lpDigits[ l_compareDigits++ ] == 0 )
			continue;

		b_isInteger = false;
		break;
	}

	return b_isInteger;
}


bool CLZNumber::IsNegative()
{
	bool b_isNegative = true;

	if ( m_lSign != kNegative )
		b_isNegative = false;

	return b_isNegative;
}


bool CLZNumber::IsPositive()
{
	bool b_isPositive = true;

	if ( m_lSign != kPositive )
		b_isPositive = false;

	return b_isPositive;
}


bool CLZNumber::IsZero()
{
	return ( (m_lFlags & kfZero) == kAllFlagsOff?false:true );
}


long CLZNumber::PadLeft( long lPadDigits )
{
	long*	lp_newStorage = nullptr;

	if ( lPadDigits > 0 )
	{
		if ( (m_lDigitsLength + lPadDigits) > m_lDigitsArrayLength )
		{
			m_lDigitsArrayLength = m_lDigitsLength + lPadDigits;

			// Create new storage & copy in digits...
			lp_newStorage = static_cast<long*>(::calloc(m_lDigitsArrayLength, sizeof(long)) );
			::memcpy( lp_newStorage, m_lpDigits, m_lDigitsLength * sizeof(long) );

			// Free the old array...
			::free( m_lpDigits );

			// Save the new info...
			m_lpDigits = lp_newStorage;
		}

		else
		{
			::memset( m_lpDigits + m_lDigitsLength, 0, lPadDigits * sizeof(long) );
		}

		m_lDigitsLength += lPadDigits;
	}

	return kNoError;
}


long CLZNumber::PadRight( long lPadDigits )
{
	if ( ShiftLeft( lPadDigits ) == kNoError )
	{
		m_lRadix += lPadDigits;
	}

	return kNoError;
}


void CLZNumber::EnterInternalWork()
{
	ms_lInternalWork++;
}


void CLZNumber::LeaveInternalWork()
{
	if ( ms_lInternalWork > 0 )
		ms_lInternalWork--;
}


// Recording counts its requests, as internal work does, so a function
//	that asks inside another leaves it on for the outer one -- see
//	LastDiscardedTail in LZNumber.h...
void CLZNumber::RecordRoundings( bool bRecord )
{
	if ( bRecord )
		ms_lRecordRoundings++;
	else if ( ms_lRecordRoundings > 0 )
		ms_lRecordRoundings--;
}


void CLZNumber::ForgetLastRounding()
{
	ms_dLastTail = -1.0;
	ms_lLastPlaces = 0;
}


double CLZNumber::LastDiscardedTail()
{
	return ms_dLastTail;
}


long CLZNumber::LastRoundedPlaces()
{
	return ms_lLastPlaces;
}


// The mode and the digit a Round that is not told them uses: the
//	number's own, except inside the engine's own work -- see
//	EnterInternalWork...
RoundingMode CLZNumber::OwnMode() const
{
	return (ms_lInternalWork > 0) ? kDefaultRoundingMode : m_rmMode;
}


RoundingDigit CLZNumber::OwnDigit() const
{
	return (ms_lInternalWork > 0) ? RD_digit_precision : m_rdDigit;
}


long CLZNumber::Round()
{
	return (Round( OwnMode(), OwnDigit() ));
}


long CLZNumber::Round( RoundingMode rm_mode )
{
	return (Round( rm_mode, OwnDigit() ));
}


long CLZNumber::Round( long lPrecision )
{
	return (Round( OwnMode(), lPrecision ));
}


long CLZNumber::Round( RoundingDigit rd_digit )
{
	return (Round( OwnMode(), rd_digit ));
}


long CLZNumber::Round( RoundingMode rm_mode, RoundingDigit rd_digit )
{
	return (RoundAt( rm_mode, PrecisionFor( rd_digit ), false ));
}


long CLZNumber::Round( RoundingMode rm_mode, long lPrecision )
{
	return (RoundAt( rm_mode, lPrecision, false ));
}


long CLZNumber::RoundMiddleRepeating()
{
	return (RoundAt( OwnMode(), PrecisionFor( OwnDigit() ), true ));
}


long CLZNumber::RoundMiddleRepeating( RoundingMode rm_mode, long lPrecision )
{
	return (RoundAt( rm_mode, lPrecision, true ));
}


long CLZNumber::PrecisionFor( RoundingDigit rd_digit )
{
	long			lPrecision = 0;
	RoundingDigit	rd_targetDigit = RD_digit_default;

	// Set target precision...
	while ( rd_targetDigit == RD_digit_default )
	{
		switch( rd_digit )
		{
		case RD_digit_zero:
			rd_targetDigit = RD_digit_zero;
			lPrecision = 0;
			break;

		case RD_digit_last:
			rd_targetDigit = RD_digit_last;
			if ( m_lRadix > 0 )
				lPrecision = m_lRadix - 1;
			else
				lPrecision = 0;
			break;

		case RD_digit_precision:
			rd_targetDigit = RD_digit_precision;
			lPrecision = m_lPrecision;
			break;

		case RD_digit_default:
		default:
			rd_digit = OwnDigit();	// Set default, go around again...
			break;
		}
	}

	return lPrecision;
}


/*
	**Where the discarded digits stand against half a place is worked
	out once, before the modes**, as lTail: below it, exactly it, or
	above it. Every half mode then asks lTail and nothing else.

	In an even base that is the comparison the modes always made: a
	guard digit above half is above, one below is below, and one equal
	to half is exactly half only if everything after it is nought.

	**An odd base has no half digit**, and until 2026-09-24 the modes
	treated its middle digit as one -- m_lHalfDigit is base / 2, three
	in base seven, and .3 there is 3/7, less than a half. So 1.33 base
	seven, 1.49 in decimal, rounded to 2 in every half mode. The rule
	there is the one for comparing .3333... with any other string: the
	first discarded digit that is not the middle digit decides, above
	it or below it; if every one is the middle digit the value is
	below half, a finite run of them always falling short. Nothing
	finite in an odd base is ever exactly half -- only a value whose
	middle digit repeats forever is, which is what bMiddleRepeats says
	for a caller that knows it; see RoundMiddleRepeating.
*/

long CLZNumber::RoundAt( RoundingMode rm_mode, long lPrecision, bool bMiddleRepeats )
{
	bool	bCheckCarry = false;
	bool	bZeroDiscard = true;
	long	lDigitIndex = 0;
	long	lBeyondGuard = 0;	// odd bases: the first discarded digit past the guard that is not the middle one, against it
	long	lTail = 0;			// the discarded digits against half a place: -1, 0 or 1

	// OK to continue?
	if ( m_bConstant )
		return kTargetIsConstant;

	/*
		**A negative count rounds to the left of the point** -- since
		2026-10-15, found by that day's review: round( 1234.5, -2 ) is
		1200. It discarded whole digits and moved the radix past them,
		answering a wrong magnitude; before that it shifted forever. The
		value over b^k, rounded at the point, then times b^k: the first a
		move of the radix, with noughts above where the digits run out, the
		last noughts below...
	*/
	if ( (lPrecision < 0) && (rm_mode != RM_round_noop) )
	{
		long	l_places = -lPrecision;

		m_lRadix += l_places;

		if ( m_lDigitsLength <= m_lRadix )
			PadLeft( m_lRadix + 1 - m_lDigitsLength );

		long	l_result = RoundAt( rm_mode, 0, bMiddleRepeats );

		CheckForZero();

		if ( !IsZero() )
			ShiftLeft( l_places );

		return l_result;
	}

	// Perform rounding...
	if ( rm_mode != RM_round_noop )
	{
		if ( m_lRadix > lPrecision )
		{
			/*
				Discard digits beyond lPrecision + 1: read them where they
				lie, then shift once.

				**This shifted one digit at a time**, and each shift moved
				and cleared the whole array, so discarding n digits of an
				m-digit number cost n times m. Rounding a product back to
				its width discards half of it: e^50000, worked at 21,750
				places, spent 98% of its 1.7 s here and 2% multiplying.
				Reading from the least significant up and shifting by the
				count is the same walk -- ShiftRight runs off the end
				exactly as that many single shifts do.
			*/

			long	lDiscard = m_lRadix - (lPrecision + 1);

			// The guard digit and the digits after it, read before they go,
			//	for a function deciding whether its answer was in doubt --
			//	see LastDiscardedTail...
			if ( ms_lRecordRoundings > 0 )
			{
				double	d_tail = 0.0;
				double	d_scale = 1.0;
				double	d_step = 1.0 / m_lBase;

				for ( long k = lDiscard, n = 0; (k >= 0) && (n < kRecordedTailDigits); k--, n++ )
				{
					d_scale *= d_step;

					if ( k < m_lDigitsLength )
						d_tail += m_lpDigits[ k ] * d_scale;

					// Clear of every window whatever follows: stop reading...
					if ( (d_tail > 0.5 + kRecordedTailReach) || (d_tail + d_scale < 0.5 - kRecordedTailReach) )
						break;
				}

				ms_dLastTail = d_tail;
				ms_lLastPlaces = lPrecision;
			}

			if ( lDiscard > 0 )
			{
				long	lScan = (lDiscard < m_lDigitsLength) ? lDiscard : m_lDigitsLength;

				// Its own index: lDigitIndex has to reach the carry below
				//	still at nought...
				for ( long lScanIndex = 0; lScanIndex < lScan; lScanIndex++ )
				{
					if ( m_lpDigits[ lScanIndex ] != 0 )
						bZeroDiscard = false;

					// Each is more significant than the last, so the last to differ wins...
					if ( m_lpDigits[ lScanIndex ] != m_lHalfDigit )
						lBeyondGuard = (m_lpDigits[ lScanIndex ] > m_lHalfDigit) ? 1 : -1;
				}

				ShiftRight( lDiscard );
			}

			if ( m_lpDigits[ 0 ] > m_lHalfDigit )
				lTail = 1;
			else if ( m_lpDigits[ 0 ] < m_lHalfDigit )
				lTail = -1;
			else if ( (m_lBase % 2) == 0 )
				lTail = bZeroDiscard ? 0 : 1;
			else if ( lBeyondGuard != 0 )
				lTail = lBeyondGuard;
			else
				lTail = bMiddleRepeats ? 0 : -1;

			// The default mode if the default is asked for, as
			//	SetRoundingMode has it -- until 2026-10-15 this truncated...
			if ( rm_mode == RM_round_default )
				rm_mode = kDefaultRoundingMode;

			// Now round the number off...
			switch( rm_mode )
			{
			case RM_round_down:
				ShiftRight( 1 );

				break;

			case RM_round_half_up:
				if ( lTail >= 0 )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				ShiftRight( 1 );

				break;

			case RM_round_half_even:

				// If greater than half, round up...
				if ( lTail > 0 )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				// If exactly half, the mode breaks the tie...
				else if ( lTail == 0 )
				{
					// If the digit is odd, increment...
					if ( m_lpDigits[ 1 ] % 2 )
					{
						(m_lpDigits[ 1 ])++;
						bCheckCarry = true;
					}
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				// In all cases, drop the leftmost digit...
				ShiftRight( 1 );

				break;

			case RM_round_ceiling:
				if ( (!bZeroDiscard || (m_lpDigits[ 0 ]) ) && (m_lSign == kPositive) )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				ShiftRight( 1 );

				break;

			case RM_round_floor:
				if ( (!bZeroDiscard || (m_lpDigits[ 0 ]) ) && (m_lSign == kNegative) )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				ShiftRight( 1 );

				break;

			case RM_round_half_down:
				if ( lTail > 0 )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				ShiftRight( 1 );

				break;

			case RM_round_half_odd:

				// If greater than half, round up...
				if ( lTail > 0 )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				// If exactly half, the mode breaks the tie...
				else if ( lTail == 0 )
				{
					// If the digit is even, increment...
					if ( !(m_lpDigits[ 1 ] % 2) )
					{
						(m_lpDigits[ 1 ])++;
						bCheckCarry = true;
					}
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				// In all cases, drop the leftmost digit...
				ShiftRight( 1 );

				break;

			case RM_round_up:
				if ( m_lpDigits[ 0 ] || !bZeroDiscard )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				ShiftRight( 1 );

				break;

			case RM_round_alternate:

				// If greater than half, round up...
				if ( lTail > 0 )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				// If exactly half, the mode breaks the tie...
				else if ( lTail == 0 )
				{
					// If flag is set, increment...
					if ( ms_bRoundAlternate )
					{
						(m_lpDigits[ 1 ])++;
						bCheckCarry = true;
						ms_bRoundAlternate = false;
					}
					else
						ms_bRoundAlternate = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				// In all cases, drop the leftmost digit...
				ShiftRight( 1 );

				break;

			case RM_round_random:

				// If greater than half, round up...
				if ( lTail > 0 )
				{
					(m_lpDigits[ 1 ])++;
					bCheckCarry = true;
				}

				// If exactly half, the mode breaks the tie...
				else if ( lTail == 0 )
				{
					// If flag is set, increment...
					if ( ms_bRoundAlternate )
					{
						(m_lpDigits[ 1 ])++;
						bCheckCarry = true;
						ms_bRoundAlternate = false;
					}
					else
						ms_bRoundAlternate = true;
				}

				if ( bCheckCarry && (m_lDigitsLength == 1) )
					m_lDigitsLength++;

				// In all cases, drop the leftmost digit...
				ShiftRight( 1 );

				break;

			case RM_round_truncate:

				// Drop the leftmost digit...
				ShiftRight( 1 );

				break;

			default:
				break;
			}

			// 'Normalize' the remaining digits if we incremented...
			if ( bCheckCarry )
			{
				while ( lDigitIndex < m_lDigitsLength )
				{
					if ( m_lpDigits[ lDigitIndex ] > m_lMaxDigit )
					{
						// Make room before carrying past the last digit...
						if ( (lDigitIndex + 1) >= m_lDigitsArrayLength )
							SetDigitsArraySize( m_lDigitsArrayLength + kDigitsLengthIncrement );

						m_lpDigits[ lDigitIndex ] = 0;
						(m_lpDigits[ lDigitIndex + 1 ])++;

						lDigitIndex++;
					}
					else
						break;
				}

				/*
					Where every digit was a nine the carry runs off the
					top, landing one place above the digits the value
					claims to have. Unless the value is made a digit
					longer that carry cannot be seen, and a figure that
					should have become one reads as zero instead.
				*/

				if ( (lDigitIndex == m_lDigitsLength) && (m_lpDigits[ lDigitIndex ] != 0) )
					m_lDigitsLength++;
			}
		}
	}

	/*
		Rounding is how a value most often becomes zero - anything
		smaller than the precision being kept rounds away to nothing.
		Since IsZero() reports a flag rather than reading the digits,
		that flag has to be brought up to date here, or a value that
		is now zero goes on claiming otherwise.
	*/

	CheckForZero();

	return kNoError;
}


long CLZNumber::SetBase( long lNewBase, bool bConvertValue )
{
	if ( (lNewBase > 1) && (lNewBase < 17) && (lNewBase != m_lBase) )
	{
		long	l_oldBase = m_lBase;

		// Set new number base parameters...
		m_lBase = lNewBase;
		m_lMaxDigit = lNewBase - 1;
		m_lHalfDigit = lNewBase / 2;

		if ( bConvertValue )
		{
			BaseConvertValue( l_oldBase );
		}
		else
		{
			// Clear the decks...
			::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );
		}
	}

	return kNoError;
}


/*
	**Convert as many places as the work in the other base needs, where
	the value holds them.** SetBase's budget is the value's own
	precision scaled by the logs, and the engine's constants are all made
	at the default 36 while holding fifty digits or more -- so a table
	entry met by a hex operand arrived with thirty hex places, 36 digits
	of its fifty, and every table-driven function in bases eleven to
	sixteen stopped there: the arcsine's coefficients, the reciprocal
	factorials, pi.

	The first version converted every place held, and that was too
	blunt: pi's 120 places became some 400 binary ones for work needing
	forty, and exp in base two ran twice as slow for nothing. So the
	budget is the working precision in the destination base, re-counted
	in the value's own and two places over, never below its own
	precision. For base two and eight that is the precision, as before;
	for twelve and sixteen, the few more places their work reaches.

	**Not capped at the places the value holds**, which it was at
	first, on the reasoning that more places than a value has add
	nothing. An exact short decimal is the case that breaks it: the
	arcsine's second coefficient is .075, three places held and exact,
	and in hex it never ends. Capped, it converted to thirty places,
	and its term stopped the arcsine at 42 digits in every base but
	ten, at any width.
*/

long CLZNumber::SetBaseForWork( long lNewBase )
{
	return SetBaseForWork( lNewBase, glWorkingPrecision );
}


long CLZNumber::SetBaseForWork( long lNewBase, long lPlaces )
{
	// Nothing to convert, and nothing to widen: a value in its own base
	//	keeps the precision it had, so base ten is untouched by any of
	//	this...
	if ( lNewBase == m_lBase )
		return kNoError;

	long	l_wanted = EquivalentPlaces( lPlaces, lNewBase, m_lBase ) + 2;

	if ( l_wanted > m_lPrecision )
		SetPrecision( l_wanted );

	return SetBase( lNewBase, true );
}


long CLZNumber::SetConstant( bool bConstant )
{
	m_bConstant = bConstant;

	return kNoError;
}


long CLZNumber::SetDigitsArraySize( long lDigitsArraySize, bool bPreserve )
{
	long*	lp_newStorage = nullptr;

	if ( lDigitsArraySize > m_lDigitsArrayLength )
	{
		m_lDigitsArrayLength = lDigitsArraySize;

		// Create new storage & copy in digits...
		lp_newStorage = static_cast<long*>(::calloc(m_lDigitsArrayLength, sizeof(long)) );

		if ( bPreserve )
			::memcpy( lp_newStorage, m_lpDigits, m_lDigitsLength * sizeof(long) );

		// Free the old array...
		::free( m_lpDigits );

		// Save the new info...
		m_lpDigits = lp_newStorage;
	}

	return kNoError;
}


long CLZNumber::SetFlag( unsigned long lFlag )
{
	m_lFlags |= lFlag;

	return kNoError;
}


long CLZNumber::SetPrecision( long lPrecision )
{
	long*	lp_newStorage = nullptr;

	m_lPrecision = lPrecision;

	if ( m_lPrecision > m_lDigitsArrayLength )
	{
		m_lDigitsArrayLength = m_lPrecision;

		// Create new storage & copy in digits...
		lp_newStorage = static_cast<long*>(::calloc(m_lDigitsArrayLength, sizeof(long)) );
		::memcpy( lp_newStorage, m_lpDigits, m_lDigitsLength * sizeof(long) );

		// Free the old array...
		::free( m_lpDigits );

		// Save the new info...
		m_lpDigits = lp_newStorage;
	}

	return kNoError;
}


long CLZNumber::SetRoundingDigit( RoundingDigit rd_digit )
{
	switch( rd_digit )
	{
	case RD_digit_zero:
	case RD_digit_last:
	case RD_digit_precision:
		m_rdDigit = rd_digit;
		break;

	case RD_digit_default:
	default:
		m_rdDigit = RD_digit_precision;
		break;
	}

	return kNoError;
}


long CLZNumber::SetRoundingMode( RoundingMode rm_mode )
{
	m_rmMode = rm_mode;

	switch( rm_mode )
	{
	case RM_round_down:
	case RM_round_half_up:
	case RM_round_half_even:
	case RM_round_ceiling:
	case RM_round_floor:
	case RM_round_half_down:
	case RM_round_half_odd:
	case RM_round_up:
	case RM_round_alternate:
	case RM_round_random:
	case RM_round_noop:
		m_rmMode = rm_mode;
		break;

	case RM_round_default:
	default:
		m_rmMode = kDefaultRoundingMode;
		break;
	}

	return kNoError;
}


long CLZNumber::ShiftLeft( long lShiftDigits )
{
	long*	lp_newStorage = nullptr;

	if ( lShiftDigits > 0 )
	{
		if ( (m_lDigitsLength + lShiftDigits) >= (m_lDigitsArrayLength - 8) )
		{
			m_lDigitsArrayLength = m_lDigitsLength + lShiftDigits + 16;

			// Create new storage & copy in digits...
			lp_newStorage = static_cast<long*>(::calloc(m_lDigitsArrayLength, sizeof(long)) );
			::memcpy( lp_newStorage + lShiftDigits, m_lpDigits, m_lDigitsLength * sizeof(long) );

			// Free the old array...
			::free( m_lpDigits );

			// Save the new info...
			m_lpDigits = lp_newStorage;
		}

		else
		{
			::memmove( m_lpDigits + lShiftDigits, m_lpDigits, m_lDigitsLength * sizeof(long) );
			::memset( m_lpDigits, 0, lShiftDigits * sizeof(long) );
		}

		m_lDigitsLength += lShiftDigits;
	}

	return kNoError;
}


long CLZNumber::ShiftRight( long lShiftDigits )
{
	if ( lShiftDigits > 0 )
	{
		if ( lShiftDigits < m_lDigitsLength )
		{
			m_lDigitsLength -= lShiftDigits;
			::memmove( m_lpDigits, m_lpDigits + lShiftDigits, m_lDigitsLength * sizeof(long) );
			::memset( m_lpDigits + m_lDigitsLength, 0, (m_lDigitsArrayLength - m_lDigitsLength) * sizeof(long) );

			if ( lShiftDigits > m_lRadix )
				m_lRadix = 0;
			else
				m_lRadix -= lShiftDigits;
		}

		else
		{
			m_lDigitsLength = 1;
			::memset( m_lpDigits, 0, m_lDigitsArrayLength * sizeof(long) );
			m_lRadix = 0;
		}
	}

	return kNoError;
}


long CLZNumber::StripFraction()
{
	long	l_Result = kNoError;

	Round( RM_round_truncate, RD_digit_zero );

	TrimTrailing();

	CheckForZero();

	return l_Result;
}


long CLZNumber::StripInteger()
{
	long	l_Result = kNoError;

	// Clear the integer digits, m_lRadix .. m_lDigitsLength-1, working
	// downward; starting on m_lDigitsLength wrote past the last digit...
	long	l_compareDigits = m_lDigitsLength - 1;

	while ( l_compareDigits >= m_lRadix )
	{
		m_lpDigits[ l_compareDigits ] = 0;

		if ( l_compareDigits == 0 )
			break;

		l_compareDigits--;
	}

	TrimLeading();

	CheckForZero();

	return l_Result;
}


long CLZNumber::TrimAll()
{
	long	l_Result = kNoError;

	if ( (l_Result = TrimLeading()) == kNoError )
		l_Result = TrimTrailing();

	return l_Result;
}


long CLZNumber::TrimLeading()
{
	long	l_leadingZeroes = 0;
	long	l_testIndex = 0;

	// Reset...
	l_testIndex = m_lDigitsLength - 1;

	// Check for leading zeroes...
	while ( l_testIndex >= 0 )
	{
		if ( m_lpDigits[ l_testIndex ] > 0 )
			break;

		l_leadingZeroes++;
		l_testIndex--;
	}

	// Only trim the whole number part...
	if ( l_leadingZeroes > ( m_lDigitsLength - m_lRadix ) )
		l_leadingZeroes = m_lDigitsLength - m_lRadix;

	// Trim...
	if ( l_leadingZeroes > 0 )
	{
		// Set new digits length...
		m_lDigitsLength = m_lDigitsLength - l_leadingZeroes;

		// Clear digits...
		::memset( m_lpDigits + m_lDigitsLength, 0, (m_lDigitsArrayLength - m_lDigitsLength) * sizeof(long) );
	}

	return kNoError;
}


long CLZNumber::TrimTrailing()
{
	long	l_trailingZeroes = 0;
	long	l_testIndex = 0;

	// Check for trailing zeroes...
	while ( l_testIndex < m_lDigitsLength )
	{
		if ( m_lpDigits[ l_testIndex ] > 0 )
			break;

		l_trailingZeroes++;
		l_testIndex++;
	}

	// Only trim the fractional part...
	if ( l_trailingZeroes > m_lRadix )
		l_trailingZeroes = m_lRadix;

	// Trim...
	if ( l_trailingZeroes > 0 )
	{
		// Set new digits length...
		m_lDigitsLength -= l_trailingZeroes;
		m_lRadix -= l_trailingZeroes;

		// Copy digits...
		::memmove( m_lpDigits, m_lpDigits + l_trailingZeroes, m_lDigitsLength * sizeof(long) );
		::memset( m_lpDigits + m_lDigitsLength, 0, (m_lDigitsArrayLength - m_lDigitsLength) * sizeof(long) );
	}

	return kNoError;
}
