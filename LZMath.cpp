/////////////////////////////////////////////////////////////////
//
//	LZMath.cpp
//
//	Part of Addaptic
//
//	Copyright (c) 2006-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#include <memory.h>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <array>
#include <climits>
#include <new>

#include "LZMath.h"

/////////////////////////////////////////////////////////////////
//	Globals

CLZNumber*	gRecipFact[ kRecipFactMax + 1 ];
CLZNumber*	gInvTrigA[ 23 ];
CLZNumber*	gBernoulli[ 30 ];

// ln of each base and the base itself, indexed by base, 2 through 16.
//	They were five named pairs and two five-way switches, whose default
//	sent the other ten bases to ln 10...
CLZNumber*	gLogBaseE[ 17 ];
CLZNumber*	gBaseN[ 17 ];

// ...and the text gLogBaseE is made from, whose length tells
//	lzi_log_base_in_base how far gLogBaseE serves...
static const char* const	kLogBaseEText[ 17 ] =
{
	nullptr,	nullptr,	kLog2BaseE,		kLog3BaseE,
	kLog4BaseE,	kLog5BaseE,	kLog6BaseE,		kLog7BaseE,
	kLog8BaseE,	kLog9BaseE,	kLog10BaseE,	kLog11BaseE,
	kLog12BaseE, kLog13BaseE, kLog14BaseE,	kLog15BaseE,
	kLog16BaseE
};

// Every entry of the tables above, and what each table entry or text
//	constant has been converted to in other bases -- keyed by the source,
//	the base, the working precision and the places asked for. See
//	lzi_table_in_base and lzi_constant_in_base...
static std::set<const CLZNumber*>	gTabled;

// ...one cache of conversions per thread, so no two threads ever share
//	one, and nothing has to be locked to read it. A thread's
//	conversions are deleted when it ends, and by
//	lzu_takedown for the thread that calls it...
class CLZConversions : public std::map< std::tuple<const void*, long, long, long>, CLZNumber* >
{
public:

	~CLZConversions()
	{
		Empty();
	}

	void Empty()
	{
		for ( auto& converted : *this )
			delete converted.second;

		clear();
	}
};

static thread_local CLZConversions	gTabledInBase;

/*
	**A computed constant this wide or wider**, for the three that are
	worth it: pi, ln(base) and, since 2026-10-15, e -- each past its text
	a computation as long as the width. Asked at a width this thread has not computed, they round down
	the narrowest wider one they hold rather than computing again -- see
	lzi_cached_rounded. The entry for pSource in lBase at lPlaces or more,
	its fourth field lCode or, for lCode < 0, anything; nullptr if none.
	The map orders on source, base, places and code, so the first match
	from the lower bound is the narrowest.
*/

static CLZNumber* lzi_cached_wider( const void* pSource, long lBase, long lPlaces, long lCode )
{
	for ( auto it = gTabledInBase.lower_bound( std::make_tuple( pSource, lBase, lPlaces, LONG_MIN ) );
		  (it != gTabledInBase.end()) && (std::get<0>( it->first ) == pSource) && (std::get<1>( it->first ) == lBase);
		  ++it )
	{
		if ( (lCode < 0) || (std::get<3>( it->first ) == lCode) )
			return it->second;
	}

	return nullptr;
}

// ...and that one rounded to lPlaces and kept under key, so the width is
//	held exactly next time. In the default mode, whatever the caller's: a
//	cached constant is every caller's...
static CLZNumber* lzi_cached_rounded( const std::tuple<const void*, long, long, long>& key, CLZNumber* pWider, long lPlaces )
{
	auto	p_value = new CLZNumber;

	p_value->SetBase( pWider->m_lBase, false );
	p_value->Assign( pWider );
	p_value->SetPrecision( lPlaces );
	p_value->Round( kDefaultRoundingMode );

	return gTabledInBase.emplace( key, p_value ).first->second;
}

// The places after a constant's point, which say how far its text
//	serves. Each caller keeps its answer in a local static, whose
//	initialiser runs once however many threads arrive...
static long lzi_text_places( const char* szText )
{
	const char*	sz_point = ::strchr( szText, '.' );

	return (sz_point != nullptr) ? static_cast<long>( ::strlen( sz_point + 1 ) ) : 0;
}


/*
	**A guard is a count of decimal digits, and the work is in places of
	the operand's base.** Every guard in this file was measured in base
	ten and added to the working precision as places, so four places of
	guard were four bits in base two and a digit and a quarter in base
	seven. At the default precision, across fourteen functions and 50
	random arguments each, bases two to nine missed correct rounding 193
	times where they miss 60 now; arcsine in base two went from 18 of 50
	to none.

	So a guard becomes the places that hold at least that many digits.
	Never fewer than the count itself: in bases eleven to sixteen a place
	is already more than a digit, and those bases, ten with them, keep
	exactly the guards they had.
*/

static long lzi_guard_places( long lDigits, long lBase )
{
	long	l_places = CLZNumber::EquivalentPlaces( lDigits, 10, lBase );

	return (l_places > lDigits) ? l_places : lDigits;
}

// The means and the AGM work in their first value's base, and raise the
//	width before the list has been checked -- see lzi_mean_prepare...
static long lzi_mean_base( CLZNumber** values, long count )
{
	if ( (values == nullptr) || (count < 1) || (values[ 0 ] == nullptr) )
		return 10;

	return values[ 0 ]->m_lBase;
}


/////////////////////////////////////////////////////////////////
//	Utility functions

long lzu_setup()
{
	gRecipFact[ 0 ] = new CLZNumber( kRecip00Fact );
	gRecipFact[ 1 ] = new CLZNumber( kRecip01Fact );
	gRecipFact[ 2 ] = new CLZNumber( kRecip02Fact );
	gRecipFact[ 3 ] = new CLZNumber( kRecip03Fact );
	gRecipFact[ 4 ] = new CLZNumber( kRecip04Fact );
	gRecipFact[ 5 ] = new CLZNumber( kRecip05Fact );
	gRecipFact[ 6 ] = new CLZNumber( kRecip06Fact );
	gRecipFact[ 7 ] = new CLZNumber( kRecip07Fact );
	gRecipFact[ 8 ] = new CLZNumber( kRecip08Fact );
	gRecipFact[ 9 ] = new CLZNumber( kRecip09Fact );
	gRecipFact[ 10 ] = new CLZNumber( kRecip10Fact );
	gRecipFact[ 11 ] = new CLZNumber( kRecip11Fact );
	gRecipFact[ 12 ] = new CLZNumber( kRecip12Fact );
	gRecipFact[ 13 ] = new CLZNumber( kRecip13Fact );
	gRecipFact[ 14 ] = new CLZNumber( kRecip14Fact );
	gRecipFact[ 15 ] = new CLZNumber( kRecip15Fact );
	gRecipFact[ 16 ] = new CLZNumber( kRecip16Fact );
	gRecipFact[ 17 ] = new CLZNumber( kRecip17Fact );
	gRecipFact[ 18 ] = new CLZNumber( kRecip18Fact );
	gRecipFact[ 19 ] = new CLZNumber( kRecip19Fact );
	gRecipFact[ 20 ] = new CLZNumber( kRecip20Fact );
	gRecipFact[ 21 ] = new CLZNumber( kRecip21Fact );
	gRecipFact[ 22 ] = new CLZNumber( kRecip22Fact );
	gRecipFact[ 23 ] = new CLZNumber( kRecip23Fact );
	gRecipFact[ 24 ] = new CLZNumber( kRecip24Fact );
	gRecipFact[ 25 ] = new CLZNumber( kRecip25Fact );
	gRecipFact[ 26 ] = new CLZNumber( kRecip26Fact );
	gRecipFact[ 27 ] = new CLZNumber( kRecip27Fact );
	gRecipFact[ 28 ] = new CLZNumber( kRecip28Fact );
	gRecipFact[ 29 ] = new CLZNumber( kRecip29Fact );
	gRecipFact[ 30 ] = new CLZNumber( kRecip30Fact );
	gRecipFact[ 31 ] = new CLZNumber( kRecip31Fact );
	gRecipFact[ 32 ] = new CLZNumber( kRecip32Fact );
	gRecipFact[ 33 ] = new CLZNumber( kRecip33Fact );
	gRecipFact[ 34 ] = new CLZNumber( kRecip34Fact );
	gRecipFact[ 35 ] = new CLZNumber( kRecip35Fact );
	gRecipFact[ 36 ] = new CLZNumber( kRecip36Fact );
	gRecipFact[ 37 ] = new CLZNumber( kRecip37Fact );
	gRecipFact[ 38 ] = new CLZNumber( kRecip38Fact );
	gRecipFact[ 39 ] = new CLZNumber( kRecip39Fact );
	gRecipFact[ 40 ] = new CLZNumber( kRecip40Fact );

	gInvTrigA[ 0 ] = new CLZNumber( kInvTrigA00 );
	gInvTrigA[ 1 ] = new CLZNumber( kInvTrigA01 );
	gInvTrigA[ 2 ] = new CLZNumber( kInvTrigA02 );
	gInvTrigA[ 3 ] = new CLZNumber( kInvTrigA03 );
	gInvTrigA[ 4 ] = new CLZNumber( kInvTrigA04 );
	gInvTrigA[ 5 ] = new CLZNumber( kInvTrigA05 );
	gInvTrigA[ 6 ] = new CLZNumber( kInvTrigA06 );
	gInvTrigA[ 7 ] = new CLZNumber( kInvTrigA07 );
	gInvTrigA[ 8 ] = new CLZNumber( kInvTrigA08 );
	gInvTrigA[ 9 ] = new CLZNumber( kInvTrigA09 );
	gInvTrigA[ 10 ] = new CLZNumber( kInvTrigA10 );
	gInvTrigA[ 11 ] = new CLZNumber( kInvTrigA11 );
	gInvTrigA[ 12 ] = new CLZNumber( kInvTrigA12 );
	gInvTrigA[ 13 ] = new CLZNumber( kInvTrigA13 );
	gInvTrigA[ 14 ] = new CLZNumber( kInvTrigA14 );
	gInvTrigA[ 15 ] = new CLZNumber( kInvTrigA15 );
	gInvTrigA[ 16 ] = new CLZNumber( kInvTrigA16 );
	gInvTrigA[ 17 ] = new CLZNumber( kInvTrigA17 );
	gInvTrigA[ 18 ] = new CLZNumber( kInvTrigA18 );
	gInvTrigA[ 19 ] = new CLZNumber( kInvTrigA19 );
	gInvTrigA[ 20 ] = new CLZNumber( kInvTrigA20 );
	gInvTrigA[ 21 ] = new CLZNumber( kInvTrigA21 );
	gInvTrigA[ 22 ] = new CLZNumber( kInvTrigA22 );

	gBernoulli[ 0 ] = new CLZNumber( kBernoulli00 );
	gBernoulli[ 1 ] = new CLZNumber( kBernoulli01 );
	gBernoulli[ 2 ] = new CLZNumber( kBernoulli02 );
	gBernoulli[ 3 ] = new CLZNumber( kBernoulli03 );
	gBernoulli[ 4 ] = new CLZNumber( kBernoulli04 );
	gBernoulli[ 5 ] = new CLZNumber( kBernoulli05 );
	gBernoulli[ 6 ] = new CLZNumber( kBernoulli06 );
	gBernoulli[ 7 ] = new CLZNumber( kBernoulli07 );
	gBernoulli[ 8 ] = new CLZNumber( kBernoulli08 );
	gBernoulli[ 9 ] = new CLZNumber( kBernoulli09 );
	gBernoulli[ 10 ] = new CLZNumber( kBernoulli10 );
	gBernoulli[ 11 ] = new CLZNumber( kBernoulli11 );
	gBernoulli[ 12 ] = new CLZNumber( kBernoulli12 );
	gBernoulli[ 13 ] = new CLZNumber( kBernoulli13 );
	gBernoulli[ 14 ] = new CLZNumber( kBernoulli14 );
	gBernoulli[ 15 ] = new CLZNumber( kBernoulli15 );
	gBernoulli[ 16 ] = new CLZNumber( kBernoulli16 );
	gBernoulli[ 17 ] = new CLZNumber( kBernoulli17 );
	gBernoulli[ 18 ] = new CLZNumber( kBernoulli18 );
	gBernoulli[ 19 ] = new CLZNumber( kBernoulli19 );
	gBernoulli[ 20 ] = new CLZNumber( kBernoulli20 );
	gBernoulli[ 21 ] = new CLZNumber( kBernoulli21 );
	gBernoulli[ 22 ] = new CLZNumber( kBernoulli22 );
	gBernoulli[ 23 ] = new CLZNumber( kBernoulli23 );
	gBernoulli[ 24 ] = new CLZNumber( kBernoulli24 );
	gBernoulli[ 25 ] = new CLZNumber( kBernoulli25 );
	gBernoulli[ 26 ] = new CLZNumber( kBernoulli26 );
	gBernoulli[ 27 ] = new CLZNumber( kBernoulli27 );
	gBernoulli[ 28 ] = new CLZNumber( kBernoulli28 );
	gBernoulli[ 29 ] = new CLZNumber( kBernoulli29 );

	for ( long l_base = 2; l_base <= 16; l_base++ )
	{
		gLogBaseE[ l_base ] = new CLZNumber( kLogBaseEText[ l_base ] );
		gBaseN[ l_base ] = new CLZNumber( l_base );
	}

	for ( CLZNumber* p_entry : gRecipFact )
		gTabled.insert( p_entry );

	for ( CLZNumber* p_entry : gInvTrigA )
		gTabled.insert( p_entry );

	for ( CLZNumber* p_entry : gBernoulli )
		gTabled.insert( p_entry );

	for ( long l_base = 2; l_base <= 16; l_base++ )
	{
		gTabled.insert( gLogBaseE[ l_base ] );
		gTabled.insert( gBaseN[ l_base ] );
	}

	return kNoError;
}


static void lzi_log_agm_forget();

long lzu_takedown()
{
	short i;

	for ( i=0; i<=kRecipFactMax; i++ )
	{
		delete gRecipFact[ i ];
	}

	for ( i=0; i<23; i++ )
	{
		delete gInvTrigA[ i ];
	}

	for ( i=0; i<30; i++ )
	{
		delete gBernoulli[ i ];
	}

	for ( i=2; i<17; i++ )
	{
		delete gLogBaseE[ i ];
		delete gBaseN[ i ];
	}

	gTabledInBase.Empty();
	gTabled.clear();

	// What lz_log_e knew of the AGM's constants went with them...
	lzi_log_agm_forget();

	return kNoError;
}


/////////////////////////////////////////////////////////////////
//	Internal functions

/*
	Keep the first failure of a sequence rather than the last. These
	functions chain several operations into one answer, and assigning
	each return over the previous one throws away everything but the
	final step -- lz_modulo caught kDivideByZero and then overwrote it
	with the success of the lz_integer that followed.
*/

/*
	Report an argument the function has no answer for, the way
	lz_divide reports division by zero: zero the result, flag what the
	answer would have been, and return a code as well, so a caller
	sees it whether it reads flags or return values.

	kFailed follows lz_arcsine, which already reports an argument
	outside a function's domain that way. An argument on the edge of
	the domain, where the limit exists but cannot be written down --
	the logarithm of zero -- is kUndefinedResult with kfInfinity
	instead, following lz_divide.
*/

static long lzi_domain_failure( CLZNumber* leftOperand, CLZNumber* targetNumber,
								unsigned long lFlag, long lCode )
{
	CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( static_cast<long>(0) );
		p_target->SetFlag( lFlag );
	}

	return lCode;
}


/*
	Carry infinity and not-a-number through an operation, so a failure
	deep inside a calculation is still visible in its answer. Before
	this, infinity plus five was five, unflagged, and every trace of
	the division that produced the infinity was gone.

	An infinity here is stored as zero with a flag, and so carries no
	sign. That decides two of the rules below: a sum or difference of
	two infinities cannot be told from a cancelling one, and NaN is
	the honest answer rather than a guess at which.

	Propagation is quiet. The operation that first failed returned a
	code at the time; an operation merely handling the result of one
	is not itself a failure, which is how a quiet NaN behaves
	elsewhere. The flag is the carrier.
*/

/*
	The same for an operation with one operand. A not-a-number always
	begets one. An infinity depends on the function: this engine holds
	an infinity as zero with a flag and so has no sign for it, which
	only matters where the two signs would lead to different answers.

	The hyperbolic sine and cosine, and the cube root, are infinite
	whichever sign went in, so infinity is a true answer. The
	exponential is not -- e to the plus infinity is infinite and e to
	the minus infinity is nought -- and neither is the hyperbolic
	tangent, which is plus or minus one. Those, and the square root,
	the factorial and the logarithm, have no answer that does not
	depend on the sign, so they give a not-a-number rather than a
	guess. The same reasoning as one infinity added to another.
*/

/*
	A comparison has no answer when either side is a not-a-number, and
	none either when a side is infinite, because an infinity is held
	here as zero with a flag and so has no sign: it might be above
	everything or below it. kUnordered says so, rather than picking.

	This applies to the lz_ entry points only. CLZNumber::Compare and
	CompareABS are left alone -- the engine's own algorithms lean on
	them constantly and expect three answers, not four.
*/

static bool lzi_is_unordered( CLZNumber* leftOperand, CLZNumber* rightOperand )
{
	return leftOperand->FlagIsOn( kfQNaN ) || rightOperand->FlagIsOn( kfQNaN ) ||
		   leftOperand->FlagIsOn( kfInfinity ) || rightOperand->FlagIsOn( kfInfinity );
}


/*
	Two operands, and the answer is always a not-a-number: quantize
	sets one value's exponent from another's, and neither an infinity
	nor a not-a-number has an exponent to give or to take.
*/

static bool lzi_propagate_as_nan( CLZNumber* leftOperand, CLZNumber* rightOperand,
								  CLZNumber* targetNumber )
{
	if ( !(leftOperand->FlagIsOn( kfQNaN ) || rightOperand->FlagIsOn( kfQNaN ) ||
		   leftOperand->FlagIsOn( kfInfinity ) || rightOperand->FlagIsOn( kfInfinity )) )
		return false;

	bool	b_cameFromDivideByZero = leftOperand->FlagIsOn( kfDiv0 ) ||
									 rightOperand->FlagIsOn( kfDiv0 );

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( static_cast<long>(0) );
		p_target->SetFlag( kfQNaN );

		if ( b_cameFromDivideByZero )
			p_target->SetFlag( kfDiv0 );
	}

	return true;
}


/*
	Raising to a power, where more survives than elsewhere. An
	infinite base still has no sign, but the exponent does, and that
	is enough: an infinity to any positive power is infinite whichever
	sign it had, to any negative power is nought, and to the nought
	power is one by the usual convention.

	An infinite exponent is the other way about and has no answer at
	all -- two to the plus infinity is infinite, and to the minus
	infinity is nought.
*/

static bool lzi_propagate_power( CLZNumber* leftOperand, CLZNumber* rightOperand,
								 CLZNumber* targetNumber )
{
	bool	b_leftNaN = leftOperand->FlagIsOn( kfQNaN );
	bool	b_rightNaN = rightOperand->FlagIsOn( kfQNaN );
	bool	b_leftInfinite = leftOperand->FlagIsOn( kfInfinity );
	bool	b_rightInfinite = rightOperand->FlagIsOn( kfInfinity );

	if ( !(b_leftNaN || b_rightNaN || b_leftInfinite || b_rightInfinite) )
		return false;

	unsigned long	l_resultFlag = kfQNaN;
	long			l_value = 0;

	// Only the base is infinite, so the exponent can be read...
	if ( !(b_leftNaN || b_rightNaN || b_rightInfinite) )
	{
		rightOperand->CheckForZero();

		if ( rightOperand->IsZero() )
		{
			l_resultFlag = 0;
			l_value = 1;
		}
		else if ( rightOperand->IsPositive() )
			l_resultFlag = kfInfinity;
		else
			l_resultFlag = 0;
	}

	bool	b_cameFromDivideByZero = leftOperand->FlagIsOn( kfDiv0 ) ||
									 rightOperand->FlagIsOn( kfDiv0 );

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( l_value );

		if ( l_resultFlag != 0 )
			p_target->SetFlag( l_resultFlag );

		if ( b_cameFromDivideByZero )
			p_target->SetFlag( kfDiv0 );
	}

	return true;
}


/*
	The edges x^y has wherever it is computed. lz_power, lz_exp_n and
	lz_raise_n each met them their own way and agreed at none: a
	negative base to a fraction was a flagged not-a-number reporting
	success from one, a refusal from another and a silently truncated
	exponent from the third; nought to the minus two was infinite,
	refused, or infinite again but reporting success; and lz_exp_n
	refused every negative and every zero base, (-8)^3 and 0^2
	included, because it took a logarithm of the base first.

	So there is one set of rules, here, and all three ask it first:

	- a not-a-number or an infinity is carried by lzi_propagate_power's
	  rules, which also keeps a flagged zero from being read below as
	  an honest nought;
	- a whole-exponent function refuses a fraction, as lz_factorial
	  refuses 2.5 -- lz_power is the function that takes one;
	- nought to a positive power is nought, and to the nought power
	  one, by the convention lz_power already kept;
	- nought to a negative power is the reciprocal's pole, and is
	  answered by dividing one by nought, so its code and flags are
	  lz_divide's and cannot drift from them;
	- a negative base to a fraction has no real answer, and is refused
	  the way every other domain failure is. To a whole power it has
	  one, signed by the exponent's parity, and that is left to the
	  caller to compute.

	Answers true when it has answered, with the code in *plResult.
*/

static bool lzi_power_edges( CLZNumber* pBase, CLZNumber* pExponent, CLZNumber* pTarget,
							 bool bWholeExponent, long* plResult )
{
	*plResult = kNoError;

	if ( lzi_propagate_power( pBase, pExponent, pTarget ) )
		return true;

	if ( bWholeExponent && !(pExponent->IsInteger()) )
	{
		*plResult = lzi_domain_failure( pBase, pTarget, kfQNaN, kFailed );
		return true;
	}

	pBase->CheckForZero();
	pExponent->CheckForZero();

	CLZNumber*	p_target = (pTarget != nullptr) ? pTarget : pBase;

	if ( pBase->IsZero() )
	{
		if ( pExponent->IsNegative() )
		{
			CLZNumber	one;
			CLZNumber	zero;

			one.SetBase( pBase->m_lBase );
			one.Assign( static_cast<long>(1) );
			zero.Assign( pBase );

			*plResult = lz_divide( &one, &zero, p_target );
		}
		else if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>( pExponent->IsZero() ? 1 : 0 ) );

		return true;
	}

	if ( pBase->IsNegative() && !(pExponent->IsInteger()) )
	{
		*plResult = lzi_domain_failure( pBase, pTarget, kfQNaN, kFailed );
		return true;
	}

	return false;
}


const long	kUnaryStaysInfinite = 0;
const long	kUnaryNeedsTheSign = 1;

static bool lzi_propagate_unary( CLZNumber* leftOperand, CLZNumber* targetNumber, long lRule )
{
	bool	b_notANumber = leftOperand->FlagIsOn( kfQNaN );
	bool	b_infinite = leftOperand->FlagIsOn( kfInfinity );

	if ( !(b_notANumber || b_infinite) )
		return false;

	unsigned long	l_resultFlag = (b_notANumber || (lRule == kUnaryNeedsTheSign))
									? kfQNaN : kfInfinity;

	bool	b_cameFromDivideByZero = leftOperand->FlagIsOn( kfDiv0 );

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( static_cast<long>(0) );
		p_target->SetFlag( l_resultFlag );

		if ( b_cameFromDivideByZero )
			p_target->SetFlag( kfDiv0 );
	}

	return true;
}


const long	kOpAdd = 0;
const long	kOpSubtract = 1;
const long	kOpMultiply = 2;
const long	kOpDivide = 3;

static bool lzi_propagate_flags( CLZNumber* leftOperand, CLZNumber* rightOperand,
								 CLZNumber* targetNumber, long lOperation )
{
	bool	b_leftNaN = leftOperand->FlagIsOn( kfQNaN );
	bool	b_rightNaN = rightOperand->FlagIsOn( kfQNaN );
	bool	b_leftInfinite = leftOperand->FlagIsOn( kfInfinity );
	bool	b_rightInfinite = rightOperand->FlagIsOn( kfInfinity );

	// Nothing to carry: let the operation proceed as usual...
	if ( !(b_leftNaN || b_rightNaN || b_leftInfinite || b_rightInfinite) )
		return false;

	// IsZero reports a flag; ask for it to be brought up to date
	//	before the two rules below lean on it...
	leftOperand->CheckForZero();
	rightOperand->CheckForZero();

	unsigned long	l_resultFlag = 0;

	if ( b_leftNaN || b_rightNaN )
		l_resultFlag = kfQNaN;

	else if ( b_leftInfinite && b_rightInfinite )
	{
		// Unsigned, so only the product is knowable...
		l_resultFlag = (lOperation == kOpMultiply) ? kfInfinity : kfQNaN;
	}

	else if ( b_leftInfinite )
	{
		if ( (lOperation == kOpMultiply) && rightOperand->IsZero() )
			l_resultFlag = kfQNaN;
		else
			l_resultFlag = kfInfinity;
	}

	else	// only the right operand is infinite
	{
		if ( lOperation == kOpDivide )
			l_resultFlag = 0;				// something finite over infinity is zero
		else if ( (lOperation == kOpMultiply) && leftOperand->IsZero() )
			l_resultFlag = kfQNaN;
		else
			l_resultFlag = kfInfinity;
	}

	bool	b_cameFromDivideByZero = leftOperand->FlagIsOn( kfDiv0 ) ||
									 rightOperand->FlagIsOn( kfDiv0 );

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( static_cast<long>(0) );

		if ( l_resultFlag != 0 )
			p_target->SetFlag( l_resultFlag );

		// Where the trouble started travels with it...
		if ( b_cameFromDivideByZero )
			p_target->SetFlag( kfDiv0 );
	}

	return true;
}


/*
	**A table entry is converted into another base once per width, and
	kept.** The tables are base ten, and every use of an entry by a
	value in another base converted it afresh -- 33 conversions of a
	fifty-digit constant for each exponential in base two, and in base
	two SetBaseForWork keeps the entry's own 36 places, so each came out
	120 binary places wide for work needing forty. That, and not the
	arithmetic, is where the time went: base-two exp took 870 us against
	76 with the tables converted once, arcsine 387 against 51. It had
	always been so for the exponential and the sine; working every
	function in its operand's base (see CLZAnswerBase) brought the
	rest of the engine to the tables and made it show.

	**Keyed by the width as well as the base, so a kept conversion is
	exactly the one that would have been made** -- the same entry, base
	and working precision go through the same SetBaseForWork. It can
	change how long an answer takes and never what it is, and the whole
	suite and the hex map are byte-identical with it and without.
	Nothing writes to an entry: they are right operands, or sources of
	a copy. The cache is emptied by lzu_takedown.
*/

static CLZNumber* lzi_table_in_base( CLZNumber* pEntry, long lBase )
{
	if ( pEntry->m_lBase == lBase )
		return pEntry;

	auto	key = std::make_tuple( static_cast<const void*>( pEntry ), lBase, glWorkingPrecision, glWorkingPrecision );
	auto	found = gTabledInBase.find( key );

	if ( found != gTabledInBase.end() )
		return found->second;

	auto	p_converted = new CLZNumber;

	p_converted->Assign( pEntry );
	p_converted->SetBaseForWork( lBase );

	gTabledInBase.emplace( key, p_converted );

	return p_converted;
}


/*
	**The constants a function makes from text are kept the same way.**
	pi and its fractions, e and ln(2 pi)/2 are written as text in
	LZMath.h and made afresh by each call, and a call in another base
	converted each one afresh as well -- once a division in lz_log_e's
	reduction loop. This replaces the value, made from szConstant in
	base ten at the current working precision, with its conversion into
	lBase at lPlaces places of that base, kept from the first time. In
	base ten it does nothing, so a function's base-ten path is exactly
	as it was written.

	**Convert before narrowing, never after.** Several functions narrow
	a 120-place constant to the working precision for speed. Narrowed
	first, it held the working precision counted as decimal places, and
	a hex operand met a pi of forty digits where it needed 48.
*/

static void lzi_constant_in_base( CLZNumber* pValue, const char* szConstant, long lBase, long lPlaces )
{
	if ( pValue->m_lBase == lBase )
		return;

	auto	key = std::make_tuple( static_cast<const void*>( szConstant ), lBase, glWorkingPrecision, lPlaces );
	auto	found = gTabledInBase.find( key );

	if ( found == gTabledInBase.end() )
	{
		auto	p_converted = new CLZNumber( szConstant );

		p_converted->SetBaseForWork( lBase, lPlaces );

		found = gTabledInBase.emplace( key, p_converted ).first;
	}

	pValue->Assign( found->second );
}


/*
	**e as wide as it is asked for.** kBaseE holds 120 places, and a
	caller asking for more met them as a ceiling. e^n wants its own
	size plus 36 places, so from about e^200 its answer came out short
	-- 118 digits right of the 123 e^200 asks for, 116 of e^50000's
	21,751 -- while lz_raise_e's widening worked at every place and
	filled none past the constant's.

	**Past the constant's width in lBase, e is computed**, by binary
	splitting the series: with Q(a,b) = (a+1)(a+2)...b and P(a,b) the
	sum over k in (a,b] of Q(k,b), the sum of 1/k! for k from 1 to K is
	P(0,K)/Q(0,K). Both are exact integers built by halving the range,
	so the work is a tree of multiplications and one division at the
	end, and all of it in lBase, so nothing is converted. K is taken
	where 1/K! falls under a unit of the guard place.

	A computation is kept in the constants' cache, keyed on the base
	and the places, so the cache stays transparent. **Within the
	constant's width nothing changes**: e comes from the text exactly
	as it did.
*/

static const char	kComputedE[] = "e, computed";

const long	kComputedEGuard = 4;

static void lzi_e_split( long lA, long lB, long lBase, CLZNumber* pP, CLZNumber* pQ )
{
	if ( (lB - lA) == 1 )
	{
		pP->SetBase( lBase, false );
		pP->Assign( 1L );
		pQ->SetBase( lBase, false );
		pQ->Assign( lB );

		return;
	}

	long		l_middle = (lA + lB) / 2;
	CLZNumber	p_left, q_left, p_right, q_right;

	lzi_e_split( lA, l_middle, lBase, &p_left, &q_left );
	lzi_e_split( l_middle, lB, lBase, &p_right, &q_right );

	// P(a,b) = P(a,m) Q(m,b) + P(m,b), and Q(a,b) = Q(a,m) Q(m,b). All
	//	integers, so the products are exact at any working precision...
	lz_multiply( &p_left, &q_right, pP );
	lz_add( pP, &p_right );
	lz_multiply( &q_left, &q_right, pQ );
}

static void lzi_e_in_base( CLZNumber* pValue, long lBase, long lPlaces )
{
	static const long	sl_constantPlaces = lzi_text_places( kBaseE );

	if ( lPlaces <= CLZNumber::EquivalentPlaces( sl_constantPlaces, 10, lBase ) )
	{
		pValue->SetBase( 10, false );
		pValue->Assign( kBaseE );
		lzi_constant_in_base( pValue, kBaseE, lBase, lPlaces );

		return;
	}

	auto	key = std::make_tuple( static_cast<const void*>( kComputedE ), lBase, lPlaces, lPlaces );
	auto	found = gTabledInBase.find( key );

	// A wider one is rounded, not computed again, as pi's and ln(base)'s
	//	are -- see lzi_cached_wider. Since 2026-10-15, performance item 8 of
	//	that day's review: e was kept only at the width it was made at, so a
	//	narrower width computed it afresh -- 2.8 times a first call in base
	//	two...
	if ( found == gTabledInBase.end() )
	{
		if ( CLZNumber* p_wider = lzi_cached_wider( kComputedE, lBase, lPlaces, -1 ) )
		{
			pValue->Assign( lzi_cached_rounded( key, p_wider, lPlaces ) );
			return;
		}
	}

	if ( found == gTabledInBase.end() )
	{
		double	d_wanted = static_cast<double>( lPlaces + kComputedEGuard ) * std::log( static_cast<double>( lBase ) );
		double	d_logFactorial = 0.0;
		long	l_terms = 1;

		while ( d_logFactorial <= d_wanted )
		{
			l_terms++;
			d_logFactorial += std::log( static_cast<double>( l_terms ) );
		}

		CLZNumber	p, q, one;
		auto		p_e = new CLZNumber;

		lzi_e_split( 0, l_terms, lBase, &p, &q );

		// lz_divide takes its width from its target...
		p_e->SetBase( lBase, false );
		p_e->SetPrecision( lPlaces + kComputedEGuard );
		lz_divide( &p, &q, p_e );

		one.SetBase( lBase, false );
		one.Assign( 1L );
		lz_add( p_e, &one );

		p_e->SetPrecision( lPlaces );
		p_e->Round();

		found = gTabledInBase.emplace( key, p_e ).first;
	}

	pValue->Assign( found->second );
}


/*
	**Pi as wide as it is asked for**, as e is. The text constants hold
	120 places, and every circular function reads pi, pi/2, pi/4 or
	2*pi: past that width the tail was nought, so until 2026-09-26 the
	whole family went wrong past 120 places in base ten -- cos at 150
	places missed correct rounding for two arguments in three even
	below pi/4, every large argument's reduction missed, and base two
	was spared only because 300 of its places are 90 digits.

	Within the constant's width nothing changes: the text serves, as it
	always did. Past it, pi comes from the Chudnovsky series by binary
	splitting,

		pi = 426880 sqrt(10005) Q(0,N) / T(0,N)

	P, Q and T exact integers built by halving the range, all in the
	operand's base, about 14.18 digits a term, then one root and one
	division. A computation is kept in the constants' cache, keyed on
	the base and the places, and 2*pi, pi/2 and pi/4 are taken from it
	at the same width.
*/

static const char	kComputedPi[] = "pi, computed";

const long	kComputedPiGuard = 4;

/*
	**Nothing here may need more than a 32-bit long**, which is what a
	long is on Windows. Until 2026-09-29 two things did, and between them
	made every pi computed past its text wrong there, without a word:
	640320^3 / 24 = 10939058860032000 was a long constant, which truncates
	to 495419392, and a leaf's 545140134 a overflows from a = 4. The cube
	is now built once, in the operand's base, from factors that fit, and
	handed to every leaf; the leaf's sum is made in the engine's own
	arithmetic. Found by building the engine with long rewritten to 32
	bits.
*/

const long	kChudnovskyCubeFactors[] = { 26680, 640320, 640320 };		// 640320^3 / 24

static void lzi_pi_split( long lA, long lB, long lBase, CLZNumber* pCube,
						  CLZNumber* pP, CLZNumber* pQ, CLZNumber* pT )
{
	if ( (lB - lA) == 1 )
	{
		CLZNumber	factor;
		CLZNumber	term;

		pP->SetBase( lBase, false );
		pQ->SetBase( lBase, false );
		pT->SetBase( lBase, false );
		factor.SetBase( lBase, false );
		term.SetBase( lBase, false );

		if ( lA == 0 )
		{
			pP->Assign( 1L );
			pQ->Assign( 1L );
		}
		else
		{
			// P = (6a-5)(2a-1)(6a-1), Q = a^3 640320^3 / 24...
			pP->Assign( (6 * lA) - 5 );
			factor.Assign( (2 * lA) - 1 );
			lz_multiply( pP, &factor );
			factor.Assign( (6 * lA) - 1 );
			lz_multiply( pP, &factor );

			pQ->Assign( lA );
			factor.Assign( lA );
			lz_multiply( pQ, &factor );
			lz_multiply( pQ, &factor );
			lz_multiply( pQ, pCube );
		}

		// T = (-1)^a P (13591409 + 545140134 a)...
		factor.Assign( 545140134L );
		term.Assign( lA );
		lz_multiply( &factor, &term );
		term.Assign( 13591409L );
		lz_add( &factor, &term );
		lz_multiply( pP, &factor, pT );

		if ( (lA % 2) != 0 )
			lz_negate( pT );

		return;
	}

	long		l_middle = (lA + lB) / 2;
	CLZNumber	p_left, q_left, t_left, p_right, q_right, t_right, product;

	lzi_pi_split( lA, l_middle, lBase, pCube, &p_left, &q_left, &t_left );
	lzi_pi_split( l_middle, lB, lBase, pCube, &p_right, &q_right, &t_right );

	// P = P1 P2, Q = Q1 Q2, T = T1 Q2 + P1 T2 -- integers, exact...
	lz_multiply( &p_left, &p_right, pP );
	lz_multiply( &q_left, &q_right, pQ );
	lz_multiply( &t_left, &q_right, pT );
	lz_multiply( &p_left, &t_right, &product );
	lz_add( pT, &product );
}

static CLZNumber* lzi_computed_pi( long lBase, long lPlaces )
{
	auto	key = std::make_tuple( static_cast<const void*>( kComputedPi ), lBase, lPlaces, 0L );
	auto	found = gTabledInBase.find( key );

	if ( found != gTabledInBase.end() )
		return found->second;

	// A wider one is rounded, not computed again -- see lzi_cached_wider...
	if ( CLZNumber* p_wider = lzi_cached_wider( kComputedPi, lBase, lPlaces, 0 ) )
		return lzi_cached_rounded( key, p_wider, lPlaces );

	long	l_saved = glWorkingPrecision;
	long	l_width = lPlaces + kComputedPiGuard;

	glWorkingPrecision = l_width;

	double	d_digits = static_cast<double>( l_width ) * std::log10( static_cast<double>( lBase ) ) + 2.0;
	long	l_terms = static_cast<long>( d_digits / 14.18 ) + 2;

	CLZNumber	p, q, t, root, scale, cube;
	auto		p_pi = new CLZNumber;

	// 640320^3 / 24, whole, in lBase, from factors a 32-bit long holds...
	cube.SetBase( lBase, false );
	cube.Assign( kChudnovskyCubeFactors[ 0 ] );
	scale.SetBase( lBase, false );

	for ( long i = 1; i < static_cast<long>( sizeof( kChudnovskyCubeFactors ) / sizeof( kChudnovskyCubeFactors[ 0 ] ) ); i++ )
	{
		scale.Assign( kChudnovskyCubeFactors[ i ] );
		lz_multiply( &cube, &scale );
	}

	lzi_pi_split( 0, l_terms, lBase, &cube, &p, &q, &t );

	root.SetBase( lBase, false );
	root.Assign( 10005L );
	lz_root_2( &root );						// sqrt(10005), to the width

	scale.SetBase( lBase, false );
	scale.Assign( 426880L );
	lz_multiply( &root, &scale );
	lz_multiply( &root, &q );

	// lz_divide takes its width from its target...
	p_pi->SetBase( lBase, false );
	p_pi->SetPrecision( l_width );
	lz_divide( &root, &t, p_pi );

	p_pi->SetPrecision( lPlaces );
	p_pi->Round();

	glWorkingPrecision = l_saved;

	return gTabledInBase.emplace( key, p_pi ).first->second;
}

/*
	pi, 2*pi, pi/2 or pi/4 -- szConstant names which -- in lBase to
	lPlaces places: the text within its width, computed past it. A
	drop-in for lzi_constant_in_base on the four, which is what every
	reader of them called.
*/

static void lzi_pi_family_in_base( CLZNumber* pValue, const char* szConstant, long lBase, long lPlaces )
{
	static const long	sl_constantPlaces = lzi_text_places( kPi );

	// Within the text's width, exactly the call every reader made
	//	before: each makes its value from the same text first...
	if ( lPlaces <= CLZNumber::EquivalentPlaces( sl_constantPlaces, 10, lBase ) )
	{
		lzi_constant_in_base( pValue, szConstant, lBase, lPlaces );

		return;
	}

	long	l_code = (szConstant == kTwoPi) ? 2 : (szConstant == kPiOverTwo) ? 3 : (szConstant == kPiOverFour) ? 4 : 1;
	auto	key = std::make_tuple( static_cast<const void*>( kComputedPi ), lBase, lPlaces, l_code );
	auto	found = gTabledInBase.find( key );

	if ( found == gTabledInBase.end() )
	{
		long		l_saved = glWorkingPrecision;
		long		l_width = lPlaces + kComputedPiGuard;
		CLZNumber*	p_pi = lzi_computed_pi( lBase, l_width );
		auto		p_value = new CLZNumber;
		CLZNumber	factor;

		glWorkingPrecision = l_width;

		p_value->SetBase( lBase, false );
		p_value->SetPrecision( l_width );
		factor.SetBase( lBase, false );

		if ( l_code == 2 )
		{
			factor.Assign( 2L );
			lz_multiply( p_pi, &factor, p_value );
		}
		else if ( l_code == 1 )
			p_value->Assign( p_pi );
		else
		{
			factor.Assign( (l_code == 3) ? 2L : 4L );
			lz_divide( p_pi, &factor, p_value );
		}

		p_value->SetPrecision( lPlaces );
		p_value->Round();

		glWorkingPrecision = l_saved;

		found = gTabledInBase.emplace( key, p_value ).first;
	}

	pValue->Assign( found->second );
}

/*
	Bring two operands into one base, so an operation on a mixed pair
	means something.

	A digit array carries no base, so an operand in another one would
	be read as though it were in this one -- 0x15 and 21 are the same
	number and different digits. Every arithmetic entry point used to
	refuse the pair outright with kBaseMismatch, which was safe and
	unhelpful: Elzed accepts "numbers in any base, or combination of
	bases, at any time", and once Addaptic learned to read 0x15 from
	text the refusal became reachable from ordinary input.

	**The left operand's base wins**, which is what lzi_divide_op has
	always done with its quotient. So 0x15 + 0b1010 is 0x1F and
	0b1010 + 0x15 is 0b11111 -- the same thirty-one, written the way
	the left-hand side was written.

	The right operand is copied rather than converted in place: an
	operand is not the operator's to modify, and callers pass
	constants. The copy is only made when the bases actually differ,
	so the ordinary same-base call pays one comparison.

	It matters that this runs before lzi_propagate_flags, for the same
	reason the refusal did.

	**Nothing is allocated unless the bases actually differ.** The
	first cut took a CLZNumber by value for the copy, whose
	constructor callocs whether it is used or not: that cost the
	addition 41 per cent and the comparison seven times over, on the
	same-base call that is all but every call. Hence the holder --
	two pointers on the stack, and a conversion only when there is
	one to do.
*/

class CLZMatchedBase
{
public:

	CLZMatchedBase( CLZNumber* leftOperand, CLZNumber* rightOperand )
		: m_pConverted( nullptr ), m_pRight( rightOperand )
	{
		if ( leftOperand->m_lBase != rightOperand->m_lBase )
		{
			// A table entry's conversion is kept -- see lzi_table_in_base...
			if ( gTabled.count( rightOperand ) != 0 )
			{
				m_pRight = lzi_table_in_base( rightOperand, leftOperand->m_lBase );
				return;
			}

			m_pConverted = new CLZNumber;

			m_pConverted->Assign( rightOperand );
			m_pConverted->SetBaseForWork( leftOperand->m_lBase );

			m_pRight = m_pConverted;
		}
	}

	~CLZMatchedBase()
	{
		delete m_pConverted;
	}

	CLZNumber*	Operand() const		{ return m_pRight; }

private:

	CLZMatchedBase( const CLZMatchedBase& );
	CLZMatchedBase& operator=( const CLZMatchedBase& );

	CLZNumber*	m_pConverted;
	CLZNumber*	m_pRight;
};


/*
	**Every answer comes back in the left operand's base** -- Rich's
	rule, 2026-09-24, and the one the twelve entry points above already
	kept by construction.

	The others did not. A function builds its answer from temporaries,
	and a fresh CLZNumber is base ten, as are the constants; Assign
	copies the base along with the value, so the answer took whichever
	base its last working value happened to have. The value was right
	and the base was not: for a hexadecimal operand, sine answered in
	base ten, lz_log_n in its right operand's base, lz_power in either
	depending on whether the exponent had a whole part, and a refusal's
	zero in whatever base the target had before. 354 of 1126 answers
	across 29 functions, measured with operands in bases sixteen and
	seven and a target set to base two, in place and not.

	So this is declared first in every function that answers a number
	and is not one of the twelve: it notes the left operand's base on
	the way in and converts the answer to it on the way out, by any
	return, early ones included. An answer already there costs a
	comparison. A constant was not written to, so it is left alone.

	**Declare it before anything else**, or an early return above it
	escapes. The rule is mechanical: every public function taking a
	target carries either this or CLZMatchedBase, and
	TestAnswerBases checks every one of them.
*/

/*
	**A fresh CLZNumber is base ten, and an operation answers in its
	left operand's base** -- so a temporary made fresh and used on the
	left, a one to divide by the operand or a sum to accumulate into,
	took the work into base ten without a word, and so did a base-ten
	constant written on the left. Everything built on that answer ran
	there too, at the working precision counted as decimal places, and
	this guard converted the answer back on the way out, which is what
	hid it: hex answers carried thirty places, base twelve's 34, and
	base two's 120, the decimal's precision re-counted. A census found
	25 such sites in 19 functions.

	So a temporary takes the operand's base before it is used, and a
	constant on the left is converted first -- by SetBaseForWork, or by
	lzi_table_in_base for a table entry. Both are nothing in base ten.
*/

/*
	**A target's rounding settings are its own**, and an operation leaves
	them as it found them. Assign copies the mode and the digit along
	with the value, so until 2026-09-26 a target came back from every
	operation set to whatever its temporary had been -- half-even,
	whatever the caller had asked for -- and a division or a
	multiplication into it never rounded in the target's mode at all.
	This takes them at the start and puts them back at the end.
	RestoreMode puts the mode back early, for the two operations that
	round in it: see lzi_divide_op and lz_multiply.

	**The digit is honoured at the units only**, since 2026-10-01: a
	target set to RD_digit_zero gets a whole number from the four plain
	operations, rounded once from the exact value in its own mode -- see
	lzi_asks_whole. 'Last' and the default are not: the last place is
	one a number already holds, and the default is the width the
	operation answers at anyway.
*/

/*
	**Whether an operation's answer should be a whole number**: its
	target set to round at the units, and no function's work around it
	(OwnDigit answers the default while that is up -- see
	CLZNumber::EnterInternalWork). Taken before anything is assigned to
	the target, since Assign brings the temporary's digit with it.

	Rounding afterwards is no substitute, being a second rounding: a
	product or quotient that lands exactly on a half at its own width
	rounds from there. 1.581138830084189665999446772216359267 squared is
	2.5 + 4.4e-37, 2.5 at 36 places, and a whole number taken from that
	in half-even is 2 where the nearest is 3.

	A target set not to round at all is not rounded at any digit: the
	whole-number division below would otherwise hand back its working
	digits unrounded -- 22/7 came back 3.1421.
*/

static bool lzi_asks_whole( const CLZNumber* pAnswer )
{
	return (pAnswer != nullptr) && (pAnswer->OwnDigit() == RD_digit_zero) &&
		   (pAnswer->OwnMode() != RM_round_noop);
}


class CLZTargetRounding
{
public:

	explicit CLZTargetRounding( CLZNumber* pAnswer )
		: m_pAnswer( pAnswer ), m_rmMode( kDefaultRoundingMode ), m_rdDigit( RD_digit_precision ),
		  m_bWhole( lzi_asks_whole( pAnswer ) )
	{
		if ( m_pAnswer != nullptr )
		{
			m_rmMode = m_pAnswer->m_rmMode;
			m_rdDigit = m_pAnswer->m_rdDigit;
		}
	}

	~CLZTargetRounding()
	{
		if ( m_pAnswer != nullptr )
		{
			m_pAnswer->m_rmMode = m_rmMode;
			m_pAnswer->m_rdDigit = m_rdDigit;
		}
	}

	void RestoreMode()
	{
		if ( m_pAnswer != nullptr )
			m_pAnswer->m_rmMode = m_rmMode;
	}

	// The target asked for a whole number -- see lzi_asks_whole...
	bool Whole() const
	{
		return m_bWhole;
	}

	// An exact answer, just assigned, rounded to a whole number in the
	//	target's mode if the target asked for one. Addition and
	//	subtraction round nowhere else...
	void RoundExactIfWhole()
	{
		if ( m_bWhole && (m_pAnswer != nullptr) )
		{
			RestoreMode();
			m_pAnswer->Round( RD_digit_zero );
		}
	}

private:

	CLZTargetRounding( const CLZTargetRounding& );
	CLZTargetRounding& operator=( const CLZTargetRounding& );

	CLZNumber*		m_pAnswer;
	RoundingMode	m_rmMode;
	RoundingDigit	m_rdDigit;
	bool			m_bWhole;
};


/*
	**The engine's own work, which no caller's rounding settings reach**
	-- see CLZNumber::EnterInternalWork. Every function raises it
	through CLZAnswerBase, and the compound operations raise it
	themselves. The three means need it for the target's settings, the
	caller's values reaching them only as right operands.

	The integer quotient and the three remainders needed it when they
	divided into their target and took a whole number from what the
	target's mode had rounded. Since 2026-09-26 they take it straight
	from the raw quotient in a mode of their own -- see
	lzi_divide_work -- and answer through lzi_divide_work or
	lz_subtract, which keep their target's settings themselves. So
	theirs are defensive now: taking them out fails nothing.

	Only the four plain operations leave it down, and round in their
	target's own mode. It keeps the answer's settings as well.
*/

class CLZInternalWork
{
public:

	explicit CLZInternalWork( CLZNumber* pAnswer )
		: m_target( pAnswer )
	{
		CLZNumber::EnterInternalWork();
	}

	~CLZInternalWork()
	{
		CLZNumber::LeaveInternalWork();
	}

private:

	CLZInternalWork( const CLZInternalWork& );
	CLZInternalWork& operator=( const CLZInternalWork& );

	CLZTargetRounding	m_target;
};


static void lzi_write_to_width( CLZNumber* pAnswer, bool bQuotientWidth );
static bool lzi_exactly( CLZNumber* pValue, long lValue );
static bool lzi_short_of_width( CLZNumber* pAnswer );
static bool lzi_power_is_exactly( CLZNumber* pValue, CLZNumber* pBase, CLZNumber* pExponent );
static void lzi_exact_power( CLZNumber* pBase, long lPower, CLZNumber* pResult );
static void lzi_root_decide( CLZNumber* pRoot, const CLZNumber* pOperand, long lPower, long lRound );
static bool lzi_rounding_in_doubt( const CLZNumber* pValue, long lRound, long lWidth, long lUnits );
static long lzi_leading_exponent( CLZNumber* value );
static bool lzi_is_power_of_base( CLZNumber* pValue );
static void lzi_trim_after_carry( CLZNumber* pAnswer );

/*
	**And written out to its width unless it is exact**, since
	2026-10-01 -- see lzi_write_to_width. A function marks an answer it
	knows to be exact with Exact(), which leaves its spelling as it is,
	or with ExactAndShort(), which takes its trailing noughts off -- for
	an answer reached through other functions, whose own guards will
	have written it out. Any other answer shorter than its width gets
	its trailing noughts as the function returns. QuotientWidth() gives
	a small answer the width the rule for a quotient gives it.
*/

typedef long (*LZUnaryFunction)( CLZNumber*, CLZNumber* );
typedef long (*LZBinaryFunction)( CLZNumber*, CLZNumber*, CLZNumber* );

// Whether pResult, a function's answer worked wider and in doubt at
//	lPlaces, can be rounded there exactly from the operands it was taken
//	of -- a tie to the even kept digit -- with numbers of no more than
//	lMostDigits digits, and if so pResult so rounded. A function of one
//	operand is passed no right one...
typedef bool (*LZTieDecider)( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* pRight, long lPlaces, long lMostDigits );

// How large a decider's numbers may grow: kTieEarlyDigits answer's widths
//	while the retake has rounds left -- cheaper than those rounds, measured
//	-- and kTieMostPlaces once it has none -- see lzi_take_again_wider...
const long	kTieEarlyDigits = 512;
const long	kTieMostPlaces = 2000000;

/*
	**Taken again wider when the guard leaves the rounding in doubt**, since
	2026-10-03. A function that rounds its answer from guard digits cannot
	tell an answer within its own error of a half from the half, and
	rounds it either way: the near-half census found every one of them
	doing it, at depths from 10^-3 of a unit down. One that calls Retake
	right after its CLZAnswerBase, at the outermost call, has the digits
	its final rounding discarded read as it returns -- see
	CLZNumber::LastDiscardedTail -- and when they lie within its doubt
	window of the half it is taken again kRetakeGuardDigits wider. That
	answer decides the rounding when it lies more than a unit of its own
	last place clear of the half; otherwise it is taken wider again, the
	extra doubling each round, and once it has reached the answer's own
	places -- see kRetakeGuardDigits -- the widest is rounded as it stands,
	which is right for an answer that is exactly a half.

	**The window is the function's, in units of the answer's last place**,
	about ten times the largest error measured against exact values over
	the census's near-half and random answers, in four bases at 36 and 80
	places -- see the kRetakeDoubt constants in LZMath.h. The retake
	itself relies on nothing about a guard: only that a function's answer
	is right to about a unit at the width it is asked for.

	A nested call never retakes -- the depth is counted here -- and an
	answer marked exact, flagged, or reached with an error is left alone.
*/

/*
	**How far the retake widens: the extra doubling each round, to the
	answer's own places**, since 2026-10-14 -- Ziv's strategy, with a cap.
	It took three rounds of five digits until then, and an answer closer to
	the half than fifteen digits went on its internal error alone. A tiny
	argument in an odd base lands there by construction: x = d / b^k, d
	odd, puts x^2 / 2 on a half from place 2k -- which in an odd base never
	ends -- and only the next term moves it off, x^3 / 6 at 3k for e^x,
	x^4 / 24 at 4k for cos and cosh. Of 79 such near ties each, cos and
	cosh misrounded 37 and 32. With the width at least 2k, the depth they
	need is at most the width again, so the cap is the answer's places and
	a round's step over: an answer exactly on a half -- only a rational one
	can be -- costs a few rounds more before it is rounded as it stands.
*/

const long	kRetakeGuardDigits = 5;		// the first round's extra, in digits

// The calls a retake counts its depth across -- CLZAnswerBase and
//	CLZListRetake both -- whether one is under way, and how many there
//	have been. Per thread, like all of a calculation's own state...
static thread_local long	glRetakeDepth = 0;
static thread_local bool	gbRetaking = false;
static thread_local long	glRetakes = 0;

// The tail lzi_mean_record's final rounding discarded, kept before the
//	spelling check's own arithmetic rounds over the record -- see
//	CLZListRetake...
static thread_local double	gdMeanTail = -1.0;
static thread_local long		glMeanPlaces = 0;

/*
	The retake itself, for both kinds of caller: fnCall( pResult ) takes
	the function again into pResult at the working precision it is given,
	kRetakeGuardDigits wider and then the extra doubling, until the answer
	lies more than a unit of its own last place clear of the half at lPlaces
	or ends before them; at the cap the widest is used. It is rounded at lPlaces
	in lBase. Answers false, leaving pResult to be ignored, if a call fails
	or flags its answer.
*/

/*
	**Whether rounding at lRound places is in doubt, read exactly however
	long the tail**: the digits from there to the value's end within lUnits
	of its last place of the half. lzi_rounding_in_doubt reads its guard as
	a 64-bit number and calls anything longer doubt, which a retake that
	widens by hundreds of places would meet every round. Here the tail T is
	doubled digit by digit, carry and all, and 2T - 1 read off it: with a
	carry out it is the doubled digits themselves, without one their
	complement less a unit, and either way in doubt when that comes to
	2 lUnits or less. An odd base's half, which no digits reach, needs
	nothing special.
*/

static bool lzi_tail_in_doubt( const CLZNumber* pValue, long lRound, long lUnits )
{
	long	l_base = pValue->m_lBase;
	long	l_places = pValue->m_lRadix - lRound;

	if ( l_places <= 0 )
		return false;

	// 2T, least significant place first: place m_lRadix is index 0...
	std::vector<long>	doubled( l_places );
	long				l_carry = 0;

	for ( long i = 0; i < l_places; i++ )
	{
		long	l_digit = (i < pValue->m_lDigitsLength) ? pValue->m_lpDigits[ i ] : 0;
		long	l_twice = (2 * l_digit) + l_carry;

		doubled[ i ] = l_twice % l_base;
		l_carry = l_twice / l_base;
	}

	// |2T - 1| in units of the last place: the doubled digits with a carry
	//	out, their complement plus one without, read low places first and
	//	given up once past 2 lUnits...
	int64_t	l_bound = 2 * static_cast<int64_t>( lUnits );
	int64_t	l_value = 0;
	int64_t	l_weight = 1;
	bool	b_borrow = (l_carry == 0);

	for ( long i = 0; i < l_places; i++ )
	{
		long	l_digit = (l_carry != 0) ? doubled[ i ] : ((l_base - 1) - doubled[ i ]);

		if ( l_digit != 0 )
		{
			if ( l_weight > l_bound )
				return false;

			l_value += l_digit * l_weight;

			// Past 2 lUnits already, and the complement only adds one...
			if ( l_value > l_bound )
				return false;
		}

		if ( l_weight <= l_bound )
			l_weight *= l_base;
	}

	// The complement is b^n - 1 - D, so b^n - D is one more...
	if ( b_borrow )
		l_value += 1;

	return l_value <= l_bound;
}

/*
	**And an answer still in doubt is decided exactly**, when a function
	has a way to -- see CLZAnswerBase::DecideTies. In an odd base a half
	never ends, so an exact tie stays in doubt however wide the answer is
	taken, and was rounded by its internal error at the cap; and a near miss
	can lie deeper than the cap. A tieCall answering true has rounded
	pResult at lPlaces itself.

	**It is asked twice at most, and what it may spend differs.** At the
	first round in doubt, with numbers of up to kTieEarlyDigits answer's
	widths: a decision that small costs less than even the next round, on a
	near miss that round would settle -- a root of order 128 in base ten
	decided in 425 microseconds at 36 places against 597 retaken, and in 3.1
	ms at 300 against 4.1 -- and an order of 512, just past the limit, about
	breaks even, 822 against 588 at 36 places. Run to the cap, where a tie
	or a deep near miss goes, the early decision is two to five times the
	cheaper. A larger one waits: if the cap is reached still in doubt it is
	asked again, with up to kTieMostPlaces digits, since the answer has no
	other road to being right.
*/

template <class FnCall, class TieCall>
static bool lzi_take_again_wider( FnCall fnCall, TieCall tieCall, long lBase, long lPlaces, CLZNumber* pResult )
{
	long	l_saved = glWorkingPrecision;
	long	l_step = lzi_guard_places( kRetakeGuardDigits, lBase );
	long	l_cap = lPlaces + l_step;
	bool	b_answered = true;
	bool	b_tieAsked = false;
	long	l_early = (kTieEarlyDigits * (lPlaces + 2));

	glRetakes++;
	gbRetaking = true;

	// The extra doubling each round, the last at the cap -- see
	//	kRetakeGuardDigits...
	for ( long l_extra = l_step; ; l_extra = ((2 * l_extra) < l_cap) ? (2 * l_extra) : l_cap )
	{
		glWorkingPrecision = l_saved + l_extra;
		pResult->SetBase( lBase, false );

		long	lResult = fnCall( pResult );

		if ( (lResult != kNoError) || pResult->FlagIsOn( kfInfinity ) || pResult->FlagIsOn( kfQNaN ) )
		{
			b_answered = false;
			break;
		}

		if ( (pResult->m_lRadix <= lPlaces) || !lzi_tail_in_doubt( pResult, lPlaces, 1 ) )
			break;

		// Still in doubt: perhaps exactly on the half, which no width settles,
		//	or nearer it than the cap reaches -- asked first while cheap...
		if ( !b_tieAsked )
		{
			b_tieAsked = true;

			if ( tieCall( pResult, lPlaces, l_early ) )
				break;
		}

		// ...and at the cap, at any size it may reach...
		if ( l_extra >= l_cap )
		{
			if ( l_early < kTieMostPlaces )
				tieCall( pResult, lPlaces, kTieMostPlaces );

			break;
		}
	}

	glWorkingPrecision = l_saved;
	gbRetaking = false;

	if ( !b_answered )
		return false;

	if ( pResult->m_lBase != lBase )
		pResult->SetBase( lBase, true );

	pResult->SetPrecision( lPlaces );
	pResult->Round();

	return true;
}

/*
	**Whether an answer just rounded at its width leaves that rounding in
	doubt**: its guard, the digits the rounding discarded -- see
	CLZNumber::LastDiscardedTail -- lying within dDoubt of the half, in
	units of the answer's last place. CLZAnswerBase asks it as a function
	returns; the two pairs, lz_sine_cosine and lz_hypersine_hypercosine,
	ask it of each answer, which is why it is out of the class -- since
	stage 3 of the complex layer, 2026-10-22. A flagged answer is never in
	doubt.
*/

static bool lzi_rounding_recorded_in_doubt( CLZNumber* pAnswer, bool bQuotientWidth, double dDoubt )
{
	if ( pAnswer->FlagIsOn( kfInfinity ) || pAnswer->FlagIsOn( kfQNaN ) || pAnswer->FlagIsOn( kfDiv0 ) )
		return false;

	/*
		**Only a rounding at the answer's own width is the answer's**:
		the working precision, and for a small quotient its noughts as
		well. A function whose answer took no final rounding -- e^x
		of 1e-20 is a sum at the width -- leaves an inner one as the
		last, at a guard's width, and a retake taken from that rounded
		at it: forty places at the default.
	*/

	long	l_width = glWorkingPrecision;

	if ( bQuotientWidth )
	{
		long	l_lead = lzi_leading_exponent( pAnswer );

		if ( l_lead < -1 )
			l_width += -l_lead - 1;
	}

	/*
		**A carry to a power of the base leaves the rounding a place
		wider than the answer**: 0.00999...97 rounds at the places
		0.00999... takes and carries to 0.01, which the rule gives one
		fewer -- see lzi_round_quotient_answer. Until 2026-10-15 the
		nought trimmed after it was a second rounding, recorded over the
		first with a tail of nought, and the doubt was never seen: a
		tangent 1e-12 of a unit under the half answered 0.01 where
		0.00999...9 is right...
	*/
	long	l_rounded = CLZNumber::LastRoundedPlaces();

	if ( (l_rounded != l_width) &&
		 !(bQuotientWidth && (l_rounded == (l_width + 1)) && lzi_is_power_of_base( pAnswer )) )
		return false;

	double	d_tail = CLZNumber::LastDiscardedTail();

	return (d_tail >= 0.0) && (std::fabs( d_tail - 0.5 ) <= dDoubt);
}

class CLZAnswerBase
{
public:

	// A missing left operand is refused by the function itself, as
	//	lz_agm does, and leaves no base to answer in. The function is the
	//	engine's own work from here to its end -- see CLZInternalWork...
	CLZAnswerBase( CLZNumber* leftOperand, CLZNumber* targetNumber )
		: m_pAnswer( nullptr ), m_lBase( 0 ), m_bExact( false ), m_bShort( false ), m_bQuotientWidth( false ),
		  m_pLeft( leftOperand ), m_pRight( nullptr ), m_pLeftCopy( nullptr ), m_pRightCopy( nullptr ),
		  m_fnUnary( nullptr ), m_fnBinary( nullptr ), m_dDoubt( 0.0 ),
		  m_work( (leftOperand == nullptr) ? nullptr : ((targetNumber != nullptr) ? targetNumber : leftOperand) )
	{
		glRetakeDepth++;

		if ( leftOperand != nullptr )
		{
			m_pAnswer = (targetNumber != nullptr) ? targetNumber : leftOperand;
			m_lBase = leftOperand->m_lBase;
		}
	}

	~CLZAnswerBase()
	{
		glRetakeDepth--;

		if ( (m_fnUnary != nullptr) || (m_fnBinary != nullptr) )
		{
			CLZNumber::RecordRoundings( false );

			if ( InDoubt() )
				TakeAgainWider();

			delete m_pLeftCopy;
			delete m_pRightCopy;
		}

		if ( (m_pAnswer == nullptr) || m_pAnswer->IsConstant() )
			return;

		if ( m_pAnswer->m_lBase != m_lBase )
			m_pAnswer->SetBase( m_lBase, true );

		if ( !m_bExact )
			lzi_write_to_width( m_pAnswer, m_bQuotientWidth );
		else if ( m_bShort && m_pAnswer->FlagIsOff( kfInfinity ) && m_pAnswer->FlagIsOff( kfQNaN ) )
			m_pAnswer->TrimTrailing();
	}

	void Exact( bool bExact = true )
	{
		if ( bExact )
			m_bExact = true;
	}

	void ExactAndShort( bool bExact = true )
	{
		if ( bExact )
		{
			m_bExact = true;
			m_bShort = true;
		}
	}

	void QuotientWidth()
	{
		m_bQuotientWidth = true;
	}

	// The function and its doubt window, at the outermost call only --
	//	see kRetakeGuardDigits above. The operands are copied now, before the
	//	answer can be written over the left one...
	void Retake( LZUnaryFunction fn, double dDoubt )
	{
		if ( Arm( dDoubt ) )
			m_fnUnary = fn;
	}

	/*
		**A way to round exactly, tie or near miss, for a retake to ask**
		-- see lzi_take_again_wider. Called after Retake, and only an
		armed retake keeps it: a nested call never retakes.
	*/

	void DecideTies( LZTieDecider fn )
	{
		if ( (m_fnBinary != nullptr) || (m_fnUnary != nullptr) )
			m_fnTie = fn;
	}

	void Retake( LZBinaryFunction fn, CLZNumber* rightOperand, double dDoubt )
	{
		if ( (rightOperand != nullptr) && Arm( dDoubt ) )
		{
			m_fnBinary = fn;
			m_pRight = rightOperand;

			if ( m_pAnswer == rightOperand )
			{
				m_pRightCopy = new CLZNumber;
				m_pRightCopy->SetBase( rightOperand->m_lBase, false );
				m_pRightCopy->Assign( rightOperand );
			}
		}
	}

	// The retakes this process has made, for the tests...
	static long RetakeCount()
	{
		return glRetakes;
	}

private:

	CLZAnswerBase( const CLZAnswerBase& );
	CLZAnswerBase& operator=( const CLZAnswerBase& );

	bool Arm( double dDoubt )
	{
		if ( (glRetakeDepth != 1) || gbRetaking || (m_pAnswer == nullptr) || m_pAnswer->IsConstant() )
			return false;

		m_dDoubt = dDoubt;

		// Copied only if the answer is to be written over it; otherwise it
		//	is still there for a retake, which is rare...
		if ( m_pAnswer == m_pLeft )
		{
			m_pLeftCopy = new CLZNumber;
			m_pLeftCopy->SetBase( m_pLeft->m_lBase, false );
			m_pLeftCopy->Assign( m_pLeft );
		}

		CLZNumber::ForgetLastRounding();
		CLZNumber::RecordRoundings( true );

		return true;
	}

	// The final rounding's guard within the window of the half -- see
	//	lzi_rounding_recorded_in_doubt...
	bool InDoubt() const
	{
		if ( m_bExact || (m_pAnswer == nullptr) || m_pAnswer->IsConstant() )
			return false;

		return lzi_rounding_recorded_in_doubt( m_pAnswer, m_bQuotientWidth, m_dDoubt );
	}

	void TakeAgainWider()
	{
		long		l_saved = glWorkingPrecision;
		CLZNumber	result;

		// Fresh operands each round, from the copies where there are any...
		auto	fnCall = [ this ]( CLZNumber* pResult ) -> long
		{
			CLZNumber	left;
			CLZNumber	right;
			CLZNumber*	p_left = (m_pLeftCopy != nullptr) ? m_pLeftCopy : m_pLeft;

			left.SetBase( p_left->m_lBase, false );
			left.Assign( p_left );

			if ( m_fnUnary != nullptr )
				return m_fnUnary( &left, pResult );

			CLZNumber*	p_right = (m_pRightCopy != nullptr) ? m_pRightCopy : m_pRight;

			right.SetBase( p_right->m_lBase, false );
			right.Assign( p_right );

			return m_fnBinary( &left, &right, pResult );
		};

		// The tie decider on the same fresh operands...
		auto	tieCall = [ this ]( CLZNumber* pResult, long lPlaces, long lMostDigits ) -> bool
		{
			if ( m_fnTie == nullptr )
				return false;

			CLZNumber	left;
			CLZNumber	right;
			CLZNumber*	p_left = (m_pLeftCopy != nullptr) ? m_pLeftCopy : m_pLeft;

			left.SetBase( p_left->m_lBase, false );
			left.Assign( p_left );

			// A function of one operand has no right one to pass...
			if ( m_fnBinary == nullptr )
				return m_fnTie( pResult, &left, nullptr, lPlaces, lMostDigits );

			CLZNumber*	p_right = (m_pRightCopy != nullptr) ? m_pRightCopy : m_pRight;

			right.SetBase( p_right->m_lBase, false );
			right.Assign( p_right );

			return m_fnTie( pResult, &left, &right, lPlaces, lMostDigits );
		};

		if ( lzi_take_again_wider( fnCall, tieCall, m_lBase, CLZNumber::LastRoundedPlaces(), &result ) )
		{
			m_pAnswer->Assign( &result );
			m_pAnswer->SetPrecision( l_saved );

			// ...and one decided up to a power of the base is a place too
			//	wide, as the first rounding's was -- 21 of 72 such were
			//	written so until 2026-10-15...
			if ( m_bQuotientWidth )
				lzi_trim_after_carry( m_pAnswer );
		}
	}

	CLZNumber*		m_pAnswer;
	long			m_lBase;
	bool			m_bExact;
	bool			m_bShort;
	bool			m_bQuotientWidth;
	CLZNumber*		m_pLeft;
	CLZNumber*		m_pRight;
	CLZNumber*		m_pLeftCopy;
	CLZNumber*		m_pRightCopy;
	LZUnaryFunction		m_fnUnary;
	LZBinaryFunction	m_fnBinary;
	LZTieDecider		m_fnTie = nullptr;
	double			m_dDoubt;
	CLZInternalWork	m_work;

};


/*
	**The retake for a function of a list**, lz_gmean and lz_hmean, which
	answer through lzi_mean_record rather than a CLZAnswerBase: the same
	window test and the same lzi_take_again_wider, on the tail
	lzi_mean_record keeps. Constructed first in the function, it counts
	the depth as CLZAnswerBase does, so the roots and powers a mean calls
	never retake on their own. The values are copied only if the target is
	one of them.
*/

typedef long (*LZListFunction)( CLZNumber**, long, CLZNumber* );

class CLZListRetake
{
public:

	CLZListRetake( LZListFunction fn, CLZNumber** values, long count, CLZNumber* targetNumber, double dDoubt )
		: m_fn( fn ), m_pValues( values ), m_lCount( count ), m_pTarget( targetNumber ), m_pCopies( nullptr ),
		  m_dDoubt( dDoubt ), m_bArmed( false )
	{
		glRetakeDepth++;

		if ( (glRetakeDepth != 1) || gbRetaking || (targetNumber == nullptr) || targetNumber->IsConstant() ||
			 (values == nullptr) || (count < 1) )
			return;

		bool	b_overwritten = false;

		for ( long i = 0; i < count; i++ )
		{
			if ( values[ i ] == nullptr )
				return;

			if ( values[ i ] == targetNumber )
				b_overwritten = true;
		}

		if ( b_overwritten )
		{
			m_pCopies = new CLZNumber[ count ];

			for ( long i = 0; i < count; i++ )
			{
				m_pCopies[ i ].SetBase( values[ i ]->m_lBase, false );
				m_pCopies[ i ].Assign( values[ i ] );
			}
		}

		gdMeanTail = -1.0;
		glMeanPlaces = 0;
		CLZNumber::ForgetLastRounding();
		CLZNumber::RecordRoundings( true );
		m_bArmed = true;
	}

	~CLZListRetake()
	{
		glRetakeDepth--;

		if ( !m_bArmed )
			return;

		CLZNumber::RecordRoundings( false );

		if ( m_pTarget->FlagIsOn( kfInfinity ) || m_pTarget->FlagIsOn( kfQNaN ) || m_pTarget->FlagIsOn( kfDiv0 ) ||
			 (glMeanPlaces != glWorkingPrecision) || (gdMeanTail < 0.0) || (std::fabs( gdMeanTail - 0.5 ) > m_dDoubt) )
		{
			delete[] m_pCopies;
			return;
		}

		long		l_saved = glWorkingPrecision;
		long		l_base = m_pTarget->m_lBase;
		CLZNumber	result;

		auto	fnCall = [ this ]( CLZNumber* pResult ) -> long
		{
			CLZNumber*	p_values = new CLZNumber[ m_lCount ];
			CLZNumber**	pp_values = new CLZNumber*[ m_lCount ];

			for ( long i = 0; i < m_lCount; i++ )
			{
				CLZNumber*	p_from = (m_pCopies != nullptr) ? &m_pCopies[ i ] : m_pValues[ i ];

				p_values[ i ].SetBase( p_from->m_lBase, false );
				p_values[ i ].Assign( p_from );
				pp_values[ i ] = &p_values[ i ];
			}

			long	lResult = m_fn( pp_values, m_lCount, pResult );

			delete[] pp_values;
			delete[] p_values;

			return lResult;
		};

		if ( lzi_take_again_wider( fnCall, []( CLZNumber*, long, long ) { return false; }, l_base, l_saved, &result ) )
		{
			m_pTarget->Assign( &result );
			m_pTarget->SetPrecision( l_saved );
		}

		delete[] m_pCopies;
	}

private:

	CLZListRetake( const CLZListRetake& );
	CLZListRetake& operator=( const CLZListRetake& );

	LZListFunction	m_fn;
	CLZNumber**		m_pValues;
	long			m_lCount;
	CLZNumber*		m_pTarget;
	CLZNumber*		m_pCopies;
	double			m_dDoubt;
	bool			m_bArmed;
};


static long lzi_first_failure( long lSoFar, long lNext )
{
	return (lSoFar != kNoError) ? lSoFar : lNext;
}


/*
	**An answer's own size comes on top of the places asked for**, and a
	function that forms it at relative precision does not deliver them
	unless it works wider. 36 places of 3.2e12 are 49 significant digits;
	lz_power( 10, 12.5 ) gave 37, 2.8e11 units short of its 36th place,
	because it multiplies an exact whole power by a fractional one that
	was right to 36 places -- relative -- and the product is no better.
	The same shape in lz_exp_b and lz_exp_n, in lz_root_n's exponential
	of a logarithm, lz_log_n's quotient of two, and lz_raise_n's
	squarings. None was in the sweep or the audit until 2026-09-24, which
	is how six functions came to lack what lz_raise_e had long had.

	So: learn the answer's magnitude, work at the working precision plus
	that magnitude plus a guard, and round back. Every function
	using this keeps its algorithm as it was, as a worker; only the width
	it runs at changes, which is why its behaviour at the edges --
	refusals, flags, codes -- is exactly what it was.

	**The guard is four places, and it is needed.** Without it 1.0001 to
	the 100000th -- seventeen squarings -- came out 550 units out in its
	36th place, 1.000001 to the 5000000th 17290, and the audit lost 17
	correctly rounded answers. A further guard of as many places as a
	power's exponent has whole digits was written in first, on the
	reasoning that a squaring chain's error grows with the exponent, and
	measured it bought nothing: every one of those cases and every audit
	count came out the same without it. So it is not here. **The
	reasoning was right and the measure was not**: correctly rounded
	random answers say only that the error is under half a unit, and the
	retake needs it under its window. A thousandth power's error is 1e-2
	of a unit, so since 2026-10-10 lzi_raise_n works as many places wider
	as its exponent has -- there, so the squaring chain carries it for
	every caller.

	lz_raise_e calls the unwidened lzi_raise_n, since it widens by its
	own answer already and would otherwise do it twice.

	**The magnitude is estimated, not probed**, and this is the one place
	LZMath touches floating point. lz_raise_e learns its answer's size by
	running the work at a dozen places, and that was tried here first --
	but a logarithm or an exponential has fixed costs even at twelve
	places, and the probe came to half of lz_root_n's old time. Every one
	of these answers has a closed-form size instead: 10^x has x digits, a
	power y log|x|, a root log|x| / n, a logarithm log|ln x / ln n|. A
	logarithm read off an operand's leading digits and where they sit
	sizes the width to within a place, and the guard covers that place --
	once it holds as many digits as a double does and takes the
	difference from one near one, which it did not until 2026-10-12: see
	lzi_ln_estimate. The doubles decide how wide to work and
	never reach an answer. Where no estimate can be formed -- a zero, a
	flagged value, a logarithm to base one -- the probe is still there,
	and the edges behave exactly as they did.

	**The guard is a count of digits**, like every other -- see
	lzi_guard_places. It was added as places until 2026-09-26, the one
	guard the census of that day missed, being written as a sum rather
	than an increment: four places were four bits in base two, and
	log_n missed there where it missed nowhere else.

	**The worker is handed copies of its operands, made at the width of
	the pass.** A CLZNumber takes its precision when it is made and
	Assign carries it along -- the trap lzi_log_gamma fell into -- so a
	worker that copies its operand, or splits it with lz_integer and
	lz_fraction, gets working values pinned at the caller's 36 however
	wide glWorkingPrecision has been set. The first version of this
	widened the precision and not the operands, and lz_power, lz_exp_n
	and lz_root_n came out at 37 to 38 significant digits as if nothing
	had happened.
*/

const long kAnswerGuardPlaces = 4;

typedef long (*LZWidenedWork)( CLZNumber*, CLZNumber*, CLZNumber* );

static long lzi_raise_n( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber );

// Copies of the operands at the width of the pass, so that nothing a
//	worker derives from them is pinned at the caller's precision...
static long lzi_work_at_width( LZWidenedWork fnWork, CLZNumber* pLeft, CLZNumber* pRight,
							   CLZNumber* pResult )
{
	CLZNumber	left;
	CLZNumber	right;

	left.Assign( pLeft );
	left.SetPrecision( glWorkingPrecision );

	if ( pRight != nullptr )
	{
		right.Assign( pRight );
		right.SetPrecision( glWorkingPrecision );
	}

	return fnWork( &left, (pRight != nullptr) ? &right : nullptr, pResult );
}

/*
	The natural logarithm of |value|, from as many leading digits as a
	double holds and the place the first of them sits in, for sizing a
	width or an answer and nothing else. False where there is no
	logarithm to estimate.

	**Good to about sixteen figures of ln|value| itself, near one as
	well.** A power's size is its exponent times this, so an error here
	is multiplied by the exponent, and until 2026-10-12 it had two:

	- **Sixteen digits, in every base.** In base two that is sixteen
	  bits, an error of 3E-5 in the logarithm: times an exponent of a
	  million, 43 places of the answer's size, against a margin of two
	  for calling it nought.
	- **Near one, the logarithm of a truncated value is mostly
	  cancellation.** 1 + 1E-20 read as 1.000000000000000, whose
	  logarithm and the place's cancel to a rounding residue of 1E-15;
	  to the 1E19th, that called e^0.1 nought, and answered 0. Within
	  kLnEstimateNearOne of one the logarithm is log1p of the difference
	  from one instead, which a subtraction gives exactly, past the
	  width as well.
*/

const double	kLnEstimateNearOne = 1.0E-3;

static bool lzi_ln_estimate( CLZNumber* pValue, double* pdLn )
{
	if ( pValue->FlagIsOn( kfQNaN ) || pValue->FlagIsOn( kfInfinity ) )
		return false;

	long	l_top = pValue->m_lDigitsLength - 1;

	while ( (l_top >= 0) && (pValue->m_lpDigits[ l_top ] == 0) )
		l_top--;

	if ( l_top < 0 )
		return false;

	// As many digits as a double's 53 bits hold exactly: 53 in base two,
	//	15 in ten, 13 in sixteen, counted once...
	static const std::array<long, 17>	sl_digits = []
	{
		std::array<long, 17>	digits = {};

		for ( long l_base = 2; l_base <= 16; l_base++ )
			digits[ l_base ] = static_cast<long>( 53.0 * std::log( 2.0 ) / std::log( static_cast<double>( l_base ) ) );

		return digits;
	}();

	long	l_wanted = sl_digits[ pValue->m_lBase ];
	double	d_mantissa = 0.0;
	long	l_index = l_top;

	for ( long l_taken = 0; (l_taken < l_wanted) && (l_index >= 0); l_taken++, l_index-- )
		d_mantissa = (d_mantissa * pValue->m_lBase) + pValue->m_lpDigits[ l_index ];

	// The last digit taken sits at place l_index + 1, counted from the radix...
	*pdLn = std::log( d_mantissa ) +
			(static_cast<double>( (l_index + 1) - pValue->m_lRadix ) * std::log( static_cast<double>( pValue->m_lBase ) ));

	if ( std::fabs( *pdLn ) < kLnEstimateNearOne )
	{
		// One itself, which every exponent of one asks about, without
		//	building anything: its only digit is the units', and a one...
		bool	b_one = (l_top == pValue->m_lRadix) && (pValue->m_lpDigits[ l_top ] == 1);

		for ( long i = l_top - 1; b_one && (i >= 0); i-- )
			b_one = (pValue->m_lpDigits[ i ] == 0);

		if ( b_one )
		{
			*pdLn = 0.0;

			return true;
		}

		CLZNumber	magnitude;
		CLZNumber	one;
		CLZNumber	difference;

		// A subtraction keeps every place, whatever the width...
		lz_abs( pValue, &magnitude );

		one.SetBase( magnitude.m_lBase, false );
		one.Assign( 1L );
		lz_subtract( &magnitude, &one, &difference );

		double	d_lnDifference = 0.0;

		difference.CheckForZero();

		if ( difference.IsZero() )
			*pdLn = 0.0;
		else if ( lzi_ln_estimate( &difference, &d_lnDifference ) )
			*pdLn = std::log1p( (difference.IsNegative() ? -1.0 : 1.0) * std::exp( d_lnDifference ) );
	}

	return true;
}

// The value itself, for an exponent or an order: its sign on the exponential
//	of the estimate above, and nought for nought...
static bool lzi_value_estimate( CLZNumber* pValue, double* pdValue )
{
	double	d_ln = 0.0;

	pValue->CheckForZero();

	if ( pValue->IsZero() && pValue->FlagIsOff( kfQNaN ) && pValue->FlagIsOff( kfInfinity ) )
	{
		*pdValue = 0.0;
		return true;
	}

	if ( !lzi_ln_estimate( pValue, &d_ln ) || (d_ln > 700.0) )
		return false;

	*pdValue = (pValue->IsNegative() ? -1.0 : 1.0) * std::exp( d_ln );

	return true;
}

// Whole places of an answer whose natural logarithm is about dLnAnswer, in
//	base lBase, plus one for the estimate's own uncertainty; -1 when the
//	estimate is not a number, which sends the caller back to the probe...
static long lzi_places_for( double dLnAnswer, long lBase )
{
	if ( !std::isfinite( dLnAnswer ) )
		return -1;

	double	d_places = dLnAnswer / std::log( static_cast<double>( lBase ) );

	return (d_places < 0.0) ? 0 : static_cast<long>( d_places ) + 2;
}

/*
	**An answer too large to hold is refused, and one too small to show
	is nought, both before any of it is worked out**. Until 2026-10-11
	neither was: e^(1E24)
	multiplied ever wider numbers until memory ran out, and so did
	e^(-1E24), whose answer at any width is nought.

	The limit is on the answer's order of magnitude in places of its
	base: the whole places of a large answer, or, for a quotient, which
	keeps its significant digits however small, the noughts after its
	point. It is the process's, like the settings the tests change, so
	it wants setting before a second thread starts.

	**The estimate is good to about sixteen figures of the logarithm**,
	so an answer within a place of the limit may fall either side of it.
	This said a dozen, of the estimate as it was before 2026-10-12, which
	was not true near one or in base two -- see lzi_ln_estimate: near one
	it called e^0.1 nought. And nought needs a
	margin below the last place, not a guess at it: an answer two places
	under the width cannot round to anything but nought.
*/

static long	gAnswerLimitPlaces = kAnswerLimitPlaces;

const long	kNoughtMarginPlaces = 2;

long lzu_set_answer_limit( long lPlaces )
{
	long	l_previous = gAnswerLimitPlaces;

	gAnswerLimitPlaces = (lPlaces > 0) ? lPlaces : kAnswerLimitPlaces;

	return l_previous;
}

// The value as a double, or plus or minus infinity past what a double
//	holds: an exponent's size only has to be compared, and an exponent of
//	1E400 is past the limit whatever it raises...
static bool lzi_value_or_infinity( CLZNumber* pValue, double* pdValue )
{
	if ( lzi_value_estimate( pValue, pdValue ) )
		return true;

	double	d_ln = 0.0;

	if ( !lzi_ln_estimate( pValue, &d_ln ) )
		return false;

	*pdValue = pValue->IsNegative() ? -HUGE_VAL : HUGE_VAL;

	return true;
}

enum ELZAnswerSize
{
	AS_fits,
	AS_too_large,
	AS_nought
};

/*
	What an answer whose natural logarithm is about dLnAnswer comes to in
	lBase: one past the limit, one that rounds to nought at the working
	precision, or one to work out. No estimate is no refusal.

	**A quotient's small answer needs nothing here.** sech and csch keep
	their significant digits however small, and are too small to hold
	exactly when cosh and sinh are too large -- whose exponentials refuse,
	and one over an infinity is an unflagged nought, which with kTooLarge
	is the refusal a quotient too small wants. They had a check of their
	own, which mutation showed could not change an answer.
*/

static ELZAnswerSize lzi_answer_size( double dLnAnswer, long lBase )
{
	if ( std::isnan( dLnAnswer ) )
		return AS_fits;

	double	d_places = dLnAnswer / std::log( static_cast<double>( lBase ) );

	if ( d_places > static_cast<double>( gAnswerLimitPlaces ) )
		return AS_too_large;

	if ( d_places < -static_cast<double>( glWorkingPrecision + kNoughtMarginPlaces ) )
		return AS_nought;

	return AS_fits;
}

// Refuse the answer or write nought, as eSize says, answering whether it
//	did, with the function's return in *plResult. A refusal zeroes the
//	result and flags kfInfinity -- see lzi_domain_failure. Nought is
//	written at the working precision, so it is spelled exactly as worked
//	out it would have been...
static bool lzi_answer_beyond( ELZAnswerSize eSize, CLZNumber* leftOperand, CLZNumber* targetNumber,
							   long* plResult )
{
	if ( eSize == AS_fits )
		return false;

	if ( eSize == AS_too_large )
		*plResult = lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kTooLarge );
	else
	{
		CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		*plResult = kNoError;
	}

	return true;
}

// ln Gamma(z) for z above nought, from Stirling's leading terms once z is
//	shifted past ten -- a few good figures, for sizing an answer and
//	nothing else. std::lgamma would do, but it writes signgam, which every
//	thread shares...
static double lzi_ln_gamma_estimate( double dZ )
{
	if ( !std::isfinite( dZ ) )
		return HUGE_VAL;

	double	d_shift = 0.0;

	while ( dZ < 10.0 )
	{
		d_shift -= std::log( dZ );
		dZ += 1.0;
	}

	return d_shift + ((dZ - 0.5) * std::log( dZ )) - dZ + 0.918938533204673;
}

// tanh and coth are plus or minus one past where 2 e^-2|x| falls two places
//	under the width -- tanh(x) is 1 - 2/(e^2x + 1) -- and answering them
//	from sinh and cosh would work out e^|x| first, which may be too large to
//	hold. Answers whether it wrote the answer...
static bool lzi_unit_beyond_width( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	double	d_x = 0.0;

	if ( !lzi_value_or_infinity( leftOperand, &d_x ) )
		return false;

	double	d_reach = (static_cast<double>( glWorkingPrecision + kNoughtMarginPlaces ) *
					   std::log( static_cast<double>( leftOperand->m_lBase ) )) + std::log( 4.0 );

	if ( 2.0 * std::fabs( d_x ) <= d_reach )
		return false;

	CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_target->IsConstant()) )
		p_target->Assign( (d_x < 0.0) ? -1L : 1L );

	return true;
}

static long lzi_widened_by_answer( LZWidenedWork fnWork, CLZNumber* pLeft, CLZNumber* pRight,
								   CLZNumber* pTarget, long lMagnitude )
{
	long	l_saved = glWorkingPrecision;
	long	l_magnitude = lMagnitude;

	// No estimate: the probe, at a dozen places, only to be measured...
	if ( l_magnitude < 0 )
	{
		glWorkingPrecision = kRaiseEProbePlaces;

		CLZNumber	probe;

		lzi_work_at_width( fnWork, pLeft, pRight, &probe );
		l_magnitude = probe.GetOrderOfMagnitude();
	}

	// ...and the real pass, as wide as that asks...
	glWorkingPrecision = l_saved + lzi_guard_places( kAnswerGuardPlaces, pLeft->m_lBase ) + l_magnitude;

	CLZNumber	result;
	long		lResult = lzi_work_at_width( fnWork, pLeft, pRight, &result );

	glWorkingPrecision = l_saved;

	// The guard places never leave...
	result.SetPrecision( l_saved );
	result.Round();

	CLZNumber*	p_answer = (pTarget != nullptr) ? pTarget : pLeft;

	if ( !(p_answer->IsConstant()) )
		p_answer->Assign( &result );

	return lResult;
}


long lzi_add_op( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber, long resultSign )
{
	if ( leftOperand->m_lBase != rightOperand->m_lBase )
		return kBaseMismatch;

	long	l_lengthRight = 0;
	long	l_adjustRight = 0;
	long*	lp_rightDigits = nullptr;

	long	l_lengthResult = 0;
	long*	lp_resultDigits = nullptr;

	long	l_carryDigit = 0;

	// Set up result...
	l_lengthResult = targetNumber->m_lDigitsLength;

	// Get right radix adjust and digits length...
	l_adjustRight = targetNumber->m_lRadix - rightOperand->m_lRadix;
	l_lengthRight = rightOperand->m_lDigitsLength + l_adjustRight;

	// Get copies of the operand pointers...
	lp_resultDigits = targetNumber->m_lpDigits;
	lp_rightDigits = rightOperand->m_lpDigits;

	// Add in the right side...
	for ( long i=0; i<l_lengthResult; i++ )
	{
		if ( l_carryDigit )
			*lp_resultDigits += l_carryDigit;

		if ( i < l_lengthRight )
		{
			if ( l_adjustRight )
				l_adjustRight--;
			else
				*lp_resultDigits += *(lp_rightDigits++);
		}

		// Carry over as required...
		if ( *lp_resultDigits > leftOperand->m_lMaxDigit )
		{
			*lp_resultDigits -= leftOperand->m_lBase;
			l_carryDigit = 1;
		}
		else
			l_carryDigit = 0;

		// Next result digit...
		lp_resultDigits++;
	}

	// Add in the remaining carry digit, if any...
	if ( l_carryDigit )
	{
		targetNumber->AlignStorage( targetNumber->m_lDigitsLength + 1, targetNumber->m_lRadix );
		targetNumber->m_lpDigits[ targetNumber->m_lDigitsLength - 1 ] += l_carryDigit;
	}

	// Set the target sign...
	targetNumber->m_lSign = resultSign;

	return kNoError;
}



/*
	Where division stops dividing and starts multiplying by a
	reciprocal. When this was set, long division cost the square of
	the quotient's length and a reciprocal a handful of
	multiplications, and the crossover measured at a few hundred
	places -- eight times in the reciprocal's favour at three hundred
	places, twenty-two at ten thousand.

	**That stopped being true, and the floor is now out of reach.**
	Long division was reworked on 2026-09-20 and beats the reciprocal
	at every size measured, up to 16,000 digits over a 16,000-digit
	divisor; on 2026-09-21 it gained a sticky digit and rounds
	correctly, which had been the reciprocal's one remaining claim.
	The reciprocal also scales quadratically where it should not, so
	nothing measured suggests it overtakes further out. A million
	places is past anything this engine has been asked for, which
	retires the path in practice while keeping it -- fixed, and tested
	through lzi_set_reciprocal_floor -- for the day it is made faster.

	**That day was 2026-09-26.** The quadratic was RoundAt discarding
	a digit at a time, and once it was not the reciprocal beat long
	division for a long divisor. It was made exact then too, and now
	takes long divisors by their shape -- see
	lzi_reciprocal_by_shape. This floor stays for the tests.

	The seed is what the engine's own division supplies as a starting
	estimate, and the guard is how far past the answer the iteration
	works so the last place is not itself the rounding.
*/

const long kReciprocalFloor = 1000000;
const long kReciprocalSeedPlaces = 12;
const long kReciprocalGuardPlaces = 24;

// The places past a pass's width that the divisor keeps for it, and the
//	widest divisor the last reciprocal's passes multiplied by, which the
//	tests read -- see lzi_reciprocal_newton...
const long kReciprocalDivisorGuard = 4;

static thread_local long gReciprocalWidest = 0;

long lzi_reciprocal_widest_divisor()
{
	return gReciprocalWidest;
}

/*
	The floor in force, which only the tests change. With the floor
	where it is, nothing an affordable test divides would ever reach
	the reciprocal, and a path no test reaches rots: the base fault
	below lived in it for as long as it existed. The suite lowers the
	floor to exercise the path at ordinary sizes in every base, and
	puts it back.
*/

static long gReciprocalFloor = kReciprocalFloor;

long lzi_set_reciprocal_floor( long lPlaces )
{
	long	l_previous = gReciprocalFloor;

	gReciprocalFloor = lPlaces;

	return l_previous;
}

/*
	**Where the reciprocal beats long division: a long divisor, from a
	width that depends on the base.** Measured 2026-09-26, once the
	rounding stopped making it quadratic and the reciprocal was made
	exact. Long division costs the quotient's length times the
	divisor's, and the reciprocal a few multiplies as long as the
	quotient whatever the divisor, so a short divisor goes to long
	division at every width and only a long one gains.

	The reciprocal first won, with the divisor at three quarters of the
	quotient, in places of the base: twelve from 4,000, sixteen from
	6,000, ten from 8,000 and eight from 16,000. But at those widths it
	won by 7 to 18%, and exp(1) at 9,000 places, whose one long division
	sits right there, came out 6% slower through it. **So each floor is
	half as far again past where it first won**, where every shape
	measured gains 15% or more -- in base ten 25% at 16,000 and 33% at
	32,000 for a short dividend, 25% at 32,000 for a long one.

	Base two packs 31 digits a limb and its long division stays ahead
	to 131,072 places, where the two tie, so it never takes this. A
	base with no limb table multiplied a digit at a time, and the
	reciprocal lost there by three times at every width measured.

	**Measured again 2026-10-02, after Karatsuba had made the
	reciprocal's products faster**, by the same rule and with both
	dividend shapes, 1/d and a dividend as long as the quotient, over a
	divisor three quarters and fully as long. Every shape now wins
	from 6,000 places in base ten, 4,000 in hex and in base twelve, and
	8,000 in base eight, the long dividend over the shorter divisor
	always last to; base two still never does, to 32,000. Half as far
	again puts the floors at 9,000, 6,000, 6,000 and 12,000, where every
	shape gains 30% or more. Base twelve's did not move. e at 9,000
	places, the 6% loss that set the rule, is 16% faster through it now,
	and every workload whose division fell between the old floors and
	these gained 16 to 30%.

	**Bases three to fifteen, by the same rule, 2026-10-12**, once they
	had limbs and their long division and multiplication both packed.
	Every shape first wins from 24,000 places in
	base three, 16,000 in four, 12,000 in five and six, 8,000 in seven,
	nine and eleven, 6,000 in fourteen and fifteen and 4,000 in thirteen;
	the floors are half as far again. A smaller base wants more places,
	each being fewer bits. Base two still never wins, to 32,000.

	**A staircase since 2026-10-16** -- performance item 7 of the review
	of 2026-10-15, on Rich's choice of a measurement in two dimensions.
	The rule above counted the
	quotient's places, noughts before its leading digit and all: 1 over a
	16,000-digit whole number at 16,000 places looked 32,001 places wide
	and went the long way, 14.0 ms against 5.0, while 1 over the same
	digits as a fraction, the same work exactly, took the reciprocal. And
	it was measured only with the divisor at three quarters of the
	quotient and more: once each pass took the divisor to its own width,
	the reciprocal won over a far wider region -- a short quotient over a
	long divisor by as much as six times.

	**So it goes by the quotient's significant digits and the divisor's.**
	Per base, the divisor at least rho of the quotient -- past that the
	gain fades -- and, for some step, the quotient at least the step's
	digits and the divisor at least its floor; the floors never rise as
	the quotient grows, a short quotient wanting a long divisor. Measured
	over a grid of quotient and divisor digits in every base, twice, with
	1 over the divisor and with a dividend as long as the quotient, and
	chosen by dynamic programming over the grid's quotient widths: every
	point the old rule took and won by 2% stays, on both families it read
	differently -- a quotient near one, and 1 over a whole number, its
	noughts counted -- and every point added gains 15% or more in both
	shapes and both passes; of the staircases that do, these leave the
	least time. **On the grid no point is slower than the old rule on
	either family, in any base**; 28 to 67 points a base are faster, by as
	much as 6.4 times, and the grid's time falls by 26 to 53%.
	Base two takes the reciprocal now, for a long divisor.

	Simpler shapes were tried first and left points slower than the old
	rule: one floor for each, which lost a divisor just short of its floor
	by as much as 1.75 times; two floors beside the old rule, which lost
	the old rule's 1 over a whole number; a product of the two. The
	staircase is what the data asked for.
*/

static bool lzi_reciprocal_by_shape( long lBase, int64_t lQuotientDigits, int64_t lDivisorDigits )
{
	// Base by base: rho in eighths, then up to four steps, each the
	//	quotient's significant digits and the divisor's floor...
	static const long	kSteps[ 17 ][ 9 ] =
	{
		{ 0,	0, 0,	0, 0,	0, 0,	0, 0 },
		{ 0,	0, 0,	0, 0,	0, 0,	0, 0 },
		{ 4,	6644, 119589,	9966, 79726,	119589, 1661,	0, 0 },
		{ 2,	3144, 33534,	6288, 25151,	75453, 1048,	0, 0 },
		{ 3,	830, 59795,	1661, 19932,	59795, 830,	0, 0 },
		{ 2,	1431, 17168,	2146, 12876,	51504, 715,	0, 0 },
		{ 2,	643, 20562,	1285, 11566,	46263, 643,	0, 0 },
		{ 3,	592, 56798,	1183, 10650,	28399, 592,	0, 0 },
		{ 3,	554, 13288,	1107, 9966,	26575, 554,	0, 0 },
		{ 3,	524, 37726,	1048, 9432,	25151, 524,	0, 0 },
		{ 3,	250, 48000,	500, 9000,	1000, 6000,	12000, 500 },
		{ 2,	480, 15364,	960, 8642,	34569, 480,	0, 0 },
		{ 3,	232, 22239,	463, 5560,	11120, 463,	0, 0 },
		{ 2,	449, 8079,	898, 5386,	21545, 449,	0, 0 },
		{ 2,	218, 13960,	436, 7853,	31410, 436,	0, 0 },
		{ 2,	425, 10203,	850, 7652,	30610, 425,	0, 0 },
		{ 3,	208, 7474,	415, 4983,	7474, 415,	0, 0 },
	};

	if ( (lBase < 2) || (lBase > 16) || (kSteps[ lBase ][ 0 ] == 0) )
		return false;

	const long*	p_steps = kSteps[ lBase ];

	if ( (8 * lDivisorDigits) < (p_steps[ 0 ] * lQuotientDigits) )
		return false;

	for ( int i = 1; (i < 9) && (p_steps[ i ] > 0); i += 2 )
	{
		if ( (lQuotientDigits >= p_steps[ i ]) && (lDivisorDigits >= p_steps[ i + 1 ]) )
			return true;
	}

	return false;
}

// Whether the shape above chooses the path, which only the tests turn
//	off -- to compare the reciprocal with long division on the very
//	shapes it now takes...
static bool gReciprocalByShape = true;

bool lzi_set_reciprocal_by_shape( bool bByShape )
{
	bool	b_previous = gReciprocalByShape;

	gReciprocalByShape = bByShape;

	return b_previous;
}

// Whether a wide dividend is cut, which only the tests turn off -- see
//	the cut in lzi_divide_work -- and the quotient digits worked, which
//	the tests count it by...
static bool gDivideCuts = true;
static thread_local long gQuotientDigitsWorked = 0;

long lzi_quotient_digits_worked()
{
	return gQuotientDigitsWorked;
}

bool lzi_set_divide_cuts( bool bCuts )
{
	bool	b_previous = gDivideCuts;

	gDivideCuts = bCuts;

	return b_previous;
}

/*
	Raised while the reciprocal is finding its starting estimate, which
	it does with an ordinary division. Without it that division would
	take the reciprocal path in turn and call itself for ever. Setting
	the floor above the seed width hides the problem rather than
	fixing it, and a later hand lowering the floor would find out the
	hard way -- which is how this was found, by lowering it to one.
	Per thread: another thread's reciprocal is no reason to divide the
	long way.
*/

static thread_local bool gSeekingReciprocal = false;

/*
	Both live further down, beside the roots that also use them.
*/

static long lzi_leading_exponent( CLZNumber* value );
static void lzi_round_quotient_answer( CLZNumber* pAnswer );
static long lzi_small_divisor_widening( CLZNumber* pValue );
static void lzi_sine_and_cosine( CLZNumber* pArgument, CLZNumber* pSine, CLZNumber* pCosine );
static long lzi_trig_fold( CLZNumber* pArgument, CLZNumber* pReduced,
						   bool* pbNegateSine, bool* pbNegateCosine, bool* pbSwapped );
long lzi_sincos_split_floor( long lBase );
static long lzi_widen_until_seen( long (*fnTrig)( CLZNumber*, CLZNumber* ),
								  CLZNumber* pArgument, CLZNumber* pValue );
static void lzi_scale_by_base_power( CLZNumber* value, long power );
static CLZNumber* lzi_log_base_in_base( long lBase, long lPlaces, bool bByAgm = false );


/////////////////////////////////////////////////////////////////
//	Two answers from one calculation -- stage 3 of the complex layer,
//	2026-10-22: lz_sine_cosine and lz_hypersine_hypercosine

/*
	**A pair's answer rounded at its width, and whether that rounding can
	stand**: rounded as its own function rounds it -- at the working
	precision, or by the rule for a quotient -- and clear of the window
	that function retakes within, read from the rounding recorded; see
	lzi_rounding_recorded_in_doubt. One that cannot stand, or is flagged,
	goes to its own function instead, which takes it again wider and
	decides it as it always does.

	**Why one clear of the window can stand.** Worked under a guard whose
	error lies well inside the window, an answer clear of it is on the
	same side of the half as the exact value, so it rounds as a correctly
	rounded function rounds it -- and the function's own answer is
	correctly rounded, taken again or not. So the pair answers as the two
	functions would, and only an answer within the window needs its
	function: exactly the ones the function would have taken again.
*/
// The answers the pairs have handed to their own functions, this thread's,
//	for the tests -- a pair that handed every one on would answer the same,
//	slower...
static thread_local long	glPairHandoffs = 0;

long lzi_pair_handoffs()
{
	return glPairHandoffs;
}

static bool lzi_rounded_clear_of_doubt( CLZNumber* pValue, double dDoubt, bool bQuotientWidth )
{
	CLZNumber::ForgetLastRounding();
	CLZNumber::RecordRoundings( true );

	if ( bQuotientWidth )
		lzi_round_quotient_answer( pValue );
	else
	{
		pValue->SetPrecision( glWorkingPrecision );
		pValue->Round();
	}

	CLZNumber::RecordRoundings( false );

	if ( pValue->FlagIsOn( kfInfinity ) || pValue->FlagIsOn( kfQNaN ) )
		return false;

	return !lzi_rounding_recorded_in_doubt( pValue, bQuotientWidth, dDoubt );
}

/*
	**An answer worked out here, written as a function writes one** -- in
	the operand's base, to its width unless exact, short where exact and
	asked to be, with the target's rounding settings left as they were.
	CLZAnswerBase does all of it, armed for nothing.
*/
static long lzi_answer_as_function( CLZNumber* pValue, CLZNumber* pOperand, CLZNumber* pTarget,
									bool bQuotientWidth, bool bExact = false, bool bShort = false )
{
	CLZAnswerBase	answer( pOperand, pTarget );

	if ( bQuotientWidth )
		answer.QuotientWidth();

	if ( bShort )
		answer.ExactAndShort( bExact );
	else
		answer.Exact( bExact );

	if ( !(pTarget->IsConstant()) )
		pTarget->Assign( pValue );

	return kNoError;
}

// A pair refused still answers, each target it was given a not-a-number --
//	see lzi_domain_failure...
static long lzi_pair_refused( CLZNumber* pFirst, CLZNumber* pSecond )
{
	for ( CLZNumber* p_target : { pFirst, pSecond } )
	{
		if ( (p_target != nullptr) && !(p_target->IsConstant()) )
		{
			p_target->Assign( 0L );
			p_target->SetFlag( kfQNaN );
		}
	}

	return kFailed;
}


/////////////////////////////////////////////////////////////////
//	The angle unit -- stage 2b of the complex layer, 2026-10-22

/*
	**The unit an angle is read or written in, one per thread**, set
	directly as the working precision is -- Rich's call, 2026-10-22.
	Radians unless a host says otherwise, and a new thread starts there.
*/
thread_local AngleUnit	geAngleUnit = AU_radians;

// A whole turn in a unit, nought for the radian or for a value that names
//	no unit at all...
static long lzi_turn_of( AngleUnit eUnit )
{
	switch ( eUnit )
	{
	case AU_degrees:	return 360;
	case AU_grads:		return 400;
	case AU_millicles:	return 1000;
	default:			return 0;
	}
}

/*
	**The unit applies to the host's own call only.** The engine calls its
	angle functions itself -- the tangent its sine and cosine, the
	arccosine its arcsine, the cotangent of an arccotangent, gamma's
	reflection the sine -- and each such call must be in radians whatever
	the host has set. So the unit is honoured only by the outermost
	function, the one whose answer guard has just counted one: glRetakeDepth
	counts every CLZAnswerBase alive. A retake re-enters at that depth, so
	an answer taken again wider is taken in the same unit; and lzc_arg,
	which keeps no guard of its own, hands lz_arctangent2 the host's call.

	Asked once, right after the answer guard, and kept: the tie deciders
	are set by it too, before anything nests.
*/
static AngleUnit lzi_angle_unit_here()
{
	return ((glRetakeDepth == 1) && (lzi_turn_of( geAngleUnit ) != 0)) ? geAngleUnit : AU_radians;
}

// The circular functions by index, each pair a function and its
//	complement, so a complement is the index with its low bit flipped...
const long	kAngleSine = 0;
const long	kAngleCosine = 1;
const long	kAngleTangent = 2;
const long	kAngleCotangent = 3;
const long	kAngleSecant = 4;
const long	kAngleCosecant = 5;

// ...and their inverses...
const long	kAngleArcSine = 0;
const long	kAngleArcCosine = 1;
const long	kAngleArcTangent = 2;
const long	kAngleArcCotangent = 3;
const long	kAngleArcSecant = 4;
const long	kAngleArcCosecant = 5;

// x times a whole number, exactly: a whole number adds no places...
static void lzi_times_whole( CLZNumber* pValue, long lWhole, CLZNumber* pResult )
{
	long		l_saved = glWorkingPrecision;
	CLZNumber	whole;

	whole.SetBase( pValue->m_lBase, false );
	whole.Assign( lWhole );

	glWorkingPrecision = pValue->m_lRadix + 2;
	lz_multiply( pValue, &whole, pResult );
	glWorkingPrecision = l_saved;
}

/*
	A whole number of a turn's parts, turn a / b, correctly rounded at the
	working precision, half-even: short where it ends, written to the
	width where it does not.

	**Rounded exactly, in whole numbers.** With P the places, q and r the
	whole quotient and remainder of turn a b^P by b, the answer is q / b^P
	or a unit more, as 2r is above b or, on a tie, as q's last digit is
	odd -- the kept digit, which in an odd base is not the number's
	parity. A division would not do: it writes a place or two past its
	width, and rounding that again goes wrong on the one tie these parts
	have, a sine of 1/2 in an odd base, wherever the kept digit is not the
	repeating one -- at nought places in base three, sin 30 degrees came
	back 1 for 0, from a division twelve places wider; rounding a plain
	division again got every one of base two's 480 thirds wrong, and none
	in six other small bases.
*/
static void lzi_part_of_turn( long lTurn, long lA, long lB, long lBase, CLZNumber* pResult )
{
	long	l_places = glWorkingPrecision;
	long	l_numerator = lTurn * lA;
	bool	b_negative = (l_numerator < 0);

	CLZNumber	scaled;
	CLZNumber	divisor;
	CLZNumber	whole;
	CLZNumber	product;
	CLZNumber	remainder;
	CLZNumber	twice;

	// turn a b^P, a whole number...
	scaled.SetBase( lBase, false );
	scaled.Assign( b_negative ? -l_numerator : l_numerator );
	lzi_scale_by_base_power( &scaled, -l_places );

	divisor.SetBase( lBase, false );
	divisor.Assign( lB );

	// ...its whole quotient and remainder by b, exactly...
	lz_divide_integer( &scaled, &divisor, &whole );
	lzi_times_whole( &whole, lB, &product );
	lz_subtract( &scaled, &product, &remainder );
	lzi_times_whole( &remainder, 2, &twice );

	long	l_order = lz_compare( &twice, &divisor );
	bool	b_up = (l_order == kGreaterThan);

	// ...a tie to the even kept digit...
	if ( l_order == kEqualTo )
	{
		CLZNumber	base;
		CLZNumber	digit;
		CLZNumber	two;
		CLZNumber	parity;

		base.SetBase( lBase, false );
		base.Assign( lBase );
		two.SetBase( lBase, false );
		two.Assign( 2L );

		lz_remainder( &whole, &base, &digit );
		lz_remainder( &digit, &two, &parity );
		parity.CheckForZero();
		b_up = !(parity.IsZero());
	}

	if ( b_up )
	{
		CLZNumber	one;

		one.SetBase( lBase, false );
		one.Assign( 1L );
		lz_add( &whole, &one );
	}

	// ...back over b^P, short only where nothing was left over...
	lzi_scale_by_base_power( &whole, l_places );

	remainder.CheckForZero();

	if ( remainder.IsZero() )
		whole.TrimTrailing();

	if ( b_negative )
		lz_negate( &whole );

	pResult->Assign( &whole );
}

/*
	**x = q quarters + phi: a circular function of x is that function (0)
	or its complement (1) of phi, negated or not, by quadrant** --

		sin	  sin,  cos, -sin, -cos		tan	  tan, -cot,  tan, -cot
		cos	  cos, -sin, -cos,  sin		cot	  cot, -tan,  cot, -tan
		sec	  sec, -csc, -sec,  csc		csc	  csc,  sec, -csc, -sec

	Out of lzi_circular_in_unit since stage 3, 2026-10-22, which shares it
	with the sine and cosine pair.
*/
static const struct { long m_lComplement; bool m_bNegate; } sl_angleQuadrant[ 6 ][ 4 ] =
{
	{ { 0, false }, { 1, false }, { 0, true  }, { 1, true  } },
	{ { 0, false }, { 1, true  }, { 0, true  }, { 1, false } },
	{ { 0, false }, { 1, true  }, { 0, false }, { 1, true  } },
	{ { 0, false }, { 1, true  }, { 0, false }, { 1, true  } },
	{ { 0, false }, { 1, true  }, { 0, true  }, { 1, false } },
	{ { 0, false }, { 1, false }, { 0, true  }, { 1, true  } },
};

/*
	**An angle in a unit brought exactly to theta in [0, an eighth of a
	turn]**: whole turns off by an exact remainder, then whole quarters,
	their count the quadrant, then past an eighth the complement of the
	rest -- a quarter and an eighth of a turn being whole in every unit.
	lzi_function_of_theta then says which function of theta a function of
	the angle is. Out of lzi_circular_in_unit since stage 3, 2026-10-22.
*/
static long lzi_reduce_in_unit( CLZNumber* pAngle, AngleUnit eUnit, CLZNumber* pTheta,
								long* plQuadrant, bool* pbComplement )
{
	long	l_base = pAngle->m_lBase;
	long	l_turn = lzi_turn_of( eUnit );

	CLZNumber	turn;
	CLZNumber	quarter;
	CLZNumber	eighth;

	turn.SetBase( l_base, false );
	turn.Assign( l_turn );
	quarter.SetBase( l_base, false );
	quarter.Assign( l_turn / 4 );
	eighth.SetBase( l_base, false );
	eighth.Assign( l_turn / 8 );

	// Whole turns off, exactly, into [0, a turn)...
	long	lResult = lz_remainder( pAngle, &turn, pTheta );

	pTheta->CheckForZero();

	if ( pTheta->IsNegative() )
		lz_add( pTheta, &turn );

	// ...whole quarters off, their count the quadrant...
	long	l_quadrant = 0;

	while ( (l_quadrant < 3) && (lz_compare( pTheta, &quarter ) != kLessThan) )
	{
		lz_subtract( pTheta, &quarter );
		l_quadrant++;
	}

	// ...and past an eighth, the complement of the rest...
	bool	b_complement = (lz_compare( pTheta, &eighth ) == kGreaterThan);

	if ( b_complement )
	{
		CLZNumber	complement;

		lz_subtract( &quarter, pTheta, &complement );
		pTheta->Assign( &complement );
	}

	pTheta->CheckForZero();

	*plQuadrant = l_quadrant;
	*pbComplement = b_complement;

	return lResult;
}

// Which circular function of theta lFunction of the angle is, by its
//	quadrant and complement, and whether negated...
static long lzi_function_of_theta( long lFunction, long lQuadrant, bool bComplement, bool* pbNegate )
{
	long	l_function = lFunction ^ sl_angleQuadrant[ lFunction ][ lQuadrant ].m_lComplement;

	*pbNegate = sl_angleQuadrant[ lFunction ][ lQuadrant ].m_bNegate;

	return bComplement ? (l_function ^ 1) : l_function;
}

/*
	**Theta in a unit as radians, theta 2pi / turn, worked to lWidth
	places**: 2pi two places wider, the product, and the quotient by the
	turn at lWidth. The caller sizes lWidth -- theta's own noughts past
	the answer it wants, so a small angle keeps its digits. Out of
	lzi_circular_in_unit since stage 3, 2026-10-22.
*/
static void lzi_unit_to_radians( CLZNumber* pTheta, AngleUnit eUnit, long lWidth, CLZNumber* pRadians )
{
	long	l_saved = glWorkingPrecision;
	long	l_base = pTheta->m_lBase;

	glWorkingPrecision = lWidth + 2;

	CLZNumber	twoPi( kTwoPi );
	CLZNumber	scaled;

	lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, glWorkingPrecision );
	lz_multiply( pTheta, &twoPi, &scaled );

	// A division answers at its target's precision, so the target is made
	//	at the width...
	glWorkingPrecision = lWidth;

	CLZNumber	turn;
	CLZNumber	radians;

	turn.SetBase( l_base, false );
	turn.Assign( lzi_turn_of( eUnit ) );
	lz_divide( &scaled, &turn, &radians );

	pRadians->Assign( &radians );

	glWorkingPrecision = l_saved;
}

/*
	**A circular function of an angle in a unit**: sin, cos, tan, cot, sec
	or csc of x degrees, grads or millicles. Called by each of the six,
	in the host's unit, from just after its answer guard.

	**Whole turns come off exactly, and so does everything else until pi
	is needed.** A remainder is exact, and a quarter and an eighth of a
	turn are whole in each unit, so x is brought to q quarters and phi in
	[0, a quarter), then phi to theta in [0, an eighth] by the complement
	where it is past one, all exactly. The function of x is then a sign
	times the function or its complement of theta -- the table below --
	and only theta is converted, which keeps every function away from its
	poles but cot and csc near nought.

	**Exact where the mathematics is**: at theta of nought, the radian
	function at an exact nought -- 0, 1, or the pole cot(0) and csc(0)
	answer -- and at a twelfth of a turn the sine 1/2 and the cosecant 2,
	at an eighth the tangent and cotangent 1. Every other rational
	multiple of a turn has an irrational sine (Niven's theorem), so
	nothing else is exact and nothing else can lie on a half.

	**Otherwise theta becomes theta 2pi / turn, wide enough not to
	matter**: kAngleUnitGuardDigits past the answer, theta's own leading
	noughts again so a small angle keeps its digits, and those twice more
	for cot and csc, whose answer is one over it. The radian function of
	that is taken kAngleUnitGuardDigits wider than the answer and rounded
	once, by the function's own rule: places for sin, cos, sec and csc,
	the quotient's for tan and cot.
*/

static long lzi_circular_in_unit( long lFunction, CLZNumber* leftOperand, CLZNumber* targetNumber,
								  CLZAnswerBase& answer, AngleUnit eUnit )
{
	static const LZUnaryFunction	sl_radians[] = { lz_sine, lz_cosine, lz_tangent, lz_cotangent, lz_secant, lz_cosecant };

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	// A not-a-number or an infinity is the radian function's to answer, the
	//	same in any unit...
	if ( leftOperand->FlagIsOn( kfQNaN ) || leftOperand->FlagIsOn( kfInfinity ) )
		return sl_radians[ lFunction ]( leftOperand, targetNumber );

	long	l_base = leftOperand->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_turn = lzi_turn_of( eUnit );
	long	lResult = kNoError;

	CLZNumber	turn;
	CLZNumber	eighth;
	CLZNumber	theta;
	long		l_quadrant = 0;
	bool		b_complement = false;

	turn.SetBase( l_base, false );
	turn.Assign( l_turn );
	eighth.SetBase( l_base, false );
	eighth.Assign( l_turn / 8 );

	// Into [0, an eighth], exactly -- see lzi_reduce_in_unit...
	lResult = lzi_first_failure( lResult, lzi_reduce_in_unit( leftOperand, eUnit, &theta, &l_quadrant, &b_complement ) );

	bool	b_negate = false;
	long	l_function = lzi_function_of_theta( lFunction, l_quadrant, b_complement, &b_negate );

	CLZNumber	value;
	bool		b_exact = true;

	value.SetBase( l_base, false );

	if ( theta.IsZero() )
	{
		// 0, 1, or a pole, as the radian function answers an exact nought...
		CLZNumber	nought;

		nought.SetBase( l_base, false );
		nought.Assign( 0L );
		lResult = lzi_first_failure( lResult, sl_radians[ l_function ]( &nought, &value ) );
		answer.ExactAndShort();
	}
	else
	{
		CLZNumber	twelve;

		lzi_times_whole( &theta, 12, &twelve );

		bool	b_twelfth = (lz_compare( &twelve, &turn ) == kEqualTo);
		bool	b_eighth = (lz_compare( &theta, &eighth ) == kEqualTo);

		if ( b_twelfth && (l_function == kAngleSine) )
		{
			// A half, which ends only in an even base, divided so it rounds
			//	correctly where it does not...
			lzi_part_of_turn( 1, 1, 2, l_base, &value );
			answer.Exact();
		}
		else if ( (b_twelfth && (l_function == kAngleCosecant)) ||
				  (b_eighth && ((l_function == kAngleTangent) || (l_function == kAngleCotangent))) )
		{
			value.Assign( (l_function == kAngleCosecant) ? 2L : 1L );
			answer.ExactAndShort();
		}
		else
			b_exact = false;
	}

	if ( !b_exact )
	{
		long	l_guard = lzi_guard_places( kAngleUnitGuardDigits, l_base );

		// Theta's noughts again, and the factor's three at most -- 2pi over a
		//	thousand is 0.00628 -- and twice for cot and csc...
		long	l_lead = lzi_leading_exponent( &theta );
		long	l_small = ((l_lead < 0) ? -l_lead : 0) + 3;
		long	l_width = l_saved + l_guard + l_small;

		if ( (l_function == kAngleCotangent) || (l_function == kAngleCosecant) )
			l_width += l_small;

		CLZNumber	radians;

		lzi_unit_to_radians( &theta, eUnit, l_width, &radians );

		// ...and the radian function of it, wider than the answer...
		glWorkingPrecision = l_saved + l_guard;
		lResult = lzi_first_failure( lResult, sl_radians[ l_function ]( &radians, &value ) );
	}

	if ( b_negate && value.FlagIsOff( kfInfinity ) && value.FlagIsOff( kfQNaN ) )
		lz_negate( &value );

	glWorkingPrecision = l_saved;

	// ...rounded once by the function's own rule...
	if ( !b_exact )
	{
		if ( (lFunction == kAngleTangent) || (lFunction == kAngleCotangent) )
			lzi_round_quotient_answer( &value );
		else
		{
			value.SetPrecision( glWorkingPrecision );
			value.Round();
		}
	}

	if ( !(p_target->IsConstant()) )
		p_target->Assign( &value );

	return lResult;
}

/*
	**An inverse circular function, answering in a unit**: arcsin, arccos,
	arctan, arccot, arcsec or arccsc in degrees, grads or millicles.

	**Exact where the mathematics is**: at the arguments whose angle is a
	whole number of twelfths or eighths of a turn -- nought, a half, one
	and two, either sign, as each function's domain allows -- the angle is
	that part of a turn, divided out and so correctly rounded, short where
	it ends: arcsin(1/2) is 30 degrees, 33.3... grads. These are the only
	rational arguments with a rational angle (Niven's theorem again), so
	nothing else can lie on a half.

	**Otherwise the radian function, kAngleUnitGuardDigits and three wider
	-- the turn over 2pi reaches 159 -- times the turn over 2pi**, rounded
	once at the places. A refusal is the radian function's, answered as it
	answers it.
*/

static long lzi_inverse_in_unit( long lFunction, CLZNumber* leftOperand, CLZNumber* targetNumber,
								 CLZAnswerBase& answer, AngleUnit eUnit )
{
	static const LZUnaryFunction	sl_radians[] = { lz_arcsine, lz_arccosine, lz_arctangent, lz_arccotangent,
													 lz_arcsecant, lz_arccosecant };

	// x = n / d answers turn a / b...
	struct SLZExactAngle { long m_lN; long m_lD; long m_lA; long m_lB; };

	static const SLZExactAngle	sl_arcsine[] = { { 0, 1, 0, 1 }, { 1, 2, 1, 12 }, { -1, 2, -1, 12 }, { 1, 1, 1, 4 }, { -1, 1, -1, 4 } };
	static const SLZExactAngle	sl_arccosine[] = { { 1, 1, 0, 1 }, { 1, 2, 1, 6 }, { 0, 1, 1, 4 }, { -1, 2, 1, 3 }, { -1, 1, 1, 2 } };
	static const SLZExactAngle	sl_arctangent[] = { { 0, 1, 0, 1 }, { 1, 1, 1, 8 }, { -1, 1, -1, 8 } };
	static const SLZExactAngle	sl_arccotangent[] = { { 0, 1, 1, 4 }, { 1, 1, 1, 8 }, { -1, 1, -1, 8 } };
	static const SLZExactAngle	sl_arcsecant[] = { { 1, 1, 0, 1 }, { 2, 1, 1, 6 }, { -2, 1, 1, 3 }, { -1, 1, 1, 2 } };
	static const SLZExactAngle	sl_arccosecant[] = { { 1, 1, 1, 4 }, { 2, 1, 1, 12 }, { -2, 1, -1, 12 }, { -1, 1, -1, 4 } };

	static const struct { const SLZExactAngle* m_pAngles; long m_lCount; } sl_exact[] =
	{
		{ sl_arcsine, 5 }, { sl_arccosine, 5 }, { sl_arctangent, 3 },
		{ sl_arccotangent, 3 }, { sl_arcsecant, 4 }, { sl_arccosecant, 4 },
	};

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	// A not-a-number or an infinity is the radian function's to answer...
	if ( leftOperand->FlagIsOn( kfQNaN ) || leftOperand->FlagIsOn( kfInfinity ) )
		return sl_radians[ lFunction ]( leftOperand, targetNumber );

	long	l_base = leftOperand->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_turn = lzi_turn_of( eUnit );

	// The arguments with an exact angle, x d = n worked exactly...
	for ( long i = 0; i < sl_exact[ lFunction ].m_lCount; i++ )
	{
		const SLZExactAngle&	exact = sl_exact[ lFunction ].m_pAngles[ i ];
		CLZNumber				times;
		CLZNumber				n;

		lzi_times_whole( leftOperand, exact.m_lD, &times );
		n.SetBase( l_base, false );
		n.Assign( exact.m_lN );

		if ( lz_compare( &times, &n ) == kEqualTo )
		{
			CLZNumber	value;

			lzi_part_of_turn( l_turn, exact.m_lA, exact.m_lB, l_base, &value );
			answer.Exact();

			if ( !(p_target->IsConstant()) )
				p_target->Assign( &value );

			return kNoError;
		}
	}

	// ...and the rest through radians, wider...
	glWorkingPrecision = l_saved + lzi_guard_places( kAngleUnitGuardDigits + 3, l_base );

	CLZNumber	theta;
	long		lResult = sl_radians[ lFunction ]( leftOperand, &theta );

	if ( (lResult == kNoError) && theta.FlagIsOff( kfQNaN ) && theta.FlagIsOff( kfInfinity ) )
	{
		CLZNumber	twoPi( kTwoPi );
		CLZNumber	scaled;
		CLZNumber	angle;

		lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, glWorkingPrecision + 2 );
		lzi_times_whole( &theta, l_turn, &scaled );
		lz_divide( &scaled, &twoPi, &angle );

		theta.Assign( &angle );
		glWorkingPrecision = l_saved;
		theta.SetPrecision( glWorkingPrecision );
		theta.Round();
	}

	glWorkingPrecision = l_saved;

	if ( !(p_target->IsConstant()) )
		p_target->Assign( &theta );

	return lResult;
}


/*
	The reciprocal of a value, by Newton, to a given number of places.

		y  <-  y * (2 - d*y)

	converges on 1/d and doubles its correct digits each pass, using
	nothing but multiplication and subtraction. That is the whole
	point: multiplication has a transform tier and long division does
	not, so past a few hundred places a reciprocal and a multiply beat
	dividing digit by digit by an order of magnitude and more.

	Three things about the shape of it, each measured rather than
	assumed:

	  - The divisor is scaled into [1, base) by a whole power of the
	    base first, which is exact -- it moves the radix, not the
	    digits -- and keeps the iteration's arithmetic near one.

	  - The starting estimate comes from the engine's own division at
	    a dozen places. Long division is only slow when the quotient is
	    long, so a short one is the cheapest good estimate available,
	    and it saves the half-dozen passes that climbing from a crude
	    guess would cost.

	  - The working precision doubles with the convergence rather than
	    sitting at the final width. A pass that knows twelve digits has
	    no business computing ten thousand, and letting the width track
	    the accuracy made the whole thing three times quicker.

	The pass count is arithmetic rather than a convergence test.
	Testing for two equal passes keeps iterating after it has arrived,
	because the last place of a value carrying spare digits oscillates,
	and a wasted pass at the final width costs as much as everything
	before it -- that was worth another eight times.
*/

static void lzi_reciprocal_newton( CLZNumber* divisor, long places, CLZNumber* target )
{
	long	l_exponent = lzi_leading_exponent( divisor );

	CLZNumber	l_scaled;
	l_scaled.Assign( divisor );
	lzi_scale_by_base_power( &l_scaled, l_exponent );

	/*
		**In the divisor's base, or the iteration runs in base ten.**
		These were made in the default base, and every mixed-base
		operation takes its left operand's base, so the seed division
		and each pass after it quietly worked in base ten -- and the
		answer was then scaled back by the divisor's exponent in its own
		base, multiplying a base-ten number by eight to the twenty-first
		where it wanted ten to the twenty-first. A base-eight quotient
		came out 108.4 times too small, which is (10/8)^21 exactly, and
		silently: no flag, no code. Every base but ten was wrong this way
		from the day the reciprocal went in, until 2026-09-21.
	*/

	CLZNumber	one;
	CLZNumber	two;

	one.SetBase( divisor->m_lBase, false );
	two.SetBase( divisor->m_lBase, false );
	one.Assign( 1L );
	two.Assign( 2L );

	CLZNumber	y;
	y.SetPrecision( kReciprocalSeedPlaces );

	gSeekingReciprocal = true;
	lz_divide( &one, &l_scaled, &y );
	gSeekingReciprocal = false;

	long	l_saved = glWorkingPrecision;
	long	l_target = places + kReciprocalGuardPlaces;
	long	l_width = kReciprocalSeedPlaces + 4;

	CLZNumber	product;
	CLZNumber	correction;
	CLZNumber	next;

	/*
		**Each pass multiplies by the divisor to its own width**, since
		2026-10-16 -- part (c) of performance item 7 of the review of
		2026-10-15. A pass knowing w places of 1/d needs d to w places and
		a guard, not all of it: y converges on the reciprocal of the
		divisor it is shown, and a divisor right to w and a guard has a
		reciprocal right to w. It multiplied by the whole divisor every
		pass, and a long divisor over a short quotient paid for its length
		at every width. The division corrects the quotient by its exact
		remainder afterwards, so nothing about the answer depends on this;
		a divisor no longer than the pass is used as it stands.
	*/

	CLZNumber	l_near;

	auto	divisor_at = [ &l_scaled, &l_near ]( long lWidth ) -> CLZNumber*
	{
		if ( l_scaled.m_lRadix <= (lWidth + kReciprocalDivisorGuard) )
			return &l_scaled;

		l_near.Assign( &l_scaled );
		l_near.Round( kDefaultRoundingMode, lWidth + kReciprocalDivisorGuard );

		return &l_near;
	};

	// Read as each pass multiplies...
	auto	widest = []( CLZNumber* pDivisor ) -> CLZNumber*
	{
		if ( pDivisor->m_lDigitsLength > gReciprocalWidest )
			gReciprocalWidest = pDivisor->m_lDigitsLength;

		return pDivisor;
	};

	gReciprocalWidest = 0;

	while ( l_width < l_target )
	{
		l_width = ((l_width * 2) < l_target) ? (l_width * 2) : l_target;

		glWorkingPrecision = l_width;

		lz_multiply( widest( divisor_at( l_width ) ), &y, &product );
		lz_subtract( &two, &product, &correction );
		lz_multiply( &y, &correction, &next );

		next.SetPrecision( l_width );
		next.Round();

		y.Assign( &next );
	}

	// Two settling passes at the full width, one more than the
	//	arithmetic needs...
	for ( int i = 0; i < 2; i++ )
	{
		glWorkingPrecision = l_target;

		lz_multiply( widest( divisor_at( l_target ) ), &y, &product );
		lz_subtract( &two, &product, &correction );
		lz_multiply( &y, &correction, &next );

		next.SetPrecision( l_target );
		next.Round();

		y.Assign( &next );
	}

	glWorkingPrecision = l_saved;

	// 1/d = (1/scaled) divided by the same power of the base...
	lzi_scale_by_base_power( &y, l_exponent );

	target->Assign( &y );
}


/*
	One digit of long division, found by estimate rather than by trial.

	The window is the part of the dividend the divisor is being taken
	from, one digit longer than the divisor; the digit wanted is how
	many whole divisors it holds, and the remainder is left in it.

	It was found by subtracting the divisor until the window went
	negative and then undoing the last subtraction -- about five and a
	half passes over the divisor per digit in base ten, up to sixteen in
	base sixteen, and nearly four fifths of the time a division took.

	Now the digit is estimated from the leading digits, and that many
	divisors are taken away in one pass. **The estimate is never low and
	at most one too high**, in every base, provided it reads
	kDivideLeadDigits of the divisor against one more of the window:
	the divisor's own leading digit is never nought -- lzi_divide_op
	trims it -- so the error is under the base over the base squared,
	which is under one. A divisor no longer than that is read whole and
	the estimate is exact.

	**The two correction loops are what make that bound safe to lean
	on.** An overshoot is taken back one divisor at a time, and whatever
	is left is compared against the divisor and taken from again while
	it still holds one. The second should never run. If the bound is
	ever wrong the answer is still exact and only the time moves, which
	is why every quotient came back byte-identical to the trial loop's
	rather than merely close.

	Measured over 266,614 quotient digits across bases 2, 8, 10 and 16:
	3,039 overshoots, 1.1%, and no shortfall at all. 2.7 times faster
	than the trial loop at 36 places, 5.4 at 390.

	**No correctness test can see how good the estimate is**, and that
	is by design rather than a gap. Making it deliberately one low, or
	one high, or reading a single lead digit, each passes every test and
	the sweep; so does removing the shortfall loop. Only removing the
	overshoot loop is caught -- the suite hangs. The estimate is a cost,
	and a timing is the only thing that measures it.

	Three lead digits, and two would do: the bound holds from two in
	every base. One is where it goes wrong, and only for some divisors
	-- a leading 1 followed by nines, where reading "1" for a value
	nearly twice that overshoots by about as much. At 36 places, dividing
	by 1999999.87 took 3,776ns with one lead digit against 2,427 with two
	and 2,438 with three; four measured the same again. Three is kept
	because it overshoots less often and costs nothing measurable.
*/

const long kDivideLeadDigits = 3;

// Take the divisor from the window once. The window must hold it...
static void lzi_window_subtract( long* lpWindow, const long* lpDivisor, long lLength, long lBase )
{
	long	l_borrow = 0;

	for ( long i = 0; i <= lLength; i++ )
	{
		long	l_value = lpWindow[ i ] - l_borrow - ((i < lLength) ? lpDivisor[ i ] : 0);

		if ( l_value < 0 )
		{
			l_value += lBase;
			l_borrow = 1;
		}
		else
			l_borrow = 0;

		lpWindow[ i ] = l_value;
	}
}

// Add the divisor back over the window, answering the carry out of its top...
static long lzi_window_add( long* lpWindow, const long* lpDivisor, long lLength, long lBase )
{
	long	l_carry = 0;

	for ( long i = 0; i <= lLength; i++ )
	{
		long	l_value = lpWindow[ i ] + l_carry + ((i < lLength) ? lpDivisor[ i ] : 0);

		if ( l_value >= lBase )
		{
			l_value -= lBase;
			l_carry = 1;
		}
		else
			l_carry = 0;

		lpWindow[ i ] = l_value;
	}

	return l_carry;
}

// Whether the window, read as a number, still holds the divisor. Settled
//	by the first digit the two differ in, which is nearly always the top...
static bool lzi_window_holds( const long* lpWindow, const long* lpDivisor, long lLength )
{
	if ( lpWindow[ lLength ] != 0 )
		return true;

	for ( long i = lLength - 1; i >= 0; i-- )
	{
		if ( lpWindow[ i ] != lpDivisor[ i ] )
			return (lpWindow[ i ] > lpDivisor[ i ]);
	}

	return true;
}

static long lzi_quotient_digit( long* lpWindow, const long* lpDivisor, long lLength, long lBase )
{
	long	l_lead = (lLength < kDivideLeadDigits) ? lLength : kDivideLeadDigits;
	long	l_divisorLead = 0;
	long	l_windowLead = lpWindow[ lLength ];

	for ( long j = 1; j <= l_lead; j++ )
	{
		l_divisorLead = (l_divisorLead * lBase) + lpDivisor[ lLength - j ];
		l_windowLead = (l_windowLead * lBase) + lpWindow[ lLength - j ];
	}

	long	l_digit = l_windowLead / l_divisorLead;

	if ( l_digit > (lBase - 1) )
		l_digit = lBase - 1;

	// Take that many divisors away in one pass...
	long	l_borrow = 0;

	if ( l_digit > 0 )
	{
		for ( long i = 0; i <= lLength; i++ )
		{
			long	l_product = l_borrow + ((i < lLength) ? (lpDivisor[ i ] * l_digit) : 0);
			long	l_value = lpWindow[ i ] - (l_product % lBase);

			l_borrow = l_product / lBase;

			if ( l_value < 0 )
			{
				l_value += lBase;
				l_borrow++;
			}

			lpWindow[ i ] = l_value;
		}
	}

	// ...give one back for each the estimate overshot...
	while ( l_borrow > 0 )
	{
		l_borrow -= lzi_window_add( lpWindow, lpDivisor, lLength, lBase );
		l_digit--;
	}

	// ...and take one more for each it fell short, which it should not...
	while ( lzi_window_holds( lpWindow, lpDivisor, lLength ) )
	{
		lzi_window_subtract( lpWindow, lpDivisor, lLength, lBase );
		l_digit++;
	}

	return l_digit;
}


/*
	Long division on packed limbs.

	The loop above runs once per digit of the divisor for every digit of
	the quotient, because every digit lives in its own long. Packing k
	digits of base b into one limb makes it one digit of base b^k -- a
	regrouping, not a conversion, so the base is respected at the
	boundaries and nowhere else -- and the loop runs k times fewer passes
	of k times fewer steps.

	**Every limb is below 2^32**, so a product of two, plus a carry, fits
	in 64 bits and the estimate is a 64-by-32 division. Nothing here
	needs __uint128_t, which MSVC does not have -- and since 2026-09-29
	nothing in the transform tier does either. That is what sets the table: 31, 10, 9, 8 and 7 digits for
	bases 2, 8, 10, 12 and 16. Base twelve is the awkward one -- 12^9 is
	over the line, so it takes eight where base ten takes nine.

	**It is Knuth's Algorithm D, normalised**, not the estimate the loop
	above uses. That one reads three lead digits unnormalised, and its
	error grows with the base it works in: fine at sixteen, useless at
	a billion. Normalised, Knuth's estimate is at most two too high
	whatever the limb base, two refinements take it to at most one, and
	the add-back takes that one back.

	**The add-back is rarer than it looks and commoner than the
	textbook says.** About twice in every limb base on random operands,
	which at a billion is never: stage one's 5,943 random long divisions
	reached it not once. But a divisor a whisker off a power of the base
	-- 1.0000000000000000005, 0.9999999999999999999999995 -- normalises
	to limbs of all nines, and that is the shape the estimate overshoots
	on. The engine makes it whenever it divides values near one: 25 times
	in the suite and 198 in the sweep. The tests reach it both ways, by
	that shape and by shrinking the limb; see lzi_set_packed_divide_factor.

	Proved before it came here, the way the transform kernel was: over
	12,660 divisions against Python's own integers, in five bases, at the
	table's limbs and at one digit per limb, and by mutation. The
	harness drove this function rather than a copy.

	**Where lzi_divide_op calls it, it answers false only when it cannot
	allocate**, and the loop above is taken instead, needing nothing.
*/

/*
	Where packing starts paying: the product of quotient digits and
	divisor digits, since that is the work the loop above does. Packing
	is a linear pass, so small divisions lose by it.

	**Measured in context**, lz_divide with the path forced each way,
	because in isolation the allocation this path
	makes is cheap and inside lzi_divide_op an allocation need not be.
	The floor is set by two things: it sits past the random shapes'
	crossing, and the worst shape measured loses no more than about a
	tenth.

	**The shape to know is a divisor that is a power of the base**,
	whose quotient is all noughts but one -- the loop's best case,
	since a nought digit skips the subtraction. In bases 8 to 16 from
	32, random shapes pack 1.03 to 3 times faster and that one loses up
	to 6 per cent; the worst single shape measured, all the largest
	digit over the same, lost 10 once. Below 32 random shapes lose too,
	as low as 0.90.

	**Base two wants its own floor**, because a quotient digit there is
	nought or one, so about half of them are the loop's cheap case. At
	256, random shapes pack at least 1.5 times faster and three times
	typically, and the power of two loses at most 3 per cent; 128 would
	sit on the crossing, the random ones' worst there being 0.99. It was
	512 until that was measured in context -- the same trade the other
	bases were already making at 32, and a milder one.

	Neither floor is a serious gate. Over the suite and the sweep, 99.9
	and 100 per cent of the loop's work was in divisions past it.

	**Measured again 2026-10-17, and both came down: 23, and 51 in base
	two** -- Rich's choice of two floors over one per base. The packed
	path had grown cheaper since: pack and unpack compiled for their base
	(2026-09-29) and every base given limbs (2026-10-12). lz_divide timed
	both ways over 1,950 shapes a base in all fifteen, divisors of one to
	sixteen digits, random and the four awkward shapes above, twice, each
	shape's ratio the two passes' geometric mean. From 23, every random
	shape packs at least as fast in every base -- 1.17 to 1.50 times at
	the median of the band it moves -- and the worst shape loses 7 per
	cent; in base two from 51, random shapes 1.01 to 3.9 times, 1.67 at
	the median, the worst losing 2. A one-digit divisor was what showed
	it: a whole quotient of 20 to 26 digits over 7 packed faster and sat
	below 32. A floor per base would move half again as many shapes, but
	below about 15 the per-base figures move with the noise allowed.
*/

const long kPackedDivideFloor = 23;
const long kPackedDivideFloorBinary = 51;

// The floor and the factor in force, which only the tests change, and
//	what the tests read back to know the paths were taken. The settings
//	are the process's and the counts each thread's own...
static long gPackedDivideFloor = kPackedFloorMeasured;
static long gPackedDivideFactor = 0;
static thread_local long gPackedDivisions = 0;
static thread_local long gPackedAddBacks = 0;

long lzi_set_packed_divide_floor( long lProduct )
{
	long	l_previous = gPackedDivideFloor;

	gPackedDivideFloor = lProduct;

	return l_previous;
}

long lzi_set_packed_divide_factor( long lFactor )
{
	long	l_previous = gPackedDivideFactor;

	gPackedDivideFactor = lFactor;

	return l_previous;
}

void lzi_packed_divide_counts( long* plDivisions, long* plAddBacks )
{
	if ( plDivisions != nullptr )
		*plDivisions = gPackedDivisions;

	if ( plAddBacks != nullptr )
		*plAddBacks = gPackedAddBacks;
}

/*
	Digits per limb in the packed divider's table and the packed
	multiply's, which differ only in base two. Named because pack and
	unpack below are compiled for exactly these: one name for each
	number, so a table cannot change without the code compiled for it.
*/

constexpr long kPackedDigitsBinaryDivide = 31;
constexpr long kPackedDigitsBinaryMultiply = 20;
constexpr long kPackedDigitsOctal = 10;
constexpr long kPackedDigitsDecimal = 9;
constexpr long kPackedDigitsDuodecimal = 8;
constexpr long kPackedDigitsHex = 7;

/*
	**Digits per limb in any base**: the most k for which b^k is below
	2^32, which keeps every product of two limbs inside 64 bits. It gives
	31, 10, 9, 8 and 7 for the five bases above -- their tables are this
	rule, written out so that pack and unpack can be compiled for them --
	and 20, 15, 13, 12, 11, 10, 9, 8, 8 and 8 for bases three to fifteen,
	which had no table until 2026-10-12. A multiplication stops lower, at
	2^30 -- see lzi_multiply_limb_digits_for. Both are worked out at
	compile time.
*/

constexpr long lzi_limb_digits_for( uint64_t lBase )
{
	long		k = 0;
	uint64_t	v = 1;

	while ( (v * lBase) < (1ULL << 32) )
	{
		v *= lBase;
		k++;
	}

	return k;
}

/*
	**A multiplication's limbs stop at 2^30, not 2^32.** Its column loop
	sums products and folds the sum before it can pass 2^64, and a fold
	is a division: a limb near 2^32 leaves room for one to four products
	between folds, and so it multiplied slower than a digit fewer. Base
	three at 4,000 digits took 110 us at twenty trits a limb and 30 at
	eighteen, base nine 307 against 64, base fifteen 273 against 86. A
	limb within 2^30 leaves room for sixteen, and that rule picked the
	fastest factor measured in every base from three to fifteen but five,
	whose factors differed by six per cent; it is also what the tables
	give for eight, ten, twelve and sixteen. Base two packs twenty bits,
	measured faster again, since packing it is a bit at a time.
*/

constexpr long lzi_multiply_limb_digits_for( uint64_t lBase )
{
	long		k = 0;
	uint64_t	v = 1;

	while ( (v * lBase) <= (1ULL << 30) )
	{
		v *= lBase;
		k++;
	}

	return k;
}

// The named tables are these rules, bar base two's multiply, measured
//	faster at twenty bits...
static_assert( kPackedDigitsBinaryDivide == lzi_limb_digits_for( 2 ), "base two's limbs" );
static_assert( kPackedDigitsOctal == lzi_limb_digits_for( 8 ), "base eight's limbs" );
static_assert( kPackedDigitsDecimal == lzi_limb_digits_for( 10 ), "base ten's limbs" );
static_assert( kPackedDigitsDuodecimal == lzi_limb_digits_for( 12 ), "base twelve's limbs" );
static_assert( kPackedDigitsHex == lzi_limb_digits_for( 16 ), "base sixteen's limbs" );
static_assert( (kPackedDigitsOctal == lzi_multiply_limb_digits_for( 8 )) &&
			   (kPackedDigitsDecimal == lzi_multiply_limb_digits_for( 10 )) &&
			   (kPackedDigitsDuodecimal == lzi_multiply_limb_digits_for( 12 )) &&
			   (kPackedDigitsHex == lzi_multiply_limb_digits_for( 16 )), "a multiplication's limbs" );

static long lzi_limb_digits( long lBase )
{
	static constexpr std::array<long, 17>	sl_digits =
	{
		0, 0, lzi_limb_digits_for( 2 ), lzi_limb_digits_for( 3 ), lzi_limb_digits_for( 4 ),
		lzi_limb_digits_for( 5 ), lzi_limb_digits_for( 6 ), lzi_limb_digits_for( 7 ),
		lzi_limb_digits_for( 8 ), lzi_limb_digits_for( 9 ), lzi_limb_digits_for( 10 ),
		lzi_limb_digits_for( 11 ), lzi_limb_digits_for( 12 ), lzi_limb_digits_for( 13 ),
		lzi_limb_digits_for( 14 ), lzi_limb_digits_for( 15 ), lzi_limb_digits_for( 16 )
	};

	return ((lBase >= 2) && (lBase <= 16)) ? sl_digits[ lBase ] : 0;
}

static long lzi_multiply_limb_digits( long lBase )
{
	static constexpr std::array<long, 17>	sl_digits =
	{
		0, 0, kPackedDigitsBinaryMultiply, lzi_multiply_limb_digits_for( 3 ), lzi_multiply_limb_digits_for( 4 ),
		lzi_multiply_limb_digits_for( 5 ), lzi_multiply_limb_digits_for( 6 ), lzi_multiply_limb_digits_for( 7 ),
		lzi_multiply_limb_digits_for( 8 ), lzi_multiply_limb_digits_for( 9 ), lzi_multiply_limb_digits_for( 10 ),
		lzi_multiply_limb_digits_for( 11 ), lzi_multiply_limb_digits_for( 12 ), lzi_multiply_limb_digits_for( 13 ),
		lzi_multiply_limb_digits_for( 14 ), lzi_multiply_limb_digits_for( 15 ), lzi_multiply_limb_digits_for( 16 )
	};

	return ((lBase >= 2) && (lBase <= 16)) ? sl_digits[ lBase ] : 0;
}

// Digits per limb: the table, or fewer if the tests asked. Never more,
//	which is what keeps a limb below 2^32 whatever they ask...
static long lzi_packed_factor( long lBase )
{
	long	l_factor = 0;

	switch ( lBase )
	{
	case 2:		l_factor = kPackedDigitsBinaryDivide;	break;
	case 8:		l_factor = kPackedDigitsOctal;			break;
	case 10:	l_factor = kPackedDigitsDecimal;		break;
	case 12:	l_factor = kPackedDigitsDuodecimal;		break;
	case 16:	l_factor = kPackedDigitsHex;			break;
	default:	l_factor = lzi_limb_digits( lBase );	break;
	}

	if ( (gPackedDivideFactor > 0) && (gPackedDivideFactor < l_factor) )
		l_factor = gPackedDivideFactor;

	return l_factor;
}

static bool lzi_packed_divide_wanted( long lBase, long lQuotientDigits, long lDivisorDigits )
{
	if ( lzi_packed_factor( lBase ) == 0 )
		return false;

	long	l_floor = gPackedDivideFloor;

	if ( l_floor == kPackedFloorMeasured )
		l_floor = (lBase == 2) ? kPackedDivideFloorBinary : kPackedDivideFloor;

	// In 64 bits: a long is 32 on Windows, and this product need not fit...
	return ((static_cast<int64_t>( lQuotientDigits ) * lDivisorDigits) >= l_floor);
}

// The fewest quotient digits that pack over a divisor of lDivisorDigits at
//	the measured floor, nought in a base that never packs. A cut dividend
//	keeps at least these, or a division the packed path took would fall
//	back to the loop and be slower for doing less -- see the cut in
//	lzi_divide_work. The measured floor and not the one in force, so a
//	test raising the floor still cuts, and reaches the loop with a cut...
static long lzi_packed_divide_digits( long lBase, long lDivisorDigits )
{
	if ( (lzi_packed_factor( lBase ) == 0) || (lDivisorDigits < 1) )
		return 0;

	long	l_floor = (lBase == 2) ? kPackedDivideFloorBinary : kPackedDivideFloor;

	return (l_floor + lDivisorDigits - 1) / lDivisorDigits;
}

/*
	**Pack and unpack are compiled for each tabled base and its table's
	digits per limb**, since 2026-09-29. With both known, the multiply by
	the base going in and the % and / taking digits out become shifts and
	masks in bases two, eight and sixteen and a multiply-high in ten and
	twelve, and each limb's loop unrolls. Moving digits in and out is most
	of a lopsided product, so it shows: 5 x 1000 packed 3.1 times faster in
	base two, 1.9 in sixteen and 1.3 in ten, and the sweep went 208 to
	197 ms.

	One implementation, three ways in: B and F are 0 for "take it at run
	time". A factor the tests force comes in with the base known and the
	factor not; a base outside two to sixteen, which SetBase refuses, with
	neither. Whole limbs go first with nothing to test inside the loop,
	then the part-filled one.
*/

// Digits to limbs, least significant first, the most significant digit
//	of each group highest...
template <uint64_t B, long F>
static long lzi_packed_pack_in( const long* lpDigits, long lLength, long lBase, long lFactor,
								uint64_t* lpLimbs )
{
	const uint64_t	l_base = (B != 0) ? B : static_cast<uint64_t>( lBase );
	const long		l_factor = (F != 0) ? F : lFactor;
	long			l_whole = lLength / l_factor;
	long			l_limbs = ((lLength + l_factor) - 1) / l_factor;

	for ( long i = 0; i < l_whole; i++ )
	{
		const long*	lp_group = lpDigits + (i * l_factor);
		uint64_t	l_value = 0;

		for ( long j = l_factor - 1; j >= 0; j-- )
			l_value = (l_value * l_base) + static_cast<uint64_t>( lp_group[ j ] );

		lpLimbs[ i ] = l_value;
	}

	if ( l_whole < l_limbs )
	{
		const long*	lp_group = lpDigits + (l_whole * l_factor);
		uint64_t	l_value = 0;

		for ( long j = (lLength - (l_whole * l_factor)) - 1; j >= 0; j-- )
			l_value = (l_value * l_base) + static_cast<uint64_t>( lp_group[ j ] );

		lpLimbs[ l_whole ] = l_value;
	}

	// ...answering how many are in use, leading noughts not counted...
	while ( (l_limbs > 0) && (lpLimbs[ l_limbs - 1 ] == 0) )
		l_limbs--;

	return l_limbs;
}

// ...and back again, into exactly the number of digits asked for...
template <uint64_t B, long F>
static void lzi_packed_unpack_in( const uint64_t* lpLimbs, long lLimbs, long lBase, long lFactor,
								  long* lpDigits, long lLength )
{
	const uint64_t	l_base = (B != 0) ? B : static_cast<uint64_t>( lBase );
	const long		l_factor = (F != 0) ? F : lFactor;
	long			i = 0;
	long			l_index = 0;

	for ( ; (i < lLimbs) && ((l_index + l_factor) <= lLength); i++, l_index += l_factor )
	{
		uint64_t	l_value = lpLimbs[ i ];

		for ( long j = 0; j < l_factor; j++ )
		{
			lpDigits[ l_index + j ] = static_cast<long>( l_value % l_base );
			l_value /= l_base;
		}
	}

	// A limb reaching past the digits asked for gives only those...
	if ( (i < lLimbs) && (l_index < lLength) )
	{
		uint64_t	l_value = lpLimbs[ i ];

		for ( ; l_index < lLength; l_index++ )
		{
			lpDigits[ l_index ] = static_cast<long>( l_value % l_base );
			l_value /= l_base;
		}
	}

	// ...and any past the last limb are noughts...
	for ( ; l_index < lLength; l_index++ )
		lpDigits[ l_index ] = 0;
}

/*
	**Every other base is compiled at its own digits per limb too**, since
	2026-10-12, when bases three to fifteen were given limbs: at a
	multiplication's factor and at the divider's, both worked out at
	compile time, and a forced factor takes the base alone. Measured on small products, packing went 3 to 9
	per cent faster in bases three, seven and thirteen than with both at
	run time, a division by a constant being a multiply.
*/

template <uint64_t B>
static long lzi_packed_pack_for( const long* lpDigits, long lLength, long lFactor, uint64_t* lpLimbs )
{
	if ( lFactor == lzi_multiply_limb_digits_for( B ) )
		return lzi_packed_pack_in<B, lzi_multiply_limb_digits_for( B )>( lpDigits, lLength, static_cast<long>( B ), lFactor, lpLimbs );

	if ( lFactor == lzi_limb_digits_for( B ) )
		return lzi_packed_pack_in<B, lzi_limb_digits_for( B )>( lpDigits, lLength, static_cast<long>( B ), lFactor, lpLimbs );

	return lzi_packed_pack_in<B, 0>( lpDigits, lLength, static_cast<long>( B ), lFactor, lpLimbs );
}

template <uint64_t B>
static void lzi_packed_unpack_for( const uint64_t* lpLimbs, long lLimbs, long lFactor, long* lpDigits, long lLength )
{
	if ( lFactor == lzi_multiply_limb_digits_for( B ) )
		lzi_packed_unpack_in<B, lzi_multiply_limb_digits_for( B )>( lpLimbs, lLimbs, static_cast<long>( B ), lFactor, lpDigits, lLength );
	else if ( lFactor == lzi_limb_digits_for( B ) )
		lzi_packed_unpack_in<B, lzi_limb_digits_for( B )>( lpLimbs, lLimbs, static_cast<long>( B ), lFactor, lpDigits, lLength );
	else
		lzi_packed_unpack_in<B, 0>( lpLimbs, lLimbs, static_cast<long>( B ), lFactor, lpDigits, lLength );
}

static long lzi_packed_pack( const long* lpDigits, long lLength, long lBase, long lFactor,
							 uint64_t* lpLimbs )
{
	switch ( lBase )
	{
	case 2:
		if ( lFactor == kPackedDigitsBinaryMultiply )
			return lzi_packed_pack_in<2, kPackedDigitsBinaryMultiply>( lpDigits, lLength, lBase, lFactor, lpLimbs );
		if ( lFactor == kPackedDigitsBinaryDivide )
			return lzi_packed_pack_in<2, kPackedDigitsBinaryDivide>( lpDigits, lLength, lBase, lFactor, lpLimbs );
		return lzi_packed_pack_in<2, 0>( lpDigits, lLength, lBase, lFactor, lpLimbs );

	case 8:
		if ( lFactor == kPackedDigitsOctal )
			return lzi_packed_pack_in<8, kPackedDigitsOctal>( lpDigits, lLength, lBase, lFactor, lpLimbs );
		return lzi_packed_pack_in<8, 0>( lpDigits, lLength, lBase, lFactor, lpLimbs );

	case 10:
		if ( lFactor == kPackedDigitsDecimal )
			return lzi_packed_pack_in<10, kPackedDigitsDecimal>( lpDigits, lLength, lBase, lFactor, lpLimbs );
		return lzi_packed_pack_in<10, 0>( lpDigits, lLength, lBase, lFactor, lpLimbs );

	case 12:
		if ( lFactor == kPackedDigitsDuodecimal )
			return lzi_packed_pack_in<12, kPackedDigitsDuodecimal>( lpDigits, lLength, lBase, lFactor, lpLimbs );
		return lzi_packed_pack_in<12, 0>( lpDigits, lLength, lBase, lFactor, lpLimbs );

	case 16:
		if ( lFactor == kPackedDigitsHex )
			return lzi_packed_pack_in<16, kPackedDigitsHex>( lpDigits, lLength, lBase, lFactor, lpLimbs );
		return lzi_packed_pack_in<16, 0>( lpDigits, lLength, lBase, lFactor, lpLimbs );

	case 3:		return lzi_packed_pack_for<3>( lpDigits, lLength, lFactor, lpLimbs );
	case 4:		return lzi_packed_pack_for<4>( lpDigits, lLength, lFactor, lpLimbs );
	case 5:		return lzi_packed_pack_for<5>( lpDigits, lLength, lFactor, lpLimbs );
	case 6:		return lzi_packed_pack_for<6>( lpDigits, lLength, lFactor, lpLimbs );
	case 7:		return lzi_packed_pack_for<7>( lpDigits, lLength, lFactor, lpLimbs );
	case 9:		return lzi_packed_pack_for<9>( lpDigits, lLength, lFactor, lpLimbs );
	case 11:	return lzi_packed_pack_for<11>( lpDigits, lLength, lFactor, lpLimbs );
	case 13:	return lzi_packed_pack_for<13>( lpDigits, lLength, lFactor, lpLimbs );
	case 14:	return lzi_packed_pack_for<14>( lpDigits, lLength, lFactor, lpLimbs );
	case 15:	return lzi_packed_pack_for<15>( lpDigits, lLength, lFactor, lpLimbs );

	default:
		return lzi_packed_pack_in<0, 0>( lpDigits, lLength, lBase, lFactor, lpLimbs );
	}
}

static void lzi_packed_unpack( const uint64_t* lpLimbs, long lLimbs, long lBase, long lFactor,
							   long* lpDigits, long lLength )
{
	switch ( lBase )
	{
	case 2:
		if ( lFactor == kPackedDigitsBinaryMultiply )
			lzi_packed_unpack_in<2, kPackedDigitsBinaryMultiply>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		else if ( lFactor == kPackedDigitsBinaryDivide )
			lzi_packed_unpack_in<2, kPackedDigitsBinaryDivide>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		else
			lzi_packed_unpack_in<2, 0>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		break;

	case 8:
		if ( lFactor == kPackedDigitsOctal )
			lzi_packed_unpack_in<8, kPackedDigitsOctal>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		else
			lzi_packed_unpack_in<8, 0>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		break;

	case 10:
		if ( lFactor == kPackedDigitsDecimal )
			lzi_packed_unpack_in<10, kPackedDigitsDecimal>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		else
			lzi_packed_unpack_in<10, 0>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		break;

	case 12:
		if ( lFactor == kPackedDigitsDuodecimal )
			lzi_packed_unpack_in<12, kPackedDigitsDuodecimal>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		else
			lzi_packed_unpack_in<12, 0>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		break;

	case 16:
		if ( lFactor == kPackedDigitsHex )
			lzi_packed_unpack_in<16, kPackedDigitsHex>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		else
			lzi_packed_unpack_in<16, 0>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		break;

	case 3:		lzi_packed_unpack_for<3>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 4:		lzi_packed_unpack_for<4>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 5:		lzi_packed_unpack_for<5>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 6:		lzi_packed_unpack_for<6>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 7:		lzi_packed_unpack_for<7>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 9:		lzi_packed_unpack_for<9>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 11:	lzi_packed_unpack_for<11>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 13:	lzi_packed_unpack_for<13>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 14:	lzi_packed_unpack_for<14>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;
	case 15:	lzi_packed_unpack_for<15>( lpLimbs, lLimbs, lFactor, lpDigits, lLength );		break;

	default:
		lzi_packed_unpack_in<0, 0>( lpLimbs, lLimbs, lBase, lFactor, lpDigits, lLength );
		break;
	}
}

/*
	A divisor of one limb: short division, with nothing to estimate and
	nothing to normalise. Two fifths of the engine's divisions have a
	divisor this short once it is packed.
*/

static bool lzi_packed_divide_short( const uint64_t* lpDividend, long lDividendLimbs,
									 uint64_t lDivisor, uint64_t lLimbBase, uint64_t* lpQuotient,
									 uint64_t* plRemainder )
{
	uint64_t	l_remainder = 0;

	for ( long i = lDividendLimbs - 1; i >= 0; i-- )
	{
		uint64_t	l_current = (l_remainder * lLimbBase) + lpDividend[ i ];

		lpQuotient[ i ] = l_current / lDivisor;
		l_remainder = l_current % lDivisor;
	}

	*plRemainder = l_remainder;

	return (l_remainder != 0);
}

/*
	**Twice a remainder against the divisor, on limbs** -- lzi_twice_against
	for the packed path, which an odd base needs since 2026-10-12, when
	bases three to fifteen were given limbs: -1 below, 0 equal, 1 above.
	Algorithm D leaves the remainder scaled by the factor that normalised
	the divisor, and the divisor scaled by it too, so the two compare as
	they stand. lPlus is added to twice the remainder first, scaled as
	the remainder is: lzi_twice_against's lPlus, for a cut dividend.
*/

static long lzi_packed_twice_against( const uint64_t* lpRemainder, long lRemainderLimbs,
									  const uint64_t* lpDivisor, long lDivisorLimbs, uint64_t lLimbBase,
									  uint64_t lPlus )
{
	// Twice the remainder, carried, compared from the top down against a
	//	divisor no longer than it plus a limb...
	long		l_limbs = (lRemainderLimbs > lDivisorLimbs) ? lRemainderLimbs : lDivisorLimbs;
	uint64_t	l_carry = lPlus;
	long		l_answer = 0;

	std::vector<uint64_t>	twice( static_cast<size_t>( l_limbs ) + 1, 0 );

	for ( long i = 0; i < lRemainderLimbs; i++ )
	{
		uint64_t	l_value = (2 * lpRemainder[ i ]) + l_carry;

		twice[ static_cast<size_t>( i ) ] = l_value % lLimbBase;
		l_carry = l_value / lLimbBase;
	}

	twice[ static_cast<size_t>( lRemainderLimbs ) ] += l_carry;

	for ( long i = l_limbs; (i >= 0) && (l_answer == 0); i-- )
	{
		uint64_t	l_divisor = (i < lDivisorLimbs) ? lpDivisor[ i ] : 0;

		if ( twice[ static_cast<size_t>( i ) ] != l_divisor )
			l_answer = (twice[ static_cast<size_t>( i ) ] > l_divisor) ? 1 : -1;
	}

	return l_answer;
}

/*
	Algorithm D. Both are normalised in place by the same factor, which
	leaves the quotient alone and makes the divisor's leading limb at
	least half the limb base -- which is what bounds the estimate. The
	dividend needs one limb of headroom for it. The factor is answered
	in *plScale, the remainder left in u being scaled by it.
*/

static bool lzi_packed_divide_long( uint64_t* u, long lDividendLimbs, uint64_t* v, long lDivisorLimbs,
									uint64_t lLimbBase, uint64_t* lpQuotient, uint64_t* plScale )
{
	long	l_n = lDivisorLimbs;
	long	l_m = lDividendLimbs - l_n;

	uint64_t	l_scale = lLimbBase / (v[ l_n - 1 ] + 1);

	*plScale = l_scale;

	{
		uint64_t	l_carry = 0;

		for ( long i = 0; i < lDividendLimbs; i++ )
		{
			uint64_t	l_value = (u[ i ] * l_scale) + l_carry;

			u[ i ] = l_value % lLimbBase;
			l_carry = l_value / lLimbBase;
		}

		u[ lDividendLimbs ] = l_carry;

		l_carry = 0;

		for ( long i = 0; i < l_n; i++ )
		{
			uint64_t	l_value = (v[ i ] * l_scale) + l_carry;

			v[ i ] = l_value % lLimbBase;
			l_carry = l_value / lLimbBase;
		}
	}

	for ( long j = l_m; j >= 0; j-- )
	{
		// The estimate, from the two leading limbs over the one...
		uint64_t	l_numerator = (u[ j + l_n ] * lLimbBase) + u[ j + l_n - 1 ];
		uint64_t	l_digit = l_numerator / v[ l_n - 1 ];
		uint64_t	l_rest = l_numerator % v[ l_n - 1 ];

		// ...brought to at most one too high by the next limb of each...
		while ( (l_digit >= lLimbBase) ||
				((l_digit * v[ l_n - 2 ]) > ((l_rest * lLimbBase) + u[ j + l_n - 2 ])) )
		{
			l_digit--;
			l_rest += v[ l_n - 1 ];

			if ( l_rest >= lLimbBase )
				break;
		}

		// ...that many divisors taken away in one pass...
		uint64_t	l_carry = 0;
		int64_t		l_borrow = 0;

		for ( long i = 0; i < l_n; i++ )
		{
			uint64_t	l_product = (l_digit * v[ i ]) + l_carry;
			int64_t		l_value = static_cast<int64_t>( u[ i + j ] ) -
								  static_cast<int64_t>( l_product % lLimbBase ) - l_borrow;

			l_carry = l_product / lLimbBase;

			if ( l_value < 0 )
			{
				l_value += static_cast<int64_t>( lLimbBase );
				l_borrow = 1;
			}
			else
				l_borrow = 0;

			u[ i + j ] = static_cast<uint64_t>( l_value );
		}

		{
			int64_t		l_value = static_cast<int64_t>( u[ j + l_n ] ) -
								  static_cast<int64_t>( l_carry ) - l_borrow;

			if ( l_value < 0 )
			{
				l_value += static_cast<int64_t>( lLimbBase );
				l_borrow = 1;
			}
			else
				l_borrow = 0;

			u[ j + l_n ] = static_cast<uint64_t>( l_value );
		}

		// ...and one given back if that was one too many...
		if ( l_borrow )
		{
			uint64_t	l_add = 0;

			l_digit--;
			gPackedAddBacks++;

			for ( long i = 0; i < l_n; i++ )
			{
				uint64_t	l_value = u[ i + j ] + v[ i ] + l_add;

				l_add = (l_value >= lLimbBase) ? 1 : 0;
				u[ i + j ] = l_value - (l_add ? lLimbBase : 0);
			}

			u[ j + l_n ] = (u[ j + l_n ] + l_add) % lLimbBase;
		}

		lpQuotient[ j ] = l_digit;
	}

	// The remainder is what is left in u, scaled, and scaling cannot turn
	//	a remainder that is not nought into one that is...
	for ( long i = 0; i < l_n; i++ )
	{
		if ( u[ i ] != 0 )
			return true;
	}

	return false;
}

bool lzi_packed_divide( const long* lpDividend, long lDividendLength,
						const long* lpDivisor, long lDivisorLength, long lBase,
						long* lpQuotient, long lQuotientLength, bool* pbRemainder, long* plHalf,
						long lHalfPlus )
{
	long	l_factor = lzi_packed_factor( lBase );

	if ( l_factor == 0 )
		return false;

	uint64_t	l_limbBase = 1;

	for ( long i = 0; i < l_factor; i++ )
		l_limbBase *= static_cast<uint64_t>( lBase );

	long	l_dividendLimbs = ((lDividendLength + l_factor) - 1) / l_factor;
	long	l_divisorLimbs = ((lDivisorLength + l_factor) - 1) / l_factor;

	// One block for the lot: the dividend with a limb of headroom, the
	//	divisor, and the quotient, which is never longer than the dividend...
	uint64_t*	lp_block = static_cast<uint64_t*>(
					::calloc( (2 * l_dividendLimbs) + l_divisorLimbs + 2, sizeof(uint64_t) ) );

	if ( lp_block == nullptr )
		return false;

	uint64_t*	lp_dividend = lp_block;
	uint64_t*	lp_divisor = lp_dividend + l_dividendLimbs + 1;
	uint64_t*	lp_quotient = lp_divisor + l_divisorLimbs;

	l_dividendLimbs = lzi_packed_pack( lpDividend, lDividendLength, lBase, l_factor, lp_dividend );
	l_divisorLimbs = lzi_packed_pack( lpDivisor, lDivisorLength, lBase, l_factor, lp_divisor );

	// lzi_divide_op never asks this, but a direct caller can...
	if ( l_divisorLimbs == 0 )
	{
		::free( lp_block );

		return false;
	}

	long	l_quotientLimbs = 0;
	bool	b_remainder = false;
	long	l_half = 0;

	gPackedDivisions++;

	// A dividend under the divisor divides to nought, and is its own remainder...
	if ( l_dividendLimbs < l_divisorLimbs )
	{
		b_remainder = (l_dividendLimbs > 0);

		if ( b_remainder && (plHalf != nullptr) )
			l_half = lzi_packed_twice_against( lp_dividend, l_dividendLimbs, lp_divisor, l_divisorLimbs, l_limbBase,
											   static_cast<uint64_t>( lHalfPlus ) );
	}
	else if ( l_divisorLimbs == 1 )
	{
		uint64_t	l_remainder = 0;

		l_quotientLimbs = l_dividendLimbs;
		b_remainder = lzi_packed_divide_short( lp_dividend, l_dividendLimbs, lp_divisor[ 0 ],
											   l_limbBase, lp_quotient, &l_remainder );

		// A remainder is under the divisor, itself under 2^32, so twice it fits...
		if ( b_remainder && (plHalf != nullptr) )
		{
			uint64_t	l_twice = (2 * l_remainder) + static_cast<uint64_t>( lHalfPlus );

			l_half = (l_twice > lp_divisor[ 0 ]) ? 1 : ((l_twice == lp_divisor[ 0 ]) ? 0 : -1);
		}
	}
	else
	{
		uint64_t	l_scale = 1;

		l_quotientLimbs = (l_dividendLimbs - l_divisorLimbs) + 1;
		b_remainder = lzi_packed_divide_long( lp_dividend, l_dividendLimbs, lp_divisor, l_divisorLimbs,
											  l_limbBase, lp_quotient, &l_scale );

		// The remainder in the dividend's low limbs, scaled as the divisor
		//	is, and what is added to it scaled the same way...
		if ( b_remainder && (plHalf != nullptr) )
			l_half = lzi_packed_twice_against( lp_dividend, l_divisorLimbs, lp_divisor, l_divisorLimbs, l_limbBase,
											   static_cast<uint64_t>( lHalfPlus ) * l_scale );
	}

	if ( plHalf != nullptr )
		*plHalf = l_half;

	lzi_packed_unpack( lp_quotient, l_quotientLimbs, lBase, l_factor, lpQuotient, lQuotientLength );

	::free( lp_block );

	*pbRemainder = b_remainder;

	return true;
}


/*
	Twice a remainder against the divisor, both whole numbers held as
	digits in one base, least significant first: -1 below, 0 equal, 1
	above. That is where the remainder stands against half the
	divisor, which an odd base needs to round a quotient; see the
	sticky digit in lzi_divide_op. lPlus, nought or one, is added to
	twice the remainder first: one when a cut dividend dropped more than
	half a unit of its last place -- see the cut in lzi_divide_work.
*/

static long lzi_twice_against( const long* lpRemainder, long lRemainderLength,
							   const long* lpDivisor, long lDivisorLength, long lBase, long lPlus )
{
	long*	lp_twice = static_cast<long*>(::calloc( lRemainderLength + 1, sizeof(long) ));

	if ( lp_twice == nullptr )
		return 1;

	long	l_carry = lPlus;

	for ( long i = 0; i < lRemainderLength; i++ )
	{
		long	l_digit = (2 * lpRemainder[ i ]) + l_carry;

		lp_twice[ i ] = l_digit % lBase;
		l_carry = l_digit / lBase;
	}

	lp_twice[ lRemainderLength ] = l_carry;

	long	l_twiceLength = lRemainderLength + 1;

	while ( (l_twiceLength > 0) && (lp_twice[ l_twiceLength - 1 ] == 0) )
		l_twiceLength--;

	while ( (lDivisorLength > 0) && (lpDivisor[ lDivisorLength - 1 ] == 0) )
		lDivisorLength--;

	long	l_answer = 0;

	if ( l_twiceLength != lDivisorLength )
		l_answer = (l_twiceLength > lDivisorLength) ? 1 : -1;
	else
	{
		for ( long i = l_twiceLength - 1; (i >= 0) && (l_answer == 0); i-- )
		{
			if ( lp_twice[ i ] != lpDivisor[ i ] )
				l_answer = (lp_twice[ i ] > lpDivisor[ i ]) ? 1 : -1;
		}
	}

	::free( lp_twice );

	return l_answer;
}


/*
	Where a quotient's leading digit lies, as a power of the base: 0 for a
	quotient from one up to the base, -1 for one from a base-th up to one,
	and so on down. The dividend's leading position less the divisor's,
	less one more when the dividend's digits from there compare below the
	divisor's. The divisor is whole and trimmed by the time this is asked;
	the dividend carries the radix the two had between them. A nought
	dividend answers 0.
*/

static long lzi_quotient_lead( const CLZNumber* pDividend, const CLZNumber* pDivisor )
{
	long	l_top = pDividend->m_lDigitsLength - 1;

	while ( (l_top >= 0) && (pDividend->m_lpDigits[ l_top ] == 0) )
		l_top--;

	if ( l_top < 0 )
		return 0;

	long	i = l_top;
	long	j = pDivisor->m_lDigitsLength - 1;
	long	l_compare = 0;

	for ( ; (i >= 0) && (j >= 0); i--, j-- )
	{
		if ( pDividend->m_lpDigits[ i ] != pDivisor->m_lpDigits[ j ] )
		{
			l_compare = (pDividend->m_lpDigits[ i ] < pDivisor->m_lpDigits[ j ]) ? -1 : 1;
			break;
		}
	}

	// The dividend's digits ran out first: smaller if any of the divisor's are left...
	for ( ; (l_compare == 0) && (j >= 0); j-- )
	{
		if ( pDivisor->m_lpDigits[ j ] != 0 )
			l_compare = -1;
	}

	return (l_top - pDividend->m_lRadix) - (pDivisor->m_lDigitsLength - 1) - ((l_compare < 0) ? 1 : 0);
}


/*
	**Division, and division to a whole number.** lzi_divide_op below is
	the first; lzi_divide_whole the second, for the integer quotient and
	the remainders.

	Those used to divide to the precision, rounding there, and take a
	whole number from the rounded quotient -- which a quotient just
	short of a whole number, or just past a half, rounds across:
	1 / 0.333...334 to fifty places is 2.999...994, rounded to 3, and
	its integer quotient came back 3. The raw quotient already carries
	what a whole number needs, its digits down to the first fractional
	one and the sticky digit after them saying whether anything was
	left and, in an odd base, where it stands against half. So bWhole
	rounds that straight to a whole number in rmWhole, once, and works
	to only a place past the point to get there.

	It never takes the reciprocal path: that rounds its quotient on the
	way, which is exactly the fault being removed. Past a million places
	of whole number, long division is the price of being right.

	**bAnswerWidth rounds the raw quotient once at the width a function's
	answer is rounded to** -- the precision, and for a quotient under a
	tenth its noughts as well; see lzi_round_quotient_answer -- in place
	of the quotient's own, a place or two wider. A function whose answer
	is a single exact division takes it, lz_reciprocal since 2026-10-03:
	dividing wider and rounding again rounded twice, and an answer within
	the extra places' last unit of a half could round the wrong way.
*/

static long lzi_divide_work( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber,
							 long precision, bool bWhole, RoundingMode rmWhole, bool bAnswerWidth = false )
{
	if ( leftOperand->m_lBase != rightOperand->m_lBase )
		return kBaseMismatch;

	// The target's rounding settings are its own -- see CLZTargetRounding
	//	and the rounding below...
	CLZTargetRounding	targetRounding( (targetNumber != nullptr) ? targetNumber : leftOperand );

	long	l_dividendOffset = 0;

	CLZNumber*	p_dividend = nullptr;
	CLZNumber*	p_divisor = nullptr;
	CLZNumber*	p_quotient = nullptr;

	long		l_precision = bWhole ? 1 : precision;
	long		l_adjustPrecision = 0;
	long		l_target = 0;			// the places the answer is rounded to; see below
	long		l_lead = 0;				// where the quotient's leading digit lies; see lzi_quotient_lead

	bool		b_middleRepeats = false;

	// A dividend cut to the places worked, whether what it dropped was not
	//	all nought, and whether it was over half a unit of the last place
	//	kept -- see the cut below...
	bool		b_cut = false;
	bool		b_cutLeft = false;
	long		l_cutHalf = 0;

	// The width the rule gives a quotient whose leading digit lies at lLead...
	auto	width_for = [ precision, bAnswerWidth ]( long lLead ) -> long
	{
		if ( bAnswerWidth )
			return precision + ((lLead < -1) ? (-lLead - 1) : 0);

		return precision + 1 + ((lLead < 0) ? -lLead : 0);
	};

	long		lResult = kNoError;

	p_quotient = new CLZNumber;

	// Built digit by digit in the operands' base, so it has to carry
	//	that base or the target inherits base ten from it...
	p_quotient->SetBase( leftOperand->m_lBase, true );

	/*
		Nothing can be divided by zero, and the loop below would never
		finish trying: it subtracts the divisor until the dividend goes
		negative, which subtracting zero never does.

		IsZero() reports a flag rather than reading the digits, so ask
		for that flag to be brought up to date before trusting it - a
		divisor that arrived here as the result of some earlier
		operation may have rounded away to nothing since it was last
		examined.
	*/

	rightOperand->CheckForZero();

	if ( rightOperand->IsZero() )
	{
		leftOperand->CheckForZero();

		CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		// Asked before the target is written: in place the target is the
		//	dividend, and until 2026-10-15 x / 0 in place zeroed it first and
		//	reported nought over nought -- NaN and kUndefinedResult...
		bool	b_noughtOverNought = leftOperand->IsZero();

		if ( !(p_target->IsConstant()) )
		{
			p_target->Assign( static_cast<long>(0) );
			p_target->SetFlag( kfDiv0 );

			// Zero divided by zero is not infinite, merely undefined...
			if ( b_noughtOverNought )
				p_target->SetFlag( kfQNaN );
			else
				p_target->SetFlag( kfInfinity );
		}

		lResult = b_noughtOverNought ? kUndefinedResult : kDivideByZero;

		goto DivideExit;
	}

	// If we can divide, proceed...
	if ( !(rightOperand->IsInfinity()) )
	{
		// Create working objects...
		p_dividend = new CLZNumber;
		p_divisor = new CLZNumber;

		// Set up divisor and eliminate the fractional part (if any)...
		p_divisor->Assign( rightOperand );
		p_divisor->m_lRadix = 0;

		/*
			And drop its leading zeros, because from here on the
			divisor's length is taken as its size.

			A fraction carries a zero in front of the radix mark -- .1
			occupies two digit slots, not one -- and setting the radix
			to nought above turns that into a leading zero rather than
			removing it. The offset below is then computed as the
			dividend's length less the divisor's, comes out one too
			small, and the first quotient digit has to absorb a whole
			factor of the base to compensate.

			It cannot: a quotient digit is incremented one at a time by
			the loop further down and nothing normalises it, so it
			simply holds a value the base has no digit for. That is why
			1 / .1 answered "A" -- a single digit holding ten -- and
			why 1.99 / .1 answered .9, the nineteen having overflowed
			the only slot it had. Anything whose divisor was below one
			and whose quotient reached ten was wrong.
		*/

		p_divisor->TrimLeading();

		// Set up dividend and shift by the change in the divisor's radix...
		p_dividend->Assign( leftOperand );

		if ( p_dividend->m_lRadix >= rightOperand->m_lRadix )
		{
			p_dividend->m_lRadix -= rightOperand->m_lRadix;
		}
		else
		{
			p_dividend->m_lRadix = 0;
			p_dividend->ShiftLeft( rightOperand->m_lRadix - leftOperand->m_lRadix );
		}

		// Set quotient's exponent to match dividend's...
		p_quotient->m_lRadix = p_dividend->m_lRadix;

		// Where the quotient's leading digit lies, which sets its width below...
		l_lead = lzi_quotient_lead( p_dividend, p_divisor );

		/*
			**A dividend wider than the answer is cut to the places the
			quotient is worked to**, since 2026-10-16 -- performance item 7
			of the review of 2026-10-15. The quotient was worked to every
			place of the dividend and all but the answer's dropped: 125
			times the work for a 16,000-place dividend at 36 places.

			Cutting is exact. With the divisor whole, D, and the dividend
			cut to the G places worked, A = A' + E with E under a unit of
			place G, and A'b^G = QD + R with R under D: R + Eb^G is then
			under D, so the quotient's digits to place G are Q, the very
			digits the whole dividend gives. What was left is not nought
			when R or E is not, which is the sticky digit in an even base.
			An odd base wants where it stands against half the divisor,
			2R + 2Eb^G against D, and with R whole that turns on one more
			fact about E: whether Eb^G is over a half, which no finite
			string of digits in an odd base equals. So 2R plus that, one or
			nought, against D settles it, the two being equal meaning over,
			since E is not nought; see lzi_twice_against.

			**It keeps kRecordedTailDigits places past the answer**, the
			digits a rounding records for a function deciding whether its
			answer was in doubt, so the quotient, its rounding and the tail
			it records are all what the whole dividend gave. Division to a
			whole number rounds at the units, and so keeps its places past
			those. The leading digit was found above, from the whole
			dividend; the cut cannot move it, the digits being the same.

			**Two limits on how far.** It leaves the dividend a digit longer
			than the divisor, so nothing below pads it: a pad is a place
			more, and the quotient would be worked a place past the cut.
			And it leaves a quotient that packed enough digits to pack, or
			the loop would take it and be slower for doing less -- see
			lzi_packed_divide_digits.
		*/

		if ( gDivideCuts )
		{
			long	l_keep = (bWhole ? 1 : width_for( l_lead )) + CLZNumber::kRecordedTailDigits;
			long	l_cut = p_dividend->m_lRadix - l_keep;
			long	l_most = (p_dividend->m_lDigitsLength - p_divisor->m_lDigitsLength) - 1;

			// ...and a quotient that packed enough digits to pack...
			long	l_packs = lzi_packed_divide_digits( leftOperand->m_lBase, p_divisor->m_lDigitsLength );
			long	l_packedMost = (p_dividend->m_lDigitsLength - p_divisor->m_lDigitsLength + 1) - l_packs;

			if ( l_most > l_packedMost )
				l_most = l_packedMost;

			if ( l_cut > l_most )
				l_cut = l_most;

			if ( l_cut > 0 )
			{
				long	l_middle = leftOperand->m_lBase / 2;
				bool	b_decided = false;

				// From the top: the first digit off the middle says which side
				//	of half it lies, and any digit that is not nought that it
				//	was not all nought...
				for ( long i = l_cut - 1; (i >= 0) && !(b_cutLeft && b_decided); i-- )
				{
					long	l_digit = p_dividend->m_lpDigits[ i ];

					if ( l_digit != 0 )
						b_cutLeft = true;

					if ( !b_decided && (l_digit != l_middle) )
					{
						l_cutHalf = (l_digit > l_middle) ? 1 : 0;
						b_decided = true;
					}
				}

				p_dividend->ShiftRight( l_cut );
				p_quotient->m_lRadix = p_dividend->m_lRadix;

				b_cut = true;
			}
		}

		/*
			Make sure the dividend is 'larger' than the divisor.

			**All the noughts at once.** This padded one at a time, and each
			PadRight moves every digit up, so a dividend k digits short of
			the divisor cost k moves of a growing array -- quadratic in the
			divisor's length. 1 over a 16,900-digit integer spent 99% of its
			79 ms here, before a single quotient digit. One PadRight of k
			leaves the same digits, length and radix; only the array's spare
			room differs, which the headroom check below sees to.
		*/

		bool	b_padded = false;

		if ( p_dividend->m_lDigitsLength < p_divisor->m_lDigitsLength )
		{
			long	l_short = p_divisor->m_lDigitsLength - p_dividend->m_lDigitsLength;

			p_quotient->PadRight( l_short );
			p_dividend->PadRight( l_short );

			b_padded = true;
		}

		// Set the initial dividend offset...
		l_dividendOffset = p_dividend->m_lDigitsLength - p_divisor->m_lDigitsLength;

		// Shift one more digit, if required...
		if ( p_dividend->CompareRaw( p_divisor ) == kLessThan )
		{
			p_quotient->PadRight( 1 );
			p_dividend->PadRight( 1 );

			b_padded = true;
		}

		/*
			**A quotient's width is set by its value, since 2026-09-30**:
			the precision plus one place from one up, and the precision
			plus two significant digits below one -- one place more for
			each nought after the point. Nothing about how the dividend
			was written enters into it.

			Until then the width was counted from the bookkeeping: the
			radix the dividend brought, the pads, and (until 2026-09-29)
			the pads again. So at 36 places a quotient below one carried
			anywhere from 32 significant digits to 43 -- 0.001/300 had 34,
			1/7 38, 0.12345678/7 43 -- and one above carried the
			dividend's places on top of the precision, so a chain of
			divisions grew by the precision at every step: 72, 108, 144
			places.

			The digits are still computed to at least the dividend's own
			places, which the loop below needs, and rounded once, at the
			width, after it -- RoundAt reads every digit it discards, so a
			longer tail rounds exactly as a guard digit and a sticky one
			would. **Division to a whole number keeps the count it always
			had**: it rounds the raw quotient to a whole number, and the
			width there is only how far past the point it looks.
		*/

		if ( bWhole )
		{
			l_adjustPrecision = p_quotient->m_lRadix + (b_padded ? 1 : 0);

			if ( l_adjustPrecision > 0 )
				l_precision += (l_adjustPrecision - 1);

			// Do one more than we need...
			l_precision++;

			l_target = l_precision;
		}
		else
		{
			l_target = width_for( l_lead );

			l_precision = (l_target > p_quotient->m_lRadix) ? l_target : p_quotient->m_lRadix;
		}

		// A cut dividend is worked to its own places and not one past them,
		//	or the digits are no longer the whole dividend's -- see the cut...
		if ( b_cut )
		{
			l_precision = p_quotient->m_lRadix - 1;

			if ( bWhole )
				l_target = l_precision;
		}

		// Set quotient precision...
		p_quotient->SetPrecision( l_precision );

		/*
			Increase the dividend digits array size, if required.

			The room needed is the divisor's length and not a fixed
			eight. The subtraction loop below walks from the offset to
			one past the divisor's last digit -- the extra element
			carries the final borrow -- so a divisor longer than the
			headroom walks off the end of the allocation.

			Eight was enough for every divisor anyone had divided by
			until the circular functions reduced a huge argument modulo
			2*pi, which is a divisor of forty-nine digits: cos(1e40)
			read eight bytes past a 560-byte array. Nothing failed, the
			slack the array is calloc'd with having absorbed it, which
			is exactly the class of fault the notes say only ASan will
			show.
		*/

		while ( (p_dividend->m_lDigitsArrayLength - l_dividendOffset) <
					(p_divisor->m_lDigitsLength + 8) )
		{
			p_dividend->SetDigitsArraySize( p_dividend->m_lDigitsArrayLength + 16 );
		}

		/*
			Past the floor, stop dividing and multiply by a reciprocal
			instead -- which, with the floor at a million places, means
			in practice only when the tests have lowered it. See
			kReciprocalFloor for why.

			Everything above has already settled the precision, the
			radix and the sign, so only the digits are wanted here.

			The reciprocal is taken of the magnitude and multiplied by
			the magnitude: the sign is applied below, as it is for the
			loop, and giving either operand a sign here would apply it
			twice.
		*/

		/*
			When the floor was four hundred, a divisor-length condition
			was written here and taken out again once the two paths
			were timed by one harness -- the comparison that motivated
			it had mixed two. **It is back, measured by one harness this
			time**: since 2026-09-26 the reciprocal wins for a long
			divisor from a width that depends on the base, and loses
			for a short one at every width -- see
			lzi_reciprocal_by_shape, which since 2026-10-16 counts the
			quotient's significant digits and the divisor's. The floor on
			the quotient alone stays, out of reach, for the tests to force
			the path.
		*/

		// Past the floor, or a quotient and divisor of a shape the
		//	reciprocal wins -- see lzi_reciprocal_by_shape. The quotient's
		//	significant digits are its places from its leading digit...
		int64_t	l_quotientDigits = static_cast<int64_t>( l_precision ) + l_lead + 1;
		bool	b_byShape = gReciprocalByShape &&
							lzi_reciprocal_by_shape( leftOperand->m_lBase, l_quotientDigits, p_divisor->m_lDigitsLength );

		// Not for a cut dividend: it multiplies the operands as they came...
		if ( ((l_precision >= gReciprocalFloor) || b_byShape) && !gSeekingReciprocal && !bWhole && !b_cut )
		{
			CLZNumber	l_reciprocal;
			CLZNumber	l_dividendMagnitude;
			CLZNumber	l_divisorMagnitude;
			CLZNumber	l_quotient;

			lz_abs( leftOperand, &l_dividendMagnitude );
			lz_abs( rightOperand, &l_divisorMagnitude );

			/*
				The reciprocal has to be carried past the answer by as
				much as the dividend will magnify it. An error of
				base^-p in 1/d becomes base^(m-p) in a*(1/d) when the
				dividend is about base^m, so a thirty-digit dividend
				eats thirty places of the reciprocal's accuracy before
				the multiplication is done.

				Without this the quotient came back with its last
				thirty-odd digits noise -- correct to four hundred
				places when four hundred and two had been asked for,
				and only for an integer dividend over a fractional
				divisor, which is the shape that magnifies most.
			*/

			/*
				**Sized by the quotient's significant digits, not the
				dividend's size alone.** That was right while the divisor
				sat near one, which is all the reciprocal was tried on. A
				divisor of 24,000 digits -- lzi_e_in_base's, for e^50000 --
				then asked for twice the digits the quotient has, four
				times the work, and a long divisor went slower through the
				reciprocal than through long division. The quotient's
				leading digit sits at the dividend's exponent less the
				divisor's, so its digits are the places asked for plus
				that, and one for the estimate. Its places also count a
				small quotient's leading noughts, which this takes back
				off.
			*/

			long	l_quotientExponent = lzi_leading_exponent( &l_dividendMagnitude ) -
										 lzi_leading_exponent( &l_divisorMagnitude ) + 2;

			/*
				**Exact, as long division is.** Until 2026-09-26 this
				rounded its quotient half-even at the width, and the
				target's rounding then had nothing left to discard, so
				past the floor a division ignored its target's mode: 2/3
				in floor mode gave ...667. So it now produces what long
				division does -- the quotient's digits exactly, one place
				past its precision, and a sticky digit after them -- and
				the rounding below treats the two the same way.

				The reciprocal gives a quotient good to a unit of that
				place or better. Truncated there, it is corrected by the
				exact remainder |a| - q|b| until 0 <= r < |b| units, and
				the sticky digit is read off r: whether anything was left
				and, in an odd base, where it stands against half. One
				more multiply, at full width, buys every mode.
			*/

			long	l_places = l_precision + 1;
			long	l_significant = l_places + l_quotientExponent;

			if ( l_significant < 1 )
				l_significant = 1;

			lzi_reciprocal_newton( &l_divisorMagnitude, l_significant, &l_reciprocal );

			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = l_places + kReciprocalGuardPlaces;

			lz_multiply( &l_dividendMagnitude, &l_reciprocal, &l_quotient );

			l_quotient.Round( RM_round_down, l_places );

			// The product at its full width, so the remainder is exact...
			glWorkingPrecision = l_places + l_divisorMagnitude.m_lRadix + 1;

			CLZNumber	l_product;
			CLZNumber	l_remainder;
			CLZNumber	l_unit;
			CLZNumber	l_step;

			lz_multiply( &l_quotient, &l_divisorMagnitude, &l_product );
			lz_subtract( &l_dividendMagnitude, &l_product, &l_remainder );

			// A unit of the last place, and the divisor's worth of it...
			l_unit.SetBase( leftOperand->m_lBase, false );
			l_unit.Assign( 1L );
			lzi_scale_by_base_power( &l_unit, l_places );

			l_step.Assign( &l_divisorMagnitude );
			lzi_scale_by_base_power( &l_step, l_places );

			l_remainder.CheckForZero();

			while ( l_remainder.IsNegative() )
			{
				lz_subtract( &l_quotient, &l_unit );
				lz_add( &l_remainder, &l_step );
				l_remainder.CheckForZero();
			}

			while ( lz_compare( &l_remainder, &l_step ) != kLessThan )
			{
				lz_add( &l_quotient, &l_unit );
				lz_subtract( &l_remainder, &l_step );
			}

			l_remainder.CheckForZero();

			long	l_sticky = l_remainder.IsZero() ? 0 : 1;

			if ( (l_sticky != 0) && ((leftOperand->m_lBase % 2) != 0) )
			{
				CLZNumber	l_twice;

				lz_add( &l_remainder, &l_remainder, &l_twice );

				long	l_against = lz_compare( &l_twice, &l_step );

				if ( l_against == kGreaterThan )
					l_sticky = leftOperand->m_lBase - 1;
				else if ( l_against == kEqualTo )
				{
					l_sticky = leftOperand->m_lBase / 2;
					b_middleRepeats = true;
				}
			}

			glWorkingPrecision = l_saved;

			// Laid out as long division lays its quotient out: the digits
			//	to l_places, and the sticky digit a place beyond...
			long	l_base = p_quotient->m_lBase;

			p_quotient->Assign( &l_quotient );
			p_quotient->m_lBase = l_base;

			if ( p_quotient->m_lRadix < l_places )
				p_quotient->PadRight( l_places - p_quotient->m_lRadix );

			p_quotient->PadRight( 1 );
			p_quotient->m_lpDigits[ 0 ] = l_sticky;
			p_quotient->SetPrecision( l_precision );
		}
		else
		{
			/*
				One digit per pass, and nothing shifted to make room for it.

				The loop used to make room for each digit by moving every
				digit of the quotient up one place, and every digit of the
				dividend too once the dividend's own digits were used up
				and noughts had to be brought down. Both grow as the
				division proceeds, so that was quadratic: a third of the
				time at 36 places and more than three quarters of it at
				390, once the digits themselves had stopped being found
				by trial.

				**How many digits, and how many of those need a nought
				brought down, is known before the first one.** The
				quotient's radix only moves on a digit that pads, and the
				loop runs until it passes the precision, so the padding is
				the precision less that radix, plus one, and the digits are
				that plus the dividend's starting offset, plus one. Checked
				against the shifting loop over 88,082 divisions -- the
				corpus, the suite and the sweep -- and it never disagreed.

				So the dividend is padded once with every nought it would
				have been given, and the window walks down through them;
				the quotient is widened once, and each digit is written
				where the shifting would have left it. The first lands on
				the nought already at the bottom, which the rest of the
				quotient is carried up with.
			*/

			// A cut dividend has a place more than the precision, and no padding...
			long	l_radix = p_quotient->m_lRadix;
			long	l_padding = (l_radix <= (l_precision + 1)) ? (l_precision - l_radix + 1) : 0;
			long	l_digits = (l_radix <= (l_precision + 1)) ? (1 + l_dividendOffset + l_padding) : 0;

			if ( l_digits > 0 )
			{
				p_dividend->PadRight( l_padding );
				l_dividendOffset += l_padding;

				/*
					The headroom check again, for the window's new
					position. **It cannot fire today**: ShiftLeft, under
					PadRight, always leaves eight elements spare past the
					digits, and the window's top never passes them. A
					mutation removing it passed the suite and the sweep
					under ASan. It is kept because that guarantee lives
					in another file, and a change there would otherwise
					walk this window off the end of its array silently.
				*/

				while ( (p_dividend->m_lDigitsArrayLength - l_dividendOffset) <
							(p_divisor->m_lDigitsLength + 8) )
				{
					p_dividend->SetDigitsArraySize( p_dividend->m_lDigitsArrayLength + 16 );
				}

				// One place more than the digits, for the sticky digit below...
				p_quotient->ShiftLeft( l_digits );
				p_quotient->m_lRadix += l_padding + 1;

				/*
					**Packed, past the floor, and it is all of them in
					practice**; see kPackedDivideFloor. The window the
					loop walks down is the whole dividend by the end of
					it, so the packed path divides the whole dividend at
					once and writes the same digits to the same places.
					It answers the remainder rather than leaving it in the
					dividend.
				*/

				bool	b_remainder = false;
				long	l_packedHalf = 0;
				bool	b_oddBase = ((leftOperand->m_lBase % 2) != 0);

				gQuotientDigitsWorked += l_digits;
				bool	b_packed = lzi_packed_divide_wanted( leftOperand->m_lBase, l_digits,
															 p_divisor->m_lDigitsLength ) &&
								   lzi_packed_divide( p_dividend->m_lpDigits, p_dividend->m_lDigitsLength,
													  p_divisor->m_lpDigits, p_divisor->m_lDigitsLength,
													  leftOperand->m_lBase,
													  p_quotient->m_lpDigits + 1, l_digits, &b_remainder,
													  b_oddBase ? &l_packedHalf : nullptr, l_cutHalf );

				if ( !b_packed )
				{
					for ( long l_digit = 1; l_digit <= l_digits; l_digit++ )
					{
						if ( l_digit > 1 )
							l_dividendOffset--;

						// This digit, leaving the remainder in the window...
						p_quotient->m_lpDigits[ l_digits - l_digit + 1 ] +=
							lzi_quotient_digit( p_dividend->m_lpDigits + l_dividendOffset,
												p_divisor->m_lpDigits,
												p_divisor->m_lDigitsLength,
												leftOperand->m_lBase );
					}

					for ( long i = 0; (i < p_dividend->m_lDigitsLength) && !b_remainder; i++ )
					{
						if ( p_dividend->m_lpDigits[ i ] != 0 )
							b_remainder = true;
					}
				}

				/*
					**The sticky digit: whether anything was left over.**

					The quotient carries one digit past its precision and is
					rounded from it, and that digit is truncated. So a
					guard digit of exactly half could not tell a quotient
					that stops there from one that goes on, and rounding
					half to even took some of the second kind down --
					"double-rounds and lands one place low now and again".
					Measured before this: 30 of 400 quotients at 400 to 800
					places wrong in their last place at full width.

					One more digit, 1 if the remainder is not nought and 0
					if it is, gives Round everything it needs in every mode:
					each that cares about the discarded tail -- half-even,
					half-down, half-odd, up, ceiling, floor -- asks only
					whether any of it is non-zero. It is discarded in the
					rounding, so the answer's width does not change, and a
					quotient that divides exactly gets a nought there and
					comes out as it always did.

					A single 1 is the right value in every base, including
					two, where bumping the guard digit instead would have
					nowhere to go from 1.
				*/

				// What the division left, and what a cut dropped...
				bool	b_dividedLeft = b_remainder;

				if ( b_cutLeft )
					b_remainder = true;

				p_quotient->m_lpDigits[ 0 ] = b_remainder ? 1 : 0;

				/*
					**An odd base needs more than whether anything was
					left: it needs where the remainder stands against
					half.** Its middle digit is under half a place, so a
					guard of 3 in base seven with anything at all after it
					used to round up, and 1/2 there is .3333... forever --
					exactly half a place past every digit, which no finite
					string of digits can say.

					So the sticky digit carries it: the largest digit when
					the remainder is over half the divisor, a 1 when it is
					under, which is below the middle digit or, in base
					three, equal to it and so still short of half; and the
					middle digit itself when it is exactly half, with the
					rounding told that the middle digit goes on. The digit
					loop leaves the remainder in the dividend, beside the
					divisor, both whole numbers; the packed path, which
					odd bases have taken since 2026-10-12, answers where it
					stands itself -- see lzi_packed_twice_against.
				*/

				if ( b_remainder && b_oddBase )
				{
					long	l_half = 0;

					// Nothing left by the division, only by the cut: twice nought
					//	and what the cut adds, against the divisor...
					if ( !b_dividedLeft )
						l_half = ((l_cutHalf != 0) && (p_divisor->m_lDigitsLength == 1) &&
								  (p_divisor->m_lpDigits[ 0 ] == 1)) ? 0 : -1;
					else
						l_half = b_packed ? l_packedHalf :
								 lzi_twice_against( p_dividend->m_lpDigits, p_dividend->m_lDigitsLength,
													p_divisor->m_lpDigits, p_divisor->m_lDigitsLength,
													leftOperand->m_lBase, l_cutHalf );

					// Past a cut, level is over: what it dropped was not nought...
					if ( b_cutLeft && (l_half == 0) )
						l_half = 1;

					if ( l_half > 0 )
						p_quotient->m_lpDigits[ 0 ] = leftOperand->m_lBase - 1;
					else if ( l_half == 0 )
					{
						p_quotient->m_lpDigits[ 0 ] = leftOperand->m_lBase / 2;
						b_middleRepeats = true;
					}
				}
			}
		}

		/*
			**Rounded at the width the value asks for, from everything
			computed**, when the dividend's places took the digits past
			it. Both paths lay the quotient out alike -- the digits, a
			guard place, a sticky place -- so rounding at the width reads
			everything after it.
		*/

		if ( l_target < l_precision )
			p_quotient->SetPrecision( l_target );

		// Determine result sign...
		if ( leftOperand->m_lSign == rightOperand->m_lSign )
			p_quotient->m_lSign = kPositive;
		else
			p_quotient->m_lSign = kNegative;

		// Normalize the exponent, as required...
		if ( p_quotient->m_lRadix > p_quotient->m_lDigitsLength )
		{
			p_quotient->PadLeft( p_quotient->m_lRadix - p_quotient->m_lDigitsLength );
		}
	}

	/*
		Assign the result, then put the target's precision back to what
		it came in with.

		The value is left exactly as it was computed. The quotient is
		deliberately wider than the caller's precision: l_adjustPrecision
		above counts the places a small quotient loses to its own
		leading zeros, and adding them back is what lets 1 / a 35-digit
		divisor carry 38 significant digits rather than the two that 36
		places would leave it. Rounding those away is not an option --
		it was tried while the noughts were still counted twice, and it
		took that division from 57 correct digits to 8.

		What must not persist is the wider figure in m_lPrecision.
		Assign copies precision along with the value, so a target came
		back wider than it went in, and the next division into it
		started from the wider figure and widened again. Dividing a
		hundred and eighty times into one target -- exactly what the
		logarithm's series does -- took it from 36 places to 1116, and
		the cost of the division with it. It also meant an answer's
		width depended on how many times its target had been used
		before: 1/3 came back with 38 places on the first call and 46
		on the fifth.

		So the value keeps its places and the field forgets them. The
		two disagreeing is deliberate: the places beyond the working
		precision are real digits of a small number, and the precision
		is what the next operation should work to.
	*/

	{
		CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
		{
			/*
				**In the target's own mode**, which Assign has just
				replaced with the quotient's half-even. The sticky digit
				above is right for every mode, so nothing else is needed:
				until 2026-09-26 this rounded half to even whatever the
				target asked for.
			*/

			p_target->Assign( p_quotient );
			targetRounding.RestoreMode();
			p_target->TrimAll();

			// Straight to a whole number, from the raw quotient -- see
			//	lzi_divide_work above...
			if ( bWhole && b_middleRepeats )
				p_target->RoundMiddleRepeating( rmWhole, 0L );
			else if ( bWhole )
				p_target->Round( rmWhole, 0L );
			else if ( b_middleRepeats )
				p_target->RoundMiddleRepeating();
			else
				p_target->Round();

			/*
				**A carry to a power of the base moves the leading digit up
				a place**, and the width set from where it lay is then a
				place wider than the rule gives the power: (10^42 - 1) /
				10^42 at 36 places rounds to 1, which takes 37 places, and
				until 2026-10-10 was written to 38 -- in every base, and
				0.0999...9 to 0.1 alike. Rounding at the narrower width
				carries the same way, so this only trims a nought.
			*/

			if ( !bWhole )
			{
				p_target->CheckForZero();

				if ( !(p_target->IsZero()) && (lzi_leading_exponent( p_target ) > l_lead) )
				{
					long	l_carried = width_for( l_lead + 1 );

					if ( l_carried < l_target )
					{
						p_target->SetPrecision( l_carried );
						p_target->Round();
					}
				}
			}

			/*
				**Trimmed before rounding and never after**, since
				2026-10-01, so a short quotient is an exact one. The trim
				above takes the noughts off a quotient that ended -- the
				guard and sticky digits nought with it -- and one that did
				not keeps every place of its width, trailing noughts and
				all. It used to trim again after rounding whenever the
				guard digit was nought, to keep the spelling the old
				truncating division gave, and so wrote a quotient a hair
				from 22.5 as 22.5.
			*/

			p_target->SetPrecision( precision );
		}
	}

	// Clean up...
DivideExit:

    delete p_dividend;
    delete p_divisor;
    delete p_quotient;

	return lResult;
}


long lzi_divide_op( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber, long precision )
{
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	// A target at the units gets a whole number, straight from the raw
	//	quotient in its own mode, as the integer quotient takes one --
	//	see lzi_asks_whole...
	if ( lzi_asks_whole( p_answer ) )
		return lzi_divide_work( leftOperand, rightOperand, targetNumber, precision, true, p_answer->OwnMode() );

	return lzi_divide_work( leftOperand, rightOperand, targetNumber, precision, false, RM_round_default );
}


// The quotient as a whole number, rounded to it once in rm_mode -- see
//	lzi_divide_work. The target keeps its precision...
static long lzi_divide_whole( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber,
							  RoundingMode rm_mode )
{
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	return lzi_divide_work( leftOperand, rightOperand, targetNumber, p_answer->m_lPrecision, true, rm_mode );
}


/*
	Sums the Taylor series for sine or cosine.

	Each term is reached from the one before it by multiplying by the
	square of the argument, so the whole series costs about as many
	multiplications as it has terms - the previous approach raised the
	argument to each power from scratch, at roughly sixteen times the
	work.

	The series is only summed for arguments already reduced into
	[0, pi/4] by lz_sine, where the term dropped after the reciprocal
	factorial table runs out was smaller than 1e-40, and since the table
	was lengthened on 2026-09-26 is smaller than 1e-48.

	**That 1e-40 is the whole reason there are two loops below.**
	kRecipFactMax is 32, so the first term the table cannot supply is
	x^33/33! for the sine and x^34/34! for the cosine. At pi/4 -- the
	largest argument the series is ever handed, everything above having
	been folded down to it -- those come to 5.6e-41 and 1.3e-42
	relative to the answers there. The table is good for a little over
	40 digits and not one more, which is exactly the default working
	width plus its guard. That is not a coincidence: it was sized for
	it, and it had no headroom at all.

	**It has eight more entries now, to 1/40!**, and serves 48 digits.
	Base sixteen showed why: 36 hex places are 44 digits, so hex's
	default width was past the table, and hex sine at 0.78125 was 740
	units out of its 36th place -- bases fourteen and fifteen, at 41 and
	42 digits, 6 and 73. Lengthened, the entries stay fifty digits wide,
	so nothing in the fast path grows (see below on widening); and the
	tabled loop now stops at the first term that rounds to nought, so
	the default path sums no more terms than it did. Base ten's answers
	did not move, the audit's and the sweep's byte for byte.

	A caller who raised the precision therefore hit a wall the table
	put there, and it was the table's **length** that did it, not the
	width of its entries -- at 100 places the sine held 42 to 50
	digits, and those numbers are the truncation above, argument by
	argument, not the 50 significant digits each entry carries.

	So past a threshold the table is set aside and every term is
	computed from the recurrence the factorials already imply:

		term(k) = term(k-2) * x^2 / (k * (k-1))

	which is one multiply and one divide where a tabled term is one
	multiply and a lookup. The divide by a small integer measures 4.5
	times the tabled multiply, which is why the table is kept for the
	precision it was sized for rather than being deleted outright.

	**Widening the table instead was tried and is the wrong answer.**
	Taking its entries to 130 significant digits works -- it reaches
	the same accuracy -- but every one of the sixteen multiplications
	in the default path is then against a 130-digit constant, and the
	circular family measured 1.71 to 2.33 times its cost for a gain no
	caller at the default precision can see. The recurrence divides by
	small integers instead, so nothing in the fast path grows at all.

	The threshold is the awkward part, because the series cannot tell a
	guard from a request. lz_cosine, lz_tangent and the three
	reciprocal functions each add kTrigGuardDigits and then call down
	through lz_sine, which adds it again, so a caller at the default 36
	arrives here as 40 by the bare sine, 44 through lz_cosine and 48
	through lz_secant. Those places are slack to be rounded away.
	Entering the computed loop on the term's own size instead -- which
	needs no threshold and has no blind spot -- cost 1.05 to 1.30 times
	at the default precision, filling guard places nobody reads. So the
	threshold allowed for the deepest stack there is, and
	kRecipFactServes + 3 * kTrigGuardDigits was what that came to.

	Its one soft edge was a bare lz_sine at 41 to 48 places, which
	stayed under the threshold and kept the table's 40 digits -- and
	hex's default width, 44 digits, landed in it. **With the table
	serving 48 the edge is gone**: the threshold is unchanged at 52,
	now kRecipFactServes plus one guard, and every call sent to the
	table gets what its width asks for, however many guards are stacked
	on it. The price is that the stacked calls now fill their guard
	places -- the cosine, secant and tangent near pi/4 cost 7 to 23 per
	cent more in bases ten and twelve -- which the early stop pays back
	elsewhere: over a full turn the family moves between -11 and +5 per
	cent.

	Since 2026-10-03 the four quotients add kTrigQuotientGuardDigits,
	six, so lz_secant arrives as 50 -- still under the threshold at the
	default. Only where a quotient takes its cosine again, near pi/2,
	does that sine run past it, which is most of what the wider guard
	costs there.
*/

/*
	How many digits of *answer* the tabled part of the series is good
	for, measured above -- 40 until the table was lengthened to 1/40! --
	and the width at which the computed tail takes over: one guard past
	what the table serves. The note in lzi_trig_series has the reasoning
	and the measurement.
*/

const long	kRecipFactServes = 48;
const long	kRecipFactThreshold = kRecipFactServes + kTrigGuardDigits;

/*
	A bound on the computed tail of the series, so that one which
	somehow failed to converge would stop rather than hang.

	It is not a precision limit. The terms fall away by a growing
	factorial against an argument below one, so the loop leaves on its
	own long before this -- around 120 terms at 200 places, 450 or so
	at a thousand. This engine has twice had a fault whose symptom was
	a hang rather than a wrong answer, and a bound costs one comparison
	a term.

	**It has to grow with the width, or it is a precision limit after
	all.** Fixed at 1024 it became one past about 2,600 places: 1/1024!
	is 1e-2640, so the exponential's fraction stopped there -- exp(0.5)
	at 4,400 places had 2,951 digits right -- and exp(10000.5), which
	works its fraction that wide, lost some two thousand places. So
	lzi_series_bound is the factorial index whose factorial passes the
	working width, with a margin, and never under 1024. A loop that
	left on its own before 1024 leaves exactly where it did.
*/

const long	kTrigSeriesMax = 1024;

static long lzi_series_bound( long lBase )
{
	double	d_wanted = static_cast<double>( glWorkingPrecision + kTrigGuardDigits ) *
					   std::log( static_cast<double>( lBase ) );
	double	d_logFactorial = 0.0;
	long	l_index = 1;

	while ( d_logFactorial <= d_wanted )
	{
		l_index++;
		d_logFactorial += std::log( static_cast<double>( l_index ) );
	}

	return ((l_index + 16) > kTrigSeriesMax) ? (l_index + 16) : kTrigSeriesMax;
}

// Both live further down, beside the exponential's kernel they follow...
long lzi_sincos_split_floor( long lBase );
static long lzi_sincos_split_work( CLZNumber* pArgument, CLZNumber* pSine, CLZNumber* pCosine, long lGuardDigits );

long lzi_trig_series( CLZNumber* argument, bool useCosineSeries, CLZNumber* targetNumber )
{
	/*
		**Past its base's floor, the kernel in place of the series** --
		stage 2 of performance item 1 of the review of 2026-10-15; see
		lzi_sincos_by_splitting and lzi_sincos_split_floor. Under the
		caller's guard and none of its own, as the exponential's kernel is
		under lz_exp_e's. Not for a flagged argument, which the table's
		arithmetic carries through to the sum -- see b_tabled below.
	*/

	if ( (glWorkingPrecision >= lzi_sincos_split_floor( argument->m_lBase )) &&
		 !(argument->FlagIsOn( kfInfinity )) && !(argument->FlagIsOn( kfQNaN )) )
	{
		return useCosineSeries ? lzi_sincos_split_work( argument, nullptr, targetNumber, 0 )
							   : lzi_sincos_split_work( argument, targetNumber, nullptr, 0 );
	}

	auto	sum = new CLZNumber;
	auto	power = new CLZNumber;
	auto	term = new CLZNumber;
	auto	argument_squared = new CLZNumber;

	// sum and power are on the left, so in the argument's base -- see
	//	CLZAnswerBase...
	sum->SetBase( argument->m_lBase );
	power->SetBase( argument->m_lBase );

	bool	b_subtract = false;
	long	l_index;

	lz_multiply( argument, argument, argument_squared );

	// Cosine runs over the even powers, sine over the odd...
	if ( useCosineSeries )
	{
		power->Assign( 1L );
		l_index = 0;
	}
	else
	{
		power->Assign( argument );
		l_index = 1;
	}

	/*
		**The threshold is in decimal digits and the working width in
		places of the argument's base**, so the width is converted
		before it is compared. Compared as they stood, a hex argument at
		44 places -- 53 digits -- stayed on a table good for 40.
	*/

	long	l_digits = CLZNumber::EquivalentPlaces( glWorkingPrecision, argument->m_lBase, 10 );

	bool	b_tabled = (l_digits <= kRecipFactThreshold) ||
					   argument->FlagIsOn( kfInfinity ) ||
					   argument->FlagIsOn( kfQNaN );

	if ( b_tabled )
	{
		/*
			The table, for as far as it goes. Every term is a lookup
			and a multiply, and each is independent of the last, so
			nothing accumulates.
		*/

		/*
			**The loop stops at the first term that rounds to nought**,
			every term after it being smaller still, the argument being
			under one. Adding a nought changes nothing, so this stops
			nothing that mattered -- and it is what lets the table be
			longer without the default path paying for it: at 40
			places the terms past x^31 round away, and a small argument
			stops sooner still. **No test catches its removal** -- every
			answer is the same without it; it is speed.

			**Not for a flagged argument.** An infinity or a
			not-a-number is held as a nought carrying a flag, so every
			term reads as nought; the flag reaches the sum only by the
			loop running on through the arithmetic that carries it. See
			b_tabled.
		*/

		bool	b_flagged = argument->FlagIsOn( kfInfinity ) || argument->FlagIsOn( kfQNaN );

		while ( l_index <= kRecipFactMax )
		{
			lz_multiply( power, gRecipFact[ l_index ], term );

			term->CheckForZero();

			if ( !b_flagged && term->IsZero() )
				break;

			if ( b_subtract )
				lz_subtract( sum, term );
			else
				lz_add( sum, term );

			b_subtract = !b_subtract;

			// Step up two powers for the next term...
			lz_multiply( power, argument_squared );

			l_index += 2;
		}
	}

	/*
		Then as many more as the working width still wants, from the
		recurrence the factorials already imply:

			term(k) = term(k-2) * x^2 / (k * (k-1))

		term still holds the last one the table supplied, so this picks
		up from it and needs neither the power nor the factorial.

		Every term computed, the table set aside. The threshold and
		why it is a threshold at all are covered at the head of this
		function; what matters here is the two things this branch has
		to get right.

		**It must not be taken for a flagged argument**, which is why
		b_tabled tests for one. An infinity or a not-a-number is held
		as a zero carrying a flag, and the loop below stops when a term
		has faded to nothing -- it reads that zero, leaves, and the
		flag never reaches sum. The tabled loop propagates it through
		the arithmetic instead, and its behaviour is what the suite
		pins: without the test, sec(infinity) comes back as 1 rather
		than a not-a-number. lz_factorial has the same trap written up
		at its head, and it is worth reading both together.

		**Each term carries the rounding of the one before**, where a
		tabled term carried none, so error compounds along the run:
		some seventy terms at 100 places, each rounded at the working
		width, costs about two digits. Two against the fifty the table
		would have capped it at is the trade, and only a caller who
		asked for more than fifty ever makes it.

		The loop ends when a term falls below the working width, which
		it must: the argument is under one and the divisors only grow.
		lzi_series_bound is a guard against a hang, not the stopping
		rule.

		The divisor stays in base ten and is matched by lz_divide,
		which costs nothing for a base-ten argument and is correct for
		any other.
	*/

	else
	{
		auto	divisor = new CLZNumber;

		term->Assign( power );		// x for the sine, 1 for the cosine

		long	l_bound = lzi_series_bound( argument->m_lBase );

		while ( l_index <= l_bound )
		{
			if ( b_subtract )
				lz_subtract( sum, term );
			else
				lz_add( sum, term );

			b_subtract = !b_subtract;

			l_index += 2;

			lz_multiply( term, argument_squared );

			divisor->Assign( l_index * (l_index - 1) );
			lz_divide( term, divisor );

			term->CheckForZero();

			if ( term->IsZero() )
				break;
		}

		delete divisor;
	}

	targetNumber->Assign( sum );

	// Clean up...
	delete sum;
	delete power;
	delete term;
	delete argument_squared;

	return kNoError;
}


long lzi_subtract_op( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber, long resultSign )
{
	if ( leftOperand->m_lBase != rightOperand->m_lBase )
		return kBaseMismatch;

	long	l_lengthRight = 0;
	long	l_adjustRight = 0;
	long*	lp_rightDigits = nullptr;

	long	l_lengthResult = 0;
	long*	lp_resultDigits = nullptr;

	long	l_carryDigit = 0;

	// Set up the result...
	l_lengthResult = targetNumber->m_lDigitsLength;

	// Get right radix adjust and digits length...
	l_adjustRight = targetNumber->m_lRadix - rightOperand->m_lRadix;
	l_lengthRight = rightOperand->m_lDigitsLength + l_adjustRight;

	// Get copies of the operand pointers...
	lp_resultDigits = targetNumber->m_lpDigits;
	lp_rightDigits = rightOperand->m_lpDigits;

	// Subtract the right side...
	for ( long i=0; i<l_lengthResult; i++ )
	{
		if ( l_carryDigit )
			*lp_resultDigits += l_carryDigit;

		if ( i < l_lengthRight )
		{
			if ( l_adjustRight )
				l_adjustRight--;
			else
				*lp_resultDigits -= *(lp_rightDigits++);
		}

		// Carry over as required...
		if ( *lp_resultDigits < 0 )
		{
			*lp_resultDigits += leftOperand->m_lBase;
			l_carryDigit = -1;
		}
		else
			l_carryDigit = 0;

		// Next result digit...
		lp_resultDigits++;
	}

	// Set the target...
	targetNumber->m_lSign = resultSign;

	return kNoError;
}


/////////////////////////////////////////////////////////////////
//	Math API functions

/*
	Whether every number written in lFrom ends in lInto: every prime of
	the one divides the other, as two's divides ten.

	**The primes, not every divisor** -- since 2026-10-21. It tested every
	divisor of lFrom, so 22 pairs read as not ending -- 4 into 2, 6, 10 or
	14; 8 into 2, 4, 6, 10, 12 or 14; 9 into 3, 6, 12 or 15; 12 into 6;
	and 16 into every even base below it -- and a mixed pair in them took
	the whole-number road and was written by the quotient's rule though it
	converts exactly: 0.333... times a quarter came back rounded to the
	width with the quarter written in base two, and 38 places long with it
	written in base four. Found building the complex layer.

	Dividing out what the two share, again and again, leaves nothing of
	lFrom exactly when each of its primes divides lInto...
*/
static bool lzi_ends_in_base( long lFrom, long lInto )
{
	long	l_from = lFrom;

	for ( long l_shared = std::gcd( l_from, lInto ); l_shared > 1; l_shared = std::gcd( l_from, lInto ) )
		l_from /= l_shared;

	return l_from == 1;
}

/*
	**A value brought into lBase keeping its relative precision**, not a
	fixed count of places -- since 2026-10-15. Exact where it ends there:
	its k places become at most 4k, four being the largest power of a
	prime among bases two to sixteen. Otherwise lSignificant places past
	its leading digit. SetBaseForWork's fixed count lost all of a small
	value: 2^-50 in base two reached base ten at ten places as nought, so
	a geometric mean of it and 1 answered nought, and an exponent of 0.1
	in base three gave 10^0.1 twenty right digits of thirty-six. The value
	keeps the precision it came with.
*/
static void lzi_into_base_keeping( CLZNumber* pValue, long lBase, long lSignificant )
{
	if ( pValue->m_lBase == lBase )
		return;

	long	l_precision = pValue->m_lPrecision;

	pValue->CheckForZero();
	pValue->TrimTrailing();

	if ( pValue->IsZero() || lzi_ends_in_base( pValue->m_lBase, lBase ) )
	{
		pValue->SetPrecision( (4 * pValue->m_lRadix) + 8 );
		pValue->SetBase( lBase, true );
		pValue->TrimTrailing();
	}
	else
	{
		long	l_lead = lzi_leading_exponent( pValue );
		long	l_noughts = (l_lead < 0) ? (CLZNumber::EquivalentPlaces( -l_lead, pValue->m_lBase, lBase ) + 1) : 0;

		pValue->SetBaseForWork( lBase, lSignificant + l_noughts );
	}

	pValue->SetPrecision( l_precision );
}

/*
	**A host's mixed pair, worked exactly** -- since 2026-10-15.
	CLZMatchedBase converts the right operand at a fixed count of the
	left base's places and the arithmetic then treats it as exact, so a
	small right operand lost its relative precision -- 1 over 2^-34
	written in base two refused at five places as a division by nought,
	and at ten answered 16666666666.67 -- and two
	values that differ past those places compared equal or the wrong
	way round.

	**Where y ends in the left base -- every prime of its base dividing
	the left one, as base two in base ten -- it is simply converted
	exactly**, and the ordinary same-base operation follows, spelled as
	it always is.

	**Where it does not, every number a base can write is a whole number
	over a power of that base**, y = Y / D with D = b^k, and whole
	numbers convert between bases exactly. So each operation becomes one
	on whole numbers and x in the left base, and the only rounding left
	is one correctly rounded division:

		compare( x, y )		compare( x D, Y )
		x / y, whole x / y	( x D ) / Y
		x + y, x - y		( x D +- Y ) / D
		x y					( x Y ) / D
		the remainders		remainder( x D, Y ) / D

	the last because scaling both operands by D leaves the quotient
	alone. **Such a sum, product or remainder is written by the
	quotient's rule**, being one, and is short only when it ends --
	Rich's call: rounded again to the width, it would round twice. Same-
	base pairs never come here, and nor does a
	registered table entry (see lzi_table_in_base), whose kept
	conversion is the engine's own working with its guards. An operand
	carrying an infinity or a NaN goes the ordinary way, to have it
	propagated.

	**Built only once Wanted() has said so**: it holds three numbers,
	and a number callocs when it is made -- the cost that took 41 per
	cent off the same-base addition when CLZMatchedBase first held its
	copy by value. The same-base call pays one comparison of bases.
*/

class CLZMixedPair
{
public:

	static bool Wanted( CLZNumber* leftOperand, CLZNumber* rightOperand )
	{
		return (leftOperand->m_lBase != rightOperand->m_lBase) && (gTabled.count( rightOperand ) == 0) &&
			   !Flagged( leftOperand ) && !Flagged( rightOperand );
	}

	CLZMixedPair( CLZNumber* leftOperand, CLZNumber* rightOperand )
		: m_pLeft( leftOperand ), m_bTerminates( lzi_ends_in_base( rightOperand->m_lBase, leftOperand->m_lBase ) )
	{
		/*
			**Where y ends in the left base, it is converted exactly** --
			every prime of its base divides the left one's, so its k places
			become at most 4k there, four being the largest power of a prime
			among bases two to sixteen -- and the ordinary same-base
			operation follows, spelled as it always is...
		*/
		if ( m_bTerminates )
		{
			m_exact.Assign( rightOperand );
			lzi_into_base_keeping( &m_exact, leftOperand->m_lBase, 0 );
			return;
		}

		long	l_saved = glWorkingPrecision;

		// ...and where it does not, y is a whole number over a power of its
		//	base, both of which convert exactly...
		m_whole.Assign( rightOperand );
		m_whole.TrimTrailing();

		long	l_places = m_whole.m_lRadix;

		lzi_scale_by_base_power( &m_whole, -l_places );
		m_scale.SetBase( rightOperand->m_lBase, false );
		m_scale.Assign( 1L );
		lzi_scale_by_base_power( &m_scale, -l_places );
		m_whole.SetBase( leftOperand->m_lBase, true );
		m_scale.SetBase( leftOperand->m_lBase, true );

		// ...and x times that power, exactly: a whole number adds no places...
		glWorkingPrecision = Exact();
		lz_multiply( leftOperand, &m_scale, &m_scaledLeft );
		glWorkingPrecision = l_saved;
	}

	long	Compare()
	{
		return m_bTerminates ? m_pLeft->Compare( &m_exact ) : m_scaledLeft.Compare( &m_whole );
	}

	long	Divide( CLZNumber* targetNumber )
	{
		return m_bTerminates ? lz_divide( m_pLeft, &m_exact, targetNumber ) :
							   lz_divide( &m_scaledLeft, &m_whole, Answer( targetNumber ) );
	}

	long	DividePrecision( long precision )
	{
		return m_bTerminates ? lz_divide( m_pLeft, &m_exact, precision ) :
							   lzi_divide_op( &m_scaledLeft, &m_whole, m_pLeft, precision );
	}

	long	DivideInteger( CLZNumber* targetNumber )
	{
		return m_bTerminates ? lz_divide_integer( m_pLeft, &m_exact, targetNumber ) :
							   lz_divide_integer( &m_scaledLeft, &m_whole, Answer( targetNumber ) );
	}

	long	Add( CLZNumber* targetNumber, bool bSubtract )
	{
		if ( m_bTerminates )
			return bSubtract ? lz_subtract( m_pLeft, &m_exact, targetNumber ) : lz_add( m_pLeft, &m_exact, targetNumber );

		CLZNumber	sum;
		long		l_saved = glWorkingPrecision;

		glWorkingPrecision = Exact();
		long	lResult = bSubtract ? lz_subtract( &m_scaledLeft, &m_whole, &sum ) : lz_add( &m_scaledLeft, &m_whole, &sum );
		glWorkingPrecision = l_saved;

		return lzi_first_failure( lResult, lz_divide( &sum, &m_scale, Answer( targetNumber ) ) );
	}

	long	Multiply( CLZNumber* targetNumber )
	{
		if ( m_bTerminates )
			return lz_multiply( m_pLeft, &m_exact, targetNumber );

		CLZNumber	product;
		long		l_saved = glWorkingPrecision;

		glWorkingPrecision = Exact();
		long	lResult = lz_multiply( m_pLeft, &m_whole, &product );
		glWorkingPrecision = l_saved;

		return lzi_first_failure( lResult, lz_divide( &product, &m_scale, Answer( targetNumber ) ) );
	}

	// lz_remainder, lz_remainder_near or lz_modulo: scaling both operands
	//	by the same power leaves the quotient alone...
	long	Remainder( LZBinaryFunction fn, CLZNumber* targetNumber )
	{
		if ( m_bTerminates )
			return fn( m_pLeft, &m_exact, targetNumber );

		CLZNumber	scaled;
		long		l_saved = glWorkingPrecision;

		glWorkingPrecision = Exact();
		long	lResult = fn( &m_scaledLeft, &m_whole, &scaled );
		glWorkingPrecision = l_saved;

		return lzi_first_failure( lResult, lz_divide( &scaled, &m_scale, Answer( targetNumber ) ) );
	}

private:

	CLZMixedPair( const CLZMixedPair& );
	CLZMixedPair& operator=( const CLZMixedPair& );

	static bool Flagged( CLZNumber* p )
	{
		return p->FlagIsOn( kfInfinity ) || p->FlagIsOn( kfQNaN ) || p->FlagIsOn( kfSNaN );
	}

	// Places enough that no intermediate rounds: x's own, and the width...
	long	Exact() const
	{
		long	l_places = (m_pLeft->m_lRadix > glWorkingPrecision) ? m_pLeft->m_lRadix : glWorkingPrecision;

		return l_places + m_scaledLeft.m_lRadix + 2;
	}

	// Written where the operation writes: the target, or the left operand...
	CLZNumber*	Answer( CLZNumber* targetNumber ) const
	{
		return (targetNumber != nullptr) ? targetNumber : m_pLeft;
	}

	CLZNumber*	m_pLeft;
	bool		m_bTerminates;
	CLZNumber	m_exact;
	CLZNumber	m_whole;
	CLZNumber	m_scale;
	CLZNumber	m_scaledLeft;
};

/*
	**An odd function of a tiny argument, decided by its side of x** --
	since 2026-10-15, item 6 of that day's review. sin, tan, arctan,
	arcsin and the four hyperbolic and inverse hyperbolic kin are x plus
	a cubic term and smaller, so when x itself lies on a half of the
	answer's places -- its last digit half an even base, one place past
	them, or for an answer written to W significant digits one digit past
	those -- only that cubic term moves the answer off it, at three times
	the argument's depth: past the retake's cap, where the half read
	exact and went by the kept digit. sin of 0.00000000055 at ten places
	answered ...6 where sin x < x makes it ...5.

	**|f(x)| - |x| has a known sign and is smaller than |x|^3** wherever
	|x| is below one over the base, at most a half, for all eight: the
	cubic coefficients are a sixth and a third, and the rest of each
	series adds less than the cubic's own size again. So with y the answer
	cut at its places and u a unit there, 2|f(x)| against 2y + u is
	2|x| - (2y + u), worked exactly, plus twice a term of that sign and
	under 2|x|^3: decided where the two agree in sign, where the exact
	part is nought, or where it outweighs the bound. Anything else
	declines, and the retake goes on as before.
*/

static bool lzi_beside_x_decide( CLZNumber* pResult, CLZNumber* pX, long lPlaces, long lSide )
{
	if ( (lPlaces < 0) || pX->FlagIsOn( kfInfinity ) || pX->FlagIsOn( kfQNaN ) )
		return false;

	long		l_base = pResult->m_lBase;
	CLZNumber	x;

	x.Assign( pX );
	x.m_lSign = kPositive;
	x.CheckForZero();
	x.TrimTrailing();

	if ( x.IsZero() || (x.m_lBase != l_base) )
		return false;

	// |x| < base^(lead + 1), and below one over the base...
	long	l_lead = lzi_leading_exponent( &x );

	if ( l_lead > -2 )
		return false;

	long	l_saved = glWorkingPrecision;
	long	l_places = (x.m_lRadix > lPlaces) ? x.m_lRadix : lPlaces;

	glWorkingPrecision = l_places + 4;

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	twice_mid;
	CLZNumber	twice_x;
	CLZNumber	gap;
	CLZNumber	bound;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &twice_mid );
	lz_add( &twice_mid, &unit );
	lz_add( &x, &x, &twice_x );
	lz_subtract( &twice_x, &twice_mid, &gap );
	gap.CheckForZero();

	// 2 base^(3 (lead + 1)), over twice the cubic term...
	bound.SetBase( l_base, false );
	bound.Assign( 2L );
	lzi_scale_by_base_power( &bound, -3 * (l_lead + 1) );

	long	l_gapSide = gap.IsZero() ? 0 : (gap.IsNegative() ? -1 : 1);
	bool	b_up = false;

	if ( l_gapSide == 0 )
		b_up = (lSide > 0);
	else if ( l_gapSide == lSide )
		b_up = (l_gapSide > 0);
	else
	{
		CLZNumber	size;

		size.Assign( &gap );
		size.m_lSign = kPositive;

		if ( lz_compare( &size, &bound ) != kGreaterThan )
		{
			glWorkingPrecision = l_saved;
			return false;
		}

		b_up = (l_gapSide > 0);
	}

	glWorkingPrecision = lPlaces + 2;

	if ( b_up )
		lz_add( &y, &unit );

	y.m_lSign = pResult->m_lSign;
	y.CheckForZero();
	pResult->Assign( &y );
	pResult->SetPrecision( lPlaces );

	glWorkingPrecision = l_saved;

	return true;
}

// sin, arctan, arcsinh and tanh lie inside x...
static bool lzi_inside_x_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long /* lMostDigits */ )
{
	return lzi_beside_x_decide( pResult, pLeft, lPlaces, -1 );
}

// ...and tan, sinh, arcsin and arctanh outside it...
static bool lzi_outside_x_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long /* lMostDigits */ )
{
	return lzi_beside_x_decide( pResult, pLeft, lPlaces, 1 );
}

/*
	**sin and cos in a unit, a hair from a half** -- since 2026-10-23,
	found building complex stage 6b. sin 30 degrees and cos 60 are 1/2,
	which in an odd base is a half at every width -- 0.333... in base
	seven -- and in an even base a half at no places. An angle a hair from
	such a one lies a hair from the half past any retake, whose widest
	answer still looks exactly like the half, and half-even then took a
	side whichever side the angle lay: sin(30 - 7^-60) degrees at ten
	places in base seven answered 0.3333333334 for ...333.

	**The side is exact.** Where the half in doubt is +-1/2, the angle r in
	[0, a turn) lies a hair from k twelfths of a turn, and the function is
	that half there, the answer lies on the side of the half that r lies of
	k/12 of a turn, times the sign of the slope there -- cos 30k degrees
	for the sine, -sin 30k for the cosine, neither nought at those k. 12r
	against k turns is exact, the remainder that makes r being exact.
	Nothing else in a unit is a half: the tangent and cotangent at an
	eighth are one, the secant and cosecant at a twelfth two. pScale is
	lzc_polar's modulus r, exact, whose part r cos theta lies on the half
	r/2 a hair from the same angles: the half is then |r| times one half,
	and the side times r's sign.
*/

// Whether the half in doubt at lPlaces is |r| / 2, r one where pScale is
//	null -- 2 Y + 1 against |r| b^P, Y = y b^P -- with y, |the answer| cut at
//	the places, into pY...
static bool lzi_answer_is_half( CLZNumber* pResult, long lPlaces, CLZNumber* pScale, CLZNumber* pY )
{
	long		l_base = pResult->m_lBase;
	CLZNumber	Y;
	CLZNumber	twice;
	CLZNumber	one;
	CLZNumber	b_power;

	pY->Assign( pResult );
	pY->m_lSign = kPositive;
	pY->CheckForZero();
	pY->Round( RM_round_down, lPlaces );
	Y.Assign( pY );
	lzi_scale_by_base_power( &Y, -lPlaces );
	lzi_times_whole( &Y, 2, &twice );
	one.SetBase( l_base, false );
	one.Assign( 1L );
	lz_add( &twice, &one, &twice );
	b_power.Assign( (pScale != nullptr) ? pScale : &one );
	b_power.m_lSign = b_power.IsZero() ? kZero : kPositive;
	lzi_scale_by_base_power( &b_power, -lPlaces );

	return lz_compare( &twice, &b_power ) == kEqualTo;
}

// An angle x in a unit of lTurn a hair from k twelfths of a turn: the sine's
//	or cosine's value there as +-1, a half of that sign, and the sign of its
//	slope, and 12 x - k turns, exact, x taken into [0, a turn) -- false
//	where the value there is no half...
static bool lzi_near_twelfth( CLZNumber* pX, long lFunction, long lTurn, long* plValue, long* plSlope, CLZNumber* pTwelveDx )
{
	long		l_base = pX->m_lBase;
	CLZNumber	turn;
	CLZNumber	r;

	turn.SetBase( l_base, false );
	turn.Assign( lTurn );
	lz_remainder( pX, &turn, &r );
	r.CheckForZero();

	if ( r.IsNegative() )
		lz_add( &r, &turn );

	double	d_r = 0.0;

	lzi_value_or_infinity( &r, &d_r );

	long	l_k = static_cast<long>( std::floor( d_r * 12.0 / static_cast<double>( lTurn ) + 0.5 ) );

	*plValue = 0;
	*plSlope = 0;

	if ( lFunction == kAngleSine )
	{
		switch ( l_k )
		{
		case 1:		*plValue = 1;	*plSlope = 1;	break;
		case 5:		*plValue = 1;	*plSlope = -1;	break;
		case 7:		*plValue = -1;	*plSlope = -1;	break;
		case 11:	*plValue = -1;	*plSlope = 1;	break;
		default:	break;
		}
	}
	else
	{
		switch ( l_k )
		{
		case 2:		*plValue = 1;	*plSlope = -1;	break;
		case 4:		*plValue = -1;	*plSlope = -1;	break;
		case 8:		*plValue = -1;	*plSlope = 1;	break;
		case 10:	*plValue = 1;	*plSlope = 1;	break;
		default:	break;
		}
	}

	if ( *plValue == 0 )
		return false;

	CLZNumber	twelve_r;
	CLZNumber	k_turns;

	lzi_times_whole( &r, 12, &twelve_r );
	k_turns.SetBase( l_base, false );
	k_turns.Assign( l_k * lTurn );
	lz_subtract( &twelve_r, &k_turns, pTwelveDx );
	pTwelveDx->CheckForZero();

	return true;
}

static bool lzi_unit_half_decide( CLZNumber* pResult, CLZNumber* pX, long lPlaces, long lFunction, CLZNumber* pScale = nullptr )
{
	long	l_turn = lzi_turn_of( geAngleUnit );

	if ( (lPlaces < 0) || (l_turn == 0) || pX->FlagIsOn( kfInfinity ) || pX->FlagIsOn( kfQNaN ) )
		return false;

	long	l_base = pResult->m_lBase;

	if ( pX->m_lBase != l_base )
		return false;

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pX->m_lRadix > lPlaces) ? pX->m_lRadix : lPlaces) + 4;

	if ( (pScale != nullptr) && (pScale->m_lRadix + lPlaces + 4 > glWorkingPrecision) )
		glWorkingPrecision = pScale->m_lRadix + lPlaces + 4;

	CLZNumber	y;
	long		l_value = 0;
	long		l_slope = 0;
	CLZNumber	twelve_dx;
	long		l_scale = ((pScale != nullptr) && pScale->IsNegative()) ? -1 : 1;
	long		l_sign_h = pResult->IsNegative() ? -1 : 1;

	// The half in doubt |r| / 2, the angle a hair from a twelfth where the
	//	function is a half of the answer's sign, and not on it...
	if ( !lzi_answer_is_half( pResult, lPlaces, pScale, &y ) ||
		 !lzi_near_twelfth( pX, lFunction, l_turn, &l_value, &l_slope, &twelve_dx ) ||
		 (l_value * l_scale != l_sign_h) || twelve_dx.IsZero() )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// The answer's side of the half: the angle's side of k twelfths, by the
	//	slope and the modulus's sign...
	long	l_side = (twelve_dx.IsPositive() ? 1 : -1) * l_slope * l_scale;
	bool	b_up = (l_side * l_sign_h > 0);

	CLZNumber	unit;

	glWorkingPrecision = lPlaces + 2;
	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	if ( b_up )
		lz_add( &y, &unit );

	y.m_lSign = pResult->m_lSign;
	y.CheckForZero();
	pResult->Assign( &y );
	pResult->SetPrecision( lPlaces );

	glWorkingPrecision = l_saved;

	return true;
}

// The sine's and the cosine's, as an answer guard's tie decider...
static bool lzi_sine_unit_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long /* lMostDigits */ )
{
	return lzi_unit_half_decide( pResult, pLeft, lPlaces, kAngleSine );
}

static bool lzi_cosine_unit_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long /* lMostDigits */ )
{
	return lzi_unit_half_decide( pResult, pLeft, lPlaces, kAngleCosine );
}

/*
	**arccos and arccosh a hair from one, decided by their square** --
	since 2026-10-15, item 9 of that day's review. arccos(1 - d) is
	sqrt(2d) (1 + d/12 + ...) and arccosh(1 + d) is sqrt(2d) (1 - d/12
	+ ...), so where sqrt(2d) lies on a half of the answer's places --
	an operand some twice the places wide -- only the d/12 term moves it
	off, past the retake's cap: arccos of 0.99999999999999999921875 at
	ten places answered ...12 for ...13.

	**Squared, it is exact.** With y the answer cut at its places and u a
	unit there, the answer against the half is 2 f against M = 2y + u,
	and f = sqrt(2d) r with r strictly between 1 + d/12 and 1 + d/12 +
	d^2/32 for arccos, 1 - d/12 and 1 - d/12 + d^2/32 for arccosh,
	wherever d is below one over the base -- the series' later terms add
	less than d^2/32 there. So 8d r^2 against M^2 at the lower r says
	above if it is not below, and at the upper says below if it is not
	above, all in finite numbers; between the two it declines. A first
	version bounded the correction by d/6 alone and left four of the
	census's near ties, where the moved operand and the d/12 term nearly
	cancel.
*/

static bool lzi_near_one_decide( CLZNumber* pResult, CLZNumber* pX, long lPlaces, long lSide )
{
	if ( (lPlaces < 0) || pX->FlagIsOn( kfInfinity ) || pX->FlagIsOn( kfQNaN ) || pResult->IsNegative() )
		return false;

	long		l_base = pResult->m_lBase;
	long		l_saved = glWorkingPrecision;
	CLZNumber	one;
	CLZNumber	d;

	if ( pX->m_lBase != l_base )
		return false;

	glWorkingPrecision = pX->m_lRadix + 4;

	one.SetBase( l_base, false );
	one.Assign( 1L );

	// d = 1 - x for arccos, x - 1 for arccosh; a hair, positive...
	if ( lSide > 0 )
		lz_subtract( &one, pX, &d );
	else
		lz_subtract( pX, &one, &d );

	d.CheckForZero();
	d.TrimTrailing();

	if ( d.IsZero() || d.IsNegative() || (lzi_leading_exponent( &d ) > -2) )
	{
		glWorkingPrecision = l_saved;
		return false;
	}

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	m;
	CLZNumber	m2;
	CLZNumber	eight;
	CLZNumber	e;

	glWorkingPrecision = lPlaces + 4;

	y.Assign( pResult );
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &m );
	lz_add( &m, &unit );

	// Products' places are their operands' together...
	glWorkingPrecision = (2 * m.m_lRadix) + d.m_lRadix + 8;

	lz_multiply( &m, &m, &m2 );
	eight.SetBase( l_base, false );
	eight.Assign( 8L );
	lz_multiply( &eight, &d, &e );

	// The answer over sqrt(2d) lies strictly between 1 -+ d/12 and that
	//	plus d^2/32: P/Q of 12 -+ d over 12, and of 96 -+ 8d + 3d^2 over
	//	96. Above the half if 8d P^2 >= M^2 Q^2 at the lower, below it if
	//	8d P^2 <= M^2 Q^2 at the upper...
	auto	against = [ & ]( bool bUpper ) -> long
	{
		CLZNumber	p;
		CLZNumber	p2;
		CLZNumber	q;
		CLZNumber	left;
		CLZNumber	right;
		CLZNumber	q2;

		glWorkingPrecision = (4 * m.m_lRadix) + (6 * d.m_lRadix) + 16;

		if ( !bUpper )
		{
			p.SetBase( l_base, false );
			p.Assign( 12L );
			q.SetBase( l_base, false );
			q.Assign( 12L );

			if ( lSide > 0 )
				lz_add( &p, &d );
			else
				lz_subtract( &p, &d );
		}
		else
		{
			CLZNumber	term;
			CLZNumber	d2;

			p.SetBase( l_base, false );
			p.Assign( 96L );
			q.SetBase( l_base, false );
			q.Assign( 96L );
			term.SetBase( l_base, false );
			term.Assign( 8L );
			lz_multiply( &term, &d );

			if ( lSide > 0 )
				lz_add( &p, &term );
			else
				lz_subtract( &p, &term );

			lz_multiply( &d, &d, &d2 );
			term.Assign( 3L );
			lz_multiply( &term, &d2 );
			lz_add( &p, &term );
		}

		lz_multiply( &p, &p, &p2 );
		lz_multiply( &e, &p2, &left );
		lz_multiply( &q, &q, &q2 );
		lz_multiply( &m2, &q2, &right );

		return lz_compare( &left, &right );
	};

	bool	b_up = false;

	if ( against( false ) != kLessThan )
		b_up = true;
	else if ( against( true ) != kGreaterThan )
		b_up = false;
	else
	{
		glWorkingPrecision = l_saved;
		return false;
	}

	glWorkingPrecision = lPlaces + 2;

	if ( b_up )
		lz_add( &y, &unit );

	pResult->Assign( &y );
	pResult->SetPrecision( lPlaces );

	glWorkingPrecision = l_saved;

	return true;
}

// arccos lies above sqrt( 2 (1 - x) )...
static bool lzi_arccosine_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long /* lMostDigits */ )
{
	return lzi_near_one_decide( pResult, pLeft, lPlaces, 1 );
}

// ...and arccosh below sqrt( 2 (x - 1) )...
static bool lzi_archypercosine_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long /* lMostDigits */ )
{
	return lzi_near_one_decide( pResult, pLeft, lPlaces, -1 );
}

long lz_abs( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	// A constant is not written, in place either -- until 2026-10-15 an
	//	operand was made positive whether it was constant or not...
	if ( p_answer->IsConstant() )
		return kNoError;

	if ( targetNumber != nullptr )
		targetNumber->Assign( leftOperand );

	/*
		**A nought keeps the sign a nought has**, kZero -- since 2026-10-15,
		item 10 of that day's review. Every answer was made kPositive, and
		Compare reads the sign first, so abs(0) ranked above nought by
		lz_compare, Compare, lz_max and lz_min alike, however the nought
		was made. Only a negative is turned: every nought the engine makes
		carries kZero already, which twelve ways of making one confirmed...
	*/
	if ( p_answer->m_lSign == kNegative )
		p_answer->m_lSign = kPositive;

	return kNoError;
}


long lz_add( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The target's rounding settings are its own -- see CLZTargetRounding...
	CLZTargetRounding	targetRounding( (targetNumber != nullptr) ? targetNumber : leftOperand );

	// Before anything -- see CLZMatchedBase.
	//	A digit array carries no base, so a mixed pair has to be
	//	brought into one before any of it is read...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Add( targetNumber, false );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	// An infinity or a NaN carries through rather than being lost...
	if ( lzi_propagate_flags( leftOperand, rightOperand, targetNumber, kOpAdd ) )
		return kNoError;

	CLZNumber*	maxOperand = nullptr;
	CLZNumber*	minOperand = nullptr;

	auto result = new CLZNumber;

	unsigned long	maxDigits = 0;
	unsigned long	maxExponent = 0;

	// Set up the result as required...
	if ( (leftOperand->m_lDigitsLength - leftOperand->m_lRadix) >
		 (rightOperand->m_lDigitsLength - rightOperand->m_lRadix) )

		maxDigits = leftOperand->m_lDigitsLength - leftOperand->m_lRadix;

	else
		maxDigits = rightOperand->m_lDigitsLength - rightOperand->m_lRadix;

	if ( leftOperand->m_lRadix > rightOperand->m_lRadix )
		maxExponent = leftOperand->m_lRadix;

	else
		maxExponent = rightOperand->m_lRadix;

	// Identify the 'big' operand...
	if ( leftOperand->CompareABS( rightOperand ) == kLessThan )
	{
		maxOperand = rightOperand;
		minOperand = leftOperand;
	}
	else
	{
		maxOperand = leftOperand;
		minOperand = rightOperand;
	}

	// Align target storage...
	result->Assign( maxOperand );
	result->AlignStorage( maxDigits, maxExponent, kPreserve );

	// Perform the operation...
	if ( (leftOperand->m_lSign == rightOperand->m_lSign) && (leftOperand->m_lSign != kZero) )
		lzi_add_op( maxOperand, minOperand, result, leftOperand->m_lSign );			// Cases I and IV

	else if ( (leftOperand->m_lSign == kZero) && (rightOperand->m_lSign == kZero) )
		result->Assign( static_cast<long>(0) );

	else if ( (leftOperand->m_lSign == kZero) || (rightOperand->m_lSign == kZero) )
		result->Assign( maxOperand );

	else
		lzi_subtract_op( maxOperand, minOperand, result, maxOperand->m_lSign );		// Cases II and III

    result->TrimLeading();

	// Record the result, a whole number if the target asks for one...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetRounding.RoundExactIfWhole();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			targetRounding.RoundExactIfWhole();
		}
	}

	// Clean up...
	delete result;

	return kNoError;
}


/////////////////////////////////////////////////////////////////
//	The means of a list: arithmetic, geometric and harmonic -- and the
//	arithmetic-geometric mean of two, which shares their preparation
//
//	The only functions in the engine that take a list rather than one
//	or two operands. The target is required -- there is no left operand
//	to fall back on, and writing into the first value would be a
//	surprise nobody asked for. The values are never modified, so the
//	target may be one of them. lz_agm takes two operands like any
//	other function, and hands them to the same preparation as a list.
//
//	**Each one is a small calculation with a trap in it**, and the
//	traps are the three this engine keeps meeting, because precision
//	counts places after the radix mark:
//
//	  - the arithmetic mean sums exactly and divides once, and needs
//	    nothing but a guard;
//
//	  - the geometric mean multiplies, and a product of small values
//	    underflows: 1e-21 * 4e-21 is 4e-42, which 36 places holds as
//	    nought, so the naive geometric mean came back nought -- and
//	    unlike division, multiplication really does round to the
//	    working places, which a mutation confirmed. The exponents are taken out and summed as
//	    integers instead, and only the mantissas are multiplied;
//
//	  - the harmonic mean takes reciprocals, which never underflow: a
//	    small quotient keeps the width plus two significant digits. But
//	    a mean of large values has many whole places to fill from them,
//	    so it is widened by those; see lz_hmean.
/////////////////////////////////////////////////////////////////

/*
	What every one of the three does before its own work: check the
	list, answer a flagged value, and copy the values into the first
	one's base at the working width.

	**An infinity or a not-a-number makes the mean a not-a-number.** The
	engine holds an infinity as a nought carrying a flag, often with no
	sign, and every one of these needs the sign -- infinity plus minus
	infinity has no average, and the other two are defined only where
	nothing is negative. It is the same reading lzi_propagate_unary
	takes as kUnaryNeedsTheSign, and it has to come before any test of
	a value, which would read the infinity as nought.

	Answers kNoError with the result written when a flag decided it,
	kFailed for a list that is not one, and -1 to mean "carry on" with
	the copies made. The copies belong to the caller, freed with
	delete[].
*/

static const long	kMeanCarryOn = -1;

static long lzi_mean_prepare( CLZNumber** values, long count, CLZNumber* targetNumber,
							  CLZNumber** ppCopies )
{
	*ppCopies = nullptr;

	if ( (values == nullptr) || (count < 1) || (targetNumber == nullptr) )
		return kFailed;

	bool	b_flagged = false;
	bool	b_fromDivideByZero = false;

	for ( long i = 0; i < count; i++ )
	{
		if ( values[ i ] == nullptr )
			return kFailed;

		if ( values[ i ]->FlagIsOn( kfQNaN ) || values[ i ]->FlagIsOn( kfInfinity ) )
			b_flagged = true;

		if ( values[ i ]->FlagIsOn( kfDiv0 ) )
			b_fromDivideByZero = true;
	}

	if ( b_flagged )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( static_cast<long>(0) );
			targetNumber->SetFlag( kfQNaN );

			if ( b_fromDivideByZero )
				targetNumber->SetFlag( kfDiv0 );
		}

		return kNoError;
	}

	/*
		Copies, in the first value's base -- the rule every binary
		operator follows, the left operand's base winning.

		**SetPrecision after Assign, not before.** Assign carries the
		source's precision along with its digits, so a copy of a number
		made at 36 places stays at 36 however wide this is running --
		the trap lzi_log_gamma fell into.

		**No test can see this line today, and that is recorded rather
		than hidden.** Every copy is only ever the right operand of an
		operation, and an operation takes its width from its left one, so
		a copy at the caller's precision changes no answer; a mutation
		removing it was caught by nothing. It is kept because the first
		change that puts a copy on the left would otherwise inherit that
		trap silently.
	*/

	CLZNumber*	p_copies = new CLZNumber[ count ];
	long		l_base = values[ 0 ]->m_lBase;

	for ( long i = 0; i < count; i++ )
	{
		p_copies[ i ].Assign( values[ i ] );
		p_copies[ i ].SetPrecision( glWorkingPrecision );

		// ...keeping each value's relative precision -- see
		//	lzi_into_base_keeping...
		lzi_into_base_keeping( &p_copies[ i ], l_base, glWorkingPrecision + 4 );

		p_copies[ i ].CheckForZero();
	}

	*ppCopies = p_copies;

	return kMeanCarryOn;
}


/*
	Whether every value in the list is the same one.

	**The mean of identical values is that value, and every one of the
	three should say so exactly.** Only the arithmetic mean did without
	help. The harmonic mean of 7 alone went through 1/7 and back and
	answered 7.000000000000000000000000000000000000 -- the right value,
	written as a rounded one -- and the geometric mean of four 5s went
	through lz_root_n. A list of one value is the commonest case of
	this, and the check costs a comparison per value.

	Asked after the domain checks, so a list of identical negatives is
	still refused by the geometric and harmonic means rather than
	answered.
*/

static bool lzi_mean_all_equal( CLZNumber* pCopies, long count )
{
	for ( long i = 1; i < count; i++ )
	{
		if ( lz_compare( &pCopies[ 0 ], &pCopies[ i ] ) != kEqualTo )
			return false;
	}

	return true;
}


/*
	Write the answer, rounded back to the width the caller is working
	at. Shared by all three so the constant-target rule lives once.
*/

// The product of the values, all but the one at lSkip (or none, for
//	-1), worked wide enough that nothing is rounded...
static void lzi_exact_product( CLZNumber** values, long count, long lSkip, CLZNumber* pResult )
{
	long	l_saved = glWorkingPrecision;
	long	l_places = 1;

	for ( long i = 0; i < count; i++ )
		l_places += values[ i ]->m_lRadix;

	glWorkingPrecision = l_places;

	CLZNumber	product;

	product.SetBase( values[ 0 ]->m_lBase, false );
	product.Assign( 1L );

	for ( long i = 0; i < count; i++ )
		if ( i != lSkip )
			lz_multiply( &product, values[ i ] );

	glWorkingPrecision = l_saved;

	pResult->Assign( &product );
}


/*
	**A short geometric or harmonic mean is checked against the values
	themselves** before it stands -- see lzi_write_to_width. Both work
	through rounded values, a product or a sum of reciprocals, and an
	answer can check out exactly against those: gmean(2, 8 + 1e-42)
	came back 4. So the geometric mean's power has to be the values'
	own product, and the harmonic mean times the sum of the products of
	all but one value has to be the count times the product of all,
	each worked exactly. Past sixteen values it is not asked, and a
	short answer is written out to its width.
*/

static void lzi_mean_spelling( CLZNumber** values, long count, CLZNumber* pAnswer, bool bHarmonic )
{
	const long	kMostCheckedValues = 16;

	if ( pAnswer->IsConstant() || !lzi_short_of_width( pAnswer ) )
		return;

	bool	b_exact = false;

	if ( count <= kMostCheckedValues )
	{
		CLZNumber	left;
		CLZNumber	right;

		if ( bHarmonic )
		{
			CLZNumber	all;
			CLZNumber	sum;
			CLZNumber	n;
			CLZNumber	answer;
			long		l_saved = glWorkingPrecision;

			lzi_exact_product( values, count, -1, &all );
			sum.SetBase( values[ 0 ]->m_lBase, false );
			sum.Assign( 0L );

			for ( long i = 0; i < count; i++ )
			{
				CLZNumber	others;

				lzi_exact_product( values, count, i, &others );
				lz_add( &sum, &others );
			}

			sum.TrimTrailing();
			all.TrimTrailing();
			answer.Assign( pAnswer );
			answer.TrimTrailing();

			glWorkingPrecision = answer.m_lRadix + sum.m_lRadix + all.m_lRadix + 1;

			n.SetBase( values[ 0 ]->m_lBase, false );
			n.Assign( count );
			lz_multiply( &answer, &sum, &left );
			lz_multiply( &n, &all, &right );

			glWorkingPrecision = l_saved;
		}
		else
		{
			lzi_exact_product( values, count, -1, &right );
			lzi_exact_power( pAnswer, count, &left );
		}

		b_exact = (lz_compare( &left, &right ) == kEqualTo);
	}

	if ( b_exact )
		pAnswer->TrimTrailing();
	else
		lzi_write_to_width( pAnswer, false );
}


// Which mean is being recorded, for the check on a short one...
const long	kMeanUnchecked = 0;
const long	kMeanGeometric = 1;
const long	kMeanHarmonic = 2;

static void lzi_mean_record( CLZNumber* result, CLZNumber* targetNumber, long lSavedPrecision,
							 CLZNumber** values = nullptr, long count = 0, long lMean = kMeanUnchecked )
{
	glWorkingPrecision = lSavedPrecision;

	// The final rounding's tail, kept before the spelling check's own
	//	arithmetic records over it -- see CLZListRetake. Forgotten first,
	//	so a result with nothing to discard leaves none...
	CLZNumber::ForgetLastRounding();

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	gdMeanTail = CLZNumber::LastDiscardedTail();
	glMeanPlaces = CLZNumber::LastRoundedPlaces();

	// Checked before it is recorded: the target may be one of the values...
	if ( lMean != kMeanUnchecked )
		lzi_mean_spelling( values, count, result, lMean == kMeanHarmonic );

	if ( !(targetNumber->IsConstant()) )
		targetNumber->Assign( result );
}


/*
	The arithmetic-geometric mean of two values, for values neither of
	which is negative.

		a  <-  (a + g) / 2,		g  <-  sqrt( a * g )

	The two close on each other and meet at the AGM, doubling the
	correct digits every cycle once they are near. It is not the
	average of anything -- that is lz_amean -- but it sits here because
	it takes its operands through the same preparation as the means of
	a list, as a list of two.

	**A negative is refused, as MPFR refuses it**, and a nought anywhere
	makes the answer nought, exactly. Before, a negative answered nought
	with a not-a-number flag and kNoError -- the failure was the square
	root's, and nothing passed its code on -- and a nought answered
	5 / 2^15 for agm(0, 5): the geometric mean stays nought, the
	arithmetic one only halves, and sixteen fixed cycles stopped there.

	**The iteration works between one and the base, and the answer is
	scaled back.** agm(k*x, k*y) = k * agm(x, y), so both are divided
	by b^e, e being the larger one's exponent, which is exact. Then the
	two traps the geometric mean has, for the same reason:

	  - **the first product underflows** where the smaller value is
	    small beside the larger. agm(1e-40, 1) multiplied them to
	    1e-40, nought at 36 places, and answered 3.05e-5 for 0.0168.
	    The width is raised by the geometric mean's exponent, which is
	    what makes that product hold its digits -- the smaller value's
	    at first, and less every cycle as the two close;

	  - **the scale back up costs places**, b^e moving the radix, so
	    the width is raised by e as well where e is positive -- the
	    cube root's band.

	And it stops when it has converged rather than after a count. The
	cycles needed are those that bring the two together, which grow
	only with the log of the log of their ratio, plus those that fill
	the width, which grow with its log. Sixteen were enough at 36
	places even for a ratio of 1e1000, and not for agm(1e-200, 1e200),
	whose answer needs some 240 places: stopped there it had half its
	digits. agm(1, 2) now takes six.

	Worked to four guard digits, which is where it stops helping: over
	800 random pairs, half of them close in size, 679 came back
	correctly rounded at none, 786 at one, 798 at two, and all of them
	at three or more.

	The answer goes to the target, or to the left operand where there is
	none; neither operand is otherwise modified, so the target may be
	the right one.
*/

static const long	kAgmMaxCycles = 256;

// For the tests: the most places past the first cycle's width of the means
//	going into it, since the last read, this thread's, which the read sets
//	back -- nought or less, LONG_MIN for none...
static thread_local long gAgmFirstOverrun = LONG_MIN;
static thread_local long gAgmCycles = 0;

long lzi_agm_first_overrun()
{
	long	l_overrun = gAgmFirstOverrun;

	gAgmFirstOverrun = LONG_MIN;

	return l_overrun;
}

// ...and the cycles this thread's AGMs have taken...
long lzi_agm_cycles()
{
	return gAgmCycles;
}

long lz_agm( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_agm, rightOperand, kRetakeDoubtWide );

	long	l_savedPrecision = glWorkingPrecision;

	CLZNumber*	values[] = { leftOperand, rightOperand };

	glWorkingPrecision += lzi_guard_places( kAgmGuardDigits, lzi_mean_base( values, 2 ) );
	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	CLZNumber*	p_copies = nullptr;
	long		lResult = lzi_mean_prepare( values, 2, p_target, &p_copies );

	if ( lResult != kMeanCarryOn )
	{
		glWorkingPrecision = l_savedPrecision;
		return lResult;
	}

	lResult = kNoError;

	if ( p_copies[ 0 ].IsNegative() || p_copies[ 1 ].IsNegative() )
	{
		delete[] p_copies;
		glWorkingPrecision = l_savedPrecision;

		return lzi_domain_failure( p_target, p_target, kfQNaN, kFailed );
	}

	CLZNumber	result;

	// Nought, and identical values below, answer exactly -- see
	//	lzi_write_to_width...
	if ( p_copies[ 0 ].IsZero() || p_copies[ 1 ].IsZero() )
	{
		result.Assign( static_cast<long>(0) );
		answer.Exact();
	}

	/*
		The AGM of identical values is that value, and the iteration
		reaches it anyway -- a mutation removing this changed no value
		over 300 identical pairs, only whether a trailing nought the
		value was written with survives. Kept for that, for the cycles it
		saves, and to answer as the means beside it do.
	*/

	else if ( lzi_mean_all_equal( p_copies, 2 ) )
	{
		result.Assign( &p_copies[ 0 ] );
		answer.Exact();
	}

	else
	{
		CLZNumber*	p_larger = &p_copies[ 0 ];
		CLZNumber*	p_smaller = &p_copies[ 1 ];

		if ( lz_compare( p_larger, p_smaller ) == kLessThan )
		{
			p_larger = &p_copies[ 1 ];
			p_smaller = &p_copies[ 0 ];
		}

		// Both divided by b^e, so the larger lands in [1, b)...
		long	l_scale = lzi_leading_exponent( p_larger );

		lzi_scale_by_base_power( p_larger, l_scale );
		lzi_scale_by_base_power( p_smaller, l_scale );

		// The width every cycle needs whatever the values are...
		long	l_width = glWorkingPrecision;

		if ( l_scale > 0 )
			l_width += l_scale;

		CLZNumber	a_value;
		CLZNumber	g_value;
		CLZNumber	half;
		CLZNumber	two;

		a_value.SetBase( p_copies[ 0 ].m_lBase, false );
		g_value.SetBase( p_copies[ 0 ].m_lBase, false );
		half.SetBase( p_copies[ 0 ].m_lBase, false );
		two.SetBase( p_copies[ 0 ].m_lBase, false );

		/*
			**In an odd base the arithmetic mean halves by a division**,
			since 2026-10-03. It multiplied by the reciprocal of two, which
			there is .333... cut at the working width and not a half, and
			every cycle carried the shortfall on: 14 of 800 random pairs in
			base seven came back a unit out, by up to 0.024 of one past the
			half. For an answer far from one it was far worse, the
			shortfall having been fixed at the width before the answer's
			places were added and then scaled back up with it: every one
			of 300 random pairs with the first operand near 7^8 or 7^30
			came back wrong, by up to 6e22 units of the 36th place. The
			square root's half was the same fault until 2026-10-01.

			**In an even base the reciprocal is exact, and is kept**: a
			division by two measured 1.4 to 1.6 times the whole AGM in
			base ten, against a multiplication that was never wrong there.
			Base seven, by contrast, got 39 times faster: with the inexact
			half no two cycles ever agreed, and every pair ran to the cap
			of 256 cycles, where it now settles in 7.

			The division is in place, so it works to the sum's precision,
			which is set to the cycle's width first: a multiplication
			rounds to glWorkingPrecision whatever its target holds, and a
			division does not.
		*/

		bool	b_halveByDivision = (p_copies[ 0 ].m_lBase % 2) != 0;

		half.Assign( static_cast<long>(2) );
		two.Assign( static_cast<long>(2) );

		if ( !b_halveByDivision )
			lz_reciprocal( &half );

		a_value.Assign( p_larger );
		g_value.Assign( p_smaller );

		CLZNumber	a_next;
		CLZNumber	g_next;

		for ( long lCycle = 0; lCycle < kAgmMaxCycles; lCycle++ )
		{
			/*
				...and the width the smaller value needs on top, taken
				afresh each cycle. It is the geometric mean that is
				small, and its exponent halves every cycle, so the
				width falls back to l_width within a few. Holding the
				first cycle's width throughout took agm(1, 1e-300) to
				650ms, where narrowing takes 15. The arithmetic mean is
				the larger, so a width that holds the geometric mean
				holds it too.

				**The operands are narrowed to the first cycle's width, and
				only there** -- since 2026-10-17, the Newton audit. This
				said the working width was all there was to narrowing, the
				operations rounding to it whatever their operands carry,
				which holds for every cycle's results and not for what goes
				in: a product works every digit it is handed, and the first
				cycle multiplied the operands at their own length, 66,426
				places over for one three widths long at 10,000 digits.
				Narrowed, such an AGM is up to 30% faster, and an ordinary
				one level; now and then one of base two's takes a cycle
				more and is up to 10% slower. Every cycle after the first
				takes means its own operations rounded.

				**Rounding the means every cycle instead is wrong**: a_value
				rounded never equals the a_next the division leaves a place
				longer, so in an odd base the second settling test below
				could not fire and pairs ran to the cap, 17 to 24 times
				slower.
			*/

			long	l_gExponent = lzi_leading_exponent( &g_value );

			glWorkingPrecision = l_width + ((l_gExponent < 0) ? -l_gExponent : 0);
			gAgmCycles++;

			if ( lCycle == 0 )
			{
				if ( a_value.m_lRadix > glWorkingPrecision )
					a_value.Round( kDefaultRoundingMode, glWorkingPrecision );

				if ( g_value.m_lRadix > glWorkingPrecision )
					g_value.Round( kDefaultRoundingMode, glWorkingPrecision );

				// For the tests, read by lzi_agm_first_overrun...
				long	l_over = std::max( a_value.m_lRadix, g_value.m_lRadix ) - glWorkingPrecision;

				if ( l_over > gAgmFirstOverrun )
					gAgmFirstOverrun = l_over;
			}

			lResult = lzi_first_failure( lResult, lz_add( &a_value, &g_value, &a_next ) );
			if ( b_halveByDivision )
			{
				a_next.SetPrecision( glWorkingPrecision );
				lResult = lzi_first_failure( lResult, lz_divide( &a_next, &two ) );
			}
			else
				lResult = lzi_first_failure( lResult, lz_multiply( &a_next, &half ) );

			lResult = lzi_first_failure( lResult, lz_multiply( &a_value, &g_value, &g_next ) );
			lResult = lzi_first_failure( lResult, lz_root_2( &g_next ) );

			/*
				Settled when the two agree, or when the arithmetic mean
				stops moving -- they can end a place apart, each rounding
				to its own side, and never meet.

				**No test can see the second test, and it is not
				optional.** Without it agm(1e-200, 1e200) never meets,
				runs all kAgmMaxCycles and answers the same thing in
				3.8s rather than 0.39. The cap is only there against
				that failing too.
			*/

			bool	b_settled = ( lz_compare( &a_next, &g_next ) == kEqualTo ) ||
								( lz_compare( &a_next, &a_value ) == kEqualTo );

			a_value.Assign( &a_next );
			g_value.Assign( &g_next );

			if ( b_settled || (lResult != kNoError) )
				break;
		}

		result.Assign( &a_value );

		// ...and times b^e, which is why the width went up by e...
		lzi_scale_by_base_power( &result, -l_scale );
	}

	lzi_mean_record( &result, p_target, l_savedPrecision );

	delete[] p_copies;

	return lResult;
}


/*
	(x1 + x2 + ... + xn) / n.

	**Rounded once, at the width** -- since 2026-10-15, item 8 of that
	day's review. The sum is exact -- addition in this engine loses
	nothing -- and so is the sum times b^W, so that over n rounded to a
	whole number and moved back is the answer correctly rounded, in one
	step. Until then the division rounded at the width and its guard and
	the record rounded again at the width: amean( 0.1234567892,
	0.1234567893000000002 ) at ten places answered ...892 for ...893, and
	a list of 100,001 values misrounded the same way. That step is taken
	only where the guard lies within a unit of the half, the one place two
	roundings can differ; worked every time it cost a quarter of the call.
*/

long lz_amean( CLZNumber** values, long count, CLZNumber* targetNumber )
{
	// The engine's own work -- see CLZInternalWork...
	CLZInternalWork	work( targetNumber );

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kMeanGuardDigits, lzi_mean_base( values, count ) );

	CLZNumber*	p_copies = nullptr;
	long		lResult = lzi_mean_prepare( values, count, targetNumber, &p_copies );

	if ( lResult != kMeanCarryOn )
	{
		glWorkingPrecision = l_savedPrecision;
		return lResult;
	}

	lResult = kNoError;

	CLZNumber	sum;
	CLZNumber	n;
	CLZNumber	result;

	if ( lzi_mean_all_equal( p_copies, count ) )
	{
		result.Assign( &p_copies[ 0 ] );
		lzi_mean_record( &result, targetNumber, l_savedPrecision );

		delete[] p_copies;

		return kNoError;
	}

	sum.SetBase( p_copies[ 0 ].m_lBase, false );
	sum.Assign( static_cast<long>(0) );
	n.SetBase( p_copies[ 0 ].m_lBase, false );
	n.Assign( count );

	for ( long i = 0; i < count; i++ )
		lResult = lzi_first_failure( lResult, lz_add( &sum, &p_copies[ i ] ) );

	lResult = lzi_first_failure( lResult, lz_divide( &sum, &n, &result ) );

	/*
		The quotient is correctly rounded a guard wider, so rounding it
		again at the width is right wherever its tail lies more than a unit
		of its last place clear of the half: it cannot have been rounded
		across. Within that unit the answer is S b^W over n rounded to a
		whole number once and moved back. It is never exact there -- one
		that ended within the width would have had no tail to doubt -- so
		it is written to the width...
	*/
	if ( (result.m_lRadix > l_savedPrecision) && lzi_tail_in_doubt( &result, l_savedPrecision, 1 ) )
	{
		lzi_scale_by_base_power( &sum, -l_savedPrecision );
		lResult = lzi_first_failure( lResult, lzi_divide_whole( &sum, &n, &result, kDefaultRoundingMode ) );
		lzi_scale_by_base_power( &result, l_savedPrecision );
	}

	lzi_mean_record( &result, targetNumber, l_savedPrecision );

	delete[] p_copies;

	return lResult;
}


/*
	The nth root of (x1 * x2 * ... * xn), for values none of which is
	negative.

	**A negative value is refused rather than answered**, even where
	the product comes out positive: an even number of negatives
	multiplies to something with a real root, but that root is not the
	geometric mean of anything. A nought anywhere makes the answer
	nought, exactly.

	**The exponents are taken out and summed as integers, and only the
	mantissas are multiplied.** Each value is m * b^e with m in [1, b),
	so the product is (m1 * m2 * ... * mn) * b^E with E the sum of the
	exponents. The mantissa product sits in [1, b^n) whatever the
	values were, so it cannot underflow and keeps its relative
	precision; multiplying the values directly took the geometric mean
	of 1e-20 with itself to nought, the product 1e-40 being below the
	working width.

	Then E = q*n + r, so the root splits cleanly:

		(M * b^E)^(1/n)  =  (M * b^r)^(1/n) * b^q

	**And the b^q at the end has to be paid for in places.** Scaling up
	by a power of the base moves the radix, so a root worked to W places
	comes out with W - q once it is multiplied back -- the cube root's
	band, again. The width is raised by q before any of the work is
	done, which is possible here because q is known from the exponents
	alone, before a single multiplication.

	Two, three and one value go to lz_root_2, lz_root_3 and nothing, so
	the geometric mean of 4 and 9 is 6 exactly and not 5.999...; more
	values go through lz_root_n, which is right to within its last place
	and so right once the guard is rounded away.

	lzi_scale_by_base_power( x, p ) DIVIDES by b^p for a positive p --
	it moves the radix left. Easy to get backwards.
*/

long lz_gmean( CLZNumber** values, long count, CLZNumber* targetNumber )
{
	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZListRetake...
	CLZListRetake	retake( lz_gmean, values, count, targetNumber, kRetakeDoubt );

	// The engine's own work -- see CLZInternalWork...
	CLZInternalWork	work( targetNumber );

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kMeanGuardDigits, lzi_mean_base( values, count ) );

	CLZNumber*	p_copies = nullptr;
	long		lResult = lzi_mean_prepare( values, count, targetNumber, &p_copies );

	if ( lResult != kMeanCarryOn )
	{
		glWorkingPrecision = l_savedPrecision;
		return lResult;
	}

	lResult = kNoError;

	bool	b_zero = false;

	for ( long i = 0; i < count; i++ )
	{
		if ( p_copies[ i ].IsNegative() )
		{
			delete[] p_copies;
			glWorkingPrecision = l_savedPrecision;

			return lzi_domain_failure( targetNumber, targetNumber, kfQNaN, kFailed );
		}

		if ( p_copies[ i ].IsZero() )
			b_zero = true;
	}

	CLZNumber*	p_result = new CLZNumber;

	if ( b_zero )
	{
		p_result->Assign( static_cast<long>(0) );
	}

	else if ( lzi_mean_all_equal( p_copies, count ) )
	{
		p_result->Assign( &p_copies[ 0 ] );
	}

	else
	{
		long	l_exponentSum = 0;

		for ( long i = 0; i < count; i++ )
		{
			long	l_exponent = lzi_leading_exponent( &p_copies[ i ] );

			l_exponentSum += l_exponent;

			// Down to its mantissa, in [1, b)...
			lzi_scale_by_base_power( &p_copies[ i ], l_exponent );
		}

		/*
			E = q*n + r. Truncating toward zero leaves r negative for a
			negative sum, and that is fine: the split is exact either way,
			b^r only moves M's radix before the root, and the root of a
			value under one is as good as the root of one over it. Written
			as a floor first, and a mutation back to plain division changed
			no answer.
		*/

		long	l_quotient = l_exponentSum / count;
		long	l_remainder = l_exponentSum % count;

		if ( l_quotient > 0 )
			glWorkingPrecision += l_quotient;

		delete p_result;
		p_result = new CLZNumber;

		CLZNumber	product;

		product.SetBase( p_copies[ 0 ].m_lBase, false );
		product.Assign( static_cast<long>(1) );

		for ( long i = 0; i < count; i++ )
		{
			p_copies[ i ].SetPrecision( glWorkingPrecision );
			lResult = lzi_first_failure( lResult, lz_multiply( &product, &p_copies[ i ] ) );
		}

		// Times b^r...
		lzi_scale_by_base_power( &product, -l_remainder );

		if ( count == 1 )
		{
			p_result->Assign( &product );
		}
		else if ( count == 2 )
		{
			lResult = lzi_first_failure( lResult, lz_root_2( &product, p_result ) );
		}
		else if ( count == 3 )
		{
			lResult = lzi_first_failure( lResult, lz_root_3( &product, p_result ) );
		}
		else
		{
			CLZNumber	degree;

			degree.SetBase( p_copies[ 0 ].m_lBase, false );
			degree.Assign( count );

			lResult = lzi_first_failure( lResult, lz_root_n( &product, &degree, p_result ) );
		}

		// ...and times b^q, which is why the width went up by q...
		lzi_scale_by_base_power( p_result, -l_quotient );
	}

	lzi_mean_record( p_result, targetNumber, l_savedPrecision, values, count, kMeanGeometric );

	delete p_result;
	delete[] p_copies;

	return lResult;
}


/*
	n / (1/x1 + 1/x2 + ... + 1/xn), for values none of which is
	negative.

	**A nought anywhere makes the answer nought.** Its reciprocal is
	infinite, so the sum is, and n over it is nought -- a finite limit
	approached from the positive side, and the only one the definition
	allows, so it is answered rather than refused. A negative is
	refused, as for the geometric mean.

	**Widened by the answer's whole places**, which are at most the
	smallest value's plus the count's: the mean lies between the
	smallest value and n times it. Each reciprocal keeps the width plus
	two significant digits, so a mean with forty whole places needs
	forty more of them for n over their sum to be right to its places.

	This was first written with a widening by the largest value's
	magnitude, and a mutation removing it was caught by no test: until
	2026-09-29 division counted a small quotient's noughts twice, so
	1/3e40 at 36 places came back holding 118 significant digits and
	the widening was dead. Counting them once took hmean(1e40, 3e40)
	to 15000000000000000000000000000000000000000.00375, and the
	widening is back, sized by the answer rather than the largest value.
*/

long lz_hmean( CLZNumber** values, long count, CLZNumber* targetNumber )
{
	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZListRetake...
	CLZListRetake	retake( lz_hmean, values, count, targetNumber, kRetakeDoubt );

	// The engine's own work -- see CLZInternalWork...
	CLZInternalWork	work( targetNumber );

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kMeanGuardDigits, lzi_mean_base( values, count ) );

	CLZNumber*	p_copies = nullptr;
	long		lResult = lzi_mean_prepare( values, count, targetNumber, &p_copies );

	if ( lResult != kMeanCarryOn )
	{
		glWorkingPrecision = l_savedPrecision;
		return lResult;
	}

	lResult = kNoError;

	bool	b_zero = false;

	for ( long i = 0; i < count; i++ )
	{
		if ( p_copies[ i ].IsNegative() )
		{
			delete[] p_copies;
			glWorkingPrecision = l_savedPrecision;

			return lzi_domain_failure( targetNumber, targetNumber, kfQNaN, kFailed );
		}

		if ( p_copies[ i ].IsZero() )
			b_zero = true;
	}

	CLZNumber	result;

	if ( b_zero )
	{
		result.Assign( static_cast<long>(0) );
	}

	else if ( lzi_mean_all_equal( p_copies, count ) )
	{
		result.Assign( &p_copies[ 0 ] );
	}

	else
	{
		// The answer's whole places, at most the smallest value's plus
		//	the count's -- see above. Before the numbers below are made,
		//	since each takes the working precision it is made at...
		long	l_least = p_copies[ 0 ].GetOrderOfMagnitude();
		long	l_countPlaces = 0;

		for ( long i = 1; i < count; i++ )
		{
			long	l_magnitude = p_copies[ i ].GetOrderOfMagnitude();

			if ( l_magnitude < l_least )
				l_least = l_magnitude;
		}

		for ( long c = count; c > 0; c /= p_copies[ 0 ].m_lBase )
			l_countPlaces++;

		glWorkingPrecision += l_least + l_countPlaces;

		CLZNumber	sum;
		CLZNumber	one;
		CLZNumber	reciprocal;
		CLZNumber	n;

		sum.SetBase( p_copies[ 0 ].m_lBase, false );
		sum.Assign( static_cast<long>(0) );
		one.SetBase( p_copies[ 0 ].m_lBase, false );
		one.Assign( static_cast<long>(1) );
		n.SetBase( p_copies[ 0 ].m_lBase, false );
		n.Assign( count );

		for ( long i = 0; i < count; i++ )
		{
			lResult = lzi_first_failure( lResult, lz_divide( &one, &p_copies[ i ], &reciprocal ) );
			lResult = lzi_first_failure( lResult, lz_add( &sum, &reciprocal ) );
		}

		lResult = lzi_first_failure( lResult, lz_divide( &n, &sum, &result ) );
	}

	lzi_mean_record( &result, targetNumber, l_savedPrecision, values, count, kMeanHarmonic );

	delete[] p_copies;

	return lResult;
}


/*
	The three inverse circular reciprocals -- lz_arccosecant,
	lz_arccotangent and lz_arcsecant -- each go through the inverse
	they are the co-function of, and none has a series of its own.
	They are filed alphabetically like everything else here, so the
	other arc functions sit between them.

	Two of them share an awkwardness worth stating once. arccsc and
	arcsec of an infinity are nought and pi/2, from either direction:
	the sign does not matter, but the answer is finite. This engine's
	unary propagation offers two readings, carry the infinity or give
	a not-a-number, and a finite limit is neither. They give a
	not-a-number rather than a number invented for an argument that
	has none. See lzi_propagate_unary.
*/

/*
	**How much wider a reciprocal has to be worked near |x| = 1**, for the
	four co-functions whose inverse has a slope without bound there:
	arccsc and arcsec take the arcsine and arccosine on 1/x, arcsech the
	hyperbolic arccosine and arccoth the hyperbolic arctangent. 1/x is
	rounded before the inverse sees it, and near one that rounding
	arrives multiplied by the slope -- 1/sqrt(x^2 - 1) for the first
	three, 1/(x^2 - 1) for arccoth. A fixed guard cannot cover a slope
	with no bound: until 2026-09-26 arccoth missed correct rounding by
	up to 600 million units in base ten, arcsech by 95 and arccsc and
	arcsec by 9 -- and those after their inverses were mended.

	So as many places as ||x| - 1| has noughts after the point, which
	covers arccoth's slope and more than covers the square-root ones.
	Nothing away from one: a gap of a tenth or more has none.
*/

static long lzi_near_one_places( CLZNumber* pValue )
{
	CLZNumber	one;
	CLZNumber	magnitude;
	CLZNumber	gap;

	one.SetBase( pValue->m_lBase );
	one.Assign( 1L );

	lz_abs( pValue, &magnitude );
	lz_subtract( &magnitude, &one, &gap );		// exact

	gap.CheckForZero();

	if ( gap.IsZero() )
		return 0;

	long	l_exponent = lzi_leading_exponent( &gap );

	return (l_exponent < 0) ? -(l_exponent + 1) : 0;
}


long lz_arccosecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arccosecant, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// In another unit, exact where the angle is a part of a turn, and
	//	otherwise converted from radians -- see lzi_inverse_in_unit...
	if ( e_unit != AU_radians )
		return lzi_inverse_in_unit( kAngleArcCosecant, leftOperand, targetNumber, answer, e_unit );

	/*
		Before any test of the value, an infinity being held as zero
		-- which would read as nought and be refused as outside the
		domain, which is the right refusal for the wrong reason.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber magnitude;
	CLZNumber reciprocal;
	CLZNumber result;

	// Set up before the domain test, which compares against one...
	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	/*
		Defined only from one outwards, which is exactly where the
		reciprocal lands inside the arcsine's own domain. Nought is
		refused here too: it is the far end of the same exclusion, not
		a pole, since the cosecant never returns a value below one.
	*/

	if ( lz_compare( &magnitude, &one ) == kLessThan )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		arccsc(x) = arcsin(1/x), reciprocal and arcsine both taken at
		guard width -- see kArcCoFunctionGuardDigits. Without it, one
		part in 1e10 outside the domain's edge costs two digits.

		The reciprocal's own precision is set explicitly rather than
		left to its construction: lz_divide takes its width from the
		target it is handed, which was built before the working
		precision went up.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase ) +
						  lzi_near_one_places( leftOperand );

	reciprocal.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( &one, leftOperand, &reciprocal ) );
	lResult = lzi_first_failure( lResult, lz_arcsine( &reciprocal, &result ) );

	glWorkingPrecision = l_savedPrecision;

	// The guard places never leave the function...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_arccosine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arccosine, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and an operand a hair from one decided by the answer's square --
	//	see lzi_near_one_decide...
	//	In radians only: a decider reads the argument as the angle itself.
	if ( e_unit == AU_radians )
		answer.DecideTies( lzi_arccosine_decide );

	// Exact at one, its nought spelled short -- until 2026-10-15 it was
	//	written to the full width, though nothing about it is in doubt...
	answer.ExactAndShort( lzi_exactly( leftOperand, 1 ) );

	// In another unit, exact where the angle is a part of a turn, and
	//	otherwise converted from radians -- see lzi_inverse_in_unit...
	if ( e_unit != AU_radians )
		return lzi_inverse_in_unit( kAngleArcCosine, leftOperand, targetNumber, answer, e_unit );

	/*
		arccos of an infinity is outside the domain whichever sign it
		had, and an infinity here is held as zero with no sign anyway.
		A not-a-number, then, and before any test of the value, which
		would read the infinity as nought and answer pi/2.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber magnitude;
	CLZNumber arcsine;
	CLZNumber result;

	// Set up before the domain test, which compares against one...
	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	// Outside [-1, 1] there is no angle to find...
	if ( lz_compare( &magnitude, &one ) == kGreaterThan )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		arccos(x) = pi/2 - arcsin(x), and nothing cleverer.

		The received advice is that this cancels as x approaches 1 --
		the arcsine climbs to pi/2, the difference shrinks to nought,
		and the leading digits go with it -- and that the half-angle
		form, 2*arcsin(sqrt((1-x)/2)), avoids it by computing a small
		answer from a small argument. Both were built and scored side
		by side, and the half-angle form is never better here and
		sometimes worse: 35 digits against 38 at x = .99, 31 against
		32 at nine nines.

		The advice assumes a representation carrying a fixed number of
		significant digits, where a small answer keeps its precision.
		This one counts places after the radix mark, so a small answer
		has few significant digits whatever produced it -- both forms
		land on the same absolute error and the same floor of roughly
		36 + log10 of the answer. Measured against that floor, this is
		already sitting on it.

		**Guard digits were once tried here and removed, and that was
		wrong** -- see kArcComplementGuardDigits. They were measured in
		digits of agreement, which cannot see the fault: the arcsine
		rounds its answer to the working width, and the subtraction
		rounds again, so the error reaches nearly a unit while never
		costing a whole digit. A hex pi/2 is a rounded conversion too.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcComplementGuardDigits, leftOperand->m_lBase );

	arcsine.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	CLZNumber	halfPi( kPiOverTwo );

	// It is on the left, so in the operand's base -- see CLZAnswerBase...
	lzi_pi_family_in_base( &halfPi, kPiOverTwo, leftOperand->m_lBase, glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_arcsine( leftOperand, &arcsine ) );
	lResult = lzi_first_failure( lResult, lz_subtract( &halfPi, &arcsine, &result ) );

	glWorkingPrecision = l_savedPrecision;

	// kPiOverTwo is 120 places wide and the arcsine taken off it has
	//	the working precision, so round back to what everything else
	//	here carries rather than returning pi/2's tail. Measured at
	//	1.00x with the constant wide, so it is not narrowed first --
	//	one subtraction cannot feel the difference the way lz_sine's
	//	comparisons and lz_cosine's addition did...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}

long lz_arccotangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arccotangent, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// In another unit, exact where the angle is a part of a turn, and
	//	otherwise converted from radians -- see lzi_inverse_in_unit...
	if ( e_unit != AU_radians )
		return lzi_inverse_in_unit( kAngleArcCotangent, leftOperand, targetNumber, answer, e_unit );

	/*
		arccot is nought at an infinity from either direction, a finite
		limit -- the arccosecant's awkwardness, and its answer: a
		not-a-number, since this engine's unary propagation carries the
		infinity or gives a not-a-number and a finite limit is neither.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber arctangent;
	CLZNumber result;
	CLZNumber magnitude;

	/*
		arccot(x) = arctan(1/x), pi/2 at nought: the range (-pi/2,
		pi/2], odd, and a jump of pi at nought -- **Rich's call,
		2026-10-23, the convention Elzed answers, and DLMF's and
		Mathematica's**, so a host moving to the engine meets the
		answers it had. Until then the engine took the continuous
		convention, pi/2 - arctan(x), range (0, pi), as lz_arccosine is
		pi/2 - arcsin; the two differ only on the negatives, where
		arccot(-1) was 3pi/4 and is -pi/4 now.

		Worked as pi/2 - arctan(|x|), and negated for a negative x:
		the subtraction is the one the continuous form took on the
		positives, and a correctly rounded answer negated is correctly
		rounded, half-even rounding being symmetric. Nought is not
		negative, so it answers pi/2.

		No cancellation needs guarding. The received worry is the cancellation
		for large x, where the answer shrinks to nought while both
		terms sit at pi/2 -- but the same reasoning as the arccosine's
		applies, and for the same reason: this representation counts
		places after the radix mark, so a small answer has few
		significant digits whatever produced it. Both forms land on the
		same absolute error. Built and scored against arctan(1/x) both
		sides of the reflection, and there is nothing between them.

		**The double rounding does need guarding**, as the arccosine's
		does -- see kArcComplementGuardDigits.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcComplementGuardDigits, leftOperand->m_lBase );

	arctangent.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	CLZNumber	halfPi( kPiOverTwo );

	// It is on the left, so in the operand's base -- see CLZAnswerBase...
	lzi_pi_family_in_base( &halfPi, kPiOverTwo, leftOperand->m_lBase, glWorkingPrecision );

	leftOperand->CheckForZero();

	bool	b_negative = leftOperand->IsNegative();

	magnitude.Assign( leftOperand );

	if ( b_negative )
		lz_negate( &magnitude );

	lResult = lzi_first_failure( lResult, lz_arctangent( &magnitude, &arctangent ) );
	lResult = lzi_first_failure( lResult, lz_subtract( &halfPi, &arctangent, &result ) );

	glWorkingPrecision = l_savedPrecision;

	// kPiOverTwo is 120 places wide and the arctangent taken off it
	//	has the working precision, so round back to what everything
	//	else here carries rather than returning pi/2's tail. Not
	//	narrowed first, for the reason given in lz_arccosine...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	if ( b_negative )
		lz_negate( &result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


/*
	The three inverse hyperbolic reciprocals -- lz_archypercosecant,
	lz_archypercotangent and lz_archypersecant -- go through the
	inverse each is the co-function of, on the reciprocal of the
	argument, exactly as the circular three do. Alphabetical order
	again puts the other arc functions between them.

	All three are nought at an infinity, from either direction. The
	same awkwardness as the arccosecant's, and the same answer: a
	not-a-number, because this engine's unary propagation says carry
	the infinity or give a not-a-number and a finite limit is neither.
*/

long lz_archypercosecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_archypercosecant, kRetakeDoubt );

	// Before any test of the value, an infinity being held as zero
	//	-- which would read as nought and be reported as the pole...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber reciprocal;
	CLZNumber result;

	/*
		Defined for every real but nought, where the limit is an
		infinity that cannot be written down. kUndefinedResult with
		kfInfinity, following lz_log_e at zero rather than the flat
		refusal an argument outside the domain gets.
	*/

	leftOperand->CheckForZero();

	if ( leftOperand->IsZero() )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kUndefinedResult );

	/*
		arccsch(x) = arcsinh(1/x), guarded -- see
		kArcCoFunctionGuardDigits -- though for a different reason from
		the other four, and it took three readings to settle which.

		There is no slope to guard against here: the hyperbolic
		arcsine's never exceeds one, so nothing amplifies the
		reciprocal's last place, and the argument for going without was
		sound as far as it went.

		The first reading watched whether the answer *changed* as the
		guard widened. It does, at x = 1.000001 and x = 1000, and that
		was nearly written down as the guard being needed. The second
		scored those same arguments against ground truth: identical to
		the digit, guarded or not. A changed answer is not a better
		one.

		The third let the magnitude sweep ask, and it is the sweep that
		is right: the guard is worth exactly one digit at four points
		neither hand-picked list contained -- +-0.1 and +-1e-40. Small,
		real, and cheap, so it stays. What it buys is not slope but the
		reciprocal's own last place, which an argument with no exact
		reciprocal hands straight on.

		It is odd, and comes out odd on its own: the reciprocal keeps
		the sign and the arcsinh carries it through.
	*/

	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase );

	reciprocal.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( &one, leftOperand, &reciprocal ) );
	lResult = lzi_first_failure( lResult, lz_archypersine( &reciprocal, &result ) );

	glWorkingPrecision = l_savedPrecision;

	// The guard places never leave the function...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_archypercosine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_archypercosine, kRetakeDoubt );

	// ...and an operand a hair from one decided by the answer's square --
	//	see lzi_near_one_decide...
	answer.DecideTies( lzi_archypercosine_decide );

	// Exact at one and nowhere else, asked before anything is done to the
	//	operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 1 ) );

	/*
		arccosh of a positive infinity is infinite, of a negative one
		outside the domain. Two different answers, so a not-a-number
		rather than a guess -- and before any test of the value, an
		infinity being held as zero.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber inner;
	CLZNumber result;

	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	// Defined only at or above one...
	if ( lz_compare( leftOperand, &one ) == kLessThan )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		arccosh(x) = ln( x + sqrt(x^2 - 1) ).

		**x^2 - 1 as (x - 1)(x + 1), kept exact, as the arcsine's
		1 - sin^2 t is** -- see lz_arcsine. It was squared, rounded and
		one taken off, which near one kept only the rounding's worth of
		digits, and the root carried that straight into the answer:
		until 2026-09-26 more than half of all arguments within 1e-2 to
		1e-15 of one missed correct rounding, in base ten by up to four
		million units and in hex by 240 million. Worked wider too -- see
		kArcHyperCosineGuardDigits -- and rounded once.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcHyperCosineGuardDigits, leftOperand->m_lBase );

	inner.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	{
		CLZNumber	below;
		CLZNumber	above;
		long		l_width = glWorkingPrecision;

		lResult = lzi_first_failure( lResult, lz_subtract( leftOperand, &one, &below ) );
		lResult = lzi_first_failure( lResult, lz_add( leftOperand, &one, &above ) );

		// Multiplication rounds to the working places, so the width
		//	holds the product, for the one multiply...
		glWorkingPrecision = below.m_lRadix + above.m_lRadix + 1;

		if ( glWorkingPrecision < l_width )
			glWorkingPrecision = l_width;

		lResult = lzi_first_failure( lResult, lz_multiply( &below, &above, &inner ) );

		glWorkingPrecision = l_width;
	}

	lResult = lzi_first_failure( lResult, lz_root_2( &inner ) );
	lResult = lzi_first_failure( lResult, lz_add( &inner, leftOperand ) );
	lResult = lzi_first_failure( lResult, lz_log_e( &inner, &result ) );

	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_archypercotangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_archypercotangent, kRetakeDoubt );

	// Before any test of the value, an infinity being held as zero...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber magnitude;
	CLZNumber reciprocal;
	CLZNumber result;

	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	/*
		Defined strictly beyond one, the mirror of the hyperbolic
		arctangent's strict inside. The two ends are told apart:
		outside the domain there is no answer at all, while at exactly
		one the limit is an infinity that cannot be written down, so
		that one is kUndefinedResult with kfInfinity by the rule
		lz_log_e follows at zero.

		lz_archypertangent refuses its own ends the other way, with a
		flat kFailed. Noted rather than copied.
	*/

	if ( lz_compare( &magnitude, &one ) == kEqualTo )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kUndefinedResult );

	if ( lz_compare( &magnitude, &one ) == kLessThan )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		arccoth(x) = arctanh(1/x), guarded for the reason set out at
		kArcCoFunctionGuardDigits: just outside one the reciprocal is
		just inside it, which is where the hyperbolic arctangent is at
		its worst.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase ) +
						  lzi_near_one_places( leftOperand );

	reciprocal.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( &one, leftOperand, &reciprocal ) );
	lResult = lzi_first_failure( lResult, lz_archypertangent( &reciprocal, &result ) );

	glWorkingPrecision = l_savedPrecision;

	// The guard places never leave the function...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


/*
	The hyperbolic arcsecant's own series, for a small argument:

		arcsech(x) = ln(2/x) - sum from n=1 of c_n * x^(2n)

		c_1 = 1/4,  c_n = c_(n-1) * (2n-1)(2n-2) / (4 n^2)

	which is the large-argument expansion of the hyperbolic arccosine
	written on 1/x. It buys one thing: no square root. The route
	through lz_archypercosine spends about a tenth of its time in
	lz_root_2 at the default precision and better than a third of it
	at three hundred places, and this pays none of that.

	It buys nothing else, and that is worth stating plainly. Roughly
	four fifths of either route is one call to lz_log_e -- both need
	exactly one logarithm and no rearrangement of the algebra avoids
	it -- so the square root is the whole of the prize. Measured
	against the arccosine route it is 1.1 to 1.2 times faster where it
	is used, never more.

	Three things about summing it here were each found the hard way.

	One, the term is rounded back to the working precision at the end
	of every pass. Multiplication is exact, so an unrounded term gains
	places every time round and the loop's cost climbs with it.

	Two, that rounding goes after the whole update and not in the
	middle of it. Rounding between the multiply and the divide leaves
	the term stalled at the last place it can hold instead of falling
	through it.

	Three, and the reason this stops on a comparison rather than on
	reaching nought: the term never reaches nought. Each pass
	multiplies by x^2 and rounds back, and once x^2 is above a half
	that is more than half a unit in the last place, which half-even
	rounding returns whole. The sequence arrives at a fixed point and
	sits on it forever -- a decaying series that cannot decay. The
	fixed point is not one unit either, but depends on x^2: at x = .9
	the term settles on two, because twice .81 rounds back to two, so
	a test against a single unit never fires. Asking whether the term
	still shrinks is the one test that does not care where it lands.
*/

static long lzi_archypersecant_series( CLZNumber* argument, CLZNumber* targetNumber )
{
	long	lResult = kNoError;

	CLZNumber	two;
	CLZNumber	quotient;
	CLZNumber	logarithm;
	CLZNumber	argument_squared;
	CLZNumber	term;
	CLZNumber	previous;
	CLZNumber	sum;
	CLZNumber	scratch;
	CLZNumber	factor;

	quotient.SetPrecision( glWorkingPrecision );
	logarithm.SetPrecision( glWorkingPrecision );
	argument_squared.SetPrecision( glWorkingPrecision );
	term.SetPrecision( glWorkingPrecision );
	sum.SetPrecision( glWorkingPrecision );
	scratch.SetPrecision( glWorkingPrecision );

	// two and sum are on the left and factor divides every term, so all
	//	three are in the argument's base -- see CLZAnswerBase...
	two.SetBase( argument->m_lBase );
	sum.SetBase( argument->m_lBase );
	factor.SetBase( argument->m_lBase );

	two.Assign( 2L );

	lResult = lzi_first_failure( lResult, lz_divide( &two, argument, &quotient ) );
	lResult = lzi_first_failure( lResult, lz_log_e( &quotient, &logarithm ) );

	lResult = lzi_first_failure( lResult,
					lz_multiply( argument, argument, &argument_squared ) );

	// The first term is x^2 / 4...
	term.Assign( &argument_squared );
	factor.Assign( 4L );
	lResult = lzi_first_failure( lResult, lz_divide( &term, &factor ) );

	long	l_index = 1;

	while ( true )
	{
		lz_add( &sum, &term );

		if ( (l_index > 1) && (lz_compare( &term, &previous ) != kLessThan) )
			break;

		previous.Assign( &term );

		l_index++;

		// term *= x^2 * (2n-1)(2n-2) / (4 n^2)...
		lz_multiply( &term, &argument_squared, &scratch );
		term.Assign( &scratch );

		factor.Assign( (2 * l_index - 1) * (2 * l_index - 2) );
		lz_multiply( &term, &factor );

		factor.Assign( 4 * l_index * l_index );
		lResult = lzi_first_failure( lResult, lz_divide( &term, &factor ) );

		term.SetPrecision( glWorkingPrecision );
		term.Round();
	}

	lResult = lzi_first_failure( lResult, lz_subtract( &logarithm, &sum, targetNumber ) );

	return lResult;
}


long lz_archypersecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_archypersecant, kRetakeDoubt );

	// Exact at one and nowhere else, asked before anything is done to the
	//	operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 1 ) );

	// Before any test of the value, an infinity being held as zero...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber reciprocal;
	CLZNumber result;

	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	/*
		Defined on (0, 1]. Nought is the pole -- the limit there is an
		infinity -- and everything below it, along with everything
		above one, is outside the domain altogether.
	*/

	leftOperand->CheckForZero();

	if ( leftOperand->IsZero() )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kUndefinedResult );

	if ( leftOperand->IsNegative() ||
		 (lz_compare( leftOperand, &one ) == kGreaterThan) )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		Two ways to the same answer, and the argument decides.

		Below kArcHyperSecantSeriesCeiling the function's own series
		is summed -- see lzi_archypersecant_series -- which skips the
		square root the other route pays for.

		At or above it, arcsech(x) = arccosh(1/x), guarded for the
		reason set out at kArcCoFunctionGuardDigits: as x climbs to
		one the reciprocal falls to one, and the hyperbolic
		arccosine's slope there is unbounded.

		The two agree to the digit everywhere both were run, so the
		threshold is a question of speed alone. See the ceiling's own
		note for where it came from.
	*/

	/*
		**And only up to sixty digits of width** -- since 2026-10-15,
		performance item 8 of that day's review. The series grows faster
		with the width than the route does: at 0.001 it led at 36 digits,
		74 us against 82, drew level at 60 and lost from 70 or 80, 1.7
		times at 600; base two the same, level near 60 digits.
	*/
	const long	kArcHyperSecantSeriesMostDigits = 60;

	CLZNumber	seriesCeiling( kArcHyperSecantSeriesCeiling );

	bool	b_useSeries = (lz_compare( leftOperand, &seriesCeiling ) != kGreaterThan) &&
						  (CLZNumber::EquivalentPlaces( glWorkingPrecision, leftOperand->m_lBase, 10 ) <= kArcHyperSecantSeriesMostDigits);

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase ) +
						  lzi_near_one_places( leftOperand );

	reciprocal.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	if ( b_useSeries )
	{
		lResult = lzi_first_failure( lResult,
						lzi_archypersecant_series( leftOperand, &result ) );
	}
	else
	{
		lResult = lzi_first_failure( lResult, lz_divide( &one, leftOperand, &reciprocal ) );
		lResult = lzi_first_failure( lResult, lz_archypercosine( &reciprocal, &result ) );
	}

	glWorkingPrecision = l_savedPrecision;

	// The guard places never leave the function...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_archypersine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_archypersine, kRetakeDoubt );

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	answer.DecideTies( lzi_inside_x_decide );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	/*
		arcsinh of an infinity is infinite whichever sign went in, so
		the missing sign costs nothing and the infinity carries -- the
		same reading the hyperbolic sine itself gets.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryStaysInfinite ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber magnitude;
	CLZNumber inner;
	CLZNumber result;

	bool b_negateResult = false;

	/*
		arcsinh(x) = ln( x + sqrt(x^2 + 1) ), over the whole real line
		and with no domain to refuse.
	*/

	one.Assign( 1L );

	/*
		**Guarded, as the arc co-functions are** -- see
		kArcCoFunctionGuardDigits. It had no guard, and inner came into
		the logarithm already rounded at the answer's own width: 8 of the
		audit's 76 cases missed correct rounding. lz_log_e's reduction
		running at its operand's width rather than its own had hidden
		some of it, the two errors cancelling by chance -- arcsinh(1.5)
		lies 0.036 of a unit past a rounding boundary -- and giving the
		logarithm its own width took the misses to ten. One guard place
		clears all 76 in base ten; four is the family's, and base two
		counts its places in bits.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase );

	inner.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	// It is odd, so work with the magnitude and put the sign back...
	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	if ( leftOperand->IsNegative() )
		b_negateResult = true;

	lResult = lzi_first_failure( lResult, lz_square( &magnitude, &inner ) );
	lResult = lzi_first_failure( lResult, lz_add( &inner, &one ) );
	lResult = lzi_first_failure( lResult, lz_root_2( &inner ) );
	lResult = lzi_first_failure( lResult, lz_add( &inner, &magnitude ) );
	lResult = lzi_first_failure( lResult, lz_log_e( &inner, &result ) );

	if ( b_negateResult )
		lResult = lzi_first_failure( lResult, lz_negate( &result ) );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_archypertangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_archypertangent, kRetakeDoubt );

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	answer.DecideTies( lzi_outside_x_decide );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	/*
		arctanh of an infinity is outside the domain whichever sign it
		had. A not-a-number, and before any test of the value.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber two;
	CLZNumber magnitude;
	CLZNumber numerator;
	CLZNumber denominator;
	CLZNumber quotient;
	CLZNumber result;

	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	two.SetBase( leftOperand->m_lBase );

	one.Assign( 1L );
	two.Assign( 2L );

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	/*
		Defined strictly inside [-1, 1], and the two ends are told
		apart from the outside.

		Beyond one there is no answer at all: kFailed with kfQNaN.

		At exactly one the quotient below is a division by zero and the
		limit is an infinity that cannot be written down, which is
		kUndefinedResult with kfInfinity -- the reading lz_log_e gives
		zero, and the one lzi_domain_failure documents for an argument
		on the edge of a domain rather than outside it.

		This used to refuse both the same way, with a flat kFailed, and
		its own comment said the limit was an infinity while the code
		reported otherwise. lz_archypercotangent mirrors this function
		exactly -- it is arctanh(1/x), so its edge at |x| = 1 is this
		one -- and was written to the documented rule, which left the
		pair disagreeing about the same pole. This is the side that was
		wrong.
	*/

	long l_order = lz_compare( &magnitude, &one );

	if ( l_order == kEqualTo )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kUndefinedResult );

	// Everything that is not strictly inside, which keeps kUnordered
	//	on the refusing side where the old single test had it...
	if ( l_order != kLessThan )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		arctanh(x) = ln( (1 + x) / (1 - x) ) / 2.
	*/

	/*
		**Guarded, as lz_archypersine and the arc co-functions are** --
		see kArcCoFunctionGuardDigits. With no guard the quotient came
		into the logarithm already rounded at the answer's own width,
		and 12 of the audit's 33 cases missed correct rounding, all by
		one unit, spread across the range. Two places clear all 33; four
		is the family's.

		The targets are widened by hand, as the arc co-functions do:
		lz_divide takes its width from the target it is handed, and these
		were built at the caller's precision before the guard went on, so
		without it the quotient would miss the guard. **On the audit's
		ladder that changes no value** -- put back, it alters only how
		two zeros and +-1e-20 are written -- so it is here for the
		principle, and costs nothing.

		The final rounding is not optional: the last step halves a
		logarithm in place, and half of a number rounded to 36 places can
		need a 37th. It used to hand that back unrounded.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase );

	numerator.SetPrecision( glWorkingPrecision );
	denominator.SetPrecision( glWorkingPrecision );
	quotient.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	/*
		**On |x|, the sign put back after**: arctanh(-x) = -arctanh(x).
		For x near -1 the quotient (1 + x)/(1 - x) is tiny -- 2e-15 at
		fourteen nines -- and at the working places it held that many
		fewer significant digits, which the logarithm then carried: until
		2026-09-26 negative arguments near -1 missed correct rounding by
		up to 680,000 units in base ten and 14 million in base twelve.
		On |x| the quotient is at least one.
	*/

	lResult = lzi_first_failure( lResult, lz_add( &one, &magnitude, &numerator ) );
	lResult = lzi_first_failure( lResult, lz_subtract( &one, &magnitude, &denominator ) );
	lResult = lzi_first_failure( lResult, lz_divide( &numerator, &denominator, &quotient ) );
	lResult = lzi_first_failure( lResult, lz_log_e( &quotient, &result ) );
	lResult = lzi_first_failure( lResult, lz_divide( &result, &two ) );

	if ( leftOperand->IsNegative() )
		lResult = lzi_first_failure( lResult, lz_negate( &result ) );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_arcsecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arcsecant, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// Exact at one, its nought spelled short -- until 2026-10-15 it was
	//	written to the full width, though nothing about it is in doubt...
	answer.ExactAndShort( lzi_exactly( leftOperand, 1 ) );

	// In another unit, exact where the angle is a part of a turn, and
	//	otherwise converted from radians -- see lzi_inverse_in_unit...
	if ( e_unit != AU_radians )
		return lzi_inverse_in_unit( kAngleArcSecant, leftOperand, targetNumber, answer, e_unit );

	// Before any test of the value, an infinity being held as zero...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber magnitude;
	CLZNumber reciprocal;
	CLZNumber result;

	// Set up before the domain test, which compares against one...
	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	// Defined only from one outwards, where the reciprocal lands
	//	inside the arccosine's domain. Nought is refused with the
	//	rest: the secant never returns a value below one...
	if ( lz_compare( &magnitude, &one ) == kLessThan )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		arcsec(x) = arccos(1/x), guarded for the reason set out at
		kArcCoFunctionGuardDigits.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcCoFunctionGuardDigits, leftOperand->m_lBase ) +
						  lzi_near_one_places( leftOperand );

	reciprocal.SetPrecision( glWorkingPrecision );
	result.SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( &one, leftOperand, &reciprocal ) );
	lResult = lzi_first_failure( lResult, lz_arccosine( &reciprocal, &result ) );

	glWorkingPrecision = l_savedPrecision;

	// The guard places never leave the function...
	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}

/*
long lz_arcsine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	long lResult = kNoError;
	long lCycle = 0;
	long lCycleValue = 0;

	CLZNumber* result = new CLZNumber;

	CLZNumber* two = new CLZNumber;
	CLZNumber* n_value = new CLZNumber;
	CLZNumber* two_n_value = new CLZNumber;
	CLZNumber* two_x_value = new CLZNumber;
	CLZNumber* cycle_result = new CLZNumber;

	// Set up...
	two->Assign( 2L );
	lz_multiply( two, leftOperand, two_x_value );

	// Calculate...
	for ( lCycle = 1; lCycle <= 32; lCycle++)
	{
		lCycleValue = 2 * lCycle;

		n_value->Assign( lCycle );
		two_n_value->Assign( lCycleValue );

		lz_combinations( two_n_value, n_value, cycle_result );
		lz_square( n_value );
		lz_multiply( n_value, cycle_result, cycle_result );
		lz_raise_n( two_x_value, two_n_value, two_n_value );
		lResult = lzi_first_failure( lResult, lz_divide( two_n_value, cycle_result, cycle_result ) );

		lz_add( result, cycle_result );
	}

	lResult = lzi_first_failure( lResult, lz_divide( result, two ) );
	lz_root_2( result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
	delete result;

	delete two;

	delete n_value;

	delete two_n_value;

	delete two_x_value;

	delete cycle_result;

	return lResult;
}
*/


// Both live further down, beside the sine/cosine kernel they ride on...
long lzi_arcsine_newton_floor( long lBase );
static void lzi_angle_by_newton( CLZNumber* pS, CLZNumber* pC, CLZNumber* pTarget );

long lz_arcsine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arcsine, kRetakeDoubtArcSine );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	//	In radians only: a decider reads the argument as the angle itself.
	if ( e_unit == AU_radians )
		answer.DecideTies( lzi_outside_x_decide );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// In another unit, exact where the angle is a part of a turn, and
	//	otherwise converted from radians -- see lzi_inverse_in_unit...
	if ( e_unit != AU_radians )
		return lzi_inverse_in_unit( kAngleArcSine, leftOperand, targetNumber, answer, e_unit );

	/*
		An infinity or a not-a-number has no arcsine: a not-a-number, as
		the arctangent answers. This used to be left to the flags spreading
		through the table's multiplications, which the computed series
		past the table does not do -- an infinity, held as nought, gave a
		first term of nought and came straight back out. Before anything
		reads the value...
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	long lCycle;
	long l_halvings = 0;

	/*
		The precision to answer in, and the one to put back on the way
		out. It is the engine's own rather than the operand's: every
		other operation builds its result in a number made fresh at the
		working precision, and an arcsine that instead honoured whatever
		its operand happened to carry would be the one function doing
		something different. Nesting composes, since this saves what it
		found rather than assuming kInitialPrecision.

		Widen before the working values are made, because that is the
		only moment their precision is settled: each takes it from the
		engine as it is made, and an operation writing to one hands it
		the precision of its own fresh result. See kArcSineGuardDigits
		for what the extra places buy.
	*/

	long l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision = l_savedPrecision + lzi_guard_places( kArcSineGuardDigits, leftOperand->m_lBase );

	CLZNumber* result = new CLZNumber;
	CLZNumber* sine = new CLZNumber;
	CLZNumber* cosine = new CLZNumber;
	CLZNumber* power = new CLZNumber;
	CLZNumber* squared = new CLZNumber;
	CLZNumber* cycle_result = new CLZNumber;
	CLZNumber* scratch = new CLZNumber;
	CLZNumber* one = new CLZNumber;
	CLZNumber* two = new CLZNumber;
	CLZNumber* threshold = new CLZNumber;

	bool b_negateResult = false;
	bool b_newton = false;

	// Set up, one and two in the operand's base -- see CLZAnswerBase.
	//	threshold stays base ten: Assign reads text in the number's own
	//	base, and a comparison converts it...
	one->SetBase( leftOperand->m_lBase );
	two->SetBase( leftOperand->m_lBase );

	one->Assign( 1L );
	two->Assign( 2L );
	threshold->Assign( kArcSineThreshold );

	sine->Assign( leftOperand );

	/*
		That assignment brought the caller's precision across with the
		value. It matters: the divide in the halving loop below has no
		target of its own and so works to the precision of this, its
		left operand. Put the guard digits back.
	*/

	sine->SetPrecision( glWorkingPrecision );

	// arcsin(-x) = -arcsin(x), so work with a positive argument...
	if ( sine->IsNegative() )
	{
		lz_negate( sine );
		b_negateResult = true;
	}

	/*
		Outside [-1, 1] there is no angle to find -- refused the way
		every other function refuses, the target a nought flagged
		not-a-number. This one alone used to return kFailed and write
		nothing, leaving the target holding whatever it held before, and
		in place leaving the operand as it was: the one refusal in the
		engine a caller could mistake for an answer.
	*/

	if ( lz_compare( sine, one ) == kGreaterThan )
	{
		lResult = lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );
		goto ArcSineExit;
	}

	/*
		This series converges far too slowly to be summed directly:
		its error is the first term dropped, and against an argument
		of .5 that leaves only sixteen correct digits.

		Halving the angle shrinks the argument geometrically, and the
		angle is recovered afterwards by doubling the result. Each
		halving costs a single square root:

			cos(t/2) = sqrt( (1 + cos t) / 2 )
			sin(t/2) = sin t / (2 cos(t/2))
	*/

	/*
		**cos t from an exact 1 - sin^2 t.** Near one, 1 - sin^2 t is tiny
		-- 2e-10 for ten nines -- and it was formed by squaring, rounding
		to the working places, and subtracting from one, which kept only
		that rounding's worth of digits: its root came out five digits
		short. Until 2026-09-26 half of all arguments within 1e-2 to 1e-15
		of one missed correct rounding, in base ten by up to 81,000 units
		and in hex by a million.

		It is formed as (1 - sin t)(1 + sin t) now, whose factors are
		exact -- addition in this engine loses nothing -- and the product
		is kept whole. **Multiplication rounds to the working places, not
		to its target's**, so it is the width that has to hold it, for
		the one multiply. The root of the exact value, taken at the
		working width, is then as good as the answer needs: widening the
		root too was tried and changed no answer at 36, 60 or 80 places.
	*/

	lz_subtract( one, sine, scratch );		// 1 - sin t
	lz_add( one, sine, cosine );			// 1 + sin t

	{
		long	l_width = glWorkingPrecision;

		glWorkingPrecision = cosine->m_lRadix + scratch->m_lRadix + 1;

		if ( glWorkingPrecision < l_width )
			glWorkingPrecision = l_width;

		lz_multiply( cosine, scratch );		// 1 - sin^2 t, exactly

		glWorkingPrecision = l_width;
	}

	lz_root_2( cosine );					// cos t

	/*
		**Past its base's floor, the angle of (sin t, cos t) by Newton on the
		sine/cosine kernel** in place of the halvings and the series below --
		stage 3 of performance item 1 of the review of 2026-10-15; see
		lzi_angle_by_newton. The pair is the one just formed, cos t from the
		exact (1 - sin t)(1 + sin t); there are no halvings to recover.
	*/

	b_newton = (glWorkingPrecision >= lzi_arcsine_newton_floor( leftOperand->m_lBase ));

	if ( b_newton )
		lzi_angle_by_newton( sine, cosine, result );

	while ( !b_newton && (lz_compare( sine, threshold ) == kGreaterThan) && (l_halvings < kArcSineMaxHalvings) )
	{
		lz_add( one, cosine, scratch );			// 1 + cos t
		lz_divide( scratch, two, cosine );		// (1 + cos t) / 2
		lz_root_2( cosine );					// cos(t/2)

		lz_multiply( cosine, two, scratch );	// 2 cos(t/2)
		lz_divide( sine, scratch );				// sin(t/2)

		l_halvings++;
	}

	/*
		**The table while it can serve the width, the recurrence past it.**
		gInvTrigA holds 22 coefficients of 50 digits, and against a sine
		under a tenth the term it drops is about 1e-48: so this, and
		everything built on it -- arccos, arctan, arccot, arcsec, arccsc
		-- stopped at 50 digits whatever was asked, until 2026-09-26. At
		60, 100 or 200 places every argument missed, right to 49 to 52
		places in base ten, 41 in hex.

		Past the table's reach the terms are computed, as the sine's and
		the exponential's are past theirs:

			p_k = p_(k-1) x^2 (2k - 1) / (2k),   term_k = p_k / (2k + 1)

		from p_0 = x, each step a multiply by the square and small whole
		numbers, until a term falls under the width. Within the reach
		nothing changes -- see kArcSineTableReach -- so every answer at
		the default is exactly as it was.
	*/

	if ( !b_newton )
	{
		result->Assign( sine );		// the '0' cycle...
		lz_square( sine, squared );
	}

	if ( b_newton )
		;
	else if ( CLZNumber::EquivalentPlaces( glWorkingPrecision, leftOperand->m_lBase, 10 ) <= kArcSineTableReach )
	{
		power->Assign( sine );

		for ( lCycle = 1; lCycle <= kInvTrigMax; lCycle++ )
		{
			lz_multiply( power, squared );
			lz_multiply( power, gInvTrigA[ lCycle ], cycle_result );
			lz_add( result, cycle_result );
		}
	}
	else
	{
		CLZNumber	factor;

		// The halvings leave a sine under a tenth -- pi/2 takes four of the
		//	eight -- so each term gains two digits at least, and the width
		//	in digits bounds the loop, which ends on a term under the width
		//	long before. lzi_series_bound is sized for 1/k!, which falls
		//	faster: its floor of 1024 terms runs short past about 2,050
		//	digits, where a sine near a tenth wants more...
		long		l_bound = CLZNumber::EquivalentPlaces( glWorkingPrecision, leftOperand->m_lBase, 10 ) + 8;

		factor.SetBase( leftOperand->m_lBase, false );

		power->Assign( sine );
		power->SetPrecision( glWorkingPrecision );
		cycle_result->SetPrecision( glWorkingPrecision );

		for ( lCycle = 1; lCycle <= l_bound; lCycle++ )
		{
			lz_multiply( power, squared );
			factor.Assign( (2 * lCycle) - 1 );
			lz_multiply( power, &factor );
			factor.Assign( 2 * lCycle );
			lz_divide( power, &factor );

			factor.Assign( (2 * lCycle) + 1 );
			lz_divide( power, &factor, cycle_result );

			cycle_result->CheckForZero();

			if ( cycle_result->IsZero() )
				break;

			lz_add( result, cycle_result );
		}
	}

	// Recover the angle that was halved away...
	while ( l_halvings-- > 0 )
		lz_multiply( result, two );

	if ( b_negateResult )
		lz_negate( result );

	// The guard digits were scaffolding: the answer goes back at the
	//	precision everything else here works to...
	result->SetPrecision( l_savedPrecision );
	result->Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

ArcSineExit:

	glWorkingPrecision = l_savedPrecision;

	// Clean up...
	delete result;
	delete sine;
	delete cosine;
	delete power;
	delete squared;
	delete cycle_result;
	delete scratch;
	delete one;
	delete two;
	delete threshold;

	return lResult;
}


/*
long lz_arcsine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	long lResult = kNoError;
	long lCycle = 0;
	long lCycleValue = 0;

    auto result = new CLZNumber;

    auto multiplier = new CLZNumber;
    auto numerator = new CLZNumber;
    auto denominator = new CLZNumber;
    auto power = new CLZNumber;
    auto x_squared = new CLZNumber;
    auto cycle_result = new CLZNumber;

	// Set up...
	result->Assign( leftOperand );

	multiplier->Assign( 1L );
	power->Assign( leftOperand );
	lz_square( leftOperand, x_squared );

	// Calculate...
	for ( lCycle = 1; lCycle <= 256; lCycle++)
	{
		lCycleValue = 2 * lCycle + 1;
		cycle_result->Assign( lCycleValue );
		lCycleValue--;
		denominator->Assign( lCycleValue );
		lCycleValue--;
		numerator->Assign( lCycleValue );

		lz_multiply( multiplier, numerator );
		lResult = lzi_first_failure( lResult, lz_divide( multiplier, denominator ) );
		lz_multiply( power, x_squared );
		lResult = lzi_first_failure( lResult, lz_divide( power, cycle_result, cycle_result ) );
		lz_multiply( cycle_result, multiplier );

		lz_add( result, cycle_result );
	}

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete multiplier;
    delete numerator;
    delete denominator;
    delete power;
    delete x_squared;
    delete cycle_result;

	return lResult;
}
*/


long lz_arctangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arctangent, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	//	In radians only: a decider reads the argument as the angle itself.
	if ( e_unit == AU_radians )
		answer.DecideTies( lzi_inside_x_decide );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// In another unit, exact where the angle is a part of a turn, and
	//	otherwise converted from radians -- see lzi_inverse_in_unit...
	if ( e_unit != AU_radians )
		return lzi_inverse_in_unit( kAngleArcTangent, leftOperand, targetNumber, answer, e_unit );

	/*
		arctan(+inf) is pi/2 and arctan(-inf) is -pi/2. An infinity is
		held here as zero with no sign to tell those apart, so this
		gives a not-a-number rather than guessing between them, by the
		rule the exponential and the hyperbolic tangent already follow.

		Before any test of the value, which would read the infinity as
		nought.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	/*
		**Worked wider than the answer, rounded back once at the end** --
		see kArcTanGuardDigits. Before any working value exists, a
		CLZNumber taking its precision at construction.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kArcTanGuardDigits, leftOperand->m_lBase );

	CLZNumber one;
	CLZNumber magnitude;
	CLZNumber argument;
	CLZNumber divisor;
	CLZNumber quotient;
	CLZNumber result;

	bool b_reflected = false;

	/*
		arctan(x) = arcsin( x / sqrt(1 + x^2) ).

		Through the arcsine rather than the arctangent's own series,
		which diverges past an argument of 1 and would need a reduction
		written to go with it. There is no domain to refuse: the
		quotient lies strictly inside (-1, 1) for every real x.

		Past |x| = 1 the quotient crowds up against 1 though, and that
		is where the arcsine is worst conditioned: its slope there is
		1/sqrt(1 - z^2), which climbs without bound, so an error in the
		quotient arrives in the answer multiplied by roughly x. Measured
		rather than argued -- arctan of a million came back with 33
		digits where the same code gives 37 either side of 1.

		So reflect instead, and keep the argument where the arcsine is
		well behaved:

			arctan(x) = sign(x) * pi/2 - arctan(1/x)

		The slope over [-1, 1] is at most sqrt(2). Nothing cancels in
		that subtraction: for large x the second term is small against
		pi/2, and just past 1 both terms are near pi/4 and so is the
		answer.
	*/

	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	argument.Assign( leftOperand );

	/*
		**The reflection divides into argument, and lz_divide takes its
		width from its target** -- which the Assign had just given the
		operand's precision. So 1/x was formed to whatever the operand
		carried: 20 places from a 20-place operand, 2e13 units out of the
		36th place. The working precision, as everything else here uses.
	*/

	argument.SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	if ( lz_compare( &magnitude, &one ) == kGreaterThan )
	{
		lResult = lzi_first_failure( lResult, lz_divide( &one, leftOperand, &argument ) );
		b_reflected = true;
	}

	lResult = lzi_first_failure( lResult, lz_square( &argument, &divisor ) );
	lResult = lzi_first_failure( lResult, lz_add( &divisor, &one ) );
	lResult = lzi_first_failure( lResult, lz_root_2( &divisor ) );
	lResult = lzi_first_failure( lResult, lz_divide( &argument, &divisor, &quotient ) );
	lResult = lzi_first_failure( lResult, lz_arcsine( &quotient, &result ) );

	// Undo the reflection, taking pi/2's sign from the operand's...
	if ( b_reflected )
	{
		CLZNumber	halfPi( kPiOverTwo );
		CLZNumber	reflected;

		// It is on the left, so in the operand's base -- see
		//	CLZAnswerBase...
		lzi_pi_family_in_base( &halfPi, kPiOverTwo, leftOperand->m_lBase, glWorkingPrecision );

		if ( leftOperand->IsNegative() )
			lResult = lzi_first_failure( lResult, lz_negate( &halfPi ) );

		lResult = lzi_first_failure( lResult, lz_subtract( &halfPi, &result, &reflected ) );

		result.Assign( &reflected );
	}

	/*
		The guard places never leave the function, and neither does the
		tail of pi/2: kPiOverTwo is carried to 120 places, so a reflected
		difference comes out wearing places it has not earned. Worth
		keeping the wide constant for the subtraction and rounding after,
		rather than narrowing it first: those places are guard digits
		while the arithmetic is happening, and one subtraction against a
		wide constant measures at 1.00x. lz_sine and lz_cosine do narrow
		theirs, because three comparisons and an addition a call did not.
	*/

	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


/*
	**A tiny angle on a half, decided by its side of the quotient** -- the
	arctangent's lzi_inside_x_decide, asked of y / x. Where x is positive
	and |y| at most |x| the answer is arctan( y / x ), which lies inside
	y / x; where that quotient ends exactly on a half of the answer's
	places, only the cubic term moves the answer off it, past the
	retake's cap -- lz_arctangent2( 1.5e-36, 1 ) at 36 places is
	arctan( 1.5e-36 ) and wants ...1, not the ...2 a half read exact and
	taken to the even digit would give. A quotient that does not end at
	those places cannot lie on a half, and declines with everything else,
	so the retake goes on as before.
*/

static bool lzi_arctangent2_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* pRight, long lPlaces, long /* lMostDigits */ )
{
	if ( (pRight == nullptr) || (lPlaces < 0) || (pLeft->m_lBase != pRight->m_lBase) )
		return false;

	if ( pLeft->FlagIsOn( kfInfinity ) || pLeft->FlagIsOn( kfQNaN ) ||
		 pRight->FlagIsOn( kfInfinity ) || pRight->FlagIsOn( kfQNaN ) )
		return false;

	pRight->CheckForZero();

	if ( pRight->IsZero() || pRight->IsNegative() )
		return false;

	long	l_saved = glWorkingPrecision;

	// On a half, the quotient ends a place past the answer's...
	glWorkingPrecision = lPlaces + 2;

	CLZNumber	ratio;
	CLZNumber	back;

	lz_divide( pLeft, pRight, &ratio );

	// ...and is exact only if it gives y back...
	glWorkingPrecision = ratio.m_lRadix + pRight->m_lRadix + 2;
	lz_multiply( &ratio, pRight, &back );
	glWorkingPrecision = l_saved;

	if ( lz_compare( &back, pLeft ) != kEqualTo )
		return false;

	return lzi_beside_x_decide( pResult, &ratio, lPlaces, -1 );
}

/*
	**atan2 answering in a unit** -- the host's own call, from
	lz_arctangent2 once the origin and the positive real axis have
	answered nought. The operands are its copies, both in the left base.

	**The axes and the diagonals are exact**: the negative real axis half
	a turn, the imaginary axis a quarter, |y| = |x| an eighth or three --
	each divided out of the turn, short where it ends, with y's sign.
	Nothing else is a rational part of a turn, |y / x| being rational.
	**The rest is the radian angle**, kAngleUnitGuardDigits and three
	wider, times the turn over 2pi, rounded once -- as
	lzi_inverse_in_unit takes the other inverses.
*/

static long lzi_arctangent2_in_unit( CLZNumber* pY, CLZNumber* pX, CLZNumber* pTarget, CLZAnswerBase& answer,
									 AngleUnit eUnit )
{
	long	l_base = pY->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_turn = lzi_turn_of( eUnit );
	long	l_a = 0;
	long	l_b = 0;

	CLZNumber	magnitudeY;
	CLZNumber	magnitudeX;

	lz_abs( pY, &magnitudeY );
	lz_abs( pX, &magnitudeX );

	if ( pY->IsZero() )
	{
		l_a = 1;
		l_b = 2;
	}
	else if ( pX->IsZero() )
	{
		l_a = 1;
		l_b = 4;
	}
	else if ( lz_compare( &magnitudeY, &magnitudeX ) == kEqualTo )
	{
		l_a = pX->IsNegative() ? 3 : 1;
		l_b = 8;
	}

	if ( l_b != 0 )
	{
		CLZNumber	value;

		lzi_part_of_turn( l_turn, pY->IsNegative() ? -l_a : l_a, l_b, l_base, &value );
		answer.Exact();

		if ( !(pTarget->IsConstant()) )
			pTarget->Assign( &value );

		return kNoError;
	}

	glWorkingPrecision = l_saved + lzi_guard_places( kAngleUnitGuardDigits + 3, l_base );

	CLZNumber	theta;
	CLZNumber	twoPi( kTwoPi );
	CLZNumber	scaled;
	CLZNumber	angle;

	long	lResult = lz_arctangent2( pY, pX, &theta );

	lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, glWorkingPrecision + 2 );
	lzi_times_whole( &theta, l_turn, &scaled );
	lz_divide( &scaled, &twoPi, &angle );

	glWorkingPrecision = l_saved;

	angle.SetPrecision( glWorkingPrecision );
	angle.Round();

	if ( !(pTarget->IsConstant()) )
		pTarget->Assign( &angle );

	return lResult;
}

/*
	**The angle of the point (x, y)** -- atan2, y first as C's is: the
	two-argument arctangent, in (-pi, pi], its quadrant taken from both
	signs. Stage 2a of the complex layer, 2026-10-22, and all of lzc_arg's
	work. **Radians only for now**: the angle unit is stage 2b's, and comes
	to every function that takes or answers an angle at once -- Rich's
	call, 2026-10-22, so no release has a unit some of them ignore.

	**On the axes it answers by the convention std::atan2 and std::arg keep
	for +0**, this engine having no negative nought: 0 along the positive
	real axis, pi along the negative one -- the side the real-to-complex
	switch arrives from, Rich's decision of 2026-10-21 -- and +-pi/2 on
	the imaginary axis. **The origin answers 0**, Rich's call, 2026-10-22:
	no angle is right there, but std::atan2( +0, +0 ) and std::arg( 0 )
	answer 0, and so has Elzed's arg(0) always.

	**An infinity or a not-a-number in either operand answers a
	not-a-number**: an infinity here carries no sign, so it names no
	direction -- the arctangent's own reading of one.

	**Worked kArcTan2GuardDigits wider and rounded once.** The smaller
	magnitude over the larger keeps the arctangent on [0, 1], where its
	slope is at most one, and never forms a quotient past one:

		|y| <= |x|:	t = arctan( |y| / |x| )
		|y| >  |x|:	t = pi/2 - arctan( |x| / |y| )

	then pi - t where x is negative, and y's sign. Each step rounds at the
	guarded width -- the quotient, the arctangent's own answer, pi or pi/2
	-- so the error before the last rounding is under two units of the
	guard's last place; nothing cancels, places being counted from the
	radix mark (see lz_arccotangent).

	**A right operand in another base is brought over exactly.** Where its
	base ends in the left's it converts exactly. Where it does not, it is a
	whole number X over b^k, and the angle of (X / b^k, y) is the angle of
	(X, y b^k) -- scaling both by the same positive number moves no angle
	-- so both are worked in the left base with no rounding added. The
	answer is in the left operand's base, as every function's is.
*/

long lz_arctangent2( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// Both are needed, and a refusal still answers -- see lzi_domain_failure...
	if ( (leftOperand == nullptr) || (rightOperand == nullptr) )
	{
		if ( (targetNumber != nullptr) && !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( 0L );
			targetNumber->SetFlag( kfQNaN );
		}

		return kFailed;
	}

	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_arctangent2, rightOperand, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and a tiny angle on a half decided by its side of y / x -- see
	//	lzi_arctangent2_decide. In radians only, where y / x is the angle...
	if ( e_unit == AU_radians )
		answer.DecideTies( lzi_arctangent2_decide );

	CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( leftOperand->FlagIsOn( kfQNaN ) || leftOperand->FlagIsOn( kfInfinity ) ||
		 rightOperand->FlagIsOn( kfQNaN ) || rightOperand->FlagIsOn( kfInfinity ) )
	{
		if ( !(p_target->IsConstant()) )
		{
			p_target->Assign( 0L );
			p_target->SetFlag( kfQNaN );

			// Where the trouble started travels with it...
			if ( leftOperand->FlagIsOn( kfDiv0 ) || rightOperand->FlagIsOn( kfDiv0 ) )
				p_target->SetFlag( kfDiv0 );
		}

		return kNoError;
	}

	long	l_base = leftOperand->m_lBase;
	long	l_savedPrecision = glWorkingPrecision;
	long	lResult = kNoError;

	// Copies, so an answer written over an operand cannot reach the work...
	CLZNumber	y;
	CLZNumber	x;

	y.Assign( leftOperand );
	x.Assign( rightOperand );

	// ...and x in y's base, exactly, scaling y with it where it must...
	if ( x.m_lBase != l_base )
	{
		x.CheckForZero();
		x.TrimTrailing();

		if ( x.IsZero() || lzi_ends_in_base( x.m_lBase, l_base ) )
			lzi_into_base_keeping( &x, l_base, 0 );
		else
		{
			long		l_places = x.m_lRadix;
			CLZNumber	scale;

			lzi_scale_by_base_power( &x, -l_places );
			scale.SetBase( x.m_lBase, false );
			scale.Assign( 1L );
			lzi_scale_by_base_power( &scale, -l_places );
			x.SetBase( l_base, true );
			scale.SetBase( l_base, true );

			// A whole number adds no places, so this product is exact...
			glWorkingPrecision = ((y.m_lRadix > glWorkingPrecision) ? y.m_lRadix : glWorkingPrecision) + 2;
			lResult = lzi_first_failure( lResult, lz_multiply( &y, &scale ) );
			glWorkingPrecision = l_savedPrecision;
		}
	}

	y.CheckForZero();
	x.CheckForZero();

	// The origin and the positive real axis answer nought, exactly...
	if ( y.IsZero() && !(x.IsNegative()) )
	{
		CLZNumber	nought;

		nought.SetBase( l_base, false );
		nought.Assign( 0L );
		answer.ExactAndShort();

		if ( !(p_target->IsConstant()) )
			p_target->Assign( &nought );

		return lResult;
	}

	// In another unit, the axes and diagonals exact and the rest converted
	//	from radians -- see lzi_arctangent2_in_unit...
	if ( e_unit != AU_radians )
		return lzi_first_failure( lResult, lzi_arctangent2_in_unit( &y, &x, p_target, answer, e_unit ) );

	/*
		**Worked wider than the answer, rounded back once at the end** --
		see kArcTan2GuardDigits. Before any working value exists, a
		CLZNumber taking its precision at construction.
	*/

	glWorkingPrecision += lzi_guard_places( kArcTan2GuardDigits, l_base );

	CLZNumber	magnitudeY;
	CLZNumber	magnitudeX;
	CLZNumber	ratio;
	CLZNumber	angle;
	CLZNumber	result;

	result.SetBase( l_base, false );

	if ( y.IsZero() )
	{
		// The negative real axis: pi, from the +0 side...
		CLZNumber	pi( kPi );

		lzi_pi_family_in_base( &pi, kPi, l_base, glWorkingPrecision );
		result.Assign( &pi );
	}
	else if ( x.IsZero() )
	{
		// The imaginary axis: pi/2, with y's sign below...
		CLZNumber	halfPi( kPiOverTwo );

		lzi_pi_family_in_base( &halfPi, kPiOverTwo, l_base, glWorkingPrecision );
		result.Assign( &halfPi );
	}
	else
	{
		lResult = lzi_first_failure( lResult, lz_abs( &y, &magnitudeY ) );
		lResult = lzi_first_failure( lResult, lz_abs( &x, &magnitudeX ) );

		// The smaller over the larger, so the arctangent stays on [0, 1]...
		bool	b_steep = (lz_compare( &magnitudeY, &magnitudeX ) == kGreaterThan);

		if ( b_steep )
			lResult = lzi_first_failure( lResult, lz_divide( &magnitudeX, &magnitudeY, &ratio ) );
		else
			lResult = lzi_first_failure( lResult, lz_divide( &magnitudeY, &magnitudeX, &ratio ) );

		lResult = lzi_first_failure( lResult, lz_arctangent( &ratio, &angle ) );

		if ( b_steep )
		{
			CLZNumber	halfPi( kPiOverTwo );

			lzi_pi_family_in_base( &halfPi, kPiOverTwo, l_base, glWorkingPrecision );
			lResult = lzi_first_failure( lResult, lz_subtract( &halfPi, &angle, &result ) );
		}
		else
			result.Assign( &angle );

		// ...into the left half-plane where x is negative...
		if ( x.IsNegative() )
		{
			CLZNumber	pi( kPi );
			CLZNumber	reflected;

			lzi_pi_family_in_base( &pi, kPi, l_base, glWorkingPrecision );
			lResult = lzi_first_failure( lResult, lz_subtract( &pi, &result, &reflected ) );
			result.Assign( &reflected );
		}
	}

	// ...and below the real axis where y is...
	if ( y.IsNegative() )
		lResult = lzi_first_failure( lResult, lz_negate( &result ) );

	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	if ( !(p_target->IsConstant()) )
		p_target->Assign( &result );

	return lResult;
}


/*
	**The ways to choose k of n, on the whole numbers, and nought past n**
	-- Rich's call, 2026-10-11, as factorial is defined: a fraction or a
	negative is refused, with kFailed and a not-a-number. It used to be
	n! / (k! (n - k)!) on anything, checking none of the three factorials
	it took, so a fraction, a negative, or k past n each came back a
	not-a-number reported as success; and a large n either never finished
	or, past the answer limit, was refused into the same not-a-number.
	C(100000, 3) answered NaN.

	**Now a product over the smaller of k and n - k**, in one of two
	ways, sized first as every function is -- see lzi_answer_size:

	- **By its primes**, where n is small enough beside the count for the
	  sieve to cost less than the product: Legendre's formula gives each prime's exponent in C(n, k)
	  without dividing anything, every prime power is at most n, so it
	  fits in a long, and their product, multiplied out as a balanced
	  tree, holds nothing larger than the answer.
	- **Otherwise the k whole numbers that end at n, over k!**, both
	  multiplied out as balanced trees and divided once, exactly. With a
	  huge n the count is small, and so is that division.

	The second alone took 7.5 s for C(332000, 166000), 99,940 digits, all
	but a second of it dividing a 900,000-digit numerator by a
	794,000-digit k!. Stepping C(n, i) = C(n, i - 1) (n - k + i) / i
	instead keeps every number under the answer and took 6.6 s for
	C(100000, 50000), multiplying a long number by a short one 50,000
	times.
*/

/*
	**Which way is cheaper, measured 2026-10-11.** The sieve and its
	product cost about 0.0075 us per unit of n; the product and division
	about 0.05 count^1.56 us, whatever n is. So the primes are taken
	while n is under seven times count^1.56 -- C(1000, 10), C(100000,
	500), C(1E7, 7000) are about where they cross -- and only within the
	sieve's reach, a byte of memory per unit of n...
*/

const long kCombinationsSieveReach = 10000000;

static bool lzi_binomial_by_primes_is_cheaper( long lN, long lCount )
{
	return (lN <= kCombinationsSieveReach) &&
		   (static_cast<double>( lN ) < (7.0 * std::pow( static_cast<double>( lCount ), 1.56 )));
}

// The power of the prime p dividing C(n, k), by Legendre's formula: the
//	sum over powers of p of n/p^i - k/p^i - (n - k)/p^i, whole quotients.
//	Never more than one power of p past n, so p^e fits where n does...
static long lzi_binomial_prime_power( long lN, long lK, long lP )
{
	long	l_power = 1;
	long	l_n = lN;
	long	l_k = lK;
	long	l_rest = lN - lK;

	while ( l_n > 0 )
	{
		l_n /= lP;
		l_k /= lP;
		l_rest /= lP;

		if ( (l_n - l_k - l_rest) > 0 )
			l_power *= lP;
	}

	return l_power;
}

// The product of factors[ lFrom ] up to factors[ lTo ], in lBase, split in
//	halves as lzi_range_product is...
static void lzi_factor_product( const std::vector<long>& factors, long lFrom, long lTo, long lBase,
								CLZNumber* pProduct )
{
	if ( (lTo - lFrom) < 16 )
	{
		CLZNumber	term;

		term.SetBase( lBase );
		pProduct->SetBase( lBase );
		pProduct->Assign( factors[ lFrom ] );

		for ( long l = lFrom + 1; l <= lTo; l++ )
		{
			term.Assign( factors[ l ] );
			lz_multiply( pProduct, &term );
		}

		return;
	}

	long		l_middle = lFrom + ((lTo - lFrom) / 2);
	CLZNumber	lower;
	CLZNumber	upper;

	lzi_factor_product( factors, lFrom, l_middle, lBase, &lower );
	lzi_factor_product( factors, l_middle + 1, lTo, lBase, &upper );
	lz_multiply( &lower, &upper, pProduct );
}

// C(n, k) by its primes, for n within the sieve's reach...
static void lzi_binomial_by_primes( long lN, long lK, long lBase, CLZNumber* pResult )
{
	std::vector<bool>	composite( static_cast<size_t>( lN ) + 1, false );
	std::vector<long>	factors;

	for ( long p = 2; p <= lN; p++ )
	{
		if ( composite[ static_cast<size_t>( p ) ] )
			continue;

		for ( long m = p * 2; m <= lN; m += p )
			composite[ static_cast<size_t>( m ) ] = true;

		long	l_power = lzi_binomial_prime_power( lN, lK, p );

		if ( l_power > 1 )
			factors.push_back( l_power );
	}

	pResult->SetBase( lBase );

	if ( factors.empty() )
		pResult->Assign( 1L );
	else
		lzi_factor_product( factors, 0, static_cast<long>( factors.size() ) - 1, lBase, pResult );
}

// A whole number's value as a long, exactly, answering false if it does not
//	fit. The size estimate will not do for a count: it is exp(ln k), and
//	took 500 for 499.99999999999994...
static bool lzi_whole_as_long( CLZNumber* pValue, long* plValue )
{
	long	l_value = 0;
	long	l_base = pValue->m_lBase;

	for ( long i = pValue->m_lDigitsLength - 1; i >= pValue->m_lRadix; i-- )
	{
		long	l_digit = pValue->m_lpDigits[ i ];

		if ( l_value > ((LONG_MAX - l_digit) / l_base) )
			return false;

		l_value = (l_value * l_base) + l_digit;
	}

	*plValue = l_value;

	return true;
}

// The product of start + lFrom up to start + lTo, whole numbers in start's
//	base, split in halves until a run is short enough to multiply along...
static void lzi_range_product( CLZNumber* pStart, long lFrom, long lTo, CLZNumber* pProduct )
{
	if ( (lTo - lFrom) < 16 )
	{
		CLZNumber	offset;
		CLZNumber	term;

		offset.SetBase( pStart->m_lBase );
		offset.Assign( lFrom );
		lz_add( pStart, &offset, pProduct );

		for ( long l = lFrom + 1; l <= lTo; l++ )
		{
			offset.Assign( l );
			lz_add( pStart, &offset, &term );
			lz_multiply( pProduct, &term );
		}

		return;
	}

	long		l_middle = lFrom + ((lTo - lFrom) / 2);
	CLZNumber	lower;
	CLZNumber	upper;

	lzi_range_product( pStart, lFrom, l_middle, &lower );
	lzi_range_product( pStart, l_middle + 1, lTo, &upper );
	lz_multiply( &lower, &upper, pProduct );
}

// ln P(n, k), the k whole numbers that end at n: Stirling's leading terms
//	for ln Gamma(n + 1) - ln Gamma(n - k + 1), taken as one difference so a
//	huge n and a small k do not cancel. dLower is n - k + 1, which the
//	caller supplies: from the estimates of n and k only when k is at most
//	half of n, since two near estimates may cross. For sizing an answer
//	and nothing else.
//
//	Under ten the lower term is lzi_ln_gamma_estimate's, shifted as it
//	shifts: the leading terms alone read ln Gamma(1) as -0.08, which at
//	k = n refused P(1516, 1516) in base 11 a thirtieth of a place under
//	the limit. Nothing cancels there -- n - k is small -- so the two terms
//	are taken apart...
static double lzi_ln_falling_estimate( double dN, double dLnN, double dK, double dLower )
{
	if ( !std::isfinite( dN ) )
		return dK * dLnN;

	if ( dLower < 10.0 )
		return lzi_ln_gamma_estimate( dN + 1.0 ) - lzi_ln_gamma_estimate( dLower );

	return ((dLower - 0.5) * std::log1p( dK / dLower )) + (dK * std::log( dN + 1.0 )) - dK;
}

// ln C(n, k), k the smaller of k and n - k: ln P(n, k) less
//	ln Gamma(k + 1)...
static double lzi_ln_binomial_estimate( double dN, double dLnN, double dK )
{
	if ( dK < 1.0 )
		return 0.0;

	if ( !std::isfinite( dK ) )
		return HUGE_VAL;

	return lzi_ln_falling_estimate( dN, dLnN, dK, dN - dK + 1.0 ) - lzi_ln_gamma_estimate( dK + 1.0 );
}

long lz_combinations( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	// An infinity or a not-a-number in either is a not-a-number...
	if ( lzi_propagate_as_nan( leftOperand, rightOperand, targetNumber ) )
		return kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	leftOperand->CheckForZero();
	rightOperand->CheckForZero();

	if ( leftOperand->IsNegative() || rightOperand->IsNegative() ||
		 !(leftOperand->IsInteger()) || !(rightOperand->IsInteger()) )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	// More than there are can be chosen no ways at all...
	if ( rightOperand->Compare( leftOperand ) == kGreaterThan )
	{
		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		return kNoError;
	}

	// C(n, k) is C(n, n - k): the product runs over the smaller...
	CLZNumber	count;
	CLZNumber	factor;

	lz_subtract( leftOperand, rightOperand, &factor );

	if ( factor.Compare( rightOperand ) == kLessThan )
	{
		count.Assign( &factor );
		factor.Assign( rightOperand );
	}
	else
		count.Assign( rightOperand );

	// ...sized before it is begun -- see lzi_answer_size...
	double	d_n = 0.0;
	double	d_lnN = 0.0;
	double	d_k = 0.0;
	long	l_beyond = kNoError;

	if ( lzi_value_or_infinity( leftOperand, &d_n ) && lzi_value_or_infinity( &count, &d_k ) &&
		 lzi_ln_estimate( leftOperand, &d_lnN ) &&
		 lzi_answer_beyond( lzi_answer_size( lzi_ln_binomial_estimate( d_n, d_lnN, d_k ), leftOperand->m_lBase ),
							leftOperand, targetNumber, &l_beyond ) )
		return l_beyond;

	// ...and taken in the operand's base -- see CLZAnswerBase. Under the
	//	limit the count is far inside a long; one that is not was past it...
	CLZNumber	result;
	long		l_count = 0;

	if ( !lzi_whole_as_long( &count, &l_count ) )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kTooLarge );

	result.SetBase( leftOperand->m_lBase );
	result.Assign( 1L );

	long	l_n = 0;

	if ( (l_count > 0) && lzi_whole_as_long( leftOperand, &l_n ) && lzi_binomial_by_primes_is_cheaper( l_n, l_count ) )
		lzi_binomial_by_primes( l_n, l_count, leftOperand->m_lBase, &result );
	else if ( l_count > 0 )
	{
		CLZNumber	zero;
		CLZNumber	falling;
		CLZNumber	factorial;

		zero.SetBase( leftOperand->m_lBase );
		zero.Assign( 0L );

		// (n - count + 1) ... n, over 1 ... count...
		lzi_range_product( &factor, 1, l_count, &falling );
		lzi_range_product( &zero, 1, l_count, &factorial );
		lz_divide_integer( &falling, &factorial, &result );
	}

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( &result );
		p_target->TrimAll();
		p_target->Round();
	}

	return kNoError;
}


long lz_compare( CLZNumber* leftOperand, CLZNumber* rightOperand )
{
	if ( lzi_is_unordered( leftOperand, rightOperand ) )
		return kUnordered;

	long lResult = kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Compare();
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();


	lResult = leftOperand->Compare( rightOperand );

	return lResult;
}


long lz_cosecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_cosecant, kRetakeDoubtCosecant );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// In another unit, the angle reduced exactly first -- see
	//	lzi_circular_in_unit...
	if ( e_unit != AU_radians )
		return lzi_circular_in_unit( kAngleCosecant, leftOperand, targetNumber, answer, e_unit );

	long lResult = kNoError;

	// Guarded across the division -- see kTrigQuotientGuardDigits, and the
	//	note in lz_cotangent...
	long l_savedPrecision = glWorkingPrecision;

	long	l_guard = lzi_guard_places( kTrigQuotientGuardDigits, leftOperand->m_lBase );

	glWorkingPrecision += l_guard;

    auto result = new CLZNumber;
    auto one = new CLZNumber;
    auto sine_result = new CLZNumber;

	// Calculate...
	// In the operand's base -- see CLZAnswerBase...
	one->SetBase( leftOperand->m_lBase );
	one->Assign( 1L );
	/*
		**Wider by twice the answer's size where the sine is small** --
		see lzi_small_divisor_widening. Near nought the argument says
		how small, before the sine is taken. Near any other multiple of
		pi only the sine can, so it is taken, and taken again wider when
		its size would spend more than half the guard -- which leaves
		every ordinary argument one sine.

		**Half, not all of it.** Until 2026-09-26 this waited for the
		size to outgrow the whole guard, the threshold exp and e^n had
		too, and a divisor between a hundredth and a tenth asked for
		exactly the guard and was refused it: 3 to 8% of such arguments
		missed correct rounding in base ten, by up to 0.88 of a unit, in
		all four quotients. The guard is in the operand's places, as
		every guard is -- see lzi_guard_places.
	*/

	long	l_widened = lzi_small_divisor_widening( leftOperand );

	glWorkingPrecision += l_widened;

	lz_sine( leftOperand, sine_result );

	// A sine too small to show at this width is taken wider until it
	//	does -- see lzi_widen_until_seen...
	l_widened += lzi_widen_until_seen( lz_sine, leftOperand, sine_result );

	long	l_wanted = lzi_small_divisor_widening( sine_result );

	if ( (l_wanted - l_widened) > (l_guard / 2) )
	{
		glWorkingPrecision += l_wanted - l_widened;

		lz_sine( leftOperand, sine_result );
	}

	// lz_divide takes its width from its target, made before any of
	//	this...
	result->SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( one, sine_result, result ) );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete one;
    delete sine_result;

	return lResult;
}


long lz_cosine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_cosine, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and in a unit an angle a hair from one whose cosine is a half
	//	decided by its side of it -- see lzi_unit_half_decide...
	if ( e_unit != AU_radians )
		answer.DecideTies( lzi_cosine_unit_decide );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// In another unit, the angle reduced exactly first -- see
	//	lzi_circular_in_unit...
	if ( e_unit != AU_radians )
		return lzi_circular_in_unit( kAngleCosine, leftOperand, targetNumber, answer, e_unit );

	long lResult = kNoError;

	/*
		Guarded in its own right and not only through lz_sine -- see
		kTrigGuardDigits. The addition below is where the cosine would
		otherwise lose its place: x + pi/2 rounds to the working
		precision before the sine ever sees it, and no width inside the
		sine puts back what was dropped on the way in.
	*/

	long l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kTrigGuardDigits, leftOperand->m_lBase );

    auto result = new CLZNumber;
    auto piOverTwo = new CLZNumber( kPiOverTwo);

	/*
		To the working width before the addition, for the reason given
		at kTwoPi in lz_sine. kPiOverTwo carries 120 places for the
		series that runs underneath; the addition here rounds to the
		working precision regardless -- the comment above says as much
		-- so the places past it are paid for and dropped. Leaving them
		on made the cosine 1.61 times its cost, and everything reached
		through it, which is the tangent, cotangent and secant.

		No test catches its removal -- it is speed, not accuracy. See
		the same note at the 2*pi narrowing in lz_sine.
	*/

	lzi_pi_family_in_base( piOverTwo, kPiOverTwo, leftOperand->m_lBase, glWorkingPrecision );

	if ( piOverTwo->m_lRadix > glWorkingPrecision )
	{
		piOverTwo->SetPrecision( glWorkingPrecision );
		piOverTwo->Round();
	}

	// Calculate...
    lResult = lzi_first_failure( lResult, lz_add(leftOperand, piOverTwo, result) );
    lResult = lzi_first_failure( lResult, lz_sine(result) );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete piOverTwo;

	return lResult;
}


long lz_cotangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_cotangent, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// In another unit, the angle reduced exactly first -- see
	//	lzi_circular_in_unit...
	if ( e_unit != AU_radians )
		return lzi_circular_in_unit( kAngleCotangent, leftOperand, targetNumber, answer, e_unit );

	long lResult = kNoError;

	/*
		Guarded across the division, not merely inside the sine and
		cosine being divided -- see kTrigQuotientGuardDigits.

		Widening those two alone is not enough here, and can cost a
		digit rather than buy one. Each comes back correctly rounded to
		the working precision, and two roundings that had been
		cancelling in the quotient by luck stop doing so. Three sweep
		cases dropped exactly that way -- cot at 1e8, -1e8 and
		.9999999999 -- until the width was carried across the divide
		as well.
	*/

	long l_savedPrecision = glWorkingPrecision;

	long	l_guard = lzi_guard_places( kTrigQuotientGuardDigits, leftOperand->m_lBase );

	glWorkingPrecision += l_guard;

    auto result = new CLZNumber;
    auto sine_result = new CLZNumber;
    auto cosine_result = new CLZNumber;

	// Calculate...
	// Wider by twice the answer's size where the sine is small, as the
	//	cosecant does and for the same reasons -- see there...
	/*
		**One fold and one kernel run past the kernel's floor** -- see
		lzi_sine_and_cosine, and the same in lz_tangent -- wherever the
		widening below would take only the argument's own: a sine and a
		cosine both seen, and neither so small that it asks for more than
		half the guard past it. Both are taken at the width the argument
		asks, which is where the two below would be taken; anything else puts
		the width back and takes them as it always did.
	*/

	bool	b_onePair = false;
	long	l_argument = lzi_small_divisor_widening( leftOperand );

	if ( ((glWorkingPrecision + l_argument) >= lzi_sincos_split_floor( leftOperand->m_lBase )) &&
		 leftOperand->FlagIsOff( kfInfinity ) && leftOperand->FlagIsOff( kfQNaN ) )
	{
		glWorkingPrecision += l_argument;

		lzi_sine_and_cosine( leftOperand, sine_result, cosine_result );

		sine_result->CheckForZero();
		cosine_result->CheckForZero();

		b_onePair = !(sine_result->IsZero()) && !(cosine_result->IsZero()) &&
					((lzi_small_divisor_widening( sine_result ) - l_argument) <= (l_guard / 2)) &&
					(((lzi_small_divisor_widening( cosine_result ) / 2) - l_argument) <= (l_guard / 2));

		if ( !b_onePair )
			glWorkingPrecision -= l_argument;
	}

	if ( !b_onePair )
	{
		long	l_widened = lzi_small_divisor_widening( leftOperand );

		glWorkingPrecision += l_widened;

		lz_sine( leftOperand, sine_result );

		// A sine or cosine too small to show at this width is taken wider
		//	until it does -- see lzi_widen_until_seen...
		l_widened += lzi_widen_until_seen( lz_sine, leftOperand, sine_result );

		long	l_wanted = lzi_small_divisor_widening( sine_result );

		if ( (l_wanted - l_widened) > (l_guard / 2) )
		{
			glWorkingPrecision += l_wanted - l_widened;
			l_widened = l_wanted;

			lz_sine( leftOperand, sine_result );
		}

		lz_cosine( leftOperand, cosine_result );

		l_widened += lzi_widen_until_seen( lz_cosine, leftOperand, cosine_result );

		// A small cosine, near the odd multiples of pi/2, is a small
		//	numerator and a small answer, and wants its size once -- see
		//	lz_tangent, which is this with the two exchanged...
		l_wanted = lzi_small_divisor_widening( cosine_result ) / 2;

		if ( (l_wanted - l_widened) > (l_guard / 2) )
		{
			glWorkingPrecision += l_wanted - l_widened;

			lz_cosine( leftOperand, cosine_result );
		}

		// lz_divide takes its width from its target, made before any of
		//	this...
	}

	result->SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( cosine_result, sine_result, result ) );

	// The guard places never leave the function, and the answer is
	//	rounded to the rule for a quotient -- see
	//	lzi_round_quotient_answer...
	glWorkingPrecision = l_savedPrecision;

	lzi_round_quotient_answer( result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete sine_result;
    delete cosine_result;

	return lResult;
}


long lz_cube( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// The operand as it came, for the check below; the answer may be
	//	written over it...
	CLZNumber	operand;
	CLZNumber	three;

	operand.Assign( leftOperand );
	three.SetBase( leftOperand->m_lBase, false );
	three.Assign( 3L );

	CLZNumber* result = nullptr;

	result = new CLZNumber;

	/*
		**The square exactly, so the product is the only rounding.** A
		square rounded to the width carried its error, times the operand,
		into the cube: 123.4567... at 36 places came out 42 units of its
		last place wrong, and -98765.43... 39,550, until 2026-10-15. A
		product's places are its operands' together, so the square is exact
		at twice the operand's...
	*/
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = (2 * leftOperand->m_lRadix > l_saved) ? (2 * leftOperand->m_lRadix) : l_saved;
	lz_multiply( leftOperand, leftOperand, result );
	glWorkingPrecision = l_saved;
	lz_multiply( result, leftOperand  );

	/*
		**Checked, where the square need not be.** The product rounds, and
		can round to nought: 3e-20 cubed is 2.7e-59, nought at 36 places,
		which would write the cube as a short 0. So a short cube is checked
		against the operand before it is called exact -- see
		lzi_write_to_width...
	*/

	if ( lzi_short_of_width( result ) && lzi_power_is_exactly( result, &operand, &three ) )
	{
		result->TrimTrailing();
		answer.Exact();
	}

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;

	return kNoError;
}


long lz_divide( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The target's rounding settings are its own -- see CLZTargetRounding.
	//	lzi_divide_op keeps them too; this covers the flags' answer...
	CLZTargetRounding	targetRounding( (targetNumber != nullptr) ? targetNumber : leftOperand );

	// Before anything -- see CLZMatchedBase.
	//	A digit array carries no base, so a mixed pair has to be
	//	brought into one before any of it is read...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Divide( targetNumber );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	// An infinity or a NaN carries through rather than being lost...
	if ( lzi_propagate_flags( leftOperand, rightOperand, targetNumber, kOpDivide ) )
		return kNoError;

	long lResult = kNoError;

	if ( targetNumber != nullptr )
		lResult = lzi_divide_op( leftOperand, rightOperand, targetNumber, targetNumber->m_lPrecision );
	else
		lResult = lzi_divide_op( leftOperand, rightOperand, nullptr, leftOperand->m_lPrecision );

	return lResult;
}


long lz_divide( CLZNumber* leftOperand, CLZNumber* rightOperand, long precision )
{
	// Before anything -- see CLZMatchedBase.
	//	A digit array carries no base, so a mixed pair has to be
	//	brought into one before any of it is read...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.DividePrecision( precision );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	// An infinity or a NaN carries through, as through the other division;
	//	until 2026-10-15 a NaN dividend was lost and a NaN divisor read as a
	//	division by nought...
	if ( lzi_propagate_flags( leftOperand, rightOperand, nullptr, kOpDivide ) )
		return kNoError;

	long lResult = kNoError;

	lResult = lzi_divide_op( leftOperand, rightOperand, nullptr, precision );

	return lResult;
}


long lz_divide_integer( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The engine's own work: its divisions are not the caller's -- see
	//	CLZInternalWork...
	CLZInternalWork	work( (targetNumber != nullptr) ? targetNumber : leftOperand );

	long lResult = kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.DivideInteger( targetNumber );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	// An infinity or a NaN carries through as through division -- infinity
	//	over a finite value infinite, a finite one over infinity nought, as
	//	the decimal arithmetic's divide-integer has them. Until 2026-10-15
	//	nothing was propagated...
	if ( lzi_propagate_flags( leftOperand, rightOperand, targetNumber, kOpDivide ) )
		return kNoError;

	// Truncated straight from the raw quotient -- see lzi_divide_work.
	//	Until 2026-09-26 this truncated a quotient already rounded at the
	//	precision, and 1 / 0.333...334 to fifty places came back 3...
	lResult = lzi_first_failure( lResult, lzi_divide_whole( leftOperand, rightOperand, targetNumber, RM_round_down ) );

	return lResult;
}


static long lzi_exp_b_work( CLZNumber* leftOperand, CLZNumber*, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	long	lResult = kNoError;

	CLZNumber*	natural_log = nullptr;
	CLZNumber*	integer_base = nullptr;
	CLZNumber*	fractional_base = nullptr;
	CLZNumber*	result = nullptr;
	CLZNumber*	working = nullptr;

	// Create working objects...
	natural_log = new CLZNumber;
	integer_base = new CLZNumber;
	fractional_base = new CLZNumber;
	result = new CLZNumber;
	working = new CLZNumber;

	/*
		gLogBaseE and gBaseN are held in base ten and are the left
		operands below, so each is brought into the operand's base
		first -- see CLZAnswerBase and lzi_table_in_base. This used to convert
		the operand to base ten instead and the answer back, from before
		mixed operands were matched, when a right operand's digits were
		read in the left one's base and the square root of eight came
		back 2.297. Worked in base ten, a hex answer carried the
		decimal's precision re-counted, and 16 to the 12.5 was 6300
		units out of its 36th place.
	*/

	/*
		**ln(base) as wide as the work.** It came from gLogBaseE, 120
		places of text, whatever the width: 10^x was out by 1e20 units at
		120 places and 1e50 at 150, and a large hex power missed from 80
		places, where its own size had widened the work past the text.
		See lzi_log_base_in_base.
	*/

	long		l_operandBase = leftOperand->m_lBase;
	CLZNumber*	p_baseN = lzi_table_in_base( gBaseN[ l_operandBase ], l_operandBase );

	working->Assign( leftOperand );

	// Initial values...
	integer_base->Assign( working );
	fractional_base->Assign( working );

	integer_base->StripFraction();
	fractional_base->StripInteger();

	/*
		**A whole operand is the base to a whole power, and needs no
		logarithm.** ln(base) was worked out at the full width before
		anything looked at the fraction -- a width that grows with the
		answer -- and then multiplied by nought: b to the 99,999 spent 96
		per cent of 22 seconds on it. The squaring alone gives the same
		answer, since e to the nought is one; and squaring a power of the
		base is cheap in every base, the digit loop skipping noughts.
	*/

	if ( fractional_base->IsZero() )
		lzi_raise_n( p_baseN, integer_base, result );
	else
	{
		CLZNumber*	p_logBase = lzi_log_base_in_base( l_operandBase, glWorkingPrecision );

		// Raise the integer & fractional parts...
		lz_multiply( p_logBase, fractional_base, natural_log );
		lzi_raise_n( p_baseN, integer_base, integer_base );

		// Calculate the result...
		lz_exp_e( natural_log, fractional_base );
		lz_multiply( integer_base, fractional_base, result );
	}

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete natural_log;
    delete integer_base;
    delete fractional_base;
    delete result;
    delete working;

	return lResult;
}


long lz_exp_b( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_exp_b, kRetakeDoubt );

	// A NaN carries through, and an infinity -- b^x past either end -- needs
	//	the sign this engine's infinities do not carry, as e^x's does. Until
	//	2026-10-15 nothing was propagated and b^NaN answered 1 with success...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	// b^x has about x whole places in its own base...
	double	d_x = 0.0;
	long	l_magnitude = -1;

	if ( lzi_value_or_infinity( leftOperand, &d_x ) )
	{
		double	d_lnAnswer = d_x * std::log( static_cast<double>( leftOperand->m_lBase ) );
		long	l_beyond = kNoError;

		// Refused past the limit, nought below the width -- see
		//	lzi_answer_size...
		if ( lzi_answer_beyond( lzi_answer_size( d_lnAnswer, leftOperand->m_lBase ),
								leftOperand, targetNumber, &l_beyond ) )
			return l_beyond;

		l_magnitude = lzi_places_for( d_lnAnswer, leftOperand->m_lBase );
	}

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	// The operand as it came, and the base as a number, for the check
	//	below...
	CLZNumber	exponent;
	CLZNumber	base;

	exponent.Assign( leftOperand );
	base.SetBase( leftOperand->m_lBase, false );
	base.Assign( leftOperand->m_lBase );

	long	lResult = lzi_widened_by_answer( lzi_exp_b_work, leftOperand, nullptr, targetNumber, l_magnitude );

	// A short answer is exact only if it is the base to the operand --
	//	see lzi_write_to_width...
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( lzi_short_of_width( p_answer ) && lzi_power_is_exactly( p_answer, &base, &exponent ) )
	{
		p_answer->TrimTrailing();
		answer.Exact();
	}

	return lResult;
}


/*
	The computed series for e^f, f a fraction: halved first past a
	thousand digits where two divides the base, summed term by term, and
	squared back -- see kExpComputedHalvings and the note in lz_exp_e.
	bMayHalve is false for the tail of a split fraction, which is small
	enough already. pSum takes the fraction's base.
*/

static void lzi_exp_fraction_computed( CLZNumber* pFraction, long lDigits, bool bMayHalve, CLZNumber* pSum )
{
	long	l_base = pFraction->m_lBase;
	long	l_halvings = (bMayHalve && (lDigits > kExpComputedHalvedFrom) && ((l_base % 2) == 0))
							 ? kExpComputedHalvings : 0;
	long	l_extra = 0;

	if ( l_halvings > 0 )
		l_extra = static_cast<long>( std::ceil( static_cast<double>( l_halvings ) * std::log( 2.0 ) /
												std::log( static_cast<double>( l_base ) ) ) ) + 2;

	glWorkingPrecision += l_extra;

	CLZNumber	argument;
	CLZNumber	term;
	CLZNumber	divisor;

	// Assign carries the source's precision, and the halving
	//	divides in place...
	argument.Assign( pFraction );
	argument.SetPrecision( glWorkingPrecision );

	if ( l_halvings > 0 )
	{
		CLZNumber	two;

		two.Assign( 2L );

		for ( long i = 0; i < l_halvings; i++ )
			lz_divide( &argument, &two );
	}

	// The sum in the fraction's base, or a fresh one sums in base ten...
	term.SetBase( l_base );
	term.Assign( 1L );
	pSum->SetBase( l_base, false );
	pSum->Assign( 1L );

	long	l_bound = lzi_series_bound( l_base );

	for ( long k = 1; k <= l_bound; k++ )
	{
		divisor.Assign( k );

		lz_multiply( &term, &argument );
		lz_divide( &term, &divisor );

		term.CheckForZero();

		if ( term.IsZero() )
			break;

		lz_add( pSum, &term );
	}

	if ( l_halvings > 0 )
	{
		CLZNumber	squared;

		for ( long i = 0; i < l_halvings; i++ )
		{
			lz_multiply( pSum, pSum, &squared );
			pSum->Assign( &squared );
		}
	}

	glWorkingPrecision -= l_extra;
}


/*
	**e to a fraction by binary splitting, a chunk at a time: the bit-burst
	method** -- part B of the slowness item, 2026-10-13. lz_exp_e takes it
	past lzi_exp_split_floor, in place of the series above.

	The series above sums a term at a time, and past a few thousand
	places that is tens of thousands of terms, each a full-width addition
	and a division by k -- quadratic in the width. Here the fraction f is
	cut into chunks of doubling length, f = c0 + c1 + ..., so that
	e^f = e^c0 e^c1 ..., and each chunk is p / b^e, a whole p of the
	chunk's places over the base to the chunk's last place. Its series is
	summed by binary splitting: the terms from a to b - 1 of x^k / k!
	come to T / (Q b^S), where

		a leaf, b = a + 1:	P = p, Q = a, T = p, S = e
		a join:				P = P1 P2, Q = Q1 Q2, S = S1 + S2,
							T = T1 Q2 b^S2 + P1 T2

	every number whole and exact, b^S2 a shift. So the only rounding is
	one division a chunk and the products of the chunks. A chunk past the
	first is under b to the minus its first place, so it needs about the
	width over that many terms, and its numbers stay about the width long
	whichever chunk it is.
*/

// Terms a to b - 1 of e^(p / b^lPlaces) - 1, as T / (Q b^S)...
static void lzi_exp_split( CLZNumber* pP, long lPlaces, long lFrom, long lTo,
						   CLZNumber* pPOut, CLZNumber* pQOut, CLZNumber* pTOut, long* plShift )
{
	long	l_base = pP->m_lBase;

	if ( (lTo - lFrom) == 1 )
	{
		pPOut->SetBase( l_base, false );
		pPOut->Assign( pP );
		pQOut->SetBase( l_base, false );
		pQOut->Assign( lFrom );
		pTOut->SetBase( l_base, false );
		pTOut->Assign( pP );
		*plShift = lPlaces;

		return;
	}

	long		l_middle = (lFrom + lTo) / 2;
	long		l_rightShift = 0;
	CLZNumber	right_p, right_q, right_t, product;

	lzi_exp_split( pP, lPlaces, lFrom, l_middle, pPOut, pQOut, pTOut, plShift );
	lzi_exp_split( pP, lPlaces, l_middle, lTo, &right_p, &right_q, &right_t, &l_rightShift );

	// T = T1 Q2 b^S2 + P1 T2, before P1 becomes P1 P2...
	lz_multiply( pTOut, &right_q );
	lzi_scale_by_base_power( pTOut, -l_rightShift );
	lz_multiply( pPOut, &right_t, &product );
	lz_add( pTOut, &product );

	lz_multiply( pPOut, &right_p );
	lz_multiply( pQOut, &right_q );
	*plShift += l_rightShift;
}

/*
	**The first chunk is sixteen bits, in as many places as the base needs**
	-- the fewest whose power of the base reaches 2^16: sixteen in base two,
	eight in four, five in ten, four in sixteen. Measured 2026-10-13 against
	first chunks of two to twenty-four places: four places was best in bases
	seven and ten, eleven and thirteen bits, and 21 to 43% slow in base two,
	where it is four bits; base two was best at twelve to sixteen places,
	base four at eight to twelve, sixteen at three or four. Counted in whole
	numbers, not a logarithm, which can put four places in base sixteen at
	4.0000001 and take five.
*/
const long	kExpSplitFirstBits = 16;
const long	kExpSplitGuardDigits = 6;

// The floor in force, which only tests and measurements change, and the
//	count of the kernel's runs they read back. The setting is the process's
//	and the count each thread's own...
static long gExpSplitFloor = kExpSplitFloorMeasured;
static thread_local long gExpSplittings = 0;

long lzi_set_exp_split_floor( long lPlaces )
{
	long	l_previous = gExpSplitFloor;

	gExpSplitFloor = lPlaces;

	return l_previous;
}

long lzi_exp_splittings()
{
	return gExpSplittings;
}

/*
	**Guarded on its own and not inside lz_exp_e.** Called directly it works
	kExpSplitGuardDigits wider, two digits more than it needs: at four it
	is clean over 1,500 random fractions, at three it differs from
	lz_exp_e, at one it is 34 wrong. lz_exp_e already works kExpGuardDigits
	wider and retakes a doubtful answer, and the two guards nest -- the
	reason kExpGuardDigits gives for being four -- so it takes the kernel
	with none of its own: 3,000 random fractions through lz_exp_e, every
	width past the table path in every base, none wrong either way, and
	3 to 7% faster near an odd base's floor, under 2% by a thousand places.
*/

static long lzi_exp_split_work( CLZNumber* pFraction, CLZNumber* pTarget, long lFirstChunk, long lGuardDigits )
{
	gExpSplittings++;

	long	l_base = pFraction->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_width = l_saved + ((lGuardDigits > 0) ? lzi_guard_places( lGuardDigits, l_base ) : 0);
	double	d_lnBase = std::log( static_cast<double>( l_base ) );

	glWorkingPrecision = l_width;

	// The fraction to the working places; what is past them moves e^f by
	//	under a unit of the last...
	CLZNumber	fraction;

	fraction.Assign( pFraction );
	fraction.SetPrecision( l_width );
	fraction.Round( RM_round_down );
	fraction.TrimTrailing();

	long		l_places = fraction.m_lRadix;
	long		l_start = 0;
	long		l_length = lFirstChunk;

	if ( l_length <= 0 )
	{
		long	l_reach = 1;

		for ( l_length = 0; l_reach < (1L << kExpSplitFirstBits); l_length++ )
			l_reach *= l_base;
	}
	CLZNumber	product;

	product.SetBase( l_base, false );
	product.Assign( 1L );

	while ( l_start < l_places )
	{
		long	l_end = ((l_start + l_length) < l_places) ? (l_start + l_length) : l_places;

		// The chunk's places as a whole number: fraction place i is digit
		//	m_lRadix - i, least significant first...
		CLZNumber	p;

		p.SetBase( l_base, false );
		p.SetDigitsArraySize( (l_end - l_start) + 8, false );

		for ( long i = 0; i < (l_end - l_start); i++ )
			p.m_lpDigits[ i ] = fraction.m_lpDigits[ fraction.m_lRadix - l_end + i ];

		// A negative fraction -- lz_exp_e's for a negative argument -- gives
		//	every chunk its sign, and the splitting is the same signed:
		//	e^c = 1 + T / (Q b^S) still, its terms alternating...
		p.m_lDigitsLength = l_end - l_start;
		p.m_lRadix = 0;
		p.m_lSign = fraction.m_lSign;
		p.TrimLeading();
		p.CheckForZero();

		if ( !(p.IsZero()) )
		{
			// The terms: the fewest K for which the chunk to the K over K!,
			//	the chunk being under b^-start, falls past the width...
			double	d_need = static_cast<double>( l_width + 2 ) * d_lnBase;
			long	l_terms = 1;

			while ( ((static_cast<double>( l_terms * l_start ) * d_lnBase) +
					 lzi_ln_gamma_estimate( static_cast<double>( l_terms ) + 1.0 )) < d_need )
				l_terms++;

			CLZNumber	split_p, split_q, split_t, chunk;
			long		l_shift = 0;

			lzi_exp_split( &p, l_end, 1, l_terms + 1, &split_p, &split_q, &split_t, &l_shift );

			// e^c = 1 + T / (Q b^S): T moved down S places first, which is
			//	exact, so the quotient is taken to the width's places...
			lzi_scale_by_base_power( &split_t, l_shift );
			chunk.SetBase( l_base, false );
			chunk.SetPrecision( l_width );
			lz_divide( &split_t, &split_q, &chunk );

			CLZNumber	one;

			one.SetBase( l_base, false );
			one.Assign( 1L );
			lz_add( &chunk, &one );
			lz_multiply( &product, &chunk );
		}

		l_start = l_end;
		l_length *= 2;
	}

	glWorkingPrecision = l_saved;
	product.SetPrecision( l_saved );
	product.Round();

	if ( !(pTarget->IsConstant()) )
		pTarget->Assign( &product );

	return kNoError;
}

long lzi_exp_by_splitting( CLZNumber* pFraction, CLZNumber* pTarget, long lFirstChunk )
{
	return lzi_exp_split_work( pFraction, pTarget, lFirstChunk, kExpSplitGuardDigits );
}


/*
	**Where lz_exp_e takes the splitting kernel: from a width that depends
	on the base**, in places of the base, against the width the series
	works at -- the caller's, with lz_exp_e's guard and the answer's
	magnitude on top, a few places that sit inside the measurement's steps.
	Below it the series stays, and below 48 digits the table path.

	Measured 2026-10-13: lz_exp_e
	with the kernel forced on and forced off through
	lzi_set_exp_split_floor, at eighteen widths from 50 to 3,000 places and
	three fractions -- twelve places, 64, and as long as the width, as
	the logarithm's Newton steps hand it -- two of each. Each floor is
	the width at which the worst loss over every width and fraction is
	least, the rule the transform's floors use.

	**An odd base takes it from 75 to 150 places**, the worst loss 0.5%
	to 17% by base: the series cannot halve its argument there. **An even
	base halves, and holds the kernel off to 500 to 800 places**, the worst
	loss 8% to 18%: a twelve-place fraction gains from the kernel from 150
	to 600 places, by base, and a long one loses to the halving series
	until 700 to 1,000, and one floor has to choose between them. **Base two takes it wherever the
	series would run**: its halvings are cheapest of all, but its
	sixteen-place first chunk made the kernel faster still, by 2% at the
	least from the first width past the table path.

	At 3,000 places the kernel is 1.8 to 3.7 times faster in an even base
	and 6.8 to 12.2 times in an odd one.
*/

long lzi_exp_split_floor( long lBase )
{
	if ( gExpSplitFloor != kExpSplitFloorMeasured )
		return gExpSplitFloor;

	switch ( lBase )
	{
	case 2:		return 0;
	case 3:		return 150;
	case 4:		return 500;
	case 5:		return 150;
	case 6:		return 600;
	case 7:		return 150;
	case 8:		return 800;
	case 9:		return 100;
	case 10:	return 600;
	case 11:	return 150;
	case 12:	return 700;
	case 13:	return 100;
	case 14:	return 600;
	case 15:	return 75;
	case 16:	return 700;
	default:	return LONG_MAX;
	}
}

/*
	**Cosine and sine together by binary splitting, a chunk at a time** --
	the bit-burst method of lzi_exp_split_work, for e^{iy} = cos y +
	i sin y. Performance item 1 of the review of 2026-10-15, on Rich's
	choice of a kernel with Newton for the inverses: the series sums a term
	at a time, and its terms grow with the width and each is a pass over
	it, so sine and cosine grew about as the square of the width.

	y is cut into chunks of doubling length as e^f is, so that e^{iy} is
	the product of each chunk's, and a chunk is p / b^e. Its terms from a
	to b - 1 of e^{iy} - 1 come to (Tre + i Tim) / (Q b^S), where P is
	(i p)^(b - a) and is held as its magnitude p^(b - a), sign and all,
	and its phase, b - a counted in quarter turns:

		a leaf, b = a + 1:	P = p, phase 1, Q = a, Tre = 0, Tim = p, S = e
		a join:				T = T1 Q2 b^S2 + i^phase1 |P1| T2, the
							turn taking (u, v) to (-v, u); P, Q and S as
							for e^f, the phases adding

	every number whole and exact, as there. So a chunk gives cos c = 1 +
	Tre / (Q b^S) and sin c = Tim / (Q b^S), and the chunks combine by
	cos(a + c) = cos a cos c - sin a sin c and sin(a + c) = sin a cos c +
	cos a sin c.

	**That is the whole part's road, and a chunk after it takes a cheaper
	one.** Every chunk past the whole part is under one, so its cosine is
	positive and sqrt(1 - sin^2 c) is well conditioned: its sine alone is
	summed, sin c = c (1 + T / (Q b^S)) over the series in c^2 = -p^2 / b^2e
	-- lzi_sine_split, the exponential's splitting with half the terms and
	four products a join for the pair's six -- and its cosine is one root.
	A whole part, which can pass pi/2 where the cosine's sign matters,
	keeps the pair. Measured against the pair for every chunk: 1.3 to 1.4
	times faster for an argument as long as the width, 1.8 for a short
	one, the census unchanged. And the rightmost joins skip P, which
	nothing reads: 3 to 4% more.

	It works to the places, as sine and cosine are answered: every value it
	handles is under one or about it, so a unit of the last place is the
	error that matters, whatever the answer's size.
*/

// Terms a to b - 1 of e^{i p / b^lPlaces} - 1, as (Tre + i Tim) / (Q b^S)...
static void lzi_sincos_split( CLZNumber* pP, long lPlaces, long lFrom, long lTo,
							  CLZNumber* pPOut, long* plPhase, CLZNumber* pQOut,
							  CLZNumber* pTre, CLZNumber* pTim, long* plShift )
{
	long	l_base = pP->m_lBase;

	if ( (lTo - lFrom) == 1 )
	{
		pPOut->SetBase( l_base, false );
		pPOut->Assign( pP );
		*plPhase = 1;
		pQOut->SetBase( l_base, false );
		pQOut->Assign( lFrom );
		pTre->SetBase( l_base, false );
		pTre->Assign( 0L );
		pTim->SetBase( l_base, false );
		pTim->Assign( pP );
		*plShift = lPlaces;

		return;
	}

	long		l_middle = (lFrom + lTo) / 2;
	long		l_rightShift = 0;
	long		l_rightPhase = 0;
	CLZNumber	right_p, right_q, right_re, right_im, turned_re, turned_im;

	lzi_sincos_split( pP, lPlaces, lFrom, l_middle, pPOut, plPhase, pQOut, pTre, pTim, plShift );
	lzi_sincos_split( pP, lPlaces, l_middle, lTo, &right_p, &l_rightPhase, &right_q, &right_re, &right_im, &l_rightShift );

	// T1 Q2 b^S2...
	lz_multiply( pTre, &right_q );
	lz_multiply( pTim, &right_q );
	lzi_scale_by_base_power( pTre, -l_rightShift );
	lzi_scale_by_base_power( pTim, -l_rightShift );

	// ...and |P1| T2 turned by P1's phase, before P1 becomes P1 P2...
	lz_multiply( pPOut, &right_re, &turned_re );
	lz_multiply( pPOut, &right_im, &turned_im );

	switch ( *plPhase & 3 )
	{
	case 0:		lz_add( pTre, &turned_re );			lz_add( pTim, &turned_im );			break;
	case 1:		lz_subtract( pTre, &turned_im );	lz_add( pTim, &turned_re );			break;
	case 2:		lz_subtract( pTre, &turned_re );	lz_subtract( pTim, &turned_im );	break;
	default:	lz_add( pTre, &turned_im );			lz_subtract( pTim, &turned_re );	break;
	}

	lz_multiply( pPOut, &right_p );
	lz_multiply( pQOut, &right_q );
	*plPhase += l_rightPhase;
	*plShift += l_rightShift;
}

// Terms a to b - 1 of sin(c) / c - 1, c = p / b^lPlaces, as T / (Q b^S):
//	the series in pSquare = -p^2, so a leaf is P = T = -p^2, Q = 2a (2a + 1),
//	S = 2 lPlaces. bNeedP is false where nothing reads P: the root, and every
//	right child of a node whose own P is not read...
static void lzi_sine_split( CLZNumber* pSquare, long lPlaces, long lFrom, long lTo,
							CLZNumber* pPOut, CLZNumber* pQOut, CLZNumber* pTOut, long* plShift, bool bNeedP )
{
	long	l_base = pSquare->m_lBase;

	if ( (lTo - lFrom) == 1 )
	{
		pPOut->SetBase( l_base, false );
		pPOut->Assign( pSquare );
		pQOut->SetBase( l_base, false );
		pQOut->Assign( (2 * lFrom) * ((2 * lFrom) + 1) );
		pTOut->SetBase( l_base, false );
		pTOut->Assign( pSquare );
		*plShift = 2 * lPlaces;

		return;
	}

	long		l_middle = (lFrom + lTo) / 2;
	long		l_rightShift = 0;
	CLZNumber	right_p, right_q, right_t, product;

	lzi_sine_split( pSquare, lPlaces, lFrom, l_middle, pPOut, pQOut, pTOut, plShift, true );
	lzi_sine_split( pSquare, lPlaces, l_middle, lTo, &right_p, &right_q, &right_t, &l_rightShift, bNeedP );

	// T = T1 Q2 b^S2 + P1 T2, before P1 becomes P1 P2...
	lz_multiply( pTOut, &right_q );
	lzi_scale_by_base_power( pTOut, -l_rightShift );
	lz_multiply( pPOut, &right_t, &product );
	lz_add( pTOut, &product );

	if ( bNeedP )
		lz_multiply( pPOut, &right_p );

	lz_multiply( pQOut, &right_q );
	*plShift += l_rightShift;
}

// The guard the kernel takes called directly, and the floor and count the
//	tests and measurements read -- see lzi_sincos_split_floor...
const long	kSinCosSplitGuardDigits = 6;

static long gSinCosSplitFloor = kSinCosSplitFloorMeasured;
static thread_local long gSinCosSplittings = 0;
static thread_local long gSinCosLastChunks = 0;

long lzi_set_sincos_split_floor( long lPlaces )
{
	long	l_previous = gSinCosSplitFloor;

	gSinCosSplitFloor = lPlaces;

	return l_previous;
}

/*
	**Where sine and cosine take the kernel: from a width that depends on the
	base**, in places, against the working width lzi_trig_series sees, the
	guard included -- a few places inside the measurement's steps. Rich's
	choice, 2026-10-16, of one floor a base.

	Measured: lz_sine and
	lz_cosine with the kernel forced on and off, at eighteen widths from 50
	to 3,000 digits and three arguments -- twelve places, 64, and as long as
	the width -- each inside pi/4 and past it, where the fold leaves an
	argument as long as the width. **Each floor is the least width from
	which the kernel is never more than 5% slower than the series** at any
	width and argument measured: 250 to 600 digits by base. The fold decides
	the crossover more than the base does: a short argument needs a chunk or
	two and the kernel wins from 100 to 300 digits, a long one needs about
	six and wins from 250 to 600, and one floor has to serve both. At 3,000
	digits the kernel is 3.0 times faster in base ten and 4.0 to 5.0 in the
	others, at the least over every argument.
*/

long lzi_sincos_split_floor( long lBase )
{
	if ( gSinCosSplitFloor != kSinCosSplitFloorMeasured )
		return gSinCosSplitFloor;

	switch ( lBase )
	{
	case 2:		return 997;		// 300 digits
	case 3:		return 629;		// 300
	case 4:		return 498;		// 300
	case 5:		return 358;		// 250
	case 6:		return 643;		// 500
	case 7:		return 473;		// 400
	case 8:		return 332;		// 300
	case 9:		return 629;		// 600
	case 10:	return 500;		// 500
	case 11:	return 384;		// 400
	case 12:	return 371;		// 400
	case 13:	return 449;		// 500
	case 14:	return 262;		// 300
	case 15:	return 340;		// 400
	case 16:	return 415;		// 500
	default:	return LONG_MAX;
	}
}

long lzi_sincos_splittings()
{
	return gSinCosSplittings;
}

long lzi_sincos_last_chunks()
{
	return gSinCosLastChunks;
}

static long lzi_sincos_split_work( CLZNumber* pArgument, CLZNumber* pSine, CLZNumber* pCosine, long lGuardDigits )
{
	gSinCosSplittings++;
	gSinCosLastChunks = 0;

	long	l_base = pArgument->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_width = l_saved + ((lGuardDigits > 0) ? lzi_guard_places( lGuardDigits, l_base ) : 0);
	double	d_lnBase = std::log( static_cast<double>( l_base ) );

	glWorkingPrecision = l_width;

	// The argument to the working places; what is past them moves both by
	//	under a unit of the last...
	CLZNumber	argument;

	argument.Assign( pArgument );
	argument.SetPrecision( l_width );
	argument.Round( RM_round_down );
	argument.TrimTrailing();

	// Its whole part, the first chunk when there is one...
	CLZNumber	whole;

	whole.Assign( &argument );
	whole.Round( RM_round_down, 0L );
	whole.TrimAll();
	whole.CheckForZero();

	long	l_places = argument.m_lRadix;
	long	l_start = 0;
	long	l_length = 0;

	for ( long l_reach = 1; l_reach < (1L << kExpSplitFirstBits); l_length++ )
		l_reach *= l_base;

	CLZNumber	cosine, sine;
	bool		b_first = true;
	bool		b_whole = !(whole.IsZero());

	cosine.SetBase( l_base, false );
	cosine.Assign( 1L );
	sine.SetBase( l_base, false );
	sine.Assign( 0L );

	/*
		**The fraction's chunks are taken to the nearest, not cut**, and what
		each leaves -- either sign -- goes on to the next. Cut, 0.6999...9xyz
		had a digit in every chunk and paid for all of them, where 0.7000...0xyz
		skipped the noughts; and the reductions leave either, by the sign of
		a rounding residue: cosine of 0.7 at 3,000 places reached the kernel
		with five digits that were not nought, the same cosine inside secant
		and tangent with 3,014, and cost four times as much. To the nearest,
		both are 0.7 and a tail. A chunk so taken is at most a unit, and each
		after the first under half a unit of the place before it, so the
		bounds the term counts read are as they were, and its cosine is still
		positive.
	*/

	CLZNumber	remaining;

	lz_subtract( &argument, &whole, &remaining );

	while ( b_whole || (l_start < l_places) )
	{
		CLZNumber	p;
		long		l_end = 0;
		double		d_lnChunk = 0.0;		// the chunk is under e^this

		p.SetBase( l_base, false );

		if ( b_whole )
		{
			p.Assign( &whole );
			p.m_lRadix = 0;
			d_lnChunk = static_cast<double>( p.m_lDigitsLength ) * d_lnBase;
		}
		else
		{
			// The remainder to the chunk's last place, as a whole number, and
			//	what is left of it for the chunks after...
			CLZNumber	chunk;

			l_end = ((l_start + l_length) < l_places) ? (l_start + l_length) : l_places;

			chunk.Assign( &remaining );
			chunk.Round( RM_round_half_even, l_end );
			lz_subtract( &remaining, &chunk );

			p.Assign( &chunk );
			lzi_scale_by_base_power( &p, -l_end );
			p.TrimAll();
			d_lnChunk = -static_cast<double>( l_start ) * d_lnBase;
		}

		// Each chunk carries its own sign -- the whole part the argument's, a
		//	remainder its own; sine follows it and cosine does not, the
		//	splitting signed as e^f's is...
		p.TrimLeading();
		p.CheckForZero();

		if ( !(p.IsZero()) )
		{
			gSinCosLastChunks++;

			// The fewest terms K for which the chunk to the K over K! falls
			//	past the width: K! over the chunk to the K past b^width...
			double	d_need = static_cast<double>( l_width + 2 ) * d_lnBase;
			long	l_terms = 1;

			while ( (lzi_ln_gamma_estimate( static_cast<double>( l_terms ) + 1.0 ) -
					 (static_cast<double>( l_terms ) * d_lnChunk)) < d_need )
				l_terms++;

			CLZNumber	split_p, split_q, split_re, split_im, chunk_cos, chunk_sin, one;
			long		l_shift = 0;

			chunk_cos.SetBase( l_base, false );
			chunk_cos.SetPrecision( l_width );
			chunk_sin.SetBase( l_base, false );
			chunk_sin.SetPrecision( l_width );
			one.SetBase( l_base, false );
			one.Assign( 1L );

			if ( b_whole )
			{
				// The pair: T moved down S places, exactly, then divided to the
				//	width...
				long	l_phase = 0;

				lzi_sincos_split( &p, l_end, 1, l_terms + 1, &split_p, &l_phase, &split_q,
								  &split_re, &split_im, &l_shift );

				lzi_scale_by_base_power( &split_re, l_shift );
				lzi_scale_by_base_power( &split_im, l_shift );

				lz_divide( &split_re, &split_q, &chunk_cos );
				lz_divide( &split_im, &split_q, &chunk_sin );
				lz_add( &chunk_cos, &one );
			}
			else
			{
				// Sine alone, half the terms of the pair's, then cosine by a root:
				//	sin c = c (1 + T / (Q b^S)), cos c = sqrt(1 - sin^2 c)...
				CLZNumber	square, ratio, chunk, sine_squared;

				lz_multiply( &p, &p, &square );
				lz_negate( &square );
				lzi_sine_split( &square, l_end, 1, (l_terms / 2) + 2, &split_p, &split_q, &split_re, &l_shift, false );
				lzi_scale_by_base_power( &split_re, l_shift );

				ratio.SetBase( l_base, false );
				ratio.SetPrecision( l_width );
				lz_divide( &split_re, &split_q, &ratio );
				lz_add( &ratio, &one );

				chunk.Assign( &p );
				lzi_scale_by_base_power( &chunk, l_end );
				lz_multiply( &ratio, &chunk, &chunk_sin );

				lz_multiply( &chunk_sin, &chunk_sin, &sine_squared );
				lz_subtract( &one, &sine_squared, &chunk_cos );
				lz_root_2( &chunk_cos );
			}

			// The first chunk is the pair as it stands; each after it turns
			//	the pair by its own angle...
			if ( b_first )
			{
				cosine.Assign( &chunk_cos );
				sine.Assign( &chunk_sin );
				b_first = false;
			}
			else
			{
				CLZNumber	cos_cos, sin_sin, sin_cos, cos_sin;

				lz_multiply( &cosine, &chunk_cos, &cos_cos );
				lz_multiply( &sine, &chunk_sin, &sin_sin );
				lz_multiply( &sine, &chunk_cos, &sin_cos );
				lz_multiply( &cosine, &chunk_sin, &cos_sin );
				lz_subtract( &cos_cos, &sin_sin, &cosine );
				lz_add( &sin_cos, &cos_sin, &sine );
			}
		}

		if ( b_whole )
		{
			b_whole = false;
			continue;
		}

		l_start = l_end;
		l_length *= 2;
	}

	glWorkingPrecision = l_saved;

	cosine.SetPrecision( l_saved );
	cosine.Round();
	sine.SetPrecision( l_saved );
	sine.Round();

	if ( (pCosine != nullptr) && !(pCosine->IsConstant()) )
		pCosine->Assign( &cosine );

	if ( (pSine != nullptr) && !(pSine->IsConstant()) )
		pSine->Assign( &sine );

	return kNoError;
}

long lzi_sincos_by_splitting( CLZNumber* pArgument, CLZNumber* pSine, CLZNumber* pCosine )
{
	return lzi_sincos_split_work( pArgument, pSine, pCosine, kSinCosSplitGuardDigits );
}


/*
	**Sine and cosine of one argument from one fold and one kernel run** --
	for tangent and cotangent past the kernel's floor, since 2026-10-16.
	They took the two by separate calls, a cosine as the sine of x + pi/2
	and then the sine, so the kernel ran twice where one run gives both,
	and the cosine's fold went through pi/2 and back. Here lzi_trig_fold
	folds x once, the kernel takes the reduced argument, and the fold's
	answers put each back: which is which, and the signs.

	Under the trig guard, as lz_sine works, and rounded back to the width
	it was asked at -- a few units at most from what lz_sine and lz_cosine
	answer there, which the quotient's own guard and retake absorb.
*/

static void lzi_sine_and_cosine( CLZNumber* pArgument, CLZNumber* pSine, CLZNumber* pCosine )
{
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kTrigGuardDigits, pArgument->m_lBase );

	CLZNumber	reduced;
	CLZNumber	sine;
	CLZNumber	cosine;
	bool		b_negateSine = false;
	bool		b_negateCosine = false;
	bool		b_swapped = false;

	lzi_trig_fold( pArgument, &reduced, &b_negateSine, &b_negateCosine, &b_swapped );
	lzi_sincos_split_work( &reduced, &sine, &cosine, 0 );

	// Folded past pi/4 the reduced argument's cosine is the sine...
	CLZNumber*	p_sine = b_swapped ? &cosine : &sine;
	CLZNumber*	p_cosine = b_swapped ? &sine : &cosine;

	if ( b_negateSine )
		lz_negate( p_sine );

	if ( b_negateCosine )
		lz_negate( p_cosine );

	glWorkingPrecision = l_saved;

	p_sine->SetPrecision( l_saved );
	p_sine->Round();
	p_cosine->SetPrecision( l_saved );
	p_cosine->Round();

	pSine->Assign( p_sine );
	pCosine->Assign( p_cosine );
}


/*
	**An angle by Newton on the sine/cosine kernel** -- stage 3 of
	performance item 1 of the review of 2026-10-15, for the inverse
	functions, which all come down to lz_arcsine. Given s and c, not both
	nought, this finds theta = atan2(s, c) by

		t = (s cos y - c sin y) / (c cos y + s sin y),   y <- y + t

	t being tan(theta - y) exactly, whatever the length of (s, c): so the
	error after a step is the tangent's cube over three, the digits about
	tripling, and the denominator is |(s, c)| cos(theta - y), near its
	length, wherever theta lies -- no reflection, no trouble near a right
	angle. One kernel run a step gives both sine and cosine. The widths
	are planned back from the full one, each a third of the next and two
	places more, from the 14 digits a double's atan2 starts it with.

	It works at the working precision under the caller's guard, as the
	kernel does under the series'.

	**Each stage but the last at its own width, the seed at its own** --
	since 2026-10-17, the Newton audit. Three things were wide, every one
	carried from stage to stage:

	- **the seed**, about 3.3 times the full width in every base but ten.
	  It was made at the full width, set to base ten unconverted and then
	  converted, and a precision is a count, so it wrote the full width's
	  worth of decimal places in the operand's base -- 110,422 places in
	  base two at 10,000 digits, for a seed good to fourteen;
	- **the iterate**, the step being a small quotient that answers to the
	  width's significant digits and its noughts, added unrounded: a third
	  again past its stage;
	- **s and c**, at their full length in all four products every stage.

	Now the seed is converted at its text's places and rounded to the first
	stage, and each stage but the last narrows s and c to the stage and two
	places and rounds the iterate to the stage. The last takes s and c
	whole, and the caller rounds the answer. Over all fifteen bases 12 to
	14% faster at 600 and 3,000 digits, base ten the least at 2 to 6% and
	base two the most at up to 25%; base two at 10,000 digits 140 ms to
	118. Of 15,900 answers censused none moved at the natural floors, nor
	any near tie within the retake's reach.
*/

const long	kAngleNewtonSeedDigits = 14;

static long gArcSineNewtonFloor = kArcSineNewtonFloorMeasured;
static thread_local long gArcSineNewtons = 0;

long lzi_set_arcsine_newton_floor( long lPlaces )
{
	long	l_previous = gArcSineNewtonFloor;

	gArcSineNewtonFloor = lPlaces;

	return l_previous;
}

/*
	**Where lz_arcsine takes Newton: from a width that depends on the base**,
	in places, against its working width, the guard included. Measured as
	sine and cosine's floors were:
	lz_arcsine forced each way at eighteen widths from 50 to 3,000 digits,
	arguments twelve places long, 64 and as long as the width, each under a
	tenth -- where the series needs no halving -- about a half, and a hair
	from one. Each floor is the least width from which Newton is never more
	than 5% slower than the series at any width and argument measured:
	250 to 600 digits by base. At 3,000 digits Newton is 2.9 to 6.7 times
	faster, at the least over every argument.
*/

long lzi_arcsine_newton_floor( long lBase )
{
	if ( gArcSineNewtonFloor != kArcSineNewtonFloorMeasured )
		return gArcSineNewtonFloor;

	switch ( lBase )
	{
	case 2:		return 830;		// 250 digits
	case 3:		return 524;		// 250
	case 4:		return 498;		// 300
	case 5:		return 429;		// 300
	case 6:		return 514;		// 400
	case 7:		return 355;		// 300
	case 8:		return 443;		// 400
	case 9:		return 419;		// 400
	case 10:	return 300;		// 300
	case 11:	return 576;		// 600
	case 12:	return 463;		// 500
	case 13:	return 449;		// 500
	case 14:	return 436;		// 500
	case 15:	return 425;		// 500
	case 16:	return 415;		// 500
	default:	return LONG_MAX;
	}
}

long lzi_arcsine_newtons()
{
	return gArcSineNewtons;
}

// For the tests: the most places past its width of the seed as converted
//	-- the larger of its text's places and the first stage -- or of the
//	iterate or s and c in any stage but the last, since the last read, this
//	thread's, which the read sets back -- nought or less, LONG_MIN for none...
static thread_local long gAngleNewtonOverrun = LONG_MIN;

long lzi_angle_newton_overrun()
{
	long	l_overrun = gAngleNewtonOverrun;

	gAngleNewtonOverrun = LONG_MIN;

	return l_overrun;
}

static void lzi_note_angle_overrun( long lOver )
{
	if ( lOver > gAngleNewtonOverrun )
		gAngleNewtonOverrun = lOver;
}

static void lzi_angle_by_newton( CLZNumber* pS, CLZNumber* pC, CLZNumber* pTarget )
{
	gArcSineNewtons++;

	long	l_base = pS->m_lBase;
	long	l_width = glWorkingPrecision;
	double	d_s = 0.0;
	double	d_c = 0.0;

	lzi_value_estimate( pS, &d_s );
	lzi_value_estimate( pC, &d_c );

	// The widths, the full one last, each a third of the next and two more...
	std::vector<long>	widths;
	long				l_seed = static_cast<long>( static_cast<double>( kAngleNewtonSeedDigits ) *
														std::log( 10.0 ) / std::log( static_cast<double>( l_base ) ) );

	for ( long w = l_width; ; w = (w / 3) + 2 )
	{
		widths.push_back( w );

		if ( ((w / 3) + 2) <= l_seed )
			break;
	}

	// From a double, by its text -- Assign( double ) writes "%f" -- converted
	//	at the text's own places, not the full width's, then to the first
	//	stage. The conversion writes the precision's worth, and Assign
	//	resets the precision, so it is set between the two...
	char		sz_seed[ 40 ];
	CLZNumber	y;

	::snprintf( sz_seed, sizeof( sz_seed ), "%.17f", std::atan2( d_s, d_c ) );
	y.SetBase( 10, false );
	y.Assign( sz_seed );
	y.SetPrecision( kAngleNewtonSeedDigits + 4 );
	y.SetBase( l_base, true );

	// As converted, past the larger of its text's places and the first stage...
	lzi_note_angle_overrun( y.m_lRadix - std::max( widths.back(), CLZNumber::EquivalentPlaces( kAngleNewtonSeedDigits + 4, 10, l_base ) ) );

	y.SetPrecision( widths.back() );
	y.Round();

	CLZNumber	s_narrowed, c_narrowed;

	for ( size_t i = widths.size(); i-- > 0; )
	{
		glWorkingPrecision = widths[ i ];

		CLZNumber	sine, cosine, s_cos, c_sin, c_cos, s_sin, numerator, denominator, step;
		CLZNumber*	p_s = pS;
		CLZNumber*	p_c = pC;

		// s and c to the stage and two, but at the last...
		if ( i > 0 )
		{
			if ( pS->m_lRadix > (widths[ i ] + 2) )
			{
				s_narrowed.Assign( pS );
				s_narrowed.Round( kDefaultRoundingMode, widths[ i ] + 2 );
				p_s = &s_narrowed;
			}

			if ( pC->m_lRadix > (widths[ i ] + 2) )
			{
				c_narrowed.Assign( pC );
				c_narrowed.Round( kDefaultRoundingMode, widths[ i ] + 2 );
				p_c = &c_narrowed;
			}
		}

		lzi_sincos_split_work( &y, &sine, &cosine, 0 );

		lz_multiply( p_s, &cosine, &s_cos );
		lz_multiply( p_c, &sine, &c_sin );
		lz_multiply( p_c, &cosine, &c_cos );
		lz_multiply( p_s, &sine, &s_sin );
		lz_subtract( &s_cos, &c_sin, &numerator );
		lz_add( &c_cos, &s_sin, &denominator );

		step.SetBase( l_base, false );
		step.SetPrecision( widths[ i ] );
		lz_divide( &numerator, &denominator, &step );
		lz_add( &y, &step );

		// ...and the iterate to the stage, but at the last...
		if ( i > 0 )
		{
			y.SetPrecision( widths[ i ] );
			y.Round();

			lzi_note_angle_overrun( y.m_lRadix - widths[ i ] );
			lzi_note_angle_overrun( p_s->m_lRadix - (widths[ i ] + 2) );
			lzi_note_angle_overrun( p_c->m_lRadix - (widths[ i ] + 2) );
		}
	}

	glWorkingPrecision = l_width;

	if ( !(pTarget->IsConstant()) )
		pTarget->Assign( &y );
}


// The exponentials this thread has taken, for the tests: sinh, cosh,
//	tanh and coth take one -- see lzi_exponential_count...
static thread_local long gExponentials = 0;

long lzi_exponential_count()
{
	return gExponentials;
}

long lz_exp_e( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	gExponentials++;

	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_exp_e, kRetakeDoubtExp );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	// Refused past the limit, nought below the width, before anything is
	//	worked out -- see lzi_answer_size...
	double	d_exponent = 0.0;
	long	l_beyond = kNoError;

	if ( lzi_value_or_infinity( leftOperand, &d_exponent ) &&
		 lzi_answer_beyond( lzi_answer_size( d_exponent, leftOperand->m_lBase ),
							leftOperand, targetNumber, &l_beyond ) )
		return l_beyond;

	long	lResult = kNoError;

	CLZNumber*	series_fract = nullptr;
	CLZNumber*	series_sum = nullptr;

	CLZNumber*	integer_base = nullptr;
	CLZNumber*	fractional_base = nullptr;

	CLZNumber*	result = nullptr;

	/*
		Wider than the answer throughout, rounded back at the end --
		see kExpGuardDigits. The raise comes before the working values
		exist, a CLZNumber taking its precision at construction.

		It has to wrap this whole function and not merely lz_raise_e,
		even though lz_raise_e is where the digits are actually lost.
		The integer and fractional parts are multiplied together at the
		end, so a fractional part computed at the narrow width caps the
		product however good the integer part was: exp(50.5) reaches 47
		digits with a guard here and 37 with it only in lz_raise_e.

		Four is enough, though, and not the sixteen lz_raise_e wants.
		The two nest, so anything more here runs the squaring chain
		wider for nothing -- see kExpGuardDigits.
	*/

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kExpGuardDigits, leftOperand->m_lBase );

	// Create working objects. series_sum and result are made further
	//	down, after the width the fractional part needs is known...
	series_fract = new CLZNumber;

	integer_base = new CLZNumber;
	fractional_base = new CLZNumber;

	// Set bases...
	series_fract->SetBase( leftOperand->m_lBase );

	// Initial values...
	integer_base->Assign( leftOperand );
	fractional_base->Assign( leftOperand );

	integer_base->StripFraction();
	fractional_base->StripInteger();

	/*
		**The integer part is taken first, because its size decides how
		wide the fractional part has to be.**

		The two are multiplied together at the end, so the product can
		be no better, relatively, than the worse of them -- and the
		fractional factor sits between one and e, where the working
		width buys it about that many significant digits and no more.
		Whatever lz_raise_e manages for the integer part, a fractional
		part computed at the narrow width caps the answer there.

		exp(13.94) is the case: e^13 is right to eighty digits now, and
		the answer still stopped at 38 where 36 places at 1.1e6 has
		room for 43, because e^0.94 had only been worked to 40 places.

		So the answer's magnitude is read off the integer part -- which
		carries essentially all of it, the fractional factor being
		worth at most one digit -- and the fractional series is then
		worked that much wider. The values it needs are made after the
		raise, a CLZNumber taking its precision at construction and not
		looking again.
	*/

	if ( !(integer_base->IsZero()) )
		lz_raise_e( integer_base, series_fract );
	else
		series_fract->Assign( static_cast<long>(1) );

	long	l_magnitude = series_fract->GetOrderOfMagnitude();

	/*
		**From two integer places, not from past the guard.** This waited
		for the answer to outgrow the guard, on the reasoning that the
		guard already covered anything smaller. It covered nothing: the
		guard is for the series and the product, and an answer of two to
		four integer places spent it on its own size. At the default, 40
		and 43 of 300 random arguments between 6.9 and 11.5 missed correct
		rounding, by up to 1.79 units; the audit never samples that band
		with a fraction. GetOrderOfMagnitude counts integer places, so one
		is anything under the base, and that needs nothing.
	*/

	if ( l_magnitude > 1 )
	{
		glWorkingPrecision += l_magnitude;

		/*
			fractional_base was made before the raise and Assign
			carried the caller's precision into it, so anything taking
			its width from it would be pinned there -- the trap
			lzi_log_gamma fell into.
		*/

		fractional_base->SetPrecision( glWorkingPrecision );
		series_fract->SetPrecision( glWorkingPrecision );
	}

	/*
		Handle the fractional part.

		Horner over the reciprocal factorial table while the table can
		serve the width, and the recurrence past that -- the same
		either/or, for the same reason, as lzi_trig_series. See the
		note at its head for why the threshold allows for stacked
		guards.

		**The table's length caps this exactly as it capped the sine.**
		The dropped term is f^33/33! for a fractional part f, which is
		1.3e-47 at f = 0.5 and 1.5e-38 at f = 0.94 -- and those were
		precisely where exp(0.5) and exp(13.94) stopped, at 47 and 38
		digits, at every precision from 36 to 100. Whole arguments were
		unaffected, having no fractional part to sum, which is why
		exp(1), exp(20) and exp(50) all scaled with the precision and
		hid it.
	*/

	series_sum = new CLZNumber;
	result = new CLZNumber;

	series_sum->SetBase( leftOperand->m_lBase );

	if ( !(fractional_base->IsZero()) )
	{
		/*
			**Halve the fraction into the table's reach, then square
			the answer back.**

			The table stops at f^32/32!, so what it can give depends
			entirely on how large f is: 70 digits at 0.1, 47 at 0.5,
			and only 38 as f approaches one. Each halving divides the
			dropped term by 2^33 and buys **exactly ten digits** --
			measured at 38, 48, 58, 68, 78 for f = 0.94 -- for one
			division and one squaring.

			That is much the cheaper way to widen this. Summing the
			terms from the recurrence instead needs no table and has no
			ceiling either, but a divide by a small integer measures
			4.5 times a tabled multiply and there are forty of them:
			exp(13.94) went to 8.9 times its cost that way, against
			1.0 with a single halving. The recurrence this replaces
			was put in one session earlier for the same fault, before
			the magnitude widening made the cost visible.

			Squaring back multiplies the relative error by two each
			time, so seven halvings cost about two digits -- which the
			guard covers, and the halvings themselves buy seventy.
		*/

		/*
			**The reaches and the threshold are decimal digits, and the
			working width is in places of the fraction's own base**, so
			it is converted first: 48 hex places are 58 digits, and
			compared as places they took the table with its halvings
			sized for 48.
		*/

		long	l_digits = CLZNumber::EquivalentPlaces( glWorkingPrecision, fractional_base->m_lBase, 10 );

		if ( l_digits <= kExpFractionThreshold )
		{
			/*
				What the table already reaches depends on how large f
				is, so ask before halving. At or under a half the
				dropped term is 8.1e-48 and the table gives 47 digits
				with no help at all; above it, as f approaches one, it
				gives 38. Sizing every call from the worst case made
				exp(0.5) pay 1.4 times for a halving it did not need.
			*/

			CLZNumber	two;
			CLZNumber	half;
			CLZNumber	squared;

			two.Assign( 2L );
			half.Assign( "0.5" );

			long	l_reach = (fractional_base->CompareABS( &half ) != kGreaterThan)
							  ? kExpTableReachHalf : kExpTableReach;

			long	l_halvings = 0;

			while ( ((l_reach + (kExpHalvingGain * l_halvings)) < l_digits) &&
					(l_halvings < kExpMaxHalvings) )
			{
				l_halvings++;
			}

			/*
				**At the working width, as lzi_exp_fraction_computed's
				halving is.** Assign gave the fraction the operand's
				precision, and an in-place division works to its left
				operand's, so until 2026-10-03 each halving cut it there,
				wherever a halving ran. In an even base that cost an
				ordinary operand nothing, halving one needing exactly the
				one place past its precision that a division keeps. In an
				odd base a half never terminates, and the cut fell inside
				the guard: at 36 places e^x itself misrounded 15 of 1,500
				random operands in bases nine, eleven and thirteen, and
				the hyperbolic functions, which work wider, 21 of 4,800
				in base seven -- by up to 0.07 of a unit past the half.
				And an operand written wider than the working width lost
				its tail: 72 of 5,400 e^x in base ten, by up to 0.14,
				every one with its fraction past a half, the only
				fractions halved at that width.
			*/

			fractional_base->SetPrecision( glWorkingPrecision );

			for ( long i = 0; i < l_halvings; i++ )
				lz_divide( fractional_base, &two );

			series_sum->Assign( static_cast<long>(0) );

			// Taylor series for the fractional part...
			for ( short i=32; i>0; i-- )
			{
				lz_add( series_sum, gRecipFact[ i ] );
				lz_multiply( series_sum, fractional_base );
			}

			// One last add...
			lz_add( series_sum, gRecipFact[ 0 ] );	// Done here to avoid the multiply
													//  operation at loop index 0

			for ( long i = 0; i < l_halvings; i++ )
			{
				lz_multiply( series_sum, series_sum, &squared );
				series_sum->Assign( &squared );
			}
		}

		// Past the base's floor, by binary splitting, under this function's
		//	guard and not its own -- see lzi_exp_split_floor and
		//	lzi_exp_split_work...
		else if ( glWorkingPrecision >= lzi_exp_split_floor( fractional_base->m_lBase ) )
			lzi_exp_split_work( fractional_base, series_sum, 0, 0 );

		else
		{
			/*
				Past what any halving can reach, because the limit is
				no longer the table's length but its entries' width:
				gRecipFact carries 50 significant digits, and halving
				an argument cannot add any. Seven halvings at 100
				places bought the length back and the answer still
				stopped at 55.

				So the terms are computed instead, forward from
				term(k) = term(k-1) * f / k, dividing by exact small
				integers with neither a length nor a width behind them.
				It costs 4.5 times a tabled multiply a term, which is
				why it is not simply how this always works.

				**term takes the fraction's base before anything else.**
				Made fresh it is decimal, and as the left operand it
				kept every product decimal too -- so in hex the whole
				recurrence ran in base ten, at the working width counted
				as decimal places, converting the fraction in and each
				term out again. That was 51 digits of a 53-digit answer
				at 44 hex places, and 5.4 times the time. The divisor
				may stay decimal, as the sine's does: an integer
				converts exactly.
			*/

			/*
				Halved first past a hundred digits, and worked wider for
				the squaring back -- see kExpComputedHalvings. Below it
				there are no halvings and no extra places.

				**Only where two divides the base**, which is where halving
				is exact and the argument stays short. In an odd base it is
				neither. A division in place keeps every place its dividend
				had and the precision's worth more, so in base seven
				sixty-four halvings at 1,300 places left 39,871 places and
				exp(3/7) took 15.5 s; rounded back each time, the argument
				still filled the width, every term became a full-width
				multiply, and it took 855 ms where the plain series takes
				17. So an odd base sums the series as it stands, which the
				bound keeps right at any width -- see lzi_series_bound.
			*/

			/*
				**Split a long fraction into a short head and a tiny
				tail**, wherever the series is computed: e^f = e^a e^b, a
				the first places of f and b the rest. Every term of the series
				multiplies by its argument, and a fraction as long as the
				width -- as the logarithm's Newton step hands it -- made
				every term a full-width multiply: at 9,000 places some
				470 of them. The head is short, so its terms are cheap;
				the tail is under a unit of the head's last place, so its
				series ends in about the width over the head's length --
				thirty terms or so. A fraction no longer than the head,
				which is every one written by hand, has no tail and goes
				exactly as it did.

				**It does not wait for the halvings.** Gated on them while
				they began at a thousand digits, it left Newton's steps
				under that summing full-width terms: exp(ln 2) at 960
				places took 12.6 ms, and 2.4 at 1,008 where both began.
			*/

			long	l_headPlaces = glWorkingPrecision / kExpSplitShare;

			if ( l_headPlaces < kExpSplitLeast )
				l_headPlaces = kExpSplitLeast;

			if ( fractional_base->m_lRadix > l_headPlaces )
			{
				CLZNumber	head;
				CLZNumber	tail;
				CLZNumber	tail_sum;
				CLZNumber	product;

				head.Assign( fractional_base );
				head.Round( RM_round_down, l_headPlaces );
				lz_subtract( fractional_base, &head, &tail );

				lzi_exp_fraction_computed( &head, l_digits, true, series_sum );
				lzi_exp_fraction_computed( &tail, l_digits, false, &tail_sum );

				lz_multiply( series_sum, &tail_sum, &product );
				series_sum->Assign( &product );
			}
			else
				lzi_exp_fraction_computed( fractional_base, l_digits, true, series_sum );
		}
	}
	else
		series_sum->Assign( static_cast<long>(1) );

	// The integer part was taken above, before the width was settled.
	//	Calculate the result...
	lz_multiply( series_sum, series_fract, result );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete series_fract;
    delete series_sum;
    delete integer_base;
    delete fractional_base;
    delete result;

	return lResult;
}


static long lzi_exp_n_work( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	/*
		lz_exp_n has settled the edges before this -- see
		lzi_power_edges -- so a negative base arrives here only with a
		whole exponent. Its answer is real, but the logarithm below has
		none for it, so the squaring takes it, sign and all. Before
		this, (-8)^3 was refused.
	*/

	if ( leftOperand->IsNegative() )
		return lzi_raise_n( leftOperand, rightOperand, targetNumber );

	/*
		**A whole exponent takes no logarithm either.** The fraction's
		power below is e to the nought times ln x, which is one, but the
		logarithm was worked out first at the full width -- and that
		width grows with the answer, so 1.2345 to the millionth spent 17
		of its 17 seconds on a logarithm it multiplied by nought. lz_power
		has always sent a whole exponent straight to the squaring. The
		answer is the same, digit for digit: the squaring's, times one.
	*/

	if ( rightOperand->IsInteger() )
		return lzi_raise_n( leftOperand, rightOperand, targetNumber );

	long	lResult = kNoError;

	CLZNumber*	natural_log = nullptr;
	CLZNumber*	integer_base = nullptr;
	CLZNumber*	fractional_base = nullptr;
	CLZNumber*	result = nullptr;

	// Create working objects...
	natural_log = new CLZNumber;
	integer_base = new CLZNumber;
	fractional_base = new CLZNumber;
	result = new CLZNumber;

	// Set bases...
	natural_log->SetBase( leftOperand->m_lBase );

	// Initial values...
	integer_base->Assign( rightOperand );
	fractional_base->Assign( rightOperand );

	integer_base->StripFraction();
	fractional_base->StripInteger();

	// Raise the fractional part...
	lResult = lzi_first_failure( lResult, lz_log_e( leftOperand, natural_log ) );

	lz_multiply( natural_log, fractional_base );
	lz_exp_e( natural_log, fractional_base );

	// Raise the integer part...
	lzi_raise_n( leftOperand, integer_base, integer_base );

	// Calculate the result...
	lz_multiply( integer_base, fractional_base, result );

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete natural_log;
    delete integer_base;
    delete fractional_base;
    delete result;

	return lResult;
}


/*
	**A root of whole order, its rounding decided exactly when the retake
	is in doubt** -- the open item of 2026-10-12, done 2026-10-14.

	**The exact tie it began as.** 8^(-1/3) is a half, and in base three,
	five or seven a half never ends, so at every width the root lies on the
	half itself and the logarithm road rounded it whichever way its error
	fell: at f9eff0f 56 of 126 such roots built in seven odd bases came out
	a unit wrong, and 42 of 84 each of power and exp_n at -1/Q. Only a root
	of negative whole order can be one there: a tie is c / (2 b^j), c odd,
	and its p-th power for a positive p has 2^p below the line, which no
	number finite in an odd base has. Nor can any other exponent's: x^(-P/Q)
	on a tie makes x^P = (2 b^j / c)^Q, whose twos number Q, while x^P's are
	P times x's -- so P divides Q, and P is one. An even base's half is a
	digit, and a tie there was already spelled and rounded exactly.

	**And what it became.** The comparison that tells a tie tells the side
	of one: for y, the root cut at lPlaces, and u a unit there, the root of
	order p is past the half when 2^p x > (2y + u)^p, and of order -p when
	x (2y + u)^p < 2^p -- whole numbers, as lzi_root_decide compares for the
	square and cube roots. So a root of whole order in doubt is decided by
	the retake -- at its first round in doubt while that is cheap, at the
	cap otherwise; see lzi_take_again_wider -- tie or not, in any base: an
	operand a hair from a tie's lies as deep as it likes, past any width a
	retake could reach -- one built a unit of the 245th place off, order
	-121 at ten places in base eleven, lay 3.6e-30 of a unit from the half,
	and the retake stops fifteen places past the width. Its numbers are
	bounded by lMostDigits, past which it declines and the retake goes on as
	before. Since 2026-10-15 power and exp_n bring any exponent P/Q here, as
	the root of order Q of x^P -- see lzi_power_decide.
*/

static bool lzi_whole_order_root_decide( CLZNumber* pResult, const CLZNumber* pMagnitude, long lOrder, long lPlaces, long lMostDigits )
{
	long	l_base = pResult->m_lBase;
	long	l_power = (lOrder < 0) ? -lOrder : lOrder;

	if ( (l_power < 2) || (lPlaces < 0) )
		return false;

	long	l_top = lzi_leading_exponent( pResult );
	long	l_digits = lPlaces + 2 + ((l_top > 0) ? l_top : 0);

	if ( (l_power > (lMostDigits / l_digits)) || (pMagnitude->m_lDigitsLength > (lMostDigits - (l_power * l_digits))) )
		return false;

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = lPlaces + 2;

	// y, |pResult| cut at lPlaces, and 2y + u...
	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	twice;
	CLZNumber	twice_power;
	CLZNumber	two;
	CLZNumber	two_power;
	CLZNumber	scaled;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &twice );
	lz_add( &twice, &unit );

	lzi_exact_power( &twice, l_power, &twice_power );

	two.SetBase( l_base, false );
	two.Assign( 2L );
	lzi_exact_power( &two, l_power, &two_power );

	// 2^p x against (2y + u)^p, or x (2y + u)^p against 2^p, exactly: a
	//	product's places are its operands' together...
	long	l_against;

	if ( lOrder > 0 )
	{
		glWorkingPrecision = pMagnitude->m_lRadix + 2;
		lz_multiply( const_cast<CLZNumber*>( pMagnitude ), &two_power, &scaled );
		l_against = lz_compare( &scaled, &twice_power );
	}
	else
	{
		glWorkingPrecision = twice_power.m_lRadix + pMagnitude->m_lRadix + 2;
		lz_multiply( const_cast<CLZNumber*>( pMagnitude ), &twice_power, &scaled );
		l_against = lz_compare( &two_power, &scaled );
	}

	bool	b_up = (l_against == kGreaterThan);

	// Exactly on the half: to the even kept digit...
	if ( l_against == kEqualTo )
	{
		long	l_index = y.m_lRadix - lPlaces;
		long	l_kept = ((l_index >= 0) && (l_index < y.m_lDigitsLength)) ? y.m_lpDigits[ l_index ] : 0;

		b_up = (l_kept % 2) != 0;
	}

	glWorkingPrecision = lPlaces + 2;

	if ( b_up )
		lz_add( &y, &unit );

	y.m_lSign = pResult->m_lSign;
	pResult->Assign( &y );
	pResult->SetPrecision( lPlaces );

	glWorkingPrecision = l_saved;

	return true;
}

// root_n's: an order that is a whole number...
static bool lzi_root_n_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* pRight, long lPlaces, long lMostDigits )
{
	long	l_power = 0;

	pRight->CheckForZero();

	if ( !(pRight->IsInteger()) || !lzi_whole_as_long( pRight, &l_power ) )
		return false;

	CLZNumber	magnitude;

	lz_abs( pLeft, &magnitude );

	return lzi_whole_order_root_decide( pResult, &magnitude, pRight->IsNegative() ? -l_power : l_power, lPlaces, lMostDigits );
}

// An exponent as a fraction P/Q in its lowest terms, exactly: Q is the
//	least whole number made of the base's primes that makes the exponent
//	times Q whole, tried in order up to lMostQ, and P that product. Every
//	fraction a base can write ends in a denominator of its own primes, so
//	none is missed below lMostQ; answers false past it, or for a P too big
//	for a long...
static bool lzi_exponent_as_fraction( CLZNumber* pExponent, long lMostQ, long* plP, long* plQ )
{
	long				l_base = pExponent->m_lBase;
	std::vector<long>	v_primes;

	for ( long l_rest = l_base, l_prime = 2; l_rest > 1; l_prime++ )
	{
		if ( (l_rest % l_prime) != 0 )
			continue;

		v_primes.push_back( l_prime );

		while ( (l_rest % l_prime) == 0 )
			l_rest /= l_prime;
	}

	std::vector<long>	v_denominators( 1, 1L );

	for ( long l_prime : v_primes )
	{
		size_t	n_before = v_denominators.size();

		for ( size_t n = 0; n < n_before; n++ )
			for ( long l_q = v_denominators[ n ]; l_q <= (lMostQ / l_prime); )
			{
				l_q *= l_prime;
				v_denominators.push_back( l_q );
			}
	}

	std::sort( v_denominators.begin(), v_denominators.end() );

	long		l_saved = glWorkingPrecision;
	bool		b_found = false;
	CLZNumber	q;
	CLZNumber	product;

	glWorkingPrecision = pExponent->m_lRadix + 24;

	for ( long l_q : v_denominators )
	{
		q.SetBase( l_base, false );
		q.Assign( l_q );
		lz_multiply( pExponent, &q, &product );

		if ( !product.IsInteger() )
			continue;

		b_found = lzi_whole_as_long( &product, plP );
		*plQ = l_q;
		break;
	}

	glWorkingPrecision = l_saved;

	return b_found;
}

// x^y's, for power and exp_n: an exponent of P/Q, read exactly, Q two or
//	more and the base positive. x^(P/Q) is the root of order Q of x^P, and
//	x^P is raised exactly...
static bool lzi_power_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* pRight, long lPlaces, long lMostDigits )
{
	pRight->CheckForZero();

	if ( pLeft->IsNegative() || pRight->IsInteger() || (lPlaces < 0) )
		return false;

	long	l_p = 0;
	long	l_q = 0;

	if ( !lzi_exponent_as_fraction( pRight, lMostDigits / (lPlaces + 2), &l_p, &l_q ) || (l_p < 1) || (l_q < 2) )
		return false;

	CLZNumber	magnitude;

	magnitude.Assign( pLeft );
	magnitude.TrimTrailing();

	if ( l_p > (lMostDigits / ((magnitude.m_lDigitsLength > 0) ? magnitude.m_lDigitsLength : 1)) )
		return false;

	if ( l_p > 1 )
		lzi_exact_power( pLeft, l_p, &magnitude );

	return lzi_whole_order_root_decide( pResult, &magnitude, pRight->IsNegative() ? -l_q : l_q, lPlaces, lMostDigits );
}


/*
	**A logarithm's rounding decided exactly in doubt** -- since
	2026-10-15, item 5 of that day's review. log_a x is p/q exactly when
	a^p = x^q, and only such an answer can lie on a half: a tie at n
	places is c / (2 B^n), c odd, so its q is made of two and the base's
	own primes. In an odd base that half never ends, and the retake ran
	to its cap and rounded the tie by its error -- log_9 3 at one place
	answered 0.5 for 0.4, and 783 of 2,044 such ties in the review's
	sample came out wrong.

	**So the answer in doubt is tried against each such q in turn**,
	smallest first, up to what its own places can tell: where q times it
	lies within its error of a whole number p, a^p = x^q is put to whole
	numbers -- a and x each a whole number over a power of its base,
	which convert between bases exactly, so x may be in another -- and
	if it holds, p/q is rounded at lPlaces exactly, a tie to the even
	kept digit. An answer that is no such fraction is irrational, and
	rounds by its digits as before. Declines past lMostDigits, and the
	retake goes on as it would have.
*/

// A positive value as a whole number over a power of its base, both
//	brought into lBase -- exactly, being whole...
static void lzi_whole_over_base_power( CLZNumber* pValue, long lBase, CLZNumber* pWhole, CLZNumber* pScale )
{
	pWhole->Assign( pValue );
	pWhole->m_lSign = kPositive;
	pWhole->CheckForZero();
	pWhole->TrimTrailing();

	long	l_places = pWhole->m_lRadix;

	lzi_scale_by_base_power( pWhole, -l_places );
	pScale->SetBase( pValue->m_lBase, false );
	pScale->Assign( 1L );
	lzi_scale_by_base_power( pScale, -l_places );
	pWhole->SetBase( lBase, true );
	pScale->SetBase( lBase, true );
}

// log_n's, a the base and x the value; and log_b's, through the one below...
static bool lzi_log_decide( CLZNumber* pResult, CLZNumber* pBase, CLZNumber* pValue, long lPlaces, long lMostDigits )
{
	if ( (lPlaces < 0) || pBase->IsNegative() || pValue->IsNegative() ||
		 pBase->FlagIsOn( kfInfinity ) || pBase->FlagIsOn( kfQNaN ) ||
		 pValue->FlagIsOn( kfInfinity ) || pValue->FlagIsOn( kfQNaN ) )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_saved = glWorkingPrecision;

	CLZNumber	a_whole;
	CLZNumber	a_scale;
	CLZNumber	x_whole;
	CLZNumber	x_scale;

	lzi_whole_over_base_power( pBase, l_base, &a_whole, &a_scale );
	lzi_whole_over_base_power( pValue, l_base, &x_whole, &x_scale );

	// A base of nought or one has no logarithm, and nor does nought...
	if ( a_whole.IsZero() || x_whole.IsZero() || (lz_compare( &a_whole, &a_scale ) == kEqualTo) )
		return false;

	long	l_aDigits = a_whole.m_lDigitsLength + a_scale.m_lDigitsLength;
	long	l_xDigits = x_whole.m_lDigitsLength + x_scale.m_lDigitsLength;

	// The answer's places tell q times it to within a unit of the last
	//	place times q, a few units allowed for its error; q no larger than
	//	a quarter of that, nor than x^q's digits allow...
	CLZNumber	answer;

	answer.Assign( pResult );
	answer.m_lSign = kPositive;

	long	l_places = answer.m_lRadix;
	double	d_mostQ = std::pow( static_cast<double>( l_base ), static_cast<double>( l_places - 2 ) ) / 16.0;
	long	l_mostQ = lMostDigits / l_xDigits;

	if ( d_mostQ < static_cast<double>( l_mostQ ) )
		l_mostQ = static_cast<long>( d_mostQ );

	// ...and four times it inside a 32-bit long -- see "A 32-bit long"...
	if ( l_mostQ > 100000000L )
		l_mostQ = 100000000L;

	// Two and the base's primes, and every product of them to the most...
	std::vector<long>	v_primes( 1, 2L );

	for ( long l_rest = l_base, l_prime = 2; l_rest > 1; l_prime++ )
	{
		if ( (l_rest % l_prime) != 0 )
			continue;

		if ( l_prime != 2 )
			v_primes.push_back( l_prime );

		while ( (l_rest % l_prime) == 0 )
			l_rest /= l_prime;
	}

	std::vector<long>	v_denominators( 1, 1L );

	for ( long l_prime : v_primes )
	{
		size_t	n_before = v_denominators.size();

		for ( size_t n = 0; n < n_before; n++ )
			for ( long l_q = v_denominators[ n ]; l_q <= (l_mostQ / l_prime); )
			{
				l_q *= l_prime;
				v_denominators.push_back( l_q );
			}
	}

	std::sort( v_denominators.begin(), v_denominators.end() );

	long				l_p = 0;
	long				l_q = 0;
	std::vector<long>	v_refuted;

	glWorkingPrecision = l_places + 4;

	for ( long l_try : v_denominators )
	{
		if ( l_try < 2 )
			continue;

		/*
			**A q refuted refutes its multiples**: a^(kp) = x^(kq) holds only
			if a^p = x^q does, and kp/kq looks as whole as p/q. Without this
			a near miss of a half tried every even q, each power larger than
			the last -- log_9 of 3.00000001 at one place took 590 ms against
			125 microseconds before the decider...
		*/
		bool	b_refuted = false;

		for ( long l_no : v_refuted )
			if ( (l_try % l_no) == 0 )
				b_refuted = true;

		if ( b_refuted )
			continue;

		CLZNumber	q;
		CLZNumber	product;
		CLZNumber	nearest;
		CLZNumber	gap;
		CLZNumber	allowed;

		q.SetBase( l_base, false );
		q.Assign( l_try );
		lz_multiply( &answer, &q, &product );

		nearest.Assign( &product );
		nearest.Round( RM_round_half_even, 0 );
		lz_subtract( &product, &nearest, &gap );
		gap.m_lSign = kPositive;

		// Four units of the answer's last place, times q...
		allowed.SetBase( l_base, false );
		allowed.Assign( 4L * l_try );
		lzi_scale_by_base_power( &allowed, l_places );

		if ( lz_compare( &gap, &allowed ) == kGreaterThan )
			continue;

		long	l_whole = 0;

		nearest.CheckForZero();

		if ( nearest.IsZero() || !lzi_whole_as_long( &nearest, &l_whole ) || (l_whole > (lMostDigits / l_aDigits)) )
			break;

		// a^p = x^q in whole numbers: A^p Dx^q against X^q Da^p, and for
		//	a negative p Da^m Dx^q against X^q A^m...
		bool		b_negative = (pResult->m_lSign == kNegative);
		CLZNumber	a_left;
		CLZNumber	x_left;
		CLZNumber	x_right;
		CLZNumber	a_right;
		CLZNumber	left;
		CLZNumber	right;

		lzi_exact_power( b_negative ? &a_scale : &a_whole, l_whole, &a_left );
		lzi_exact_power( &x_scale, l_try, &x_left );
		lzi_exact_power( &x_whole, l_try, &x_right );
		lzi_exact_power( b_negative ? &a_whole : &a_scale, l_whole, &a_right );
		lz_multiply( &a_left, &x_left, &left );
		lz_multiply( &x_right, &a_right, &right );

		if ( lz_compare( &left, &right ) == kEqualTo )
		{
			l_p = b_negative ? -l_whole : l_whole;
			l_q = l_try;
			break;
		}

		v_refuted.push_back( l_try );
	}

	if ( l_q == 0 )
	{
		glWorkingPrecision = l_saved;
		return false;
	}

	// p/q at lPlaces exactly: N = |p| B^n, k the answer cut there as a
	//	whole, and r = N - k q. The answer lies within units of its last
	//	place of a half, so the cut is p/q's own and r lies in [0, q);
	//	were it not, the answer is not what it seemed, and the retake goes
	//	on...
	CLZNumber	q;
	CLZNumber	n;
	CLZNumber	k;
	CLZNumber	kq;
	CLZNumber	r;
	CLZNumber	one;

	q.SetBase( l_base, false );
	q.Assign( l_q );
	n.SetBase( l_base, false );
	n.Assign( (l_p < 0) ? -l_p : l_p );
	lzi_scale_by_base_power( &n, -lPlaces );
	one.SetBase( l_base, false );
	one.Assign( 1L );

	k.Assign( &answer );
	k.Round( RM_round_down, lPlaces );
	lzi_scale_by_base_power( &k, -lPlaces );
	k.TrimTrailing();

	lz_multiply( &k, &q, &kq );
	lz_subtract( &n, &kq, &r );

	r.CheckForZero();

	if ( (r.IsNegative() && !r.IsZero()) || (lz_compare( &r, &q ) != kLessThan) )
	{
		glWorkingPrecision = l_saved;
		return false;
	}

	// Past the half, or on it with an odd kept digit...
	CLZNumber	twice;

	lz_add( &r, &r, &twice );

	long	l_against = lz_compare( &twice, &q );
	long	l_kept = (k.m_lRadix < k.m_lDigitsLength) ? k.m_lpDigits[ k.m_lRadix ] : 0;

	if ( (l_against == kGreaterThan) || ((l_against == kEqualTo) && ((l_kept % 2) != 0)) )
		lz_add( &k, &one );

	lzi_scale_by_base_power( &k, lPlaces );
	k.m_lSign = (l_p < 0) ? kNegative : kPositive;
	k.CheckForZero();

	pResult->Assign( &k );
	pResult->SetPrecision( lPlaces );

	glWorkingPrecision = l_saved;

	return true;
}

// log_b's: the base the value is written in...
static bool lzi_log_b_decide( CLZNumber* pResult, CLZNumber* pLeft, CLZNumber* /* pRight */, long lPlaces, long lMostDigits )
{
	CLZNumber	base;

	base.SetBase( pLeft->m_lBase, false );
	base.Assign( pLeft->m_lBase );

	return lzi_log_decide( pResult, &base, pLeft, lPlaces, lMostDigits );
}

long lz_exp_n( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_exp_n, rightOperand, kRetakeDoubt );

	// ...and a root of whole order decided exactly in doubt -- see
	//	lzi_whole_order_root_decide...
	answer.DecideTies( lzi_power_decide );

	// The edges every x^y shares, answered before any work -- see
	//	lzi_power_edges...
	long	l_edge = kNoError;

	// An edge answers exactly, or flagged -- see lzi_write_to_width...
	if ( lzi_power_edges( leftOperand, rightOperand, targetNumber, false, &l_edge ) )
	{
		answer.Exact();
		return l_edge;
	}

	// x^y has about y ln|x| of logarithm...
	double	d_ln = 0.0;
	double	d_y = 0.0;
	long	l_magnitude = -1;

	if ( lzi_ln_estimate( leftOperand, &d_ln ) && lzi_value_or_infinity( rightOperand, &d_y ) )
	{
		long	l_beyond = kNoError;

		// Refused past the limit, nought below the width -- see
		//	lzi_answer_size...
		if ( lzi_answer_beyond( lzi_answer_size( d_y * d_ln, leftOperand->m_lBase ),
								leftOperand, targetNumber, &l_beyond ) )
			return l_beyond;

		l_magnitude = lzi_places_for( d_y * d_ln, leftOperand->m_lBase );
	}

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	// The operands as they came, for the check below; the answer may be
	//	written over the left...
	CLZNumber	base;
	CLZNumber	exponent;

	base.Assign( leftOperand );
	exponent.Assign( rightOperand );

	long	lResult = lzi_widened_by_answer( lzi_exp_n_work, leftOperand, rightOperand, targetNumber, l_magnitude );

	// A short answer is checked against the operands before it is called
	//	exact, and spelled without its noughts if it is -- see
	//	lzi_write_to_width...
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( lzi_short_of_width( p_answer ) && lzi_power_is_exactly( p_answer, &base, &exponent ) )
	{
		p_answer->TrimTrailing();
		answer.Exact();
	}

	return lResult;
}


long lz_factorial( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	/*
		Defined here on the non-negative integers, and refused
		everywhere else.

		The loop below counts down to nought by ones, so anything that
		does not reach nought by that route never stops. A negative
		counted away from zero forever; so did a fraction, stepping
		2.5, 1.5, .5, -.5 and past it. Both spun with no answer and no
		error, which is worse than either a wrong number or a refusal.

		The mathematics has answers this does not, and **lz_gamma
		provides them** -- the continuation to every argument but the
		poles at nought and the negative whole numbers.

		It stays a separate function on purpose. This refusal is part
		of a published entry point and something may be relying on it
		to catch a bad argument, so the wider domain is asked for by
		name rather than arriving unannounced.
	*/

	leftOperand->CheckForZero();

	if ( leftOperand->IsNegative() || !(leftOperand->IsInteger()) )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	// n! is Gamma(n + 1), refused past the limit before the product is
	//	begun -- see lzi_answer_size...
	double	d_n = 0.0;
	long	l_beyond = kNoError;

	if ( lzi_value_or_infinity( leftOperand, &d_n ) && (d_n > 2.0) &&
		 lzi_answer_beyond( lzi_answer_size( lzi_ln_gamma_estimate( d_n + 1.0 ), leftOperand->m_lBase ),
							leftOperand, targetNumber, &l_beyond ) )
		return l_beyond;

	CLZNumber*	index = nullptr;
	CLZNumber*	one	= nullptr;
	CLZNumber*	result = nullptr;

	index = new CLZNumber;
	index->Assign( leftOperand );

	one = new CLZNumber( static_cast<long>(1) );
	result = new CLZNumber( static_cast<long>(1) );

	// In the operand's base -- see CLZAnswerBase...
	one->SetBase( leftOperand->m_lBase );
	result->SetBase( leftOperand->m_lBase );

	/*
		**n! by halves**, since 2026-10-15, performance item 2 of that
		day's review: the product of 1 to n split as lzi_range_product
		splits it, where one factor at a time multiplied a growing result
		by a small number n times -- 25,000! took 2.76 s. Whole numbers,
		so exact either way. A count past a long keeps the old walk...
	*/
	long	l_count = 0;

	if ( lzi_whole_as_long( index, &l_count ) )
	{
		CLZNumber	zero;

		zero.SetBase( leftOperand->m_lBase, false );
		zero.Assign( 0L );
		lzi_range_product( &zero, 1, l_count, result );
	}
	else
	{
		while ( !(index->IsZero()) )
		{
			lz_multiply( result, index );
			lz_subtract( index, one );
		}
	}

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->TrimAll();
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->TrimAll();
			leftOperand->Round();
		}
	}

	// Clean up...
    delete index;
    delete one;
    delete result;

	return kNoError;
}




/*
	Stirling's term for B_2n, n = lIndex + 1, as a numerator over
	2n(2n-1) times its denominator, both whole and exact in lBase,
	made once and kept -- see kBernoulliNumerators.

	The denominator is kept whole rather than multiplied in per term:
	that way hex gamma at the default, where the exact rationals serve,
	was 4 to 8% slower than the table, and is 1 to 5% kept. It is built
	here rather than as a long because 60 * 59 times B_60's denominator
	is 2e11: past a long where a long is 32 bits.
*/

static void lzi_bernoulli_rational( long lIndex, long lBase, CLZNumber** ppNumerator, CLZNumber** ppDenominator )
{
	auto	numeratorKey = std::make_tuple( static_cast<const void*>( kBernoulliNumerators[ lIndex ] ), lBase, -1L, -1L );
	auto	numerator = gTabledInBase.find( numeratorKey );

	if ( numerator == gTabledInBase.end() )
	{
		auto	p_numerator = new CLZNumber;

		p_numerator->Assign( kBernoulliNumerators[ lIndex ] );		// base ten text
		p_numerator->SetBase( lBase, true );						// a whole number converts exactly

		numerator = gTabledInBase.emplace( numeratorKey, p_numerator ).first;
	}

	auto	denominatorKey = std::make_tuple( static_cast<const void*>( &kBernoulliDenominators[ lIndex ] ), lBase, -1L, -1L );
	auto	denominator = gTabledInBase.find( denominatorKey );

	if ( denominator == gTabledInBase.end() )
	{
		long		l_n = lIndex + 1;
		auto		p_denominator = new CLZNumber;
		CLZNumber	stirling;

		p_denominator->SetBase( lBase );
		p_denominator->Assign( kBernoulliDenominators[ lIndex ] );
		stirling.SetBase( lBase );
		stirling.Assign( 2 * l_n * ((2 * l_n) - 1) );
		lz_multiply( p_denominator, &stirling );					// whole, so exact

		denominator = gTabledInBase.emplace( denominatorKey, p_denominator ).first;
	}

	*ppNumerator = numerator->second;
	*ppDenominator = denominator->second;
}

/*
	**Stirling's term lTerm, as a whole numerator over a whole denominator,
	for any term**: the table's exact rationals to the thirtieth, and past
	it the tangent numbers, made in lBase as far as asked and kept.

	Brent and Harvey's recurrence builds T_1 .. T_n with nothing but
	multiplications by small whole numbers and additions:

		T_k = (k - 1) T_(k-1), then for each k from 2, for j from k to n,
		T_j = (j - k) T_(j-1) + (j - k + 2) T_j

	and B_2n = (-1)^(n-1) 2n T_n / (4^n (4^n - 1)), so the term for B_2n is
	(-1)^(n-1) T_n / ((2n - 1) 4^n (4^n - 1) w^(2n-1)). Checked against the
	Bernoulli numbers exactly, to n = 80.

	The recurrence is not incremental, so a longer run is made afresh, and a
	quarter longer than asked and eight more: lz_gamma works a large answer
	a second time, up to thirty digits wider, and that asks about eight
	terms more. Its cost grows as the square of the count, so half as long
	again made the first gamma at 1000 places 0.96 s where this takes less.
*/

static const char	kTangentNumerator[] = "Stirling numerator, tangent";
static const char	kTangentDenominator[] = "Stirling denominator, tangent";

static void lzi_stirling_rational( long lTerm, long lBase, CLZNumber** ppNumerator, CLZNumber** ppDenominator )
{
	if ( lTerm <= (kBernoulliMax + 1) )
	{
		lzi_bernoulli_rational( lTerm - 1, lBase, ppNumerator, ppDenominator );
		return;
	}

	auto	numeratorKey = std::make_tuple( static_cast<const void*>( kTangentNumerator ), lBase, lTerm, 0L );
	auto	denominatorKey = std::make_tuple( static_cast<const void*>( kTangentDenominator ), lBase, lTerm, 0L );
	auto	found = gTabledInBase.find( numeratorKey );

	if ( found == gTabledInBase.end() )
	{
		long					l_count = lTerm + (lTerm / 4) + 8;
		std::vector<CLZNumber>	tangent( l_count + 1 );
		CLZNumber				factor;
		CLZNumber				scratch;

		for ( auto& t_k : tangent )
			t_k.SetBase( lBase, false );

		factor.SetBase( lBase, false );

		// All whole numbers, so every product is exact at any width...
		tangent[ 1 ].Assign( 1L );

		for ( long k = 2; k <= l_count; k++ )
		{
			factor.Assign( k - 1 );
			lz_multiply( &tangent[ k - 1 ], &factor, &tangent[ k ] );
		}

		for ( long k = 2; k <= l_count; k++ )
		{
			for ( long j = k; j <= l_count; j++ )
			{
				factor.Assign( j - k );
				lz_multiply( &tangent[ j - 1 ], &factor, &scratch );
				factor.Assign( j - k + 2 );
				lz_multiply( &tangent[ j ], &factor );
				lz_add( &tangent[ j ], &scratch );
			}
		}

		CLZNumber	power;					// 4^n
		CLZNumber	four;
		CLZNumber	one;

		power.SetBase( lBase, false );
		power.Assign( 1L );
		four.SetBase( lBase, false );
		four.Assign( 4L );
		one.SetBase( lBase, false );
		one.Assign( 1L );

		for ( long n = 1; n <= l_count; n++ )
		{
			lz_multiply( &power, &four );

			auto	key = std::make_tuple( static_cast<const void*>( kTangentNumerator ), lBase, n, 0L );

			if ( (n <= (kBernoulliMax + 1)) || (gTabledInBase.find( key ) != gTabledInBase.end()) )
				continue;

			auto	p_numerator = new CLZNumber;
			auto	p_denominator = new CLZNumber;

			p_numerator->Assign( &tangent[ n ] );

			if ( (n % 2) == 0 )
				lz_negate( p_numerator );

			p_denominator->SetBase( lBase, false );
			lz_subtract( &power, &one, p_denominator );
			lz_multiply( p_denominator, &power );
			factor.Assign( (2 * n) - 1 );
			lz_multiply( p_denominator, &factor );

			gTabledInBase.emplace( key, p_numerator );
			gTabledInBase.emplace( std::make_tuple( static_cast<const void*>( kTangentDenominator ), lBase, n, 0L ), p_denominator );
		}

		found = gTabledInBase.find( numeratorKey );
	}

	*ppNumerator = found->second;
	*ppDenominator = gTabledInBase.find( denominatorKey )->second;
}

/*
	ln(2*pi)/2 in lBase to lPlaces: the text within its 120 places, as
	before, and past them half the logarithm of 2*pi computed at the
	width -- see lzi_pi_family_in_base -- kept once made.
*/

static const char	kComputedHalfLogTwoPi[] = "ln(2 pi)/2, computed";

static void lzi_half_log_two_pi( CLZNumber* pValue, long lBase, long lPlaces )
{
	static const long	sl_constantPlaces = lzi_text_places( kHalfLogTwoPi );

	if ( lPlaces <= CLZNumber::EquivalentPlaces( sl_constantPlaces, 10, lBase ) )
	{
		lzi_constant_in_base( pValue, kHalfLogTwoPi, lBase, lPlaces );
		return;
	}

	auto	key = std::make_tuple( static_cast<const void*>( kComputedHalfLogTwoPi ), lBase, lPlaces, 0L );
	auto	found = gTabledInBase.find( key );

	if ( found == gTabledInBase.end() )
	{
		long		l_saved = glWorkingPrecision;
		auto		p_value = new CLZNumber;
		CLZNumber	twoPi( kTwoPi );
		CLZNumber	two;

		glWorkingPrecision = lPlaces + kComputedPiGuard;

		lzi_pi_family_in_base( &twoPi, kTwoPi, lBase, glWorkingPrecision );

		p_value->SetBase( lBase, false );
		p_value->SetPrecision( glWorkingPrecision );
		lz_log_e( &twoPi, p_value );

		two.SetBase( lBase, false );
		two.Assign( 2L );
		lz_divide( p_value, &two );

		p_value->SetPrecision( lPlaces );
		p_value->Round();

		glWorkingPrecision = l_saved;

		found = gTabledInBase.emplace( key, p_value ).first;
	}

	pValue->Assign( found->second );
}


/*
	Stirling's series for the logarithm of the gamma function, for a
	positive argument.

		lnGamma(w) = (w - 1/2) ln w  -  w  +  ln(2*pi)/2
		             + sum B(2n) / ( 2n (2n-1) w^(2n-1) )

	**The series is asymptotic, not convergent.** At a fixed argument
	its terms shrink for a while and then grow without bound, so past
	a certain point adding terms makes the answer worse. Accuracy comes
	instead from the argument: every term carries w^(2n-1) underneath,
	so a larger w makes all of them smaller at once.

	That is what the recurrence is for. lnGamma(z) = lnGamma(z + m) -
	ln( z (z+1) ... (z+m-1) ) moves the argument up for the price of m
	multiplications and one logarithm, and multiplications are the
	cheapest thing in this engine. So the shift buys the accuracy and
	the table only sets the price -- see the note at kBernoulliMax for
	what each table length costs in shift.

	The subtraction at the end is where the places go. For a small
	argument the two sides are close -- lnGamma(0.5) is about 0.57
	while both terms are near 60 -- so the guard digits are not
	optional here. See kGammaGuardDigits.
*/

static long lzi_log_gamma( CLZNumber* argument, CLZNumber* targetNumber )
{
	long	lResult = kNoError;

	/*
		Where the argument has to reach before the series is summed.

		**A target, not a count.** The series settles once its argument
		is large enough, so an argument that is already large needs no
		shift at all -- and shifting it anyway is not merely wasted, it
		is ruinous: the rising product grows like a factorial, so a
		thousand needless steps means multiplying thousand-digit
		numbers a thousand times. Written as a count first, this took
		Gamma(100.5) from milliseconds to not finishing.

		Wider precision needs a larger target; see kGammaBaseShift for
		the measurements behind the rate.

		**Those measurements are in decimal digits, and the working
		precision is in places of the argument's base**, so it is
		converted before it is compared: 44 places are 44 digits in base
		ten, 53 in base sixteen and 13 in base two. Compared as they
		stood, base sixteen asked for a series good to 44 digits when it
		needed 53.
	*/

	long	l_target = kGammaBaseShift;
	long	l_digits = CLZNumber::EquivalentPlaces( glWorkingPrecision, argument->m_lBase, 10 );

	/*
		**Past the default, as far as the first term left out needs.**
		The series stops at B_60, so what it leaves is about
		B_62 / (62 * 61 * w^61), and w is shifted until that falls a
		margin under the working width. |B_2n| is 2 (2n)! / (2 pi)^2n
		times zeta(2n), which at 62 is one to a part in 2^62, so the
		bound is computed rather than tabled. The rate this replaced
		grew half as much again every ten digits, which is 1.4 to 2.3
		times what this asks for from 50 digits to 200 -- 3079 at 150
		places where this asks 1445 -- and the rising product it paid
		for was most of a wide gamma's time once its logarithm stopped
		being quadratic. Never less than kGammaBaseShift, so the default
		is as it was.
	*/

	/*
		**How many terms: thirty to 120 digits, and a share of the width
		past them** -- kGammaTermsPer100Digits. Thirty terms make the
		shift grow tenfold every 61 digits, so past about 225 places a
		gamma took seconds and past 300 the shift ran to hundreds of
		thousands; see lzi_stirling_rational for the rest of the terms.
	*/

	long	l_terms = kBernoulliMax + 1;
	long	l_share = ((l_digits * kGammaTermsPer100Digits) + 99) / 100;

	if ( (l_digits > kBernoulliTableReach) && (l_share > l_terms) )
		l_terms = l_share;

	if ( l_digits > kGammaShiftFrom )
	{
		double	d_omitted = 2.0 * (l_terms + 1);
		double	d_logOmitted = (std::lgamma( d_omitted + 1.0 ) + std::log( 2.0 ) -
								d_omitted * std::log( 2.0 * std::acos( -1.0 ) ) -
								std::log( d_omitted * (d_omitted - 1.0) )) / std::log( 10.0 );
		double	d_logTarget = (static_cast<double>( l_digits + kGammaShiftMargin ) + d_logOmitted) / (d_omitted - 1.0);
		long	l_needed = static_cast<long>( std::ceil( std::pow( 10.0, d_logTarget ) ) );

		if ( l_needed > l_target )
			l_target = l_needed;
	}

	/*
		**A tiny argument costs places at the other end.**

		The rising product below is z(z+1)(z+2)..., and for a very
		small z every one of those factors contributes a correction the
		size of z itself: z*(1+z) is z + z^2, and z^2 sits twice as far
		down as z does. At 1e-30 the whole correction is z*H(m), about
		3.7e-30, and dropping it left lnGamma(1e-30) with 31 correct
		digits where 36 places at 69.08 affords 38 -- which is exactly
		what 3.7e-30 against 69 predicts.

		So the argument's own smallness is added to the width, the
		mirror of the magnitude widening lz_gamma does at the top end.
		It is the argument that decides here, not the answer, so no
		probe is needed: lzi_leading_exponent already says how far down
		the value starts.
	*/

	long	l_leading = lzi_leading_exponent( argument );
	long	l_savedPrecision = glWorkingPrecision;

	if ( l_leading < 0 )
		glWorkingPrecision -= l_leading;

	CLZNumber	product;
	CLZNumber	shifted;
	CLZNumber	step;
	CLZNumber	term;
	CLZNumber	power;
	CLZNumber	square;
	CLZNumber	scratch;
	CLZNumber	total;

	// The rising product z (z+1) ... (z+shift-1), and the shifted
	//	argument the series is summed at...
	// In the operand's base -- see CLZAnswerBase...
	product.SetBase( argument->m_lBase );
	step.SetBase( argument->m_lBase );

	product.Assign( 1L );
	shifted.Assign( argument );
	step.Assign( 1L );

	/*
		**Assign carries the source's precision, and the caller's
		number was not made at ours.**

		A CLZNumber takes its precision at construction and Assign
		copies it over along with the digits, so shifted came out of
		this at the argument's width -- the default 36 -- however wide
		lzi_log_gamma was running. Every operation that then took its
		width from shifted was pinned there too, and the one that
		mattered is the logarithm below: ln(36.5) came back with 37
		correct digits inside this function while the same call from
		outside, on a number made at the working precision, gave 51.

		Gamma(10.5) stopped at 37 digits for it. The symptom looked
		exactly like a series that had run out of terms, which is what
		sent me through the Bernoulli table twice before looking here.
	*/

	shifted.SetPrecision( glWorkingPrecision );

	CLZNumber	target;

	target.Assign( l_target );

	/*
		**The rising product by halves**, since 2026-10-15, performance
		item 2 of that day's review. It grows like a factorial -- 9,400
		digits at 1,000 places -- and was multiplied by one factor of the
		width at a time, 42 ms at 1,000 places; split as lzi_range_product
		splits it the halves stay balanced. Each multiplication rounds at
		the working width either way. The count of factors is the steps
		the walk took: the least k with z + k at the target...
	*/
	CLZNumber	gap;
	CLZNumber	target_here;
	long		l_factors = 0;

	// ...in the argument's base, where the gap is exact...
	target_here.SetBase( argument->m_lBase, false );
	target_here.Assign( l_target );
	lz_subtract( &target_here, &shifted, &gap );
	gap.CheckForZero();

	if ( !(gap.IsZero()) && !(gap.IsNegative()) )
	{
		gap.Round( RM_round_ceiling, 0 );

		if ( !lzi_whole_as_long( &gap, &l_factors ) )
			l_factors = -1;
	}

	if ( l_factors > 0 )
	{
		CLZNumber	count;

		lzi_range_product( &shifted, 0, l_factors - 1, &product );
		count.SetBase( argument->m_lBase, false );
		count.Assign( l_factors );
		lResult = lzi_first_failure( lResult, lz_add( &shifted, &count ) );
	}
	else if ( l_factors < 0 )
	{
		while ( shifted.Compare( &target ) == kLessThan )
		{
			lResult = lzi_first_failure( lResult, lz_multiply( &product, &shifted ) );
			lResult = lzi_first_failure( lResult, lz_add( &shifted, &step ) );
		}
	}

	// (w - 1/2) ln w - w + ln(2*pi)/2...
	CLZNumber	half;
	CLZNumber	logShifted;

	half.Assign( "0.5" );

	/*
		**ln w is worked as many places wider as w has whole places**,
		because w - 1/2 multiplies it next and carries its last place up
		by that many -- the rule lzi_log_scaled keeps for k ln(base).
		Rounded at the width, ln w was half a unit of the guard's last
		place out, and the product w/2 units. Until 2026-10-07 that came
		out of the guard:

		  - **as the width grew**, since the shift's target does: w is 24
		    at the default and about 2,940 at 1000 digits, where gamma and
		    log-gamma had used up three of their eight guard digits and
		    the near-half census found log-gamma rounded wrong;
		  - **for a large argument at any width**, which is not shifted
		    at all: w is the argument, so from about 1e8 at the default
		    nearly every log-gamma was wrong, by 4 units at 1e8 and
		    2.6e22 at 1e30.
	*/

	long	l_wholePlaces = lzi_leading_exponent( &shifted ) + 1;

	glWorkingPrecision += l_wholePlaces;
	lResult = lzi_first_failure( lResult, lz_log_e( &shifted, &logShifted ) );
	glWorkingPrecision -= l_wholePlaces;

	lResult = lzi_first_failure( lResult, lz_subtract( &shifted, &half, &scratch ) );
	lResult = lzi_first_failure( lResult, lz_multiply( &scratch, &logShifted, &total ) );
	lResult = lzi_first_failure( lResult, lz_subtract( &total, &shifted ) );

	CLZNumber	halfLogTwoPi( kHalfLogTwoPi );

	lzi_half_log_two_pi( &halfLogTwoPi, argument->m_lBase, glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_add( &total, &halfLogTwoPi ) );

	/*
		The Bernoulli terms, over the odd powers of w. power walks
		w, w^3, w^5, ... by one multiplication by w^2 a term, the way
		lzi_trig_series walks its own.

		The loop leaves early when a term has fallen below the working
		width, which is the ordinary case: at the default precision the
		terms run out well before the table does. It has to be able to
		leave, because carrying on past that point is where an
		asymptotic series starts undoing itself.
	*/

	lResult = lzi_first_failure( lResult, lz_square( &shifted, &square ) );

	power.Assign( &shifted );

	bool	b_exactBernoulli = (CLZNumber::EquivalentPlaces( glWorkingPrecision, argument->m_lBase, 10 ) > kBernoulliTableReach);

	for ( long n = 1; n <= l_terms; n++ )
	{
		/*
			**Past the table's fifty digits, the exact rational** -- see
			kBernoulliNumerators. The fifty-digit entries moved lnGamma by
			1e-53 to 7e-53 whatever the width, so every large gamma
			stopped at 52 significant digits in every base until
			2026-09-26 -- at the default too, where lz_gamma works wider
			for a large answer's size. Within the reach the table serves,
			as it always did.
		*/

		CLZNumber	denominator;
		CLZNumber*	p_denominator = &denominator;
		CLZNumber*	p_bernoulli = nullptr;

		if ( b_exactBernoulli )
			lzi_stirling_rational( n, argument->m_lBase, &p_bernoulli, &p_denominator );
		else
		{
			// Both on the left, so in the argument's base -- see
			//	CLZAnswerBase and lzi_table_in_base...
			denominator.SetBase( argument->m_lBase );
			denominator.Assign( 2 * n * ((2 * n) - 1) );

			p_bernoulli = lzi_table_in_base( gBernoulli[ n - 1 ], argument->m_lBase );
		}

		lResult = lzi_first_failure( lResult, lz_multiply( p_denominator, &power, &scratch ) );
		lResult = lzi_first_failure( lResult, lz_divide( p_bernoulli, &scratch, &term ) );

		term.CheckForZero();

		if ( term.IsZero() )
			break;

		lResult = lzi_first_failure( lResult, lz_add( &total, &term ) );

		lResult = lzi_first_failure( lResult, lz_multiply( &power, &square ) );
	}

	// ...less the logarithm of everything the shift multiplied in...
	CLZNumber	logProduct;

	lResult = lzi_first_failure( lResult, lz_log_e( &product, &logProduct ) );
	lResult = lzi_first_failure( lResult, lz_subtract( &total, &logProduct ) );

	/*
		The extra places for a small argument never leave this
		function. **They have to be given back here rather than by the
		caller**: lz_gamma calls this twice, once to measure the
		answer's magnitude and once for real, and a widening left in
		place compounds across the two -- the second call widened on
		top of the first and the whole thing stopped finishing.
	*/

	glWorkingPrecision = l_savedPrecision;

	targetNumber->Assign( &total );

	return lResult;
}


/*
	The logarithm of the gamma function, for a positive argument.

	Refused at nought and below rather than answering ln|Gamma|. Gamma
	alternates sign below nought, so a real logarithm of it does not
	exist there; the absolute value is a different function and
	pretending otherwise is the sort of quiet substitution this engine
	avoids elsewhere. lz_gamma answers the whole domain, and a caller
	who wants the log of a negative gamma can take it from there and
	handle the sign in the open.
*/

long lz_log_gamma( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_log_gamma, kRetakeDoubtLogGamma );

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	leftOperand->CheckForZero();

	if ( leftOperand->IsZero() || leftOperand->IsNegative() )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	long	lResult = kNoError;

	// Guarded across the shift and the subtraction that undoes it --
	//	see kGammaGuardDigits...
	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kGammaGuardDigits, leftOperand->m_lBase );

	CLZNumber	result;

	lResult = lzi_first_failure( lResult, lzi_log_gamma( leftOperand, &result ) );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


/*
	The gamma function, which continues the factorial to every argument
	but its poles: Gamma(n) = (n-1)! for a whole number.

	lz_factorial is unchanged and still refuses anything that is not a
	non-negative whole number. The continuation is a different function
	and is asked for by name; a caller relying on that refusal to catch
	a bad argument keeps it.

	Four cases:

	  - **A whole number one or more** is (n-1)! exactly, and goes to
	    lz_factorial rather than through the series. An exact answer is
	    worth more than a well-rounded one, and this is the same
	    reasoning that has lz_root_3 answer the cube root of 8 as 2.

	  - **Nought and the negative whole numbers are poles**, and the
	    refusal there is not the usual "outside the domain". Gamma runs
	    to plus infinity on one side of each and minus infinity on the
	    other, so there is no answer, not even an infinity with a sign
	    -- kUndefinedResult and a not-a-number rather than kfInfinity.

	  - **Anything else above nought** is exp of the logarithm.

	  - **Anything else below nought** comes from the reflection

	        Gamma(z) = pi / ( sin(pi z) Gamma(1 - z) )

	    which turns a negative argument into a positive one. Near nought
	    and near the poles sin(pi z) is small and the answer large, and
	    this said the accuracy lost there was honest, the answer being
	    that badly conditioned. It is conditioned badly only against a
	    change in z, and z is taken as exact here as everywhere: the
	    work is widened for both since 2026-09-26 -- see below.
*/

long lz_gamma( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_gamma, kRetakeDoubtGamma );

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	leftOperand->CheckForZero();

	// The poles: nought and every negative whole number...
	if ( leftOperand->IsInteger() && !(leftOperand->IsPositive()) )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kUndefinedResult );

	// A whole number one or more is a factorial, exactly...
	if ( leftOperand->IsInteger() )
	{
		CLZNumber	one;
		CLZNumber	lowered;

		one.Assign( 1L );

		lz_subtract( leftOperand, &one, &lowered );

		// ...and exactly -- see lzi_write_to_width. Written where the
		//	operation writes: in place it went into the copy until
		//	2026-10-15, and gamma(5) in place was left 5...
		answer.Exact();
		return lz_factorial( &lowered, (targetNumber != nullptr) ? targetNumber : leftOperand );
	}

	/*
		Sized before anything is worked out -- see lzi_answer_size. Past
		the limit, refused. A negative argument's answer is pi over
		sin(pi x) Gamma(1 - x), and sin(pi x) can be no smaller than the
		argument's own places allow, so that bounds the answer from
		above: nought if even the bound is two places under the width,
		and otherwise worked out as before. A bound from above says
		nothing about too large...
	*/
	double	d_x = 0.0;
	double	d_lnX = 0.0;
	long	l_beyond = kNoError;

	if ( lzi_value_or_infinity( leftOperand, &d_x ) && lzi_ln_estimate( leftOperand, &d_lnX ) )
	{
		double	d_lnBase = std::log( static_cast<double>( leftOperand->m_lBase ) );
		double	d_lnAnswer = 0.0;

		// Gamma(x) is about 1/x for a tiny x, which a double would round
		//	to nought...
		if ( d_x > 0.0 )
			d_lnAnswer = (d_lnX < -20.0) ? -d_lnX : lzi_ln_gamma_estimate( d_x );
		else
			d_lnAnswer = 1.144729885849400 + (static_cast<double>( leftOperand->m_lRadix + 1 ) * d_lnBase) -
						 lzi_ln_gamma_estimate( 1.0 - d_x );

		ELZAnswerSize	e_size = lzi_answer_size( d_lnAnswer, leftOperand->m_lBase );

		if ( (d_x < 0.0) && (e_size != AS_nought) )
			e_size = AS_fits;

		if ( lzi_answer_beyond( e_size, leftOperand, targetNumber, &l_beyond ) )
			return l_beyond;
	}

	long	lResult = kNoError;

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kGammaGuardDigits, leftOperand->m_lBase );

	CLZNumber	result;

	if ( leftOperand->IsPositive() )
	{
		lResult = lzi_first_failure( lResult, lzi_log_gamma( leftOperand, &result ) );
		lResult = lzi_first_failure( lResult, lz_exp_e( &result ) );

		/*
			**The logarithm's places are the answer's significant
			digits, whatever size the answer is**, and for a large
			gamma that is not nearly enough.

			exp turns an absolute error into a relative one: if
			lnGamma is right to 44 places then Gamma is right to 44
			significant digits, and no more, because 44 places of a
			number near fourteen is all the logarithm was carrying.
			Meanwhile 36 places of an answer near 1e6 has room for 43
			digits and 1e17 has room for 54. Gamma(20.5) came out at
			37 digits against the 54 its size affords.

			So the answer's own magnitude says how much wider the
			logarithm has to be, and the cheapest way to learn it is
			to work the thing out once and look. Only an argument
			whose gamma is large pays for the second pass, and it
			pays for it in the places it gets back.

			The same shape as lz_root_3's band: a magnitude that
			appears at the end has to be carried through the middle.
		*/

		long	l_magnitude = result.GetOrderOfMagnitude();

		/*
			Capped, and the cap is the same judgement the rounding
			audit already makes: past about 1e30 an answer has more
			integer digits than 36 places of fraction, and "correctly
			rounded to 36 places" stops being a reasonable thing to ask
			of a series. auditrounding.py sets those cases aside for
			the same reason. Without a cap, Gamma(100.5) asks for 159
			extra places and the shift that needs makes the rising
			product a number with thousands of digits.
		*/

		/*
			**The cap is thirty decimal places, and the magnitude is
			counted in places of the answer's own base**, so the cap is
			converted to that base before it is applied. Taken as it
			stood it allowed thirty binary places -- nine digits -- and
			a large gamma in base two stopped at twenty digits where
			base ten reached fifty-two.
		*/

		long	l_cap = CLZNumber::EquivalentPlaces( kGammaMaxWidening, 10, result.m_lBase );

		if ( l_magnitude > l_cap )
			l_magnitude = l_cap;

		if ( l_magnitude > 1 )
		{
			glWorkingPrecision += l_magnitude;

			CLZNumber	wider;

			lResult = lzi_first_failure( lResult, lzi_log_gamma( leftOperand, &wider ) );
			lResult = lzi_first_failure( lResult, lz_exp_e( &wider ) );

			result.Assign( &wider );
		}
	}

	else
	{
		/*
			Gamma(z) = pi / ( sin(pi z) Gamma(1 - z) ), which needs
			the gamma of a positive argument and gets it from the
			series directly rather than by calling back in.

			**Worked wider by two things the guard does not cover**,
			both known before any work is done, until 2026-09-26 neither
			taken: 64 of 179 negative arguments missed at the
			default, gamma(-1e-20) by 2e31 units, where gamma(1e-20) was
			right.

			  - **A small sine.** sin(pi z) comes back to so many places,
			    so its relative accuracy falls with its size, and near
			    nought or a pole it is small: as many places as it has
			    noughts. z's distance to the nearest whole number says
			    how small, exactly.
			  - **The answer's own size**, as for a positive argument:
			    capped the same way, past about 1e30.

			The argument is taken as exact, as the circular quotients
			near their poles take theirs.

			**z's own size needs nothing**, though pi's last place is
			multiplied by it: the guard absorbs that below about 1e5,
			and past it Gamma(1 - z) has grown so fast that the answer is
			under the width unless z sits within 1/Gamma(|z|) of a pole
			-- a distance of tens of thousands of places. Widening for it
			cost 1 to 3% on every z of 10 or more and moved no answer in
			the census, nor did setting the quotient's width, which
			rounds in places after the point already.
		*/

		long	l_base = leftOperand->m_lBase;
		long	l_widening = 0;

		{
			CLZNumber	fraction;
			CLZNumber	twice;
			CLZNumber	whole;
			double		d_lnDistance = 0.0;
			double		d_z = 0.0;

			// The distance to the nearest whole number, exactly: the
			//	fraction, or one less it past a half. Doubled rather than
			//	compared with a half, which an odd base cannot write...
			lz_abs( leftOperand, &fraction );
			fraction.StripInteger();

			whole.SetBase( l_base, false );
			whole.Assign( 1L );

			lz_add( &fraction, &fraction, &twice );

			if ( lz_compare( &twice, &whole ) == kGreaterThan )
				lz_subtract( &whole, &fraction, &fraction );

			if ( lzi_ln_estimate( &fraction, &d_lnDistance ) )
			{
				// ln |sin(pi z)|, which is ln(pi d) for a small distance d...
				double	d_lnSine = (d_lnDistance < -20.0)
									? (std::log( std::acos( -1.0 ) ) + d_lnDistance)
									: std::log( std::sin( std::acos( -1.0 ) * std::exp( d_lnDistance ) ) );
				double	d_lnBase = std::log( static_cast<double>( l_base ) );

				// The sine's noughts after the point...
				if ( -d_lnSine > d_lnBase )
					l_widening += static_cast<long>( -d_lnSine / d_lnBase );

				if ( lzi_value_estimate( leftOperand, &d_z ) )
				{
					// ...and the answer's whole places, as GetOrderOfMagnitude
					//	would count them...
					double	d_lnAnswer = std::log( std::acos( -1.0 ) ) - d_lnSine - std::lgamma( 1.0 - d_z );
					long	l_magnitude = (d_lnAnswer < 0.0) ? 0 : (static_cast<long>( d_lnAnswer / d_lnBase ) + 1);
					long	l_cap = CLZNumber::EquivalentPlaces( kGammaMaxWidening, 10, l_base );

					if ( l_magnitude > l_cap )
						l_magnitude = l_cap;

					// ...from two whole places, as for a positive argument...
					if ( l_magnitude > 1 )
						l_widening += l_magnitude;
				}
			}
		}

		glWorkingPrecision += l_widening;

		CLZNumber	one;
		CLZNumber	pi( kPi );
		CLZNumber	reflected;
		CLZNumber	sine;
		CLZNumber	logGamma;

		// Both on the left, so in the operand's base -- see CLZAnswerBase...
		one.SetBase( leftOperand->m_lBase );
		lzi_pi_family_in_base( &pi, kPi, leftOperand->m_lBase, glWorkingPrecision );

		one.Assign( 1L );

		lResult = lzi_first_failure( lResult, lz_subtract( &one, leftOperand, &reflected ) );

		lResult = lzi_first_failure( lResult, lzi_log_gamma( &reflected, &logGamma ) );
		lResult = lzi_first_failure( lResult, lz_exp_e( &logGamma ) );

		lResult = lzi_first_failure( lResult, lz_multiply( &pi, leftOperand, &sine ) );
		lResult = lzi_first_failure( lResult, lz_sine( &sine ) );

		lResult = lzi_first_failure( lResult, lz_multiply( &sine, &logGamma ) );
		lResult = lzi_first_failure( lResult, lz_divide( &pi, &sine, &result ) );
	}

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result.SetPrecision( glWorkingPrecision );
	result.Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}

long lz_fraction( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( leftOperand );
			targetNumber->StripInteger();
		}
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->StripInteger();

	return kNoError;
}


/*
	**sinh and cosh from one exponential** -- since 2026-10-15,
	performance item 5 of that day's review. e^-|x| is the reciprocal of
	e^|x|, a division where it was a second exponential, and tanh and coth
	take both from the one: tanh at 300 places cost five and a half
	exponentials, and one exponential and a division gave the same digits
	three to seven times faster. |x|, so the reciprocal is of the larger:
	e^x of a large negative x rounds to nought at the width.

	Worked at the width given. **A small argument widens the sine by its
	own noughts**, e^x - e^-x cancelling to 2x + x^3/3 + ...: worked to W
	places it held only W - k significant digits of an answer about
	base^-k, where the rule for a quotient writes W of them --
	sinh(-0.00000000002332) in base four at three digits answered
	-0.000000000100 until 2026-10-15. **Where x^3/6 lies past the width
	beside x, the sine is x and the cosine one**, taken as they stand
	rather than widened by noughts past any width the answer needs; an x
	near a half of its own places is left to the retake and decided by its
	side of x -- see lzi_beside_x_decide. Either answer may be wanted
	alone.
*/

static long lzi_sinh_cosh( CLZNumber* pX, CLZNumber* pSine, CLZNumber* pCosine )
{
	long	l_saved = glWorkingPrecision;
	long	l_lead = lzi_leading_exponent( pX );
	long	l_small = (l_lead < 0) ? -l_lead : 0;

	pX->CheckForZero();

	if ( (pSine != nullptr) && (pX->IsZero() || ((l_small > 0) && ((2 * l_small) >= (l_saved + 2)))) )
	{
		pSine->Assign( pX );

		if ( pCosine != nullptr )
		{
			pCosine->SetBase( pX->m_lBase, false );
			pCosine->Assign( 1L );
		}

		return kNoError;
	}

	if ( pSine != nullptr )
		glWorkingPrecision += l_small;

	CLZNumber	magnitude;
	CLZNumber	larger;
	CLZNumber	smaller;
	CLZNumber	one;
	CLZNumber	two;

	magnitude.Assign( pX );

	if ( magnitude.IsNegative() )
		magnitude.m_lSign = kPositive;

	one.SetBase( pX->m_lBase, false );
	one.Assign( 1L );
	two.SetBase( pX->m_lBase, false );
	two.Assign( 2L );

	long	lResult = lz_exp_e( &magnitude, &larger );

	lResult = lzi_first_failure( lResult, lz_divide( &one, &larger, &smaller ) );

	if ( pSine != nullptr )
	{
		pSine->SetPrecision( glWorkingPrecision );
		lResult = lzi_first_failure( lResult, lz_subtract( &larger, &smaller, pSine ) );
		lResult = lzi_first_failure( lResult, lz_divide( pSine, &two ) );

		if ( pX->IsNegative() )
			lz_negate( pSine );
	}

	if ( pCosine != nullptr )
	{
		pCosine->SetPrecision( glWorkingPrecision );
		lResult = lzi_first_failure( lResult, lz_add( &larger, &smaller, pCosine ) );
		lResult = lzi_first_failure( lResult, lz_divide( pCosine, &two ) );
	}

	glWorkingPrecision = l_saved;

	return lResult;
}

long lz_hypersine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_hypersine, kRetakeDoubt );

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	answer.DecideTies( lzi_outside_x_decide );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryStaysInfinite ) )
		return kNoError;

    long lResult = kNoError;

    CLZNumber result;

    // Wider than the answer, rounded once at the end -- see
    //	kQuotientGuardDigits...
    long	l_savedPrecision = glWorkingPrecision;

    glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, leftOperand->m_lBase );

    result.SetPrecision( glWorkingPrecision );

    /*
        sinh(x) = (e^x - e^-x) / 2, from one exponential and its
        reciprocal, a small argument widened by its noughts -- see
        lzi_sinh_cosh.

        Two things used to sit in front of that and both were wrong.
        The argument was first taken modulo 2*pi, which is a reduction
        for the circular functions and means nothing here -- these are
        not periodic, and sinh(10) came back as sinh(10 - 2*pi).
        Then the exponentials went through lz_raise_e, which raises e
        to an integer power by squaring and simply drops a fraction:
        sinh(0.5) answered sinh(0), which is nought, and sinh(1.1)
        answered sinh(1).

        Between them the three hyperbolic functions were right only
        for whole-numbered arguments below 2*pi, which is to say for
        0 through 6. The suite tested them at 1.
    */
    lResult = lzi_sinh_cosh( leftOperand, &result, nullptr );

    // The guard places never leave, and the answer is rounded to
    //	the rule for a quotient -- see lzi_round_quotient_answer...
    glWorkingPrecision = l_savedPrecision;

    lzi_round_quotient_answer( &result );

    // Record the result...
    if ( targetNumber != nullptr )
    {
        if ( !(targetNumber->IsConstant()) )
            targetNumber->Assign( &result );
    }
    else if ( !(leftOperand->IsConstant()) )
        leftOperand->Assign( &result );

    return lResult;
}


long lz_hypercosine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_hypercosine, kRetakeDoubt );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryStaysInfinite ) )
		return kNoError;

    long lResult = kNoError;

    CLZNumber result;

    // Wider than the answer, rounded once at the end -- see
    //	kQuotientGuardDigits...
    long	l_savedPrecision = glWorkingPrecision;

    glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, leftOperand->m_lBase );

    result.SetPrecision( glWorkingPrecision );

    // cosh(x) = (e^x + e^-x) / 2, from one exponential and its
    //	reciprocal -- see lzi_sinh_cosh...
    lResult = lzi_sinh_cosh( leftOperand, nullptr, &result );

    // The guard places never leave, and the answer is rounded to
    //	the rule for a quotient -- see lzi_round_quotient_answer...
    glWorkingPrecision = l_savedPrecision;

    lzi_round_quotient_answer( &result );

    // Record the result...
    if ( targetNumber != nullptr )
    {
        if ( !(targetNumber->IsConstant()) )
            targetNumber->Assign( &result );
    }
    else if ( !(leftOperand->IsConstant()) )
        leftOperand->Assign( &result );

    return lResult;
}


/*
	**sinh and cosh of one argument, each answered exactly as lz_hypersine
	and lz_hypercosine answer it** -- stage 3 of the complex layer,
	2026-10-22, on Rich's calls that day: this name, and both targets
	required and different, either of them the operand if wanted.

	**Both from one exponential and its reciprocal** -- see lzi_sinh_cosh
	-- under the guard the two functions work at, each rounded once by the
	rule for a quotient. An answer within its function's doubt window, a
	flagged one, nought's exact pair, or any where the work failed is that
	function's to answer instead, as the host's call -- see
	lzi_rounded_clear_of_doubt and lz_sine_cosine. So every answer is the
	function's to the character.
*/

long lz_hypersine_hypercosine( CLZNumber* leftOperand, CLZNumber* sineTarget, CLZNumber* cosineTarget )
{
	if ( (leftOperand == nullptr) || (sineTarget == nullptr) || (cosineTarget == nullptr) ||
		 (sineTarget == cosineTarget) )
		return lzi_pair_refused( sineTarget, cosineTarget );

	CLZNumber	x;
	CLZNumber	sine;
	CLZNumber	cosine;
	bool		b_sineOwn = true;
	bool		b_cosineOwn = true;

	x.Assign( leftOperand );

	if ( x.FlagIsOff( kfQNaN ) && x.FlagIsOff( kfInfinity ) && !lzi_exactly( &x, 0 ) )
	{
		// The work under a guard of its own, so every call inside nests...
		CLZNumber		scratch;
		CLZAnswerBase	work( &x, &scratch );

		work.Exact();

		long	l_saved = glWorkingPrecision;

		glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, x.m_lBase );

		sine.SetPrecision( glWorkingPrecision );
		cosine.SetPrecision( glWorkingPrecision );

		long	lResult = lzi_sinh_cosh( &x, &sine, &cosine );

		glWorkingPrecision = l_saved;

		if ( lResult == kNoError )
		{
			b_sineOwn = !lzi_rounded_clear_of_doubt( &sine, kRetakeDoubt, true );
			b_cosineOwn = !lzi_rounded_clear_of_doubt( &cosine, kRetakeDoubt, true );
		}
	}

	// ...and an answer that cannot stand its own function's, as the host's
	//	call...
	glPairHandoffs += (b_sineOwn ? 1 : 0) + (b_cosineOwn ? 1 : 0);

	long	l_sineResult = b_sineOwn ? lz_hypersine( &x, sineTarget )
									 : lzi_answer_as_function( &sine, &x, sineTarget, true );
	long	l_cosineResult = b_cosineOwn ? lz_hypercosine( &x, cosineTarget )
										 : lzi_answer_as_function( &cosine, &x, cosineTarget, true );

	return lzi_first_failure( l_sineResult, l_cosineResult );
}


long lz_hypertangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_hypertangent, kRetakeDoubt );

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	answer.DecideTies( lzi_inside_x_decide );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	// Plus or minus one past the width, without sinh or cosh -- see
	//	lzi_unit_beyond_width...
	if ( lzi_unit_beyond_width( leftOperand, targetNumber ) )
		return kNoError;

    long lResult = kNoError;

    CLZNumber part1;
    CLZNumber part2;
    CLZNumber result;

    // Wider than the answer, rounded once at the end -- see
    //	kQuotientGuardDigits. The division's target is widened by
    //	hand, since lz_divide takes its width from it...
    long	l_savedPrecision = glWorkingPrecision;

    glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, leftOperand->m_lBase );

    result.SetPrecision( glWorkingPrecision );

    // Calculate...
    // sinh over cosh, both from one exponential -- see lzi_sinh_cosh...
    part1.SetPrecision( glWorkingPrecision );
    part2.SetPrecision( glWorkingPrecision );
    lResult = lzi_first_failure( lResult, lzi_sinh_cosh( leftOperand, &part1, &part2 ) );
    lResult = lzi_first_failure( lResult, lz_divide( &part1, &part2, &result ) );

    // The guard places never leave, and the answer is rounded to
    //	the rule for a quotient -- see lzi_round_quotient_answer...
    glWorkingPrecision = l_savedPrecision;

    lzi_round_quotient_answer( &result );

    // Record the result...
    if ( targetNumber != nullptr )
    {
        if ( !(targetNumber->IsConstant()) )
            targetNumber->Assign( &result );
    }
    else if ( !(leftOperand->IsConstant()) )
        leftOperand->Assign( &result );

    return lResult;
}


/*
	The three hyperbolic reciprocals, kept with the family they belong
	to rather than filed alphabetically, the way the sine, cosine and
	tangent above them are.

	Each is a quotient, and each leaves its pole to lz_divide to report
	rather than testing for it first -- which is what lz_cotangent and
	lz_cosecant already do on the circular side, and it gives the
	fuller answer: kDivideByZero with kfInfinity and kfDiv0 both set,
	rather than a bare refusal. The hyperbolic secant has no pole at
	all, the hyperbolic cosine never falling below one.
*/

long lz_hypercosecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_hypercosecant, kRetakeDoubt );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// An infinity is held as zero here, so this has to come before any
	//	test of the value. csch is nought at either infinity -- a finite
	//	limit, which this engine's propagation cannot return, so a
	//	not-a-number rather than an invented answer...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber hyperbolicSine;
	CLZNumber result;

	// Wider than the answer, rounded once at the end -- see
	//	kQuotientGuardDigits. The division's target is widened by
	//	hand, since lz_divide takes its width from it...
	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, leftOperand->m_lBase );

	// ...and by the answer's size below one, read off the argument -- see
	//	lzi_small_divisor_widening. Twice that until 2026-10-15, once for
	//	the digits e^x - e^-x lost and once for those the answer adds in
	//	front; sinh keeps its own digits since, see lz_hypersine...
	glWorkingPrecision += lzi_small_divisor_widening( leftOperand ) / 2;

	result.SetPrecision( glWorkingPrecision );

	// csch(x) = 1 / sinh(x)...
	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	lResult = lzi_first_failure( lResult, lz_hypersine( leftOperand, &hyperbolicSine ) );
	lResult = lzi_first_failure( lResult, lz_divide( &one, &hyperbolicSine, &result ) );

	// The guard places never leave, and the answer is rounded to
	//	the rule for a quotient -- see lzi_round_quotient_answer...
	glWorkingPrecision = l_savedPrecision;

	lzi_round_quotient_answer( &result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_hypercotangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_hypercotangent, kRetakeDoubt );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// coth is plus one at a positive infinity and minus one at a
	//	negative, and an infinity here carries no sign to tell them
	//	apart -- the same reading the hyperbolic tangent takes...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	// Plus or minus one past the width, without sinh or cosh -- see
	//	lzi_unit_beyond_width...
	if ( lzi_unit_beyond_width( leftOperand, targetNumber ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber hyperbolicSine;
	CLZNumber hyperbolicCosine;
	CLZNumber result;

	// Wider than the answer, rounded once at the end -- see
	//	kQuotientGuardDigits. The division's target is widened by
	//	hand, since lz_divide takes its width from it...
	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, leftOperand->m_lBase );

	// ...and by the answer's size below one, read off the argument -- see
	//	lzi_small_divisor_widening. Twice that until 2026-10-15, once for
	//	the digits e^x - e^-x lost and once for those the answer adds in
	//	front; sinh keeps its own digits since, see lz_hypersine...
	glWorkingPrecision += lzi_small_divisor_widening( leftOperand ) / 2;

	result.SetPrecision( glWorkingPrecision );

	/*
		coth(x) = cosh(x) / sinh(x).

		Taken from the two directly rather than as one over the
		hyperbolic tangent, which would divide twice and round twice
		for the same answer -- and which is how lz_cotangent takes the
		circular one.
	*/

	// Both from one exponential -- see lzi_sinh_cosh...
	hyperbolicSine.SetPrecision( glWorkingPrecision );
	hyperbolicCosine.SetPrecision( glWorkingPrecision );
	lResult = lzi_first_failure( lResult, lzi_sinh_cosh( leftOperand, &hyperbolicSine, &hyperbolicCosine ) );
	lResult = lzi_first_failure( lResult,
					lz_divide( &hyperbolicCosine, &hyperbolicSine, &result ) );

	// The guard places never leave, and the answer is rounded to
	//	the rule for a quotient -- see lzi_round_quotient_answer...
	glWorkingPrecision = l_savedPrecision;

	lzi_round_quotient_answer( &result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_hypersecant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_hypersecant, kRetakeDoubt );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// sech is nought at either infinity. A finite limit, which this
	//	engine's propagation cannot return, so a not-a-number...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long lResult = kNoError;

	CLZNumber one;
	CLZNumber hyperbolicCosine;
	CLZNumber result;

	// Wider than the answer, rounded once at the end -- see
	//	kQuotientGuardDigits. The division's target is widened by
	//	hand, since lz_divide takes its width from it...
	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kQuotientGuardDigits, leftOperand->m_lBase );

	result.SetPrecision( glWorkingPrecision );

	// sech(x) = 1 / cosh(x), and no pole: cosh never falls below one...
	// In the operand's base -- see CLZAnswerBase...
	one.SetBase( leftOperand->m_lBase );
	one.Assign( 1L );

	lResult = lzi_first_failure( lResult, lz_hypercosine( leftOperand, &hyperbolicCosine ) );
	lResult = lzi_first_failure( lResult, lz_divide( &one, &hyperbolicCosine, &result ) );

	// The guard places never leave, and the answer is rounded to
	//	the rule for a quotient -- see lzi_round_quotient_answer...
	glWorkingPrecision = l_savedPrecision;

	lzi_round_quotient_answer( &result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_integer( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( leftOperand );
			targetNumber->StripFraction();
		}
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->StripFraction();

	return kNoError;
}


static long lzi_log_b_work( CLZNumber* leftOperand, CLZNumber* /* rightOperand */, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	long lResult = kNoError;

	auto logLeftBase_e = new CLZNumber;
	auto result = new CLZNumber;

	lResult = lzi_first_failure( lResult, lz_log_e( leftOperand, logLeftBase_e ) );

	// ln(base) as wide as the work, where it was gLogBaseE's 120 places
	//	at any width -- see lzi_exp_b_work and lzi_log_base_in_base...
	lResult = lzi_first_failure( lResult, lz_divide( logLeftBase_e, lzi_log_base_in_base( leftOperand->m_lBase, glWorkingPrecision ), result ) );

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete logLeftBase_e;
    delete result;

	return lResult;
}

/*
	**The logarithm in the operand's own base: ln x over ln b, worked wider
	and rounded once**, as lz_log_n beside it always was. Until 2026-09-26
	it took ln x rounded to the working width and divided that by the
	table's ln b, and the division rounded again: 600 of 3,500 random
	arguments missed correct rounding at the default, in every base and
	at every size -- 10 to 18% in base ten, 40% in base two, by up to
	1.36 units. The audit's last six near misses were all this.
*/

long lz_log_b( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_log_b, kRetakeDoubt );

	// ...and an answer that is a fraction decided exactly in doubt -- see
	//	lzi_log_decide...
	answer.DecideTies( lzi_log_b_decide );

	// Exact at one and nowhere else, asked before anything is done to the
	//	operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 1 ) );

	// ln x / ln b, sized by its own logarithm, as lz_log_n's is...
	double	d_lnX = 0.0;
	long	l_magnitude = -1;

	if ( lzi_ln_estimate( leftOperand, &d_lnX ) && (d_lnX != 0.0) )
		l_magnitude = lzi_places_for( std::log( std::fabs( d_lnX / std::log( static_cast<double>( leftOperand->m_lBase ) ) ) ),
									  leftOperand->m_lBase );

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	long	lResult = lzi_widened_by_answer( lzi_log_b_work, leftOperand, nullptr, targetNumber, l_magnitude );

	/*
		**The spelling it always had: every place written.** The answer
		used to come straight from a division at the working width, which
		writes its quotient out to the full places. Worked wider, a wide
		quotient that lands exactly on a whole number comes back bare, so
		log10(100) became "2" while log10(1000), rounded up from
		2.999..., kept its noughts. Zero, and a flagged answer held as
		zero, stay as they were.
	*/

	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( !(p_answer->IsConstant()) )
	{
		p_answer->CheckForZero();

		if ( !(p_answer->IsZero()) && (p_answer->m_lRadix < glWorkingPrecision) )
			p_answer->PadRight( glWorkingPrecision - p_answer->m_lRadix );
	}

	return lResult;
}


/*
	**Past a few hundred digits the logarithm inverts the exponential by
	Newton**, rather than summing its series. The series gains about
	0.67 digits a term at the worst of its band and every term is a
	multiply at the full width, so its cost grows as the cube of the
	width: ln(1e40) at 9,000 places took 7.4 s, some 13,500 full-width
	multiplies. The exponential has been right at any width since
	lzi_e_in_base, and

		y  <-  y + 2 (x - e^y) / (x + e^y)

	converges on ln(x) cubically. So it starts from the series at a
	narrow width and doubles the width each step, and the whole costs
	about two exponentials at the full width.

	**It takes over at 150 digits, where the series would be long** --
	more than 128 terms, sized from its y; see lz_log_e. Until 2026-10-15
	it took over at 200, half again past where it first won, about 130,
	measured once the exponential halved from a hundred digits; the
	binary-splitting exponential since moved Newton's own cost, and the
	series sized from y moved the series'.

	**An argument under one is inverted first**, and the answer negated,
	as the series path does for a small one: e^y of a large negative y
	is tiny, and a fixed number of places would leave it few significant
	digits.

	**The guard is eight digits, converted to places** -- see
	lzi_guard_places. Until 2026-10-07 it was added as eight places, so
	in base two it was eight bits, two and a half digits, and from 680
	binary places, where this takes over, a base-two logarithm came
	out as far as 1.5e-3 of a unit of its last place against 5e-8 just
	below the switch, and 2 of the 10 base-two logarithms the near-half
	census built at 1000 places rounded wrong. Bases ten to sixteen kept
	the guard they had.
*/

const long kLogNewtonFrom	= 150;		// decimal digits of working width
const long kLogNewtonTerms	= 128;		// and the series longer than this -- see lz_log_e
const long kLogNewtonSeed	= 60;		// decimal digits the series starts it at
const long kLogNewtonGuard	= 8;		// decimal digits

/*
	The width in force, which only the tests move -- to force the series
	past it and compare the two, which are independent. Never below the
	seed, and lz_log_e refuses Newton at the seed's own width as well --
	see there. A width set here forces Newton wherever it reaches, however
	short the series would be.
*/

static long gLogNewtonFrom = kLogNewtonFrom;

long lzi_set_log_newton_from( long lDigits )
{
	long	l_previous = gLogNewtonFrom;

	gLogNewtonFrom = (lDigits > kLogNewtonSeed) ? lDigits : (kLogNewtonSeed + 1);

	return l_previous;
}

// For the tests: the most places past its stage of an iterate or operand
//	in any stage but the last since the last read, this thread's, which the
//	read sets back -- nought or less, and LONG_MIN for none...
static thread_local long gLogNewtonOverrun = LONG_MIN;

long lzi_log_newton_overrun()
{
	long	l_overrun = gLogNewtonOverrun;

	gLogNewtonOverrun = LONG_MIN;

	return l_overrun;
}

static long lzi_log_newton( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	long	lResult = kNoError;
	long	l_saved = glWorkingPrecision;
	long	l_base = leftOperand->m_lBase;
	long	l_target = l_saved + lzi_guard_places( kLogNewtonGuard, l_base );

	CLZNumber	x;
	CLZNumber	one;
	bool		b_negate = false;

	x.Assign( leftOperand );
	one.SetBase( l_base, false );
	one.Assign( 1L );

	if ( lz_compare( &x, &one ) == kLessThan )
	{
		CLZNumber	inverted;

		// lz_divide takes its width from its target...
		glWorkingPrecision = l_target;
		inverted.SetPrecision( l_target );
		lResult = lzi_first_failure( lResult, lz_divide( &one, &x, &inverted ) );

		x.Assign( &inverted );
		b_negate = true;
	}

	CLZNumber	y;
	long		l_width = CLZNumber::EquivalentPlaces( kLogNewtonSeed, 10, l_base );

	if ( l_width > l_target )
		l_width = l_target;

	glWorkingPrecision = l_width;
	lResult = lzi_first_failure( lResult, lz_log_e( &x, &y ) );

	CLZNumber	e_y;
	CLZNumber	difference;
	CLZNumber	sum;
	CLZNumber	step;
	CLZNumber	narrowed;

	/*
		**Each stage but the last at its own width, operand and iterate
		both** -- since 2026-10-17, the Newton audit. The operand went into
		every stage at its full length, and the step is a small quotient,
		which answers to the width's worth of significant digits and its
		noughts as well, so the iterate came out half again as wide as its
		stage and the next exponential was handed all of it. x is narrowed
		to the stage and two places, as the cube root narrows its operand,
		and the iterate rounded to the stage; the last stage takes x whole
		and is rounded below, so the answer is worked from all of it.
		Over all fifteen bases an ordinary operand is up to 6% faster and
		one three widths long 8 to 11% at 3,000 digits; no answer moved in
		6,360 censused.
	*/

	while ( l_width < l_target )
	{
		l_width = ((l_width * 2) < l_target) ? (l_width * 2) : l_target;
		glWorkingPrecision = l_width;

		bool		b_last = (l_width >= l_target);
		CLZNumber*	p_x = &x;

		if ( !b_last && (x.m_lRadix > (l_width + 2)) )
		{
			narrowed.Assign( &x );
			narrowed.Round( kDefaultRoundingMode, l_width + 2 );
			p_x = &narrowed;
		}

		lResult = lzi_first_failure( lResult, lz_exp_e( &y, &e_y ) );

		lz_subtract( p_x, &e_y, &difference );
		lz_add( p_x, &e_y, &sum );

		step.SetPrecision( l_width );
		lResult = lzi_first_failure( lResult, lz_divide( &difference, &sum, &step ) );

		lz_add( &step, &step );
		lz_add( &y, &step );

		if ( !b_last )
		{
			y.SetPrecision( l_width );
			y.Round();

			// For the tests: the most past its stage of the iterate or the
			//	operand, read by lzi_log_newton_overrun...
			long	l_over = (y.m_lRadix - l_width > p_x->m_lRadix - (l_width + 2)) ? (y.m_lRadix - l_width)
																					: (p_x->m_lRadix - (l_width + 2));

			if ( l_over > gLogNewtonOverrun )
				gLogNewtonOverrun = l_over;
		}
	}

	glWorkingPrecision = l_saved;

	if ( b_negate )
		lz_negate( &y );

	y.SetPrecision( l_saved );
	y.Round();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &y );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &y );

	return lResult;
}


/*
	**The logarithm by the arithmetic-geometric mean** -- the second half
	of part B of the slowness item, 2026-10-14. lz_log_e takes it by the
	rule at lzi_log_agm_floor. For s large,

		ln s  =  pi / (2 AGM(1, 4/s))

	to within 4 ln(s) / s^2, so with s = x b^m past b^(W/2) the logarithm
	of x is that less m ln b. The AGM settles in about twice log2 of the
	width cycles, each a multiplication and a square root, and pi and
	ln b are the constants' cache's, computed once a width.

	**It needs a square root of its own.** lz_root_2 decides its rounding
	exactly and costs 15 to 29 multiplications at three to thirty thousand
	places; on it the AGM took 2 to 5 times the logarithm by Newton on the
	exponential. lzi_root_for_agm answers a root good to the width, which
	is all the AGM wants, by Newton on the reciprocal square root --
	multiplications only -- with two things that made the difference,
	measured on a prototype in base ten at 10,000 places:

	- **the argument rounded to each step's width.** The engine multiplies
	  every digit it is handed and rounds afterwards, so the early, narrow
	  steps were each a full-width product: 0.098 s to 0.071.
	- **Karp and Markstein's last step**, the root itself at half the width
	  corrected once, s + r (m - s^2) / 2, in place of a last reciprocal
	  step at the whole width and a multiplication: 0.071 s to 0.051.

	And it is seeded from the cycle before once the two means are close:
	the next product moves by about their difference, so the reciprocal
	root already held is good to that many places.

	**lz_root_2 takes it too**, past its floors, in place of Heron's
	iteration -- see lzi_root_core_floor.
*/

const long	kLogAgmGuardDigits	= 8;		// decimal digits past the answer
const long	kLogAgmMaxCycles	= 300;		// against a failure only; widths to 100,000 places take 44
const long	kLogAgmSlackDigits	= 2;		// decimal digits each working width carries past its need

/*
	**Every working width's margin is digits, converted to places** -- see
	lzi_guard_places. Written as two places, it was two digits in base ten
	and two bits in base two, and base two's error grew with the width and
	the forty-odd cycles' roundings: 11 units of the last place at 3,000
	places, 72 at 10,000 and 153 at 30,000, where base ten stayed under a
	unit.
*/

static long lzi_agm_slack( long lBase )
{
	return lzi_guard_places( kLogAgmSlackDigits, lBase );
}

// The reciprocal square root held between an AGM's cycles: r is about
//	1 / sqrt( p / b^k ), good to good places...
struct SLZRootSeed
{
	CLZNumber	r;
	long		k = LONG_MIN;
	long		good = 0;
};

// Half of pValue in place, at the working width: exact by the reciprocal
//	in an even base, by a division in an odd one, where a half never
//	ends -- the AGM's own rule, see lz_agm...
static void lzi_halve_at_width( CLZNumber* pValue )
{
	CLZNumber	two;

	two.SetBase( pValue->m_lBase, false );
	two.Assign( 2L );

	if ( (pValue->m_lBase % 2) == 0 )
	{
		lz_reciprocal( &two );
		lz_multiply( pValue, &two );
	}
	else
	{
		pValue->SetPrecision( glWorkingPrecision );
		lz_divide( pValue, &two );
	}
}

/*
	sqrt( pValue ), positive, into pRoot to lRelPlaces places relative.
	pSeed carries the reciprocal root between calls; lSeedGood is how many
	places the caller knows the argument has moved by no more than, nought
	for none.
*/

// For the tests: the most places past its working width of any root
//	answered since the last read, this thread's, which the read sets back --
//	nought or less since each is rounded, and LONG_MIN for none...
static thread_local long gRootCoreOverrun = LONG_MIN;

long lzi_root_core_overrun()
{
	long	l_overrun = gRootCoreOverrun;

	gRootCoreOverrun = LONG_MIN;

	return l_overrun;
}

static void lzi_root_for_agm( CLZNumber* pValue, CLZNumber* pRoot, long lRelPlaces, SLZRootSeed* pSeed, long lSeedGood )
{
	long	l_base = pValue->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_exponent = lzi_leading_exponent( pValue );
	long	l_even = l_exponent - (((l_exponent % 2) + 2) % 2);
	long	l_seedPlaces = CLZNumber::EquivalentPlaces( 14, 10, l_base );

	// m = value / b^k, k even, in [1, b^2): its root is the value's moved
	//	back by k/2, exactly...
	CLZNumber	m;

	m.Assign( pValue );
	lzi_scale_by_base_power( &m, l_even );

	long	l_good;

	if ( (pSeed->k != l_even) || (lSeedGood < l_seedPlaces) )
	{
		// From a double, by its text: Assign( double ) writes "%f", six
		//	places after the point, and a six-place seed loses a step...
		double	d_value = 1.0;
		char	sz_seed[ 40 ];

		lzi_value_estimate( &m, &d_value );
		::snprintf( sz_seed, sizeof( sz_seed ), "%.17f", 1.0 / std::sqrt( d_value ) );

		glWorkingPrecision = l_seedPlaces + 4;
		pSeed->r.SetBase( 10, false );
		pSeed->r.Assign( sz_seed );
		pSeed->r.SetBase( l_base, true );
		l_good = l_seedPlaces;
	}
	else
		l_good = (lSeedGood < pSeed->good) ? lSeedGood : pSeed->good;

	pSeed->k = l_even;

	// 1/sqrt(m) to half the width: r <- r + r (1 - m r^2) / 2, each step
	//	at its own width with m rounded to it...
	long		l_slack = lzi_agm_slack( l_base );
	long		l_half = ((lRelPlaces + 1) / 2) + l_slack;
	CLZNumber	m_step, square, residual, one;

	one.SetBase( l_base, false );
	one.Assign( 1L );

	while ( l_good < l_half )
	{
		l_good = (((2 * l_good) - l_slack) < l_half) ? ((2 * l_good) - l_slack) : l_half;
		glWorkingPrecision = l_good + l_slack;

		m_step.Assign( &m );
		m_step.SetPrecision( l_good + (2 * l_slack) );
		m_step.Round( kDefaultRoundingMode );
		pSeed->r.SetPrecision( l_good + (2 * l_slack) );
		pSeed->r.Round( kDefaultRoundingMode );

		lz_multiply( &pSeed->r, &pSeed->r, &square );
		lz_multiply( &square, &m_step );
		lz_subtract( &one, &square, &residual );
		lz_multiply( &residual, &pSeed->r );
		lzi_halve_at_width( &residual );

		pSeed->r.SetPrecision( glWorkingPrecision );
		lz_add( &pSeed->r, &residual );
	}

	pSeed->good = l_half;

	// Karp and Markstein: s = m r at half the width, then
	//	s + r (m - s^2) / 2 at the whole...
	CLZNumber	r_half;

	glWorkingPrecision = l_half + l_slack;
	m_step.Assign( &m );
	m_step.SetPrecision( l_half + (2 * l_slack) );
	m_step.Round( kDefaultRoundingMode );
	r_half.Assign( &pSeed->r );
	r_half.SetPrecision( l_half + (2 * l_slack) );
	r_half.Round( kDefaultRoundingMode );
	lz_multiply( &m_step, &r_half, pRoot );
	pRoot->SetPrecision( l_half + l_slack );
	pRoot->Round( kDefaultRoundingMode );

	glWorkingPrecision = lRelPlaces + l_slack;
	lz_multiply( pRoot, pRoot, &square );
	lz_subtract( &m, &square, &residual );
	lz_multiply( &residual, &r_half );
	lzi_halve_at_width( &residual );

	pRoot->SetPrecision( glWorkingPrecision );
	lz_add( pRoot, &residual );

	/*
		**Rounded to the working width** -- since 2026-10-17. In an odd base
		the halving above is a division, a half never ending there, and a
		division answers wider than it was asked; the sum kept every place,
		so the root came back as much as half again as long as asked --
		31,472 places for 20,970 in base three, 925 of 2,017 odd-base roots
		in a scan of 4,680. lz_root_2 rounds its own and paid nothing, but
		the logarithm by the AGM multiplies this root into its next cycle,
		and a product works every digit it is handed: lz_log_e in odd bases
		3 to 10% faster rounded, base nine at 10,480 places 47.1 ms to 43.0,
		even bases level. Nothing past the width was ever relied on, and no
		answer moves...
	*/
	pRoot->Round( kDefaultRoundingMode );

	if ( (pRoot->m_lRadix - glWorkingPrecision) > gRootCoreOverrun )
		gRootCoreOverrun = pRoot->m_lRadix - glWorkingPrecision;

	lzi_scale_by_base_power( pRoot, -(l_even / 2) );

	glWorkingPrecision = l_saved;
}

/*
	AGM( 1, z ), z small and positive, to lRelPlaces places relative. Each
	cycle works wide enough to hold the geometric mean that many places,
	as lz_agm does, and seeds the root from the cycle before by the places
	the two means agree to. Answers the cycles it took.
*/

static thread_local long gAgmLogCycles = 0;

static long lzi_agm_for_log( CLZNumber* pSmall, CLZNumber* pMean, long lRelPlaces )
{
	long		l_base = pSmall->m_lBase;
	long		l_saved = glWorkingPrecision;
	long		l_slack = lzi_agm_slack( l_base );
	long		l_seedGood = 0;
	long		l_cycles = 0;
	SLZRootSeed	seed;
	CLZNumber	a_value, g_value, a_next, g_next, product, difference;

	a_value.SetBase( l_base, false );
	a_value.Assign( 1L );
	g_value.Assign( pSmall );

	for ( ; l_cycles < kLogAgmMaxCycles; l_cycles++ )
	{
		long	l_gExponent = lzi_leading_exponent( &g_value );

		glWorkingPrecision = lRelPlaces + l_slack + ((l_gExponent < 0) ? -l_gExponent : 0);

		lz_add( &a_value, &g_value, &a_next );
		lzi_halve_at_width( &a_next );
		lz_multiply( &a_value, &g_value, &product );
		lzi_root_for_agm( &product, &g_next, lRelPlaces + l_slack, &seed, l_seedGood );

		glWorkingPrecision = lRelPlaces + l_slack + ((l_gExponent < 0) ? -l_gExponent : 0);
		lz_subtract( &a_next, &g_next, &difference );
		difference.CheckForZero();

		long	l_agree = difference.IsZero() ? (lRelPlaces + l_slack)
											  : (lzi_leading_exponent( &a_next ) - lzi_leading_exponent( &difference ));

		a_value.Assign( &a_next );
		g_value.Assign( &g_next );

		if ( l_agree > (lRelPlaces + 1) )
			break;

		// The next product moves by about the two means' difference...
		l_seedGood = l_agree - l_slack;
	}

	glWorkingPrecision = lRelPlaces + l_slack;
	lz_add( &a_value, &g_value, pMean );
	lzi_halve_at_width( pMean );

	glWorkingPrecision = l_saved;
	gAgmLogCycles += l_cycles + 1;

	return l_cycles + 1;
}

// The cycles this thread's AGM logarithms have taken, for the tests and
//	for measuring...
long lzi_agm_log_cycles()
{
	return gAgmLogCycles;
}

// The places of a whole number's magnitude in lBase, at least one...
static long lzi_whole_places( long lValue, long lBase )
{
	long	l_places = 1;

	for ( long l_rest = (lValue < 0) ? -lValue : lValue; l_rest >= lBase; l_rest /= lBase )
		l_places++;

	return l_places;
}

/*
	ln( b^n ) = pi / (2 AGM(1, 4 / b^n)), to lPlaces places, for n large
	enough that the formula's 4 n ln(b) / b^(2n) is past them. No
	logarithm is asked of anything, so ln(base) can be built from it
	without the logarithm asking for ln(base).
*/

static void lzi_log_power_by_agm( long lBase, long lPower, CLZNumber* pLog, long lPlaces )
{
	long	l_saved = glWorkingPrecision;
	double	d_logS = static_cast<double>( lPower ) * std::log( static_cast<double>( lBase ) );
	long	l_slack = lzi_agm_slack( lBase );
	long	l_extra = lzi_whole_places( static_cast<long>( d_logS ) + 1, lBase ) + l_slack;
	long	l_rel = lPlaces + l_extra;

	CLZNumber	small, mean, half_pi( kPiOverTwo );

	small.SetBase( lBase, false );
	small.Assign( 4L );
	lzi_scale_by_base_power( &small, lPower );

	lzi_agm_for_log( &small, &mean, l_rel );
	lzi_pi_family_in_base( &half_pi, kPiOverTwo, lBase, l_rel );

	// (pi/2) / AGM, never pi / (2 AGM): the AGM is about 1 / ln s, and a
	//	product rounds to places, so 2 AGM lost as many significant ones as
	//	the AGM has noughts after its point -- three at a thousand places,
	//	times ln s, a thousand units...
	glWorkingPrecision = lPlaces + l_slack;
	pLog->SetBase( lBase, false );
	pLog->SetPrecision( lPlaces + l_slack );
	lz_divide( &half_pi, &mean, pLog );

	glWorkingPrecision = l_saved;
}

// The power of the base an AGM logarithm at lPlaces places takes s to:
//	s^2 past 4 ln(s) b^lPlaces, with a place to spare...
static long lzi_agm_log_reach( long lPlaces, long lBase )
{
	double	d_logS = (static_cast<double>( lPlaces ) / 2.0 + 2.0) * std::log( static_cast<double>( lBase ) );

	return (lPlaces / 2) + 1 + lzi_whole_places( static_cast<long>( 4.0 * d_logS ) + 1, lBase ) + 1;
}


/*
	ln of lBase in lBase, good to lPlaces: gLogBaseE within its 120
	places, converted and kept as every table entry is -- see
	lzi_table_in_base -- and past them the logarithm of lBase itself,
	worked at the width and kept, as e and pi are past their text.
	lBase written in its own base is "10", so that logarithm never takes
	lzi_log_scaled, which asks for this. The value is the cache's own:
	read it, or copy it before working on it.
*/

static const char	kComputedLogBase[] = "ln(base), computed";

static CLZNumber* lzi_log_base_in_base( long lBase, long lPlaces, bool bByAgm )
{
	// Each text's places, counted once: lz_log_b and lz_exp_b ask on
	//	every call...
	static const std::array<long, 17>	sl_textPlaces = []
	{
		std::array<long, 17>	places = {};

		for ( long l_base = 2; l_base <= 16; l_base++ )
			places[ l_base ] = lzi_text_places( kLogBaseEText[ l_base ] );

		return places;
	}();

	if ( lPlaces <= CLZNumber::EquivalentPlaces( sl_textPlaces[ lBase ], 10, lBase ) )
		return lzi_table_in_base( gLogBaseE[ lBase ], lBase );

	auto	key = std::make_tuple( static_cast<const void*>( kComputedLogBase ), lBase, lPlaces, lPlaces );
	auto	found = gTabledInBase.find( key );

	// A wider one is rounded, not computed again -- see lzi_cached_wider...
	if ( found == gTabledInBase.end() )
	{
		if ( CLZNumber* p_wider = lzi_cached_wider( kComputedLogBase, lBase, lPlaces, -1 ) )
			return lzi_cached_rounded( key, p_wider, lPlaces );
	}

	if ( found == gTabledInBase.end() )
	{
		long		l_saved = glWorkingPrecision;
		auto		p_log = new CLZNumber;
		CLZNumber	base;

		glWorkingPrecision = lPlaces;
		p_log->SetBase( lBase, false );

		if ( bByAgm )
		{
			// ln(b^n) / n, n as large as an AGM logarithm at this width
			//	takes s, and the places n's magnitude costs on top...
			long		l_power = lzi_agm_log_reach( lPlaces, lBase );
			long		l_places = lPlaces + lzi_whole_places( l_power, lBase ) + lzi_agm_slack( lBase );
			CLZNumber	log_power, count;

			lzi_log_power_by_agm( lBase, l_power, &log_power, l_places );

			glWorkingPrecision = l_places;
			count.SetBase( lBase, false );
			count.Assign( l_power );
			p_log->SetPrecision( l_places );
			lz_divide( &log_power, &count, p_log );

			glWorkingPrecision = lPlaces;
			p_log->SetPrecision( lPlaces );
			p_log->Round( kDefaultRoundingMode );
		}
		else
		{
			base.SetBase( lBase, false );
			base.Assign( lBase );
			lz_log_e( &base, p_log );
		}

		glWorkingPrecision = l_saved;

		found = gTabledInBase.emplace( key, p_log ).first;
	}

	return found->second;
}


/*
	**ln x by the AGM** -- see above lzi_root_for_agm. x positive; the
	working precision is the answer's places. An exact power of the base,
	b^j, is j ln b; anything else is pi / (2 AGM(1, 4/s)) - m ln b with
	s = x b^m past b^(W/2), worked kLogAgmGuardDigits past the answer.
*/

/*
	**The guard in force, which only the tests lower.** At eight digits it
	covers every margin above it whether that margin is there or not: the
	places for ln s, for m on ln b, the reach's log term and the slack each
	survive their removal at eight, and at none fail 83 to 1,097 of 1,500
	arguments, by up to 37,000 units -- an error the guard then hides only
	down to some 10^-4 of a unit, far outside lz_log_e's retake window. So
	a test needs a way to take the guard away and see each margin work.
*/

static long gLogAgmGuard = -1;

long lzi_set_log_agm_guard( long lDigits )
{
	long	l_previous = gLogAgmGuard;

	gLogAgmGuard = lDigits;

	return l_previous;
}

long lzi_log_by_agm( CLZNumber* pValue, CLZNumber* pTarget )
{
	long	l_base = pValue->m_lBase;
	long	l_saved = glWorkingPrecision;
	long	l_guard = (gLogAgmGuard >= 0) ? gLogAgmGuard : kLogAgmGuardDigits;
	long	l_width = l_saved + ((l_guard > 0) ? lzi_guard_places( l_guard, l_base ) : 0);
	long	l_exponent = lzi_leading_exponent( pValue );

	// x = x' b^e, x' in [1, b)...
	CLZNumber	mantissa, one, answer;

	mantissa.Assign( pValue );
	lzi_scale_by_base_power( &mantissa, l_exponent );
	one.SetBase( l_base, false );
	one.Assign( 1L );
	answer.SetBase( l_base, false );

	if ( lz_compare( &mantissa, &one ) == kEqualTo )
	{
		// b^j: j ln b, ln b by AGM when it has to be computed...
		CLZNumber	count;

		glWorkingPrecision = l_width;
		answer.Assign( lzi_log_base_in_base( l_base, l_width + lzi_whole_places( l_exponent, l_base ), true ) );
		count.SetBase( l_base, false );
		count.Assign( l_exponent );
		lz_multiply( &answer, &count );
	}
	else
	{
		// s = x' b^n, n the reach at this width, so m = n - e...
		long	l_reach = lzi_agm_log_reach( l_width, l_base );
		long	l_power = l_reach - l_exponent;
		double	d_logS = static_cast<double>( l_reach ) * std::log( static_cast<double>( l_base ) );
		long	l_slack = lzi_agm_slack( l_base );
		long	l_extra = lzi_whole_places( static_cast<long>( d_logS ) + 1, l_base ) + l_slack;
		long	l_rel = l_width + l_extra;

		// 4/s = (4 / x') b^-n, the division at the AGM's width...
		CLZNumber	small, four, mean, half_pi( kPiOverTwo ), multiple, count;

		four.SetBase( l_base, false );
		four.Assign( 4L );
		glWorkingPrecision = l_rel + l_slack;
		small.SetBase( l_base, false );
		small.SetPrecision( l_rel + l_slack );
		lz_divide( &four, &mantissa, &small );
		lzi_scale_by_base_power( &small, l_reach );

		lzi_agm_for_log( &small, &mean, l_rel );
		lzi_pi_family_in_base( &half_pi, kPiOverTwo, l_base, l_rel );

		// ln s = (pi/2) / AGM -- see lzi_log_power_by_agm for why not
		//	pi / (2 AGM)...
		glWorkingPrecision = l_width + l_slack;
		answer.SetPrecision( l_width + l_slack );
		lz_divide( &half_pi, &mean, &answer );

		// ...less m ln b, ln b carrying m's places in front of its error...
		multiple.Assign( lzi_log_base_in_base( l_base, l_width + lzi_whole_places( l_power, l_base ) + l_slack, true ) );
		count.SetBase( l_base, false );
		count.Assign( l_power );
		lz_multiply( &multiple, &count );
		lz_subtract( &answer, &multiple );
	}

	glWorkingPrecision = l_saved;
	answer.SetPrecision( l_saved );
	answer.Round();

	if ( !(pTarget->IsConstant()) )
		pTarget->Assign( &answer );

	return kNoError;
}


/*
	**A large or a small argument is brought to [1, base) by moving its
	radix, which is exact**, and its logarithm is ln m + k ln(base).

	The series path reduces by dividing by e, once for every unit of
	the answer, so ln x cost ln x divisions of a number as long as x:
	a thousand-digit argument took 11 ms and three thousand 80. That is
	quadratic, and lzi_log_gamma's rising product is 9,400 digits
	at 150 places. A small argument went the same way after
	being inverted.

	**k ln(base) carries k's own places in front of ln(base)'s error**,
	so ln(base) is taken that many places wider. The two logarithms are
	added kLogScaleGuardDigits past the width and rounded once.

	**From eight digits of exponent, where it starts to pay.** Moving
	the radix costs a fixed 10 to 20 microseconds at the default, the
	division path about 1.5 a digit, and the two cross at 8 digits in
	base ten, 6 hex places and some 30 binary ones -- about eight
	digits in each. At 158 places they cross at four. Below it the
	division path is left as it was.
*/

const long kLogScaleFrom		= 8;	// decimal digits of exponent, either side of one
const long kLogScaleGuardDigits	= 4;

// The threshold in force, which only the tests move. Never below three:
//	ln(base) is the logarithm of "10", a place of any base and never more
//	than two digits, which must not take this path...
static long gLogScaleFrom = kLogScaleFrom;

long lzi_set_log_scale_from( long lDigits )
{
	long	l_previous = gLogScaleFrom;

	gLogScaleFrom = (lDigits > 3) ? lDigits : 3;

	return l_previous;
}

static long lzi_log_scaled( CLZNumber* leftOperand, CLZNumber* targetNumber, long lExponent )
{
	long	lResult = kNoError;
	long	l_saved = glWorkingPrecision;
	long	l_base = leftOperand->m_lBase;
	long	l_exponentPlaces = 1;

	for ( long l_rest = (lExponent < 0) ? -lExponent : lExponent; l_rest >= l_base; l_rest /= l_base )
		l_exponentPlaces++;

	glWorkingPrecision = l_saved + lzi_guard_places( kLogScaleGuardDigits, l_base );

	// m = x / base^k, in [1, base)...
	CLZNumber	mantissa;

	mantissa.Assign( leftOperand );
	lzi_scale_by_base_power( &mantissa, lExponent );
	mantissa.SetPrecision( glWorkingPrecision );
	mantissa.Round();

	CLZNumber	logMantissa;

	logMantissa.SetBase( l_base, false );
	lResult = lzi_first_failure( lResult, lz_log_e( &mantissa, &logMantissa ) );

	// ...plus k ln(base)...
	CLZNumber	logBase;
	CLZNumber	count;

	logBase.Assign( lzi_log_base_in_base( l_base, glWorkingPrecision + l_exponentPlaces ) );

	count.SetBase( l_base, false );
	count.Assign( lExponent );

	lResult = lzi_first_failure( lResult, lz_multiply( &logBase, &count ) );
	lResult = lzi_first_failure( lResult, lz_add( &logMantissa, &logBase ) );

	glWorkingPrecision = l_saved;

	logMantissa.SetPrecision( glWorkingPrecision );
	logMantissa.Round();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &logMantissa );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &logMantissa );

	return lResult;
}


/*
	**Where lz_log_e takes the AGM**, measured 2026-10-14. Rich's rule,
	three conditions:

	- **Past the base's floor**, from lz_log_e timed with the AGM forced on
	  and off, warm, on ordinary arguments: past it the AGM never lost by
	  more than 0.2%, and gained up to 2.7 times. An odd base and the bases
	  from eleven take it from 50 places. **Base two never does**: warm it
	  gained 1.33 to 1.42 times at 1,000 places but 1.00 to 1.06 at 3,000,
	  so a run of twelve logarithms there came out 9% slower, the warm-up
	  call costing more than the rest repaid.
	- **Not within 2^-16 of one.** There Newton's exponential skips its
	  first chunk of sixteen bits, and gets cheap: the AGM led by 1.3 to 2.2
	  times out to 16 bits from one in every base and width measured, and by
	  24 bits was even or behind, Newton six times faster at 10^-1000.
	- **Only once its constants are held, or paid for once.** Past their
	  text a first call at a width builds pi and ln(base), and lost to
	  Newton in every base. So a width's first logarithm goes to Newton and
	  its second builds them; within the text, or at or below a width the AGM
	  has already run at in this base -- whose constants serve it, rounded --
	  every call takes the AGM. A first call is never slower than before.

	**A floor set through the hook forces the width rule alone**, for the
	tests and for measuring: nought takes the AGM everywhere, a huge one never.
*/

static long gLogAgmFloor = kLogAgmFloorMeasured;

long lzi_set_log_agm_floor( long lPlaces )
{
	long	l_previous = gLogAgmFloor;

	gLogAgmFloor = lPlaces;

	return l_previous;
}

long lzi_log_agm_floor( long lBase )
{
	if ( gLogAgmFloor != kLogAgmFloorMeasured )
		return gLogAgmFloor;

	switch ( lBase )
	{
	case 3:		return 150;
	case 4:		return 100;
	case 6:		return 75;
	case 7:		return 75;
	case 8:		return 75;
	case 10:	return 75;
	case 5:
	case 9:
	case 11:
	case 12:
	case 13:
	case 14:
	case 15:
	case 16:	return 50;
	default:	return LONG_MAX;
	}
}

// pi and ln(base) come from their text to here, with the AGM's guard and
//	margins to spare -- past it a first call builds them...
const long	kLogAgmTextDigits = 100;

// This thread's logarithms past the floor, by base and width, and the
//	widest width the AGM has run at in each base -- its constants held to
//	there. Both go with the constants, in lzu_takedown...
static thread_local std::map<std::pair<long, long>, long>	gLogAgmSeen;
static thread_local long									gLogAgmHeld[ 17 ] = {};

static void lzi_log_agm_forget()
{
	gLogAgmSeen.clear();

	for ( long& l_held : gLogAgmHeld )
		l_held = 0;
}

static bool lzi_log_agm_wanted( CLZNumber* pValue )
{
	long	l_base = pValue->m_lBase;
	long	l_width = glWorkingPrecision;

	if ( l_width < lzi_log_agm_floor( l_base ) )
		return false;

	if ( gLogAgmFloor != kLogAgmFloorMeasured )
		return true;

	double	d_value = 0.0;

	if ( lzi_value_estimate( pValue, &d_value ) &&
		 (std::fabs( d_value - 1.0 ) < std::ldexp( 1.0, -kExpSplitFirstBits )) )
		return false;

	// Held, or this width's second logarithm, which builds them -- the
	//	count is reached only when neither of the others answers...
	bool	b_held = (l_width <= CLZNumber::EquivalentPlaces( kLogAgmTextDigits, 10, l_base )) ||
					 (l_width <= gLogAgmHeld[ l_base ]) ||
					 (++gLogAgmSeen[ std::make_pair( l_base, l_width ) ] >= 2);

	if ( b_held && (l_width > gLogAgmHeld[ l_base ]) )
		gLogAgmHeld[ l_base ] = l_width;

	return b_held;
}

// The natural logarithms this thread has taken, for the tests: a whole
//	power takes none -- see lzi_logarithm_count...
static thread_local long gLogarithms = 0;

long lzi_logarithm_count()
{
	return gLogarithms;
}

long lz_log_e( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	gLogarithms++;

	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_log_e, kRetakeDoubtLog );

	// Exact at one and nowhere else, asked before anything is done to the
	//	operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 1 ) );

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long	lResult = kNoError;

	// The logarithm of zero is minus infinity; of a negative, nothing
	//	real at all. Neither was noticed before, and ln(-5) answered
	//	1659318788593735693051314843533547208079...
	leftOperand->CheckForZero();

	if ( leftOperand->IsZero() )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kUndefinedResult );

	if ( leftOperand->IsNegative() )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	// Far from one, by the argument's exponent -- see lzi_log_scaled...
	long	l_exponent = lzi_leading_exponent( leftOperand );

	if ( CLZNumber::EquivalentPlaces( (l_exponent < 0) ? -l_exponent : l_exponent, leftOperand->m_lBase, 10 ) >= gLogScaleFrom )
		return lzi_log_scaled( leftOperand, targetNumber, l_exponent );

	// By the AGM, past the base's floor, away from one, its constants held
	//	or worth building -- see lzi_log_agm_floor...
	if ( lzi_log_agm_wanted( leftOperand ) )
		return lzi_log_by_agm( leftOperand, (targetNumber != nullptr) ? targetNumber : leftOperand );

	/*
		Past a hundred and fifty digits, Newton on the exponential rather
		than the series -- see lzi_log_newton -- **where the series would
		be long**. Since 2026-10-15 the series is sized from its y, so near
		one it is short and beats Newton at any width measured: 844 us
		against 1,788 at 800 digits for 1 + 1.2E-7. Measured against it,
		Newton wins once the series wants more than about 128 terms,
		whatever the width -- 962 us against 1,235 for 2.5 in base two at
		150 digits, while 1.5625 there, at 114 terms, is the series' until
		200. So y is estimated as the reduction below would leave it, from
		the double estimate of ln x: folded into [-1, 1] and halved through
		tanh. kLogNewtonFrom was 200 until then, set before the exponential
		took binary splitting...

		**Never at the seed's own width**, which is Newton's own call for its
		start -- since 2026-10-17. The floor the hook keeps, a digit past the
		seed's sixty, was all that stood between the two, and it is digits
		where the seed is places: in base fifteen sixty digits are 52 places,
		and 52 places read back as 62 digits. So with the hook at its floor
		the seed asked Newton for itself there, and lz_log_e and
		lzi_log_newton called each other until the stack ran out -- no other
		base, and only through the hook, the floor being 150 digits...
	*/
	if ( (CLZNumber::EquivalentPlaces( glWorkingPrecision, leftOperand->m_lBase, 10 ) > gLogNewtonFrom) &&
		 (glWorkingPrecision > CLZNumber::EquivalentPlaces( kLogNewtonSeed, 10, leftOperand->m_lBase )) )
	{
		double	d_ln = 0.0;
		bool	b_long = true;

		if ( lzi_ln_estimate( leftOperand, &d_ln ) )
		{
			double	d_folded = std::fabs( d_ln );

			if ( d_folded > 1.0 )
				d_folded -= std::ceil( d_folded - 1.0 );

			double	d_y = std::tanh( d_folded / 2.0 ) * (1.0 + 1e-9);

			if ( d_y < 1e-300 )
				b_long = false;
			else
			{
				double	d_terms = ((static_cast<double>( glWorkingPrecision ) * std::log( static_cast<double>( leftOperand->m_lBase ) )) /
								   -std::log( d_y ) - 1.0) / 2.0;

				b_long = (d_terms > kLogNewtonTerms);
			}
		}

		// ...unless a test has set the width, to compare the two roads,
		//	which wants Newton everywhere it reaches -- see
		//	lzi_set_log_newton_from...
		if ( b_long || (gLogNewtonFrom != kLogNewtonFrom) )
			return lzi_log_newton( leftOperand, targetNumber );
	}

	/*
		How many terms of the series to sum.

		This used to be a flat 180, with a comment beside it saying it
		produced 32 digits for ln(10). It did, and that was the fault:
		32 digits where the precision holds 36. See the reduction
		below for why the band the argument is folded into decides
		this, and why it is now the natural base rather than the
		representation's.

		The series is in y = (x-1)/(x+1) over powers of y squared, so
		N terms leave about y^(2N+1). Folded into [1/e, e] the worst y
		is (e-1)/(e+1), a shade under .4622, and log10 of one over
		that is .3354 -- so reaching P places wants P / (2 * .3354)
		terms, which is P * 1.491. Three halves of P plus a little is
		that with room to spare.

		A flat count was wrong a second way as well: nothing tied it
		to the working precision, so a caller who raised the precision
		got the same 180 terms and a quietly worse answer. It scales
		now.
	*/

	/*
		The whole of the work below is done wider than the answer, and
		rounded back at the end.

		The reduction divides by e until the argument is under e, and
		every one of those divisions rounds. Folding into [1/e, e]
		rather than [1/base, base] costs an extra division or two for
		the same argument -- that being the point, since it is what
		makes y small -- so it also costs an extra rounding or two.
		Without guard places that shows: arccosh(10) gave up a last
		digit to the tighter band even while nine genuine faults were
		being fixed around it. With them it does not.
	*/

	const long	kLogGuardDigits = 8;

	long	l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kLogGuardDigits, leftOperand->m_lBase );

	/*
		**Three halves of the width in decimal digits, not in places.**
		The 1.491 above is a count of decimal digits a term gains, and it
		was applied to places of the operand's base: a hex place is 1.2
		digits, so in hex the series ran some 20% short. ln(1000.5) at
		420 hex places came out 8e19 units out of its last place --
		found against Newton, which does not count terms. Base ten is
		unchanged, and bases two and eight were counting over.
	*/

	long	series_depth = ((CLZNumber::EquivalentPlaces( glWorkingPrecision, leftOperand->m_lBase, 10 ) * 3) / 2) + 8;
	long	series_index = 0;

	if ( series_depth < 24 )
		series_depth = 24;

	bool	b_negateResult = false;

	CLZNumber*	one = nullptr;
	CLZNumber*	two = nullptr;
	CLZNumber*	base_reciprocal = nullptr;
	CLZNumber*	series_denom = nullptr;
	CLZNumber*	series_frac = nullptr;
	CLZNumber*	base_e = nullptr;

	CLZNumber*	y = nullptr;
	CLZNumber*	y_squared = nullptr;
	CLZNumber*	natural_log = nullptr;
	CLZNumber*	integer_log = nullptr;

	// Create working objects...
	one = new CLZNumber;
	two = new CLZNumber;
	base_reciprocal = new CLZNumber;
	series_denom = new CLZNumber;
	series_frac = new CLZNumber;
	base_e = new CLZNumber;

	y = new CLZNumber;
	y_squared = new CLZNumber;
	natural_log = new CLZNumber;
	integer_log = new CLZNumber;

	// Set bases...
	one->SetBase( leftOperand->m_lBase );
	two->SetBase( leftOperand->m_lBase );
	base_reciprocal->SetBase( leftOperand->m_lBase );
	series_denom->SetBase( leftOperand->m_lBase );
	series_frac->SetBase( leftOperand->m_lBase );

	y->SetBase( leftOperand->m_lBase );
	y_squared->SetBase( leftOperand->m_lBase );
	natural_log->SetBase( leftOperand->m_lBase );
	integer_log->SetBase( leftOperand->m_lBase );

	/*
		Whole numbers, never doubles. Assign( double ) goes by way of base
		ten and re-counts the precision into the operand's base, so a hex
		one came out at 37 places of 44 -- and x - 1, y and the whole
		series took the shorter precision from it. From 60 hex places up,
		ln of anything between 1/e and 1 fell short: 17 units of the last
		place at 60, 4.6e9 at 100, 2e19 at 150; base twelve from 150.
		These were the only two doubles assigned in the engine.
	*/

	one->Assign( 1L );
	two->Assign( 2L );

	series_denom->Assign( leftOperand );

	/*
		**Assign brought the operand's precision across with its value**,
		and the reduction below divides series_denom in place, which
		works to the precision of its left operand. So the reduction ran
		at whatever the operand happened to carry: a 20-place operand gave
		a 20-place logarithm, 4e14 units out of its 36th place, while the
		answer reported 36. The arcsine and lzi_log_gamma put the working
		precision back the same way, and say why.
	*/

	series_denom->SetPrecision( glWorkingPrecision );

	/*
		The divisor of the reduction loop, into the operand's base once,
		and as wide as the working precision -- computed past the
		constant's width, as lz_raise_e's is. See lzi_e_in_base.

		**The two have to agree.** power goes through a logarithm and an
		exponential, and while both took e from the same 120 places its
		error cancelled in the round trip: power(1e40, 2.5) came out
		exactly 1e100. With an accurate e in the exponential alone, the
		logarithm's own shortfall past 120 digits showed as 3e16 units
		of that answer's 36th place.
	*/

	lzi_e_in_base( base_e, leftOperand->m_lBase, glWorkingPrecision );

	/*
		Narrow e to the width being worked at before the reduction
		loop starts dividing by it.

		kBaseE carries 120 places for the sake of a caller who has
		raised the precision, and at the default width all but a few
		dozen of those are discarded by the first division anyway.
		Leaving them on made ln of a million take 2.95 ms against
		1.65, the reduction dividing by a 120-place divisor once per
		power of e it strips off.

		This costs no accuracy: the answer is bounded by the working
		width, and where that width is the wider one -- a caller at
		200 places -- this leaves the constant untouched.
	*/

	if ( base_e->m_lRadix > glWorkingPrecision )
	{
		base_e->SetPrecision( glWorkingPrecision );
		base_e->Round();
	}

	/*
		Reduce the argument into [1/e, e] before summing.

		The series below is in y = (x-1)/(x+1), and its terms carry
		powers of y squared, so it converges quickly while y is small
		and hardly at all as y approaches one. The band the argument
		is folded into therefore sets the worst y, and the worst y
		sets how many terms are needed.

		The band used to be [1/base, base] -- the representation's
		base, ten -- while the loop that does the folding divides by
		e. Those two do not match: dividing by e can only bring a
		value down to just above ten, not to ten, so the reduced
		argument landed anywhere in (10/e, 10] and the worst y was
		9/11, a shade over .818. The tail left after 180 terms at that
		y is about 1.9e-34, which is exactly the two or three digits
		ln(10), ln(200), ln(4000) and their neighbours were losing.
		Arguments that happened to reduce into the middle of the band
		were unaffected, which is why ln(3) and ln(100) were always
		right and the fault looked scattered.

		Folding into [1/e, e] instead matches the band to the divisor
		actually used. The worst y becomes (e-1)/(e+1), just under
		.4622, and the tail after the same 180 terms would be 5e-124.
		It is cheaper as well as better: 67 terms now reach what 258
		would have needed before, so the series does well under half
		the work even after the count was raised to track the
		precision.

		Both edges need doing, and for a long time only the upper one
		did. Nothing brought a small argument up, so ln(.01) summed
		against y = -.98 -- where the term dropped after 180 of them
		is still around 1e-3 -- and answered 4 correct digits.
		ln(.001) managed one. That was fixed by inverting rather than
		climbing; the threshold for it moves here from 1/base to 1/e,
		to keep the lower edge as tight as the upper. At the old
		threshold ln(.1) summed against y = -.818 and lost the same
		digits the top of the band did.
	*/

	lResult = lzi_first_failure( lResult, lz_divide( one, base_e, base_reciprocal ) );

	/*
		Below the band, invert rather than climb.

		The reduction used to multiply a small argument by e until it
		came up into range, which is the mirror of the loop below and
		looks right. It is not: the multiplication rounds to the
		working precision, and once the argument is smaller than that
		precision can hold, multiplying it by e and rounding gives
		nought. The loop then compares nought against a tenth for ever.
		The cliff sat one place past the working precision -- ln(1e-36)
		answered, ln(1e-37) hung -- and lz_root_3 went through it, so
		the cube root of anything that small hung too.

		ln(x) = -ln(1/x), and the reciprocal of something very small is
		something very large, which the loop below already reduces
		correctly. So invert, note that the answer wants negating, and
		let the existing path do the work.
	*/

	if ( series_denom->Compare( base_reciprocal ) == kLessThan )
	{
		lResult = lzi_first_failure( lResult, lz_reciprocal( series_denom ) );
		b_negateResult = true;
	}

	// Get down to e, which is what the loop divides by...
	while ( series_denom->Compare( base_e ) == kGreaterThan )
	{
		lResult = lzi_first_failure( lResult, lz_divide( series_denom, base_e ) );
		lz_add( integer_log, one );
	}

	lz_subtract( series_denom, one, y );
	lz_add( series_denom, one, y_squared );
	lResult = lzi_first_failure( lResult, lz_divide( y, y_squared ) );
	lz_multiply( y, y, y_squared );

	/*
		**Sized from y, not from the band's worst** -- since 2026-10-15,
		performance item 6 of that day's review. N terms leave about
		|y|^(2N+1), so the width P in places of base b wants N of about
		(P ln b / -ln|y| - 1) / 2; the count above is that for the worst
		y, .4622. An argument near one has a small y: ln(1.0000001234567)
		at 36 places summed 62 terms for the seven it needed, and ln(3.2)
		in base seven, folded to y = .081, three times its share. Two
		terms over, and never more than the count above...
	*/
	double	d_y = 0.0;

	if ( lzi_value_or_infinity( y, &d_y ) )
	{
		d_y = std::fabs( d_y ) * (1.0 + 1e-9);

		if ( d_y < 1e-300 )
			series_depth = 2;
		else if ( d_y < 1.0 )
		{
			double	d_wanted = ((static_cast<double>( glWorkingPrecision ) * std::log( static_cast<double>( leftOperand->m_lBase ) )) /
								-std::log( d_y )) - 1.0;
			long	l_wanted = static_cast<long>( std::ceil( d_wanted / 2.0 ) ) + 2;

			if ( l_wanted < 2 )
				l_wanted = 2;

			if ( l_wanted < series_depth )
				series_depth = l_wanted;
		}
	}

	// Now do the series...
	for ( series_index = series_depth; series_index >= 2; series_index-- )
	{
		series_denom->Assign( (( 2 * series_index ) - 1) );
		lResult = lzi_first_failure( lResult, lz_divide( one, series_denom, series_frac ) );
		lz_add( series_frac, natural_log, series_denom );
		lz_multiply( y_squared, series_denom, natural_log );
	}

	lz_add( natural_log, one );
	lz_multiply( natural_log, two );
	lz_multiply( natural_log, y );

	// Combine the results...
	lz_add( natural_log, integer_log );

	// ...and undo the inversion, if the argument came in below the band...
	if ( b_negateResult )
		lz_negate( natural_log );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	natural_log->SetPrecision( glWorkingPrecision );
	natural_log->Round();

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( natural_log );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( natural_log );
			leftOperand->Round();
		}
	}

	// Clean up -- base_e among them, which this list left out until
	//	2026-10-15, so every logarithm down the series path leaked a number...
    delete one;
    delete two;
    delete base_reciprocal;
    delete series_denom;
    delete series_frac;
    delete base_e;
    delete y;
    delete y_squared;
    delete natural_log;
    delete integer_log;

	return lResult;
}


static long lzi_log_n_work( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	long lResult = kNoError;

	auto logLeftBase_e = new CLZNumber;
	auto logRightBase_e = new CLZNumber;
	auto result = new CLZNumber;

	/*
		**A base near one divides by a small logarithm, so both are
		worked wider by the places one over it has.** A logarithm is good
		to its last place, absolutely, and divided by ln n that error
		grows by 1 / |ln n|: to a base of 1 + 1E-20 it is twenty places
		larger, and the answer's own size, which the width already
		allows for, does not cover it. Until 2026-10-12 nothing did, and
		HEAD was right there only by luck -- its estimate of ln n came out
		nought, and the probe that took over widened by far too much. When
		the estimate became exact the luck ran out: 238 of 400 logarithms
		to bases a hair from one were wrong, 19 of them right before.
		Now none is. The quotient itself needs only the answer's width,
		which its target already has.
	*/

	long	l_saved = glWorkingPrecision;
	double	d_lnN = 0.0;

	// Only where a whole place is lost, so an ordinary base works as it
	//	always did: lzi_places_for's two places' margin cost a base of
	//	two 16 per cent for an error ln 2 grows by 1.4...
	if ( lzi_ln_estimate( leftOperand, &d_lnN ) && (d_lnN != 0.0) )
	{
		double	d_lost = -std::log( std::fabs( d_lnN ) ) / std::log( static_cast<double>( leftOperand->m_lBase ) );

		if ( d_lost >= 1.0 )
			glWorkingPrecision += static_cast<long>( std::ceil( d_lost ) );
	}

	// A value in another base, brought into this one first, or its
	//	logarithm and the quotient are worked in its base at its width --
	//	see lzi_into_base_keeping. ln x needs x only relatively, which is
	//	what that keeps...
	CLZNumber	value;

	if ( rightOperand->m_lBase != leftOperand->m_lBase )
	{
		value.Assign( rightOperand );
		lzi_into_base_keeping( &value, leftOperand->m_lBase, glWorkingPrecision + 4 );
		rightOperand = &value;
	}

	lResult = lzi_first_failure( lResult, lz_log_e( leftOperand, logLeftBase_e ) );
	lResult = lzi_first_failure( lResult, lz_log_e( rightOperand, logRightBase_e ) );
	lResult = lzi_first_failure( lResult, lz_divide( logRightBase_e, logLeftBase_e, result ) );

	glWorkingPrecision = l_saved;

	// Assign result to target...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete logLeftBase_e;
    delete logRightBase_e;
    delete result;

	return lResult;
}


long lz_log_n( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_log_n, rightOperand, kRetakeDoubt );

	// ...and an answer that is a fraction decided exactly in doubt -- see
	//	lzi_log_decide...
	answer.DecideTies( lzi_log_decide );

	// The log of x to base n is ln x / ln n; its own logarithm sizes it...
	double	d_lnX = 0.0;
	double	d_lnN = 0.0;
	long	l_magnitude = -1;

	if ( lzi_ln_estimate( rightOperand, &d_lnX ) && lzi_ln_estimate( leftOperand, &d_lnN ) &&
		 (d_lnN != 0.0) && (d_lnX != 0.0) )
		l_magnitude = lzi_places_for( std::log( std::fabs( d_lnX / d_lnN ) ), leftOperand->m_lBase );

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	// The operands as they came, for the check below...
	CLZNumber	base;
	CLZNumber	value;

	base.Assign( leftOperand );
	value.Assign( rightOperand );

	long	lResult = lzi_widened_by_answer( lzi_log_n_work, leftOperand, rightOperand, targetNumber, l_magnitude );

	// A short logarithm is exact only if the base to it is the value --
	//	see lzi_write_to_width...
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( lzi_short_of_width( p_answer ) && lzi_power_is_exactly( &value, &base, p_answer ) )
	{
		p_answer->TrimTrailing();
		answer.Exact();
	}

	return lResult;
}


long lz_max( CLZNumber* leftOperand, CLZNumber* rightOperand )
{
	if ( lzi_is_unordered( leftOperand, rightOperand ) )
		return kUnordered;

	long lResult = 0;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		long	l_order = mixed.Compare();

		return (l_order == kLessThan) ? kRight : ((l_order == kGreaterThan) ? kLeft : kEqual);
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();


	switch( leftOperand->Compare( rightOperand ) )
	{
	case kLessThan:
		lResult = kRight;
		break;

	case kEqualTo:
		lResult = kEqual;
		break;

	case kGreaterThan:
		lResult = kLeft;
		break;
	}

	return lResult;
}


long lz_min( CLZNumber* leftOperand, CLZNumber* rightOperand )
{
	if ( lzi_is_unordered( leftOperand, rightOperand ) )
		return kUnordered;

	long lResult = 0;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		long	l_order = mixed.Compare();

		return (l_order == kLessThan) ? kLeft : ((l_order == kGreaterThan) ? kRight : kEqual);
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();


	switch( leftOperand->Compare( rightOperand ) )
	{
	case kLessThan:
		lResult = kLeft;
		break;

	case kEqualTo:
		lResult = kEqual;
		break;

	case kGreaterThan:
		lResult = kRight;
		break;
	}

	return lResult;
}


/*
	**A remainder is exact**: left less a whole multiple of right, the
	whole number taken from the raw quotient in rm_mode -- see
	lzi_divide_work. Truncated, it is the modulo and the remainder;
	nearest, ties to even, the remainder-near.

	Two things were wrong, and both are why 1 over 0.333...334 to fifty
	places had a modulo of 0 where it is 0.333...332:

	- the whole number was taken from a quotient already rounded, which
	  took 2.999...994 to 3;
	- the product of it and right was rounded at the working precision,
	  so a right operand wider than that lost places the subtraction
	  then could not give back.

	So the product is taken at right's own width, which a whole number
	times right needs and cannot exceed, and the subtraction is exact
	as it always was. The remainder keeps every place it has, as a sum
	or a difference does.
*/

static long lzi_remainder_by( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber,
							  RoundingMode rm_mode )
{
	/*
		**An infinity or a NaN carries through**, as the decimal arithmetic
		has it -- since 2026-10-15, found beside the division faults of
		that day's review. A NaN on either side, or an infinite left, is a
		NaN; a finite left beside an infinite right is the left itself, no
		whole multiple of infinity fitting in it. Until then the division
		below read a NaN or an infinite divisor as nought: remainder( 7,
		NaN ) reported a division by nought, remainder( 7, infinity )
		answered infinity, and remainder( infinity, 7 ) infinity too.
	*/
	bool	b_leftNaN = leftOperand->FlagIsOn( kfQNaN );
	bool	b_rightNaN = rightOperand->FlagIsOn( kfQNaN );
	bool	b_leftInfinite = leftOperand->FlagIsOn( kfInfinity );
	bool	b_rightInfinite = rightOperand->FlagIsOn( kfInfinity );

	if ( b_leftNaN || b_rightNaN || b_leftInfinite || b_rightInfinite )
	{
		CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
		{
			if ( !(b_leftNaN || b_rightNaN || b_leftInfinite) )
			{
				if ( p_target != leftOperand )
					p_target->Assign( leftOperand );
			}
			else
			{
				bool	b_cameFromDivideByZero = leftOperand->FlagIsOn( kfDiv0 ) || rightOperand->FlagIsOn( kfDiv0 );

				p_target->Assign( static_cast<long>(0) );
				p_target->SetFlag( kfQNaN );

				// Where the trouble started travels with it...
				if ( b_cameFromDivideByZero )
					p_target->SetFlag( kfDiv0 );
			}
		}

		return kNoError;
	}

	CLZNumber	whole;
	CLZNumber	product;

	long	lResult = lzi_divide_whole( leftOperand, rightOperand, &whole, rm_mode );
	long	l_savedPrecision = glWorkingPrecision;

	if ( rightOperand->m_lRadix > glWorkingPrecision )
		glWorkingPrecision = rightOperand->m_lRadix;

	lResult = lzi_first_failure( lResult, lz_multiply( &whole, rightOperand, &product ) );

	glWorkingPrecision = l_savedPrecision;

	return lzi_first_failure( lResult, lz_subtract( leftOperand, &product, targetNumber ) );
}


long lz_modulo( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The engine's own work: its divisions are not the caller's -- see
	//	CLZInternalWork...
	CLZInternalWork	work( (targetNumber != nullptr) ? targetNumber : leftOperand );

    long lResult = kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Remainder( lz_modulo, targetNumber );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();


	// Truncated -- see lzi_remainder_by...
	lResult = lzi_first_failure( lResult, lzi_remainder_by( leftOperand, rightOperand, targetNumber, RM_round_down ) );

    return lResult;
}


/*
	Multiplication of two digit arrays, least significant digit first,
	each digit below the base, in four tiers: the plain method a digit
	at a time for the smallest products, the plain method on packed
	limbs from there, Karatsuba on packed limbs from about 96 limbs, and
	the transform for the largest. See lzi_multiply_tier_for for where
	each takes over, and lzi_multiply_digits for the dispatch.

	**Every floor is measured in context** -- lz_multiply itself, the
	tiers forced each way -- over about 1,800 shapes, and scored the
	divider's way: past the random shapes' crossing, and asking what the
	digit loop's best cases lose.

	**Packed from a product of 48 digit-pairs, the shorter operand two
	digits or more; 96 in base two.** Random shapes from there pack at
	least 1.18 times faster (1.32 in base two), and 4.4 at the median in
	bases 8 to 16. A one-digit operand stays on the digit loop, whose
	divisions come out quicker when the quotients are tiny: a multiplier
	of exactly one is its best case, and packing loses 20-30% on it at
	any length when it is the left operand.

	**The one case the floor does not protect**: a power of the base as
	the right operand, both operands under 32 digits. The digit loop
	skips noughts in its right operand, so it does a single pass, and
	packing loses up to about a fifth. That is past the divider's "about
	a tenth", and accepted: it is under half a per cent of the
	multiplications the gates make, and the rest of the same band packs
	1.2 to 10 times faster. Base two is noisier, half of a random binary
	operand being noughts.

	**There was a fourth tier, and packing left it nothing.** The split
	kernel -- Karatsuba, three half-products in place of four, n^1.585 --
	ran from 32 digits to the transform. Over packed leaves it lost to
	packed alone at every size and floor measured, by 2 to 47 times: its
	additions and subtractions went a digit at a time with a division
	apiece, which cost more than the quarter of the products it saved.
	Retired 2026-09-24. It was also the only sub-quadratic path for a
	base with no table -- SetBase accepts 2 to 16 and only five had one
	-- so such a base multiplied a digit at a time at every length.

	**It is back on packed limbs, since 2026-09-29**, as
	lzi_karatsuba_multiply: a limb's carry is a comparison, not a
	division, and it beats packed from about 96 limbs and the transform
	to beyond 100,000 digits -- see lzi_transform_floor_for. Neither is
	reached by anything the engine does at its default: the largest
	multiplication the gates make is under 800 digits.

	**Every base has all four tiers since 2026-10-12**: three to fifteen
	have limbs -- see lzi_multiply_limb_digits_for -- and so Karatsuba and
	the transform too. The floor of 48 suits them, measured in bases
	three, seven and thirteen: packing wins by 1.15 to 2 times there, and
	a power of the base on the right loses as it does in base ten.
*/

const long kPackedMultiplyFloor = 48;
const long kPackedMultiplyFloorBinary = 96;

/*
	**Where the transform takes over, per base, on the shorter operand** --
	measured again for every base on 2026-10-12, when bases three to
	fifteen were given limbs, from 150,000 to 1.5 million bits of a square
	in steps of three per cent.

	**The crossing is a staircase, and no single floor serves all of it.**
	The transform's cost is flat between power-of-two lengths and doubles
	at each, so just past a step Karatsuba wins again, in bands that
	shrink as the lengths grow. Each floor is the one whose worst loss
	over the whole scan is least: the worst cases are those bands, 18 to
	39 per cent, and the mean loss is one to two. Base eleven's lands a
	step later than the rest, its transform element being the smallest,
	seventeen bits.

	The five old floors were set on 2026-09-29 at the first crossing, from
	a scan that ended just below the next step; past it Karatsuba won
	again, by a third at 100,000 digits in base ten. Their worst losses
	were 50 per cent in base two, 91 in eight and 48 in ten, and their
	means 3.2, 9.0 and 2.8. **A pair four times as long as it is wide goes
	to the transform from half the floor**: that halved the lopsided shapes
	that lost more than one per cent, when it was measured.
*/

static long lzi_transform_floor_for( long lBase )
{
	switch ( lBase )
	{
	case 2:		return 220280;
	case 3:		return 290995;
	case 4:		return 230608;
	case 5:		return 198635;
	case 6:		return 178423;
	case 7:		return 164288;
	case 8:		return 149261;
	case 9:		return 145497;
	case 10:	return 130870;
	case 11:	return 207710;
	case 12:	return 128653;
	case 13:	return 121008;
	case 14:	return 114184;
	case 15:	return 114613;
	case 16:	return 118763;
	default:	return 0x7FFFFFFFL;
	}
}

/*
	The most Newton passes the square root will take before giving up
	on settling. With a starting estimate within a factor of the base
	it converges in under ten; this is only here so an estimate that
	oscillates between two adjacent values cannot spin.
*/

const long kRootMaxCycles = 64;

/*
	**The widths a root's Newton iteration climbs through**, since
	2026-10-02: from about kRootSeedPlaces up to the full width, each a
	little over half the next. Newton doubles the correct digits each
	pass, so one pass at each width brings the iterate to the next
	from half of it, two places to spare. The square and cube roots
	used to work at full width every pass: about fifteen full-width
	divisions for the square root of 2 at 10,000 places, where these
	cost about two. Answers the count and fills plWidths in rising
	order, the last being lTarget itself.

	**The seed width is 64, measured.** The square root of 2 and of a
	36-digit value, best of three: from 16 or 32 the root at 36 places
	cost 8% more than working at full width, the extra passes at small
	widths costing more than they save; from 64 it is level there,
	1.2 times faster at 100 places, 3.2 at 500, 4.4 at 1,000 and 5.9 at
	3,000. Below 64 places the climb is a single width, the old way.
*/

const long kRootSeedPlaces = 64;
const long kRootMostWidths = 64;

/*
	**Whether rounding at lRound places is in doubt**: the digits from
	there to lWidth -- the guard -- within lUnits of the last place of
	the half, the one boundary rounding to nearest can cross. Since
	2026-10-02 a root's last width takes one pass and stops, and a root
	whose guard says doubt has its rounding decided exactly instead --
	see lzi_root_decide. The pass from half the width lands within a
	unit of where settling would have: the square root's matched exactly
	in 4,733 climbs of 4,736 and was a unit out in the rest, and the
	cube root's matched in 5,431 of 5,640 and was a unit out in the
	others, those being its two-pass cycles. Four units leaves room, so
	a value clear of the half rounds as the true root does. Compared as
	twice the guard against the base's power, so an odd base's half,
	which no finite digits reach, needs nothing special. A guard too
	long to read, or no guard at all, is doubt.
*/

const long kRootDoubtUnits = 4;

static bool lzi_rounding_in_doubt( const CLZNumber* pValue, long lRound, long lWidth, long lUnits )
{
	long	l_guard = lWidth - lRound;

	if ( (lRound < 0) || (l_guard <= 0) )
		return true;

	int64_t	l_digits = 0;
	int64_t	l_power = 1;

	for ( long l_place = lRound + 1; l_place <= lWidth; l_place++ )
	{
		long	l_index = pValue->m_lRadix - l_place;
		long	l_digit = ((l_index >= 0) && (l_index < pValue->m_lDigitsLength)) ? pValue->m_lpDigits[ l_index ] : 0;

		if ( l_power > (static_cast<int64_t>( 1 ) << 40) )
			return true;

		l_digits = (l_digits * pValue->m_lBase) + l_digit;
		l_power *= pValue->m_lBase;
	}

	int64_t	l_fromHalf = (2 * l_digits) - l_power;

	if ( l_fromHalf < 0 )
		l_fromHalf = -l_fromHalf;

	return l_fromHalf <= (2 * static_cast<int64_t>( lUnits ));
}


/*
	**A root's rounding decided exactly**: pRoot, near enough the root
	of pOperand to the lPower that only the side of the half at lRound
	places is in question, is cut there to y and rounded by comparing
	the half itself, raised to the power, with the operand. Both sides
	are doubled to stay finite in an odd base: with t = 2y + a unit, the
	root is past the half when 2^lPower times the operand is more than
	t^lPower, short of it when less, and on it exactly when equal, which
	goes to the even kept digit as RoundAt's half-even does. Since
	2026-10-02, when roots built a hair from a half showed that a guard
	cannot tell: 241 of 480 came out a unit wrong, settled or not. The
	operand and pRoot are positive here, a cube root's sign going on
	after.
*/

static void lzi_root_decide( CLZNumber* pRoot, const CLZNumber* pOperand, long lPower, long lRound )
{
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = pOperand->m_lRadix + lRound + 8;

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	twice;
	CLZNumber	half_power;
	CLZNumber	operand_scaled;
	CLZNumber	factor;

	y.Assign( pRoot );
	y.Round( RM_round_down, lRound );

	unit.SetBase( pRoot->m_lBase, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lRound );

	lz_add( &y, &y, &twice );
	lz_add( &twice, &unit );
	lzi_exact_power( &twice, lPower, &half_power );

	factor.SetBase( pRoot->m_lBase, false );
	factor.Assign( (lPower == 2) ? 4L : 8L );
	lz_multiply( const_cast<CLZNumber*>( pOperand ), &factor, &operand_scaled );

	long	l_against = lz_compare( &operand_scaled, &half_power );
	bool	b_up = (l_against == kGreaterThan);

	// Exactly on the half: to the even kept digit...
	if ( l_against == kEqualTo )
	{
		long	l_index = y.m_lRadix - lRound;
		long	l_kept = ((l_index >= 0) && (l_index < y.m_lDigitsLength)) ? y.m_lpDigits[ l_index ] : 0;

		b_up = (l_kept % 2) != 0;
	}

	if ( b_up )
		lz_add( &y, &unit );

	glWorkingPrecision = l_saved;

	pRoot->Assign( &y );
	pRoot->SetPrecision( l_saved );
}


static long lzi_root_widths( long lTarget, long* plWidths )
{
	long	la_falling[ kRootMostWidths ];
	long	l_count = 0;
	long	l_width = lTarget;

	la_falling[ l_count++ ] = l_width;

	while ( (l_width > kRootSeedPlaces) && (l_count < kRootMostWidths) )
	{
		l_width = (l_width / 2) + 2;
		la_falling[ l_count++ ] = l_width;
	}

	for ( long i = 0; i < l_count; i++ )
		plWidths[ i ] = la_falling[ l_count - 1 - i ];

	return l_count;
}

/*
	The transform runs over the field of integers modulo

		p = 2^64 - 2^32 + 1 = 18446744069414584321

	which is prime, and whose multiplicative group has order
	2^32 * 3 * 5 * 17 * 257 * 65537. The 2-adic part of that is what
	matters: p - 1 is divisible by 2^32, so a primitive root of unity
	exists for every power-of-two transform length up to 2^32, and 7
	is a primitive root of the whole group.

	Working in a finite field rather than the complex numbers is the
	whole point. Every step is exact modular integer arithmetic, so
	the product comes back exact with nothing to round and no error
	bound to argue. A floating-point transform would be quicker and
	would owe a proof.
*/

const uint64_t kTransformPrime = 18446744069414584321ULL;
const uint64_t kTransformGenerator = 7;

/*
	The largest transform this will attempt, and so the largest
	operand. Beyond it the caller falls through to Karatsuba on packed
	limbs.

	It is not an arbitrary ceiling: the packing factors below are
	chosen against it, and raising it would shrink them. At 2^22
	elements the factors are 20, 6, 6, 5 and 5 with at least fourfold
	headroom on the coefficient bound; at 2^24 base ten would drop to
	5 digits per element. 2^22 elements is twenty-five million
	base-ten digits, which is three orders of magnitude past anything
	this engine is asked for, so the trade goes to the packing.
*/

const long kTransformMaxLength = 1L << 22;

// The reach in force, which only the tests lower...
static long gTransformReach = kTransformMaxLength;

long lzi_set_transform_reach( long lLength )
{
	long	l_previous = gTransformReach;

	gTransformReach = ((lLength <= 0) || (lLength > kTransformMaxLength)) ? kTransformMaxLength : lLength;

	return l_previous;
}


/////////////////////////////////////////////////////////////////
//	Modular arithmetic

static uint64_t lzi_ntt_add( uint64_t a, uint64_t b )
{
	uint64_t	l_sum = a + b;

	// Wrapped, or landed at or above the modulus...
	if ( (l_sum < a) || (l_sum >= kTransformPrime) )
		l_sum -= kTransformPrime;

	return l_sum;
}

static uint64_t lzi_ntt_subtract( uint64_t a, uint64_t b )
{
	if ( a >= b )
		return a - b;

	return (a + kTransformPrime) - b;
}

/*
	**64-bit words only, so the transform builds where there is no 128-bit
	integer.** MSVC has none, and Elzed builds with it. Until 2026-09-29
	this tier used __uint128_t at three sites with no guard; each is now
	exact in 64 bits, on every compiler, so the path MSVC takes is the one
	every test here runs.

	The whole product of two words, as its high and low words, from their
	32-bit halves...
*/

static void lzi_multiply_words( uint64_t a, uint64_t b, uint64_t* pHigh, uint64_t* pLow )
{
	const uint64_t	kHalf = 0xFFFFFFFFULL;

	uint64_t	l_lowLow = (a & kHalf) * (b & kHalf);
	uint64_t	l_lowHigh = (a & kHalf) * (b >> 32);
	uint64_t	l_highLow = (a >> 32) * (b & kHalf);
	uint64_t	l_highHigh = (a >> 32) * (b >> 32);

	// The middle column, where the two cross products meet the carry out
	//	of the lowest; three numbers under 2^32 each, so it cannot wrap...
	uint64_t	l_middle = (l_lowLow >> 32) + (l_lowHigh & kHalf) + (l_highLow & kHalf);

	*pLow = (l_middle << 32) | (l_lowLow & kHalf);
	*pHigh = l_highHigh + (l_lowHigh >> 32) + (l_highLow >> 32) + (l_middle >> 32);
}

/*
	a * b modulo p, reduced by p's own shape rather than divided.

	p = 2^64 - 2^32 + 1, so 2^64 is 2^32 - 1 modulo p -- call it e -- and
	2^96 is -1. A product high * 2^64 + low, with high split as
	hh * 2^32 + hl, is therefore low - hh + hl * e, and that takes two
	corrections and a final subtraction:

	  - low - hh borrowing has added 2^64, which is e too many;
	  - adding hl * e carrying has dropped 2^64, which is e too few;
	  - what is left is under 2^64, so one subtraction of p brings it in.

	Neither correction can overflow in turn: a borrow leaves at least
	2^64 - 2^32 to take e from, and a carry leaves under e^2 to add e to.
	The remainder it replaces was a call into the compiler's runtime for
	a 128-bit division.
*/

static uint64_t lzi_ntt_multiply( uint64_t a, uint64_t b )
{
	const uint64_t	kEpsilon = 0xFFFFFFFFULL;

	uint64_t	l_high = 0;
	uint64_t	l_low = 0;

	lzi_multiply_words( a, b, &l_high, &l_low );

	uint64_t	l_highHigh = l_high >> 32;
	uint64_t	l_highLow = l_high & kEpsilon;
	uint64_t	l_value = l_low - l_highHigh;

	if ( l_low < l_highHigh )
		l_value -= kEpsilon;

	uint64_t	l_added = l_highLow * kEpsilon;
	uint64_t	l_sum = l_value + l_added;

	if ( l_sum < l_value )
		l_sum += kEpsilon;

	if ( l_sum >= kTransformPrime )
		l_sum -= kTransformPrime;

	return l_sum;
}

static uint64_t lzi_ntt_power( uint64_t base, uint64_t exponent )
{
	uint64_t	l_result = 1;

	while ( exponent )
	{
		if ( exponent & 1 )
			l_result = lzi_ntt_multiply( l_result, base );

		base = lzi_ntt_multiply( base, base );
		exponent >>= 1;
	}

	return l_result;
}

static uint64_t lzi_ntt_inverse( uint64_t a )
{
	// The field is prime, so Fermat gives the inverse directly...
	return lzi_ntt_power( a, kTransformPrime - 2 );
}

/////////////////////////////////////////////////////////////////
//	The transform itself

/*
	Iterative Cooley-Tukey, decimation in time, with the input
	permuted into bit-reversed order first. Length must be a power
	of two.
*/

static void lzi_ntt_transform( uint64_t* values, long length, bool bInverse )
{
	// Bit-reversal permutation...
	for ( long i = 1, j = 0; i < length; i++ )
	{
		long	l_bit = length >> 1;

		for ( ; j & l_bit; l_bit >>= 1 )
			j ^= l_bit;

		j ^= l_bit;

		if ( i < j )
		{
			uint64_t l_swap = values[ i ];
			values[ i ] = values[ j ];
			values[ j ] = l_swap;
		}
	}

	/*
		**Each span's powers of its root are worked once, into a table** --
		since 2026-10-15, performance item 8 of that day's review. They
		were stepped along afresh in every block, a second multiplication
		in every butterfly; tabled, the butterflies take one. The same
		products in the same field, so the same answers. Without the
		memory for the table it steps them as it did...
	*/
	uint64_t*	lp_powers = static_cast<uint64_t*>( ::malloc( static_cast<size_t>( (length > 1) ? (length / 2) : 1 ) * sizeof( uint64_t ) ) );

	for ( long l_span = 2; l_span <= length; l_span <<= 1 )
	{
		/*
			The root for this span. p - 1 is divisible by 2^32, so
			raising the generator to (p-1)/span lands on a primitive
			span-th root of unity for every power of two up to that.
		*/

		uint64_t	l_root = lzi_ntt_power( kTransformGenerator,
										 (kTransformPrime - 1) / static_cast<uint64_t>( l_span ) );

		if ( bInverse )
			l_root = lzi_ntt_inverse( l_root );

		long	l_half = l_span / 2;

		if ( lp_powers != nullptr )
		{
			lp_powers[ 0 ] = 1;

			for ( long j = 1; j < l_half; j++ )
				lp_powers[ j ] = lzi_ntt_multiply( lp_powers[ j - 1 ], l_root );

			for ( long i = 0; i < length; i += l_span )
			{
				for ( long j = 0; j < l_half; j++ )
				{
					uint64_t	l_even = values[ i + j ];
					uint64_t	l_odd = lzi_ntt_multiply( values[ i + j + l_half ], lp_powers[ j ] );

					values[ i + j ] = lzi_ntt_add( l_even, l_odd );
					values[ i + j + l_half ] = lzi_ntt_subtract( l_even, l_odd );
				}
			}

			continue;
		}

		for ( long i = 0; i < length; i += l_span )
		{
			uint64_t	l_power = 1;

			for ( long j = 0; j < l_half; j++ )
			{
				uint64_t	l_even = values[ i + j ];
				uint64_t	l_odd = lzi_ntt_multiply( values[ i + j + l_half ], l_power );

				values[ i + j ] = lzi_ntt_add( l_even, l_odd );
				values[ i + j + l_half ] = lzi_ntt_subtract( l_even, l_odd );

				l_power = lzi_ntt_multiply( l_power, l_root );
			}
		}
	}

	::free( lp_powers );

	if ( bInverse )
	{
		uint64_t	l_scale = lzi_ntt_inverse( static_cast<uint64_t>( length ) );

		for ( long i = 0; i < length; i++ )
			values[ i ] = lzi_ntt_multiply( values[ i ], l_scale );
	}
}

/////////////////////////////////////////////////////////////////
//	Packing

/*
	How many base-b digits go into one transform element.

	Two things bound it. The packed radix b^k must fit a 64-bit word,
	and -- the binding one -- every coefficient of the convolution
	must stay below the modulus, or it wraps and the answer is
	silently wrong with nothing to notice it. A coefficient is at
	most n * (b^k - 1)^2 for a transform of length n, so with the
	transform capped at kTransformMaxLength and a factor of two kept in
	hand,

		(b^k - 1)^2 <= p / 2 / kTransformMaxLength

	and k is the largest integer satisfying it: 20, 6, 6, 5 and 5 for
	bases two, eight, ten, twelve and sixteen, with fourfold headroom at
	worst. The margin is not decoration -- the arithmetic is close
	enough that a cap of 2^24 instead of 2^22 would push base ten down
	to five digits per element, and taking the bound as b^(2k) rather
	than (b^k - 1)^2 would overflow base two outright.

	**Worked out from that bound for every base**, since 2026-10-12, when
	bases three to fifteen were given the transform: 12, 10, 8, 7, 7, 6, 5,
	5, 5 and 5 for those, and the five values above, which were a table
	until then and are the same. lzi_ntt_bound_holds below re-derives the
	bound at run time from the actual lengths rather than trusting this.
*/

long lzi_ntt_pack_factor( long base )
{
	static const std::array<long, 17>	sl_factors = []
	{
		std::array<long, 17>	factors = {};
		const uint64_t			l_bound = (kTransformPrime - 1) / 2 / static_cast<uint64_t>( kTransformMaxLength );

		for ( long l_base = 2; l_base <= 16; l_base++ )
		{
			uint64_t	l_radix = 1;

			while ( ((l_radix * static_cast<uint64_t>( l_base )) - 1) * ((l_radix * static_cast<uint64_t>( l_base )) - 1) <= l_bound )
			{
				l_radix *= static_cast<uint64_t>( l_base );
				factors[ l_base ]++;
			}
		}

		return factors;
	}();

	return ((base >= 2) && (base <= 16)) ? sl_factors[ base ] : 0;
}

static uint64_t lzi_ntt_packed_radix( long base, long factor )
{
	uint64_t	l_radix = 1;

	for ( long i = 0; i < factor; i++ )
		l_radix *= static_cast<uint64_t>( base );

	return l_radix;
}

/*
	Re-derive the coefficient bound from the lengths actually in hand. A
	wrap here would not raise anything -- it would return a plausible
	wrong product -- so this is worth the handful of instructions it costs
	once per multiplication.

	largest^2 * length < p is asked as largest^2 <= (p - 1) / length,
	which is the same question exactly and has no product to overflow.
	It was asked in 128 bits until 2026-09-29. A largest digit past 2^32
	squares past p, so it fails outright.
*/

static bool lzi_ntt_bound_holds( uint64_t packedRadix, long length )
{
	uint64_t	l_largestDigit = packedRadix - 1;

	if ( (length <= 0) || (l_largestDigit > 0xFFFFFFFFULL) )
		return false;

	return (l_largestDigit * l_largestDigit) <= ((kTransformPrime - 1) / static_cast<uint64_t>( length ));
}

/////////////////////////////////////////////////////////////////
//	The kernel

static bool lzi_multiply_transform( const long* a, long lenA, const long* b, long lenB,
				   long* out, long base )
{
	if ( (lenA <= 0) || (lenB <= 0) )
		return false;

	long	l_factor = lzi_ntt_pack_factor( base );

	if ( l_factor == 0 )
		return false;

	uint64_t	l_packedRadix = lzi_ntt_packed_radix( base, l_factor );

	long	l_packedA = ((lenA + l_factor) - 1) / l_factor;
	long	l_packedB = ((lenB + l_factor) - 1) / l_factor;

	// The transform length: a power of two holding the whole product...
	long	l_length = 1;

	while ( l_length < (l_packedA + l_packedB) )
		l_length <<= 1;

	if ( l_length > gTransformReach )
		return false;

	// Belt and braces against the one failure that would be silent...
	if ( !lzi_ntt_bound_holds( l_packedRadix, l_length ) )
		return false;

	uint64_t*	lp_a = static_cast<uint64_t*>( ::calloc( l_length, sizeof(uint64_t) ) );
	uint64_t*	lp_b = static_cast<uint64_t*>( ::calloc( l_length, sizeof(uint64_t) ) );

	if ( (lp_a == nullptr) || (lp_b == nullptr) )
	{
		::free( lp_a );
		::free( lp_b );

		return false;
	}

	// Pack, most significant digit of each group highest...
	for ( long i = 0; i < l_packedA; i++ )
	{
		uint64_t	l_value = 0;

		for ( long j = l_factor - 1; j >= 0; j-- )
		{
			long	l_index = (i * l_factor) + j;

			l_value *= static_cast<uint64_t>( base );

			if ( l_index < lenA )
				l_value += static_cast<uint64_t>( a[ l_index ] );
		}

		lp_a[ i ] = l_value;
	}

	for ( long i = 0; i < l_packedB; i++ )
	{
		uint64_t	l_value = 0;

		for ( long j = l_factor - 1; j >= 0; j-- )
		{
			long	l_index = (i * l_factor) + j;

			l_value *= static_cast<uint64_t>( base );

			if ( l_index < lenB )
				l_value += static_cast<uint64_t>( b[ l_index ] );
		}

		lp_b[ i ] = l_value;
	}

	// Transform, multiply pointwise, transform back...
	lzi_ntt_transform( lp_a, l_length, false );
	lzi_ntt_transform( lp_b, l_length, false );

	for ( long i = 0; i < l_length; i++ )
		lp_a[ i ] = lzi_ntt_multiply( lp_a[ i ], lp_b[ i ] );

	lzi_ntt_transform( lp_a, l_length, true );

	/*
		The coefficients are now the exact convolution, each below the
		modulus. Carrying them into the packed radix and then writing
		the digits out is the last step. A coefficient can approach p
		and the carry from below can push the sum past what a 64-bit
		word holds, so the two are divided by the radix separately and
		put back together: the remainders add to under twice the radix,
		and at most one more goes into the quotient. This was a 128-bit
		sum until 2026-09-29.
	*/

	uint64_t	l_carry = 0;
	long		l_outLength = lenA + lenB;

	for ( long i = 0; i < l_length; i++ )
	{
		uint64_t	l_quotient = (lp_a[ i ] / l_packedRadix) + (l_carry / l_packedRadix);
		uint64_t	l_digitGroup = (lp_a[ i ] % l_packedRadix) + (l_carry % l_packedRadix);

		if ( l_digitGroup >= l_packedRadix )
		{
			l_digitGroup -= l_packedRadix;
			l_quotient++;
		}

		l_carry = l_quotient;

		// Unpack the group into base-b digits...
		for ( long j = 0; j < l_factor; j++ )
		{
			long	l_index = (i * l_factor) + j;

			if ( l_index < l_outLength )
				out[ l_index ] = static_cast<long>( l_digitGroup % static_cast<uint64_t>( base ) );

			l_digitGroup /= static_cast<uint64_t>( base );
		}
	}

	::free( lp_a );
	::free( lp_b );

	return true;
}


// The plain method: every digit of one against every digit of the other...
static void lzi_multiply_school( const long* a, long lenA, const long* b, long lenB,
								 long* out, long lBase )
{
	for ( long j = 0; j < lenB; j++ )
	{
		if ( b[ j ] == 0 )
			continue;

		long	l_carry = 0;

		for ( long i = 0; i < lenA; i++ )
		{
			long	l_current = out[ i + j ] + (a[ i ] * b[ j ]) + l_carry;

			out[ i + j ] = l_current % lBase;
			l_carry = l_current / lBase;
		}

		for ( long k = lenA; l_carry != 0; k++ )
		{
			long	l_current = out[ k + j ] + l_carry;

			out[ k + j ] = l_current % lBase;
			l_carry = l_current / lBase;
		}
	}
}


/*
	The plain method on packed limbs.

	The loop above works a digit at a time, with a remainder and a
	division by the base for every pair of digits -- the shape long
	division had before it was packed. This packs k digits into a limb
	with the divider's own helpers, lzi_packed_pack and
	lzi_packed_unpack, and multiplies limbs. Measured on the same arrays
	before it was written: ten digits by nine in 30 ns against 259, and
	forty by forty in 98 against the split kernel's 2,470.

	**The middle tier of three**, from kPackedMultiplyFloor to the
	transform; see lzi_multiply_digits. It came in four stages, as the
	divider did -- proved against Python's integers and lz_multiply by
	a harness, measured in context against every tier,
	wired in, the figures re-recorded -- and it retired the split kernel
	on the way.

	**One carry per column, not one per product.** Column c of the
	product is the sum of left[i] * right[c - i]; it is summed whole and
	split into a limb and a carry once, where the loop above splits after
	every product. Worth 1.65 times at forty digits, and nothing at ten.

	**And no column can overflow, at any length.** Every product is
	below (L - 1)^2 for a limb base L, and 64 bits hold only so many of
	them: sixteen in base eight, eighteen in base ten, ninety-nine in
	base twelve, 256 in base sixteen. So the sum is folded after every
	so many products -- its carry set aside and only its limb kept --
	and how many is worked out from L, so that no sum between two folds
	can reach 2^64. The carry into a column is split the same way before
	the column starts, since the bound assumes a sum below one limb. The
	transform tier's notes call its coefficient bound the thing that
	fails silently; here there is no bound to fail.

	**The limbs are the divider's, except in base two.** 31 bits there
	leaves room for four products between folds, and a fold is a
	division. At 20 bits no column a test can afford ever folds, and it
	measured 10 to 20 per cent faster at every size tried, from ten
	digits to six hundred; 16 to 28 bits were all within a few per cent
	of it. Pack and unpack, a bit at a time, are most of base two's cost
	at any limb size. Elsewhere the table was fastest or level: a digit
	fewer per limb was slower in bases eight, twelve and sixteen, and
	level in base ten.

	**Short operands stay on the stack.** A buffer of
	kPackedMultiplyStackLimbs saves a constant 20 to 35 ns against an
	allocation: 40 to 60 per cent of a product of a few digits, 13 per
	cent at forty, 3 at 150. The divider's allocation measured 0 to 6
	per cent; this kernel is quick enough for the same cost to matter.

	**It writes exactly lLeftLength + lRightLength digits** rather than
	adding into what is there, and answers false for a base it has no
	table for, an empty operand, or a failed allocation.
*/

// The factor and the fold in force, which only the tests change, and
//	what the tests read back to know the paths were taken. The settings
//	are the process's and the counts each thread's own...
static long gPackedMultiplyFactor = 0;
static long gPackedMultiplyFold = 0;
static thread_local long gPackedMultiplications = 0;
static thread_local long gPackedMultiplyFolds = 0;

// How many limbs the kernel keeps on the stack before it allocates...
const long kPackedMultiplyStackLimbs = 256;

long lzi_set_packed_multiply_factor( long lFactor )
{
	long	l_previous = gPackedMultiplyFactor;

	gPackedMultiplyFactor = lFactor;

	return l_previous;
}

long lzi_set_packed_multiply_fold( long lProducts )
{
	long	l_previous = gPackedMultiplyFold;

	gPackedMultiplyFold = lProducts;

	return l_previous;
}

void lzi_packed_multiply_counts( long* plMultiplications, long* plFolds )
{
	if ( plMultiplications != nullptr )
		*plMultiplications = gPackedMultiplications;

	if ( plFolds != nullptr )
		*plFolds = gPackedMultiplyFolds;
}

// Digits per limb: the table, or fewer if the tests asked. Never more,
//	which keeps a limb below 2^32 and a product inside 64 bits. Base two
//	is 20 rather than the divider's 31; see above...
static long lzi_packed_multiply_factor( long lBase )
{
	long	l_factor = 0;

	switch ( lBase )
	{
	case 2:		l_factor = kPackedDigitsBinaryMultiply;	break;
	case 8:		l_factor = kPackedDigitsOctal;			break;
	case 10:	l_factor = kPackedDigitsDecimal;		break;
	case 12:	l_factor = kPackedDigitsDuodecimal;		break;
	case 16:	l_factor = kPackedDigitsHex;			break;
	default:	l_factor = lzi_multiply_limb_digits( lBase );	break;
	}

	if ( (gPackedMultiplyFactor > 0) && (gPackedMultiplyFactor < l_factor) )
		l_factor = gPackedMultiplyFactor;

	return l_factor;
}

/*
	The column loop on limbs already packed: writes exactly lLeft +
	lRight limbs of the product, least significant first. Shared by
	lzi_packed_multiply and, below its floor, the Karatsuba kernel.
*/

static void lzi_packed_multiply_limbs( const uint64_t* lpLeft, long lLeftLimbs,
									   const uint64_t* lpRight, long lRightLimbs,
									   uint64_t* lpProduct, uint64_t lLimbBase, long lFold )
{
	long		l_productLimbs = lLeftLimbs + lRightLimbs;
	uint64_t	l_carry = 0;

	for ( long c = 0; c < l_productLimbs; c++ )
	{
		long	l_first = ((c - lRightLimbs) + 1 > 0) ? (c - lRightLimbs) + 1 : 0;
		long	l_last = (c < lLeftLimbs - 1) ? c : lLeftLimbs - 1;

		// The carry in, split so the sum starts below a limb...
		uint64_t	l_sum = l_carry % lLimbBase;
		uint64_t	l_over = l_carry / lLimbBase;

		for ( long i = l_first; i <= l_last; )
		{
			long	l_stop = ((i + lFold) - 1 < l_last) ? (i + lFold) - 1 : l_last;

			for ( ; i <= l_stop; i++ )
				l_sum += lpLeft[ i ] * lpRight[ c - i ];

			// ...and folded before the next run could overflow it...
			if ( i <= l_last )
			{
				l_over += l_sum / lLimbBase;
				l_sum %= lLimbBase;
				gPackedMultiplyFolds++;
			}
		}

		lpProduct[ c ] = l_sum % lLimbBase;
		l_carry = l_over + (l_sum / lLimbBase);
	}
}

bool lzi_packed_multiply( const long* lpLeft, long lLeftLength,
						  const long* lpRight, long lRightLength,
						  long lBase, long* lpProduct )
{
	long	l_factor = lzi_packed_multiply_factor( lBase );

	if ( (l_factor == 0) || (lLeftLength <= 0) || (lRightLength <= 0) )
		return false;

	uint64_t	l_limbBase = 1;

	for ( long i = 0; i < l_factor; i++ )
		l_limbBase *= static_cast<uint64_t>( lBase );

	/*
		The fold: as many products as a column can take on top of a
		limb's worth of sum without reaching 2^64. Capped, because at the
		tests' smallest limbs it is astronomical and a long has to hold
		it -- 32 bits of one on Windows.
	*/

	uint64_t	l_largest = (l_limbBase - 1) * (l_limbBase - 1);
	uint64_t	l_room = (UINT64_MAX - l_limbBase) / l_largest;
	long		l_fold = (l_room > (1UL << 30)) ? (1L << 30) : static_cast<long>( l_room );

	if ( (gPackedMultiplyFold > 0) && (gPackedMultiplyFold < l_fold) )
		l_fold = gPackedMultiplyFold;

	long	l_leftLimbs = ((lLeftLength + l_factor) - 1) / l_factor;
	long	l_rightLimbs = ((lRightLength + l_factor) - 1) / l_factor;

	// Both operands and the product, on the stack while they fit...
	uint64_t	l_stack[ kPackedMultiplyStackLimbs ];
	long		l_need = 2 * (l_leftLimbs + l_rightLimbs);
	uint64_t*	lp_block = l_stack;

	if ( l_need > kPackedMultiplyStackLimbs )
	{
		lp_block = static_cast<uint64_t*>( ::calloc( l_need, sizeof(uint64_t) ) );

		if ( lp_block == nullptr )
			return false;
	}

	uint64_t*	lp_left = lp_block;
	uint64_t*	lp_right = lp_left + l_leftLimbs;
	uint64_t*	lp_product = lp_right + l_rightLimbs;

	l_leftLimbs = lzi_packed_pack( lpLeft, lLeftLength, lBase, l_factor, lp_left );
	l_rightLimbs = lzi_packed_pack( lpRight, lRightLength, lBase, l_factor, lp_right );

	long		l_productLimbs = l_leftLimbs + l_rightLimbs;

	gPackedMultiplications++;

	lzi_packed_multiply_limbs( lp_left, l_leftLimbs, lp_right, l_rightLimbs, lp_product, l_limbBase, l_fold );

	lzi_packed_unpack( lp_product, l_productLimbs, lBase, l_factor,
					   lpProduct, lLeftLength + lRightLength );

	if ( lp_block != l_stack )
		::free( lp_block );

	return true;
}

/*
	**Karatsuba on packed limbs**, lz_multiply's third tier since 2026-09-29,
	between packed schoolbook and the transform -- see lzi_multiply_tier_for
	for where.

	The digit-level split kernel was removed on 2026-09-24 because packing
	left it no band: its additions went a digit at a time with a division
	apiece. This one splits the packed limbs, so its additions and
	subtractions are a limb at a time with a comparison apiece -- a limb is
	below 2^32, so any sum of two is below 2^33, and a carry is "take L
	off". Below its floor it hands the limbs to lzi_packed_multiply_limbs,
	the packed kernel's own column loop, so the leaves are exactly the
	packed multiply.

	For x = x1 L^h + x0 and y = y1 L^h + y0, three products instead of
	four: z0 = x0 y0, z2 = x1 y1, and (x0 + x1)(y0 + y1) - z0 - z2 for the
	middle. A lopsided pair -- the shorter no longer than half the longer,
	rounded up -- is multiplied a chunk of the longer at a time instead,
	each chunk the shorter's length and added in at its place.

	**The floor must be at least four limbs**, or it never stops: the middle
	product's operands are h + 1 limbs long, and only from four is that
	shorter than the operands it came from.
*/

const long kKaratsubaLimbFloor = 96;		// measured: within 2% of every base's best
const long kKaratsubaLeastFloor = 4;

// The floor is the process's, the split count each thread's own...
static long gKaratsubaFloor = kKaratsubaLimbFloor;
static thread_local long gKaratsubaSplits = 0;

long lzi_set_karatsuba_floor( long lLimbs )
{
	long	l_previous = gKaratsubaFloor;

	if ( lLimbs <= 0 )
		gKaratsubaFloor = kKaratsubaLimbFloor;
	else
		gKaratsubaFloor = (lLimbs < kKaratsubaLeastFloor) ? kKaratsubaLeastFloor : lLimbs;

	return l_previous;
}

long lzi_karatsuba_splits()
{
	return gKaratsubaSplits;
}

// Into lpSum, lpA plus lpB, limbs in base lLimbBase; answers the sum's
//	length, which is the longer's or one more...
static long lzi_limbs_add( const uint64_t* lpA, long lA, const uint64_t* lpB, long lB,
						   uint64_t* lpSum, uint64_t lLimbBase )
{
	if ( lA < lB )
	{
		const uint64_t*	lp_swap = lpA;
		long			l_swap = lA;

		lpA = lpB;
		lA = lB;
		lpB = lp_swap;
		lB = l_swap;
	}

	uint64_t	l_carry = 0;

	for ( long i = 0; i < lA; i++ )
	{
		uint64_t	l_value = lpA[ i ] + ((i < lB) ? lpB[ i ] : 0) + l_carry;

		l_carry = (l_value >= lLimbBase) ? 1 : 0;
		lpSum[ i ] = l_value - (l_carry * lLimbBase);
	}

	if ( l_carry == 0 )
		return lA;

	lpSum[ lA ] = 1;

	return lA + 1;
}

// lpA less lpB, in place; lpA is known to be the larger...
static void lzi_limbs_subtract( uint64_t* lpA, long lA, const uint64_t* lpB, long lB,
								uint64_t lLimbBase )
{
	uint64_t	l_borrow = 0;

	for ( long i = 0; (i < lA) && ((i < lB) || (l_borrow != 0)); i++ )
	{
		uint64_t	l_take = ((i < lB) ? lpB[ i ] : 0) + l_borrow;

		l_borrow = (lpA[ i ] < l_take) ? 1 : 0;
		lpA[ i ] = (lpA[ i ] + (l_borrow * lLimbBase)) - l_take;
	}
}

// lpInto from limb lAt upward, plus lpAdd, carried within lInto limbs.
//	Anything of lpAdd past the end is nought, the product's true length
//	being shorter than the arrays that hold it...
static void lzi_limbs_add_at( uint64_t* lpInto, long lInto, long lAt, const uint64_t* lpAdd, long lAdd,
							  uint64_t lLimbBase )
{
	uint64_t	l_carry = 0;
	long		i = 0;

	for ( ; (i < lAdd) && ((lAt + i) < lInto); i++ )
	{
		uint64_t	l_value = lpInto[ lAt + i ] + lpAdd[ i ] + l_carry;

		l_carry = (l_value >= lLimbBase) ? 1 : 0;
		lpInto[ lAt + i ] = l_value - (l_carry * lLimbBase);
	}

	for ( long k = lAt + i; (l_carry != 0) && (k < lInto); k++ )
	{
		uint64_t	l_value = lpInto[ k ] + l_carry;

		l_carry = (l_value >= lLimbBase) ? 1 : 0;
		lpInto[ k ] = l_value - (l_carry * lLimbBase);
	}
}

// lpProduct[ 0 .. lX + lY ) = lpX times lpY, limbs in base lLimbBase...
static void lzi_karatsuba_limbs( const uint64_t* lpX, long lX, const uint64_t* lpY, long lY,
								 uint64_t* lpProduct, uint64_t lLimbBase, long lFold )
{
	if ( lX < lY )
	{
		const uint64_t*	lp_swap = lpX;
		long			l_swap = lX;

		lpX = lpY;
		lX = lY;
		lpY = lp_swap;
		lY = l_swap;
	}

	if ( lY < gKaratsubaFloor )
	{
		lzi_packed_multiply_limbs( lpX, lX, lpY, lY, lpProduct, lLimbBase, lFold );
		return;
	}

	long	l_half = (lX + 1) / 2;

	// Lopsided: the shorter against the longer a chunk of its own length at a time...
	if ( lY <= l_half )
	{
		std::vector<uint64_t>	l_part( 2 * lY );

		std::fill( lpProduct, lpProduct + lX + lY, 0 );

		for ( long l_at = 0; l_at < lX; l_at += lY )
		{
			long	l_chunk = ((lX - l_at) < lY) ? (lX - l_at) : lY;

			lzi_karatsuba_limbs( lpX + l_at, l_chunk, lpY, lY, l_part.data(), lLimbBase, lFold );
			lzi_limbs_add_at( lpProduct, lX + lY, l_at, l_part.data(), l_chunk + lY, lLimbBase );
		}

		return;
	}

	gKaratsubaSplits++;

	long	l_highX = lX - l_half;
	long	l_highY = lY - l_half;

	std::vector<uint64_t>	l_sumX( l_half + 1 );
	std::vector<uint64_t>	l_sumY( l_half + 1 );
	std::vector<uint64_t>	l_middle( (2 * l_half) + 2 );

	long	l_lengthX = lzi_limbs_add( lpX, l_half, lpX + l_half, l_highX, l_sumX.data(), lLimbBase );
	long	l_lengthY = lzi_limbs_add( lpY, l_half, lpY + l_half, l_highY, l_sumY.data(), lLimbBase );

	// z0 in the low 2h limbs, z2 above it: the product's own array...
	lzi_karatsuba_limbs( lpX, l_half, lpY, l_half, lpProduct, lLimbBase, lFold );
	lzi_karatsuba_limbs( lpX + l_half, l_highX, lpY + l_half, l_highY, lpProduct + (2 * l_half), lLimbBase, lFold );

	// ...and the middle, less both, added in h limbs up...
	lzi_karatsuba_limbs( l_sumX.data(), l_lengthX, l_sumY.data(), l_lengthY, l_middle.data(), lLimbBase, lFold );

	lzi_limbs_subtract( l_middle.data(), l_lengthX + l_lengthY, lpProduct, 2 * l_half, lLimbBase );
	lzi_limbs_subtract( l_middle.data(), l_lengthX + l_lengthY, lpProduct + (2 * l_half), l_highX + l_highY, lLimbBase );

	lzi_limbs_add_at( lpProduct, lX + lY, l_half, l_middle.data(), l_lengthX + l_lengthY, lLimbBase );
}

bool lzi_karatsuba_multiply( const long* lpLeft, long lLeftLength,
							 const long* lpRight, long lRightLength,
							 long lBase, long* lpProduct )
{
	long	l_factor = lzi_packed_multiply_factor( lBase );

	if ( (l_factor == 0) || (lLeftLength <= 0) || (lRightLength <= 0) )
		return false;

	uint64_t	l_limbBase = 1;

	for ( long i = 0; i < l_factor; i++ )
		l_limbBase *= static_cast<uint64_t>( lBase );

	// The packed kernel's fold, worked out the same way...
	uint64_t	l_largest = (l_limbBase - 1) * (l_limbBase - 1);
	uint64_t	l_room = (UINT64_MAX - l_limbBase) / l_largest;
	long		l_fold = (l_room > (1UL << 30)) ? (1L << 30) : static_cast<long>( l_room );

	if ( (gPackedMultiplyFold > 0) && (gPackedMultiplyFold < l_fold) )
		l_fold = gPackedMultiplyFold;

	long	l_leftLimbs = ((lLeftLength + l_factor) - 1) / l_factor;
	long	l_rightLimbs = ((lRightLength + l_factor) - 1) / l_factor;

	// Its vectors throw where the packed kernel's calloc answers null, and
	//	both have to decline the same way, so the dispatch can hand on...
	try
	{
		std::vector<uint64_t>	l_left( l_leftLimbs );
		std::vector<uint64_t>	l_right( l_rightLimbs );

		l_leftLimbs = lzi_packed_pack( lpLeft, lLeftLength, lBase, l_factor, l_left.data() );
		l_rightLimbs = lzi_packed_pack( lpRight, lRightLength, lBase, l_factor, l_right.data() );

		std::vector<uint64_t>	l_product( l_leftLimbs + l_rightLimbs );

		lzi_karatsuba_limbs( l_left.data(), l_leftLimbs, l_right.data(), l_rightLimbs, l_product.data(), l_limbBase, l_fold );

		lzi_packed_unpack( l_product.data(), l_leftLimbs + l_rightLimbs, lBase, l_factor,
						   lpProduct, lLeftLength + lRightLength );
	}
	catch ( const std::bad_alloc& )
	{
		return false;
	}

	return true;
}

// The floors in force, which only the tests change, and what the tests
//	read back to know which tier did the work. The settings are the
//	process's and the counts each thread's own...
static long gPackedMultiplyFloor = kPackedFloorMeasured;
static long gTransformFloor = kPackedFloorMeasured;
static thread_local long gMultiplyTiers[ 4 ] = { 0, 0, 0, 0 };

long lzi_set_packed_multiply_floor( long lProduct )
{
	long	l_previous = gPackedMultiplyFloor;

	gPackedMultiplyFloor = lProduct;

	return l_previous;
}

long lzi_set_transform_floor( long lDigits )
{
	long	l_previous = gTransformFloor;

	gTransformFloor = lDigits;

	return l_previous;
}

long lzi_retake_count()
{
	return CLZAnswerBase::RetakeCount();
}


void lzi_multiply_tier_counts( long* plDigits, long* plPacked, long* plTransform, long* plKaratsuba )
{
	if ( plKaratsuba != nullptr )
		*plKaratsuba = gMultiplyTiers[ kMultiplyKaratsuba ];

	if ( plDigits != nullptr )
		*plDigits = gMultiplyTiers[ kMultiplyByDigits ];

	if ( plPacked != nullptr )
		*plPacked = gMultiplyTiers[ kMultiplyPacked ];

	if ( plTransform != nullptr )
		*plTransform = gMultiplyTiers[ kMultiplyByTransform ];
}

/*
	Which tier a product of these lengths belongs to. One function, so
	that the dispatch below and the tests pinning the floors ask the
	same question -- the floors are measured, and a test that re-derived
	them would pin its own copy rather than the engine's.

	A tier that has no table for the base is passed over. An overridden
	floor is a plain product for packing and a plain length for the
	transform: the tests use it to force a tier at any size, the
	transform's included, since nothing the tests can afford reaches it
	otherwise.
*/

long lzi_multiply_tier_for( long lBase, long lLeftLength, long lRightLength )
{
	long	l_shorter = (lLeftLength < lRightLength) ? lLeftLength : lRightLength;
	long	l_longer = (lLeftLength < lRightLength) ? lRightLength : lLeftLength;

	// In 64 bits: a long is 32 on Windows, and this product need not fit...
	int64_t	l_product = static_cast<int64_t>( lLeftLength ) * lRightLength;

	if ( lzi_ntt_pack_factor( lBase ) > 0 )
	{
		if ( gTransformFloor != kPackedFloorMeasured )
		{
			if ( l_shorter >= gTransformFloor )
				return kMultiplyByTransform;
		}
		else
		{
			long	l_floor = lzi_transform_floor_for( lBase );

			if ( (l_shorter >= l_floor) || ((l_shorter >= (l_floor / 2)) && ((l_longer / 4) >= l_shorter)) )
				return kMultiplyByTransform;
		}
	}

	long	l_factor = lzi_packed_multiply_factor( lBase );

	/*
		**Karatsuba from its floor, 96 limbs, when the pair is near square,
		and from twice it at any shape.** A lopsided pair goes a chunk of
		the longer at a time, each chunk the shorter's length, and a chunk
		that splits only once does not pay for itself: just past the floor
		a lopsided pair lost up to 9 per cent, and from twice it gained
		15 to 32. The floor is the kernel's own, so the tests' hook moves
		both at once.
	*/

	if ( l_factor > 0 )
	{
		int64_t	l_shortLimbs = (static_cast<int64_t>( l_shorter ) + l_factor - 1) / l_factor;
		int64_t	l_longLimbs = (static_cast<int64_t>( l_longer ) + l_factor - 1) / l_factor;
		int64_t	l_floor = gKaratsubaFloor;

		if ( (l_shortLimbs >= l_floor) && ((l_longLimbs < (2 * l_shortLimbs)) || (l_shortLimbs >= (2 * l_floor))) )
			return kMultiplyKaratsuba;
	}

	if ( l_factor > 0 )
	{
		if ( gPackedMultiplyFloor != kPackedFloorMeasured )
		{
			if ( l_product >= gPackedMultiplyFloor )
				return kMultiplyPacked;
		}
		else if ( (l_shorter >= 2) &&
				  (l_product >= ((lBase == 2) ? kPackedMultiplyFloorBinary : kPackedMultiplyFloor)) )
			return kMultiplyPacked;
	}

	return kMultiplyByDigits;
}

/*
	The dispatch. out must hold lenA + lenB zeroed digits: the digit
	loop adds into it, where the other three write it.

	A tier that declines -- the transform past its reach or short of
	memory, Karatsuba or the packed kernel short of memory -- hands on to
	the one below, and the digit loop needs nothing, so every product is
	made. The transform's reach is 2^22 elements, some 25 million digits
	between the two operands in base ten.
*/

static void lzi_multiply_digits( const long* a, long lenA, const long* b, long lenB,
								 long* out, long lBase )
{
	if ( (lenA <= 0) || (lenB <= 0) )
		return;

	long	l_tier = lzi_multiply_tier_for( lBase, lenA, lenB );

	/*
		**A one-digit operand packs too, but for a left multiplier of
		exactly one** -- since 2026-10-15, performance item 8 of that day's
		review. The rule above keeps every one-digit operand on the digit
		loop, after a multiplier of one, the loop's best case, lost 20-30%
		packed on the left. Measured again across lengths 50 to 2,000 and
		bases two, ten and sixteen: every other one-digit operand packs 1.4
		to 5.3 times faster, on either side, and a left one still loses 2
		to 18% in bases ten and sixteen. Past the same floor, and only
		while the floor is the measured one, so the tests' hook still moves
		it alone...
	*/
	if ( (l_tier == kMultiplyByDigits) && (gPackedMultiplyFloor == kPackedFloorMeasured) && ((lenA == 1) || (lenB == 1)) &&
		 (lzi_packed_multiply_factor( lBase ) > 0) &&
		 ((static_cast<int64_t>( lenA ) * lenB) >= ((lBase == 2) ? kPackedMultiplyFloorBinary : kPackedMultiplyFloor)) &&
		 !((lenA == 1) && (a[ 0 ] == 1)) )
		l_tier = kMultiplyPacked;

	if ( (l_tier == kMultiplyByTransform) && lzi_multiply_transform( a, lenA, b, lenB, out, lBase ) )
	{
		gMultiplyTiers[ kMultiplyByTransform ]++;
		return;
	}

	// Beneath a declined transform as well: past its reach, packed alone
	//	would be quadratic on tens of millions of digits...
	if ( ((l_tier == kMultiplyKaratsuba) || (l_tier == kMultiplyByTransform)) &&
		 lzi_karatsuba_multiply( a, lenA, b, lenB, lBase, out ) )
	{
		gMultiplyTiers[ kMultiplyKaratsuba ]++;
		return;
	}

	if ( (l_tier != kMultiplyByDigits) && lzi_packed_multiply( a, lenA, b, lenB, lBase, out ) )
	{
		gMultiplyTiers[ kMultiplyPacked ]++;
		return;
	}

	gMultiplyTiers[ kMultiplyByDigits ]++;
	lzi_multiply_school( a, lenA, b, lenB, out, lBase );
}

// The dispatch for CLZNumber's base conversion, which holds digit arrays
//	rather than numbers to hand lz_multiply...
void lzi_multiply_digit_arrays( const long* lpLeft, long lLeftLength, const long* lpRight, long lRightLength,
								long* lpProduct, long lBase )
{
	lzi_multiply_digits( lpLeft, lLeftLength, lpRight, lRightLength, lpProduct, lBase );
}


long lz_multiply( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The target's rounding settings are its own -- see CLZTargetRounding...
	CLZTargetRounding	targetRounding( (targetNumber != nullptr) ? targetNumber : leftOperand );

	// Before anything -- see CLZMatchedBase.
	//	A digit array carries no base, so a mixed pair has to be
	//	brought into one before any of it is read...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Multiply( targetNumber );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	// An infinity or a NaN carries through rather than being lost...
	if ( lzi_propagate_flags( leftOperand, rightOperand, targetNumber, kOpMultiply ) )
		return kNoError;

	long	l_lengthLeft = 0;

	long	l_lengthRight = 0;



	CLZNumber*	p_result = nullptr;

	// Create the result...
	p_result = new CLZNumber;

	/*
		The digits below are worked out in the operands' base, but a
		fresh number is base ten, and the target takes its base from
		whatever is assigned to it. Without this the answer came back
		with the right digits wearing the wrong base, and the next
		operation on it either refused as a mismatch or read the
		digits as decimal. Both are zero here, so this only labels
		them -- there is nothing yet to convert.
	*/

	p_result->SetBase( leftOperand->m_lBase, true );

	l_lengthLeft = leftOperand->m_lDigitsLength;
	l_lengthRight = rightOperand->m_lDigitsLength;

	/*
		The whole product at once, rather than a partial product per
		digit of the right operand added in one at a time. The digits
		go straight into the result's own array, which has to be big
		enough and empty first.
	*/

	long	l_productLength = l_lengthLeft + l_lengthRight;

	if ( l_productLength > p_result->m_lDigitsArrayLength )
		p_result->SetDigitsArraySize( l_productLength, kNoPreserve );

	::memset( p_result->m_lpDigits, 0, p_result->m_lDigitsArrayLength * sizeof(long) );

	lzi_multiply_digits( leftOperand->m_lpDigits, l_lengthLeft,
						 rightOperand->m_lpDigits, l_lengthRight,
						 p_result->m_lpDigits, leftOperand->m_lBase );

	p_result->m_lDigitsLength = l_productLength;

	// Determine result sign...
	if ( leftOperand->m_lSign == rightOperand->m_lSign )
		p_result->m_lSign = kPositive;
	else
		p_result->m_lSign = kNegative;

	// Calculate result radix...
	p_result->m_lRadix = leftOperand->m_lRadix + rightOperand->m_lRadix;

	// Assign to target, and round in the target's own mode -- the product
	//	is exact, so every mode rounds it correctly, at the precision or
	//	at the units if the target asks for a whole number. See
	//	CLZTargetRounding...
	{
		CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
		{
			p_target->Assign( p_result );
			targetRounding.RestoreMode();
			p_target->TrimAll();

			if ( targetRounding.Whole() )
				p_target->Round( RD_digit_zero );
			else
				p_target->Round();
		}
	}

	// Clean up...
    delete p_result;

	return kNoError;
}


long lz_negate( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

    // Record the result...
    if ( targetNumber != nullptr )
    {
        if ( !(targetNumber->IsConstant()) )
        {
            targetNumber->Assign( leftOperand );

            switch(targetNumber->m_lSign)
            {
                case kPositive:
                    targetNumber->m_lSign = kNegative;
                    break;

                case kNegative:
                    targetNumber->m_lSign = kPositive;
                    break;

                case kZero:
                default:
                    break;
            }
        }
    }
    else
    {
        if ( !(leftOperand->IsConstant()) )
        {
            leftOperand->CheckForZero();

            switch(leftOperand->m_lSign)
            {
                case kPositive:
                    leftOperand->m_lSign = kNegative;
                    break;

                case kNegative:
                    leftOperand->m_lSign = kPositive;
                    break;

                case kZero:
                default:
                    break;
            }
        }
    }

    return kNoError;
}


long lz_normalize( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( leftOperand );
			targetNumber->TrimAll();
		}
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->TrimAll();

	return kNoError;
}


/*
	**permutations, P(n, k) = n! / (n - k)!**, the ways of arranging k of
	n things in order -- added 2026-10-19, on combinations' terms: whole
	numbers only, a fraction or a negative refused as factorial refuses
	it, and arranging more than there are is nought.

	**The k whole numbers that end at n, multiplied out as a balanced
	tree** -- see lzi_range_product -- with nothing divided. There is no
	symmetry to take, as combinations takes C(n, n - k), and none is
	wanted: the product is the answer, so it holds nothing larger, where
	n! over (n - k)! would work out (n - k)! twice over only to divide it
	away. Nor is there a way by primes worth having: combinations gains
	from one because the alternative divides by k!, and here nothing does.
*/
long lz_permutations( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	// An infinity or a not-a-number in either is a not-a-number...
	if ( lzi_propagate_as_nan( leftOperand, rightOperand, targetNumber ) )
		return kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	leftOperand->CheckForZero();
	rightOperand->CheckForZero();

	if ( leftOperand->IsNegative() || rightOperand->IsNegative() ||
		 !(leftOperand->IsInteger()) || !(rightOperand->IsInteger()) )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

	// More than there are can be arranged no ways at all...
	if ( rightOperand->Compare( leftOperand ) == kGreaterThan )
	{
		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		return kNoError;
	}

	/*
		Sized before it is begun -- see lzi_answer_size -- with n - k
		taken exactly. k may be as large as n here, and nothing keeps the
		estimates of two near values in order where their digit counts
		differ, each adding its own exponent: n - k + 1 worked from them
		could fall far below nought, and lzi_ln_gamma_estimate counts up
		from there to ten by ones. Measured over 8,443 near pairs, 3,408
		of them across a digit count, they never crossed; this does not
		rest on it...
	*/
	CLZNumber	before;
	double		d_n = 0.0;
	double		d_lnN = 0.0;
	double		d_k = 0.0;
	double		d_before = 0.0;
	long		l_beyond = kNoError;

	lz_subtract( leftOperand, rightOperand, &before );

	if ( lzi_value_or_infinity( leftOperand, &d_n ) && lzi_value_or_infinity( rightOperand, &d_k ) &&
		 lzi_value_or_infinity( &before, &d_before ) && lzi_ln_estimate( leftOperand, &d_lnN ) &&
		 lzi_answer_beyond( lzi_answer_size( lzi_ln_falling_estimate( d_n, d_lnN, d_k, d_before + 1.0 ),
											 leftOperand->m_lBase ),
							leftOperand, targetNumber, &l_beyond ) )
		return l_beyond;

	// ...and taken in the operand's base -- see CLZAnswerBase. Under the
	//	limit the count is far inside a long, P(n, k) being at least k!;
	//	one that is not was past it...
	CLZNumber	result;
	long		l_count = 0;

	if ( !lzi_whole_as_long( rightOperand, &l_count ) )
		return lzi_domain_failure( leftOperand, targetNumber, kfInfinity, kTooLarge );

	result.SetBase( leftOperand->m_lBase );
	result.Assign( 1L );

	// (n - k + 1) ... n...
	if ( l_count > 0 )
		lzi_range_product( &before, 1, l_count, &result );

	if ( !(p_target->IsConstant()) )
	{
		p_target->Assign( &result );
		p_target->TrimAll();
		p_target->Round();
	}

	return kNoError;
}


long lz_plus( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	// While this operation is a no-op, we still have
	//   to copy the "result" into the target operand
	//   if a separate target is supplied...

	// Record the result...
	if ( targetNumber != nullptr )
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( leftOperand );

	return kNoError;
}


static long lzi_power_work( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Flags, noughts and negative bases were settled by lz_power before
	//	this was called -- see lzi_power_edges...

	/*
		An exponent in another base, brought into this one first, or the
		work below goes on in its base at its width: 10^0.1 with 0.1 in
		base three was right to 20 digits of 36. y ln x wants this width
		past the point, so y's error is multiplied by ln x -- which the
		answer's guard and four places more absorb until ln x passes
		base^8, 2^18 in base two: an x of a quarter of a million bits and
		an answer as long. Adding ln x's own places was tried and no case
		within reach could tell -- see lzi_into_base_keeping...
	*/
	CLZNumber	exponent;

	if ( rightOperand->m_lBase != leftOperand->m_lBase )
	{
		exponent.Assign( rightOperand );
		lzi_into_base_keeping( &exponent, leftOperand->m_lBase, glWorkingPrecision + 4 );
		rightOperand = &exponent;
	}

	CLZNumber* integer = nullptr;
	CLZNumber* fraction = nullptr;
	CLZNumber* result = nullptr;

	integer = new CLZNumber;
	fraction = new CLZNumber;
	result = new CLZNumber;

	lz_integer( rightOperand, integer );
	lz_fraction( rightOperand, fraction );

	// Calculate the integer power...
	if ( integer->IsZero() )
		integer->Assign( static_cast<long>(1) );
	else
		lzi_raise_n( leftOperand, integer, integer );

	// Calculate the fractional power...
	if ( fraction->IsZero() )
		fraction->Assign( static_cast<long>(1) );
	else
	{
		lz_log_e( leftOperand, result );
		lz_multiply( fraction, result );
		lz_exp_e( fraction );
	}

	// Calculate the result...
	lz_multiply( integer, fraction, result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete integer;
    delete fraction;
    delete result;

	return kNoError;
}


long lz_power( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// x to a half is root_2's, for x past nought -- see lz_root_n. Twice
	//	the exponent is one exactly: a half is finite in an even base only...
	if ( leftOperand->IsPositive() && !(leftOperand->FlagIsOn( kfInfinity ) || leftOperand->FlagIsOn( kfQNaN )) &&
		 !(rightOperand->FlagIsOn( kfInfinity ) || rightOperand->FlagIsOn( kfQNaN )) && ((rightOperand->m_lBase % 2) == 0) )
	{
		CLZNumber	twice;
		long		l_saved = glWorkingPrecision;

		glWorkingPrecision = rightOperand->m_lRadix + 2;
		lz_add( rightOperand, rightOperand, &twice );
		glWorkingPrecision = l_saved;

		if ( lzi_exactly( &twice, 1 ) )
			return lz_root_2( leftOperand, targetNumber );
	}

	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_power, rightOperand, kRetakeDoubt );

	// ...and a root of whole order decided exactly in doubt -- see
	//	lzi_whole_order_root_decide...
	answer.DecideTies( lzi_power_decide );

	// The edges every x^y shares, answered before any work -- see
	//	lzi_power_edges...
	long	l_edge = kNoError;

	// An edge answers exactly, or flagged -- see lzi_write_to_width...
	if ( lzi_power_edges( leftOperand, rightOperand, targetNumber, false, &l_edge ) )
	{
		answer.Exact();
		return l_edge;
	}

	// x^y has about y ln|x| of logarithm...
	double	d_ln = 0.0;
	double	d_y = 0.0;
	long	l_magnitude = -1;

	if ( lzi_ln_estimate( leftOperand, &d_ln ) && lzi_value_or_infinity( rightOperand, &d_y ) )
	{
		long	l_beyond = kNoError;

		// Refused past the limit, nought below the width -- see
		//	lzi_answer_size...
		if ( lzi_answer_beyond( lzi_answer_size( d_y * d_ln, leftOperand->m_lBase ),
								leftOperand, targetNumber, &l_beyond ) )
			return l_beyond;

		l_magnitude = lzi_places_for( d_y * d_ln, leftOperand->m_lBase );
	}

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	// The operands as they came, for the check below; the answer may be
	//	written over the left...
	CLZNumber	base;
	CLZNumber	exponent;

	base.Assign( leftOperand );
	exponent.Assign( rightOperand );

	long	lResult = lzi_widened_by_answer( lzi_power_work, leftOperand, rightOperand, targetNumber, l_magnitude );

	// A short answer is checked against the operands before it is called
	//	exact, and spelled without its noughts if it is -- see
	//	lzi_write_to_width...
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( lzi_short_of_width( p_answer ) && lzi_power_is_exactly( p_answer, &base, &exponent ) )
	{
		p_answer->TrimTrailing();
		answer.Exact();
	}

	return lResult;
}


long lz_quantize( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	if ( lzi_propagate_as_nan( leftOperand, rightOperand, targetNumber ) )
		return kNoError;

	CLZNumber* p_quantized = nullptr;

	// The quantum's places, read before anything is written: the target
	//	may be the quantum itself, and until 2026-10-15 quantize( 0.7, 3 )
	//	into the 3 read the places of the 0.7 just copied over it...
	long	l_places = rightOperand->m_lRadix;

	if ( targetNumber != nullptr )
	{
		p_quantized = targetNumber;
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( leftOperand );
	}
	else
		p_quantized = leftOperand;

	if ( !(p_quantized->IsConstant()) )
	{
		// If the left exponent is less than the right, pad the left...
		if ( p_quantized->m_lRadix < l_places )
			p_quantized->PadRight( l_places - p_quantized->m_lRadix );

		// If the left exponent is greater than the right, round the left...
		else if ( p_quantized->m_lRadix > l_places )
			p_quantized->Round( l_places );
	}

	return kNoError;
}


long lz_raise_e( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_raise_e, kRetakeDoubtRaiseE );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	/*
		Carried and refused as lz_exp_e and lz_raise_n are. This took
		neither: an infinity or a not-a-number answered 1, unflagged,
		and a fraction was stripped with success reported, so e^2.5
		answered e^2. A whole exponent is what this function is for;
		lz_exp_e takes the rest.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	if ( !(leftOperand->IsInteger()) )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	// Refused past the limit, nought below the width, before anything is
	//	worked out -- see lzi_answer_size...
	double	d_exponent = 0.0;
	long	l_beyond = kNoError;

	if ( lzi_value_or_infinity( leftOperand, &d_exponent ) &&
		 lzi_answer_beyond( lzi_answer_size( d_exponent, leftOperand->m_lBase ),
							leftOperand, targetNumber, &l_beyond ) )
		return l_beyond;

	CLZNumber* base_e = nullptr;
	CLZNumber* exponent = nullptr;
	CLZNumber* one = nullptr;
	CLZNumber* result = nullptr;

	// Wider than the answer, rounded back at the end. This is where
	//	kBaseE's thirteen unused places were being thrown away on the
	//	first squaring -- see kRaiseEGuardDigits...
	long	l_savedPrecision = glWorkingPrecision;

	long	l_guard = lzi_guard_places( kRaiseEGuardDigits, leftOperand->m_lBase );

	glWorkingPrecision += l_guard;

	base_e = new CLZNumber;
	exponent = new CLZNumber;
	one = new CLZNumber( static_cast<long>(1) );
	result = new CLZNumber;

	/*
		e is the left operand of the whole squaring chain, so it works in
		the operand's base -- see CLZAnswerBase. Left in base ten, the
		chain ran there, and the magnitude the probe below reads off its
		answer was a count of decimal digits applied as places of the
		operand's base. As wide as the working precision, computed past
		the constant's width -- see lzi_e_in_base...
	*/

	lzi_e_in_base( base_e, leftOperand->m_lBase, glWorkingPrecision );
	one->SetBase( leftOperand->m_lBase );

	/*
		Narrow e to the width actually being worked at before the
		squaring starts.

		kBaseE carries 120 places and the squaring runs at the working
		precision, so without this the first multiplication is 120
		digits by 120 and everything past the working width is thrown
		away immediately anyway -- paid for and discarded. Rounding it
		here costs nothing in accuracy, the answer being limited by
		the working width in either case, and took sinh back from
		0.13 ms to 0.09 when the constant was widened.
	*/

	if ( base_e->m_lRadix > glWorkingPrecision )
	{
		base_e->SetPrecision( glWorkingPrecision );
		base_e->Round();
	}

	exponent->Assign( leftOperand );
	exponent->StripFraction();

	/*
		This carried its own copy of the counting loop that lz_raise_n
		used to have -- one multiplication for every unit of the
		exponent, which is what made a large one impossible rather
		than merely slow. e to the n is a power like any other, so it
		hands the work to lz_raise_n and gets the squaring with it.
	*/

	/*
		**The magnitude is worth learning before the real work, not
		after.**

		How big e^n comes out decides how wide this has to run -- see
		below -- but the size of an answer costs nothing to find. So
		the chain is run once at a dozen places purely to be measured,
		and then once properly at the width that measurement asks for.
		Squaring at twelve places against fifty or a hundred and fifty
		is close to free, where computing the whole thing twice at full
		width cost 2.5 times at e^50 and 3.4 at e^100.

		Reading the magnitude off the answer rather than working it out
		as n/ln(base) also keeps this base-agnostic, GetOrderOfMagnitude
		counting digits in the number's own base.
	*/

	long	l_magnitude = 0;

	/*
		Not for a small exponent, which cannot produce an answer the
		guard does not already cover. e^n has magnitude n*log(e) in the
		number's own base, and the base that reaches sixteen soonest is
		two, at n = 11 -- so below that the probe could only ever say
		"no widening needed" and exp(1) paid 1.8 times its cost to hear
		it.
	*/

	CLZNumber	probeFrom;

	probeFrom.Assign( kRaiseEProbeFrom );

	if ( exponent->CompareABS( &probeFrom ) != kLessThan )
	{
		long	l_probePrecision = glWorkingPrecision;

		glWorkingPrecision = kRaiseEProbePlaces;

		CLZNumber*	probe_e = new CLZNumber( kBaseE );
		CLZNumber*	probe = new CLZNumber;

		lzi_constant_in_base( probe_e, kBaseE, leftOperand->m_lBase, glWorkingPrecision );
		probe_e->SetPrecision( glWorkingPrecision );
		probe_e->Round();

		lzi_raise_n( probe_e, exponent, probe );

		l_magnitude = probe->GetOrderOfMagnitude();

		delete probe_e;
		delete probe;

		glWorkingPrecision = l_probePrecision;
	}

	/*
		**A large answer needs more places than the guard provides,
		and the squaring chain cannot make them.**

		Narrowing e above fixes its *relative* accuracy: at
		glWorkingPrecision places a value near 2.7 carries about that
		many significant digits and no more. Squaring preserves
		relative error, so the answer comes out with the same relative
		accuracy however large it grows -- about the caller's precision
		plus fifteen, measured, whether the answer is 1e2 or 1e87.

		Meanwhile the places the caller asked for afford the answer its
		own magnitude on top. e^100 is 2.7e43, so 36 places has room
		for 80 significant digits and the chain delivered 50. The
		shortfall was the magnitude less what the guard already
		covered: 7 digits at e^50, 30 at e^100, 73 at e^200, and the
		same at every precision -- the signature of a fixed relative
		width rather than a working one.

		So the magnitude is read off the first answer and the whole
		thing done again wider. Same shape as lz_root_3's band and
		lz_gamma's logarithm.

		**From half the guard, not from past it.** This used to wait for
		the answer to outgrow the guard, as though the guard were there
		to hold the answer's size. It is there for the squaring chain,
		and an answer as large as the guard spent all of it: e^35 and
		e^36 missed by up to 2.27 units in base ten, e^43 and e^44 by
		5.65 in hex. Below half, what is left is well over what the
		chain spends, a few digits at those powers; widening from the
		first integer place instead would take e again, wider, for
		every small power -- 17% on e^30.
	*/

	if ( l_magnitude > (l_guard / 2) )
	{
		glWorkingPrecision += l_magnitude;

		// base_e was narrowed to the old width above, so it has to be
		//	taken again to gain by this -- from the constant, or computed
		//	past its width. See lzi_e_in_base...
		lzi_e_in_base( base_e, leftOperand->m_lBase, glWorkingPrecision );

		if ( base_e->m_lRadix > glWorkingPrecision )
		{
			base_e->SetPrecision( glWorkingPrecision );
			base_e->Round();
		}
	}

	lzi_raise_n( base_e, exponent, result );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete base_e;
    delete exponent;
    delete one;
    delete result;

	return kNoError;
}


static long lzi_raise_n( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	CLZNumber* exponent = nullptr;
	CLZNumber* result = nullptr;

	exponent = new CLZNumber;
	result = new CLZNumber;

	exponent->Assign( rightOperand );
	exponent->StripFraction();

	/*
		The loop below pairs the running product with the operand, and
		an operation reads its right operand's digits in the left's
		base, so the product is put into the operand's base.
	*/

	result->SetBase( leftOperand->m_lBase, true );
	result->Assign( "1" );

	// Raise n...
	if ( exponent->IsZero() )
	{
		result->Assign( static_cast<long>(1) );
	}

	else
	{
		/*
			Squaring, rather than multiplying the operand by itself
			the exponent's own number of times. That loop cost one
			multiplication per unit of exponent -- e to the two
			thousandth took a hundred thousand operations, and a
			larger exponent simply never finished. This takes one
			squaring per digit of the exponent in binary, so the same
			answer costs about eleven.

			**The parity comes from the halving, not the last digit.**
			This said every base the engine allows is even, so the
			exponent is odd exactly when its last digit is -- and
			SetBase has always taken three, five, seven and the rest.
			Seven is 10 in base seven, whose last digit is even, so e
			to the 7.25 there answered e to the 6.25 and 2 to the 7
			answered 64. The remainder of the halving is the parity in
			every base, and the halving was being done anyway.
		*/

		bool	b_reciprocal = exponent->IsNegative();

		if ( b_reciprocal )
			lz_negate( exponent );

		/*
			**Worked as many places wider as the exponent has**, and
			rounded back to the width it came in at. Each squaring and
			each multiplication rounds, and the chain carries those
			roundings into the answer about as many times over as the
			exponent: x^1000 came out as far as 1.75e-2 of a unit from
			its value, at the default and at 200 places alike, where
			x^25 came out 6e-4. Until 2026-10-10 nothing here allowed for
			it, and near a half that was past the retake's window -- 19
			of 470 hundredth and thousandth powers built there rounded
			wrong, and lz_power and lz_exp_n, which come here for their
			exponent's whole part, reached 1.5e-2. Now x^1000 is good to
			5e-5 of a unit, like x^25. See the note at
			lzi_widened_by_answer for why this was once left out.
		*/

		long	l_entryWidth = glWorkingPrecision;
		double	d_lnExponent = 0.0;

		if ( lzi_ln_estimate( exponent, &d_lnExponent ) )
			glWorkingPrecision += lzi_places_for( d_lnExponent, leftOperand->m_lBase );

		// A multiplication rounds to the working precision whatever its
		//	operands carry, so raising that is all the chain needs...
		auto	running = new CLZNumber;
		auto	scratch = new CLZNumber;

		running->SetBase( leftOperand->m_lBase, true );
		running->Assign( leftOperand );

		/*
			**A negative power of a fraction takes the reciprocal first,
			and powers that.** It used to build the positive power and
			divide one by it at the end, which loses digits for a
			fraction: the power is then tiny, and held to the working
			places it keeps fewer significant digits the smaller it gets
			-- as many fewer as the answer has whole digits. The
			reciprocal had nothing to give them back, so 0.0012345 to the
			-20th, 59 whole digits, came out wrong from its 47th digit,
			off by about 10^16; 0.5 to the -100th, which is 2^100
			exactly, was wrong in its 18th place; and 0.001 to the -2nd
			in base two came out a unit low -- 1,570 of 2,715 random
			cases, by up to past 10^400 units. Turned round, a fraction's
			reciprocal grows as it is powered, and a growing value keeps
			every place it is worked to.

			**A base past one keeps the old order, for speed alone.** It
			is right either way -- its reciprocal's power is a small
			answer whose places are all that count, and the two orders
			answered identically on 2,959 such bases, the gate
			byte-identical with either -- but a short base like 3.7
			stays short while it is squared, and its reciprocal does not:
			taken first, 3.7 to the -37th cost 30 per cent more.
		*/

		bool		b_reciprocalFirst = false;
		CLZNumber*	p_unit = nullptr;

		if ( b_reciprocal )
		{
			// One, built once for the comparison and whichever division
			//	follows, and only for a negative exponent...
			p_unit = new CLZNumber;
			p_unit->SetBase( leftOperand->m_lBase, false );
			p_unit->Assign( 1L );

			if ( leftOperand->CompareABS( p_unit ) == kLessThan )
			{
				// **A division rounds to its target's precision**, not
				//	the working one as a multiplication does, and the
				//	running value took the operand's with it. Held
				//	there, a reciprocal's error is multiplied by the
				//	exponent: 1 - 1E-12 to the -10^12th came out 3 units
				//	low...
				running->SetPrecision( glWorkingPrecision );
				lz_divide( p_unit, leftOperand, running );
				b_reciprocalFirst = true;
			}
		}

		scratch->SetBase( leftOperand->m_lBase, true );

		while ( !(exponent->IsZero()) )
		{
			/*
				Short division in the exponent's own base, done here
				rather than through lz_divide_integer: this runs once
				per binary digit of the exponent, and that function
				allocates working values and performs a full division
				for what is a halving. It has no fraction to worry
				about either -- StripFraction saw to that.

				It also used to be that lz_divide_integer could not be
				used at all, handing back a base-ten answer whatever
				base it was given. That is fixed; the reason to keep
				this here is now only the cost.
			*/

			long	l_carry = 0;

			for ( long i = exponent->m_lDigitsLength - 1; i >= 0; i-- )
			{
				long	l_current = (l_carry * exponent->m_lBase) + exponent->m_lpDigits[ i ];

				exponent->m_lpDigits[ i ] = l_current / 2;
				l_carry = l_current % 2;
			}

			// What is left over from halving is the parity...
			if ( l_carry != 0 )
				lz_multiply( result, running );

			while ( (exponent->m_lDigitsLength > 1) &&
					(exponent->m_lpDigits[ exponent->m_lDigitsLength - 1 ] == 0) )
				exponent->m_lDigitsLength--;

			exponent->CheckForZero();

			// One squaring too many would double the work of the last
			//	step for nothing...
			if ( !(exponent->IsZero()) )
			{
				lz_square( running, scratch );
				running->Assign( scratch );
			}
		}

		delete running;
		delete scratch;

		if ( b_reciprocal && !b_reciprocalFirst )
		{
			// The counting is over; the reciprocal is taken against
			//	the product, in the operand's base as one is. Its answer
			//	is under one, so only its places count, and the
			//	product's precision holds them...
			lz_divide( p_unit, result, result );
		}

		delete p_unit;

		// The exponent's places never leave...
		glWorkingPrecision = l_entryWidth;

		result->SetPrecision( glWorkingPrecision );
		result->Round();
	}

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetNumber->Round();
		}
	}
	else
	{
		if ( !(leftOperand->IsConstant()) )
		{
			leftOperand->Assign( result );
			leftOperand->Round();
		}
	}

	// Clean up...
    delete exponent;
    delete result;

	return kNoError;
}


long lz_raise_n( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_raise_n, rightOperand, kRetakeDoubt );

	/*
		The edges every x^y shares, answered before any work -- see
		lzi_power_edges. This one takes a whole exponent and so refuses
		a fraction: it used to strip it and report success, answering 4
		for 2^2.5 and 1 for (-8)^0.5.

		The refusal is here and not in lzi_raise_n, which the engine's
		own functions call with exponents they have already made whole.
	*/

	long	l_edge = kNoError;

	// An edge answers exactly, or flagged -- see lzi_write_to_width...
	if ( lzi_power_edges( leftOperand, rightOperand, targetNumber, true, &l_edge ) )
	{
		answer.Exact();
		return l_edge;
	}

	// x^n, n whole, has about n ln|x| of logarithm...
	double	d_ln = 0.0;
	double	d_n = 0.0;
	long	l_magnitude = -1;

	if ( lzi_ln_estimate( leftOperand, &d_ln ) && lzi_value_or_infinity( rightOperand, &d_n ) )
	{
		long	l_beyond = kNoError;

		// Refused past the limit, nought below the width -- see
		//	lzi_answer_size...
		if ( lzi_answer_beyond( lzi_answer_size( std::trunc( d_n ) * d_ln, leftOperand->m_lBase ),
								leftOperand, targetNumber, &l_beyond ) )
			return l_beyond;

		l_magnitude = lzi_places_for( std::trunc( d_n ) * d_ln, leftOperand->m_lBase );
	}

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	// The operands as they came, for the check below; the answer may be
	//	written over the left...
	CLZNumber	base;
	CLZNumber	exponent;

	base.Assign( leftOperand );
	exponent.Assign( rightOperand );

	long	lResult = lzi_widened_by_answer( lzi_raise_n, leftOperand, rightOperand, targetNumber, l_magnitude );

	// A short answer is checked against the operands before it is called
	//	exact, and spelled without its noughts if it is -- see
	//	lzi_write_to_width...
	CLZNumber*	p_answer = (targetNumber != nullptr) ? targetNumber : leftOperand;

	if ( lzi_short_of_width( p_answer ) && lzi_power_is_exactly( p_answer, &base, &exponent ) )
	{
		p_answer->TrimTrailing();
		answer.Exact();
	}

	return lResult;
}


long lz_reciprocal( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// One division of the operand, which writes short only what is exact,
	//	and a small answer as wide as a quotient's -- see lzi_write_to_width...
	answer.Exact();
	answer.QuotientWidth();

	long lResult = kNoError;

	auto result = new CLZNumber;
	auto one = new CLZNumber;

	// In the operand's base -- see CLZAnswerBase...
	one->SetBase( leftOperand->m_lBase );
	one->Assign( static_cast<long>(1) );

	/*
		**One division, rounded once, at the answer's width** -- see
		bAnswerWidth at lzi_divide_work. Until 2026-10-03 this divided
		kQuotientGuardDigits wider, which the division answered a place
		wider still, and rounded the quotient again to the answer's width:
		two roundings, and 10 of 200 answers built within a millionth of a
		unit of a half came back a unit wrong. The raw quotient carries
		every digit the rounding needs, so it is rounded at the width the
		answer has, from those.

		The rule for a quotient is applied again after. Until 2026-10-10
		that trimmed a nought: the width was read off the quotient before
		rounding, so a quotient a unit short of a power of the base rounded
		up to it a place wider than the rule then asks. The division reads
		the width again itself now, so this changes nothing, and is kept
		as every function ending in a quotient keeps it.
	*/

	if ( !lzi_propagate_flags( one, leftOperand, result, kOpDivide ) )
		lResult = lzi_divide_work( one, leftOperand, result, glWorkingPrecision, false, RM_round_default, true );

	lzi_round_quotient_answer( result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete one;

	return lResult;
}


long lz_remainder( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The engine's own work: its divisions are not the caller's -- see
	//	CLZInternalWork...
	CLZInternalWork	work( (targetNumber != nullptr) ? targetNumber : leftOperand );

	long lResult = kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Remainder( lz_remainder, targetNumber );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();


	// Truncated -- see lzi_remainder_by...
	lResult = lzi_first_failure( lResult, lzi_remainder_by( leftOperand, rightOperand, targetNumber, RM_round_down ) );

	return lResult;
}


long lz_remainder_near( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The engine's own work: its divisions are not the caller's -- see
	//	CLZInternalWork...
	CLZInternalWork	work( (targetNumber != nullptr) ? targetNumber : leftOperand );

	long lResult = kNoError;

	// A mixed pair comes into one base first -- see CLZMatchedBase...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Remainder( lz_remainder_near, targetNumber );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();


	// To the nearest, ties to even -- see lzi_remainder_by...
	lResult = lzi_first_failure( lResult, lzi_remainder_by( leftOperand, rightOperand, targetNumber, RM_round_half_even ) );

	return lResult;
}


/*
	Where a value's leading digit sits, as a power of the base: the
	value is about base to this. Positive for anything at or above
	one, negative below it. GetOrderOfMagnitude clamps at nought and
	so cannot say how far below one a value is, which is exactly what
	the square root's starting estimate needs to know.
*/

static long lzi_leading_exponent( CLZNumber* value )
{
	long	l_top = value->m_lDigitsLength - 1;

	while ( (l_top >= 0) && (value->m_lpDigits[ l_top ] == 0) )
		l_top--;

	return l_top - value->m_lRadix;
}


/*
	**Rounding an answer that is a quotient.** Division answers wider
	than it is asked -- one place past the precision, and two
	significant digits past it for a small quotient -- and a function
	that handed a quotient straight
	back handed that back too: a 37th place on an ordinary answer,
	unrounded, and on a small one places that are real accuracy. The
	cube root's note records what rounding those away cost.

	Rich's rule, 2026-09-25: **36 places for an answer of 0.1 or more**,
	as everywhere else, **and 36 significant digits below it** -- the
	working precision plus the noughts after the point. A small answer
	keeps up to 36 significant digits and no answer carries an
	unrounded tail. It is in the working precision's own base, the
	noughts counted there too. A zero, an infinity and a not-a-number
	are left alone.

	The square root ends in a division as well and does not come here,
	which is Rich's call too (2026-10-01) -- see the end of lz_root_2.
*/

// The places the rule above gives a non-zero answer: the working
//	precision, and below a base-th its noughts as well...
static long lzi_quotient_answer_places( CLZNumber* pAnswer )
{
	long	l_places = glWorkingPrecision;
	long	l_exponent = lzi_leading_exponent( pAnswer );

	if ( l_exponent < -1 )
		l_places += -l_exponent - 1;

	return l_places;
}

static void lzi_round_quotient_answer( CLZNumber* pAnswer )
{
	if ( pAnswer->FlagIsOn( kfInfinity ) || pAnswer->FlagIsOn( kfQNaN ) )
		return;

	pAnswer->CheckForZero();

	if ( pAnswer->IsZero() )
		return;

	long	l_places = lzi_quotient_answer_places( pAnswer );

	pAnswer->SetPrecision( l_places );
	pAnswer->Round();

	/*
		**A carry to a power of the base moves the leading digit up a
		place, and the rule then gives one place fewer.** The width above
		is read before rounding, so 0.0999...9 rounding to 0.1 was written
		to 37 places where 0.1 takes 36 -- by every function ending here,
		in every base, until 2026-10-10, though never in value: rounding
		at the narrower width carries the same way, so reading the width
		again only trims a nought. lz_reciprocal applied the rule a second
		time for this from 2026-10-03, which this now does for everyone.
	*/

	lzi_trim_after_carry( pAnswer );
}

/*
	**The trim takes noughts off; it does not round**, since 2026-10-15.
	An answer has places past the rule's only after a carry, when
	everything past its leading digit is nought, so trimming and writing
	back to the rule's places gives what rounding there gave -- without a
	second rounding recorded over the first, whose tail is the one a
	retake must read. See CLZAnswerBase::InDoubt.
*/

static void lzi_trim_after_carry( CLZNumber* pAnswer )
{
	if ( pAnswer->FlagIsOn( kfInfinity ) || pAnswer->FlagIsOn( kfQNaN ) )
		return;

	pAnswer->CheckForZero();

	if ( pAnswer->IsZero() )
		return;

	long	l_places = lzi_quotient_answer_places( pAnswer );

	if ( pAnswer->m_lRadix <= l_places )
		return;

	pAnswer->SetPrecision( l_places );
	pAnswer->TrimTrailing();

	if ( pAnswer->m_lRadix < l_places )
		pAnswer->PadRight( l_places - pAnswer->m_lRadix );
}

// One digit of one and noughts besides: a power of the base, either sign...
static bool lzi_is_power_of_base( CLZNumber* pValue )
{
	long	l_ones = 0;

	for ( long i = 0; i < pValue->m_lDigitsLength; i++ )
	{
		long	l_digit = pValue->m_lpDigits[ i ];

		if ( l_digit == 0 )
			continue;

		if ( (l_digit != 1) || (++l_ones > 1) )
			return false;
	}

	return l_ones == 1;
}


/*
	**A short answer is an exact one**, since 2026-10-01 -- Rich's call,
	the decimal arithmetic standards' way. An answer written with fewer
	places than its width ended there; one that did not end is written
	out to its width, trailing noughts and all. Until then the spelling
	followed the path the answer took: sin(1e-42) was 0 and arcsin(1e-42)
	0.000...0, both rounded noughts, exp(1e-42) was 1, and a quotient a
	hair from 22.5 was 22.5.

	CLZAnswerBase calls this as every function returns, unless the
	function marked its answer exact. The width is the working
	precision, or for a small answer of a function ending in a quotient
	the rule for one -- see lzi_round_quotient_answer. A flagged answer
	is left alone.
*/

static void lzi_write_to_width( CLZNumber* pAnswer, bool bQuotientWidth )
{
	if ( pAnswer->FlagIsOn( kfInfinity ) || pAnswer->FlagIsOn( kfQNaN ) )
		return;

	long	l_places = glWorkingPrecision;

	if ( bQuotientWidth )
	{
		pAnswer->CheckForZero();

		if ( !(pAnswer->IsZero()) )
			l_places = lzi_quotient_answer_places( pAnswer );
	}

	if ( pAnswer->m_lRadix < l_places )
		pAnswer->PadRight( l_places - pAnswer->m_lRadix );
}


// The value is exactly lValue as it was given: nothing rounded, nothing
//	flagged. Asked of an operand before anything is done to it...
static bool lzi_exactly( CLZNumber* pValue, long lValue )
{
	if ( (pValue == nullptr) || pValue->FlagIsOn( kfInfinity ) || pValue->FlagIsOn( kfQNaN ) )
		return false;

	CLZNumber	value;

	value.SetBase( pValue->m_lBase, false );
	value.Assign( lValue );

	return lz_compare( pValue, &value ) == kEqualTo;
}


// Fewer places than the working precision, trailing noughts aside: an
//	answer that could be exact, and is worth asking...
static bool lzi_short_of_width( CLZNumber* pAnswer )
{
	if ( pAnswer->FlagIsOn( kfInfinity ) || pAnswer->FlagIsOn( kfQNaN ) )
		return false;

	long	l_places = pAnswer->m_lRadix;

	for ( long i = 0; (i < pAnswer->m_lRadix) && (i < pAnswer->m_lDigitsLength) && (pAnswer->m_lpDigits[ i ] == 0); i++ )
		l_places--;

	return l_places < glWorkingPrecision;
}


// pBase to the whole lPower, nought or more, worked wide enough that
//	nothing is rounded...
static void lzi_exact_power( CLZNumber* pBase, long lPower, CLZNumber* pResult )
{
	CLZNumber	base;

	base.Assign( pBase );
	base.TrimTrailing();

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = (base.m_lRadix * lPower) + 1;

	CLZNumber	power;
	CLZNumber	square;

	power.SetBase( base.m_lBase, false );
	power.Assign( 1L );
	square.Assign( &base );
	square.SetPrecision( glWorkingPrecision );

	for ( long n = lPower; n > 0; n >>= 1 )
	{
		if ( n & 1 )
			lz_multiply( &power, &square );

		if ( n > 1 )
			lz_multiply( &square, &square );
	}

	glWorkingPrecision = l_saved;

	pResult->Assign( &power );
}


/*
	**Whether pValue is exactly pBase to the power pExponent**, for an
	exponent p/q with q no more than 64 and p no more than 4096 either
	way: pValue^q and pBase^p compared, both worked exactly, with a
	negative p moved across as pValue^q * pBase^-p = 1. The roots, the
	powers, the logarithm to a base and the base's exponential ask it of
	a short answer before calling it exact; any other exponent is taken
	as inexact, which costs only a written-out spelling.
*/

/*
	**The places of y^k, or of 1 / y^k, from y's own** -- since
	2026-10-15, performance item 3 of that day's review: the check on a
	short answer raised its operands exactly to as much as the 4096th
	power, and raise_n( x, 4000 ) with an answer ending in a nought took
	4.1 ms against 112 microseconds beside it. Unequal places rule the
	power out without working it.

	y = Y / b^t, Y with its own trailing noughts taken off, t shifted to
	match. b is the product of l^e over its primes, and Y^k ends in the
	least over them of floor( k v_l(Y) / e ) noughts, so y^k has k t less
	that many places, or none. 1 / y^k is b^(kt) / Y^k, finite only when Y
	is made of the base's primes -- and when it is, it has the most over
	them of ceil( k v_l(Y) / e ) less k t places; when it is not, the
	power is not exact anyway. v_l(Y) is read off Y's last digits, b^J
	being a multiple of l^j once J e reaches j. Answers -1 when it cannot
	say: a nought, or a valuation past what 64 bits hold.
*/

// The exact powers worked out to check a short answer, for the tests --
//	see lzi_exactness_proofs...
static thread_local long gExactnessProofs = 0;

long lzi_exactness_proofs()
{
	return gExactnessProofs;
}

static long lzi_places_of_power( CLZNumber* pValue, long lPower, bool bReciprocal )
{
	pValue->CheckForZero();

	if ( pValue->IsZero() || pValue->FlagIsOn( kfInfinity ) || pValue->FlagIsOn( kfQNaN ) || (lPower < 0) )
		return -1;

	long	l_base = pValue->m_lBase;
	long	l_length = pValue->m_lDigitsLength;
	long	l_low = 0;

	while ( (l_low < l_length) && (pValue->m_lpDigits[ l_low ] == 0) )
		l_low++;

	if ( l_low >= l_length )
		return -1;

	long long	ll_places = static_cast<long long>( pValue->m_lRadix ) - l_low;
	long long	ll_noughts = -1;		// k = 1: the least floor( k v / e )
	long long	ll_most = 0;			// k = 1: the most ceil( k v / e )

	for ( long l_prime = 2, l_rest = l_base; l_rest > 1; l_prime++ )
	{
		if ( (l_rest % l_prime) != 0 )
			continue;

		long	l_exponent = 0;

		while ( (l_rest % l_prime) == 0 )
		{
			l_rest /= l_prime;
			l_exponent++;
		}

		// v_l(Y), from the last digits, while l^j stays well inside 64 bits...
		long long	ll_valuation = 0;

		for ( ; ; )
		{
			long long	ll_modulus = 1;
			long		l_power = static_cast<long>( ll_valuation + 1 );

			for ( long i = 0; i < l_power; i++ )
			{
				if ( ll_modulus > (static_cast<long long>( 1 ) << 50) )
					return -1;

				ll_modulus *= l_prime;
			}

			long	l_digits = (l_power + l_exponent - 1) / l_exponent;
			long long	ll_remainder = 0;

			for ( long i = l_low + l_digits - 1; i >= l_low; i-- )
				ll_remainder = ((ll_remainder * l_base) + ((i < l_length) ? pValue->m_lpDigits[ i ] : 0)) % ll_modulus;

			if ( ll_remainder != 0 )
				break;

			ll_valuation++;
		}

		long long	ll_floor = (static_cast<long long>( lPower ) * ll_valuation) / l_exponent;
		long long	ll_ceiling = ((static_cast<long long>( lPower ) * ll_valuation) + l_exponent - 1) / l_exponent;

		if ( (ll_noughts < 0) || (ll_floor < ll_noughts) )
			ll_noughts = ll_floor;

		if ( ll_ceiling > ll_most )
			ll_most = ll_ceiling;
	}

	long long	ll_answer = bReciprocal ? (ll_most - (static_cast<long long>( lPower ) * ll_places))
										: ((static_cast<long long>( lPower ) * ll_places) - ll_noughts);

	if ( ll_answer < 0 )
		ll_answer = 0;

	return (ll_answer > LONG_MAX) ? -1 : static_cast<long>( ll_answer );
}

static bool lzi_power_is_exactly( CLZNumber* pValue, CLZNumber* pBase, CLZNumber* pExponent )
{
	const long	kMostDenominator = 64;
	const long	kMostNumerator = 4096;

	if ( pValue->FlagIsOn( kfInfinity ) || pValue->FlagIsOn( kfQNaN ) ||
		 pBase->FlagIsOn( kfInfinity ) || pBase->FlagIsOn( kfQNaN ) ||
		 pExponent->FlagIsOn( kfInfinity ) || pExponent->FlagIsOn( kfQNaN ) )
		return false;

	long	l_p = 0;
	long	l_q = 0;

	{
		CLZNumber	exponent;

		exponent.Assign( pExponent );
		exponent.TrimTrailing();

		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = exponent.m_lRadix + 1;

		for ( long q = 1; (q <= kMostDenominator) && (l_q == 0); q++ )
		{
			CLZNumber	scaled;
			CLZNumber	factor;

			factor.SetBase( exponent.m_lBase, false );
			factor.Assign( q );
			lz_multiply( &exponent, &factor, &scaled );
			scaled.TrimTrailing();

			if ( scaled.m_lRadix == 0 )
			{
				long	l_value = 0;

				// A whole number: small enough to work with, or not worth it...
				for ( long i = scaled.m_lDigitsLength - 1; (i >= 0) && (l_value <= kMostNumerator); i-- )
					l_value = (l_value * scaled.m_lBase) + scaled.m_lpDigits[ i ];

				if ( l_value > kMostNumerator )
					break;

				l_p = (scaled.m_lSign == kNegative) ? -l_value : l_value;
				l_q = q;
			}
		}

		glWorkingPrecision = l_saved;
	}

	if ( l_q == 0 )
		return false;

	// Places first: unequal places rule an exact power out without working
	//	it -- see lzi_places_of_power. Asked only of two in one base...
	if ( pValue->m_lBase == pBase->m_lBase )
	{
		long	l_left = lzi_places_of_power( pValue, l_q, false );
		long	l_right = lzi_places_of_power( pBase, (l_p < 0) ? -l_p : l_p, l_p < 0 );

		if ( (l_left >= 0) && (l_right >= 0) && (l_left != l_right) )
			return false;
	}

	CLZNumber	left;
	CLZNumber	right;

	gExactnessProofs++;

	lzi_exact_power( pValue, l_q, &left );
	lzi_exact_power( pBase, (l_p < 0) ? -l_p : l_p, &right );

	if ( l_p < 0 )
	{
		CLZNumber	product;
		CLZNumber	one;
		long		l_saved = glWorkingPrecision;

		left.TrimTrailing();
		right.TrimTrailing();
		glWorkingPrecision = left.m_lRadix + right.m_lRadix + 1;
		lz_multiply( &left, &right, &product );
		glWorkingPrecision = l_saved;

		one.SetBase( product.m_lBase, false );
		one.Assign( 1L );

		return lz_compare( &product, &one ) == kEqualTo;
	}

	return lz_compare( &left, &right ) == kEqualTo;
}


/*
	**How much wider a division by a small value has to work: twice the
	answer's own size.** A divisor about base^-k, worked to W places,
	carries only W - k significant digits -- and the answer, about
	base^k, wants k + 36 of them, 36 places after its point. So the
	width wants the size once for the digits the divisor lacks and once
	for the digits the answer adds in front: 2k. csch(1e-12) held 32
	digits of its 48 without it, 1e16 units out of its 36th place, and
	csc(1e-12) 28.

	Handed the divisor, or the argument that predicts it: the
	hyperbolic sine is small only where its argument is, and is read
	off that; the circular sine and cosine are small away from nought
	as well, so the circular functions read the value they computed.
	Half of it is a small numerator's widening -- see lz_tangent.
	The size is the leading digit's place, exactly and in the value's
	own base, where the powers' widening has to estimate. Nothing for
	a value of one or more, and nothing for nought, which the division
	then refuses as it always did.
*/

static long lzi_small_divisor_widening( CLZNumber* pValue )
{
	long	l_top = pValue->m_lDigitsLength - 1;

	while ( (l_top >= 0) && (pValue->m_lpDigits[ l_top ] == 0) )
		l_top--;

	if ( l_top < 0 )
		return 0;

	long	l_exponent = l_top - pValue->m_lRadix;

	return (l_exponent < 0) ? (2 * -l_exponent) : 0;
}


/*
	**A sine or cosine that reads nought has not been seen, unless its
	argument is nought.** No finite argument is a zero of either but
	nought itself, the sine's; the value is only smaller than the width
	it was taken at can show, and then there is no size to widen by.
	Base ten's first take at forty places sees down to 1e-40, but base
	two's forty binary places see only to about 1e-12, so there pi/2 to
	twenty digits had a cosine of nought: tan and sec refused it as a
	division by zero and cot answered nought, where its answer is
	1.3e-21.

	Taken again the argument's own places wider -- about as near as a
	value with that many places comes to a multiple of pi/2 -- doubling
	until it shows, a few times at most. Answers the places it added;
	the caller's widening then sizes the value as it would any other.
	Nothing is paid unless the value reads nought.

	**Shows means stands clear of the width it was taken at**, by the
	guard. Until 2026-09-26 it was also capped at 2*pi's 120 places,
	since the reduction subtracts that constant and past its width the
	difference was the constant's own rounding: pi/2 written to 130
	places had a cosine of 1e-131 and read 1e-121 at any width. Pi is
	computed as wide as it is asked for now -- see
	lzi_pi_family_in_base -- so the retakes alone bound the widening. A
	value still unseen after them is left nought, and the division
	refuses, as it always did.
*/

const long	kUnseenRetakes = 4;

static long lzi_widen_until_seen( long (*fnTrig)( CLZNumber*, CLZNumber* ),
								  CLZNumber* pArgument, CLZNumber* pValue )
{
	if ( pValue->FlagIsOn( kfQNaN ) || pValue->FlagIsOn( kfInfinity ) )
		return 0;

	pValue->CheckForZero();
	pArgument->CheckForZero();

	if ( !(pValue->IsZero()) || pArgument->IsZero() )
		return 0;

	auto	Unseen = [ pValue ]() -> bool
	{
		pValue->CheckForZero();

		return pValue->IsZero() ||
			   ((lzi_small_divisor_widening( pValue ) / 2) > (glWorkingPrecision - kTrigGuardDigits));
	};

	long	l_step = pArgument->m_lRadix + kTrigGuardDigits;
	long	l_added = 0;

	for ( long l_try = 0; (l_try < kUnseenRetakes) && Unseen(); l_try++ )
	{
		long	l_this = l_step;

		glWorkingPrecision += l_this;
		l_added += l_this;
		l_step *= 2;

		fnTrig( pArgument, pValue );
	}

	if ( Unseen() )
		pValue->Assign( static_cast<long>(0) );

	return l_added;
}


/*
	Multiply or divide by a whole power of the base, exactly.

	A positive power divides, a negative one multiplies. Nothing moves
	but the radix where it can, which makes this free and exact -- the
	digits are untouched, so no rounding enters. Both roots build their
	starting estimate with it.
*/

static void lzi_scale_by_base_power( CLZNumber* value, long power )
{
	if ( power > 0 )
	{
		value->m_lRadix += power;

		if ( value->m_lRadix > value->m_lDigitsLength )
			value->PadLeft( value->m_lRadix - value->m_lDigitsLength );
	}
	else if ( power < 0 )
	{
		long	l_up = -power;

		if ( value->m_lRadix >= l_up )
		{
			value->m_lRadix -= l_up;
		}
		else
		{
			value->ShiftLeft( l_up - value->m_lRadix );
			value->m_lRadix = 0;
		}
	}
}


/*
	**Where the square root takes the AGM's core** -- lzi_root_for_agm,
	Newton on the reciprocal root with multiplications only and Karp and
	Markstein's last step -- in place of Heron's iteration, whose every
	pass is a division. Measured 2026-10-14 against the root's significant
	places, the answer's and the
	root's own with the guard apart:

	- **From the base's floor**: 5 places, so the default 36 too, where it
	  is 1.2 to 1.8 times faster in bases five to sixteen; 10 in bases
	  three and four; base two not until 5,000, short of which Heron won.
	- **But a short argument keeps Heron's iteration within its base's
	  band**, Rich's call. In a small base at middling widths Heron divides
	  a short argument cheaply, and the core took up to 1.35 times its time;
	  past the band the core's multiplications reach Karatsuba and it wins
	  again.
	  Short is holding fewer digits than seven-eighths of the root's
	  significant places: the core lost at every length to half the width
	  and won from 67 to 85% of it, by base. **Held, not significant**: two
	  written out with the width's noughts lost Heron its advantage -- the
	  core 1.05 to 1.15 times faster on it in bases three and six -- where
	  a bare 2 kept it, so it is the digits the division is handed that
	  count.

	The floor and the bands in force are the process's, which only the
	tests and measurements move -- a floor set through the hook applies alone,
	bands and all set aside; the count is of the core's runs, this
	thread's.
*/

static long gRootCoreFloor = kRootCoreFloorMeasured;
static thread_local long gRootCoreRuns = 0;

long lzi_set_root_core_floor( long lPlaces )
{
	long	l_previous = gRootCoreFloor;

	gRootCoreFloor = lPlaces;

	return l_previous;
}

long lzi_root_core_floor( long lBase )
{
	if ( gRootCoreFloor != kRootCoreFloorMeasured )
		return gRootCoreFloor;

	switch ( lBase )
	{
	case 2:		return 5000;
	case 3:
	case 4:		return 10;
	default:	return 5;
	}
}

// The band, in the root's significant places, within which a short argument
//	keeps Heron's iteration -- each drawn to take in every width at which the
//	core lost to it by more than 3%, measured at 22 widths from 50 to 3,000.
//	False for a base with none...
static bool lzi_root_short_band( long lBase, long* plFrom, long* plTo )
{
	switch ( lBase )
	{
	case 3:		*plFrom = 100;	*plTo = 2000;	return true;
	case 4:		*plFrom = 80;	*plTo = 1000;	return true;
	case 5:		*plFrom = 150;	*plTo = 700;	return true;
	case 6:		*plFrom = 64;	*plTo = 1250;	return true;
	case 8:		*plFrom = 125;	*plTo = 600;	return true;
	case 12:	*plFrom = 175;	*plTo = 250;	return true;
	case 14:	*plFrom = 175;	*plTo = 300;	return true;
	default:	return false;
	}
}

static bool lzi_root_core_wanted( const CLZNumber* pOperand, long lRootPlaces )
{
	long	l_base = pOperand->m_lBase;

	if ( lRootPlaces < lzi_root_core_floor( l_base ) )
		return false;

	if ( gRootCoreFloor != kRootCoreFloorMeasured )
		return true;

	long	l_from = 0;
	long	l_to = 0;

	if ( lzi_root_short_band( l_base, &l_from, &l_to ) && (lRootPlaces >= l_from) && (lRootPlaces < l_to) &&
		 ((8 * pOperand->m_lDigitsLength) < (7 * lRootPlaces)) )
		return false;

	return true;
}

long lzi_root_core_runs()
{
	return gRootCoreRuns;
}

long lz_root_2( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Nought roots to nought exactly -- see lzi_write_to_width...
	answer.Exact( lzi_exactly( leftOperand, 0 ) );

	// An infinity is held as zero here, so this has to come before any
	//	test of the value, which would read it as nought...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	// No real square root of a negative. Newton's method converged on
	//	whatever it liked and reported success -- -9 gave 5.311...
	if ( leftOperand->IsNegative() )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	/*
		Zero roots to zero, and has to be taken before the iteration
		rather than through it. Newton divides by the running
		estimate, which starts at zero for a zero argument, so the
		loop performed a nought-over-nought sixty-four times and came
		out with the right answer carrying a not-a-number flag. That
		was invisible while nothing propagated flags.
	*/

	leftOperand->CheckForZero();

	if ( leftOperand->IsZero() )
	{
		CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		return kNoError;
	}

	/*
		**Worked with a guard, halved by a division, and rounded once at
		the end**, since 2026-10-01. It worked at the answer's width
		with no guard, multiplied by a half made with lz_reciprocal, and
		answered whatever the last multiplication rounded to. So 867 of
		9,100 random square roots in bases two to sixteen were not
		correctly rounded at their places -- 509 by more than a unit --
		and 551 of 1,300 in base seven, where a half is .333... and so
		every pass multiplied by a value slightly off one.

		Raised before the working values are made: a CLZNumber takes its
		precision from glWorkingPrecision at construction.
	*/

	long	l_savedPrecision = glWorkingPrecision;
	long	l_guard = lzi_guard_places( kSquareRootGuardDigits, leftOperand->m_lBase );

	// The root's size, the place of its leading digit, from the
	//	operand's: floor rather than truncate, so negatives halve the
	//	same way...
	long	l_exponent = lzi_leading_exponent( leftOperand );
	long	l_half = (l_exponent >= 0) ? (l_exponent / 2) : -(((-l_exponent) + 1) / 2);

	/*
		**A root under a place past the answer's width is a rounded
		nought**, and is answered before anything is iterated. The root
		is under base^(l_half + 1), so from l_half + 1 <= -(W + 1) it is
		under half a place at W, whatever the base. The guard writes the
		nought out to its width. Until 2026-10-02 the iteration took
		these, rounded its iterate to nought at the working width and
		divided by it on the next pass: the root of anything below
		about 1e-83 answered 0 flagged infinite, and arcsin(1 - 1e-100),
		which takes the root of 1 - x^2, answered 0 flagged
		not-a-number in place of pi/2.
	*/

	if ( l_half < -(l_savedPrecision + 1) )
	{
		CLZNumber* p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		return kNoError;
	}

	/*
		**Worked in a band, at widths that double.** The operand is
		scaled by an even power of the base into [1, base^2), so its
		root lies in [1, base) and its places are its significant
		digits less one; the root is scaled back at the end. Scaling
		moves the radix and is exact. The answer's W places are
		W + l_half places of the band's root, and the guard rides on
		top. Newton climbs to that width a pass at each of
		lzi_root_widths, settling at the first and at the last, which
		is where the iteration always settled.

		Made at the full width: a CLZNumber takes its precision from
		glWorkingPrecision at construction, and each pass sets it.
	*/

	long	l_target = l_savedPrecision + l_half + l_guard;

	if ( l_target < 1 )
		l_target = 1;

	glWorkingPrecision = l_target;

	auto two = new CLZNumber;
	auto band = new CLZNumber;
	auto x_current = new CLZNumber;
	auto x_next = new CLZNumber;
    CLZNumber* swap = nullptr;	// Don't assign new!

	two->SetBase( leftOperand->m_lBase, false );
	two->Assign( static_cast<long>(2) );

	band->Assign( leftOperand );
	band->SetPrecision( l_target );
	lzi_scale_by_base_power( band, 2 * l_half );

	/*
		Start near the answer: the band's operand itself, which lies
		within a factor of the base of its root, whichever side.

		The estimate used to be half the operand, which is a very long
		way from the root of anything large, and Newton only halves
		while it is far above. That was not merely slow, it was wrong.
		Sixty-three halvings from 5e39 land on 5.4e20, and that is
		precisely what sqrt(1e40) used to return in place of 1e20 -- no
		flag, no code, a plausible number of the wrong magnitude.
	*/

	x_current->Assign( band );

	/*
		The AGM's core to the band's width, which works two digits past it
		and lands well inside the window below, or Heron's iteration a pass
		at each width -- see lzi_root_core_floor...
	*/

	long	la_widths[ kRootMostWidths ];
	long	l_widths = 0;

	if ( lzi_root_core_wanted( leftOperand, l_savedPrecision + l_half ) )
	{
		SLZRootSeed	seed;

		gRootCoreRuns++;
		lzi_root_for_agm( band, x_current, l_target, &seed, 0 );
		glWorkingPrecision = l_target;
	}
	else
		l_widths = lzi_root_widths( l_target, la_widths );

	for ( long l_stage = 0; l_stage < l_widths; l_stage++ )
	{
		long	l_width = la_widths[ l_stage ];
		bool	b_settle = (l_stage == 0);

		// Each width is set on the working values, as a division takes
		//	its width from its target...
		glWorkingPrecision = l_width;
		x_current->SetPrecision( l_width );
		x_next->SetPrecision( l_width );

		/*
			Settled at the first width, from the estimate, and a pass at
			each after it, the last included: that pass lands within a
			unit of the root, and the guard decides the rest -- see
			lzi_rounding_in_doubt. The cap stays as a guard against an
			estimate that oscillates between two adjacent values rather
			than settling on one.
		*/

		for ( long lCycle = 1; lCycle < kRootMaxCycles; lCycle++ )
		{
			// Calculate this cycle: (x/y + y) / 2, brought back to the
			//	working width before it is compared...
			lz_divide( band, x_current, x_next );
			lz_add( x_next, x_current );
			lz_divide( x_next, two );
			x_next->Round();

			bool	b_settled = ( x_next->Compare( x_current ) == kEqualTo );

			// Set up the next cycle...
			swap = x_current;
			x_current = x_next;
			x_next = swap;

			if ( b_settled || !b_settle )
				break;
		}
	}

	// Out of the band...
	lzi_scale_by_base_power( x_current, -l_half );

	/*
		**Rounded at the answer's places by the guard, or decided exactly
		when the guard is too near the half to say** -- see
		lzi_rounding_in_doubt and lzi_root_decide. The root is W places
		and the guard here, out of the band.
	*/

	if ( lzi_rounding_in_doubt( x_current, l_savedPrecision, l_savedPrecision + l_guard, kRootDoubtUnits ) )
		lzi_root_decide( x_current, leftOperand, 2, l_savedPrecision );

	/*
		The guard off, once, at the answer's places -- **W places for a
		small answer too**, where the cube root keeps W significant
		digits by the rule for a quotient. Rich's call, 2026-10-01, and
		not an oversight: root_n, power and exp_n taken at a half answer
		the same value at W places, and the rule keeps digits a division
		has already made, which this iteration does not make for a
		small answer unless it is widened by the answer's noughts.
	*/

	glWorkingPrecision = l_savedPrecision;
	x_current->SetPrecision( l_savedPrecision );
	x_current->Round();

	// A short root is checked against the operand before it is called
	//	exact, and spelled without its noughts if it is: the iterates
	//	pass through inexact values on the way to an exact one -- see
	//	lzi_write_to_width...
	if ( lzi_short_of_width( x_current ) && lzi_power_is_exactly( leftOperand, x_current, two ) )
	{
		x_current->TrimTrailing();
		answer.Exact();
	}

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( x_current );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( x_current );

	// Clean up...
    delete two;
    delete band;
    delete x_current;
    delete x_next;

	return kNoError;
}


// The cube root's Newton passes this thread has made, for the tests...
static thread_local long gRoot3Passes = 0;

long lzi_root_3_passes()
{
	return gRoot3Passes;
}

/*
	**Past a floor, the cube root takes its stages after the first by its
	reciprocal** -- since 2026-10-17, the remainder under performance item 8
	of the review of 2026-10-15. Each pass of the root's own Newton divides
	the operand by the iterate's square, and one division costs two products
	at 100 digits, six at 1,000 and ten at 3,000 in base ten: 71% of a cube
	root at 1,000 digits, 89% at 3,000. Newton on r = x^(-1/3) multiplies
	instead,

		r <- r + r (1 - x r^3) / 3

	a square and three products a pass and a division by three, converging
	as the root's own does, one pass to each doubled width. The first
	stage settles the root at its narrow width as before and one division
	there seeds r; after the last, the root is x r^2.

	**Two things it needs that the root's own pass does not.** r lies
	between 1/b and one, its cube as low as 1/b^3, so a product counting
	places loses up to two of its digits to leading noughts:
	kRoot3InverseGuard places more on every pass and on x r^2 cover that and
	the rounding -- four, where two left 55 of 1,800 near ties wrong and
	three 2 of 7,200. And every number in it is in the operand's base: a
	reciprocal of a default number came back in base ten, whose places hold
	fewer digits than twelve's, and 752 of 3,000 cube roots went wrong in
	their last places. The guard, the doubt test and the exact decision
	after are the root's own, unchanged.

	**The floors are each base's**, measured as the sine kernel's were:
	the cube root forced each way at 36 to 10,000 digits, a short operand,
	one as wide as the width and one three times as wide, two passes, and
	the floor the least width from which the reciprocal is never more than
	5% slower. At 10,000 digits it is 2.8 to 4.8 times faster on a short
	operand, 1.7 to 2.8 on a full-width one, 1.4 to 2.1 on a long one. Base
	two's floor is the highest, then three's and ten's; ten's old road was
	the quick one, the three buffers' fault never having reached it.
*/

constexpr long	kRoot3InverseGuard = 4;

static long gRoot3InverseFloor = kRoot3InverseFloorMeasured;
static long gRoot3InverseGuard = kRoot3InverseGuard;

long lzi_set_root_3_inverse_guard( long lPlaces )
{
	long	l_previous = gRoot3InverseGuard;

	gRoot3InverseGuard = (lPlaces < 0) ? kRoot3InverseGuard : lPlaces;

	return l_previous;
}
static thread_local long gRoot3InversePasses = 0;

long lzi_set_root_3_inverse_floor( long lPlaces )
{
	long	l_previous = gRoot3InverseFloor;

	gRoot3InverseFloor = lPlaces;

	return l_previous;
}

long lzi_root_3_inverse_floor( long lBase )
{
	if ( gRoot3InverseFloor != kRoot3InverseFloorMeasured )
		return gRoot3InverseFloor;

	// Each base's, in places: 1,500 digits in base two, 600 in three and
	//	ten, 300 in four, 200 in five, 100 in six, 150 in seven to fifteen,
	//	and sixteen from its first width with a second stage...
	static const long	kFloors[ 17 ] =
	{
		0, 0, 4983, 1258, 499, 287, 129, 178, 167, 158, 600, 145, 139, 135, 131, 128, 30
	};

	return ((lBase >= 2) && (lBase <= 16)) ? kFloors[ lBase ] : LONG_MAX;
}

long lzi_root_3_inverse_passes()
{
	return gRoot3InversePasses;
}

// One pass of Newton on r = x^(-1/3) at lWidth places, multiplying where the
//	cube root's own pass divides: r <- r + r (1 - x r^3) / 3. pSeed, on the
//	first, is the cube root the narrower stages settled at lSeedWidth, whose
//	reciprocal starts r. glWorkingPrecision is lWidth already...
static long lzi_root_3_inverse_pass( CLZNumber* pOperand, CLZNumber* pInverse, long lWidth,
									 CLZNumber* pSeed, long lSeedWidth )
{
	long	lResult = kNoError;
	long	l_saved = glWorkingPrecision;

	lWidth += gRoot3InverseGuard;
	glWorkingPrecision = lWidth;

	long	l_base = pOperand->m_lBase;

	if ( pSeed != nullptr )
	{
		CLZNumber	one;

		one.SetBase( l_base, false );
		one.Assign( 1L );
		pInverse->SetBase( l_base, false );
		pInverse->SetPrecision( lSeedWidth );
		lResult = lzi_first_failure( lResult, lz_divide( &one, pSeed, pInverse ) );
	}

	pInverse->SetPrecision( lWidth );

	CLZNumber	square;
	CLZNumber	cube;
	CLZNumber	product;
	CLZNumber	error;
	CLZNumber	step;
	CLZNumber	one;
	CLZNumber	three;

	one.SetBase( l_base, false );
	one.Assign( 1L );
	three.SetBase( l_base, false );
	three.Assign( 3L );

	lResult = lzi_first_failure( lResult, lz_square( pInverse, &square ) );
	lResult = lzi_first_failure( lResult, lz_multiply( &square, pInverse, &cube ) );
	lResult = lzi_first_failure( lResult, lz_multiply( pOperand, &cube, &product ) );
	lResult = lzi_first_failure( lResult, lz_subtract( &one, &product, &error ) );
	lResult = lzi_first_failure( lResult, lz_multiply( pInverse, &error, &step ) );
	lResult = lzi_first_failure( lResult, lz_divide( &step, &three ) );
	lResult = lzi_first_failure( lResult, lz_add( pInverse, &step ) );

	pInverse->Round();
	gRoot3InversePasses++;
	glWorkingPrecision = l_saved;

	return lResult;
}

long lz_root_3( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Nought roots to nought exactly -- see lzi_write_to_width...
	answer.Exact( lzi_exactly( leftOperand, 0 ) );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// Infinite whichever sign went in, so the missing sign costs
	//	nothing and the infinity carries...
	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryStaysInfinite ) )
		return kNoError;

	long	lResult = kNoError;

	/*
		Zero roots to zero, and has to be taken before the iteration
		rather than through it: the estimate would be nought and the
		division by its square a nought over nought.

		It used to go the other way and through lz_log_e, which
		refuses zero, so the cube root of nought came back correct but
		carrying the not-a-number flag that refusal set.
	*/

	leftOperand->CheckForZero();

	if ( leftOperand->IsZero() )
	{
		CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		return kNoError;
	}

	/*
		The band reduction below takes the argument's magnitude out
		before iterating and puts it back afterwards, and **the places
		have to be taken out and put back with it**.

		Precision counts places after the radix mark. In the band the
		iterate sits between one and the base, where the working width
		is 36 places and about 37 significant digits; scaling it back
		up by base^block moves the radix, so those 37 digits become 37
		digits of a much larger number -- 36 - block places, not 36.
		The cube root of 1e40 came out with 37 significant digits where
		36 places at 1e13 has room for 50, and the shortfall was
		exactly the block count: 1 at 1e8, 5 at 1e20, 13 at 1e40.

		So the iteration runs block places wider and the answer is
		rounded back at the end. The scaling itself then lands it on 36
		places, the width having been carried through the band with it.

		lz_root_2 has no such trouble and the difference is worth
		knowing: it scales only its starting estimate and divides the
		unscaled operand, so its iterate is at full magnitude
		throughout and never loses the places. The cube root cannot do
		that -- squaring a small iterate is what drops it under the
		floor, which is the whole reason for the band.

		Raised before the working values are made, a CLZNumber taking
		its precision from glWorkingPrecision at construction and not
		looking again.
	*/

	long	l_exponent = lzi_leading_exponent( leftOperand );

	// Floor rather than truncate, so a negative exponent blocks the
	//	same way a positive one does...
	long	l_block = (l_exponent >= 0) ? (l_exponent / 3)
										: -(((-l_exponent) + 2) / 3);

	long	l_savedPrecision = glWorkingPrecision;

	if ( l_block > 0 )
		glWorkingPrecision += l_block;

	// ...and the iteration's guard, which the answer's rounding below takes off...
	long	l_guard = lzi_guard_places( kCubeRootGuardDigits, leftOperand->m_lBase );

	glWorkingPrecision += l_guard;

	CLZNumber	magnitude;
	CLZNumber	two;
	CLZNumber	three;
	CLZNumber	x_current;
	CLZNumber	x_next;
	CLZNumber	x_previous;
	CLZNumber	squared;
	CLZNumber	quotient;

	bool	b_negateResult = leftOperand->IsNegative();

	two.Assign( 2L );
	three.Assign( 3L );

	/*
		Newton on y^3 = x, which is

			y <- (2y + x/y^2) / 3

		rather than exp(ln(x)/3), which is what this was.

		Two things were wrong with going through the logarithm. It has
		no answer for a negative, so the cube root of one -- which is
		perfectly well defined, -8 rooting to -2 -- came back as nought
		carrying a not-a-number flag, and with kNoError, because the
		refusal was never looked at. And the round trip through two
		series cost digits that an iteration does not: the cube root of
		27 came back as 2.999999999999999999999999999999999651.

		The cube root is odd, so the sign comes off the front and goes
		back on the answer.
	*/

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	/*
		Bring the argument into [1, base^3) before iterating, and put
		the answer back afterwards.

		Precision counts places after the radix mark, so squaring a
		small value drops it under the floor: the cube root of 1e-40 is
		about 4.6e-14, whose square is 2.2e-27, and 36 places leaves
		that nine significant digits. The iteration then divides by it
		and inherits the loss -- eleven correct digits where the answer
		had room for twenty-two.

		Scaling by a power of base^3 costs nothing and is exact, being
		a move of the radix, and it commutes with the root:
		cbrt(x * base^3k) is cbrt(x) * base^k. Inside the band the
		square of the estimate stays near one and no places are lost.
	*/

	x_current.Assign( &magnitude );
	lzi_scale_by_base_power( &magnitude, 3 * l_block );
	lzi_scale_by_base_power( &x_current, 3 * l_block );

	// ...and a starting estimate within the band, the scaled value
	//	having an exponent of 0, 1 or 2...
	{
		long	l_scaledExponent = l_exponent - (3 * l_block);
		long	l_twice = 2 * l_scaledExponent;

		lzi_scale_by_base_power( &x_current, l_twice / 3 );
	}

	/*
		**At widths that double**, since 2026-10-02 -- see
		lzi_root_widths. The band's root lies in [1, base), so its
		places are its significant digits less one, and the full width
		is the one raised above, block and guard. Settled at the first
		width, and a pass at each after it, the last included.
	*/

	long	l_target = glWorkingPrecision;

	long	la_widths[ kRootMostWidths ];
	long	l_widths = lzi_root_widths( l_target, la_widths );

	// The iterate, the next, and the one before, turned by pointer so
	//	that no pass copies a value...
	CLZNumber*	p_current = &x_current;
	CLZNumber*	p_next = &x_next;
	CLZNumber*	p_previous = &x_previous;

	/*
		**The operand narrowed to each width, two places over** -- since
		2026-10-15, performance item 8 of that day's review. It was divided
		at its own width in every pass, however narrow the stage, and a
		division works every quotient digit to its dividend's places before
		dropping them: a cube root of an operand far wider than the answer
		paid for all of it every pass. Newton corrects what the narrowing
		leaves, and the exact decision at the end reads the operand
		itself...
	*/
	CLZNumber	narrowed;

	// Past its floor, the stages after the first by the reciprocal...
	CLZNumber	inverse;
	bool		b_inverse = (l_widths > 1) && (l_target >= lzi_root_3_inverse_floor( leftOperand->m_lBase ));

	for ( long l_stage = 0; l_stage < l_widths; l_stage++ )
	{
		long	l_width = la_widths[ l_stage ];
		bool	b_settle = (l_stage == 0);

		/*
			Each width is set on the working values, as a division and a
			product take theirs from their targets -- **all three of them**,
			since 2026-10-17. The one before was left at the full width it
			was made with, and the turn of the three hands it on as the next:
			the product and the division by three then answered at the full
			width, and the rounding after kept it, so at 50 places one pass
			in three wrote an iterate 11,842 digits long and the pass after
			it squared that. Such an iterate never equals a narrow one, so
			the stage took the three for a cycle it could not see and ran to
			kRootMaxCycles: 71 passes for the cube root of 2 at 10,000 digits
			in bases five and seven, 13.0 ms against 4.8 set right. Answers
			were right, the later stages correcting them; 47 of 72 widths and
			bases measured were more than 5% slower, base ten untouched...
		*/
		glWorkingPrecision = l_width;
		p_current->SetPrecision( l_width );
		p_next->SetPrecision( l_width );
		p_previous->SetPrecision( l_width );
		squared.SetPrecision( l_width );
		quotient.SetPrecision( l_width );

		CLZNumber*	p_dividend = &magnitude;

		if ( magnitude.m_lRadix > (l_width + 2) )
		{
			narrowed.Assign( &magnitude );
			narrowed.Round( kDefaultRoundingMode, l_width + 2 );
			p_dividend = &narrowed;
		}

		if ( b_inverse && (l_stage > 0) )
		{
			lResult = lzi_first_failure( lResult, lzi_root_3_inverse_pass( p_dividend, &inverse, l_width,
												(l_stage == 1) ? p_current : nullptr, la_widths[ 0 ] ) );
			continue;
		}

		for ( long lCycle = 1; lCycle < kRootMaxCycles; lCycle++ )
		{
			gRoot3Passes++;

			lResult = lzi_first_failure( lResult, lz_square( p_current, &squared ) );
			lResult = lzi_first_failure( lResult, lz_divide( p_dividend, &squared, &quotient ) );
			lResult = lzi_first_failure( lResult, lz_multiply( p_current, &two, p_next ) );
			lResult = lzi_first_failure( lResult, lz_add( p_next, &quotient ) );
			lResult = lzi_first_failure( lResult, lz_divide( p_next, &three ) );

			/*
				Bring the iterate back to the working precision before
				comparing it. Division deliberately answers wider than it
				was asked -- see the sharp edge on that -- and the last
				step here is a division, so without this the estimate
				carries places below the precision, two passes never look
				equal, and the loop runs to its cap and returns an answer
				with noise on the end. The square root does the same since
				2026-10-01, when its half became a division.
			*/

			p_next->Round();

			bool	b_settled = ( p_next->Compare( p_current ) == kEqualTo );

			/*
				**Two passes that swap one value for its neighbour have
				settled too.** The iterate can land a unit either side
				of the root at this width and alternate between the two
				for ever, and until 2026-10-02 it ran to the cap: 63
				passes in place of eleven, on 301 of 13,600 random cube
				roots, cbrt(1.2345...) at 200 places among them. Either
				value rounds alike past the guard, so the smaller is
				kept, whichever pass the cycle is noticed on.
			*/

			bool	b_cycle = !b_settled && (lCycle > 1) && ( p_next->Compare( p_previous ) == kEqualTo );

			if ( b_cycle && ( p_next->Compare( p_current ) == kLessThan ) )
			{
				CLZNumber*	p_turn = p_current;

				p_current = p_next;
				p_next = p_turn;
			}
			else if ( !b_cycle )
			{
				CLZNumber*	p_turn = p_previous;

				p_previous = p_current;
				p_current = p_next;
				p_next = p_turn;
			}

			if ( b_settled || b_cycle || !b_settle )
				break;
		}
	}

	// The reciprocal's root: x r^2, at the full width...
	if ( b_inverse )
	{
		glWorkingPrecision = l_target + gRoot3InverseGuard;

		CLZNumber	square;
		CLZNumber*	p_operand = &magnitude;

		if ( magnitude.m_lRadix > (l_target + 2) )
		{
			narrowed.Assign( &magnitude );
			narrowed.Round( kDefaultRoundingMode, l_target + 2 );
			p_operand = &narrowed;
		}

		square.SetPrecision( l_target + gRoot3InverseGuard );
		x_current.SetPrecision( l_target + gRoot3InverseGuard );
		lResult = lzi_first_failure( lResult, lz_square( &inverse, &square ) );
		lResult = lzi_first_failure( lResult, lz_multiply( p_operand, &square, &x_current ) );
		glWorkingPrecision = l_target;
		x_current.SetPrecision( l_target );
		x_current.Round();
		p_current = &x_current;
	}

	if ( p_current != &x_current )
		x_current.Assign( p_current );

	// Put back the power of the base taken out above...
	lzi_scale_by_base_power( &x_current, -l_block );
	lzi_scale_by_base_power( &magnitude, -3 * l_block );

	glWorkingPrecision = l_savedPrecision;

	/*
		**Rounded at the answer's places by the guard, or decided exactly
		when the guard is too near the half to say** -- see
		lzi_rounding_in_doubt and lzi_root_decide. The places are the
		rule for a quotient's, as lzi_round_quotient_answer will take
		them: W, and for an answer under a tenth its noughts as well.
		The magnitude is the operand's again, scaled back. The guard
		ends where the iteration's width does, out of the band: an
		iterate that came out short, as an exact root does, has noughts
		to there and is clear of the half. Read to its own last digit
		instead, it had no guard at all, and 62 of 13,600 random cube
		roots were decided for nothing.
	*/

	{
		long	l_places = l_savedPrecision;
		long	l_lead = lzi_leading_exponent( &x_current );

		if ( l_lead < -1 )
			l_places += -l_lead - 1;

		if ( lzi_rounding_in_doubt( &x_current, l_places, l_target - l_block, kRootDoubtUnits ) )
			lzi_root_decide( &x_current, &magnitude, 3, l_places );
	}

	/*
		The extra places never leave the function, and the answer is
		rounded to the rule for a quotient -- see
		lzi_round_quotient_answer.

		Rounding to 36 places unconditionally cost the small arguments
		badly: the iteration's last step is a division, which answers
		wider than it was asked, and for an argument under one those
		spare places are real accuracy the function has always
		returned. Cutting them took the cube root of 1e-40 from 38
		significant digits to 24. So this rounded only where a block
		had been taken out, until 2026-09-25, and handed the division's
		extra places back everywhere else -- a 37th on the cube root of
		0.3. The rule keeps a small answer's 36 significant digits and
		takes the tail off the rest.
	*/

	lzi_round_quotient_answer( &x_current );

	if ( b_negateResult )
		lResult = lzi_first_failure( lResult, lz_negate( &x_current ) );

	// A short root is checked against the operand before it is called
	//	exact, and spelled without its noughts if it is -- see
	//	lzi_write_to_width...
	if ( lzi_short_of_width( &x_current ) && lzi_power_is_exactly( leftOperand, &x_current, &three ) )
	{
		x_current.TrimTrailing();
		answer.Exact();
	}

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &x_current );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &x_current );

	return lResult;
}


/*
	**Whether one over pValue is a whole number, and if so which.** A
	value with r places is m over b^r, so its reciprocal is b^r over m,
	whole exactly when m divides b^r: the whole quotient, multiplied back
	exactly, has to come to one. Both are worked at the value's own
	places, so nothing rounds. A whole value has none but one and minus
	one; a half, a fifth, 1E-6 in base ten, and a third in base twelve
	all do.
*/

static bool lzi_reciprocal_is_whole( CLZNumber* pValue, CLZNumber* pReciprocal )
{
	CLZNumber	value;
	CLZNumber	one;
	CLZNumber	product;

	value.Assign( pValue );
	value.TrimTrailing();

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = value.m_lRadix + 2;

	one.SetBase( value.m_lBase, false );
	one.Assign( 1L );

	lz_divide_integer( &one, &value, pReciprocal );
	lz_multiply( pReciprocal, &value, &product );

	glWorkingPrecision = l_saved;

	return (lz_compare( &product, &one ) == kEqualTo);
}


static long lzi_root_n_work( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	/*
		The nth root of an infinity is infinite for a positive n and
		nought for a negative one. Two answers, and an infinity here
		carries no sign to choose between them, so a not-a-number --
		and before any test of the value, which would read the
		infinity as nought.
	*/

	if ( lzi_propagate_unary( leftOperand, targetNumber, kUnaryNeedsTheSign ) )
		return kNoError;

	long	lResult = kNoError;

	/*
		This function used to be three lines with no guards at all:
		the logarithm, a division and the exponential, every return
		code dropped and kNoError returned whatever happened. A
		negative operand went into lz_log_e, which refuses it, and the
		refusal's nought came back out as the answer -- so the cube
		root of -8 was 0, reported as success, while lz_root_3 beside
		it gave -2. That was the same fault lz_root_3 had been fixed
		for, still living here because nothing tested it.
	*/

	CLZNumber	order;
	CLZNumber	two;
	CLZNumber	remainder;

	order.Assign( rightOperand );
	order.CheckForZero();

	// There is no root of order nought...
	if ( order.IsZero() )
		return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

	leftOperand->CheckForZero();

	/*
		Nought roots to nought for a positive order, and has to be
		taken before lz_log_e, which refuses it. For a negative order
		it is one over nought, an infinity that cannot be written.
	*/

	if ( leftOperand->IsZero() )
	{
		if ( order.IsNegative() )
			return lzi_domain_failure( leftOperand, targetNumber,
									   kfInfinity, kUndefinedResult );

		CLZNumber*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;

		if ( !(p_target->IsConstant()) )
			p_target->Assign( static_cast<long>(0) );

		return kNoError;
	}

	/*
		A negative has a real root only when the order is an odd
		integer -- the cube root of -8 is -2, the square root of -8 is
		nothing real. Take it on the magnitude and put the sign back,
		the way lz_root_3 does.
	*/

	bool	b_negateResult = false;

	if ( leftOperand->IsNegative() )
	{
		CLZNumber	orderMagnitude;

		lz_abs( &order, &orderMagnitude );

		two.Assign( 2L );

		lz_modulo( &orderMagnitude, &two, &remainder );
		remainder.CheckForZero();

		if ( !(order.IsInteger()) || remainder.IsZero() )
			return lzi_domain_failure( leftOperand, targetNumber, kfQNaN, kFailed );

		b_negateResult = true;
	}

	CLZNumber	magnitude;
	CLZNumber	result;

	lResult = lzi_first_failure( lResult, lz_abs( leftOperand, &magnitude ) );

	/*
		**An order whose reciprocal is whole is a power, and takes no
		logarithm.** The root of order 1E-6 is the millionth power, and
		the logarithm below is worked at the answer's full width, which
		for 1.2345 is 91,552 digits: 66 seconds, where the squaring takes
		under one. The sign was settled above, by the same rule as
		before -- a negative takes a real root only of an odd whole
		order, and one and minus one are the only whole orders with a
		whole reciprocal.
	*/

	CLZNumber	power;

	if ( lzi_reciprocal_is_whole( &order, &power ) )
		lResult = lzi_first_failure( lResult, lzi_raise_n( &magnitude, &power, &result ) );
	else
	{
		/*
			**A small order divides the logarithm by a small number, so the
			two are worked wider by the places one over the order has.** The
			logarithm is good to its last place, absolutely; divided by
			1.5E-12 that error is twelve places larger, and e to the
			quotient carries it as relative error into every digit of the
			answer. The width allowed for the answer's own size and nothing
			more, so near one, where a tiny order still gives an answer of
			ordinary size, a third of such roots were wrong -- 115 of 300 in
			HEAD, by up to 10^17 units. The same fault as lz_log_n's
			to a base near one, and the same remedy: only where a whole
			place is lost, so an ordinary order works as it always did.
		*/

		long	l_saved = glWorkingPrecision;
		double	d_lnOrder = 0.0;

		if ( lzi_ln_estimate( &order, &d_lnOrder ) )
		{
			double	d_lost = -d_lnOrder / std::log( static_cast<double>( order.m_lBase ) );

			if ( d_lost >= 1.0 )
				glWorkingPrecision += static_cast<long>( std::ceil( d_lost ) );
		}

		lResult = lzi_first_failure( lResult, lz_log_e( &magnitude, &result ) );
		lResult = lzi_first_failure( lResult, lz_divide( &result, &order ) );

		glWorkingPrecision = l_saved;

		lResult = lzi_first_failure( lResult, lz_exp_e( &result ) );
	}

	if ( b_negateResult )
		lResult = lzi_first_failure( lResult, lz_negate( &result ) );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( &result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( &result );

	return lResult;
}


long lz_root_n( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	/*
		**A square or cube root goes to root_2 or root_3** -- since
		2026-10-15, performance item 4 of that day's review: the logarithm
		and the exponential took 166 microseconds against 6.9 at 36 places,
		about 143 times at 10,000. Only where the two spell alike: root_2
		writes W places as this does, but root_3 writes a small answer to
		W significant digits, the rule for a quotient, so only a cube root
		of a base-th or more goes there -- one of x from base^-3 up. An
		exact root comes back short, as everywhere else; this wrote 8's
		cube root as 2.000... to the width.
	*/
	if ( !(leftOperand->FlagIsOn( kfInfinity ) || leftOperand->FlagIsOn( kfQNaN )) )
	{
		if ( lzi_exactly( rightOperand, 2 ) && !(leftOperand->IsNegative()) )
			return lz_root_2( leftOperand, targetNumber );

		leftOperand->CheckForZero();

		if ( lzi_exactly( rightOperand, 3 ) && !(leftOperand->IsZero()) && (lzi_leading_exponent( leftOperand ) >= -3) )
			return lz_root_3( leftOperand, targetNumber );
	}

	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_root_n, rightOperand, kRetakeDoubt );

	// ...and a root of whole order decided exactly in doubt -- see
	//	lzi_whole_order_root_decide...
	answer.DecideTies( lzi_root_n_decide );

	// The nth root of x has about ln|x| / n of logarithm...
	double	d_ln = 0.0;
	double	d_n = 0.0;
	long	l_magnitude = -1;

	if ( lzi_ln_estimate( leftOperand, &d_ln ) && lzi_value_or_infinity( rightOperand, &d_n ) &&
		 (d_n != 0.0) && !(leftOperand->IsNegative() && !(rightOperand->IsInteger())) )
	{
		long	l_beyond = kNoError;

		// Refused past the limit, nought below the width -- see
		//	lzi_answer_size. A tiny index is a huge power. A negative has
		//	no real root but of an odd whole order, and that refusal comes
		//	first...
		if ( lzi_answer_beyond( lzi_answer_size( d_ln / d_n, leftOperand->m_lBase ),
								leftOperand, targetNumber, &l_beyond ) )
			return l_beyond;

		l_magnitude = lzi_places_for( d_ln / d_n, leftOperand->m_lBase );
	}

	// Worked wide enough for the answer's own size -- see
	//	lzi_widened_by_answer...
	return lzi_widened_by_answer( lzi_root_n_work, leftOperand, rightOperand, targetNumber, l_magnitude );
}


long lz_round( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( leftOperand );
			targetNumber->Round();
		}
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Round();

	return kNoError;
}


long lz_round_integral( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( leftOperand );
			targetNumber->Round( RD_digit_zero );
		}
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Round( RD_digit_zero );

	return kNoError;
}


long lz_round_last_digit( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Exact by definition, so its spelling stands -- see lzi_write_to_width...
	answer.Exact();

	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( leftOperand );
			targetNumber->Round( RD_digit_last );
		}
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Round( RD_digit_last );

	return kNoError;
}


long lz_same_quantum( CLZNumber* leftOperand, CLZNumber* rightOperand )
{
	if ( lzi_is_unordered( leftOperand, rightOperand ) )
		return kUnordered;

	long lResult = kNoError;

	if ( leftOperand->m_lRadix == rightOperand->m_lRadix )
		lResult = kEqual;
	else
		lResult = kNotEqual;

	return lResult;
}


long lz_secant( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_secant, kRetakeDoubtWide );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// In another unit, the angle reduced exactly first -- see
	//	lzi_circular_in_unit...
	if ( e_unit != AU_radians )
		return lzi_circular_in_unit( kAngleSecant, leftOperand, targetNumber, answer, e_unit );

	long lResult = kNoError;

	// Guarded across the division -- see kTrigQuotientGuardDigits, and the
	//	note in lz_cotangent...
	long l_savedPrecision = glWorkingPrecision;

	long	l_guard = lzi_guard_places( kTrigQuotientGuardDigits, leftOperand->m_lBase );

	glWorkingPrecision += l_guard;

	auto result = new CLZNumber;
	auto one = new CLZNumber;
	auto cosine_result = new CLZNumber;

	// Calculate...
	// In the operand's base -- see CLZAnswerBase...
	one->SetBase( leftOperand->m_lBase );
	one->Assign( 1L );

	/*
		**Wider by twice the answer's size where the cosine is small** --
		see lzi_small_divisor_widening, and lz_cosecant, whose sine this
		is the mirror of. The cosine is small near the odd multiples of
		pi/2, which no argument's own size predicts, so it is taken, and
		taken again wider when its size would spend more than half the
		guard -- see lz_cosecant. Until 2026-09-26 it was not taken
		again at all: sec(1.5707963268) was 7e17 units out of its 36th
		place.
	*/

	lz_cosine( leftOperand, cosine_result );

	// A cosine too small to show at this width is taken wider until it
	//	does -- see lzi_widen_until_seen...
	long	l_widened = lzi_widen_until_seen( lz_cosine, leftOperand, cosine_result );
	long	l_wanted = lzi_small_divisor_widening( cosine_result );

	if ( (l_wanted - l_widened) > (l_guard / 2) )
	{
		glWorkingPrecision += l_wanted - l_widened;

		lz_cosine( leftOperand, cosine_result );
	}

	// lz_divide takes its width from its target, made before any of
	//	this...
	result->SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( one, cosine_result, result ) );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete one;
    delete cosine_result;

	return lResult;
}


/*
	**The circular functions' fold, out of lz_sine** -- shared since
	2026-10-16 with lzi_sine_and_cosine, which wants the cosine of the same
	argument from the same fold: into [0, pi/4], answering what the sine
	needs -- whether to negate it and whether the reduced argument's cosine
	is the sine -- and what the cosine needs besides, whether to negate it.
	It works at the working precision as it stands, the caller having raised
	it, every working value made after.
*/

static long lzi_trig_fold( CLZNumber* pArgument, CLZNumber* pReduced,
						   bool* pbNegateSine, bool* pbNegateCosine, bool* pbSwapped )
{
	long		lResult = kNoError;
	CLZNumber*	reduced = pReduced;
	CLZNumber	scratch_value;
	CLZNumber*	scratch = &scratch_value;

	CLZNumber	two_pi_value( kTwoPi );
	CLZNumber	pi_value( kPi );
	CLZNumber	pi_over_two_value( kPiOverTwo );
	CLZNumber	pi_over_four_value( kPiOverFour );

	CLZNumber*	two_pi = &two_pi_value;
	CLZNumber*	pi = &pi_value;
	CLZNumber*	pi_over_two = &pi_over_two_value;
	CLZNumber*	pi_over_four = &pi_over_four_value;

	bool b_negateResult = false;
	bool b_useCosineSeries = false;
	bool b_negateCosine = false;

	/*
		The series converges quickly only for small arguments: its
		error is the first term dropped, and that grows with the
		argument. Folding the argument down into [0, pi/4] before
		summing keeps the dropped term below 1e-40 no matter how
		large the argument started out.
	*/

	reduced->Assign( pArgument );

	// sin(-x) = -sin(x)...
	if ( reduced->IsNegative() )
	{
		lz_negate( reduced );
		b_negateResult = !b_negateResult;
	}

	/*
		Narrow 2*pi to the width this reduction needs, before the
		modulo divides by it.

		kTwoPi carries 120 places, and unlike e in lz_raise_e the right
		width here is not the working precision. The modulo multiplies
		2*pi by a quotient the size of the argument, so **every digit
		of magnitude the argument carries costs a digit of the
		constant** -- that is the whole reason it is 120 wide. Working
		precision plus that magnitude is exactly enough.

		Which also makes the wide constant free for the arguments that
		do not need it. Left at 120 the whole circular family paid
		1.35 to 1.96 times its old cost, worst at sin(0.5), which has
		no magnitude at all and was buying places the modulo then threw
		away. GetOrderOfMagnitude is 0 below one, so those arguments
		get the constant exactly as narrow as it always was.

		**Nothing in either suite will catch you removing this.** It
		costs no accuracy -- the calibration table is byte for byte the
		same with it and without -- so both suites pass either way, and
		two mutations that disable it were caught by nothing. What it
		costs is speed. Get the formula wrong rather than removing it and
		the six reduction-at-scale tests do fire.
	*/

	long l_neededPlaces = glWorkingPrecision + reduced->GetOrderOfMagnitude();

	/*
		All four constants are base ten and are subtracted from or
		divided into the argument, so they come into its base first --
		see CLZAnswerBase -- and are narrowed after, in places of that
		base. Narrowed first, they held the working precision counted as
		decimal places, and a hex argument met a pi of forty digits
		where it wanted 48.
	*/

	lzi_pi_family_in_base( two_pi, kTwoPi, reduced->m_lBase, l_neededPlaces );

	if ( two_pi->m_lRadix > l_neededPlaces )
	{
		two_pi->SetPrecision( l_neededPlaces );
		two_pi->Round();
	}

	/*
		And the three the folds use, to the plain working width.

		These carry 120 places too, for the series that runs after
		them, but nothing here multiplies their error: each is
		subtracted once from a value the modulo has already brought
		into [0, 2*pi). So the width the answer is going to be rounded
		to is the width worth keeping, with no magnitude term of the
		kind 2*pi needs.

		Skipping this cost 1.10 to 1.61 times across the circular
		family -- three comparisons and up to two subtractions a call,
		every one of them against a constant more than twice as wide as
		the answer. That is the third time in this engine a widened
		constant has had to be narrowed again at the point of use, and
		it will not be the last.

		Like the 2*pi narrowing above, **no test catches you removing
		this**: it costs no accuracy, so both suites pass either way
		and a mutation disabling it was caught by nothing.
	*/

	lzi_pi_family_in_base( pi, kPi, reduced->m_lBase, glWorkingPrecision );
	lzi_pi_family_in_base( pi_over_two, kPiOverTwo, reduced->m_lBase, glWorkingPrecision );
	lzi_pi_family_in_base( pi_over_four, kPiOverFour, reduced->m_lBase, glWorkingPrecision );

	if ( pi->m_lRadix > glWorkingPrecision )
	{
		pi->SetPrecision( glWorkingPrecision );
		pi->Round();

		pi_over_two->SetPrecision( glWorkingPrecision );
		pi_over_two->Round();

		pi_over_four->SetPrecision( glWorkingPrecision );
		pi_over_four->Round();
	}

	// Into [0, 2pi)...
	lResult = lzi_first_failure( lResult, lz_modulo( reduced, two_pi, scratch ) );
	reduced->Assign( scratch );

	// sin(x) = -sin(x - pi), into [0, pi)...
	if ( lz_compare( reduced, pi ) == kGreaterThan )
	{
		lz_subtract( reduced, pi );
		b_negateResult = !b_negateResult;
		b_negateCosine = !b_negateCosine;		// cos(x) = -cos(x - pi)
	}

	// sin(x) = sin(pi - x), into [0, pi/2]...
	if ( lz_compare( reduced, pi_over_two ) == kGreaterThan )
	{
		lz_subtract( pi, reduced, scratch );
		reduced->Assign( scratch );
		b_negateCosine = !b_negateCosine;		// cos(x) = -cos(pi - x)
	}

	// sin(x) = cos(pi/2 - x), into [0, pi/4]...
	if ( lz_compare( reduced, pi_over_four ) == kGreaterThan )
	{
		lz_subtract( pi_over_two, reduced, scratch );
		reduced->Assign( scratch );
		b_useCosineSeries = true;
	}

	*pbNegateSine = b_negateResult;
	*pbNegateCosine = b_negateCosine;
	*pbSwapped = b_useCosineSeries;

	return lResult;
}


long lz_sine( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_sine, kRetakeDoubtSine );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	//	In radians a decider reads the argument as the angle itself; in a
	//	unit an angle a hair from one whose sine is a half is decided by its
	//	side of it -- see lzi_unit_half_decide.
	if ( e_unit == AU_radians )
		answer.DecideTies( lzi_inside_x_decide );
	else
		answer.DecideTies( lzi_sine_unit_decide );

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// In another unit, the angle reduced exactly first -- see
	//	lzi_circular_in_unit...
	if ( e_unit != AU_radians )
		return lzi_circular_in_unit( kAngleSine, leftOperand, targetNumber, answer, e_unit );

	long lResult = kNoError;

	/*
		Everything below runs wider than the answer and is rounded back
		at the end -- see kTrigGuardDigits. The raise has to come
		before the working values are created, because a CLZNumber
		takes its precision from glWorkingPrecision at construction and
		does not look again.
	*/

	long l_savedPrecision = glWorkingPrecision;

	glWorkingPrecision += lzi_guard_places( kTrigGuardDigits, leftOperand->m_lBase );

	auto reduced = new CLZNumber;
	auto result = new CLZNumber;

	bool b_negateResult = false;
	bool b_useCosineSeries = false;
	bool b_negateCosine = false;		// the sine needs only its own sign

	// Folded into [0, pi/4] -- see lzi_trig_fold...
	lResult = lzi_first_failure( lResult, lzi_trig_fold( leftOperand, reduced, &b_negateResult,
														 &b_negateCosine, &b_useCosineSeries ) );

	// Calculate...
	lzi_trig_series( reduced, b_useCosineSeries, result );

	if ( b_negateResult )
		lz_negate( result );

	// The guard places never leave the function...
	glWorkingPrecision = l_savedPrecision;

	result->SetPrecision( glWorkingPrecision );
	result->Round();

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
	delete reduced;
	delete result;

	return lResult;
}


/*
	**A sine or cosine the mathematics makes exact**: n / d, with d one or
	two. Nought, plus or minus one and plus or minus a half are the only
	rational sines and cosines of a rational part of a turn (Niven's
	theorem), and in radians only an angle of nought has one.
*/
struct SLZExactPart
{
	bool	m_bExact;
	long	m_lNumerator;
	long	m_lDenominator;
};

/*
	**Sine and cosine of one angle, from one reduction and one run** --
	stage 3 of the complex layer, 2026-10-22, for lz_sine_cosine and
	lzc_polar. At the working precision as it stands, the caller having
	raised it by its guard, and nothing rounded back.

	**In radians** the angle is folded into [0, pi/4] as lz_sine folds it
	-- see lzi_trig_fold. Past the kernel's floor the kernel gives both
	from one run, as lzi_sine_and_cosine has them for the tangent. Below
	it the series gives sin theta, and the cosine is the square root of
	1 - sin^2 theta: theta being at most pi/4, the cosine is at least
	0.707, so nothing cancels, and the root's slope there, sin over cos,
	is at most one, so the sine's error is not magnified. Measured at 36
	places, 6 to 9 microseconds against a cosine of its own at 17 to 48;
	at 300, 19 to 26 against 225 to 710.

	**In a unit** the angle comes to theta in [0, an eighth] exactly --
	see lzi_reduce_in_unit -- and only theta is converted, as
	lzi_circular_in_unit converts it. Each answer is plus or minus sin
	theta or cos theta by its quadrant.

	**Exact parts come back marked, their value meaningless**: an angle of
	nought answers 0 and 1, and in a unit theta of nought does, and a
	twelfth of a turn a sine of one half, which does not end in an odd
	base. bExactAllowed false works every part out instead, for a caller
	whose angle stands in for one given in another base -- see lzc_polar.
*/

static long lzi_sine_cosine_work( CLZNumber* pAngle, AngleUnit eUnit, bool bExactAllowed,
								  CLZNumber* pSine, CLZNumber* pCosine,
								  SLZExactPart* pSineExact, SLZExactPart* pCosineExact )
{
	long	l_base = pAngle->m_lBase;
	long	lResult = kNoError;

	*pSineExact = { false, 0, 1 };
	*pCosineExact = { false, 0, 1 };

	if ( eUnit == AU_radians )
	{
		if ( bExactAllowed && lzi_exactly( pAngle, 0 ) )
		{
			*pSineExact = { true, 0, 1 };
			*pCosineExact = { true, 1, 1 };

			return kNoError;
		}

		CLZNumber	reduced;
		CLZNumber	sine;
		CLZNumber	cosine;
		bool		b_negateSine = false;
		bool		b_negateCosine = false;
		bool		b_swapped = false;

		lResult = lzi_trig_fold( pAngle, &reduced, &b_negateSine, &b_negateCosine, &b_swapped );

		if ( glWorkingPrecision >= lzi_sincos_split_floor( l_base ) )
			lzi_sincos_split_work( &reduced, &sine, &cosine, 0 );
		else
		{
			CLZNumber	one;
			CLZNumber	square;
			CLZNumber	difference;

			lzi_trig_series( &reduced, false, &sine );

			one.SetBase( l_base, false );
			one.Assign( 1L );
			lz_multiply( &sine, &sine, &square );
			lz_subtract( &one, &square, &difference );
			lResult = lzi_first_failure( lResult, lz_root_2( &difference, &cosine ) );
		}

		// Folded past pi/4 the reduced argument's cosine is the sine...
		CLZNumber*	p_sine = b_swapped ? &cosine : &sine;
		CLZNumber*	p_cosine = b_swapped ? &sine : &cosine;

		if ( b_negateSine )
			lz_negate( p_sine );

		if ( b_negateCosine )
			lz_negate( p_cosine );

		pSine->Assign( p_sine );
		pCosine->Assign( p_cosine );

		return lResult;
	}

	CLZNumber	theta;
	long		l_quadrant = 0;
	bool		b_complement = false;

	lResult = lzi_reduce_in_unit( pAngle, eUnit, &theta, &l_quadrant, &b_complement );

	// sin theta and cos theta, exact where they are...
	CLZNumber		sinTheta;
	CLZNumber		cosTheta;
	SLZExactPart	sin_exact = { false, 0, 1 };
	SLZExactPart	cos_exact = { false, 0, 1 };
	bool			b_twelfth = false;

	if ( bExactAllowed && !(theta.IsZero()) )
	{
		CLZNumber	twelve;
		CLZNumber	turn;

		turn.SetBase( l_base, false );
		turn.Assign( lzi_turn_of( eUnit ) );
		lzi_times_whole( &theta, 12, &twelve );
		b_twelfth = (lz_compare( &twelve, &turn ) == kEqualTo);
	}

	if ( bExactAllowed && theta.IsZero() )
	{
		sin_exact = { true, 0, 1 };
		cos_exact = { true, 1, 1 };
	}
	else if ( theta.IsZero() )
	{
		sinTheta.SetBase( l_base, false );
		sinTheta.Assign( 0L );
		cosTheta.SetBase( l_base, false );
		cosTheta.Assign( 1L );
	}
	else
	{
		// Theta's noughts again, and the factor's three at most -- see
		//	lzi_circular_in_unit...
		long	l_lead = lzi_leading_exponent( &theta );
		long	l_small = ((l_lead < 0) ? -l_lead : 0) + 3;

		CLZNumber		radians;
		SLZExactPart	unused_sine;
		SLZExactPart	unused_cosine;

		lzi_unit_to_radians( &theta, eUnit, glWorkingPrecision + l_small, &radians );
		lResult = lzi_first_failure( lResult, lzi_sine_cosine_work( &radians, AU_radians, false, &sinTheta, &cosTheta,
																	&unused_sine, &unused_cosine ) );

		if ( b_twelfth )
			sin_exact = { true, 1, 2 };
	}

	// ...each answer plus or minus one of them, by quadrant -- see
	//	lzi_function_of_theta...
	struct { CLZNumber* m_pValue; SLZExactPart* m_pExact; long m_lFunction; } sl_parts[] =
	{
		{ pSine, pSineExact, kAngleSine },
		{ pCosine, pCosineExact, kAngleCosine },
	};

	for ( auto& part : sl_parts )
	{
		bool	b_negate = false;
		long	l_of = lzi_function_of_theta( part.m_lFunction, l_quadrant, b_complement, &b_negate );

		*part.m_pExact = (l_of == kAngleSine) ? sin_exact : cos_exact;

		if ( part.m_pExact->m_bExact )
		{
			if ( b_negate )
				part.m_pExact->m_lNumerator = -part.m_pExact->m_lNumerator;

			continue;
		}

		part.m_pValue->Assign( (l_of == kAngleSine) ? &sinTheta : &cosTheta );

		if ( b_negate )
			lz_negate( part.m_pValue );
	}

	return lResult;
}

/*
	**Sine and cosine of one angle, each answered exactly as lz_sine and
	lz_cosine answer it** -- stage 3 of the complex layer, 2026-10-22, on
	Rich's calls that day: this name, and both targets required and
	different, either of them the operand if wanted.

	**One reduction and one run give both** -- see lzi_sine_cosine_work --
	and each is rounded once, at the working precision. An answer within
	its function's doubt window, an exact one, a flagged one, or any where
	the work failed is that function's to answer instead -- see
	lzi_rounded_clear_of_doubt -- outside this function's own guard, so as
	the host's call: in its unit, taken again wider, ties decided. So every
	answer is the function's to the character, flags and codes included,
	in every unit and base.

	The angle is copied first, so a target that is the operand cannot reach
	the other answer's work.
*/

long lz_sine_cosine( CLZNumber* leftOperand, CLZNumber* sineTarget, CLZNumber* cosineTarget )
{
	if ( (leftOperand == nullptr) || (sineTarget == nullptr) || (cosineTarget == nullptr) ||
		 (sineTarget == cosineTarget) )
		return lzi_pair_refused( sineTarget, cosineTarget );

	CLZNumber	angle;
	CLZNumber	sine;
	CLZNumber	cosine;
	bool		b_sineOwn = true;
	bool		b_cosineOwn = true;

	angle.Assign( leftOperand );

	if ( angle.FlagIsOff( kfQNaN ) && angle.FlagIsOff( kfInfinity ) && !lzi_exactly( &angle, 0 ) )
	{
		// The work under a guard of its own, so every call inside nests and
		//	the host's unit is read here -- see lzi_angle_unit_here...
		CLZNumber		scratch;
		CLZAnswerBase	work( &angle, &scratch );

		work.Exact();

		AngleUnit		e_unit = lzi_angle_unit_here();
		long			l_saved = glWorkingPrecision;
		long			l_guard = lzi_guard_places( kSineCosineGuardDigits, angle.m_lBase );
		SLZExactPart	sine_exact;
		SLZExactPart	cosine_exact;

		glWorkingPrecision += l_guard;

		long	lResult = lzi_sine_cosine_work( &angle, e_unit, true, &sine, &cosine, &sine_exact, &cosine_exact );

		glWorkingPrecision = l_saved;

		if ( lResult == kNoError )
		{
			b_sineOwn = sine_exact.m_bExact || !lzi_rounded_clear_of_doubt( &sine, kRetakeDoubtSine, false );
			b_cosineOwn = cosine_exact.m_bExact || !lzi_rounded_clear_of_doubt( &cosine, kRetakeDoubt, false );
		}
	}

	// ...and an answer that cannot stand its own function's, as the host's
	//	call...
	glPairHandoffs += (b_sineOwn ? 1 : 0) + (b_cosineOwn ? 1 : 0);

	long	l_sineResult = b_sineOwn ? lz_sine( &angle, sineTarget )
									 : lzi_answer_as_function( &sine, &angle, sineTarget, false );
	long	l_cosineResult = b_cosineOwn ? lz_cosine( &angle, cosineTarget )
										 : lzi_answer_as_function( &cosine, &angle, cosineTarget, false );

	return lzi_first_failure( l_sineResult, l_cosineResult );
}


long lz_square( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// A product of the operand, which writes short only what is exact --
	//	see lzi_write_to_width...
	answer.Exact();

	CLZNumber* result = nullptr;

	result = new CLZNumber;

	lz_multiply( leftOperand, leftOperand, result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;

	return kNoError;
}


long lz_subtract( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber )
{
	// The target's rounding settings are its own -- see CLZTargetRounding...
	CLZTargetRounding	targetRounding( (targetNumber != nullptr) ? targetNumber : leftOperand );

	// Before anything -- see CLZMatchedBase.
	//	A digit array carries no base, so a mixed pair has to be
	//	brought into one before any of it is read...
	// A host's mixed pair, worked exactly -- see CLZMixedPair...
	if ( CLZMixedPair::Wanted( leftOperand, rightOperand ) )
	{
		CLZMixedPair	mixed( leftOperand, rightOperand );

		return mixed.Add( targetNumber, true );
	}

	CLZMatchedBase	matched( leftOperand, rightOperand );

	rightOperand = matched.Operand();

	// An infinity or a NaN carries through rather than being lost...
	if ( lzi_propagate_flags( leftOperand, rightOperand, targetNumber, kOpSubtract ) )
		return kNoError;

	CLZNumber*	maxOperand = nullptr;
	CLZNumber*	minOperand = nullptr;

	auto result = new CLZNumber;

	unsigned long	maxDigits = 0;
	unsigned long	maxExponent = 0;

	// Set up the result as required...
	if ( (leftOperand->m_lDigitsLength - leftOperand->m_lRadix) >
		 (rightOperand->m_lDigitsLength - rightOperand->m_lRadix) )
		maxDigits = leftOperand->m_lDigitsLength - leftOperand->m_lRadix;

	else
		maxDigits = rightOperand->m_lDigitsLength - rightOperand->m_lRadix;

	if ( leftOperand->m_lRadix > rightOperand->m_lRadix )
		maxExponent = leftOperand->m_lRadix;

	else
		maxExponent = rightOperand->m_lRadix;

	// Identify the 'big' operand...
	if ( leftOperand->CompareABS( rightOperand ) == kLessThan )
	{
		maxOperand = rightOperand;
		minOperand = leftOperand;
	}
	else
	{
		maxOperand = leftOperand;
		minOperand = rightOperand;
	}

	// Align target storage...
	result->Assign( maxOperand );
	result->AlignStorage( maxDigits, maxExponent, kPreserve );

	// Perform the operation...
	if ( leftOperand->m_lSign == rightOperand->m_lSign )
	{
		if ( leftOperand == maxOperand )
			lzi_subtract_op( maxOperand, minOperand, result, maxOperand->m_lSign );
		else
			lzi_subtract_op( maxOperand, minOperand, result, (maxOperand->m_lSign)==kPositive?kNegative:kPositive );
	}

	else if ( leftOperand->m_lSign == kZero )
	{
		if ( rightOperand->m_lSign == kZero )
			lzi_add_op( maxOperand, minOperand, result, kZero );
		else
			lzi_add_op( maxOperand, minOperand, result, (rightOperand->m_lSign)==kPositive?kNegative:kPositive );
	}

	else
		lzi_add_op( maxOperand, minOperand, result, leftOperand->m_lSign );

	// Record the result, a whole number if the target asks for one...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( result );
			targetRounding.RoundExactIfWhole();
		}
	}
	else if ( !(leftOperand->IsConstant()) )
	{
		leftOperand->Assign( result );
		targetRounding.RoundExactIfWhole();
	}

	// Clean up...
	delete result;

	return kNoError;
}


long lz_tangent( CLZNumber* leftOperand, CLZNumber* targetNumber )
{
	// The answer in the left operand's base -- see CLZAnswerBase...
	CLZAnswerBase	answer( leftOperand, targetNumber );

	// Taken again wider if the guard leaves the rounding in doubt -- see
	//	CLZAnswerBase::Retake...
	answer.Retake( lz_tangent, kRetakeDoubt );

	// The host's unit, if not radians: asked here, before anything nests
	//	-- see lzi_angle_unit_here...
	AngleUnit	e_unit = lzi_angle_unit_here();

	// ...and a tiny argument on a half decided by its side of x -- see
	//	lzi_beside_x_decide...
	//	In radians only: a decider reads the argument as the angle itself.
	if ( e_unit == AU_radians )
		answer.DecideTies( lzi_outside_x_decide );

	// A small answer is as wide as a quotient's -- see lzi_write_to_width...
	answer.QuotientWidth();

	// Exact at nought and nowhere else, asked before anything is done to
	//	the operand -- see lzi_write_to_width...
	answer.ExactAndShort( lzi_exactly( leftOperand, 0 ) );

	// In another unit, the angle reduced exactly first -- see
	//	lzi_circular_in_unit...
	if ( e_unit != AU_radians )
		return lzi_circular_in_unit( kAngleTangent, leftOperand, targetNumber, answer, e_unit );

	long lResult = kNoError;

	// Guarded across the division -- see kTrigQuotientGuardDigits, and the
	//	note in lz_cotangent on why widening the sine and cosine
	//	themselves is not enough...
	long l_savedPrecision = glWorkingPrecision;

	long	l_guard = lzi_guard_places( kTrigQuotientGuardDigits, leftOperand->m_lBase );

	glWorkingPrecision += l_guard;

	auto result = new CLZNumber;
	auto sine_result = new CLZNumber;
	auto cosine_result = new CLZNumber;

	/*
		**Wider by the answer's size wherever the sine or the cosine is
		small**, and never both are. Until 2026-09-26 this was not
		widened at all.

		- **A small cosine is a small divisor**, near the odd multiples
		  of pi/2, and wants twice its size -- see
		  lzi_small_divisor_widening. The answer is large, and its
		  numerator needs the width as well as its divisor, so the
		  cosine goes first and the sine is taken after it at the width
		  it asked for. tan(1.5707963268) was 7e17 units out of its 36th
		  place, and at pi/2 to twenty digits 9e36.
		- **A small sine is a small numerator and a small answer**, near
		  the multiples of pi, which Rich's rule gives 36 significant
		  digits -- see lzi_round_quotient_answer. A numerator about
		  base^-k worked to W places carries W - k of them, so it wants
		  its size once. Near nought the argument says how small, before
		  the sine is taken, as in lz_cosecant.

		Each is taken again wider only when its size would spend more
		than half the guard -- see lz_cosecant -- so an ordinary argument
		takes one sine and one cosine as it always did. lz_cotangent is the same with the two exchanged,
		and lz_secant is this function's divisor alone.
	*/

	/*
		**One fold and one kernel run past the kernel's floor** -- see
		lzi_sine_and_cosine -- wherever the widening below would take only
		the argument's own: a cosine and a sine both seen, and neither so
		small that it asks for more than half the guard past it. The two
		below would then be taken a place or so apart -- the sine as wide as
		a small numerator's argument asks -- and these are taken at once at
		the sine's width, the cosine a little wider than it needs. Anything
		else puts the width back and takes them as it always did.
	*/

	bool	b_onePair = false;
	long	l_numerator = lzi_small_divisor_widening( leftOperand ) / 2;

	if ( ((glWorkingPrecision + l_numerator) >= lzi_sincos_split_floor( leftOperand->m_lBase )) &&
		 leftOperand->FlagIsOff( kfInfinity ) && leftOperand->FlagIsOff( kfQNaN ) )
	{
		glWorkingPrecision += l_numerator;

		lzi_sine_and_cosine( leftOperand, sine_result, cosine_result );

		sine_result->CheckForZero();
		cosine_result->CheckForZero();

		b_onePair = !(sine_result->IsZero()) && !(cosine_result->IsZero()) &&
					(lzi_small_divisor_widening( cosine_result ) <= (l_guard / 2)) &&
					(((lzi_small_divisor_widening( sine_result ) / 2) - l_numerator) <= (l_guard / 2));

		if ( !b_onePair )
			glWorkingPrecision -= l_numerator;
	}

	if ( !b_onePair )
	{
		lz_cosine( leftOperand, cosine_result );

		// A cosine too small to show at this width is taken wider until it
		//	does -- see lzi_widen_until_seen...
		long	l_widened = lzi_widen_until_seen( lz_cosine, leftOperand, cosine_result );
		long	l_wanted = lzi_small_divisor_widening( cosine_result );

		if ( (l_wanted - l_widened) > (l_guard / 2) )
		{
			glWorkingPrecision += l_wanted - l_widened;
			l_widened = l_wanted;

			lz_cosine( leftOperand, cosine_result );
		}

		l_wanted = lzi_small_divisor_widening( leftOperand ) / 2;

		if ( l_wanted > l_widened )
		{
			glWorkingPrecision += l_wanted - l_widened;
			l_widened = l_wanted;
		}

		lz_sine( leftOperand, sine_result );

		l_widened += lzi_widen_until_seen( lz_sine, leftOperand, sine_result );
		l_wanted = lzi_small_divisor_widening( sine_result ) / 2;

		if ( (l_wanted - l_widened) > (l_guard / 2) )
		{
			glWorkingPrecision += l_wanted - l_widened;

			lz_sine( leftOperand, sine_result );
		}
	}

	// lz_divide takes its width from its target, made before any of
	//	this...
	result->SetPrecision( glWorkingPrecision );

	lResult = lzi_first_failure( lResult, lz_divide( sine_result, cosine_result, result ) );

	// The guard places never leave the function, and the answer is
	//	rounded to the rule for a quotient -- see
	//	lzi_round_quotient_answer...
	glWorkingPrecision = l_savedPrecision;

	lzi_round_quotient_answer( result );

	// Record the result...
	if ( targetNumber != nullptr )
	{
		if ( !(targetNumber->IsConstant()) )
			targetNumber->Assign( result );
	}
	else if ( !(leftOperand->IsConstant()) )
		leftOperand->Assign( result );

	// Clean up...
    delete result;
    delete sine_result;
    delete cosine_result;

	return lResult;
}

/////////////////////////////////////////////////////////////////
//	Complex numbers -- stage 1, 2026-10-21

/*
	**Each part is rounded to the working precision in places, in its own
	rounding mode** -- Rich's call, 2026-10-21, with the layer here rather
	than in a third source file, and stage 1 the operations that are
	correctly rounded from exact intermediates.

	**Nothing here rounds a part itself.** Each part of an answer is
	worked out exactly -- products at a width that holds every place,
	sums and differences exact as they always are -- and then handed to
	one real operation writing into the target's own part: a
	multiplication, a division or a square root. So a part is rounded,
	flagged and spelled exactly as that operation rounds, flags and
	spells, and the rules for doing so live in one place. Rich's contract
	for the real ones carries over whole: the four plain operations round
	in the target part's own mode and digit, and the square and the
	modulus answer as functions do, the same whatever the target's
	settings. With both imaginary parts nought, every answer here is the
	real operation's, character for character.

	**A value's two parts are in one base**, the value's own: a digit
	array carries no base, and neither does a pair of them. A value whose
	parts differ is refused, kBaseMismatch. A right operand in another
	base is brought into the left's exactly -- see CLZComplexInBase.
*/

/*
	**An infinity anywhere is the point at infinity.** This engine holds
	an infinity as nought with a flag, unsigned, so a complex one has no
	direction either: a value with an infinite part, and no not-a-number,
	is the one point at infinity, and comes back with both parts
	infinite. The rules are lzi_propagate_flags', lifted: a not-a-number
	anywhere begets one; infinity plus or minus infinity is a
	not-a-number, as are infinity times nought and infinity over
	infinity; anything else involving one is infinite, but for something
	finite over infinity, which is nought.
*/

const long	kComplexAdd = 0;
const long	kComplexSubtract = 1;
const long	kComplexMultiply = 2;
const long	kComplexDivide = 3;
const long	kComplexUnary = 4;

static bool lzi_complex_has( CLZComplex* pValue, unsigned long lFlag )
{
	return pValue->m_real.FlagIsOn( lFlag ) || pValue->m_imag.FlagIsOn( lFlag );
}

// IsZero reports a flag, brought up to date first...
static bool lzi_complex_is_zero( CLZComplex* pValue )
{
	pValue->m_real.CheckForZero();
	pValue->m_imag.CheckForZero();

	return pValue->m_real.IsZero() && pValue->m_imag.IsZero();
}

// Both parts of the answer nought in lBase, flagged as lFlag says -- none
//	for a plain nought -- and carrying kfDiv0 if the trouble began there...
static void lzi_complex_flagged( CLZComplex* pTarget, long lBase, unsigned long lFlag, bool bFromDivideByZero )
{
	for ( CLZNumber* p_part : { &pTarget->m_real, &pTarget->m_imag } )
	{
		if ( p_part->IsConstant() )
			continue;

		p_part->Assign( 0L );
		p_part->SetBase( lBase, false );

		if ( lFlag != 0 )
			p_part->SetFlag( lFlag );

		if ( bFromDivideByZero )
			p_part->SetFlag( kfDiv0 );
	}
}

// A value whose two parts are in different bases, refused with both parts
//	of the answer a not-a-number, answering whether it was...
static bool lzi_complex_refused( CLZComplex* pLeft, CLZComplex* pRight, CLZComplex* pTarget, long* plResult )
{
	if ( (pLeft->m_real.m_lBase == pLeft->m_imag.m_lBase) &&
		 ((pRight == nullptr) || (pRight->m_real.m_lBase == pRight->m_imag.m_lBase)) )
		return false;

	lzi_complex_flagged( pTarget, pLeft->m_real.m_lBase, kfQNaN, false );
	*plResult = kBaseMismatch;

	return true;
}

// The point at infinity and the not-a-number carried through, answering
//	whether either was -- pRight null for an operation with one operand...
static bool lzi_complex_propagate( CLZComplex* pLeft, CLZComplex* pRight, CLZComplex* pTarget, long lOperation )
{
	bool	b_leftNaN = lzi_complex_has( pLeft, kfQNaN );
	bool	b_rightNaN = (pRight != nullptr) && lzi_complex_has( pRight, kfQNaN );
	bool	b_leftInfinite = lzi_complex_has( pLeft, kfInfinity );
	bool	b_rightInfinite = (pRight != nullptr) && lzi_complex_has( pRight, kfInfinity );

	if ( !(b_leftNaN || b_rightNaN || b_leftInfinite || b_rightInfinite) )
		return false;

	unsigned long	l_resultFlag = kfInfinity;

	if ( b_leftNaN || b_rightNaN )
		l_resultFlag = kfQNaN;

	else if ( lOperation == kComplexUnary )
		l_resultFlag = kfInfinity;

	else if ( b_leftInfinite && b_rightInfinite )
		l_resultFlag = (lOperation == kComplexMultiply) ? kfInfinity : kfQNaN;

	else if ( b_leftInfinite )
		l_resultFlag = ((lOperation == kComplexMultiply) && lzi_complex_is_zero( pRight )) ? kfQNaN : kfInfinity;

	else if ( lOperation == kComplexDivide )
		l_resultFlag = 0;				// something finite over infinity is nought

	else
		l_resultFlag = ((lOperation == kComplexMultiply) && lzi_complex_is_zero( pLeft )) ? kfQNaN : kfInfinity;

	bool	b_fromDivideByZero = lzi_complex_has( pLeft, kfDiv0 ) ||
								 ((pRight != nullptr) && lzi_complex_has( pRight, kfDiv0 ));

	lzi_complex_flagged( pTarget, pLeft->m_real.m_lBase, l_resultFlag, b_fromDivideByZero );

	return true;
}

/*
	**A right operand in the left's base, exactly** -- CLZMixedPair's way
	with a complex value. Where every prime of its base divides the
	left's, both parts convert exactly. Where not, the value is held as a
	whole number over a power of its own base -- (Yr + i Yi) / D, with
	D = b^k for the more places of its two parts -- and whole numbers
	convert exactly. An answer over D is a quotient and is written by the
	quotient's rule, as a real mixed product is.

	**Nothing is made when the bases agree**: a number callocs when it is
	made, which is what cost the same-base addition 41 per cent when
	CLZMatchedBase first held its copy by value.
*/

class CLZComplexInBase
{
public:

	CLZComplexInBase( CLZComplex* pValue, long lBase )
		: m_pOwned( nullptr ), m_pReal( &pValue->m_real ), m_pImag( &pValue->m_imag ), m_bScaled( false )
	{
		long	l_from = pValue->m_real.m_lBase;

		if ( l_from == lBase )
			return;

		m_pOwned = new CLZNumber[ 3 ];
		m_pOwned[ 0 ].Assign( &pValue->m_real );
		m_pOwned[ 1 ].Assign( &pValue->m_imag );
		m_pReal = &m_pOwned[ 0 ];
		m_pImag = &m_pOwned[ 1 ];

		if ( lzi_ends_in_base( l_from, lBase ) )
		{
			lzi_into_base_keeping( m_pReal, lBase, 0 );
			lzi_into_base_keeping( m_pImag, lBase, 0 );
			return;
		}

		m_bScaled = true;
		m_pReal->TrimTrailing();
		m_pImag->TrimTrailing();

		long	l_places = (m_pReal->m_lRadix > m_pImag->m_lRadix) ? m_pReal->m_lRadix : m_pImag->m_lRadix;
		CLZNumber*	p_scale = &m_pOwned[ 2 ];

		lzi_scale_by_base_power( m_pReal, -l_places );
		lzi_scale_by_base_power( m_pImag, -l_places );
		p_scale->SetBase( l_from, false );
		p_scale->Assign( 1L );
		lzi_scale_by_base_power( p_scale, -l_places );

		m_pReal->SetBase( lBase, true );
		m_pImag->SetBase( lBase, true );
		p_scale->SetBase( lBase, true );
	}

	~CLZComplexInBase()
	{
		delete [] m_pOwned;
	}

	CLZNumber*	Real() const	{ return m_pReal; }
	CLZNumber*	Imag() const	{ return m_pImag; }
	bool		Scaled() const	{ return m_bScaled; }
	CLZNumber*	Scale() const	{ return &m_pOwned[ 2 ]; }

private:

	CLZComplexInBase( const CLZComplexInBase& );
	CLZComplexInBase& operator=( const CLZComplexInBase& );

	CLZNumber*	m_pOwned;
	CLZNumber*	m_pReal;
	CLZNumber*	m_pImag;
	bool		m_bScaled;
};

// x times y exactly, both in one base: the working precision widened to
//	hold every place of the product...
static void lzi_complex_product( CLZNumber* x, CLZNumber* y, CLZNumber* pProduct )
{
	long	l_saved = glWorkingPrecision;
	long	l_places = x->m_lRadix + y->m_lRadix;

	if ( l_places > glWorkingPrecision )
		glWorkingPrecision = l_places;

	lz_multiply( x, y, pProduct );

	glWorkingPrecision = l_saved;
}

/*
	**An exact value rounded as a product is**, into a target part --
	through lz_multiply by one, so the rule for rounding a product (the
	part's own mode and digit, a whole number if it asks for one, trimmed
	where it ends and the full width where it does not) is written once,
	there, and this cannot drift from it. One is a single digit beside
	the products that made the value.
*/
static long lzi_complex_round_as_product( CLZNumber* pExact, CLZNumber* pPart )
{
	CLZNumber	one;

	one.SetBase( pExact->m_lBase, false );
	one.Assign( 1L );

	return lz_multiply( pExact, &one, pPart );
}

/*
	**An exact value rounded as lz_square rounds its product** -- a
	function's answer, the same whatever the target's settings, which is
	Rich's contract for every operation but the four plain ones. Through
	CLZAnswerBase as lz_square goes, the product into a fresh number and
	that into the target, so the two cannot differ...
*/
static long lzi_complex_round_as_square( CLZNumber* pExact, CLZNumber* pPart )
{
	CLZAnswerBase	answer( pExact, pPart );

	answer.Exact();

	CLZNumber	result;
	long		lResult = lzi_complex_round_as_product( pExact, &result );

	if ( !(pPart->IsConstant()) )
		pPart->Assign( &result );

	return lResult;
}

// Two exact parts into the target: rounded as products, or divided once by
//	the right operand's scale where it was held over one...
static long lzi_complex_finish( CLZNumber* pReal, CLZNumber* pImag, const CLZComplexInBase& right,
								CLZComplex* pTarget )
{
	if ( right.Scaled() )
		return lzi_first_failure( lz_divide( pReal, right.Scale(), &pTarget->m_real ),
								  lz_divide( pImag, right.Scale(), &pTarget->m_imag ) );

	return lzi_first_failure( lzi_complex_round_as_product( pReal, &pTarget->m_real ),
							  lzi_complex_round_as_product( pImag, &pTarget->m_imag ) );
}


long lzc_abs( CLZComplex* leftOperand, CLZNumber* targetNumber )
{
	// Required, as the means' is: the answer is real, so there is no
	//	complex operand to answer in place...
	if ( targetNumber == nullptr )
		return kFailed;

	long	l_base = leftOperand->m_real.m_lBase;
	bool	b_flagged = (leftOperand->m_imag.m_lBase != l_base) || lzi_complex_has( leftOperand, kfQNaN ) ||
						lzi_complex_has( leftOperand, kfInfinity );

	if ( b_flagged )
	{
		// The modulus of the point at infinity is infinite...
		long			l_result = (leftOperand->m_imag.m_lBase != l_base) ? kBaseMismatch : kNoError;
		unsigned long	l_flag = ((l_result == kNoError) && !lzi_complex_has( leftOperand, kfQNaN )) ? kfInfinity : kfQNaN;

		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( 0L );
			targetNumber->SetBase( l_base, false );
			targetNumber->SetFlag( l_flag );

			if ( lzi_complex_has( leftOperand, kfDiv0 ) )
				targetNumber->SetFlag( kfDiv0 );
		}

		return l_result;
	}

	leftOperand->m_real.CheckForZero();
	leftOperand->m_imag.CheckForZero();

	// One part nought, and the modulus is the other's, exactly...
	if ( leftOperand->m_imag.IsZero() )
		return lz_abs( &leftOperand->m_real, targetNumber );

	if ( leftOperand->m_real.IsZero() )
		return lz_abs( &leftOperand->m_imag, targetNumber );

	// ...and otherwise the square root of an exact sum, rounded once...
	CLZNumber	aa;
	CLZNumber	bb;
	CLZNumber	sum;

	lzi_complex_product( &leftOperand->m_real, &leftOperand->m_real, &aa );
	lzi_complex_product( &leftOperand->m_imag, &leftOperand->m_imag, &bb );
	lz_add( &aa, &bb, &sum );

	return lz_root_2( &sum, targetNumber );
}


/*
	**The argument, the angle of the point the value names** -- stage 2a,
	2026-10-22: lz_arctangent2 of the imaginary part over the real one,
	which already answers everything this has to -- the quadrant, the
	axes by the +0 convention, nought for nought, and a not-a-number for
	a not-a-number or an infinity, the point at infinity having no
	direction. So this does nothing of its own but refuse what the layer
	refuses: a missing target, the answer being real, and a value whose
	two parts are in different bases.
*/

long lzc_arg( CLZComplex* leftOperand, CLZNumber* targetNumber )
{
	// Required, as lzc_abs's is...
	if ( targetNumber == nullptr )
		return kFailed;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( leftOperand->m_imag.m_lBase != l_base )
	{
		if ( !(targetNumber->IsConstant()) )
		{
			targetNumber->Assign( 0L );
			targetNumber->SetBase( l_base, false );
			targetNumber->SetFlag( kfQNaN );
		}

		return kBaseMismatch;
	}

	return lz_arctangent2( &leftOperand->m_imag, &leftOperand->m_real, targetNumber );
}


/*
	**An exact value times n / d, correctly rounded at lPlaces places**,
	half to the even kept digit -- lzc_polar's exact parts, r, -r and
	r / 2, n being nought or plus or minus one. In whole numbers, as
	lzi_part_of_turn rounds a part of a turn, since a half does not end in
	an odd base: the value is V / b^k, V whole and k its places, so the
	answer times b^P is V n b^(P - k) / d, or V n over d b^(k - P), and
	the whole quotient and remainder of that decide. Answers whether
	nothing was left over, the answer then exact and short.
*/
static bool lzi_ratio_at_places( CLZNumber* pValue, long lNumerator, long lDenominator, long lPlaces,
								 CLZNumber* pResult )
{
	long	l_base = pValue->m_lBase;

	pResult->SetBase( l_base, false );

	if ( lNumerator == 0 )
	{
		pResult->Assign( 0L );

		return true;
	}

	bool	b_negative = (pValue->IsNegative() != (lNumerator < 0));

	CLZNumber	scaled;
	CLZNumber	divisor;
	CLZNumber	whole;
	CLZNumber	product;
	CLZNumber	remainder;
	CLZNumber	twice;

	scaled.Assign( pValue );
	scaled.TrimTrailing();

	if ( scaled.IsNegative() )
		lz_negate( &scaled );

	long	l_places = scaled.m_lRadix;

	divisor.SetBase( l_base, false );
	divisor.Assign( lDenominator );

	// V, a whole number, and the power of the base onto one side or the
	//	other...
	lzi_scale_by_base_power( &scaled, -l_places );

	if ( lPlaces >= l_places )
		lzi_scale_by_base_power( &scaled, -(lPlaces - l_places) );
	else
		lzi_scale_by_base_power( &divisor, -(l_places - lPlaces) );

	// ...the whole quotient and remainder, exactly...
	lz_divide_integer( &scaled, &divisor, &whole );
	lz_multiply( &whole, &divisor, &product );
	lz_subtract( &scaled, &product, &remainder );
	lzi_times_whole( &remainder, 2, &twice );

	long	l_order = lz_compare( &twice, &divisor );
	bool	b_up = (l_order == kGreaterThan);

	// ...a tie to the even kept digit...
	if ( l_order == kEqualTo )
	{
		CLZNumber	base;
		CLZNumber	digit;
		CLZNumber	two;
		CLZNumber	parity;

		base.SetBase( l_base, false );
		base.Assign( l_base );
		two.SetBase( l_base, false );
		two.Assign( 2L );

		lz_remainder( &whole, &base, &digit );
		lz_remainder( &digit, &two, &parity );
		parity.CheckForZero();
		b_up = !(parity.IsZero());
	}

	if ( b_up )
	{
		CLZNumber	one;

		one.SetBase( l_base, false );
		one.Assign( 1L );
		lz_add( &whole, &one );
	}

	// ...back over b^P, short only where nothing was left over...
	lzi_scale_by_base_power( &whole, lPlaces );

	remainder.CheckForZero();

	bool	b_exact = remainder.IsZero();

	if ( b_exact )
		whole.TrimTrailing();

	whole.CheckForZero();

	if ( b_negative && !(whole.IsZero()) )
		lz_negate( &whole );

	pResult->Assign( &whole );

	return b_exact;
}

/*
	**A part of lzc_polar on a half, decided by a bound it lies under** --
	stage 3, 2026-10-22. Wherever neither part is exact, |r cos theta| and
	|r sin theta| are under |r|, and in radians |r sin theta| is under
	|r theta| too. So with y the part's magnitude cut at its places and u a
	unit there, a bound whose double is at most 2y + u -- on the half or
	below it -- puts the part below the half, and its magnitude rounds down
	to y. Anything else declines, and the retake goes on as before.

	**It is the tie a polar value meets by construction**: r with a place
	past the answer's, a half there, and theta tiny, so that the cosine is
	short of one by theta^2 / 2 -- at a depth the retake never reaches once
	theta is small enough. The census met it at random, three times in
	eleven thousand: 1747.875 |/ -1.8e-7 degrees at one place in base four
	answered the half's upper neighbour. As lzi_beside_x_decide, which
	decides the sine's half by its side of x.
*/
static bool lzi_polar_decide( CLZNumber* pResult, CLZNumber* pR, CLZNumber* pTheta, bool bSine, bool bRadians,
							  long lPlaces )
{
	if ( lPlaces < 0 )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_saved = glWorkingPrecision;

	CLZNumber	bounds[ 2 ];
	long		l_bounds = 1;

	bounds[ 0 ].Assign( pR );
	bounds[ 0 ].m_lSign = kPositive;
	bounds[ 0 ].CheckForZero();

	if ( bSine && bRadians && (pTheta->m_lBase == l_base) )
	{
		lzi_complex_product( pR, pTheta, &bounds[ 1 ] );
		bounds[ 1 ].m_lSign = kPositive;
		bounds[ 1 ].CheckForZero();
		l_bounds = 2;
	}

	long	l_places = lPlaces;

	for ( long i = 0; i < l_bounds; i++ )
	{
		if ( bounds[ i ].m_lRadix > l_places )
			l_places = bounds[ i ].m_lRadix;
	}

	glWorkingPrecision = l_places + 4;

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	twice_mid;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &twice_mid );
	lz_add( &twice_mid, &unit );

	bool	b_down = false;

	for ( long i = 0; (i < l_bounds) && !b_down; i++ )
	{
		CLZNumber	twice_bound;

		lz_add( &bounds[ i ], &bounds[ i ], &twice_bound );
		b_down = (lz_compare( &twice_bound, &twice_mid ) != kGreaterThan);
	}

	glWorkingPrecision = l_saved;

	if ( !b_down )
		return false;

	bool	b_negative = pResult->IsNegative();

	y.CheckForZero();

	if ( b_negative && !(y.IsZero()) )
		lz_negate( &y );

	pResult->Assign( &y );

	return true;
}

/*
	**A complex value from its modulus and its angle, r cos theta + i r sin
	theta** -- stage 3 of the complex layer, 2026-10-22: the constructor
	Elzed's r |/ theta wants. The
	angle is read in the host's unit, as every angle function reads one.

	**Each part is correctly rounded, half-even, at the working precision
	in places**, as a function's answer is. The sine and cosine come from
	one run -- see lzi_sine_cosine_work -- as wide again as r's whole
	digits, so the product with r, exact, keeps the guard; rounded once,
	and a part within kRetakeDoubt of the half taken again wider, as
	lzi_take_again_wider takes a function's answer.

	**Exact where the mathematics is**: at a whole number of quarter turns
	the parts are r, -r and an exact nought, and at a twelfth one part is
	r / 2 -- rounded in whole numbers where it does not end, by
	lzi_ratio_at_places. So in degrees 1 |/ 180 is -1 + 0i exactly, where
	Elzed's doubles give an imaginary part of 1.2e-16. Short where the part
	ends, as every exact answer is.

	**A negative modulus is taken as it is**, as Elzed's r |/ theta takes
	it: the point -r cos theta + i -r sin theta. Nought at any finite angle
	is nought, exactly. **The answer is in r's base**; an angle in another
	is brought into it as power brings its exponent -- see
	lzi_into_base_keeping -- with its exact parts read from the angle as
	given, since a twelfth of a turn in grads ends in base three and not in
	base ten. A not-a-number in either, or an infinite angle, which has no
	direction, answers a not-a-number; an infinite modulus at a finite
	angle the point at infinity. The target is required, kFailed if null;
	a missing operand is refused the same way, the target a not-a-number.
*/

long lzc_polar( CLZNumber* modulus, CLZNumber* angle, CLZComplex* targetNumber )
{
	// Required, as lzc_abs's is: no operand is complex...
	if ( targetNumber == nullptr )
		return kFailed;

	if ( (modulus == nullptr) || (angle == nullptr) )
	{
		lzi_complex_flagged( targetNumber, (modulus != nullptr) ? modulus->m_lBase :
										   ((angle != nullptr) ? angle->m_lBase : 10L), kfQNaN, false );

		return kFailed;
	}

	long	l_base = modulus->m_lBase;
	bool	b_fromDivideByZero = modulus->FlagIsOn( kfDiv0 ) || angle->FlagIsOn( kfDiv0 );

	if ( modulus->FlagIsOn( kfQNaN ) || angle->FlagIsOn( kfQNaN ) || angle->FlagIsOn( kfInfinity ) )
	{
		lzi_complex_flagged( targetNumber, l_base, kfQNaN, b_fromDivideByZero );

		return kNoError;
	}

	if ( modulus->FlagIsOn( kfInfinity ) )
	{
		lzi_complex_flagged( targetNumber, l_base, kfInfinity, b_fromDivideByZero );

		return kNoError;
	}

	// Copies, so a part of the target that is an operand cannot reach the
	//	work...
	CLZNumber	r;
	CLZNumber	theta;
	CLZNumber	parts[ 2 ];			// real, imaginary
	bool		b_short[ 2 ] = { true, true };
	long		lResult = kNoError;

	r.Assign( modulus );
	theta.Assign( angle );
	r.CheckForZero();

	parts[ 0 ].SetBase( l_base, false );
	parts[ 0 ].Assign( 0L );
	parts[ 1 ].SetBase( l_base, false );
	parts[ 1 ].Assign( 0L );

	if ( !(r.IsZero()) )
	{
		// The work under a guard of its own, so every call inside nests and
		//	the host's unit is read here -- see lzi_angle_unit_here...
		CLZNumber		scratch;
		CLZAnswerBase	work( &r, &scratch );

		work.Exact();

		AngleUnit	e_unit = lzi_angle_unit_here();
		long		l_saved = glWorkingPrecision;
		long		l_lead = lzi_leading_exponent( &r );
		long		l_magnitude = (l_lead >= 0) ? (l_lead + 1) : 0;
		long		l_guard = lzi_guard_places( kSineCosineGuardDigits, l_base );

		// An angle in another base that does not end in this one stands in
		//	for itself only approximately, so its exact parts are read from it
		//	as given -- the work at one place, only its marks wanted...
		bool			b_approximate = (theta.m_lBase != l_base) && !lzi_ends_in_base( theta.m_lBase, l_base );
		SLZExactPart	given[ 2 ];			// cosine, sine

		if ( b_approximate )
		{
			CLZNumber	sine;
			CLZNumber	cosine;

			glWorkingPrecision = 1;
			lResult = lzi_first_failure( lResult, lzi_sine_cosine_work( &theta, e_unit, true, &sine, &cosine,
																		&given[ 1 ], &given[ 0 ] ) );
			glWorkingPrecision = l_saved;
		}

		// The cosine and the sine at lWidth places of r's base, and which
		//	are exact...
		auto	fnValues = [ & ]( long lWidth, CLZNumber* pCosine, CLZNumber* pSine, SLZExactPart* pExact ) -> long
		{
			long	l_was = glWorkingPrecision;
			long	l_result = kNoError;

			glWorkingPrecision = lWidth;

			CLZNumber		sine;
			CLZNumber		cosine;
			CLZNumber		within;
			SLZExactPart	exact[ 2 ];

			within.Assign( &theta );

			if ( within.m_lBase != l_base )
			{
				long	l_thetaLead = lzi_leading_exponent( &within );

				lzi_into_base_keeping( &within, l_base, lWidth + ((l_thetaLead > 0) ? l_thetaLead : 0) + 4 );
			}

			l_result = lzi_sine_cosine_work( &within, e_unit, !b_approximate, &sine, &cosine, &exact[ 1 ], &exact[ 0 ] );

			if ( b_approximate )
			{
				exact[ 0 ] = given[ 0 ];
				exact[ 1 ] = given[ 1 ];
			}

			pCosine->Assign( &cosine );
			pSine->Assign( &sine );
			pExact[ 0 ] = exact[ 0 ];
			pExact[ 1 ] = exact[ 1 ];

			glWorkingPrecision = l_was;

			return l_result;
		};

		CLZNumber		values[ 2 ];
		SLZExactPart	exact[ 2 ];

		lResult = lzi_first_failure( lResult, fnValues( l_saved + l_guard + l_magnitude, &values[ 0 ], &values[ 1 ], exact ) );

		for ( long l_part = 0; l_part < 2; l_part++ )
		{
			if ( exact[ l_part ].m_bExact )
			{
				b_short[ l_part ] = lzi_ratio_at_places( &r, exact[ l_part ].m_lNumerator, exact[ l_part ].m_lDenominator,
														 l_saved, &parts[ l_part ] );
				continue;
			}

			b_short[ l_part ] = false;

			lzi_complex_product( &r, &values[ l_part ], &parts[ l_part ] );

			if ( lzi_rounded_clear_of_doubt( &parts[ l_part ], kRetakeDoubt, false ) )
				continue;

			// Within the window of the half: taken again wider, the extra
			//	doubling each round -- see lzi_take_again_wider...
			auto	fnCall = [ & ]( CLZNumber* pResult ) -> long
			{
				CLZNumber		cosine;
				CLZNumber		sine;
				SLZExactPart	unused[ 2 ];

				long	l_result = fnValues( glWorkingPrecision + l_guard + l_magnitude, &cosine, &sine, unused );

				lzi_complex_product( &r, (l_part == 0) ? &cosine : &sine, pResult );
				pResult->SetPrecision( glWorkingPrecision );
				pResult->Round();

				return l_result;
			};

			// ...and a half decided by a bound the part lies under -- see
			//	lzi_polar_decide -- or, in a unit, r/2 a hair from an angle whose
			//	cosine or sine is a half by the angle's side of it -- see
			//	lzi_unit_half_decide...
			auto	tieCall = [ & ]( CLZNumber* pResult, long lPlaces, long ) -> bool
			{
				if ( (e_unit != AU_radians) &&
					 lzi_unit_half_decide( pResult, &theta, lPlaces, (l_part == 0) ? kAngleCosine : kAngleSine, &r ) )
					return true;

				return lzi_polar_decide( pResult, &r, &theta, (l_part == 1), (e_unit == AU_radians), lPlaces );
			};

			CLZNumber	retaken;

			if ( lzi_take_again_wider( fnCall, tieCall, l_base, l_saved, &retaken ) )
				parts[ l_part ].Assign( &retaken );
		}
	}

	lResult = lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 0 ], &r, &targetNumber->m_real, false,
																  b_short[ 0 ], b_short[ 0 ] ) );

	return lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 1 ], &r, &targetNumber->m_imag, false,
															   b_short[ 1 ], b_short[ 1 ] ) );
}


/////////////////////////////////////////////////////////////////
//	Complex numbers, stage 4 -- the logarithm, the exponential and the
//	square root, 2026-10-22

/*
	**The angle unit held at radians while a call lasts** -- for the
	logarithm's argument and the exponential's angle, radians in every
	unit, Rich's call, 2026-10-22: Elzed's ln and power answer so today,
	through std::complex, and exp(ln z) is z only so. A function called as
	the host's with this held is taken again wider and decided as the
	host's call is, in radians.
*/
class CLZRadiansHere
{
public:

	CLZRadiansHere()
		: m_eSaved( geAngleUnit )
	{
		geAngleUnit = AU_radians;
	}

	~CLZRadiansHere()
	{
		geAngleUnit = m_eSaved;
	}

private:

	CLZRadiansHere( const CLZRadiansHere& );
	CLZRadiansHere& operator=( const CLZRadiansHere& );

	AngleUnit	m_eSaved;
};

/*
	**A part worked under a guard, rounded at the width, and taken again
	wider while it is in doubt** -- the logarithm's real part and the
	exponential's two. pValue holds the part as fnPart( pValue, lWidth )
	worked it, lGuard places past the width; it is rounded at the working
	precision and, within kRetakeDoubt of a half, taken again by
	lzi_take_again_wider, as a function's answer is. None of the stage 4
	parts can lie exactly on a half -- each is transcendental -- so they
	decide no tie: one nearer a half than the retake reaches is rounded as
	the widest stands, as every function's is. A part with a bound to
	decide by passes tieCall -- see lzi_complex_trig_decide.
*/
template <class FnPart, class TieCall>
static void lzi_complex_part_rounded( FnPart fnPart, TieCall tieCall, long lBase, long lGuard, CLZNumber* pValue )
{
	long	l_saved = glWorkingPrecision;

	if ( lzi_rounded_clear_of_doubt( pValue, kRetakeDoubt, false ) )
		return;

	auto	fnCall = [ & ]( CLZNumber* pWider ) -> long
	{
		long	l_result = fnPart( pWider, glWorkingPrecision + lGuard );

		pWider->SetPrecision( glWorkingPrecision );
		pWider->Round();

		return l_result;
	};

	CLZNumber	retaken;

	if ( lzi_take_again_wider( fnCall, tieCall, lBase, l_saved, &retaken ) )
		pValue->Assign( &retaken );
}

// ...and with no tie to decide: the logarithm's and the exponential's parts...
template <class FnPart>
static void lzi_complex_part_rounded( FnPart fnPart, long lBase, long lGuard, CLZNumber* pValue )
{
	lzi_complex_part_rounded( fnPart, []( CLZNumber*, long, long ) -> bool { return false; }, lBase, lGuard, pValue );
}

/*
	**The natural logarithm of a complex value, ln |z| + i arg z** --
	stage 4, 2026-10-22: the principal value, its imaginary part in (-pi,
	pi], the cut on the negative real axis taken from the +0 side, so
	ln(-1) is pi i; in radians whatever the angle unit, Rich's
	call. Each part correctly rounded, half-even, in places.

	- **The real part is half ln(a^2 + b^2)**: the sum exact, the half
	  taken before the one rounding, under kComplexGuardDigits and taken
	  again wider in doubt -- not ln |z|, which would round the modulus
	  first. It is exact only where a^2 + b^2 is one, as at 0.6 + 0.8i: a
	  rational's logarithm is irrational but at one.
	- **The imaginary part is lz_arctangent2( b, a )**, called as the
	  host's with radians held, so taken again and decided as that function
	  is -- an exact, short nought on the positive real axis.
	- **On an axis the real part is a real logarithm's**: ln |a| beside
	  nought or pi, ln |b| beside +-pi/2. On the positive real axis the
	  answer is lz_log_e's to the character, with an exact nought beside it.

	ln 0 is refused as lz_log_e( 0 ) is, kUndefinedResult, the answer the
	point at infinity. A not-a-number, or the point at infinity, which has
	no direction, answers a not-a-number.
*/

long lzc_log_e( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( lzi_complex_has( leftOperand, kfQNaN ) || lzi_complex_has( leftOperand, kfInfinity ) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, lzi_complex_has( leftOperand, kfDiv0 ) );

		return kNoError;
	}

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	a;
	CLZNumber	b;

	a.Assign( &leftOperand->m_real );
	b.Assign( &leftOperand->m_imag );
	a.CheckForZero();
	b.CheckForZero();

	if ( a.IsZero() && b.IsZero() )
	{
		lzi_complex_flagged( p_target, l_base, kfInfinity, false );

		return kUndefinedResult;
	}

	// The argument, as the host's call in radians...
	{
		CLZRadiansHere	radians;

		lResult = lz_arctangent2( &b, &a, &p_target->m_imag );
	}

	// On an axis, ln |a| or ln |b|...
	if ( a.IsZero() || b.IsZero() )
	{
		CLZNumber	magnitude;

		magnitude.Assign( b.IsZero() ? &a : &b );
		magnitude.m_lSign = kPositive;

		return lzi_first_failure( lResult, lz_log_e( &magnitude, &p_target->m_real ) );
	}

	CLZNumber	aa;
	CLZNumber	bb;
	CLZNumber	sum;
	CLZNumber	one;
	CLZNumber	real;

	lzi_complex_product( &a, &a, &aa );
	lzi_complex_product( &b, &b, &bb );
	lz_add( &aa, &bb, &sum );
	one.SetBase( l_base, false );
	one.Assign( 1L );

	// ...exact on the unit circle...
	if ( lz_compare( &sum, &one ) == kEqualTo )
	{
		real.SetBase( l_base, false );
		real.Assign( 0L );

		return lzi_first_failure( lResult, lzi_answer_as_function( &real, &a, &p_target->m_real, false, true, true ) );
	}

	// ...and otherwise half ln( a^2 + b^2 ), under a guard of its own...
	{
		CLZNumber		scratch;
		CLZAnswerBase	work( &a, &scratch );

		work.Exact();

		auto	fnHalfLog = [ & ]( CLZNumber* pValue, long lWidth ) -> long
		{
			long	l_was = glWorkingPrecision;

			glWorkingPrecision = lWidth;

			CLZNumber	logarithm;
			CLZNumber	two;
			CLZNumber	half;

			two.SetBase( l_base, false );
			two.Assign( 2L );

			long	l_result = lz_log_e( &sum, &logarithm );

			lz_divide( &logarithm, &two, &half );
			pValue->Assign( &half );

			glWorkingPrecision = l_was;

			return l_result;
		};

		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base );

		lResult = lzi_first_failure( lResult, fnHalfLog( &real, glWorkingPrecision + l_guard ) );
		lzi_complex_part_rounded( fnHalfLog, l_base, l_guard, &real );
	}

	return lzi_first_failure( lResult, lzi_answer_as_function( &real, &a, &p_target->m_real, false ) );
}

/*
	**e to a complex power, e^x (cos y + i sin y)** -- stage 4, 2026-10-22,
	y in radians whatever the angle unit, Rich's call. Each part correctly
	rounded, half-even, in places: e^x, and the sine and cosine of y from
	one run as wide again as e^x's whole digits, the products exact and
	rounded once, a part in doubt taken again wider. Off the real axis
	neither part is rational (Lindemann-Weierstrass), so neither lies on a
	half.

	- **On the real axis, lz_exp_e's answer to the character**, with an
	  exact nought beside it; **on the imaginary axis, cos y and sin y as
	  lz_sine_cosine answers them**, in radians.
	- **e^x past the answer limit is refused as lz_exp_e refuses it**,
	  kTooLarge, the answer the point at infinity; under the width, both
	  parts are noughts rounded there.

	A not-a-number, or the point at infinity, answers a not-a-number.
*/

long lzc_exp_e( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( lzi_complex_has( leftOperand, kfQNaN ) || lzi_complex_has( leftOperand, kfInfinity ) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, lzi_complex_has( leftOperand, kfDiv0 ) );

		return kNoError;
	}

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	x;
	CLZNumber	y;
	CLZNumber	nought;

	x.Assign( &leftOperand->m_real );
	y.Assign( &leftOperand->m_imag );
	x.CheckForZero();
	y.CheckForZero();
	nought.SetBase( l_base, false );
	nought.Assign( 0L );

	// On the real axis, e^x, and the point at infinity where it is refused...
	if ( y.IsZero() )
	{
		lResult = lz_exp_e( &x, &p_target->m_real );

		if ( p_target->m_real.FlagIsOn( kfInfinity ) )
		{
			lzi_complex_flagged( p_target, l_base, kfInfinity, false );

			return lResult;
		}

		return lzi_first_failure( lResult, lzi_answer_as_function( &nought, &x, &p_target->m_imag, false, true, true ) );
	}

	// ...on the imaginary axis, cos y + i sin y, in radians...
	if ( x.IsZero() )
	{
		CLZRadiansHere	radians;

		return lz_sine_cosine( &y, &p_target->m_imag, &p_target->m_real );
	}

	// ...refused past the limit, nought under the width -- see
	//	lzi_answer_size...
	double	d_exponent = 0.0;

	if ( lzi_value_or_infinity( &x, &d_exponent ) )
	{
		ELZAnswerSize	e_size = lzi_answer_size( d_exponent, l_base );

		if ( e_size == AS_too_large )
		{
			lzi_complex_flagged( p_target, l_base, kfInfinity, false );

			return kTooLarge;
		}

		if ( e_size == AS_nought )
		{
			lResult = lzi_answer_as_function( &nought, &x, &p_target->m_real, false );

			return lzi_first_failure( lResult, lzi_answer_as_function( &nought, &x, &p_target->m_imag, false ) );
		}
	}

	CLZNumber	parts[ 2 ];			// real, imaginary

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( &x, &scratch );

		work.Exact();

		// Both parts to lWidth places: e^x, then the pair wider by e^x's
		//	whole digits, so the products keep the width...
		auto	fnParts = [ & ]( long lWidth, CLZNumber* pReal, CLZNumber* pImag ) -> long
		{
			long	l_was = glWorkingPrecision;

			glWorkingPrecision = lWidth;

			CLZNumber	modulus;
			long		l_result = lz_exp_e( &x, &modulus );
			long		l_lead = lzi_leading_exponent( &modulus );

			glWorkingPrecision = lWidth + ((l_lead >= 0) ? (l_lead + 1) : 0);

			CLZNumber		sine;
			CLZNumber		cosine;
			SLZExactPart	exact[ 2 ];

			l_result = lzi_first_failure( l_result, lzi_sine_cosine_work( &y, AU_radians, false, &sine, &cosine,
																		  &exact[ 1 ], &exact[ 0 ] ) );

			if ( pReal != nullptr )
				lzi_complex_product( &modulus, &cosine, pReal );

			if ( pImag != nullptr )
				lzi_complex_product( &modulus, &sine, pImag );

			glWorkingPrecision = l_was;

			return l_result;
		};

		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base );

		lResult = fnParts( glWorkingPrecision + l_guard, &parts[ 0 ], &parts[ 1 ] );

		lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth ) { return fnParts( lWidth, pValue, nullptr ); },
								  l_base, l_guard, &parts[ 0 ] );
		lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth ) { return fnParts( lWidth, nullptr, pValue ); },
								  l_base, l_guard, &parts[ 1 ] );
	}

	lResult = lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 0 ], &x, &p_target->m_real, false ) );

	return lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 1 ], &x, &p_target->m_imag, false ) );
}

/*
	**Which side of the half h a part of a square root lies, exactly** --
	-1 under it, 0 on it, 1 over it -- with A = |a|, B = |b|, S = a^2 + b^2
	and H2 = 2h, all exact, and m = sqrt S. For t = sqrt((A + m) / 2), t
	against h is 2m against V = H2^2 - 2A; for |b| / 2t it is 2B^2 - A H2^2
	against m H2^2. Each rational side's sign first, then both sides
	squared, so m is never worked out: the comparison is exact whether the
	part is rational or not, and equal only where it is exactly the half.
*/
static long lzi_root_part_side( bool bOther, CLZNumber* pA, CLZNumber* pB, CLZNumber* pS, CLZNumber* pH2 )
{
	CLZNumber	h2_squared;
	CLZNumber	rational;
	CLZNumber	left;
	CLZNumber	right;

	lzi_complex_product( pH2, pH2, &h2_squared );

	if ( !bOther )
	{
		CLZNumber	twice_a;
		CLZNumber	four;

		lz_add( pA, pA, &twice_a );
		lz_subtract( &h2_squared, &twice_a, &rational );
		rational.CheckForZero();

		// 2m is over nought, b being nought nowhere here...
		if ( rational.IsZero() || rational.IsNegative() )
			return 1;

		four.SetBase( pS->m_lBase, false );
		four.Assign( 4L );
		lzi_complex_product( pS, &four, &left );
		lzi_complex_product( &rational, &rational, &right );
	}
	else
	{
		CLZNumber	bb;
		CLZNumber	twice_bb;
		CLZNumber	a_h2;
		CLZNumber	h2_fourth;

		lzi_complex_product( pB, pB, &bb );
		lz_add( &bb, &bb, &twice_bb );
		lzi_complex_product( pA, &h2_squared, &a_h2 );
		lz_subtract( &twice_bb, &a_h2, &rational );
		rational.CheckForZero();

		// m H2^2 is over nought...
		if ( rational.IsZero() || rational.IsNegative() )
			return -1;

		lzi_complex_product( &rational, &rational, &right );
		lzi_complex_product( &h2_squared, &h2_squared, &h2_fourth );
		lzi_complex_product( pS, &h2_fourth, &left );

		// rational over m H2^2 is rational^2 over S H2^4...
		long	l_order = lz_compare( &right, &left );

		return (l_order == kGreaterThan) ? 1 : ((l_order == kEqualTo) ? 0 : -1);
	}

	long	l_order = lz_compare( &left, &right );

	return (l_order == kGreaterThan) ? 1 : ((l_order == kEqualTo) ? 0 : -1);
}

/*
	**A part of a square root in doubt, decided exactly** at lPlaces: y
	its magnitude cut there, u a unit there, the side of y + u/2 settled by
	lzi_root_part_side, and a part exactly on it taken to the even kept
	digit. pValue holds the part as worked, wider than lPlaces; the
	magnitude decided comes back in pResult.
*/
static void lzi_root_part_decide( bool bOther, CLZNumber* pValue, CLZNumber* pA, CLZNumber* pB, CLZNumber* pS,
								  long lPlaces, CLZNumber* pResult )
{
	long	l_base = pValue->m_lBase;
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pValue->m_lRadix > lPlaces) ? pValue->m_lRadix : lPlaces) + 4;

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	h2;

	y.Assign( pValue );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &h2 );
	lz_add( &h2, &unit );

	long	l_side = lzi_root_part_side( bOther, pA, pB, pS, &h2 );
	bool	b_up = (l_side > 0);

	// ...a part on the half to the even kept digit...
	if ( l_side == 0 )
	{
		CLZNumber	whole;
		CLZNumber	base;
		CLZNumber	digit;
		CLZNumber	two;
		CLZNumber	parity;

		whole.Assign( &y );
		lzi_scale_by_base_power( &whole, -lPlaces );
		base.SetBase( l_base, false );
		base.Assign( l_base );
		two.SetBase( l_base, false );
		two.Assign( 2L );
		lz_remainder( &whole, &base, &digit );
		lz_remainder( &digit, &two, &parity );
		parity.CheckForZero();
		b_up = !(parity.IsZero());
	}

	if ( b_up )
		lz_add( &y, &unit );

	glWorkingPrecision = l_saved;

	y.SetPrecision( lPlaces );
	pResult->Assign( &y );
}

/*
	**Whether a part of a square root, as rounded, is the root exactly** --
	squared back with no root taken: t is exact when m = 2t^2 - A is not
	below nought and m^2 is S, and |b| / 2t when B^2 - 2A o^2 is not below
	nought and its square is 4 o^4 S, o the part.
*/
static bool lzi_root_part_exact( bool bOther, CLZNumber* pPart, CLZNumber* pA, CLZNumber* pB, CLZNumber* pS )
{
	CLZNumber	part;
	CLZNumber	part_squared;
	CLZNumber	rational;
	CLZNumber	left;
	CLZNumber	right;

	part.Assign( pPart );
	part.m_lSign = kPositive;
	part.CheckForZero();

	if ( part.IsZero() )
		return false;

	lzi_complex_product( &part, &part, &part_squared );

	if ( !bOther )
	{
		CLZNumber	twice;

		lz_add( &part_squared, &part_squared, &twice );
		lz_subtract( &twice, pA, &rational );
		left.Assign( pS );
	}
	else
	{
		CLZNumber	bb;
		CLZNumber	a_o2;
		CLZNumber	twice;
		CLZNumber	o4;
		CLZNumber	four;

		lzi_complex_product( pB, pB, &bb );
		lzi_complex_product( pA, &part_squared, &a_o2 );
		lz_add( &a_o2, &a_o2, &twice );
		lz_subtract( &bb, &twice, &rational );
		lzi_complex_product( &part_squared, &part_squared, &o4 );
		four.SetBase( pS->m_lBase, false );
		four.Assign( 4L );
		lzi_complex_product( &o4, &four, &right );
		lzi_complex_product( &right, pS, &left );
	}

	rational.CheckForZero();

	if ( rational.IsNegative() )
		return false;

	lzi_complex_product( &rational, &rational, &right );

	return lz_compare( &left, &right ) == kEqualTo;
}

/*
	**The principal square root of a complex value** -- stage 4,
	2026-10-22: the root whose real part is not below nought, the cut on
	the negative real axis taken from the +0 side, so sqrt(-4) is 0 + 2i
	exactly. Each part correctly rounded, half-even, in
	places, and exact where the mathematics is: sqrt(3 + 4i) is 2 + i.

	**From |z| without cancelling**: with m = |z| and t = sqrt((|a| + m) /
	2), the root is t + i b / 2t for a at or above nought and |b| / 2t + i
	sign(b) t below -- |a| + m adds two magnitudes, and b / 2t divides by
	the larger part. Worked under kComplexGuardDigits and the places a
	small t would lose.

	**A part in doubt is decided exactly**, rational or not, as lz_root_2
	decides one -- see lzi_root_part_side -- and one exactly on the half
	goes to the even kept digit; **a part short of the width is squared
	back** before it is spelled short -- see lzi_root_part_exact.

	**On the real axis each part is lz_root_2's**: sqrt a beside an exact
	nought, or an exact nought beside sqrt -a. The point at infinity roots
	to itself; a not-a-number answers a not-a-number.
*/

long lzc_root_2( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, nullptr, p_target, kComplexUnary ) )
		return kNoError;

	long	l_base = leftOperand->m_real.m_lBase;

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	a;
	CLZNumber	b;
	CLZNumber	nought;

	a.Assign( &leftOperand->m_real );
	b.Assign( &leftOperand->m_imag );
	a.CheckForZero();
	b.CheckForZero();
	nought.SetBase( l_base, false );
	nought.Assign( 0L );

	// On the real axis, sqrt a or i sqrt -a...
	if ( b.IsZero() )
	{
		CLZNumber	magnitude;
		bool		b_negative = a.IsNegative();

		magnitude.Assign( &a );
		magnitude.m_lSign = a.IsZero() ? kZero : kPositive;

		lResult = lz_root_2( &magnitude, b_negative ? &p_target->m_imag : &p_target->m_real );

		return lzi_first_failure( lResult, lzi_answer_as_function( &nought, &a, b_negative ? &p_target->m_real
																						: &p_target->m_imag,
																	false, true, true ) );
	}

	CLZNumber	A;
	CLZNumber	B;
	CLZNumber	aa;
	CLZNumber	bb;
	CLZNumber	S;
	CLZNumber	t;
	CLZNumber	o;

	A.Assign( &a );
	A.m_lSign = a.IsZero() ? kZero : kPositive;
	B.Assign( &b );
	B.m_lSign = kPositive;
	lzi_complex_product( &a, &a, &aa );
	lzi_complex_product( &b, &b, &bb );
	lz_add( &aa, &bb, &S );

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( &a, &scratch );

		work.Exact();

		// t and |b| / 2t to lWidth places...
		auto	fnParts = [ & ]( long lWidth, CLZNumber* pT, CLZNumber* pO ) -> long
		{
			long	l_was = glWorkingPrecision;

			glWorkingPrecision = lWidth;

			CLZNumber	two;
			CLZNumber	modulus;
			CLZNumber	sum;
			CLZNumber	u;
			CLZNumber	root;
			CLZNumber	twice;
			CLZNumber	other;

			two.SetBase( l_base, false );
			two.Assign( 2L );

			long	l_result = lz_root_2( &S, &modulus );

			lz_add( &A, &modulus, &sum );
			lz_divide( &sum, &two, &u );
			l_result = lzi_first_failure( l_result, lz_root_2( &u, &root ) );
			lz_add( &root, &root, &twice );
			lz_divide( &B, &twice, &other );

			pT->Assign( &root );
			pO->Assign( &other );

			glWorkingPrecision = l_was;

			return l_result;
		};

		// A small z's modulus, sqrt S, has half of S's noughts, and is a
		//	nought at the width without them -- |b| / 2t then divides by
		//	one. The census found it: sqrt( -2.6e-36 i ) at a place in base
		//	thirteen answered an infinite imaginary part...
		long	l_lead = lzi_leading_exponent( &S );
		long	l_extra = (l_lead < 0) ? (((-l_lead) / 2) + 2) : 0;
		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base ) + l_extra;

		lResult = fnParts( glWorkingPrecision + l_guard, &t, &o );

		for ( bool b_other : { false, true } )
		{
			CLZNumber*	p_part = b_other ? &o : &t;
			CLZNumber	worked;

			worked.Assign( p_part );

			if ( !lzi_rounded_clear_of_doubt( p_part, kRetakeDoubt, false ) )
				lzi_root_part_decide( b_other, &worked, &A, &B, &S, glWorkingPrecision, p_part );
		}
	}

	// The parts in place: t and |b| / 2t, by the sign of a; b's sign on the
	//	imaginary part; short where squared back exact...
	bool		b_left = a.IsNegative();
	CLZNumber*	p_real = b_left ? &o : &t;
	CLZNumber*	p_imag = b_left ? &t : &o;
	bool		b_realExact = lzi_short_of_width( p_real ) && lzi_root_part_exact( b_left, p_real, &A, &B, &S );
	bool		b_imagExact = lzi_short_of_width( p_imag ) && lzi_root_part_exact( !b_left, p_imag, &A, &B, &S );

	p_real->m_lSign = kPositive;
	p_imag->m_lSign = b.IsNegative() ? kNegative : kPositive;
	p_real->CheckForZero();
	p_imag->CheckForZero();

	lResult = lzi_first_failure( lResult, lzi_answer_as_function( p_real, &a, &p_target->m_real, false,
																  b_realExact, b_realExact ) );

	return lzi_first_failure( lResult, lzi_answer_as_function( p_imag, &a, &p_target->m_imag, false,
															   b_imagExact, b_imagExact ) );
}


/////////////////////////////////////////////////////////////////
//	Complex numbers, stage 5a -- the circular and hyperbolic functions,
//	2026-10-22

// The twelve by index: the circular six in kAngle order, the hyperbolic
//	six after them in the same order...
const long	kComplexHyperbolic = 6;

/*
	**Every part of the twelve is a sign and a product, over a denominator
	for the quotients.** With s and c the sine and cosine of the circular
	argument and sh, ch the sinh and cosh of the hyperbolic one -- for sin z
	the real part and the imaginary, for sinh z the imaginary part and the
	real -- sin z is s ch + i c sh. **The quotients are divided through by
	ch^2**, so t = tanh and u = sech stand for sh and ch, both bounded
	however large the argument: tan z is (s c u^2 + i t) / (c^2 u^2 + t^2),
	where the textbook form needs cosh 2y, past the answer limit long
	before tan z is anything but i. The two denominators, Dc = c^2 u^2 +
	t^2 and Ds = s^2 u^2 + t^2, are |cos z|^2 and |sin z|^2 over ch^2: sums
	of squares, so nothing cancels near a pole.
*/
enum ELZTrigFactor { TF_none = 0, TF_s, TF_c, TF_sh, TF_ch, TF_t, TF_u, TF_count };

/*
	**A quotient starts two places wider than the guard**: its denominator,
	from one over the base up to two, magnifies the error by up to the base
	and the answer by up to its root -- b^1.5, two places in every base.
	Measured: without them the worst quotient error at the guard is 1.4e-4
	of a unit, seven times inside kRetakeDoubt; they had been there by
	accident until 2026-10-22, a denominator of 0.98 read as one nought
	and the whole run taken again six places wider.
*/
const long	kTrigQuotientPlaces = 2;

const long	kTrigNoDenominator = 0;
const long	kTrigDc = 1;
const long	kTrigDs = 2;

struct SLZTrigPart
{
	long			m_lSign;
	ELZTrigFactor	m_eFactors[ 4 ];
};

struct SLZTrigFunction
{
	SLZTrigPart		m_parts[ 2 ];			// real, imaginary
	long			m_lDenominator;
};

static const SLZTrigFunction	sl_trigFunctions[ 12 ] =
{
	{ { { 1, { TF_s, TF_ch } },				{ 1, { TF_c, TF_sh } } },				kTrigNoDenominator },	// sin
	{ { { 1, { TF_c, TF_ch } },				{ -1, { TF_s, TF_sh } } },				kTrigNoDenominator },	// cos
	{ { { 1, { TF_s, TF_c, TF_u, TF_u } },	{ 1, { TF_t } } },						kTrigDc },				// tan
	{ { { 1, { TF_s, TF_c, TF_u, TF_u } },	{ -1, { TF_t } } },						kTrigDs },				// cot
	{ { { 1, { TF_c, TF_u } },				{ 1, { TF_s, TF_t, TF_u } } },			kTrigDc },				// sec
	{ { { 1, { TF_s, TF_u } },				{ -1, { TF_c, TF_t, TF_u } } },			kTrigDs },				// csc
	{ { { 1, { TF_sh, TF_c } },				{ 1, { TF_ch, TF_s } } },				kTrigNoDenominator },	// sinh
	{ { { 1, { TF_ch, TF_c } },				{ 1, { TF_sh, TF_s } } },				kTrigNoDenominator },	// cosh
	{ { { 1, { TF_t } },					{ 1, { TF_s, TF_c, TF_u, TF_u } } },	kTrigDc },				// tanh
	{ { { 1, { TF_t } },					{ -1, { TF_s, TF_c, TF_u, TF_u } } },	kTrigDs },				// coth
	{ { { 1, { TF_c, TF_u } },				{ -1, { TF_s, TF_t, TF_u } } },			kTrigDc },				// sech
	{ { { 1, { TF_c, TF_t, TF_u } },		{ -1, { TF_s, TF_u } } },				kTrigDs },				// csch
};

// The real functions, for the real axis...
static const LZUnaryFunction	sl_realTrig[ 12 ] =
{
	lz_sine, lz_cosine, lz_tangent, lz_cotangent, lz_secant, lz_cosecant,
	lz_hypersine, lz_hypercosine, lz_hypertangent, lz_hypercotangent, lz_hypersecant, lz_hypercosecant,
};

/*
	**On the imaginary axis each is i, -i or one times a real function of
	the other family**, in radians -- sin iy = i sinh y, cot iy = -i coth y,
	cosh iy = cos y: which part it is, the function, and its sign.
*/
struct SLZTrigAxis
{
	long			m_lPart;
	LZUnaryFunction	m_fn;
	long			m_lSign;
};

static const SLZTrigAxis	sl_trigOnImaginaryAxis[ 12 ] =
{
	{ 1, lz_hypersine, 1 }, { 0, lz_hypercosine, 1 }, { 1, lz_hypertangent, 1 },
	{ 1, lz_hypercotangent, -1 }, { 0, lz_hypersecant, 1 }, { 1, lz_hypercosecant, -1 },
	{ 1, lz_sine, 1 }, { 0, lz_cosine, 1 }, { 1, lz_tangent, 1 },
	{ 1, lz_cotangent, -1 }, { 0, lz_secant, 1 }, { 1, lz_cosecant, -1 },
};

/*
	**Each part near nought as a polynomial over a denominator, and a bound
	on the rest** -- for the ties a small z makes by construction. There
	sin z is z + z^3/6 ..., cot z is 1/z - z/3 ..., cos z is 1 - z^2/2 ...,
	so a part is a component of z, or of 1/z, or one, or the product xy,
	moved by terms at twice its own depth and more: where that leading term
	lies on a half of the answer's places -- tan(2.2e-9 - 2^-37 i) in base
	two at 36 places has y on one -- only those terms move the part off
	it, past any width the retake reaches. The census found two.

	**S part = c poly(x^2, y^2) / N + R, |R| at most K |c| (r^2)^m / N**,
	for |z| under one over the base: S one or r^2 = x^2 + y^2, c one, x, y
	or xy, the coefficients over 1, x^2, y^2, x^4, x^2 y^2, y^4, x^6, x^4
	y^2, x^2 y^4 and y^6. **Generated, not typed**: from the series in
	exact Fractions, the remainder bounded by each tail at |z| = 1/2 with
	|Re z^k / x|, |Im z^k / y| at most k r^(k-1), |Im z^k / xy| at most
	k r^(k-2) and |Re z^k| at most r^k, and the bound checked against
	Decimal at random small z -- the worst remainder 0.78 of its bound.
	See lzi_complex_trig_near_decide.
*/
enum ELZNearFactor { TN_one = 0, TN_x, TN_y, TN_xy };

struct SLZTrigNearNought
{
	bool			m_bByModulusSquared;		// S is r^2
	ELZNearFactor	m_eFactor;					// c
	long			m_lScale;					// N
	long			m_lCoefficients[ 10 ];
	long			m_lRemainder;				// K
	long			m_lRemainderPower;			// m
};

static const SLZTrigNearNought	sl_trigNearNought[ 12 ][ 2 ] =
{
	{ { false, TN_x, 120, { 120, -20, 60, 1, -10, 5, 0, 0, 0, 0 }, 1, 3 },
	  { false, TN_y, 120, { 120, -60, 20, 5, -10, 1, 0, 0, 0, 0 }, 1, 3 } },		// sin
	{ { false, TN_one, 720, { 720, -360, 360, 30, -180, 30, -1, 15, -15, 1 }, 1, 4 },
	  { false, TN_xy, 360, { -360, 60, -60, -3, 10, -3, 0, 0, 0, 0 }, 1, 3 } },		// cos
	{ { false, TN_x, 15, { 15, 5, -15, 2, -20, 10, 0, 0, 0, 0 }, 9, 3 },
	  { false, TN_y, 15, { 15, 15, -5, 10, -20, 2, 0, 0, 0, 0 }, 9, 3 } },		// tan
	{ { true, TN_x, 45, { 45, -15, -15, -1, 2, 3, 0, 0, 0, 0 }, 1, 3 },
	  { true, TN_y, 45, { -45, -15, -15, -3, -2, 1, 0, 0, 0, 0 }, 1, 3 } },		// cot
	{ { false, TN_one, 720, { 720, 360, -360, 150, -900, 150, 61, -915, 915, -61 }, 35, 4 },
	  { false, TN_xy, 360, { 360, 300, -300, 183, -610, 183, 0, 0, 0, 0 }, 142, 3 } },		// sec
	{ { true, TN_x, 360, { 360, 60, 60, 7, -14, -21, 0, 0, 0, 0 }, 5, 3 },
	  { true, TN_y, 360, { -360, 60, 60, 21, 14, -7, 0, 0, 0, 0 }, 5, 3 } },		// csc
	{ { false, TN_x, 120, { 120, 20, -60, 1, -10, 5, 0, 0, 0, 0 }, 1, 3 },
	  { false, TN_y, 120, { 120, 60, -20, 5, -10, 1, 0, 0, 0, 0 }, 1, 3 } },		// sinh
	{ { false, TN_one, 720, { 720, 360, -360, 30, -180, 30, 1, -15, 15, -1 }, 1, 4 },
	  { false, TN_xy, 360, { 360, 60, -60, 3, -10, 3, 0, 0, 0, 0 }, 1, 3 } },		// cosh
	{ { false, TN_x, 15, { 15, -5, 15, 2, -20, 10, 0, 0, 0, 0 }, 9, 3 },
	  { false, TN_y, 15, { 15, -15, 5, 10, -20, 2, 0, 0, 0, 0 }, 9, 3 } },		// tanh
	{ { true, TN_x, 45, { 45, 15, 15, -1, 2, 3, 0, 0, 0, 0 }, 1, 3 },
	  { true, TN_y, 45, { -45, 15, 15, -3, -2, 1, 0, 0, 0, 0 }, 1, 3 } },		// coth
	{ { false, TN_one, 720, { 720, -360, 360, 150, -900, 150, -61, 915, -915, 61 }, 35, 4 },
	  { false, TN_xy, 360, { -360, 300, -300, -183, 610, -183, 0, 0, 0, 0 }, 142, 3 } },		// sech
	{ { true, TN_x, 360, { 360, -60, -60, 7, -14, -21, 0, 0, 0, 0 }, 5, 3 },
	  { true, TN_y, 360, { -360, -60, -60, 21, 14, -7, 0, 0, 0, 0 }, 5, 3 } },		// csch
};

// The values the parts are made of, at one width, and which of the sine
//	and cosine are exact...
struct SLZTrigValues
{
	CLZNumber		m_values[ TF_count ];
	SLZExactPart	m_sineExact;
	SLZExactPart	m_cosineExact;
};

/*
	**s and c of alpha, and sh, ch or t, u of beta**, at the working
	precision as it stands. An exact s or c is written out at the width --
	a half, in an odd base, does not end there -- and marked. For a
	quotient t and u come from sinh and cosh, and **bFar** -- e^-|beta| past
	the width -- from F = e^-|beta|, as (1 - F^2) / (1 + F^2) and 2F / (1 +
	F^2), so an argument whose cosh is past the answer limit still answers.
*/
static long lzi_complex_trig_values( CLZNumber* pAlpha, AngleUnit eUnit, CLZNumber* pBeta, bool bQuotient, bool bFar,
									 SLZTrigValues* pValues )
{
	long	l_base = pAlpha->m_lBase;
	long	lResult = lzi_sine_cosine_work( pAlpha, eUnit, true, &pValues->m_values[ TF_s ], &pValues->m_values[ TF_c ],
											&pValues->m_sineExact, &pValues->m_cosineExact );

	for ( ELZTrigFactor e_factor : { TF_s, TF_c } )
	{
		SLZExactPart&	exact = (e_factor == TF_s) ? pValues->m_sineExact : pValues->m_cosineExact;

		if ( !exact.m_bExact )
			continue;

		CLZNumber	numerator;
		CLZNumber	denominator;

		numerator.SetBase( l_base, false );
		numerator.Assign( exact.m_lNumerator );
		denominator.SetBase( l_base, false );
		denominator.Assign( exact.m_lDenominator );
		lz_divide( &numerator, &denominator, &pValues->m_values[ e_factor ] );
	}

	CLZNumber	one;

	one.SetBase( l_base, false );
	one.Assign( 1L );

	if ( !bFar )
	{
		lResult = lzi_first_failure( lResult, lzi_sinh_cosh( pBeta, &pValues->m_values[ TF_sh ], &pValues->m_values[ TF_ch ] ) );

		if ( bQuotient )
		{
			lz_divide( &pValues->m_values[ TF_sh ], &pValues->m_values[ TF_ch ], &pValues->m_values[ TF_t ] );
			lz_divide( &one, &pValues->m_values[ TF_ch ], &pValues->m_values[ TF_u ] );
		}

		return lResult;
	}

	CLZNumber	minus;
	CLZNumber	f;
	CLZNumber	ff;
	CLZNumber	sum;
	CLZNumber	difference;
	CLZNumber	twice;

	minus.Assign( pBeta );
	minus.m_lSign = kNegative;
	lResult = lzi_first_failure( lResult, lz_exp_e( &minus, &f ) );
	lz_multiply( &f, &f, &ff );
	lz_add( &one, &ff, &sum );
	lz_subtract( &one, &ff, &difference );
	lz_divide( &difference, &sum, &pValues->m_values[ TF_t ] );
	lz_add( &f, &f, &twice );
	lz_divide( &twice, &sum, &pValues->m_values[ TF_u ] );

	if ( pBeta->IsNegative() )
		lz_negate( &pValues->m_values[ TF_t ] );

	return lResult;
}

// Whether a part is exactly nought -- a sine or cosine in it exactly
//	nought, which only a quarter turn in a unit gives...
static bool lzi_trig_part_nought( const SLZTrigPart& part, const SLZTrigValues& values )
{
	for ( ELZTrigFactor e_factor : part.m_eFactors )
	{
		if ( ((e_factor == TF_s) && values.m_sineExact.m_bExact && (values.m_sineExact.m_lNumerator == 0)) ||
			 ((e_factor == TF_c) && values.m_cosineExact.m_bExact && (values.m_cosineExact.m_lNumerator == 0)) )
			return true;
	}

	return false;
}

// Whether a part holds a factor, and the part's sine or cosine...
static bool lzi_trig_part_has( const SLZTrigPart& part, ELZTrigFactor eFactor )
{
	for ( ELZTrigFactor e_factor : part.m_eFactors )
	{
		if ( e_factor == eFactor )
			return true;
	}

	return false;
}

/*
	**A part on a half, decided by the bound it lies over** -- the real
	part of sin z or cos z in a unit, e cosh beta with e an exact sine or
	cosine, one or a half. bOver false is the bound it lies under, for a
	power's e M with M under one -- see lzi_complex_power_general. cosh beta is over one wherever beta is not
	nought, so the part's magnitude is over |e|, and where |e| is at or
	over the half in doubt the magnitude rounds up. **It is the tie a unit
	makes by construction**: sin(30 + iy) degrees is (1 + y^2 pi^2 /
	64800 ...) / 2, a half in an odd base and at no places in any, short of
	it by a depth no retake reaches once y is small. The polar counterpart
	decides by a bound it lies under -- see lzi_polar_decide. An x a hair
	from such an angle is lzi_complex_trig_half_decide's.
*/
static bool lzi_complex_trig_decide( CLZNumber* pResult, const SLZExactPart& exact, long lPlaces, bool bOver = true )
{
	if ( (lPlaces < 0) || !exact.m_bExact || (exact.m_lNumerator == 0) )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pResult->m_lRadix > lPlaces) ? pResult->m_lRadix : lPlaces) + 4;

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	twice_mid;
	CLZNumber	twice_bound;
	CLZNumber	denominator;
	CLZNumber	scaled;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &twice_mid );
	lz_add( &twice_mid, &unit );

	// 2|e| against 2y + u, with |e| = |n| / d: 2|n| against (2y + u) d...
	twice_bound.SetBase( l_base, false );
	twice_bound.Assign( 2 * ((exact.m_lNumerator < 0) ? -exact.m_lNumerator : exact.m_lNumerator) );
	denominator.SetBase( l_base, false );
	denominator.Assign( exact.m_lDenominator );
	lzi_complex_product( &twice_mid, &denominator, &scaled );

	// Over |e| and |e| at or over the half: up. Under |e| and |e| at or
	//	under it: down, to y...
	long	l_order = lz_compare( &twice_bound, &scaled );
	bool	b_decided = bOver ? (l_order != kLessThan) : (l_order != kGreaterThan);

	if ( b_decided )
	{
		if ( bOver )
			lz_add( &y, &unit );

		if ( pResult->IsNegative() )
			lz_negate( &y );

		y.CheckForZero();
		pResult->Assign( &y );
	}

	glWorkingPrecision = l_saved;

	return b_decided;
}

/*
	**The real part of sin z or cos z in a unit, x a hair from an angle
	whose sine or cosine is a half and y small** -- since 2026-10-23, found
	with the real sine's and cosine's halves (see lzi_unit_half_decide).
	sin(30 + dx + iy) degrees is sin(30 + dx) cosh y, a hair from the half
	1/2 in both coordinates at once, and lzi_complex_trig_decide covers only
	an exact x: sin(30 + 7^-60 + 7^-70 i) at ten places in base seven
	answered ...333 for ...334.

	With a the angle, d and b the hair dx and y in radians, s = f(a) = +-1/2
	and c = f'(a) = +-sqrt(3)/2, the part less s is

		c d + s (b^2 - d^2) / 2 + R,

	|R| at most (d^2 + b^2)^2 / 8 + 0.44 |d| b^2 + 0.15 |d|^3 while |d| and
	|b| are at most a hundredth -- Taylor's remainder for cos d cosh b - 1,
	the cross term d^2 b^2 / 4 among it, and for sin d cosh b - d. So where
	c d + s (b^2 - d^2) / 2, worked wide enough, lies further from nought
	than that bound, its sign is the part's side of the half; where it does
	not, the two leading terms nearly cancel, and it declines.
*/
static bool lzi_complex_trig_half_decide( CLZNumber* pResult, long lFunction, CLZNumber* pX, CLZNumber* pY, long lPlaces )
{
	long	l_turn = lzi_turn_of( geAngleUnit );

	if ( (lPlaces < 0) || (l_turn == 0) )
		return false;

	long	l_base = pResult->m_lBase;

	if ( (pX->m_lBase != l_base) || (pY->m_lBase != l_base) || pX->FlagIsOn( kfInfinity ) || pX->FlagIsOn( kfQNaN ) ||
		 pY->FlagIsOn( kfInfinity ) || pY->FlagIsOn( kfQNaN ) )
		return false;

	CLZNumber	y_part;

	y_part.Assign( pY );
	y_part.CheckForZero();

	if ( y_part.IsZero() )
		return false;

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pX->m_lRadix > lPlaces) ? pX->m_lRadix : lPlaces) + 4;

	if ( y_part.m_lRadix + 4 > glWorkingPrecision )
		glWorkingPrecision = y_part.m_lRadix + 4;

	CLZNumber	y;
	long		l_value = 0;
	long		l_slope = 0;
	CLZNumber	twelve_dx;

	if ( !lzi_answer_is_half( pResult, lPlaces, nullptr, &y ) ||
		 !lzi_near_twelfth( pX, lFunction, l_turn, &l_value, &l_slope, &twelve_dx ) )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	long	l_sign_h = pResult->IsNegative() ? -1 : 1;

	twelve_dx.CheckForZero();

	if ( (l_value != l_sign_h) || twelve_dx.IsZero() )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// Wide enough to hold the smaller of d and b^2 to thirty places of its
	//	own: twice the larger lead's noughts, and forty...
	long	l_lead_d = lzi_leading_exponent( &twelve_dx );
	long	l_lead_b = lzi_leading_exponent( &y_part );
	long	l_noughts = 0;

	if ( -l_lead_d > l_noughts )
		l_noughts = -l_lead_d;

	if ( -l_lead_b > l_noughts )
		l_noughts = -l_lead_b;

	long	l_wide = 2 * l_noughts + 40;

	if ( l_wide < glWorkingPrecision )
		l_wide = glWorkingPrecision;

	glWorkingPrecision = l_wide;

	// d = 12 dx (2 pi / turn) / 12 and b = y (2 pi / turn), in radians...
	CLZNumber	twoPi( kTwoPi );
	CLZNumber	per_unit;
	CLZNumber	turn;
	CLZNumber	twelve;
	CLZNumber	d;
	CLZNumber	b;
	CLZNumber	work;

	lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, l_wide + 4 );
	turn.SetBase( l_base, false );
	turn.Assign( l_turn );
	lz_divide( &twoPi, &turn, &per_unit );
	twelve.SetBase( l_base, false );
	twelve.Assign( 12L );
	lz_multiply( &twelve_dx, &per_unit, &work );
	lz_divide( &work, &twelve, &d );
	lz_multiply( &y_part, &per_unit, &b );

	// Both a hair: past a hundredth the bound below does not hold...
	CLZNumber	hundredth;
	CLZNumber	d_abs;
	CLZNumber	b_abs;

	hundredth.SetBase( l_base, false );
	hundredth.Assign( 1L );
	work.Assign( 100L );
	lz_divide( &hundredth, &work, &hundredth );
	d_abs.Assign( &d );
	d_abs.m_lSign = d.IsZero() ? kZero : kPositive;
	b_abs.Assign( &b );
	b_abs.m_lSign = b.IsZero() ? kZero : kPositive;

	if ( (lz_compare( &d_abs, &hundredth ) == kGreaterThan) || (lz_compare( &b_abs, &hundredth ) == kGreaterThan) )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// c d + s (b^2 - d^2) / 2, c = slope sqrt(3) / 2 and s = value / 2...
	CLZNumber	three;
	CLZNumber	root3;
	CLZNumber	cd;
	CLZNumber	dd;
	CLZNumber	bb;
	CLZNumber	second;
	CLZNumber	sum;

	three.SetBase( l_base, false );
	three.Assign( 3L );
	lz_root_2( &three, &root3 );
	lz_multiply( &root3, &d, &work );
	lz_divide( &work, &twelve, &cd );					// sqrt(3) d / 12...
	lzi_times_whole( &cd, 6, &cd );						// ...times six, sqrt(3) d / 2

	if ( l_slope < 0 )
		lz_negate( &cd );

	lz_multiply( &d, &d, &dd );
	lz_multiply( &b, &b, &bb );
	lz_subtract( &bb, &dd, &second );
	lz_divide( &second, &twelve, &work );				// (b^2 - d^2) / 12...
	lzi_times_whole( &work, 3, &second );				// ...times three, (b^2 - d^2) / 4

	if ( l_value < 0 )
		lz_negate( &second );

	lz_add( &cd, &second, &sum );
	sum.CheckForZero();

	// ...against (d^2 + b^2)^2 / 8 + 0.44 |d| b^2 + 0.15 |d|^3, and a unit of
	//	the work's last place for its rounding...
	CLZNumber	bound;
	CLZNumber	term;
	CLZNumber	square_sum;
	CLZNumber	scale;

	lz_add( &dd, &bb, &square_sum );
	lz_multiply( &square_sum, &square_sum, &term );
	scale.SetBase( l_base, false );
	scale.Assign( 8L );
	lz_divide( &term, &scale, &bound );
	lz_multiply( &d_abs, &bb, &term );
	scale.Assign( 44L );
	lz_multiply( &term, &scale, &work );
	scale.Assign( 100L );
	lz_divide( &work, &scale, &term );
	lz_add( &bound, &term, &bound );
	lz_multiply( &d_abs, &dd, &term );
	scale.Assign( 15L );
	lz_multiply( &term, &scale, &work );
	scale.Assign( 100L );
	lz_divide( &work, &scale, &term );
	lz_add( &bound, &term, &bound );

	CLZNumber	ulp;

	ulp.SetBase( l_base, false );
	ulp.Assign( 16L );
	lzi_scale_by_base_power( &ulp, l_wide - 2 );
	lz_add( &bound, &ulp, &bound );

	CLZNumber	sum_abs;

	sum_abs.Assign( &sum );
	sum_abs.m_lSign = sum.IsZero() ? kZero : kPositive;

	if ( lz_compare( &sum_abs, &bound ) != kGreaterThan )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// The magnitude goes up where the part lies beyond the half...
	bool		b_up = ((sum.IsPositive() ? 1 : -1) * l_sign_h > 0);
	CLZNumber	unit;

	glWorkingPrecision = lPlaces + 2;
	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	if ( b_up )
		lz_add( &y, &unit );

	y.m_lSign = pResult->m_lSign;
	y.CheckForZero();
	pResult->Assign( &y );
	pResult->SetPrecision( lPlaces );

	glWorkingPrecision = l_saved;

	return true;
}

/*
	**A part of a small z on a half, decided by its polynomial** -- see
	sl_trigNearNought. With y the part's magnitude cut at its places and u a
	unit there, 2 N S |part| against N S (2y + u) is 2 |W| against it, W =
	c poly worked exactly, give or take twice the remainder's bound: above
	where the difference outweighs it, below where it falls short of minus
	it, and otherwise declined, the retake going on as before. Nothing here
	rounds: x and y are as given, so every product and sum is exact.
*/
static bool lzi_complex_trig_near_decide( CLZNumber* pResult, const SLZTrigNearNought& near, CLZNumber* pX, CLZNumber* pY,
										  long lPlaces )
{
	long	l_base = pResult->m_lBase;

	if ( (lPlaces < 0) || (pX->m_lBase != l_base) || (pY->m_lBase != l_base) )
		return false;

	long	l_saved = glWorkingPrecision;
	long	l_digits = (pX->m_lRadix > pY->m_lRadix) ? pX->m_lRadix : pY->m_lRadix;
	long	l_places = (pResult->m_lRadix > lPlaces) ? pResult->m_lRadix : lPlaces;

	// Ten of x's or y's places at most in any product, and the answer's...
	glWorkingPrecision = (10 * l_digits) + l_places + 8;

	CLZNumber	one;
	CLZNumber	x2;
	CLZNumber	y2;
	CLZNumber	r2;
	CLZNumber	base_squared;
	CLZNumber	scaled;

	one.SetBase( l_base, false );
	one.Assign( 1L );
	lzi_complex_product( pX, pX, &x2 );
	lzi_complex_product( pY, pY, &y2 );
	lz_add( &x2, &y2, &r2 );

	// |z| under one over the base, where the bounds hold...
	base_squared.SetBase( l_base, false );
	base_squared.Assign( l_base * l_base );
	lzi_complex_product( &r2, &base_squared, &scaled );

	if ( lz_compare( &scaled, &one ) != kLessThan )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// 1, x^2, y^2, x^4, x^2 y^2, y^4, x^6, x^4 y^2, x^2 y^4, y^6...
	CLZNumber	monomials[ 10 ];

	monomials[ 0 ].Assign( &one );
	monomials[ 1 ].Assign( &x2 );
	monomials[ 2 ].Assign( &y2 );
	lzi_complex_product( &x2, &x2, &monomials[ 3 ] );
	lzi_complex_product( &x2, &y2, &monomials[ 4 ] );
	lzi_complex_product( &y2, &y2, &monomials[ 5 ] );
	lzi_complex_product( &monomials[ 3 ], &x2, &monomials[ 6 ] );
	lzi_complex_product( &monomials[ 3 ], &y2, &monomials[ 7 ] );
	lzi_complex_product( &x2, &monomials[ 5 ], &monomials[ 8 ] );
	lzi_complex_product( &monomials[ 5 ], &y2, &monomials[ 9 ] );

	CLZNumber	poly;

	poly.SetBase( l_base, false );
	poly.Assign( 0L );

	for ( long i = 0; i < 10; i++ )
	{
		if ( near.m_lCoefficients[ i ] == 0 )
			continue;

		CLZNumber	coefficient;
		CLZNumber	term;
		CLZNumber	sum;

		coefficient.SetBase( l_base, false );
		coefficient.Assign( near.m_lCoefficients[ i ] );
		lzi_complex_product( &coefficient, &monomials[ i ], &term );
		lz_add( &poly, &term, &sum );
		poly.Assign( &sum );
	}

	// c, W = c poly, and N S...
	CLZNumber	c;
	CLZNumber	w;
	CLZNumber	scale;
	CLZNumber	ns;

	if ( near.m_eFactor == TN_one )
		c.Assign( &one );
	else if ( near.m_eFactor == TN_x )
		c.Assign( pX );
	else if ( near.m_eFactor == TN_y )
		c.Assign( pY );
	else
		lzi_complex_product( pX, pY, &c );

	lzi_complex_product( &c, &poly, &w );

	scale.SetBase( l_base, false );
	scale.Assign( near.m_lScale );

	if ( near.m_bByModulusSquared )
		lzi_complex_product( &scale, &r2, &ns );
	else
		ns.Assign( &scale );

	// ...and K |c| (r^2)^m, the remainder's bound...
	CLZNumber	bound;
	CLZNumber	power;

	power.Assign( &one );

	for ( long i = 0; i < near.m_lRemainderPower; i++ )
	{
		CLZNumber	next;

		lzi_complex_product( &power, &r2, &next );
		power.Assign( &next );
	}

	{
		CLZNumber	k;
		CLZNumber	kc;

		c.m_lSign = c.IsZero() ? kZero : kPositive;
		k.SetBase( l_base, false );
		k.Assign( near.m_lRemainder );
		lzi_complex_product( &k, &c, &kc );
		lzi_complex_product( &kc, &power, &bound );
	}

	// |W| past the bound, so the part has W's sign...
	w.m_lSign = w.IsZero() ? kZero : kPositive;
	w.CheckForZero();

	if ( lz_compare( &w, &bound ) != kGreaterThan )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	twice_mid;
	CLZNumber	scaled_mid;
	CLZNumber	twice_w;
	CLZNumber	gap;
	CLZNumber	twice_bound;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	lz_add( &y, &y, &twice_mid );
	lz_add( &twice_mid, &unit );
	lzi_complex_product( &ns, &twice_mid, &scaled_mid );
	lz_add( &w, &w, &twice_w );
	lz_subtract( &twice_w, &scaled_mid, &gap );
	gap.CheckForZero();
	lz_add( &bound, &bound, &twice_bound );

	CLZNumber	size;

	size.Assign( &gap );
	size.m_lSign = size.IsZero() ? kZero : kPositive;

	if ( lz_compare( &size, &twice_bound ) != kGreaterThan )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	if ( !gap.IsNegative() )
		lz_add( &y, &unit );

	if ( pResult->IsNegative() )
		lz_negate( &y );

	pResult->Assign( &y );
	glWorkingPrecision = l_saved;

	return true;
}

/*
	**The twelve complex circular and hyperbolic functions** -- stage 5a,
	2026-10-22. The circular six read the whole complex angle in the host's
	unit, as Elzed's do: z degrees is z pi / 180, both parts.
	The hyperbolic six are in radians whatever the unit, as Elzed's are.
	Each part correctly rounded, half-even, in places.

	- **On the real axis, the real function's answer to the character** --
	  its rule for a small quotient included, Rich's call -- beside an exact
	  nought; a pole there is the point at infinity, with the real
	  function's code. **On the imaginary axis in radians**, the real
	  function of the other family, times i, -i or one, beside an exact
	  nought -- see sl_trigOnImaginaryAxis.
	- **Off the axes both parts in places**, from lzi_complex_trig_values
	  by sl_trigFunctions, worked under kComplexGuardDigits: wider by cosh's
	  whole digits for sin, cos, sinh and cosh, and near a pole by twice
	  the denominator's noughts, and taken again wider in doubt. **A part
	  is exactly nought** where its sine or cosine is, which only a quarter
	  turn in a unit gives: cos(90 + iy) degrees is exactly -i sinh(y pi /
	  180). sin, cos, sinh and cosh past the answer limit are refused,
	  kTooLarge, the point at infinity, and so is a quotient's near a pole.

	A not-a-number, or the point at infinity, answers a not-a-number.
*/

static long lzi_complex_trig( long lFunction, CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( lzi_complex_has( leftOperand, kfQNaN ) || lzi_complex_has( leftOperand, kfInfinity ) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, lzi_complex_has( leftOperand, kfDiv0 ) );

		return kNoError;
	}

	bool					b_hyperbolic = (lFunction >= kComplexHyperbolic);
	const SLZTrigFunction&	function = sl_trigFunctions[ lFunction ];
	bool					b_quotient = (function.m_lDenominator != kTrigNoDenominator);

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	x;
	CLZNumber	y;
	CLZNumber	nought;

	x.Assign( &leftOperand->m_real );
	y.Assign( &leftOperand->m_imag );
	x.CheckForZero();
	y.CheckForZero();
	nought.SetBase( l_base, false );
	nought.Assign( 0L );

	// The host's unit, read as an answer guard reads it -- see
	//	lzi_angle_unit_here...
	AngleUnit	e_unit = AU_radians;

	if ( !b_hyperbolic )
	{
		CLZNumber		scratch;
		CLZAnswerBase	probe( &x, &scratch );

		probe.Exact();
		e_unit = lzi_angle_unit_here();
	}

	// On the real axis, the real function as the host's call...
	if ( y.IsZero() )
	{
		lResult = sl_realTrig[ lFunction ]( &x, &p_target->m_real );

		if ( p_target->m_real.FlagIsOn( kfInfinity ) )
		{
			lzi_complex_flagged( p_target, l_base, kfInfinity, p_target->m_real.FlagIsOn( kfDiv0 ) );

			return lResult;
		}

		return lzi_first_failure( lResult, lzi_answer_as_function( &nought, &x, &p_target->m_imag, false, true, true ) );
	}

	// ...and on the imaginary axis in radians, one of the other family's...
	if ( x.IsZero() && (e_unit == AU_radians) )
	{
		const SLZTrigAxis&	axis = sl_trigOnImaginaryAxis[ lFunction ];
		CLZNumber*			p_part = (axis.m_lPart == 0) ? &p_target->m_real : &p_target->m_imag;
		CLZNumber*			p_other = (axis.m_lPart == 0) ? &p_target->m_imag : &p_target->m_real;

		{
			CLZRadiansHere	radians;

			lResult = axis.m_fn( &y, p_part );
		}

		if ( p_part->FlagIsOn( kfInfinity ) )
		{
			lzi_complex_flagged( p_target, l_base, kfInfinity, p_part->FlagIsOn( kfDiv0 ) );

			return lResult;
		}

		if ( axis.m_lSign < 0 )
			lz_negate( p_part );

		return lzi_first_failure( lResult, lzi_answer_as_function( &nought, &x, p_other, false, true, true ) );
	}

	// The circular argument and the hyperbolic one, its size in radians as
	//	a double -- an infinity past what one holds...
	CLZNumber*	p_alpha = b_hyperbolic ? &y : &x;
	CLZNumber*	p_beta = b_hyperbolic ? &x : &y;
	bool		b_inUnit = !b_hyperbolic && (e_unit != AU_radians);
	double		d_lnBase = std::log( static_cast<double>( l_base ) );
	double		d_beta = 0.0;
	CLZNumber	beta_magnitude;

	beta_magnitude.Assign( p_beta );
	beta_magnitude.m_lSign = kPositive;

	if ( !lzi_value_or_infinity( &beta_magnitude, &d_beta ) )
		d_beta = 0.0;

	if ( b_inUnit )
		d_beta *= 6.283185307179586 / static_cast<double>( lzi_turn_of( e_unit ) );

	// sin, cos, sinh and cosh past the answer limit refused, as cosh is --
	//	see lzi_answer_size -- and worked wider by its whole digits...
	long	l_magnitude = 0;

	if ( !b_quotient )
	{
		if ( lzi_answer_size( d_beta - std::log( 2.0 ), l_base ) == AS_too_large )
		{
			lzi_complex_flagged( p_target, l_base, kfInfinity, false );

			return kTooLarge;
		}

		l_magnitude = static_cast<long>( d_beta / d_lnBase ) + 2;
	}

	CLZNumber	parts[ 2 ];						// real, imaginary
	bool		b_nought[ 2 ] = { false, false };
	bool		b_refused = false;

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( &x, &scratch );

		work.Exact();

		long			l_guard = lzi_guard_places( kComplexGuardDigits, l_base );
		long			l_pole = b_quotient ? kTrigQuotientPlaces : 0;
		SLZExactPart	exact[ 2 ];				// sine, cosine

		// Both parts at lWidth places, cosh's digits and a small
		//	denominator's noughts on top...
		auto	fnParts = [ & ]( long lWidth, CLZNumber* pReal, CLZNumber* pImag, CLZNumber* pDenominator ) -> long
		{
			long	l_was = glWorkingPrecision;

			glWorkingPrecision = lWidth + l_magnitude + l_pole;

			// The hyperbolic argument in radians -- in a unit beta 2pi / turn,
			//	wide by beta's whole digits, or for a quotient with e^-|beta|
			//	past the width only near enough...
			bool		b_far = b_quotient && ((d_beta / d_lnBase) > static_cast<double>( glWorkingPrecision + 10 ));
			CLZNumber	beta;

			if ( b_inUnit )
			{
				long	l_lead = lzi_leading_exponent( p_beta );

				lzi_unit_to_radians( p_beta, e_unit, b_far ? 20 : (glWorkingPrecision + ((l_lead > 0) ? l_lead : 0) + 2), &beta );
			}
			else
				beta.Assign( p_beta );

			SLZTrigValues	values;
			long			l_result = lzi_complex_trig_values( p_alpha, b_inUnit ? e_unit : AU_radians, &beta, b_quotient,
															   b_far, &values );
			CLZNumber		denominator;

			exact[ 0 ] = values.m_sineExact;
			exact[ 1 ] = values.m_cosineExact;

			// Dc = c^2 u^2 + t^2, Ds = s^2 u^2 + t^2...
			if ( b_quotient )
			{
				ELZTrigFactor	e_circular = (function.m_lDenominator == kTrigDc) ? TF_c : TF_s;
				CLZNumber		uu;
				CLZNumber		tt;
				CLZNumber		square;
				CLZNumber		scaled;

				lz_multiply( &values.m_values[ TF_u ], &values.m_values[ TF_u ], &uu );
				lz_multiply( &values.m_values[ TF_t ], &values.m_values[ TF_t ], &tt );
				lz_multiply( &values.m_values[ e_circular ], &values.m_values[ e_circular ], &square );
				lz_multiply( &square, &uu, &scaled );
				lz_add( &scaled, &tt, &denominator );

				if ( pDenominator != nullptr )
					pDenominator->Assign( &denominator );
			}

			for ( long l_part = 0; l_part < 2; l_part++ )
			{
				CLZNumber*	p_out = (l_part == 0) ? pReal : pImag;

				if ( p_out == nullptr )
					continue;

				const SLZTrigPart&	part = function.m_parts[ l_part ];
				CLZNumber			product;

				b_nought[ l_part ] = lzi_trig_part_nought( part, values );

				product.SetBase( l_base, false );
				product.Assign( 1L );

				for ( ELZTrigFactor e_factor : part.m_eFactors )
				{
					if ( e_factor == TF_none )
						continue;

					CLZNumber	next;

					lz_multiply( &product, &values.m_values[ e_factor ], &next );
					product.Assign( &next );
				}

				if ( b_quotient )
				{
					CLZNumber	quotient;

					lz_divide( &product, &denominator, &quotient );
					product.Assign( &quotient );
				}

				if ( part.m_lSign < 0 )
					lz_negate( &product );

				p_out->Assign( &product );
			}

			glWorkingPrecision = l_was;

			return l_result;
		};

		/*
			**Near a pole the denominator is small**, and both parts worked at
			the width lose what the quotient magnifies: the work is taken
			again as wide as twice its noughts and four, and wider until it is
			seen at all. A denominator past twice the answer limit is a quotient past
			the limit -- |tan z| is about one over its square root -- and is
			refused, as the real function's is.
		*/
		CLZNumber	denominator;

		for ( ; ; )
		{
			lResult = fnParts( glWorkingPrecision + l_guard, &parts[ 0 ], &parts[ 1 ], b_quotient ? &denominator : nullptr );

			if ( !b_quotient )
				break;

			long	l_need = 0;

			denominator.CheckForZero();

			if ( denominator.IsZero() )
				l_need = (l_pole == 0) ? (glWorkingPrecision + l_guard) : (2 * l_pole);
			else
			{
				// ...its noughts after the point, which a lead of -1 has none
				//	of -- see kTrigQuotientPlaces...
				long	l_noughts = -lzi_leading_exponent( &denominator ) - 1;

				l_need = (l_noughts > 0) ? ((2 * l_noughts) + 4) : 0;
			}

			if ( l_need <= l_pole )
				break;

			if ( l_need > ((4 * gAnswerLimitPlaces) + 8) )
			{
				b_refused = true;
				break;
			}

			l_pole = l_need;
		}

		for ( long l_part = 0; (l_part < 2) && !b_refused; l_part++ )
		{
			if ( b_nought[ l_part ] )
				continue;

			// The real part of sin or cos in a unit is decided by the bound it
			//	lies over -- see lzi_complex_trig_decide -- and in radians a part of
			//	a small z by its polynomial -- see lzi_complex_trig_near_decide...
			const SLZTrigPart&			part = function.m_parts[ l_part ];
			bool						b_bounded = !b_hyperbolic && !b_quotient && lzi_trig_part_has( part, TF_ch );
			SLZExactPart				circular = lzi_trig_part_has( part, TF_s ) ? exact[ 0 ] : exact[ 1 ];
			const SLZTrigNearNought&	near = sl_trigNearNought[ lFunction ][ l_part ];

			auto	tieCall = [ & ]( CLZNumber* pResult, long lPlaces, long ) -> bool
			{
				if ( b_bounded && lzi_complex_trig_decide( pResult, circular, lPlaces ) )
					return true;

				// ...and in a unit with x a hair from an angle whose sine or
				//	cosine is a half, by the leading terms -- see
				//	lzi_complex_trig_half_decide...
				if ( b_bounded && b_inUnit &&
					 lzi_complex_trig_half_decide( pResult, lzi_trig_part_has( part, TF_s ) ? kAngleSine : kAngleCosine, &x, &y, lPlaces ) )
					return true;

				return !b_inUnit && lzi_complex_trig_near_decide( pResult, near, &x, &y, lPlaces );
			};

			lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth )
									  {
										  return fnParts( lWidth, (l_part == 0) ? pValue : nullptr,
														  (l_part == 1) ? pValue : nullptr, nullptr );
									  },
									  tieCall, l_base, l_guard, &parts[ l_part ] );
		}
	}

	if ( b_refused )
	{
		lzi_complex_flagged( p_target, l_base, kfInfinity, false );

		return kTooLarge;
	}

	for ( long l_part = 0; l_part < 2; l_part++ )
	{
		CLZNumber*	p_out = (l_part == 0) ? &p_target->m_real : &p_target->m_imag;

		lResult = lzi_first_failure( lResult, lzi_answer_as_function( b_nought[ l_part ] ? &nought : &parts[ l_part ], &x, p_out,
																	  false, b_nought[ l_part ], b_nought[ l_part ] ) );
	}

	return lResult;
}

long lzc_sine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kAngleSine, leftOperand, targetNumber );
}

long lzc_cosine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kAngleCosine, leftOperand, targetNumber );
}

long lzc_tangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kAngleTangent, leftOperand, targetNumber );
}

long lzc_cotangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kAngleCotangent, leftOperand, targetNumber );
}

long lzc_secant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kAngleSecant, leftOperand, targetNumber );
}

long lzc_cosecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kAngleCosecant, leftOperand, targetNumber );
}

long lzc_hypersine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kComplexHyperbolic + kAngleSine, leftOperand, targetNumber );
}

long lzc_hypercosine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kComplexHyperbolic + kAngleCosine, leftOperand, targetNumber );
}

long lzc_hypertangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kComplexHyperbolic + kAngleTangent, leftOperand, targetNumber );
}

long lzc_hypercotangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kComplexHyperbolic + kAngleCotangent, leftOperand, targetNumber );
}

long lzc_hypersecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kComplexHyperbolic + kAngleSecant, leftOperand, targetNumber );
}

long lzc_hypercosecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_trig( kComplexHyperbolic + kAngleCosecant, leftOperand, targetNumber );
}


/////////////////////////////////////////////////////////////////
//	Complex numbers, stage 5b -- the power, the roots and the other
//	bases, 2026-10-22

// A value as N / D, two whole numbers in lBase, exactly: N its digits with
//	its sign, D its own base to the power of its places, each converted as
//	a whole, which converts exactly...
static void lzi_as_whole_ratio( CLZNumber* pValue, long lBase, CLZNumber* pN, CLZNumber* pD )
{
	long	l_from = pValue->m_lBase;

	pN->Assign( pValue );
	pN->CheckForZero();
	pN->TrimTrailing();

	long	l_places = pN->m_lRadix;

	lzi_scale_by_base_power( pN, -l_places );
	pD->SetBase( l_from, false );
	pD->Assign( 1L );
	lzi_scale_by_base_power( pD, -l_places );

	if ( l_from != lBase )
	{
		pN->SetBase( lBase, true );
		pD->SetBase( lBase, true );
	}
}

/*
	**An exact quotient rounded at lPlaces, half-even by its magnitude's
	kept digit** -- pNumerator / pDenominator, both exact, in one base, the
	denominator not nought. Answers whether nothing was left over, the
	answer then exact and short. As lzi_ratio_at_places, for a quotient of
	two numbers rather than a number times a small fraction.
*/
static bool lzi_exact_quotient_at_places( CLZNumber* pNumerator, CLZNumber* pDenominator, long lPlaces, CLZNumber* pResult )
{
	long		l_base = pNumerator->m_lBase;
	CLZNumber	scaled;
	CLZNumber	divisor;

	scaled.Assign( pNumerator );
	scaled.CheckForZero();
	pResult->SetBase( l_base, false );

	if ( scaled.IsZero() )
	{
		pResult->Assign( 0L );

		return true;
	}

	bool	b_negative = (pNumerator->IsNegative() != pDenominator->IsNegative());

	scaled.m_lSign = kPositive;
	scaled.TrimTrailing();
	divisor.Assign( pDenominator );
	divisor.m_lSign = kPositive;
	divisor.TrimTrailing();

	// Both whole, and the base's power onto one side or the other...
	long	l_shift = divisor.m_lRadix - scaled.m_lRadix + lPlaces;

	lzi_scale_by_base_power( &scaled, -scaled.m_lRadix );
	lzi_scale_by_base_power( &divisor, -divisor.m_lRadix );

	if ( l_shift >= 0 )
		lzi_scale_by_base_power( &scaled, -l_shift );
	else
		lzi_scale_by_base_power( &divisor, l_shift );

	CLZNumber	whole;
	CLZNumber	product;
	CLZNumber	remainder;
	CLZNumber	twice;

	lz_divide_integer( &scaled, &divisor, &whole );
	lzi_complex_product( &whole, &divisor, &product );
	lz_subtract( &scaled, &product, &remainder );
	lzi_times_whole( &remainder, 2, &twice );

	long	l_order = lz_compare( &twice, &divisor );
	bool	b_up = (l_order == kGreaterThan);

	// ...a tie to the even kept digit...
	if ( l_order == kEqualTo )
	{
		CLZNumber	base;
		CLZNumber	digit;
		CLZNumber	two;
		CLZNumber	parity;

		base.SetBase( l_base, false );
		base.Assign( l_base );
		two.SetBase( l_base, false );
		two.Assign( 2L );
		lz_remainder( &whole, &base, &digit );
		lz_remainder( &digit, &two, &parity );
		parity.CheckForZero();
		b_up = !(parity.IsZero());
	}

	if ( b_up )
	{
		CLZNumber	one;

		one.SetBase( l_base, false );
		one.Assign( 1L );
		lz_add( &whole, &one );
	}

	lzi_scale_by_base_power( &whole, lPlaces );
	remainder.CheckForZero();

	bool	b_exact = remainder.IsZero();

	if ( b_exact )
		whole.TrimTrailing();

	whole.CheckForZero();

	if ( b_negative && !(whole.IsZero()) )
		lz_negate( &whole );

	pResult->Assign( &whole );

	return b_exact;
}

// (r + i m) times (c + i d), exactly, into r and m; either pair may be the
//	other...
static void lzi_complex_times( CLZNumber* pR, CLZNumber* pM, CLZNumber* pC, CLZNumber* pD )
{
	CLZNumber	rc;
	CLZNumber	md;
	CLZNumber	rd;
	CLZNumber	mc;

	lzi_complex_product( pR, pC, &rc );
	lzi_complex_product( pM, pD, &md );
	lzi_complex_product( pR, pD, &rd );
	lzi_complex_product( pM, pC, &mc );
	lz_subtract( &rc, &md, pR );
	lz_add( &rd, &mc, pM );
	pR->TrimTrailing();
	pM->TrimTrailing();
}

// (a + bi)^n exactly, n nought or more: binary powering, every product
//	exact...
static void lzi_complex_exact_power( CLZNumber* pA, CLZNumber* pB, long lPower, CLZNumber* pReal, CLZNumber* pImag )
{
	long		l_base = pA->m_lBase;
	CLZNumber	real;
	CLZNumber	imag;
	CLZNumber	square_real;
	CLZNumber	square_imag;

	real.SetBase( l_base, false );
	real.Assign( 1L );
	imag.SetBase( l_base, false );
	imag.Assign( 0L );
	square_real.Assign( pA );
	square_imag.Assign( pB );
	square_real.TrimTrailing();
	square_imag.TrimTrailing();

	for ( long n = lPower; n > 0; n >>= 1 )
	{
		if ( n & 1 )
			lzi_complex_times( &real, &imag, &square_real, &square_imag );

		if ( n > 1 )
			lzi_complex_times( &square_real, &square_imag, &square_real, &square_imag );
	}

	pReal->Assign( &real );
	pImag->Assign( &imag );
}

// |a + bi|^2, exactly...
static void lzi_complex_norm( CLZNumber* pA, CLZNumber* pB, CLZNumber* pNorm )
{
	CLZNumber	aa;
	CLZNumber	bb;

	lzi_complex_product( pA, pA, &aa );
	lzi_complex_product( pB, pB, &bb );
	lz_add( &aa, &bb, pNorm );
}

// An exact part written as a function's answer: short where it is exact,
//	rounded half-even at the width where not...
static long lzi_complex_exact_part( CLZNumber* pExact, CLZNumber* pOperand, CLZNumber* pTarget )
{
	CLZNumber	one;
	CLZNumber	rounded;

	one.SetBase( pExact->m_lBase, false );
	one.Assign( 1L );

	bool	b_exact = lzi_exact_quotient_at_places( pExact, &one, glWorkingPrecision, &rounded );

	return lzi_answer_as_function( &rounded, pOperand, pTarget, false, b_exact, b_exact );
}

/*
	**The parts of z^n for a whole n, exactly** -- worked by binary powering
	with every product exact and each part rounded once, so an answer that
	ends within the width is short and exact and one on a half goes to the
	even kept digit. A negative power is the conjugate over |z^n|^2, two
	exact quotients. Past kMostExactPowerPlaces of exact result the general
	road takes it -- see lzi_complex_power_general -- since an answer that
	long can be neither short nor exactly on a half.
*/
const long	kMostExactPowerPlaces = 20000;

static long lzi_complex_power_general( CLZNumber* pA, CLZNumber* pB, CLZComplex* pW, bool bReciprocal, CLZComplex* pTarget );

static long lzi_complex_whole_power( CLZNumber* pA, CLZNumber* pB, long lPower, CLZComplex* pTarget )
{
	long	l_base = pA->m_lBase;
	long	l_magnitude = (lPower < 0) ? -lPower : lPower;

	CLZNumber	a;
	CLZNumber	b;
	CLZNumber	norm;

	a.Assign( pA );
	b.Assign( pB );
	a.TrimTrailing();
	b.TrimTrailing();
	lzi_complex_norm( &a, &b, &norm );

	long	l_places = (a.m_lRadix > b.m_lRadix) ? a.m_lRadix : b.m_lRadix;
	bool	b_exactRoad = (static_cast<double>( l_places ) * static_cast<double>( l_magnitude )) <=
						  static_cast<double>( kMostExactPowerPlaces );

	// Past the exact road, the general one, which keeps a part that is
	//	exactly nought by its phase -- (1.0001i)^-563625's real part is, its
	//	imaginary part under the width...
	if ( !b_exactRoad )
	{
		CLZComplex	w;

		w.m_real.SetBase( l_base, false );
		w.m_real.Assign( lPower );
		w.m_imag.SetBase( l_base, false );
		w.m_imag.Assign( 0L );

		return lzi_complex_power_general( &a, &b, &w, false, pTarget );
	}

	// Refused past the limit, from |z|^n's size -- see lzi_answer_size --
	//	and worked exactly below the width, which keeps a part that is
	//	exactly nought, as (bi)^-4's imaginary part is...
	double	d_lnNorm = 0.0;

	if ( lzi_ln_estimate( &norm, &d_lnNorm ) &&
		 (lzi_answer_size( 0.5 * d_lnNorm * static_cast<double>( lPower ), l_base ) == AS_too_large) )
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, false );

		return kTooLarge;
	}

	CLZNumber	real;
	CLZNumber	imag;

	lzi_complex_exact_power( &a, &b, l_magnitude, &real, &imag );

	if ( lPower >= 0 )
	{
		long	lResult = lzi_complex_exact_part( &real, &a, &pTarget->m_real );

		return lzi_first_failure( lResult, lzi_complex_exact_part( &imag, &a, &pTarget->m_imag ) );
	}

	// One over z^n: its conjugate over |z^n|^2, two exact quotients...
	CLZNumber	divisor;
	CLZNumber	parts[ 2 ];
	bool		b_exact[ 2 ];

	lzi_complex_norm( &real, &imag, &divisor );
	lz_negate( &imag );
	b_exact[ 0 ] = lzi_exact_quotient_at_places( &real, &divisor, glWorkingPrecision, &parts[ 0 ] );
	b_exact[ 1 ] = lzi_exact_quotient_at_places( &imag, &divisor, glWorkingPrecision, &parts[ 1 ] );

	long	lResult = lzi_answer_as_function( &parts[ 0 ], &a, &pTarget->m_real, false, b_exact[ 0 ], b_exact[ 0 ] );

	return lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 1 ], &a, &pTarget->m_imag, false, b_exact[ 1 ], b_exact[ 1 ] ) );
}

/*
	**A part of a rational power in doubt, decided by a candidate root and
	its exact residue** -- for z^(p/q), the root of order q of Z = z^|p|,
	of order -q where p is negative. The part in doubt is pinned to the
	half h it lies near, the other part taken from the work, and with that
	candidate r every term is exact: for a positive order the root is
	r (1 + t)^(1/q), t = (Z - r^q) / r^q, and for a negative one r (1 +
	t)^(-1/q), t = Z r^q - 1. The first-order term is exact -- (Z - r^q) /
	(q r^(q-1)) and -r t / q -- and the rest is under 2 |r| |t|^2 wherever
	|t| is a quarter or less. So the part's side of h is the first-order
	term's sign wherever that term outweighs the bound, all compared in
	squares of exact numbers; and the other part is refined by the same
	term, Newton's step, until it does or the places run out.

	**Worked in 2r**, which ends where r does not: a half of the places
	never ends in an odd base, and twice it does. The positive order's
	tests keep their form there with Z taken as 2^q Z; the negative's
	carry the powers of two.

	**A residue of nought is the root itself**: the part is exactly the
	half, and goes to the even kept digit. That is the tie a rational root
	makes -- (0.5 + 1.5i)^3 rooted at no places, or a half in base three --
	and the near tie beside it, a residue as small as one likes, is the
	one no retake reaches. The real roots decide theirs by comparing whole
	numbers -- see lzi_whole_order_root_decide -- which a complex root has
	no order for.
*/
static bool lzi_complex_root_decide( CLZNumber* pResult, long lPart, CLZNumber* pOther, CLZNumber* pZr, CLZNumber* pZi,
									 long lOrder, long lPlaces, long lMostPlaces )
{
	if ( (lPlaces < 0) || (lOrder == 0) )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_q = (lOrder < 0) ? -lOrder : lOrder;
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pResult->m_lRadix > lPlaces) ? pResult->m_lRadix : lPlaces) + 4;

	CLZNumber	y;
	CLZNumber	unit;
	CLZNumber	two;
	CLZNumber	one;
	CLZNumber	q;
	CLZNumber	powerOfTwo;				// 2^q

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );
	two.SetBase( l_base, false );
	two.Assign( 2L );
	one.SetBase( l_base, false );
	one.Assign( 1L );
	q.SetBase( l_base, false );
	q.Assign( l_q );
	lzi_exact_power( &two, l_q, &powerOfTwo );

	bool	b_negative = pResult->IsNegative();

	// 2r: the part twice the half, 2y + u, and the other part doubled...
	CLZNumber	r[ 2 ];

	lz_add( &y, &y, &r[ lPart ] );
	lz_add( &r[ lPart ], &unit );

	if ( b_negative )
		lz_negate( &r[ lPart ] );

	r[ 1 - lPart ].Assign( pOther );
	r[ 1 - lPart ].TrimTrailing();
	lz_add( &r[ 1 - lPart ], &r[ 1 - lPart ] );

	// ...and Z, times 2^q for a positive order...
	CLZNumber	z_real;
	CLZNumber	z_imag;

	z_real.Assign( pZr );
	z_imag.Assign( pZi );

	if ( lOrder > 0 )
	{
		CLZNumber	scaled;

		lzi_complex_product( &z_real, &powerOfTwo, &scaled );
		z_real.Assign( &scaled );
		lzi_complex_product( &z_imag, &powerOfTwo, &scaled );
		z_imag.Assign( &scaled );
	}

	long	l_most = lMostPlaces / (4 * l_q);
	long	l_places = r[ 1 - lPart ].m_lRadix + 2;
	long	l_side = 0;
	bool	b_tie = false;

	for ( ; ; )
	{
		CLZNumber	rq_real;
		CLZNumber	rq_imag;
		CLZNumber	t_real;
		CLZNumber	t_imag;
		CLZNumber	step_real;
		CLZNumber	step_imag;
		CLZNumber	norm_t;
		CLZNumber	norm_r;
		CLZNumber	left;
		CLZNumber	right;
		CLZNumber	sixteen;
		CLZNumber	xx;
		CLZNumber	tt;
		CLZNumber	factor;

		sixteen.SetBase( l_base, false );
		sixteen.Assign( 16L );
		lzi_complex_norm( &r[ 0 ], &r[ 1 ], &norm_r );

		if ( lOrder > 0 )
		{
			// delta = 2^q Z - (2r)^q, T = delta conj((2r)^(q-1)), K = q |(2r)^(q-1)|^2...
			CLZNumber	rq1_real;
			CLZNumber	rq1_imag;
			CLZNumber	norm_rq1;
			CLZNumber	norm_rq;
			CLZNumber	conj_imag;
			CLZNumber	k;

			lzi_complex_exact_power( &r[ 0 ], &r[ 1 ], l_q - 1, &rq1_real, &rq1_imag );
			rq_real.Assign( &rq1_real );
			rq_imag.Assign( &rq1_imag );
			lzi_complex_times( &rq_real, &rq_imag, &r[ 0 ], &r[ 1 ] );
			lz_subtract( &z_real, &rq_real, &t_real );
			lz_subtract( &z_imag, &rq_imag, &t_imag );
			t_real.CheckForZero();
			t_imag.CheckForZero();

			if ( t_real.IsZero() && t_imag.IsZero() )
			{
				b_tie = true;
				break;
			}

			// |t| a quarter or less: 16 |delta|^2 <= |(2r)^q|^2...
			lzi_complex_norm( &t_real, &t_imag, &norm_t );
			lzi_complex_norm( &rq_real, &rq_imag, &norm_rq );
			lzi_complex_product( &norm_t, &sixteen, &left );

			if ( lz_compare( &left, &norm_rq ) == kGreaterThan )
				break;

			conj_imag.Assign( &rq1_imag );
			lz_negate( &conj_imag );
			step_real.Assign( &t_real );
			step_imag.Assign( &t_imag );
			lzi_complex_times( &step_real, &step_imag, &rq1_real, &conj_imag );

			// Decided where X^2 |2r|^2 > 4 q^2 |delta|^4...
			CLZNumber*	p_x = (lPart == 0) ? &step_real : &step_imag;

			lzi_complex_product( p_x, p_x, &xx );
			lzi_complex_product( &xx, &norm_r, &left );
			lzi_complex_product( &norm_t, &norm_t, &tt );
			factor.SetBase( l_base, false );
			factor.Assign( 4 * l_q * l_q );
			lzi_complex_product( &tt, &factor, &right );
			p_x->CheckForZero();

			if ( !(p_x->IsZero()) && (lz_compare( &left, &right ) == kGreaterThan) )
			{
				l_side = p_x->IsNegative() ? -1 : 1;
				break;
			}

			// ...else Newton's step on the other part, its component of T over K...
			lzi_complex_norm( &rq1_real, &rq1_imag, &norm_rq1 );
			lzi_complex_product( &norm_rq1, &q, &k );
			glWorkingPrecision = (2 * l_places) + 4;

			if ( glWorkingPrecision > l_most )
				break;

			CLZNumber	increment;

			lz_divide( (lPart == 0) ? &step_imag : &step_real, &k, &increment );
			lz_add( &r[ 1 - lPart ], &increment );
			r[ 1 - lPart ].TrimTrailing();
			l_places = glWorkingPrecision;
		}
		else
		{
			// t' = Z (2r)^q - 2^q, and T = 2r t'...
			lzi_complex_exact_power( &r[ 0 ], &r[ 1 ], l_q, &rq_real, &rq_imag );
			t_real.Assign( &z_real );
			t_imag.Assign( &z_imag );
			lzi_complex_times( &t_real, &t_imag, &rq_real, &rq_imag );
			lz_subtract( &t_real, &powerOfTwo );
			t_real.CheckForZero();
			t_imag.CheckForZero();

			if ( t_real.IsZero() && t_imag.IsZero() )
			{
				b_tie = true;
				break;
			}

			// |t| a quarter or less: 16 |t'|^2 <= 2^2q...
			CLZNumber	four_q;

			lzi_complex_norm( &t_real, &t_imag, &norm_t );
			lzi_complex_product( &norm_t, &sixteen, &left );
			lzi_complex_product( &powerOfTwo, &powerOfTwo, &four_q );

			if ( lz_compare( &left, &four_q ) == kGreaterThan )
				break;

			step_real.Assign( &r[ 0 ] );
			step_imag.Assign( &r[ 1 ] );
			lzi_complex_times( &step_real, &step_imag, &t_real, &t_imag );

			// Decided where X^2 2^(2q-2) > q^2 |2r|^2 |t'|^4...
			CLZNumber*	p_x = (lPart == 0) ? &step_real : &step_imag;
			CLZNumber	scaled;
			CLZNumber	quarter;

			lzi_complex_product( p_x, p_x, &xx );
			lzi_complex_product( &xx, &four_q, &left );
			lzi_complex_product( &norm_t, &norm_t, &tt );
			lzi_complex_product( &tt, &norm_r, &scaled );
			factor.SetBase( l_base, false );
			factor.Assign( 4 * l_q * l_q );
			lzi_complex_product( &scaled, &factor, &right );
			p_x->CheckForZero();

			if ( !(p_x->IsZero()) && (lz_compare( &left, &right ) == kGreaterThan) )
			{
				l_side = p_x->IsNegative() ? 1 : -1;			// the step is -r t / q
				break;
			}

			// ...else Newton's step: the other part less its component of T
			//	over q 2^q...
			CLZNumber	k;

			lzi_complex_product( &q, &powerOfTwo, &k );
			glWorkingPrecision = (2 * l_places) + 4;

			if ( glWorkingPrecision > l_most )
				break;

			CLZNumber	increment;

			lz_divide( (lPart == 0) ? &step_imag : &step_real, &k, &increment );
			lz_subtract( &r[ 1 - lPart ], &increment );
			r[ 1 - lPart ].TrimTrailing();
			l_places = glWorkingPrecision;
		}
	}

	if ( !b_tie && (l_side == 0) )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// The part above the half or below it, signed; on it, the even kept
	//	digit of the magnitude...
	bool	b_up = false;

	if ( b_tie )
	{
		CLZNumber	whole;
		CLZNumber	base;
		CLZNumber	digit;
		CLZNumber	parity;

		whole.Assign( &y );
		lzi_scale_by_base_power( &whole, -lPlaces );
		base.SetBase( l_base, false );
		base.Assign( l_base );
		lz_remainder( &whole, &base, &digit );
		lz_remainder( &digit, &two, &parity );
		parity.CheckForZero();
		b_up = !(parity.IsZero());
	}
	else
		b_up = ((l_side > 0) != b_negative);

	glWorkingPrecision = lPlaces + 2;

	if ( b_up )
		lz_add( &y, &unit );

	if ( b_negative )
		lz_negate( &y );

	y.CheckForZero();
	pResult->Assign( &y );
	glWorkingPrecision = l_saved;

	return true;
}

// cos^2 of a whole number of degrees, as n/d, for a multiple of 30 or 45
//	-- the only whole degrees whose squared cosine is rational -- and false
//	for any other...
static bool lzi_cosine_squared( long lDegrees, long* plN, long* plD )
{
	static const struct { long m_lDegrees; long m_lN; long m_lD; } kSquares[] =
	{
		{ 0, 1, 1 }, { 30, 3, 4 }, { 45, 1, 2 }, { 60, 1, 4 }, { 90, 0, 1 }, { 120, 1, 4 }, { 135, 1, 2 }, { 150, 3, 4 },
	};

	long	l_mod = ((lDegrees % 180) + 180) % 180;

	for ( auto& square : kSquares )
	{
		if ( square.m_lDegrees == l_mod )
		{
			*plN = square.m_lN;
			*plD = square.m_lD;

			return true;
		}
	}

	return false;
}

/*
	**z^E for E = w or E = 1/w, e^(E ln z)** -- the general road of the
	power, the roots and the base's exponential. Each part correctly
	rounded, half-even, in places, retaken in doubt.

	- **ln z is ½ ln(a^2 + b^2) + i atan2(b, a)**, worked as wide again as
	  the answer's whole digits and E's, so E ln z keeps the width the
	  exponential needs; then e^Re times the cosine and sine of Im, the
	  pair wider by e^Re's digits, as lzc_exp_e works.
	- **On an axis with a real exponent, or at a unit on an axis, the
	  phase is exact in degrees** where it is whole -- E times 0, 180 or
	  plus or minus 90 -- and its sine and cosine come back marked exact
	  where they are, as a unit's do: a sine or cosine of nought is an exact
	  nought, and a part of rational square is proved rational -- so
	  (-8)^(1/3) is 2 (1/2 + i sqrt(3)/2), its real part exactly 1, and
	  (-4)^(1/2) is exactly 2i.
	- **A rational exponent's rational answer is proved** by raising it
	  back -- r^q against z^p, as the real power proves one, see
	  lzi_power_is_exactly -- and then rounded exactly, ties included.
	- **Ties no retake reaches** are decided by a candidate root's exact
	  residue -- see lzi_complex_root_decide -- and, on an axis, by e M's
	  side of e, M being |z|^E and its side of one known -- see
	  lzi_complex_trig_decide.

	z is not nought and neither is E; the edges are the callers'.
*/

static long lzi_complex_power_general( CLZNumber* pA, CLZNumber* pB, CLZComplex* pW, bool bReciprocal, CLZComplex* pTarget )
{
	long	l_base = pA->m_lBase;
	long	lResult = kNoError;

	CLZNumber	a;
	CLZNumber	b;
	CLZNumber	norm;

	a.Assign( pA );
	b.Assign( pB );
	a.CheckForZero();
	b.CheckForZero();
	a.TrimTrailing();
	b.TrimTrailing();
	lzi_complex_norm( &a, &b, &norm );

	// The exponent's parts as whole ratios in the value's base, exactly...
	CLZNumber	c_n;
	CLZNumber	c_d;
	CLZNumber	d_n;
	CLZNumber	d_d;

	lzi_as_whole_ratio( &pW->m_real, l_base, &c_n, &c_d );
	lzi_as_whole_ratio( &pW->m_imag, l_base, &d_n, &d_d );
	c_n.CheckForZero();
	d_n.CheckForZero();

	bool	b_realExponent = d_n.IsZero();

	// ...and a real one as p/q, as the real power reads one: q to 64, p to
	//	4096 either way...
	long	l_p = 0;
	long	l_q = 0;

	if ( b_realExponent )
	{
		long	l_cn = 0;
		long	l_cd = 0;

		if ( lzi_whole_as_long( &c_n, &l_cn ) && lzi_whole_as_long( &c_d, &l_cd ) && (l_cn != 0) )
		{
			long	l_gcd = std::gcd( l_cn, l_cd );
			long	l_num = bReciprocal ? (l_cd / l_gcd) : (l_cn / l_gcd);
			long	l_den = bReciprocal ? (l_cn / l_gcd) : (l_cd / l_gcd);

			if ( c_n.IsNegative() )
				l_num = -l_num;

			if ( (l_den <= 64) && (l_num <= 4096) && (l_num >= -4096) )
			{
				l_p = l_num;
				l_q = l_den;
			}
		}
	}

	// On an axis, the angle in degrees...
	bool	b_onReal = b.IsZero();
	bool	b_onImag = a.IsZero();
	long	l_theta = b_onReal ? (a.IsNegative() ? 180 : 0) : (b.IsNegative() ? -90 : 90);
	bool	b_unitModulus = lzi_exactly( &norm, 1 );
	bool	b_degreePhase = false;
	CLZNumber	phase;					// E's phase in degrees, where whole

	if ( (b_onReal || b_onImag) && (b_realExponent || b_unitModulus) )
	{
		// theta c for a power, theta c / (c^2 + d^2) for a root, as whole
		//	numbers over whole numbers...
		CLZNumber	theta;
		CLZNumber	numerator;
		CLZNumber	denominator;

		theta.SetBase( l_base, false );
		theta.Assign( l_theta );

		if ( !bReciprocal )
		{
			lzi_complex_product( &theta, &c_n, &numerator );
			denominator.Assign( &c_d );
		}
		else
		{
			CLZNumber	dd2;
			CLZNumber	step;
			CLZNumber	cc;
			CLZNumber	dn2;
			CLZNumber	cd2;
			CLZNumber	first;
			CLZNumber	second;

			lzi_complex_product( &d_d, &d_d, &dd2 );
			lzi_complex_product( &theta, &c_n, &step );
			lzi_complex_product( &step, &c_d, &cc );
			lzi_complex_product( &cc, &dd2, &numerator );
			lzi_complex_product( &c_n, &c_n, &step );
			lzi_complex_product( &step, &dd2, &first );
			lzi_complex_product( &d_n, &d_n, &dn2 );
			lzi_complex_product( &c_d, &c_d, &cd2 );
			lzi_complex_product( &dn2, &cd2, &second );
			lz_add( &first, &second, &denominator );
		}

		denominator.CheckForZero();

		if ( !(denominator.IsZero()) )
		{
			CLZNumber	quotient;
			CLZNumber	product;
			CLZNumber	rest;

			lz_divide_integer( &numerator, &denominator, &quotient );
			lzi_complex_product( &quotient, &denominator, &product );
			lz_subtract( &numerator, &product, &rest );
			rest.CheckForZero();

			if ( rest.IsZero() )
			{
				b_degreePhase = true;
				phase.Assign( &quotient );
			}
		}
	}

	/*
		**A root of a unit, of complex order c + id with d small** --
		(-1)^(1/(3 - i 3^-34)). Its phase theta c / (c^2 + d^2) is theta /
		c moved by d^2, while its modulus e^u, u = theta d / (c^2 + d^2),
		moves by d itself: so where theta / c is whole degrees with an
		exact sine or cosine e, the part is e^u e plus a change under
		|theta d^2 / c|, and wherever |d| < |e c| / 4 and 7 |d| <= c^2 --
		e^u then under 1.65 -- that is outweighed, and the part lies on u's
		side of e: over |e| for u past nought, under it below. A tie a hair from a half by construction, past any
		retake; decided as the unit's half is -- see lzi_complex_trig_decide.
	*/
	SLZExactPart	unit_root[ 2 ];				// cosine, sine
	long			l_unitRootSide = 0;

	unit_root[ 0 ] = { false, 0, 1 };
	unit_root[ 1 ] = { false, 0, 1 };

	if ( b_unitModulus && (b_onReal || b_onImag) && bReciprocal && !b_realExponent && !b_degreePhase && (l_theta != 0) )
	{
		CLZNumber	theta;
		CLZNumber	numerator;
		CLZNumber	quotient;
		CLZNumber	product;
		CLZNumber	rest;

		theta.SetBase( l_base, false );
		theta.Assign( l_theta );
		lzi_complex_product( &theta, &c_d, &numerator );
		c_n.CheckForZero();

		if ( !(c_n.IsZero()) )
		{
			lz_divide_integer( &numerator, &c_n, &quotient );
			lzi_complex_product( &quotient, &c_n, &product );
			lz_subtract( &numerator, &product, &rest );
			rest.CheckForZero();

			if ( rest.IsZero() )
			{
				CLZNumber	sine;
				CLZNumber	cosine;
				long		l_was = glWorkingPrecision;

				glWorkingPrecision = 4;
				lzi_sine_cosine_work( &quotient, AU_degrees, true, &sine, &cosine, &unit_root[ 1 ], &unit_root[ 0 ] );
				glWorkingPrecision = l_was;

				// |d| < |e c| / 4, as whole numbers: 4 |d_n| den c_d < |n| |c_n| d_d...
				for ( long l_part = 0; l_part < 2; l_part++ )
				{
					SLZExactPart&	e = unit_root[ l_part ];

					if ( !e.m_bExact || (e.m_lNumerator == 0) )
					{
						e.m_bExact = false;
						continue;
					}

					CLZNumber	left;
					CLZNumber	right;
					CLZNumber	step;
					CLZNumber	factor;
					CLZNumber	dn;
					CLZNumber	cn;

					dn.Assign( &d_n );
					dn.m_lSign = kPositive;
					cn.Assign( &c_n );
					cn.m_lSign = kPositive;
					factor.SetBase( l_base, false );
					factor.Assign( 4 * e.m_lDenominator );
					lzi_complex_product( &dn, &factor, &step );
					lzi_complex_product( &step, &c_d, &left );
					factor.Assign( (e.m_lNumerator < 0) ? -e.m_lNumerator : e.m_lNumerator );
					lzi_complex_product( &cn, &factor, &step );
					lzi_complex_product( &step, &d_d, &right );

					if ( lz_compare( &left, &right ) != kLessThan )
						e.m_bExact = false;
				}

				// ...|u| at most a half, so e^u is under 1.65: 7 |d| <= c^2, 2 pi
				//	being under seven...
				{
					CLZNumber	dn;
					CLZNumber	seven;
					CLZNumber	step;
					CLZNumber	left;
					CLZNumber	cc;
					CLZNumber	right;
					CLZNumber	cd2;

					dn.Assign( &d_n );
					dn.m_lSign = kPositive;
					seven.SetBase( l_base, false );
					seven.Assign( 7L );
					lzi_complex_product( &dn, &seven, &step );
					lzi_complex_product( &c_d, &c_d, &cd2 );
					lzi_complex_product( &step, &cd2, &left );
					lzi_complex_product( &c_n, &c_n, &cc );
					lzi_complex_product( &cc, &d_d, &right );

					if ( lz_compare( &left, &right ) == kGreaterThan )
					{
						unit_root[ 0 ].m_bExact = false;
						unit_root[ 1 ].m_bExact = false;
					}
				}

				// ...and u's side: the sign of theta d...
				l_unitRootSide = ((d_n.IsNegative() == (l_theta < 0)) ? 1 : -1);
			}
		}
	}

	// Sizes, as doubles: E ln z = u + iv, the answer about e^u...
	double	d_lnNorm = 0.0;
	double	d_cValue = 0.0;
	double	d_dValue = 0.0;
	double	d_a = 0.0;
	double	d_b = 0.0;

	lzi_ln_estimate( &norm, &d_lnNorm );
	lzi_value_or_infinity( &pW->m_real, &d_cValue );
	lzi_value_or_infinity( &pW->m_imag, &d_dValue );
	lzi_value_or_infinity( &a, &d_a );
	lzi_value_or_infinity( &b, &d_b );

	double	d_lnr = 0.5 * d_lnNorm;
	double	d_theta = std::atan2( d_b, d_a );
	double	d_cd = (d_cValue * d_cValue) + (d_dValue * d_dValue);
	double	d_u = bReciprocal ? (((d_cValue * d_lnr) + (d_dValue * d_theta)) / d_cd) : ((d_cValue * d_lnr) - (d_dValue * d_theta));
	double	d_size = bReciprocal ? (1.0 / std::sqrt( d_cd )) : std::sqrt( d_cd );

	ELZAnswerSize	e_size = lzi_answer_size( d_u, l_base );

	if ( e_size == AS_too_large )
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, false );

		return kTooLarge;
	}

	// ...nought below the width, a part whose sine or cosine is exactly
	//	nought an exact one: the phase whole degrees, a cosine of nought at
	//	90 more than a multiple of 180 and a sine at a multiple...
	if ( e_size == AS_nought )
	{
		CLZNumber	nought;
		bool		b_exactNought[ 2 ] = { false, false };

		nought.SetBase( l_base, false );
		nought.Assign( 0L );

		if ( b_degreePhase )
		{
			CLZNumber	turn;
			CLZNumber	rest;

			turn.SetBase( l_base, false );
			turn.Assign( 180L );
			lz_remainder( &phase, &turn, &rest );
			rest.CheckForZero();
			b_exactNought[ 1 ] = rest.IsZero();

			if ( !b_exactNought[ 1 ] )
			{
				CLZNumber	ninety;

				ninety.SetBase( l_base, false );
				ninety.Assign( 90L );

				if ( rest.IsNegative() )
					lz_add( &rest, &turn );

				b_exactNought[ 0 ] = (lz_compare( &rest, &ninety ) == kEqualTo);
			}
		}

		lResult = lzi_answer_as_function( &nought, &a, &pTarget->m_real, false, b_exactNought[ 0 ], b_exactNought[ 0 ] );

		return lzi_first_failure( lResult, lzi_answer_as_function( &nought, &a, &pTarget->m_imag, false,
																	b_exactNought[ 1 ], b_exactNought[ 1 ] ) );
	}

	long	l_magnitude = lzi_places_for( d_u, l_base );
	long	l_exponent = lzi_places_for( std::log( (d_size > 1.0) ? d_size : 1.0 ) + std::log( std::fabs( d_lnr ) + std::fabs( d_theta ) + 1.0 ),
										 l_base );

	if ( l_magnitude < 0 )
		l_magnitude = 0;

	if ( l_exponent < 0 )
		l_exponent = 0;

	// The modulus's side of one, for the axis's decider: |z| against one
	//	and E's sign, or at a unit the sign of u...
	long	l_side = 0;

	if ( b_degreePhase )
	{
		// u is -theta d for a power at a unit, theta d / (c^2 + d^2) for a
		//	root; elsewhere ln |z| times c's sign, |z| against one exactly...
		if ( b_unitModulus )
		{
			bool	b_up = ((d_n.IsNegative() != (l_theta < 0)) == bReciprocal);

			l_side = b_realExponent ? 0 : (b_up ? -1 : 1);
		}
		else
		{
			CLZNumber	one;

			one.SetBase( l_base, false );
			one.Assign( 1L );

			bool	b_large = (lz_compare( &norm, &one ) == kGreaterThan);

			l_side = ((b_large == !(c_n.IsNegative())) ? 1 : -1);
		}
	}

	CLZNumber		parts[ 2 ];
	bool			b_exact[ 2 ] = { false, false };
	bool			b_nought[ 2 ] = { false, false };

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( &a, &scratch );

		work.Exact();

		long			l_guard = lzi_guard_places( kComplexGuardDigits, l_base );
		SLZExactPart	exact[ 2 ];				// cosine, sine

		// Both parts at lWidth places...
		auto	fnParts = [ & ]( long lWidth, CLZNumber* pReal, CLZNumber* pImag ) -> long
		{
			long	l_was = glWorkingPrecision;
			long	l_pairWidth = lWidth + l_magnitude;
			long	l_logWidth = l_pairWidth + l_exponent + 4;
			long	l_result = kNoError;

			glWorkingPrecision = l_logWidth;

			CLZNumber	c;
			CLZNumber	d;
			CLZNumber	u;
			CLZNumber	v;

			lz_divide( &c_n, &c_d, &c );
			lz_divide( &d_n, &d_d, &d );

			if ( b_degreePhase )
			{
				// u: |z|^E's logarithm, or at a unit -theta d, theta d / (c^2 + d^2)...
				{
					if ( b_unitModulus )
					{
						CLZNumber	theta;
						CLZNumber	radians;

						theta.SetBase( l_base, false );
						theta.Assign( l_theta );
						lzi_unit_to_radians( &theta, AU_degrees, l_logWidth, &radians );

						if ( bReciprocal )
						{
							CLZNumber	cc;
							CLZNumber	dd;
							CLZNumber	sum;
							CLZNumber	product;

							lzi_complex_product( &c, &c, &cc );
							lzi_complex_product( &d, &d, &dd );
							lz_add( &cc, &dd, &sum );
							lz_multiply( &radians, &d, &product );
							lz_divide( &product, &sum, &u );
						}
						else
						{
							lz_multiply( &radians, &d, &u );
							lz_negate( &u );
						}
					}
					else
					{
						CLZNumber	lnr;
						CLZNumber	magnitude;

						magnitude.Assign( b_onReal ? &a : &b );
						magnitude.m_lSign = kPositive;
						l_result = lzi_first_failure( l_result, lz_log_e( &magnitude, &lnr ) );

						if ( bReciprocal )
							lz_divide( &lnr, &c, &u );
						else
							lz_multiply( &lnr, &c, &u );
					}
				}
			}
			else
			{
				CLZNumber	lnr;
				CLZNumber	theta;
				CLZNumber	two;

				two.SetBase( l_base, false );
				two.Assign( 2L );
				l_result = lzi_first_failure( l_result, lz_log_e( &norm, &lnr ) );
				lz_divide( &lnr, &two );

				{
					CLZRadiansHere	radians;

					l_result = lzi_first_failure( l_result, lz_arctangent2( &b, &a, &theta ) );
				}

				CLZNumber	cl;
				CLZNumber	dt;
				CLZNumber	ct;
				CLZNumber	dl;

				lz_multiply( &c, &lnr, &cl );
				lz_multiply( &d, &theta, &dt );
				lz_multiply( &c, &theta, &ct );
				lz_multiply( &d, &lnr, &dl );

				if ( bReciprocal )
				{
					CLZNumber	cc;
					CLZNumber	dd;
					CLZNumber	sum;
					CLZNumber	top;

					lzi_complex_product( &c, &c, &cc );
					lzi_complex_product( &d, &d, &dd );
					lz_add( &cc, &dd, &sum );
					lz_add( &cl, &dt, &top );
					lz_divide( &top, &sum, &u );
					lz_subtract( &ct, &dl, &top );
					lz_divide( &top, &sum, &v );
				}
				else
				{
					lz_subtract( &cl, &dt, &u );
					lz_add( &dl, &ct, &v );
				}
			}

			// The modulus, e^u...
			CLZNumber	m;

			glWorkingPrecision = lWidth + 2;
			l_result = lzi_first_failure( l_result, lz_exp_e( &u, &m ) );

			// ...and the pair, in degrees where the phase is exact...
			glWorkingPrecision = l_pairWidth;

			CLZNumber	sine;
			CLZNumber	cosine;

			if ( b_degreePhase )
				l_result = lzi_first_failure( l_result, lzi_sine_cosine_work( &phase, AU_degrees, true, &sine, &cosine,
																			  &exact[ 1 ], &exact[ 0 ] ) );
			else
				l_result = lzi_first_failure( l_result, lzi_sine_cosine_work( &v, AU_radians, false, &sine, &cosine,
																			  &exact[ 1 ], &exact[ 0 ] ) );

			for ( long l_part = 0; l_part < 2; l_part++ )
			{
				CLZNumber*		p_out = (l_part == 0) ? pReal : pImag;
				CLZNumber*		p_factor = (l_part == 0) ? &cosine : &sine;
				SLZExactPart&	e = exact[ l_part ];

				if ( p_out == nullptr )
					continue;

				if ( e.m_bExact )
				{
					CLZNumber	numerator;
					CLZNumber	denominator;

					numerator.SetBase( l_base, false );
					numerator.Assign( e.m_lNumerator );
					denominator.SetBase( l_base, false );
					denominator.Assign( e.m_lDenominator );
					lz_divide( &numerator, &denominator, p_factor );
				}

				CLZNumber	product;

				b_nought[ l_part ] = e.m_bExact && (e.m_lNumerator == 0);
				lz_multiply( &m, p_factor, &product );
				p_out->Assign( &product );
			}

			glWorkingPrecision = l_was;

			return l_result;
		};

		lResult = fnParts( glWorkingPrecision + l_guard, &parts[ 0 ], &parts[ 1 ] );

		/*
			**A rational exponent's rational answer, proved** -- in 2r,
			whose parts end where a half of them never does: a root of
			(2/9)i to -1/2 is 3/2 - 3/2 i, a half at every width in base
			three, and twice it is 10 - 10i. The work doubled and rounded to
			the places 2r could have -- |p| times z's over q and two, the
			reciprocal of |z|^2|p| besides for p below nought, or all the
			work has where that cannot be told -- and raised back: (2r)^q
			against 2^q z^p, or (2r)^q z^|p| against 2^q. Where it holds each
			part is exact, half of a part of 2r, rounded so: short where it
			ends, and on a half to the even kept digit. Asked before the
			retake, so an exact tie never reaches it.
		*/
		bool	b_proved = false;

		if ( l_q >= 1 )
		{
			long	l_zPlaces = (a.m_lRadix > b.m_lRadix) ? a.m_lRadix : b.m_lRadix;
			long	l_powerPlaces = l_zPlaces * ((l_p < 0) ? -l_p : l_p);

			if ( l_p < 0 )
			{
				long	l_reciprocal = lzi_places_of_power( &norm, -l_p, true );

				l_powerPlaces = (l_reciprocal >= 0) ? (l_powerPlaces + l_reciprocal) : -1;
			}

			long	l_rootPlaces = (l_powerPlaces >= 0) ? ((l_powerPlaces / l_q) + 2) : (glWorkingPrecision + l_guard - 4);

			if ( l_rootPlaces <= ((4 * glWorkingPrecision) + 400) )
			{
				CLZNumber	wide[ 2 ];

				if ( l_rootPlaces + 2 > glWorkingPrecision + l_guard )
					fnParts( l_rootPlaces + l_guard, &wide[ 0 ], &wide[ 1 ] );
				else
				{
					wide[ 0 ].Assign( &parts[ 0 ] );
					wide[ 1 ].Assign( &parts[ 1 ] );
				}

				long	l_was = glWorkingPrecision;

				glWorkingPrecision = l_rootPlaces + 4;

				for ( long l_part = 0; l_part < 2; l_part++ )
				{
					CLZNumber	twice;

					lz_add( &wide[ l_part ], &wide[ l_part ], &twice );
					twice.SetPrecision( l_rootPlaces );
					twice.Round();
					twice.TrimTrailing();
					wide[ l_part ].Assign( &twice );
				}

				glWorkingPrecision = l_was;

				CLZNumber	r_real;
				CLZNumber	r_imag;
				CLZNumber	z_real;
				CLZNumber	z_imag;
				CLZNumber	two;
				CLZNumber	power_of_two;
				CLZNumber	nought;

				two.SetBase( l_base, false );
				two.Assign( 2L );
				nought.SetBase( l_base, false );
				nought.Assign( 0L );
				lzi_exact_power( &two, l_q, &power_of_two );
				lzi_complex_exact_power( &wide[ 0 ], &wide[ 1 ], l_q, &r_real, &r_imag );
				lzi_complex_exact_power( &a, &b, (l_p < 0) ? -l_p : l_p, &z_real, &z_imag );

				bool	b_holds = false;

				if ( l_p < 0 )
				{
					lzi_complex_times( &r_real, &r_imag, &z_real, &z_imag );
					b_holds = (lz_compare( &r_real, &power_of_two ) == kEqualTo) && (lz_compare( &r_imag, &nought ) == kEqualTo);
				}
				else
				{
					CLZNumber	want_real;
					CLZNumber	want_imag;

					lzi_complex_product( &z_real, &power_of_two, &want_real );
					lzi_complex_product( &z_imag, &power_of_two, &want_imag );
					b_holds = (lz_compare( &r_real, &want_real ) == kEqualTo) && (lz_compare( &r_imag, &want_imag ) == kEqualTo);
				}

				if ( b_holds )
				{
					b_proved = true;

					for ( long l_part = 0; l_part < 2; l_part++ )
						b_exact[ l_part ] = lzi_exact_quotient_at_places( &wide[ l_part ], &two, glWorkingPrecision, &parts[ l_part ] );
				}
			}
		}

		// Z = z^|p|, for the root's decider...
		CLZNumber	z_real;
		CLZNumber	z_imag;
		bool		b_rootDecider = (l_q >= 2);

		if ( b_rootDecider )
			lzi_complex_exact_power( &a, &b, (l_p < 0) ? -l_p : l_p, &z_real, &z_imag );

		for ( long l_part = 0; (l_part < 2) && !b_proved; l_part++ )
		{
			if ( b_nought[ l_part ] )
			{
				parts[ l_part ].SetBase( l_base, false );
				parts[ l_part ].Assign( 0L );
				b_exact[ l_part ] = true;
				continue;
			}

			/*
				**A part on an axis, rational where the mathematics makes it
				so**: |z|^E times an exact sine or cosine -- (-8)^(1/3)'s real
				part is 2 times a half -- or two irrationals -- the sixth root of
				-27 is sqrt(3) (cos 30 + i sin 30), its real part exactly 3/2, a
				half in base three. Where the squared cosine or sine of the
				phase is rational -- a multiple of 30 or 45 degrees, which every
				exact sine and cosine is -- a part c = M cos phase is rational
				exactly when, 2c taken from the work, (2c)^2q is 2^2q |z|^2p
				(cos^2)^q; and then it is rounded exactly.
			*/
			if ( b_degreePhase && b_realExponent && (l_q >= 1) && !(parts[ l_part ].IsZero()) )
			{
				long		l_degrees = 0;
				long		l_n = 0;
				long		l_d = 0;
				CLZNumber	turn;
				CLZNumber	reduced;

				turn.SetBase( l_base, false );
				turn.Assign( 360L );
				lz_remainder( &phase, &turn, &reduced );

				if ( lzi_whole_as_long( &reduced, &l_degrees ) )
				{
					if ( reduced.IsNegative() )
						l_degrees = -l_degrees;

					if ( l_part == 1 )
						l_degrees -= 90;			// sin^2 x is cos^2 (x - 90)

					if ( lzi_cosine_squared( l_degrees, &l_n, &l_d ) && (l_n != 0) )
					{
						CLZNumber	magnitude;
						CLZNumber	twice;
						CLZNumber	worked;
						long		l_places = glWorkingPrecision + l_guard - 4;
						long		l_was = glWorkingPrecision;

						magnitude.Assign( b_onReal ? &a : &b );
						magnitude.m_lSign = kPositive;
						worked.Assign( &parts[ l_part ] );
						glWorkingPrecision = l_places + 4;
						lz_add( &worked, &worked, &twice );
						twice.SetPrecision( l_places );
						twice.Round();
						twice.TrimTrailing();
						glWorkingPrecision = l_was;

						// (2c)^2q d^q [|z|^2|p| for p below nought] against
						//	2^2q n^q [|z|^2p otherwise]...
						CLZNumber	left;
						CLZNumber	right;
						CLZNumber	factor;
						CLZNumber	power;
						CLZNumber	product;
						CLZNumber	two;

						lzi_exact_power( &twice, 2 * l_q, &left );
						factor.SetBase( l_base, false );
						factor.Assign( l_d );
						lzi_exact_power( &factor, l_q, &power );
						lzi_complex_product( &left, &power, &product );
						left.Assign( &product );
						two.SetBase( l_base, false );
						two.Assign( 2L );
						lzi_exact_power( &two, 2 * l_q, &right );
						factor.Assign( l_n );
						lzi_exact_power( &factor, l_q, &power );
						lzi_complex_product( &right, &power, &product );
						right.Assign( &product );
						lzi_exact_power( &magnitude, 2 * ((l_p < 0) ? -l_p : l_p), &power );
						lzi_complex_product( (l_p < 0) ? &left : &right, &power, &product );
						((l_p < 0) ? left : right).Assign( &product );

						if ( lz_compare( &left, &right ) == kEqualTo )
						{
							b_exact[ l_part ] = lzi_exact_quotient_at_places( &twice, &two, glWorkingPrecision, &parts[ l_part ] );
							continue;
						}
					}
				}
			}

			SLZExactPart	e = exact[ l_part ];
			long			l_this = l_part;

			auto	tieCall = [ & ]( CLZNumber* pResult, long lPlaces, long lMostDigits ) -> bool
			{
				if ( b_degreePhase && e.m_bExact && (l_side != 0) && lzi_complex_trig_decide( pResult, e, lPlaces, l_side > 0 ) )
					return true;

				if ( unit_root[ l_this ].m_bExact && lzi_complex_trig_decide( pResult, unit_root[ l_this ], lPlaces, l_unitRootSide > 0 ) )
					return true;

				if ( !b_rootDecider )
					return false;

				CLZNumber	other;

				fnParts( glWorkingPrecision + l_guard, (l_this == 0) ? nullptr : &other, (l_this == 0) ? &other : nullptr );

				return lzi_complex_root_decide( pResult, l_this, &other, &z_real, &z_imag, (l_p < 0) ? -l_q : l_q, lPlaces,
												lMostDigits );
			};

			lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth )
									  {
										  return fnParts( lWidth, (l_this == 0) ? pValue : nullptr, (l_this == 1) ? pValue : nullptr );
									  },
									  tieCall, l_base, l_guard, &parts[ l_part ] );
		}
	}

	lResult = lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 0 ], &a, &pTarget->m_real, false, b_exact[ 0 ], b_exact[ 0 ] ) );

	return lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 1 ], &a, &pTarget->m_imag, false, b_exact[ 1 ], b_exact[ 1 ] ) );
}

// A complex value's two parts, refused where they or the operand's are in
//	two bases, and a flagged one carried: a not-a-number in either answers a
//	not-a-number, and so does an infinity, the point at infinity having no
//	direction to raise...
static bool lzi_complex_power_refused( CLZComplex* pLeft, CLZComplex* pRight, CLZComplex* pTarget, long* plResult )
{
	if ( lzi_complex_refused( pLeft, pRight, pTarget, plResult ) )
		return true;

	bool	b_flagged = lzi_complex_has( pLeft, kfQNaN ) || lzi_complex_has( pLeft, kfInfinity ) ||
						((pRight != nullptr) && (lzi_complex_has( pRight, kfQNaN ) || lzi_complex_has( pRight, kfInfinity )));

	if ( !b_flagged )
		return false;

	bool	b_divide = lzi_complex_has( pLeft, kfDiv0 ) || ((pRight != nullptr) && lzi_complex_has( pRight, kfDiv0 ));

	lzi_complex_flagged( pTarget, pLeft->m_real.m_lBase, kfQNaN, b_divide );
	*plResult = kNoError;

	return true;
}

// A real answer into the real part, an exact nought beside it; refused as
//	infinite, the point at infinity with the real function's code...
static long lzi_complex_real_answer( long lResult, CLZNumber* pOperand, CLZComplex* pTarget )
{
	if ( pTarget->m_real.FlagIsOn( kfInfinity ) )
	{
		lzi_complex_flagged( pTarget, pOperand->m_lBase, kfInfinity, pTarget->m_real.FlagIsOn( kfDiv0 ) );

		return lResult;
	}

	if ( pTarget->m_real.FlagIsOn( kfQNaN ) )
	{
		lzi_complex_flagged( pTarget, pOperand->m_lBase, kfQNaN, pTarget->m_real.FlagIsOn( kfDiv0 ) );

		return lResult;
	}

	CLZNumber	nought;

	nought.SetBase( pOperand->m_lBase, false );
	nought.Assign( 0L );

	return lzi_first_failure( lResult, lzi_answer_as_function( &nought, pOperand, &pTarget->m_imag, false, true, true ) );
}

// One, exactly, an exact nought beside it...
static long lzi_complex_one( CLZNumber* pOperand, CLZComplex* pTarget )
{
	CLZNumber	one;
	CLZNumber	nought;

	one.SetBase( pOperand->m_lBase, false );
	one.Assign( 1L );
	nought.SetBase( pOperand->m_lBase, false );
	nought.Assign( 0L );

	long	lResult = lzi_answer_as_function( &one, pOperand, &pTarget->m_real, false, true, true );

	return lzi_first_failure( lResult, lzi_answer_as_function( &nought, pOperand, &pTarget->m_imag, false, true, true ) );
}

// Nought to a complex power E: nought where Re E is past nought, the point
//	at infinity where it is below -- one over nought -- and a not-a-number
//	where it is nought and E is not, which has no limit...
static long lzi_complex_nought_to( CLZNumber* pReal, CLZNumber* pOperand, CLZComplex* pTarget )
{
	long	l_base = pOperand->m_lBase;

	pReal->CheckForZero();

	if ( pReal->IsZero() )
	{
		lzi_complex_flagged( pTarget, l_base, kfQNaN, false );

		return kFailed;
	}

	if ( pReal->IsNegative() )
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, true );

		return kDivideByZero;
	}

	CLZNumber	nought;

	nought.SetBase( l_base, false );
	nought.Assign( 0L );

	long	lResult = lzi_answer_as_function( &nought, pOperand, &pTarget->m_real, false, true, true );

	return lzi_first_failure( lResult, lzi_answer_as_function( &nought, pOperand, &pTarget->m_imag, false, true, true ) );
}

/*
	**z^w, the principal value e^(w ln z)** -- stage 5b, 2026-10-22, w complex,
	Rich's call. Each part correctly rounded, half-even, in places, and
	exact where the mathematics is. On the real axis, a base not below
	nought or a whole exponent, the real power's answer to the character,
	an exact nought beside it; a negative base to any other real power the
	principal value, as Elzed's ^ answers it. A whole w is worked exactly --
	see lzi_complex_whole_power. Nought to w is nought where Re w is past
	nought, the point at infinity where it is below, and refused where it
	is nought; anything to nought is one.
*/

long lzc_power( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_power_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	a;
	CLZNumber	b;
	CLZComplex	w;

	(void)l_base;
	a.Assign( &leftOperand->m_real );
	b.Assign( &leftOperand->m_imag );
	w.m_real.Assign( &rightOperand->m_real );
	w.m_imag.Assign( &rightOperand->m_imag );
	a.CheckForZero();
	b.CheckForZero();
	w.m_real.CheckForZero();
	w.m_imag.CheckForZero();

	// Anything to nought is one...
	if ( w.m_real.IsZero() && w.m_imag.IsZero() )
		return lzi_complex_one( &a, p_target );

	// On the real axis, the real power where it answers...
	if ( b.IsZero() && w.m_imag.IsZero() && (!a.IsNegative() || w.m_real.IsInteger()) )
		return lzi_complex_real_answer( lz_power( &a, &w.m_real, &p_target->m_real ), &a, p_target );

	if ( a.IsZero() && b.IsZero() )
		return lzi_complex_nought_to( &w.m_real, &a, p_target );

	// One to any power is one...
	if ( b.IsZero() && lzi_exactly( &a, 1 ) )
		return lzi_complex_one( &a, p_target );

	// A whole power, exactly...
	long	l_power = 0;

	if ( w.m_imag.IsZero() && w.m_real.IsInteger() && lzi_whole_as_long( &w.m_real, &l_power ) )
	{
		if ( w.m_real.IsNegative() )
			l_power = -l_power;

		return lzi_complex_whole_power( &a, &b, l_power, p_target );
	}

	return lzi_complex_power_general( &a, &b, &w, false, p_target );
}

/*
	**z^n for a whole n**, as lz_raise_n takes one -- stage 5b. Worked
	exactly, each part rounded once, short where exact; a fraction refused,
	kFailed, as the real one refuses it. On the real axis the real
	function's answer to the character, an exact nought beside it.
*/

long lzc_raise_n( CLZComplex* leftOperand, CLZNumber* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( rightOperand == nullptr )
		return kFailed;

	if ( lzi_complex_power_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( rightOperand->FlagIsOn( kfQNaN ) || rightOperand->FlagIsOn( kfInfinity ) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, rightOperand->FlagIsOn( kfDiv0 ) );

		return kNoError;
	}

	if ( !(rightOperand->IsInteger()) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, false );

		return kFailed;
	}

	CLZNumber	a;
	CLZNumber	b;
	CLZNumber	n;

	a.Assign( &leftOperand->m_real );
	b.Assign( &leftOperand->m_imag );
	n.Assign( rightOperand );
	a.CheckForZero();
	b.CheckForZero();
	n.CheckForZero();

	if ( b.IsZero() )
		return lzi_complex_real_answer( lz_raise_n( &a, &n, &p_target->m_real ), &a, p_target );

	long	l_power = 0;

	if ( lzi_whole_as_long( &n, &l_power ) )
	{
		if ( n.IsNegative() )
			l_power = -l_power;

		return lzi_complex_whole_power( &a, &b, l_power, p_target );
	}

	CLZComplex	w;

	w.m_real.Assign( &n );
	w.m_imag.SetBase( n.m_lBase, false );
	w.m_imag.Assign( 0L );

	return lzi_complex_power_general( &a, &b, &w, false, p_target );
}

/*
	**The principal root of complex order w, z^(1/w)** -- stage 5b. On the
	real axis with a real order, a value not below nought or an odd whole
	order, the real root's answer to the character -- the cube root of -8
	is -2, Rich's call, as Elzed's root answers -- and the principal value
	otherwise, so the square root of -4 is exactly 2i. An order whose
	reciprocal is whole is a whole power, worked exactly. Order nought is
	refused, kFailed, as the real root refuses it.
*/

static long lzi_complex_root( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* p_target, bool bCube )
{
	long	lResult = kNoError;
	long	l_base = leftOperand->m_real.m_lBase;

	CLZNumber	a;
	CLZNumber	b;
	CLZComplex	w;

	a.Assign( &leftOperand->m_real );
	b.Assign( &leftOperand->m_imag );
	w.m_real.Assign( &rightOperand->m_real );
	w.m_imag.Assign( &rightOperand->m_imag );
	a.CheckForZero();
	b.CheckForZero();
	w.m_real.CheckForZero();
	w.m_imag.CheckForZero();

	if ( w.m_real.IsZero() && w.m_imag.IsZero() )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, false );

		return kFailed;
	}

	// On the real axis, the real root where it answers...
	if ( b.IsZero() && w.m_imag.IsZero() )
	{
		bool	b_odd = false;

		if ( w.m_real.IsInteger() )
		{
			CLZNumber	two;
			CLZNumber	rest;

			two.SetBase( w.m_real.m_lBase, false );
			two.Assign( 2L );
			lz_remainder( &w.m_real, &two, &rest );
			rest.CheckForZero();
			b_odd = !(rest.IsZero());
		}

		if ( !a.IsNegative() || b_odd )
		{
			if ( bCube )
				lResult = lz_root_3( &a, &p_target->m_real );
			else
				lResult = lz_root_n( &a, &w.m_real, &p_target->m_real );

			return lzi_complex_real_answer( lResult, &a, p_target );
		}
	}

	if ( a.IsZero() && b.IsZero() )
		return lzi_complex_nought_to( &w.m_real, &a, p_target );

	if ( b.IsZero() && lzi_exactly( &a, 1 ) )
		return lzi_complex_one( &a, p_target );

	// An order of plus or minus one over a whole number is a whole power...
	if ( w.m_imag.IsZero() )
	{
		CLZNumber	n;
		CLZNumber	d;
		long		l_n = 0;
		long		l_d = 0;

		lzi_as_whole_ratio( &w.m_real, l_base, &n, &d );

		if ( lzi_whole_as_long( &n, &l_n ) && lzi_whole_as_long( &d, &l_d ) && (l_n != 0) && ((l_d % l_n) == 0) )
		{
			long	l_power = l_d / l_n;

			if ( n.IsNegative() )
				l_power = -l_power;

			return lzi_complex_whole_power( &a, &b, l_power, p_target );
		}
	}

	return lzi_complex_power_general( &a, &b, &w, true, p_target );
}

long lzc_root_n( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_power_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	return lzi_complex_root( leftOperand, rightOperand, p_target, false );
}

/*
	**The cube root**: on the real axis lz_root_3's answer to the character,
	the cube root of -8 being -2, and the principal value off it, each part
	in places -- stage 5b.
*/

long lzc_root_3( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_power_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	CLZComplex	order;

	order.m_real.SetBase( leftOperand->m_real.m_lBase, false );
	order.m_real.Assign( 3L );
	order.m_imag.SetBase( leftOperand->m_real.m_lBase, false );
	order.m_imag.Assign( 0L );

	return lzi_complex_root( leftOperand, &order, p_target, true );
}

/*
	**The base's exponential, b^z** -- b the base z is written in, as lz_exp_b
	raises it -- stage 5b. On the real axis lz_exp_b's answer to the
	character; off it b^x (cos(y ln b) + i sin(y ln b)), neither part
	rational (Gelfond-Schneider), so none lies on a half.
*/

long lzc_exp_b( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_power_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long		l_base = leftOperand->m_real.m_lBase;
	CLZComplex	z;

	z.m_real.Assign( &leftOperand->m_real );
	z.m_imag.Assign( &leftOperand->m_imag );
	z.m_real.CheckForZero();
	z.m_imag.CheckForZero();

	if ( z.m_imag.IsZero() )
		return lzi_complex_real_answer( lz_exp_b( &z.m_real, &p_target->m_real ), &z.m_real, p_target );

	CLZNumber	base;
	CLZNumber	nought;

	base.SetBase( l_base, false );
	base.Assign( l_base );
	nought.SetBase( l_base, false );
	nought.Assign( 0L );

	return lzi_complex_power_general( &base, &nought, &z, false, p_target );
}


/*
	**z^q against n^p, exactly** -- whether p/q is the logarithm of z to the
	base n, both parts compared: z^q n^-p against one where p is below
	nought. The branch is the caller's: an answer that rounds to p/q lies
	within its error of the principal one, and another branch lies 2 pi / |q
	ln n| away.
*/
static bool lzi_complex_power_is_exactly( CLZNumber* pZr, CLZNumber* pZi, CLZNumber* pNr, CLZNumber* pNi, long lP, long lQ )
{
	CLZNumber	left_real;
	CLZNumber	left_imag;
	CLZNumber	right_real;
	CLZNumber	right_imag;

	lzi_complex_exact_power( pZr, pZi, lQ, &left_real, &left_imag );
	lzi_complex_exact_power( pNr, pNi, (lP < 0) ? -lP : lP, &right_real, &right_imag );

	if ( lP < 0 )
	{
		lzi_complex_times( &left_real, &left_imag, &right_real, &right_imag );
		right_real.SetBase( pZr->m_lBase, false );
		right_real.Assign( 1L );
		right_imag.SetBase( pZr->m_lBase, false );
		right_imag.Assign( 0L );
	}

	return (lz_compare( &left_real, &right_real ) == kEqualTo) && (lz_compare( &left_imag, &right_imag ) == kEqualTo);
}

/*
	**The logarithm of z to the base n, ln z / ln n** -- stage 5b,
	2026-10-22, n complex, Rich's call: the principal value, each part
	correctly rounded, half-even, in places. bBaseOfValue is log_b's, n
	the base the value is written in and the real axis lz_log_b's.

	- **On the real axis, n and z both past nought, the real function's
	  answer to the character**, an exact nought beside it.
	- **A real base past nought**: the real part is half the logarithm of
	  |z|^2 to n -- the real lz_log_n to the base n^2, so it is exact where
	  that is and its ties are decided as that decides them -- and on an
	  axis the real logarithm of |z|, lz_log_b's for log_b, its spelling
	  included. The imaginary part is arg z / ln n, which is no fraction:
	  nought on the positive real axis, transcendental elsewhere.
	- **Any other base**: ln z times the conjugate of ln n over |ln n|^2,
	  wider by the quotient's whole digits and by |ln n|^2's noughts for a
	  base near one. An answer short of the width is raised back -- z^q
	  against n^p -- before it is spelled short: log_i(-1) is exactly 2.
	  For |n| one a part is exactly nought where ln z makes it so -- the
	  real for z past nought on the real axis, the imaginary for |z| one
	  -- and spelled short: log_i 9 is 0 - 1.398...i.

	Nought is refused as ln 0 is, kUndefinedResult, the point at infinity,
	and so is base nought; base one is one over nought, kDivideByZero.
*/

static long lzi_complex_log( CLZComplex* pBase, CLZComplex* pValue, bool bBaseOfValue, CLZComplex* pTarget )
{
	long	l_base = pBase->m_real.m_lBase;
	long	lResult = kNoError;

	CLZNumber	nr;
	CLZNumber	ni;
	CLZNumber	zr;
	CLZNumber	zi;

	nr.Assign( &pBase->m_real );
	ni.Assign( &pBase->m_imag );
	zr.Assign( &pValue->m_real );
	zi.Assign( &pValue->m_imag );

	// A value in another base is brought into the base's, exactly where it
	//	ends there...
	if ( zr.m_lBase != l_base )
	{
		long	l_keep = glWorkingPrecision + lzi_guard_places( kComplexGuardDigits, l_base ) + 8;

		lzi_into_base_keeping( &zr, l_base, l_keep );
		lzi_into_base_keeping( &zi, l_base, l_keep );
	}

	nr.CheckForZero();
	ni.CheckForZero();
	zr.CheckForZero();
	zi.CheckForZero();

	bool	b_realBase = ni.IsZero() && nr.IsPositive();

	// On the real axis, the real function...
	if ( zi.IsZero() && b_realBase && !zr.IsNegative() && !lzi_exactly( &nr, 1 ) )
	{
		if ( bBaseOfValue )
			lResult = lz_log_b( &zr, &pTarget->m_real );
		else
			lResult = lz_log_n( &nr, &zr, &pTarget->m_real );

		return lzi_complex_real_answer( lResult, &zr, pTarget );
	}

	// Nought, base nought and base one...
	if ( zr.IsZero() && zi.IsZero() )
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, false );

		return kUndefinedResult;
	}

	if ( nr.IsZero() && ni.IsZero() )
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, false );

		return kUndefinedResult;
	}

	if ( ni.IsZero() && lzi_exactly( &nr, 1 ) )
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, true );

		return kDivideByZero;
	}

	CLZNumber	norm;
	CLZNumber	parts[ 2 ];
	bool		b_exact[ 2 ] = { false, false };

	lzi_complex_norm( &zr, &zi, &norm );

	// With a real base the real part is a real logarithm, called as the
	//	host's so it is retaken and decided as that function is, and
	//	written as it writes it -- the operands are copies...
	if ( b_realBase )
	{
		if ( zr.IsZero() || zi.IsZero() )
		{
			CLZNumber	magnitude;

			magnitude.Assign( zr.IsZero() ? &zi : &zr );
			magnitude.m_lSign = kPositive;

			if ( bBaseOfValue )
				lResult = lz_log_b( &magnitude, &pTarget->m_real );
			else
				lResult = lz_log_n( &nr, &magnitude, &pTarget->m_real );
		}
		else
		{
			CLZNumber	square;

			lzi_complex_product( &nr, &nr, &square );
			lResult = lz_log_n( &square, &norm, &pTarget->m_real );
		}
	}

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( &nr, &scratch );

		work.Exact();

		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base );

		if ( b_realBase )
		{
			// ...and the imaginary, arg z over ln n...
			auto	fnImag = [ & ]( CLZNumber* pValue, long lWidth ) -> long
			{
				long	l_saved = glWorkingPrecision;
				long	l_result = kNoError;

				glWorkingPrecision = lWidth + 4;

				CLZNumber	angle;
				CLZNumber	logarithm;

				{
					CLZRadiansHere	radians;

					l_result = lz_arctangent2( &zi, &zr, &angle );
				}

				l_result = lzi_first_failure( l_result, lz_log_e( &nr, &logarithm ) );

				// A division answers at its target's precision, so the
				//	quotient is made at the width -- pValue was made at the
				//	answer's, and at no places divided to none...
				glWorkingPrecision = lWidth;

				CLZNumber	quotient;

				l_result = lzi_first_failure( l_result, lz_divide( &angle, &logarithm, &quotient ) );
				pValue->Assign( &quotient );
				glWorkingPrecision = l_saved;

				return l_result;
			};

			lResult = lzi_first_failure( lResult, fnImag( &parts[ 1 ], glWorkingPrecision + l_guard ) );
			lzi_complex_part_rounded( fnImag, l_base, l_guard, &parts[ 1 ] );
		}
		else
		{
			// ln z conj(ln n) / |ln n|^2, sized by doubles...
			double	d_zr = 0.0;
			double	d_zi = 0.0;
			double	d_nr = 0.0;
			double	d_ni = 0.0;
			double	d_lnNorm = 0.0;
			double	d_lnBase = 0.0;
			CLZNumber	base_norm;

			lzi_complex_norm( &nr, &ni, &base_norm );
			lzi_value_or_infinity( &zr, &d_zr );
			lzi_value_or_infinity( &zi, &d_zi );
			lzi_value_or_infinity( &nr, &d_nr );
			lzi_value_or_infinity( &ni, &d_ni );
			lzi_ln_estimate( &norm, &d_lnNorm );
			lzi_ln_estimate( &base_norm, &d_lnBase );

			double	d_lzr = 0.5 * d_lnNorm;
			double	d_lzi = std::atan2( d_zi, d_zr );
			double	d_lnr = 0.5 * d_lnBase;
			double	d_lni = std::atan2( d_ni, d_nr );
			double	d_den = (d_lnr * d_lnr) + (d_lni * d_lni);
			double	d_num = std::sqrt( (d_lzr * d_lzr) + (d_lzi * d_lzi) );
			long	l_answer = lzi_places_for( std::log( d_num / std::sqrt( d_den ) + 1.0 ), l_base );
			long	l_noughts = lzi_places_for( -std::log( d_den ), l_base );

			if ( l_answer < 0 )
				l_answer = 0;

			if ( l_noughts < 0 )
				l_noughts = 0;

			auto	fnParts = [ & ]( long lWidth, CLZNumber* pReal, CLZNumber* pImag ) -> long
			{
				long	l_saved = glWorkingPrecision;
				long	l_result = kNoError;

				glWorkingPrecision = lWidth + l_answer + (2 * l_noughts) + 4;

				CLZNumber	two;
				CLZNumber	lzr;
				CLZNumber	lzi;
				CLZNumber	lnr;
				CLZNumber	lni;

				two.SetBase( l_base, false );
				two.Assign( 2L );
				l_result = lzi_first_failure( l_result, lz_log_e( &norm, &lzr ) );
				lz_divide( &lzr, &two );
				l_result = lzi_first_failure( l_result, lz_log_e( &base_norm, &lnr ) );
				lz_divide( &lnr, &two );

				{
					CLZRadiansHere	radians;

					l_result = lzi_first_failure( l_result, lz_arctangent2( &zi, &zr, &lzi ) );
					l_result = lzi_first_failure( l_result, lz_arctangent2( &ni, &nr, &lni ) );
				}

				CLZNumber	den;
				CLZNumber	square;
				CLZNumber	ac;
				CLZNumber	bd;
				CLZNumber	bc;
				CLZNumber	ad;
				CLZNumber	top;

				lz_multiply( &lnr, &lnr, &den );
				lz_multiply( &lni, &lni, &square );
				lz_add( &den, &square );
				lz_multiply( &lzr, &lnr, &ac );
				lz_multiply( &lzi, &lni, &bd );
				lz_multiply( &lzi, &lnr, &bc );
				lz_multiply( &lzr, &lni, &ad );

				glWorkingPrecision = lWidth;

				CLZNumber	quotient;

				if ( pReal != nullptr )
				{
					lz_add( &ac, &bd, &top );
					lz_divide( &top, &den, &quotient );
					pReal->Assign( &quotient );
				}

				if ( pImag != nullptr )
				{
					lz_subtract( &bc, &ad, &top );
					lz_divide( &top, &den, &quotient );
					pImag->Assign( &quotient );
				}

				glWorkingPrecision = l_saved;

				return l_result;
			};

			lResult = fnParts( glWorkingPrecision + l_guard, &parts[ 0 ], &parts[ 1 ] );

			CLZNumber	wide;			// the real part as worked, for the proof below

			wide.Assign( &parts[ 0 ] );

			// A part exactly nought where ln z and ln n make it so: |n| one
			//	puts ln n on the imaginary axis, so the real part is nought
			//	where z is past nought on the real axis, ln z real, and the
			//	imaginary part is nought where |z| is one too, ln z
			//	imaginary -- log_i 9 is 0 - 1.398...i, the nought short...
			CLZNumber	one;

			one.SetBase( l_base, false );
			one.Assign( 1L );

			bool	b_unitBase = (lz_compare( &base_norm, &one ) == kEqualTo);
			bool	b_nought[ 2 ] = { b_unitBase && zi.IsZero() && zr.IsPositive(), b_unitBase && (lz_compare( &norm, &one ) == kEqualTo) };

			for ( long l_part = 0; l_part < 2; l_part++ )
			{
				long	l_this = l_part;

				if ( b_nought[ l_part ] )
				{
					parts[ l_part ].Assign( 0L );
					b_exact[ l_part ] = true;

					continue;
				}

				// A real answer exactly on a half is settled by the proof below,
				//	which rounds the fraction itself...
				lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth )
										  {
											  return fnParts( lWidth, (l_this == 0) ? pValue : nullptr, (l_this == 1) ? pValue : nullptr );
										  },
										  l_base, l_guard, &parts[ l_part ] );
			}

			/*
				**A real answer p/q, proved**: where the imaginary part rounds
				to nought, the work tried against each q to 64 -- p the
				nearest whole number to q times it, a candidate only where q
				times it lies that near -- and z^q against n^p. Where it holds
				the imaginary part is exactly nought and the real part p/q,
				rounded exactly: short where it ends, a half to the even kept
				digit, and in an odd base a half that never ends --
				log_{rho^2}(rho) is one half at every width there.
			*/
			parts[ 1 ].CheckForZero();

			if ( parts[ 1 ].IsZero() )
			{
				long	l_was = glWorkingPrecision;
				bool	b_proved = false;

				glWorkingPrecision = wide.m_lRadix + 4;

				for ( long q = 1; (q <= 64) && !b_proved; q++ )
				{
					CLZNumber	factor;
					CLZNumber	scaled;
					CLZNumber	nearest;
					CLZNumber	gap;
					CLZNumber	limit;
					long		l_p = 0;

					factor.SetBase( l_base, false );
					factor.Assign( q );
					lzi_complex_product( &wide, &factor, &scaled );
					nearest.Assign( &scaled );
					nearest.Round( RM_round_half_even, 0 );
					lz_subtract( &scaled, &nearest, &gap );
					gap.m_lSign = gap.IsZero() ? kZero : kPositive;
					limit.SetBase( l_base, false );
					limit.Assign( 1L );
					lzi_scale_by_base_power( &limit, l_was );

					if ( (lz_compare( &gap, &limit ) == kGreaterThan) || !lzi_whole_as_long( &nearest, &l_p ) || (l_p > 4096) )
						continue;

					if ( nearest.IsNegative() )
						l_p = -l_p;

					if ( lzi_complex_power_is_exactly( &zr, &zi, &nr, &ni, l_p, q ) )
					{
						CLZNumber	numerator;
						CLZNumber	denominator;

						numerator.SetBase( l_base, false );
						numerator.Assign( l_p );
						denominator.SetBase( l_base, false );
						denominator.Assign( q );
						glWorkingPrecision = l_was;
						b_exact[ 0 ] = lzi_exact_quotient_at_places( &numerator, &denominator, l_was, &parts[ 0 ] );
						b_exact[ 1 ] = true;
						parts[ 1 ].Assign( 0L );
						b_proved = true;
					}
				}

				glWorkingPrecision = l_was;
			}
		}
	}

	if ( !b_realBase )
		lResult = lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 0 ], &nr, &pTarget->m_real, false, b_exact[ 0 ], b_exact[ 0 ] ) );

	return lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 1 ], &nr, &pTarget->m_imag, false, b_exact[ 1 ], b_exact[ 1 ] ) );
}

long lzc_log_n( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_power_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	return lzi_complex_log( leftOperand, rightOperand, false, p_target );
}

/*
	**The logarithm to the base the value is written in** -- lzc_log_n with n
	that base, and on the axes lz_log_b's answers, spelling and all.
*/

long lzc_log_b( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_power_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	CLZComplex	base;
	CLZComplex	value;

	base.m_real.SetBase( leftOperand->m_real.m_lBase, false );
	base.m_real.Assign( leftOperand->m_real.m_lBase );
	base.m_imag.SetBase( leftOperand->m_real.m_lBase, false );
	base.m_imag.Assign( 0L );
	value.m_real.Assign( &leftOperand->m_real );
	value.m_imag.Assign( &leftOperand->m_imag );

	return lzi_complex_log( &base, &value, true, p_target );
}


/////////////////////////////////////////////////////////////////
//	Complex numbers, stage 6a -- the inverses that take square roots,
//	2026-10-23

// The eight by index: the circular four, then their hyperbolic twins...
const long	kInverseArcSine = 0;
const long	kInverseArcCosine = 1;
const long	kInverseArcSecant = 2;
const long	kInverseArcCosecant = 3;
const long	kInverseArcHyperSine = 4;
const long	kInverseArcHyperCosine = 5;
const long	kInverseArcHyperSecant = 6;
const long	kInverseArcHyperCosecant = 7;
const long	kInverseHyperbolic = 4;

/*
	**The principal square root of an exact a + ib, to lWidth places** -- two
	for each part of an inverse on the general road, rounded by nothing
	after: m = |z|, t = sqrt((|a| + m)/2), and the root t + i b/2t for a at
	or above nought, |b|/2t + i sign(b) t below, as lzc_root_2 works it, so
	nothing cancels. A root of a small value keeps only half its places, so
	the caller widens lWidth by half its noughts. On the negative real axis
	the +0 side: sqrt(-4) is 2i.

	**Its own work is wider again by half the noughts of |z|.** The root
	of a value k places down is k/2 down, so lWidth places of it ask for
	lWidth + k/2 of the modulus and of (|a| + m)/2 -- without them,
	arccsch(b^-28 + i) in base fourteen took the root of (x^2, -x) with
	its modulus rounded to nought, and divided by the root that made.
*/
static void lzi_complex_root_at( CLZNumber* pA, CLZNumber* pB, long lWidth, CLZNumber* pReal, CLZNumber* pImag )
{
	long	l_saved = glWorkingPrecision;
	long	l_base = pA->m_lBase;
	double	d_a = 0.0;
	double	d_b = 0.0;
	long	l_noughts = 0;

	lzi_value_or_infinity( pA, &d_a );
	lzi_value_or_infinity( pB, &d_b );

	double	d_size = std::hypot( d_a, d_b );

	if ( (d_size > 0.0) && (d_size < 1.0) )
		l_noughts = lzi_places_for( -std::log( d_size ), l_base );

	glWorkingPrecision = lWidth + (l_noughts + 1) / 2;

	CLZNumber	a;
	CLZNumber	b;
	CLZNumber	A;
	CLZNumber	two;
	CLZNumber	t;

	a.Assign( pA );
	b.Assign( pB );
	a.CheckForZero();
	b.CheckForZero();
	A.Assign( &a );
	A.m_lSign = a.IsZero() ? kZero : kPositive;
	two.SetBase( l_base, false );
	two.Assign( 2L );

	if ( b.IsZero() )
	{
		CLZNumber	nought;

		nought.SetBase( l_base, false );
		nought.Assign( 0L );
		lz_root_2( &A, &t );
		pReal->Assign( a.IsNegative() ? &nought : &t );
		pImag->Assign( a.IsNegative() ? &t : &nought );
		glWorkingPrecision = l_saved;

		return;
	}

	CLZNumber	aa;
	CLZNumber	bb;
	CLZNumber	S;
	CLZNumber	modulus;
	CLZNumber	sum;
	CLZNumber	half;
	CLZNumber	twice;
	CLZNumber	B;
	CLZNumber	other;

	lzi_complex_product( &a, &a, &aa );
	lzi_complex_product( &b, &b, &bb );
	lz_add( &aa, &bb, &S );
	lz_root_2( &S, &modulus );
	lz_add( &A, &modulus, &sum );
	lz_divide( &sum, &two, &half );
	lz_root_2( &half, &t );
	lz_add( &t, &t, &twice );
	B.Assign( &b );
	B.m_lSign = kPositive;
	lz_divide( &B, &twice, &other );

	if ( b.IsNegative() )
		lz_negate( a.IsNegative() ? &t : &other );

	pReal->Assign( a.IsNegative() ? &other : &t );
	pImag->Assign( a.IsNegative() ? &t : &other );
	glWorkingPrecision = l_saved;
}

// An exact sum or difference, x +- y, into pResult...
static void lzi_exact_sum( CLZNumber* pX, CLZNumber* pY, bool bSubtract, CLZNumber* pResult )
{
	if ( bSubtract )
		lz_subtract( pX, pY, pResult );
	else
		lz_add( pX, pY, pResult );

	pResult->CheckForZero();
}

/*
	**A real inverse's radians written in the host's unit** -- the imaginary
	part of arcsin 2 in degrees is arccosh 2 times 180 / pi. fnRadians at
	the width, the conversion's three places -- the turn over 2pi reaches
	159 -- and the guard, times the turn over 2pi, rounded once and taken
	again wider in doubt. In radians no conversion is made: the part is the
	host's call, the real function's answer to the character.
*/
static long lzi_inverse_part( LZUnaryFunction fnRadians, CLZNumber* pValue, bool bNegate, AngleUnit eUnit, CLZNumber* pTarget )
{
	long	lResult = kNoError;
	long	l_turn = lzi_turn_of( eUnit );

	if ( l_turn == 0 )
	{
		{
			CLZRadiansHere	radians;

			lResult = fnRadians( pValue, pTarget );
		}

		if ( bNegate )
			lz_negate( pTarget );

		return lResult;
	}

	long		l_base = pValue->m_lBase;
	CLZNumber	part;

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( pValue, &scratch );

		work.Exact();

		auto	fnPart = [ & ]( CLZNumber* pOut, long lWidth ) -> long
		{
			long	l_was = glWorkingPrecision;

			glWorkingPrecision = lWidth + lzi_guard_places( 3, l_base );

			CLZNumber	theta;
			CLZNumber	twoPi( kTwoPi );
			CLZNumber	scaled;
			long		l_result = fnRadians( pValue, &theta );

			lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, glWorkingPrecision + 2 );
			lzi_times_whole( &theta, l_turn, &scaled );

			// A division answers at its target's precision...
			glWorkingPrecision = lWidth;

			CLZNumber	angle;

			lz_divide( &scaled, &twoPi, &angle );

			if ( bNegate )
				lz_negate( &angle );

			pOut->Assign( &angle );
			glWorkingPrecision = l_was;

			return l_result;
		};

		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base );

		lResult = fnPart( &part, glWorkingPrecision + l_guard );
		lzi_complex_part_rounded( fnPart, l_base, l_guard, &part );
	}

	return lzi_first_failure( lResult, lzi_answer_as_function( &part, pValue, pTarget, false ) );
}

/*
	**On an axis each part is a real function's** -- a table checked
	against cmath, and against what Elzed answers today on every cut, from
	the operand's +0 side: arcsin 2 = pi/2 + 1.317i, arcsec 0.5 = 1.317i,
	arccsc 0.5 = pi/2 - 1.317i, arcsech 2 = -1.047i, arccsch(-0.5i) = 1.317 +
	(pi/2)i. A reciprocal's cut on the real axis is taken from below in 1/z,
	z's +0 side being 1/z's -0. The circular four's angles are in the host's
	unit, exact at its quarter and half turns; their other parts are a real
	function's radians written in the unit -- see lzi_inverse_part. The
	hyperbolic four are in radians.
*/
static long lzi_complex_inverse_on_axis( long lFunction, CLZNumber* pX, CLZNumber* pY, AngleUnit eUnit, CLZComplex* pTarget )
{
	long	l_base = pX->m_lBase;
	long	lResult = kNoError;
	bool	b_real = pY->IsZero();

	CLZNumber	v;				// the axis coordinate
	CLZNumber	magnitude;
	CLZNumber	sign_one;		// +-1, its sign
	CLZNumber	one;
	CLZNumber	nought;
	CLZNumber	minus_one;

	v.Assign( b_real ? pX : pY );
	magnitude.Assign( &v );
	magnitude.m_lSign = v.IsZero() ? kZero : kPositive;
	one.SetBase( l_base, false );
	one.Assign( 1L );
	minus_one.SetBase( l_base, false );
	minus_one.Assign( -1L );
	nought.SetBase( l_base, false );
	nought.Assign( 0L );
	sign_one.Assign( v.IsNegative() ? &minus_one : &one );

	long	l_against_one = lz_compare( &magnitude, &one );
	bool	b_inside = (l_against_one != kGreaterThan);		// |v| <= 1
	bool	b_within = (l_against_one != kLessThan);		// |v| >= 1

	CLZNumber*	p_real = &pTarget->m_real;
	CLZNumber*	p_imag = &pTarget->m_imag;

	// A part that is the host's call -- in the host's unit for a circular
	//	angle -- or, hyperbolic, the same call in radians...
	auto	Host = [ & ]( LZUnaryFunction fn, CLZNumber* pArg, CLZNumber* pOut, bool bNegate ) -> long
	{
		long	l_result = fn( pArg, pOut );

		if ( bNegate )
			lz_negate( pOut );

		return l_result;
	};

	auto	HostRadians = [ & ]( LZUnaryFunction fn, CLZNumber* pArg, CLZNumber* pOut, bool bNegate ) -> long
	{
		CLZRadiansHere	radians;

		return Host( fn, pArg, pOut, bNegate );
	};

	auto	Nought = [ & ]( CLZNumber* pOut ) -> long
	{
		return lzi_answer_as_function( &nought, pX, pOut, false, true, true );
	};

	// A pole, or a real function's refusal, is the point at infinity with
	//	its code...
	auto	Refused = [ & ]( long lCode ) -> long
	{
		if ( p_real->FlagIsOn( kfInfinity ) || p_real->FlagIsOn( kfQNaN ) )
			lzi_complex_flagged( pTarget, l_base, p_real->FlagIsOn( kfInfinity ) ? kfInfinity : kfQNaN, p_real->FlagIsOn( kfDiv0 ) );

		return lCode;
	};

	switch ( lFunction )
	{
	case kInverseArcSine:
		if ( b_real )
		{
			if ( b_inside )
				return Refused( lzi_first_failure( Host( lz_arcsine, &v, p_real, false ), Nought( p_imag ) ) );

			lResult = Host( lz_arcsine, &sign_one, p_real, false );

			return lzi_first_failure( lResult, lzi_inverse_part( lz_archypercosine, &magnitude, false, eUnit, p_imag ) );
		}

		lResult = Nought( p_real );

		return lzi_first_failure( lResult, lzi_inverse_part( lz_archypersine, &v, false, eUnit, p_imag ) );

	case kInverseArcCosine:
		if ( b_real )
		{
			if ( b_inside )
				return Refused( lzi_first_failure( Host( lz_arccosine, &v, p_real, false ), Nought( p_imag ) ) );

			lResult = Host( lz_arccosine, &sign_one, p_real, false );

			return lzi_first_failure( lResult, lzi_inverse_part( lz_archypercosine, &magnitude, true, eUnit, p_imag ) );
		}

		lResult = Host( lz_arccosine, &nought, p_real, false );

		return lzi_first_failure( lResult, lzi_inverse_part( lz_archypersine, &v, true, eUnit, p_imag ) );

	case kInverseArcSecant:
	case kInverseArcCosecant:
	{
		bool	b_secant = (lFunction == kInverseArcSecant);

		if ( b_real )
		{
			if ( v.IsZero() )
			{
				lzi_complex_flagged( pTarget, l_base, kfInfinity, false );

				return kUndefinedResult;
			}

			if ( b_within )
				return Refused( lzi_first_failure( Host( b_secant ? lz_arcsecant : lz_arccosecant, &v, p_real, false ), Nought( p_imag ) ) );

			lResult = Host( b_secant ? lz_arcsecant : lz_arccosecant, &sign_one, p_real, false );

			return lzi_first_failure( lResult, lzi_inverse_part( lz_archypersecant, &magnitude, !b_secant, eUnit, p_imag ) );
		}

		lResult = b_secant ? Host( lz_arccosine, &nought, p_real, false ) : Nought( p_real );

		return lzi_first_failure( lResult, lzi_inverse_part( lz_archypercosecant, &v, !b_secant, eUnit, p_imag ) );
	}

	case kInverseArcHyperSine:
		if ( b_real )
			return Refused( lzi_first_failure( Host( lz_archypersine, &v, p_real, false ), Nought( p_imag ) ) );

		if ( b_inside )
			return lzi_first_failure( Nought( p_real ), HostRadians( lz_arcsine, &v, p_imag, false ) );

		lResult = Host( lz_archypercosine, &magnitude, p_real, false );

		return lzi_first_failure( lResult, HostRadians( lz_arcsine, &sign_one, p_imag, false ) );

	case kInverseArcHyperCosine:
		if ( b_real )
		{
			if ( lz_compare( &v, &one ) != kLessThan )
				return Refused( lzi_first_failure( Host( lz_archypercosine, &v, p_real, false ), Nought( p_imag ) ) );

			if ( lz_compare( &v, &minus_one ) != kLessThan )
				return lzi_first_failure( Nought( p_real ), HostRadians( lz_arccosine, &v, p_imag, false ) );

			lResult = Host( lz_archypercosine, &magnitude, p_real, false );

			return lzi_first_failure( lResult, HostRadians( lz_arccosine, &minus_one, p_imag, false ) );
		}

		lResult = Host( lz_archypersine, &magnitude, p_real, false );

		return lzi_first_failure( lResult, HostRadians( lz_arcsine, &sign_one, p_imag, false ) );

	case kInverseArcHyperSecant:
		if ( b_real )
		{
			if ( v.IsZero() || (!v.IsNegative() && b_inside) )
				return Refused( lzi_first_failure( Host( lz_archypersecant, &v, p_real, false ), Nought( p_imag ) ) );

			if ( v.IsNegative() && b_inside )
			{
				lResult = Host( lz_archypersecant, &magnitude, p_real, false );

				return lzi_first_failure( lResult, HostRadians( lz_arccosine, &minus_one, p_imag, true ) );
			}

			return lzi_first_failure( Nought( p_real ), HostRadians( lz_arcsecant, &v, p_imag, true ) );
		}

		lResult = Host( lz_archypercosecant, &magnitude, p_real, false );

		return lzi_first_failure( lResult, HostRadians( lz_arcsine, &sign_one, p_imag, true ) );

	default:	// kInverseArcHyperCosecant
		if ( b_real )
			return Refused( lzi_first_failure( Host( lz_archypercosecant, &v, p_real, false ), Nought( p_imag ) ) );

		if ( b_within )
			return lzi_first_failure( Nought( p_real ), HostRadians( lz_arccosecant, &v, p_imag, true ) );

		lResult = Host( lz_archypersecant, &magnitude, p_real, false );

		return lzi_first_failure( lResult, HostRadians( lz_arcsine, &sign_one, p_imag, true ) );
	}
}

/*
	**A part as a real inverse of an algebraic number** -- for the tie
	decider below. Off the axes every part of the eight is sigma F(gamma),
	F one of asin, acos and acosh and gamma = (sqrt p +- sqrt q) / (2 sqrt s),
	p, q and s exact: Re arcsin z is asin of beta and Im arcsin z sign(y)
	acosh of alpha, beta and alpha those with p and q the squares of |z +- 1|
	(Abramowitz and Stegun 4.4.37), and the rest follow -- 1/z puts s = |z|^2,
	iz and i/z new p and q. Checked against cmath, sixteen parts.
*/
struct SLZInversePart
{
	long		m_lFunction;		// 0 asin, 1 acos, 2 acosh
	long		m_lSign;			// sigma, +-1
	bool		m_bPlus;			// sqrt p + sqrt q, or minus
	CLZNumber	m_p;
	CLZNumber	m_q;
	CLZNumber	m_s;
};

// The sign of (sqrt A +- sqrt B) - (rhs sign) sqrt C, A, B and C exact and
//	not below nought, by squaring: nought only where they are equal...
static long lzi_root_pair_sign( CLZNumber* pA, CLZNumber* pB, bool bPlus, CLZNumber* pC, long lRhsSign )
{
	CLZNumber	A;
	CLZNumber	B;

	A.Assign( pA );
	B.Assign( pB );

	long	l_left = bPlus ? 1 : lz_compare( &A, &B );			// the sign of the left side

	if ( l_left == kGreaterThan )
		l_left = 1;
	else if ( l_left == kLessThan )
		l_left = -1;
	else if ( !bPlus )
		l_left = 0;

	pC->CheckForZero();

	if ( pC->IsZero() )
		lRhsSign = 0;

	// Opposite signs, or one side nought, settle it at once...
	if ( (l_left >= 0) && (lRhsSign <= 0) )
		return ((l_left == 0) && (lRhsSign == 0)) ? 0 : 1;

	if ( (l_left <= 0) && (lRhsSign >= 0) )
		return ((l_left == 0) && (lRhsSign == 0)) ? 0 : -1;

	// ...both positive, or both negative and compared as magnitudes...
	bool	b_negative = (l_left < 0);

	if ( b_negative )
	{
		CLZNumber	swap;

		swap.Assign( &A );
		A.Assign( &B );
		B.Assign( &swap );
	}

	// |left|^2 = A + B +- 2 sqrt(AB) against C: D = A + B - C...
	CLZNumber	sum;
	CLZNumber	D;
	CLZNumber	AB;
	CLZNumber	four_AB;
	CLZNumber	four;
	CLZNumber	DD;

	lz_add( &A, &B, &sum );
	lz_subtract( &sum, pC, &D );
	D.CheckForZero();
	lzi_complex_product( &A, &B, &AB );
	four.SetBase( A.m_lBase, false );
	four.Assign( 4L );
	lzi_complex_product( &AB, &four, &four_AB );
	lzi_complex_product( &D, &D, &DD );

	long	l_magnitude = 0;

	if ( bPlus )
	{
		// D + 2 sqrt(AB): positive where D is not below nought...
		if ( !D.IsNegative() )
			l_magnitude = (D.IsZero() && AB.IsZero()) ? 0 : 1;
		else
		{
			long	l_order = lz_compare( &four_AB, &DD );

			l_magnitude = (l_order == kGreaterThan) ? 1 : ((l_order == kLessThan) ? -1 : 0);
		}
	}
	else
	{
		// D - 2 sqrt(AB): negative where D is not above nought...
		if ( !D.IsPositive() )
			l_magnitude = (D.IsZero() && AB.IsZero()) ? 0 : -1;
		else
		{
			long	l_order = lz_compare( &DD, &four_AB );

			l_magnitude = (l_order == kGreaterThan) ? 1 : ((l_order == kLessThan) ? -1 : 0);
		}
	}

	return b_negative ? -l_magnitude : l_magnitude;
}

// x^n exactly, n nought or more...
static void lzi_whole_power( CLZNumber* pX, long lPower, CLZNumber* pResult )
{
	CLZNumber	result;

	result.SetBase( pX->m_lBase, false );
	result.Assign( 1L );

	for ( long i = 0; i < lPower; i++ )
		lz_multiply( &result, pX );

	pResult->Assign( &result );
}

/*
	**sin, cos, cosh or exp of a rational q = Qn / Qd, as a rational
	interval** -- lKind 0, 1, 2, 3 -- from lTerms terms of the series over
	one denominator, Qd^N N!, N the remainder's power: the sum, and the
	bound on the rest, the first term left out. For sin and cos, whose
	terms alternate, that holds wherever the terms shrink from it on -- |q|
	below N; cosh's and exp's do not alternate, so their rest is bounded by
	twice that term, which holds only where each term is at most half the
	one before -- q^2 at most (N+1)(N+2)/2 for cosh, which the caller
	assures by N above 1.5 |q|, and q at most (N+1)/2 for exp, by N at
	least 2|q|. exp takes every power, N = lTerms, and wants q past
	nought. The answers are whole numbers over pDen.
*/
static void lzi_series_interval( long lKind, CLZNumber* pQn, CLZNumber* pQd, long lTerms,
								 CLZNumber* pLower, CLZNumber* pUpper, CLZNumber* pDen )
{
	long	l_base = pQn->m_lBase;
	bool	b_every = (lKind == 3);						// exp takes every power
	long	l_first = (lKind == 0) ? 1 : 0;				// the series' first power
	long	l_rest = b_every ? lTerms : l_first + 2 * lTerms;	// the remainder's power, N

	CLZNumber	sum;
	CLZNumber	factor;					// N! / m!, m running down
	CLZNumber	term;
	CLZNumber	qn_power;
	CLZNumber	qd_power;
	CLZNumber	rest;

	sum.SetBase( l_base, false );
	sum.Assign( 0L );
	factor.SetBase( l_base, false );
	factor.Assign( 1L );

	CLZNumber	Qn_abs;

	Qn_abs.Assign( pQn );
	Qn_abs.m_lSign = pQn->IsZero() ? kZero : kPositive;

	for ( long m = l_rest; m >= l_first; m-- )
	{
		if ( b_every || ((m - l_first) % 2 == 0) )
		{
			// Qn^m Qd^(N - m) N! / m!...
			lzi_whole_power( (m == l_rest) ? &Qn_abs : pQn, m, &qn_power );
			lzi_whole_power( pQd, l_rest - m, &qd_power );
			lzi_complex_product( &qn_power, &qd_power, &term );
			lzi_complex_product( &term, &factor, &term );

			if ( m == l_rest )
				rest.Assign( &term );
			else
			{
				bool	b_subtract = (lKind < 2) && (((m - l_first) / 2) % 2 == 1);

				if ( b_subtract )
					lz_subtract( &sum, &term, &sum );
				else
					lz_add( &sum, &term, &sum );
			}
		}

		CLZNumber	m_number;

		m_number.SetBase( l_base, false );
		m_number.Assign( (m > 1) ? m : 1L );
		lzi_complex_product( &factor, &m_number, &factor );
	}

	if ( lKind >= 2 )
		lz_add( &rest, &rest, &rest );

	lzi_whole_power( pQd, l_rest, &qd_power );
	lzi_complex_product( &qd_power, &factor, pDen );			// Qd^N N!
	lz_subtract( &sum, &rest, pLower );
	lz_add( &sum, &rest, pUpper );

	if ( lKind >= 2 )
		pLower->Assign( &sum );
}

/*
	**A part of the eight in doubt, decided exactly** -- in radians, where a
	part can lie a hair from a rational: arcsin z is z and a little near
	nought, arccsc z 1/z and a little far out, arcsinh's imaginary part y/x
	and a little along the real axis, and a leading term on a half puts the
	part nearer it than any retake reaches. Measured before this, 186 of
	1,021 such were rounded the wrong way.

	With h the half the part lies near, the part is above it exactly when
	F(gamma) is above h / sigma, and since F is monotone that is gamma
	against sin, cos or cosh of h / sigma -- a rational interval from the
	series, as narrow as its terms make it -- compared exactly, by
	lzi_root_pair_sign. A rational h other than nought has a
	transcendental sine, cosine and cosh, so gamma never equals one, and
	more terms always separate them; the budget lMostPlaces is the limit.
	Near the end of asin's or acos's range the interval's function stops
	being monotone, so a half there is left to the retake.
*/
static bool lzi_inverse_decide( CLZNumber* pResult, long lPlaces, long lMostPlaces, const SLZInversePart& part )
{
	if ( lPlaces < 0 )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pResult->m_lRadix > lPlaces) ? pResult->m_lRadix : lPlaces) + 4;

	// The half: |part| rounded down at the places, y, and h = y + u/2,
	//	Qn / Qd = sign (2 Y + 1) / (2 b^P), with Y = y b^P...
	CLZNumber	y;
	CLZNumber	Y;
	CLZNumber	Qn;
	CLZNumber	Qd;
	CLZNumber	one;
	CLZNumber	two;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );
	Y.Assign( &y );
	lzi_scale_by_base_power( &Y, -lPlaces );
	one.SetBase( l_base, false );
	one.Assign( 1L );
	two.SetBase( l_base, false );
	two.Assign( 2L );
	lzi_complex_product( &Y, &two, &Qn );
	lz_add( &Qn, &one, &Qn );
	Qd.Assign( &two );
	lzi_scale_by_base_power( &Qd, -lPlaces );

	long	l_sign_h = pResult->IsNegative() ? -1 : 1;

	// q = h / sigma...
	if ( l_sign_h * part.m_lSign < 0 )
		lz_negate( &Qn );

	double	d_qn = 0.0;
	double	d_qd = 1.0;

	lzi_value_or_infinity( &Qn, &d_qn );
	lzi_value_or_infinity( &Qd, &d_qd );

	double	d_q = d_qn / d_qd;
	long	l_versus = 0;				// the sign of F(gamma) - q, once known

	switch ( part.m_lFunction )
	{
	case 0:		// asin, increasing, (-pi/2, pi/2]
		if ( std::fabs( d_q ) > 1.581 )
			l_versus = (d_q > 0) ? -1 : 1;
		else if ( std::fabs( d_q ) > 1.56 )
		{
			glWorkingPrecision = l_saved;

			return false;
		}
		break;

	case 1:		// acos, decreasing, [0, pi]
		if ( d_q <= 0.0 )
			l_versus = 1;
		else if ( d_q > 3.152 )
			l_versus = -1;
		else if ( d_q > 3.13 )
		{
			glWorkingPrecision = l_saved;

			return false;
		}
		break;

	default:	// acosh, increasing, at or above nought
		if ( d_q <= 0.0 )
			l_versus = 1;
		break;
	}

	// The series, more terms each round, until gamma clears its interval...
	for ( long l_terms = 2; l_versus == 0; l_terms *= 2 )
	{
		double	d_digits = (2.0 * l_terms + 1.0) * (std::log( d_qd > 1.0 ? d_qd : 2.0 ) / std::log( (double)l_base ) + 1.0)
						   + std::lgamma( 2.0 * l_terms + 2.0 ) / std::log( (double)l_base );

		// The remainder's bound wants N above |q|, and for cosh above 1.5 |q|...
		double	d_need = (part.m_lFunction == 2) ? 1.5 * std::fabs( d_q ) + 1.0 : std::fabs( d_q );

		if ( (d_digits > (double)lMostPlaces) || ((2.0 * l_terms) <= d_need) )
		{
			if ( d_digits > (double)lMostPlaces )
			{
				glWorkingPrecision = l_saved;

				return false;
			}

			continue;
		}

		CLZNumber	lower;
		CLZNumber	upper;
		CLZNumber	den;

		lzi_series_interval( part.m_lFunction == 0 ? 0 : (part.m_lFunction == 1 ? 1 : 2), &Qn, &Qd, l_terms, &lower, &upper, &den );

		// gamma against c = num / den: (sqrt p +- sqrt q) den against
		//	2 num sqrt s, so A = den^2 p, B = den^2 q, C = 4 num^2 s...
		auto	Versus = [ & ]( CLZNumber* pNum ) -> long
		{
			CLZNumber	den2;
			CLZNumber	A;
			CLZNumber	B;
			CLZNumber	num2;
			CLZNumber	C;
			CLZNumber	four;
			CLZNumber	num_abs;

			lzi_complex_product( &den, &den, &den2 );
			lzi_complex_product( &den2, const_cast<CLZNumber*>( &part.m_p ), &A );
			lzi_complex_product( &den2, const_cast<CLZNumber*>( &part.m_q ), &B );
			lzi_complex_product( pNum, pNum, &num2 );
			four.SetBase( l_base, false );
			four.Assign( 4L );
			lzi_complex_product( &num2, &four, &C );
			lzi_complex_product( &C, const_cast<CLZNumber*>( &part.m_s ), &C );
			pNum->CheckForZero();

			return lzi_root_pair_sign( &A, &B, part.m_bPlus, &C, pNum->IsZero() ? 0 : (pNum->IsNegative() ? -1 : 1) );
		};

		long	l_low = Versus( &lower );
		long	l_high = Versus( &upper );

		// gamma above the whole interval, or below it...
		long	l_gamma = ((l_low > 0) && (l_high > 0)) ? 1 : (((l_low < 0) && (l_high < 0)) ? -1 : 0);

		if ( l_gamma != 0 )
			l_versus = (part.m_lFunction == 1) ? -l_gamma : l_gamma;
	}

	// sign(part - h) = sigma sign(F - q); the magnitude goes up where the
	//	part lies beyond h...
	bool		b_up = (part.m_lSign * l_versus * l_sign_h > 0);
	CLZNumber	unit;

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	if ( b_up )
		lz_add( &y, &unit );

	if ( l_sign_h < 0 )
		lz_negate( &y );

	y.CheckForZero();
	pResult->Assign( &y );
	glWorkingPrecision = l_saved;

	return true;
}

/*
	**Off both axes, Kahan's forms over two square roots** -- W. Kahan,
	"Branch Cuts for Complex Elementary Functions" (1987), checked against
	cmath at 20,000 points each. With S1 and S2 the principal roots of two
	exact values, D an exact divisor and N an exact numerator:

	| form | real part | imaginary part |
	|---|---|---|
	| arcsin | atan2(N, Re S1 S2) | arcsinh(Im(conj(S1) S2) / D) |
	| arccos | 2 atan2(Re S1, Re S2) | arcsinh(Im(conj(S2) S1) / D) |
	| arccosh | arcsinh(Re(conj(S1) S2) / D) | 2 atan2(Im S1, Re S2) |

	arcsin z takes 1 -+ z, N = x; arccos the same; arccosh z -+ 1. **The
	reciprocals never form 1/z**, which need not end in the base:
	arcsec and arccsc take (z -+ 1) conj z and D = |z|^2, arcsech (1 -+ z)
	conj z, each 1 -+ 1/z times |z|^2, the positive scale cancelling in
	atan2 and dividing out of arcsinh's argument. arcsinh is the arcsin
	form at iz turned back, arccsch at i/z, (|z|^2 -+ i conj z). Off the
	axes no root's argument lies on its cut, and none of these cancels.

	**Widened by half the noughts of the smallest root's argument**, which
	a root of it loses, and **by half those of D**, whose division magnifies
	what is left -- a reciprocal near nought. Each part taken again wider
	in doubt.

	**In a unit the circular four write both parts in it**, and a real part
	off the axes is exact at 30 and 45 degrees and nowhere else: it is the
	arcsine of beta, beta^2 = (p + q - 2 sqrt(pq)) / 4 with p and q the
	squares of |z +- 1| (1/z's for the reciprocals), in Q(sqrt(pq)); a
	Gaussian rational z makes beta^2 and its conjugate the roots of B^2 -
	(1 + |z|^2) B + x^2, which rules out every whole-degree sine of degree
	two, 15, 18, 22.5 and 36 among them, and 60 by a sum of two squares. So
	beta^2 is a quarter or a half, tested exactly.
*/
static long lzi_complex_inverse_general( long lFunction, CLZNumber* pX, CLZNumber* pY, AngleUnit eUnit, CLZComplex* pTarget )
{
	long	l_base = pX->m_lBase;
	long	lResult = kNoError;
	bool	b_hyperbolic = (lFunction >= kInverseHyperbolic);
	long	l_turn = b_hyperbolic ? 0 : lzi_turn_of( eUnit );

	CLZNumber	one;
	CLZNumber	xx;
	CLZNumber	yy;
	CLZNumber	n2;

	one.SetBase( l_base, false );
	one.Assign( 1L );
	lzi_complex_product( pX, pX, &xx );
	lzi_complex_product( pY, pY, &yy );
	lz_add( &xx, &yy, &n2 );

	// The form, the roots' exact arguments, N and D...
	enum ELZInverseForm { kFormArcSine, kFormArcCosine, kFormArcHyperCosine };

	ELZInverseForm	e_form = kFormArcSine;
	bool			b_turned = false;			// arcsinh and arccsch: the arcsin form turned back
	bool			b_divided = false;			// D = |z|^2 rather than one
	CLZNumber		u1r;
	CLZNumber		u1i;
	CLZNumber		u2r;
	CLZNumber		u2i;
	CLZNumber		N;
	CLZNumber		negated;

	switch ( lFunction )
	{
	case kInverseArcSine:
	case kInverseArcCosine:
		e_form = (lFunction == kInverseArcSine) ? kFormArcSine : kFormArcCosine;
		lzi_exact_sum( &one, pX, true, &u1r );			// 1 - z
		u1i.Assign( pY );
		lz_negate( &u1i );
		lzi_exact_sum( &one, pX, false, &u2r );			// 1 + z
		u2i.Assign( pY );
		N.Assign( pX );
		break;

	case kInverseArcSecant:
	case kInverseArcCosecant:
		e_form = (lFunction == kInverseArcCosecant) ? kFormArcSine : kFormArcCosine;
		b_divided = true;
		lzi_exact_sum( &n2, pX, true, &u1r );			// (z - 1) conj z
		u1i.Assign( pY );
		lzi_exact_sum( &n2, pX, false, &u2r );			// (z + 1) conj z
		u2i.Assign( pY );
		lz_negate( &u2i );
		N.Assign( pX );
		break;

	case kInverseArcHyperSine:
		b_turned = true;
		lzi_exact_sum( &one, pY, false, &u1r );			// 1 - iz
		u1i.Assign( pX );
		lz_negate( &u1i );
		lzi_exact_sum( &one, pY, true, &u2r );			// 1 + iz
		u2i.Assign( pX );
		N.Assign( pY );
		lz_negate( &N );
		break;

	case kInverseArcHyperCosecant:
		b_turned = true;
		b_divided = true;
		lzi_exact_sum( &n2, pY, true, &u1r );			// |z|^2 - i conj z
		u1i.Assign( pX );
		lz_negate( &u1i );
		lzi_exact_sum( &n2, pY, false, &u2r );			// |z|^2 + i conj z
		u2i.Assign( pX );
		N.Assign( pY );
		break;

	case kInverseArcHyperCosine:
		e_form = kFormArcHyperCosine;
		lzi_exact_sum( pX, &one, true, &u1r );			// z - 1
		u1i.Assign( pY );
		lzi_exact_sum( pX, &one, false, &u2r );			// z + 1
		u2i.Assign( pY );
		break;

	default:	// kInverseArcHyperSecant
		e_form = kFormArcHyperCosine;
		b_divided = true;
		lzi_exact_sum( pX, &n2, true, &u1r );			// (1 - z) conj z
		u1i.Assign( pY );
		lz_negate( &u1i );
		lzi_exact_sum( pX, &n2, false, &u2r );			// (1 + z) conj z
		u2i.Assign( pY );
		lz_negate( &u2i );
		break;
	}

	// The widening, from doubles: half the noughts of the smaller root's
	//	argument and of D...
	auto	Noughts = [ & ]( double dValue ) -> long
	{
		if ( !(dValue > 0.0) || (dValue >= 1.0) )
			return 0;

		return lzi_places_for( -std::log( dValue ), l_base );
	};

	double	d_u1r = 0.0;
	double	d_u1i = 0.0;
	double	d_u2r = 0.0;
	double	d_u2i = 0.0;
	double	d_n2 = 0.0;

	lzi_value_or_infinity( &u1r, &d_u1r );
	lzi_value_or_infinity( &u1i, &d_u1i );
	lzi_value_or_infinity( &u2r, &d_u2r );
	lzi_value_or_infinity( &u2i, &d_u2i );
	lzi_value_or_infinity( &n2, &d_n2 );

	double	d_first = std::hypot( d_u1r, d_u1i );
	double	d_second = std::hypot( d_u2r, d_u2i );
	double	d_smaller = (d_first < d_second) ? d_first : d_second;
	long	l_extra = (Noughts( d_smaller ) + 1) / 2 + (b_divided ? (Noughts( d_n2 ) + 1) / 2 : 0);

	if ( l_turn != 0 )
		l_extra += lzi_guard_places( 3, l_base );

	// In a unit, a real part off the axes exact at 30 or 45 degrees: beta^2
	//	a quarter where p + q - s = 2 sqrt(pq), a half where p + q - 2s does,
	//	s one or, for a reciprocal, |z|^2...
	bool		b_exact_real = false;
	CLZNumber	exact_real;

	if ( l_turn != 0 )
	{
		CLZNumber	xp1;
		CLZNumber	xm1;
		CLZNumber	p;
		CLZNumber	q;
		CLZNumber	square;
		CLZNumber	pq;
		CLZNumber	four_pq;
		CLZNumber	four;

		lzi_exact_sum( pX, &one, false, &xp1 );
		lzi_exact_sum( pX, &one, true, &xm1 );
		lzi_complex_product( &xp1, &xp1, &square );
		lz_add( &square, &yy, &p );
		lzi_complex_product( &xm1, &xm1, &square );
		lz_add( &square, &yy, &q );
		lzi_complex_product( &p, &q, &pq );
		four.SetBase( l_base, false );
		four.Assign( 4L );
		lzi_complex_product( &pq, &four, &four_pq );

		for ( long l_times = 1; (l_times <= 2) && !b_exact_real; l_times++ )
		{
			CLZNumber	s;
			CLZNumber	times;
			CLZNumber	sum;
			CLZNumber	left;

			times.SetBase( l_base, false );
			times.Assign( l_times );
			lzi_complex_product( b_divided ? &n2 : &one, &times, &s );
			lz_add( &p, &q, &sum );
			lz_subtract( &sum, &s, &left );
			left.CheckForZero();

			if ( left.IsNegative() )
				continue;

			lzi_complex_product( &left, &left, &square );

			if ( lz_compare( &square, &four_pq ) != kEqualTo )
				continue;

			// beta^2 a quarter: 30 degrees; a half: 45. beta's sign is x's...
			bool	b_arcsine = (lFunction == kInverseArcSine) || (lFunction == kInverseArcCosecant);
			long	l_a = 1;
			long	l_b = (l_times == 1) ? 12 : 8;

			if ( !b_arcsine )
			{
				// the arccosine of beta, a quarter turn less the arcsine...
				if ( l_times == 1 )
					l_b = pX->IsNegative() ? 3 : 6;
				else
					l_a = pX->IsNegative() ? 3 : 1;
			}
			else if ( pX->IsNegative() )
				l_a = -1;

			lzi_part_of_turn( l_turn, l_a, l_b, l_base, &exact_real );
			b_exact_real = true;
		}
	}

	// Each part as sigma F(gamma), for the tie decider, in radians -- see
	//	lzi_inverse_decide and SLZInversePart...
	SLZInversePart	spec[ 2 ];

	if ( l_turn == 0 )
	{
		CLZNumber	xp1;
		CLZNumber	xm1;
		CLZNumber	square;
		CLZNumber	p;
		CLZNumber	q;
		CLZNumber	n4;
		CLZNumber	shifted;

		lzi_exact_sum( pX, &one, false, &xp1 );
		lzi_exact_sum( pX, &one, true, &xm1 );
		lzi_complex_product( &xp1, &xp1, &square );
		lz_add( &square, &yy, &p );
		lzi_complex_product( &xm1, &xm1, &square );
		lz_add( &square, &yy, &q );
		lzi_complex_product( &n2, &n2, &n4 );

		long	l_sx = pX->IsNegative() ? -1 : 1;
		long	l_sy = pY->IsNegative() ? -1 : 1;

		// (F, sigma, plus) for the real part and the imaginary...
		struct SLZPartForm { long m_lFunction; long m_lSign; bool m_bPlus; };

		SLZPartForm	forms[ 2 ] = { { 0, 1, false }, { 2, 1, true } };
		CLZNumber*	p_s = &one;

		switch ( lFunction )
		{
		case kInverseArcSine:			forms[ 1 ].m_lSign = l_sy; break;
		case kInverseArcCosine:			forms[ 0 ].m_lFunction = 1; forms[ 1 ].m_lSign = -l_sy; break;
		case kInverseArcSecant:			forms[ 0 ].m_lFunction = 1; forms[ 1 ].m_lSign = l_sy; p_s = &n2; break;
		case kInverseArcCosecant:		forms[ 1 ].m_lSign = -l_sy; p_s = &n2; break;
		case kInverseArcHyperCosine:	forms[ 0 ] = { 2, 1, true }; forms[ 1 ] = { 1, l_sy, false }; break;
		case kInverseArcHyperSecant:	forms[ 0 ] = { 2, 1, true }; forms[ 1 ] = { 1, -l_sy, false }; p_s = &n2; break;

		case kInverseArcHyperSine:
		case kInverseArcHyperCosecant:
		{
			// at iz, (1 -+ y)^2 + x^2; at i/z, (y +- |z|^2)^2 + x^2 over |z|^4...
			bool		b_cosecant = (lFunction == kInverseArcHyperCosecant);
			CLZNumber*	p_shift = b_cosecant ? &n2 : &one;

			lzi_exact_sum( b_cosecant ? pY : &one, b_cosecant ? p_shift : pY, !b_cosecant ? true : false, &shifted );
			lzi_complex_product( &shifted, &shifted, &square );
			lz_add( &square, &xx, &p );
			lzi_exact_sum( pY, p_shift, b_cosecant ? true : false, &shifted );
			lzi_complex_product( &shifted, &shifted, &square );
			lz_add( &square, &xx, &q );
			forms[ 0 ] = { 2, l_sx, true };
			forms[ 1 ] = { 0, -1, false };
			p_s = b_cosecant ? &n4 : &one;
			break;
		}

		default:
			break;
		}

		for ( long l_part = 0; l_part < 2; l_part++ )
		{
			spec[ l_part ].m_lFunction = forms[ l_part ].m_lFunction;
			spec[ l_part ].m_lSign = forms[ l_part ].m_lSign;
			spec[ l_part ].m_bPlus = forms[ l_part ].m_bPlus;
			spec[ l_part ].m_p.Assign( &p );
			spec[ l_part ].m_q.Assign( &q );
			spec[ l_part ].m_s.Assign( p_s );
		}
	}

	CLZNumber	parts[ 2 ];

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( pX, &scratch );

		work.Exact();

		// Part lPart at lWidth places: both roots, then the form's quantity...
		auto	fnPart = [ & ]( long lPart, CLZNumber* pOut, long lWidth ) -> long
		{
			long	l_was = glWorkingPrecision;
			long	l_result = kNoError;
			long	l_wide = lWidth + l_extra;

			glWorkingPrecision = l_wide;

			CLZNumber	s1r;
			CLZNumber	s1i;
			CLZNumber	s2r;
			CLZNumber	s2i;
			CLZNumber	value;
			CLZNumber	left;
			CLZNumber	right;
			CLZNumber	difference;
			CLZNumber	two;

			two.SetBase( l_base, false );
			two.Assign( 2L );
			lzi_complex_root_at( &u1r, &u1i, l_wide, &s1r, &s1i );
			lzi_complex_root_at( &u2r, &u2i, l_wide, &s2r, &s2i );

			// The turned forms answer (eta, -xi)...
			long	l_which = b_turned ? (1 - lPart) : lPart;

			// arcsinh of a quantity, over D where it is |z|^2...
			auto	ArcSinhOf = [ & ]( CLZNumber* pQuantity ) -> long
			{
				if ( !b_divided )
					return lz_archypersine( pQuantity, &value );

				CLZNumber	quotient;

				lz_divide( pQuantity, &n2, &quotient );

				return lz_archypersine( &quotient, &value );
			};

			switch ( e_form )
			{
			case kFormArcSine:
				if ( l_which == 0 )
				{
					lz_multiply( &s1r, &s2r, &left );
					lz_multiply( &s1i, &s2i, &right );
					lz_subtract( &left, &right, &difference );
					l_result = lz_arctangent2( &N, &difference, &value );
				}
				else
				{
					lz_multiply( &s1r, &s2i, &left );
					lz_multiply( &s1i, &s2r, &right );
					lz_subtract( &left, &right, &difference );
					l_result = ArcSinhOf( &difference );
				}
				break;

			case kFormArcCosine:
				if ( l_which == 0 )
				{
					l_result = lz_arctangent2( &s1r, &s2r, &difference );
					lz_multiply( &difference, &two, &value );
				}
				else
				{
					lz_multiply( &s2r, &s1i, &left );
					lz_multiply( &s2i, &s1r, &right );
					lz_subtract( &left, &right, &difference );
					l_result = ArcSinhOf( &difference );
				}
				break;

			default:	// kFormArcHyperCosine
				if ( l_which == 0 )
				{
					lz_multiply( &s1r, &s2r, &left );
					lz_multiply( &s1i, &s2i, &right );
					lz_add( &left, &right, &difference );
					l_result = ArcSinhOf( &difference );
				}
				else
				{
					l_result = lz_arctangent2( &s1i, &s2r, &difference );
					lz_multiply( &difference, &two, &value );
				}
				break;
			}

			if ( b_turned && (lPart == 1) )
				lz_negate( &value );

			// ...in the host's unit, for a circular angle...
			if ( l_turn != 0 )
			{
				CLZNumber	twoPi( kTwoPi );
				CLZNumber	scaled;

				CLZNumber	angle;

				lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, glWorkingPrecision + 2 );
				lzi_times_whole( &value, l_turn, &scaled );
				lz_divide( &scaled, &twoPi, &angle );
				value.Assign( &angle );
			}

			// ...written at the width, a division answering at its target's
			//	precision...
			glWorkingPrecision = lWidth;

			CLZNumber	out;

			out.Assign( &value );
			out.SetPrecision( lWidth );
			out.Round();
			pOut->Assign( &out );
			glWorkingPrecision = l_was;

			return l_result;
		};

		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base );

		for ( long l_part = 0; l_part < 2; l_part++ )
		{
			if ( (l_part == 0) && b_exact_real )
				continue;

			long	l_this = l_part;

			lResult = lzi_first_failure( lResult, fnPart( l_this, &parts[ l_part ], glWorkingPrecision + l_guard ) );
			lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth ) { return fnPart( l_this, pValue, lWidth ); },
									  [ & ]( CLZNumber* pResult, long lPlaces, long lMostPlaces ) -> bool
									  {
										  return (l_turn == 0) && lzi_inverse_decide( pResult, lPlaces, lMostPlaces, spec[ l_this ] );
									  },
									  l_base, l_guard, &parts[ l_part ] );
		}
	}

	if ( b_exact_real )
		lResult = lzi_first_failure( lResult, lzi_answer_as_function( &exact_real, pX, &pTarget->m_real, false, true ) );
	else
		lResult = lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 0 ], pX, &pTarget->m_real, false ) );

	return lzi_first_failure( lResult, lzi_answer_as_function( &parts[ 1 ], pX, &pTarget->m_imag, false ) );
}

/*
	**The inverse functions that take square roots** -- stage 6a, 2026-10-23,
	on Rich's calls: arcsin, arccos, arcsec and arccsc and their hyperbolic
	twins first, the four that take logarithms as 6b. The principal values,
	each part correctly rounded, half-even, in places; on a cut, the
	operand's +0 side, as std::complex answers and Elzed does today; on the
	real axis inside its domain each is its real function's answer to the
	character, an exact nought beside it. The circular four write the
	whole complex angle in the host's unit, the hyperbolic four in radians.
	Nought is a pole of the four reciprocals, the point at infinity,
	kUndefinedResult, as lz_archypersecant( 0 ) answers. A not-a-number, or
	the point at infinity, answers a not-a-number.
*/
static long lzi_complex_inverse( long lFunction, CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( lzi_complex_has( leftOperand, kfQNaN ) || lzi_complex_has( leftOperand, kfInfinity ) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, lzi_complex_has( leftOperand, kfDiv0 ) );

		return kNoError;
	}

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	x;
	CLZNumber	y;

	x.Assign( &leftOperand->m_real );
	y.Assign( &leftOperand->m_imag );
	x.CheckForZero();
	y.CheckForZero();

	// The host's unit, read as an answer guard reads it, for the circular
	//	four -- see lzi_angle_unit_here...
	AngleUnit	e_unit = AU_radians;

	if ( lFunction < kInverseHyperbolic )
	{
		CLZNumber		scratch;
		CLZAnswerBase	probe( &x, &scratch );

		probe.Exact();
		e_unit = lzi_angle_unit_here();
	}

	if ( x.IsZero() || y.IsZero() )
		return lzi_complex_inverse_on_axis( lFunction, &x, &y, e_unit, p_target );

	return lzi_complex_inverse_general( lFunction, &x, &y, e_unit, p_target );
}

long lzc_arcsine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcSine, leftOperand, targetNumber );
}

long lzc_arccosine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcCosine, leftOperand, targetNumber );
}

long lzc_arcsecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcSecant, leftOperand, targetNumber );
}

long lzc_arccosecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcCosecant, leftOperand, targetNumber );
}

long lzc_archypersine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcHyperSine, leftOperand, targetNumber );
}

long lzc_archypercosine( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcHyperCosine, leftOperand, targetNumber );
}

long lzc_archypersecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcHyperSecant, leftOperand, targetNumber );
}

long lzc_archypercosecant( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_inverse( kInverseArcHyperCosecant, leftOperand, targetNumber );
}

/////////////////////////////////////////////////////////////////
//	Complex numbers, stage 6b -- the inverses that take logarithms,
//	2026-10-23

// The four by index: the circular two, then their hyperbolic twins...
const long	kLogInverseArcTangent = 0;
const long	kLogInverseArcCotangent = 1;
const long	kLogInverseArcHyperTangent = 2;
const long	kLogInverseArcHyperCotangent = 3;
const long	kLogInverseHyperbolic = 2;

/*
	**A part of the four, as the decider below compares it**: off the
	axes every part is 1/4 ln(P/Q) or 1/2 atan2(Y, X), with P, Q, Y and X
	exact -- see lzi_complex_log_inverse_general.
*/
struct SLZLogInversePart
{
	bool		m_bLogarithm;		// 1/4 ln(P/Q), or 1/2 atan2(Y, X)
	CLZNumber	m_a;				// P, or Y
	CLZNumber	m_b;				// Q, or X
};

/*
	**A part of the four in doubt, decided exactly** -- in radians, where a
	part can lie a hair from a rational past any retake: arctan z is z and
	a little near nought, arccot z is 1/z and a little far out, and along
	an axis a small part is its leading term, whose derivative is rational
	at every rational point -- arctanh(x + i eps) has the imaginary part
	eps/(1 - x^2) and a little.

	With h the half the part lies near, 1/4 ln(P/Q) is past h exactly when
	P/Q is past e^(4h), compared as P against Q e^(4h), or P e^(-4h)
	against Q; and 1/2 atan2(Y, X) is past h exactly when its angle theta
	is past 2h, which is the sign of sin(theta - 2h), so of Y cos 2h - X
	sin 2h -- the two angles so near that their difference is far below
	pi. The exponential, sine and cosine of the rational come from their
	series as rational intervals, and every comparison is exact; more terms
	each round until one decides. A nonzero rational has a transcendental
	exponential, sine and cosine, so neither comparison is ever an
	equality, and the budget lMostPlaces is the limit.
*/
static bool lzi_log_inverse_decide( CLZNumber* pResult, long lPlaces, long lMostPlaces, const SLZLogInversePart& part )
{
	if ( lPlaces < 0 )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pResult->m_lRadix > lPlaces) ? pResult->m_lRadix : lPlaces) + 4;

	// The half: |part| rounded down at the places, y, and h = y + u/2,
	//	Qn / Qd = sign (2 Y + 1) / (2 b^P), with Y = y b^P...
	CLZNumber	y;
	CLZNumber	Y;
	CLZNumber	Qn;
	CLZNumber	Qd;
	CLZNumber	one;
	CLZNumber	two;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );
	Y.Assign( &y );
	lzi_scale_by_base_power( &Y, -lPlaces );
	one.SetBase( l_base, false );
	one.Assign( 1L );
	two.SetBase( l_base, false );
	two.Assign( 2L );
	lzi_complex_product( &Y, &two, &Qn );
	lz_add( &Qn, &one, &Qn );
	Qd.Assign( &two );
	lzi_scale_by_base_power( &Qd, -lPlaces );

	long	l_sign_h = pResult->IsNegative() ? -1 : 1;

	if ( l_sign_h < 0 )
		lz_negate( &Qn );

	// q = 4h for the logarithm, 2h for the angle...
	CLZNumber	q;
	CLZNumber	q_abs;
	CLZNumber	times;

	times.SetBase( l_base, false );
	times.Assign( part.m_bLogarithm ? 4L : 2L );
	lzi_complex_product( &Qn, &times, &q );
	q_abs.Assign( &q );
	q_abs.m_lSign = kPositive;

	double	d_qn = 0.0;
	double	d_qd = 1.0;

	lzi_value_or_infinity( &q, &d_qn );
	lzi_value_or_infinity( &Qd, &d_qd );

	double	d_q = d_qn / d_qd;
	long	l_versus = 0;				// the sign of the part - h, once known

	auto	Above = []( CLZNumber* pLeft, CLZNumber* pRight ) -> bool { return lz_compare( pLeft, pRight ) == kGreaterThan; };

	for ( long l_terms = 2; l_versus == 0; l_terms *= 2 )
	{
		long	l_power = part.m_bLogarithm ? l_terms : 2 * l_terms + 1;	// the largest power summed
		double	d_digits = (l_power + 1.0) * (std::log( d_qd > 1.0 ? d_qd : 2.0 ) / std::log( (double)l_base ) + 1.0)
						   + std::lgamma( l_power + 2.0 ) / std::log( (double)l_base );

		if ( d_digits > (double)lMostPlaces )
		{
			glWorkingPrecision = l_saved;

			return false;
		}

		// The remainder's bound wants exp's terms halving, N at least
		//	2|q|, and sin's and cos's shrinking, N above |q|...
		if ( part.m_bLogarithm ? ((double)l_terms < 2.0 * std::fabs( d_q ) + 1.0) : ((2.0 * l_terms) <= std::fabs( d_q )) )
			continue;

		if ( part.m_bLogarithm )
		{
			// P against Q e^(4h): e^|q| in [lower, upper] / den...
			CLZNumber	lower;
			CLZNumber	upper;
			CLZNumber	den;
			CLZNumber	left_low;
			CLZNumber	left_high;
			CLZNumber	right_low;
			CLZNumber	right_high;

			lzi_series_interval( 3, &q_abs, &Qd, l_terms, &lower, &upper, &den );

			if ( l_sign_h > 0 )
			{
				// P den against Q lower .. Q upper...
				lzi_complex_product( const_cast<CLZNumber*>( &part.m_a ), &den, &left_low );
				lzi_complex_product( const_cast<CLZNumber*>( &part.m_b ), &lower, &right_low );
				lzi_complex_product( const_cast<CLZNumber*>( &part.m_b ), &upper, &right_high );

				if ( Above( &left_low, &right_high ) )
					l_versus = 1;
				else if ( Above( &right_low, &left_low ) )
					l_versus = -1;
			}
			else
			{
				// P lower .. P upper against Q den...
				lzi_complex_product( const_cast<CLZNumber*>( &part.m_a ), &lower, &left_low );
				lzi_complex_product( const_cast<CLZNumber*>( &part.m_a ), &upper, &left_high );
				lzi_complex_product( const_cast<CLZNumber*>( &part.m_b ), &den, &right_low );

				if ( Above( &left_low, &right_low ) )
					l_versus = 1;
				else if ( Above( &right_low, &left_high ) )
					l_versus = -1;
			}
		}
		else
		{
			// Y cos 2h - X sin 2h over the box of the two intervals, each
			//	end over its own denominator...
			CLZNumber	s_low;
			CLZNumber	s_high;
			CLZNumber	s_den;
			CLZNumber	c_low;
			CLZNumber	c_high;
			CLZNumber	c_den;

			lzi_series_interval( 0, &q, &Qd, l_terms, &s_low, &s_high, &s_den );
			lzi_series_interval( 1, &q, &Qd, l_terms, &c_low, &c_high, &c_den );

			CLZNumber*	p_Y = const_cast<CLZNumber*>( &part.m_a );
			CLZNumber*	p_X = const_cast<CLZNumber*>( &part.m_b );
			bool				b_y = !p_Y->IsNegative();
			bool				b_x = !p_X->IsNegative();

			// The least: Y at the smaller cos, less X at the larger sin;
			//	the greatest the other way...
			auto	Edge = [ & ]( bool bLeast, CLZNumber* pOut )
			{
				CLZNumber	yc;
				CLZNumber	xs;
				CLZNumber	ycd;
				CLZNumber	xsd;

				lzi_complex_product( p_Y, (b_y == bLeast) ? &c_low : &c_high, &yc );
				lzi_complex_product( p_X, (b_x == bLeast) ? &s_high : &s_low, &xs );
				lzi_complex_product( &yc, &s_den, &ycd );
				lzi_complex_product( &xs, &c_den, &xsd );
				lz_subtract( &ycd, &xsd, pOut );
				pOut->CheckForZero();
			};

			CLZNumber	least;
			CLZNumber	greatest;

			Edge( true, &least );
			Edge( false, &greatest );

			if ( least.IsPositive() )
				l_versus = 1;
			else if ( greatest.IsNegative() )
				l_versus = -1;
		}
	}

	// The magnitude goes up where the part lies beyond h...
	bool		b_up = (l_versus * l_sign_h > 0);
	CLZNumber	unit;

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	if ( b_up )
		lz_add( &y, &unit );

	if ( l_sign_h < 0 )
		lz_negate( &y );

	y.CheckForZero();
	pResult->Assign( &y );
	glWorkingPrecision = l_saved;

	return true;
}

/*
	**The angle part in a unit, a hair from an exact one** -- the half of a
	direction that is a multiple of 45 degrees, k/16 of a turn, is 22.5,
	67.5 or their kin in degrees, 62.5 or 187.5 millicles: a half at no
	places, and in an odd base a half at every width, since half a unit
	never ends there -- 62.5 is 42.777... in base fifteen. A point a hair
	from such a direction then lies a hair from a half past any retake,
	and the census found it rounded the wrong way: arctan at 15^-24 from a
	point where 2|x| = 1 - |z|^2, in millicles at one place in base
	fifteen, answered -42.7 for -42.8 -- the exact value there -62.5,
	-42.777... in base fifteen.

	So where the half in doubt is k/16 of a turn exactly, the part lies on
	the side of it that theta, the angle of (Y, X), lies of k 45 degrees --
	the sign of Y cos 45k - X sin 45k, the cosine and sine nought or +-1,
	or both +-1 times sqrt(2)/2, which does not move the sign: exact. A
	half that is not such a value lies at least 1/32 of a unit from one,
	which the retake reaches, so it is left to it.
*/
static bool lzi_log_inverse_unit_decide( CLZNumber* pResult, long lPlaces, long lTurn, const SLZLogInversePart& angle )
{
	if ( lPlaces < 0 )
		return false;

	long	l_base = pResult->m_lBase;
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = ((pResult->m_lRadix > lPlaces) ? pResult->m_lRadix : lPlaces) + 4;

	// The half: |part| rounded down at the places, y, and h = y + u/2,
	//	Qn / Qd = sign (2 Y + 1) / (2 b^P)...
	CLZNumber	y;
	CLZNumber	Y;
	CLZNumber	Qn;
	CLZNumber	Qd;
	CLZNumber	one;
	CLZNumber	two;

	y.Assign( pResult );
	y.m_lSign = kPositive;
	y.CheckForZero();
	y.Round( RM_round_down, lPlaces );
	Y.Assign( &y );
	lzi_scale_by_base_power( &Y, -lPlaces );
	one.SetBase( l_base, false );
	one.Assign( 1L );
	two.SetBase( l_base, false );
	two.Assign( 2L );
	lzi_complex_product( &Y, &two, &Qn );
	lz_add( &Qn, &one, &Qn );
	Qd.Assign( &two );
	lzi_scale_by_base_power( &Qd, -lPlaces );

	long	l_sign_h = pResult->IsNegative() ? -1 : 1;

	if ( l_sign_h < 0 )
		lz_negate( &Qn );

	// The nearest k/16 of a turn, and whether h is it: 16 Qn = k turn Qd...
	double	d_part = 0.0;

	lzi_value_or_infinity( pResult, &d_part );

	long		l_k = static_cast<long>( std::floor( d_part * 16.0 / static_cast<double>( lTurn ) + 0.5 ) );
	CLZNumber	left;
	CLZNumber	right;
	CLZNumber	factor;

	factor.SetBase( l_base, false );
	factor.Assign( 16L );
	lzi_complex_product( &Qn, &factor, &left );
	factor.Assign( l_k * lTurn );
	lzi_complex_product( &Qd, &factor, &right );

	if ( lz_compare( &left, &right ) != kEqualTo )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// ...and theta's side of k 45 degrees, from the sign of Y cos - X sin,
	//	the eight directions' cosines and sines as signs...
	static const int	sl_cos[ 8 ] = { 1, 1, 0, -1, -1, -1, 0, 1 };
	static const int	sl_sin[ 8 ] = { 0, 1, 1, 1, 0, -1, -1, -1 };

	long		l_d = ((l_k % 8) + 8) % 8;
	CLZNumber	yc;
	CLZNumber	xs;
	CLZNumber	E;

	factor.Assign( static_cast<long>( sl_cos[ l_d ] ) );
	lzi_complex_product( const_cast<CLZNumber*>( &angle.m_a ), &factor, &yc );
	factor.Assign( static_cast<long>( sl_sin[ l_d ] ) );
	lzi_complex_product( const_cast<CLZNumber*>( &angle.m_b ), &factor, &xs );
	lz_subtract( &yc, &xs, &E );
	E.CheckForZero();

	if ( E.IsZero() )
	{
		glWorkingPrecision = l_saved;

		return false;
	}

	// The magnitude goes up where the part lies beyond h...
	bool		b_up = ((E.IsPositive() ? 1 : -1) * l_sign_h > 0);
	CLZNumber	unit;

	unit.SetBase( l_base, false );
	unit.Assign( 1L );
	lzi_scale_by_base_power( &unit, lPlaces );

	if ( b_up )
		lz_add( &y, &unit );

	if ( l_sign_h < 0 )
		lz_negate( &y );

	y.CheckForZero();
	pResult->Assign( &y );
	glWorkingPrecision = l_saved;

	return true;
}

/*
	**On an axis each part is a real function's** -- a table checked
	against cmath, and against what Elzed answers today on every cut, from
	the operand's +0 side: arctanh 2 = 0.549 + (pi/2)i, arccoth 0.5 = 0.549
	- (pi/2)i, arctan 2i = pi/2 + 0.549i, arctan -2i = pi/2 - 0.549i,
	arccot 0.5i = pi/2 - 0.549i. arccoth 0 is -(pi/2)i, the cut's value
	there from either side, where Elzed answers a not-a-number for
	forming 1/0. arctanh and arccoth of +-1, and arctan and arccot of
	+-i, are poles: the point at infinity, kUndefinedResult, as
	lz_archypertangent( 1 ) answers. The circular two's angles are in the
	host's unit, a quarter turn exact; their other parts are a real
	function's radians written in the unit -- see lzi_inverse_part. The
	hyperbolic two are in radians.
*/
static long lzi_complex_log_inverse_on_axis( long lFunction, CLZNumber* pX, CLZNumber* pY, AngleUnit eUnit, CLZComplex* pTarget )
{
	long	l_base = pX->m_lBase;
	long	lResult = kNoError;
	bool	b_real = pY->IsZero();

	CLZNumber	v;				// the axis coordinate
	CLZNumber	magnitude;
	CLZNumber	one;
	CLZNumber	nought;

	v.Assign( b_real ? pX : pY );
	magnitude.Assign( &v );
	magnitude.m_lSign = v.IsZero() ? kZero : kPositive;
	one.SetBase( l_base, false );
	one.Assign( 1L );
	nought.SetBase( l_base, false );
	nought.Assign( 0L );

	long	l_against_one = lz_compare( &magnitude, &one );
	bool	b_inside = (l_against_one == kLessThan);		// |v| < 1
	bool	b_beyond = (l_against_one == kGreaterThan);		// |v| > 1

	CLZNumber*	p_real = &pTarget->m_real;
	CLZNumber*	p_imag = &pTarget->m_imag;

	// A part that is the host's call -- in the host's unit for a circular
	//	angle -- or, inside a hyperbolic answer, the same call in radians...
	auto	Host = [ & ]( LZUnaryFunction fn, CLZNumber* pArg, CLZNumber* pOut, bool bNegate ) -> long
	{
		long	l_result = fn( pArg, pOut );

		if ( bNegate )
			lz_negate( pOut );

		return l_result;
	};

	auto	HostRadians = [ & ]( LZUnaryFunction fn, CLZNumber* pArg, CLZNumber* pOut, bool bNegate ) -> long
	{
		CLZRadiansHere	radians;

		return Host( fn, pArg, pOut, bNegate );
	};

	auto	Nought = [ & ]( CLZNumber* pOut ) -> long
	{
		return lzi_answer_as_function( &nought, pX, pOut, false, true, true );
	};

	// A real function's refusal is the point at infinity with its code...
	auto	Refused = [ & ]( long lCode ) -> long
	{
		if ( p_real->FlagIsOn( kfInfinity ) || p_real->FlagIsOn( kfQNaN ) )
			lzi_complex_flagged( pTarget, l_base, p_real->FlagIsOn( kfInfinity ) ? kfInfinity : kfQNaN, p_real->FlagIsOn( kfDiv0 ) );

		return lCode;
	};

	auto	Pole = [ & ]() -> long
	{
		lzi_complex_flagged( pTarget, l_base, kfInfinity, false );

		return kUndefinedResult;
	};

	switch ( lFunction )
	{
	case kLogInverseArcTangent:
		if ( b_real )
			return Refused( lzi_first_failure( Host( lz_arctangent, &v, p_real, false ), Nought( p_imag ) ) );

		if ( b_inside )
			return lzi_first_failure( Nought( p_real ), lzi_inverse_part( lz_archypertangent, &v, false, eUnit, p_imag ) );

		if ( !b_beyond )
			return Pole();

		lResult = Host( lz_arcsine, &one, p_real, false );

		return lzi_first_failure( lResult, lzi_inverse_part( lz_archypercotangent, &v, false, eUnit, p_imag ) );

	case kLogInverseArcCotangent:
		if ( b_real )
			return Refused( lzi_first_failure( Host( lz_arccotangent, &v, p_real, false ), Nought( p_imag ) ) );

		if ( b_beyond )
			return lzi_first_failure( Nought( p_real ), lzi_inverse_part( lz_archypercotangent, &v, true, eUnit, p_imag ) );

		if ( !b_inside )
			return Pole();

		lResult = Host( lz_arcsine, &one, p_real, false );

		return lzi_first_failure( lResult, lzi_inverse_part( lz_archypertangent, &v, true, eUnit, p_imag ) );

	case kLogInverseArcHyperTangent:
		if ( !b_real )
			return lzi_first_failure( Nought( p_real ), HostRadians( lz_arctangent, &v, p_imag, false ) );

		if ( !b_beyond )
			return Refused( lzi_first_failure( Host( lz_archypertangent, &v, p_real, false ), Nought( p_imag ) ) );

		lResult = Host( lz_archypercotangent, &v, p_real, false );

		return lzi_first_failure( lResult, HostRadians( lz_arcsine, &one, p_imag, false ) );

	default:	// kLogInverseArcHyperCotangent
		if ( !b_real )
			return lzi_first_failure( Nought( p_real ), HostRadians( lz_arccotangent, &v, p_imag, true ) );

		if ( !b_inside )
			return Refused( lzi_first_failure( Host( lz_archypercotangent, &v, p_real, false ), Nought( p_imag ) ) );

		lResult = Host( lz_archypertangent, &v, p_real, false );

		return lzi_first_failure( lResult, HostRadians( lz_arcsine, &one, p_imag, true ) );
	}
}

/*
	**Off both axes, a logarithm and an angle of exact values**, with
	x and y exact and |z|^2 = x^2 + y^2:

	| | real part | imaginary part |
	|---|---|---|
	| arctanh | 1/4 ln(|z + 1|^2 / |z - 1|^2) | 1/2 atan2(2y, 1 - |z|^2) |
	| arccoth | the same | 1/2 atan2(-2y, |z|^2 - 1) |
	| arctan | 1/2 atan2(2x, 1 - |z|^2) | 1/4 ln(|z + i|^2 / |z - i|^2) |
	| arccot | 1/2 atan2(2x, |z|^2 - 1) | 1/4 ln(|z - i|^2 / |z + i|^2) |

	arccoth z is arctanh(1/z) and arccot z arctan(1/z), each worked
	without forming 1/z, which need not end in the base. Checked against
	cmath at 20,000 points each. **The logarithm is ln P - ln Q**, each of
	an exact value and correctly rounded in places, so nothing cancels
	that places can feel; the angle is lz_arctangent2 of exact arguments.
	Each part taken again wider in doubt, with the tie deciders above --
	in radians lzi_log_inverse_decide, in a unit for the angle part
	lzi_log_inverse_unit_decide.

	**In a unit the circular two write both parts in it**, and the angle
	part is exact where the angle's direction is a multiple of 45
	degrees, the only rational tangents of a rational part of a turn being
	nought, one and infinity: at |z| = 1 it is a quarter turn's half, 45
	degrees, and where 2|x| = |1 - |z|^2| it is 22.5 or 67.5. The
	logarithm part is never exact: ln of a rational other than one,
	over pi, is transcendental.
*/
static long lzi_complex_log_inverse_general( long lFunction, CLZNumber* pX, CLZNumber* pY, AngleUnit eUnit, CLZComplex* pTarget )
{
	long	l_base = pX->m_lBase;
	long	lResult = kNoError;
	bool	b_hyperbolic = (lFunction >= kLogInverseHyperbolic);
	long	l_turn = b_hyperbolic ? 0 : lzi_turn_of( eUnit );

	CLZNumber	one;
	CLZNumber	two;
	CLZNumber	xx;
	CLZNumber	yy;
	CLZNumber	n2;
	CLZNumber	one_less;			// 1 - |z|^2
	CLZNumber	shifted;
	CLZNumber	square;

	one.SetBase( l_base, false );
	one.Assign( 1L );
	two.SetBase( l_base, false );
	two.Assign( 2L );
	lzi_complex_product( pX, pX, &xx );
	lzi_complex_product( pY, pY, &yy );
	lz_add( &xx, &yy, &n2 );
	lzi_exact_sum( &one, &n2, true, &one_less );

	// The real part's form and the imaginary's...
	SLZLogInversePart	forms[ 2 ];
	long				l_angle = b_hyperbolic ? 1 : 0;		// which part is the angle
	SLZLogInversePart&	angle = forms[ l_angle ];
	SLZLogInversePart&	logarithm = forms[ 1 - l_angle ];

	angle.m_bLogarithm = false;
	logarithm.m_bLogarithm = true;

	// The angle's Y and X...
	lzi_complex_product( b_hyperbolic ? pY : pX, &two, &angle.m_a );
	angle.m_b.Assign( &one_less );

	if ( lFunction == kLogInverseArcHyperCotangent )
		lz_negate( &angle.m_a );

	if ( (lFunction == kLogInverseArcHyperCotangent) || (lFunction == kLogInverseArcCotangent) )
		lz_negate( &angle.m_b );

	angle.m_a.CheckForZero();
	angle.m_b.CheckForZero();

	// ...and the logarithm's P and Q, |z +- 1|^2 or |z +- i|^2, the other
	//	way about for arccot...
	CLZNumber*	p_along = b_hyperbolic ? pX : pY;
	CLZNumber*	p_square = b_hyperbolic ? &yy : &xx;

	lzi_exact_sum( p_along, &one, false, &shifted );
	lzi_complex_product( &shifted, &shifted, &square );
	lz_add( &square, p_square, &logarithm.m_a );
	lzi_exact_sum( p_along, &one, true, &shifted );
	lzi_complex_product( &shifted, &shifted, &square );
	lz_add( &square, p_square, &logarithm.m_b );

	if ( lFunction == kLogInverseArcCotangent )
	{
		CLZNumber	swap;

		swap.Assign( &logarithm.m_a );
		logarithm.m_a.Assign( &logarithm.m_b );
		logarithm.m_b.Assign( &swap );
	}

	// In a unit, the angle exact where its direction is a multiple of 45
	//	degrees: half of a quarter, an eighth or three eighths of a turn...
	bool		b_exact_angle = false;
	CLZNumber	exact_angle;

	if ( l_turn != 0 )
	{
		long	l_y = angle.m_a.IsNegative() ? -1 : 1;
		long	l_a = 0;
		long	l_b = 0;

		if ( angle.m_b.IsZero() )
		{
			l_a = l_y;
			l_b = 8;
		}
		else
		{
			CLZNumber	y_abs;
			CLZNumber	x_abs;

			y_abs.Assign( &angle.m_a );
			y_abs.m_lSign = kPositive;
			x_abs.Assign( &angle.m_b );
			x_abs.m_lSign = kPositive;

			if ( lz_compare( &y_abs, &x_abs ) == kEqualTo )
			{
				l_a = angle.m_b.IsNegative() ? 3 * l_y : l_y;
				l_b = 16;
			}
		}

		if ( l_b != 0 )
		{
			lzi_part_of_turn( l_turn, l_a, l_b, l_base, &exact_angle );
			b_exact_angle = true;
		}
	}

	CLZNumber	parts[ 2 ];

	{
		CLZNumber		scratch;
		CLZAnswerBase	work( pX, &scratch );

		work.Exact();

		// Part lPart at lWidth places, and three more in a unit, whose turn
		//	over 2pi reaches 159...
		auto	fnPart = [ & ]( long lPart, CLZNumber* pOut, long lWidth ) -> long
		{
			long	l_was = glWorkingPrecision;
			long	l_result = kNoError;
			long	l_wide = lWidth + ((l_turn != 0) ? lzi_guard_places( 3, l_base ) : 0);

			glWorkingPrecision = l_wide;

			const SLZLogInversePart&	form = forms[ lPart ];
			CLZNumber					value;
			CLZNumber					work_a;
			CLZNumber					work_b;
			CLZNumber					difference;

			if ( form.m_bLogarithm )
			{
				CLZNumber	four;

				four.SetBase( l_base, false );
				four.Assign( 4L );
				l_result = lz_log_e( const_cast<CLZNumber*>( &form.m_a ), &work_a );
				l_result = lzi_first_failure( l_result, lz_log_e( const_cast<CLZNumber*>( &form.m_b ), &work_b ) );
				lz_subtract( &work_a, &work_b, &difference );
				lz_divide( &difference, &four, &value );
			}
			else
			{
				l_result = lz_arctangent2( const_cast<CLZNumber*>( &form.m_a ), const_cast<CLZNumber*>( &form.m_b ), &difference );
				lz_divide( &difference, &two, &value );
			}

			// ...in the host's unit, for a circular part...
			if ( l_turn != 0 )
			{
				CLZNumber	twoPi( kTwoPi );
				CLZNumber	scaled;
				CLZNumber	turned;

				lzi_pi_family_in_base( &twoPi, kTwoPi, l_base, glWorkingPrecision + 2 );
				lzi_times_whole( &value, l_turn, &scaled );
				lz_divide( &scaled, &twoPi, &turned );
				value.Assign( &turned );
			}

			// ...written at the width, a division answering at its target's
			//	precision...
			glWorkingPrecision = lWidth;

			CLZNumber	out;

			out.Assign( &value );
			out.SetPrecision( lWidth );
			out.Round();
			pOut->Assign( &out );
			glWorkingPrecision = l_was;

			return l_result;
		};

		long	l_guard = lzi_guard_places( kComplexGuardDigits, l_base );

		for ( long l_part = 0; l_part < 2; l_part++ )
		{
			if ( (l_part == l_angle) && b_exact_angle )
				continue;

			long	l_this = l_part;

			lResult = lzi_first_failure( lResult, fnPart( l_this, &parts[ l_part ], glWorkingPrecision + l_guard ) );
			lzi_complex_part_rounded( [ & ]( CLZNumber* pValue, long lWidth ) { return fnPart( l_this, pValue, lWidth ); },
									  [ & ]( CLZNumber* pResult, long lPlaces, long lMostPlaces ) -> bool
									  {
										  if ( l_turn != 0 )
											  return (l_this == l_angle) && lzi_log_inverse_unit_decide( pResult, lPlaces, l_turn, forms[ l_this ] );

										  return lzi_log_inverse_decide( pResult, lPlaces, lMostPlaces, forms[ l_this ] );
									  },
									  l_base, l_guard, &parts[ l_part ] );
		}
	}

	for ( long l_part = 0; l_part < 2; l_part++ )
	{
		CLZNumber*	p_out = (l_part == 0) ? &pTarget->m_real : &pTarget->m_imag;

		if ( (l_part == l_angle) && b_exact_angle )
			lResult = lzi_first_failure( lResult, lzi_answer_as_function( &exact_angle, pX, p_out, false, true ) );
		else
			lResult = lzi_first_failure( lResult, lzi_answer_as_function( &parts[ l_part ], pX, p_out, false ) );
	}

	return lResult;
}

/*
	**The inverse functions that take logarithms** -- stage 6b, 2026-10-23,
	on Rich's calls: arctan, arccot, arctanh and arccoth, arccot being
	arctan(1/z), as the real lz_arccotangent is arctan(1/x). The principal
	values, each part correctly rounded, half-even, in places; on a cut,
	the operand's +0 side, as std::complex answers and Elzed does today; on
	the real axis inside its domain each is its real function's answer to
	the character, an exact nought beside it. The circular two write the
	whole complex angle in the host's unit, the hyperbolic two in radians.
	arctanh and arccoth of +-1, and arctan and arccot of +-i, are poles,
	kUndefinedResult, the point at infinity. A not-a-number, or the point
	at infinity, answers a not-a-number.
*/
static long lzi_complex_log_inverse( long lFunction, CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	long	l_base = leftOperand->m_real.m_lBase;

	if ( lzi_complex_has( leftOperand, kfQNaN ) || lzi_complex_has( leftOperand, kfInfinity ) )
	{
		lzi_complex_flagged( p_target, l_base, kfQNaN, lzi_complex_has( leftOperand, kfDiv0 ) );

		return kNoError;
	}

	// Copies, so an answer written in place cannot reach the work...
	CLZNumber	x;
	CLZNumber	y;

	x.Assign( &leftOperand->m_real );
	y.Assign( &leftOperand->m_imag );
	x.CheckForZero();
	y.CheckForZero();

	// The host's unit, read as an answer guard reads it, for the circular
	//	two -- see lzi_angle_unit_here...
	AngleUnit	e_unit = AU_radians;

	if ( lFunction < kLogInverseHyperbolic )
	{
		CLZNumber		scratch;
		CLZAnswerBase	probe( &x, &scratch );

		probe.Exact();
		e_unit = lzi_angle_unit_here();
	}

	if ( x.IsZero() || y.IsZero() )
		return lzi_complex_log_inverse_on_axis( lFunction, &x, &y, e_unit, p_target );

	return lzi_complex_log_inverse_general( lFunction, &x, &y, e_unit, p_target );
}

long lzc_arctangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_log_inverse( kLogInverseArcTangent, leftOperand, targetNumber );
}

long lzc_arccotangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_log_inverse( kLogInverseArcCotangent, leftOperand, targetNumber );
}

long lzc_archypertangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_log_inverse( kLogInverseArcHyperTangent, leftOperand, targetNumber );
}

long lzc_archypercotangent( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	return lzi_complex_log_inverse( kLogInverseArcHyperCotangent, leftOperand, targetNumber );
}


long lzc_add( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, rightOperand, p_target, kComplexAdd ) )
		return kNoError;

	// Part by part, each an ordinary sum -- exact, and a mixed pair worked
	//	as lz_add works one...
	lResult = lz_add( &leftOperand->m_real, &rightOperand->m_real, &p_target->m_real );

	return lzi_first_failure( lResult, lz_add( &leftOperand->m_imag, &rightOperand->m_imag, &p_target->m_imag ) );
}


long lzc_conjugate( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, nullptr, p_target, kComplexUnary ) )
		return kNoError;

	// In place the real part stands as it is...
	if ( p_target == leftOperand )
		return lz_negate( &leftOperand->m_imag );

	lResult = lz_plus( &leftOperand->m_real, &p_target->m_real );

	return lzi_first_failure( lResult, lz_negate( &leftOperand->m_imag, &p_target->m_imag ) );
}


long lzc_divide( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, rightOperand, p_target, kComplexDivide ) )
		return kNoError;

	long				l_base = leftOperand->m_real.m_lBase;
	CLZComplexInBase	right( rightOperand, l_base );

	right.Real()->CheckForZero();
	right.Imag()->CheckForZero();

	/*
		**Over nought, both parts are lz_divide's answer for one over
		nought** -- the point at infinity, flagged as a division by
		nought -- or for nought over nought, a not-a-number, where the
		dividend is nought too. Asked before anything is written: in
		place the target is the dividend...
	*/
	if ( right.Real()->IsZero() && right.Imag()->IsZero() )
	{
		CLZNumber	dividend;
		CLZNumber	zero;

		dividend.SetBase( l_base, false );
		dividend.Assign( lzi_complex_is_zero( leftOperand ) ? 0L : 1L );
		zero.SetBase( l_base, false );
		zero.Assign( 0L );

		lResult = lz_divide( &dividend, &zero, &p_target->m_real );

		return lzi_first_failure( lResult, lz_divide( &dividend, &zero, &p_target->m_imag ) );
	}

	/*
		(a + bi) / (c + di) is ((ac + bd) + (bc - ad) i) / (c^2 + d^2):
		three exact numbers and one correctly rounded division a part. A
		right operand held over D is (Yr + i Yi) / D, and dividing by it
		multiplies the numerators by D...
	*/
	CLZNumber	ac;
	CLZNumber	bd;
	CLZNumber	bc;
	CLZNumber	ad;
	CLZNumber	cc;
	CLZNumber	dd;
	CLZNumber	real;
	CLZNumber	imag;
	CLZNumber	divisor;

	lzi_complex_product( &leftOperand->m_real, right.Real(), &ac );
	lzi_complex_product( &leftOperand->m_imag, right.Imag(), &bd );
	lzi_complex_product( &leftOperand->m_imag, right.Real(), &bc );
	lzi_complex_product( &leftOperand->m_real, right.Imag(), &ad );
	lzi_complex_product( right.Real(), right.Real(), &cc );
	lzi_complex_product( right.Imag(), right.Imag(), &dd );
	lz_add( &ac, &bd, &real );
	lz_subtract( &bc, &ad, &imag );
	lz_add( &cc, &dd, &divisor );

	if ( right.Scaled() )
	{
		CLZNumber	scaled;

		lzi_complex_product( &real, right.Scale(), &scaled );
		real.Assign( &scaled );
		lzi_complex_product( &imag, right.Scale(), &scaled );
		imag.Assign( &scaled );
	}

	lResult = lz_divide( &real, &divisor, &p_target->m_real );

	return lzi_first_failure( lResult, lz_divide( &imag, &divisor, &p_target->m_imag ) );
}


long lzc_multiply( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, rightOperand, p_target, kComplexMultiply ) )
		return kNoError;

	// (a + bi)(c + di) is (ac - bd) + (ad + bc) i, both exact before either
	//	part of the target is written -- in place it is an operand...
	CLZComplexInBase	right( rightOperand, leftOperand->m_real.m_lBase );

	CLZNumber	ac;
	CLZNumber	bd;
	CLZNumber	ad;
	CLZNumber	bc;
	CLZNumber	real;
	CLZNumber	imag;

	lzi_complex_product( &leftOperand->m_real, right.Real(), &ac );
	lzi_complex_product( &leftOperand->m_imag, right.Imag(), &bd );
	lzi_complex_product( &leftOperand->m_real, right.Imag(), &ad );
	lzi_complex_product( &leftOperand->m_imag, right.Real(), &bc );
	lz_subtract( &ac, &bd, &real );
	lz_add( &ad, &bc, &imag );

	return lzi_complex_finish( &real, &imag, right, p_target );
}


long lzc_negate( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, nullptr, p_target, kComplexUnary ) )
		return kNoError;

	if ( p_target == leftOperand )
		return lzi_first_failure( lz_negate( &leftOperand->m_real ), lz_negate( &leftOperand->m_imag ) );

	lResult = lz_negate( &leftOperand->m_real, &p_target->m_real );

	return lzi_first_failure( lResult, lz_negate( &leftOperand->m_imag, &p_target->m_imag ) );
}


long lzc_square( CLZComplex* leftOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, nullptr, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, nullptr, p_target, kComplexUnary ) )
		return kNoError;

	// (a + bi)^2 is (a^2 - b^2) + 2ab i...
	CLZNumber	aa;
	CLZNumber	bb;
	CLZNumber	ab;
	CLZNumber	real;
	CLZNumber	imag;

	lzi_complex_product( &leftOperand->m_real, &leftOperand->m_real, &aa );
	lzi_complex_product( &leftOperand->m_imag, &leftOperand->m_imag, &bb );
	lzi_complex_product( &leftOperand->m_real, &leftOperand->m_imag, &ab );
	lz_subtract( &aa, &bb, &real );
	lz_add( &ab, &ab, &imag );

	// ...each part answered as lz_square answers, a function's answer...
	lResult = lzi_complex_round_as_square( &real, &p_target->m_real );

	return lzi_first_failure( lResult, lzi_complex_round_as_square( &imag, &p_target->m_imag ) );
}


long lzc_subtract( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber )
{
	CLZComplex*	p_target = (targetNumber != nullptr) ? targetNumber : leftOperand;
	long		lResult = kNoError;

	if ( lzi_complex_refused( leftOperand, rightOperand, p_target, &lResult ) )
		return lResult;

	if ( lzi_complex_propagate( leftOperand, rightOperand, p_target, kComplexSubtract ) )
		return kNoError;

	// Part by part, as lzc_add...
	lResult = lz_subtract( &leftOperand->m_real, &rightOperand->m_real, &p_target->m_real );

	return lzi_first_failure( lResult, lz_subtract( &leftOperand->m_imag, &rightOperand->m_imag, &p_target->m_imag ) );
}
