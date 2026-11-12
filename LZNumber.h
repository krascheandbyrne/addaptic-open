/////////////////////////////////////////////////////////////////
//
//	LZNumber.h
//
//	Part of Addaptic
//
//	Copyright (c) 2006-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#pragma once

#include "LZAddaptic.h"

// Error codes...
const long	kFailed = 1;
const long	kTargetIsConstant = 2;
const long	kDivideByZero = 3;
const long	kUndefinedResult = 4;
const long	kBaseMismatch = 5;
const long	kTooLarge = 6;		// an answer past the limit -- see lzu_set_answer_limit

// Options and indicators...
const long	kPositive = 1;
const long	kZero = 0;
const long	kNegative = -1;

const long	kLessThan = -1;
const long	kEqualTo = 0;
const long	kGreaterThan = 1;

// Neither less, greater, nor equal: what a comparison says when one
//	side is a not-a-number, or an infinity whose sign is not held...
const long	kUnordered = 2;

const long	kLeft = -1;
const long	kEqual = 0;
const long	kRight = 1;

const long	kNotEqual = -1;

const bool	kPreserve = true;
const bool	kNoPreserve = false;

const long	kDigitsLengthIncrement = 16;
const long	kInitialPrecision = 36;
const unsigned long	kInitialDigitsLength = 64;

/*
	The precision a freshly made number starts with, and so the one
	every intermediate result is rounded to. kInitialPrecision unless
	an algorithm has widened it for the span of its own working.

	Widening this is the only way to give an algorithm guard digits.
	Raising a single value's precision does not stick: every operation
	builds its answer in a number made fresh - at this precision - and
	then assigns that over the target, precision and all, before
	rounding. So the working precision is a property of the engine
	while an algorithm runs, not of the values it runs on.

	An algorithm that widens it owns restoring it, on every path out,
	and owes the caller an answer rounded back to the precision it
	would have had. lz_arcsine is the worked example.

	**One per thread.** A thread starts at kInitialPrecision whatever
	another has set, so set it in the thread that computes.
*/
extern thread_local long	glWorkingPrecision;

// Values for the flags member field...
const unsigned long kAllFlagsOff	= 0x00000000;
const unsigned long kAllFlagsOn		= 0xFFFFFFFF;

const unsigned long	kfZero			= 0x00000001;
const unsigned long kfInfinity		= 0x00000002;
const unsigned long kfQNaN			= 0x00000004;
const unsigned long	kfSNaN			= 0x00000008;
const unsigned long kfDiv0			= 0x00000010;

// Rounding modes and digits...
enum RoundingMode
{
	RM_round_down,
	RM_round_half_up,
	RM_round_half_even,
	RM_round_ceiling,
	RM_round_floor,
	RM_round_half_down,
	RM_round_half_odd,
	RM_round_up,
	RM_round_alternate,
	RM_round_random,
	RM_round_truncate,
	RM_round_noop,

	RM_round_default
};

enum RoundingDigit
{
	RD_digit_zero,
	RD_digit_last,
	RD_digit_precision,

	RD_digit_default
};

/*
	The rounding applied to a number that has not asked for any in
	particular.

	Half-even is what the decimal arithmetic standards settle on, and
	for the same reason: of the modes in common use it is the only one
	that does not push a long series of results steadily in one
	direction. Rounding away from zero, which this was, cost a whole
	digit on any value whose discarded part was small - a third came
	back ending in a four.
*/

const RoundingMode	kDefaultRoundingMode = RM_round_half_even;

/*
	How a value is written out.

	RS_positional is what GetStringRep has always done and stays its
	default: every digit in its place, however many noughts that takes.

	The other two write an exponent, to Elzed's convention -- see
	lzi_is_exponent_mark. The mark is E where the base leaves E free
	and X where it does not, and the base raised is the number's own,
	so a base-sixteen value reads 1.0X5 for sixteen to the fifth.

	The same three styles serve both renderers, and the difference
	between the renderers is the digit count -- see GetDisplayRep.
*/
enum RenderStyle
{
	RS_positional,		// 602000000000000000000000
	RS_scientific,		// 6.02E23
	RS_engineering		// 602E21 -- the exponent a multiple of three
};

// The class, its own self...
class CLZNumber
{

public:

	CLZNumber();
	explicit CLZNumber( const char *cpTextNumber );
    explicit CLZNumber( long lLongNumber );
    explicit CLZNumber( double dFloatNumber );
    CLZNumber( const CLZNumber& sourceNumber );

	// Not assignable: the compiler's own assignment would leave two
	// numbers holding one digits array, and both would free it...
	CLZNumber& operator=( const CLZNumber& ) = delete;

	~CLZNumber();


public:

	long	AlignStorage		( unsigned long lDigitsLength, unsigned long lExponent = 0, bool bPreserve = kPreserve );

	/*
		Text is read in the number's own base unless it begins with a
		prefix -- 0b, 0o, 0d, 0t, 0x or 0h, the letter lower case -- which
		reads it in that base and gives the number that base. **A prefix
		wins at the front, in every base**, Rich's call of 2026-10-15: in
		bases twelve to sixteen, where b or d is also a digit, "0b1" is
		one in base two, not 0xB1. Write such digits upper case, as the
		engine writes them -- "0B1" -- or without the leading nought.

		A double is written "%f" and read as decimal: every whole digit,
		six places after the point. An infinity is flagged as one, and a
		NaN as one, refused.
	*/
	long	Assign				( const char *cpTextNumber );
	long	Assign				( long lLongNumber );
	long	Assign				( double dFloatNumber );
	long	Assign				( CLZNumber* p_sourceNumber );

    void	CheckForZero		();
    long	Clear				( unsigned long lDigitsLength = kInitialDigitsLength );
	long	ClearFlag			( unsigned long lFlag );
	// Compare and CompareABS bring a value in another base into this one
	//	first, as lz_compare does; CompareRaw compares the digits as they
	//	stand, and is for two numbers already in one base...
	long	Compare				( CLZNumber* p_compareNumber );
	long	CompareABS			( CLZNumber* p_compareNumber );
	long	CompareRaw			( CLZNumber* p_compareNumber );

    bool	FlagIsOff			( unsigned long lFlag );
    bool	FlagIsOn			( unsigned long lFlag );

	long	GetOrderOfMagnitude	();

	/*
		Two renderers, and the whole difference between them is who
		decides how many digits appear.

		GetStringRep writes the value: every digit it holds, and no
		digit it does not. Nothing is rounded and nothing is padded, so
		what comes back is the value itself and reads back as the value
		itself. That is what the engine's own comparisons want, and what
		anything storing or transmitting a number wants.

		GetDisplayRep writes a display of the value, at a width the
		caller names, in the three modes a calculator has: lDigits is
		places after the radix mark for RS_positional (FIX), and digits
		after the leading one for RS_scientific and RS_engineering
		(SCI and ENG), which is lDigits+1 figures either way. Short
		values are padded out to the width, long ones rounded down to
		it, and the mark carries a nought in front of it where the
		value has no integer part -- a display shows 0.50, not .5.

		**The rounding is done to a copy.** A calculator in FIX 2 reads
		2.35 and goes on multiplying by 2.3456, and keeping the shown
		value and the held value apart is the entire point of a display
		mode. GetDisplayRep therefore changes nothing, and calling it
		on a constant is allowed.

		Neither renderer reads the flags. An infinity or a not-a-number
		has digits like anything else and they will be written; a
		display that means to show something other than digits for
		those asks IsInfinity() or FlagIsOn( kfQNaN ) first. The engine
		has no spelling of its own for either and this is not the place
		to invent one.

		Both return a calloc'd string the caller frees.
	*/

	char*	GetDisplayRep		( RenderStyle rsStyle, long lDigits );
	char*	GetStringRep		();
	char*	GetStringRep		( RenderStyle rsStyle );

	bool	IsConstant			();
	bool	IsFraction			();
	bool	IsInfinity			();
	bool	IsInteger			();
	bool	IsNegative			();
	bool	IsPositive			();
	bool	IsZero				();

	long	PadLeft				( long lPadDigits );
	long	PadRight			( long lPadDigits );

	/*
		**Inside the engine's own work a number's rounding settings are
		not consulted.** While any function is running, Round() and the
		other forms that take the mode or the digit from the number use
		the defaults instead. A temporary copied from a caller's number
		carries that number's mode and digit with it, and until
		2026-09-26 arcsin(1e-8) came out a unit lower in its last place
		when its operand was set to round down. CLZInternalWork in
		LZMath.cpp raises this for every function; the four plain
		operations -- add, subtract, multiply, divide -- do not, and
		round in their target's own mode, to a whole number when its
		digit is the units.

		OwnMode and OwnDigit answer the settings in effect: the
		number's own, or the defaults while this is raised.
	*/

	static void	EnterInternalWork	();
	static void	LeaveInternalWork	();

	RoundingMode	OwnMode		() const;
	RoundingDigit	OwnDigit	() const;

	/*
		**The digits the last rounding discarded**, while recording is on:
		the guard digit and those after it, read as a fraction of the
		place kept, and the places kept. A function that rounds its answer
		from a guard learns from this how near the half the guard lay, and
		decides whether to take the answer again wider -- see
		CLZAnswerBase::Retake in LZMath.cpp. Recording is on only inside a
		function that asked for it, so arithmetic called on its own pays a
		comparison. LastDiscardedTail answers -1 when nothing has been
		rounded since ForgetLastRounding.
	*/

	static void		RecordRoundings		( bool bRecord );
	static void		ForgetLastRounding	();
	static double	LastDiscardedTail	();
	static long		LastRoundedPlaces	();

	// How many discarded digits a recorded tail reads at most, the guard
	//	digit first -- enough for a double in base two. A division cutting
	//	a wide dividend keeps this many places past its answer, so the
	//	tail it records reads the digits it always did...
	static constexpr long	kRecordedTailDigits = 24;

	// A count of places after the point; a negative count rounds to the
	//	left of it, -2 to the hundreds. RM_round_default is the default
	//	mode, half-even by the kept digit...
	long	Round				();
	long	Round				( RoundingMode rm_mode );
	long	Round				( long lPrecision );
	long	Round				( RoundingDigit rd_digit );
	long	Round				( RoundingMode rm_mode, RoundingDigit rd_digit );
	long	Round				( RoundingMode rm_mode, long lPrecision );

	/*
		Round() for a value that goes on past its last digit with the
		middle digit of an odd base repeated forever -- which is
		exactly half a place, and which no finite string of digits in
		an odd base can write. A half in base seven is .3333...

		Long division is what needs it: a quotient whose remainder is
		exactly half the divisor has that tail, and Round() reading the
		digits alone would take it for less than half. In an even base
		half a place is a digit and a run of noughts, so this is plain
		Round() there.
	*/

	long	RoundMiddleRepeating();
	long	RoundMiddleRepeating( RoundingMode rm_mode, long lPrecision );

	long	SetBase				( long lNewBase, bool bConvertValue = true );

	// SetBase for a copy standing in for another value in work done in
	//	lNewBase: it first takes as many of the places it holds as the
	//	working precision needs there, which SetBase alone would drop.
	//	The arithmetic and comparison entry points convert a right
	//	operand this way, and so does CompareABS...
	long	SetBaseForWork		( long lNewBase );

	// The same, for work needing lPlaces places of lNewBase rather than
	//	the working precision -- a modulus sized for its quotient, say...
	long	SetBaseForWork		( long lNewBase, long lPlaces );

	// Places of lToBase worth lPlaces places of lFromBase, rounded up.
	//	A count measured in decimal goes through this before it sizes
	//	anything in another base...
	static long	EquivalentPlaces	( long lPlaces, long lFromBase, long lToBase );

	// For the tests and for measuring, answering the previous setting: the floor
	//	from which SetBase converts by products rather than a digit at a time,
	//	in the digit-at-a-time work it would save, and
	//	kConversionFloorMeasured puts back the measured one; and the count of
	//	conversions made by products, this thread's...
	static constexpr long	kConversionFloorMeasured = -1;

	static long	SetConversionFloor		( long lWork );
	static long	ConversionsByProducts	();

	long	SetConstant			( bool bConstant = true );
	long	SetDigitsArraySize	( long lDigitsArraySize, bool bPreserve = kPreserve );
	long	SetFlag				( unsigned long lFlag );
	long	SetPrecision		( long lPrecision );
	long	SetRoundingDigit	( RoundingDigit rd_digit );
	long	SetRoundingMode		( RoundingMode rm_mode );
	long	ShiftLeft			( long lShiftDigits );
	long	ShiftRight			( long lShiftDigits );
	long	StripFraction		();
	long	StripInteger		();

	long	TrimAll				();
	long	TrimLeading			();
	long	TrimTrailing		();


protected:

	void	BaseConvertValue	( long lOldBase );
	bool	BaseConvertByProducts( long lOldBase, long lFractionBudget, long lIntegerLength );
	void	Initialize			( const char *cpTextNumber = nullptr, unsigned long lDigitsLength = kInitialDigitsLength );
	long	PrecisionFor		( RoundingDigit rd_digit );
	long	RoundAt				( RoundingMode rm_mode, long lPrecision, bool bMiddleRepeats );


public:

	bool			m_bInitialized;
	bool			m_bConstant;

	unsigned long	m_lFlags;
	long			m_lSign;
	long			m_lBase;
	long			m_lMaxDigit;
	long			m_lHalfDigit;
	long			m_lAdjust;
	long			m_lPrecision;
	RoundingMode	m_rmMode;
	RoundingDigit	m_rdDigit;

	long	m_lDigitsArrayLength;
	long	m_lDigitsLength;
	long	m_lRadix;
	long*	m_lpDigits;


protected:

	// What rounding keeps between calls, per thread like the working
	//	precision...
	static thread_local bool		ms_bRoundAlternate;
	static thread_local long		ms_lInternalWork;
	static thread_local long		ms_lRecordRoundings;
	static thread_local double	ms_dLastTail;
	static thread_local long		ms_lLastPlaces;
	static thread_local long		ms_lConversionsByProducts;

};
