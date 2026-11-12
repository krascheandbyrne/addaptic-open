/////////////////////////////////////////////////////////////////
//
//	AddapticThreads.cpp
//
//	The threads gate: the engine's functions computed by several
//	threads at once, every answer compared with the one a single
//	thread gives.
//
//	Copyright (c) 2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

/*
	Why this exists.

	Until 2026-10-10 the engine kept its working precision, its retake
	state, its rounding records and its constant cache in process-wide
	variables, so two threads computing at once overwrote each other.
	Eight threads gave wrong answers, crashed or hung on every run. Each
	thread now keeps its own, and this is
	what holds it there.

	Every function in the census runs on ten arguments, in four bases
	at four widths. Eight threads, two to each base and width, work the
	answers out at once, and every answer is compared with the one the
	main thread gives:

	  own numbers      each thread makes its own operands
	  shared numbers   the two threads of a pair read the same operand
	                   objects, made once by the main thread. A number
	                   may be read by several threads at once, so the
	                   engine must never write to one it only reads
	  reciprocal path  division and the reciprocal again, with the
	                   reciprocal's floor lowered so that they take it
	  precision        a new thread starts at kInitialPrecision whatever
	                   the main thread has set, and setting its own
	                   leaves the main thread's alone

	**The threads go first and the main thread after**, and every cache
	is emptied between phases, so each phase starts cold. That is not a
	nicety. Worked out first, the main thread's answers filled the
	cache with every conversion the threads would need and counted the
	places in every constant's text, and a cache or a lazy static shared
	by mistake was then only ever read: ThreadSanitizer saw nothing.

	**What it cannot see.** A race that changes an answer, crashes or
	hangs shows here. A race that leaves every answer right -- a write
	of the value already there -- does not, and nor does one that
	struck on no run. ThreadSanitizer sees those: configure with
	-DADDAPTIC_TSAN=ON and run this.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

#include "LZNumber.h"
#include "LZMath.h"


typedef long (*Unary)( CLZNumber*, CLZNumber* );
typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );
typedef long (*OfList)( CLZNumber**, long, CLZNumber* );

struct SUnary	{ const char* szName; Unary fnCall; };
struct SBinary	{ const char* szName; Binary fnCall; };
struct SOfList	{ const char* szName; OfList fnCall; };

static const SUnary	kUnaries[] =
{
	{ "abs",               lz_abs },
	{ "arccosecant",       lz_arccosecant },
	{ "arccosine",         lz_arccosine },
	{ "arccotangent",      lz_arccotangent },
	{ "archypercosecant",  lz_archypercosecant },
	{ "archypercosine",    lz_archypercosine },
	{ "archypercotangent", lz_archypercotangent },
	{ "archypersecant",    lz_archypersecant },
	{ "archypersine",      lz_archypersine },
	{ "archypertangent",   lz_archypertangent },
	{ "arcsecant",         lz_arcsecant },
	{ "arcsine",           lz_arcsine },
	{ "arctangent",        lz_arctangent },
	{ "cosecant",          lz_cosecant },
	{ "cosine",            lz_cosine },
	{ "cotangent",         lz_cotangent },
	{ "cube",              lz_cube },
	{ "exp_b",             lz_exp_b },
	{ "exp_e",             lz_exp_e },
	{ "factorial",         lz_factorial },
	{ "gamma",             lz_gamma },
	{ "log_gamma",         lz_log_gamma },
	{ "fraction",          lz_fraction },
	{ "hypercosecant",     lz_hypercosecant },
	{ "hypercosine",       lz_hypercosine },
	{ "hypercotangent",    lz_hypercotangent },
	{ "hypersecant",       lz_hypersecant },
	{ "hypersine",         lz_hypersine },
	{ "hypertangent",      lz_hypertangent },
	{ "integer",           lz_integer },
	{ "log_b",             lz_log_b },
	{ "log_e",             lz_log_e },
	{ "negate",            lz_negate },
	{ "normalize",         lz_normalize },
	{ "plus",              lz_plus },
	{ "raise_e",           lz_raise_e },
	{ "reciprocal",        lz_reciprocal },
	{ "root_2",            lz_root_2 },
	{ "root_3",            lz_root_3 },
	{ "round",             lz_round },
	{ "round_integral",    lz_round_integral },
	{ "round_last_digit",  lz_round_last_digit },
	{ "secant",            lz_secant },
	{ "sine",              lz_sine },
	{ "square",            lz_square },
	{ "tangent",           lz_tangent },
};

static const SBinary	kBinaries[] =
{
	{ "add",               static_cast<Binary>( lz_add ) },
	{ "agm",               static_cast<Binary>( lz_agm ) },
	{ "combinations",      static_cast<Binary>( lz_combinations ) },
	{ "divide",            static_cast<Binary>( lz_divide ) },
	{ "divide_integer",    static_cast<Binary>( lz_divide_integer ) },
	{ "exp_n",             static_cast<Binary>( lz_exp_n ) },
	{ "log_n",             static_cast<Binary>( lz_log_n ) },
	{ "modulo",            static_cast<Binary>( lz_modulo ) },
	{ "multiply",          static_cast<Binary>( lz_multiply ) },
	{ "permutations",      static_cast<Binary>( lz_permutations ) },
	{ "power",             static_cast<Binary>( lz_power ) },
	{ "quantize",          static_cast<Binary>( lz_quantize ) },
	{ "raise_n",           static_cast<Binary>( lz_raise_n ) },
	{ "remainder",         static_cast<Binary>( lz_remainder ) },
	{ "remainder_near",    static_cast<Binary>( lz_remainder_near ) },
	{ "root_n",            static_cast<Binary>( lz_root_n ) },
	{ "subtract",          static_cast<Binary>( lz_subtract ) },
};

static const SOfList	kOfLists[] =
{
	{ "amean",             lz_amean },
	{ "gmean",             lz_gmean },
	{ "hmean",             lz_hmean },
};

// The bases and widths, each worked by two threads at once. Base two,
//	seven and sixteen convert their constants, so each thread fills its
//	own cache...
struct SConfig { long lBase; long lPlaces; };

static const SConfig	kConfigs[] = { { 10, 36 }, { 16, 80 }, { 2, 50 }, { 7, 60 } };
const int		kConfigCount = 4;
const int		kThreadCount = 2 * kConfigCount;

// A binary function takes an even-numbered argument on the left and an
//	odd-numbered one on the right, so 3E-25 is never a root's index or a
//	power's exponent: x to the 3E+24 is past the answer limit, and every
//	answer compared here is one worked out...
static const char*	kArguments[] = { "1.2345678901234567", "7.5", "-0.31415", "12.375", "0.0009765625",
									 "2", "0", "-2", "3E-25", "250.5" };
const int		kArgumentCount = 10;

typedef std::vector<CLZNumber*>		Operands;
typedef std::vector<std::string>	Answers;


// The arguments in a configuration's base, made at its width...
static Operands MakeOperands( const SConfig& config )
{
	Operands	operands;

	glWorkingPrecision = config.lPlaces;

	for ( const char* szArgument : kArguments )
	{
		auto	p_operand = new CLZNumber( szArgument );

		if ( config.lBase != 10 )
			p_operand->SetBase( config.lBase );

		operands.push_back( p_operand );
	}

	return operands;
}

static void DeleteOperands( Operands& operands )
{
	for ( CLZNumber* p_operand : operands )
		delete p_operand;

	operands.clear();
}

static std::string Answer( long lError, CLZNumber& result )
{
	char*		sz_result = result.GetStringRep();
	std::string	s_answer = std::to_string( lError ) + " " + ((sz_result != nullptr) ? sz_result : "(none)");

	free( sz_result );

	return s_answer;
}

// Every answer a configuration gives, always in the same order...
static Answers Work( const SConfig& config, const Operands& operands )
{
	Answers		answers;
	CLZNumber*	const* a = operands.data();

	for ( const SUnary& unary : kUnaries )
		for ( int i = 0; i < kArgumentCount; i++ )
		{
			glWorkingPrecision = config.lPlaces;

			CLZNumber	result;
			long		l_error = unary.fnCall( a[ i ], &result );

			answers.push_back( Answer( l_error, result ) );
		}

	for ( const SBinary& binary : kBinaries )
		for ( int i = 0; i < kArgumentCount; i += 2 )
			for ( int j = 1; j < kArgumentCount; j += 2 )
			{
				glWorkingPrecision = config.lPlaces;

				CLZNumber	result;
				long		l_error = binary.fnCall( a[ i ], a[ j ], &result );

				answers.push_back( Answer( l_error, result ) );
			}

	for ( const SOfList& ofList : kOfLists )
		for ( int i = 0; i + 2 < kArgumentCount; i++ )
		{
			glWorkingPrecision = config.lPlaces;

			CLZNumber	result;
			CLZNumber*	values[] = { a[ i ], a[ i + 1 ], a[ i + 2 ] };
			long		l_error = ofList.fnCall( values, 3, &result );

			answers.push_back( Answer( l_error, result ) );
		}

	return answers;
}

// Division and the reciprocal alone, for the reciprocal path...
static Answers WorkDivisions( const SConfig& config, const Operands& operands )
{
	Answers		answers;
	CLZNumber*	const* a = operands.data();

	for ( int i = 0; i < kArgumentCount; i++ )
	{
		glWorkingPrecision = config.lPlaces;

		CLZNumber	quotient;
		CLZNumber	reciprocal;
		long		l_error = lz_divide( a[ (i + 1) % kArgumentCount ], a[ i ], &quotient );

		answers.push_back( Answer( l_error, quotient ) );

		l_error = lz_reciprocal( a[ i ], &reciprocal );

		answers.push_back( Answer( l_error, reciprocal ) );
	}

	return answers;
}

typedef Answers (*WorkFunction)( const SConfig&, const Operands& );

// The main thread's answers, which every other thread must match...
static std::vector<Answers> Reference( WorkFunction fnWork )
{
	std::vector<Answers>	reference;

	for ( const SConfig& config : kConfigs )
	{
		Operands	operands = MakeOperands( config );

		reference.push_back( fnWork( config, operands ) );

		DeleteOperands( operands );
	}

	return reference;
}

static long CountDifferences( const Answers& answers, const Answers& reference )
{
	long	l_differences = (answers.size() == reference.size()) ? 0 : 1;

	for ( size_t i = 0; (i < answers.size()) && (i < reference.size()); i++ )
		if ( answers[ i ] != reference[ i ] )
			l_differences++;

	return l_differences;
}

// Eight threads at once, each answering for its configuration. With
//	shared, a pair reads one set of operands; otherwise each thread makes
//	its own...
static std::vector<Answers> RunThreads( WorkFunction fnWork, std::vector<Operands>* pShared )
{
	std::vector<Answers>		answers( kThreadCount );
	std::vector<std::thread>	threads;

	for ( int t = 0; t < kThreadCount; t++ )
	{
		threads.emplace_back( [ &, t ]()
		{
			int			c = t % kConfigCount;
			Operands	own;

			if ( pShared == nullptr )
				own = MakeOperands( kConfigs[ c ] );

			answers[ t ] = fnWork( kConfigs[ c ], (pShared != nullptr) ? (*pShared)[ c ] : own );

			DeleteOperands( own );
		} );
	}

	for ( std::thread& thread : threads )
		thread.join();

	return answers;
}

// Each thread's answers against the main thread's for its configuration,
//	reported, answering whether all agreed...
static bool Compare( const char* szPhase, const std::vector<Answers>& answers, const std::vector<Answers>& reference )
{
	long	l_differences = 0;
	long	l_answers = 0;

	for ( int t = 0; t < kThreadCount; t++ )
	{
		l_differences += CountDifferences( answers[ t ], reference[ t % kConfigCount ] );
		l_answers += static_cast<long>( answers[ t ].size() );
	}

	printf( "%-16s %ld of %ld answers differ\n", szPhase, l_differences, l_answers );

	return l_differences == 0;
}

// Every cache emptied, so the next phase starts cold...
static void StartCold()
{
	lzu_takedown();
	lzu_setup();
}

int main()
{
	lzu_setup();

	long	l_failures = 0;

	// The threads first, while everything is cold, and then the main
	//	thread's answers to compare them with...
	std::vector<Answers>	answers = RunThreads( Work, nullptr );
	std::vector<Answers>	reference = Reference( Work );

	l_failures += Compare( "own numbers", answers, reference ) ? 0 : 1;

	StartCold();

	std::vector<Operands>	shared;

	for ( const SConfig& config : kConfigs )
		shared.push_back( MakeOperands( config ) );

	answers = RunThreads( Work, &shared );

	l_failures += Compare( "shared numbers", answers, reference ) ? 0 : 1;

	for ( Operands& operands : shared )
		DeleteOperands( operands );

	// The reciprocal's floor as the tests lower it, so that a division
	//	takes it; a setting of the process's, so made before the threads
	//	start...
	StartCold();

	long	l_floor = lzi_set_reciprocal_floor( 2 );

	answers = RunThreads( WorkDivisions, nullptr );
	reference = Reference( WorkDivisions );

	l_failures += Compare( "reciprocal path", answers, reference ) ? 0 : 1;

	lzi_set_reciprocal_floor( l_floor );

	// A new thread starts at kInitialPrecision, and its own setting
	//	stays its own...
	glWorkingPrecision = 50;

	long		l_started = 0;
	std::thread	thread( [ & ]()
	{
		l_started = glWorkingPrecision;
		glWorkingPrecision = 80;
	} );

	thread.join();

	bool	b_precision = (l_started == kInitialPrecision) && (glWorkingPrecision == 50);

	printf( "precision        a new thread started at %ld, the main thread kept %ld: %s\n",
			l_started, glWorkingPrecision, b_precision ? "pass" : "FAIL" );

	l_failures += b_precision ? 0 : 1;

	glWorkingPrecision = kInitialPrecision;

	// ...and so does the angle unit: a new thread reads in radians whatever
	//	another has set -- its sine of 30 is this thread's radian one, to the
	//	character -- while this thread keeps degrees...
	CLZNumber	thirty;
	CLZNumber	radianSine;

	thirty.Assign( "30" );
	lz_sine( &thirty, &radianSine );

	char*	sz_reference = radianSine.GetStringRep();

	geAngleUnit = AU_degrees;

	AngleUnit	e_started = AU_degrees;
	bool		b_radianSine = false;
	std::thread	unitThread( [ & ]()
	{
		e_started = geAngleUnit;

		CLZNumber	x;
		CLZNumber	sine;

		x.Assign( "30" );
		lz_sine( &x, &sine );

		char*	sz_sine = sine.GetStringRep();

		b_radianSine = (sz_sine != nullptr) && (sz_reference != nullptr) && (::strcmp( sz_sine, sz_reference ) == 0);
		::free( sz_sine );

		geAngleUnit = AU_grads;
	} );

	unitThread.join();
	::free( sz_reference );

	bool	b_unit = (e_started == AU_radians) && b_radianSine && (geAngleUnit == AU_degrees);

	printf( "angle unit       a new thread started in radians, the main thread kept degrees: %s\n",
			b_unit ? "pass" : "FAIL" );

	l_failures += b_unit ? 0 : 1;

	geAngleUnit = AU_radians;

	lzu_takedown();

	printf( "%s\n", (l_failures == 0) ? "passed" : "FAILED" );

	return (l_failures == 0) ? 0 : 1;
}
