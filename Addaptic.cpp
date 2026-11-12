/////////////////////////////////////////////////////////////////
//
//	Addaptic.cpp
//
//	Console bench for the Addaptic engine
//
//	Copyright (c) 2006-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#include <cstdio>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include "LZNumber.h"
#include "LZMath.h"

using namespace std;

///////////////////////////////////////////////////////////
//	Global variables
///////////////////////////////////////////////////////////

long glNewBase = 0;

///////////////////////////////////////////////////////////
//	Function prototypes
///////////////////////////////////////////////////////////

bool NeedsSecondOperand( char c_operand );
void DoOperation( char* firstOperand, char* operation, char* secondOperand, CLZNumber* result );
bool ReadLine( char* szBuffer, long lSize );

const long	kOperandSize = 2048;
const long	kOperationSize = 16;

///////////////////////////////////////////////////////////
//	The Main Line
///////////////////////////////////////////////////////////

int main(int /* argc */, char* /* argv */[])
{
	char*	theAnswer = nullptr;

	char*	firstOperand;
	char*	secondOperand;
	char	operation[ kOperationSize ];
	long	lFirstSize = kOperandSize;

	long	lOperandCount = 0;

	CLZNumber result;

	firstOperand = static_cast<char*>(::calloc( kOperandSize, sizeof(char) ));
	secondOperand = static_cast<char*>(::calloc( kOperandSize, sizeof(char) ));

	bool bKeepTesting = true;

    lzu_setup();

	do
	{
		// Reset...
		lOperandCount = 0;

		// Get the first operand...
		cout << "first operand  >> ";
		cout.flush();

		// Null first operand signals exit, and so does the end of the
		//	input, which a piped script reaches without a blank line...
		if ( !ReadLine( firstOperand, lFirstSize ) || (firstOperand[ 0 ] == '\0') )
		{
			bKeepTesting = false;
			break;
		}

		// Use the last result in a chain calculation, the operand grown to
		//	hold it: a product of two long operands is longer than either...
		if ( (firstOperand[ 0 ] == 'r') && (theAnswer != nullptr) )
		{
			long	l_needed = static_cast<long>( ::strlen( theAnswer ) ) + 1;

			if ( l_needed > lFirstSize )
			{
				char*	p_grown = static_cast<char*>(::realloc( firstOperand, l_needed ));

				if ( p_grown == nullptr )
				{
					bKeepTesting = false;
					break;
				}

				firstOperand = p_grown;
				lFirstSize = l_needed;
			}

			::memcpy( firstOperand, theAnswer, l_needed );
		}

		lOperandCount++;

		// Get the operator...
		cout << "operator       >> ";
		cout.flush();

		if ( !ReadLine( operation, kOperationSize ) )
		{
			bKeepTesting = false;
			break;
		}

		// Get the second operand...
		if ( NeedsSecondOperand( operation[ 0 ] ) )
		{
			cout << "second operand >> ";
			cout.flush();

			if ( !ReadLine( secondOperand, kOperandSize ) )
			{
				bKeepTesting = false;
				break;
			}
		}
		else
		{
			secondOperand[ 0 ] = '\0';
		}

		if ( secondOperand[ 0 ] != '\0' )
		{
			lOperandCount++;
		}

		// Clear the result number...
		result.Clear();

		// Do the operation...
		DoOperation( firstOperand, operation, secondOperand, &result );

		::free( theAnswer );
		theAnswer = nullptr;

		// Show results...
		result.TrimTrailing();
		theAnswer = result.GetStringRep();

		if ( theAnswer != nullptr )
		{
			if ( lOperandCount == 1 )
				printf( "\n%s( %s ) = %s\n\n", operation, firstOperand, theAnswer );

			if ( lOperandCount == 2 )
				printf( "\n%s %s %s = %s\n\n", firstOperand, operation, secondOperand, theAnswer );

			cout << "Digits array length ... " << result.m_lDigitsArrayLength << '\n';
			cout << "Digits used ........... " << result.m_lDigitsLength << '\n';
			cout << "Radix position ........ " << result.m_lRadix << "\n\n";
		}

	} while ( bKeepTesting );

	::free( firstOperand );
	firstOperand = nullptr;

	::free( secondOperand );
	secondOperand = nullptr;

	::free( theAnswer );
	theAnswer = nullptr;

    lzu_takedown();

	return 0;
}

/*
	One line of input into szBuffer, no more of it than lSize - 1
	characters, and false at the end of the input with nothing read.
	The loops this replaced stored getchar() in a char with no bound and
	no test for the end: an input ending without a newline wrote past the
	operator's sixteen bytes until it crashed, and a line longer than an
	operand's buffer overran that. A longer line is cut, and says so.

	Trailing white space is dropped, a carriage return with it: a script
	saved on Windows ends every line in one when it is piped in here, and
	the engine reads "2\r", like "2 ", as malformed -- nought, flagged --
	which the bench would show as a plain 0.
*/

bool ReadLine( char* szBuffer, long lSize )
{
	long	l_index = 0;
	long	l_dropped = 0;
	int		i_char = EOF;

	while ( ((i_char = getchar()) != EOF) && (i_char != '\n') )
	{
		if ( l_index < (lSize - 1) )
			szBuffer[ l_index++ ] = static_cast<char>( i_char );
		else
			l_dropped++;
	}

	while ( (l_index > 0) && ((szBuffer[ l_index - 1 ] == ' ') || (szBuffer[ l_index - 1 ] == '\t') || (szBuffer[ l_index - 1 ] == '\r')) )
		l_index--;

	szBuffer[ l_index ] = '\0';

	if ( l_dropped > 0 )
		fprintf( stderr, "(%ld characters past the first %ld were dropped)\n", l_dropped, lSize - 1 );

	return (i_char != EOF) || (l_index > 0) || (l_dropped > 0);
}

bool NeedsSecondOperand( char c_operand )
{
    bool b_needsSecondOperand = true;

	// Check the operation...
	switch ( c_operand )
	{
	case '|':	// absolute value
	case 'L':	// log via Taylor series
	case 'l':	// log via table
	case 'n':	// negate
	case '~':	// normalize
	case 'p':	// positive
	case 'B':	// set base
	case 'F':	// is fraction
	case 'I':	// is integer
	case 'M':	// order of magnitude
	case 'N':	// is negative
	case 'Z':	// is zero
	case '[':	// shift left
	case ']':	// shift right
	case '(':	// pad left
	case ')':	// pad right
    case 's':   // sine
    case 'c':   // cosine
    case 't':   // tangent
    case 'S':   // hyperbolic sine
    case 'C':   // hyperbolic cosine
    case 'T':   // hyperbolic tan
    case 'x':   // arcsin
    case 'X':   // arctan
    case 'y':   // arccos
    case 'j':   // arcsinh
    case 'h':   // arccosh
    case 'g':   // arctanh
    case 'o':   // cotangent
    case 'e':   // secant
    case 'k':   // cosecant
    case 'O':   // hyperbolic cotangent
    case 'E':   // hyperbolic secant
    case 'K':   // hyperbolic cosecant
    case 'v':   // arccot
    case 'V':   // arcsec
    case 'w':   // arccsc
    case 'W':   // arccoth
    case 'u':   // arcsech
    case 'U':   // arccsch
		b_needsSecondOperand = false;
		break;

	default:
		break;
	}

	return b_needsSecondOperand;
}

void DoOperation( char* firstOperand, char* operation, char* secondOperand, CLZNumber* result )
{
	CLZNumber	first;
	CLZNumber	second;

	long		lNewBase = 0;

	if ( glNewBase )
	{
		first.m_lBase = glNewBase;
		second.m_lBase = glNewBase;
	}

	first.Assign( firstOperand );

	if ( secondOperand[ 0 ] != '\0' )
	{
		second.Assign( secondOperand );
	}

	// Do the operation...
	switch ( operation[ 0 ] )
	{
	case '|':	// absolute value
		lz_abs( &first, result );
		break;

	case '+':	// add
		lz_add( &first, &second, result );
		break;

	case '?':	// compare
		result->Assign( lz_compare( &first, &second ) );
		break;

	case '/':	// divide
		lz_divide( &first, &second, result );
		break;

	case '>':	// max
		result->Assign( lz_max( &first, &second ) );
		break;

	case '<':	// min
		result->Assign( lz_min( &first, &second ) );
		break;

	case 'L':	// log e
		lz_log_e( &first, result );
		break;

	case 'l':	// log current base
		lz_log_b( &first, result );
		break;

	case 'n':	// negate
        lz_negate(&first, result);
		break;

	case '*':	// multiply
		lz_multiply( &first, &second, result );
		break;

    case '^':   // raise to power
        lz_exp_n( &first, &second, result );
        break;

	case '~':	// normalize
		lz_normalize( &first, result );
		break;

    case '%':	// modulo
        lz_modulo( &first, &second, result );
        break;

	case 'p':	// plus
		lz_plus( &first, result );
		break;

	case 'q':	// quantize
		lz_quantize( &first, &second, result );
		break;

	case 'R':	// round
		if ( secondOperand[ 0 ] == 'D' )
			first.m_rmMode = RM_round_down;
		else if ( secondOperand[ 0 ] == 'u' )
			first.m_rmMode = RM_round_half_up;
		else if ( secondOperand[ 0 ] == 'e' )
			first.m_rmMode = RM_round_half_even;
		else if ( secondOperand[ 0 ] == 'C' )
			first.m_rmMode = RM_round_ceiling;
		else if ( secondOperand[ 0 ] == 'F' )
			first.m_rmMode = RM_round_floor;
		else if ( secondOperand[ 0 ] == 'd' )
			first.m_rmMode = RM_round_half_down;
		else if ( secondOperand[ 0 ] == 'U' )
			first.m_rmMode = RM_round_up;
		else if ( secondOperand[ 0 ] == 'A' )
			first.m_rmMode = RM_round_alternate;
		else if ( secondOperand[ 0 ] == 'o' )
			first.m_rmMode = RM_round_half_odd;
		else if ( secondOperand[ 0 ] == 'R' )
			first.m_rmMode = RM_round_random;
		else
			first.m_rmMode = RM_round_default;

		lz_round( &first, result );
		break;

    case 's':   // sin
        lz_sine( &first, result );
        break;

    case 'c':   // cosine
        lz_cosine( &first, result );
        break;

    case 't':   // tangent
        lz_tangent( &first, result );
        break;

    case 'x':   // arcsin
        lz_arcsine( &first, result );
        break;

    case 'X':   // arctan
        lz_arctangent( &first, result );
        break;

    case 'y':   // arccos
        lz_arccosine( &first, result );
        break;

    case 'j':   // arcsinh
        lz_archypersine( &first, result );
        break;

    case 'h':   // arccosh
        lz_archypercosine( &first, result );
        break;

    case 'g':   // arctanh
        lz_archypertangent( &first, result );
        break;

    case 'S':   // hyperbolic sin
        lz_hypersine( &first, result );
        break;

    case 'C':   // hyperbolic cosine
        lz_hypercosine( &first, result );
        break;

    case 'T':   // hyperbolic tangent
        lz_hypertangent( &first, result );
        break;

    /*
        The reciprocal family. The cotangent, secant and cosecant were
        in the engine all along but had no key here; the other nine are
        new with them.
    */

    case 'o':   // cotangent
        lz_cotangent( &first, result );
        break;

    case 'e':   // secant
        lz_secant( &first, result );
        break;

    case 'k':   // cosecant
        lz_cosecant( &first, result );
        break;

    case 'O':   // hyperbolic cotangent
        lz_hypercotangent( &first, result );
        break;

    case 'E':   // hyperbolic secant
        lz_hypersecant( &first, result );
        break;

    case 'K':   // hyperbolic cosecant
        lz_hypercosecant( &first, result );
        break;

    case 'v':   // arccot
        lz_arccotangent( &first, result );
        break;

    case 'V':   // arcsec
        lz_arcsecant( &first, result );
        break;

    case 'w':   // arccsc
        lz_arccosecant( &first, result );
        break;

    case 'W':   // arccoth
        lz_archypercotangent( &first, result );
        break;

    case 'u':   // arcsech
        lz_archypersecant( &first, result );
        break;

    case 'U':   // arccsch
        lz_archypercosecant( &first, result );
        break;

	case '0':	// round integral
		if ( secondOperand[ 0 ] == 'D' )
			first.m_rmMode = RM_round_down;
		else if ( secondOperand[ 0 ] == 'u' )
			first.m_rmMode = RM_round_half_up;
		else if ( secondOperand[ 0 ] == 'e' )
			first.m_rmMode = RM_round_half_even;
		else if ( secondOperand[ 0 ] == 'C' )
			first.m_rmMode = RM_round_ceiling;
		else if ( secondOperand[ 0 ] == 'F' )
			first.m_rmMode = RM_round_floor;
		else if ( secondOperand[ 0 ] == 'd' )
			first.m_rmMode = RM_round_half_down;
		else if ( secondOperand[ 0 ] == 'U' )
			first.m_rmMode = RM_round_up;
		else if ( secondOperand[ 0 ] == 'A' )
			first.m_rmMode = RM_round_alternate;
		else if ( secondOperand[ 0 ] == 'o' )
			first.m_rmMode = RM_round_half_odd;
		else if ( secondOperand[ 0 ] == 'R' )
			first.m_rmMode = RM_round_random;
		else
			first.m_rmMode = RM_round_default;

		lz_round_integral( &first, result );
		break;

	case '9':	// round last digit
		if ( secondOperand[ 0 ] == 'D' )
			first.m_rmMode = RM_round_down;
		else if ( secondOperand[ 0 ] == 'u' )
			first.m_rmMode = RM_round_half_up;
		else if ( secondOperand[ 0 ] == 'e' )
			first.m_rmMode = RM_round_half_even;
		else if ( secondOperand[ 0 ] == 'C' )
			first.m_rmMode = RM_round_ceiling;
		else if ( secondOperand[ 0 ] == 'F' )
			first.m_rmMode = RM_round_floor;
		else if ( secondOperand[ 0 ] == 'd' )
			first.m_rmMode = RM_round_half_down;
		else if ( secondOperand[ 0 ] == 'U' )
			first.m_rmMode = RM_round_up;
		else if ( secondOperand[ 0 ] == 'A' )
			first.m_rmMode = RM_round_alternate;
		else if ( secondOperand[ 0 ] == 'o' )
			first.m_rmMode = RM_round_half_odd;
		else if ( secondOperand[ 0 ] == 'R' )
			first.m_rmMode = RM_round_random;
		else
			first.m_rmMode = RM_round_default;

		lz_round_last_digit( &first, result );
		break;

	case '=':	// assign
		result->Assign( lz_same_quantum( &first, &second ) );
		break;

	case '-':	// subtract
		lz_subtract( &first, &second, result );
		break;

	case 'B':	// set base
		if ( firstOperand[ 0 ] == '2' )
			lNewBase = 2;
		else if ( firstOperand[ 0 ] == '8' )
			lNewBase = 8;
		else if ( firstOperand[ 0 ] == '1' )
		{
			if ( firstOperand[ 1 ] == '0' )
				lNewBase = 10;
			else if ( firstOperand[ 1 ] == '2' )
				lNewBase = 12;
			else if ( firstOperand[ 1 ] == '6' )
				lNewBase = 16;
		}

		if ( lNewBase )
		{
			first.SetBase( lNewBase );
			second.SetBase( lNewBase );
			result->SetBase( lNewBase );

			lNewBase = 0;
			glNewBase = 0;
		}

		break;

	case 'F':	// is fraction
		result->Assign( &first );

		if ( first.IsFraction() )
			printf( "%s", "is fraction:  " );
		else
			printf( "%s", "is not fraction:  " );

		break;

	case 'I':	// is integer
		result->Assign( &first );

		if ( first.IsInteger() )
			printf( "%s", "is integer:  " );
		else
			printf( "%s", "is not integer:  " );

		break;

	case 'M':	// order of magnitude
		result->Assign( first.GetOrderOfMagnitude() );
		printf( "%s", "order of magnitude:  " );
		break;

	case 'N':	// is negative
		result->Assign( &first );

		if ( first.IsNegative() )
			printf( "%s", "is negative:  " );
		else
			printf( "%s", "is not negative:  " );

		break;

	case ':':	// trim
		result->Assign( &first );

		if ( secondOperand[ 0 ] == '1' )
			result->TrimLeading();
		else if ( secondOperand[ 0 ] == '2' )
			result->TrimTrailing();
		else
			result->TrimAll();

		break;

	case 'Z':	// is zero
		result->Assign( &first );

		if ( first.IsZero() )
			printf( "%s", "is zero:  " );
		else
			printf( "%s", "is not zero:  " );

		break;

	case '[':	// shift left
		result->Assign( &first );
		result->ShiftLeft( 1 );

		break;

	case ']':	// shift right
		result->Assign( &first );
		result->ShiftRight( 1 );

		break;

	case '(':	// pad left
		result->Assign( &first );
		result->PadLeft( 1 );

		break;

	case ')':	// pad right
		result->Assign( &first );
		result->PadRight( 1 );

		break;
	}
}