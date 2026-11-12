/////////////////////////////////////////////////////////////////
//
//	LZMath.h
//
//	Part of Addaptic
//
//	Copyright (c) 2006-2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#pragma once

#include "LZNumber.h"

/////////////////////////////////////////////////////////////////
//  Some Useful Constants

/*
	The trigonometric argument reduction subtracts these, so their
	width sets the ceiling on how accurate a reduced argument - and
	therefore a sine - can be.

	kTwoPi carries 120 places for the reason kBaseE does, and it is a
	sharper reason. The reduction takes x modulo 2*pi, which multiplies
	2*pi by a quotient the size of the argument, so **every digit of
	the argument's magnitude costs a digit of the constant**: at 50
	places the sine of 1e20 held 30 digits, 1e30 held 19, and 1e40 held
	10, falling exactly one for one. Those were the thirty cases the
	sweep set aside as constant-width, and the only thirty it had.

	120 places carries the sweep's whole ladder to full width with room
	over, and lz_sine narrows it per call to the width that call needs
	-- see the note there, without which the whole circular family
	slowed by a third or more.

	**The other three carry 120 now as well, and the order that
	happened in is the lesson.** They are subtracted once each from a
	value the modulo has already brought into [0, 2*pi), so nothing
	multiplies their error the way the reduction multiplies 2*pi's --
	which is why widening them was measured first and found to buy
	nothing: at 100 places the sine held 42 to 50 digits with them at
	50, and 42 to 51 with them at 120.

	That measurement was right and its conclusion was wrong, because
	something else was binding in front of them. The series stopped at
	x^32/32! and the truncation, not any constant, was what held those
	42 to 50. Once lzi_trig_series could carry the series past the
	table, the same test came back 50 flat -- and then widening these
	took it to the full working width. **A constant that measures as
	not binding may only be standing behind one that is.**
*/

const char kPiOverFour[] =
		"0.7853981633974483096156608458198757210492923498437764552437"
		"361480769541015715522496570087063355292669955370216283205766"
		"62";
const char kPiOverTwo[] =
		"1.5707963267948966192313216916397514420985846996875529104874"
		"722961539082031431044993140174126710585339910740432566411533"
		"24";
const char kPi[] =
		"3.1415926535897932384626433832795028841971693993751058209749"
		"445923078164062862089986280348253421170679821480865132823066"
		"47";
const char kPiSquared[]			= "9.8696044010893586188344909998761511353136994072408";

const char kTwoPi[] =
		"6.2831853071795864769252867665590057683943387987502116419498"
		"891846156328125724179972560696506842341359642961730265646132"
		"94";

/////////////////////////////////////////////////////////////////
//  Constants for Log Calculations

/*
	ln of every base SetBase accepts, for lz_exp_b and lz_log_b.

	**120 places, where they were fifty until 2026-09-24.** lz_exp_b
	multiplies the fractional part of its exponent by one of these and
	takes the exponential, then multiplies that by a whole power of the
	base -- so the constant's relative error becomes an absolute error
	the size of the answer. At fifty digits 10^20.5 stopped at 49
	significant digits where 36 places of it afford 57, and 10^29.5 at
	49 where it affords 66, once lz_exp_b had begun widening by its
	answer's size and every other limit had gone. The same width kBaseE
	and the pi family carry, and for the same reason.

	**And a ceiling at 120 until 2026-09-26**, as kBaseE and pi had
	been: past them lz_exp_b and lz_log_b read these regardless, and
	were out by 1e20 units at 120 places and 1e50 at 150. Past their
	width lzi_log_base_in_base computes ln(base) now, and within it
	hands these back exactly as before.

	Until 2026-09-24 only 2, 8, 10, 12 and 16 had one, and the other ten
	bases fell through to ln 10. All fifteen were generated with Python's
	decimal and checked against an independent atanh series to 1e-190
	and against ln(ab) = ln a + ln b, as every one of these should be.
*/

const char kLog2BaseE[] =
		".69314718055994530941723212145817656807550013436025525412068"
		"000949339362196969471560586332699641868754200148102057068573"
		"4";

const char kLog3BaseE[] =
		"1.0986122886681096913952452369225257046474905578227494517346"
		"943336374942932186089668736157548137320887879700290659578657"
		"42";

const char kLog4BaseE[] =
		"1.3862943611198906188344642429163531361510002687205105082413"
		"600189867872439393894312117266539928373750840029620411413714"
		"67";

const char kLog5BaseE[] =
		"1.6094379124341003746007593332261876395256013542685177219126"
		"478914741789877076577646301338780931796107999663030217155629"
		"00";

const char kLog6BaseE[] =
		"1.7917594692280550008124773583807022727229906921830047058553"
		"743431308879151883036824794790818101507763299715100865285514"
		"76";

const char kLog7BaseE[] =
		"1.9459101490553133051053527434431797296370847295818611884593"
		"901499375798627520692677876584985878715269930616942058511409"
		"12";

const char kLog8BaseE[] =
		"2.0794415416798359282516963643745297042265004030807657623620"
		"400284801808659090841468175899809892560626260044430617120572"
		"01";

const char kLog9BaseE[] =
		"2.1972245773362193827904904738450514092949811156454989034693"
		"886672749885864372179337472315096274641775759400581319157314"
		"85";

const char kLog10BaseE[] =
		"2.3025850929940456840179914546843642076011014886287729760333"
		"279009675726096773524802359972050895982983419677840422862486"
		"33";

const char kLog11BaseE[] =
		"2.3978952727983705440619435779651292998217068539374171752185"
		"677091305736239132367130750547080026347914147157258881379985"
		"22";

const char kLog12BaseE[] =
		"2.4849066497880003102297094798388788407984908265432599599760"
		"543526242815371579983980853424088065694638719729911070992372"
		"10";

const char kLog13BaseE[] =
		"2.5649493574615367360534874415653186048052679447602071164190"
		"455106634646673244101793995746634404894887692581920927627216"
		"32";

const char kLog14BaseE[] =
		"2.6390573296152586145225848649013562977125848639421164425800"
		"701594309734847217639833935218255842902145350631752264218266"
		"45";

const char kLog15BaseE[] =
		"2.7080502011022100659960045701487133441730919120912671736473"
		"422251116732809262667315037496329069116995879363320876734286"
		"42";

const char kLog16BaseE[] =
		"2.7725887222397812376689284858327062723020005374410210164827"
		"200379735744878787788624234533079856747501680059240822827429"
		"35";


/////////////////////////////////////////////////////////////////
//  Constants for Exponent Calculations

/*
	e to 120 places, where every other constant here carries 50.

	It is wider than its neighbours because it was measurably the
	limit and they are not. lz_raise_e squares this to reach e^n, and
	the answer can be no better than the constant it started from:
	at 49 places, exp(100) stopped at 47 significant digits where the
	representation could hold 80, and ln of anything needing reduction
	stopped at 49 however wide the caller set the precision.

	Widening it alone is not enough and not the point -- the squaring
	runs at the working precision, so it discards whatever the
	constant carries past that. See kRaiseEGuardDigits, which had to
	move with it.

	Generated and checked two ways before being written down: ln(e)
	back to 1 to 170 places, and against the series summed to 140
	terms.
*/
const char kBaseE[] =
		"2.7182818284590452353602874713526624977572470936999595749669"
		"676277240766303535475945713821785251664274274663919320030599"
		"22";

// The highest index held in gRecipFact, and so the highest power any
//	series summed against that table can reach. 32 until 2026-09-26,
//	when eight more took the circular series from 40 digits to 48 --
//	see lzi_trig_series. The exponential's Horner loop still starts at
//	32, the reaches it was measured with...
const long kRecipFactMax = 40;

const char kRecip00Fact[]	= "1";	//	0
const char kRecip01Fact[]	= "1";	//	1
const char kRecip02Fact[]	= ".5";	//	2
const char kRecip03Fact[]	= ".16666666666666666666666666666666666666666666666667";	//	3
const char kRecip04Fact[]	= ".041666666666666666666666666666666666666666666666667";	//	4
const char kRecip05Fact[]	= ".0083333333333333333333333333333333333333333333333333";	//	5
const char kRecip06Fact[]	= ".0013888888888888888888888888888888888888888888888889";	//	6
const char kRecip07Fact[]	= ".00019841269841269841269841269841269841269841269841270";	//	7
const char kRecip08Fact[]	= ".000024801587301587301587301587301587301587301587301587";	//	8
const char kRecip09Fact[]	= ".0000027557319223985890652557319223985890652557319223986";	//	9
const char kRecip10Fact[]	= ".00000027557319223985890652557319223985890652557319223986";	//	10
const char kRecip11Fact[]	= ".000000025052108385441718775052108385441718775052108385442";	//	11
const char kRecip12Fact[]	= ".0000000020876756987868098979210090321201432312543423654535";	//	12
const char kRecip13Fact[]	= ".00000000016059043836821614599392377170154947932725710503488";	//	13
const char kRecip14Fact[]	= ".000000000011470745597729724713851697978682105666232650359634";	//	14
const char kRecip15Fact[]	= ".00000000000076471637318198164759011319857880704441551002397563";	//	15
const char kRecip16Fact[]	= ".000000000000047794773323873852974382074911175440275969376498477";	//	16
const char kRecip17Fact[]	= ".0000000000000028114572543455207631989455830103200162334927352045";	//	17
const char kRecip18Fact[]	= ".00000000000000015619206968586226462216364350057333423519404084470";	//	18
const char kRecip19Fact[]	= ".0000000000000000082206352466243297169559812368722807492207389918261";	//	19
const char kRecip20Fact[]	= ".00000000000000000041103176233121648584779906184361403746103694959131";	//	20
const char kRecip21Fact[]	= ".000000000000000000019572941063391261230847574373505430355287473790062";	//	21
const char kRecip22Fact[]	= ".00000000000000000000088967913924505732867488974425024683433124880863919";	//	22
const char kRecip23Fact[]	= ".000000000000000000000038681701706306840377169119315228123231793426462573";	//	23
const char kRecip24Fact[]	= ".0000000000000000000000016117375710961183490487133048011718013247261026072";	//	24
const char kRecip25Fact[]	= ".000000000000000000000000064469502843844733961948532192046872052989044104289";	//	25
const char kRecip26Fact[]	= ".0000000000000000000000000024795962632247974600749435458479566174226555424727";	//	26
const char kRecip27Fact[]	= ".000000000000000000000000000091836898637955461484257168364739133978616871943432";	//	27
const char kRecip28Fact[]	= ".0000000000000000000000000000032798892370698379101520417273121119278077454265511";	//	28
const char kRecip29Fact[]	= ".00000000000000000000000000000011309962886447716931558764576938316992440501470866";	//	29
const char kRecip30Fact[]	= ".0000000000000000000000000000000037699876288159056438529215256461056641468338236220";	//	30
const char kRecip31Fact[]	= ".00000000000000000000000000000000012161250415535179496299746856922921497247851043942";	//	31
const char kRecip32Fact[]	= ".0000000000000000000000000000000000038003907548547435925936708927884129678899534512318";	//	32
const char kRecip33Fact[]	= ".00000000000000000000000000000000000011516335620771950280586881493298221114818040761309";	//	33
const char kRecip34Fact[]	= ".0000000000000000000000000000000000000033871575355211618472314357333230062102406002239143";	//	34
const char kRecip35Fact[]	= ".000000000000000000000000000000000000000096775929586318909920898163809228748864017149254694";	//	35
const char kRecip36Fact[]	= ".0000000000000000000000000000000000000000026882202662866363866916156613674652462226985904082";	//	36
const char kRecip37Fact[]	= ".000000000000000000000000000000000000000000072654601791530713153827450307228790438451313254275";	//	37
const char kRecip38Fact[]	= ".0000000000000000000000000000000000000000000019119632050402819251007223765060208010118766645862";	//	38
const char kRecip39Fact[]	= ".000000000000000000000000000000000000000000000049024697565135433976941599397590276949022478579133";	//	39
const char kRecip40Fact[]	= ".0000000000000000000000000000000000000000000000012256174391283858494235399849397569237255619644783";	//	40


/////////////////////////////////////////////////////////////////
//	Constants for Inverse Trig Calculations

/////////////////////////////////////////////////////////////////
//  Constants for the Gamma Function

/*
	The even Bernoulli numbers B_2 through B_60, for Stirling's series
	in lzi_log_gamma. The odd ones past B_1 are all nought and are not
	held.

	**The count is not a ceiling on the precision, which is why there
	is a table here at all.** Stirling is asymptotic: at a fixed
	argument, adding terms stops helping and eventually hurts. But the
	recurrence lnGamma(z) = lnGamma(z+m) - ln(z(z+1)...(z+m-1)) moves
	the argument up as far as wanted for the price of m multiplications
	and one logarithm, and a larger argument makes every term smaller.
	So accuracy is bought with the shift, and the table only decides
	how dear the shift is:

	    terms   shift for 40 digits   for 60   for 80   for 100
	      12            55              345      --        --
	      20            25               70      205       --
	      30            20               35       75       150

	Thirty is where the shift stops running away. Twelve cannot reach
	80 digits at any shift worth paying for, and twenty gives out
	before a hundred. Past thirty the shift barely moves and the table
	is just longer.

	Generated by two independent exact derivations -- Akiyama-Tanigawa
	and the binomial recurrence -- agreeing on all thirty as exact
	rationals before being rounded to fifty digits.

	**The count is no ceiling, but the width of the entries was, until
	2026-09-26.** Held to fifty digits, they move lnGamma by 1e-53 to
	7e-53 whatever the shift, so it stopped at 52 or 53 places however
	wide the caller asked -- measured 2026-09-24 at precisions 60, 80
	and 100 -- and every large gamma at 52 significant digits in every
	base, the default included. Past kBernoulliTableReach lzi_log_gamma
	sums over the exact rationals below instead, and the table above,
	which assumes exact numbers, is reached: correctly rounded to 150
	places, measured. What was left was the table's last row: thirty
	terms need the shift to grow tenfold every 61 digits, so past about
	225 places a gamma took seconds, and past 300 the shift ran to
	hundreds of thousands. Past 120 digits lzi_log_gamma now sums more
	terms than these, from the tangent numbers -- see
	lzi_stirling_rational and kGammaTermsPer100Digits.
*/

// The highest index held in gBernoulli, so the highest term Stirling's
//	series can reach from the table...
const long kBernoulliMax = 29;

/*
	**The same thirty exactly**, as numerator over denominator, for
	lzi_log_gamma past kBernoulliTableReach -- where the fifty-digit
	entries above were the wall the note describes. Generated by the
	same two derivations, Akiyama-Tanigawa and the binomial recurrence,
	which agree on all thirty as exact rationals, and within 3.5e-50 of
	every entry above.

	Whole numbers convert between bases exactly, so these serve any
	width in any base. The largest denominator fits a 32-bit long;
	multiplied by Stirling's 2n(2n-1) it does not, which is why
	lzi_bernoulli_rational builds that product as a number.
*/

const char* const kBernoulliNumerators[ kBernoulliMax + 1 ] =
{
	"1",	// B_2
	"-1",	// B_4
	"1",	// B_6
	"-1",	// B_8
	"5",	// B_10
	"-691",	// B_12
	"7",	// B_14
	"-3617",	// B_16
	"43867",	// B_18
	"-174611",	// B_20
	"854513",	// B_22
	"-236364091",	// B_24
	"8553103",	// B_26
	"-23749461029",	// B_28
	"8615841276005",	// B_30
	"-7709321041217",	// B_32
	"2577687858367",	// B_34
	"-26315271553053477373",	// B_36
	"2929993913841559",	// B_38
	"-261082718496449122051",	// B_40
	"1520097643918070802691",	// B_42
	"-27833269579301024235023",	// B_44
	"596451111593912163277961",	// B_46
	"-5609403368997817686249127547",	// B_48
	"495057205241079648212477525",	// B_50
	"-801165718135489957347924991853",	// B_52
	"29149963634884862421418123812691",	// B_54
	"-2479392929313226753685415739663229",	// B_56
	"84483613348880041862046775994036021",	// B_58
	"-1215233140483755572040304994079820246041491",	// B_60
};

const long kBernoulliDenominators[ kBernoulliMax + 1 ] =
{
	6L,
	30L,
	42L,
	30L,
	66L,
	2730L,
	6L,
	510L,
	798L,
	330L,
	138L,
	2730L,
	6L,
	870L,
	14322L,
	510L,
	6L,
	1919190L,
	6L,
	13530L,
	1806L,
	690L,
	282L,
	46410L,
	66L,
	1590L,
	798L,
	870L,
	354L,
	56786730L
};

// The widest the fifty-digit entries serve, in decimal digits. Past it
//	lzi_log_gamma sums over the exact ones...
const long kBernoulliTableReach = 50;

/*
	Half the natural logarithm of 2*pi, the constant term of Stirling's
	series. 120 places, for a caller who raises the precision; it is
	added once per call, which is not enough to be worth narrowing --
	measured, as the three pi subtractions in the inverse functions
	were. Past its 120 digits lzi_half_log_two_pi computes it.

	Checked by exponentiating twice back to 2*pi.
*/

const char kHalfLogTwoPi[] =
		"0.9189385332046727417803297364056176398613974736377834128171"
		"515404827656959272603976947432986359541976220056466246343374"
		"46";

/*
	How far lzi_log_gamma shifts its argument up before summing
	Stirling's series at the default working width, and the width past
	which it asks the first term left out instead.

	Twenty reaches 40 digits with the full table; twenty-four is that
	with a little room. Above the default the shift has to grow. It
	grew half as much again per ten places until 2026-09-26, which is
	1.3 to 3 times what the Bernoulli note's measurements need -- 3079
	at 150 places where thirty terms need 1340 -- and erring high was
	said to cost little because the shift is paid for in
	multiplications. Below a hundred places it did; past them the
	product runs to thousands of digits and a gamma took seconds. It
	is now the bound from the first term left out -- see
	kGammaShiftMargin and lzi_log_gamma -- which reproduces those
	measurements, and never less than twenty-four.

	**Both are decimal digits.** lzi_log_gamma converts the working
	precision, which is in places of the argument's base, before
	comparing it with kGammaShiftFrom: 44 hex places are 53 digits and
	44 binary ones 13. Compared as places, base sixteen got a series
	short of what it needed, and base two, once lz_gamma had widened
	it by a hundred binary places, a target of over a thousand --
	gamma(30.5) there took 9.5 seconds rather than 0.28.
*/

const long kGammaBaseShift	= 24;
const long kGammaShiftFrom	= 44;

// How far under the working width, in decimal digits, lzi_log_gamma
//	puts the first term it leaves out...
const long kGammaShiftMargin	= 2;

/*
	Past 120 digits, how many of Stirling's terms lzi_log_gamma sums for
	every hundred digits of working width, the terms past the thirtieth
	coming from the tangent numbers. Fewer terms need a larger shift, and
	more cost a longer tangent run and a longer sum. Measured for
	log-gamma at 300, 500 and 1000 places, a first call and then a warm
	one: 15 took 0.42, 1.3 and 6.3 s warm, the shift growing again; 25
	took 8, 24 and 122 ms; 35 was as fast warm and slower first; 50
	slower both ways.
*/

const long kGammaTermsPer100Digits	= 25;

/*
	The guard digits lzi_log_gamma works to.

	Stirling's series is summed at the shifted argument and then a
	logarithm of the shifted product is taken off it, and that
	subtraction is where a gamma loses its places: the two sides are
	close for a small argument. Measured at 4, 8 and 16.
*/

const long kGammaGuardDigits	= 8;

/*
	The most lz_gamma will widen itself by for a large answer.

	Thirty, which is where auditrounding.py stops asking for 36 places
	as well: past 1e30 an answer has more integer digits than the
	fraction is worth, and the shift needed to serve them makes the
	rising product a number of thousands of digits.

	**Thirty decimal places, converted to the answer's base** --
	about 100 binary places, 34 octal, 25 hex. Applied as it stood it
	allowed thirty binary places, nine digits, and a large gamma in
	base two stopped at 20 digits where base ten reached 52.

	This said Gamma(20.5) gets its full 54 digits and Gamma(100.5) 66.
	Measured 2026-09-24 both got 52: the Bernoulli numbers' fifty-digit
	width stopped lnGamma first -- see the note on the table. Since
	2026-09-26 Gamma(20.5) is correctly rounded and Gamma(100.5) gets
	72, the thirty and lzi_log_gamma's guard digits.
*/

const long kGammaMaxWidening	= 30;

const char kBernoulli00[]	= ".16666666666666666666666666666666666666666666666667";	//	B_2
const char kBernoulli01[]	= "-.033333333333333333333333333333333333333333333333333";	//	B_4
const char kBernoulli02[]	= ".023809523809523809523809523809523809523809523809524";	//	B_6
const char kBernoulli03[]	= "-.033333333333333333333333333333333333333333333333333";	//	B_8
const char kBernoulli04[]	= ".075757575757575757575757575757575757575757575757576";	//	B_10
const char kBernoulli05[]	= "-.25311355311355311355311355311355311355311355311355";	//	B_12
const char kBernoulli06[]	= "1.1666666666666666666666666666666666666666666666667";	//	B_14
const char kBernoulli07[]	= "-7.0921568627450980392156862745098039215686274509804";	//	B_16
const char kBernoulli08[]	= "54.971177944862155388471177944862155388471177944862";	//	B_18
const char kBernoulli09[]	= "-529.12424242424242424242424242424242424242424242424";	//	B_20
const char kBernoulli10[]	= "6192.1231884057971014492753623188405797101449275362";	//	B_22
const char kBernoulli11[]	= "-86580.253113553113553113553113553113553113553113553";	//	B_24
const char kBernoulli12[]	= "1425517.1666666666666666666666666666666666666666667";	//	B_26
const char kBernoulli13[]	= "-27298231.067816091954022988505747126436781609195402";	//	B_28
const char kBernoulli14[]	= "601580873.90064236838430386817483591677140064236838";	//	B_30
const char kBernoulli15[]	= "-15116315767.092156862745098039215686274509803921569";	//	B_32
const char kBernoulli16[]	= "429614643061.16666666666666666666666666666666666667";	//	B_34
const char kBernoulli17[]	= "-13711655205088.332772159087948561632772159087948562";	//	B_36
const char kBernoulli18[]	= "488332318973593.16666666666666666666666666666666667";	//	B_38
const char kBernoulli19[]	= "-19296579341940068.148632668144863266814486326681449";	//	B_40
const char kBernoulli20[]	= "841693047573682615.00055370985603543743078626799557";	//	B_42
const char kBernoulli21[]	= "-40338071854059455413.076811594202898550724637681159";	//	B_44
const char kBernoulli22[]	= "2115074863808199160560.1453900709219858156028368794";	//	B_46
const char kBernoulli23[]	= "-120866265222965259346027.31193708252531781943546649";	//	B_48
const char kBernoulli24[]	= "7500866746076964366855720.0757575757575757575757576";	//	B_50
const char kBernoulli25[]	= "-503877810148106891413789303.05220125786163522012579";	//	B_52
const char kBernoulli26[]	= "36528776484818123335110430842.971177944862155388471";	//	B_54
const char kBernoulli27[]	= "-2849876930245088222626914643291.0678160919540229885";	//	B_56
const char kBernoulli28[]	= "238654274996836276446459819192192.14971751412429379";	//	B_58
const char kBernoulli29[]	= "-21399949257225333665810744765191097.392674151161724";	//	B_60


const char kInvTrigA00[]	= "1";
const char kInvTrigA01[]	= ".16666666666666666666666666666666666666666666666667";
const char kInvTrigA02[]	= ".075";
const char kInvTrigA03[]	= ".044642857142857142857142857142857142857142857142857";
const char kInvTrigA04[]	= ".030381944444444444444444444444444444444444444444444";
const char kInvTrigA05[]	= ".022372159090909090909090909090909090909090909090909";
const char kInvTrigA06[]	= ".017352764423076923076923076923076923076923076923077";
const char kInvTrigA07[]	= ".01396484375";
const char kInvTrigA08[]	= ".011551800896139705882352941176470588235294117647059";
const char kInvTrigA09[]	= ".0097616095291940789473684210526315789473684210526316";
const char kInvTrigA10[]	= ".0083903358096168154761904761904761904761904761904762";
const char kInvTrigA11[]	= ".0073125258735988451086956521739130434782608695652174";
const char kInvTrigA12[]	= ".0064472103118896484375";
const char kInvTrigA13[]	= ".0057400376708419234664351851851851851851851851851852";
const char kInvTrigA14[]	= ".0051533096823199041958512931034482758620689655172414";
const char kInvTrigA15[]	= ".0046601434869150961599042338709677419354838709677419";
const char kInvTrigA16[]	= ".0042409070936793630773370916193181818181818181818182";
const char kInvTrigA17[]	= ".0038809645588376692363194056919642857142857142857143";
const char kInvTrigA18[]	= ".0035692053938259345454138678473395270270270270270270";
const char kInvTrigA19[]	= ".0032970595034734847453924325796274038461538461538462";
const char kInvTrigA20[]	= ".0030578216492580306693548109473251714939024390243902";
const char kInvTrigA21[]	= ".0028461784011089421678767647854117460029069767441860";
const char kInvTrigA22[]	= ".0026578706382072899335374434789021809895833333333333";
//const char kInvTrigA22[]	= ".12050598141624290916106116114875";

// The highest index held in gInvTrigA, and so the highest power the
//	inverse sine series can reach...
const long kInvTrigMax = 22;

// How small the sine has to be before that series is summed, and how
//	many half-angle steps are allowed to get it there. Below a tenth,
//	the first term dropped is smaller than 1e-45...
const char kArcSineThreshold[]	= ".1";
const long kArcSineMaxHalvings	= 8;

// The widest the coefficient table serves, in decimal digits: its 50
//	digits, and a first dropped term of about 1e-48 against a sine under a
//	tenth. Past it lz_arcsine computes its terms. At the default a base-
//	ten arcsine works to 41 digits and a hex one to 49.4, both inside, so
//	nothing at the default takes the new path...
const long kArcSineTableReach	= 50;

/*
	Guard digits for the arcsine's working.

	Precision here counts places after the radix mark, and the halving
	sums the series against a value below a tenth before doubling the
	answer back up. The series result is therefore held to the same
	number of places as an answer 2^h times its size, and undoing the
	halving multiplies its error by that 2^h - h*log10(2) places, which
	is the two digits every argument that halves at all used to give
	up. The answer is rounded back to the working precision afterwards,
	so these places never leave the function.

	**Five, and two was wrong.** Two was chosen because 2, 3, 4 and 6
	"all score identically across the domain" -- in digits of agreement,
	which cannot see an error under a unit. Up to eight halvings
	multiply the series' error by 256, two and a half digits, which two
	guard digits cannot hold. Scored for correct rounding over 3,500
	random arguments in seven bases at the default, misses went 126 at
	two, 9 at three, 3 at four and 1 at five, six and eight -- that one
	3.5e-14 of a unit from its rounding boundary, which no guard sees.
	Four still missed in base twelve at 40 and 44 places, where five
	missed nothing short of the series' own ceiling. It costs 5 to 8%
	in base ten and 20% in base two over two. (Those counts are after
	the fix to cos t near one -- see lz_arcsine -- without which half
	the arguments near one missed whatever the guard.)
*/
const long kArcSineGuardDigits	= 5;

/*
	The guard digits lz_sine and lz_cosine work at.

	The series rounds every term at the working precision and sums
	thirty of them, and the reduction subtracts pi/2 and pi/4 before
	that, rounding again each time. None of it is large, but it
	accumulated: measured against the correctly rounded answer, the
	circular functions were landing a digit or two short of what their
	own answers could hold, and everything built on them inherited it.

	Four, because four is where it stops. Counting how many of the
	sweep's 1542 scoreable cases come back correctly rounded to 36
	places: no guard 1127, two 1324, three 1336, four 1338, six 1338.
	Six buys nothing at all, and the whole family costs 6 to 24 per
	cent more time for it.

	Three is the better bargain and is worth knowing about -- it takes
	all but two of the cases at essentially two's cost, the step from
	three to four being 4 per cent of the time for those two. Four is
	kept because it is where the curve flattens, which is the rule the
	other guard constants here were chosen by.

	The answer is rounded back afterwards, so these places never leave
	the function.

	This does nothing for an argument already inside [0, pi/4] and
	answering at its ceiling -- sin(.2) and cos(.2) do not move. What
	it recovers is the reduction's roundings and the series' own, which
	is why the gain shows on the reduced arguments and, most of all, on
	the reciprocals that divide two of these together.
*/
const long kTrigGuardDigits		= 4;

/*
	The guard digits lz_tangent, lz_cotangent, lz_secant and lz_cosecant
	work at: two more than the sine's, since 2026-10-03.

	A quotient spends its guard twice over. With the sine and cosine each
	right to half a unit of the guard's last place, the quotient's error
	is about (1 + |answer|) / |divisor| times that -- 110 for a tangent
	near ten, whose cosine is near a tenth -- and the small-divisor
	retake does not begin until the divisor is under a hundredth. Four
	digits left an ordinary tangent, cotangent or secant blind within
	about 10^-3 of a unit of a half: 11 of 30,000 random operands in
	bases ten and sixteen misrounded, by up to 0.0025 of a unit.

	Measured, by guard: five took that to 1 of 30,000, six to none, the
	blind zone to about 4e-6 of a unit, level with the sine's -- and
	near the poles from 51 wrong of 960 answers a hair from a half to 7.
	Seven buys little more and cost the cotangent a third at the default,
	the series crossing its table threshold. Six costs 2 to 8 per cent
	on an ordinary argument at 36 places and nothing at 200, and 17 to
	24 where the cosine is already taken again, near pi/2, whose sine
	then runs past the table.

	**The sweep cannot see any of this**, which is why "six buys nothing"
	above was a fair reading of it: its arguments never land within a
	thousandth of a unit of a half.
*/
const long kTrigQuotientGuardDigits	= 6;

/*
	**The doubt windows of the retake**, in units of the answer's last
	place: a function whose final rounding discards a guard lying within
	its window of the half is taken again wider -- see
	CLZAnswerBase::Retake in LZMath.cpp. Since 2026-10-03.

	Each is about ten times the largest error measured in the guard's
	digits, against exact values, over the near-half census's answers and
	random ones -- some 30,000 in all, in bases ten, sixteen, two and
	seven, at 36 places and again at 80, where the errors were no larger:

	  e to a whole power			5.9e-8
	  the logarithm					1.6e-8
	  the cosecant					6.9e-7
	  log-gamma						4.9e-7
	  gamma							2.9e-6
	  the harmonic mean				1.7e-5
	  most of the rest				2.5e-5 to 9.4e-5,
	    the geometric mean and x to a whole power among them
	  the secant					1.1e-4
	  the AGM						1.5e-4
	  the sine						2.8e-4
	  the arcsine					4.3e-4
	  the exponential				7.6e-4

	A wider window costs only retakes: about twice its width of a
	function's calls are taken again. A narrower one than the error
	misses roundings, as every function did before there was a window.
*/
const double kRetakeDoubtLog		= 1.0e-6;
const double kRetakeDoubtRaiseE	= 1.0e-6;
const double kRetakeDoubtCosecant	= 1.0e-5;
const double kRetakeDoubtLogGamma	= 5.0e-6;
const double kRetakeDoubtGamma		= 3.0e-5;
const double kRetakeDoubt			= 1.0e-3;
const double kRetakeDoubtWide		= 2.0e-3;
const double kRetakeDoubtSine		= 3.0e-3;
const double kRetakeDoubtArcSine	= 5.0e-3;
const double kRetakeDoubtExp		= 8.0e-3;

/*
	The guard digits lz_raise_e works at, and why this one is four
	times the trigonometric guard.

	kBaseE was written to 49 places when this was found. Squaring it at
	36 threw thirteen of them away on the very first step -- e squared
	came back 37 digits wide where its operand was 50 -- and every
	squaring after doubled whatever relative error was already there.
	Six squarings reach e^64, so the thirteen places lost at the start
	arrived at the answer multiplied by about sixty. That was the whole
	of it: e^50 scored 36 significant digits where the representation
	could hold 58.

	Sixteen, because thirteen was what kBaseE's width demanded and the
	squaring wants a few more on top. Measured from .5 to 200, 12 and
	16 both reach the ceiling and 20 and 24 add nothing.

	kBaseE carries 120 places now, so the constant is no longer what
	bounds this at the default precision -- the working width is, and
	e^50 and e^200 sit at 50 or 51 digits where they sat at 48 and 47.
	Widening the guard past sixteen used to chase five audit cases --
	exp(50), sinh and cosh at plus and minus 50 -- at the price of
	half as much again on every sinh and cosh. Sixteen stayed, and the
	trade was written up here as settled.

	**It is not a trade any more, and the reason is worth keeping.**
	Those five were short because the guard is a fixed number of
	places while the answers are large, so their *relative* accuracy
	was capped wherever the guard sat; widening it lifted them by
	brute force and charged everyone for it. lz_raise_e now widens by
	the answer's own magnitude instead, which lifts the same five and
	costs sinh(1) nothing. Re-measured after that change, 32 buys
	**nothing at all** -- identical audit, identical sweep. Sixteen
	stays for the opposite reason it used to.

	The constant does still bound a caller who raises the precision a
	long way -- at 200 places the answer stops near 121 -- and that
	remains a constant-width limit of the same kind kTwoPi used to
	place on the sine before it was widened too.
*/
const long kRaiseEGuardDigits	= 16;

/*
	The width lz_raise_e runs its throwaway magnitude probe at.

	Only the size of e^n is wanted, not its digits, and a dozen places
	is ample for that at any exponent the engine will finish. Running
	the real chain twice instead cost 2.5x at e^50 and 3.4x at e^100.
*/

const long kRaiseEProbePlaces	= 12;

/*
	And the exponent below which the probe is not run at all. e^n has
	magnitude n*log(e) in the number's own base, and base two reaches
	kRaiseEGuardDigits soonest, at n = 11. Below that the probe could
	only ever report that no widening is needed, and exp(1) paid 1.8
	times its cost to be told so.
*/

const long kRaiseEProbeFrom	= 11;

/*
	The guard digits lz_exp_e adds on top of that.

	Only four, because lz_raise_e is already doing the heavy work for
	the integer part. What is left for this one is the Taylor series
	over the fractional part and the multiplication that puts the two
	halves together -- a fractional part computed narrow caps the
	product however good the integer part was.

	Measured, holding lz_raise_e at sixteen: no guard here scores 1355
	of the sweep's cases correctly rounded, four scores 1374, and
	eight and sixteen score 1374 as well. Four is where it stops.

	It matters that it stops there rather than matching the sixteen
	below it. The two nest -- lz_exp_e raises, then lz_raise_e raises
	again on top -- so sixteen here would run the squaring chain at
	sixty-eight places for no gain at all, and sinh, which calls the
	exponential twice, measured 0.12 ms that way against 0.089 at
	four.
*/
const long kExpGuardDigits		= 4;

/*
	How lz_exp_e sizes the halving of its fractional part.

	The reciprocal factorial table stops at f^32/32!, so what it can
	give depends on how large f is -- 70 digits at 0.1, 47 at 0.5, and
	only 38 as f approaches one, which is the worst case that matters.
	Each halving divides the dropped term by 2^33 and buys exactly ten
	digits: 38, 48, 58, 68, 78, measured for f = 0.94.

	Two reaches are kept rather than one, because sizing every call
	from the worst case made exp(0.5) pay for a halving it did not
	need -- at or under a half the table already gives 47.

	So the halvings needed are (width wanted - 38) / 10, and the answer
	is squared back the same number of times. kExpMaxHalvings only
	stops a runaway; at a thousand places the count is under a hundred.

	**Halving is why exp reaches so much further on the cheap path
	than lzi_trig_series can.** The sine's argument is already folded
	as far as it goes; exp can fold its own further whenever it likes,
	and one division beats forty. It still needs the recurrence
	eventually -- see kExpFractionThreshold -- because halving cannot
	widen the entries, only shorten the argument.
*/

const long kExpTableReach		= 38;
const long kExpTableReachHalf	= 47;
const long kExpHalvingGain	= 10;
const long kExpMaxHalvings	= 32;

/*
	And the width past which no halving helps, so the terms are
	computed instead.

	Halving fixes the table's *length*; nothing fixes its entries'
	**width**, which is 50 significant digits. Seven halvings at 100
	places restored the length and the answer still stopped at 55 --
	the same wall lzi_trig_series met from the other side. Forty-eight
	leaves room under the fifty and keeps the default, and every
	widening the magnitude asks for at the default, on the cheap path.
*/

const long kExpFractionThreshold	= 48;

/*
	**Past a hundred digits the computed terms are halved first too.**
	The recurrence costs a multiply, a divide and an add a term, each as
	long as the width, and at 21,750 places e^0.5 took some 6,000 terms
	and 3.4 s. Sixty-four halvings leave about a thousand, for
	sixty-four squarings; each squaring doubles the relative error, so
	the series and the squarings work 64 log(2) places wider -- twenty
	digits in base ten. An odd base never halves: there the halved
	argument cannot terminate, fills the width, and makes every term a
	full-width multiply -- see lz_exp_e.

	**The hundred is measured.** It was a thousand until 2026-09-26, put
	there so nothing narrower moved, and exp at 300 places took 1.65 ms
	where it takes 0.50 halved. Halving first starts to win between 64
	and 84 digits; under that the sixty-four squarings cost more than
	the terms they save, 10 to 30% at 45 to 60.
*/

const long kExpComputedHalvedFrom	= 100;
const long kExpComputedHalvings		= 64;

/*
	A long fraction is split wherever the series is computed, e^f =
	e^a e^b, a its first places and b the rest, so the head's terms are
	cheap and the tail's series is short -- see lz_exp_e. The head takes
	a thirty-second of the width, and never fewer than 64 places, so a
	fraction written by hand is never split.
*/

const long kExpSplitShare			= 32;
const long kExpSplitLeast			= 64;

/*
	The guard digits the three means of a list work to.

	Each rounds only a few times -- the arithmetic mean once, at its
	division; the geometric mean at each mantissa product and its root;
	the harmonic mean at each reciprocal and its quotient -- so a small
	guard is enough. What the geometric mean needs beyond it is not a
	guard at all but width for the answer's magnitude, and it takes that
	separately; see lz_gmean.
*/
const long kMeanGuardDigits = 4;

/*
	The guard digits lz_agm works to. Each cycle rounds a sum, a halving,
	a product and a square root, and three is the least that rounds all
	of 800 random pairs correctly; see lz_agm for the measurement. The
	width it needs for its arguments' magnitudes is separate, as for
	lz_gmean.
*/
const long kAgmGuardDigits = 4;

/*
	The guard digits lz_root_3's Newton iteration works to. It had none
	until 2026-09-30, and needed none while each division carried its
	dividend's places on top of the width: those were dozens of guard
	digits nobody had asked for. With a quotient's width set by its value
	they went, and the iteration settled a unit off the correctly rounded
	root -- 16 of the cube roots in the census, in bases two, ten and
	sixteen. They had never been enough: before any of this the cube
	root missed correct rounding in 393 of 9,020 random arguments, a
	unit off in every base.

	Measured on 9,020 random cube roots in bases two to sixteen at 36
	and 80 places, each against the exact root: guards of nought to
	three digits missed 203, 66, 8 and 1 of the first 3,000; four and
	five missed none of any. Five, a digit past the least that never
	missed, as kAgmGuardDigits is.
*/
const long kCubeRootGuardDigits = 5;

/*
	The guard digits lz_root_2's Newton iteration works to, which it had
	none of until 2026-10-01: it worked at the answer's places and
	answered whatever its last multiplication rounded to. Measured as
	the cube root's was, on 9,100 random square roots in bases two to
	sixteen at 36 and 80 places, each against the exact root: with the
	half made a division, guards of nought to three digits missed 166,
	81, 5 and 1 of the first 3,080; four, five and six missed none of
	any. Five, a digit past the least that never missed.
*/
const long kSquareRootGuardDigits = 5;

/*
	The guard digits the five reciprocal-then-inverse functions work
	at: the arccosecant and arcsecant, and the hyperbolic arccosecant,
	arcsecant and arccotangent. Each is its co-function's inverse taken
	on 1/x.

	Four of the five are worst just outside |x| = 1, where the
	reciprocal lands just inside it and the inverse's slope -- 1/sqrt(1
	- z^2), or its hyperbolic counterpart -- climbs without bound, so a
	reciprocal rounded to the working precision before it is handed
	over arrives multiplied by that slope. Measured against the floor
	each answer can hold, going without costs two digits on the
	arcsecant and arccosecant, two on the hyperbolic arcsecant, and six
	on the hyperbolic arccotangent. The fifth, the hyperbolic
	arccosecant, has no slope to fear and gains only one digit at four
	points -- see the note there, and the three readings it took.

	Four, because four is where it stops helping: two is enough for the
	arcsecant and arccosecant, the rest need four, and 6, 8, 12 and 24
	all score identically across the ladder.

	The width has to wrap the inverse call and not only the division.
	Widening the division alone changes nothing, because lz_arcsine
	narrows its own argument to glWorkingPrecision on the way in -- so
	the extra places are handed over and immediately discarded. That
	was measured the wrong way round first, and read as the guard being
	worthless. The answer is rounded back afterwards, so these places
	never leave the function.
*/
const long kArcCoFunctionGuardDigits = 4;

/*
	The guard digits the arctangent works at. It had none until
	2026-09-26, and every step it takes rounds: 1/x for the reflection,
	the square, the sum, the root and the quotient all at the working
	width, so the arcsine was handed an argument up to half a unit out
	and the reflection added another rounding. At the default it missed
	correct rounding everywhere but at tiny arguments, below one as well
	as above: in base ten 8 of 60 below one and 10 of 60 between one and
	two; in bases eight, twelve and sixteen about 40% of the arguments
	from one up, by up to 1.15 units.

	Over 3,500 random arguments in seven bases, misses went 718 at no
	guard, 63 at one, 11 at two, and none at three, four, six or eight.
	**Four, not three**: three still missed in base two at 40 and 48
	places, where four missed nothing up to the arcsine's own ceiling --
	which in hex begins at 44 places, and no guard reaches.
*/
const long kArcTanGuardDigits = 4;

/*
	The guard digits the two complements work at, arccos and arccot:
	each is pi/2 less an inverse. Neither had any until 2026-09-26, and
	both said they needed none -- the arccosine because guards had been
	"measured at not one digit anywhere". Digits cannot see this. The
	inverse rounds its answer to the working width and the subtraction
	rounds it again, so the error reaches nearly a unit, and outside
	base ten pi/2 is itself a rounded conversion.

	At the default, over random arguments in seven bases, arccos missed
	correct rounding 334 times in 1,400 -- 12% below one in base ten,
	38% in hex -- and arccot 820 in 3,500, by up to 1.63 units. Making
	the arctangent correct took arccot's worst under a unit and its
	count up, to 1,001: its errors no longer cancelled arccot's own.

	Misses went 1,335 at no guard, 62 at one, 5 at two, none at three or
	four. **Four, not three**, as the arctangent's: three still missed
	once in base two at 40 places.
*/
const long kArcComplementGuardDigits = 4;

/*
	The guard digits the two-argument arctangent works at, lz_arctangent2,
	since 2026-10-22. Up to four steps round at the guarded width -- the
	quotient of the smaller magnitude by the larger, the arctangent's own
	answer, and pi/2 and pi where the angle is reflected -- and then the
	answer rounds once. Measured before that last rounding, over 5,800
	random pairs in every base at 20 to 80 places, the worst error was
	5e-4 of a unit at a guard of three, 5e-5 at four and 5e-6 at five.
	**Four**, as the arctangent and the complements work at: twenty times
	inside the retake's window, kRetakeDoubt, where three would leave it
	only twice.
*/
const long kArcTan2GuardDigits = 4;

/*
	The guard digits the hyperbolic arccosine works at. It had none until
	2026-09-26. Over 1,050 arguments near one and 1,050 across its
	ordinary range, misses went 260 and 0 at none, 31 and 24 at one, 5
	and 1 at two, 0 and 1 at three, and none at four. **One is worse than
	none away from one**: without a guard the logarithm's rounding is the
	only one, and with one the answer is rounded twice -- by the
	logarithm at the wider width, then on the way out. Four makes that
	negligible; three still missed once, by 0.5002, in base two.
*/
const long kArcHyperCosineGuardDigits = 4;

/*
	The guard for a function whose answer is a quotient: the hyperbolic
	functions, which end in a division by two or by one another, their
	reciprocals, and the reciprocal. Worked this much wider and rounded
	once -- see lzi_round_quotient_answer.

	Measured on the rounding audit at 1, 2, 4 and 8. One place clears
	sinh, cosh, tanh and sech of their misses, 6 to 15 each without a
	guard. The cosecant and cotangent go on improving past four, from 8
	misses each to 6 and 4 at eight, and what is left of them is the
	small-divisor floor, which a guard only nibbles. Four is the
	family's.
*/
const long kQuotientGuardDigits = 4;

/*
	Where lz_archypersecant stops summing its own series and goes
	through the hyperbolic arccosine instead.

	The series costs about p * ln(10) / (2 * ln(1/x)) terms to reach p
	places, and measurement puts the crossover at a term count of
	roughly p/4 -- past that the terms cost more than the square root
	they save. Setting those equal:

		p * ln(10) / (2 * ln(1/x))  =  p / 4
		ln(1/x) = 2 * ln(10),  so  x = 1/100

	The precision cancels, which is why one number serves at every
	width, and measurement agrees: timed at 36, 80, 160 and 320
	places, the series is 1.12, 1.21, 1.16 and 0.99 times the
	arccosine route's speed at exactly this argument. Above it the
	series loses badly -- 0.18 at x = .5 and 320 places -- and below
	it wins by 1.1 to 1.2.

	Held as text rather than a double so it converts exactly, the way
	every other constant here does.
*/
const char kArcHyperSecantSeriesCeiling[] = ".01";


/////////////////////////////////////////////////////////////////
//  Utility functions

long	lzu_setup			();
long	lzu_takedown		();

// The most places an answer's order of magnitude may reach beyond the
//	working precision: the whole places of a large answer, or the noughts
//	after the point of a small one that keeps its significant digits, as
//	a quotient does. An answer past it is refused before any of it is
//	worked out: the result zeroed, kTooLarge returned, and kfInfinity set
//	if it was too large. One too small that would round to nought anyway
//	is nought, promptly. The default is the most digits number entry
//	accepts; a host may raise or lower it, for the whole process, before
//	a second thread starts. Answers the previous limit, and a limit of
//	nought or less puts back the default...
constexpr long	kAnswerLimitPlaces = 100000;
long	lzu_set_answer_limit	( long lPlaces );


/////////////////////////////////////////////////////////////////
//	Internal functions

long	lzi_add_op			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber, long resultSign );
long	lzi_divide_op		( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber, long precision );
long	lzi_subtract_op		( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber, long resultSign );
long	lzi_trig_series		( CLZNumber* argument, bool useCosineSeries, CLZNumber* targetNumber );

// Sets the quotient length past which division takes the Newton reciprocal,
//	answering the previous one. For the tests: the floor is out of reach
//	otherwise, and a path nothing reaches is a path nothing tests...
long	lzi_set_reciprocal_floor( long lPlaces );

// Whether a long divisor takes the reciprocal from its base's width,
//	answering the previous setting. For the tests, to compare the two
//	paths on the shapes the reciprocal takes...
bool	lzi_set_reciprocal_by_shape( bool bByShape );

// Whether division cuts a dividend wider than its answer needs, answering
//	the previous setting. For the tests, to compare a cut quotient with
//	the whole dividend's...
bool	lzi_set_divide_cuts		( bool bCuts );

// The width, in decimal digits, past which the logarithm takes Newton on
//	the exponential rather than its series, answering the previous one.
//	For the tests, to compare the two...
long	lzi_set_log_newton_from( long lDigits );

// For the tests: the most places past its stage's width of the iterate,
//	or past the stage and two of the operand, in any stage of the
//	logarithm by Newton but the last since the last read, this thread's;
//	the read sets it back. Nought or less, and LONG_MIN for none...
long	lzi_log_newton_overrun	();

// The exponent, in decimal digits either side of one, from which the
//	logarithm moves the radix rather than dividing by e, answering the
//	previous one; never below three. For the tests, to compare the two...
long	lzi_set_log_scale_from( long lDigits );

// Long division on packed limbs: the quotient's digits, least significant
//	first, and in *pbRemainder whether anything was left over -- and, when
//	plHalf is given, where twice the remainder stands against the divisor,
//	-1, 0 or 1, which an odd base needs to round -- lHalfPlus, nought or
//	one, added to twice the remainder first, for a cut dividend. Answers
//	false for a base outside two to sixteen, a divisor of nought, or if it
//	cannot allocate...
bool	lzi_packed_divide	( const long* lpDividend, long lDividendLength,
							  const long* lpDivisor, long lDivisorLength, long lBase,
							  long* lpQuotient, long lQuotientLength, bool* pbRemainder,
							  long* plHalf = nullptr, long lHalfPlus = 0 );

// For the tests, each answering the previous setting. The floor is the
//	product of quotient and divisor digits from which division packs, and
//	kPackedFloorMeasured puts back the measured one for each base. The
//	factor is digits per limb, and can only be lowered: at the table's
//	limbs the add-back never fires in anything a test can afford...
constexpr long	kPackedFloorMeasured = -1;

long	lzi_set_packed_divide_floor	( long lProduct );
long	lzi_set_packed_divide_factor( long lFactor );
void	lzi_packed_divide_counts	( long* plDivisions, long* plAddBacks );

// The plain multiplication on packed limbs: writes exactly lLeftLength +
//	lRightLength digits of the product, least significant first. Answers
//	false for a base it has no table for, an empty operand, or if it
//	cannot allocate. lz_multiply's middle tier...
bool	lzi_packed_multiply	( const long* lpLeft, long lLeftLength,
							  const long* lpRight, long lRightLength,
							  long lBase, long* lpProduct );

// Which of lz_multiply's four tiers a product of these lengths takes, and
//	how many products each has made. The dispatch asks the same function,
//	so a test pinning a floor pins the engine's own...
constexpr long	kMultiplyByDigits = 0;
constexpr long	kMultiplyPacked = 1;
constexpr long	kMultiplyByTransform = 2;
constexpr long	kMultiplyKaratsuba = 3;

long	lzi_multiply_tier_for	( long lBase, long lLeftLength, long lRightLength );
void	lzi_multiply_tier_counts( long* plDigits, long* plPacked, long* plTransform, long* plKaratsuba = nullptr );

// For the tests: the answers taken again wider since the process began,
//	because the guard left their rounding in doubt -- see CLZAnswerBase...
long	lzi_retake_count		();

// For the tests: the natural logarithms this thread has taken since it
//	began, counting each lz_log_e however it is reached. A whole power, a
//	whole exponent of the base and a root of whole reciprocal order take
//	none -- they go to the squaring...
long	lzi_logarithm_count		();

// For the tests: the calls of lz_exp_e this thread has made since it began,
//	a retake's included. sinh, cosh, tanh and coth take one each, the other
//	power of e being its reciprocal...
long	lzi_exponential_count	();

// For the tests: the quotient digits this thread's long divisions have
//	worked since it began, packed or a digit at a time. A dividend wider
//	than its answer is cut first, so it works the answer's places and not
//	its own...
long	lzi_quotient_digits_worked	();

// For the tests: the digits of the widest divisor the last reciprocal this
//	thread worked out multiplied by. Each pass takes the divisor to its own
//	width and a guard, so a long divisor over a short quotient is not
//	multiplied whole...
long	lzi_reciprocal_widest_divisor	();

// For the tests: the short answers this thread has checked by working
//	their powers out exactly -- an answer whose places rule it out is not
//	worked, see lzi_places_of_power in LZMath.cpp...
long	lzi_exactness_proofs	();

// e to a fraction, minus one to one, by binary splitting a chunk at a time
//	-- the bit-burst method -- at the working precision, under a guard of
//	its own; lz_exp_e takes the same kernel past lzi_exp_split_floor, under
//	lz_exp_e's guard instead. lFirstChunk sets the first chunk's places,
//	nought the default, for the tests and for measuring...
long	lzi_exp_by_splitting	( CLZNumber* pFraction, CLZNumber* pTarget, long lFirstChunk = 0 );

// For the tests and for measuring, answering the previous setting: the width
//	in places from which lz_exp_e takes the kernel above, for every base,
//	and kExpSplitFloorMeasured puts back each base's measured floor; the
//	floor in force for a base; and the count of the kernel's runs, this
//	thread's...
// ln x by the arithmetic-geometric mean, x positive, at the working
//	precision; lz_log_e takes it by lzi_log_agm_floor's rule. The cycle
//	count is this thread's, for the tests and for measuring...
long	lzi_log_by_agm			( CLZNumber* pValue, CLZNumber* pTarget );
long	lzi_agm_log_cycles		();

// For the tests and for measuring, answering the previous setting: the width
//	in places from which lz_log_e takes the AGM, for every base, and
//	kLogAgmFloorMeasured puts back each base's measured floor; and the
//	floor in force for a base...
constexpr long	kLogAgmFloorMeasured = -1;

long	lzi_set_log_agm_floor	( long lPlaces );
long	lzi_log_agm_floor		( long lBase );

// For the tests, answering the previous setting: the AGM logarithm's guard
//	in decimal digits, so that the margins it would otherwise cover can be
//	seen to work; a negative one puts back the guard it has...
long	lzi_set_log_agm_guard	( long lDigits );

// For the tests and for measuring, answering the previous setting: the band's
//	width from which lz_root_2 takes the AGM's core rather than Heron's
//	iteration, for every base, and kRootCoreFloorMeasured puts back each
//	base's measured floor; the floor in force for a base; and the count of
//	the core's runs, this thread's...
constexpr long	kRootCoreFloorMeasured = -1;

long	lzi_set_root_core_floor	( long lPlaces );
long	lzi_root_core_floor		( long lBase );
long	lzi_root_core_runs		();

// For the tests: the most places past its working width of any root the
//	AGM's square root has answered since the last read, which sets it back,
//	this thread's -- the core of lz_root_2 past its floors and of the
//	logarithm by the AGM; nought or less, each being rounded, and LONG_MIN
//	when none was taken...
long	lzi_root_core_overrun	();

// For the tests: the cube root's Newton passes this thread has made, every
//	width's together, so a stage running to its cap can be seen...
long	lzi_root_3_passes		();

// For the tests and for measuring, answering the previous setting: the width in
//	places from which the cube root takes its stages after the first by its
//	reciprocal, every base, and kRoot3InverseFloorMeasured puts back each
//	base's floor; the floor in force for a base; and the reciprocal's passes,
//	this thread's...
constexpr long	kRoot3InverseFloorMeasured = -1;

long	lzi_set_root_3_inverse_floor	( long lPlaces );
long	lzi_root_3_inverse_floor		( long lBase );
long	lzi_root_3_inverse_passes		();

// For the tests, answering the previous setting: the reciprocal's guard in
//	places, so the margin it covers can be seen to matter; a negative one
//	puts back the four it has...
long	lzi_set_root_3_inverse_guard	( long lPlaces );

constexpr long	kExpSplitFloorMeasured = -1;

long	lzi_set_exp_split_floor	( long lPlaces );
long	lzi_exp_split_floor		( long lBase );
long	lzi_exp_splittings		();

// Sine and cosine of y together by binary splitting e^{iy} a chunk at a
//	time -- the bit-burst method -- at the working precision in places,
//	under a guard of its own; either target may be null. lz_sine and
//	lz_cosine take the same kernel past lzi_sincos_split_floor, under their
//	own guard instead...
long	lzi_sincos_by_splitting	( CLZNumber* pArgument, CLZNumber* pSine, CLZNumber* pCosine );

// For the tests and for measuring, answering the previous setting: the width in
//	places from which sine and cosine take the kernel above, for every base,
//	and kSinCosSplitFloorMeasured puts back each base's floor; the floor in
//	force for a base; and the count of the kernel's runs, this thread's...
constexpr long	kSinCosSplitFloorMeasured = -1;

long	lzi_set_sincos_split_floor	( long lPlaces );
long	lzi_sincos_split_floor		( long lBase );
long	lzi_sincos_splittings		();

// For the tests: the answers lz_sine_cosine and lz_hypersine_hypercosine
//	have handed to their own functions, this thread's -- one in doubt, exact
//	or flagged. A pair that handed every one on would answer the same...
long	lzi_pair_handoffs			();

// For the tests: the chunks that were not nought in the kernel's last run,
//	this thread's -- each is taken to the nearest, so a remainder of either
//	sign costs what its size does...
long	lzi_sincos_last_chunks		();

// For the tests and for measuring, answering the previous setting: the width in
//	places from which lz_arcsine -- and every inverse circular function,
//	which all come down to it -- takes the angle by Newton on the kernel
//	above, for every base, and kArcSineNewtonFloorMeasured puts back each
//	base's floor; the floor in force for a base; and the count of Newton's
//	runs, this thread's...
constexpr long	kArcSineNewtonFloorMeasured = -1;

long	lzi_set_arcsine_newton_floor	( long lPlaces );
long	lzi_arcsine_newton_floor		( long lBase );
long	lzi_arcsine_newtons				();

// For the tests: the most places past its width of Newton's seed as it
//	is converted, the larger of its text's places and the first stage; or
//	of the iterate or s and c in any stage but the last, since the last
//	read, this thread's; the read sets it back. Nought or less, LONG_MIN
//	for none...
long	lzi_angle_newton_overrun		();

// For the tests: the most places past the first cycle's width of the two
//	means lz_agm takes into it, since the last read, this thread's; the read
//	sets it back. Nought or less, LONG_MIN for none...
long	lzi_agm_first_overrun			();

// For the tests: the cycles this thread's calls of lz_agm have taken...
long	lzi_agm_cycles					();

// For the tests, each answering the previous setting; kPackedFloorMeasured
//	puts back the measured floor. Packing's is a product of lengths, the
//	transform's a length of the shorter operand -- which nothing the tests
//	can afford reaches at its measured value...
long	lzi_set_packed_multiply_floor	( long lProduct );
long	lzi_set_transform_floor			( long lDigits );

// For the tests, answering the previous setting: the longest transform
//	attempted, in elements. It can only be lowered, and 0 puts back 2^22.
//	Past it the transform declines and Karatsuba takes the product, which
//	nothing the tests can afford would otherwise show...
long	lzi_set_transform_reach			( long lLength );

// For the tests, each answering the previous setting. The factor is
//	digits per limb and the fold is products per column between folds;
//	either can only be lowered, and 0 puts back the worked-out value...
long	lzi_set_packed_multiply_factor	( long lFactor );
long	lzi_set_packed_multiply_fold	( long lProducts );
void	lzi_packed_multiply_counts		( long* plMultiplications, long* plFolds );

// Karatsuba on packed limbs, lz_multiply's third tier: the same contract as
//	lzi_packed_multiply. The floor is in limbs, never below four, and 0 puts
//	back the measured 96; it is the dispatch's threshold as well as the
//	kernel's, so the tests can force the tier, or switch it off with a huge
//	floor. The split count is for the tests, to know the recursion ran...
// The product of two digit arrays by lz_multiply's own dispatch, every tier
//	included: lpProduct holds lLeftLength + lRightLength digits, zeroed, and
//	neither length is nought. CLZNumber's base conversion multiplies this
//	way past its floor...
void	lzi_multiply_digit_arrays	( const long* lpLeft, long lLeftLength,
									  const long* lpRight, long lRightLength,
									  long* lpProduct, long lBase );

bool	lzi_karatsuba_multiply	( const long* lpLeft, long lLeftLength,
								  const long* lpRight, long lRightLength,
								  long lBase, long* lpProduct );
long	lzi_set_karatsuba_floor	( long lLimbs );
long	lzi_karatsuba_splits	();


/////////////////////////////////////////////////////////////////
//  The angle unit -- stage 2b of the complex layer, 2026-10-22

/*
	The unit the angle functions read an angle in or write one in: the
	six circular functions read it, their six inverses, lz_arctangent2
	and lzc_arg write it. **One per thread, set directly, as the working
	precision is** -- Rich's call, 2026-10-22. Radians unless a host says
	otherwise, and a new thread starts at radians whatever another has set.

	A degree is a 360th of a turn, a grad a 400th and a millicle a
	thousandth, Elzed's four units. In each but the radian a quarter and
	an eighth of a turn are whole, so an angle is reduced by whole turns
	exactly, before pi is used at all, and an angle that is a whole
	number of twelfths or eighths of a turn answers exactly where the
	mathematics does: sin(30) in degrees is 0.5, arctan(1) is 45, tan(90)
	a pole, answered as cot(0) is in radians.

	**It applies to the host's own call only.** The engine's own uses of
	these functions -- a tangent's sine and cosine, an arccosine's
	arcsine, gamma's reflection -- stay in radians whatever is set here.

	A value that names none of the four reads as radians. **The enum has a
	fixed underlying type so that such a value is defined at all**: a host
	casting a setting of its own into it, as the suite does with 7, would
	otherwise hold a value outside the enumeration's range, which UBSan
	reported the first time the sanitized suite ran.
*/
enum AngleUnit : int
{
	AU_radians = 0,
	AU_degrees,
	AU_grads,
	AU_millicles
};

extern thread_local AngleUnit	geAngleUnit;

/*
	The guard digits an angle function works at in a unit other than the
	radian, past whatever its radian twin takes for itself: the angle's
	conversion with pi, and the turn over two pi on the way back, which
	takes three more. Measured before the final rounding, over 4,000
	random arguments to all twelve in every base at 20 to 80 places, the
	worst error is half a unit of the guard's last place -- 5e-5 at four,
	5e-6 at five, 5e-7 at six -- for the forward functions, and a sixth of
	that for the inverses. **Six**: twenty times inside the cosecant's
	retake window, 1e-5, the narrowest, where four would sit outside it.
*/
const long kAngleUnitGuardDigits = 6;

/*
	The guard digits lz_sine_cosine and lzc_polar work at, in radians and
	in a unit alike -- stage 3 of the complex layer, 2026-10-22. Measured
	before the final rounding, over 1,500 arguments in every base at 20,
	36 and 80 places, against lz_sine and lz_cosine thirty places wider:
	the worst error in units of the answer's last place --

	  guard		radians		in a unit
	  3			4.2e-3		4.0e-3
	  4			2.4e-4		2.8e-4
	  5			3.4e-5		2.5e-5
	  6			2.6e-6		2.6e-6

	**Five**: thirty times inside the cosine's doubt window, kRetakeDoubt,
	where four sits four times inside it. An answer in doubt goes to its
	own function, which retakes it, so the guard decides only how often.
*/
const long kSineCosineGuardDigits = 5;

/*
	The guard digits lzc_log_e, lzc_exp_e and lzc_root_2 work their parts
	at -- stage 4 of the complex layer, 2026-10-22. Measured before the
	final rounding, over 800 operands for each in every base at 20, 36 and
	80 places, against the same function thirty places wider: the worst
	error in units of the answer's last place --

	  guard		log			exp			root
	  3			2.4e-4		2.7e-3		6.4e-4
	  4			2.5e-5		2.1e-4		5.1e-5
	  5			2.5e-6		1.4e-5		5.0e-6
	  6			2.6e-7		1.5e-6		6.0e-7

	**Five**: the exponential, the worst, seventy times inside
	kRetakeDoubt, where four sits five times inside it.
*/
const long kComplexGuardDigits = 5;


/////////////////////////////////////////////////////////////////
//  Math API functions

long	lz_abs				( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_add				( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_agm				( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

/*
	The means of a list -- the only functions here that take one. The
	target is required and is refused if null; the values are never
	modified, so the target may be one of them. An infinity or a
	not-a-number among the values makes the answer a not-a-number,
	since each of the three needs a sign the engine may not hold. The
	geometric and harmonic means refuse a negative value (kFailed) and
	answer nought for a nought.
*/
long	lz_amean			( CLZNumber** values, long count, CLZNumber* targetNumber );
long	lz_arccosecant			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_arccosine			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_arccotangent			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_archypercosecant		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_archypercosine		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_archypercotangent	( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_archypersecant		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_archypersine			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_archypertangent		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_arcsecant			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_arcsine				( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_arctangent			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

/*
	The angle of the point (x, y), atan2: the LEFT operand is y and the
	right x, as C's atan2( y, x ) takes them. In (-pi, pi], in radians for
	now, correctly rounded. On the axes it answers as std::atan2 does for
	+0: 0 along the positive real axis, pi along the negative one, and
	+-pi/2 on the imaginary axis; the origin answers 0. An infinity or a
	not-a-number in either answers a not-a-number. A missing operand is
	refused, kFailed. A right operand in another base is brought into the
	left's exactly.
*/
long	lz_arctangent2			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_combinations		( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_compare			( CLZNumber* leftOperand, CLZNumber* rightOperand );
long	lz_cosecant			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_cosine			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_cotangent		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_cube				( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

long	lz_divide			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_divide			( CLZNumber* leftOperand, CLZNumber* rightOperand, long precision );
long	lz_divide_integer	( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_exp_b			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_exp_e			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_exp_n			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_factorial		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_gamma			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_gmean			( CLZNumber** values, long count, CLZNumber* targetNumber );
long	lz_log_gamma		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_fraction			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

long	lz_hmean			( CLZNumber** values, long count, CLZNumber* targetNumber );
long	lz_hypercosecant	( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_hypercosine	    ( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_hypercotangent	( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_hypersecant		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_hypersine	    ( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

/*
	sinh and cosh of one argument, from one exponential: each answer
	exactly what lz_hypersine and lz_hypercosine answer, flags and codes
	included, for about the cost of one of them. Both targets are
	required and must be different numbers -- either may be the operand
	-- or the call is refused, kFailed, each target given a not-a-number.
	Stage 3 of the complex layer, 2026-10-22.
*/
long	lz_hypersine_hypercosine	( CLZNumber* leftOperand, CLZNumber* sineTarget, CLZNumber* cosineTarget );
long	lz_hypertangent	    ( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

long	lz_integer			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

long	lz_log_b			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_log_e			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_log_n			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_max				( CLZNumber* leftOperand, CLZNumber* rightOperand );
long	lz_min				( CLZNumber* leftOperand, CLZNumber* rightOperand );
long	lz_modulo			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_multiply			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_negate			(CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_normalize		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

long	lz_permutations		( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_plus				( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_power			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_quantize			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_raise_e			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_raise_n			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_reciprocal		( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_remainder		( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_remainder_near	( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_root_n			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );
long	lz_root_2			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_root_3			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_round			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_round_integral	( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_round_last_digit	( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

long	lz_same_quantum		( CLZNumber* leftOperand, CLZNumber* rightOperand );
long	lz_secant			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_sine				( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );

/*
	Sine and cosine of one angle, from one reduction and one run: each
	answer exactly what lz_sine and lz_cosine answer, flags and codes
	included, in the host's angle unit as they are, for about half the
	cost of the two -- 1.3 to 3.4 times less, measured. Both targets are required and must be different
	numbers -- either may be the operand -- or the call is refused,
	kFailed, each target given a not-a-number. Stage 3 of the complex
	layer, 2026-10-22.
*/
long	lz_sine_cosine		( CLZNumber* leftOperand, CLZNumber* sineTarget, CLZNumber* cosineTarget );
long	lz_square			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );
long	lz_subtract			( CLZNumber* leftOperand, CLZNumber* rightOperand, CLZNumber* targetNumber = nullptr );

long	lz_tangent			( CLZNumber* leftOperand, CLZNumber* targetNumber = nullptr );


/////////////////////////////////////////////////////////////////
//  Complex numbers -- stage 1, 2026-10-21.

/*
	A complex value: its real part and its imaginary part, in one base,
	which is the value's base. Each part of an answer is rounded to the
	working precision in places, exactly as the real operation that
	writes it rounds: the four plain operations in that part's own
	rounding mode and digit -- so a target's two parts may be set to
	round differently -- and lzc_square and lzc_abs as functions do, the
	same whatever the target's settings. Every answer is
	in its left operand's base; a right operand in another base is
	brought into it exactly, as a real pair is. A value whose two parts
	are in different bases is refused, kBaseMismatch, with both parts of
	the answer a not-a-number.

	An infinity in either part, and no not-a-number, is the one point at
	infinity -- this engine's infinities carry no sign, so a complex one
	has no direction -- and comes back with both parts infinite. A
	not-a-number in either part makes both parts one.

	Stage 1 is the operations correctly rounded from exact intermediates.
	The target is optional, as for the real functions, but for lzc_abs,
	whose answer is real: it is required and refused, kFailed, if null.

	Stage 2a, 2026-10-22, adds lzc_arg, the argument: lz_arctangent2 of
	the imaginary part over the real one, in (-pi, pi] and in radians for
	now, 0 for nought. Its answer is real, so its target is required as
	lzc_abs's is; a not-a-number or the point at infinity has no angle,
	and answers a not-a-number.
*/
struct CLZComplex
{
	CLZNumber	m_real;
	CLZNumber	m_imag;
};

long	lzc_abs				( CLZComplex* leftOperand, CLZNumber* targetNumber );
long	lzc_add				( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );
long	lzc_arg				( CLZComplex* leftOperand, CLZNumber* targetNumber );
long	lzc_conjugate		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_divide			( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );
long	lzc_multiply		( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );
long	lzc_negate			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );

/*
	A complex value from its modulus and its angle, r cos theta + i r sin
	theta, the angle read in the host's unit -- stage 3, 2026-10-22. Each
	part correctly rounded, half-even, at the working precision in places,
	and exact where the mathematics is: at whole quarter turns the parts
	are r, -r and an exact nought, so 1 |/ 180 degrees is -1 + 0i. A
	negative modulus is taken as it is. The answer is in the modulus's
	base, an angle in another brought into it. The target is required,
	kFailed if null; a missing operand is refused the same way. A
	not-a-number in either, or an infinite angle, answers a not-a-number,
	and an infinite modulus at a finite angle the point at infinity.
*/
long	lzc_polar			( CLZNumber* modulus, CLZNumber* angle, CLZComplex* targetNumber );

/*
	The principal natural logarithm, e to a complex power, and the
	principal square root -- stage 4, 2026-10-22. Each part correctly
	rounded, half-even, at the working precision in places, and exact
	where the mathematics is; the logarithm's imaginary part and the
	exponential's angle in radians whatever the angle unit. The cut is
	the negative real axis, taken from the +0 side: ln(-1) is pi i and
	sqrt(-4) is 0 + 2i, exactly. On the real axis each answers as the
	real function does, with an exact nought beside it. ln 0 is refused,
	kUndefinedResult, and an exponential past the answer limit, kTooLarge,
	each answering the point at infinity. The target is optional.
*/
long	lzc_exp_e			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_log_e			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_root_2			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );

/*
	The twelve circular and hyperbolic functions -- stage 5a, 2026-10-22.
	Each part correctly rounded, half-even, at the working precision in
	places. The circular six read the whole complex angle in the host's
	unit, so z degrees is z pi / 180, both parts; the hyperbolic six are in
	radians whatever the unit. On the real axis each answers as the real
	function does, with an exact nought beside it, and a pole there is the
	point at infinity; in a unit, a part whose sine or cosine is exactly
	nought is an exact nought, so cos(90 + i) degrees is 0 - 0.01745...i.
	An answer past the answer limit is refused, kTooLarge, the point at
	infinity. The target is optional.
*/
long	lzc_sine			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_cosine			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_tangent			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_cotangent		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_secant			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_cosecant		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_hypersine		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_hypercosine		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_hypertangent	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_hypercotangent	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_hypersecant		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_hypercosecant	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );

/*
	The power, the roots and the base's exponential -- stage 5b,
	2026-10-22. lzc_power is z^w, the principal value e^(w ln z), w
	complex; lzc_raise_n z^n for a whole n, a fraction refused; lzc_root_n
	the principal root of complex order w, z^(1/w), order nought refused;
	lzc_root_3 the cube root; lzc_exp_b b^z, b the base z is written in.
	Each part correctly rounded, half-even, at the working precision in
	places, and exact where the mathematics is: (2 + i)^3 is 2 + 11i, its
	cube root 2 + i, the square root of -4 exactly 2i, and (-8)^(1/3) is
	1 + 1.732...i, its real part exactly one, in a base that holds a
	third. On the real axis each answers as its real function does, to
	the character, an exact nought beside -- so an odd root of a negative
	real is the real root, the cube root of -8 being -2 -- and a negative
	base to any other real power is the principal value. Nought to w is
	nought where Re w is past nought, the point at infinity where it is
	below, and refused where it is nought; anything to nought is one. An
	answer past the answer limit is refused, kTooLarge, the point at
	infinity. The target is optional.
*/
long	lzc_power			( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );
long	lzc_raise_n			( CLZComplex* leftOperand, CLZNumber* rightOperand, CLZComplex* targetNumber = nullptr );
long	lzc_root_n			( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );
long	lzc_root_3			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_exp_b			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );

/*
	The logarithms to another base -- stage 5b, 2026-10-22. lzc_log_n is
	the logarithm of z to the base n, ln z / ln n, n complex and first, as
	lz_log_n takes them; lzc_log_b the logarithm to the base z is written
	in. The principal value, each part correctly rounded, half-even, at
	the working precision in places, exact where the mathematics is:
	log_i(-1) is exactly 2. On the real axis, base and value past nought,
	each answers as its real function does, to the character, and with a
	real base its real part is a real logarithm's. Nought is refused,
	kUndefinedResult, the point at infinity, as is base nought; base one is
	kDivideByZero. The answer is in the base's base. The target is
	optional.
*/
long	lzc_log_b			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_log_n			( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );

/*
	The inverse functions that take square roots -- stage 6a, 2026-10-23:
	the arcsine, arccosine, arcsecant and arccosecant and their hyperbolic
	twins, the principal values, each part correctly rounded, half-even,
	at the working precision in places. On a cut the operand's +0 side, as
	std::complex answers: arcsin 2 is pi/2 + 1.317...i, arcsec 0.5 is
	1.317...i. On the real axis inside its domain each answers as its real
	function does, to the character, an exact nought beside. The circular
	four write the whole complex angle in the thread's unit -- arcsin 2 in
	degrees is 90 + 75.456...i -- the hyperbolic four in radians. Nought is
	a pole of the four reciprocals, kUndefinedResult, the point at
	infinity. The target is optional.
*/
long	lzc_arcsine			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_arccosine		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_arcsecant		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_arccosecant		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_archypersine	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_archypercosine	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_archypersecant	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_archypercosecant( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );

/*
	The inverse functions that take logarithms -- stage 6b, 2026-10-23:
	the arctangent and arccotangent and their hyperbolic twins, the
	principal values, each part correctly rounded, half-even, at the
	working precision in places. arccot z is arctan(1/z), as lz_arccotangent
	is arctan(1/x). On a cut the operand's +0 side, as std::complex
	answers: arctanh 2 is 0.549... + (pi/2)i, arctan 2i is pi/2 + 0.549...i,
	arccoth 0 is -(pi/2)i. On the real axis inside its domain each answers
	as its real function does, to the character, an exact nought beside.
	The circular two write the whole complex angle in the thread's unit --
	in degrees arctan(0.6 + 0.8i) has a real part of exactly 45 -- the
	hyperbolic two in radians. arctanh and arccoth of +-1 and arctan and
	arccot of +-i are poles, kUndefinedResult, the point at infinity. The
	target is optional.
*/
long	lzc_arctangent		( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_arccotangent	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_archypertangent	( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_archypercotangent( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_square			( CLZComplex* leftOperand, CLZComplex* targetNumber = nullptr );
long	lzc_subtract		( CLZComplex* leftOperand, CLZComplex* rightOperand, CLZComplex* targetNumber = nullptr );
