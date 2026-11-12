/////////////////////////////////////////////////////////////////
//
//	AddapticTests.cpp
//
//	Accuracy and precision tests for the Addaptic engine.
//
//	Each check reports how many leading significant digits the
//	engine agrees with a known-correct reference to. Reference
//	values are correctly rounded to 50 digits (generated with
//	Python's decimal module; trigonometric values from Taylor
//	series verified against sin^2 + cos^2 = 1).
//
//	Copyright (c) 2026 by Krasche & Byrne LLC
//	Released under the MIT License; see LICENSE.
//	SPDX-License-Identifier: MIT
//
/////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <cstdint>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>

#include "LZNumber.h"
#include "LZMath.h"

#include <atomic>
#include <new>

/*
	**Every allocation through new is counted, live**, so a test can ask
	whether a function hands back everything it made -- see
	TestNothingLeftBehind. The replacements are ordinary C++, so the
	count works under every compiler the engine builds with; the digits
	themselves are calloc'd and not counted, but every number made with
	new is. lz_log_e left one behind on every call until 2026-10-15.
*/

static std::atomic<long>	g_liveNews( 0 );

/*
	**GCC can warn -Wmismatched-new-delete on these, and is wrong to.**
	Once it inlines a library allocation and this delete into one function,
	it sees free() given what the library's operator new returned, and does
	not see that the operator new it returned from is the one below, which
	mallocs. It said so once, in TestArcSineByNewton, when an unrelated
	change to this file moved its inlining. CI fails on any warning, so the
	pair is told, here alone and for GCC alone -- clang and MSVC say
	nothing.
*/

#if defined( __GNUC__ ) && !defined( __clang__ )
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

void* operator new( size_t n )
{
	void*	p = ::malloc( (n > 0) ? n : 1 );

	if ( p == nullptr )
		throw std::bad_alloc();

	g_liveNews++;
	return p;
}

void* operator new[]( size_t n )						{ return ::operator new( n ); }
void operator delete( void* p ) noexcept				{ if ( p != nullptr ) { g_liveNews--; ::free( p ); } }
void operator delete[]( void* p ) noexcept				{ ::operator delete( p ); }
void operator delete( void* p, size_t ) noexcept		{ ::operator delete( p ); }
void operator delete[]( void* p, size_t ) noexcept		{ ::operator delete( p ); }

#if defined( __GNUC__ ) && !defined( __clang__ )
	#pragma GCC diagnostic pop
#endif

/////////////////////////////////////////////////////////////////
//	Reference values, correctly rounded to 50 digits
/////////////////////////////////////////////////////////////////

const char kRefSqrt2[]			= "1.41421356237309504880168872420969807856967187537694";
const char kRefSqrt3[]			= "1.73205080756887729352744634150587236694280525381038";
const char kRefLn2[]			= "0.69314718055994530941723212145817656807550013436025";

// ln(10) far wider than the default precision, for the test that
//	kBaseE is wide enough to serve a caller who has raised it. Ten is
//	used *because* it reduces: the reduction divides by e, so this is
//	the one that the constant's width bounds. It stopped at 49 digits
//	however wide the precision was set until kBaseE was widened.
//	Generated and round-tripped through exp to 160 places...
const char kRefLn10Wide[] =
		"2.3025850929940456840179914546843642076011014886287729760333"
		"27900967572609677352480235997205089598298341967784042";

// ln(2) far wider than the default precision, for the test that the
//	logarithm's series depth tracks the working precision. Generated
//	and round-tripped through exp to 200 places before being written
//	down. Two is used because it needs no argument reduction, and so
//	was never capped by kBaseE the way ln(10) was...
const char kRefLn2Wide[] =
		"0.693147180559945309417232121458176568075500134360255254120680"
		"00949339362196969471560586332699641868754200148102057068573368"
		"55202357581305570326707516350759619307275708283714351903070386"
		"2389167347112335011536449795524";
const char kRefLn10[]			= "2.30258509299404568401799145468436420760110148862877";

// Logarithms of small arguments, generated and round-tripped through
//	exp() -- never typed by hand...
const char kRefLnTwentieth[]   = "-2.9957322735539909934352235761425407756766016229890";
const char kRefLnHundredth[]   = "-4.6051701859880913680359829093687284152022029772575";
const char kRefLnThousandth[]  = "-6.9077552789821370520539743640530926228033044658863";
const char kRefLnTenThousandth[] = "-9.2103403719761827360719658187374568304044059545151";
const char kRefLnTenMillionth[] = "-16.118095650958319788125940182790549453207710420401";
const char kRefLnTenth[]       = "-2.3025850929940456840179914546843642076011014886288";
// The exponential and logarithm within a whisker of one, where the
//	old rounding threw the answer away...
// Cube roots, checked by cubing them back...
const char kRefCbrtTwo[]     = "1.2599210498948731647672106072782283505702514647015";
const char kRefCbrtTiny[]    = "0.000000000000046415888336127788924100763509194465765513491250112";
const char kRefCbrtBig[]     = "10000000000000.000000000000000000000000000000000000";

// The exponential at a size where 36 places affords far more than 36
//	significant digits. See "the exponential's own magnitude"...
const char kRefExp50[] =
		"5184705528587072464087.453322933485384827469100583846401904056933806857";
const char kRefExp100[] =
		"26881171418161354484126255515800135873611118.7737419224151916086152802870349095649141588711";
const char kRefExp1394[] =
		"1132570.062762860298253376972614711189783306013524945553";
const char kRefSinh50[] =
		"2592352764293536232043.726661466742692413734453854430802756139316036288";

const char kRefExp605[] =
		"188285187488585157559895080.92490219528064434389478598502508"
		"622005680567804797944";

// The means of a list, generated by Python's decimal at 200 digits and
//	checked against a closed form or an identity where one exists...
const char kRefGmean1234[]	= "2.2133638394006431848175805468806632971078363650598";	// gm(1, 2, 3, 4)
const char kRefGmean123[]	= "1.8171205928321396588912117563272605024282104631412";	// gm(1, 2, 3)
const char kRefHmean49[]	= "5.5384615384615384615384615384615384615384615384615";	// hm(4, 9)
const char kRefHmean123[]	= "1.6363636363636363636363636363636363636363636363636";	// hm(1, 2, 3)
const char kRefAmeanMixed[]	= "54.360800000000000000000000000000000000000000000000";	// am(3.7, 1.1, 250, 0.004, 17)
const char kRefGmeanMixed[]	= "2.3335046256885920877510194827538930201978890736592";	// gm(3.7, 1.1, 250, 0.004, 17)
const char kRefHmeanMixed[]	= "0.019901116549339608642722359624841166999955767115083";	// hm(3.7, 1.1, 250, 0.004, 17)

// The geometric mean at a scale where the width has to be paid for...
const char kRefGmeanBigEven[] =
		"14142135623730950488016887242096980785696.71875376948073176679737990732478462107038850";
const char kRefGmeanBigOdd[] =
		"44721359549995793928183473374625524708812.36719223051448541794490821041851275609798829";

// exp of a fractional argument at a raised precision, which the
//	fixed 32-term series used to cap...
const char kRefExpHalfWide[] =
		"1.6487212707001281468486507878141635716537761007101480115750"
		"7931164066102119421560863277652005636664300286663776";

// The gamma function. Generated by the recurrence from Gamma(1/2) =
//	sqrt(pi), and the negatives by reflection off it -- no series
//	anywhere, so nothing in common with what they are checking...
const char kRefGammaHalf[] =
		"1.77245385090551602729816748334114518279754945612238712821381";
const char kRefGammaOneHalf[] =
		"0.886226925452758013649083741670572591398774728061193564106904";
const char kRefGammaTwoHalf[] =
		"1.32934038817913702047362561250585888709816209209179034616036";
const char kRefGammaTenHalf[] =
		"1133278.38894878556733457416558889247556029830827515977660872";
const char kRefGammaTwentyHalf[] =
		"540624298233507504.473687364780822141227316639578171478960195";
const char kRefGammaMinusHalf[] =
		"-3.544907701811032054596334966682290365595098912244774256427616";
const char kRefGammaMinusOneHalf[] =
		"2.363271801207354703064223311121526910396732608163182837618410";

// Large arguments whose cube root is not itself a whole number. See
//	"the magnitude the band takes out" in TestTranscendentals...
const char kRefCbrt4000[]    = "15.87401051968199474751705639272308260391493327899853010";
const char kRefCbrt1e20[]    = "4641588.833612778892410076350919446576551349125011243638";
const char kRefCbrt1e40[]    = "21544346900318.83721759293566519350495259344942192108582";
const char kRefExpE18[]      = "1.0000000000000000010000000000000000005000000000000";
const char kRefExpE20[]      = "1.0000000000000000000100000000000000000000500000000";
const char kRefExpNegE20[]   = "0.99999999999999999999000000000000000000005000000000";
const char kRefLnNear1[]  = "0.0000000000000000000099999999999999999999500000000000000000003333333333";
const char kRefLnE20[]     = "-46.051701859880913680359829093687284152022029772575";
const char kRefLnE30[]     = "-69.077552789821370520539743640530926228033044658863";
const char kRefLnE40[]     = "-92.103403719761827360719658187374568304044059545151";
const char kRefLnE100[]    = "-230.25850929940456840179914546843642076011014886288";
const char kRefLn12345[]       = "9.4210064017792798779058775355940914883016248742817";
const char kRefE[]				= "2.71828182845904523536028747135266249775724709369995";
const char kRefSin1[]			= "0.84147098480789650665250232163029899962256306079837";
const char kRefCos1[]			= "0.54030230586813971740093660744297660373231042061792";

// The series at a raised precision, where the table used to stop it.
//	See "The series past the table" below...
const char kRefSin1Wide[] =
		"0.8414709848078965066525023216302989996225630607983710656727"
		"5170999191040439123966894863974354305269585434903791";
const char kRefCos1Wide[] =
		"0.5403023058681397174009366074429766037323104206179222276700"
		"9725538110039477447176451795185608718308934357173116";
const char kRefSin4Wide[] =
		"-0.756802495307928251372639094511829094135912887336472571485"
		"41677340131049361917941642357281056242274808769969717";

// The reduction's own references: what kTwoPi's width decides. See
//	the note in TestArgumentMagnitude...
const char kRefSin1e20[]		= "-0.64525128526578084420581171131252300740690419668690";
const char kRefCos1e20[]		= "0.76397040444172830040014680273788112283447344174700";
const char kRefSin1e30[]		= "-0.09011690191213805803038642895298733027439633299304";
const char kRefCos1e30[]		= "-0.99593119440539570239424858799704864113024773495505";
const char kRefSin1e40[]		= "-0.56963340095363632730803418157356872313292131914787";
const char kRefCos1e40[]		= "-0.82189889190702392144480253644325575173454251195073";

const char kRefTan1[]			= "1.55740772465490223050697480745836017308725077238152";
const char kRefSinh1[]			= "1.17520119364380145688238185059560081515571798133409";
const char kRefCosh1[]			= "1.54308063481524377847790562075706168260152911236586";
const char kRefTanh1[]			= "0.76159415595576488811945828260479359041276859725793";
const char kRefAsinHalf[]		= "0.52359877559829887307710723054658381403286156656251";
const char kRefAsinMinusHalf[]	= "-0.52359877559829887307710723054658381403286156656251";
const char kRefAsinTwentieth[]	= "0.050020856805770014662744386820464114977806080494689";
const char kRefAsinNinth[]		= "1.1197695149986341866866770558453996158951621864033";
const char kRefAsinThreeQtrs[]	= "0.84806207898148100805294433899841808007336621326311";

// Arctangent, both sides of the reflection at |x| = 1. Each was
//	cross-checked against arctan(x) + arctan(1/x) = pi/2 before being
//	written down...
const char kRefAtanTenth[]		= "0.099668652491162027378446119878020590243278322504315";
const char kRefAtanHalf[]		= "0.46364760900080611621425623146121440202853705428612";
const char kRefAtanOne[]		= "0.78539816339744830961566084581987572104929234984378";
const char kRefAtanTwo[]		= "1.1071487177940905030170654601785370400700476454014";
const char kRefAtanMinusTwo[]	= "-1.1071487177940905030170654601785370400700476454014";
const char kRefAtanMillion[]	= "1.5707953267948966195646550249728847754319181758780";

// Arccosine. Cross-checked against arccos(x) + arcsin(x) = pi/2 before
//	being written down...
const char kRefAcosMinusOne[]	= "3.1415926535897932384626433832795028841971693993751";
const char kRefAcosMinusHalf[]	= "2.0943951023931954923084289221863352561314462662501";
const char kRefAcosHalf[]		= "1.0471975511965977461542144610931676280657231331250";
const char kRefAcosNineTenths[]	= "0.45102681179626243254464463579435182620342251328425";
const char kRefAcosNearOne[]	= "0.00044721359922673798563403587961532526145234241692909";

// Inverse hyperbolics. Generated, and every one round-tripped
//	through its forward function before being written down...
// Hyperbolics at arguments that are not whole numbers, and past
//	2*pi -- generated and round-tripped through the inverses...
const char kRefSinhHalf[]      = "0.52109530549374736162242562641149155910592898261148";
const char kRefSinhOnePt1[]    = "1.3356474701241767793847805235786784358351007522235";
const char kRefSinhTen[]       = "11013.232874703393377236524554846364402901451190319";
const char kRefSinhMinusHalf[] = "-0.52109530549374736162242562641149155910592898261148";
const char kRefCoshHalf[]      = "1.1276259652063807852262251614026720125478471180987";
const char kRefCoshTen[]       = "11013.232920103323139721376090437879963452061428237";
const char kRefTanhHalf[]      = "0.46211715726000975850231848364367254873028928033011";
const char kRefTanhTen[]       = "0.99999999587769276361959283713827574105081461849502";
const char kRefAsinhHalf[]     = "0.48121182505960344749775891342436842313518433438566";
const char kRefAsinhOne[]      = "0.88137358701954302523260932497979230902816032826164";
const char kRefAsinhTwo[]      = "1.4436354751788103424932767402731052694055530031570";
const char kRefAsinhMinus2[]   = "-1.4436354751788103424932767402731052694055530031570";
const char kRefAsinhFifty[]    = "4.6052701709914238266212392672083064130089022456922";
const char kRefAsinhMinus1k[]  = "-7.6009027095419886115232897846649396335683391204172";
const char kRefAsinhTiny[]     = "0.00099999983333340833328869050657239826277889676200226";
const char kRefAcoshMilli[]    = "0.044717633608309309571769903613797397347973180574293";
const char kRefAcoshTwo[]      = "1.3169578969248167086250463473079684440269819714675";
const char kRefAcoshHundred[]  = "5.2982923656104845907016668349432471689372518831401";
const char kRefAtanhHalf[]     = "0.54930614433405484569762261846126285232374527891137";
const char kRefAtanhMinusH[]   = "-0.54930614433405484569762261846126285232374527891137";
const char kRefAtanhNineT[]    = "1.4722194895832202300045137159439267686186896306496";
const char kRefAtanh99[]       = "2.6466524123622461977050606459342686009455526402847";
const char kRefAtanhNearOne[]  = "8.4056213907591319237715653187902484272770219427454";
// The reciprocal trigonometric family -- cotangent, secant and
//	cosecant, their hyperbolic counterparts, and the six inverses.
//	Generated, and every one checked by an identity or a round trip
//	through its forward function before being written down. Several
//	land byte-identical on constants already above, which is the
//	check arriving for free: arccot(2) is arctan(1/2), arccsc(2) is
//	arcsin(1/2), and arccsch(1) is arcsinh(1)...
const char kRefCot1[]                  = "0.64209261593433070300641998659426562023027811391817";
const char kRefCotHalf[]               = "1.8304877217124519192680194389688166237581079480161";
const char kRefCotTwo[]                = "-0.45765755436028576375027741043204727642848632923167";
const char kRefCotMinusOne[]           = "-0.64209261593433070300641998659426562023027811391817";
const char kRefSec1[]                  = "1.8508157176809256179117532413986501934703966550940";
const char kRefSecHalf[]               = "1.1394939273245491223133277682049499284237252460490";
const char kRefSecTwo[]                = "-2.4029979617223809897546004014200662262451210931545";
const char kRefCsc1[]                  = "1.1883951057781212162615994523745510035278298340980";
const char kRefCscHalf[]               = "2.0858296429334881857725016754592903019623095868170";
const char kRefCscTwo[]                = "1.0997501702946164667566973970263128966587644431498";
const char kRefCoth1[]                 = "1.3130352854993313036361612469308478329120139412405";
const char kRefCothHalf[]              = "2.1639534137386528487700040102180231170937386021508";
const char kRefCothTen[]               = "1.0000000041223072533738241840280803124601800267562";
const char kRefCothMinusOne[]          = "-1.3130352854993313036361612469308478329120139412405";
const char kRefSech1[]                 = "0.64805427366388539957497735322615032310848931207194";
const char kRefSechHalf[]              = "0.88681888397007390865889779778340856253408908871261";
const char kRefSechTwo[]               = "0.26580222883407969212086273981988897153078265443227";
const char kRefCsch1[]                 = "0.85091812823932154513384276328717528418172466091034";
const char kRefCschHalf[]              = "1.9190347513349437194922028787270061595871798715330";
const char kRefCschMinusOne[]          = "-0.85091812823932154513384276328717528418172466091034";
const char kRefAcot1[]                 = "0.78539816339744830961566084581987572104929234984378";
const char kRefAcotHalf[]              = "1.1071487177940905030170654601785370400700476454014";
const char kRefAcotTwo[]               = "0.46364760900080611621425623146121440202853705428612";
const char kRefAcotMinusOne[]          = "-0.78539816339744830961566084581987572104929234984378";
const char kRefAcotMinusTwo[]          = "-0.46364760900080611621425623146121440202853705428612";
const char kRefAcotMillion[]           = "0.00000099999999999966666666666686666666666652380952380963";
const char kRefAsecTwo[]               = "1.0471975511965977461542144610931676280657231331250";
const char kRefAsecMinusTwo[]          = "2.0943951023931954923084289221863352561314462662501";
const char kRefAsecElevenTenths[]      = "0.42969966615142451742718719630135821756460507096529";
const char kRefAsecMillion[]           = "1.5707953267948966190646550249730097754319179883780";
const char kRefAcscTwo[]               = "0.52359877559829887307710723054658381403286156656252";
const char kRefAcscMinusTwo[]          = "-0.52359877559829887307710723054658381403286156656252";
const char kRefAcscElevenTenths[]      = "1.1410966606434721018041344953383932245339796287223";
const char kRefAcscMillion[]           = "0.0000010000000000001666666666667416666666667113095238096";
const char kRefAcothTwo[]              = "0.54930614433405484569762261846126285232374527891137";
const char kRefAcothMinusTwo[]         = "-0.54930614433405484569762261846126285232374527891137";
const char kRefAcothElevenTenths[]     = "1.5222612188617114982502989901828527171422876437023";
const char kRefAcothMillion[]          = "0.0000010000000000003333333333335333333333334761904761906";
const char kRefAsechHalf[]             = "1.3169578969248167086250463473079684440269819714675";
const char kRefAsechTenth[]            = "2.9932228461263808979126677137741829130836604511810";
const char kRefAsechNineTenths[]       = "0.46714530810326201812838149166592294470091920423633";
// The hyperbolic arcsecant on both sides of the threshold at which it
//	stops summing its own series -- kArcHyperSecantSeriesCeiling...
const char kRefAsechElevenThousandths[] = "5.2029769353710257912798216459701206794025773254287";
const char kRefAsechHundredth[]        = "5.2982923656104845907016668349432471689372518831401";
const char kRefAsechThousandth[]       = "7.6009022095419886114191231179982237481516724238259";
const char kRefAsechHundredThousandth[] = "12.206072645505173729506251894879945522747670826202";

const char kRefAcsch1[]                = "0.88137358701954302523260932497979230902816032826164";
const char kRefAcschMinusOne[]         = "-0.88137358701954302523260932497979230902816032826164";
const char kRefAcschHalf[]             = "1.4436354751788103424932767402731052694055530031570";
const char kRefAcschTwo[]              = "0.48121182505960344749775891342436842313518433438566";
const char kRefAcschMillion[]          = "0.00000099999999999983333333333340833333333328869047619051";

const char kRefPiOverTwo[]		= "1.5707963267948966192313216916397514420985846996876";
const char kRefMinusPiOverTwo[]	= "-1.5707963267948966192313216916397514420985846996876";
const char kRefThird[]			= "0.33333333333333333333333333333333333333333333333333";
const char kRefSeventh[]		= "0.14285714285714285714285714285714285714285714285714";
const char kRefSin2[]			= "0.90929742682568169539601986591174484270225497144789";
const char kRefSin3[]			= "0.14112000805986722210074480280811027984693326425226";
const char kRefSin4[]			= "-0.7568024953079282513726390945118290941359128873364";
const char kRefSin6[]			= "-0.2794154981989258728115554466118947596279948643182";
const char kRef22Over7[]		= "3.14285714285714285714285714285714285714285714285714";

// 1 / 12345678901234567890123456789012345, to 50 significant digits
const char kRefLongDivisor[] = "0.000000000000000000000000000000000081000000729000006"
							   "6339000603684905538075350792915085823480";

/////////////////////////////////////////////////////////////////
//	Test scaffolding
/////////////////////////////////////////////////////////////////

static long	g_passed = 0;
static long	g_failed = 0;
static long	g_known	= 0;

// Every measured accuracy, so the weak spots can be summarized at the end...
const int	kMaxMeasurements = 64;
const int	kAccuracyTarget = 30;

struct SMeasurement
{
	char	szLabel[ 48 ];		// copied: callers often pass a reused buffer
	int		iDigits;
};

static SMeasurement	g_measurements[ kMaxMeasurements ];
static int			g_measurementCount = 0;

static void Record( const char* szLabel, int iDigits )
{
	if ( g_measurementCount < kMaxMeasurements )
	{
		::snprintf( g_measurements[ g_measurementCount ].szLabel, sizeof( g_measurements[ 0 ].szLabel ), "%s", szLabel );
		g_measurements[ g_measurementCount ].iDigits = iDigits;
		g_measurementCount++;
	}
}

// Reduce a decimal string to a sign, a power-of-ten exponent, and a
//	run of significant digits, so two values can be compared without
//	regard to how they were formatted...
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

	// A value of zero...
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

// A reciprocal floor no test division reaches, to force long division...
const long kReciprocalOutOfReach = 100000000;

// And a width no test logarithm reaches, to force the series...
const long kLogNewtonOutOfReach = 100000000;

// How many significant digits of the reference does this value get
//	right? Measured as the ratio of the reference to the difference,
//	so that values straddling a power of ten (.99999 against 1) are
//	scored on how far apart they actually are, not on how their
//	leading digits happen to line up.
//
//	**As long as the two strings are.** Until 2026-09-26 this normalised
//	into 256 characters and lined up in 320, so it compared the first
//	255 digits and called everything after them agreement: CheckExact
//	never saw the last 47 of 2^1000, and a check at 1,200 places passed
//	with the answer wrong from its 985th...
static int AgreementDigits( const char* szGot, const char* szReference )
{
	int	iSize = static_cast<int>( ::strlen( szGot ) + ::strlen( szReference ) ) + 16;

	std::vector<char>	gotDigits( static_cast<size_t>( iSize ) );
	std::vector<char>	refDigits( static_cast<size_t>( iSize ) );

	char*	szGotDigits = gotDigits.data();
	char*	szRefDigits = refDigits.data();

	int		iGotSign, iGotExp;
	int		iRefSign, iRefExp;

	Normalize( szGot, iGotSign, iGotExp, szGotDigits, iSize );
	Normalize( szReference, iRefSign, iRefExp, szRefDigits, iSize );

	// Both zero is perfect agreement; one zero, or opposing signs, is none...
	if ( (iGotSign == 0) && (iRefSign == 0) )
		return kExactAgreement;

	if ( (iGotSign == 0) || (iRefSign == 0) || (iGotSign != iRefSign) )
		return 0;

	// Line both values up on a common power of ten, wide enough for both...
	int		iAlignExp = (iGotExp > iRefExp) ? iGotExp : iRefExp;
	int		iLowExp = (iGotExp < iRefExp) ? iGotExp : iRefExp;
	int		iLongest = static_cast<int>( std::max( ::strlen( szGotDigits ), ::strlen( szRefDigits ) ) );
	int		kWidth = (iAlignExp - iLowExp) + iLongest + 8;

	std::vector<char>	gotAligned( static_cast<size_t>( kWidth ) + 1, '0' );
	std::vector<char>	refAligned( static_cast<size_t>( kWidth ) + 1, '0' );
	std::vector<char>	diff( static_cast<size_t>( kWidth ) + 1, '0' );

	char*	szGotAligned = gotAligned.data();
	char*	szRefAligned = refAligned.data();
	char*	szDiff = diff.data();

	szGotAligned[ kWidth ] = '\0';
	szRefAligned[ kWidth ] = '\0';

	for ( int i = 0; szGotDigits[ i ]; i++ )
	{
		int iPos = (iAlignExp - iGotExp) + i;
		if ( iPos < kWidth ) szGotAligned[ iPos ] = szGotDigits[ i ];
	}

	for ( int i = 0; szRefDigits[ i ]; i++ )
	{
		int iPos = (iAlignExp - iRefExp) + i;
		if ( iPos < kWidth ) szRefAligned[ iPos ] = szRefDigits[ i ];
	}

	// Identical to every digit compared...
	if ( ::memcmp( szGotAligned, szRefAligned, static_cast<size_t>( kWidth ) ) == 0 )
		return kExactAgreement;

	// Subtract the smaller from the larger...
	const char*	pHigh = szGotAligned;
	const char*	pLow = szRefAligned;

	if ( ::memcmp( szGotAligned, szRefAligned, static_cast<size_t>( kWidth ) ) < 0 )
	{
		pHigh = szRefAligned;
		pLow = szGotAligned;
	}

	int		iBorrow = 0;

	szDiff[ kWidth ] = '\0';

	for ( int i = kWidth - 1; i >= 0; i-- )
	{
		int d = (pHigh[ i ] - '0') - (pLow[ i ] - '0') - iBorrow;

		if ( d < 0 )	{ d += 10; iBorrow = 1; }
		else			{ iBorrow = 0; }

		szDiff[ i ] = static_cast<char>('0' + d);
	}

	// The position of the difference's leading digit gives its magnitude...
	for ( int i = 0; i < kWidth; i++ )
	{
		if ( szDiff[ i ] != '0' )
			return iRefExp - (iAlignExp - i);
	}

	return kExactAgreement;
}

static char* Render( CLZNumber* pNumber )
{
	pNumber->TrimTrailing();
	return pNumber->GetStringRep();
}

// The engine agrees with the reference to at least iRequired digits...
static void CheckDigits( const char* szLabel, CLZNumber* pGot, const char* szReference, int iRequired )
{
	char*	szGot = Render( pGot );
	int		iAgree = AgreementDigits( (szGot != nullptr) ? szGot : "", szReference );

	Record( szLabel, iAgree );

	if ( iAgree >= iRequired )
	{
		g_passed++;
		if ( iAgree == kExactAgreement )
			printf( "  pass  %-34s exact\n", szLabel );
		else
			printf( "  pass  %-34s %3d digits\n", szLabel, iAgree );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s %3d digits (wanted %d)\n", szLabel, iAgree, iRequired );
		printf( "          got %s\n          ref %s\n", (szGot != nullptr) ? szGot : "(null)", szReference );
	}

	::free( szGot );
}

// The engine reproduces this value exactly, digit for digit...
static void CheckExact( const char* szLabel, CLZNumber* pGot, const char* szExpected )
{
	char*	szGot = Render( pGot );

	if ( (szGot != nullptr) && (AgreementDigits( szGot, szExpected ) == kExactAgreement) )
	{
		g_passed++;
		printf( "  pass  %-34s exact\n", szLabel );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s not exact\n", szLabel );
		printf( "          got %s\n          expected %s\n", (szGot != nullptr) ? szGot : "(null)", szExpected );
	}

	::free( szGot );
}

// Compare the rendering character for character. CheckExact cannot
//	be used for this: it normalises both sides, so a leading zero
//	disappears before the comparison and the fault reads as passing...
static void CheckRendering( const char* szLabel, CLZNumber* pNumber, const char* szExpected )
{
	char*	szGot = pNumber->GetStringRep();

	if ( (szGot != nullptr) && (::strcmp( szGot, szExpected ) == 0) )
	{
		g_passed++;
		printf( "  pass  %-34s %s\n", szLabel, szGot );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s got %s, wanted %s\n", szLabel,
				szGot ? szGot : "(null)", szExpected );
	}

	if ( szGot )
		::free( szGot );
}

// The same, for a value written with an exponent...
static void CheckStyled( const char* szLabel, CLZNumber* pNumber,
						 RenderStyle rsStyle, const char* szExpected )
{
	char*	szGot = pNumber->GetStringRep( rsStyle );

	if ( (szGot != nullptr) && (::strcmp( szGot, szExpected ) == 0) )
	{
		g_passed++;
		printf( "  pass  %-34s %s\n", szLabel, szGot );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s got %s, wanted %s\n", szLabel,
				szGot ? szGot : "(null)", szExpected );
	}

	if ( szGot )
		::free( szGot );
}

// And for a display of a value at a named width...
static void CheckDisplay( const char* szLabel, CLZNumber* pNumber,
						  RenderStyle rsStyle, long lDigits, const char* szExpected )
{
	char*	szGot = pNumber->GetDisplayRep( rsStyle, lDigits );

	if ( (szGot != nullptr) && (::strcmp( szGot, szExpected ) == 0) )
	{
		g_passed++;
		printf( "  pass  %-34s %s\n", szLabel, szGot );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s got %s, wanted %s\n", szLabel,
				szGot ? szGot : "(null)", szExpected );
	}

	if ( szGot )
		::free( szGot );
}

// A plain assertion, for behaviour that is not a digit count...
static void CheckTrue( const char* szLabel, bool bCondition, const char* szDetail )
{
	if ( bCondition )
	{
		g_passed++;
		printf( "  pass  %-34s %s\n", szLabel, szDetail );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s %s\n", szLabel, szDetail );
	}
}

// Behavior we know to be wrong, recorded so a change is noticed...
[[maybe_unused]] static void CheckKnown( const char* szLabel, CLZNumber* pGot, const char* szCurrent, const char* szNote )
{
	char*	szGot = Render( pGot );
	bool	bStillTrue = szGot && (::strcmp( szGot, szCurrent ) == 0);

	g_known++;
	printf( "  %s  %-34s %s\n", bStillTrue ? "KNOWN" : "CHNGD", szLabel, szNote );

	if ( !bStillTrue )
		printf( "          behavior changed: got %s, previously %s\n", (szGot != nullptr) ? szGot : "(null)", szCurrent );

	::free( szGot );
}

static void Banner( const char* szTitle )
{
	printf( "\n%s\n", szTitle );
}

static void Set( CLZNumber& theNumber, const char* szValue )
{
	theNumber.Assign( szValue );
}

/////////////////////////////////////////////////////////////////
//	Exact arithmetic - the digit primitives must not approximate
/////////////////////////////////////////////////////////////////

static void TestExactArithmetic()
{
	CLZNumber	a, b, result;

	Banner( "Exact arithmetic" );

	Set( a, "2" ); Set( b, "2" ); result.Clear();
	lz_add( &a, &b, &result );
	CheckExact( "2 + 2", &result, "4" );

	Set( a, "99999999999999999999" ); Set( b, "99999999999999999999" ); result.Clear();
	lz_multiply( &a, &b, &result );
	CheckExact( "20-digit square", &result, "9999999999999999999800000000000000000001" );

	Set( a, "123456789" ); Set( b, "987654321" ); result.Clear();
	lz_multiply( &a, &b, &result );
	CheckExact( "123456789 * 987654321", &result, "121932631112635269" );

	Set( a, "1000000000000000000000000000001" ); Set( b, "1" ); result.Clear();
	lz_subtract( &a, &b, &result );
	CheckExact( "carry across 30 digits", &result, "1000000000000000000000000000000" );

	Set( a, "20" ); result.Clear();
	lz_factorial( &a, &result );
	CheckExact( "20!", &result, "2432902008176640000" );

	Set( a, "25" ); Set( b, "12" ); result.Clear();
	lz_combinations( &a, &b, &result );
	CheckExact( "C(25,12)", &result, "5200300" );
}

/////////////////////////////////////////////////////////////////
//	Division - exactness, accuracy, and edge cases
/////////////////////////////////////////////////////////////////

static void TestDivisionExact()
{
	CLZNumber	a, b, result;

	Banner( "Division - terminating quotients must be exact" );

	Set( a, "10" ); Set( b, "2" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "10 / 2", &result, "5" );

	Set( a, "1" ); Set( b, ".5" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / .5", &result, "2" );

	Set( a, "3" ); Set( b, "8" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "3 / 8", &result, ".375" );

	Set( a, "1" ); Set( b, "1024" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / 1024", &result, ".0009765625" );

	Set( a, "0" ); Set( b, "5" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "0 / 5", &result, "0" );

	/*
		A divisor below one, with a quotient that reaches ten.

		This is the shape the suite had no test for, and it was wrong
		for years: a fraction carries a leading zero, that zero was
		counted as part of the divisor's length, the dividend offset
		came out one too small, and the leading quotient digit had to
		hold a value the base has no digit for. 1 / .1 answered "A".
		1 / .03 answered .333 where the answer is 33.3 -- wrong by two
		orders of magnitude and perfectly plausible-looking.

		Every case below has a divisor under one. The first group has
		quotients of ten and up, which is where it broke; the rest walk
		the quotient up through the powers of the base.
	*/

	/*
		Both spellings of the same divisor, because they are not the
		same value in storage and only one of them used to fail. ".1"
		occupies a single digit; "0.1" keeps the zero that was written
		and occupies two, and it was the leading zero that broke it. A
		test written the tidy way would have missed this entirely.
	*/

	Set( a, "1" ); Set( b, ".1" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / .1", &result, "10" );

	/*
		Dividing into the same target repeatedly must not widen it.

		lzi_divide_op computes the quotient wider than asked, because a
		small quotient loses places to its own leading zeros and has to
		be given them back. Assign then copied that wider precision on
		to the target, so the next division into it started wider still
		and widened again. Nothing failed -- the digits were right --
		but the cost climbed with every call, and a hundred and eighty
		divisions into one target, which is what the logarithm's series
		does, carried it from 36 places to 1116 and took ln from two
		milliseconds to ninety.

		The value is allowed to be wider than the precision; the
		precision is what must not drift.
	*/

	{
		CLZNumber	one, divisor, target;

		Set( one, "1" ); Set( divisor, "359" );

		long l_entryPrecision = target.m_lPrecision;

		for ( int i = 0; i < 200; i++ )
			lz_divide( &one, &divisor, &target );

		CheckTrue( "200 divides do not widen the target",
				   target.m_lPrecision == l_entryPrecision, "" );
	}

	// ...and the answer must not depend on how used the target is...
	{
		CLZNumber	one, three, fresh, reused;

		Set( one, "1" ); Set( three, "3" );

		lz_divide( &one, &three, &fresh );

		for ( int i = 0; i < 20; i++ )
			lz_divide( &one, &three, &reused );

		char* szFresh = fresh.GetStringRep();
		char* szReused = reused.GetStringRep();

		CheckTrue( "the twentieth quotient matches the first",
				   (szFresh != nullptr) && (szReused != nullptr) &&
				   (::strcmp( szFresh, szReused ) == 0), "" );

		::free( szFresh );
		::free( szReused );
	}

	Set( a, "1" ); Set( b, "0.1" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / 0.1 (written leading zero)", &result, "10" );

	Set( a, "1.99" ); Set( b, "0.1" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1.99 / 0.1 (written leading zero)", &result, "19.9" );

	Set( a, "5" ); Set( b, "0.1" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "5 / 0.1 (written leading zero)", &result, "50" );

	Set( a, "2" ); Set( b, ".1" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "2 / .1", &result, "20" );

	Set( a, "1.99" ); Set( b, ".1" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1.99 / .1", &result, "19.9" );

	Set( a, "1" ); Set( b, ".01" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / .01", &result, "100" );

	Set( a, "1.99" ); Set( b, ".01" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1.99 / .01", &result, "199" );

	Set( a, "10" ); Set( b, ".01" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "10 / .01", &result, "1000" );

	Set( a, "1" ); Set( b, ".0001" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / .0001", &result, "10000" );

	// ...and one that does not terminate, so the magnitude is checked
	//	as well as the digits...
	Set( a, "1" ); Set( b, ".03" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / .03", &result, "33.333333333333333333333333333333333333333", 39 );

	Set( a, "1" ); Set( b, ".0003" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / .0003", &result, "3333.3333333333333333333333333333333333333", 36 );

	// The reciprocal goes the same road, and had the same fault...
	Set( a, ".01" ); result.Clear();
	lz_reciprocal( &a, &result );
	CheckExact( "reciprocal of .01", &result, "100" );

	Set( a, ".2" ); result.Clear();
	lz_reciprocal( &a, &result );
	CheckExact( "reciprocal of .2", &result, "5" );
}

static void TestDivisionAccuracy()
{
	CLZNumber	a, b, result;

	Banner( "Division - accuracy against exact references" );

	Set( a, "1" ); Set( b, "3" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / 3", &result, kRefThird, 38 );

	Set( a, "1" ); Set( b, "7" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / 7", &result, kRefSeventh, 38 );

	Set( a, "22" ); Set( b, "7" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "22 / 7", &result, kRef22Over7, 38 );

	Set( a, "1" ); Set( b, "12345678901234567890123456789012345" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / 35-digit divisor", &result, kRefLongDivisor, 38 );

	// The same quotient, reached from very different magnitudes...
	Set( a, "1" ); Set( b, "3" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / 3 (unit magnitude)", &result, kRefThird, 38 );

	Set( a, "1000000000000000000000000000000" );
	Set( b, "3000000000000000000000000000000" );
	result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1e30 / 3e30 (same quotient)", &result, kRefThird, 38 );

	Set( a, ".000000000000000000000000000001" );
	Set( b, ".000000000000000000000000000003" );
	result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1e-30 / 3e-30 (same quotient)", &result, kRefThird, 37 );
}

static void TestDivisionSigns()
{
	CLZNumber	a, b, result;

	Banner( "Division - signs" );

	Set( a, "-1" ); Set( b, "7" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "-1 / 7", &result, "-0.14285714285714285714285714285714285714", 36 );

	Set( a, "1" ); Set( b, "-7" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / -7", &result, "-0.14285714285714285714285714285714285714", 36 );

	Set( a, "-1" ); Set( b, "-7" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckDigits( "-1 / -7", &result, kRefSeventh, 38 );
}

static void TestDivisionRoundTrip()
{
	CLZNumber	a, b, quotient, product;

	Banner( "Division - round trips" );

	const char* szDivisors[] = { "3", "7", "13", "97", "12345678901234567890", nullptr };

	for ( int i = 0; szDivisors[ i ] != nullptr; i++ )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "(1 / %s) * %s", szDivisors[ i ], szDivisors[ i ] );

		Set( a, "1" ); Set( b, szDivisors[ i ] );
		quotient.Clear(); product.Clear();

		lz_divide( &a, &b, &quotient );
		lz_multiply( &quotient, &b, &product );

		CheckDigits( szLabel, &product, "1", 30 );
	}
}

static void TestDivisionEdges()
{
	CLZNumber	a, b, result;

	Banner( "Division - edge cases" );

	// Nothing divides by zero: the caller is told so, and the result
	//	carries the flags that say why...
	Set( a, "1" ); Set( b, "0" ); result.Clear();
	CheckTrue( "1 / 0 reports the error", lz_divide( &a, &b, &result ) == kDivideByZero, "" );
	CheckTrue( "1 / 0 flags divide by zero", result.FlagIsOn( kfDiv0 ), "" );
	CheckTrue( "1 / 0 flags infinity", result.FlagIsOn( kfInfinity ), "" );

	Set( a, "-5" ); Set( b, "0" ); result.Clear();
	CheckTrue( "-5 / 0 reports the error", lz_divide( &a, &b, &result ) == kDivideByZero, "" );

	// Zero over zero is undefined rather than infinite...
	Set( a, "0" ); Set( b, "0" ); result.Clear();
	CheckTrue( "0 / 0 is undefined", lz_divide( &a, &b, &result ) == kUndefinedResult, "" );
	CheckTrue( "0 / 0 flags not-a-number", result.FlagIsOn( kfQNaN ), "" );

	/*
		A divisor that rounded away to nothing after it was last
		examined is still a zero divisor. Left undetected it used to
		hang: the division subtracts the divisor until the dividend
		goes negative, and subtracting zero never gets there.
	*/

	CLZNumber tiny, one, underflowed;

	Set( tiny, "0.00000000000000000000000000000000000001" );
	one.Assign( 1L );
	lz_multiply( &tiny, &one, &underflowed );

	CheckTrue( "a rounded-away divisor is zero", underflowed.IsZero(), "" );

	Set( a, "1" ); result.Clear();
	CheckTrue( "dividing by it terminates", lz_divide( &a, &underflowed, &result ) == kDivideByZero, "" );

}

/////////////////////////////////////////////////////////////////
//	Division on packed limbs, against the digit loop it replaced
/////////////////////////////////////////////////////////////////

/*
	Past a small floor, division packs its digits into limbs -- see
	kPackedDivideFloor -- so the division tests above all run packed.
	What they cannot show is that the packed path and the digit loop
	agree wherever both are allowed, or that Algorithm D's rare corner
	is reached at all. This shows both.

	Each division is done four ways: by the digit loop alone, packed at
	the table's limbs, and packed at one and at two digits per limb,
	where the add-back fires constantly. All four must render the same.
	The operands come from a fixed generator, so the battery is the same
	on every run, and a quarter of the divisors are runs of the largest
	digit -- values just under a power of the base, which is the shape
	that makes the estimate overshoot at the table's own limbs.
*/

// Sixty-four bits whatever a long is: on Windows an unsigned long is 32,
//	the generator wrapped there and its shift by 33 was undefined...
static unsigned long long g_packedSeed = 20260923;

static long PackedRandom( long lBelow )
{
	g_packedSeed = (g_packedSeed * 6364136223846793005ULL) + 1442695040888963407ULL;

	return static_cast<long>( (g_packedSeed >> 33) % static_cast<unsigned long long>( lBelow ) );
}

// An operand of lDigits digits in the base's own spelling: random, all
//	the largest digit, or the largest digit then random...
static void PackedOperand( char* szOut, long lBase, long lDigits, long iShape )
{
	const char*	szSpelling = (lBase == 12) ? "0123456789DE" : "0123456789ABCDEF";
	long		l_radixAt = PackedRandom( lDigits + 1 );
	char*		p = szOut;

	if ( PackedRandom( 2 ) == 0 )
		*p++ = '-';

	for ( long i = 0; i < lDigits; i++ )
	{
		long	l_digit = PackedRandom( lBase );

		if ( (iShape == 1) || ((iShape == 2) && (i < (lDigits / 2))) )
			l_digit = lBase - 1;
		else if ( (i == 0) && (l_digit == 0) )
			l_digit = 1;

		if ( (i == l_radixAt) && (i > 0) )
			*p++ = '.';

		*p++ = szSpelling[ l_digit ];
	}

	*p = '\0';
}

static void TestPackedDivision()
{
	Banner( "Division - packed limbs agree with the digit loop" );

	long	l_savedFloor = lzi_set_packed_divide_floor( kPackedFloorMeasured );
	long	l_savedFactor = lzi_set_packed_divide_factor( 0 );
	long	l_savedPrecision = glWorkingPrecision;

	// The four ways, as a floor and a factor...
	struct { long lFloor; long lFactor; } ways[] =
	{
		{ LONG_MAX, 0 },	// the digit loop alone
		{ 0, 0 },			// packed, at the table's limbs
		{ 0, 1 },			// packed, one digit per limb
		{ 0, 2 },			// packed, two
	};

	const long	kWays = sizeof( ways ) / sizeof( ways[ 0 ] );

	// Every base: three to fifteen have divided on limbs since 2026-10-12...
	const long	bases[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	const long	precisions[] = { 5, 36, 100 };
	const long	kDivisions = static_cast<long>( (sizeof( bases ) / sizeof( bases[ 0 ] )) * (sizeof( precisions ) / sizeof( precisions[ 0 ] )) ) * 24;

	long	l_addBacks[ kWays ] = { 0 };
	long	l_packed[ kWays ] = { 0 };

	for ( long lBase : bases )
	{
		long	l_cases = 0;
		long	l_disagreements = 0;
		char	szFirst[ 600 ] = "";

		for ( long lPrecision : precisions )
		{
			for ( long iCase = 0; iCase < 24; iCase++ )
			{
				char	szDividend[ 160 ], szDivisor[ 160 ];

				PackedOperand( szDividend, lBase, 1 + PackedRandom( 60 ), (iCase % 8 == 7) ? 1 : 0 );
				PackedOperand( szDivisor, lBase, 1 + PackedRandom( 40 ), iCase % 4 );

				char*	szAnswers[ kWays ] = { nullptr };

				for ( long w = 0; w < kWays; w++ )
				{
					lzi_set_packed_divide_floor( ways[ w ].lFloor );
					lzi_set_packed_divide_factor( ways[ w ].lFactor );

					// Before the numbers are made, which take their precision from it...
					glWorkingPrecision = lPrecision;

					CLZNumber	dividend, divisor, quotient;
					long		l_divisionsBefore, l_addBacksBefore, l_divisionsAfter, l_addBacksAfter;

					dividend.SetBase( lBase, true );
					divisor.SetBase( lBase, true );
					quotient.SetBase( lBase, true );
					dividend.Assign( szDividend );
					divisor.Assign( szDivisor );

					lzi_packed_divide_counts( &l_divisionsBefore, &l_addBacksBefore );
					lz_divide( &dividend, &divisor, &quotient );
					lzi_packed_divide_counts( &l_divisionsAfter, &l_addBacksAfter );

					l_packed[ w ] += l_divisionsAfter - l_divisionsBefore;
					l_addBacks[ w ] += l_addBacksAfter - l_addBacksBefore;

					szAnswers[ w ] = quotient.GetStringRep();
				}

				l_cases++;

				for ( long w = 1; w < kWays; w++ )
				{
					if ( (szAnswers[ w ] == nullptr) || (szAnswers[ 0 ] == nullptr) ||
						 (::strcmp( szAnswers[ w ], szAnswers[ 0 ] ) != 0) )
					{
						if ( l_disagreements++ == 0 )
							::snprintf( szFirst, sizeof( szFirst ), "%s / %s at %ld places: loop %s, way %ld %s",
										szDividend, szDivisor, lPrecision,
										szAnswers[ 0 ] ? szAnswers[ 0 ] : "(null)", w,
										szAnswers[ w ] ? szAnswers[ w ] : "(null)" );
					}
				}

				for ( long w = 0; w < kWays; w++ )
					::free( szAnswers[ w ] );
			}
		}

		char	szLabel[ 64 ], szDetail[ 700 ];

		::snprintf( szLabel, sizeof( szLabel ), "packed agrees, base %ld", lBase );

		if ( l_disagreements == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "%ld divisions, four ways", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld disagree; first %s",
						l_disagreements, l_cases, szFirst );

		CheckTrue( szLabel, l_disagreements == 0, szDetail );
	}

	glWorkingPrecision = l_savedPrecision;

	/*
		**That the paths were taken.** Without these, a floor that never
		let anything pack, or a factor the kernel ignored, would pass
		every agreement above -- the digit loop agreeing with itself.
	*/

	{
		char	szDetail[ 120 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld packed", l_packed[ 0 ] );
		CheckTrue( "the loop alone packs nothing", l_packed[ 0 ] == 0, szDetail );

		::snprintf( szDetail, sizeof( szDetail ), "%ld, %ld and %ld packed of %ld",
					l_packed[ 1 ], l_packed[ 2 ], l_packed[ 3 ], kDivisions );
		CheckTrue( "the forced ways pack everything",
				   (l_packed[ 1 ] == kDivisions) && (l_packed[ 2 ] == kDivisions) && (l_packed[ 3 ] == kDivisions), szDetail );

		::snprintf( szDetail, sizeof( szDetail ), "%ld and %ld add-backs", l_addBacks[ 2 ], l_addBacks[ 3 ] );
		CheckTrue( "small limbs reach the add-back", (l_addBacks[ 2 ] > 0) && (l_addBacks[ 3 ] > 0), szDetail );
	}

	/*
		**The add-back at the table's own limbs.** On random operands it
		fires about twice per limb base, so at a billion it never does;
		the battery above reached it only through the small limbs. But a
		divisor a whisker off a power of the base normalises to limbs of
		all nines, and that is the shape that makes the estimate
		overshoot -- the sweep reaches it 198 times, dividing values near
		one. One such division per base, found by search at 36 places,
		each answer checked against Python's exact rationals, and each by
		the loop as well. None is exact, so each is written out to its 37
		places -- see TestShortMeansExact. Base two's rounds up to one, and
		was written to 38 until 2026-10-10, recorded here that way: the
		width came from the quotient before it carried. See
		TestCarrySpelled.
	*/

	{
		struct { long lBase; const char* szDividend; const char* szDivisor;
				 const char* szExpect; } addBack[] =
		{
			{ 2, "0.11111111111111111111111111111111111111111111111111111111111111",
			  "0.111111111111111111111111111111111111111111111111111111111111111",
			  "1.0000000000000000000000000000000000000" },
			{ 8, "1.00000000000000000001",
			  "1.000000000000000000004",
			  "1.0000000000000000000040000000000000000" },
			{ 10, "1.000000000000000001",
			  "1.0000000000000000005",
			  "1.0000000000000000004999999999999999998" },
			{ 12, "1.0000000000000001",
			  "1.00000000000000006",
			  "1.00000000000000005EEEEEEEEEEEEEEE90000" },
			{ 16, "1.00000000000001",
			  "1.000000000000008",
			  "1.000000000000007FFFFFFFFFFFFFC00000000" },
		};

		for ( auto& theCase : addBack )
		{
			char	szLabel[ 64 ];
			long	l_before, l_after;

			for ( long w = 0; w < 2; w++ )
			{
				lzi_set_packed_divide_floor( (w == 0) ? 0 : LONG_MAX );
				lzi_set_packed_divide_factor( 0 );
				glWorkingPrecision = 36;

				CLZNumber	dividend, divisor, quotient;

				dividend.SetBase( theCase.lBase, true );
				divisor.SetBase( theCase.lBase, true );
				quotient.SetBase( theCase.lBase, true );
				dividend.Assign( theCase.szDividend );
				divisor.Assign( theCase.szDivisor );

				lzi_packed_divide_counts( nullptr, &l_before );
				lz_divide( &dividend, &divisor, &quotient );
				lzi_packed_divide_counts( nullptr, &l_after );

				glWorkingPrecision = l_savedPrecision;

				::snprintf( szLabel, sizeof( szLabel ), "add-back case, base %ld, %s", theCase.lBase,
							(w == 0) ? "packed" : "loop" );
				CheckRendering( szLabel, &quotient, theCase.szExpect );

				if ( w == 0 )
				{
					::snprintf( szLabel, sizeof( szLabel ), "table limbs add back, base %ld", theCase.lBase );
					CheckTrue( szLabel, l_after > l_before, "" );
				}
			}
		}
	}

	/*
		**An exact half on the guard digit.** The remainder is what the
		sticky digit is made of, and the packed path answers it rather
		than leaving it in the dividend for a scan. Reported as not
		nought when it is, an exact quotient whose guard digit is exactly
		half rounds away from even -- and only then, since half to even
		takes an odd digit up anyway, which is why nothing else here
		noticed a remainder always reported as left over. Found by
		search: an odd numerator over twice a power of the base, at 36
		places. Each answer checked against Python's exact rationals, and
		each by the loop as well.

		**In an odd base it is a tie at every width**, a half being .333...
		in base seven, so the sticky digit is the middle one, told that it
		goes on -- which the packed path, taken by odd bases since
		2026-10-12, answers from twice its remainder against the divisor.
		The divisor is 2 b^37, so the half falls just past the 37th place,
		where the quotient stops; a kept digit of each parity, in bases
		three, five, seven and thirteen.
	*/

	{
		struct { long lBase; const char* szDividend; const char* szDivisor;
				 const char* szExpect; } half[] =
		{
			{ 2, "100010111011101001001100111111111111000100111001111000011101100110101010100001111010110101",
			  "100000000000000000000000000000000000000",
			  "1000101110111010010011001111111111110001001110011110.0001110110011010101010000111101011010" },
			{ 8, "166035211046274463072234337410232671417466310705343643471605",
			  "20000000000000000000000000000000000000",
			  "7301650442313623143511.6157604115334607633144342561721634702" },
			{ 10, "113570378579842288015964968110320493030094614712035792809781",
			  "20000000000000000000000000000000000000",
			  "5678518928992114400798.2484055160246515047307356017896404890" },
			{ 12, "247850076474519866E1796E4813D8764D46599576580468495EDE767651",
			  "20000000000000000000000000000000000000",
			  "123D2603923826DD33569D9.58407E43925232DD8992D0234248EE5993926" },
			{ 16, "246F23965A83DA16EFE00A672BE73E65BC386621C00E38D7E670ECD37981",
			  "20000000000000000000000000000000000000",
			  "123791CB2D41ED0B77F0053.395F39F32DE1C3310E0071C6BF3387669BCC0" },
			{ 3, "120012112112020021121202200210100100022201110001111111012011",
			  "20000000000000000000000000000000000000",
			  "2112102102101001021021.2211220011200011100201112020202002120" },
			{ 3, "111002002200121100121202011122202101220112011222220112021010",
			  "20000000000000000000000000000000000000",
			  "2011211221121020002210.1002022212200221202120222221202122002" },
			{ 5, "140121214202230342242121302303330221441430400422244422032000",
			  "20000000000000000000000000000000000000",
			  "4231033210111242112103.3123401412333220440200211122211013222" },
			{ 5, "421140131113030134220311324400203300311443102310413441323143",
			  "20000000000000000000000000000000000000",
			  "21032004030401231433240.3134422324122403221301130204220411322" },
			{ 7, "365240221402661631065105346015220600451436121662615153242306",
			  "20000000000000000000000000000000000000",
			  "16612011053466431403240.2523006110300224216410664642425121136" },
			{ 7, "562563164360243536462444203532501210626006536544151661011213",
			  "20000000000000000000000000000000000000",
			  "26463143216345525323122.2101614600440313003253255424330340442" },
			{ 13, "C929ACA51C212AAAC35216C2B2606358B8904473B773980A850166316025",
			  "20000000000000000000000000000000000000",
			  "647B565276107BBBC827736.15796982AC44688A1C3A1B405426733173012" },
			{ 13, "3478B9B65CC574B78C5B5054513029B59674C822B2694279099AB70A6A9A",
			  "20000000000000000000000000000000000000",
			  "18A45B599662A25A462C269.227167B5949A26411579B213B04B55A05354C" },
		};

		for ( auto& theCase : half )
		{
			for ( long w = 0; w < 2; w++ )
			{
				lzi_set_packed_divide_floor( (w == 0) ? 0 : LONG_MAX );
				lzi_set_packed_divide_factor( 0 );
				glWorkingPrecision = 36;

				CLZNumber	dividend, divisor, quotient;
				char		szLabel[ 64 ];

				dividend.SetBase( theCase.lBase, true );
				divisor.SetBase( theCase.lBase, true );
				quotient.SetBase( theCase.lBase, true );
				dividend.Assign( theCase.szDividend );
				divisor.Assign( theCase.szDivisor );
				lz_divide( &dividend, &divisor, &quotient );

				glWorkingPrecision = l_savedPrecision;

				::snprintf( szLabel, sizeof( szLabel ), "exact half to even, base %ld, %s", theCase.lBase,
							(w == 0) ? "packed" : "loop" );
				CheckRendering( szLabel, &quotient, theCase.szExpect );
			}
		}
	}

	/*
		**What lzi_divide_op never asks.** It pads the dividend past the
		divisor before dividing, and never divides by nought, so neither
		of these answers is reached through it -- only by calling the
		kernel directly.
	*/

	{
		// Ten digits, least significant first, so two limbs to the dividend's
		//	one -- a divisor of one limb would be short division instead...
		const long	three[] = { 3 };
		const long	divisor[] = { 7, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
		const long	nought[] = { 0 };
		long		quotient[ 2 ] = { 9, 9 };
		bool		b_remainder = false;
		bool		b_answered = lzi_packed_divide( three, 1, divisor, 10, 10, quotient, 2, &b_remainder );

		CheckTrue( "3 / 1000000007 is nought, remainder 3",
				   b_answered && (quotient[ 0 ] == 0) && (quotient[ 1 ] == 0) && b_remainder, "" );

		b_answered = lzi_packed_divide( nought, 1, divisor, 10, 10, quotient, 2, &b_remainder );

		CheckTrue( "0 / 1000000007, nothing left", b_answered && !b_remainder, "" );

		CheckTrue( "a nought divisor is declined",
				   !lzi_packed_divide( divisor, 10, nought, 1, 10, quotient, 2, &b_remainder ), "" );
	}

	/*
		**And that each floor is where it says, to the digit.** A floor
		that packed nothing would cost only time, and nothing else here
		could see it: the answers are the same either way. So each is
		pinned from both sides, one step of the product apart: 7 / 3 is
		22 digits of quotient over one at 19 places and 23 at 20, and in
		base two 1100111011 / 11 is 25 over two at 14 places, 50, and
		1100111011 / 101 is 17 over three at 7, 51. Moving a floor fails
		one of these, which is the point -- it is measured, and a change
		to it should say so here. They were 32 and 256 until the floors
		were measured again on 2026-10-17.
	*/

	{
		struct { const char* szLabel; long lBase; long lPrecision;
				 const char* szDividend; const char* szDivisor; bool bPacks; } gate[] =
		{
			{ "product 22 keeps the loop", 10, 19, "7", "3", false },
			{ "product 23 packs", 10, 20, "7", "3", true },
			{ "base two, product 50, the loop", 2, 14, "1100111011", "11", false },
			{ "base two, product 51, packs", 2, 7, "1100111011", "101", true },
		};

		lzi_set_packed_divide_floor( kPackedFloorMeasured );
		lzi_set_packed_divide_factor( 0 );

		for ( auto& theCase : gate )
		{
			glWorkingPrecision = theCase.lPrecision;

			CLZNumber	dividend, divisor, quotient;
			long		l_before, l_after;

			dividend.SetBase( theCase.lBase, true );
			divisor.SetBase( theCase.lBase, true );
			quotient.SetBase( theCase.lBase, true );
			dividend.Assign( theCase.szDividend );
			divisor.Assign( theCase.szDivisor );

			lzi_packed_divide_counts( &l_before, nullptr );
			lz_divide( &dividend, &divisor, &quotient );
			lzi_packed_divide_counts( &l_after, nullptr );

			glWorkingPrecision = l_savedPrecision;

			CheckTrue( theCase.szLabel, (l_after > l_before) == theCase.bPacks,
					   (l_after > l_before) ? "packed" : "digit loop" );
		}
	}

	lzi_set_packed_divide_floor( l_savedFloor );
	lzi_set_packed_divide_factor( l_savedFactor );

	CheckTrue( "the packed floor and factor are put back",
			   (lzi_set_packed_divide_floor( l_savedFloor ) == l_savedFloor) &&
			   (lzi_set_packed_divide_factor( l_savedFactor ) == l_savedFactor), "" );
}

/////////////////////////////////////////////////////////////////
//	Multiplication's four tiers, against each other
/////////////////////////////////////////////////////////////////

/*
	lz_multiply takes one of four tiers -- a digit at a time, packed limbs,
	Karatsuba on packed limbs, the transform -- by the lengths of its
	operands; see lzi_multiply_tier_for. Every product is exact, so every
	tier must answer the same digits, and this makes each product seven
	ways and requires all seven to render the same: the digit loop alone,
	packed at the table's limbs, packed at one digit per limb folding after
	every product, packed at three folding every two, the transform, and
	Karatsuba at a floor of four limbs, at the table's limbs and at one
	digit per limb folding after every product.

	**The transform and Karatsuba need forcing to be tested at all.** They
	take over at 96 limbs and tens of thousands of digits, and the gates
	make nothing that long, so without the hooks they would be tiers no
	test reaches. At four limbs Karatsuba splits nearly every product here,
	and a digit a limb makes its additions carry nearly every time.
*/

struct SMultiplyWay { const char* szName; long lPacked; long lTransform; long lFactor; long lFold; long lKaratsuba; };

static const SMultiplyWay kMultiplyWays[] =
{
	{ "digit loop",	LONG_MAX,				LONG_MAX,				0, 0, LONG_MAX },
	{ "packed",		0,						LONG_MAX,				0, 0, LONG_MAX },
	{ "one a limb",	0,						LONG_MAX,				1, 1, LONG_MAX },
	{ "three a limb", 0,					LONG_MAX,				3, 2, LONG_MAX },
	{ "transform",	LONG_MAX,				1,						0, 0, LONG_MAX },
	{ "karatsuba",	0,						LONG_MAX,				0, 0, 4 },
	{ "karatsuba, one a limb", 0,			LONG_MAX,				1, 1, 4 },
};

static char* MultiplyOneWay( const SMultiplyWay& theWay, long lBase,
							 const char* szLeft, const char* szRight )
{
	long	l_packed = lzi_set_packed_multiply_floor( theWay.lPacked );
	long	l_transform = lzi_set_transform_floor( theWay.lTransform );
	long	l_factor = lzi_set_packed_multiply_factor( theWay.lFactor );
	long	l_fold = lzi_set_packed_multiply_fold( theWay.lFold );
	long	l_karatsuba = lzi_set_karatsuba_floor( theWay.lKaratsuba );

	CLZNumber	left, right, product;

	left.SetBase( lBase, true );
	right.SetBase( lBase, true );
	product.SetBase( lBase, true );
	left.Assign( szLeft );
	right.Assign( szRight );
	lz_multiply( &left, &right, &product );

	lzi_set_packed_multiply_floor( l_packed );
	lzi_set_transform_floor( l_transform );
	lzi_set_packed_multiply_factor( l_factor );
	lzi_set_packed_multiply_fold( l_fold );
	lzi_set_karatsuba_floor( l_karatsuba );

	return product.GetStringRep();
}

// An operand of lDigits digits: random, a power of the base, all the largest
//	digit, or with a radix mark somewhere in it...
static void MultiplyOperand( char* szOut, long lBase, long lDigits, long iShape )
{
	const char*	szSpelling = (lBase == 12) ? "0123456789DE" : "0123456789ABCDEF";
	char*		p = szOut;
	long		l_radixAt = (iShape == 3) ? PackedRandom( lDigits + 1 ) : 0;

	if ( PackedRandom( 4 ) == 0 )
		*p++ = '-';

	for ( long i = 0; i < lDigits; i++ )
	{
		long	l_digit = PackedRandom( lBase );

		if ( iShape == 1 )
			l_digit = (i == 0) ? 1 : 0;
		else if ( iShape == 2 )
			l_digit = lBase - 1;
		else if ( (i == 0) && (l_digit == 0) )
			l_digit = 1;

		if ( (l_radixAt > 0) && (i == l_radixAt) )
			*p++ = '.';

		*p++ = szSpelling[ l_digit ];
	}

	*p = '\0';
}

static void TestPackedMultiply()
{
	Banner( "Multiplication - the four tiers agree" );

	// From the measured floors, whatever this started with, and back to
	//	that at the end -- a harness may have set them for the whole run...
	long	l_savedPacked = lzi_set_packed_multiply_floor( kPackedFloorMeasured );
	long	l_savedTransform = lzi_set_transform_floor( kPackedFloorMeasured );
	long	l_savedFactor = lzi_set_packed_multiply_factor( 0 );
	long	l_savedFold = lzi_set_packed_multiply_fold( 0 );
	long	l_savedKaratsuba = lzi_set_karatsuba_floor( 0 );
	long	l_savedReach = lzi_set_transform_reach( 0 );

	const long	kWays = sizeof( kMultiplyWays ) / sizeof( kMultiplyWays[ 0 ] );

	// Every base: three to fifteen have had limbs since 2026-10-12, and
	//	until then multiplied a digit at a time at every length...
	const long	bases[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };

	// Lengths from a digit to past the table's natural fold, so the fold is
	//	reached at the table's own limbs as well as when forced...
	const long	lengths[][ 2 ] =
	{
		{ 1, 1 }, { 1, 7 }, { 7, 1 }, { 2, 23 }, { 2, 24 }, { 6, 8 }, { 5, 5 },
		{ 12, 16 }, { 31, 31 }, { 40, 40 }, { 3, 300 }, { 300, 3 }, { 64, 200 },
		{ 200, 200 },
	};

	long	l_before[ 4 ], l_after[ 4 ];
	long	l_wayFolds[ kWays ] = { 0 };
	long	l_splitsBefore = lzi_karatsuba_splits();

	lzi_multiply_tier_counts( &l_before[ 0 ], &l_before[ 1 ], &l_before[ 2 ], &l_before[ 3 ] );

	for ( long lBase : bases )
	{
		long	l_cases = 0;
		long	l_disagreements = 0;
		char	szFirst[ 400 ] = "";

		for ( auto& theLengths : lengths )
		{
			for ( long iShape = 0; iShape < 4; iShape++ )
			{
				static char	szLeft[ 400 ], szRight[ 400 ];

				MultiplyOperand( szLeft, lBase, theLengths[ 0 ], (iShape == 1) ? 0 : iShape );
				MultiplyOperand( szRight, lBase, theLengths[ 1 ], iShape );

				char*	szAnswers[ kWays ];

				for ( long w = 0; w < kWays; w++ )
				{
					long	l_folds0 = 0, l_folds1 = 0;

					lzi_packed_multiply_counts( nullptr, &l_folds0 );
					szAnswers[ w ] = MultiplyOneWay( kMultiplyWays[ w ], lBase, szLeft, szRight );
					lzi_packed_multiply_counts( nullptr, &l_folds1 );

					l_wayFolds[ w ] += l_folds1 - l_folds0;
				}

				l_cases++;

				for ( long w = 1; w < kWays; w++ )
				{
					if ( (szAnswers[ w ] == nullptr) || (szAnswers[ 0 ] == nullptr) ||
						 (::strcmp( szAnswers[ w ], szAnswers[ 0 ] ) != 0) )
					{
						if ( l_disagreements++ == 0 )
							::snprintf( szFirst, sizeof( szFirst ), "%.40s x %.40s: %s against %s",
										szLeft, szRight, kMultiplyWays[ w ].szName, kMultiplyWays[ 0 ].szName );
					}
				}

				for ( long w = 0; w < kWays; w++ )
					::free( szAnswers[ w ] );
			}
		}

		char	szLabel[ 64 ], szDetail[ 500 ];

		::snprintf( szLabel, sizeof( szLabel ), "tiers agree, base %ld", lBase );

		if ( l_disagreements == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "%ld products, %ld ways", l_cases, kWays );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld disagree; first %s",
						l_disagreements, l_cases, szFirst );

		CheckTrue( szLabel, l_disagreements == 0, szDetail );
	}

	lzi_multiply_tier_counts( &l_after[ 0 ], &l_after[ 1 ], &l_after[ 2 ], &l_after[ 3 ] );

	/*
		**That each way took its tier.** A floor the hooks could not move
		would pass every agreement above -- one tier agreeing with itself.
		Every way above makes the same number of products, so the digit
		loop and the transform should each have made one way's share, and
		packed and Karatsuba between them five: three ways pack everything,
		and the two Karatsuba ways pack what is under four limbs. And the
		Karatsuba ways must have split, or the tier ran its leaf alone.
	*/

	{
		long	l_digits = l_after[ 0 ] - l_before[ 0 ];
		long	l_packed = l_after[ 1 ] - l_before[ 1 ];
		long	l_transform = l_after[ 2 ] - l_before[ 2 ];
		long	l_karatsuba = l_after[ 3 ] - l_before[ 3 ];
		long	l_splits = lzi_karatsuba_splits() - l_splitsBefore;
		char	szDetail[ 160 ];

		::snprintf( szDetail, sizeof( szDetail ), "digit loop %ld, packed %ld, transform %ld, karatsuba %ld (%ld splits)",
					l_digits, l_packed, l_transform, l_karatsuba, l_splits );
		CheckTrue( "each way took its own tier",
				   (l_digits > 0) && ((l_packed + l_karatsuba) == 5 * l_digits) && (l_transform == l_digits) &&
				   (l_karatsuba > 0) && (l_splits > 0), szDetail );

		// At one and three digits per limb a column could take millions of
		//	products before folding, so every fold there is the hook's --
		//	counted in total, the table's own natural folds hid a hook
		//	that was being ignored...
		::snprintf( szDetail, sizeof( szDetail ), "%ld and %ld folds", l_wayFolds[ 2 ], l_wayFolds[ 3 ] );
		CheckTrue( "the forced folds were made", (l_wayFolds[ 2 ] > 0) && (l_wayFolds[ 3 ] > 0), szDetail );
	}

	/*
		**A transform that declines hands on to Karatsuba**, not to packed:
		past the transform's reach, packed alone would be quadratic on tens
		of millions of digits. Nothing a test can afford gets there, so the
		reach comes down to one element with the transform forced on, and
		every product must still come out right, each of them by Karatsuba.
	*/

	{
		long	l_cases = 0, l_wrong = 0, l_transform = 0, l_karatsuba = 0;

		lzi_set_transform_reach( 1 );

		for ( long lBase : bases )
		{
			static char	szLeft[ 400 ], szRight[ 400 ];
			long		l_b[ 4 ], l_a[ 4 ];

			MultiplyOperand( szLeft, lBase, 200, 0 );
			MultiplyOperand( szRight, lBase, 64, 0 );

			char*	szWant = MultiplyOneWay( kMultiplyWays[ 0 ], lBase, szLeft, szRight );

			lzi_multiply_tier_counts( &l_b[ 0 ], &l_b[ 1 ], &l_b[ 2 ], &l_b[ 3 ] );

			char*	szGot = MultiplyOneWay( kMultiplyWays[ 4 ], lBase, szLeft, szRight );

			lzi_multiply_tier_counts( &l_a[ 0 ], &l_a[ 1 ], &l_a[ 2 ], &l_a[ 3 ] );

			l_cases++;
			l_transform += l_a[ 2 ] - l_b[ 2 ];
			l_karatsuba += l_a[ 3 ] - l_b[ 3 ];

			if ( (szWant == nullptr) || (szGot == nullptr) || (::strcmp( szWant, szGot ) != 0) )
				l_wrong++;

			::free( szWant );
			::free( szGot );
		}

		lzi_set_transform_reach( 0 );

		char	szDetail[ 120 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld products, %ld wrong; karatsuba %ld, transform %ld",
					l_cases, l_wrong, l_karatsuba, l_transform );
		CheckTrue( "a declined transform goes to Karatsuba",
				   (l_wrong == 0) && (l_karatsuba == l_cases) && (l_transform == 0), szDetail );
	}

	/*
		**The fold, at the table's own limbs, where it is load-bearing.** A
		column of all the largest digit overflows 64 bits one product past
		the fold, and a carry into a column left unsplit overflows once
		enough columns precede it -- in base eight past 31 limbs, in base
		sixteen past 512, and in the bases given limbs on 2026-10-12 past
		16 in base four, 27 in thirteen, 122 in nine and 231 in seven. So
		these are squares of all the largest digit at those lengths, and
		the answer is known without multiplying:
		(b^n - 1)^2 is n - 1 of the largest digit, then b - 2, then n - 1
		noughts, then a 1. The transform checks them too, and the digit
		loop is left out: at 3,600 digits it would cost the whole suite
		again.
	*/

	{
		struct { long lBase; long lDigits; } edges[] =
		{
			{ 8, 170 }, { 8, 330 }, { 10, 171 }, { 12, 800 }, { 16, 1799 }, { 16, 3600 },
			{ 4, 270 }, { 13, 232 }, { 9, 1116 }, { 7, 2330 },
		};

		for ( auto& theEdge : edges )
		{
			const char*	szSpelling = (theEdge.lBase == 12) ? "0123456789DE" : "0123456789ABCDEF";
			std::string	l_operand( static_cast<size_t>( theEdge.lDigits ), szSpelling[ theEdge.lBase - 1 ] );
			std::string	l_expect = std::string( static_cast<size_t>( theEdge.lDigits - 1 ), szSpelling[ theEdge.lBase - 1 ] ) +
								   szSpelling[ theEdge.lBase - 2 ] +
								   std::string( static_cast<size_t>( theEdge.lDigits - 1 ), '0' ) + "1";
			char		szLabel[ 64 ];

			for ( long w : { 1L, 4L, 5L } )
			{
				char*	szGot = MultiplyOneWay( kMultiplyWays[ w ], theEdge.lBase,
												l_operand.c_str(), l_operand.c_str() );

				::snprintf( szLabel, sizeof( szLabel ), "largest digit squared, b%ld n%ld, %s",
							theEdge.lBase, theEdge.lDigits, kMultiplyWays[ w ].szName );
				CheckTrue( szLabel, (szGot != nullptr) && (l_expect == szGot), "" );
				::free( szGot );
			}
		}
	}

	/*
		**And that each floor is where it says, to the digit**, asked of
		the function the dispatch itself asks. A floor that packed nothing
		would cost only time, and nothing else here could see it. The
		transform's is pinned the same way, without multiplying anything
		that long.
	*/

	{
		struct { const char* szLabel; long lBase; long lLeft; long lRight; long lTier; } pins[] =
		{
			{ "one digit keeps the loop",			10, 1, 48, kMultiplyByDigits },
			{ "product 46 keeps the loop",			10, 2, 23, kMultiplyByDigits },
			{ "product 48 packs",					10, 2, 24, kMultiplyPacked },
			{ "6 x 8 packs",						10, 6, 8, kMultiplyPacked },
			{ "base two, product 94, the loop",		2, 2, 47, kMultiplyByDigits },
			{ "base two, product 96, packs",		2, 2, 48, kMultiplyPacked },
			{ "95 limbs square packs",				10, 855, 855, kMultiplyPacked },
			{ "96 limbs square, Karatsuba",			10, 856, 856, kMultiplyKaratsuba },
			{ "96 by 191 limbs, Karatsuba",			10, 864, 1719, kMultiplyKaratsuba },
			{ "96 by 192 limbs packs",				10, 864, 1720, kMultiplyPacked },
			{ "191 limbs, lopsided, packs",			10, 1719, 40000, kMultiplyPacked },
			{ "192 limbs, lopsided, Karatsuba",		10, 1720, 40000, kMultiplyKaratsuba },
			{ "130,869 square, Karatsuba",			10, 130869, 130869, kMultiplyKaratsuba },
			{ "130,870 square, the transform",		10, 130870, 130870, kMultiplyByTransform },
			{ "half the floor, 4x, the transform",	10, 65435, 261740, kMultiplyByTransform },
			{ "under half, 4x, Karatsuba",			10, 65434, 261736, kMultiplyKaratsuba },
			{ "half the floor, under 4x, Karatsuba",	10, 65435, 261739, kMultiplyKaratsuba },
			{ "hex 118,762, Karatsuba",				16, 118762, 118762, kMultiplyKaratsuba },
			{ "hex 118,763, the transform",			16, 118763, 118763, kMultiplyByTransform },
			{ "base two, 95 limbs square, packs",	2, 1900, 1900, kMultiplyPacked },
			{ "base two, 96 limbs, Karatsuba",		2, 1901, 1901, kMultiplyKaratsuba },
			{ "base two, 220,279, Karatsuba",		2, 220279, 220279, kMultiplyKaratsuba },
			{ "base two, 220,280, the transform",	2, 220280, 220280, kMultiplyByTransform },
			{ "base seven, product 46, the loop",	7, 2, 23, kMultiplyByDigits },
			{ "base seven, product 48, packs",		7, 2, 24, kMultiplyPacked },
			{ "base seven, 95 limbs square, packs",	7, 950, 950, kMultiplyPacked },
			{ "base seven, 96 limbs, Karatsuba",	7, 951, 951, kMultiplyKaratsuba },
			{ "base seven, 164,287, Karatsuba",		7, 164287, 164287, kMultiplyKaratsuba },
			{ "base seven, 164,288, the transform",	7, 164288, 164288, kMultiplyByTransform },
			{ "base three, 95 limbs, packs",		3, 1710, 1710, kMultiplyPacked },
			{ "base three, 96 limbs, Karatsuba",	3, 1711, 1711, kMultiplyKaratsuba },
			{ "base eleven, 207,709, Karatsuba",	11, 207709, 207709, kMultiplyKaratsuba },
			{ "base eleven, 207,710, the transform",	11, 207710, 207710, kMultiplyByTransform },
		};

		for ( auto& thePin : pins )
		{
			long	l_tier = lzi_multiply_tier_for( thePin.lBase, thePin.lLeft, thePin.lRight );
			char	szDetail[ 32 ];

			::snprintf( szDetail, sizeof( szDetail ), "tier %ld", l_tier );
			CheckTrue( thePin.szLabel, l_tier == thePin.lTier, szDetail );
		}
	}

	lzi_set_packed_multiply_floor( l_savedPacked );
	lzi_set_transform_floor( l_savedTransform );
	lzi_set_packed_multiply_factor( l_savedFactor );
	lzi_set_packed_multiply_fold( l_savedFold );
	lzi_set_karatsuba_floor( l_savedKaratsuba );
	lzi_set_transform_reach( l_savedReach );

	CheckTrue( "the multiply floors are put back",
			   (lzi_set_packed_multiply_floor( l_savedPacked ) == l_savedPacked) &&
			   (lzi_set_transform_floor( l_savedTransform ) == l_savedTransform) &&
			   (lzi_set_packed_multiply_factor( l_savedFactor ) == l_savedFactor) &&
			   (lzi_set_packed_multiply_fold( l_savedFold ) == l_savedFold) &&
			   (lzi_set_karatsuba_floor( l_savedKaratsuba ) == l_savedKaratsuba) &&
			   (lzi_set_transform_reach( l_savedReach ) == l_savedReach), "" );
}

/////////////////////////////////////////////////////////////////
//	Rounding. The default decides the last digit of very nearly
//	every inexact result the engine produces.
/////////////////////////////////////////////////////////////////

static void TestRounding()
{
	CLZNumber	a, b, result;

	Banner( "Rounding" );

	// A discarded part below half leaves the last digit alone...
	Set( a, "1" ); Set( b, "3" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "1 / 3 rounds down", &result, ".33333333333333333333333333333333333333" );

	// Above half, it carries...
	Set( a, "2" ); Set( b, "3" ); result.Clear();
	lz_divide( &a, &b, &result );
	CheckExact( "2 / 3 rounds up", &result, ".66666666666666666666666666666666666667" );

	/*
		**The last place of a division, pinned.** Long division computes
		one digit past the precision and rounds from it. Until these
		existed, "2 / 3 rounds up" was the only test that could tell
		whether it had: a mutation computing one digit fewer changed the
		last place of 1,406 of 3,016 quotients in a random corpus and
		failed that one test alone. These five are drawn from the
		quotients it changed, across two precisions, a negative, a large
		quotient and two small ones whose extra places are part of the
		answer. Each expected string is the correctly rounded quotient,
		checked in Python's decimal rather than copied from the engine.
	*/

	{
		struct { const char* szLabel; long lPrecision; const char* szDividend;
				 const char* szDivisor; const char* szExpect; } lastPlace[] =
		{
			{ "last place: negative, 17 places", 17,
			  "1811472116958961462189", "-2536831893397623262305",
			  "-0.7140686466744255637" },
			{ "last place: large quotient", 17,
			  "1999015829655172583135912", "836621583631981332",
			  "2389390.697974764377028283" },
			{ "last place: short operands, 36 places", 36,
			  "-32.5382", "326",
			  "-0.099810429447852760736196319018404907975" },
			{ "last place: small quotient", 36,
			  "16910", "2705502033.1741427",
			  "0.0000062502263138796793926082942091202012429" },
			{ "last place: tiny quotient", 17,
			  "79", "484796326216449165759861046032",
			  "0.0000000000000000000000000001629550302423878510" },
		};

		for ( auto& theCase : lastPlace )
		{
			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lPrecision;

			CLZNumber	dividend, divisor, quotient;

			dividend.Assign( theCase.szDividend );
			divisor.Assign( theCase.szDivisor );
			lz_divide( &dividend, &divisor, &quotient );

			glWorkingPrecision = l_saved;

			CheckRendering( theCase.szLabel, &quotient, theCase.szExpect );
		}
	}

	/*
		**A dividend with more places than the answer's width.** Since
		2026-09-30 a quotient's width is set by its value, so a long
		dividend no longer widens it: the digits are computed to the
		dividend's places and rounded once at the width. The first is a
		quotient whose digit past the width is a nought; it trimmed to
		0.4 until 2026-10-01, and since a short quotient is an exact one
		it keeps its seven places, 1.2 / 3 being the one that is 0.4.
		Each expected string was computed in Python's decimal.
	*/

	{
		struct { const char* szDividend; const char* szDivisor; const char* szExpect; } longDividend[] =
		{
			{ "1.2000000000000000001", "3", "0.4000000" },
			{ "1.23456789012345678901", "7", "0.1763668" },
			{ "-98765.432109876543210987", "0.007", "-14109347.444268" },
		};

		for ( auto& theCase : longDividend )
		{
			long	l_saved = glWorkingPrecision;
			char	szLabel[ 80 ];

			// Before the numbers are made, which take the width they are made at...
			glWorkingPrecision = 5;

			CLZNumber	dividend, divisor, quotient;

			dividend.Assign( theCase.szDividend );
			divisor.Assign( theCase.szDivisor );
			lz_divide( &dividend, &divisor, &quotient );
			glWorkingPrecision = l_saved;

			::snprintf( szLabel, sizeof( szLabel ), "long dividend: %.24s / %s", theCase.szDividend, theCase.szDivisor );
			CheckRendering( szLabel, &quotient, theCase.szExpect );
		}
	}

	/*
		**Past a guard digit of exactly half, the rest decides.** Long
		division rounds from one truncated guard digit, and until it
		carried a sticky digit -- whether anything was left over -- a
		guard of exactly half could not tell a quotient ending there from
		one going on, and half-to-even took some of the second kind down.
		Each of these came out one place low; each expected string is the
		correctly rounded quotient, checked in Python's decimal.

		The octal one is the textbook case. Four is half in base eight,
		so a quotient of repeating fours always ends on a guard digit of
		exactly half with more fours behind it -- which only the sticky
		digit can see.
	*/

	{
		struct { const char* szLabel; long lBase; long lPrecision; const char* szDividend;
				 const char* szDivisor; const char* szExpect; } sticky[] =
		{
			{ "sticky: base ten, 17 places", 10, 17,
			  "340667155272307418832668", "926526086297994593181",
			  "367.682205941409521689" },
			{ "sticky: base ten, over seven", 10, 5,
			  "-33311195625953816970351824", "7",
			  "-4758742232279116710050260.571429" },
			{ "sticky: base sixteen", 16, 36,
			  "F52D71199B50A.8896C527E92D9", "271A2160F99AF43E45",
			  "0.00006452B51742BF91575A2A3846D2D85497952CC9" },
			{ "sticky: base two", 2, 36,
			  "11001110", "111",
			  "11101.0110110110110110110110110110110110111" },
			{ "sticky: base eight, repeating fours", 8, 5,
			  "20611256673", "7",
			  "2312364321.444445" },
		};

		for ( auto& theCase : sticky )
		{
			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lPrecision;

			CLZNumber	dividend, divisor, quotient;

			dividend.SetBase( theCase.lBase, true );
			divisor.SetBase( theCase.lBase, true );
			quotient.SetBase( theCase.lBase, true );
			dividend.Assign( theCase.szDividend );
			divisor.Assign( theCase.szDivisor );
			lz_divide( &dividend, &divisor, &quotient );

			glWorkingPrecision = l_saved;

			CheckRendering( theCase.szLabel, &quotient, theCase.szExpect );
		}
	}

	/*
		**The Newton reciprocal, in every base.** The floor is a million
		places, so no ordinary division reaches this path; the floor is
		lowered to twenty for these and put back. Until 2026-09-21 the
		reciprocal worked in base ten whatever base it was given -- its
		constants were made in the default base -- and then scaled its
		answer by the divisor's exponent in the divisor's own base, so
		the base-eight case here came out 108.4 times too small. Base ten
		is the control, and was always right. Each expected string was
		checked against Python's decimal, and against long division's
		answer for the same division.

		**Since 2026-09-26 the reciprocal is exact**, and writes what
		long division writes, character for character: the digits one
		place past the precision and a sticky digit, rounded by the same
		code. Bases two and eight here used to lose a trailing nought
		long division keeps -- the same value, spelt shorter. So each
		answer is now compared with long division's as well.
	*/

	{
		struct { long lBase; const char* szDividend; const char* szDivisor;
				 const char* szExpect; } reciprocal[] =
		{
			{ 2, "11101000110110001000.010010111101110101001011101101010011110001011", "1101101001000010111001",
			  "0.010001000100011011001100011110011100001" },
			{ 8, "25154642326044745607.153715473071721311235534003", "1025107637474421206514",
			  "0.024276122110753656245224640554210400656" },
			{ 10, "1234567.891", "9876.54321",
			  "124.9999989621875000129726562498378417969" },
			{ 12, "E5D2819.7DE3", "3D8E1",
			  "2E4.807E9078596E2E62136D230342D8081D8E8D7" },
			{ 16, "B3925E", "9C2A9D7865D11C2ABA.1E143F323AF529",
			  "0.000000000001265E2341B528CD8A2D7D87E4195EFDC1483C6" },
		};

		long	l_savedFloor = lzi_set_reciprocal_floor( 20 );

		for ( auto& theCase : reciprocal )
		{
			CLZNumber	dividend, divisor, quotient;
			char		szLabel[ 64 ];

			::snprintf( szLabel, sizeof( szLabel ), "reciprocal path in base %ld", theCase.lBase );

			dividend.SetBase( theCase.lBase, true );
			divisor.SetBase( theCase.lBase, true );
			quotient.SetBase( theCase.lBase, true );
			dividend.Assign( theCase.szDividend );
			divisor.Assign( theCase.szDivisor );
			lz_divide( &dividend, &divisor, &quotient );

			CheckRendering( szLabel, &quotient, theCase.szExpect );

			// ...and long division's answer, to the character...
			CLZNumber	longQuotient;
			long		l_reciprocalFloor = lzi_set_reciprocal_floor( kReciprocalOutOfReach );

			longQuotient.SetBase( theCase.lBase, true );
			lz_divide( &dividend, &divisor, &longQuotient );
			lzi_set_reciprocal_floor( l_reciprocalFloor );

			char*	szReciprocal = quotient.GetStringRep();
			char*	szLong = longQuotient.GetStringRep();

			::snprintf( szLabel, sizeof( szLabel ), "reciprocal is long division, base %ld", theCase.lBase );
			CheckTrue( szLabel, ::strcmp( szReciprocal, szLong ) == 0, "" );

			::free( szReciprocal );
			::free( szLong );
		}

		lzi_set_reciprocal_floor( l_savedFloor );

		/*
			**Every mode, every base, both paths.** Random divisions in
			bases two, three, seven, eight, ten, twelve and sixteen --
			signs, fractions, lopsided lengths, exact quotients among
			them -- each divided in all seven modes through the
			reciprocal and through long division, and the two required
			to be identical. Until the reciprocal was made exact it
			rounded half-even whatever its target asked, and the same
			comparison over 10,500 divisions found 2,642 wrong in value.
		*/

		{
			unsigned long long	l_state = 88172645463325252ULL;

			auto Next = [ &l_state ]() -> unsigned long long
			{
				l_state ^= l_state << 13;
				l_state ^= l_state >> 7;
				l_state ^= l_state << 17;

				return l_state;
			};

			auto Random = [ & ]( long lBase, int iWhole, int iFraction ) -> std::string
			{
				const char*	szDigits = "0123456789ABCDEF";
				std::string	s_value = ((Next() % 2) != 0) ? "-" : "";

				if ( iWhole == 0 )
					s_value += "0";
				else
				{
					s_value += szDigits[ 1 + (Next() % static_cast<unsigned long long>( lBase - 1 )) ];

					for ( int i = 1; i < iWhole; i++ )
						s_value += szDigits[ Next() % static_cast<unsigned long long>( lBase ) ];
				}

				if ( iFraction > 0 )
				{
					s_value += ".";

					for ( int i = 0; i < iFraction; i++ )
						s_value += szDigits[ Next() % static_cast<unsigned long long>( lBase ) ];

					if ( iWhole == 0 )
						s_value[ s_value.size() - 1 ] = szDigits[ 1 + (Next() % static_cast<unsigned long long>( lBase - 1 )) ];
				}

				return s_value;
			};

			const long			bases[] = { 2, 3, 7, 8, 10, 12, 16 };
			const RoundingMode	modes[] = { RM_round_down, RM_round_half_up, RM_round_half_even, RM_round_ceiling,
											RM_round_floor, RM_round_half_down, RM_round_up };
			long				l_divisions = 0;
			long				l_differ = 0;
			char				szWhere[ 96 ] = "";

			for ( int iCase = 0; iCase < 120; iCase++ )
			{
				long		lBase = bases[ Next() % 7 ];
				long		lPlaces = 5 + static_cast<long>( Next() % 120 );
				int			iLeftWhole = static_cast<int>( Next() % 20 ), iLeftFraction = static_cast<int>( Next() % 40 );
				int			iRightWhole = static_cast<int>( Next() % 20 ), iRightFraction = static_cast<int>( Next() % 40 );

				if ( (iLeftWhole + iLeftFraction) == 0 )
					iLeftWhole = 1;

				if ( (iRightWhole + iRightFraction) == 0 )
					iRightWhole = 1;

				std::string	s_left = Random( lBase, iLeftWhole, iLeftFraction );
				std::string	s_right = Random( lBase, iRightWhole, iRightFraction );
				bool		b_exact = (Next() % 8) == 0;

				for ( auto rmMode : modes )
				{
					std::string	s_answers[ 2 ];

					for ( int iPath = 0; iPath < 2; iPath++ )
					{
						long	l_floorWas = lzi_set_reciprocal_floor( (iPath == 0) ? 2 : kReciprocalOutOfReach );
						long	l_saved = glWorkingPrecision;

						glWorkingPrecision = lPlaces;

						CLZNumber	left, right, answer;

						left.SetBase( lBase, true );
						left.Assign( s_left.c_str() );
						right.SetBase( lBase, true );
						right.Assign( s_right.c_str() );

						if ( b_exact )
						{
							CLZNumber	product;

							lz_multiply( &left, &right, &product );
							left.Assign( &product );
						}

						answer.SetBase( lBase, true );
						answer.SetPrecision( lPlaces );
						answer.SetRoundingMode( rmMode );
						lz_divide( &left, &right, &answer );

						glWorkingPrecision = l_saved;
						lzi_set_reciprocal_floor( l_floorWas );

						char*	szAnswer = answer.GetStringRep();

						s_answers[ iPath ] = szAnswer;
						::free( szAnswer );
					}

					l_divisions++;

					if ( s_answers[ 0 ] != s_answers[ 1 ] )
					{
						l_differ++;
						::snprintf( szWhere, sizeof( szWhere ), "base %ld, mode %d", lBase, static_cast<int>( rmMode ) );
					}
				}
			}

			char	szDetail[ 96 ];

			if ( l_differ == 0 )
				::snprintf( szDetail, sizeof( szDetail ), "all %ld divisions identical", l_divisions );
			else
				::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, one in %s", l_differ, l_divisions, szWhere );

			CheckTrue( "reciprocal and long division, every mode", l_differ == 0, szDetail );
		}

		/*
			**An exact quotient whose truncation lands a unit high**, which
			120 random divisions did not meet and 10,500 met seven times:
			(a*b)/b in base eight at 19 places. The reciprocal's quotient
			truncates to a + a unit there, the remainder comes out
			negative, and it takes the correction to bring it back. Without
			that correction this is the one division here that differs.
		*/

		{
			const RoundingMode	modes[] = { RM_round_down, RM_round_half_up, RM_round_half_even, RM_round_ceiling,
											RM_round_floor, RM_round_half_down, RM_round_up };
			long				l_differ = 0;

			for ( auto rmMode : modes )
			{
				std::string	s_answers[ 2 ];

				for ( int iPath = 0; iPath < 2; iPath++ )
				{
					long	l_floorWas = lzi_set_reciprocal_floor( (iPath == 0) ? 2 : kReciprocalOutOfReach );
					long	l_saved = glWorkingPrecision;

					glWorkingPrecision = 19;

					CLZNumber	left, right, product, answer;

					left.SetBase( 8, true );
					left.Assign( "-353315640630304500447.6313562304210567" );
					right.SetBase( 8, true );
					right.Assign( "-3153544447437634757522532462.26677572544147203107430571762143723563117272312131" );
					lz_multiply( &left, &right, &product );

					answer.SetBase( 8, true );
					answer.SetPrecision( 19 );
					answer.SetRoundingMode( rmMode );
					lz_divide( &product, &right, &answer );

					glWorkingPrecision = l_saved;
					lzi_set_reciprocal_floor( l_floorWas );

					char*	szAnswer = answer.GetStringRep();

					s_answers[ iPath ] = szAnswer;
					::free( szAnswer );
				}

				if ( s_answers[ 0 ] != s_answers[ 1 ] )
					l_differ++;
			}

			CheckTrue( "base 8: an exact quotient a unit high", l_differ == 0,
					   (l_differ == 0) ? "corrected, in all seven modes" : "the reciprocal left it a unit high" );
		}

		/*
			**An exact half in an odd base**, which random divisions never
			meet: 1/2 in base seven is .333... with the middle digit
			repeating forever, exactly half a place past every digit. The
			reciprocal has to find its remainder exactly half the divisor
			and say the middle digit goes on, or the half modes break the
			tie as though it fell short. Every mode, through both paths.
		*/

		{
			const RoundingMode	modes[] = { RM_round_down, RM_round_half_up, RM_round_half_even, RM_round_ceiling,
											RM_round_floor, RM_round_half_down, RM_round_up, RM_round_half_odd };
			long				l_differ = 0;

			for ( auto rmMode : modes )
			{
				std::string	s_answers[ 2 ];

				for ( int iPath = 0; iPath < 2; iPath++ )
				{
					long	l_floorWas = lzi_set_reciprocal_floor( (iPath == 0) ? 2 : kReciprocalOutOfReach );

					CLZNumber	one, two, answer;

					one.SetBase( 7, true );
					one.Assign( "1" );
					two.SetBase( 7, true );
					two.Assign( "2" );
					answer.SetBase( 7, true );
					answer.SetRoundingMode( rmMode );
					lz_divide( &one, &two, &answer );

					lzi_set_reciprocal_floor( l_floorWas );

					char*	szAnswer = answer.GetStringRep();

					s_answers[ iPath ] = szAnswer;
					::free( szAnswer );
				}

				if ( s_answers[ 0 ] != s_answers[ 1 ] )
					l_differ++;
			}

			CheckTrue( "base 7: 1/2 through both paths, eight modes", l_differ == 0,
					   (l_differ == 0) ? "identical" : "the reciprocal broke a tie long division did not" );
		}

		/*
			**And the shapes it now takes by default.** A long divisor past
			its base's width goes to the reciprocal with no floor lowered --
			see lzi_reciprocal_floors -- so these compare that with long
			division, the shape rule turned off: 1/d over a divisor as long
			as the quotient, where each base starts taking it, in floor and
			half-even. The floors moved on 2026-10-02, and on 2026-10-16 the
			rule became a staircase over the quotient's significant digits
			and the divisor's; in this shape, the divisor a digit longer
			than the places, each base starts at its first step that both
			reach -- see lzi_reciprocal_by_shape.
		*/

		{
			struct { long lBase; long lPlaces; } shapes[] = { { 10, 5999 }, { 16, 4982 }, { 12, 5559 }, { 8, 9965 },
															  { 7, 10649 }, { 13, 5385 }, { 3, 25150 }, { 2, 79725 } };

			for ( auto& theShape : shapes )
			{
				std::string	s_divisor = "1.";

				for ( long i = 0; i < theShape.lPlaces; i++ )
					s_divisor += "0123456789ABCDEF"[ (i * 7 + 3) % theShape.lBase ];

				for ( RoundingMode rmMode : { RM_round_floor, RM_round_half_even } )
				{
					std::string	s_answers[ 2 ];

					for ( int iPath = 0; iPath < 2; iPath++ )
					{
						bool	b_shapeWas = lzi_set_reciprocal_by_shape( iPath == 0 );
						long	l_saved = glWorkingPrecision;

						glWorkingPrecision = theShape.lPlaces;

						CLZNumber	one, divisor, answer;

						one.SetBase( theShape.lBase, true );
						one.Assign( "1" );
						divisor.SetBase( theShape.lBase, true );
						divisor.Assign( s_divisor.c_str() );
						answer.SetBase( theShape.lBase, true );
						answer.SetPrecision( theShape.lPlaces );
						answer.SetRoundingMode( rmMode );
						lz_divide( &one, &divisor, &answer );

						glWorkingPrecision = l_saved;
						lzi_set_reciprocal_by_shape( b_shapeWas );

						char*	szAnswer = answer.GetStringRep();

						s_answers[ iPath ] = szAnswer;
						::free( szAnswer );
					}

					char	szLabel[ 72 ];

					::snprintf( szLabel, sizeof( szLabel ), "by shape, base %ld at %ld, mode %d",
								theShape.lBase, theShape.lPlaces, static_cast<int>( rmMode ) );
					CheckTrue( szLabel, s_answers[ 0 ] == s_answers[ 1 ], "the reciprocal's answer is long division's" );
				}
			}
		}

		/*
			**And the floors themselves.** Only the reciprocal makes
			Karatsuba products -- long division multiplies a digit at a
			time -- so 1/d, the divisor as long as the quotient, makes some
			where its base starts taking it and none a thousand places
			below. Measured 2026-10-16 over a grid of quotient and divisor
			digits; see lzi_reciprocal_by_shape.
		*/

		{
			struct { long lBase; long lFloor; } floors[] = { { 10, 5999 }, { 16, 4982 }, { 12, 5559 }, { 8, 9965 },
															 { 7, 10649 }, { 13, 5385 }, { 3, 25150 }, { 2, 79725 } };

			for ( auto& theFloor : floors )
			{
				long	la_products[ 2 ] = { 0, 0 };

				for ( int iSide = 0; iSide < 2; iSide++ )
				{
					long	l_places = (iSide == 0) ? theFloor.lFloor : (theFloor.lFloor - 1000);
					long	l_saved = glWorkingPrecision;
					long	l_digits = 0, l_packed = 0, l_transform = 0, l_before = 0, l_after = 0;
					std::string	s_divisor = "1.";

					for ( long i = 0; i < l_places; i++ )
						s_divisor += "0123456789ABCDEF"[ (i * 7 + 3) % theFloor.lBase ];

					glWorkingPrecision = l_places;

					CLZNumber	one, divisor, answer;

					one.SetBase( theFloor.lBase, true );
					one.Assign( "1" );
					divisor.SetBase( theFloor.lBase, true );
					divisor.Assign( s_divisor.c_str() );
					answer.SetBase( theFloor.lBase, true );

					lzi_multiply_tier_counts( &l_digits, &l_packed, &l_transform, &l_before );
					lz_divide( &one, &divisor, &answer );
					lzi_multiply_tier_counts( &l_digits, &l_packed, &l_transform, &l_after );
					glWorkingPrecision = l_saved;

					la_products[ iSide ] = l_after - l_before;
				}

				char	szLabel[ 72 ];
				char	szDetail[ 72 ];

				::snprintf( szLabel, sizeof( szLabel ), "base %ld: the reciprocal from %ld", theFloor.lBase, theFloor.lFloor );
				::snprintf( szDetail, sizeof( szDetail ), "%ld Karatsuba products there, %ld a thousand below",
							la_products[ 0 ], la_products[ 1 ] );
				CheckTrue( szLabel, (la_products[ 0 ] > 0) && (la_products[ 1 ] == 0), szDetail );
			}
		}

		CheckTrue( "the reciprocal floor is put back",
				   lzi_set_reciprocal_floor( l_savedFloor ) == l_savedFloor, "" );
	}

	// An exact half goes to the even digit, either way...
	Set( a, ".125" );
	a.Round( 2L );
	CheckExact( "half to even, down (.125)", &a, ".12" );

	Set( a, ".135" );
	a.Round( 2L );
	CheckExact( "half to even, up (.135)", &a, ".14" );

	// Every digit a nine: the carry has to run off the top...
	Set( a, ".99999999999999999999999999999999999" );
	a.Round( 32L );
	CheckExact( "carry out of the leading digit", &a, "1" );

	Set( a, "9.9999999999999999999999999999999999" );
	a.Round( 32L );
	CheckExact( "carry lengthening the value", &a, "10" );
}

static void TestPrecisionControl()
{
	CLZNumber	a, b, result;

	Banner( "Precision control" );

	// Asking for more precision should deliver more correct digits...
	Set( a, "1" ); Set( b, "7" );
	result.Clear();
	result.SetPrecision( 48 );
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / 7 at precision 48", &result, kRefSeventh, 44 );

	Set( a, "1" ); Set( b, "3" );
	result.Clear();
	result.SetPrecision( 48 );
	lz_divide( &a, &b, &result );
	CheckDigits( "1 / 3 at precision 48", &result, kRefThird, 44 );

	/*
		The logarithm's series depth has to track the working
		precision, and until it did it was a flat 180 terms. At the
		default width that is far more than enough, so nothing here
		could see the difference -- which is why this test exists at
		two hundred places, where a flat count is 27 digits short.

		Two, deliberately: it is under e, so it goes through no
		argument reduction and is not capped by the fifty digits
		kBaseE carries. Anything that does reduce stops at about 49
		digits however wide the precision is set, which is a limit of
		the constant and not of this.
	*/

	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = 200;

		CLZNumber	wide;
		CLZNumber	wideResult;

		wide.SetPrecision( 200 );
		wide.Assign( "2" );
		wideResult.SetPrecision( 200 );

		lz_log_e( &wide, &wideResult );

		glWorkingPrecision = l_saved;

		CheckDigits( "ln(2) at precision 200", &wideResult, kRefLn2Wide, 201 );
	}

	/*
		And the same question asked of kBaseE rather than of the series
		depth. Ten reduces where two does not -- the reduction divides
		by e -- so this one is bounded by how many digits of e the
		table carries, and stopped flat at 49 however wide the caller
		set the precision until kBaseE was widened to 120 places.

		A hundred rather than two hundred: at two hundred the constant
		is the limit again, which is a true constant-width ceiling and
		is not gated here.
	*/

	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = 100;

		CLZNumber	wide;
		CLZNumber	wideResult;

		wide.SetPrecision( 100 );
		wide.Assign( "10" );
		wideResult.SetPrecision( 100 );

		lz_log_e( &wide, &wideResult );

		glWorkingPrecision = l_saved;

		CheckDigits( "ln(10) at precision 100", &wideResult, kRefLn10Wide, 101 );
	}
}

/////////////////////////////////////////////////////////////////
//	Magnitude comparison, which decides operand order inside every
//	addition and subtraction. Comparing whole-digit counts rather
//	than significant ones made the answer depend on whether a value
//	had been written with a leading zero.
/////////////////////////////////////////////////////////////////

static void TestMagnitudeComparison()
{
	Banner( "Magnitude comparison" );

	struct { const char* szLeft; const char* szRight; long lExpected; const char* szLabel; } cases[] =
	{
		{ ".84",  "0.5",  kGreaterThan, ".84 > 0.5" },
		{ "0.84", "0.5",  kGreaterThan, "0.84 > 0.5" },
		{ ".84",  ".5",   kGreaterThan, ".84 > .5" },
		{ ".9",   "0.1",  kGreaterThan, ".9 > 0.1" },
		{ ".5",   "0.84", kLessThan,    ".5 < 0.84" },
		{ "0.5",  ".84",  kLessThan,    "0.5 < .84" },
		{ "0.50", ".5",   kEqualTo,     "0.50 == .5" },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber left, right;

		left.Assign( theCase.szLeft );
		right.Assign( theCase.szRight );

		if ( left.CompareABS( &right ) == theCase.lExpected )
		{
			g_passed++;
			printf( "  pass  %-34s\n", theCase.szLabel );
		}
		else
		{
			g_failed++;
			printf( "  FAIL  %-34s\n", theCase.szLabel );
		}
	}

	// A term too small to survive the working precision must subtract
	//	as a zero, not reverse the operand order and wreck the total...
	CLZNumber tiny, one, underflowed, running, result;

	tiny.Assign( "0.00000000000000000000000000000000000001" );
	one.Assign( 1L );
	lz_multiply( &tiny, &one, &underflowed );

	running.Assign( ".84147098480789650665250232163031" );
	result.Clear();
	lz_subtract( &running, &underflowed, &result );

	CheckDigits( "subtracting an underflowed term", &result,
				 ".84147098480789650665250232163031", 30 );
}

/////////////////////////////////////////////////////////////////
//	Transcendentals against references
/////////////////////////////////////////////////////////////////

static void TestTranscendentals()
{
	CLZNumber	a, result;

	Banner( "Transcendental functions against references" );

	Set( a, "2" ); result.Clear();
	lz_root_2( &a, &result );
	CheckDigits( "sqrt(2)", &result, kRefSqrt2, 37 );

	Set( a, "3" ); result.Clear();
	lz_root_2( &a, &result );
	CheckDigits( "sqrt(3)", &result, kRefSqrt3, 38 );

	/*
		The cube root, which is now Newton on y^3 = x rather than
		exp(ln(x)/3).

		Through the logarithm it had no answer for a negative -- the
		logarithm refuses one, and the refusal was discarded, so
		cbrt(-8) came back as nought with kNoError and a stray
		not-a-number flag. Nought itself came back flagged for the same
		reason. And the round trip through two series cost digits an
		iteration does not: 27 rooted to 2.999999999999999999999999999
		999999651.

		The argument is scaled into [1, base^3) by whole powers of
		base^3 first, which is exact and commutes with the root.
		Without it the square of a small estimate falls under the
		precision floor and takes the answer with it -- the cube root
		of 1e-40 came back with eleven correct digits instead of
		twenty-two.
	*/

	Set( a, "8" ); result.Clear();
	lz_root_3( &a, &result );
	CheckExact( "cbrt(8)", &result, "2" );

	Set( a, "27" ); result.Clear();
	lz_root_3( &a, &result );
	CheckExact( "cbrt(27)", &result, "3" );

	Set( a, "1000" ); result.Clear();
	lz_root_3( &a, &result );
	CheckExact( "cbrt(1000)", &result, "10" );

	// Negative arguments, which are perfectly well defined...
	Set( a, "-8" ); result.Clear();
	lz_root_3( &a, &result );
	CheckExact( "cbrt(-8)", &result, "-2" );

	Set( a, "-27" ); result.Clear();
	lz_root_3( &a, &result );
	CheckExact( "cbrt(-27)", &result, "-3" );

	// ...and nought, which must come back clean...
	Set( a, "0" ); result.Clear();
	lz_root_3( &a, &result );
	CheckExact( "cbrt(0)", &result, "0" );
	CheckTrue( "cbrt(0) carries no flag", !result.FlagIsOn( kfQNaN ) &&
										   !result.FlagIsOn( kfInfinity ), "" );

	Set( a, "2" ); result.Clear();
	lz_root_3( &a, &result );
	CheckDigits( "cbrt(2)", &result, kRefCbrtTwo, 37 );

	// The extremes, where the scaling earns its place...
	Set( a, ".0000000000000000000000000000000000000001" ); result.Clear();
	lz_root_3( &a, &result );
	// 36 significant digits, the rule for an answer under 0.1 -- see
	//	lzi_round_quotient_answer. It carried a 37th, the division's...
	CheckDigits( "cbrt(1e-40)", &result, kRefCbrtTiny, 36 );

	Set( a, "1000000000000000000000000000000000000000" ); result.Clear();
	lz_root_3( &a, &result );
	CheckDigits( "cbrt(1e39)", &result, kRefCbrtBig, 36 );

	/*
		The magnitude the band takes out, which has to be put back in
		places as well as in size.

		The three above are perfect cubes -- 8, 27, 1e39 -- and a
		perfect cube cannot show this: the iteration lands on the exact
		answer whatever width it works to. **That is why the shortfall
		lived here unnoticed.** The cube root of 1e40 held 37
		significant digits where 36 places at 1e13 has room for 50, and
		the loss was exactly the block count the band scaling takes
		out: 1 at 1e8, 5 at 1e20, 13 at 1e40.

		Requiring more than the old code could give is the whole point
		of these three; at 36 places the argument's own size is what
		decides how many digits that is.
	*/

	Set( a, "4000" ); result.Clear();
	lz_root_3( &a, &result );
	CheckDigits( "cbrt(4000)", &result, kRefCbrt4000, 38 );

	Set( a, "100000000000000000000" ); result.Clear();
	lz_root_3( &a, &result );
	CheckDigits( "cbrt(1e20)", &result, kRefCbrt1e20, 42 );

	Set( a, "10000000000000000000000000000000000000000" ); result.Clear();
	lz_root_3( &a, &result );
	CheckDigits( "cbrt(1e40)", &result, kRefCbrt1e40, 49 );

	/*
		**Correctly rounded, which the cube root was not**: its Newton
		iteration had no guard and settled a unit off the right answer
		in about one cube root in twenty-three, in every base -- 393 of
		9,020 random ones, base ten among them. It leaned on the extra
		places each division carried until 2026-09-30, and they were not
		enough. kCubeRootGuardDigits gives it its own. These are cube
		roots that guards of three digits or fewer, or the old code,
		round wrongly; each expected string is the exact root rounded to
		the answer's width, computed in Python and not copied from the
		engine.
	*/

	{
		struct { long lBase; long lPrecision; const char* szValue; const char* szExpect; } cbrts[] =
		{
			{ 8, 36, "753666042603336113050",
			  "7711437.300154114007652651075147053365577035" },
			{ 2, 36, "111111100110001011010011.10010100",
			  "11111111.011101011111110000101011111010010001" },
			{ 2, 80, "11100.111110011001",
			  "11.00010010010010010010010000100110010001000000110110000110000011001101100100101001" },
			{ 12, 36, "-820099E68274D255DD0D06350.53",
			  "-201ED9E77.D8117555375819647DDE7998690E272943D1" },
			{ 10, 36, "-0.000000000000000001112754714261638",
			  "-0.00000103625462275239245220859501894860505" },
			{ 10, 36, "-0.0000000000000078",
			  "-0.0000198319248268077472688430025350065727" },
			{ 12, 36, "0.0000000000000000000000005856706E02E744893226",
			  "0.00000000944810732738E4805226EE44E2E947353303" },
		};

		for ( auto& theCase : cbrts )
		{
			long		l_saved = glWorkingPrecision;
			CLZNumber	value, root;
			char		szLabel[ 64 ];

			glWorkingPrecision = theCase.lPrecision;
			value.SetBase( theCase.lBase, true );
			root.SetBase( theCase.lBase, true );
			value.Assign( theCase.szValue );
			lz_root_3( &value, &root );
			glWorkingPrecision = l_saved;

			::snprintf( szLabel, sizeof( szLabel ), "cbrt rounds right, base %ld at %ld", theCase.lBase, theCase.lPrecision );
			CheckRendering( szLabel, &root, theCase.szExpect );
		}
	}


	/*
		Square roots of exact squares, large and small.

		Newton only halves while it is far above the root, and the
		estimate used to be half the operand -- a very long way from
		the root of anything large. Sixty-three halvings from 5e39 land
		on 5.4e20, and that is what sqrt(1e40) returned in place of
		1e20: no flag, no code, a plausible number of the wrong
		magnitude. Anything past about 1e38 was wrong and nothing here
		looked. The estimate now starts within a factor of the base,
		whichever side, so these are exact.
	*/

	Set( a, "10000000000000000000000000000000000000000" ); result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt(1e40)", &result, "100000000000000000000" );

	Set( a, "1000000000000000000000000000000000000000000000000000000000000" );
	result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt(1e60)", &result, "1000000000000000000000000000000" );

	Set( a, "152415787532388367501905199875019052100" ); result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt of a 39-digit square", &result, "12345678901234567890" );

	// ...and below one, where the estimate starts under the root and
	//	the old code had the same climb to make in the other direction...
	Set( a, ".0000000000000000000000000000000000000001" ); result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt(1e-40)", &result, ".00000000000000000001" );

	Set( a, ".00000000000000000000000001" ); result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt(1e-26)", &result, ".0000000000001" );

	/*
		**Correctly rounded, which the square root was not**: its Newton
		iteration worked at the answer's places with no guard, and
		halved by multiplying by a half made with lz_reciprocal -- which
		in an odd base is .333... and not a half at all. 867 of 9,100
		random square roots were not correctly rounded at their places,
		509 of them by more than a unit, and 551 of 1,300 in base seven.
		These are square roots that guards of three digits or fewer, or
		the old code, get wrong -- the base-seven two by more than a
		unit; each expected string is the exact root rounded to the
		answer's places, computed in Python and not copied from the
		engine.
	*/

	{
		struct { long lBase; long lPrecision; const char* szValue; const char* szExpect; } sqrts[] =
		{
			{ 8, 36, "0.0000000000000000000000325652764",
			  "0.000000000005127035121213175214702425" },
			{ 2, 80, "11.1000",
			  "1.11011110111011101010000100010110100000111111010010010010000001010101010111001001" },
			{ 8, 80, "41406065567577113.312",
			  "202772257.17051415460101067245240547464103407653760344600575000146355244357630745347004651" },
			{ 7, 36, "140052566213205306232344440535322150.214251044000122535145",
			  "321426306134146050.236510554360441262103642424624350546" },
			{ 7, 36, "32135633506452154356343005656015.6046422134504040346",
			  "4550603346003030.606242250134362105260525335562562051" },
			{ 10, 36, "1304.43782468715399907216514448144177503",
			  "36.117001878438830336679760609240391709" },
			{ 2, 36, "10101111001011001110.00010110011110011011000111100101",
			  "1101001111.000100000111110000100010001001010001" },
		};

		for ( auto& theCase : sqrts )
		{
			long	l_saved = glWorkingPrecision;
			char	szLabel[ 64 ];

			// Before the numbers are made, which take the width they are made at...
			glWorkingPrecision = theCase.lPrecision;

			CLZNumber	value, root;

			value.SetBase( theCase.lBase, true );
			root.SetBase( theCase.lBase, true );
			value.Assign( theCase.szValue );
			lz_root_2( &value, &root );
			glWorkingPrecision = l_saved;

			::snprintf( szLabel, sizeof( szLabel ), "sqrt rounds right, base %ld at %ld", theCase.lBase, theCase.lPrecision );
			CheckRendering( szLabel, &root, theCase.szExpect );
		}
	}

	// The ones that were always right, held still...
	Set( a, "9" ); result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt(9)", &result, "3" );

	Set( a, "100" ); result.Clear();
	lz_root_2( &a, &result );
	CheckExact( "sqrt(100)", &result, "10" );

	Set( a, "2" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(2)", &result, kRefLn2, 37 );

	Set( a, "10" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(10)", &result, kRefLn10, 37 );

	/*
		Below a tenth, where the reduction used to give up.

		The series is in y = (x-1)/(x+1) and needs y small. Arguments
		above the base were divided down by e; nothing brought small
		ones up, so ln(.01) summed against y = -.98 and returned four
		correct digits, ln(.001) one. The suite had no logarithm of a
		small argument at all, which is why it never showed.

		Every case here is below the old reduction's reach, and the
		last of them is far enough below to need six multiplications.
	*/

	Set( a, ".05" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(.05)", &result, kRefLnTwentieth, 37 );

	Set( a, ".01" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(.01)", &result, kRefLnHundredth, 37 );

	Set( a, ".001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(.001)", &result, kRefLnThousandth, 37 );

	Set( a, ".0001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(.0001)", &result, kRefLnTenThousandth, 37 );

	Set( a, ".0000001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(.0000001)", &result, kRefLnTenMillionth, 38 );

	// Just inside the band on the small side, which never needed it...
	Set( a, ".1" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(.1)", &result, kRefLnTenth, 37 );

	// ...and a large argument, to hold the half that always worked...
	Set( a, "12345" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(12345)", &result, kRefLn12345, 37 );

	/*
		Far below the band, where the reduction used to give up
		entirely. Multiplying a small argument up by e rounds it to
		nothing once it is smaller than the working precision holds,
		and the loop then compared nought against a tenth for ever:
		ln(1e-36) answered, ln(1e-37) hung, and lz_root_3 hung with it
		for any argument that small. Below the band it now inverts
		instead -- ln(x) = -ln(1/x) -- and lets the reduction that was
		always right handle the large reciprocal.

		The climb was quietly inaccurate long before it hung, too:
		1e-20 came back with 18 correct digits, 1e-30 with 8, 1e-36
		with 2. None of that was visible because nothing here looked
		below 1e-7.
	*/

	Set( a, ".00000000000000000001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(1e-20)", &result, kRefLnE20, 39 );

	Set( a, ".000000000000000000000000000001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(1e-30)", &result, kRefLnE30, 39 );

	Set( a, ".0000000000000000000000000000000000000001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(1e-40)", &result, kRefLnE40, 39 );

	Set( a, ".0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001" );
	result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(1e-100)", &result, kRefLnE100, 39 );


	Set( a, "1" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "e", &result, kRefE, 37 );

	/*
		Within a whisker of one, where the answer IS the last digit.

		These used to round with RM_round_up at RD_digit_last, which
		takes the value to one place fewer than it has and rounds away
		from zero. Where a value had digits to spare that merely nudged
		the last one; where the perturbation was itself the last digit,
		it was discarded and replaced by one ten times larger.
		exp(1e-20) came back as 1 + 1e-19, exp(-1e-20) as exactly 1,
		and ln(1 + 1e-20) as 1e-19. All six exponential and logarithm
		entry points did it.
	*/

	Set( a, ".000000000000000001" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(1e-18)", &result, kRefExpE18, 37 );

	Set( a, ".00000000000000000001" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(1e-20)", &result, kRefExpE20, 41 );

	Set( a, "-.00000000000000000001" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(-1e-20)", &result, kRefExpNegE20, 40 );

	Set( a, "1.00000000000000000001" ); result.Clear();
	lz_log_e( &a, &result );
	CheckDigits( "ln(1 + 1e-20)", &result, kRefLnNear1, 20 );


	// The only caller of the kLogNBaseE constants, and the only base-b
	//	operation the suite reaches. 10^.5 is the square root of ten...
	Set( a, ".5" ); result.Clear();
	lz_exp_b( &a, &result );
	CheckDigits( "10^.5 (base-b exponential)", &result,
				 "3.1622776601683793319988935444327185337195551393252", 37 );

	Set( a, "1" ); result.Clear();
	lz_sine( &a, &result );
	CheckDigits( "sin(1)", &result, kRefSin1, 36 );

	Set( a, "1" ); result.Clear();
	lz_cosine( &a, &result );
	CheckDigits( "cos(1)", &result, kRefCos1, 36 );

	Set( a, "1" ); result.Clear();
	lz_tangent( &a, &result );
	CheckDigits( "tan(1)", &result, kRefTan1, 38 );

	Set( a, "1" ); result.Clear();
	lz_hypersine( &a, &result );
	CheckDigits( "sinh(1)", &result, kRefSinh1, 37 );

	Set( a, "1" ); result.Clear();
	lz_hypercosine( &a, &result );
	CheckDigits( "cosh(1)", &result, kRefCosh1, 37 );

	Set( a, "1" ); result.Clear();
	lz_hypertangent( &a, &result );
	CheckDigits( "tanh(1)", &result, kRefTanh1, 36 );

	/*
		Hyperbolics away from whole numbers, and past 2*pi.

		These took the argument modulo 2*pi first -- a reduction that
		belongs to the circular functions and means nothing for one
		that is not periodic -- and then raised e to it with
		lz_raise_e, which takes an integer exponent and then dropped
		any fraction. So they were right for whole-numbered arguments
		below 2*pi and nowhere else: sinh(.5) answered nought,
		sinh(1.1) answered sinh(1), sinh(10) answered sinh(3). The
		suite tested them at 1, which is one of the seven values that
		worked.
	*/

	Set( a, ".5" ); result.Clear();
	lz_hypersine( &a, &result );
	CheckDigits( "sinh(.5)", &result, kRefSinhHalf, 36 );

	Set( a, "1.1" ); result.Clear();
	lz_hypersine( &a, &result );
	CheckDigits( "sinh(1.1)", &result, kRefSinhOnePt1, 37 );

	Set( a, "-.5" ); result.Clear();
	lz_hypersine( &a, &result );
	CheckDigits( "sinh(-.5)", &result, kRefSinhMinusHalf, 36 );

	/*
		These read a digit higher until 2026-09-25: the answer is a
		quotient and carried the division's extra place unrounded, which
		counts as a digit here. They are correctly rounded at 36 places
		now -- see lzi_round_quotient_answer.
	*/

	// Past 2*pi, where the old reduction wrapped it round...
	Set( a, "10" ); result.Clear();
	lz_hypersine( &a, &result );
	CheckDigits( "sinh(10)", &result, kRefSinhTen, 41 );

	Set( a, ".5" ); result.Clear();
	lz_hypercosine( &a, &result );
	CheckDigits( "cosh(.5)", &result, kRefCoshHalf, 37 );

	Set( a, "10" ); result.Clear();
	lz_hypercosine( &a, &result );
	CheckDigits( "cosh(10)", &result, kRefCoshTen, 41 );

	Set( a, ".5" ); result.Clear();
	lz_hypertangent( &a, &result );
	CheckDigits( "tanh(.5)", &result, kRefTanhHalf, 36 );

	Set( a, "10" ); result.Clear();
	lz_hypertangent( &a, &result );
	CheckDigits( "tanh(10)", &result, kRefTanhTen, 37 );


	/*
		The thresholds below are what the halve-and-double achieves now
		that it sums against guard digits. They were 31, 35, 32 and 31,
		which the unguarded code met with two digits to spare on every
		argument that halves at all - so they gated nothing that
		mattered. Anything here dropping by a digit means the guard
		places have stopped arriving.
	*/

	Set( a, ".5" ); result.Clear();
	lz_arcsine( &a, &result );
	CheckDigits( "arcsin(.5)", &result, kRefAsinHalf, 37 );

	// Across the domain, and through every count of half-angle steps...
	Set( a, ".05" ); result.Clear();
	lz_arcsine( &a, &result );
	CheckDigits( "arcsin(.05)", &result, kRefAsinTwentieth, 36 );

	// The deepest reduction, and what used to be the worst answer...
	Set( a, ".75" ); result.Clear();
	lz_arcsine( &a, &result );
	CheckDigits( "arcsin(.75)", &result, kRefAsinThreeQtrs, 37 );

	Set( a, ".9" ); result.Clear();
	lz_arcsine( &a, &result );
	CheckDigits( "arcsin(.9)", &result, kRefAsinNinth, 37 );

	Set( a, "1" ); result.Clear();
	lz_arcsine( &a, &result );
	CheckDigits( "arcsin(1) = pi/2", &result, kRefPiOverTwo, 38 );

	Set( a, "-.5" ); result.Clear();
	lz_arcsine( &a, &result );
	CheckDigits( "arcsin(-.5)", &result, kRefAsinMinusHalf, 37 );

	/*
		The arctangent, which has no domain to refuse but does have a
		reflection at |x| = 1: past there it answers
		sign(x)*pi/2 - arctan(1/x), because the arcsine it routes
		through is worst conditioned exactly where the unreflected
		quotient would put it. Both sides are covered here, and a
		million to hold the far end -- 33 digits unreflected, 37 with.
	*/

	Set( a, ".1" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckDigits( "arctan(.1)", &result, kRefAtanTenth, 35 );

	Set( a, ".5" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckDigits( "arctan(.5)", &result, kRefAtanHalf, 37 );

	// The boundary itself, which does not reflect...
	Set( a, "1" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckDigits( "arctan(1) = pi/4", &result, kRefAtanOne, 37 );

	// ...and the first argument that does...
	Set( a, "2" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckDigits( "arctan(2)", &result, kRefAtanTwo, 38 );

	Set( a, "-2" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckDigits( "arctan(-2)", &result, kRefAtanMinusTwo, 38 );

	Set( a, "1000000" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckDigits( "arctan(1000000)", &result, kRefAtanMillion, 37 );

	// arctan(0) is a clean nought...
	Set( a, "0" ); result.Clear();
	lz_arctangent( &a, &result );
	CheckExact( "arctan(0)", &result, "0" );

	/*
		The arccosine, which is pi/2 less the arcsine and nothing
		cleverer. Its accuracy falls away as x approaches 1 and that is
		the representation rather than the algorithm: the answer goes
		to nought, precision counts places after the radix mark, and a
		small answer therefore has few significant digits whatever
		produced it. The half-angle form usually recommended against
		the cancellation was built and scored, and is never better
		here. So arccos(.9999999) is gated at what the floor allows,
		not at what the others manage.
	*/

	Set( a, "-1" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckDigits( "arccos(-1) = pi", &result, kRefAcosMinusOne, 37 );

	Set( a, "-.5" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckDigits( "arccos(-.5)", &result, kRefAcosMinusHalf, 37 );

	Set( a, ".5" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckDigits( "arccos(.5)", &result, kRefAcosHalf, 38 );

	Set( a, ".9" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckDigits( "arccos(.9)", &result, kRefAcosNineTenths, 36 );

	// Near the top of the domain, where the answer is small...
	Set( a, ".9999999" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckDigits( "arccos(.9999999)", &result, kRefAcosNearOne, 33 );

	// The two ends, exactly...
	Set( a, "1" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckExact( "arccos(1)", &result, "0" );

	Set( a, "0" ); result.Clear();
	lz_arccosine( &a, &result );
	CheckDigits( "arccos(0) = pi/2", &result, kRefPiOverTwo, 38 );

	/*
		The three inverse hyperbolics, each a logarithm of something
		built with a square root.

		arcsinh is odd and is folded onto its magnitude before the
		work, because x + sqrt(x^2+1) cancels for a negative x -- at
		-1000 the two terms agree to six digits and the argument to ln
		collapses from a thousand to five ten-thousandths. Folding
		keeps the argument large and the answer whole.

		arctanh was the casualty of two faults fixed just before this:
		it divides (1+x) by (1-x), which for a positive x divides by
		something below one, and it takes the logarithm of a quotient
		that grows without bound. It scored 19 digits at .9 and none
		at all at .99 until both were dealt with.
	*/

	Set( a, ".5" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(.5)", &result, kRefAsinhHalf, 36 );

	Set( a, "1" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(1)", &result, kRefAsinhOne, 37 );

	Set( a, "2" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(2)", &result, kRefAsinhTwo, 37 );

	Set( a, "50" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(50)", &result, kRefAsinhFifty, 39 );

	Set( a, ".001" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(.001)", &result, kRefAsinhTiny, 33 );

	// The odd fold: a negative must match its positive, mirrored...
	Set( a, "-2" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(-2)", &result, kRefAsinhMinus2, 37 );

	Set( a, "-1000" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckDigits( "arcsinh(-1000)", &result, kRefAsinhMinus1k, 37 );

	Set( a, "0" ); result.Clear();
	lz_archypersine( &a, &result );
	CheckExact( "arcsinh(0)", &result, "0" );

	Set( a, "1.001" ); result.Clear();
	lz_archypercosine( &a, &result );
	CheckDigits( "arccosh(1.001)", &result, kRefAcoshMilli, 35 );

	Set( a, "2" ); result.Clear();
	lz_archypercosine( &a, &result );
	CheckDigits( "arccosh(2)", &result, kRefAcoshTwo, 38 );

	Set( a, "100" ); result.Clear();
	lz_archypercosine( &a, &result );
	CheckDigits( "arccosh(100)", &result, kRefAcoshHundred, 38 );

	Set( a, "1" ); result.Clear();
	lz_archypercosine( &a, &result );
	CheckExact( "arccosh(1)", &result, "0" );

	// Below one there is no angle to find...
	Set( a, ".5" ); result.Clear();
	if ( lz_archypercosine( &a, &result ) != kNoError )
	{
		g_passed++;
		printf( "  pass  %-34s reports failure\n", "arccosh(.5) out of domain" );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s accepted an impossible argument\n", "arccosh(.5)" );
	}

	Set( a, ".5" ); result.Clear();
	lz_archypertangent( &a, &result );
	CheckDigits( "arctanh(.5)", &result, kRefAtanhHalf, 36 );

	Set( a, "-.5" ); result.Clear();
	lz_archypertangent( &a, &result );
	CheckDigits( "arctanh(-.5)", &result, kRefAtanhMinusH, 36 );

	// Toward the end of the domain, which used to return nothing right...
	Set( a, ".9" ); result.Clear();
	lz_archypertangent( &a, &result );
	CheckDigits( "arctanh(.9)", &result, kRefAtanhNineT, 37 );

	Set( a, ".99" ); result.Clear();
	lz_archypertangent( &a, &result );
	CheckDigits( "arctanh(.99)", &result, kRefAtanh99, 38 );

	Set( a, ".9999999" ); result.Clear();
	lz_archypertangent( &a, &result );
	CheckDigits( "arctanh(.9999999)", &result, kRefAtanhNearOne, 37 );

	Set( a, "0" ); result.Clear();
	lz_archypertangent( &a, &result );
	CheckExact( "arctanh(0)", &result, "0" );

	// The ends do not answer either, though they answer differently
	//	from the outside: at exactly one the quotient divides by nought
	//	and the limit is an infinity, which is kUndefinedResult with
	//	kfInfinity rather than a flat refusal. Pinned in detail among
	//	the failures that report; this only holds the line that the
	//	edge is not quietly accepted...
	Set( a, "1" ); result.Clear();
	if ( lz_archypertangent( &a, &result ) != kNoError )
	{
		g_passed++;
		printf( "  pass  %-34s reports\n", "arctanh(1) at the edge" );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s accepted the domain's edge\n", "arctanh(1)" );
	}

	// Outside [-1, 1] there is no angle to find, for either of them...
	Set( a, "2" ); result.Clear();
	if ( lz_arccosine( &a, &result ) != kNoError )
	{
		g_passed++;
		printf( "  pass  %-34s reports failure\n", "arccos(2) out of domain" );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s accepted an impossible argument\n", "arccos(2)" );
	}

	Set( a, "2" ); result.Clear();
	if ( lz_arcsine( &a, &result ) != kNoError )
	{
		g_passed++;
		printf( "  pass  %-34s reports failure\n", "arcsin(2) out of domain" );
	}
	else
	{
		g_failed++;
		printf( "  FAIL  %-34s accepted an impossible argument\n", "arcsin(2)" );
	}
}

/////////////////////////////////////////////////////////////////
//	Argument magnitude - the sine series uses a fixed 16 terms, so
//	its truncation error grows with the size of the argument. Cosine
//	is computed as sin(x + pi/2), so it inherits the loss one
//	quadrant early.
/////////////////////////////////////////////////////////////////

static void TestArgumentMagnitude()
{
	CLZNumber	a, result;

	Banner( "Sine accuracy against argument magnitude" );

	struct { const char* szArg; const char* szRef; int iRequired; } cases[] =
	{
		{ "1", kRefSin1,  33 },
		{ "2", kRefSin2,  34 },
		{ "3", kRefSin3,  34 },
		{ "4", kRefSin4,  33 },
		{ "6", kRefSin6,  35 },
	};

	for ( auto& theCase : cases )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "sin(%s)", theCase.szArg );

		Set( a, theCase.szArg );
		result.Clear();
		lz_sine( &a, &result );
		CheckDigits( szLabel, &result, theCase.szRef, theCase.iRequired );
	}

	/*
		What kTwoPi's width decides, and the reason it carries 120
		places where its three neighbours carry 50.

		The reduction takes the argument modulo 2*pi, which multiplies
		the constant by a quotient the size of the argument -- so every
		digit of magnitude costs a digit of the constant. At 50 places
		these three held 30, 19 and 10 digits, falling exactly one for
		one with the argument. **Narrowing the constant fails this at
		1e40 by twenty-six digits**, which is what the test is for; the
		sweep records the same thirty cases, but records them as a
		table of numbers rather than as a stated expectation.
	*/

	Banner( "The reduction at scale" );

	struct { const char* szArg; const char* szSin; const char* szCos; int iRequired; } wideCases[] =
	{
		{ "100000000000000000000",
		  kRefSin1e20, kRefCos1e20, 35 },
		{ "1000000000000000000000000000000",
		  kRefSin1e30, kRefCos1e30, 34 },
		{ "10000000000000000000000000000000000000000",
		  kRefSin1e40, kRefCos1e40, 35 },
	};

	for ( auto& theCase : wideCases )
	{
		char szLabel[ 72 ];

		::snprintf( szLabel, 72, "sin(1e%d)", static_cast<int>( ::strlen( theCase.szArg ) - 1 ) );
		Set( a, theCase.szArg );
		result.Clear();
		lz_sine( &a, &result );
		CheckDigits( szLabel, &result, theCase.szSin, theCase.iRequired );

		::snprintf( szLabel, 72, "cos(1e%d)", static_cast<int>( ::strlen( theCase.szArg ) - 1 ) );
		Set( a, theCase.szArg );
		result.Clear();
		lz_cosine( &a, &result );
		CheckDigits( szLabel, &result, theCase.szCos, theCase.iRequired );
	}

	/*
		The series past the table, which is the whole of what a caller
		who raises the precision used to be denied.

		gRecipFact stops at x^32/32!, and that truncation -- not the
		width of anything -- held the sine to between 42 and 50 digits
		however high the precision went. Past a threshold
		lzi_trig_series now computes every term from the recurrence
		instead, which has no length limit and no table to be bounded
		by, and the pi family behind it went to 120 places.

		Both of those had to happen, and the order in which they were
		found is the lesson: widening the pi family alone measured as
		no gain at all, because the truncation was binding in front of
		it. Removing the truncation alone took every argument to 50 and
		no further, because then the pi family was.

		A hundred rather than two hundred, for the same reason ln(10)
		is tested at a hundred: past about 120 the constants are the
		ceiling again, which is a true constant-width limit and does not
		belong in a gate.
	*/

	Banner( "The series past the table" );

	{
		long	l_saved = glWorkingPrecision;

		struct { const char* szArg; bool bCosine; const char* szRef; } deep[] =
		{
			{ "1", false, kRefSin1Wide },
			{ "1", true,  kRefCos1Wide },
			{ "4", false, kRefSin4Wide },
		};

		for ( auto& theCase : deep )
		{
			CLZNumber	wide, wideResult;
			char		szLabel[ 72 ];

			glWorkingPrecision = 100;

			wide.SetPrecision( 100 );
			wide.Assign( theCase.szArg );
			wideResult.SetPrecision( 100 );

			if ( theCase.bCosine )
				lz_cosine( &wide, &wideResult );
			else
				lz_sine( &wide, &wideResult );

			glWorkingPrecision = l_saved;

			::snprintf( szLabel, 72, "%s(%s) at precision 100",
						theCase.bCosine ? "cos" : "sin", theCase.szArg );

			CheckDigits( szLabel, &wideResult, theCase.szRef, 99 );
		}
	}

	/*
		And an infinity at that same precision, which is the one case
		the computed tail must not touch.

		An infinity is held as a zero carrying a flag. The tail stops
		when a term has faded to nothing, reads that zero, and leaves
		-- taking the flag with it, so sec(infinity) came back as 1.
		The tail is skipped for a flagged argument because of it.

		This has to be tested *at a raised precision*: at the default
		the tail does not run at all, so the guard is never reached and
		removing it breaks nothing. A mutation that dropped it was
		caught by no test until this one existed.
	*/

	{
		long	l_saved = glWorkingPrecision;

		struct { const char* szName; bool bCosine; } both[] =
		{
			{ "sin", false }, { "cos", true },
		};

		for ( auto& theCase : both )
		{
			CLZNumber	wide, wideResult;
			char		szLabel[ 72 ];

			glWorkingPrecision = 100;

			wide.SetPrecision( 100 );
			wide.SetFlag( kfInfinity );
			wideResult.SetPrecision( 100 );

			if ( theCase.bCosine )
				lz_cosine( &wide, &wideResult );
			else
				lz_sine( &wide, &wideResult );

			glWorkingPrecision = l_saved;

			::snprintf( szLabel, 72, "%s(infinity) at precision 100", theCase.szName );
			CheckTrue( szLabel, wideResult.FlagIsOn( kfQNaN ), "" );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	Identities - these need no reference values, and they expose
//	error that accumulates inside the iterative algorithms
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
//	The reciprocal trigonometric family. The cotangent, secant and
//	cosecant were in the engine all along with no test, no sweep
//	case and no mention in the notes; the other nine are new. All
//	twelve are covered here, then again in the magnitude sweep.
/////////////////////////////////////////////////////////////////

// Defined with the flag-propagation tests further down, and used here
//	as well: an infinity is stored as zero with a flag, so it cannot
//	simply be assigned...
static void MakeInfinite( CLZNumber& theNumber );
static void MakeNotANumber( CLZNumber& theNumber );

static void TestReciprocalTrig()
{
	CLZNumber	a, x, y, result;

	Banner( "Reciprocal trigonometry" );

	/*
		The three circular reciprocals, each a quotient of the sine and
		cosine it is built on and so no better than they are.
	*/

	Set( a, "1" );   result.Clear(); lz_cotangent( &a, &result );
	CheckDigits( "cot(1)", &result, kRefCot1, 36 );

	Set( a, ".5" );  result.Clear(); lz_cotangent( &a, &result );
	CheckDigits( "cot(.5)", &result, kRefCotHalf, 37 );

	Set( a, "2" );   result.Clear(); lz_cotangent( &a, &result );
	CheckDigits( "cot(2)", &result, kRefCotTwo, 36 );

	Set( a, "-1" );  result.Clear(); lz_cotangent( &a, &result );
	CheckDigits( "cot(-1)", &result, kRefCotMinusOne, 36 );

	Set( a, "1" );   result.Clear(); lz_secant( &a, &result );
	CheckDigits( "sec(1)", &result, kRefSec1, 37 );

	Set( a, ".5" );  result.Clear(); lz_secant( &a, &result );
	CheckDigits( "sec(.5)", &result, kRefSecHalf, 37 );

	Set( a, "2" );   result.Clear(); lz_secant( &a, &result );
	CheckDigits( "sec(2)", &result, kRefSecTwo, 37 );

	// sec(0) is one exactly, the cosine there being one exactly...
	Set( a, "0" );   result.Clear(); lz_secant( &a, &result );
	CheckExact( "sec(0)", &result, "1" );

	Set( a, "1" );   result.Clear(); lz_cosecant( &a, &result );
	CheckDigits( "csc(1)", &result, kRefCsc1, 37 );

	Set( a, ".5" );  result.Clear(); lz_cosecant( &a, &result );
	CheckDigits( "csc(.5)", &result, kRefCscHalf, 38 );

	Set( a, "2" );   result.Clear(); lz_cosecant( &a, &result );
	CheckDigits( "csc(2)", &result, kRefCscTwo, 37 );

	/*
		The three hyperbolic reciprocals. coth is taken from cosh over
		sinh directly rather than as one over tanh, which would divide
		and round twice for the same answer.
	*/

	Set( a, "1" );   result.Clear(); lz_hypercotangent( &a, &result );
	CheckDigits( "coth(1)", &result, kRefCoth1, 37 );

	Set( a, ".5" );  result.Clear(); lz_hypercotangent( &a, &result );
	CheckDigits( "coth(.5)", &result, kRefCothHalf, 36 );

	// coth(10) and the three sechs read a digit higher until 2026-09-25,
	//	for the extra place a quotient carried -- see sinh(10) above...

	// Past 2*pi, where the old hyperbolics used to reduce and be wrong...
	Set( a, "10" );  result.Clear(); lz_hypercotangent( &a, &result );
	CheckDigits( "coth(10)", &result, kRefCothTen, 37 );

	Set( a, "-1" );  result.Clear(); lz_hypercotangent( &a, &result );
	CheckDigits( "coth(-1)", &result, kRefCothMinusOne, 37 );

	Set( a, "1" );   result.Clear(); lz_hypersecant( &a, &result );
	CheckDigits( "sech(1)", &result, kRefSech1, 36 );

	Set( a, ".5" );  result.Clear(); lz_hypersecant( &a, &result );
	CheckDigits( "sech(.5)", &result, kRefSechHalf, 36 );

	Set( a, "2" );   result.Clear(); lz_hypersecant( &a, &result );
	CheckDigits( "sech(2)", &result, kRefSechTwo, 36 );

	// sech has no pole anywhere, cosh never falling below one...
	Set( a, "0" );   result.Clear(); lz_hypersecant( &a, &result );
	CheckExact( "sech(0)", &result, "1" );

	Set( a, "1" );   result.Clear(); lz_hypercosecant( &a, &result );
	CheckDigits( "csch(1)", &result, kRefCsch1, 36 );

	Set( a, ".5" );  result.Clear(); lz_hypercosecant( &a, &result );
	CheckDigits( "csch(.5)", &result, kRefCschHalf, 36 );

	Set( a, "-1" );  result.Clear(); lz_hypercosecant( &a, &result );
	CheckDigits( "csch(-1)", &result, kRefCschMinusOne, 36 );

	/*
		The arccotangent, arctan(1/x) and pi/2 at nought: range (-pi/2,
		pi/2], odd, as Elzed answers -- Rich's call, 2026-10-23. Until
		then it was pi/2 - arctan(x), the continuous convention, range
		(0, pi), and the two differ only on the negatives.
	*/

	Set( a, "1" );   result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(1) = pi/4", &result, kRefAcot1, 37 );

	Set( a, ".5" );  result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(.5)", &result, kRefAcotHalf, 38 );

	Set( a, "2" );   result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(2)", &result, kRefAcotTwo, 37 );

	// The negatives, which is where the two conventions part company...
	Set( a, "-1" );  result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(-1) = -pi/4", &result, kRefAcotMinusOne, 37 );

	Set( a, "-2" );  result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(-2)", &result, kRefAcotMinusTwo, 37 );

	// arccot(0) is pi/2, the positive side of the jump...
	Set( a, "0" );   result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(0) = pi/2", &result, kRefPiOverTwo, 38 );

	// Far out, where the answer is small and the floor is what limits it...
	Set( a, "1000000" ); result.Clear(); lz_arccotangent( &a, &result );
	CheckDigits( "arccot(1000000)", &result, kRefAcotMillion, 30 );

	/*
		The arcsecant and arccosecant, each its co-function's inverse
		taken on the reciprocal, and each guarded: just outside |x| = 1
		the reciprocal is just inside it, where the inverse's slope
		climbs without bound. See kArcCoFunctionGuardDigits.
	*/

	Set( a, "2" );   result.Clear(); lz_arcsecant( &a, &result );
	CheckDigits( "arcsec(2) = pi/3", &result, kRefAsecTwo, 38 );

	Set( a, "-2" );  result.Clear(); lz_arcsecant( &a, &result );
	CheckDigits( "arcsec(-2) = 2pi/3", &result, kRefAsecMinusTwo, 37 );

	// Close to the edge of the domain, where the guard digits earn their keep...
	Set( a, "1.1" ); result.Clear(); lz_arcsecant( &a, &result );
	CheckDigits( "arcsec(1.1)", &result, kRefAsecElevenTenths, 36 );

	Set( a, "1000000" ); result.Clear(); lz_arcsecant( &a, &result );
	CheckDigits( "arcsec(1000000)", &result, kRefAsecMillion, 37 );

	// The two ends of the domain, exactly...
	Set( a, "1" );   result.Clear(); lz_arcsecant( &a, &result );
	CheckExact( "arcsec(1)", &result, "0" );

	Set( a, "-1" );  result.Clear(); lz_arcsecant( &a, &result );
	CheckDigits( "arcsec(-1) = pi", &result, kRefAcosMinusOne, 37 );

	Set( a, "2" );   result.Clear(); lz_arccosecant( &a, &result );
	CheckDigits( "arccsc(2) = pi/6", &result, kRefAcscTwo, 37 );

	Set( a, "-2" );  result.Clear(); lz_arccosecant( &a, &result );
	CheckDigits( "arccsc(-2)", &result, kRefAcscMinusTwo, 37 );

	Set( a, "1.1" ); result.Clear(); lz_arccosecant( &a, &result );
	CheckDigits( "arccsc(1.1)", &result, kRefAcscElevenTenths, 37 );

	Set( a, "1000000" ); result.Clear(); lz_arccosecant( &a, &result );
	CheckDigits( "arccsc(1000000)", &result, kRefAcscMillion, 31 );

	Set( a, "1" );   result.Clear(); lz_arccosecant( &a, &result );
	CheckDigits( "arccsc(1) = pi/2", &result, kRefPiOverTwo, 38 );

	/*
		The three inverse hyperbolic reciprocals.
	*/

	Set( a, "2" );   result.Clear(); lz_archypercotangent( &a, &result );
	CheckDigits( "arccoth(2)", &result, kRefAcothTwo, 36 );

	Set( a, "-2" );  result.Clear(); lz_archypercotangent( &a, &result );
	CheckDigits( "arccoth(-2)", &result, kRefAcothMinusTwo, 36 );

	Set( a, "1.1" ); result.Clear(); lz_archypercotangent( &a, &result );
	CheckDigits( "arccoth(1.1)", &result, kRefAcothElevenTenths, 37 );

	Set( a, "1000000" ); result.Clear(); lz_archypercotangent( &a, &result );
	CheckDigits( "arccoth(1000000)", &result, kRefAcothMillion, 31 );

	Set( a, ".5" );  result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.5)", &result, kRefAsechHalf, 38 );

	Set( a, ".1" );  result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.1)", &result, kRefAsechTenth, 38 );

	Set( a, ".9" );  result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.9)", &result, kRefAsechNineTenths, 36 );

	/*
		Both sides of the threshold at which the hyperbolic arcsecant
		stops going through the arccosine and sums its own series
		instead -- kArcHyperSecantSeriesCeiling, at a hundredth. The
		pair either side of it, .011 and .01, take different routes to
		answers that agree to the last place but one.

		arcsech(.01) sat three digits below its floor when the series
		went in, and so did arcsech(.2) and arcsech(.0005) on the
		other route -- all three spending their accuracy in the one
		call to lz_log_e that every route here shares. That was the
		logarithm's reduction band, since fixed, and all three sit on
		the floor now.
	*/

	Set( a, ".011" ); result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.011), arccosine route", &result, kRefAsechElevenThousandths, 37 );

	Set( a, ".01" );  result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.01), series route", &result, kRefAsechHundredth, 38 );

	Set( a, ".001" ); result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.001)", &result, kRefAsechThousandth, 37 );

	Set( a, ".00001" ); result.Clear(); lz_archypersecant( &a, &result );
	CheckDigits( "arcsech(.00001)", &result, kRefAsechHundredThousandth, 38 );

	// The top of the domain, where the answer is nought exactly...
	Set( a, "1" );   result.Clear(); lz_archypersecant( &a, &result );
	CheckExact( "arcsech(1)", &result, "0" );

	Set( a, "1" );   result.Clear(); lz_archypercosecant( &a, &result );
	CheckDigits( "arccsch(1)", &result, kRefAcsch1, 37 );

	Set( a, "-1" );  result.Clear(); lz_archypercosecant( &a, &result );
	CheckDigits( "arccsch(-1)", &result, kRefAcschMinusOne, 37 );

	Set( a, ".5" );  result.Clear(); lz_archypercosecant( &a, &result );
	CheckDigits( "arccsch(.5)", &result, kRefAcschHalf, 37 );

	Set( a, "2" );   result.Clear(); lz_archypercosecant( &a, &result );
	CheckDigits( "arccsch(2)", &result, kRefAcschTwo, 36 );

	Set( a, "1000000" ); result.Clear(); lz_archypercosecant( &a, &result );
	CheckDigits( "arccsch(1000000)", &result, kRefAcschMillion, 30 );

	/*
		Identities. Each of the six reciprocals against the function it
		is the reciprocal of, and the two Pythagorean relations that
		tie the pairs together. These would catch a function wired to
		the wrong partner, which the value tests above would not if the
		reference had been generated the same wrong way.
	*/

	Banner( "Reciprocal identities" );

	struct { const char* szLabel;
			 long (*pForward)( CLZNumber*, CLZNumber* );
			 long (*pReciprocal)( CLZNumber*, CLZNumber* );
			 int iDigits; } reciprocals[] =
	{
		{ "tan(1) * cot(1)",   lz_tangent,       lz_cotangent,        33 },
		{ "cos(1) * sec(1)",   lz_cosine,        lz_secant,           34 },
		{ "sin(1) * csc(1)",   lz_sine,          lz_cosecant,         33 },
		{ "tanh(1) * coth(1)", lz_hypertangent,  lz_hypercotangent,   33 },
		{ "cosh(1) * sech(1)", lz_hypercosine,   lz_hypersecant,      33 },
		{ "sinh(1) * csch(1)", lz_hypersine,     lz_hypercosecant,    33 },
	};

	for ( auto& thePair : reciprocals )
	{
		Set( a, "1" );
		x.Clear(); y.Clear(); result.Clear();

		thePair.pForward( &a, &x );
		thePair.pReciprocal( &a, &y );
		lz_multiply( &x, &y, &result );

		CheckDigits( thePair.szLabel, &result, "1", thePair.iDigits );
	}

	// sec^2 - tan^2 = 1...
	Set( a, "1" );
	x.Clear(); y.Clear(); result.Clear();
	lz_secant( &a, &x );
	lz_tangent( &a, &y );
	lz_multiply( &x, &x );
	lz_multiply( &y, &y );
	lz_subtract( &x, &y, &result );
	CheckDigits( "sec(1)^2 - tan(1)^2", &result, "1", 36 );

	// csc^2 - cot^2 = 1...
	Set( a, "1" );
	x.Clear(); y.Clear(); result.Clear();
	lz_cosecant( &a, &x );
	lz_cotangent( &a, &y );
	lz_multiply( &x, &x );
	lz_multiply( &y, &y );
	lz_subtract( &x, &y, &result );
	CheckDigits( "csc(1)^2 - cot(1)^2", &result, "1", 36 );

	// tanh^2 + sech^2 = 1...
	Set( a, "1" );
	x.Clear(); y.Clear(); result.Clear();
	lz_hypertangent( &a, &x );
	lz_hypersecant( &a, &y );
	lz_multiply( &x, &x );
	lz_multiply( &y, &y );
	lz_add( &x, &y, &result );
	CheckDigits( "tanh(1)^2 + sech(1)^2", &result, "1", 33 );

	// coth^2 - csch^2 = 1...
	Set( a, "1" );
	x.Clear(); y.Clear(); result.Clear();
	lz_hypercotangent( &a, &x );
	lz_hypercosecant( &a, &y );
	lz_multiply( &x, &x );
	lz_multiply( &y, &y );
	lz_subtract( &x, &y, &result );
	CheckDigits( "coth(1)^2 - csch(1)^2", &result, "1", 36 );

	// arctan(x) + arccot(x) = pi/2 for x past nought...
	Set( a, "2" );
	x.Clear(); y.Clear(); result.Clear();
	lz_arctangent( &a, &x );
	lz_arccotangent( &a, &y );
	lz_add( &x, &y, &result );
	CheckDigits( "arctan(2) + arccot(2) = pi/2", &result, kRefPiOverTwo, 38 );

	// ...and -pi/2 below it, which pins the convention: the continuous
	//	one answered pi/2 here too...
	Set( a, "-2" );
	x.Clear(); y.Clear(); result.Clear();
	lz_arctangent( &a, &x );
	lz_arccotangent( &a, &y );
	lz_add( &x, &y, &result );
	CheckDigits( "arctan(-2) + arccot(-2) = -pi/2", &result, kRefMinusPiOverTwo, 38 );

	// arcsec(x) + arccsc(x) = pi/2...
	Set( a, "3" );
	x.Clear(); y.Clear(); result.Clear();
	lz_arcsecant( &a, &x );
	lz_arccosecant( &a, &y );
	lz_add( &x, &y, &result );
	CheckDigits( "arcsec(3) + arccsc(3) = pi/2", &result, kRefPiOverTwo, 38 );

	/*
		Round trips. Every inverse taken back through its forward
		function, which is the check that catches a wrong branch or a
		wrong sign where a value test cannot.
	*/

	struct { const char* szLabel; const char* szArgument;
			 long (*pInverse)( CLZNumber*, CLZNumber* );
			 long (*pForward)( CLZNumber*, CLZNumber* );
			 int iDigits; } trips[] =
	{
		{ "cot(arccot(3))",     "3",   lz_arccotangent,       lz_cotangent,      35 },
		{ "cot(arccot(-3))",    "-3",  lz_arccotangent,       lz_cotangent,      35 },
		{ "sec(arcsec(4))",     "4",   lz_arcsecant,          lz_secant,         35 },
		{ "sec(arcsec(-4))",    "-4",  lz_arcsecant,          lz_secant,         35 },
		{ "csc(arccsc(4))",     "4",   lz_arccosecant,        lz_cosecant,       35 },
		{ "csc(arccsc(-4))",    "-4",  lz_arccosecant,        lz_cosecant,       35 },
		{ "coth(arccoth(3))",   "3",   lz_archypercotangent,  lz_hypercotangent, 36 },
		{ "coth(arccoth(-3))",  "-3",  lz_archypercotangent,  lz_hypercotangent, 36 },
		{ "sech(arcsech(.25))", ".25", lz_archypersecant,     lz_hypersecant,    36 },
		{ "csch(arccsch(3))",   "3",   lz_archypercosecant,   lz_hypercosecant,  36 },
		{ "csch(arccsch(-3))",  "-3",  lz_archypercosecant,   lz_hypercosecant,  36 },
	};

	for ( auto& theTrip : trips )
	{
		Set( a, theTrip.szArgument );
		x.Clear(); result.Clear();

		theTrip.pInverse( &a, &x );
		theTrip.pForward( &x, &result );

		CheckDigits( theTrip.szLabel, &result, theTrip.szArgument, theTrip.iDigits );
	}

	/*
		Poles and domains. Every one of these once had no test at all,
		and four of the twelve functions have a pole where the answer
		is an infinity rather than a refusal.
	*/

	Banner( "Reciprocal poles and domains" );

	struct { const char* szLabel; const char* szArgument;
			 long (*pFunction)( CLZNumber*, CLZNumber* );
			 long lExpect; unsigned long lFlag; } edges[] =
	{
		// The four poles, reported by the division underneath...
		{ "cot(0)",     "0",   lz_cotangent,         kDivideByZero,    kfInfinity },
		{ "csc(0)",     "0",   lz_cosecant,          kDivideByZero,    kfInfinity },
		{ "coth(0)",    "0",   lz_hypercotangent,    kDivideByZero,    kfInfinity },
		{ "csch(0)",    "0",   lz_hypercosecant,     kDivideByZero,    kfInfinity },

		// Inside |x| = 1, where the arcsecant and arccosecant have no answer...
		{ "arcsec(.5)",  ".5", lz_arcsecant,         kFailed,          kfQNaN },
		{ "arcsec(0)",   "0",  lz_arcsecant,         kFailed,          kfQNaN },
		{ "arcsec(-.5)", "-.5",lz_arcsecant,         kFailed,          kfQNaN },
		{ "arccsc(.5)",  ".5", lz_arccosecant,       kFailed,          kfQNaN },
		{ "arccsc(0)",   "0",  lz_arccosecant,       kFailed,          kfQNaN },

		// The hyperbolic arccotangent: outside is a refusal, the edge
		//	itself an infinity...
		{ "arccoth(.5)", ".5", lz_archypercotangent, kFailed,          kfQNaN },
		{ "arccoth(0)",  "0",  lz_archypercotangent, kFailed,          kfQNaN },
		{ "arccoth(1)",  "1",  lz_archypercotangent, kUndefinedResult, kfInfinity },
		{ "arccoth(-1)", "-1", lz_archypercotangent, kUndefinedResult, kfInfinity },

		// The hyperbolic arcsecant, defined on (0, 1]...
		{ "arcsech(0)",  "0",  lz_archypersecant,    kUndefinedResult, kfInfinity },
		{ "arcsech(-1)", "-1", lz_archypersecant,    kFailed,          kfQNaN },
		{ "arcsech(2)",  "2",  lz_archypersecant,    kFailed,          kfQNaN },

		// The hyperbolic arccosecant, defined everywhere but nought...
		{ "arccsch(0)",  "0",  lz_archypercosecant,  kUndefinedResult, kfInfinity },
	};

	for ( auto& theEdge : edges )
	{
		CLZNumber	argument, outcome;
		char		szLabel[ 72 ];

		Set( argument, theEdge.szArgument );

		long lReturned = theEdge.pFunction( &argument, &outcome );

		::snprintf( szLabel, 72, "reports: %s", theEdge.szLabel );
		CheckTrue( szLabel, lReturned == theEdge.lExpect, "" );

		::snprintf( szLabel, 72, "flags: %s", theEdge.szLabel );
		CheckTrue( szLabel, outcome.FlagIsOn( theEdge.lFlag ), "" );
	}

	/*
		An infinity or a not-a-number going in comes out as a
		not-a-number from all twelve. Six of them have a finite limit
		at an infinity -- sech and csch are nought there, arcsec is
		pi/2, arccsc and both inverse hyperbolic reciprocals are nought
		-- but this engine's unary propagation offers only carry the
		infinity or refuse, and a finite limit is neither. They refuse
		rather than invent a number for an argument that has none.
	*/

	long (*everyone[])( CLZNumber*, CLZNumber* ) =
	{
		lz_cotangent, lz_secant, lz_cosecant,
		lz_hypercotangent, lz_hypersecant, lz_hypercosecant,
		lz_arccotangent, lz_arcsecant, lz_arccosecant,
		lz_archypercotangent, lz_archypersecant, lz_archypercosecant,
	};

	const char* theNames[] =
	{
		"cot", "sec", "csc", "coth", "sech", "csch",
		"arccot", "arcsec", "arccsc", "arccoth", "arcsech", "arccsch",
	};

	for ( size_t i = 0; i < sizeof(everyone) / sizeof(everyone[0]); i++ )
	{
		char szLabel[ 72 ];

		{
			CLZNumber argument, outcome;
			MakeInfinite( argument );
			everyone[ i ]( &argument, &outcome );
			::snprintf( szLabel, 72, "%s(infinity) is a NaN", theNames[ i ] );
			CheckTrue( szLabel, outcome.FlagIsOn( kfQNaN ), "" );
		}

		{
			CLZNumber argument, outcome;
			MakeNotANumber( argument );
			everyone[ i ]( &argument, &outcome );
			::snprintf( szLabel, 72, "%s(NaN) is a NaN", theNames[ i ] );
			CheckTrue( szLabel, outcome.FlagIsOn( kfQNaN ), "" );
		}
	}
}


static void TestIdentities()
{
	CLZNumber	a, b, x, y, result;

	Banner( "Identities" );

	// sqrt(x) squared returns x - exercises 64 Newton iterations...
	Set( a, "2" ); x.Clear();
	lz_root_2( &a, &x );
	result.Clear();
	lz_multiply( &x, &x, &result );
	CheckDigits( "sqrt(2) squared", &result, "2", 36 );

	Set( a, "1000000" ); x.Clear();
	lz_root_2( &a, &x );
	CheckDigits( "sqrt(1000000)", &x, "1000", 28 );

	// exp(ln(x)) returns x - exercises the 180-term log series...
	Set( a, "5" ); x.Clear();
	lz_log_e( &a, &x );
	result.Clear();
	lz_exp_e( &x, &result );
	CheckDigits( "exp(ln(5))", &result, "5", 36 );

	// sin^2 + cos^2 = 1...
	Set( a, "1" );
	x.Clear(); y.Clear();
	lz_sine( &a, &x );
	lz_cosine( &a, &y );
	lz_multiply( &x, &x );
	lz_multiply( &y, &y );
	result.Clear();
	lz_add( &x, &y, &result );
	CheckDigits( "sin(1)^2 + cos(1)^2", &result, "1", 36 );

	// tan = sin / cos...
	Set( a, "1" );
	x.Clear(); y.Clear();
	lz_sine( &a, &x );
	lz_cosine( &a, &y );
	result.Clear();
	lz_divide( &x, &y, &result );
	CheckDigits( "sin(1) / cos(1) = tan(1)", &result, kRefTan1, 38 );

	// cosh^2 - sinh^2 = 1...
	Set( a, "1" );
	x.Clear(); y.Clear();
	lz_hypercosine( &a, &x );
	lz_hypersine( &a, &y );
	lz_multiply( &x, &x );
	lz_multiply( &y, &y );
	result.Clear();
	lz_subtract( &x, &y, &result );
	CheckDigits( "cosh(1)^2 - sinh(1)^2", &result, "1", 36 );

	// A reciprocal, taken twice, returns the original -- to 36 digits:
	//	1/7 is a 36-place answer since 2026-09-25, where it carried a
	//	37th, and 7.000...007 is exactly what two correct roundings at 36
	//	places give...
	Set( a, "7" ); x.Clear();
	lz_reciprocal( &a, &x );
	result.Clear();
	lz_reciprocal( &x, &result );
	CheckDigits( "reciprocal twice of 7", &result, "7", 36 );
}

/////////////////////////////////////////////////////////////////
//	Integer and fractional parts - reached from the exponential
//	paths through StripInteger, but never directly until now
/////////////////////////////////////////////////////////////////

static void TestNumberParts()
{
	CLZNumber	a;

	Banner( "Integer and fractional parts" );

	// The predicates ask whether a part is absent...
	Set( a, "5" );      CheckTrue( "5 is not a fraction",     !a.IsFraction(), "" );
	Set( a, "5" );      CheckTrue( "5 is an integer",          a.IsInteger(),  "" );
	Set( a, "0.5" );    CheckTrue( ".5 is a fraction",         a.IsFraction(), "" );
	Set( a, "0.5" );    CheckTrue( ".5 is not an integer",    !a.IsInteger(),  "" );
	Set( a, "12.34" );  CheckTrue( "12.34 is not a fraction", !a.IsFraction(), "" );
	Set( a, "12.34" );  CheckTrue( "12.34 is not an integer", !a.IsInteger(),  "" );

	// A trailing zero in the integer part must not read as absent...
	Set( a, "10" );     CheckTrue( "10 is not a fraction",    !a.IsFraction(), "" );
	Set( a, "10" );     CheckTrue( "10 is an integer",         a.IsInteger(),  "" );

	// Splitting a value keeps each part whole...
	Set( a, "12.34" );  a.StripInteger();   CheckExact( "12.34 less its integer part",    &a, ".34" );
	Set( a, "12.34" );  a.StripFraction();  CheckExact( "12.34 less its fractional part", &a, "12" );

	// Stripping the part that is not there leaves the value alone...
	Set( a, "7" );      a.StripFraction();  CheckExact( "7 less its fractional part",     &a, "7" );
	Set( a, ".25" );    a.StripInteger();   CheckExact( ".25 less its integer part",      &a, ".25" );

	// Stripping the only part there is leaves zero...
	Set( a, "7" );      a.StripInteger();   CheckExact( "7 less its integer part",        &a, "0" );
	Set( a, ".25" );    a.StripFraction();  CheckExact( ".25 less its fractional part",   &a, "0" );

	// A carry-length value, where the digit count runs to the array...
	Set( a, "999.999" ); a.StripInteger();  CheckExact( "999.999 less its integer part",  &a, ".999" );
}

/////////////////////////////////////////////////////////////////
//	Number bases. Storage, parsing, rendering and carry propagation
//	are all base-aware already; this records that, so the base
//	conversion work has something underneath it. Note that base 12
//	spells ten and eleven D and E, not A and B.
/////////////////////////////////////////////////////////////////

static void SetInBase( CLZNumber& theNumber, long lBase, const char* szValue )
{
	theNumber.SetBase( lBase, true );
	theNumber.Assign( szValue );
}

static void TestNumberBases()
{
	Banner( "Number bases" );

	struct { long lBase; const char* szValue; } trips[] =
	{
		{  2, "1011" }, {  8, "777" }, { 10, "90210" },
		{ 12, "DE"   }, { 16, "FF"  }, { 16, "DEADBEEF" },
	};

	for ( auto& theTrip : trips )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "base %ld holds %s", theTrip.lBase, theTrip.szValue );

		CLZNumber a;
		SetInBase( a, theTrip.lBase, theTrip.szValue );
		CheckExact( szLabel, &a, theTrip.szValue );
	}

	// Carries and borrows have to run in the number's own base...
	struct { long lBase; const char* szLeft; char cOp; const char* szRight; const char* szExpect; } sums[] =
	{
		{  2, "1011", '+', "1",    "1100" },
		{  2, "111",  '+', "1",    "1000" },
		{  2, "1100", '-', "1",    "1011" },
		{  8, "7",    '+', "1",    "10"   },
		{  8, "777",  '+', "1",    "1000" },
		{  8, "1000", '-', "1",    "777"  },
		{ 10, "99",   '+', "1",    "100"  },
		{ 12, "E",    '+', "1",    "10"   },
		{ 12, "DE",   '+', "1",    "E0"   },
		{ 12, "100",  '-', "1",    "EE"   },
		{ 16, "FF",   '+', "1",    "100"  },
		{ 16, "100",  '-', "1",    "FF"   },
		{ 16, "F",    '*', "F",    "E1"   },
	};

	for ( auto& theSum : sums )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "base %ld: %s %c %s",
					theSum.lBase, theSum.szLeft, theSum.cOp, theSum.szRight );

		CLZNumber a, b, result;
		SetInBase( a, theSum.lBase, theSum.szLeft );
		SetInBase( b, theSum.lBase, theSum.szRight );
		result.SetBase( theSum.lBase, true );

		switch ( theSum.cOp )
		{
		case '+':	lz_add( &a, &b, &result );		break;
		case '-':	lz_subtract( &a, &b, &result );	break;
		default:	lz_multiply( &a, &b, &result );	break;
		}

		CheckExact( szLabel, &result, theSum.szExpect );
	}

	// Conversion between bases, integer values. Every pair is checked
	//	both ways, so neither direction can be wrong on its own...
	struct { const char* szValue; long lFrom; long lTo; const char* szExpect; } conversions[] =
	{
		{ "255",        10, 16, "FF"        },
		{ "11",         10,  2, "1011"      },
		{ "511",        10,  8, "777"       },
		{ "131",        10, 12, "DE"        },
		{ "FF",         16,  2, "11111111"  },
		{ "777",         8, 16, "1FF"       },
		{ "0",          10, 16, "0"         },
		{ "-255",       10, 16, "-FF"       },
		{ "3735928559", 10, 16, "DEADBEEF"  },

		// Long enough that the digit array has to grow on the way...
		{ "123456789012345678901234567890", 10, 16, "18EE90FF6C373E0EE4E3F0AD2" },
		{ "123456789012345678901234567890", 10, 12, "D94E8867D2385E6E5056131D416" },
		{ "618970019642690137449562111",    10, 16, "1FFFFFFFFFFFFFFFFFFFFFF" },
		{ "618970019642690137449562111",    10,  2,
		  "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111" },

		// Every other base SetBase accepts. Until 2026-09-24 these kept
		//	their digits and changed only the label: 21 in base seven
		//	still read 21. Expected digits generated, not typed...
		{ "21", 10,  3, "210" },
		{ "21", 10,  4, "111" },
		{ "21", 10,  5, "41" },
		{ "21", 10,  6, "33" },
		{ "21", 10,  7, "30" },
		{ "21", 10,  9, "23" },
		{ "21", 10, 11, "1A" },
		{ "21", 10, 13, "18" },
		{ "21", 10, 14, "17" },
		{ "21", 10, 15, "16" },
		{ "123456789012345678901234567890", 10,  7, "21653251153414601406403630240331250" },
		{ "123456789012345678901234567890", 10, 13, "1465A068454C58B8757720C20B0" },
		{ "FF", 16,  7, "513" },
		{ "666",  7,  3, "110200" },
		{ "-100", 10,  9, "-121" },
	};

	for ( auto& theCase : conversions )
	{
		char szLabel[ 64 ];

		::snprintf( szLabel, 64, "base %ld to %ld: %.18s",
					theCase.lFrom, theCase.lTo, theCase.szValue );

		CLZNumber a;
		SetInBase( a, theCase.lFrom, theCase.szValue );
		a.SetBase( theCase.lTo, true );
		CheckExact( szLabel, &a, theCase.szExpect );

		// ...and back again, which has to return the original...
		::snprintf( szLabel, 64, "base %ld to %ld and back: %.10s",
					theCase.lFrom, theCase.lTo, theCase.szValue );

		a.SetBase( theCase.lFrom, true );
		CheckExact( szLabel, &a, theCase.szValue );
	}

	// Fractions that land exactly in the new base. A fraction
	//	terminates only where the new base's prime factors cover the
	//	old denominator, so these are the cases with a right answer to
	//	the last digit...
	struct { const char* szValue; long lFrom; long lTo; const char* szExpect; } fractions[] =
	{
		{ ".5",     10,  2, ".1"     },
		{ ".25",    10,  2, ".01"    },
		{ ".75",    10,  2, ".11"    },
		{ ".5",     10,  8, ".4"     },
		{ ".5",     10, 16, ".8"     },
		{ ".0625",  10, 16, ".1"     },
		{ ".8",     16, 10, ".5"     },
		{ ".4",      8, 10, ".5"     },
		{ ".1",      2, 10, ".5"     },
		{ "10.5",   10,  2, "1010.1" },
		{ "-.5",    10,  2, "-.1"    },
		{ "3.75",   10, 16, "3.C"    },
		{ ".5",     10,  4, ".2"     },
		{ ".5",     10,  6, ".3"     },
		{ ".5",     10, 14, ".7"     },
		{ ".2",     10,  5, ".1"     },
		{ ".1",      5, 10, ".2"     },
		{ ".1",      3,  9, ".3"     },
		{ ".3",      9,  3, ".1"     },
		{ ".3",      6, 10, ".5"     },
		{ "2.4",     5, 10, "2.8"    },
		{ ".7",     14,  2, ".1"     },
		{ "-.25",   10,  4, "-.1"    },
	};

	for ( auto& theCase : fractions )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "base %ld to %ld: %s",
					theCase.lFrom, theCase.lTo, theCase.szValue );

		CLZNumber a;
		SetInBase( a, theCase.lFrom, theCase.szValue );
		a.SetBase( theCase.lTo, true );
		CheckExact( szLabel, &a, theCase.szExpect );
	}

	// And ones that do not land: a tenth is the standard example
	//	against base two. The budget is the working precision scaled by
	//	the ratio of the logs, so the value has to survive the trip out
	//	and back to the precision it started with, not to the digit...
	struct { const char* szValue; long lVia; } inexact[] =
	{
		{ ".1",    2 }, { ".333",  16 }, { ".7",   12 },
		{ "2.718", 16 }, { ".05",   8 },
		{ ".1",    3 }, { ".333",   7 }, { ".7",   13 },
		{ "2.718",  9 }, { ".05",  11 }, { ".9",   15 },
	};

	for ( auto& theCase : inexact )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "%s through base %ld and back",
					theCase.szValue, theCase.lVia );

		CLZNumber a;
		SetInBase( a, 10, theCase.szValue );
		a.SetBase( theCase.lVia, true );
		a.SetBase( 10, true );
		CheckDigits( szLabel, &a, theCase.szValue, 34 );
	}

	/*
		How many places of one base are worth so many of another. This
		is the conversion's own budget, and it is what a count measured
		in decimal -- gamma's widening cap, its shift threshold -- goes
		through before it sizes anything in another base. Rounded up.
	*/

	struct { long lPlaces; long lFrom; long lTo; long lExpect; } equivalents[] =
	{
		{ 36, 10,  2, 120 }, { 36, 10, 16,  30 }, { 44, 16, 10,  53 },
		{ 30, 10,  2, 100 }, { 30, 10,  8,  34 }, { 30, 10, 16,  25 },
		{ 36,  7,  7,  36 }, { 44, 10, 10,  44 },
	};

	for ( auto& theCase : equivalents )
	{
		char	szLabel[ 72 ];
		long	lGot = CLZNumber::EquivalentPlaces( theCase.lPlaces, theCase.lFrom, theCase.lTo );

		::snprintf( szLabel, 72, "%ld places of base %ld in base %ld",
					theCase.lPlaces, theCase.lFrom, theCase.lTo );
		CheckTrue( szLabel, lGot == theCase.lExpect, (lGot == theCase.lExpect) ? "" : "wrong count" );
	}

	/*
		An answer has to come back in the base it was asked in. The
		working values inside multiply and divide are made fresh, and
		a fresh number is base ten, so the digits were right and the
		label was wrong -- which is worse than either, because the
		next operation then either refuses the mismatch or reads hex
		digits as decimal.

		Checked against base ten rather than by eye: the same sum done
		in another base and converted back has to give the same
		answer.
	*/

	{
		struct { const char* szLabel; char cOp; const char* szExpect; } hexSums[] =
		{
			{ "FF + 10",  '+', "10F" },  { "FF - 10",  '-', "EF"  },
			{ "FF * 10",  '*', "FF0" },  { "FF / 10",  '/', "F.F" },
			{ "FF div 10",'d', "F"   },  { "FF mod 10",'m', "F"   },
			{ "FF rem 10",'r', "F"   },  { "FF remnear 10", 'n', "-1" },
		};

		for ( auto& theCase : hexSums )
		{
			CLZNumber a, b, result;
			SetInBase( a, 16, "FF" );
			SetInBase( b, 16, "10" );
			result.SetBase( 16, true );

			switch ( theCase.cOp )
			{
			case '+':	lz_add( &a, &b, &result );				break;
			case '-':	lz_subtract( &a, &b, &result );			break;
			case '*':	lz_multiply( &a, &b, &result );			break;
			case '/':	lz_divide( &a, &b, &result );			break;
			case 'd':	lz_divide_integer( &a, &b, &result );	break;
			case 'm':	lz_modulo( &a, &b, &result );			break;
			case 'r':	lz_remainder( &a, &b, &result );		break;
			default:	lz_remainder_near( &a, &b, &result );	break;
			}

			char szLabel[ 64 ];

			::snprintf( szLabel, 64, "base 16: %s", theCase.szLabel );
			CheckExact( szLabel, &result, theCase.szExpect );

			::snprintf( szLabel, 64, "answer stays base 16: %s", theCase.szLabel );
			CheckTrue( szLabel, result.m_lBase == 16, "" );
		}
	}

	/*
		Two faults recorded rather than fixed. Both come from the same
		place: a digit array carries no base of its own, so an operand
		in another base is read as though it were in this one. The
		engine does that to itself -- lz_raise_n decrements its
		exponent with a base-10 one, and lz_exp_b multiplies a base-10
		logarithm by a fraction in the operand's base.

		Refusing mismatched operands in the three primitives was tried
		once and could not work then: lz_raise_n's loop never
		terminated, because the engine hands its own primitives a
		base-ten counter against an operand in another base. This note
		ended "fixing it needs base conversion first", and that is
		what lzi_match_base is -- the operators convert now rather
		than refuse, and the engine doing it to itself is exactly the
		case that made conversion the right answer instead of a
		guard.
	*/

	/*
		Upper and lower case name the same digit, in every base that
		has letters for digits.

		Both the decoder and the validator have always accepted either,
		and **nothing tested it**: a mutation removing the lowercase
		arm of lzi_digit_value was caught by no test at all, found when
		the decoder's own sixteen-case copy of the mapping was
		collapsed into that function. One place deciding it is the
		reason to check it in one place too.

		The rendering is upper case whichever way it went in, which is
		what makes the second column the same for both rows.
	*/

	Banner( "Upper and lower case digits" );

	struct { long lBase; const char* szUpper; const char* szLower; const char* szExpect; } cases[] =
	{
		{ 16, "DEADBEEF", "deadbeef", "DEADBEEF" },
		{ 16, "FACE",     "face",     "FACE"     },
		{ 16, "A.B",      "a.b",      "A.B"      },
		{ 12, "DE",       "de",       "DE"       },
		{ 12, "D.E",      "d.e",      "D.E"      },
		{ 14, "CD",       "cd",       "CD"       },
		{ 11, "A",        "a",        "A"        },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber	upper, lower;
		char		szLabel[ 80 ];

		upper.SetBase( theCase.lBase, false );
		lower.SetBase( theCase.lBase, false );

		long	lUpper = upper.Assign( theCase.szUpper );
		long	lLower = lower.Assign( theCase.szLower );

		::snprintf( szLabel, 80, "base %ld %s reads", theCase.lBase, theCase.szLower );
		CheckTrue( szLabel, (lUpper == kNoError) && (lLower == kNoError) &&
						    lower.FlagIsOff( kfQNaN ), "" );

		::snprintf( szLabel, 80, "base %ld %s == %s", theCase.lBase,
					theCase.szLower, theCase.szUpper );
		CheckTrue( szLabel, lz_compare( &upper, &lower ) == kEqualTo, "" );

		::snprintf( szLabel, 80, "base %ld %s writes back upper", theCase.lBase, theCase.szLower );
		CheckRendering( szLabel, &lower, theCase.szExpect );
	}

	/*
		Operands in different bases are brought into one and answered,
		which they were not until base conversion went into the
		operators. The note above this used to end "fixing it needs
		base conversion first"; this is that.

		**The left operand's base wins**, which is what lzi_divide_op
		had always done with its quotient. So FF base sixteen plus ten
		base ten is 109 base sixteen, and the same sum the other way
		round is 265 base ten.
	*/

	struct { const char* szLabel; long lBase; const char* szLeft; char cOp;
			 const char* szRight; const char* szInBase; const char* szDecimal; } mixed[] =
	{
		{ "add",      16, "FF",  '+', "10", "109",  "265" },
		{ "subtract",  8, "70",  '-', "2",  "66",   "54"  },
		{ "multiply",  2, "101", '*', "3",  "1111", "15"  },
		{ "divide",   12, "E",   '/', "2",  "5.6",  "5.5" },

		// And bases with no entry in the conversion's log table until
		//	2026-09-24, where the right operand kept its decimal digits
		//	-- a 7 in base five, a 9 in base nine...
		{ "add base 5",      5, "2.4", '+', "7",  "14.4", "9.8" },
		{ "subtract base 7", 7, "30",  '-', "8",  "16",   "13"  },
		{ "multiply base 3", 3, "12",  '*', "4",  "202",  "20"  },
		{ "divide base 9",   9, "80",  '/', "18", "4",    "4"   },
	};

	for ( auto& theCase : mixed )
	{
		char szLabel[ 72 ];

		CLZNumber a, b, result;
		SetInBase( a, theCase.lBase, theCase.szLeft );
		b.Assign( theCase.szRight );		// left in base ten, deliberately

		long lResult = kNoError;

		switch ( theCase.cOp )
		{
		case '+':	lResult = lz_add( &a, &b, &result );		break;
		case '-':	lResult = lz_subtract( &a, &b, &result );	break;
		case '*':	lResult = lz_multiply( &a, &b, &result );	break;
		default:	lResult = lz_divide( &a, &b, &result );		break;
		}

		::snprintf( szLabel, 72, "mixed bases answered: %s", theCase.szLabel );
		CheckTrue( szLabel, lResult == kNoError, "" );

		::snprintf( szLabel, 72, "%s takes the left operand's base", theCase.szLabel );
		CheckTrue( szLabel, result.m_lBase == theCase.lBase, "" );
		CheckRendering( szLabel, &result, theCase.szInBase );

		// ...and is the same number whichever base it is written in...
		::snprintf( szLabel, 72, "%s is %s in base ten", theCase.szLabel,
					theCase.szDecimal );
		result.SetBase( 10, true );
		CheckRendering( szLabel, &result, theCase.szDecimal );
	}

	// The other way round, where the fault showed: a decimal left
	//	operand read the right one's base-seven digits as decimal, so
	//	1 + 3.3 base seven answered 4.3...
	struct { const char* szLabel; const char* szLeft; char cOp; long lBase;
			 const char* szRight; const char* szDecimal; } mixedRight[] =
	{
		{ "1 + 2.4 base 5",    "1",   '+', 5, "2.4", "3.8" },
		{ "100 - 30 base 7",   "100", '-', 7, "30",  "79"  },
		{ "5 * 12 base 3",     "5",   '*', 3, "12",  "25"  },
		{ "72 / 80 base 9",    "72",  '/', 9, "80",  "1"   },
	};

	for ( auto& theCase : mixedRight )
	{
		CLZNumber a, b, result;
		a.Assign( theCase.szLeft );
		SetInBase( b, theCase.lBase, theCase.szRight );

		switch ( theCase.cOp )
		{
		case '+':	lz_add( &a, &b, &result );		break;
		case '-':	lz_subtract( &a, &b, &result );	break;
		case '*':	lz_multiply( &a, &b, &result );	break;
		default:	lz_divide( &a, &b, &result );	break;
		}

		CheckRendering( theCase.szLabel, &result, theCase.szDecimal );
	}

	/*
		And the functions, which reach base-ten constants: an operand
		in base seven used to be combined with them digit for digit,
		so sqrt of three sevenths answered 1.015... Measured in base
		ten after converting back. Thirty-six places of base seven is
		about thirty decimal digits, hence the threshold.
	*/

	{
		struct { const char* szLabel; long (*fn)( CLZNumber*, CLZNumber* ); const char* szRef; } sevenths[] =
		{
			{ "sin of 3/7 in base 7",  lz_sine,   "0.41557185499305200807304366539942007870604329514826" },
			{ "ln of 3/7 in base 7",   lz_log_e,  "-0.84729786038720361371010750652065402498959417175911" },
			{ "sqrt of 3/7 in base 7", lz_root_2, "0.65465367070797714379829245624685835556920808239542" },
			{ "exp of 3/7 in base 7",  lz_exp_e,  "1.5350630092552098974443491684097758374476613892049" },
		};

		for ( auto& theCase : sevenths )
		{
			CLZNumber a, result;
			SetInBase( a, 7, ".3" );
			theCase.fn( &a, &result );
			result.SetBase( 10, true );
			CheckDigits( theCase.szLabel, &result, theCase.szRef, 29 );
		}
	}

	/*
		Assigning a C number to a value in another base. The text
		these build is decimal, so it has to be read as decimal and
		converted -- reading it in the value's own base put a 2 into a
		base-two number.

		This is what carried the fault furthest, because the engine
		builds its own constants that way: lz_log_e does Assign( 2L ),
		so the logarithm of a base-two value hung before the text was
		validated and answered nought after it. Nothing tested a
		transcendental outside base ten, so nothing noticed.
	*/

	struct { long lBase; long lValue; const char* szInBase; } assigned[] =
	{
		{  2,   2, "10"    },
		{  2,  10, "1010"  },
		{  8,   9, "11"    },
		{ 12,  11, "E"     },		// eleven is E in twelve
		{ 12,  12, "10"    },
		{ 16,  15, "F"     },
		{ 16, 255, "FF"    },
	};

	for ( auto& theCase : assigned )
	{
		CLZNumber	a;
		char		szLabel[ 72 ];

		a.SetBase( theCase.lBase, true );

		::snprintf( szLabel, 72, "base %ld holds Assign( %ldL )",
					theCase.lBase, theCase.lValue );

		CheckTrue( szLabel, a.Assign( theCase.lValue ) == kNoError, "" );
		CheckRendering( szLabel, &a, theCase.szInBase );
	}

	// The double form takes the same route...
	{
		CLZNumber a;
		a.SetBase( 2, true );
		CheckTrue( "base 2 holds Assign( 2.0 )", a.Assign( 2.0 ) == kNoError, "" );
		// the conversion trims the fractional noughts "%f" wrote...
		CheckRendering( "base 2 holds Assign( 2.0 )", &a, "10" );
	}

	/*
		And the case that was actually lost: a transcendental on a
		value outside base ten. ln(1010 base two) is ln(10), which is
		2.302585..., written in binary. Thirty-six places of base two
		is about eleven decimal digits, so that is what it is checked
		to -- the precision counts places of the value's own base.
	*/

	{
		CLZNumber a, result;
		a.SetBase( 2, true );
		a.Assign( "1010" );

		CheckTrue( "ln of a base-two ten", lz_log_e( &a, &result ) == kNoError, "" );

		result.SetBase( 10, true );
		CheckDigits( "ln of a base-two ten", &result, kRefLn10, 10 );
	}

	/*
		The comparisons need it too, and are easier to get wrong: they
		read the digit arrays directly, so a mismatch gave a confident
		answer rather than a refusal. 0x15 and 21 are the same number
		and were reported as different.
	*/

	{
		CLZNumber hex, ten;
		hex.Assign( "0x15" );
		ten.Assign( "0d21" );

		CheckTrue( "0x15 equals 21", lz_compare( &hex, &ten ) == kEqualTo, "" );

		CLZNumber bigger;
		bigger.Assign( "0x16" );
		CheckTrue( "0x16 exceeds 21", lz_compare( &bigger, &ten ) == kGreaterThan, "" );
		CheckTrue( "21 is under 0x16", lz_compare( &ten, &bigger ) == kLessThan, "" );
	}

	/*
		And the member comparisons, which the engine uses on itself and
		which read raw digits until 2026-09-24 -- lz_compare had been
		fixed and they had not. Gamma's shift compared a base-two
		argument with a 24 made in base ten, so 4.1, written 100.0001...,
		passed for more than 24 and the series was summed at 4. Both
		orders, both members, a fraction, a sign, and an odd base; and
		the argument must come back in the base it went in.
	*/

	{
		CLZNumber	binary, target;
		binary.Assign( "4.1" );
		binary.SetBase( 2, true );
		target.Assign( 24L );

		CheckTrue( "4.1 in base 2 is under 24", binary.Compare( &target ) == kLessThan, "" );
		CheckTrue( "24 exceeds 4.1 in base 2", target.Compare( &binary ) == kGreaterThan, "" );
		CheckTrue( "...and 24 stays in base 10", target.m_lBase == 10, "" );

		CLZNumber	hex, ten;
		hex.Assign( "0x15" );
		ten.Assign( "21" );

		CheckTrue( "0x15 Compare 21 is equal", hex.Compare( &ten ) == kEqualTo, "" );
		CheckTrue( "21 CompareABS 0x15 is equal", ten.CompareABS( &hex ) == kEqualTo, "" );

		// 21.5 is 30.333... in base seven, so its digits beat 0x1F's
		//	while its value is under thirty-one...
		CLZNumber	seven, negative, hex31;
		seven.Assign( "21.5" );
		seven.SetBase( 7, true );
		negative.Assign( "-21.5" );
		negative.SetBase( 7, true );
		hex31.Assign( "0x1F" );

		CheckTrue( "21.5 base 7 is under 0x1F", seven.Compare( &hex31 ) == kLessThan, "" );
		CheckTrue( "|-21.5| base 7 is under 0x1F", negative.CompareABS( &hex31 ) == kLessThan, "" );
		CheckTrue( "...and 0x1F stays in base 16", hex31.m_lBase == 16, "" );
	}

	/*
		And the remainder family, which reaches lzi_divide_op directly
		rather than through lz_divide.
	*/

	{
		CLZNumber a, b, result;
		a.Assign( "0x15" );					// twenty-one
		b.Assign( "0b1010" );				// ten

		CheckTrue( "modulo across bases",
				   lz_modulo( &a, &b, &result ) == kNoError, "" );
		result.SetBase( 10, true );
		CheckRendering( "21 modulo 10", &result, "1" );
	}

	// The engine used to do this to itself: lz_raise_n counted its
	//	exponent down with a base-ten one, so guarding the primitives
	//	stopped the loop terminating...
	{
		CLZNumber base, expo, result;
		base.Assign( "2" );
		SetInBase( expo, 16, "A" );			// ten, written in hex
		long lResult = lz_raise_n( &base, &expo, &result );

		CheckTrue( "raise to a hex exponent", lResult == kNoError, "" );
		CheckExact( "2 to the A base 16", &result, "1024" );
	}

	/*
		Exponentials in the operand's own base. b squared is written
		100 in base b, whatever b is, which makes the integer case
		checkable without converting anything.

		The fractional cases are measured in base ten instead, because
		AgreementDigits borrows in tens and cannot read a hexadecimal
		string. Converting the answer back before measuring exercises
		the round trip as well.
	*/

	struct { long lBase; const char* szExponent; } squares[] =
	{
		{ 2, "10" }, { 8, "2" }, { 10, "2" }, { 12, "2" }, { 16, "2" },

		// Bases with no logarithm of their own until 2026-09-24, which
		//	fell through to ln 10 and answered ten squared...
		{ 3, "2" }, { 7, "2" }, { 13, "2" },
	};

	for ( auto& theCase : squares )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "base %ld squared", theCase.lBase );

		CLZNumber a, result;
		SetInBase( a, theCase.lBase, theCase.szExponent );
		lz_exp_b( &a, &result );
		CheckExact( szLabel, &result, "100" );
	}

	{
		CLZNumber a, result;
		SetInBase( a, 2, "11" );			// three, in binary
		lz_exp_b( &a, &result );
		CheckExact( "base 2 cubed", &result, "1000" );
	}

	struct { long lBase; const char* szExponent; const char* szRef; } roots[] =
	{
		{  8, ".4", "2.8284271247461900976033774484193961571393437507539" },
		{ 16, ".8", "4" },
	};

	/*
		Base two is checked by rendering, against the correctly rounded
		binary, and the reason is worth a line. This row was the square
		root of two scored to 33 decimal digits, and lz_exp_b used to
		pass it with 37 -- because it worked in base ten and converted
		the answer back without rounding, so a base-two answer carried
		120 binary places. Since it widens by its answer's size and
		rounds, as lz_gamma and lz_root_2 always have, it answers with
		the 36 binary places a base-two value holds: about eleven
		digits, every one of them right, and the last bit correctly
		rounded. Expected strings from Python.
	*/

	struct { const char* szExponent; const char* szExpect; } binary[] =
	{
		{ ".1",   "1.011010100000100111100110011001111111"   },		// 2^(1/2)
		{ ".11",  "1.101011101000100111111001100101011011"   },		// 2^(3/4)
		{ "10.1", "101.101010000010011110011001100111111101" },		// 2^(5/2)
	};

	for ( auto& theCase : binary )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "2 to the %s, base 2, rounded", theCase.szExponent );

		CLZNumber a, result;
		SetInBase( a, 2, theCase.szExponent );
		lz_exp_b( &a, &result );
		CheckRendering( szLabel, &result, theCase.szExpect );
	}

	for ( auto& theCase : roots )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "base %ld to the half", theCase.lBase );

		CLZNumber a, result;
		SetInBase( a, theCase.lBase, theCase.szExponent );
		lz_exp_b( &a, &result );
		result.SetBase( 10, true );
		CheckDigits( szLabel, &result, theCase.szRef, 33 );
	}

	// Fractional powers in the bases that fell through to ln 10.
	//	References generated, and each raised back to its whole power
	//	to check them; the threshold is what thirty-six places of the
	//	operand's own base hold, measured in base ten...
	struct { long lBase; const char* szExponent; const char* szRef; int iDigits; } fractionals[] =
	{
		{  4, ".2", "2", 33 },								// .2 in base 4 is one half
		{  9, ".3", "2.0800838230519041145300568243578853863378053403733", 35 },	// a third
		{  7, ".3", "2.3024217239402238137470075067782911548599735484354", 31 },	// three sevenths
		{ 13, ".1", "1.2181140435607681966590257333736671593280798531997", 37 },	// a thirteenth
	};

	for ( auto& theCase : fractionals )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "base %ld to the %s", theCase.lBase, theCase.szExponent );

		CLZNumber a, result;
		SetInBase( a, theCase.lBase, theCase.szExponent );
		lz_exp_b( &a, &result );
		result.SetBase( 10, true );
		CheckDigits( szLabel, &result, theCase.szRef, theCase.iDigits );
	}

	/*
		lz_log_b, the logarithm to the operand's own base, which had no
		test anywhere until 2026-09-24 and took the same fall-through
		to ln 10 as lz_exp_b. The log of 100 is 2 in every base; the
		fractional cases are generated.
	*/

	const long kLogBases[] = { 2, 3, 7, 10, 13, 16 };

	for ( long lBase : kLogBases )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "log base %ld of 100", lBase );

		CLZNumber a, result;
		SetInBase( a, lBase, "100" );
		CheckTrue( szLabel, lz_log_b( &a, &result ) == kNoError, "" );
		result.SetBase( 10, true );
		CheckDigits( szLabel, &result, "2", 37 );
	}

	struct { long lBase; const char* szValue; const char* szRef; int iDigits; } logs[] =
	{
		{ 10, "2",  "0.30102999566398119521373889472449302676818988146211", 36 },
		{  2, "11", "1.5849625007211561814537389439478165087598144076925", 12 },	// log2 of 3
		{  7, ".3", "-0.43542496594642038619544983282509146385677208881321", 30 },	// log7 of 3/7
		{ 13, "A",  "0.89771171750262309292247769627456071248501032733323", 40 },	// log13 of 10
	};

	for ( auto& theCase : logs )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "log base %ld of %s", theCase.lBase, theCase.szValue );

		CLZNumber a, result;
		SetInBase( a, theCase.lBase, theCase.szValue );
		lz_log_b( &a, &result );
		result.SetBase( 10, true );
		CheckDigits( szLabel, &result, theCase.szRef, theCase.iDigits );
	}

	/*
		Rounding in an odd base, which has no digit for half a place.
		Its middle digit -- three in base seven -- is under half, and
		until 2026-09-24 every half mode took it for exactly half: 1.33
		base seven, 1.49 in decimal, rounded to 2. Long division and
		conversion made the same mistake at their own guard digits, and
		neither could say that a half in base seven is .3333... forever,
		a true tie that a finite string of digits cannot write.

		Every expected answer here was computed with exact rationals,
		and every odd-base row but the ties that round up anyway was
		answered wrongly before. The even-base rows are controls: the
		rounding there must not have moved. A random differential over
		30,000 cases in twelve bases found nothing.
	*/

	Banner( "Rounding in an odd base" );

	struct { long lBase; RoundingMode rmMode; long lPlaces; const char* szValue; const char* szExpect; } rounded[] =
	{
		{  7, RM_round_half_even, 0, "1.33",          "1"      },
		{  7, RM_round_half_even, 0, "1.3",           "1"      },
		{  7, RM_round_half_even, 0, "1.34",          "2"      },
		{  7, RM_round_half_up,   0, "2.3",           "2"      },
		{  7, RM_round_half_down, 0, "2.36",          "3"      },
		{  7, RM_round_half_up,   3, "6.41031025",    "6.410"  },
		{  5, RM_round_half_odd,  4, "1.0143222011",  "1.0143" },
		{ 13, RM_round_half_up,   0, "3.6",           "3"      },
		{  3, RM_round_half_even, 1, "2.11111",       "2.1"    },
		{  3, RM_round_half_even, 1, "2.1112",        "2.2"    },
		{  9, RM_round_half_even, 2, "0.4444",        "0.44"   },
		{  9, RM_round_half_even, 2, "0.4445",        "0.45"   },
		{ 15, RM_round_up,        0, "4.7",           "5"      },
		{ 15, RM_round_down,      0, "4.7",           "4"      },
		{ 10, RM_round_half_even, 0, "2.5",           "2"      },
		{ 10, RM_round_half_even, 0, "3.5",           "4"      },
		{ 10, RM_round_half_up,   0, "2.5",           "3"      },
		{ 10, RM_round_half_down, 0, "2.5",           "2"      },
		{ 16, RM_round_half_even, 0, "1.8",           "2"      },
		{  2, RM_round_half_even, 0, "10.1",          "10"     },
	};

	for ( auto& theCase : rounded )
	{
		char szLabel[ 72 ];
		::snprintf( szLabel, 72, "base %ld: %s to %ld places", theCase.lBase,
					theCase.szValue, theCase.lPlaces );

		CLZNumber a;
		a.SetBase( theCase.lBase, false );
		a.Assign( theCase.szValue );
		a.Round( theCase.rmMode, theCase.lPlaces );
		CheckRendering( szLabel, &a, theCase.szExpect );
	}

	/*
		Whole powers in an odd base. lz_raise_n squares its way along
		the exponent's binary digits and took the parity from the last
		digit, on the grounds that every base was even -- and seven is
		10 in base seven. So 2 to the 7 answered 64 there, and e to the
		7.25 answered e to the 6.25, the exponential reaching its whole
		part through this. The exponent is written so that its last
		digit is even and the number is odd.
	*/

	struct { long lBase; const char* szLeft; const char* szExponent; const char* szExpect; } powers[] =
	{
		{  7, "2", "10",  "242"    },		// 2^7
		{  5, "3", "10",  "1433"   },		// 3^5
		{  3, "2", "100", "200222" },		// 2^9
		{ 13, "2", "10",  "3962"   },		// 2^13
	};

	for ( auto& theCase : powers )
	{
		char szLabel[ 72 ];
		::snprintf( szLabel, 72, "base %ld: %s to the %s", theCase.lBase,
					theCase.szLeft, theCase.szExponent );

		CLZNumber a, n, result;
		SetInBase( a, theCase.lBase, theCase.szLeft );
		SetInBase( n, theCase.lBase, theCase.szExponent );
		lz_raise_n( &a, &n, &result );
		CheckRendering( szLabel, &result, theCase.szExpect );
	}

	/*
		The digits are what 36 places of base seven hold, about 34 after
		a four-digit whole part, and both answers are the correctly
		rounded ones, checked in base seven against Python's Fraction.
		lz_raise_e used to score 40 here: it raised a base-ten e and
		converted back, and the answer carried the decimal's width
		unrounded -- see CLZAnswerBase.
	*/

	{
		struct { const char* szLabel; long (*fn)( CLZNumber*, CLZNumber* ); const char* szValue; const char* szRef; int iDigits; const char* szBaseSeven; } whole[] =
		{
			{ "exp of 7.25 in base 7",   lz_exp_e,   "7.25", "1408.10484820469557502008632701363782242007184452645", 35, "4051.050651200150462223360120236211011440" },
			{ "raise_e of 7 in base 7",  lz_raise_e, "7",    "1096.63315842845859926372023828812143244221913483361", 35, "3124.430113312251013021151654322661541524" },
		};

		for ( auto& theCase : whole )
		{
			char		szLabel[ 72 ];
			CLZNumber	a, result;

			a.Assign( theCase.szValue );
			a.SetBase( 7, true );
			theCase.fn( &a, &result );

			::snprintf( szLabel, 72, "%s, its digits", theCase.szLabel );
			CheckRendering( szLabel, &result, theCase.szBaseSeven );

			result.SetBase( 10, true );
			CheckDigits( theCase.szLabel, &result, theCase.szRef, theCase.iDigits );
		}
	}

	// Quotients. A remainder of exactly half the divisor is the tie,
	//	and the kept digit decides it: a half in base five is .2222...
	//	and stays, in base seven it is .3333... and goes up...
	struct { long lBase; const char* szLeft; const char* szRight; const char* szExpect; } quotients[] =
	{
		{  5, "1",      "2", "0.22222222222222222222222222222222222222" },
		{  3, "1",      "2", "0.11111111111111111111111111111111111112" },
		{  7, "1",      "2", "0.33333333333333333333333333333333333334" },
		{  9, "3",      "6", "0.44444444444444444444444444444444444444" },
		{ 13, "B",      "2", "5.6666666666666666666666666666666666666" },
		{  5, "342",    "2", "143.2222222222222222222222222222222222222" },
		{  7, "335415", "2", "152542.3333333333333333333333333333333333334" },
		{ 10, "1",      "8", "0.125" },
		{ 16, "1",      "3", "0.55555555555555555555555555555555555555" },
	};

	for ( auto& theCase : quotients )
	{
		char szLabel[ 72 ];
		::snprintf( szLabel, 72, "base %ld: %s / %s", theCase.lBase,
					theCase.szLeft, theCase.szRight );

		CLZNumber a, b, result;
		a.SetBase( theCase.lBase, false );
		a.Assign( theCase.szLeft );
		b.SetBase( theCase.lBase, false );
		b.Assign( theCase.szRight );
		lz_divide( &a, &b, &result );
		CheckRendering( szLabel, &result, theCase.szExpect );
	}

	// And conversions into an odd base, where the leftover beyond the
	//	guard is still in the old base and is compared with a half
	//	written there...
	struct { long lBase; const char* szDecimal; const char* szExpect; } converted[] =
	{
		{  5, "1631.5",  "23011.2222222222222222222222222222222222222222222222222222" },
		{  7, "0.5",     "0.3333333333333333333333333333333333333333334" },
		{ 13, "42447.5", "16422.666666666666666666666666666666666" },
		{  7, "0.11231", "0.053344105413354461345553243034331564415436" },
		{ 10, "0.5",     "0.5" },
		{ 16, "0.1",     "0.19999999999999999999999999999A" },
	};

	for ( auto& theCase : converted )
	{
		char szLabel[ 72 ];
		::snprintf( szLabel, 72, "%s into base %ld", theCase.szDecimal, theCase.lBase );

		CLZNumber a;
		a.Assign( theCase.szDecimal );
		a.SetBase( theCase.lBase, true );
		CheckRendering( szLabel, &a, theCase.szExpect );
	}

	/*
		Widths in another base that the tables have to be told about.
		Their thresholds are decimal digits, and until 2026-09-24 they
		were compared with a count of hex places: 44 of them are 53
		digits, and the sine stayed on a table good for forty, the
		exponential took halvings sized for 44 digits where it needed
		53, and at 33 places (40 digits) took none where it needed one.
		The arcsine's coefficients, made in base ten, arrived in hex
		with the thirty places their own precision was worth.

		Each row fails with one of those faults put back: the sine with
		the series threshold (42 digits), the arcsine with the
		conversion budget (39, and the exponential at 33 with it, 38),
		the exponential at 44 with its threshold (52) or with its
		recurrence working in base ten (51), and at 33 with its halvings
		(39). The inputs are exact in hex; the references are the Taylor
		series in Python's Decimal at 120 digits, and its ln.

		The floors went up on 2026-09-25, when every function's work
		moved into its operand's base: the sine and arcsine from 44 and
		42 digits to 53, and the logarithm row was added for lz_log_e's
		e, which divided its reduction in base ten (53 digits without
		it, 54 with). The sine's series summing in base ten fails the
		first row, the arcsine's one in base ten the second.
	*/

	{
		struct { const char* szLabel; long lPlaces; char cFn; const char* szHex; const char* szRef; int iDigits; } wide[] =
		{
			{ "sin 0.B at 44 hex places",    44, 's', "0.B", "0.6346070800152692968503099142036714408128736065378100879387283929", 52 },
			{ "arcsin 0.4 at 44 hex places", 44, 'a', "0.4", "0.2526802551420786534856574369937109722521937330968381936339237787", 52 },
			{ "exp 0.E at 44 hex places",    44, 'e', "0.E", "2.398875293967097914691645733716296748285465871698578890902254684",  53 },
			{ "exp 0.E at 33 hex places",    33, 'e', "0.E", "2.398875293967097914691645733716296748285465871698578890902254684",  40 },
			{ "ln 3E8.8 at 44 hex places",   44, 'l', "3E8.8", "6.9082551540237880999680379796364568605088836522256740185343200532698746983715649", 54 },
		};

		for ( auto& theCase : wide )
		{
			long	l_saved = glWorkingPrecision;

			// Before the numbers are made: a CLZNumber takes its
			//	precision at construction...
			glWorkingPrecision = theCase.lPlaces;

			CLZNumber	x, answer;

			SetInBase( x, 16, theCase.szHex );

			if ( theCase.cFn == 's' )
				lz_sine( &x, &answer );
			else if ( theCase.cFn == 'a' )
				lz_arcsine( &x, &answer );
			else if ( theCase.cFn == 'l' )
				lz_log_e( &x, &answer );
			else
				lz_exp_e( &x, &answer );

			answer.SetBase( 10, true );

			glWorkingPrecision = l_saved;

			CheckDigits( theCase.szLabel, &answer, theCase.szRef, theCase.iDigits );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	Failures that used to be quiet. Every one of these once came
//	back as kNoError, several with a plausible-looking number
//	attached: ln(-5) answered 1659318788593735693051314843533547
//	208079, and the square root of -9 answered 5.311...
/////////////////////////////////////////////////////////////////

static void TestQuietFailures()
{
	Banner( "Failures that report" );

	struct { const char* szLabel; const char* szLeft; const char* szRight;
			 char cOp; long lExpect; unsigned long lFlag; } cases[] =
	{
		{ "ln(0) is infinite",      "0",  nullptr, 'L', kUndefinedResult, kfInfinity },
		{ "ln(-5) is not real",     "-5", nullptr, 'L', kFailed,          kfQNaN     },
		{ "sqrt(-9) is not real",   "-9", nullptr, 'R', kFailed,          kfQNaN     },
		{ "1/0",                    "0",  nullptr, 'V', kDivideByZero,    kfDiv0     },
		// A tan(pi/2) row stood here, pinning a refusal that was never a
		//	pole: pi/2 to 49 places until the answer width was fixed, then
		//	to 130 until pi was computed past its 120. No written argument
		//	is a pole; tan answers them -- see TestAnswerWidth and
		//	TestWidePi...
		// arcsin(2) returned kFailed and wrote nothing until 2026-09-26,
		//	the one refusal in the engine that did; the flag pins it...
		{ "arcsin(2) is out of domain", "2", nullptr, 'A', kFailed,      kfQNaN     },
		{ "arccos(2) is out of domain", "2", nullptr, 'c', kFailed,      kfQNaN     },
		{ "arccosh(.5) is out of domain", ".5", nullptr, 'h', kFailed,   kfQNaN     },
		// The edge of arctanh's domain, where the limit is an infinity,
		//	against the outside of it, where there is no answer at all.
		//	Both used to report a flat kFailed...
		{ "arctanh(1) is infinite",  "1",  nullptr, 'g', kUndefinedResult, kfInfinity },
		{ "arctanh(-1) is infinite", "-1", nullptr, 'g', kUndefinedResult, kfInfinity },
		{ "arctanh(2) is out of domain", "2", nullptr, 'g', kFailed,   kfQNaN     },
		{ "arctanh(-2) is out of domain","-2",nullptr, 'g', kFailed,   kfQNaN     },

		/*
			The factorial here is the product loop, defined on the
			non-negative integers. Everything else used to spin
			forever counting away from a nought it would never reach:
			-3 counted down past it, 2.5 stepped over it.
		*/
		{ "(-3)! has no answer",    "-3",   nullptr, '!', kFailed, kfQNaN },
		{ "(2.5)! wants gamma",     "2.5",  nullptr, '!', kFailed, kfQNaN },
		{ "(-0.5)! wants gamma",    "-0.5", nullptr, '!', kFailed, kfQNaN },

		/*
			lz_root_n had no guards at all and dropped every return
			code: a negative went into lz_log_e, which refuses it, and
			the refusal's nought came back out as the answer. The cube
			root of -8 was 0 reported as success.
		*/
		{ "even root of a negative", "-8", "2",   'N', kFailed, kfQNaN },
		{ "root of order nought",    "8",  "0",   'N', kFailed, kfQNaN },
		{ "fractional root of a negative", "-8", "2.5", 'N', kFailed, kfQNaN },
		{ "negative order of nought", "0", "-3",  'N', kUndefinedResult, kfInfinity },
		{ "5 / 0",                  "5",  "0",     '/', kDivideByZero,    kfDiv0     },
		{ "0 / 0",                  "0",  "0",     '/', kUndefinedResult, kfQNaN     },
		{ "5 modulo 0",             "5",  "0",     '%', kDivideByZero,    0          },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber	a, b, result;
		long		lReturned = kNoError;

		Set( a, theCase.szLeft );

		if ( theCase.szRight != nullptr )
			Set( b, theCase.szRight );

		switch ( theCase.cOp )
		{
		case 'L':	lReturned = lz_log_e( &a, &result );			break;
		case 'R':	lReturned = lz_root_2( &a, &result );			break;
		case 'V':	lReturned = lz_reciprocal( &a, &result );		break;
		case 'T':	lReturned = lz_tangent( &a, &result );			break;
		case 'A':	lReturned = lz_arcsine( &a, &result );			break;
		case 'c':	lReturned = lz_arccosine( &a, &result );		break;
		case 'h':	lReturned = lz_archypercosine( &a, &result );	break;
		case 'g':	lReturned = lz_archypertangent( &a, &result );	break;
		case '!':	lReturned = lz_factorial( &a, &result );		break;
		case 'N':	lReturned = lz_root_n( &a, &b, &result );		break;
		case '%':	lReturned = lz_modulo( &a, &b, &result );		break;
		default:	lReturned = lz_divide( &a, &b, &result );		break;
		}

		char szLabel[ 72 ];

		::snprintf( szLabel, 72, "reports: %s", theCase.szLabel );
		CheckTrue( szLabel, lReturned == theCase.lExpect, "" );

		if ( theCase.lFlag != 0 )
		{
			::snprintf( szLabel, 72, "flags: %s", theCase.szLabel );
			CheckTrue( szLabel, result.FlagIsOn( theCase.lFlag ), "" );
		}
	}

	// ...while the edge of the domain that does have an answer keeps it...
	{
		CLZNumber a, result;
		Set( a, "0" );
		CheckTrue( "sqrt(0) is not a failure", lz_root_2( &a, &result ) == kNoError, "" );
		CheckExact( "sqrt(0)", &result, "0" );

		Set( a, "1" );
		CheckTrue( "ln(1) is not a failure", lz_log_e( &a, &result ) == kNoError, "" );
	}

	/*
		An odd root of a negative is real and has to come back, which
		is the other half of the lz_root_n fix. lz_root_3 always got
		these right; lz_root_n answered nought for every one of them.
	*/

	Banner( "Odd roots of negatives" );

	struct { const char* szValue; const char* szOrder; const char* szExpect; } odds[] =
	{
		{ "-8",   "3", "-2" }, { "-27", "3", "-3" },
		{ "-32",  "5", "-2" }, { "-1",  "3", "-1" },
	};

	for ( auto& theOdd : odds )
	{
		CLZNumber	a, n, result;
		char		szLabel[ 64 ];

		Set( a, theOdd.szValue );
		Set( n, theOdd.szOrder );

		::snprintf( szLabel, 64, "root %s of %s", theOdd.szOrder, theOdd.szValue );

		CheckTrue( szLabel, lz_root_n( &a, &n, &result ) == kNoError, "" );
		CheckDigits( szLabel, &result, theOdd.szExpect, 35 );
	}

	// Nought still roots to nought for a positive order...
	{
		CLZNumber a, n, result;
		Set( a, "0" ); Set( n, "3" );
		CheckTrue( "root 3 of 0 is not a failure",
				   lz_root_n( &a, &n, &result ) == kNoError, "" );
		CheckExact( "root 3 of 0", &result, "0" );
	}

	/*
		Text that is not a number. Every one of these used to return
		kNoError with something plausible attached -- "12x34" became
		1234, and the minus inside "1.5E-8" made the whole value
		negative. See lzi_text_is_a_number.
	*/

	Banner( "Text that is not a number" );

	const char* rubbish[] =
	{
		"1.2.3", "", "-", ".", "1 2", "twelve", "1,000", "12y34",
		"--1", "1-2", "1.2-",

		// Exponents with nothing in them, two of them, or one too
		//	long to be anything but a typo...
		"1e", "e5", "1e+", "1e5.5", "1ee5", "1X", "X5", "1e99999999",

		// Base prefixes with nothing after them, a digit the base
		//	does not have, or an identifier that is not one...
		"0x", "0b", "0t", "0b2", "0o8", "0tC", "0z1", "0d1F",
	};

	for ( auto& szText : rubbish )
	{
		CLZNumber	a;
		char		szLabel[ 64 ];

		::snprintf( szLabel, 64, "refuses \"%s\"", szText );
		CheckTrue( szLabel, a.Assign( szText ) == kFailed, "" );

		::snprintf( szLabel, 64, "flags \"%s\"", szText );
		CheckTrue( szLabel, a.FlagIsOn( kfQNaN ), "" );
	}

	// ...and everything that is a number still reads...
	const char* proper[] =
	{
		"1", "0", "-5", "+5", ".5", "-.5", "5.", "3.14159", "000123",
		"-0.0000000001", "99999999999999999999",

		// Exponent notation, both marks and both signs...
		"6.02e23", "1.5E-8", "1e5", "1X5", "-2.5e-3", "+1.5e2",
		".5e1", "5.e2", "0e10", "1e-40", "123.456e-3",
	};

	for ( auto& szText : proper )
	{
		CLZNumber	a;
		char		szLabel[ 64 ];

		::snprintf( szLabel, 64, "reads \"%s\"", szText );
		CheckTrue( szLabel, a.Assign( szText ) == kNoError, "" );
	}

	/*
		Exponent notation, to Elzed's convention, since the two are to
		merge. X marks an exponent as well as E, because E is a digit
		in twelve, fifteen and sixteen; and the base raised is the
		number's own, so 1.0X5 in base sixteen is sixteen to the fifth.

		The value matters as much as the acceptance -- reading "1e5"
		without complaint and storing the wrong thing is the fault
		this replaced.
	*/

	Banner( "Exponent notation" );

	struct { const char* szText; const char* szPlain; } exponents[] =
	{
		{ "1e5",         "100000"    },
		{ "1E5",         "100000"    },
		{ "1X5",         "100000"    },
		{ "1e0",         "1"         },
		{ "6.02e23",     "602000000000000000000000" },
		{ "1.5e-8",      ".000000015" },
		{ "-2.5e-3",     "-.0025"    },
		{ "+1.5e2",      "150"       },
		{ ".5e1",        "5"         },
		{ "5.e2",        "500"       },
		{ "123.456e-3",  ".123456"   },
		{ "0e10",        "0"         },
		{ "1e-40",       ".0000000000000000000000000000000000000001" },
	};

	for ( auto& theCase : exponents )
	{
		CLZNumber	a;
		char		szLabel[ 64 ];

		::snprintf( szLabel, 64, "\"%s\" reads as %s", theCase.szText,
					theCase.szPlain );

		CheckTrue( szLabel, a.Assign( theCase.szText ) == kNoError, "" );
		CheckExact( szLabel, &a, theCase.szPlain );
	}

	/*
		Base prefixes, the other half of Elzed's entry convention. The
		leading character is always a nought and the identifier always
		lower case; a prefix takes the value's base with it, so after
		Assign("0x15") the number is a base-sixteen one holding
		fifteen sixteen.

		Twenty-one written every way Elzed accepts, which is the
		example its own documentation gives.
	*/

	Banner( "Base prefixes" );

	struct { const char* szText; long lBase; const char* szPlain; } prefixes[] =
	{
		{ "0b10101", 2,  "10101" },
		{ "0o25",    8,  "25"    },
		{ "21",      10, "21"    },
		{ "0d21",    10, "21"    },
		{ "0t19",    12, "19"    },
		{ "0x15",    16, "15"    },
		{ "0h15",    16, "15"    },
	};

	for ( auto& thePrefix : prefixes )
	{
		CLZNumber	a, decimal;
		char		szLabel[ 72 ];

		::snprintf( szLabel, 72, "%s is twenty-one", thePrefix.szText );

		CheckTrue( szLabel, a.Assign( thePrefix.szText ) == kNoError, "" );
		CheckTrue( szLabel, a.m_lBase == thePrefix.lBase, "wrong base" );
		CheckRendering( szLabel, &a, thePrefix.szPlain );

		// ...and every one of them is the same number underneath...
		decimal.Assign( &a );
		decimal.SetBase( 10, true );
		CheckRendering( szLabel, &decimal, "21" );
	}

	// Extended digits, signs, radix marks and exponents all compose
	//	with a prefix...
	struct { const char* szText; const char* szDecimal; } composed[] =
	{
		{ "0xFF",       "255"           },
		{ "0tDE",       "131"           },	// D is ten and E eleven in twelve
		{ "-0x15",      "-21"           },
		{ "+0b101",     "5"             },
		{ "0x1.8",      "1.5"           },
		{ "0x1Ax2",     "6656"          },	// 1A sixteen, times sixteen squared
		{ "0d1.5e2",    "150"           },
	};

	for ( auto& theCase : composed )
	{
		CLZNumber	a;
		char		szLabel[ 72 ];

		::snprintf( szLabel, 72, "%s is %s", theCase.szText, theCase.szDecimal );

		CheckTrue( szLabel, a.Assign( theCase.szText ) == kNoError, "" );

		a.SetBase( 10, true );
		CheckRendering( szLabel, &a, theCase.szDecimal );
	}

	/*
		The identifier is lower case, and that is what keeps a prefix
		apart from an exponent: a lower-case x at the front is
		hexadecimal, an upper-case one anywhere is an exponent mark.
		"0X15" is nought times ten to the fifteenth, which is nought.
	*/

	{
		CLZNumber a;
		CheckTrue( "0X15 is not a prefix", a.Assign( "0X15" ) == kNoError, "" );
		CheckTrue( "0X15 stays base ten", a.m_lBase == 10, "" );
		CheckRendering( "0X15 is nought", &a, "0" );
	}

	/*
		And across the bases, where the convention earns its keep. In
		sixteen and twelve E is a digit, so "1E5" is a plain number
		there and only X introduces an exponent.
	*/

	struct { long lBase; const char* szText; const char* szPlain; } basedExponents[] =
	{
		{ 16, "1E5",    "1E5"    },		// a plain hex number, not an exponent
		{ 16, "1.0X5",  "100000" },		// 1.0 sixteen, times sixteen to the fifth
		{ 12, "1E5",    "1E5"    },		// eleven is a digit in twelve
		{ 12, "1X3",    "1000"   },		// one, times twelve cubed
		{  2, "1.01e3", "1010"   },		// binary, times two cubed
		{  8, "1e2",    "100"    },		// octal, times sixty-four
		{ 10, "6.02X23","602000000000000000000000" },
	};

	for ( auto& theCase : basedExponents )
	{
		CLZNumber	a;
		char		szLabel[ 72 ];

		a.SetBase( theCase.lBase, true );

		::snprintf( szLabel, 72, "base %ld reads \"%s\" as %s",
					theCase.lBase, theCase.szText, theCase.szPlain );

		CheckTrue( szLabel, a.Assign( theCase.szText ) == kNoError, "" );
		CheckExact( szLabel, &a, theCase.szPlain );
	}

	/*
		Writing an exponent. The mark follows the same rule as reading
		one -- E where the base leaves it free, X where E is a digit --
		so what the engine writes, it reads.
	*/

	Banner( "Writing an exponent" );

	struct { const char* szText; const char* szSci; const char* szEng; } styles[] =
	{
		{ "12345",      "1.2345E4",   "12.345E3"  },
		{ "1",          "1E0",        "1E0"       },
		{ "0",          "0",          "0"         },
		{ "-1",         "-1E0",       "-1E0"      },
		{ ".000015",    "1.5E-5",     "15E-6"     },
		{ "-.0025",     "-2.5E-3",    "-2.5E-3"   },
		{ "999.9",      "9.999E2",    "999.9E0"   },
		{ "123456789",  "1.23456789E8", "123.456789E6" },
		{ "7",          "7E0",        "7E0"       },
	};

	for ( auto& theStyle : styles )
	{
		CLZNumber	a;
		char		szLabel[ 64 ];

		a.Assign( theStyle.szText );

		::snprintf( szLabel, 64, "%s scientific", theStyle.szText );
		CheckStyled( szLabel, &a, RS_scientific, theStyle.szSci );

		::snprintf( szLabel, 64, "%s engineering", theStyle.szText );
		CheckStyled( szLabel, &a, RS_engineering, theStyle.szEng );
	}

	/*
		The engineering form pads its mantissa to fill the group, so
		100E-3 is a tenth written with three digits where the value
		had one. The value is the same; the digits it carries are not.
	*/

	{
		CLZNumber a;
		a.Assign( ".1" );
		CheckStyled( "a tenth, scientific", &a, RS_scientific, "1E-1" );
		CheckStyled( "a tenth, engineering", &a, RS_engineering, "100E-3" );
	}

	// Where E is a digit, X does the work -- and the exponent counts
	//	the number's own base...
	struct { long lBase; const char* szText; const char* szSci; } basedStyles[] =
	{
		{ 16, "FF",   "F.FX1"    },
		{ 16, "1E5",  "1.E5X2"   },
		{ 12, "DE",   "D.EX1"    },
		{  2, "1010", "1.010E3"  },
		{  8, "777",  "7.77E2"   },
	};

	for ( auto& theStyle : basedStyles )
	{
		CLZNumber	a;
		char		szLabel[ 72 ];

		a.SetBase( theStyle.lBase, true );
		a.Assign( theStyle.szText );

		::snprintf( szLabel, 72, "base %ld %s scientific",
					theStyle.lBase, theStyle.szText );

		CheckStyled( szLabel, &a, RS_scientific, theStyle.szSci );
	}

	/*
		And the property that matters more than any particular
		spelling: what the engine writes, it reads back as the same
		value. Checked numerically, not by matching text, because the
		engineering form's padding changes the digits without changing
		the value.
	*/

	struct { long lBase; const char* szText; } trips2[] =
	{
		{ 10, "602000000000000000000000" }, { 10, "-98765.4321" },
		{ 10, "1e-40" }, { 10, ".5" }, { 10, "3.14159265358979323846" },
		{ 16, "DEADBEEF" }, { 16, "-A.B" }, { 12, "DE" },
		{  2, "1.01" }, {  8, "-7.7" },
	};

	for ( auto& theTrip : trips2 )
	{
		for ( int iStyle = 0; iStyle < 2; iStyle++ )
		{
			RenderStyle	rsStyle = iStyle ? RS_engineering : RS_scientific;

			CLZNumber	a, b;
			char		szLabel[ 72 ];

			a.SetBase( theTrip.lBase, true );
			a.Assign( theTrip.szText );

			char*	szWritten = a.GetStringRep( rsStyle );

			b.SetBase( theTrip.lBase, true );

			long	lReturned = b.Assign( szWritten );

			::snprintf( szLabel, 72, "base %ld %s round trips as %s",
						theTrip.lBase, theTrip.szText,
						iStyle ? "engineering" : "scientific" );

			CheckTrue( szLabel,
					   (lReturned == kNoError) && (lz_compare( &a, &b ) == kEqualTo), "" );

			::free( szWritten );
		}
	}

	/*
		Nought written with more than one digit used to render as the
		empty string -- every digit suppressed as a leading zero and
		nothing left. Only reachable through text, which is why it
		lasted; expanding 0e10 into eleven noughts is what found it.
	*/

	/*
		CheckRendering and not CheckExact: this is about the characters
		written, and CheckExact normalises both sides, so it reads the
		empty string and "0" as the same number -- which they are. It
		passed the whole time the bug was there, and a mutation that
		put the bug back is what showed it.
	*/

	{
		CLZNumber a;
		a.Assign( "00" );   CheckRendering( "\"00\" renders as a nought", &a, "0" );
		a.Assign( "000" );  CheckRendering( "\"000\" renders as a nought", &a, "0" );
		a.Assign( "-000" ); CheckRendering( "\"-000\" renders as a nought", &a, "0" );
		a.Assign( "0e10" ); CheckRendering( "\"0e10\" renders as a nought", &a, "0" );
	}

	/*
		A digit has to be one in the base the value is held in. Base
		twelve spells ten and eleven D and E, so C is not a digit
		there; base ten has no A through F at all, which is how
		"6.02e23" used to put a fourteen into a base-ten number.
	*/

	struct { long lBase; const char* szText; bool bValid; } digits[] =
	{
		{ 10, "A",  false }, { 10, "F",  false }, { 10, "99", true  },
		{ 12, "D",  true  }, { 12, "E",  true  }, { 12, "C",  false },
		{ 12, "A",  true  }, { 16, "F",  true  }, { 16, "G",  false },
		{  2, "1",  true  }, {  2, "2",  false }, {  8, "8",  false },
	};

	for ( auto& theDigit : digits )
	{
		CLZNumber	a;
		char		szLabel[ 72 ];

		a.SetBase( theDigit.lBase, true );

		long lReturned = a.Assign( theDigit.szText );

		::snprintf( szLabel, 72, "base %ld %s \"%s\"", theDigit.lBase,
					theDigit.bValid ? "reads" : "refuses", theDigit.szText );

		CheckTrue( szLabel,
				   theDigit.bValid ? (lReturned == kNoError) : (lReturned == kFailed), "" );
	}
}

/////////////////////////////////////////////////////////////////
//	Carrying a failure onward. An infinity is stored as zero with a
//	flag and so has no sign, which is why two of them added or
//	subtracted give a NaN: the cancelling case cannot be told from
//	the reinforcing one, and a guess would be worse than an answer
//	of "no idea".
/////////////////////////////////////////////////////////////////

static void MakeInfinite( CLZNumber& theNumber )
{
	CLZNumber	a, b;
	Set( a, "5" ); Set( b, "0" );
	lz_divide( &a, &b, &theNumber );
}

static void MakeNotANumber( CLZNumber& theNumber )
{
	CLZNumber	z;
	Set( z, "0" );
	lz_divide( &z, &z, &theNumber );
}

static void TestFlagPropagation()
{
	Banner( "Carrying a failure onward" );

	struct { const char* szLabel; bool bLeftBad; bool bLeftNaN;
			 bool bRightBad; bool bRightNaN; char cOp; unsigned long lExpect; } cases[] =
	{
		{ "inf + 5",   true,  false, false, false, '+', kfInfinity },
		{ "5 + inf",   false, false, true,  false, '+', kfInfinity },
		{ "inf - 5",   true,  false, false, false, '-', kfInfinity },
		{ "inf * 5",   true,  false, false, false, '*', kfInfinity },
		{ "inf / 5",   true,  false, false, false, '/', kfInfinity },
		{ "inf * inf", true,  false, true,  false, '*', kfInfinity },
		{ "inf + inf", true,  false, true,  false, '+', kfQNaN     },
		{ "inf - inf", true,  false, true,  false, '-', kfQNaN     },
		{ "inf / inf", true,  false, true,  false, '/', kfQNaN     },
		{ "nan + 5",   true,  true,  false, false, '+', kfQNaN     },
		{ "nan * 5",   true,  true,  false, false, '*', kfQNaN     },
		{ "5 * nan",   false, false, true,  true,  '*', kfQNaN     },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber	a, b, result;

		if ( theCase.bLeftBad )
		{
			if ( theCase.bLeftNaN )	MakeNotANumber( a );
			else					MakeInfinite( a );
		}
		else
			Set( a, "5" );

		if ( theCase.bRightBad )
		{
			if ( theCase.bRightNaN )	MakeNotANumber( b );
			else						MakeInfinite( b );
		}
		else
			Set( b, "5" );

		switch ( theCase.cOp )
		{
		case '+':	lz_add( &a, &b, &result );		break;
		case '-':	lz_subtract( &a, &b, &result );	break;
		case '*':	lz_multiply( &a, &b, &result );	break;
		default:	lz_divide( &a, &b, &result );	break;
		}

		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "carries: %s", theCase.szLabel );
		CheckTrue( szLabel, result.FlagIsOn( theCase.lExpect ), "" );
	}

	// Zero times an infinity has no answer; something finite over one is zero...
	{
		CLZNumber a, b, result;
		MakeInfinite( a ); Set( b, "0" );
		lz_multiply( &a, &b, &result );
		CheckTrue( "inf * 0 has no answer", result.FlagIsOn( kfQNaN ), "" );

		MakeInfinite( b ); Set( a, "5" );
		lz_divide( &a, &b, &result );
		CheckTrue( "5 / inf is finite", !result.FlagIsOn( kfInfinity ) &&
										!result.FlagIsOn( kfQNaN ), "" );
		CheckExact( "5 / inf is zero", &result, "0" );
	}

	// ...and it survives more than one step...
	{
		CLZNumber a, b, step, two, product, sum;
		Set( a, "1" ); Set( b, "0" ); Set( two, "2" );
		lz_divide( &a, &b, &step );
		lz_multiply( &step, &two, &product );
		lz_add( &product, &two, &sum );
		CheckTrue( "(1/0)*2 + 2 is still infinite", sum.FlagIsOn( kfInfinity ), "" );
	}

	/*
		One operand, and the rule turns on whether the two signs of an
		infinity would lead to different answers. They do not for the
		hyperbolic sine and cosine or the cube root, which are
		infinite either way. They do for the exponential -- e to the
		minus infinity is nought, not infinity -- and for the
		hyperbolic tangent, which is plus or minus one, and those give
		a not-a-number instead of choosing.
	*/

	struct { const char* szLabel; char cOp; unsigned long lExpect; } unary[] =
	{
		{ "cbrt(inf)",  '3', kfInfinity },
		{ "sinh(inf)",  'S', kfInfinity },
		{ "cosh(inf)",  'H', kfInfinity },
		{ "exp(inf)",   'E', kfQNaN     },
		{ "tanh(inf)",  'T', kfQNaN     },
		{ "sqrt(inf)",  'R', kfQNaN     },
		{ "ln(inf)",    'L', kfQNaN     },
		{ "sin(inf)",   's', kfQNaN     },
		{ "arcsin(inf)",'A', kfQNaN     },
		{ "arctan(inf)",'a', kfQNaN     },
		{ "arccos(inf)",'c', kfQNaN     },
		{ "arcsinh(inf)",'j', kfInfinity },
		{ "arccosh(inf)",'h', kfQNaN     },
		{ "arctanh(inf)",'g', kfQNaN     },
	};

	for ( auto& theCase : unary )
	{
		CLZNumber	a, result;
		MakeInfinite( a );

		switch ( theCase.cOp )
		{
		case '3':	lz_root_3( &a, &result );		break;
		case 'S':	lz_hypersine( &a, &result );	break;
		case 'H':	lz_hypercosine( &a, &result );	break;
		case 'E':	lz_exp_e( &a, &result );		break;
		case 'T':	lz_hypertangent( &a, &result );	break;
		case 'R':	lz_root_2( &a, &result );		break;
		case 'L':	lz_log_e( &a, &result );		break;
		case 's':	lz_sine( &a, &result );			break;
		case 'a':	lz_arctangent( &a, &result );	break;
		case 'c':	lz_arccosine( &a, &result );	break;
		case 'j':	lz_archypersine( &a, &result );	break;
		case 'h':	lz_archypercosine( &a, &result );	break;
		case 'g':	lz_archypertangent( &a, &result );	break;
		default:	lz_arcsine( &a, &result );		break;
		}

		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "carries: %s", theCase.szLabel );
		CheckTrue( szLabel, result.FlagIsOn( theCase.lExpect ), "" );
	}

	// A reciprocal of an infinity is a plain zero, not a failure...
	{
		CLZNumber a, result;
		MakeInfinite( a );
		lz_reciprocal( &a, &result );
		CheckTrue( "1 / inf is finite", !result.FlagIsOn( kfInfinity ) &&
										!result.FlagIsOn( kfQNaN ), "" );
	}

	// ...and a not-a-number begets one through every unary operation...
	{
		const char* szNames[] = { "sqrt", "cbrt", "exp", "ln", "sin", "cos",
								  "sinh", "cosh", "tanh", "arcsin", "arctan",
								  "arccos", "arcsinh", "arccosh",
								  "arctanh", "factorial" };

		for ( int i = 0; i < 16; i++ )
		{
			CLZNumber a, result;
			MakeNotANumber( a );

			switch ( i )
			{
			case 0:		lz_root_2( &a, &result );		break;
			case 1:		lz_root_3( &a, &result );		break;
			case 2:		lz_exp_e( &a, &result );		break;
			case 3:		lz_log_e( &a, &result );		break;
			case 4:		lz_sine( &a, &result );			break;
			case 5:		lz_cosine( &a, &result );		break;
			case 6:		lz_hypersine( &a, &result );	break;
			case 7:		lz_hypercosine( &a, &result );	break;
			case 8:		lz_hypertangent( &a, &result );	break;
			case 9:		lz_arcsine( &a, &result );		break;
			case 10:	lz_arctangent( &a, &result );	break;
			case 11:	lz_arccosine( &a, &result );	break;
			case 12:	lz_archypersine( &a, &result );	break;
			case 13:	lz_archypercosine( &a, &result );	break;
			case 14:	lz_archypertangent( &a, &result );	break;
			default:	lz_factorial( &a, &result );	break;
			}

			char szLabel[ 64 ];
			::snprintf( szLabel, 64, "nan through %s", szNames[ i ] );
			CheckTrue( szLabel, result.FlagIsOn( kfQNaN ), "" );
		}
	}

	/*
		A comparison against a flagged value has no answer. Neither
		has one against an infinity, which is held as zero with a flag
		and so might be above everything or below it.
	*/

	{
		CLZNumber a, b;

		Set( a, "5" ); Set( b, "3" );
		CheckTrue( "5 compares above 3", lz_compare( &a, &b ) == kGreaterThan, "" );
		CheckTrue( "max of 5 and 3 is the left", lz_max( &a, &b ) == kLeft, "" );
		CheckTrue( "min of 5 and 3 is the right", lz_min( &a, &b ) == kRight, "" );

		Set( a, "5" ); Set( b, "5" );
		CheckTrue( "5 compares equal to 5", lz_compare( &a, &b ) == kEqualTo, "" );
	}

	{
		CLZNumber nan, inf, five;
		MakeNotANumber( nan );
		MakeInfinite( inf );
		Set( five, "5" );

		CheckTrue( "nan against a number is unordered",
				   lz_compare( &nan, &five ) == kUnordered, "" );
		CheckTrue( "a number against nan is unordered",
				   lz_compare( &five, &nan ) == kUnordered, "" );
		CheckTrue( "nan against itself is unordered",
				   lz_compare( &nan, &nan ) == kUnordered, "" );
		CheckTrue( "an unsigned infinity is unordered",
				   lz_compare( &inf, &five ) == kUnordered, "" );
		CheckTrue( "two infinities are unordered",
				   lz_compare( &inf, &inf ) == kUnordered, "" );
		CheckTrue( "max cannot choose against a nan",
				   lz_max( &nan, &five ) == kUnordered, "" );
		CheckTrue( "min cannot choose against a nan",
				   lz_min( &nan, &five ) == kUnordered, "" );
		CheckTrue( "quantum of a nan is unordered",
				   lz_same_quantum( &nan, &five ) == kUnordered, "" );
	}

	/*
		Raising to a power keeps more than the rest, because the
		exponent's sign is readable even when the infinity's is not:
		an infinity to a positive power is infinite whichever sign it
		had, to a negative power is nought, and to the nought power is
		one. An infinite exponent is the other way about and has no
		answer -- two to the plus infinity is infinite and to the
		minus infinity is nought.

		Quantizing has none of that. It sets one value's exponent from
		another's, and neither an infinity nor a NaN has one.
	*/

	{
		struct { const char* szLabel; const char* szExponent;
				 unsigned long lWanted; bool bWantClean; const char* szValue; } powers[] =
		{
			{ "inf to the tenth",    "10", kfInfinity, false, nullptr },
			{ "inf to the second",   "2",  kfInfinity, false, nullptr },
			{ "inf to the minus 2",  "-2", 0,          true,  "0"     },
			{ "inf to the nought",   "0",  0,          true,  "1"     },
		};

		for ( auto& theCase : powers )
		{
			CLZNumber a, b, result;
			MakeInfinite( a );
			Set( b, theCase.szExponent );
			lz_power( &a, &b, &result );

			char szLabel[ 64 ];
			::snprintf( szLabel, 64, "power: %s", theCase.szLabel );

			if ( theCase.bWantClean )
			{
				CheckTrue( szLabel, !result.FlagIsOn( kfInfinity ) &&
									!result.FlagIsOn( kfQNaN ), "" );
				::snprintf( szLabel, 64, "power value: %s", theCase.szLabel );
				CheckExact( szLabel, &result, theCase.szValue );
			}
			else
				CheckTrue( szLabel, result.FlagIsOn( theCase.lWanted ), "" );
		}

		// An infinite or unreadable exponent, which has no answer...
		CLZNumber a, b, result;

		Set( a, "2" ); MakeInfinite( b );
		lz_power( &a, &b, &result );
		CheckTrue( "power: 2 to the infinite", result.FlagIsOn( kfQNaN ), "" );

		Set( a, "2" ); MakeNotANumber( b );
		lz_power( &a, &b, &result );
		CheckTrue( "power: 2 to a nan", result.FlagIsOn( kfQNaN ), "" );

		MakeNotANumber( a ); Set( b, "10" );
		lz_power( &a, &b, &result );
		CheckTrue( "power: a nan to the tenth", result.FlagIsOn( kfQNaN ), "" );

		// Quantizing, where nothing survives...
		MakeInfinite( a ); Set( b, "2" );
		lz_quantize( &a, &b, &result );
		CheckTrue( "quantize: an infinity has no exponent", result.FlagIsOn( kfQNaN ), "" );

		Set( a, "1.234" ); MakeNotANumber( b );
		lz_quantize( &a, &b, &result );
		CheckTrue( "quantize: to a nan's exponent", result.FlagIsOn( kfQNaN ), "" );

		// ...while ordinary use is untouched...
		Set( a, "2" ); Set( b, "10" );
		lz_power( &a, &b, &result );
		CheckExact( "power: 2 to the tenth", &result, "1024" );

		Set( a, "2" ); Set( b, "-3" );
		lz_power( &a, &b, &result );
		CheckExact( "power: 2 to the minus third", &result, ".125" );

		Set( a, "1.234" ); Set( b, "0.01" );
		lz_quantize( &a, &b, &result );
		CheckExact( "quantize: 1.234 to hundredths", &result, "1.23" );
	}

	/*
		The engine divides by zero itself as a matter of course, and
		before anything propagated that went unnoticed. lz_root_2
		divides by its running estimate, which starts at zero for a
		zero argument, so the square root of nought came back correct
		and flagged not-a-number. Ordinary answers must stay clean.
	*/

	struct { const char* szLabel; const char* szValue; char cOp; } clean[] =
	{
		{ "sqrt(0)",   "0", 'R' }, { "sqrt(2)", "2", 'R' },
		{ "ln(1)",     "1", 'L' }, { "sin(0)",  "0", 'S' },
		{ "cos(0)",    "0", 'C' }, { "arcsin(0)", "0", 'A' },
		{ "arcsin(1)", "1", 'A' }, { "exp(0)",  "0", 'E' },
		// Both sides of the reflection, the reciprocal being the one
		//	place the arctangent divides by something it computed...
		{ "arctan(0)", "0", 'a' }, { "arctan(1)", "1", 'a' },
		{ "arctan(5)", "5", 'a' },
		{ "arccos(0)", "0", 'c' }, { "arccos(1)", "1", 'c' },
		{ "arcsinh(0)", "0", 'j' }, { "arcsinh(2)", "2", 'j' },
		{ "arccosh(1)", "1", 'h' }, { "arccosh(5)", "5", 'h' },
		{ "arctanh(0)", "0", 'g' }, { "arctanh(.5)", ".5", 'g' },
	};

	for ( auto& theCase : clean )
	{
		CLZNumber a, result;
		Set( a, theCase.szValue );

		switch ( theCase.cOp )
		{
		case 'R':	lz_root_2( &a, &result );	break;
		case 'L':	lz_log_e( &a, &result );	break;
		case 'S':	lz_sine( &a, &result );		break;
		case 'C':	lz_cosine( &a, &result );	break;
		case 'A':	lz_arcsine( &a, &result );	break;
		case 'c':	lz_arccosine( &a, &result );	break;
		case 'a':	lz_arctangent( &a, &result );	break;
		case 'j':	lz_archypersine( &a, &result );	break;
		case 'h':	lz_archypercosine( &a, &result );	break;
		case 'g':	lz_archypertangent( &a, &result );	break;
		default:	lz_exp_e( &a, &result );	break;
		}

		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "no stray flag: %s", theCase.szLabel );
		CheckTrue( szLabel, !result.FlagIsOn( kfQNaN ) &&
							!result.FlagIsOn( kfInfinity ), "" );
	}
}

/////////////////////////////////////////////////////////////////
//	Raising to a power. lz_raise_n multiplied the operand by itself
//	once per unit of the exponent, and lz_raise_e carried its own
//	copy of the same loop, so e to the ten thousandth took half a
//	million operations and anything larger never came back. Squaring
//	makes it one per binary digit of the exponent.
/////////////////////////////////////////////////////////////////

const char kRef2To1000[] =
	"10715086071862673209484250490600018105614048117055336074437503883703510511249"
	"36122493198378815695858127594672917553146825187145285692314043598457757469857"
	"48039345677748242309854210746050623711418779541821530464749835819412673987675"
	"59165543946077062914571196477686542167660429831652624386837205668069376";

static void TestRaising()
{
	Banner( "Raising to a power" );

	// Exact powers, where every digit of the answer is checkable...
	struct { const char* szBase; const char* szExponent; const char* szExpect; } exact[] =
	{
		{ "2",  "10",   "1024"                  },
		{ "3",  "5",    "243"                   },
		{ "7",  "0",    "1"                     },
		{ "2",  "64",   "18446744073709551616"  },
		{ "3",  "40",   "12157665459056928801"  },
		{ "2",  "-10",  ".0009765625"           },
		{ "10", "20",   "100000000000000000000" },
		{ "2",  "1000", kRef2To1000             },
	};

	for ( auto& theCase : exact )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "%s to the %s", theCase.szBase, theCase.szExponent );

		CLZNumber a, b, result;
		Set( a, theCase.szBase );
		Set( b, theCase.szExponent );
		lz_raise_n( &a, &b, &result );
		CheckExact( szLabel, &result, theCase.szExpect );
	}

	/*
		Products either side of the point where multiplication stops
		doing it the plain way and starts splitting. Thirty-one digits
		is below the threshold and thirty-three above it, so these
		cover both paths and the boundary between them. All exact, so
		every digit is checked.
	*/

	struct { const char* szLabel; const char* szLeft; const char* szRight;
			 const char* szExpect; } products[] =
	{
		{ "31 nines by 31 eights",
		  "9999999999999999999999999999999",
		  "8888888888888888888888888888888",
		  "88888888888888888888888888888871111111111111111111111111111112" },

		{ "33 sevens by 33 sixes",
		  "777777777777777777777777777777777",
		  "666666666666666666666666666666666",
		  "518518518518518518518518518518517481481481481481481481481481481482" },

		{ "64 ones squared",
		  "1111111111111111111111111111111111111111111111111111111111111111",
		  "1111111111111111111111111111111111111111111111111111111111111111",
		  "1234567901234567901234567901234567901234567901234567901234567900"
		  "987654320987654320987654320987654320987654320987654320987654321" },
	};

	for ( auto& theCase : products )
	{
		CLZNumber a, b, result;
		Set( a, theCase.szLeft );
		Set( b, theCase.szRight );
		lz_multiply( &a, &b, &result );
		CheckExact( theCase.szLabel, &result, theCase.szExpect );
	}

	// The short way still has to work, and signs and fractions with it...
	{
		CLZNumber a, b, result;

		Set( a, "1.5" );   Set( b, "2.5" );
		lz_multiply( &a, &b, &result );
		CheckExact( "1.5 by 2.5", &result, "3.75" );

		Set( a, "-12" );   Set( b, "12" );
		lz_multiply( &a, &b, &result );
		CheckExact( "-12 by 12", &result, "-144" );

		Set( a, "0" );     Set( b, "999999999999999999999999999999999999" );
		lz_multiply( &a, &b, &result );
		CheckExact( "nought by a long number", &result, "0" );
	}

	// An exponent written in another base counts the same...
	{
		CLZNumber a, b, result;
		Set( a, "2" );
		SetInBase( b, 16, "40" );			// sixty-four, in hex
		lz_raise_n( &a, &b, &result );
		CheckExact( "2 to the 40 base 16", &result, "18446744073709551616" );

		SetInBase( b, 2, "1010" );			// ten, in binary
		lz_raise_n( &a, &b, &result );
		CheckExact( "2 to the 1010 base 2", &result, "1024" );
	}

	// ...and the exponential, which used to walk there one step at a
	//	time. A hundred was five thousand operations and is now five
	//	hundred; a thousand no longer has to be avoided...
	{
		CLZNumber a, result;
		Set( a, "100" );
		lz_exp_e( &a, &result );
		CheckDigits( "e to the hundredth", &result,
					 "26881171418161354484126255515800135873611118.773741", 48 );

		Set( a, "1000" );
		lz_exp_e( &a, &result );
		CheckTrue( "e to the thousandth has 435 digits",
				   result.m_lDigitsLength - result.m_lRadix == 435, "" );
	}

	/*
		Powers, roots and logarithms that widen for their answer. Six
		functions formed a large answer at relative precision and gave
		it 36 to 39 significant digits whatever its size -- lz_power( 10,
		12.5 ) was 2.8e11 units short of its 36th place -- and none was
		in the sweep or the audit to show it. Each now works out the
		answer's magnitude at a dozen places first and runs that much
		wider, as lz_raise_e does; see lzi_widened_by_answer.

		Every expected string is Python's value correctly rounded at 36
		places, and every one was wrong before by the units noted.
		CheckRendering rather than CheckDigits: the claim is the 36th
		place, not a count of digits.
	*/

	{
		typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );

		struct { const char* szLabel; Binary fn; const char* szLeft; const char* szRight; const char* szExpect; } widened[] =
		{
			{ "10 to the 12.5 by lz_power",  lz_power,   "10",  "12.5", "3162277660168.379331998893544432718533719555139325" },		// was 2.8e11 units
			{ "3 to the 20.7 by lz_power",   lz_power,   "3",   "20.7", "7523327587.932312895850435922591876293043296360" },			// 4.6e9
			{ "3 to the 29.5 by lz_exp_n",   lz_exp_n,   "3",   "29.5", "118871300538602.396939682955758993971684854288873779" },		// 6.5e13
			{ "the 2.5th root of 1e40",      lz_root_n,  "10000000000000000000000000000000000000000", "2.5",
			  "10000000000000000.000000000000000000000000000000000000" },																// 3.2e11
			{ "the 2.5th root of 123456789", lz_root_n,  "123456789", "2.5", "1724.272853629601966444928899118619877014" },			// 197
			{ "1.1 to the 300th",            lz_raise_n, "1.1", "300",  "2617010996188.399907017032528972038342491649416953" },		// 1.8e9
			{ "1.0001 to the 100000th",      lz_raise_n, "1.0001", "100000", "22015.456048552198645701456581658715520598" },		// the guard: 550 without it
			{ "log to base 3 of 1e30",       lz_log_n,   "3",   "1000000000000000000000000000000",
			  "62.877098228681538128897025660642037518" },																			// 20
		};

		for ( auto& theCase : widened )
		{
			CLZNumber	left, right, answer;

			Set( left, theCase.szLeft );
			Set( right, theCase.szRight );
			theCase.fn( &left, &right, &answer );
			CheckRendering( theCase.szLabel, &answer, theCase.szExpect );
		}

		CLZNumber	a, answer;

		Set( a, "20.5" );
		lz_exp_b( &a, &answer );
		CheckRendering( "10 to the 20.5 by lz_exp_b", &answer,
						"316227766016837933199.889354443271853371955513932521682686" );		// 49 digits of 57, the constant's width

		Set( a, "2.5" );
		lz_exp_b( &a, &answer );
		CheckRendering( "10 to the 2.5 by lz_exp_b", &answer, "316.227766016837933199889354443271853372" );	// 28 units
	}
}

/////////////////////////////////////////////////////////////////
//	The edges of x^y. lz_power, lz_exp_n and lz_raise_n met them each
//	their own way and agreed at none of them until 2026-09-26: a
//	negative base to a fraction was a flagged not-a-number reporting
//	success, a refusal, or a truncated exponent -- (-8)^0.5 was 1 from
//	lz_raise_n; nought to a negative power was infinite reporting
//	success, or refused; and lz_exp_n refused every negative and zero
//	base, (-8)^3 and 0^2 among them. lz_raise_e and lz_raise_n dropped
//	a fraction and reported success, so 2^2.5 was 4 and e^2.5 was e^2.
//	One set of rules now -- see lzi_power_edges.
/////////////////////////////////////////////////////////////////

static void TestPowerEdges()
{
	Banner( "The edges of x^y" );

	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );

	// Held as a zero with a flag, as the engine holds them...
	auto Make = []( CLZNumber& theNumber, const char* szValue )
	{
		if ( ::strcmp( szValue, "inf" ) == 0 )
		{
			theNumber.Assign( 0L );
			theNumber.SetFlag( kfInfinity );
		}
		else if ( ::strcmp( szValue, "nan" ) == 0 )
		{
			theNumber.Assign( 0L );
			theNumber.SetFlag( kfQNaN );
		}
		else
			theNumber.Assign( szValue );
	};

	auto Flags = []( CLZNumber& theNumber ) -> int
	{
		return (theNumber.FlagIsOn( kfQNaN ) ? 1 : 0) | (theNumber.FlagIsOn( kfInfinity ) ? 2 : 0) |
			   (theNumber.FlagIsOn( kfDiv0 ) ? 4 : 0);
	};

	const int	kNaN = 1;
	const int	kInf = 2;
	const int	kDiv = 4;

	/*
		Absolute answers first, since three functions agreeing proves
		nothing if they agree on the wrong thing. Every row here was
		wrong in at least one of them before.
	*/

	struct { const char* szLabel; Binary fn; const char* szBase; const char* szExponent;
			 long lCode; const char* szExpect; int iFlags; } pins[] =
	{
		{ "power (-8)^3",           lz_power,   "-8", "3",    kNoError,      "-512",   0 },
		{ "exp_n (-8)^3",           lz_exp_n,   "-8", "3",    kNoError,      "-512",   0 },		// was refused
		{ "exp_n (-2)^-3",          lz_exp_n,   "-2", "-3",   kNoError,      "-.125",  0 },		// was refused
		{ "power (-8)^2.5 refused", lz_power,   "-8", "2.5",  kFailed,       "0",      kNaN },	// was success
		{ "exp_n (-8)^2.5 refused", lz_exp_n,   "-8", "2.5",  kFailed,       "0",      kNaN },
		{ "raise_n (-8)^.5",        lz_raise_n, "-8", "0.5",  kFailed,       "0",      kNaN },	// was 1
		{ "raise_n 2^2.5",          lz_raise_n, "2",  "2.5",  kFailed,       "0",      kNaN },	// was 4
		{ "power 0^2.5",            lz_power,   "0",  "2.5",  kNoError,      "0",      0 },		// was not-a-number
		{ "exp_n 0^2",              lz_exp_n,   "0",  "2",    kNoError,      "0",      0 },		// was refused
		{ "exp_n 0^0",              lz_exp_n,   "0",  "0",    kNoError,      "1",      0 },		// was refused
		{ "power 0^-1, the pole",   lz_power,   "0",  "-1",   kDivideByZero, "0",      kInf | kDiv },	// was success
		{ "exp_n 0^-1, the pole",   lz_exp_n,   "0",  "-1",   kDivideByZero, "0",      kInf | kDiv },	// was 4, not-a-number
		{ "raise_n 0^-2, the pole", lz_raise_n, "0",  "-2",   kDivideByZero, "0",      kInf | kDiv },	// was success
		{ "raise_n 2^nan",          lz_raise_n, "2",  "nan",  kNoError,      "0",      kNaN },	// was 1
		{ "raise_n 2^inf",          lz_raise_n, "2",  "inf",  kNoError,      "0",      kNaN },	// was 1
		{ "exp_n inf^-2",           lz_exp_n,   "inf", "-2",  kNoError,      "0",      0 },		// was not-a-number
		{ "exp_n inf^0",            lz_exp_n,   "inf", "0",   kNoError,      "1",      0 },		// was not-a-number
	};

	for ( auto& thePin : pins )
	{
		CLZNumber	x, y, target;

		Make( x, thePin.szBase );
		Make( y, thePin.szExponent );
		target.Assign( 12345L );

		long	l_code = thePin.fn( &x, &y, &target );
		char*	szGot = Render( &target );
		bool	b_right = (l_code == thePin.lCode) && (Flags( target ) == thePin.iFlags) &&
						  (szGot != nullptr) && (AgreementDigits( szGot, thePin.szExpect ) == kExactAgreement);
		char	szDetail[ 96 ];

		::snprintf( szDetail, 96, "code %ld, %s, flags %d", l_code, (szGot != nullptr) ? szGot : "(null)", Flags( target ) );
		CheckTrue( thePin.szLabel, b_right, szDetail );

		::free( szGot );
	}

	// The pole is the reciprocal's, code and flags, by construction...
	{
		CLZNumber	one, zero, reciprocal, x, y, target;

		one.Assign( 1L );
		zero.Assign( 0L );

		long	l_divide = lz_divide( &one, &zero, &reciprocal );

		x.Assign( 0L );
		y.Assign( "-3" );

		long	l_power = lz_power( &x, &y, &target );

		CheckTrue( "0^-3 answers as 1/0 does", (l_power == l_divide) && (Flags( target ) == Flags( reciprocal )), "" );
	}

	// In place, where the base is the target...
	{
		CLZNumber	x, y;

		x.Assign( 0L );
		y.Assign( "-1" );

		long	l_code = lz_power( &x, &y, nullptr );

		CheckTrue( "0^-1 in place", (l_code == kDivideByZero) && x.FlagIsOn( kfInfinity ) && x.IsZero(),
				   "the base becomes a flagged nought" );

		x.Assign( 2L );
		y.Assign( "2.5" );
		l_code = lz_raise_n( &x, &y, nullptr );

		CheckTrue( "raise_n 2^2.5 in place", (l_code == kFailed) && x.FlagIsOn( kfQNaN ) && x.IsZero(),
				   "the base becomes a flagged nought" );
	}

	// A negative base in another base, and an exponent in a third...
	{
		CLZNumber	x, y, target;

		SetInBase( x, 16, "-8" );
		SetInBase( y, 3, "10" );			// three
		lz_exp_n( &x, &y, &target );
		CheckRendering( "exp_n (-8)^3 in hex", &target, "-200" );
	}

	/*
		Then agreement, over every pairing of these. lz_power and
		lz_exp_n must agree on the code, the value and the flags
		everywhere. lz_raise_n must agree with them wherever the exponent
		is whole, and refuse the rest -- unless an operand is a
		not-a-number or an infinity, which is carried before the
		exponent is looked at.
	*/

	{
		const char*	bases[] = { "-8", "-2", "-1", "-0.5", "0", "0.5", "1", "2", "inf", "nan" };
		const char*	exponents[] = { "-2.5", "-2", "-1", "0", "0.5", "1", "2", "2.5", "3", "inf", "nan" };
		long		l_disagree = 0;
		long		l_cells = 0;
		char		szWhere[ 64 ] = "";

		for ( const char* szBase : bases )
		{
			for ( const char* szExponent : exponents )
			{
				CLZNumber	answers[ 3 ];
				long		l_codes[ 3 ];
				Binary		fns[ 3 ] = { lz_power, lz_exp_n, lz_raise_n };

				for ( int i = 0; i < 3; i++ )
				{
					CLZNumber	x, y;

					Make( x, szBase );
					Make( y, szExponent );
					answers[ i ].Assign( 12345L );
					l_codes[ i ] = fns[ i ]( &x, &y, &answers[ i ] );
				}

				CLZNumber	y;

				Make( y, szExponent );

				bool	b_carried = (::strcmp( szBase, "inf" ) == 0) || (::strcmp( szBase, "nan" ) == 0) ||
									(::strcmp( szExponent, "inf" ) == 0) || (::strcmp( szExponent, "nan" ) == 0);
				bool	b_whole = b_carried || y.IsInteger();

				// A flagged answer is a flagged nought, and lz_compare rightly
				//	calls it unordered, so values are compared only unflagged...
				auto Same = [&]( int i, int j ) -> bool
				{
					return (l_codes[ i ] == l_codes[ j ]) && (Flags( answers[ i ] ) == Flags( answers[ j ] )) &&
						   ((Flags( answers[ i ] ) != 0) || (lz_compare( &answers[ i ], &answers[ j ] ) == kEqualTo));
				};

				bool	b_agree = Same( 0, 1 );

				if ( b_whole )
					b_agree = b_agree && Same( 0, 2 );
				else
					b_agree = b_agree && (l_codes[ 2 ] == kFailed) && (Flags( answers[ 2 ] ) == kNaN);

				l_cells++;

				if ( !b_agree )
				{
					l_disagree++;
					::snprintf( szWhere, 64, "%s^%s", szBase, szExponent );
				}
			}
		}

		char	szDetail[ 96 ];

		if ( l_disagree == 0 )
			::snprintf( szDetail, 96, "all %ld pairings", l_cells );
		else
			::snprintf( szDetail, 96, "%ld of %ld disagree, one at %s", l_disagree, l_cells, szWhere );

		CheckTrue( "power, exp_n and raise_n agree", l_disagree == 0, szDetail );
	}

	// lz_raise_e: a fraction refused, flags carried as lz_exp_e carries them...
	{
		struct { const char* szValue; long lCode; int iFlags; } cases[] =
		{
			{ "2.5",  kFailed,  kNaN },		// was e^2
			{ "-0.5", kFailed,  kNaN },		// was 1
			{ "inf",  kNoError, kNaN },		// was 1, unflagged
			{ "nan",  kNoError, kNaN },		// was 1, unflagged
		};

		for ( auto& theCase : cases )
		{
			CLZNumber	x, target, xExp, targetExp;

			Make( x, theCase.szValue );
			target.Assign( 12345L );

			long	l_code = lz_raise_e( &x, &target );

			char	szLabel[ 64 ];
			char	szDetail[ 64 ];

			::snprintf( szLabel, 64, "raise_e %s", theCase.szValue );
			::snprintf( szDetail, 64, "code %ld, flags %d", l_code, Flags( target ) );
			CheckTrue( szLabel, (l_code == theCase.lCode) && (Flags( target ) == theCase.iFlags) && target.IsZero(), szDetail );

			// ...and where lz_exp_e answers, a flagged operand, the two agree...
			if ( theCase.lCode == kNoError )
			{
				Make( xExp, theCase.szValue );
				lz_exp_e( &xExp, &targetExp );

				::snprintf( szLabel, 64, "raise_e %s as exp_e", theCase.szValue );
				CheckTrue( szLabel, Flags( target ) == Flags( targetExp ), "" );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////
//	How a value reads back. A leading zero is one with no digit
//	written before it, and the test for that used to ask whether
//	anything had been written at all -- which a minus sign answers,
//	so every negative value kept its leading zeros and sixteen less
//	a hundred read -084.
/////////////////////////////////////////////////////////////////

static void TestRendering()
{
	Banner( "Reading a value back" );

	struct { const char* szValue; const char* szExpect; } written[] =
	{
		{ "0",      "0"      },  { "-0",     "0"      },
		{ "1",      "1"      },  { "-1",     "-1"     },
		{ "100",    "100"    },  { "-100",   "-100"   },
		{ "0.5",    "0.5"    },  { "-0.5",   "-0.5"   },
		{ "0.05",   "0.05"   },  { "-0.05",  "-0.05"  },
		{ "10.01",  "10.01"  },  { "-10.01", "-10.01" },
	};

	for ( auto& theCase : written )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "%s reads back", theCase.szValue );

		CLZNumber a;
		Set( a, theCase.szValue );
		CheckRendering( szLabel, &a, theCase.szExpect );
	}

	// Computed negatives, where the answer is narrower than the
	//	working room it was given...
	struct { const char* szLeft; const char* szRight; const char* szExpect; } sums[] =
	{
		{ "16",  "100",  "-84"   },  { "7",    "16",   "-9"    },
		{ "1",   "1000", "-999"  },  { "0.5",  "100",  "-99.5" },
		{ "100", "200",  "-100"  },  { "1",    "1.05", "-0.05" },
		{ "5",   "5",    "0"     },  { "100",  "16",   "84"    },
	};

	for ( auto& theCase : sums )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "%s - %s", theCase.szLeft, theCase.szRight );

		CLZNumber a, b, result;
		Set( a, theCase.szLeft );
		Set( b, theCase.szRight );
		lz_subtract( &a, &b, &result );
		CheckRendering( szLabel, &result, theCase.szExpect );
	}

	// The one that led here: 255 less the nearest multiple of 16...
	{
		CLZNumber a, b, result;
		Set( a, "255" ); Set( b, "16" );
		lz_remainder_near( &a, &b, &result );
		CheckRendering( "255 remainder-near 16", &result, "-1" );
	}

	/*
		A value with nothing above the radix mark still shows the
		nought in front of it. Every digit up there was suppressed as a
		leading zero, which is right for the ones further up and wrong
		for the one against the mark.

		Most of these hold no integer digit at all -- m_lRadix equals
		m_lDigitsLength -- which is the ordinary shape of any fraction
		the engine computes, not an edge case. It is also what the
		buffer's fourth byte is for: at three the terminator goes one
		past the end, which ASan reports and the plain build does not.
	*/

	Banner( "The nought in front of the mark" );

	struct { const char* szValue; const char* szExpect; } fronted[] =
	{
		{ ".5",     "0.5"    },  { "-.5",    "-0.5"   },
		{ ".05",    "0.05"   },  { "-.05",   "-0.05"  },
		{ "0.5",    "0.5"    },  { ".000",   "0.000"  },
		{ ".123456789", "0.123456789" },
	};

	for ( auto& theCase : fronted )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "%s keeps its nought", theCase.szValue );

		CLZNumber a;
		Set( a, theCase.szValue );
		CheckRendering( szLabel, &a, theCase.szExpect );
	}

	// Computed fractions, which is where the shape actually comes from...
	{
		CLZNumber a, b, result;

		Set( a, "1" ); Set( b, "2" );
		lz_divide( &a, &b, &result );
		CheckRendering( "1 / 2", &result, "0.5" );

		Set( a, "0.3" ); Set( b, "0.3" );
		lz_multiply( &a, &b, &result );
		CheckRendering( "0.3 * 0.3", &result, "0.09" );

		Set( a, "12.34" );
		a.StripInteger();
		CheckRendering( "12.34 less its integer part", &a, "0.34" );
	}

	// And in a base where the digits are not decimal...
	for ( long lBase = 2; lBase <= 16; lBase += 2 )
	{
		CLZNumber	a;
		char		szLabel[ 64 ];

		a.SetBase( lBase, true );
		a.Assign( ".1" );

		::snprintf( szLabel, 64, "base %ld .1 keeps its nought", lBase );
		CheckRendering( szLabel, &a, "0.1" );
	}
}

static void TestDisplay()
{
	/*
		The three display modes at a named width. FIX counts places
		after the mark; SCI and ENG count digits after the leading one,
		so all three of these show a value to lDigits+1 figures in the
		exponent styles and to whatever the value's size makes of
		lDigits places in FIX.
	*/

	Banner( "FIX, SCI and ENG" );

	struct { const char* szText; long lDigits;
			 const char* szFix; const char* szSci; const char* szEng; } shown[] =
	{
		// The ordinary case, and the one that makes the point: what is
		//	shown is not what is held...
		{ "2.3456",    2, "2.35",         "2.35E0",     "2.35E0"    },

		// A value with no integer part gets the nought in front of the
		//	mark that GetStringRep leaves off...
		{ "0.5",       2, "0.50",         "5.00E-1",    "500E-3"    },
		{ "-0.5",      2, "-0.50",        "-5.00E-1",   "-500E-3"   },

		// Short of the width, and padded out to it...
		{ "1",         4, "1.0000",       "1.0000E0",   "1.0000E0"  },
		{ "12345",     2, "12345.00",     "1.23E4",     "12.3E3"    },

		// Nought, which has no leading digit to measure from...
		{ "0",         2, "0.00",         "0.00E0",     "0.00E0"    },
		{ "0",         0, "0",            "0E0",        "0E0"       },

		// Rounding that carries all the way out of the leading digit,
		//	which moves the exponent with it...
		{ "9.99",      1, "10.0",         "1.0E1",      "10E0"      },
		{ "0.999",     2, "1.00",         "9.99E-1",    "999E-3"    },
		{ "999.9",     0, "1000",         "1E3",        "1E3"       },

		// Zero width leaves the radix mark off altogether...
		{ "123.456",   0, "123",          "1E2",        "100E0"     },

		// Sizes FIX has no useful answer for, which is the caller's
		//	problem and not the renderer's -- see the note on width...
		{ "0.00012",   2, "0.00",         "1.20E-4",    "120E-6"    },
		{ "602000000000000000000000", 2,
		  "602000000000000000000000.00", "6.02E23",     "602E21"    },

		{ "-1234.5678", 3, "-1234.568",   "-1.235E3",   "-1.235E3"  },
		{ "0.0001234",  4, "0.0001",      "1.2340E-4",  "123.40E-6" },
	};

	for ( auto& theCase : shown )
	{
		CLZNumber	a;
		char		szLabel[ 80 ];

		a.Assign( theCase.szText );

		::snprintf( szLabel, 80, "%s FIX %ld", theCase.szText, theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_positional, theCase.lDigits, theCase.szFix );

		::snprintf( szLabel, 80, "%s SCI %ld", theCase.szText, theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_scientific, theCase.lDigits, theCase.szSci );

		::snprintf( szLabel, 80, "%s ENG %ld", theCase.szText, theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_engineering, theCase.lDigits, theCase.szEng );
	}

	/*
		A value too small to show at the width does not keep its sign.
		Rounding is what takes it to nothing and CheckForZero clears the
		sign as it goes, so the display reads 0.000 and never -0.000 --
		which is the answer a calculator gives and the one that reads
		back as the value shown.
	*/

	{
		CLZNumber	a;

		a.Assign( "-0.0004" );
		CheckDisplay( "a small negative shows unsigned", &a, RS_positional, 3, "0.000" );
		CheckDisplay( "and keeps its sign in SCI", &a, RS_scientific, 3, "-4.000E-4" );
	}

	/*
		The engineering form pads its group, so it shows lDigits+1
		figures where the grouping allows and more where it does not:
		one significant digit of 123456 still has to fill three places
		in front of the mark before the exponent can be a multiple of
		three. The figures beyond the first are noughts, and say so.
	*/

	Banner( "The engineering group" );

	struct { long lDigits; const char* szEng; const char* szSci; } grouped[] =
	{
		{ 0, "100E3",   "1E5"      },
		{ 1, "120E3",   "1.2E5"    },
		{ 2, "123E3",   "1.23E5"   },
		{ 3, "123.5E3", "1.235E5"  },
		{ 4, "123.46E3","1.2346E5" },
	};

	for ( auto& theCase : grouped )
	{
		CLZNumber	a;
		char		szLabel[ 72 ];

		a.Assign( "123456" );

		::snprintf( szLabel, 72, "123456 ENG %ld", theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_engineering, theCase.lDigits, theCase.szEng );

		::snprintf( szLabel, 72, "123456 SCI %ld", theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_scientific, theCase.lDigits, theCase.szSci );
	}

	/*
		Every base, with the exponent mark chosen the same way the
		reader chooses it and the exponent counting the number's own
		base -- DEADBEEF to four hexadecimal figures is D.EAEX7.
	*/

	Banner( "A display in every base" );

	struct { long lBase; const char* szText; long lDigits;
			 const char* szFix; const char* szSci; const char* szEng; } based[] =
	{
		{ 16, "DEADBEEF", 3, "DEADBEEF.000", "D.EAEX7",  "DE.AEX6"  },
		{ 16, "A.B",      2, "A.B0",         "A.B0X0",   "A.B0X0"   },
		{ 12, "DE",       2, "DE.00",        "D.E0X1",   "DE.0X0"   },
		{  2, "1010.11",  4, "1010.1100",    "1.0110E3", "1.0110E3" },
		{  8, "777",      1, "777.0",        "1.0E3",    "1.0E3"    },
	};

	for ( auto& theCase : based )
	{
		CLZNumber	a;
		char		szLabel[ 80 ];

		a.SetBase( theCase.lBase, true );
		a.Assign( theCase.szText );

		::snprintf( szLabel, 80, "base %ld %s FIX %ld",
					theCase.lBase, theCase.szText, theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_positional, theCase.lDigits, theCase.szFix );

		::snprintf( szLabel, 80, "base %ld %s SCI %ld",
					theCase.lBase, theCase.szText, theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_scientific, theCase.lDigits, theCase.szSci );

		::snprintf( szLabel, 80, "base %ld %s ENG %ld",
					theCase.lBase, theCase.szText, theCase.lDigits );
		CheckDisplay( szLabel, &a, RS_engineering, theCase.lDigits, theCase.szEng );
	}

	/*
		The property behind all of it: a display reads back as the value
		it shows. Written, read, and written again at the same width has
		to give the same characters, or the display is saying something
		the reader does not hear -- which is how the exponent mark and
		the base prefix would go wrong without anyone noticing.
	*/

	Banner( "A display reads back as itself" );

	struct { long lBase; const char* szText; } trips[] =
	{
		{ 10, "602000000000000000000000" }, { 10, "-98765.4321" },
		{ 10, "1e-40" },  { 10, ".5" },     { 10, "0" },
		{ 10, "3.14159265358979323846" },   { 10, "-0.0004" },
		{ 16, "DEADBEEF" }, { 16, "-A.B" }, { 12, "DE" },
		{  2, "1.01" },     {  8, "-7.7" },
	};

	for ( auto& theTrip : trips )
	{
		for ( int iStyle = 0; iStyle < 3; iStyle++ )
		{
			static const char*	szStyleName[] = { "FIX", "SCI", "ENG" };

			RenderStyle	rsStyle = static_cast<RenderStyle>( iStyle );

			for ( long lDigits = 0; lDigits <= 6; lDigits += 3 )
			{
				CLZNumber	a, b;
				char		szLabel[ 96 ];

				a.SetBase( theTrip.lBase, true );
				a.Assign( theTrip.szText );

				char*	szFirst = a.GetDisplayRep( rsStyle, lDigits );

				b.SetBase( theTrip.lBase, true );

				long	lReturned = b.Assign( szFirst );
				char*	szSecond = b.GetDisplayRep( rsStyle, lDigits );

				::snprintf( szLabel, 96, "base %ld %s %s %ld",
							theTrip.lBase, theTrip.szText,
							szStyleName[ iStyle ], lDigits );

				CheckTrue( szLabel,
						   (lReturned == kNoError) && (szFirst != nullptr) &&
						   (szSecond != nullptr) && (::strcmp( szFirst, szSecond ) == 0),
						   szFirst ? szFirst : "(null)" );

				::free( szFirst );
				::free( szSecond );
			}
		}
	}

	/*
		And the guarantee the whole thing rests on: displaying a value
		does not change it. A calculator in FIX 2 shows 2.35 and goes on
		multiplying by 2.3456. A constant is displayable for the same
		reason -- the rounding never touches the original, so there is
		nothing for the constant flag to protect.
	*/

	Banner( "A display does not change the value" );

	{
		CLZNumber	a;

		a.Assign( "2.3456" );
		a.SetConstant( true );

		char*	szDisplay = a.GetDisplayRep( RS_positional, 2 );

		CheckTrue( "a constant can be displayed",
				   (szDisplay != nullptr) && (::strcmp( szDisplay, "2.35" ) == 0),
				   szDisplay ? szDisplay : "(null)" );

		::free( szDisplay );

		CheckTrue( "and is still a constant", a.IsConstant(), "" );

		a.SetConstant( false );
		CheckRendering( "and still holds every digit", &a, "2.3456" );
	}

	// A width below nothing is nothing, not a walk backwards off the
	//	front of the array...
	{
		CLZNumber	a;

		a.Assign( "3.7" );
		CheckDisplay( "a negative width shows no places", &a, RS_positional, -3, "4" );
	}
}

static void TestGamma()
{
	CLZNumber	a, result;

	/*
		The gamma function, which continues the factorial to every
		argument but its poles.

		lz_factorial is deliberately unchanged and still refuses
		anything that is not a non-negative whole number; the
		continuation is a different function and is asked for by name.
		The first check here is that the refusal still stands.
	*/

	Banner( "The gamma function" );

	Set( a, "2.5" ); result.Clear();
	CheckTrue( "factorial still refuses 2.5",
			   lz_factorial( &a, &result ) == kFailed, "" );

	struct { const char* szArg; const char* szRef; int iRequired; } cases[] =
	{
		{ "0.5",   kRefGammaHalf,        35 },
		{ "1.5",   kRefGammaOneHalf,     35 },
		{ "2.5",   kRefGammaTwoHalf,     36 },
		{ "-0.5",  kRefGammaMinusHalf,   35 },
		{ "-1.5",  kRefGammaMinusOneHalf, 35 },
	};

	for ( auto& theCase : cases )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "gamma(%s)", theCase.szArg );

		Set( a, theCase.szArg );
		result.Clear();
		lz_gamma( &a, &result );
		CheckDigits( szLabel, &result, theCase.szRef, theCase.iRequired );
	}

	/*
		A whole number is a factorial and comes back exact, going to
		lz_factorial rather than through the series -- the same
		reasoning that has lz_root_3 answer the cube root of 8 as 2.

		**Checked by rendering, not by CheckExact.** The series gets
		these right to the last place and then writes them as
		24.000000000000000000000000000000000000, which CheckExact
		normalises to 24 and passes -- a mutation disabling the
		shortcut altogether was caught by nothing until these became
		CheckRendering. The same normalising blindness hid the
		multi-digit nought for as long as it existed.

		gamma(20) is in the list because it is where the series stops
		being merely untidy and starts being wrong: it ends ...004
		where the factorial ends ...000.
	*/

	struct { const char* szArg; const char* szExpect; } whole[] =
	{
		{ "1", "1" }, { "2", "1" }, { "3", "2" }, { "5", "24" },
		{ "10", "362880" }, { "13", "479001600" },
		{ "20", "121645100408832000" },
	};

	for ( auto& theCase : whole )
	{
		char szLabel[ 72 ];
		::snprintf( szLabel, 72, "gamma(%s) is (n-1)! exactly", theCase.szArg );

		Set( a, theCase.szArg );
		result.Clear();
		lz_gamma( &a, &result );
		CheckRendering( szLabel, &result, theCase.szExpect );
	}

	/*
		The poles, at nought and every negative whole number. Not the
		usual "outside the domain": gamma runs to plus infinity on one
		side of each and minus infinity on the other, so there is no
		answer at all, not even an infinity with a sign.
	*/

	struct { const char* szArg; } poles[] =
	{
		{ "0" }, { "-1" }, { "-2" }, { "-10" }, { "-100" },
	};

	for ( auto& thePole : poles )
	{
		char szLabel[ 64 ];
		::snprintf( szLabel, 64, "gamma(%s) is a pole", thePole.szArg );

		Set( a, thePole.szArg );
		result.Clear();

		long lReturned = lz_gamma( &a, &result );

		CheckTrue( szLabel, (lReturned == kUndefinedResult) &&
							result.FlagIsOn( kfQNaN ), "" );
	}

	/*
		**A large gamma needs more places than the answer's own width,
		and getting them was the whole difficulty.**

		exp turns an absolute error in the logarithm into a relative
		one in the answer, so however many places lnGamma carries are
		the significant digits gamma gets -- and 36 places of an answer
		near 1e17 has room for 54. Gamma(20.5) held 37 until lz_gamma
		learned to widen itself by the answer's own magnitude.
	*/

	Set( a, "10.5" ); result.Clear();
	lz_gamma( &a, &result );
	CheckDigits( "gamma(10.5)", &result, kRefGammaTenHalf, 41 );

	// ...and 52 from then until the Bernoulli numbers went past their
	//	fifty-digit table, on 2026-09-26; 55 since, correctly rounded at
	//	its 36 places...
	Set( a, "20.5" ); result.Clear();
	lz_gamma( &a, &result );
	CheckDigits( "gamma(20.5)", &result, kRefGammaTwentyHalf, 54 );

	/*
		The logarithm on its own, which is what a caller wants for an
		argument whose gamma would be enormous. Refused at nought and
		below rather than answering ln|Gamma|: gamma alternates sign
		there, so a real logarithm of it does not exist, and the
		absolute value is a different function.
	*/

	Set( a, "-0.5" ); result.Clear();
	CheckTrue( "log gamma refuses a negative",
			   lz_log_gamma( &a, &result ) == kFailed, "" );

	Set( a, "0" ); result.Clear();
	CheckTrue( "log gamma refuses nought",
			   lz_log_gamma( &a, &result ) == kFailed, "" );

	/*
		**The exponential's own magnitude.**

		e^100 is 2.7e43, so 36 places has room for 80 significant
		digits -- and the squaring chain delivered 50, because
		narrowing e to the working width fixes its *relative* accuracy
		and squaring preserves that however large the answer grows. The
		shortfall was the same at every precision, which is the tell:
		7 digits at e^50, 30 at e^100, 73 at e^200.

		**Neither gate could see it.** The sweep's references are
		generated to 50 significant digits, so it cannot score above
		that and sinh(-50) was already recorded at its ceiling; the
		accuracy suite had nothing above exp(1). Only auditrounding.py,
		which computes its own references, put these five in a bucket
		-- and only once that bucket started listing its contents.

		These are here so a gate in the build sees them too.
	*/

	Banner( "The exponential at scale" );

	Set( a, "50" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(50)", &result, kRefExp50, 57 );

	Set( a, "100" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(100)", &result, kRefExp100, 79 );

	Set( a, "13.94" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(13.94)", &result, kRefExp1394, 42 );

	/*
		A fraction on top of a large integer part is what separates the
		two widenings. lz_raise_e getting e^60 right is not enough --
		the fractional factor is multiplied into it and caps the
		product at its own width, so lz_exp_e has to widen for the
		answer's magnitude before summing it.

		exp(13.94) cannot show that: it sits close enough to the edge
		that the halving covers it either way, and a mutation removing
		the widening was caught by nothing. exp(60.5) is short by 23
		digits without it.
	*/

	Set( a, "60.5" ); result.Clear();
	lz_exp_e( &a, &result );
	CheckDigits( "exp(60.5)", &result, kRefExp605, 62 );

	// ...and through a caller, since the hyperbolics are where this
	//	was actually costing digits...
	Set( a, "50" ); result.Clear();
	lz_hypersine( &a, &result );
	CheckDigits( "sinh(50)", &result, kRefSinh50, 57 );

	/*
		And what the gamma work turned up in the exponential, which is
		gated here because nothing else reaches it.

		lz_exp_e summed the fractional part of its argument over a
		fixed 32 terms of gRecipFact, so exp of anything with a
		fraction stopped at the truncation whatever the precision --
		1.3e-47 at 0.5 and 1.5e-38 at 0.94, which is exactly where
		exp(0.5) and exp(13.94) stopped. Whole arguments have no
		fractional part to sum and scaled properly, which is what hid
		it: exp(1), exp(20) and exp(50) were all correct.
	*/

	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = 100;

		CLZNumber	wide, wideResult;

		wide.SetPrecision( 100 );
		wide.Assign( "0.5" );
		wideResult.SetPrecision( 100 );

		lz_exp_e( &wide, &wideResult );

		glWorkingPrecision = l_saved;

		CheckDigits( "exp(0.5) at precision 100", &wideResult, kRefExpHalfWide, 99 );
	}

	/*
		Gamma and its logarithm in other bases. In base two both were
		good to five to eight digits below about 2, where 36 binary
		places hold eleven: the shift that carries the argument up to 24
		before Stirling's series is summed compared it with a 24 made in
		base ten, digit for digit, and "100" -- four -- passed for more.
		Fixed in CLZNumber::CompareABS, which now brings a value in
		another base into its own before comparing.

		The arguments are exact in bases two, eight and sixteen, so the
		references need no rounding of the input. They come from
		Gamma(1/2) = sqrt(pi), Gamma(1/4) by the AGM formula, Gamma(3/4)
		by reflection and the rest by the recurrence; Gamma(1/4) and
		Gamma(3/4) agree with an independent Stirling sum to 1e-75. The
		thresholds are what each base achieves -- 36 places of base two
		are about eleven digits.
	*/

	{
		struct { const char* szValue; const char* szGamma; const char* szLogGamma; } exact[] =
		{
			{ "0.25", "3.625609908221908311930685155867672002995167682880065467",   "1.288022524698077457370610440219717295925377565112860550" },
			{ "0.5",  "1.772453850905516027298167483341145182797549456122387128",   "0.5723649429247000870717136756765293558236474064576557858" },
			{ "0.75", "1.225416702465177645129098303362890526851239248108070611",   "0.2032809514312953714814329718624296997596673149825786481" },
			{ "1.25", "0.9064024770554770779826712889669180007487919207200163669",  "-0.09827183642181316146385380269663584022562270360764995774" },
			{ "2.25", "1.133003096319346347478339111208647500935989900900020459",   "0.1248717148923965943024412876131986631489783819403572559" },
			{ "3.5",  "3.323350970447842551184064031264647217745405230229475865",   "1.200973602347074224816021881450712995770238915468157197" },
		};

		// Base sixteen's gamma was 38 to 40 digits until the exponential's
		//	table threshold was read in digits and base-ten constants were
		//	converted at the width the work needs, on 2026-09-24; 44 to 45
		//	since...
		struct { long lBase; int iGamma; int iLogGamma; } bases[] =
		{
			{  2, 11, 10 }, {  8, 33, 32 }, { 16, 44, 42 },
		};

		for ( auto& theBase : bases )
		{
			for ( auto& theCase : exact )
			{
				// Fresh numbers each time: Assign reads text in the
				//	number's current base, so a reused one reads "0.25"
				//	in base two...
				char		szLabel[ 72 ];
				CLZNumber	x, answer, y, logAnswer;

				x.Assign( theCase.szValue );
				x.SetBase( theBase.lBase, true );
				lz_gamma( &x, &answer );
				answer.SetBase( 10, true );
				::snprintf( szLabel, 72, "gamma(%s) in base %ld", theCase.szValue, theBase.lBase );
				CheckDigits( szLabel, &answer, theCase.szGamma, theBase.iGamma );

				y.Assign( theCase.szValue );
				y.SetBase( theBase.lBase, true );
				lz_log_gamma( &y, &logAnswer );
				logAnswer.SetBase( 10, true );
				::snprintf( szLabel, 72, "loggamma(%s) in base %ld", theCase.szValue, theBase.lBase );
				CheckDigits( szLabel, &logAnswer, theCase.szLogGamma, theBase.iLogGamma );
			}
		}
	}

	/*
		A large gamma in base two. lz_gamma widens its logarithm by the
		answer's magnitude, capped at thirty places -- meant as thirty
		decimal places, applied as thirty binary ones, nine digits. So
		base two stopped at twenty digits for every gamma past about
		1e10, where base ten went on to fifty-two. The cap is converted
		to the answer's base now: gamma(20.5) takes all 29 digits its
		size affords in base two and gamma(30.5) all 42.

		References by the recurrence from sqrt(pi), agreeing with an
		independent Stirling sum to 1e-114.
	*/

	{
		struct { const char* szValue; const char* szRef; int iDigits; } large[] =
		{
			{ "20.5", "540624298233507504.473687364780822141227316639578171478960195", 30 },
			{ "30.5", "48226969334909086010917483030261.2545074465525956899632826414", 42 },
		};

		for ( auto& theCase : large )
		{
			char		szLabel[ 72 ];
			CLZNumber	x, answer;

			x.Assign( theCase.szValue );
			x.SetBase( 2, true );
			lz_gamma( &x, &answer );
			answer.SetBase( 10, true );
			::snprintf( szLabel, 72, "gamma(%s) in base 2", theCase.szValue );
			CheckDigits( szLabel, &answer, theCase.szRef, theCase.iDigits );
		}
	}

	/*
		And the shift, whose threshold is in decimal digits while the
		working precision is in places of the argument's base: 44 hex
		places are 53 digits, and compared as places they asked for a
		series good to 44. At the default precision that showed in one
		case of thirty-six measured -- log-gamma of 2.875 in hex came
		back a unit high in its last place, 0.64 of a unit from the
		truth, where the converted threshold answers it correctly
		rounded. The expected string is Python's, rounded half-even.
	*/

	{
		CLZNumber	x, answer;

		x.Assign( "2.875" );
		x.SetBase( 16, true );
		lz_log_gamma( &x, &answer );
		CheckRendering( "loggamma(2.875) in hex, rounded", &answer, "0.94B838517FDE7E26B72A43A0C607330EF1AD" );
	}
}

// Build a list from text, in base ten unless told otherwise...
static void MakeList( CLZNumber* pValues, CLZNumber** ppList, const char** szValues, long count )
{
	for ( long i = 0; i < count; i++ )
	{
		pValues[ i ].Assign( szValues[ i ] );
		ppList[ i ] = &pValues[ i ];
	}
}

static void TestMeans()
{
	/*
		The means of a list, arithmetic, geometric and harmonic -- the
		only functions that take one. Each is a small calculation with
		one of the engine's recurring traps in it; see lz_amean.
	*/

	Banner( "The means of a list" );

	CLZNumber	values[ 8 ];
	CLZNumber*	list[ 8 ];
	CLZNumber	result;

	/*
		Exact answers, checked by rendering. CheckExact normalises, and
		reads 6.000...000 as 6 -- which is how a whole-number answer
		written as a rounded one got past it for the gamma function.
	*/

	struct { const char* szLabel; long lCount; const char* szValues[ 5 ];
			 char cMean; const char* szExpect; } exact[] =
	{
		{ "amean(4, 9)",      2, { "4", "9" },         'a', "6.5"  },
		{ "gmean(4, 9)",      2, { "4", "9" },         'g', "6"    },
		{ "hmean(2, 8)",      2, { "2", "8" },         'h', "3.2"  },
		{ "amean(2, 8, 32)",  3, { "2", "8", "32" },   'a', "14"   },
		{ "gmean(2, 8, 32)",  3, { "2", "8", "32" },   'g', "8"    },
		{ "hmean(1,2,3,4)",   4, { "1", "2", "3", "4" }, 'h', "1.92" },
		{ "gmean(10, 1000)",  2, { "10", "1000" },     'g', "100"  },
	};

	for ( auto& theCase : exact )
	{
		MakeList( values, list, theCase.szValues, theCase.lCount );
		result.Clear();

		if ( theCase.cMean == 'a' )			lz_amean( list, theCase.lCount, &result );
		else if ( theCase.cMean == 'g' )	lz_gmean( list, theCase.lCount, &result );
		else								lz_hmean( list, theCase.lCount, &result );

		CheckRendering( theCase.szLabel, &result, theCase.szExpect );
	}

	// And where the answer is irrational, against generated references...
	struct { const char* szLabel; long lCount; const char* szValues[ 5 ];
			 char cMean; const char* szRef; int iRequired; } measured[] =
	{
		{ "gmean(1,2,3,4)",   4, { "1", "2", "3", "4" },            'g', kRefGmean1234,  36 },
		{ "gmean(1,2,3)",     3, { "1", "2", "3" },                 'g', kRefGmean123,   36 },
		{ "hmean(4,9)",       2, { "4", "9" },                      'h', kRefHmean49,    36 },
		{ "hmean(1,2,3)",     3, { "1", "2", "3" },                 'h', kRefHmean123,   36 },
		{ "amean, scattered", 5, { "3.7", "1.1", "250", "0.004", "17" }, 'a', kRefAmeanMixed, 37 },
		{ "gmean, scattered", 5, { "3.7", "1.1", "250", "0.004", "17" }, 'g', kRefGmeanMixed, 36 },
		{ "hmean, scattered", 5, { "3.7", "1.1", "250", "0.004", "17" }, 'h', kRefHmeanMixed, 35 },
	};

	for ( auto& theCase : measured )
	{
		MakeList( values, list, theCase.szValues, theCase.lCount );
		result.Clear();

		if ( theCase.cMean == 'a' )			lz_amean( list, theCase.lCount, &result );
		else if ( theCase.cMean == 'g' )	lz_gmean( list, theCase.lCount, &result );
		else								lz_hmean( list, theCase.lCount, &result );

		CheckDigits( theCase.szLabel, &result, theCase.szRef, theCase.iRequired );
	}

	/*
		**The traps, with values that DIFFER.** A mean of identical
		values is answered from the first one without any arithmetic,
		so two copies of 1e-21 would never reach the code these exist
		to guard -- the obvious test for each trap tests nothing.

		The geometric mean of 1e-21 and 4e-21 is 2e-21. Multiplied
		directly the product is 4e-42, nought at 36 places, and the
		mean came back nought.

		The harmonic mean of 1e40 and 3e40 looks like the same trap and
		is not: division answers wider than it is asked, so no
		reciprocal underflows. It is kept as an accuracy check at a
		large size -- it would fail if division ever stopped doing that
		-- and not as a test of any widening, of which there is none.
	*/

	{
		const char*	szSmall[] = { "0.000000000000000000001", "0.000000000000000000004" };
		MakeList( values, list, szSmall, 2 );
		result.Clear();
		lz_gmean( list, 2, &result );
		CheckRendering( "gmean(1e-21, 4e-21)", &result, "0.000000000000000000002" );

		const char*	szLarge[] = { "10000000000000000000000000000000000000000",
								  "30000000000000000000000000000000000000000" };
		MakeList( values, list, szLarge, 2 );
		result.Clear();
		lz_hmean( list, 2, &result );
		CheckRendering( "hmean(1e40, 3e40)", &result, "15000000000000000000000000000000000000000" );

		/*
			And the width the scale-up at the end has to be paid for in.
			The root is worked with the exponents taken out and then
			multiplied back by b^q, which moves the radix, so a root
			worked to W places comes out with W - q. sqrt(2) * 1e40 at 36
			places has room for 76 digits; without the widening it got
			about 40. The second has an odd exponent sum, so the root is
			taken of 10 * M rather than M -- the remainder path.

			gmean(10, 1000) could not show this: it is 100 exactly, and
			an exact root rounds the same at any width.
		*/

		const char*	szBigEven[] = { "10000000000000000000000000000000000000000",
									"20000000000000000000000000000000000000000" };
		MakeList( values, list, szBigEven, 2 );
		result.Clear();
		lz_gmean( list, 2, &result );
		CheckDigits( "gmean(1e40, 2e40)", &result, kRefGmeanBigEven, 75 );

		const char*	szBigOdd[] = { "10000000000000000000000000000000000000000",
								   "200000000000000000000000000000000000000000" };
		MakeList( values, list, szBigOdd, 2 );
		result.Clear();
		lz_gmean( list, 2, &result );
		CheckDigits( "gmean(1e40, 2e41)", &result, kRefGmeanBigOdd, 76 );

		// The two ends of the range at once, which multiplies back to one...
		const char*	szBoth[] = { "0.0000000000000000000000000000000000000001",
								 "10000000000000000000000000000000000000000" };
		MakeList( values, list, szBoth, 2 );
		result.Clear();
		lz_gmean( list, 2, &result );
		CheckRendering( "gmean(1e-40, 1e40)", &result, "1" );
	}

	/*
		**Correctly rounded, not merely close.** These pin the guard, and
		they have to be checked by rendering against the correctly
		rounded answer: every other test here uses CheckDigits against a
		threshold, which passes a result one place off, and removing the
		guard altogether was caught by none of them.

		Measured over 400 random lists, the share correctly rounded at
		36 places: at no guard 305 geometric and 374 harmonic, at two
		397 and 400, at four all 400 of all three, and eight buys
		nothing more. These cases are drawn from the ones that fail at
		two or at none -- the first three fail at two -- and the
		expected strings are the correctly rounded values, generated.
	*/

	struct { const char* szLabel; long lCount; const char* szValues[ 7 ];
			 char cMean; const char* szExpect; } rounded[] =
	{
		{ "gmean rounds right (1)", 2, { "65205500000000", "7374280000000" }, 'g',
		  "21928146628021.256461656663365933185851371821149711" },
		{ "gmean rounds right (2)", 5, { "3198220000", "0.796392", "30404600",
		  "524381000000000", "18.6542" }, 'g',
		  "3765993.956893153622090170339396797917214813" },
		{ "gmean rounds right (3)", 7, { "1619500000000000", "674715000000000",
		  "0.00383972", "0.00265403", "4.90457", "0.000782953", "417603000" }, 'g',
		  "29143.643698509545918492774272713513808003" },
		{ "gmean rounds right (4)", 2, { "158493000000", "94880700" }, 'g',
		  "3877876581.984011517466048171466311703949445021" },
		{ "hmean rounds right (1)", 2, { "158493000000", "94880700" }, 'h',
		  "189647868.660874285836900019813430800201115242" },
		{ "hmean rounds right (2)", 2, { "5657520000000", "90479300000000000" }, 'h',
		  "11314332533622.656001441420580555299353363723610443" },
		{ "hmean rounds right (3)", 2, { "6.90505", "12650000" }, 'h',
		  "13.810092461709173697634004676072411975" },
	};

	for ( auto& theCase : rounded )
	{
		MakeList( values, list, theCase.szValues, theCase.lCount );
		result.Clear();

		if ( theCase.cMean == 'g' )	lz_gmean( list, theCase.lCount, &result );
		else						lz_hmean( list, theCase.lCount, &result );

		CheckRendering( theCase.szLabel, &result, theCase.szExpect );
	}

	// The mean of identical values is that value, exactly, in all three...
	{
		const char*	szSeven[] = { "7" };
		MakeList( values, list, szSeven, 1 );

		result.Clear(); lz_hmean( list, 1, &result );
		CheckRendering( "hmean(7)", &result, "7" );

		const char*	szFives[] = { "5", "5", "5", "5" };
		MakeList( values, list, szFives, 4 );

		result.Clear(); lz_gmean( list, 4, &result );
		CheckRendering( "gmean(5, 5, 5, 5)", &result, "5" );

		result.Clear(); lz_hmean( list, 4, &result );
		CheckRendering( "hmean(5, 5, 5, 5)", &result, "5" );
	}

	/*
		The domain. A negative is refused by the geometric and harmonic
		means -- even a list of identical negatives, which the identical
		values shortcut would otherwise have answered -- and a nought
		makes either of them nought.
	*/

	{
		const char*	szNeg[] = { "-2", "-2" };
		MakeList( values, list, szNeg, 2 );

		result.Clear();
		CheckTrue( "gmean(-2, -2) refused",
				   (lz_gmean( list, 2, &result ) == kFailed) && result.FlagIsOn( kfQNaN ), "" );

		result.Clear();
		CheckTrue( "hmean(-2, -2) refused",
				   (lz_hmean( list, 2, &result ) == kFailed) && result.FlagIsOn( kfQNaN ), "" );

		result.Clear();
		lz_amean( list, 2, &result );
		CheckRendering( "amean(-2, -2) is fine", &result, "-2" );

		const char*	szMixed[] = { "4", "-1", "9" };
		MakeList( values, list, szMixed, 3 );
		result.Clear();
		CheckTrue( "gmean refuses a negative among positives",
				   lz_gmean( list, 3, &result ) == kFailed, "" );

		const char*	szZero[] = { "4", "0", "9" };
		MakeList( values, list, szZero, 3 );

		result.Clear(); lz_gmean( list, 3, &result );
		CheckRendering( "gmean with a nought", &result, "0" );

		result.Clear();
		long	lHmeanZero = lz_hmean( list, 3, &result );
		CheckRendering( "hmean with a nought", &result, "0" );

		/*
			**Nought, and cleanly.** Without its own branch the harmonic
			mean still answers nought -- one over nought is infinite, and
			n over an infinite sum is nought -- but it gets there by a
			division by zero, and comes back with kDivideByZero and the
			flag set. The rendering alone could not tell the two apart,
			and a mutation removing the branch passed until this check
			existed.
		*/

		CheckTrue( "hmean with a nought is not an error",
				   (lHmeanZero == kNoError) && result.FlagIsOff( kfDiv0 ) &&
				   result.FlagIsOff( kfInfinity ) && result.FlagIsOff( kfQNaN ), "" );
	}

	/*
		An infinity or a not-a-number makes each of them a not-a-number:
		the engine holds an infinity as nought with a flag, often with
		no sign, and all three need one.
	*/

	{
		const char*	szTwo[] = { "4", "9" };

		for ( int iWhich = 0; iWhich < 2; iWhich++ )
		{
			MakeList( values, list, szTwo, 2 );

			if ( iWhich == 0 )	MakeInfinite( values[ 1 ] );
			else				MakeNotANumber( values[ 1 ] );

			long	lA, lG, lH;
			CLZNumber	a, g, h;

			lA = lz_amean( list, 2, &a );
			lG = lz_gmean( list, 2, &g );
			lH = lz_hmean( list, 2, &h );

			CheckTrue( iWhich ? "a NaN among the values" : "an infinity among the values",
					   (lA == kNoError) && (lG == kNoError) && (lH == kNoError) &&
					   a.FlagIsOn( kfQNaN ) && g.FlagIsOn( kfQNaN ) && h.FlagIsOn( kfQNaN ), "" );
		}
	}

	// A list that is not one...
	{
		const char*	szTwo[] = { "4", "9" };
		MakeList( values, list, szTwo, 2 );

		CheckTrue( "an empty list is refused",   lz_amean( list, 0, &result ) == kFailed, "" );
		CheckTrue( "a null list is refused",     lz_gmean( nullptr, 2, &result ) == kFailed, "" );
		CheckTrue( "a null target is refused",   lz_hmean( list, 2, nullptr ) == kFailed, "" );

		CLZNumber*	holey[] = { &values[ 0 ], nullptr };
		CheckTrue( "a null value is refused",    lz_amean( holey, 2, &result ) == kFailed, "" );
	}

	/*
		The first value's base wins, as it does for every binary
		operator. Sixteen in base sixteen with sixteen in base ten is a
		mean of sixteen, answered in base sixteen as 10.
	*/

	{
		values[ 0 ].SetBase( 16, true );
		values[ 0 ].Assign( "10" );
		values[ 1 ].SetBase( 10, true );
		values[ 1 ].Assign( "16" );
		list[ 0 ] = &values[ 0 ];
		list[ 1 ] = &values[ 1 ];

		CLZNumber	mixed;
		lz_amean( list, 2, &mixed );
		CheckTrue( "mixed bases: answered in the first value's base", mixed.m_lBase == 16, "" );
		CheckRendering( "mixed bases: 0x10 with 16", &mixed, "10" );

		/*
			The arithmetic mean cannot show whether the values were
			converted -- lz_add converts its right operand itself. The
			geometric mean can: it reads each value's exponent, and 100
			has exponent 2 in base ten but 1 in base sixteen. Left
			unconverted, 0x10 with 100 answered 0x40; the mean is
			sqrt(1600), which is 40, or 0x28.
		*/

		values[ 1 ].Assign( "100" );
		lz_gmean( list, 2, &mixed );
		CheckRendering( "mixed bases: gmean(0x10, 100)", &mixed, "28" );

		values[ 0 ].SetBase( 10, true );
	}

	// The target may be one of the values, since none of them is modified...
	{
		const char*	szTwo[] = { "4", "9" };
		MakeList( values, list, szTwo, 2 );

		lz_gmean( list, 2, &values[ 0 ] );
		CheckRendering( "target may be a value", &values[ 0 ], "6" );
	}
}

static void TestAgm()
{
	/*
		The arithmetic-geometric mean of two values. It had no test at
		all until it was brought up to the standard of the means beside
		it, and not one of 20 ordinary cases came back correctly rounded
		-- see lz_agm for the four ways it was wrong.

		**Every expected string here is the correctly rounded value at
		36 places, generated**, and checked by rendering: the answer is
		rounded now, so nothing short of the exact string is right.
	*/

	Banner( "The arithmetic-geometric mean" );

	CLZNumber	a, b, result;

	struct { const char* szLabel; const char* szLeft; const char* szRight;
			 const char* szExpect; } rounded[] =
	{
		{ "agm(1, 2)", "1", "2", "1.456791031046906869186432383265081975" },
		{ "agm(2, 1)", "2", "1", "1.456791031046906869186432383265081975" },

		/*
			**The first product underflowed.** 1e-40 * 1 is nought at 36
			places, so the geometric mean stayed nought and the answer
			was the arithmetic mean halved sixteen times: 3.05e-5 for
			0.0168. The next three are the same trap at other sizes.
		*/

		{ "agm(1e-40, 1)", "1e-40", "1", "0.016801811954039440829290115663535324" },
		{ "agm(3, 1e-10)", "3", "1e-10", "0.184721640102498661203112249796897106" },
		{ "agm(5e-30, 7e-31)", "5e-30", "7e-31", "0.000000000000000000000000000002334678" },
		{ "agm(1e-40, 1e40)", "1e-40", "1e40",
		  "84636568469133149231401786358729602483.272277540258568054543129701522205191" },

		/*
			**And the scale back up costs places.** The iteration works
			between one and the base and the answer is multiplied back
			by b^e, which moves the radix -- the cube root's band. These
			had 67 and 76 digits of the 75 and 77 on offer.
		*/

		{ "agm(1, 1e40)", "1", "1e40",
		  "168018119540394408292901156635353237918.097350825989424367069423521611820199" },
		{ "agm(1e40, 2e40)", "1e40", "2e40",
		  "14567910310469068691864323832650819749738.639432213055907941723832679264545803" },

		/*
			**Past what sixteen fixed cycles reach.** A ratio of 1e400
			takes about ten cycles to bring together, and the width the
			scale back needs about eight more to fill; stopped at
			sixteen, this had half its digits. It stops when it has
			converged now.
		*/

		{ "agm(1e-200, 1e200)", "1e-200", "1e200",
		  "170290731139946796577372574069075996958949823097659399133597592883587125333566549825135378129815991051988584796984389923854656662184723265984183733391884979628425290113424702592466484436955302862345.768896881882764800625747703632811756" },

		/*
			**These pin the guard.** Scored over 800 random pairs --
			half of them close in size, which is where the guard does
			its work -- against the correctly rounded answer: at no guard
			679 were right, at one 786, at two 798, and at three or more
			all 800. The first two fail at a guard of two, the next two
			at one, and the rest at none.
		*/

		{ "agm rounds right (1)", "991330838000790878e-14", "585878375119317650e-14",
		  "7752.966274950731638309562549281778688401" },
		{ "agm rounds right (2)", "559587194218305469e-9", "898770178753293726e-9",
		  "719146007.425016352655479059272535552947367009" },
		{ "agm rounds right (3)", "133556912321936839e-14", "731437578574800053e-15",
		  "1010.813325657115270798920567842983720275" },
		{ "agm rounds right (4)", "269548682153298641e-11", "995342114727349787e-10",
		  "31294945.243313664841890904717964844773827591" },
		{ "agm rounds right (5)", "4151782429e32", "5012770943698954e26",
		  "457213767825491989688464089686077354264721.123855325886873103086889977470896397" },
		{ "agm rounds right (6)", "6970810323e33", "3e-37",
		  "59466176274040852745506283109726328883299.015565528924989425858828250406031585" },
		{ "agm rounds right (7)", "73658073984223723e-5", "465879e-43",
		  "10089889746.491692903254529894877515307212843077" },
	};

	for ( auto& theCase : rounded )
	{
		Set( a, theCase.szLeft );
		Set( b, theCase.szRight );
		result.Clear();

		lz_agm( &a, &b, &result );

		CheckRendering( theCase.szLabel, &result, theCase.szExpect );
	}

	// At a raised precision, which the fixed count was never checked at...
	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = 100;

		CLZNumber	one, two, wideResult;

		one.SetPrecision( 100 );
		one.Assign( "1" );
		two.SetPrecision( 100 );
		two.Assign( "2" );
		wideResult.SetPrecision( 100 );

		lz_agm( &one, &two, &wideResult );

		glWorkingPrecision = l_saved;

		CheckRendering( "agm(1, 2) at precision 100", &wideResult,
						"1.4567910310469068691864323832650819749738639432213055907941723832679264545802509002574737128184484443" );
	}

	/*
		Exact answers. The AGM of identical values is that value, and a
		nought anywhere makes it nought -- which it did not: agm(0, 5)
		answered 5 / 2^15, the arithmetic mean halving on a geometric
		one stuck at nought until the sixteen cycles ran out.
	*/

	{
		Set( a, "3" ); Set( b, "3" ); result.Clear();
		lz_agm( &a, &b, &result );
		CheckRendering( "agm(3, 3)", &result, "3" );

		Set( a, "0" ); Set( b, "5" ); result.Clear();
		long	lZero = lz_agm( &a, &b, &result );
		CheckRendering( "agm(0, 5)", &result, "0" );
		CheckTrue( "agm(0, 5) is not an error",
				   (lZero == kNoError) && result.FlagIsOff( kfDiv0 ) &&
				   result.FlagIsOff( kfInfinity ) && result.FlagIsOff( kfQNaN ), "" );

		Set( a, "5" ); Set( b, "0" ); result.Clear();
		lz_agm( &a, &b, &result );
		CheckRendering( "agm(5, 0)", &result, "0" );
	}

	/*
		**A negative is refused**, as MPFR refuses it -- including two
		identical negatives, which the identical values shortcut would
		otherwise have answered. It used to answer nought with a
		not-a-number flag and kNoError: the square root failed, and
		nothing passed its code on.
	*/

	{
		Set( a, "-1" ); Set( b, "4" ); result.Clear();
		CheckTrue( "agm(-1, 4) refused",
				   (lz_agm( &a, &b, &result ) == kFailed) && result.FlagIsOn( kfQNaN ), "" );

		Set( a, "4" ); Set( b, "-1" ); result.Clear();
		CheckTrue( "agm(4, -1) refused",
				   (lz_agm( &a, &b, &result ) == kFailed) && result.FlagIsOn( kfQNaN ), "" );

		Set( a, "-4" ); Set( b, "-4" ); result.Clear();
		CheckTrue( "agm(-4, -4) refused",
				   (lz_agm( &a, &b, &result ) == kFailed) && result.FlagIsOn( kfQNaN ), "" );
	}

	// An infinity or a not-a-number makes it a not-a-number, as for the means...
	for ( int iWhich = 0; iWhich < 2; iWhich++ )
	{
		Set( a, "4" ); Set( b, "9" ); result.Clear();

		if ( iWhich == 0 )	MakeInfinite( b );
		else				MakeNotANumber( b );

		long	lFlagged = lz_agm( &a, &b, &result );

		CheckTrue( iWhich ? "agm with a NaN" : "agm with an infinity",
				   (lFlagged == kNoError) && result.FlagIsOn( kfQNaN ), "" );
	}

	// A missing operand...
	{
		Set( a, "4" ); Set( b, "9" ); result.Clear();

		CheckTrue( "agm with no right operand is refused", lz_agm( &a, nullptr, &result ) == kFailed, "" );
		CheckTrue( "agm with no left operand is refused",  lz_agm( nullptr, &b, &result ) == kFailed, "" );
	}

	/*
		With no target the answer goes to the left operand, as for every
		function with two; the right one is untouched. And a target may
		be the right operand, since neither is modified before the
		answer is written.
	*/

	{
		Set( a, "1" ); Set( b, "2" );
		lz_agm( &a, &b );
		CheckRendering( "agm with no target writes the left", &a, "1.456791031046906869186432383265081975" );
		CheckRendering( "...and leaves the right alone", &b, "2" );

		Set( a, "1" ); Set( b, "2" );
		lz_agm( &a, &b, &b );
		CheckRendering( "agm target may be the right operand", &b, "1.456791031046906869186432383265081975" );
	}

	/*
		The left operand's base wins. Sixteen in base sixteen with a
		hundred in base ten is agm(16, 100), answered in base sixteen to
		36 hexadecimal places -- the correctly rounded hexadecimal,
		generated. Left unconverted, the hundred's exponent would be read
		as two where it is one, and the scaling would be wrong.
	*/

	{
		a.SetBase( 16, true );
		a.Assign( "10" );
		b.SetBase( 10, true );
		b.Assign( "100" );

		CLZNumber	mixed;
		lz_agm( &a, &b, &mixed );

		CheckTrue( "agm mixed bases: answered in the left's base", mixed.m_lBase == 16, "" );
		CheckRendering( "agm mixed bases: (0x10, 100)", &mixed, "30.95114B7AE15C038D660DB9E026236A748B22" );

		a.SetBase( 10, true );
	}
}

/////////////////////////////////////////////////////////////////
//	Every answer in the left operand's base
/////////////////////////////////////////////////////////////////

/*
	Rich's rule, 2026-09-24: an answer comes back in the base of the
	left operand -- of the operand, for a function of one, and of the
	first value, for a mean of a list.

	Until then only the arithmetic, the division family and the
	comparisons kept it. The rest answered in whichever base their last
	working value had, and a fresh number is base ten: sine of a
	hexadecimal operand came back in base ten, lz_log_n in its right
	operand's base, lz_power in either, and a refusal's zero in the
	target's old base. 354 of the 1126 answers below were in the wrong
	base.

	**The target is set to base two first**, a base no operand is in,
	so an answer that never touched it cannot pass by accident; and
	every function is run in place as well, where the operand itself
	must keep its base. Operands in sixteen and seven, several values
	each, some of them refused. Every public function that answers a
	number is here -- the table is the list, so a new one belongs in it.
*/

/////////////////////////////////////////////////////////////////
//	Work done in the operand's own base. A fresh CLZNumber is base
//	ten, and a function whose temporary or constant sat on the left of
//	an operation did the work from there on in base ten and converted
//	the answer back. Until 2026-09-25 a hex reciprocal came back with
//	thirty places of 36, 5.6 million units out, and base two's with 120
//	unrounded.
/////////////////////////////////////////////////////////////////

// |got - ref| in units of the lPlaces'th place of lBase, got in lBase
//	and ref in base ten. Worked at 200 places, which holds a hex
//	answer's decimal expansion exactly...
static void UnitsOfThePlace( CLZNumber* pGot, CLZNumber* pRef, long lBase, long lPlaces, CLZNumber* pUnits )
{
	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = 200;

	CLZNumber	got, difference, scale, base, places;

	got.Assign( pGot );
	got.SetPrecision( 400 );
	got.SetBase( 10, true );

	lz_subtract( &got, pRef, &difference );
	lz_abs( &difference );

	base.Assign( lBase );
	places.Assign( lPlaces );
	lz_raise_n( &base, &places, &scale );

	lz_multiply( &difference, &scale, pUnits );

	glWorkingPrecision = l_saved;
}

static void UnitsOfThe36th( CLZNumber* pGot, CLZNumber* pRef, long lBase, CLZNumber* pUnits )
{
	UnitsOfThePlace( pGot, pRef, lBase, 36, pUnits );
}

static void TestWorkInOperandBase()
{
	Banner( "Work done in the operand's own base" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	/*
		Every one-operand function in bases sixteen, twelve, eight and
		two, on arguments exact in all four, against the same function
		in base ten at eighty places on the same value -- base ten being
		what the audit scores against Python. lz_exp_b and lz_log_b are
		relative to the operand's base, so theirs is lz_exp_n or
		lz_log_n with that base written out. The worst error over the
		arguments, in units of the 36th place, has to be within two:
		the arc-cotangent is pi/2 less an arctangent already rounded,
		and comes to 1.5 in base eight. Work done in base ten shows as
		hundreds of units in base twelve and millions in hex -- in
		bases two and eight it is more accurate than their own places,
		so it is hex and twelve that catch it.

		The negative arguments reach gamma's reflection, and 2^-12 the
		arc-hypersecant's series, below 0.01: the first census passed
		neither, and missed the four sites only they reach. 2^-12 is
		given to that function alone, since at it the cosecant and its
		kin answer 4096 and meet the small-divisor floor in every base,
		ten included, 1e4 units there -- the floor the audit sets aside,
		not this.

		A decimal value converted to hex keeps thirty places, and until
		TestAnswerPrecision's fix four of these answered at that width;
		the operands here were once widened to dodge it.
	*/

	#define ENTRY( f ) { #f, f }

	struct { const char* szName; Unary fn; } fns[] =
	{
		ENTRY( lz_arccosecant ),		ENTRY( lz_arccosine ),			ENTRY( lz_arccotangent ),
		ENTRY( lz_archypercosecant ),	ENTRY( lz_archypercosine ),		ENTRY( lz_archypercotangent ),
		ENTRY( lz_archypersecant ),		ENTRY( lz_archypersine ),		ENTRY( lz_archypertangent ),
		ENTRY( lz_arcsecant ),			ENTRY( lz_arcsine ),			ENTRY( lz_arctangent ),
		ENTRY( lz_cosecant ),			ENTRY( lz_cosine ),				ENTRY( lz_cotangent ),
		ENTRY( lz_exp_b ),				ENTRY( lz_exp_e ),				ENTRY( lz_gamma ),
		ENTRY( lz_hypercosecant ),		ENTRY( lz_hypercosine ),		ENTRY( lz_hypercotangent ),
		ENTRY( lz_hypersecant ),		ENTRY( lz_hypersine ),			ENTRY( lz_hypertangent ),
		ENTRY( lz_log_b ),				ENTRY( lz_log_e ),				ENTRY( lz_log_gamma ),
		ENTRY( lz_raise_e ),			ENTRY( lz_reciprocal ),			ENTRY( lz_root_2 ),
		ENTRY( lz_root_3 ),				ENTRY( lz_secant ),				ENTRY( lz_sine ),
		ENTRY( lz_tangent ),
	};

	#undef ENTRY

	/*
		**No row is allowed more than two units.** There were loose
		bounds here, each at its own argument: hex sine near pi/4, 193
		units, the series table's forty digits against hex's 44; and the
		floor inherited from a small divisor in the hyperbolic and
		circular cosecant and cotangent, up to 110 units at 1/16. All
		were fixed on 2026-09-25 and 26 -- the pairs widened by their
		answer's size, the table lengthened -- and the bounds went with
		them.
	*/

	const char*	args[] = { "0.25", "0.75", "1.5", "3.25", "0.0625", "12.5", "-0.75", "-3.25", "0.000244140625" };
	const long	bases[] = { 16, 12, 8, 2 };

	/*
		lz_raise_e takes a whole exponent and refuses a fraction, and
		every argument above is one. It was measured here at their
		truncations -- nought, one, three, twelve and minus three --
		until it began refusing on 2026-09-26, and its four rows then
		went on passing with nothing measured, both sides refusing
		alike. So it has whole arguments of its own, one for each, and
		a row that measures nothing fails.
	*/

	const char*	wholeArgs[] = { "1", "2", "5", "3", "0", "12", "-1", "-3", "0" };

	static_assert( sizeof( wholeArgs ) == sizeof( args ), "one whole argument for each" );
	const long	l_lastArg = static_cast<long>( sizeof( args ) / sizeof( args[ 0 ] ) ) - 1;

	for ( auto& theFn : fns )
	{
		bool	b_relative = (theFn.fn == lz_exp_b) || (theFn.fn == lz_log_b);

		/*
			The references, in base ten at eighty places on the same
			values, made once for all four bases -- they are most of
			this test's time. The two base-relative functions make
			theirs per base, below. Made while the width is eighty,
			since a CLZNumber takes its precision at construction...
		*/

		const long	kArgs = static_cast<long>( sizeof( args ) / sizeof( args[ 0 ] ) );
		CLZNumber*	references[ kArgs ];
		long		referenceCodes[ kArgs ];
		long		l_saved = glWorkingPrecision;

		glWorkingPrecision = 80;

		for ( long l_arg = 0; l_arg < kArgs; l_arg++ )
		{
			references[ l_arg ] = new CLZNumber;
			referenceCodes[ l_arg ] = kNoError;

			// 2^-12 serves the arc-hypersecant alone -- see above...
			if ( !b_relative && ((l_arg != l_lastArg) || (theFn.fn == lz_archypersecant)) )
			{
				CLZNumber	xTen;

				xTen.Assign( (theFn.fn == lz_raise_e) ? wholeArgs[ l_arg ] : args[ l_arg ] );
				referenceCodes[ l_arg ] = theFn.fn( &xTen, references[ l_arg ] );
			}
		}

		glWorkingPrecision = l_saved;

		for ( long lBase : bases )
		{
			CLZNumber	worst;
			long		l_codesDiffer = 0;
			long		l_over = 0;
			long		l_measured = 0;

			worst.Assign( 0L );

			for ( long l_arg = 0; l_arg <= l_lastArg; l_arg++ )
			{
				const char*	szArg = (theFn.fn == lz_raise_e) ? wholeArgs[ l_arg ] : args[ l_arg ];

				if ( (l_arg == l_lastArg) && (theFn.fn != lz_archypersecant) )
					continue;

				CLZNumber	x, answer;

				x.Assign( szArg );
				x.SetBase( lBase, true );

				long	l_code = theFn.fn( &x, &answer );

				// A base-relative function's reference is b to the x, or
				//	the log of x to base b, for this base...
				if ( b_relative )
				{
					glWorkingPrecision = 80;

					CLZNumber	xTen, base, reference;

					xTen.Assign( szArg );
					base.Assign( lBase );

					referenceCodes[ l_arg ] = (theFn.fn == lz_exp_b) ? lz_exp_n( &base, &xTen, &reference )
																	 : lz_log_n( &base, &xTen, &reference );
					references[ l_arg ]->Assign( &reference );

					glWorkingPrecision = l_saved;
				}

				CLZNumber&	reference = *references[ l_arg ];
				long		l_codeTen = referenceCodes[ l_arg ];

				if ( (l_code != kNoError) || (l_codeTen != kNoError) )
				{
					if ( (l_code == kNoError) != (l_codeTen == kNoError) )
						l_codesDiffer++;

					continue;
				}

				CLZNumber	units;
				CLZNumber	bound;

				UnitsOfThe36th( &answer, &reference, lBase, &units );
				l_measured++;

				bound.Assign( 2L );

				if ( lz_compare( &units, &bound ) == kGreaterThan )
					l_over++;

				if ( lz_compare( &units, &worst ) == kGreaterThan )
					worst.Assign( &units );
			}

			char	szLabel[ 72 ];
			char	szDetail[ 96 ];
			char*	szWorst = worst.GetDisplayRep( RS_scientific, 2 );

			::snprintf( szLabel, 72, "%s in base %ld", theFn.szName, lBase );
			::snprintf( szDetail, 96, "%s units at worst, %ld over its bound", szWorst, l_over );

			CheckTrue( szLabel, (l_codesDiffer == 0) && (l_over == 0) && (l_measured > 0),
					   (l_codesDiffer != 0) ? "refused where base ten answered, or answered where it refused"
											: ((l_measured == 0) ? "every argument refused, so nothing was measured" : szDetail) );

			::free( szWorst );
		}

		for ( CLZNumber* p_reference : references )
			delete p_reference;
	}

	/*
		**Bases thirteen to fifteen near pi/4.** 36 places of them are 40
		to 42 digits, just past the forty the series table served, and
		they were out by up to 73 units in base fifteen until it was
		lengthened on 2026-09-26. The arguments are not exact in these
		bases, so the reference is taken on the value the operand
		actually holds, written out in base ten.
	*/

	{
		struct { const char* szName; Unary fn; } circular[] =
		{
			{ "lz_sine", lz_sine }, { "lz_cosine", lz_cosine }, { "lz_tangent", lz_tangent },
		};

		const char*	nearQuarter[] = { "0.7", "0.75", "0.78125", "2.3125" };

		for ( long lBase = 13; lBase <= 15; lBase++ )
		{
			for ( auto& theFn : circular )
			{
				CLZNumber	worst;
				long		l_over = 0;

				worst.Assign( 0L );

				for ( const char* szArg : nearQuarter )
				{
					CLZNumber	x, answer;

					x.Assign( szArg );
					x.SetBase( lBase, true );

					// A refusal or a flagged answer is over the bound: it
					//	compares unordered, and would otherwise pass...
					if ( (theFn.fn( &x, &answer ) != kNoError) ||
						 answer.FlagIsOn( kfQNaN ) || answer.FlagIsOn( kfInfinity ) )
					{
						l_over++;
						continue;
					}

					long	l_saved = glWorkingPrecision;

					glWorkingPrecision = 90;

					CLZNumber	held, reference;

					held.Assign( &x );
					held.SetPrecision( 400 );
					held.SetBase( 10, true );
					held.SetPrecision( 90 );
					theFn.fn( &held, &reference );

					glWorkingPrecision = l_saved;

					CLZNumber	units, bound;

					UnitsOfThe36th( &answer, &reference, lBase, &units );
					bound.Assign( 2L );

					if ( lz_compare( &units, &bound ) == kGreaterThan )
						l_over++;

					if ( lz_compare( &units, &worst ) == kGreaterThan )
						worst.Assign( &units );
				}

				char	szLabel[ 72 ];
				char	szDetail[ 72 ];
				char*	szWorst = worst.GetDisplayRep( RS_scientific, 2 );

				::snprintf( szLabel, 72, "%s in base %ld near pi/4", theFn.szName, lBase );
				::snprintf( szDetail, 72, "%s units at worst, %ld over two", szWorst, l_over );
				CheckTrue( szLabel, l_over == 0, szDetail );

				::free( szWorst );
			}
		}
	}

	/*
		**tan, cot, sec and csc where the sine or cosine is small**, in
		bases sixteen, twelve, seven and two -- see lz_tangent. The
		widening reads each size in the value's own base, so base ten
		alone would not show it reading the wrong one. The arguments sit
		near pi/2 and pi and are not exact in these bases, so the
		reference is taken on the value held, as above. **Scored in
		units of the rule's last place**, not the 36th: a small answer's
		places past the 36th are what the rule adds, and the 36th cannot
		see them. The rule's places come from the reference's leading
		digit in this base. Until 2026-09-26 these were out by up to
		1e36 units in hex.
	*/

	{
		struct { const char* szName; Unary fn; } quotients[] =
		{
			{ "lz_tangent", lz_tangent }, { "lz_cotangent", lz_cotangent }, { "lz_secant", lz_secant },
			{ "lz_cosecant", lz_cosecant },
		};

		const char*	nearZeros[] = { "1.5707963268", "1.57079632679489661923", "4.7123889804",
									"3.1415926536", "3.14159265358979323846", "0.0000001" };

		/*
			Base two sees only forty binary places at the first take, so
			near pi/2 and pi to twenty digits its cosine and sine read
			nought: tan, sec and csc refused, and tan and cot answered
			nought, until lzi_widen_until_seen. Base two alone caught
			that.
		*/
		const long	quotientBases[] = { 16, 12, 7, 2 };

		for ( long lBase : quotientBases )
		{
			for ( auto& theFn : quotients )
			{
				CLZNumber	worst;
				long		l_over = 0;

				worst.Assign( 0L );

				for ( const char* szArg : nearZeros )
				{
					CLZNumber	x, answer;

					x.Assign( szArg );
					x.SetBase( lBase, true );

					// A refusal or a flagged answer is over the bound: it
					//	compares unordered, and would otherwise pass...
					if ( (theFn.fn( &x, &answer ) != kNoError) ||
						 answer.FlagIsOn( kfQNaN ) || answer.FlagIsOn( kfInfinity ) )
					{
						l_over++;
						continue;
					}

					long	l_saved = glWorkingPrecision;

					glWorkingPrecision = 120;

					CLZNumber	held, reference, inBase;

					held.Assign( &x );
					held.SetPrecision( 400 );
					held.SetBase( 10, true );
					held.SetPrecision( 120 );
					theFn.fn( &held, &reference );

					// The rule's places, from where the reference's leading
					//	digit falls in this base...
					inBase.Assign( &reference );
					inBase.SetPrecision( 400 );
					inBase.SetBase( lBase, true );

					long	l_top = inBase.m_lDigitsLength - 1;

					while ( (l_top >= 0) && (inBase.m_lpDigits[ l_top ] == 0) )
						l_top--;

					long	l_exponent = l_top - inBase.m_lRadix;
					long	l_places = 36 + ((l_exponent < -1) ? (-l_exponent - 1) : 0);

					glWorkingPrecision = l_saved;

					CLZNumber	units, bound;

					UnitsOfThePlace( &answer, &reference, lBase, l_places, &units );
					bound.Assign( 2L );

					if ( lz_compare( &units, &bound ) == kGreaterThan )
						l_over++;

					if ( lz_compare( &units, &worst ) == kGreaterThan )
						worst.Assign( &units );
				}

				char	szLabel[ 72 ];
				char	szDetail[ 72 ];
				char*	szWorst = worst.GetDisplayRep( RS_scientific, 2 );

				::snprintf( szLabel, 72, "%s in base %ld near its zeros", theFn.szName, lBase );
				::snprintf( szDetail, 72, "%s units at worst, %ld over two", szWorst, l_over );
				CheckTrue( szLabel, l_over == 0, szDetail );

				::free( szWorst );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////
//	An answer at the engine's precision, whatever its operand
//	carries. Until 2026-09-25 lz_log_e and lz_arctangent took their
//	width from the operand in one place each -- the logarithm's
//	reduction divided a copy of it in place, the arctangent's
//	reflection divided into one -- so a 20-place operand gave a
//	20-place answer that reported 36, and lz_log_b, lz_arccotangent
//	and lz_archypercosine, built on those two, inherited it.
/////////////////////////////////////////////////////////////////

static void TestAnswerPrecision()
{
	Banner( "An answer at the engine's precision, whatever the operand's" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	/*
		The same value at 20, 36 and 60 places must give the same
		answer, digit for digit, and needs no reference to say so. The
		arguments reach the logarithm's reduction (3.25, 5), the
		arctangent's reflection (past 1) and both signs; they stop short
		of the huge ones, which only make gamma and the exponentials
		slow. In bases ten and sixteen, the second because a decimal
		value converted to hex keeps thirty places, which is how this
		was found.
	*/

	#define ENTRY( f ) { #f, f }

	struct { const char* szName; Unary fn; } fns[] =
	{
		ENTRY( lz_arccosecant ),		ENTRY( lz_arccosine ),			ENTRY( lz_arccotangent ),
		ENTRY( lz_archypercosecant ),	ENTRY( lz_archypercosine ),		ENTRY( lz_archypercotangent ),
		ENTRY( lz_archypersecant ),		ENTRY( lz_archypersine ),		ENTRY( lz_archypertangent ),
		ENTRY( lz_arcsecant ),			ENTRY( lz_arcsine ),			ENTRY( lz_arctangent ),
		ENTRY( lz_cosecant ),			ENTRY( lz_cosine ),				ENTRY( lz_cotangent ),
		ENTRY( lz_cube ),				ENTRY( lz_exp_b ),				ENTRY( lz_exp_e ),
		ENTRY( lz_factorial ),			ENTRY( lz_gamma ),				ENTRY( lz_hypercosecant ),
		ENTRY( lz_hypercosine ),		ENTRY( lz_hypercotangent ),		ENTRY( lz_hypersecant ),
		ENTRY( lz_hypersine ),			ENTRY( lz_hypertangent ),		ENTRY( lz_log_b ),
		ENTRY( lz_log_e ),				ENTRY( lz_log_gamma ),			ENTRY( lz_raise_e ),
		ENTRY( lz_reciprocal ),			ENTRY( lz_root_2 ),				ENTRY( lz_root_3 ),
		ENTRY( lz_secant ),				ENTRY( lz_sine ),				ENTRY( lz_square ),
		ENTRY( lz_tangent ),
	};

	#undef ENTRY

	const char*	args[] = { "0.75", "1.5", "3.25", "5", "0.0625", "-0.75", "-3.25" };
	const long	bases[] = { 10, 16 };
	const long	places[] = { 20, 36, 60 };

	// lz_raise_e refuses a fraction, and was measured here at these
	//	arguments' truncations until it did -- see TestWorkInOperandBase...
	const char*	wholeArgs[] = { "1", "2", "3", "5", "0", "-1", "-3" };

	static_assert( sizeof( wholeArgs ) == sizeof( args ), "one whole argument for each" );

	const long	kArgs = static_cast<long>( sizeof( args ) / sizeof( args[ 0 ] ) );

	for ( auto& theFn : fns )
	{
		for ( long lBase : bases )
		{
			long		l_differ = 0;
			const char*	szWhere = "";

			for ( long l_arg = 0; l_arg < kArgs; l_arg++ )
			{
				const char*	szArg = (theFn.fn == lz_raise_e) ? wholeArgs[ l_arg ] : args[ l_arg ];
				char*		szAnswers[ 3 ];
				long	lCodes[ 3 ];

				for ( int i = 0; i < 3; i++ )
				{
					CLZNumber	x, answer;

					x.Assign( szArg );
					x.SetBase( lBase, true );
					x.SetPrecision( places[ i ] );

					lCodes[ i ] = theFn.fn( &x, &answer );
					szAnswers[ i ] = answer.GetStringRep();
				}

				for ( int i = 0; i < 3; i += 2 )
				{
					if ( (lCodes[ i ] != lCodes[ 1 ]) || (::strcmp( szAnswers[ i ], szAnswers[ 1 ] ) != 0) )
					{
						l_differ++;
						szWhere = szArg;
					}
				}

				for ( char* szAnswer : szAnswers )
					::free( szAnswer );
			}

			char	szLabel[ 72 ];
			char	szDetail[ 72 ];

			::snprintf( szLabel, 72, "%s in base %ld", theFn.szName, lBase );

			if ( l_differ == 0 )
				::snprintf( szDetail, 72, "the same at 20, 36 and 60 places" );
			else
				::snprintf( szDetail, 72, "%ld answers moved with the operand, one at %s", l_differ, szWhere );

			CheckTrue( szLabel, l_differ == 0, szDetail );
		}
	}

	/*
		lz_archypersine had no guard, and arcsinh(1.5) lies 0.036 of a
		unit past a rounding boundary: ...523|536... With the logarithm
		given its own width and no guard here, it came back ...523. The
		reference is ln( x + sqrt( x^2 + 1 ) ) in Python's Decimal at 80
		digits.
	*/

	{
		CLZNumber	x, answer, xNegative, answerNegative;

		x.Assign( "1.5" );
		lz_archypersine( &x, &answer );
		CheckRendering( "arcsinh 1.5, a hair past half", &answer, "1.194763217287109304111930828519090524" );

		xNegative.Assign( "-1.5" );
		lz_archypersine( &xNegative, &answerNegative );
		CheckRendering( "arcsinh -1.5", &answerNegative, "-1.194763217287109304111930828519090524" );
	}

	/*
		lz_archypertangent had no guard either: 12 of its 33 audit cases
		missed correct rounding, all by one unit. It also answered past
		its 36th place -- it halves a logarithm already rounded to 36, in
		place, and half of that can need a 37th -- which the sweep's
		digit count read as a digit more than it had. The references are
		ln( (1 + x) / (1 - x) ) / 2 in Python's Decimal at 90 digits,
		rounded half-even to 36 places; every one of these was wrong
		before.
	*/

	{
		struct { const char* szLabel; const char* szValue; const char* szExpect; } tanhs[] =
		{
			{ "arctanh 0.05",     "0.05",     "0.050041729278491268245785274238925949" },
			{ "arctanh 0.3",      "0.3",      "0.309519604203111715474067349061069438" },
			{ "arctanh 0.7",      "0.7",      "0.867300527694053194427144690475300415" },
			{ "arctanh 0.999999", "0.999999", "7.254328619262047206741757083636344449" },
			{ "arctanh -0.05",    "-0.05",    "-0.050041729278491268245785274238925949" },
			{ "arctanh -0.7",     "-0.7",     "-0.867300527694053194427144690475300415" },
		};

		for ( auto& theCase : tanhs )
		{
			CLZNumber	x, answer;

			x.Assign( theCase.szValue );
			lz_archypertangent( &x, &answer );
			CheckRendering( theCase.szLabel, &answer, theCase.szExpect );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	No answer wider than it should be. Division answers wider than it
//	is asked, and until 2026-09-25 eight functions handed a quotient
//	straight back: a 37th place on an ordinary answer, unrounded, and
//	on a small one more places than 36 significant digits. The rule,
//	Rich's: 36 places for an answer of 0.1 or more, and 36 significant
//	digits below it -- see lzi_round_quotient_answer.
/////////////////////////////////////////////////////////////////

static void TestAnswerWidth()
{
	Banner( "No answer wider than 36 places, or 36 digits below 0.1" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	/*
		Every one-operand function, in bases ten and sixteen: the places
		an answer is written with may not pass the rule, the noughts
		after the point counted in the answer's own base. It needs no
		reference. The eight quotient functions broke it 133 times on
		these arguments before; no other function ever did.
	*/

	#define ENTRY( f ) { #f, f }

	struct { const char* szName; Unary fn; } fns[] =
	{
		ENTRY( lz_arccosecant ),		ENTRY( lz_arccosine ),			ENTRY( lz_arccotangent ),
		ENTRY( lz_archypercosecant ),	ENTRY( lz_archypercosine ),		ENTRY( lz_archypercotangent ),
		ENTRY( lz_archypersecant ),		ENTRY( lz_archypersine ),		ENTRY( lz_archypertangent ),
		ENTRY( lz_arcsecant ),			ENTRY( lz_arcsine ),			ENTRY( lz_arctangent ),
		ENTRY( lz_cosecant ),			ENTRY( lz_cosine ),				ENTRY( lz_cotangent ),
		ENTRY( lz_cube ),				ENTRY( lz_exp_b ),				ENTRY( lz_exp_e ),
		ENTRY( lz_factorial ),			ENTRY( lz_gamma ),				ENTRY( lz_hypercosecant ),
		ENTRY( lz_hypercosine ),		ENTRY( lz_hypercotangent ),		ENTRY( lz_hypersecant ),
		ENTRY( lz_hypersine ),			ENTRY( lz_hypertangent ),		ENTRY( lz_log_b ),
		ENTRY( lz_log_e ),				ENTRY( lz_log_gamma ),			ENTRY( lz_raise_e ),
		ENTRY( lz_reciprocal ),			ENTRY( lz_root_2 ),				ENTRY( lz_root_3 ),
		ENTRY( lz_secant ),				ENTRY( lz_sine ),				ENTRY( lz_square ),
		ENTRY( lz_tangent ),
	};

	#undef ENTRY

	const char*	args[] = { "0.3", "0.75", "1.5", "3.25", "7", "12.5", "50", "0.0625",
						   "0.00000001", "1E-20", "-0.75", "-3.25" };
	const long	bases[] = { 10, 16 };

	// lz_raise_e refuses a fraction, and was measured here at these
	//	arguments' truncations until it did -- see TestWorkInOperandBase...
	const char*	wholeArgs[] = { "1", "2", "3", "5", "7", "12", "50", "0",
								"-1", "-2", "-5", "-3" };

	static_assert( sizeof( wholeArgs ) == sizeof( args ), "one whole argument for each" );

	const long	kArgs = static_cast<long>( sizeof( args ) / sizeof( args[ 0 ] ) );

	for ( auto& theFn : fns )
	{
		for ( long lBase : bases )
		{
			long		l_over = 0;
			const char*	szWhere = "";

			for ( long l_arg = 0; l_arg < kArgs; l_arg++ )
			{
				const char*	szArg = (theFn.fn == lz_raise_e) ? wholeArgs[ l_arg ] : args[ l_arg ];
				CLZNumber	x, answer;

				x.Assign( szArg );
				x.SetBase( lBase, true );

				if ( theFn.fn( &x, &answer ) != kNoError )
					continue;

				answer.CheckForZero();

				if ( answer.IsZero() )
					continue;

				char*		szAnswer = answer.GetStringRep();
				const char*	szPoint = ::strchr( szAnswer, '.' );
				long		l_places = (szPoint != nullptr) ? static_cast<long>( ::strlen( szPoint + 1 ) ) : 0;

				// Where the leading digit sits, as lzi_leading_exponent
				//	has it: negative below one...
				long	l_top = answer.m_lDigitsLength - 1;

				while ( (l_top >= 0) && (answer.m_lpDigits[ l_top ] == 0) )
					l_top--;

				long	l_exponent = l_top - answer.m_lRadix;
				long	l_allowed = glWorkingPrecision + ((l_exponent < -1) ? (-l_exponent - 1) : 0);

				if ( l_places > l_allowed )
				{
					l_over++;
					szWhere = szArg;
				}

				::free( szAnswer );
			}

			char	szLabel[ 72 ];
			char	szDetail[ 72 ];

			::snprintf( szLabel, 72, "%s in base %ld", theFn.szName, lBase );

			if ( l_over == 0 )
				::snprintf( szDetail, 72, "no answer past the rule" );
			else
				::snprintf( szDetail, 72, "%ld answers past the rule, one at %s", l_over, szWhere );

			CheckTrue( szLabel, l_over == 0, szDetail );
		}
	}

	/*
		And the values, correctly rounded under the rule. References are
		Python's Decimal at 250 digits, rounded half-even at 36 places,
		or at 36 significant digits below 0.1. The last is csch where
		sinh by itself has run out of places and is nought; csch works
		wider by its answer's size, and answers.
	*/

	{
		struct { const char* szLabel; Unary fn; const char* szValue; const char* szExpect; } pins[] =
		{
			{ "1/0.3, 36 places",           lz_reciprocal,    "0.3",   "3.333333333333333333333333333333333333" },
			{ "cbrt 0.3, 36 places",        lz_root_3,        "0.3",   "0.669432950082169521882659324639930793" },
			{ "cbrt 1e-40, 36 digits",      lz_root_3,        "1E-40", "0.0000000000000464158883361277889241007635091944658" },
			{ "sech 50, 36 digits",         lz_hypersecant,   "50",    "0.000000000000000000000385749969592783556603468563305402515" },
			{ "csch 1e-41, widened",        lz_hypercosecant, "1E-41", "100000000000000000000000000000000000000000.000000000000000000000000000000000000" },
		};

		for ( auto& thePin : pins )
		{
			CLZNumber	x, answer;

			x.Assign( thePin.szValue );
			thePin.fn( &x, &answer );
			CheckRendering( thePin.szLabel, &answer, thePin.szExpect );
		}
	}

	/*
		csch and coth work wider by twice their answer's size below one,
		so the sine they divide by never runs out of places, however
		small the argument -- see lzi_small_divisor_widening. Their only
		pole is nought itself, pinned with the other reciprocal poles.
		These refused as divisions by zero until 2026-09-25; the answers
		are Python's, correctly rounded at 36 places.
	*/

	{
		struct { const char* szLabel; Unary fn; const char* szValue; const char* szExpect; } tiny[] =
		{
			{ "csch 1e-45",  lz_hypercosecant,  "1E-45",  "1000000000000000000000000000000000000000000000.000000000000000000000000000000000000" },
			{ "csch -1e-45", lz_hypercosecant,  "-1E-45", "-1000000000000000000000000000000000000000000000.000000000000000000000000000000000000" },
			{ "coth 1e-45",  lz_hypercotangent, "1E-45",  "1000000000000000000000000000000000000000000000.000000000000000000000000000000000000" },
			{ "csch 1e-12",  lz_hypercosecant,  "1E-12",  "999999999999.999999999999833333333333333333333333" },
			{ "coth 1e-12",  lz_hypercotangent, "1E-12",  "1000000000000.000000000000333333333333333333333333" },
		};

		for ( auto& theCase : tiny )
		{
			CLZNumber	x, answer;

			x.Assign( theCase.szValue );
			theCase.fn( &x, &answer );
			CheckRendering( theCase.szLabel, &answer, theCase.szExpect );
		}
	}

	/*
		The circular pair, the same way, and near pi as well as near
		nought: the circular sine is small near every multiple of pi, so
		the cosecant and cotangent read how small off the sine itself
		and take it again wider -- see lz_cosecant. The rows near pi
		fail with that second sine taken out; the rows near nought
		failed without the argument's widening until 2026-09-26, when
		lzi_widen_until_seen came to cover them and the widening became
		speed alone. cot(3.1415926535897932), whose sine is 4e-17, was
		wrong from its eighth decimal place. Every
		one of these was wrong before 2026-09-25. Python's Decimal at
		300 digits, rounded half-even at 36 places.
	*/

	{
		struct { const char* szLabel; Unary fn; const char* szValue; const char* szExpect; } circular[] =
		{
			{ "csc 1e-12",       lz_cosecant,  "1E-12",     "1000000000000.000000000000166666666666666666666667" },
			{ "cot 1e-12",       lz_cotangent, "1E-12",     "999999999999.999999999999666666666666666666666667" },
			{ "csc 1e-45",       lz_cosecant,  "1E-45",     "1000000000000000000000000000000000000000000000.000000000000000000000000000000000000" },
			{ "csc 3.14159",     lz_cosecant,  "3.14159",   "376847.997587737714116236245158352008495956" },
			{ "cot 3.14159",     lz_cotangent, "3.14159",   "-376847.997586410919219616235280182705497063" },
			{ "csc -3.14159",    lz_cosecant,  "-3.14159",  "-376847.997587737714116236245158352008495956" },
			{ "csc 6.2831853",   lz_cosecant,  "6.2831853", "-139283787.891396457477386942524774905330788427" },
			{ "cot pi to 17",    lz_cotangent, "3.1415926535897932",
			  "-25999253094360135.330097370830613424622244926673382004" },
		};

		for ( auto& theCase : circular )
		{
			CLZNumber	x, answer;

			x.Assign( theCase.szValue );
			theCase.fn( &x, &answer );
			CheckRendering( theCase.szLabel, &answer, theCase.szExpect );
		}
	}

	/*
		And where the cosine is small, or the sine is small and on top.
		tan and sec divide by the cosine near the odd multiples of pi/2,
		and were not widened for it until 2026-09-26: tan and sec at
		pi/2 to 21 digits were 9e36 units out of their 36th place. A
		small numerator -- the sine for tan near the multiples of pi,
		the cosine for cot near pi/2 -- gives a small answer, which the
		rule gives 36 significant digits, and they had 36 places, as few
		as 16 digits right. See lz_tangent. Every one of these was wrong
		before. gensweep's Decimal at 260 digits, rounded half-even
		under the rule; tan(1e-40) is compared by value, the engine
		writing it without its trailing noughts.
	*/

	{
		struct { const char* szLabel; Unary fn; const char* szValue; const char* szExpect; } cosine[] =
		{
			{ "tan 1e-8",            lz_tangent,    "1E-8",                    "0.0000000100000000000000003333333333333333467" },
			{ "tan -1e-8",           lz_tangent,    "-1E-8",                   "-0.0000000100000000000000003333333333333333467" },
			{ "tan 3.1415926536",    lz_tangent,    "3.1415926536",            "0.0000000000102067615373566167204974702427673778" },
			{ "tan pi to 21",        lz_tangent,    "3.14159265358979323846",  "-0.00000000000000000000264338327950288419716939937510582097" },
			{ "tan 1.5707963268",    lz_tangent,    "1.5707963268",            "-195948537905.977872854633345673562619927674709558" },
			{ "tan pi/2 to 21",      lz_tangent,    "1.57079632679489661923",  "756606132568153667453.844815332809344259564370261015900747" },
			{ "tan 4.7123889804",    lz_tangent,    "4.7123889804",            "-65316179301.992624284873245552726492590682459945" },
			{ "cot 1.5",             lz_cotangent,  "1.5",                     "0.0709148443026524487889808929348032891" },
			{ "cot 1.5707963268",    lz_cotangent,  "1.5707963268",            "-0.00000000000510338076867830836024860220640739746" },
			{ "cot pi/2 to 21",      lz_cotangent,  "1.57079632679489661923",  "0.00000000000000000000132169163975144209858469968755291049" },
			{ "sec 1.5707963268",    lz_secant,     "1.5707963268",            "-195948537905.977872854635897363946959081854833843" },
			{ "sec pi/2 to 21",      lz_secant,     "1.57079632679489661923",  "756606132568153667453.844815332809344259565031106835776468" },
			{ "sec 4.7123889804",    lz_secant,     "4.7123889804",            "65316179301.992624284880900623879510053222832931" },
		};

		for ( auto& theCase : cosine )
		{
			CLZNumber	x, answer;

			x.Assign( theCase.szValue );
			theCase.fn( &x, &answer );
			CheckRendering( theCase.szLabel, &answer, theCase.szExpect );
		}

		CLZNumber	x, answer;

		x.Assign( "1E-40" );
		lz_tangent( &x, &answer );
		CheckExact( "tan 1e-40, which was nought", &answer, "0.0000000000000000000000000000000000000001" );

		/*
			pi/2 to 49 places is 4.7e-50 short of pi/2, and its cosine
			read nought at forty places, so tan refused it as a pole
			until 2026-09-26 -- and TestQuietFailures pinned that. Its
			cosine is seen now. It carried what 2*pi's 120 places left
			for it, 71 digits of the 86 a correct 36th place wants, until
			pi was computed past them the same day -- see TestWidePi --
			and it is correctly rounded now.
		*/

		x.Assign( "1.5707963267948966192313216916397514420985846996876" );
		lz_tangent( &x, &answer );
		CheckDigits( "tan of pi/2 to 49 places", &answer,
					 "-21236151030692384854558538473739128298113204417313.765914935043385311411384064844058191720727686816138873998447", 85 );
	}
}

/////////////////////////////////////////////////////////////////
//	Every refusal answers. A function that refuses writes a nought
//	flagged not-a-number or infinite into its target -- in place, into
//	its operand -- so a caller that ignores the code still cannot take
//	the refusal for an answer. lz_arcsine alone returned kFailed and
//	wrote nothing until 2026-09-26: the target kept whatever it held.
/////////////////////////////////////////////////////////////////

static void TestRefusalsAnswer()
{
	Banner( "Every refusal writes a flagged answer" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );
	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );
	typedef long (*OfList)( CLZNumber**, long, CLZNumber* );

	/*
		Arguments on and past every domain's edges, into a target set to
		12345 first. Any refusal must leave that target flagged, and not
		holding 12345. It needs no reference. Found by running exactly
		this over every function: the arcsine was the only one of the
		fifty that left its target as it was.
	*/

	#define ENTRY( f ) { #f, f }

	struct { const char* szName; Unary fn; } unary[] =
	{
		ENTRY( lz_arccosecant ),		ENTRY( lz_arccosine ),			ENTRY( lz_arccotangent ),
		ENTRY( lz_archypercosecant ),	ENTRY( lz_archypercosine ),		ENTRY( lz_archypercotangent ),
		ENTRY( lz_archypersecant ),		ENTRY( lz_archypersine ),		ENTRY( lz_archypertangent ),
		ENTRY( lz_arcsecant ),			ENTRY( lz_arcsine ),			ENTRY( lz_arctangent ),
		ENTRY( lz_cosecant ),			ENTRY( lz_cosine ),				ENTRY( lz_cotangent ),
		ENTRY( lz_exp_b ),				ENTRY( lz_exp_e ),				ENTRY( lz_factorial ),
		ENTRY( lz_gamma ),				ENTRY( lz_hypercosecant ),		ENTRY( lz_hypercosine ),
		ENTRY( lz_hypercotangent ),		ENTRY( lz_hypersecant ),		ENTRY( lz_hypersine ),
		ENTRY( lz_hypertangent ),		ENTRY( lz_log_b ),				ENTRY( lz_log_e ),
		ENTRY( lz_log_gamma ),			ENTRY( lz_raise_e ),			ENTRY( lz_reciprocal ),
		ENTRY( lz_root_2 ),				ENTRY( lz_root_3 ),				ENTRY( lz_secant ),
		ENTRY( lz_sine ),				ENTRY( lz_tangent ),
	};

	struct { const char* szName; Binary fn; } binary[] =
	{
		ENTRY( lz_power ),				ENTRY( lz_root_n ),				ENTRY( lz_log_n ),
		ENTRY( lz_exp_n ),				ENTRY( lz_raise_n ),			ENTRY( lz_agm ),
		ENTRY( lz_combinations ),		ENTRY( lz_divide ),				ENTRY( lz_divide_integer ),
		ENTRY( lz_modulo ),				ENTRY( lz_remainder ),			ENTRY( lz_remainder_near ),
		ENTRY( lz_permutations ),
	};

	struct { const char* szName; OfList fn; } ofList[] =
	{
		ENTRY( lz_amean ),				ENTRY( lz_gmean ),				ENTRY( lz_hmean ),
	};

	#undef ENTRY

	const char*	args[] = { "2", "-2", "1.5", "-1.5", "0", "-1", "1", "-0.5", "0.5", "-3" };
	const char*	pairs[][ 2 ] = { { "-8", "2.5" }, { "0", "-1" }, { "2", "0" }, { "-2", "-2" },
								 { "3", "-1" }, { "-1", "3" }, { "0", "0" }, { "1", "0" } };
	const char*	lists[][ 2 ] = { { "-1", "2" }, { "0", "2" }, { "1", "-2" } };

	auto	answered = []( CLZNumber* pTarget ) -> bool
	{
		char*	szTarget = pTarget->GetStringRep();
		bool	b_untouched = (::strcmp( szTarget, "12345" ) == 0);

		::free( szTarget );

		return !b_untouched && (pTarget->FlagIsOn( kfQNaN ) || pTarget->FlagIsOn( kfInfinity ));
	};

	for ( auto& theFn : unary )
	{
		long		l_refused = 0;
		long		l_silent = 0;
		const char*	szWhere = "";

		for ( const char* szArg : args )
		{
			CLZNumber	x, target;

			x.Assign( szArg );
			target.Assign( 12345L );

			if ( theFn.fn( &x, &target ) == kNoError )
				continue;

			l_refused++;

			if ( !answered( &target ) )
			{
				l_silent++;
				szWhere = szArg;
			}
		}

		// Two longs, the words and an argument: 72 bytes held what these
		//	counts ever write, and GCC 15 counts what they could...
		char	szDetail[ 128 ];

		if ( l_silent == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "%ld refusals, each a flagged answer", l_refused );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld refusals wrote nothing, one at %s", l_silent, l_refused, szWhere );

		CheckTrue( theFn.szName, l_silent == 0, szDetail );
	}

	for ( auto& theFn : binary )
	{
		long	l_refused = 0;
		long	l_silent = 0;

		for ( auto& thePair : pairs )
		{
			CLZNumber	x, y, target;

			x.Assign( thePair[ 0 ] );
			y.Assign( thePair[ 1 ] );
			target.Assign( 12345L );

			if ( theFn.fn( &x, &y, &target ) == kNoError )
				continue;

			l_refused++;

			if ( !answered( &target ) )
				l_silent++;
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, 72, "%ld of %ld refusals wrote nothing", l_silent, l_refused );
		CheckTrue( theFn.szName, l_silent == 0, szDetail );
	}

	for ( auto& theFn : ofList )
	{
		long	l_refused = 0;
		long	l_silent = 0;

		for ( auto& theList : lists )
		{
			CLZNumber	x, y, target;
			CLZNumber*	pValues[] = { &x, &y };

			x.Assign( theList[ 0 ] );
			y.Assign( theList[ 1 ] );
			target.Assign( 12345L );

			if ( theFn.fn( pValues, 2, &target ) == kNoError )
				continue;

			l_refused++;

			if ( !answered( &target ) )
				l_silent++;
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, 72, "%ld of %ld refusals wrote nothing", l_silent, l_refused );
		CheckTrue( theFn.szName, l_silent == 0, szDetail );
	}

	// And in place, where the operand is the target...
	{
		CLZNumber	x;

		x.Assign( "2" );

		long	l_code = lz_arcsine( &x, nullptr );

		CheckTrue( "arcsin(2) in place", (l_code == kFailed) && x.FlagIsOn( kfQNaN ) && x.IsZero(),
				   "the operand becomes a flagged nought" );
	}
}

/////////////////////////////////////////////////////////////////
//	Rounding modes. Until 2026-09-26 no operation honoured its target's
//	mode: Assign copies the mode and the rounding digit along with the
//	value, so every target came back set to its temporary's half-even.
//	And a caller's settings leaked the other way, into the functions --
//	arcsin(0.1) answered 0 when its operand was set to round at the
//	units digit. See CLZTargetRounding and CLZInternalWork.
/////////////////////////////////////////////////////////////////

static void TestRoundingModes()
{
	Banner( "Rounding modes: the plain operations honour them, functions never see them" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );
	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );

	const RoundingMode	modes[] = { RM_round_down, RM_round_half_up, RM_round_half_even, RM_round_ceiling,
									RM_round_floor, RM_round_half_down, RM_round_up };
	const char*			modeNames[] = { "down", "half up", "half even", "ceiling", "floor", "half down", "up" };

	/*
		Division and multiplication in each of seven modes. Both are
		exact before they round -- the product outright, the quotient
		through its sticky digit -- so each answer is the exact value
		rounded that way. The expectations are Python's Decimal on the
		exact fraction, at the width the engine answers to. 11/16 and
		13/16 at one place leave an exact half, the sticky digit nought,
		for the half modes to break.
	*/

	{
		struct { const char* szLeft; const char* szRight; long lPlaces; const char* szExpect[ 7 ]; } quotients[] =
		{
		{ "2", "3", 36,
		    { "0.66666666666666666666666666666666666666",
		      "0.66666666666666666666666666666666666667",
		      "0.66666666666666666666666666666666666667",
		      "0.66666666666666666666666666666666666667",
		      "0.66666666666666666666666666666666666666",
		      "0.66666666666666666666666666666666666667",
		      "0.66666666666666666666666666666666666667" } },
		{ "-2", "3", 36,
		    { "-0.66666666666666666666666666666666666666",
		      "-0.66666666666666666666666666666666666667",
		      "-0.66666666666666666666666666666666666667",
		      "-0.66666666666666666666666666666666666666",
		      "-0.66666666666666666666666666666666666667",
		      "-0.66666666666666666666666666666666666667",
		      "-0.66666666666666666666666666666666666667" } },
		{ "11", "16", 1, { "0.687", "0.688", "0.688", "0.688", "0.687", "0.687", "0.688" } },
		{ "13", "16", 1, { "0.812", "0.813", "0.812", "0.813", "0.812", "0.812", "0.813" } },
		{ "-13", "16", 1, { "-0.812", "-0.813", "-0.812", "-0.812", "-0.813", "-0.812", "-0.813" } },
		{ "1", "7", 5, { "0.1428571", "0.1428571", "0.1428571", "0.1428572", "0.1428571", "0.1428571", "0.1428572" } },
		{ "-1", "7", 5, { "-0.1428571", "-0.1428571", "-0.1428571", "-0.1428571", "-0.1428572", "-0.1428571", "-0.1428572" } },
		};

		for ( auto& theCase : quotients )
		{
			long		l_wrong = 0;
			long		l_lost = 0;
			const char*	szWhere = "";

			for ( int m = 0; m < 7; m++ )
			{
				CLZNumber	a, b, target;

				Set( a, theCase.szLeft );
				Set( b, theCase.szRight );
				target.SetPrecision( theCase.lPlaces );
				target.SetRoundingMode( modes[ m ] );

				lz_divide( &a, &b, &target );

				l_lost += (target.m_rmMode != modes[ m ]);

				char*	szGot = Render( &target );

				if ( (szGot == nullptr) || (AgreementDigits( szGot, theCase.szExpect[ m ] ) != kExactAgreement) )
				{
					l_wrong++;
					szWhere = modeNames[ m ];
				}

				::free( szGot );
			}

			char	szLabel[ 72 ];
			char	szDetail[ 72 ];

			::snprintf( szLabel, 72, "%s / %s at %ld, seven modes", theCase.szLeft, theCase.szRight, theCase.lPlaces );

			if ( (l_wrong == 0) && (l_lost == 0) )
				::snprintf( szDetail, 72, "each rounded its way, each mode kept" );
			else
				::snprintf( szDetail, 72, "%ld wrong, one in %s; %ld lost their mode", l_wrong, szWhere, l_lost );

			CheckTrue( szLabel, (l_wrong == 0) && (l_lost == 0), szDetail );
		}
	}

	// Products at five places, where the rounding is easy to see...
	{
		struct { const char* szLeft; const char* szRight; const char* szExpect[ 7 ]; } products[] =
		{
		{ "1.23456789", "9.87654321", { "12.19326", "12.19326", "12.19326", "12.19327", "12.19326", "12.19326", "12.19327" } },
		{ "-1.23456789", "9.87654321", { "-12.19326", "-12.19326", "-12.19326", "-12.19326", "-12.19327", "-12.19326", "-12.19327" } },
		{ "0.00125", "0.5", { "0.00062", "0.00063", "0.00062", "0.00063", "0.00062", "0.00062", "0.00063" } },
		{ "-0.00125", "0.5", { "-0.00062", "-0.00063", "-0.00062", "-0.00062", "-0.00063", "-0.00062", "-0.00063" } },
		};

		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = 5;

		for ( auto& theCase : products )
		{
			long		l_wrong = 0;
			long		l_lost = 0;
			const char*	szWhere = "";

			for ( int m = 0; m < 7; m++ )
			{
				CLZNumber	a, b, target;

				Set( a, theCase.szLeft );
				Set( b, theCase.szRight );
				target.SetRoundingMode( modes[ m ] );

				lz_multiply( &a, &b, &target );

				l_lost += (target.m_rmMode != modes[ m ]);

				char*	szGot = Render( &target );

				if ( (szGot == nullptr) || (AgreementDigits( szGot, theCase.szExpect[ m ] ) != kExactAgreement) )
				{
					l_wrong++;
					szWhere = modeNames[ m ];
				}

				::free( szGot );
			}

			char	szLabel[ 72 ];
			char	szDetail[ 72 ];

			::snprintf( szLabel, 72, "%s * %s, seven modes", theCase.szLeft, theCase.szRight );

			if ( (l_wrong == 0) && (l_lost == 0) )
				::snprintf( szDetail, 72, "each rounded its way, each mode kept" );
			else
				::snprintf( szDetail, 72, "%ld wrong, one in %s; %ld lost their mode", l_wrong, szWhere, l_lost );

			CheckTrue( szLabel, (l_wrong == 0) && (l_lost == 0), szDetail );
		}

		glWorkingPrecision = l_saved;
	}

	// In place, the left operand is the target, and its mode governs...
	{
		CLZNumber	a, three;

		Set( a, "-2" );
		Set( three, "3" );
		a.SetRoundingMode( RM_round_ceiling );
		lz_divide( &a, &three, nullptr );
		CheckExact( "-2/3 in place, ceiling", &a, "-0.66666666666666666666666666666666666666" );
		CheckTrue( "the operand keeps its mode", a.m_rmMode == RM_round_ceiling, "" );
	}

	/*
		**Every kind of operation leaves its target's settings as it found
		them**: the four plain ones, a division by nought, a function
		and one refusing, the compound operations and the three means,
		and a not-a-number carried through division and multiplication.
		Each used to hand back its temporary's half-even and digit.
	*/

	{
		long		l_changed = 0;
		const char*	szWhere = "";

		auto Check = [ & ]( CLZNumber& theTarget, const char* szName )
		{
			if ( (theTarget.m_rmMode != RM_round_ceiling) || (theTarget.m_rdDigit != RD_digit_last) )
			{
				l_changed++;
				szWhere = szName;
			}
		};

		auto Fresh = []( CLZNumber& theTarget )
		{
			theTarget.SetRoundingMode( RM_round_ceiling );
			theTarget.SetRoundingDigit( RD_digit_last );
		};

		CLZNumber	a, b, zero, notANumber, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
		CLZNumber	t12, t13, t14, t15, t16;
		CLZNumber*	pValues[] = { &a, &b };

		Set( a, "7.5" );
		Set( b, "2" );
		Set( zero, "0" );
		Set( notANumber, "0" );
		notANumber.SetFlag( kfQNaN );

		Fresh( t1 );	lz_add( &a, &b, &t1 );				Check( t1, "add" );
		Fresh( t2 );	lz_subtract( &a, &b, &t2 );			Check( t2, "subtract" );
		Fresh( t3 );	lz_multiply( &a, &b, &t3 );			Check( t3, "multiply" );
		Fresh( t4 );	lz_divide( &a, &b, &t4 );			Check( t4, "divide" );
		Fresh( t5 );	lz_divide( &a, &zero, &t5 );		Check( t5, "divide by nought" );
		Fresh( t6 );	lz_sine( &a, &t6 );					Check( t6, "sine" );
		Fresh( t7 );	lz_modulo( &a, &b, &t7 );			Check( t7, "modulo" );
		Fresh( t8 );	lz_amean( pValues, 2, &t8 );		Check( t8, "amean" );
		Fresh( t9 );	lz_arcsine( &a, &t9 );				Check( t9, "arcsine, refused" );
		Fresh( t10 );	lz_divide( &notANumber, &b, &t10 );	Check( t10, "divide, a nan carried" );
		Fresh( t11 );	lz_multiply( &notANumber, &b, &t11 );	Check( t11, "multiply, a nan carried" );
		Fresh( t12 );	lz_gmean( pValues, 2, &t12 );		Check( t12, "gmean" );
		Fresh( t13 );	lz_hmean( pValues, 2, &t13 );		Check( t13, "hmean" );
		Fresh( t14 );	lz_divide_integer( &a, &b, &t14 );	Check( t14, "divide_integer" );
		Fresh( t15 );	lz_remainder( &a, &b, &t15 );		Check( t15, "remainder" );
		Fresh( t16 );	lz_remainder_near( &a, &b, &t16 );	Check( t16, "remainder_near" );

		char	szDetail[ 72 ];

		if ( l_changed == 0 )
			::snprintf( szDetail, 72, "after all sixteen operations" );
		else
			::snprintf( szDetail, 72, "%ld of 16 changed them, one in %s", l_changed, szWhere );

		CheckTrue( "a target keeps its mode and digit", l_changed == 0, szDetail );
	}

	/*
		**A function answers the same whatever the caller's settings.**
		Every sweep case was run with its operand, then its target, in
		each mode and each rounding digit: until 2026-09-26 the
		operand's mode moved 22 answers and its digit 58 -- arcsin,
		arccos, csc, cot and tan, arcsin(0.1) answering 0 at the units
		digit -- and now none move. This is that census, smaller: every
		one-operand function on four arguments, and six two-operand ones
		and the three means on five pairs, with the operands or the
		target set eight ways. The unset answers are taken once and each setting compared
		to them.
	*/

	{
		#define ENTRY( f ) { #f, f }

		struct { const char* szName; Unary fn; } unary[] =
		{
			ENTRY( lz_arccosecant ),		ENTRY( lz_arccosine ),			ENTRY( lz_arccotangent ),
			ENTRY( lz_archypercosecant ),	ENTRY( lz_archypercosine ),		ENTRY( lz_archypercotangent ),
			ENTRY( lz_archypersecant ),		ENTRY( lz_archypersine ),		ENTRY( lz_archypertangent ),
			ENTRY( lz_arcsecant ),			ENTRY( lz_arcsine ),			ENTRY( lz_arctangent ),
			ENTRY( lz_cosecant ),			ENTRY( lz_cosine ),				ENTRY( lz_cotangent ),
			ENTRY( lz_exp_b ),				ENTRY( lz_exp_e ),				ENTRY( lz_gamma ),
			ENTRY( lz_hypercosecant ),		ENTRY( lz_hypercosine ),		ENTRY( lz_hypercotangent ),
			ENTRY( lz_hypersecant ),		ENTRY( lz_hypersine ),			ENTRY( lz_hypertangent ),
			ENTRY( lz_log_b ),				ENTRY( lz_log_e ),				ENTRY( lz_log_gamma ),
			ENTRY( lz_reciprocal ),			ENTRY( lz_root_2 ),				ENTRY( lz_root_3 ),
			ENTRY( lz_secant ),				ENTRY( lz_sine ),				ENTRY( lz_tangent ),
		};

		struct { const char* szName; Binary fn; } binary[] =
		{
			ENTRY( lz_modulo ),				ENTRY( lz_remainder ),			ENTRY( lz_remainder_near ),
			ENTRY( lz_divide_integer ),		ENTRY( lz_power ),				ENTRY( lz_log_n ),
			{ "lz_amean", []( CLZNumber* x, CLZNumber* y, CLZNumber* r ) -> long
						  { CLZNumber* v[] = { x, y }; return lz_amean( v, 2, r ); } },
			{ "lz_gmean", []( CLZNumber* x, CLZNumber* y, CLZNumber* r ) -> long
						  { CLZNumber* v[] = { x, y }; return lz_gmean( v, 2, r ); } },
			{ "lz_hmean", []( CLZNumber* x, CLZNumber* y, CLZNumber* r ) -> long
						  { CLZNumber* v[] = { x, y }; return lz_hmean( v, 2, r ); } },
		};

		#undef ENTRY

		const char*	args[] = { "0.001", "0.1", "1.5", "-0.75" };

		/*
			The two-operand ones take 7.5 over each of those, and one
			pairs more. 1 over 0.333...334 to fifty places is 2.999...994,
			just short of a whole number, and 1 over 0.399...9 is
			2.5000...006, just past a half: where a whole number taken in
			a caller's mode would land differently. It was so when these
			operations divided into their target and rounded there; since
			they take the whole number straight from the raw quotient in
			a mode of their own -- see TestWholeQuotients -- no mode can
			reach it, and the pairs stay as the proof.
		*/

		const char*	pairs[][ 2 ] = { { "7.5", "0.001" }, { "7.5", "0.1" }, { "7.5", "1.5" }, { "7.5", "-0.75" },
									 { "1", "0.33333333333333333333333333333333333333333333333334" },
									 { "1", "0.39999999999999999999999999999999999999999999999999" } };

		struct { const char* szName; bool bTarget; bool bDigit; long lSetting; } settings[] =
		{
			{ "operand rounding down",    false, false, RM_round_down },
			{ "operand ceiling",          false, false, RM_round_ceiling },
			{ "operand floor",            false, false, RM_round_floor },
			{ "operand rounding up",      false, false, RM_round_up },
			{ "operand at the units",     false, true,  RD_digit_zero },
			{ "operand at its last",      false, true,  RD_digit_last },
			{ "target floor",             true,  false, RM_round_floor },
			{ "target at the units",      true,  true,  RD_digit_zero },
		};

		auto Apply = [ & ]( CLZNumber& theNumber, int iSetting )
		{
			if ( settings[ iSetting ].bDigit )
				theNumber.SetRoundingDigit( static_cast<RoundingDigit>( settings[ iSetting ].lSetting ) );
			else
				theNumber.SetRoundingMode( static_cast<RoundingMode>( settings[ iSetting ].lSetting ) );
		};

		// The answer as a string, a code and the flags, so a refusal is
		//	compared as carefully as an answer...
		auto Answer = []( CLZNumber& theAnswer, long lCode ) -> std::string
		{
			char*		szGot = theAnswer.GetStringRep();
			std::string	s_answer = std::string( (szGot != nullptr) ? szGot : "" ) + "|" + std::to_string( lCode ) + "|" +
								   std::to_string( theAnswer.m_lFlags & (kfQNaN | kfInfinity | kfDiv0) );

			::free( szGot );

			return s_answer;
		};

		std::vector<std::string>	plain;

		for ( const char* szArg : args )
		{
			for ( auto& theFn : unary )
			{
				CLZNumber	x, answer;

				x.Assign( szArg );
				plain.push_back( Answer( answer, theFn.fn( &x, &answer ) ) );
			}
		}

		for ( auto& thePair : pairs )
		{
			for ( auto& theFn : binary )
			{
				CLZNumber	x, y, answer;

				x.Assign( thePair[ 0 ] );
				y.Assign( thePair[ 1 ] );
				plain.push_back( Answer( answer, theFn.fn( &x, &y, &answer ) ) );
			}
		}

		for ( int s = 0; s < static_cast<int>( sizeof( settings ) / sizeof( settings[ 0 ] ) ); s++ )
		{
			long	l_moved = 0;
			long	l_cases = 0;
			char	szWhere[ 64 ] = "";

			for ( const char* szArg : args )
			{
				for ( auto& theFn : unary )
				{
					CLZNumber	x, answer;

					x.Assign( szArg );

					if ( settings[ s ].bTarget )
						Apply( answer, s );
					else
						Apply( x, s );

					if ( Answer( answer, theFn.fn( &x, &answer ) ) != plain[ static_cast<size_t>( l_cases ) ] )
					{
						l_moved++;
						::snprintf( szWhere, 64, "%s(%s)", theFn.szName, szArg );
					}

					l_cases++;
				}
			}

			for ( auto& thePair : pairs )
			{
				for ( auto& theFn : binary )
				{
					CLZNumber	x, y, answer;

					x.Assign( thePair[ 0 ] );
					y.Assign( thePair[ 1 ] );

					if ( settings[ s ].bTarget )
						Apply( answer, s );
					else
					{
						Apply( x, s );
						Apply( y, s );
					}

					if ( Answer( answer, theFn.fn( &x, &y, &answer ) ) != plain[ static_cast<size_t>( l_cases ) ] )
					{
						l_moved++;
						::snprintf( szWhere, 64, "%s(%s, %.8s...)", theFn.szName, thePair[ 0 ], thePair[ 1 ] );
					}

					l_cases++;
				}
			}

			char	szLabel[ 72 ];
			char	szDetail[ 96 ];

			::snprintf( szLabel, 72, "functions ignore: %s", settings[ s ].szName );

			if ( l_moved == 0 )
				::snprintf( szDetail, 96, "all %ld answers the same", l_cases );
			else
				::snprintf( szDetail, 96, "%ld of %ld moved, one at %s", l_moved, l_cases, szWhere );

			CheckTrue( szLabel, l_moved == 0, szDetail );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	A target at the units. Since 2026-10-01 the four plain operations
//	answer a whole number into a target whose rounding digit is
//	RD_digit_zero, rounded once from the exact value in the target's
//	mode. Until then they rounded at their own width and kept the
//	digit, and a whole number taken afterwards was a second rounding.
//	See lzi_asks_whole.
/////////////////////////////////////////////////////////////////

static void TestTargetAtTheUnits()
{
	Banner( "A target at the units: a whole number, rounded once" );

	const RoundingMode	modes[] = { RM_round_down, RM_round_half_up, RM_round_half_even, RM_round_ceiling,
									RM_round_floor, RM_round_half_down, RM_round_half_odd, RM_round_up };
	const char*			modeNames[] = { "down", "half up", "half even", "ceiling", "floor", "half down", "half odd", "up" };

	/*
		Each exact value rounded to a whole number in eight modes, by
		Python on the exact fraction. A tie goes by the kept digit, as
		the standard defines half-even and half-odd; in an odd base that
		is not the whole number's parity.

		The first four were wrong by the only routes there were. Each is
		a hair from a half and rounds to exactly a half at its own width
		-- 1.5811...267 squared is 2.5 + 4.4e-37 -- so a whole number
		taken afterwards, by Round or lz_round_integral, rounded twice:
		4 for the first three in half-even, 2 for the fourth. The rest
		pin exact ties, both signs, all four operations, and bases seven,
		sixteen and two.
	*/

	struct { char cOp; long lBase; const char* szLeft; const char* szRight; const char* szExpect[ 8 ]; } cases[] =
	{
		{ '/', 10, "24.4999999999999999999999999999999999999993", "7", { "3", "3", "3", "4", "3", "3", "3", "4" } },
		{ '/', 10, "6.99999999999999999999999999999999999999998", "2", { "3", "3", "3", "4", "3", "3", "3", "4" } },
		{ '*', 10, "1.75", "1.99999999999999999999999999999999999999", { "3", "3", "3", "4", "3", "3", "3", "4" } },
		{ '*', 10, "1.581138830084189665999446772216359267", "1.581138830084189665999446772216359267", { "2", "3", "3", "3", "2", "3", "3", "3" } },
		{ '/', 10, "5", "2", { "2", "3", "2", "3", "2", "2", "3", "3" } },
		{ '/', 10, "-7", "2", { "-3", "-4", "-4", "-3", "-4", "-3", "-3", "-4" } },
		{ '+', 10, "1.25", "1.25", { "2", "3", "2", "3", "2", "2", "3", "3" } },
		{ '-', 10, "-1.25", "1.25", { "-2", "-3", "-2", "-2", "-3", "-2", "-3", "-3" } },
		{ '*', 10, "0.5", "-7", { "-3", "-4", "-4", "-3", "-4", "-3", "-3", "-4" } },
		{ '/', 10, "22", "7", { "3", "3", "3", "4", "3", "3", "3", "4" } },
		{ '/', 10, "-22", "7", { "-3", "-3", "-3", "-3", "-4", "-3", "-3", "-4" } },
		{ '/', 7, "5", "2", { "2", "3", "2", "3", "2", "2", "3", "3" } },
		{ '/', 7, "10", "2", { "3", "4", "4", "4", "3", "3", "3", "4" } },
		{ '*', 7, "1.3", "1.3", { "2", "2", "2", "3", "2", "2", "2", "3" } },
		{ '+', 16, "1.8", "-3", { "-1", "-2", "-2", "-1", "-2", "-1", "-1", "-2" } },
		{ '/', 2, "1011", "10", { "101", "110", "110", "110", "101", "101", "101", "110" } },
	};

	for ( auto& theCase : cases )
	{
		long		l_wrong = 0;
		long		l_lost = 0;
		const char*	szWhere = "";

		for ( int m = 0; m < 8; m++ )
		{
			CLZNumber	a, b, target;

			SetInBase( a, theCase.lBase, theCase.szLeft );
			SetInBase( b, theCase.lBase, theCase.szRight );
			target.SetBase( theCase.lBase, true );
			target.SetRoundingMode( modes[ m ] );
			target.SetRoundingDigit( RD_digit_zero );

			switch ( theCase.cOp )
			{
			case '+':	lz_add( &a, &b, &target );		break;
			case '-':	lz_subtract( &a, &b, &target );	break;
			case '*':	lz_multiply( &a, &b, &target );	break;
			default:	lz_divide( &a, &b, &target );	break;
			}

			l_lost += (target.m_rmMode != modes[ m ]) || (target.m_rdDigit != RD_digit_zero);

			char*	szGot = target.GetStringRep();

			if ( (szGot == nullptr) || (::strcmp( szGot, theCase.szExpect[ m ] ) != 0) )
			{
				l_wrong++;
				szWhere = modeNames[ m ];
			}

			::free( szGot );
		}

		char	szLabel[ 72 ];
		char	szDetail[ 72 ];

		::snprintf( szLabel, 72, "%.12s%s %c %.12s%s, base %ld", theCase.szLeft, (::strlen( theCase.szLeft ) > 12) ? "..." : "",
					theCase.cOp, theCase.szRight, (::strlen( theCase.szRight ) > 12) ? "..." : "", theCase.lBase );

		if ( (l_wrong == 0) && (l_lost == 0) )
			::snprintf( szDetail, 72, "whole in all eight modes, settings kept" );
		else
			::snprintf( szDetail, 72, "%ld wrong, one in %s; %ld lost a setting", l_wrong, szWhere, l_lost );

		CheckTrue( szLabel, (l_wrong == 0) && (l_lost == 0), szDetail );
	}

	// In place, the left operand is the target, and its digit governs --
	//	the explicit-precision division too...
	{
		CLZNumber	a, b, c, d, e, f, g, seven;

		Set( seven, "7" );
		Set( a, "22" );		a.SetRoundingDigit( RD_digit_zero );
		lz_divide( &a, &seven, nullptr );
		CheckExact( "22/7 in place, at the units", &a, "3" );

		Set( b, "22" );		b.SetRoundingDigit( RD_digit_zero );	b.SetRoundingMode( RM_round_ceiling );
		lz_divide( &b, &seven, 36L );
		CheckExact( "22/7 to 36 places in place, ceiling", &b, "4" );

		Set( c, "1.5" );	Set( d, "1.5" );	c.SetRoundingDigit( RD_digit_zero );
		lz_multiply( &c, &d, nullptr );
		CheckExact( "1.5 * 1.5 in place, at the units", &c, "2" );

		Set( e, "1.25" );	Set( f, "1.25" );	e.SetRoundingDigit( RD_digit_zero );	e.SetRoundingMode( RM_round_half_up );
		lz_add( &e, &f, nullptr );
		CheckExact( "1.25 + 1.25 in place, half up", &e, "3" );

		Set( g, "-1.25" );	g.SetRoundingDigit( RD_digit_zero );
		lz_subtract( &g, &f, nullptr );
		CheckExact( "-1.25 - 1.25 in place, at the units", &g, "-2" );
	}

	/*
		**Only the units.** 'Last' means a place a number already holds,
		so an operation leaves it to Round; and a target set not to round
		gets no rounding at any digit, where the whole-number division
		would have handed back its working digits. Each answers exactly
		what a target with the default digit does.
	*/

	{
		struct { const char* szName; RoundingMode rmMode; RoundingDigit rdDigit; } settings[] =
		{
			{ "a target at its last place", RM_round_half_even, RD_digit_last },
			{ "a target set not to round, at the units", RM_round_noop, RD_digit_zero },
		};

		for ( auto& theSetting : settings )
		{
			long	l_moved = 0;

			for ( char cOp : { '+', '*', '/' } )
			{
				CLZNumber	a, b, plain, target;

				Set( a, "22" );
				Set( b, (cOp == '/') ? "7" : "0.142857142857142857142857142857142857" );
				plain.SetRoundingMode( theSetting.rmMode );
				target.SetRoundingMode( theSetting.rmMode );
				target.SetRoundingDigit( theSetting.rdDigit );

				if ( cOp == '+' )		{ lz_add( &a, &b, &plain );			lz_add( &a, &b, &target ); }
				else if ( cOp == '*' )	{ lz_multiply( &a, &b, &plain );	lz_multiply( &a, &b, &target ); }
				else					{ lz_divide( &a, &b, &plain );		lz_divide( &a, &b, &target ); }

				char*	szPlain = plain.GetStringRep();
				char*	szTarget = target.GetStringRep();

				l_moved += (::strcmp( szPlain, szTarget ) != 0);

				::free( szPlain );
				::free( szTarget );
			}

			char	szDetail[ 72 ];

			::snprintf( szDetail, 72, "%ld of 3 operations moved", l_moved );
			CheckTrue( theSetting.szName, l_moved == 0, szDetail );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	A short answer is an exact one. Since 2026-10-01 an answer written
//	with fewer places than its width ended there, and one that did not
//	end is written out to its width, trailing noughts and all. Until
//	then the spelling followed the path the answer took: sin(1e-42) was
//	0 and arcsin(1e-42) 0.000...0, exp(1e-42) was 1, and a quotient a
//	hair from 22.5 was 22.5. See lzi_write_to_width.
/////////////////////////////////////////////////////////////////

static void TestShortMeansExact()
{
	Banner( "A short answer is an exact one" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );
	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );

	// The places an answer is written with, and its width: the precision,
	//	or for a small answer of a function ending in a quotient the
	//	precision plus its noughts...
	auto Places = []( CLZNumber& theAnswer ) -> long
	{
		return theAnswer.m_lRadix;
	};

	auto Width = []( CLZNumber& theAnswer, bool bQuotient ) -> long
	{
		long	l_width = glWorkingPrecision;

		theAnswer.CheckForZero();

		if ( bQuotient && !(theAnswer.IsZero()) )
		{
			long	l_top = theAnswer.m_lDigitsLength - 1;

			while ( (l_top >= 0) && (theAnswer.m_lpDigits[ l_top ] == 0) )
				l_top--;

			long	l_exponent = l_top - theAnswer.m_lRadix;

			if ( l_exponent < -1 )
				l_width += -l_exponent - 1;
		}

		return l_width;
	};

	auto Flagged = []( CLZNumber& theAnswer ) -> bool
	{
		return theAnswer.FlagIsOn( kfQNaN ) || theAnswer.FlagIsOn( kfInfinity );
	};

	/*
		**The rule itself, needing no reference.** Every answer written
		shorter than its width at 36 places is taken again at 96 and has
		to be the same value: an exact answer does not move with the
		width, and one that only rounded short does -- exp(1e-42) at 96
		places is 1.000...0001. The arguments are the shapes that wrote
		short answers which were not exact: tiny ones, ones a hair from a
		special value, and 97, where tanh is 1 to 84 places. The exact
		ones beside them make sure short answers are examined at all, and
		the count of them is required too.
	*/

	{
		#define ENTRY( f, q ) { #f, f, q }

		struct { const char* szName; Unary fn; bool bQuotient; } unary[] =
		{
			ENTRY( lz_sine, false ),			ENTRY( lz_cosine, false ),			ENTRY( lz_tangent, true ),
			ENTRY( lz_secant, false ),			ENTRY( lz_cosecant, false ),		ENTRY( lz_cotangent, true ),
			ENTRY( lz_arcsine, false ),			ENTRY( lz_arccosine, false ),		ENTRY( lz_arctangent, false ),
			ENTRY( lz_arcsecant, false ),		ENTRY( lz_arccosecant, false ),		ENTRY( lz_arccotangent, false ),
			ENTRY( lz_hypersine, true ),		ENTRY( lz_hypercosine, true ),		ENTRY( lz_hypertangent, true ),
			ENTRY( lz_hypersecant, true ),		ENTRY( lz_hypercosecant, true ),	ENTRY( lz_hypercotangent, true ),
			ENTRY( lz_archypersine, false ),	ENTRY( lz_archypercosine, false ),	ENTRY( lz_archypertangent, false ),
			ENTRY( lz_archypersecant, false ),	ENTRY( lz_archypercosecant, false ),	ENTRY( lz_archypercotangent, false ),
			ENTRY( lz_exp_e, false ),			ENTRY( lz_raise_e, false ),			ENTRY( lz_log_e, false ),
			ENTRY( lz_exp_b, false ),			ENTRY( lz_log_b, false ),			ENTRY( lz_root_2, false ),
			ENTRY( lz_root_3, true ),			ENTRY( lz_reciprocal, true ),		ENTRY( lz_square, false ),
			ENTRY( lz_cube, false ),			ENTRY( lz_gamma, false ),
		};

		#undef ENTRY

		const char*	args[] =
		{
			"1E-42", "-1E-42", "3E-20", "4.000000000000000000000000000000000000000001", "8.000000000000000000000000000000000000000001",
			"0.99999999999999999999999999999999999999999999", "1.000000000000000000000000000000000000000001",
			"97", "0", "1", "-8", "4", "6.25", "0.5", "0.001",
		};

		long		l_short = 0;
		long		l_moved = 0;
		char		szWhere[ 72 ] = "";

		for ( auto& theFn : unary )
		{
			for ( const char* szArg : args )
			{
				CLZNumber	x, answer;

				x.Assign( szArg );

				if ( (theFn.fn( &x, &answer ) != kNoError) || Flagged( answer ) ||
					 (Places( answer ) >= Width( answer, theFn.bQuotient )) )
					continue;

				l_short++;

				long	l_saved = glWorkingPrecision;

				glWorkingPrecision = l_saved + 60;

				CLZNumber	wideX, wide;

				wideX.Assign( szArg );
				theFn.fn( &wideX, &wide );
				glWorkingPrecision = l_saved;

				if ( lz_compare( &answer, &wide ) != kEqualTo )
				{
					l_moved++;
					::snprintf( szWhere, 72, "%s(%.12s)", theFn.szName, szArg );
				}
			}
		}

		char	szDetail[ 96 ];

		::snprintf( szDetail, 96, "%ld short answers, %ld moved at 96 places%s%s", l_short, l_moved,
					(l_moved > 0) ? ", one " : "", szWhere );
		CheckTrue( "one operand: short is exact", (l_moved == 0) && (l_short >= 60), szDetail );
	}

	// The same of the powers, the logarithm to a base, the roots and the
	//	plain division and product...
	{
		struct { const char* szName; Binary fn; const char* szLeft; const char* szRight; } binary[] =
		{
			{ "power", lz_power, "4.000000000000000000000000000000000000000001", "0.5" },			{ "power", lz_power, "4", "0.5" },
			{ "power", lz_power, "4", "1.5" },				{ "power", lz_power, "8", "-1" },
			{ "exp_n", lz_exp_n, "4.000000000000000000000000000000000000000001", "0.5" },			{ "exp_n", lz_exp_n, "16", "0.25" },
			{ "raise_n", lz_raise_n, "2", "-3" },			{ "raise_n", lz_raise_n, "1.5", "2" },
			{ "log_n", lz_log_n, "2", "8.000000000000000000000000000000000000000001" },				{ "log_n", lz_log_n, "4", "2" },
			{ "log_n", lz_log_n, "2", "0.125" },			{ "log_n", lz_log_n, "3", "3" },
			{ "root_n", lz_root_n, "8.000000000000000000000000000000000000000001", "3" },			{ "root_n", lz_root_n, "4", "2" },
			{ "divide", lz_divide, "83.24294", "50.09425" },	{ "divide", lz_divide, "1", "4" },
			{ "divide", lz_divide, "24.4999999999999999999999999999999999999993", "7" },
			{ "multiply", lz_multiply, "1.5", "1.5" },
			{ "multiply", lz_multiply, "1.581138830084189665999446772216359267", "1.581138830084189665999446772216359267" },
		};

		long		l_short = 0;
		long		l_moved = 0;
		char		szWhere[ 72 ] = "";

		for ( auto& theCase : binary )
		{
			CLZNumber	x, y, answer;

			x.Assign( theCase.szLeft );
			y.Assign( theCase.szRight );

			if ( (theCase.fn( &x, &y, &answer ) != kNoError) || Flagged( answer ) ||
				 (Places( answer ) >= Width( answer, false )) )
				continue;

			l_short++;

			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = l_saved + 60;

			CLZNumber	wideX, wideY, wide;

			wideX.Assign( theCase.szLeft );
			wideY.Assign( theCase.szRight );
			theCase.fn( &wideX, &wideY, &wide );
			glWorkingPrecision = l_saved;

			if ( lz_compare( &answer, &wide ) != kEqualTo )
			{
				l_moved++;
				::snprintf( szWhere, 72, "%s(%.10s, %.10s)", theCase.szName, theCase.szLeft, theCase.szRight );
			}
		}

		char	szDetail[ 96 ];

		::snprintf( szDetail, 96, "%ld short answers, %ld moved at 96 places%s%s", l_short, l_moved,
					(l_moved > 0) ? ", one " : "", szWhere );
		CheckTrue( "two operands: short is exact", (l_moved == 0) && (l_short >= 10), szDetail );
	}

	/*
		**Pins, from Python's decimal.** Answers that are not exact and
		used to be written short, now out to their width; and exact
		answers, whose spelling did not move. log10(100) is exact and
		written out: the logarithm cannot tell, and a written-out answer
		claims only to be rounded.
	*/

	{
		struct { const char* szLabel; Unary fn; const char* szArg; const char* szExpect; } pins[] =
		{
			{ "sin(1e-42), a rounded nought", lz_sine, "1E-42", "0.000000000000000000000000000000000000" },
			{ "cos(1e-42), a rounded one", lz_cosine, "1E-42", "1.000000000000000000000000000000000000" },
			{ "exp(1e-42)", lz_exp_e, "1E-42", "1.000000000000000000000000000000000000" },
			{ "tanh(97)", lz_hypertangent, "97", "1.000000000000000000000000000000000000" },
			{ "arcsin(3e-20), the operand back", lz_arcsine, "3E-20", "0.000000000000000000030000000000000000" },
			{ "tan(1e-42), a quotient's width", lz_tangent, "1E-42", "0.00000000000000000000000000000000000000000100000000000000000000000000000000000" },
			{ "sqrt(4 + 1e-42)", lz_root_2, "4.000000000000000000000000000000000000000001", "2.000000000000000000000000000000000000" },
			{ "log10(100), not proven exact", lz_log_b, "100", "2.000000000000000000000000000000000000" },
			{ "sqrt(4), exact", lz_root_2, "4", "2" },
			{ "sqrt(6.25), exact", lz_root_2, "6.25", "2.5" },
			{ "cbrt(-8), exact", lz_root_3, "-8", "-2" },
			{ "1/4, exact", lz_reciprocal, "4", "0.25" },
			{ "10^-3, exact", lz_exp_b, "-3", "0.001" },
			{ "sin(0), exact", lz_sine, "0", "0" },
			{ "cos(0), exact", lz_cosine, "0", "1" },
			{ "exp(0), exact", lz_exp_e, "0", "1" },
			{ "ln(1), exact", lz_log_e, "1", "0" },
			{ "arccosh(1), exact", lz_archypercosine, "1", "0" },
			{ "gamma(5), exact", lz_gamma, "5", "24" },
		};

		for ( auto& thePin : pins )
		{
			CLZNumber	x, answer;

			x.Assign( thePin.szArg );
			thePin.fn( &x, &answer );
			CheckRendering( thePin.szLabel, &answer, thePin.szExpect );
		}

		struct { const char* szLabel; Binary fn; const char* szLeft; const char* szRight; const char* szExpect; } pairs[] =
		{
			{ "4^0.5 a hair over, not exact", lz_power, "4.000000000000000000000000000000000000000001", "0.5", "2.000000000000000000000000000000000000" },
			{ "83.24294 / 50.09425, guard nought", lz_divide, "83.24294", "50.09425", "1.6617264456499498445430363764304286420" },
			{ "4^0.5, exact", lz_power, "4", "0.5", "2" },
			{ "4^1.5, exact", lz_power, "4", "1.5", "8" },
			{ "2^-3, exact", lz_raise_n, "2", "-3", "0.125" },
			{ "log_4 2, exact", lz_log_n, "4", "2", "0.5" },
			{ "log_2 0.125, exact", lz_log_n, "2", "0.125", "-3" },
		};

		for ( auto& thePair : pairs )
		{
			CLZNumber	x, y, answer;

			x.Assign( thePair.szLeft );
			y.Assign( thePair.szRight );
			thePair.fn( &x, &y, &answer );
			CheckRendering( thePair.szLabel, &answer, thePair.szExpect );
		}

		// The means check a short answer against their values; the target
		//	may be one of them...
		{
			CLZNumber	a, b, c, mean;
			CLZNumber*	pList[] = { &a, &b, &c };

			a.Assign( "2" );		b.Assign( "8.000000000000000000000000000000000000000001" );
			lz_gmean( pList, 2, &mean );
			CheckRendering( "gmean(2, 8 + 1e-42), not exact", &mean, "4.000000000000000000000000000000000000" );

			a.Assign( "1.2E100" );	b.Assign( "7.5" );	c.Assign( "3E100" );
			lz_hmean( pList, 3, &mean );
			CheckRendering( "hmean(1.2e100, 7.5, 3e100)", &mean, "22.500000000000000000000000000000000000" );

			a.Assign( "2" );		b.Assign( "8" );
			lz_gmean( pList, 2, &mean );
			CheckRendering( "gmean(2, 8), exact", &mean, "4" );

			a.Assign( "3" );		b.Assign( "6" );
			lz_hmean( pList, 2, &a );
			CheckRendering( "hmean(3, 6) into the first, exact", &a, "4" );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	A root under the working width. Until 2026-10-02 the square root of
//	anything below about 1e-83 answered 0 flagged infinite and divided by
//	zero, with no error: its iterate rounded to nought at the working
//	width, and the next pass divided by it. arcsin and arccos near one,
//	which take the root of 1 - x^2, answered 0 flagged not-a-number --
//	arcsin(1 - 1e-100) in place of pi/2 -- and arccosh near one carried
//	the same flag on the right value. See lz_root_2.
/////////////////////////////////////////////////////////////////

static void TestRootUnderTheWidth()
{
	Banner( "A root under the working width" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	// No infinity, no not-a-number, no division by nought: what a string
	//	comparison cannot see, and the whole of the fault for half these...
	auto Clean = []( CLZNumber& theAnswer ) -> bool
	{
		return theAnswer.FlagIsOff( kfInfinity ) && theAnswer.FlagIsOff( kfQNaN ) && theAnswer.FlagIsOff( kfDiv0 );
	};

	/*
		Each value from Python, correctly rounded at 36 places, and each
		required unflagged. The square roots are rounded noughts, written
		out, from the last that worked (1e-83) on; pi/2 and pi are what
		arcsin and arccos answer a hair inside one either way, in base
		ten and base seven.
	*/

	{
		struct { const char* szLabel; Unary fn; long lBase; const char* szArg; const char* szExpect; } pins[] =
		{
			{ "sqrt(1e-83), the last that worked", lz_root_2, 10, "1E-83", "0.000000000000000000000000000000000000" },
			{ "sqrt(1e-84)", lz_root_2, 10, "1E-84", "0.000000000000000000000000000000000000" },
			{ "sqrt(1e-200)", lz_root_2, 10, "1E-200", "0.000000000000000000000000000000000000" },
			{ "sqrt(16^-100), base 16", lz_root_2, 16, "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "0.000000000000000000000000000000000000" },
			{ "arcsin(1 - 1e-100)", lz_arcsine, 10, "0.9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999", "1.570796326794896619231321691639751442" },
			{ "arcsin(-(1 - 1e-100))", lz_arcsine, 10, "-0.9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999", "-1.570796326794896619231321691639751442" },
			{ "arccos(1 - 1e-150)", lz_arccosine, 10, "0.999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999", "0.000000000000000000000000000000000000" },
			{ "arccos(-(1 - 1e-150))", lz_arccosine, 10, "-0.999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999", "3.141592653589793238462643383279502884" },
			{ "arccosh(1 + 1e-103)", lz_archypercosine, 10, "1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "0.000000000000000000000000000000000000" },
			{ "arcsin(1 - 7^-100), base 7", lz_arcsine, 7, "0.6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666", "1.366532421435305204013153445566111332" },
		};

		for ( auto& thePin : pins )
		{
			CLZNumber	x, answer;

			x.SetBase( thePin.lBase, true );
			x.Assign( thePin.szArg );
			answer.SetBase( thePin.lBase, true );

			long	lResult = thePin.fn( &x, &answer );
			char*	szGot = answer.GetStringRep();
			bool	b_right = (szGot != nullptr) && (::strcmp( szGot, thePin.szExpect ) == 0);
			char	szDetail[ 96 ];

			::snprintf( szDetail, 96, "%s%s, error %ld, flags 0x%lx", b_right ? "exact" : "got ", b_right ? "" : szGot,
						lResult, static_cast<unsigned long>( answer.m_lFlags ) );
			CheckTrue( thePin.szLabel, b_right && (lResult == kNoError) && Clean( answer ), szDetail );

			::free( szGot );
		}
	}

	/*
		**Every root under half a place is a clean nought**, needing no
		reference: b^-k for k from 75 to 400 in bases 2, 7, 10 and 16,
		where the root is b^-k/2, under half of the 36th place.
	*/

	{
		long		l_cases = 0;
		long		l_wrong = 0;
		char		szWhere[ 48 ] = "";

		for ( long lBase : { 2L, 7L, 10L, 16L } )
		{
			for ( long k = 75; k <= 400; k += 13 )
			{
				std::string	s_value = "0." + std::string( static_cast<size_t>( k - 1 ), '0' ) + "1";
				CLZNumber	x, answer;

				x.SetBase( lBase, true );
				x.Assign( s_value.c_str() );
				answer.SetBase( lBase, true );

				long	lResult = lz_root_2( &x, &answer );

				answer.CheckForZero();
				l_cases++;

				if ( (lResult != kNoError) || !Clean( answer ) || !(answer.IsZero()) )
				{
					l_wrong++;
					::snprintf( szWhere, 48, "base %ld, k %ld", lBase, k );
				}
			}
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, 72, "%ld of %ld wrong%s%s", l_wrong, l_cases, (l_wrong > 0) ? ", one at " : "", szWhere );
		CheckTrue( "every root under half a place", l_wrong == 0, szDetail );
	}
}

/////////////////////////////////////////////////////////////////
//	The roots climb through widths that double. Since 2026-10-02 the
//	square and cube roots work in a band, their root between one and
//	the base, and take one Newton pass at each width from 64 places up
//	to the full one, where they used to work at full width every pass:
//	sqrt(2) at 10,000 places went from 69 ms to 9.5. Every check in the
//	suite before this ran at 36 places, a single width, so none of them
//	reaches the climb. See lzi_root_widths.
/////////////////////////////////////////////////////////////////

static void TestRootsByDoubling()
{
	Banner( "The roots climb through widths that double" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	/*
		Past the seed width, so the climb is taken: each from Python,
		correctly rounded at its places, or for the cube root's small
		answer at the rule for a quotient. A large root, a small one, a
		negative cube root, and bases seven and sixteen.
	*/

	{
		struct { const char* szLabel; Unary fn; long lBase; long lPlaces; const char* szArg; const char* szExpect; } pins[] =
		{
			{ "sqrt(2) at 300", lz_root_2, 10, 300, "2",
			  "1.41421356237309504880168872420969807856967187537694807317667973799073"
			  "2478462107038850387534327641572735013846230912297024924836055850737212"
			  "6441214970999358314132226659275055927557999505011527820605714701095599"
			  "7160597027453459686201472851741864088919860955232923048430871432145083"
			  "9762603627995251407990" },
			{ "sqrt(2E100) at 150", lz_root_2, 10, 150, "2E100",
			  "141421356237309504880168872420969807856967187537694.807317667973799073"
			  "2478462107038850387534327641572735013846230912297024924836055850737212"
			  "64412149709993583141322266592750559275579995050115278206057147" },
			{ "sqrt(3E-100) at 200", lz_root_2, 10, 200, "3E-100",
			  "0.00000000000000000000000000000000000000000000000001732050807568877293"
			  "5274463415058723669428052538103806280558069794519330169088000370811461"
			  "86757248575675626141415406703029969945094998952478811655512094" },
			{ "sqrt(3) at 150, base 7", lz_root_2, 7, 150, "3",
			  "1.50604402141016645614000153150653453234661520155424556352313660021662"
			  "5003643561251366033554116223663315125543541005620631256443246603564021"
			  "446324153545" },
			{ "cbrt(2) at 300", lz_root_3, 10, 300, "2",
			  "1.25992104989487316476721060727822835057025146470150798008197511215529"
			  "9676513959483729396562436255094154310256035615665259399024040613737228"
			  "4591103042693552469606426166250009774745265654803068671854055186892458"
			  "7251676419937370969509838278316139915512931369536618394746344857657030"
			  "3119095895984741105981" },
			{ "cbrt(A) at 150, base 16", lz_root_3, 16, 150, "A",
			  "2.278908270E09D951445AE49BD41233E7FA91197DB36756B15828F976EBF7B0575A59"
			  "97FDCA18E780125D496485D0C2E1DD2639B8F071A51849937FDC408045B345BAC32F89"
			  "E14BA1AD1C31" },
			{ "cbrt(-5E-90) at 200", lz_root_3, 10, 200, "-5E-90",
			  "-0.0000000000000000000000000000017099759466766969893531088725438601098"
			  "6805511054305492438286170744429592050417321625718701002018900220450328"
			  "9390454018087197576646800785374713893894415032333156031462264979150704"
			  "3419693392331556442949" },
			{ "cbrt(1.2345678...) at 200", lz_root_3, 10, 200, "1.23456789012345678901234567890123456",
			  "1.07276597967684621673776909867671517676881897476660223969600064547611"
			  "2620578996502417265258449045147046966521309224373768503543531971578325"
			  "61009564818051247519231444074754127824461730960362006668029634" },
		};

		for ( auto& thePin : pins )
		{
			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = thePin.lPlaces;

			CLZNumber	x, answer;

			x.SetBase( thePin.lBase, true );
			x.Assign( thePin.szArg );
			answer.SetBase( thePin.lBase, true );
			thePin.fn( &x, &answer );
			glWorkingPrecision = l_saved;

			CheckRendering( thePin.szLabel, &answer, thePin.szExpect );
		}
	}

	/*
		**The cube root settles a two-pass cycle.** Its iterate can land a
		unit either side of the root and alternate for ever; until
		2026-10-02 it then ran to the pass cap, 126 divisions for this
		one at 200 places, where it now makes 22. Counted with the
		packed divider's own counter, so the check needs no clock.
	*/

	{
		long	l_saved = glWorkingPrecision;
		long	l_before = 0;
		long	l_after = 0;
		long	l_addBacks = 0;

		glWorkingPrecision = 200;

		CLZNumber	x, answer;

		x.Assign( "1.23456789012345678901234567890123456" );
		lzi_packed_divide_counts( &l_before, &l_addBacks );
		lz_root_3( &x, &answer );
		lzi_packed_divide_counts( &l_after, &l_addBacks );
		glWorkingPrecision = l_saved;

		char	szDetail[ 72 ];

		::snprintf( szDetail, 72, "%ld divisions, at most 40 wanted", l_after - l_before );
		CheckTrue( "cbrt settles a cycle at 200 places", (l_after - l_before) <= 40, szDetail );
	}
}

/////////////////////////////////////////////////////////////////
//	A root near a rounding half is decided exactly. Since 2026-10-02 the
//	square and cube roots round by their guard digits only where the
//	guard is clear of the half, and elsewhere compare the half itself,
//	raised to the power, with the operand. Roots built a hair from a
//	half came out a unit wrong in 241 of 480 before, settled or not: no
//	guard can tell a root within a unit of its last place of the half.
//	See lzi_root_decide.
/////////////////////////////////////////////////////////////////

static void TestRootsDecided()
{
	Banner( "A root near a rounding half is decided exactly" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	const char*		kDigitChars = "0123456789ABCDEF";
	unsigned long long	l_seed = 20261002;

	// A digit from a fixed generator, nought allowed or not...
	auto Digit = [ & ]( long lBase, bool bNought ) -> char
	{
		l_seed = (l_seed * 6364136223846793005ULL) + 1442695040888963407ULL;

		return bNought ? kDigitChars[ (l_seed >> 33) % lBase ] : kDigitChars[ 1 + ((l_seed >> 33) % (lBase - 1)) ];
	};

	auto Digits = [ & ]( long lBase, long lCount ) -> std::string
	{
		std::string	s_digits;

		for ( long i = 0; i < lCount; i++ )
			s_digits += Digit( lBase, true );

		return s_digits;
	};

	// A run of a digit, from the place after lPlace on...
	auto Past = []( long lPlace, const std::string& sRun ) -> std::string
	{
		return "0." + std::string( static_cast<size_t>( lPlace ), '0' ) + sRun;
	};

	/*
		A kept value y with its places: the places the root rounds at,
		which for the square root are W after the point however small
		the answer, and for the cube root W after its leading noughts,
		by the rule for a quotient. iMag is the integer digits past the
		first, or less than nought for that many noughts after the point.
	*/

	auto Kept = [ & ]( long lBase, int iPower, long lPlaces, int iMag, long* plRound ) -> std::string
	{
		if ( iMag >= 0 )
		{
			*plRound = lPlaces;

			return std::string( 1, Digit( lBase, false ) ) + Digits( lBase, iMag ) + "." + Digits( lBase, lPlaces );
		}

		long	l_noughts = -iMag;

		*plRound = (iPower == 2) ? lPlaces : (lPlaces + l_noughts);

		return "0." + std::string( static_cast<size_t>( l_noughts ), '0' ) + Digit( lBase, false ) + Digits( lBase, *plRound - l_noughts - 1 );
	};

	// The root of x to iPower at lPlaces, against what it should be...
	auto RootIs = [ & ]( Unary fn, CLZNumber* pX, long lBase, long lPlaces, CLZNumber* pExpect ) -> bool
	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = lPlaces;

		CLZNumber	answer;

		answer.SetBase( lBase, true );
		fn( pX, &answer );
		glWorkingPrecision = l_saved;

		return lz_compare( &answer, pExpect ) == kEqualTo;
	};

	/*
		**A hair either side of a half.** r is y plus the half at y's
		last place, plus or less a unit twenty places past it, deep
		inside any guard; x is r squared or cubed, worked exactly. The
		root of x rounds to y plus a unit above the half and to y below
		it. An odd base has no finite half, so its half is the middle
		digit twenty-four times over, short of it by far less than the
		unit. Answers a dozen digits long, a few digits, and a dozen
		noughts after the point: the cube root decides against its
		operand scaled back from the band, at places that grow with the
		noughts.
	*/

	for ( int iPower = 2; iPower <= 3; iPower++ )
	{
		Unary	fn = (iPower == 2) ? lz_root_2 : lz_root_3;

		for ( long lBase : { 10L, 16L, 2L, 7L } )
		{
			long	l_wrong = 0;
			long	l_cases = 0;

			for ( long lPlaces : { 36L, 100L } )
			{
				for ( int iMag : { 12, 0, -12 } )
				{
					for ( int iSide = 0; iSide < 2; iSide++ )
					{
						for ( int k = 0; k < 2; k++ )
						{
							long	l_saved = glWorkingPrecision;
							long	l_round = 0;
							std::string	s_y = Kept( lBase, iPower, lPlaces, iMag, &l_round );

							glWorkingPrecision = (iPower * (l_round + 40 + 16)) + 8;

							CLZNumber	y, half, hair, r, x, unit;

							y.SetBase( lBase, true );
							y.Assign( s_y.c_str() );
							half.SetBase( lBase, true );
							half.Assign( Past( l_round, std::string( (lBase % 2 == 0) ? 1 : 24, kDigitChars[ lBase / 2 ] ) ).c_str() );
							hair.SetBase( lBase, true );
							hair.Assign( Past( l_round + 19, "1" ).c_str() );
							unit.SetBase( lBase, true );
							unit.Assign( Past( l_round - 1, "1" ).c_str() );

							lz_add( &y, &half, &r );

							if ( iSide == 0 )
								lz_add( &r, &hair );
							else
								lz_subtract( &r, &hair );

							x.Assign( &r );

							for ( int iTimes = 1; iTimes < iPower; iTimes++ )
								lz_multiply( &x, &r );

							if ( iSide == 0 )
								lz_add( &y, &unit );

							glWorkingPrecision = l_saved;

							l_cases++;

							if ( !RootIs( fn, &x, lBase, lPlaces, &y ) )
								l_wrong++;
						}
					}
				}
			}

			char	szLabel[ 72 ];
			char	szDetail[ 72 ];

			::snprintf( szLabel, sizeof( szLabel ), "%s, base %ld, a hair from a half", (iPower == 2) ? "sqrt" : "cbrt", lBase );
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld rounded the wrong way", l_wrong, l_cases );
			CheckTrue( szLabel, l_wrong == 0, szDetail );
		}
	}

	/*
		**Exactly on a half**, in an even base: r is y plus the half at
		its last place and nothing else, so the root is r exactly, and
		goes to the even kept digit as RoundAt's half-even does.
	*/

	for ( int iPower = 2; iPower <= 3; iPower++ )
	{
		Unary	fn = (iPower == 2) ? lz_root_2 : lz_root_3;

		for ( long lBase : { 10L, 16L, 2L } )
		{
			long	l_wrong = 0;
			long	l_cases = 0;

			for ( long lPlaces : { 1L, 36L, 100L } )
			{
				for ( int iMag : { 12, 0, -12 } )
				{
					if ( (iMag < 0) && (lPlaces < 36) )
						continue;

					for ( int k = 0; k < 2; k++ )
					{
						long	l_saved = glWorkingPrecision;
						long	l_round = 0;
						std::string	s_y = Kept( lBase, iPower, lPlaces, iMag, &l_round );

						glWorkingPrecision = (iPower * (l_round + 16)) + 8;

						CLZNumber	y, half, r, x, unit;

						y.SetBase( lBase, true );
						y.Assign( s_y.c_str() );
						half.SetBase( lBase, true );
						half.Assign( Past( l_round, std::string( 1, kDigitChars[ lBase / 2 ] ) ).c_str() );
						unit.SetBase( lBase, true );
						unit.Assign( Past( l_round - 1, "1" ).c_str() );

						lz_add( &y, &half, &r );
						x.Assign( &r );

						for ( int iTimes = 1; iTimes < iPower; iTimes++ )
							lz_multiply( &x, &r );

						// The kept digit odd goes up to the even one...
						if ( std::string( kDigitChars ).find( s_y.back() ) % 2 )
							lz_add( &y, &unit );

						glWorkingPrecision = l_saved;

						l_cases++;

						if ( !RootIs( fn, &x, lBase, lPlaces, &y ) )
							l_wrong++;
					}
				}
			}

			char	szLabel[ 72 ];
			char	szDetail[ 72 ];

			::snprintf( szLabel, sizeof( szLabel ), "%s, base %ld, exactly on a half", (iPower == 2) ? "sqrt" : "cbrt", lBase );
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld not to the even digit", l_wrong, l_cases );
			CheckTrue( szLabel, l_wrong == 0, szDetail );
		}
	}
}

/////////////////////////////////////////////////////////////////
//	Halving in an odd base, and a long operand's tail. Until 2026-10-03
//	lz_exp_e halved its fraction in place at the operand's precision,
//	where the width it worked at was wider: in an odd base the half
//	never terminates and the cut fell inside the guard, and an operand
//	written past the width lost the places past it. And lz_agm halved by
//	the reciprocal of two, which in an odd base is not a half. Each pin
//	is an answer HEAD's engine got wrong, found by the near-half census's
//	random operands; the expected strings are the exact values rounded to
//	36 places, computed in Python.
/////////////////////////////////////////////////////////////////

struct SOddBasePin
{
	const char*	m_szFunction;
	long		m_lBase;
	const char*	m_szOperand;
	const char*	m_szSecond;
	const char*	m_szExpect;
};

// e^x of an operand written to the width, in bases nine, eleven and thirteen:
// 14 answers HEAD's engine got wrong, each the exact value rounded to 36 places
static const SOddBasePin kOddBaseExp[] = {
	{ "exp_e", 9, "0.546440558416488311740573867005854423", nullptr,
	  "1.755754834240875251355088040500757586" },
	{ "exp_e", 9, "1.743657620444628312877843877851776581", nullptr,
	  "6.220824156266011044464872226263860576" },
	{ "exp_e", 9, "1.761603842166313651263427356266645827", nullptr,
	  "6.342510164806047516131858884678577753" },
	{ "exp_e", 9, "1.750557852146253134782710177867251870", nullptr,
	  "6.261770778058637520308237120437611151" },
	{ "exp_e", 9, "1.725317856434028214672218455383601268", nullptr,
	  "6.087627083875414032234852274517658053" },
	{ "exp_e", 9, "1.714436710771305064446846586015475376", nullptr,
	  "6.021744680300140154032607862765206051" },
	{ "exp_e", 9, "1.762278186643614304057416452302077323", nullptr,
	  "6.346632544643420205167326248612124884" },
	{ "exp_e", 9, "0.753320034668036502540116333362336137", nullptr,
	  "2.283527417420703356784455327314730573" },
	{ "exp_e", 9, "1.650247050227306260848543158421536511", nullptr,
	  "5.563021521007251855273644627053481574" },
	{ "exp_e", 11, "0.A2604733A986253376942506474A67358001", nullptr,
	  "2.597A62206246817927740224624836101172" },
	{ "exp_e", 11, "1.700620A371AA806762893408942A10047573", nullptr,
	  "5.15857451219577522585190871A85777A944" },
	{ "exp_e", 11, "1.953485767561262629A159323A99785450A0", nullptr,
	  "6.4899172173A4759545251303070706463220" },
	{ "exp_e", 13, "1.BA024B20202496A0A8C64B6006C12B504985", nullptr,
	  "6.95075686C05C609CCC8B3C7A179C7612065B" },
	{ "exp_e", 13, "0.A45CC01A3A5201C20657064A20A7097CC661", nullptr,
	  "2.2A6249733A79C5061117394047411AA866C7" },
};

// e^x of an operand written two to eight places past the width, in base ten:
// 18 answers HEAD's engine got wrong, each the exact value rounded to 36 places
static const SOddBasePin kLongOperandExp[] = {
	{ "exp_e", 10, "1.86816440780570463607245935666565773761", nullptr,
	  "6.476397457205531153263701481687429976" },
	{ "exp_e", 10, "1.60378831822853410766070379399879304663", nullptr,
	  "4.971831673701500035546696173585692697" },
	{ "exp_e", 10, "1.89596813307927867318161462229153060019", nullptr,
	  "6.658992075467048268560001121688361158" },
	{ "exp_e", 10, "1.76776495246897306785460977941864529103", nullptr,
	  "5.857746377894023670412838598345294958" },
	{ "exp_e", 10, "1.504201053282733404577512622187181093277", nullptr,
	  "4.500556488685081228530654176122124113" },
	{ "exp_e", 10, "1.929588855210285216923573461546589294309", nullptr,
	  "6.886678237563521966627452382222675010" },
	{ "exp_e", 10, "1.966236984262654037711222448759574188124", nullptr,
	  "7.143743825671475462357993505077699916" },
	{ "exp_e", 10, "1.523407783570661966313192144440434354325", nullptr,
	  "4.587832927826674248443585201012845209" },
	{ "exp_e", 10, "0.9457311985292778310946245779061909136323", nullptr,
	  "2.574695303893516068629237035616567219" },
	{ "exp_e", 10, "1.8047694318133389545801185496315358559572", nullptr,
	  "6.078569762126857138308448415802535631" },
	{ "exp_e", 10, "1.5227832922547179187037176738559349752383", nullptr,
	  "4.584968760421163322565445369334292545" },
	{ "exp_e", 10, "2.9882779344723675626252284652221860159717", nullptr,
	  "19.851467511926445237756335166095434036" },
	{ "exp_e", 10, "0.79657862630595869982065270143189471356598755", nullptr,
	  "2.217939532328180989748931934758014292" },
	{ "exp_e", 10, "1.59425169017895915648753221009847286701315749", nullptr,
	  "4.924642534581288212196350450419134034" },
	{ "exp_e", 10, "1.6142666355137517607291474129636579005825078", nullptr,
	  "5.024202000752448590207808012972312902" },
	{ "exp_e", 10, "1.83546171685036426522730280224822051127644514", nullptr,
	  "6.268027531351569020514000781364105215" },
	{ "exp_e", 10, "1.63909218685158861514929078003574830270967758", nullptr,
	  "5.150491705177818900911981892674881514" },
	{ "exp_e", 10, "2.9198500404941816544545502435728517615134983", nullptr,
	  "18.538507225908626628163119417204000938" },
};

// The six hyperbolic functions in bases five, seven, nine and eleven:
// 39 answers HEAD's engine got wrong, each the exact value rounded to 36 places
static const SOddBasePin kOddBaseHyperbolic[] = {
	{ "hypercosine", 7, "1.542045522211141511110223104056433031", nullptr,
	  "3.053630463030101350102010452131312552" },
	{ "hypercosine", 7, "1.501366651606053664520261355552155022", nullptr,
	  "2.610256321226355361331562065266264023" },
	{ "hypercosine", 7, "1.636506145130110346025454502164321345", nullptr,
	  "3.354505430525043362516052324022634242" },
	{ "hypercosine", 7, "1.434542461126541561541215016340053403", nullptr,
	  "2.456031115444666253630624163042055023" },
	{ "hypercosine", 7, "0.501514440636201042622206561415500363", nullptr,
	  "1.161435032015504603065222036340041306" },
	{ "hypersine", 7, "1.431511052326426202054164141540352412", nullptr,
	  "2.321516535310223006424054132125054322" },
	{ "hypersine", 7, "-1.511045660550030424012311010066002241", nullptr,
	  "-2.516610036340226065665641340351003231" },
	{ "hypersine", 7, "1.634345566126451263232424632341465233", nullptr,
	  "3.243356255003330015412000560515310463" },
	{ "hypersine", 7, "0.463406453046362612132416621644502364", nullptr,
	  "0.523032221241531114353455124003123361" },
	{ "hypersine", 7, "-1.345312662144013654556402661350152523", nullptr,
	  "-2.122444412016621141246110055062450410" },
	{ "hypertangent", 7, "0.501343522501165532416225015404253664", nullptr,
	  "0.421220312113612541052451643553441213" },
	{ "hypersecant", 7, "0.635526001436410252440364263516301550", nullptr,
	  "0.452213143133130465556010033634004512" },
	{ "hypersecant", 7, "0.645231662323244510553632531062540400", nullptr,
	  "0.446016051524413341150314666361600426" },
	{ "hypercosecant", 7, "0.366136461503366224104200335146366550", nullptr,
	  "1.444201050623334163305046053052651545" },
	{ "hypercosecant", 7, "0.340652441433503326605131516445535210", nullptr,
	  "1.603064015652155103006533442232545314" },
	{ "hypercosecant", 7, "-0.406425550241350112624653112612256106", nullptr,
	  "-1.412111150663305540664160603115224355" },
	{ "hypercotangent", 7, "-0.344020564445124535054346415303125031", nullptr,
	  "-2.041502535241212534225415245063231300" },
	{ "hypercotangent", 7, "-0.422064230220063446251262632305150040", nullptr,
	  "-1.550365166633655253400645061404263224" },
	{ "hypercotangent", 7, "0.603402112532540506413604241211200162", nullptr,
	  "1.266621366331162332023406311224031654" },
	{ "hypercotangent", 7, "0.446542610456561253040421013005403265", nullptr,
	  "1.463304653006664503542206336443512622" },
	{ "hypercotangent", 7, "0.353464653150153351630532640442162525", nullptr,
	  "2.011145356654553002533513412505034212" },
	{ "hypercosine", 9, "2.761456452625367182656738533516020865", nullptr,
	  "8.632181748644412000226320465707364374" },
	{ "hypercosine", 9, "2.560407417727041664105467780538656416", nullptr,
	  "6.868474025223602502076852228873384811" },
	{ "hypercosine", 9, "2.678574320877371387762638235038240225", nullptr,
	  "7.865512833341385046778041685552885667" },
	{ "hypercosine", 11, "2.66091A600861366A3017534313A56A405906", nullptr,
	  "6.81603934324152AA8A2A7825A66271233A04" },
	{ "hypercosine", 11, "2.670354453023435394A563439A0331770916", nullptr,
	  "6.87A204041211025554A65A4427757A089258" },
	{ "hypercosine", 11, "2.947A3726401357A0AA1A1248201015045825", nullptr,
	  "8.80A47630913035934A436114027947107807" },
	{ "hypersine", 9, "2.731587430220376278205816331574240072", nullptr,
	  "8.300578214617414077815134470382812334" },
	{ "hypersine", 9, "-2.613841040282364517324665446332371366", nullptr,
	  "-7.254822564046755354764155447640528600" },
	{ "hypersine", 9, "-2.750118751082053623770577178738701378", nullptr,
	  "-8.465423225060608281538261318352418151" },
	{ "hypersine", 11, "2.A5A885691472765887702A2745630080A764", nullptr,
	  "9.677385191A6A15697A04422516A623227A83" },
	{ "hypersine", 11, "2.A8625792708270779A546382777925A09542", nullptr,
	  "9.8A9A4981A9A292A39A865152799A37A02410" },
	{ "hypersine", 11, "2.60A99A65766A531451216395A89724433228", nullptr,
	  "6.42A51A2101639A132A293A270610867A5418" },
	{ "hypercosecant", 5, "-0.224324011400111042400124144330120141", nullptr,
	  "-1.410444202343004112104030120412404143" },
	{ "hypercosecant", 5, "0.331043342111321112031323333002241234", nullptr,
	  "1.112010003042304224201203010101032114" },
	{ "hypercosecant", 5, "0.311123130242111322322030120040233003", nullptr,
	  "1.204112414302130221311042233203330113" },
	{ "hypercotangent", 5, "0.223001020222042013400233334020131233", nullptr,
	  "2.033304312030300402100241121014002142" },
	{ "hypercotangent", 5, "0.311420431004123123002443232132104403", nullptr,
	  "1.332134230331433220211403202344231240" },
	{ "hypercotangent", 5, "0.312341101241403104312212131011243243", nullptr,
	  "1.330244013002213104423203413443110120" },
};

// The AGM in base seven: fourteen pairs with the second operand among 1, 17/7,
// 5/7, 4 and 3, and eight near 7^-30, 7^-8, 7^8 and 7^30.
// 22 answers HEAD's engine got wrong, each the exact value rounded to 36 places
static const SOddBasePin kOddBaseAgm[] = {
	{ "agm", 7, "3.536234062653565153123223443053422250", "2.300000000000000000000000000000000000",
	  "3.034063636206464124423460430106300044" },
	{ "agm", 7, "5.250055633560004526431023000251630540", "1.000000000000000000000000000000000000",
	  "2.511621066602615440314312010266530112" },
	{ "agm", 7, "0.653232052510135150162253066112642505", "0.500000000000000000000000000000000000",
	  "0.560004202554051232216413334101134604" },
	{ "agm", 7, "6.230110424563244416564346535461263320", "1.000000000000000000000000000000000000",
	  "3.032536105613066061562126004132646035" },
	{ "agm", 7, "5.014126050414423534423223064024612554", "0.500000000000000000000000000000000000",
	  "2.234154054320302100516565664615015510" },
	{ "agm", 7, "5.215320450302600362665610056253025201", "2.300000000000000000000000000000000000",
	  "3.506506445156263304330142352043036243" },
	{ "agm", 7, "6.046463300241324654206520616424200466", "4.000000000000000000000000000000000000",
	  "4.665232013053232506604350462203531066" },
	{ "agm", 7, "5.502126552215354560126153313356225364", "2.300000000000000000000000000000000000",
	  "3.620243562542000050530162636605366234" },
	{ "agm", 7, "10.122163203455421426524556435626403401", "3.000000000000000000000000000000000000",
	  "4.603306235242216500602510143246521050" },
	{ "agm", 7, "5.223365426535462014200510530514335022", "0.500000000000000000000000000000000000",
	  "2.313446252343321224013422133413400230" },
	{ "agm", 7, "5.464532301200316443334353412106230232", "2.300000000000000000000000000000000000",
	  "3.615363141366200014612344204306625244" },
	{ "agm", 7, "3.466301563545154256544504153622124021", "4.000000000000000000000000000000000000",
	  "3.566154253503164544650424003161624305" },
	{ "agm", 7, "10.202312221536331005634411336355363145", "1.000000000000000000000000000000000000",
	  "3.245656016036343015626545115243342525" },
	{ "agm", 7, "1.546204235612652406602310155513005551", "3.000000000000000000000000000000000000",
	  "2.240611523523222203405232416025106216" },
	{ "agm", 7, "0.00000000000000000000000000000311360660232150224532232231160333303600", "4261030136005345331653032665615.62066500000000000000000000000000000000",
	  "26030114341612036216362143500.012035056216411541314063236410141101" },
	{ "agm", 7, "0.00000000000000000000000000000313403621514662616031065220252115314200", "4143140021413413060015354056365.06414600000000000000000000000000000000",
	  "25156336516052235615101501431.533336025530063215603316200433003634" },
	{ "agm", 7, "0.0000000211065034605502462202100056334462263200", "205226660.25412441010121031600256261440000000000",
	  "4665404.135423523606363223432242164650623365" },
	{ "agm", 7, "0.0000000634102302621500031551211243622636546600", "301151023.20356323136534050134214033340000000000",
	  "10222243.315410566663350521065454611042512130" },
	{ "agm", 7, "644635420.30442265656256100660612032330000000000", "3.4514565622541566210314451445630330160000000000",
	  "41156261.425016163643215312656154513263116656" },
	{ "agm", 7, "362426360.25310423255100344511466013110000000000", "0.0000000326031015303603251520061525114354655200",
	  "12134323.144205232656452022514526566403545220" },
	{ "agm", 7, "3125304410436452125266206604513.55501400000000000000000000000000000000", "5.56645553434162422330616230015125041300000000000000000000000000000000",
	  "41065540452546345431031526633.015610156112035046602322035654443211" },
	{ "agm", 7, "3122002105342220244660432144265.23060500000000000000000000000000000000", "4401540216033210.42235341541115560003400000000000000000000000000000000000000000000000",
	  "110606450562615550335501053501.435300300145204105154155146512434211" },
};

static void TestHalvedInOddBases()
{
	Banner( "Halving in an odd base, and a long operand's tail" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	struct { const char* m_szName; Unary m_fn; } kUnary[] =
	{
		{ "exp_e", lz_exp_e }, { "hypersine", lz_hypersine }, { "hypercosine", lz_hypercosine },
		{ "hypertangent", lz_hypertangent }, { "hypersecant", lz_hypersecant },
		{ "hypercosecant", lz_hypercosecant }, { "hypercotangent", lz_hypercotangent }
	};

	struct { const char* m_szLabel; const SOddBasePin* m_pPins; size_t m_nPins; } kGroups[] =
	{
		{ "e^x at the width, bases 9, 11 and 13", kOddBaseExp, sizeof( kOddBaseExp ) / sizeof( kOddBaseExp[ 0 ] ) },
		{ "e^x of an operand past the width", kLongOperandExp, sizeof( kLongOperandExp ) / sizeof( kLongOperandExp[ 0 ] ) },
		{ "hyperbolics, bases 5, 7, 9 and 11", kOddBaseHyperbolic, sizeof( kOddBaseHyperbolic ) / sizeof( kOddBaseHyperbolic[ 0 ] ) },
		{ "the AGM in base seven", kOddBaseAgm, sizeof( kOddBaseAgm ) / sizeof( kOddBaseAgm[ 0 ] ) }
	};

	long	l_saved = glWorkingPrecision;

	for ( auto& group : kGroups )
	{
		long	l_wrong = 0;

		for ( size_t i = 0; i < group.m_nPins; i++ )
		{
			const SOddBasePin&	pin = group.m_pPins[ i ];

			glWorkingPrecision = 36;

			CLZNumber	x, y, answer, expect;

			x.SetBase( pin.m_lBase, true );
			y.SetBase( pin.m_lBase, true );
			answer.SetBase( pin.m_lBase, true );
			expect.SetBase( pin.m_lBase, true );
			x.Assign( pin.m_szOperand );
			expect.Assign( pin.m_szExpect );

			long	l_error = kFailed;

			if ( pin.m_szSecond != nullptr )
			{
				y.Assign( pin.m_szSecond );
				l_error = lz_agm( &x, &y, &answer );
			}
			else
			{
				for ( auto& u : kUnary )
					if ( ::strcmp( u.m_szName, pin.m_szFunction ) == 0 )
						l_error = u.m_fn( &x, &answer );
			}

			if ( (l_error != kNoError) || (lz_compare( &answer, &expect ) != kEqualTo) )
				l_wrong++;
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrong, group.m_nPins );
		CheckTrue( group.m_szLabel, l_wrong == 0, szDetail );
	}

	glWorkingPrecision = l_saved;
}

/////////////////////////////////////////////////////////////////
//	The reciprocal, rounded once. Until 2026-10-03 lz_reciprocal divided
//	kQuotientGuardDigits wider than its answer, which the division
//	answered a place wider still, and rounded again to the answer's
//	width: an answer within the extra places' last unit of a half could
//	round the wrong way. It rounds the raw quotient once at the answer's
//	width now. Expected strings are the exact rationals rounded, computed
//	in Python.
/////////////////////////////////////////////////////////////////

struct SReciprocalPin
{
	long		m_lBase;
	long		m_lPlaces;
	const char*	m_szOperand;
	const char*	m_szExpect;
};

// Answers within a millionth of a unit of a half, at 36 places: 1/x for x
// steered there by the near-half census. 10 HEAD's engine got wrong.
static const SReciprocalPin kReciprocalNearHalf[] = {
	{ 10, 36, "-0.100159586509810913888132106296152998",
	  "-9.984066776294520551005007803875898435" },
	{ 10, 36, "0.318107223381193664742978953180002620",
	  "3.143594129585928431680486747901939247" },
	{ 10, 36, "0.146601435145174283231846548909000878",
	  "6.821215624592862347697065447069953919" },
	{ 10, 36, "0.727978175759461379610189720593053224",
	  "1.373667553916369188371505209774034365" },
	{ 16, 36, "-0.7A608A80339D21C47C700EA55A10E2ECD61B",
	  "-2.17866474E037000000001AC1C645F8B43783" },
	{ 2, 36, "1.011111111110001111001011010111001011",
	  "0.101010101011011100110100110000101001" },
	{ 2, 36, "0.110100000001011000110010110101110010",
	  "1.001110101111001000010100011111011111" },
	{ 2, 36, "-1.001101110100101011011110001001100110",
	  "-0.110100101000011101011001110111111111" },
	{ 2, 36, "1.110101000001010011111010000000111100",
	  "0.100011000000001001111010000001010111" },
	{ 2, 36, "-0.100011000001010010000111100111000011",
	  "-1.110100111101100010100110110110000101" },
};

// Small answers, a tenth to a hundred-millionth, rounded at 36 places plus
// their noughts and built by inverting a target a hair past a half there.
// 8 HEAD's engine got wrong.
static const SReciprocalPin kReciprocalSmall[] = {
	{ 10, 36, "27436.535579999694143103187710080350554866",
	  "0.0000364477503759245083156313632963228519" },
	{ 10, 36, "-339812.351887722235719750491710039879728846",
	  "-0.00000294280062053309667329569882548963211" },
	{ 10, 36, "1279687.816846945176384863949831213545143460",
	  "0.000000781440587958338944908210980729528781" },
	{ 10, 36, "-13344801.304137612169849574870996087527569873",
	  "-0.0000000749355480991647109972741156088207565" },
	{ 16, 36, "-A3C7D.123284A9E4837454E8E816C3F19B5B7D1FD0",
	  "-0.0000190252DDD82D5307872AFA0EE1B01C1D8CBD" },
	{ 16, 36, "49B48A.BDD778CF28398557BD7604A3C9A9A7ED32D7",
	  "0.0000037929D161B36C3F5E55DC45686AA6CEA0787" },
	{ 16, 36, "15A4A72.9439E78684AFDA90C76CAF3EA1376C41297B",
	  "0.000000BD4057D2691C9ACCDC6898CF41295FE2CB83" },
	{ 16, 36, "-1DD6287E.C09E12E708BD2800D58FAD3773EE5DF8EA1F",
	  "-0.0000000894800D499C2CD9F905F429111E9F12524F1" },
};

// Exact ties, 1/2^k in base ten at a width one place short of its digits:
// each goes to its even kept digit. HEAD's engine had these right.
static const SReciprocalPin kReciprocalTie[] = {
	{ 10, 1, "4",
	  "0.2" },
	{ 10, 2, "8",
	  "0.12" },
	{ 10, 2, "16",
	  "0.062" },
	{ 10, 5, "256",
	  "0.0039062" },
	{ 10, 20, "536870912",
	  "0.0000000018626451492309570312" },
	{ 10, 20, "1073741824",
	  "0.00000000093132257461547851562" },
	{ 10, 36, "4503599627370496",
	  "0.000000000000000222044604925031308084726333618164062" },
	{ 10, 77, "2596148429267413814265248164610048",
	  "0.00000000000000000000000000000000038518598887744717061119558851698546370762032964307763904798775911331176757812" },
};

static void TestReciprocalRoundedOnce()
{
	Banner( "The reciprocal, rounded once" );

	struct { const char* m_szLabel; const SReciprocalPin* m_pPins; size_t m_nPins; } kGroups[] =
	{
		{ "1/x a hair from a half", kReciprocalNearHalf, sizeof( kReciprocalNearHalf ) / sizeof( kReciprocalNearHalf[ 0 ] ) },
		{ "small 1/x a hair from a half", kReciprocalSmall, sizeof( kReciprocalSmall ) / sizeof( kReciprocalSmall[ 0 ] ) },
		{ "1/2^k exactly on a half", kReciprocalTie, sizeof( kReciprocalTie ) / sizeof( kReciprocalTie[ 0 ] ) }
	};

	long	l_saved = glWorkingPrecision;

	for ( auto& group : kGroups )
	{
		long	l_wrong = 0;

		for ( size_t i = 0; i < group.m_nPins; i++ )
		{
			const SReciprocalPin&	pin = group.m_pPins[ i ];

			glWorkingPrecision = pin.m_lPlaces;

			CLZNumber	x, answer, expect;

			x.SetBase( pin.m_lBase, true );
			answer.SetBase( pin.m_lBase, true );
			expect.SetBase( pin.m_lBase, true );
			x.Assign( pin.m_szOperand );
			expect.Assign( pin.m_szExpect );

			if ( (lz_reciprocal( &x, &answer ) != kNoError) || (lz_compare( &answer, &expect ) != kEqualTo) )
				l_wrong++;
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrong, group.m_nPins );
		CheckTrue( group.m_szLabel, l_wrong == 0, szDetail );
	}

	/*
		**A quotient a unit short of a power of the base rounds up to it**,
		and is then written to the width the rule gives the power: 36
		places for 0.1, not the 37 the quotient's own lead asked for. HEAD
		wrote 37.
	*/

	struct { long m_lBase; const char* m_szOperand; const char* m_szExpect; } kCarry[] =
	{
		{ 10, "10.000000000000000000000000000000000001", "0.100000000000000000000000000000000000" },
		{ 16, "10.000000000000000000000000000000000001", "0.100000000000000000000000000000000000" }
	};

	for ( auto& carry : kCarry )
	{
		glWorkingPrecision = 36;

		CLZNumber	x, answer;

		x.SetBase( carry.m_lBase, true );
		answer.SetBase( carry.m_lBase, true );
		x.Assign( carry.m_szOperand );
		lz_reciprocal( &x, &answer );

		char	szLabel[ 72 ];

		::snprintf( szLabel, sizeof( szLabel ), "1/(b + a unit), base %ld", carry.m_lBase );
		CheckRendering( szLabel, &answer, carry.m_szExpect );
	}

	glWorkingPrecision = l_saved;
}

/////////////////////////////////////////////////////////////////
//	A carry to a power of the base, written to the width the rule gives
//	the power. Every function ending in the rule for a quotient, and
//	division with its own, read the width off the answer before rounding
//	it, so 0.0999...9 rounding to 0.1 came out a place wide -- in every
//	base, never in value -- until 2026-10-10.
/////////////////////////////////////////////////////////////////

struct SCarryPin
{
	const char*	m_szFunction;
	long		m_lBase;
	const char*	m_szOperand;
	const char*	m_szSecond;
	const char*	m_szExpect;
};

// Operands written past the width so the answer can be placed just under
// the power; division's are whole numbers whose quotient is exact. 12,
// all of them a place wide at c87187a.
static const SCarryPin kCarryPins[] = {
	{ "tangent", 10, "0.099668652491162027378446119878020590218525847", nullptr,
	  "0.100000000000000000000000000000000000" },
	{ "hypercosecant", 10, "5.2983423656105887573688256891129063023923710361", nullptr,
	  "0.0100000000000000000000000000000000000" },
	{ "root_3", 10, "0.000999999999999999999999999999999999999250000", nullptr,
	  "0.100000000000000000000000000000000000" },
	{ "divide", 10, "99999999999999999999999999999999999999999999999", "1000000000000000000000000000000000000000000000000",
	  "0.10000000000000000000000000000000000000" },
	{ "divide", 10, "999999999999999999999999999999999999999999999999", "1000000000000000000000000000000000000000000000000",
	  "1.0000000000000000000000000000000000000" },
	{ "cotangent", 16, "1.911FB59997F3A6B03DDB9A1B964AEA243B88839628F3D0", nullptr,
	  "0.0100000000000000000000000000000000000" },
	{ "hypersecant", 16, "3.76FA5FC9BAF1B964F7971126EF56F087DECC11DCA43CB", nullptr,
	  "0.100000000000000000000000000000000000" },
	{ "divide", 16, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", "1000000000000000000000000000000000000000000000000",
	  "0.00010000000000000000000000000000000000000" },
	{ "hypertangent", 2, "0.000100000000010101011000100010101101001100110101", nullptr,
	  "0.000100000000000000000000000000000000000" },
	{ "hypersine", 7, "0.066556231515340563410561304401243410603626455", nullptr,
	  "0.100000000000000000000000000000000000" },
	{ "root_3", 7, "0.0000006666666666666666666666666666666666661515", nullptr,
	  "0.0100000000000000000000000000000000000" },
	{ "divide", 7, "6666666666666666666666666666666666666666666666", "1000000000000000000000000000000000000000000000000",
	  "0.010000000000000000000000000000000000000" },
};

static void TestCarrySpelled()
{
	Banner( "A carry to a power of the base, spelled" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	struct { const char* m_szName; Unary m_fn; } kUnary[] =
	{
		{ "tangent", lz_tangent }, { "cotangent", lz_cotangent }, { "hypersine", lz_hypersine },
		{ "hypertangent", lz_hypertangent }, { "hypercosecant", lz_hypercosecant },
		{ "hypersecant", lz_hypersecant }, { "root_3", lz_root_3 }
	};

	long		l_saved = glWorkingPrecision;
	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kCarryPins )
	{
		glWorkingPrecision = 36;

		CLZNumber	x, y, answer;

		x.SetBase( pin.m_lBase, true );
		y.SetBase( pin.m_lBase, true );
		answer.SetBase( pin.m_lBase, true );
		x.Assign( pin.m_szOperand );

		if ( pin.m_szSecond != nullptr )
		{
			y.Assign( pin.m_szSecond );
			lz_divide( &x, &y, &answer );
		}
		else
		{
			for ( auto& u : kUnary )
				if ( ::strcmp( u.m_szName, pin.m_szFunction ) == 0 )
					u.m_fn( &x, &answer );
		}

		char*	szGot = answer.GetStringRep();

		if ( (szGot == nullptr) || (::strcmp( szGot, pin.m_szExpect ) != 0) )
		{
			if ( l_wrong++ == 0 )
				s_first = std::string( pin.m_szFunction ) + " base " + std::to_string( pin.m_lBase ) + ": " + (szGot ? szGot : "(null)");
		}

		::free( szGot );
	}

	glWorkingPrecision = l_saved;

	char	szDetail[ 160 ];
	size_t	n_pins = sizeof( kCarryPins ) / sizeof( kCarryPins[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu at the rule's width", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "a power of the base, carried to", l_wrong == 0, szDetail );
}

/////////////////////////////////////////////////////////////////
//	The circular quotients' guard. Until 2026-10-03 lz_tangent,
//	lz_cotangent, lz_secant and lz_cosecant worked four places wider than
//	their answers, the sine's guard, and a quotient spends that twice: on
//	its divisor's smallness and on its own size. They were blind within
//	about 10^-3 of a unit of a half, where the sine is blind within
//	10^-5. They work six places wider now -- see kTrigQuotientGuardDigits.
//	Expected strings are the exact values rounded to 36 places, computed
//	in Python.
/////////////////////////////////////////////////////////////////

struct SQuotientPin
{
	const char*	m_szFunction;
	long		m_lBase;
	const char*	m_szOperand;
	const char*	m_szExpect;
};

// Random operands written to the width that HEAD's engine misrounded:
// 11 of 30,000 tangents, cotangents and secants in bases ten and sixteen.
static const SQuotientPin kQuotientRandom[] = {
	{ "tangent", 10, "1.427657799879782640851545897080162359",
	  "6.938461144334330467761261486087251939" },
	{ "tangent", 10, "-1.456044806716218272530486406686178909",
	  "-8.676197382783424446533919247462178127" },
	{ "tangent", 10, "-1.352895943401395031076262151834012933",
	  "-4.516388636546414785254775257205981152" },
	{ "cotangent", 10, "3.000563381142706457808686171041098377",
	  "-7.043654295744239734256777253666695374" },
	{ "cotangent", 10, "2.963854038518862447405242714748264642",
	  "-5.566868047407817420648520184111449174" },
	{ "cotangent", 16, "3.128907F66C7B6E06DB7322100A4E1EC277E2",
	  "-E.6E17912E1C42F6D9B500E91974A21AF61CE7" },
	{ "cotangent", 16, "3.124EE5AD81936A2B7943D2CABF40C1FB4AE4",
	  "-E.3F2E24BE6E59E557EE7AA08A5BE61E690FF5" },
	{ "secant", 10, "1.452711817583573813979973346168089462",
	  "8.488224110036036233684672923518407554" },
	{ "secant", 10, "1.677823504770557943982623388449061697",
	  "-9.361282895931103621461945707591336803" },
	{ "secant", 10, "1.439208814394253662618491968255942184",
	  "7.621480952694635638709375921325451379" },
	{ "secant", 16, "1.A935059A852E706BC2219468D246E0A7EB11",
	  "-B.1AF705C84F009054E4FD275BE37012DD0C24" },
};

// Answers a hair from a half, found by the near-half census: 47 HEAD's
// engine got wrong.
static const SQuotientPin kQuotientNearHalf[] = {
	{ "tangent", 10, "1.358168897229432050022071765604072290",
	  "4.631971719028885381419513502030876823" },
	{ "tangent", 10, "1.358168897229432050022071765604453295",
	  "4.631971719028885381419513502039432353" },
	{ "tangent", 10, "1.355671428881279487559029960918098857",
	  "4.576531886622793088292837332182274974" },
	{ "tangent", 16, "-1.7248305663493A83D951CBE7BDCFB7381272",
	  "-7.FF8BA76258C4BFFFFFFF20017E8C29A57C22" },
	{ "tangent", 16, "-1.7248305663493A83D951CBE7BDCFB737F77D",
	  "-7.FF8BA76258C4BFFFFFFF20017E8C299EA4B0" },
	{ "tangent", 16, "-1.6381A6EE3320147266A47F891D0A31170123",
	  "-5.6E40D3497D4A7FFFFFD4023F3355F6A5B0C4" },
	{ "tangent", 2, "0.110001001101011011001101100110010000",
	  "0.111101111011000011001011011010100110" },
	{ "tangent", 2, "-1.000001110011100011000100101101100011",
	  "-1.101010001001010000011111110100001100" },
	{ "tangent", 2, "0.111010101101011101100101010010011101",
	  "1.010011100101101001011001110010101110" },
	{ "tangent", 7, "0.664420660433515530226113565601051423",
	  "1.351052204105234322461430636563630035" },
	{ "tangent", 7, "1.161440032134123606524113654662345154",
	  "3.140125214044060114520414630105562123" },
	{ "tangent", 7, "-1.265251122510016145625350436620331631",
	  "-6.520600356564513605243435362002615611" },
	{ "secant", 10, "1.223412110222299683674445456791138234",
	  "2.937379979750214253897411331394565983" },
	{ "secant", 10, "0.629902098698420008276942183826375774",
	  "1.237493301056535655035872163915365350" },
	{ "secant", 10, "1.808153474350061645017162432954535813",
	  "-4.252881499740056803204880343016444146" },
	{ "secant", 16, "1.6956DBBD2AF982CD89ADFAAB600F6D42ECFF",
	  "6.4DB3298187A97FFFFFFDF1A4A0F54586EBC7" },
	{ "secant", 16, "1.6956DBBD2AF982CD89ADFAAB600F6D465F89",
	  "6.4DB3298187A97FFFFFFDF1A4A0F5460E2AC6" },
	{ "secant", 16, "1.A91BD2911F9E35BC1E29D0393C7E81B9A295",
	  "-B.271B77B571C6FFFFFFFFF58319DD00DB82A8" },
	{ "secant", 2, "0.110010001000010101011111110010100101",
	  "1.011010010100011010101101110101101111" },
	{ "secant", 2, "0.101101111100000101100000110111000010",
	  "1.010100111101101101011110010011100111" },
	{ "secant", 2, "0.111001100011111001000110010100010100",
	  "1.100110111000010000110001011100011011" },
	{ "secant", 7, "2.003441542613305141555436021551625602",
	  "-2.230334204035341002402565515351501223" },
	{ "secant", 7, "1.501254101350523125034631252446506352",
	  "-6.541524564664151623630561363650314304" },
	{ "secant", 7, "1.202334312254324035605132330356311530",
	  "3.434535152002015502661261146305563633" },
	{ "cosecant", 10, "0.152824924192343575626448518220390286",
	  "6.568975670921991749651168978993308683" },
	{ "cosecant", 10, "-0.116192241713201023464335795859801901",
	  "-8.625822520774125479192661167354438633" },
	{ "cosecant", 10, "-0.214089005952038585348731869229574299",
	  "-4.706827700697332694801384989554465612" },
	{ "cosecant", 16, "0.177C57ACFE44AAEF42B72EAF49A95D9333CF",
	  "A.EA616E8A05517FFFFFFEBF00EFE8BCC76451" },
	{ "cosecant", 16, "-0.1082AFD6DE9624DD820416478FFEE707B8D6",
	  "-F.841B630E660C800000002F49C1597E939509" },
	{ "cosecant", 2, "-0.110000001000000001100110101001010111",
	  "-1.011101101100011100101100010110101110" },
	{ "cosecant", 2, "0.100100011000000110100001111101111011",
	  "1.110110111001100011000111110001010011" },
	{ "cosecant", 2, "1.000010100111111000000100101010001101",
	  "1.001010001010110001010000011010000000" },
	{ "cosecant", 7, "-0.150036313116352466545412534303455361",
	  "-4.056235115102163440423246035406223302" },
	{ "cosecant", 7, "0.330251532322263354106044431240625362",
	  "2.056350506004332231031411345150125302" },
	{ "cosecant", 7, "0.114603142565553505530140550320601222",
	  "5.444040060666050436644206132600141233" },
	{ "cotangent", 10, "0.100117712567985825012779182986166295",
	  "9.954847690202797139136236238324561291" },
	{ "cotangent", 10, "2.953079292877836639333856698622639692",
	  "-5.241676633099016768824188810179225167" },
	{ "cotangent", 10, "0.214583254430729562690985228579618092",
	  "4.588447598500684954103687187792182365" },
	{ "cotangent", 16, "0.1740A1150FA80DD52F186E7E4A07D55F7BD6",
	  "A.FAB31BDD559F0000000103D248D93A8F1C92" },
	{ "cotangent", 16, "0.1FE1B9DA6D4BE2DD3BDC4434EA09D211F36B",
	  "7.FCF554BD58D93FFFFFFE7D00AA85467732C1" },
	{ "cotangent", 16, "3.0C9777174E3057D6F1A948D663932BB4AB1E",
	  "-A.CA7A8BEFFBF87FFFFFFFB9000C57A8DF053C" },
	{ "cotangent", 2, "10.001001010001110111100010010011010100",
	  "-0.101001011001101100001110100001101011" },
	{ "cotangent", 2, "10.101010101110000010001011100101111101",
	  "-1.111100101110001111000010011000101010" },
	{ "cotangent", 2, "10.101010101110000011001101111010001000",
	  "-1.111100101110010100000000100011011001" },
	{ "cotangent", 7, "2.335613600663412631033563503543521606",
	  "-1.233613546450106535255656544516604540" },
	{ "cotangent", 7, "0.224000432454250005034146615113531616",
	  "2.562143124160043524006240261424630645" },
	{ "cotangent", 7, "2.660054422655034650364256525464353130",
	  "-6.063111426351604510122302146036406010" },
};

// The same near a pole, answers from the base to three times its square:
// 39 HEAD's engine got wrong.
static const SQuotientPin kQuotientNearPole[] = {
	{ "tangent", 10, "-1.566479686221109756208202635566181793",
	  "-231.660193099284271056603759967826554942" },
	{ "tangent", 10, "1.566325855381242471597651683922052516",
	  "223.688565647005560776960919188192443686" },
	{ "tangent", 10, "1.552262409254657453760063698799027940",
	  "53.948955647060621743094088494574282196" },
	{ "tangent", 16, "-1.911630A783FD4014FC0596621780B3100E9B",
	  "-F6.D2631CBA724FFFFFFFFF7B473465D0D6CE86" },
	{ "tangent", 16, "-1.907FB7B636E4411C2A8FB32CF2B8F2294FAE",
	  "-9D.8A3B012181500000000225C92C4D5C068492" },
	{ "tangent", 2, "-1.010101100111000110100000110010001011",
	  "-100.001101100010100001010100011100111100" },
	{ "tangent", 2, "1.011011010011010100101010111111011101",
	  "110.111000101111000101011100100101111111" },
	{ "tangent", 7, "-1.362201656405056065161203045602021510",
	  "-136.234165255653204636060406603434464102" },
	{ "tangent", 7, "1.364215363265163316012505045102530203",
	  "256.245223005160232642230041320413336665" },
	{ "tangent", 7, "1.364145423363650551134351223641562121",
	  "251.626343460052651424353130542316216555" },
	{ "cotangent", 10, "0.016618306595258883388172667633538269",
	  "60.169063348026675441815314553233591666" },
	{ "cotangent", 10, "0.008235080392481463061439381440448160",
	  "121.428977826908507609005027669041083761" },
	{ "cotangent", 10, "3.137420152867408757387803110218346350",
	  "-239.663037414095327903851285187996267489" },
	{ "cotangent", 16, "0.01AC26CF80D6EA00A2D46284264BBA3D47FB",
	  "99.10AA53FF88E0000000002F4737D515FE7613" },
	{ "cotangent", 16, "0.011EB304CCE5995BF77E59CE8E2272861D61",
	  "E4.96238ABA10C000000000676A6669D7FE8A77" },
	{ "cotangent", 16, "3.22B4E53250F44D690C2A541680DDB9A4E0CD",
	  "-A6.1D094C962B6000000000F04FE39071D2069F" },
	{ "cotangent", 2, "10.111111111010111101000000000000001100",
	  "-110.111101000011001110010011100001010101" },
	{ "cotangent", 2, "10.110101110001001010000011101000000101",
	  "-11.001101110100110011011000000100001101" },
	{ "cotangent", 2, "0.010100110001100111101011001100110100",
	  "10.111110001011101100111100100101010011" },
	{ "cotangent", 7, "0.004220203634612231503304400323633536",
	  "142.154366516225120101122445642051302205" },
	{ "cotangent", 7, "3.064031643523536024315030424424505544",
	  "-254.142414324564142033503466356664246451" },
	{ "secant", 10, "1.565204008991054971138829157696283773",
	  "178.817665131346188900351990917976345923" },
	{ "secant", 10, "1.579236290056979306859430775534885606",
	  "-118.485334735356119040261545409460120987" },
	{ "secant", 16, "1.9415E1D7A817F7CFC99AF52A00F5AE2C18FF",
	  "-82.817CD47AA1F000000001513F4F1A4379E3B4" },
	{ "secant", 2, "1.000100111101001001001010111000111110",
	  "10.000111001000101110101011000001010001" },
	{ "secant", 2, "1.101110101001101101110000001101000101",
	  "-110.010110011001101101000001101101110110" },
	{ "secant", 2, "1.101110101001101111000100000010010001",
	  "-110.010110011000111000110111101100100100" },
	{ "secant", 7, "1.414162636400034363546564265442316326",
	  "-41.566342022642441344402015332540050304" },
	{ "secant", 7, "1.402141604410605160201636016405114614",
	  "-260.120101423066400554435435506121165012" },
	{ "secant", 7, "1.344552664415542114261265311363042150",
	  "30.334015002220331642503436120432065060" },
	{ "cosecant", 10, "-0.018246597296657078887104969386295142",
	  "-54.807780081686118478501156809964015685" },
	{ "cosecant", 10, "0.017124792614965607757844737869247956",
	  "58.397721979925080404427713379917784251" },
	{ "cosecant", 10, "0.017124792614965607757844737869500896",
	  "58.397721979925080404427713379055310984" },
	{ "cosecant", 16, "0.0139051C3E70F2CADC6C2DAFAA82D2CC026C",
	  "D1.5E1D52E6A3F800000000297C45E579101505" },
	{ "cosecant", 16, "-0.01B3743D56A3A426AF39A93E45A16DF0A9FF",
	  "-96.80613B69381FFFFFFFFF3934B1481D767557" },
	{ "cosecant", 2, "0.001111110001001100110010010101001010",
	  "100.000110011001101011001101001001111110" },
	{ "cosecant", 7, "0.031156604212546462603443606404241411",
	  "21.262522560023242015464455356563535624" },
	{ "cosecant", 7, "0.002453426010001251556461162244441556",
	  "241.536025503306115231022222244661460034" },
	{ "cosecant", 7, "0.010244324114345035332661031616125245",
	  "64.322636405464456210215612050321560510" },
};

// Small answers a hair from a half, a hundredth to a ten-thousandth of the
// base, rounded at 36 places plus their noughts: 24 HEAD's engine got wrong.
static const SQuotientPin kQuotientSmall[] = {
	{ "tangent", 10, "-0.056362335588513673947232591797095358",
	  "-0.0564220938425938189810791116390067392" },
	{ "tangent", 10, "0.004813145369040359396831576612812734",
	  "0.00481318253708428223019476551677924524" },
	{ "tangent", 10, "-0.005233890839320985706359776354293285",
	  "-0.00523393863157193361996404086108279415" },
	{ "tangent", 16, "0.0038A63A61B83241BCCDEC3081C657F07937",
	  "0.0038A63B4E6EE08C0000066BF6CC96FE93DDBD" },
	{ "tangent", 16, "3.25E82F7F2BBFA9C0392470BF9B84AE7A3FDB",
	  "0.01A8C67C7867FE7FC53DF390781BE1596E0BE" },
	{ "tangent", 16, "3.161533E31E5D8AEDFC02DB9A54B4E86C3A31",
	  "-0.0E2DEB2D1E2981FFDEDA5035498D9984268DE" },
	{ "tangent", 2, "11.000011100111110011110000010100100111",
	  "-0.000101011100111111101110011100101100010" },
	{ "tangent", 2, "11.000011100111110011111001111101110010",
	  "-0.000101011100111111100100101111000010100" },
	{ "tangent", 2, "10.110011000110001000111011100000001111",
	  "-0.0101101101111100000111000100011001110" },
	{ "tangent", 7, "-0.044013661460442303505502556410561514",
	  "-0.0440610153540465333145266325034361365" },
	{ "tangent", 7, "0.003613564013425555120350030216322251",
	  "0.00361362243046230024440012156215325351" },
	{ "tangent", 7, "3.160330135351031345162053665006045242",
	  "0.0611022205316321004324035435355161633" },
	{ "cotangent", 10, "1.625286883761341566178944675151265978",
	  "-0.0545445525956857402194603476656106427" },
	{ "cotangent", 10, "1.573756990513655910671824678424172307",
	  "-0.00296067236938483136005903999458675363" },
	{ "cotangent", 10, "1.648703171753686866612536305684447916",
	  "-0.0780648464858346269432340314096285890" },
	{ "cotangent", 16, "1.9D4426D63A3A2CA5EFAD155C0399B6517ABE",
	  "-0.0B263F09259843FFF6BC5FDC388B62508DBE9" },
	{ "cotangent", 16, "1.88117491989C972B681E87C3491774D13EF4",
	  "0.0A0F93D6CDB9420017A5FFC5BEE0F06CE22B9" },
	{ "cotangent", 16, "1.88117491989C972B681E87C3491774D2AC41",
	  "0.0A0F93D6CDB9420017A5FFC5BEE0F06B744E3" },
	{ "cotangent", 2, "1.101010011001000100000100001100100010",
	  "-0.000101111000001000100011101110011000101" },
	{ "cotangent", 2, "1.101010011001000100011010110100100100",
	  "-0.000101111000001000111010100010101000010" },
	{ "cotangent", 2, "1.101101101011000011111000010011000010",
	  "-0.00100100110100010111010000110011111100" },
	{ "cotangent", 7, "1.360430005162360310205100506160114064",
	  "0.00610255363664063133621632444225541645" },
	{ "cotangent", 7, "1.400336655120110243354254602043201665",
	  "-0.000504233443435635151221566146405224462" },
	{ "cotangent", 7, "1.354454123326130065351356264644322156",
	  "0.0120461041321114136012504215133014121" },
};

static void TestQuotientGuard()
{
	Banner( "The circular quotients' guard" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	struct { const char* m_szName; Unary m_fn; } kUnary[] =
	{
		{ "tangent", lz_tangent }, { "cotangent", lz_cotangent },
		{ "secant", lz_secant }, { "cosecant", lz_cosecant }
	};

	struct { const char* m_szLabel; const SQuotientPin* m_pPins; size_t m_nPins; } kGroups[] =
	{
		{ "random operands HEAD misrounded", kQuotientRandom, sizeof( kQuotientRandom ) / sizeof( kQuotientRandom[ 0 ] ) },
		{ "a hair from a half", kQuotientNearHalf, sizeof( kQuotientNearHalf ) / sizeof( kQuotientNearHalf[ 0 ] ) },
		{ "a hair from a half, near a pole", kQuotientNearPole, sizeof( kQuotientNearPole ) / sizeof( kQuotientNearPole[ 0 ] ) },
		{ "small answers a hair from a half", kQuotientSmall, sizeof( kQuotientSmall ) / sizeof( kQuotientSmall[ 0 ] ) }
	};

	long	l_saved = glWorkingPrecision;

	for ( auto& group : kGroups )
	{
		long	l_wrong = 0;

		for ( size_t i = 0; i < group.m_nPins; i++ )
		{
			const SQuotientPin&	pin = group.m_pPins[ i ];

			glWorkingPrecision = 36;

			CLZNumber	x, answer, expect;

			x.SetBase( pin.m_lBase, true );
			answer.SetBase( pin.m_lBase, true );
			expect.SetBase( pin.m_lBase, true );
			x.Assign( pin.m_szOperand );
			expect.Assign( pin.m_szExpect );

			long	l_error = kFailed;

			for ( auto& u : kUnary )
				if ( ::strcmp( u.m_szName, pin.m_szFunction ) == 0 )
					l_error = u.m_fn( &x, &answer );

			if ( (l_error != kNoError) || (lz_compare( &answer, &expect ) != kEqualTo) )
				l_wrong++;
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrong, group.m_nPins );
		CheckTrue( group.m_szLabel, l_wrong == 0, szDetail );
	}

	glWorkingPrecision = l_saved;
}

/////////////////////////////////////////////////////////////////
//	Taken again when in doubt. Since 2026-10-03 a function whose final
//	rounding discards a guard lying within its doubt window of the half
//	is taken again wider, and that answer decides -- see
//	CLZAnswerBase::Retake. Until then every guard-rounded function
//	misrounded within its own blind spot of a half. Expected strings are
//	the exact values rounded, computed in Python.
/////////////////////////////////////////////////////////////////

struct SRetakePin
{
	const char*	m_szFunction;
	long		m_lBase;
	long		m_lPlaces;
	const char*	m_szOperand;
	const char*	m_szSecond;
	const char*	m_szExpect;
};

// Answers a hair from a half, two per function from the near-half census:
// 63 HEAD's engine got wrong. The arccotangent's negative one was replaced
// on 2026-10-23, when its convention became arctan(1/x): found by the same
// scan of x0 + j units, 9.4e-8 of a unit below a half.
static const SRetakePin kRetakeNearHalf[] = {
	{ "exp_e", 10, 36, "-0.240001887070314450438812772047895605", nullptr,
	  "0.786626376645868874426525075740029559" },
	{ "exp_e", 10, 36, "1.765643628616518195635454363241509825", nullptr,
	  "5.845333371439195069285688762074160102" },
	{ "exp_b", 10, 36, "0.787268946452098209797864126859746717", nullptr,
	  "6.127297206743289592623117240390884501" },
	{ "exp_b", 10, 36, "0.985150596827329687985173180893033322", nullptr,
	  "9.663859268216196340972601295161300783" },
	{ "log_b", 10, 36, "0.183846696969644126187183268141570122", nullptr,
	  "-0.735544168355593108721279605748596893" },
	{ "log_b", 10, 36, "45428.910097613918898492243070338406378515", nullptr,
	  "4.657332317614316430365306506253410041" },
	{ "sine", 10, 36, "-0.200458944759878472872749278366917888", nullptr,
	  "-0.199119106276614082906561783293184299" },
	{ "sine", 10, 36, "-0.200458944759878472872749278366589194", nullptr,
	  "-0.199119106276614082906561783292862187" },
	{ "cosine", 10, 36, "1.295543252405763731401218879465856792", nullptr,
	  "0.271790493426653512187846832555283101" },
	{ "cosine", 10, 36, "1.310776552879654468605306415865870818", nullptr,
	  "0.257099661155412106960939899780208919" },
	{ "tangent", 10, 36, "1.355671428881279487559029960918325066", nullptr,
	  "4.576531886622793088292837332187239049" },
	{ "tangent", 10, 36, "-1.375703497247042481433964426260606393", nullptr,
	  "-5.060568459038708333252997124792182681" },
	{ "secant", 10, 36, "1.927371818197481647552744792965288424", nullptr,
	  "-2.864778004410112544064759131214147707" },
	{ "secant", 10, 36, "1.808153474350061645017162432954877768", nullptr,
	  "-4.252881499740056803204880343010432613" },
	{ "cosecant", 2, 36, "1.010000111001000001110011011000011111", nullptr,
	  "1.000011001000101110110101000100010001" },
	{ "cotangent", 10, 36, "2.953079292877836639333856698622828566", nullptr,
	  "-5.241676633099016768824188810184603387" },
	{ "cotangent", 10, 36, "3.008331697217701383036131734615456231", nullptr,
	  "-7.459600418174738400978184489338040336" },
	{ "arcsine", 10, 36, "0.574053906063793642282651426009177695", nullptr,
	  "0.611448240120468389235953764412355605" },
	{ "arcsine", 10, 36, "-0.862863657372964372314743516546520745", nullptr,
	  "-1.040908272457002681221639689443359581" },
	{ "arccosine", 10, 36, "0.858988144265709513005573989268126238", nullptr,
	  "0.537506241103670845404849646228502815" },
	{ "arccosine", 10, 36, "-0.247452345768717994764779408242257105", nullptr,
	  "1.820846266109984012970812237182587680" },
	{ "arctangent", 10, 36, "0.680091128865406291598101110208752479", nullptr,
	  "0.597238970048342365792849496889893729" },
	{ "arctangent", 10, 36, "-1.596916768872019221454992760770695317", nullptr,
	  "-1.011329733703204558636911871771128461" },
	{ "arcsecant", 10, 36, "1.042129406430942849915667661518128774", nullptr,
	  "0.285312555784190235996790174919313611" },
	{ "arcsecant", 10, 36, "1.047787941654771461954039988141580940", nullptr,
	  "0.303181057486241134490967853158442187" },
	{ "arccosecant", 10, 36, "1.111118631647051266949631391231832692", nullptr,
	  "1.119755540127741033629143343988986567" },
	{ "arccosecant", 10, 36, "-2.491339891668487698309598571703821292", nullptr,
	  "-0.413034436936650872773667158150034319" },
	{ "arccotangent", 10, 36, "-0.918265736386524393226868592503133317", nullptr,
	  "-0.827980822121691007249719499482610727" },
	{ "arccotangent", 10, 36, "0.308828133622422297636847749283430582", nullptr,
	  "1.271260134399712971031705221322844864" },
	{ "hypersine", 10, 36, "-1.234997825005948343944512503779444583", nullptr,
	  "-1.573767823664004916750513866677361479" },
	{ "hypersine", 10, 36, "-2.424719643191765863098948491712553844", nullptr,
	  "-5.605279139243389607827339928128919577" },
	{ "hypercosine", 10, 36, "1.492756739614840393913290507498586933", nullptr,
	  "2.337048264501197625264693885857071267" },
	{ "hypercosine", 10, 36, "1.834827937352158460760487179730992739", nullptr,
	  "3.211848605127178846174727216094298287" },
	{ "hypertangent", 10, 36, "1.924615189621299590147923304097769428", nullptr,
	  "0.958295895479177817222193293306833452" },
	{ "hypertangent", 10, 36, "-0.806873136142256195271164156632929862", nullptr,
	  "-0.667861715210189554881026748932197007" },
	{ "hypersecant", 10, 36, "0.388763560436254506860638197070435294", nullptr,
	  "0.928914684377648680602133630783224685" },
	{ "hypersecant", 10, 36, "0.276060611330718342036970124825419315", nullptr,
	  "0.963068869438382546555034310822481760" },
	{ "hypercosecant", 10, 36, "-0.138069862700040353155933747034546889", nullptr,
	  "-7.219749611660446930727630491035272352" },
	{ "hypercosecant", 10, 36, "0.126290545328578580743733824783919866", nullptr,
	  "7.897239800567790624086226264511096295" },
	{ "hypercotangent", 10, 36, "0.135172834931657832727775327590978624", nullptr,
	  "7.442938956631902414073945283554925825" },
	{ "hypercotangent", 10, 36, "0.126279574238791808837264687287800234", nullptr,
	  "7.960985567052444977775852748090008054" },
	{ "archypersine", 10, 36, "21.863213275945838597960287498362419703", nullptr,
	  "3.778475248527843000090342468549803629" },
	{ "archypersine", 10, 36, "21.863213275945838597960287498362579340", nullptr,
	  "3.778475248527843000090342468549810923" },
	{ "archypercosine", 10, 36, "53.376380520134276480032584579673823557", nullptr,
	  "4.670427755850599105258423117064186919" },
	{ "archypercosine", 10, 36, "5908.034822959624720471669065188159256196", nullptr,
	  "9.377215711235926676181407499737066031" },
	{ "archypertangent", 10, 36, "-0.991400759470372092402783503167577895", nullptr,
	  "-2.722459840238766748399259075969023309" },
	{ "archypertangent", 10, 36, "0.968365270944158356994466119906531017", nullptr,
	  "2.065351563839043613768354105127161245" },
	{ "archypersecant", 10, 36, "0.765677267859548306287814325048857761", nullptr,
	  "0.763655256208763733027920916544995549" },
	{ "archypersecant", 10, 36, "0.186745479895035459399106914044197674", nullptr,
	  "2.362321092879934081310008254371874635" },
	{ "archypercosecant", 10, 36, "-0.120539121151323849861256947012493053", nullptr,
	  "-2.812540888915699088812516687549090855" },
	{ "archypercosecant", 10, 36, "-2.424385070347516174143267515632830831", nullptr,
	  "-0.401593677905549784057370762082546343" },
	{ "archypercotangent", 10, 36, "1.132797557933594409472847141304081049", nullptr,
	  "1.388181978299958041489504112793638591" },
	{ "archypercotangent", 10, 36, "1.409671847057175436590663329042247264", nullptr,
	  "0.885944693984642483997049482359397407" },
	{ "power", 10, 36, "2.968610777410386258878821608773603860", "3.250000000000000000000000000000000000",
	  "34.339826324561621257908516707492197685" },
	{ "power", 16, 36, "0.6C54BBF0CBB7DBFFFFFEB738DFA1F15F350B", "0.C00000000000000000000000000000000000",
	  "0.8650A4303F7881CB4B7109020D695E0F3E05" },
	{ "exp_n", 10, 36, "3.427172963326232668634930225969525966", "-1.500000000000000000000000000000000000",
	  "0.157614426257274965874396412035230105" },
	{ "exp_n", 10, 36, "5.865696623033923806644907746330789715", "2.500000000000000000000000000000000000",
	  "83.329541686326001317460119098925955784" },
	{ "root_n", 10, 36, "14.544676379712155878110025250114088076", "3.000000000000000000000000000000000000",
	  "2.441001357096362256832024689371850303" },
	{ "root_n", 2, 36, "1011.011010000011000011110100011111010011", "101.000000000000000000000000000000000000",
	  "1.101000001000111010001010100101110011" },
	{ "log_n", 10, 36, "3.000000000000000000000000000000000000", "21.319460840820799063521916847803636291",
	  "2.784986424824255517429294921911426315" },
	{ "log_n", 10, 36, "5.500000000000000000000000000000000000", "3.403431710656795949887510567295672839",
	  "0.718454682445195197578391807459053839" },
	{ "agm", 10, 36, "1.199848119176811778174848293065310711", "2.500000000000000000000000000000000000",
	  "1.790446682863065547839751170034804647" },
	{ "agm", 10, 36, "7.008909899155563394401723738821864930", "3.000000000000000000000000000000000000",
	  "4.792683766090288641111379306432691777" },
};

// Exact halves through power and root_n, which go to the even kept digit:
// 6 HEAD's engine got wrong.
static const SRetakePin kRetakeTie[] = {
	{ "power", 10, 5, "2.060531057025000000", "0.5",
	  "1.43546" },
	{ "power", 10, 36, "2.66221602133522557861189752480038956045134006559369815890065673656656475625000000", "0.5",
	  "1.631629866524643382059987018667703172" },
	{ "power", 10, 36, "1.04988910174393000634443587881732444291737638798578463878085328754097650625000000", "0.5",
	  "1.024640962358976748684565849281937302" },
	{ "power", 16, 5, "1.6F07D5EDEE40000000", "0.8",
	  "1.32876" },
	{ "power", 16, 5, "2.ACD9334A2440000000", "0.8",
	  "1.A2B6C" },
	{ "power", 16, 36, "1.16A97724B46E3424FAFEFF495F6EF404D9BF22F46EB4981689F2D92A7033A2F9803F8A7240000000", "0.8",
	  "1.0B173AFBB5FC65133757DC083C637E4695D2" },
};

// Random operands that landed in a blind spot -- e^x in base eleven, and
// two written past the width: 3, all HEAD's engine got wrong.
static const SRetakePin kRetakeRandom[] = {
	{ "exp_e", 11, 36, "2.2A29A422A2A936670928537815A540301406", nullptr,
	  "9.71314354811151072047703531178746273A" },
	{ "exp_e", 10, 36, "2.859514425374940965845486052012703727279", nullptr,
	  "17.453050120414775150223176150403011716" },
	{ "exp_e", 10, 36, "0.6984690882686517527233556676574122253146", nullptr,
	  "2.010672188429593474902438776578051397" },
};

// x to a whole power a hair from a half: 8 HEAD's engine got wrong.
static const SRetakePin kRetakeRaiseN[] = {
	{ "raise_n", 10, 36, "1.086508570138168261465240121084457247", "2",
	  "1.180500872983686900387058765808527943" },
	{ "raise_n", 10, 36, "1.513454621194644911241767146901325544", "5",
	  "7.940484802841039228206793834916383214" },
	{ "raise_n", 16, 36, "0.F5C21164F3F84FFFFFFF4AC01A41B0C64FCC", "2",
	  "0.EBED086999299DAFBD1E7CBBFDC81805D2DD" },
	{ "raise_n", 16, 36, "0.A1E3B0F06122AFFFFFFE7FB0673CF882DF77", "5",
	  "0.19E3D89FCB2ABDF2572F0CB6BDA48F5F049D" },
	{ "raise_n", 2, 36, "1.100010111101101011111101011011010001", "10",
	  "10.011001000001110110000101011000110011" },
	{ "raise_n", 2, 36, "0.110110111000001110100101111000001100", "10",
	  "0.101111000011101010000001100000011011" },
	{ "raise_n", 7, 36, "0.535462413404414242521423416220011200", "3",
	  "0.332325363060354556664526351550603011" },
	{ "raise_n", 7, 36, "1.635133013155046510405055221562642142", "10",
	  "203.021131233442155050404332101503441546" },
};

// x to a large whole power: the squaring chain's error grows with the
// exponent, and until 2026-10-10 nothing allowed for it. Six a hair from a
// half, then four to a hundred thousand and a million, powers be4eefc got
// as much as 169 units wrong and two places of widening for every
// exponent still gets wrong, then two random operands that power and
// exp_n, which square for their exponent's whole part, both rounded
// wrong: 14, all wrong at be4eefc, the last two in both functions.
static const SRetakePin kRetakeBigPower[] = {
	{ "raise_n", 2, 36, "1.000000000000000011000100001110010010", "11000011010100000",
	  "11.001110000111101010100000100110100100" },
	{ "raise_n", 2, 36, "1.000000000000000000110000101000001111", "11110100001001000000",
	  "10010.001001011000111110010000001111011000" },
	{ "raise_n", 10, 36, "1.000021565960373726156262819124769244", "100000",
	  "8.641470656248222453635275182693648109" },
	{ "raise_n", 10, 36, "0.999999283540370043077416539745557734", "1000000",
	  "0.488478467826911818769010018278075095" },
	{ "raise_n", 10, 36, "0.996655163382900499691174900362991213", "100",
	  "0.715306391253341095628741485061536840" },
	{ "raise_n", 10, 36, "0.999462772937801835482218848485231652", "1000",
	  "0.584282067559532703338987390639137018" },
	{ "raise_n", 10, 36, "0.999283796709503957522713335495077428", "1000",
	  "0.488478467835500216269325404108608107" },
	{ "raise_n", 16, 36, "0.FFA7EB9EF86707FFFFFFFFFFCD09DCC9D155", "3E8",
	  "0.42B492B22CC5F1C2AD3ABEAFE59DEB40780B" },
	{ "raise_n", 7, 36, "0.666212322520653544623042336013016516", "202",
	  "0.550500101561044122254046326102040622" },
	{ "raise_n", 7, 36, "1.001124316530541255443206650446543253", "2626",
	  "44.360340323243406510201360222300101313" },
	{ "power", 10, 36, "0.998507888412055622163006400136393492", "1000.753000000000000000000000000000000000",
	  "0.224394347706685391994088113129945932" },
	{ "exp_n", 10, 36, "0.998507888412055622163006400136393492", "1000.753000000000000000000000000000000000",
	  "0.224394347706685391994088113129945932" },
	{ "power", 10, 36, "1.003739859682550417409636988924269833", "1000.132000000000000000000000000000000000",
	  "41.820076996554972331127886983759512265" },
	{ "exp_n", 10, 36, "1.003739859682550417409636988924269833", "1000.132000000000000000000000000000000000",
	  "41.820076996554972331127886983759512265" },
};

// The geometric and harmonic means a hair from a half, the values written
// comma-separated: 14 HEAD's engine got wrong.
struct SRetakeListPin
{
	const char*	m_szFunction;
	long		m_lBase;
	long		m_lPlaces;
	const char*	m_szValues;
	const char*	m_szExpect;
};

static const SRetakeListPin kRetakeMean[] = {
	{ "gmean", 10, 36, "0.815167697530196910804539219643787023,5.487865290618841207503913707720321480,1.959325283400398021171896329041293910,0.532102356816443355924662088656173379,6.831866011209529965242383373320887962",
	  "1.998288982367064822898744946267731137" },
	{ "gmean", 10, 36, "2.313340007098957595843558037356076167,0.424422153628142451164959218774072189",
	  "0.990874738747073887336610424507046219" },
	{ "gmean", 16, 36, "6.4D53881C3720C0000001E03492C62E5EB905,7.8096E221E28A40000001284F5AFE824BA763",
	  "6.E0442AABA2216AC3587A7489DA49C2FDBEE5" },
	{ "gmean", 16, 36, "5.2F051DF71790C00000026E3B278BD132C6DA,6.586ACF24423EBFFFFFFF78C07FBAED24EAD6,8.56DEAE43018B80000000AC7EC4ED5098A06B",
	  "6.7F58233E9473D54F7338574A32A130BF6607" },
	{ "gmean", 2, 36, "11.111100000101100001111011100100001111,110.100110110011111111111111010101011000,10.111000110000100100110111100100100001",
	  "100.001110000010110011100101111001010001" },
	{ "gmean", 2, 36, "1.011110000100011111001000100110110010,100.010111011000001100111010111111110110,110.001001111000010100100111000101111111,10.100000000100110100101000011011010101",
	  "11.001001110000100101010000000001100111" },
	{ "gmean", 7, 36, "6.546633365504120562320562410250006225,0.542043506020500442250325342424244552",
	  "2.223651514045500326310506332120054360" },
	{ "gmean", 7, 36, "6.060042051552200441112313412404405305,11.632652266303645045352666541040343436,11.360242500656133414056466333040303331",
	  "10.521650623650312336456030604403420552" },
	{ "hmean", 10, 36, "4.169913190794079405065985437392797567,6.046542157309800380460274916607375228",
	  "4.935871599563762411925213540770254445" },
	{ "hmean", 10, 36, "6.416969617569199257900436288694290869,2.422766545391087600336277912948097326,4.738756452596948776090355355326935377",
	  "3.848058593577579752565316219626495735" },
	{ "hmean", 16, 36, "8.31E1F68410CBFFFFFFFFC95BA507D2735079,2.899DF33E3B84800000022486016BBABDEA6F",
	  "3.E00B415B73BDC6ACCDFDBC96C86A6ECD3B9F" },
	{ "hmean", 2, 36, "101.001110010001111001000100011100111110,101.000111110011001100010010000000001110",
	  "101.001011000000100000110011001100000111" },
	{ "hmean", 2, 36, "1.001000000001110000101000010000111001,10.000111100010010011000100111011010010,101.001100011110110011000000001111010010",
	  "1.111011100111001100101110000011111011" },
	{ "hmean", 16, 80, "8.31E1F68410CBFFFFFFFFC95BA507D2752A6F45E97EC609F40F62161B25A457CC1B77C3C2CE8C3682,2.899DF33E3B84800000022486016BBABDEA6E91F44E3885B8C2839AE4675204156F9ECBB9789729B2",
	  "3.E00B415B73BDC6ACCDFDBC96C86A6ECD709D08C04C0965A12288B1A13806AEC1ACA453D1CAA98CEB" },
};

typedef long (*RetakeUnary)( CLZNumber*, CLZNumber* );
typedef long (*RetakeBinary)( CLZNumber*, CLZNumber*, CLZNumber* );

static RetakeUnary RetakeUnaryNamed( const char* szName )
{
	struct { const char* m_szName; RetakeUnary m_fn; } kUnary[] =
	{
		{ "exp_e", lz_exp_e }, { "exp_b", lz_exp_b }, { "log_e", lz_log_e }, { "log_b", lz_log_b },
		{ "sine", lz_sine }, { "cosine", lz_cosine }, { "tangent", lz_tangent }, { "secant", lz_secant },
		{ "cosecant", lz_cosecant }, { "cotangent", lz_cotangent }, { "arcsine", lz_arcsine }, { "arccosine", lz_arccosine },
		{ "arctangent", lz_arctangent }, { "arcsecant", lz_arcsecant }, { "arccosecant", lz_arccosecant },
		{ "arccotangent", lz_arccotangent }, { "hypersine", lz_hypersine }, { "hypercosine", lz_hypercosine },
		{ "hypertangent", lz_hypertangent }, { "hypersecant", lz_hypersecant }, { "hypercosecant", lz_hypercosecant },
		{ "hypercotangent", lz_hypercotangent }, { "archypersine", lz_archypersine }, { "archypercosine", lz_archypercosine },
		{ "archypertangent", lz_archypertangent }, { "archypersecant", lz_archypersecant },
		{ "archypercosecant", lz_archypercosecant }, { "archypercotangent", lz_archypercotangent },
		{ "gamma", lz_gamma }, { "log_gamma", lz_log_gamma }
	};

	for ( auto& u : kUnary )
		if ( ::strcmp( u.m_szName, szName ) == 0 )
			return u.m_fn;

	return nullptr;
}

static RetakeBinary RetakeBinaryNamed( const char* szName )
{
	struct { const char* m_szName; RetakeBinary m_fn; } kBinary[] =
	{
		{ "power", lz_power }, { "exp_n", lz_exp_n }, { "root_n", lz_root_n }, { "log_n", lz_log_n }, { "agm", lz_agm },
		{ "raise_n", lz_raise_n }
	};

	for ( auto& u : kBinary )
		if ( ::strcmp( u.m_szName, szName ) == 0 )
			return u.m_fn;

	return nullptr;
}

// One pin: into a fresh target, or with none, so the answer is written
//	over the operand the retake has to take again...
static bool RetakePinRight( const SRetakePin& pin, bool bOverOperand )
{
	glWorkingPrecision = pin.m_lPlaces;

	CLZNumber	x, y, answer, expect;

	x.SetBase( pin.m_lBase, true );
	y.SetBase( pin.m_lBase, true );
	answer.SetBase( pin.m_lBase, true );
	expect.SetBase( pin.m_lBase, true );
	x.Assign( pin.m_szOperand );
	expect.Assign( pin.m_szExpect );

	CLZNumber*	p_target = bOverOperand ? nullptr : &answer;
	long		l_error = kFailed;

	if ( pin.m_szSecond != nullptr )
	{
		y.Assign( pin.m_szSecond );

		if ( RetakeBinary fn = RetakeBinaryNamed( pin.m_szFunction ) )
			l_error = fn( &x, &y, p_target );
	}
	else if ( RetakeUnary fn = RetakeUnaryNamed( pin.m_szFunction ) )
		l_error = fn( &x, p_target );

	return (l_error == kNoError) && (lz_compare( bOverOperand ? &x : &answer, &expect ) == kEqualTo);
}

// One mean: into a fresh target, or over its first value, so the retake
//	has to have kept a copy of the list...
static bool RetakeMeanRight( const SRetakeListPin& pin, bool bOverValue )
{
	glWorkingPrecision = pin.m_lPlaces;

	CLZNumber	values[ 8 ];
	CLZNumber*	p_values[ 8 ];
	long		l_count = 0;
	char		szValues[ 1024 ];

	::snprintf( szValues, sizeof( szValues ), "%s", pin.m_szValues );

	// Split at the commas by hand -- MSVC deprecates strtok...
	for ( char* p_token = szValues; (p_token != nullptr) && (*p_token != '\0') && (l_count < 8); )
	{
		char*	p_comma = ::strchr( p_token, ',' );

		if ( p_comma != nullptr )
			*p_comma = '\0';

		values[ l_count ].SetBase( pin.m_lBase, true );
		values[ l_count ].Assign( p_token );
		p_values[ l_count ] = &values[ l_count ];
		l_count++;

		p_token = (p_comma != nullptr) ? (p_comma + 1) : nullptr;
	}

	CLZNumber	answer, expect;

	answer.SetBase( pin.m_lBase, true );
	expect.SetBase( pin.m_lBase, true );
	expect.Assign( pin.m_szExpect );

	CLZNumber*	p_target = bOverValue ? &values[ 0 ] : &answer;
	long		l_error = (::strcmp( pin.m_szFunction, "gmean" ) == 0) ? lz_gmean( p_values, l_count, p_target )
																	: lz_hmean( p_values, l_count, p_target );

	return (l_error == kNoError) && (lz_compare( p_target, &expect ) == kEqualTo);
}

static void TestRetakeInDoubt()
{
	Banner( "Taken again when in doubt" );

	struct { const char* m_szLabel; const SRetakePin* m_pPins; size_t m_nPins; } kGroups[] =
	{
		{ "a hair from a half, 32 functions", kRetakeNearHalf, sizeof( kRetakeNearHalf ) / sizeof( kRetakeNearHalf[ 0 ] ) },
		{ "exact halves, power", kRetakeTie, sizeof( kRetakeTie ) / sizeof( kRetakeTie[ 0 ] ) },
		{ "random operands in a blind spot", kRetakeRandom, sizeof( kRetakeRandom ) / sizeof( kRetakeRandom[ 0 ] ) },
		{ "x to a whole power", kRetakeRaiseN, sizeof( kRetakeRaiseN ) / sizeof( kRetakeRaiseN[ 0 ] ) },
		{ "x to a large whole power", kRetakeBigPower, sizeof( kRetakeBigPower ) / sizeof( kRetakeBigPower[ 0 ] ) }
	};

	long	l_saved = glWorkingPrecision;

	for ( auto& group : kGroups )
	{
		long	l_wrong = 0;

		for ( size_t i = 0; i < group.m_nPins; i++ )
			if ( !RetakePinRight( group.m_pPins[ i ], false ) )
				l_wrong++;

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrong, group.m_nPins );
		CheckTrue( group.m_szLabel, l_wrong == 0, szDetail );
	}

	// The same answers written over their own operand, which the retake
	//	has to have kept a copy of...
	{
		long	l_wrong = 0;
		size_t	n_pins = sizeof( kRetakeNearHalf ) / sizeof( kRetakeNearHalf[ 0 ] );

		for ( size_t i = 0; i < n_pins; i++ )
			if ( !RetakePinRight( kRetakeNearHalf[ i ], true ) )
				l_wrong++;

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrong, n_pins );
		CheckTrue( "the same, over the operand", l_wrong == 0, szDetail );
	}

	// The means, into a target and over their first value, and a near-half
	//	mean retaken once...
	{
		size_t	n_pins = sizeof( kRetakeMean ) / sizeof( kRetakeMean[ 0 ] );
		long	l_wrong = 0;
		long	l_wrongOver = 0;

		for ( size_t i = 0; i < n_pins; i++ )
		{
			if ( !RetakeMeanRight( kRetakeMean[ i ], false ) )
				l_wrong++;

			if ( !RetakeMeanRight( kRetakeMean[ i ], true ) )
				l_wrongOver++;
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrong, n_pins );
		CheckTrue( "the geometric and harmonic means", l_wrong == 0, szDetail );
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong", l_wrongOver, n_pins );
		CheckTrue( "the same, over their first value", l_wrongOver == 0, szDetail );

		long	l_before = lzi_retake_count();

		RetakeMeanRight( kRetakeMean[ 0 ], false );
		::snprintf( szDetail, sizeof( szDetail ), "%ld retakes", lzi_retake_count() - l_before );
		CheckTrue( "a near-half mean retaken once", (lzi_retake_count() - l_before) == 1, szDetail );
	}

	// A near-half answer is taken again, an ordinary one is not...
	{
		long	l_before = lzi_retake_count();

		RetakePinRight( kRetakeNearHalf[ 0 ], false );

		long	l_near = lzi_retake_count() - l_before;

		glWorkingPrecision = 36;

		CLZNumber	x, answer;

		x.Assign( "0.5" );
		l_before = lzi_retake_count();
		lz_sine( &x, &answer );
		lz_exp_e( &x, &answer );
		lz_tangent( &x, &answer );

		long	l_plain = lzi_retake_count() - l_before;
		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld near the half, %ld ordinary", l_near, l_plain );
		CheckTrue( "retaken only when in doubt", (l_near == 1) && (l_plain == 0), szDetail );
	}

	glWorkingPrecision = l_saved;
}

/////////////////////////////////////////////////////////////////
//	Whole quotients and remainders. Until 2026-09-26 the integer
//	quotient and the three remainders took a whole number from a
//	quotient already rounded at the precision, so one just short of a
//	whole number or just past a half landed on the wrong side:
//	1 / 0.333...334 to fifty places, 2.999...994, had an integer
//	quotient of 3 and a modulo of 0. See lzi_divide_work and
//	lzi_remainder_by.
/////////////////////////////////////////////////////////////////

static void TestWholeQuotients()
{
	Banner( "Whole quotients and remainders" );

	/*
		Exact fractions in Python: the integer quotient truncates, the
		modulo and the remainder are left less that multiple of right,
		and the remainder-near takes the nearest multiple, ties to even.
		Every remainder here is exact, the fifty-place ones included --
		the product of the whole number and right used to be rounded at
		the precision. The first three pairs were wrong before; the ties
		and the large quotient pin what was already right.
	*/

	struct { const char* szLeft; const char* szRight; const char* szWhole; const char* szModulo; const char* szNear; } cases[] =
	{
		{ "1",  "0.33333333333333333333333333333333333333333333333334", "2",
		  "0.33333333333333333333333333333333333333333333333332", "-0.00000000000000000000000000000000000000000000000002" },
		{ "-1", "0.33333333333333333333333333333333333333333333333334", "-2",
		  "-0.33333333333333333333333333333333333333333333333332", "0.00000000000000000000000000000000000000000000000002" },
		{ "1",  "0.39999999999999999999999999999999999999999999999999", "2",
		  "0.20000000000000000000000000000000000000000000000002", "-0.19999999999999999999999999999999999999999999999997" },
		{ "5",  "2",   "2",  "1",  "1" },
		{ "7",  "2",   "3",  "1",  "-1" },
		{ "-5", "2",   "-2", "-1", "-1" },
		{ "1234567890123456789.25", "0.3", "4115226300411522630", "0.25", "-0.05" },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber	a, b, whole, modulo, remainder, nearest;

		Set( a, theCase.szLeft );
		Set( b, theCase.szRight );

		lz_divide_integer( &a, &b, &whole );
		lz_modulo( &a, &b, &modulo );
		lz_remainder( &a, &b, &remainder );
		lz_remainder_near( &a, &b, &nearest );

		struct { const char* szName; CLZNumber* pGot; const char* szWant; } answers[] =
		{
			{ "integer quotient", &whole, theCase.szWhole }, { "modulo", &modulo, theCase.szModulo },
			{ "remainder", &remainder, theCase.szModulo },   { "remainder-near", &nearest, theCase.szNear },
		};

		long		l_wrong = 0;
		const char*	szWhere = "";

		for ( auto& theAnswer : answers )
		{
			char*	szGot = Render( theAnswer.pGot );

			if ( (szGot == nullptr) || (AgreementDigits( szGot, theAnswer.szWant ) != kExactAgreement) )
			{
				l_wrong++;
				szWhere = theAnswer.szName;
			}

			::free( szGot );
		}

		char	szLabel[ 72 ];
		char	szDetail[ 72 ];

		::snprintf( szLabel, 72, "%s by %.12s%s", theCase.szLeft, theCase.szRight, (::strlen( theCase.szRight ) > 12) ? "..." : "" );

		if ( l_wrong == 0 )
			::snprintf( szDetail, 72, "all four exact" );
		else
			::snprintf( szDetail, 72, "%ld of 4 wrong, one the %s", l_wrong, szWhere );

		CheckTrue( szLabel, l_wrong == 0, szDetail );
	}

	/*
		Ties in base seven, which has no half digit: 3/2 there is
		1.3333... with the middle digit repeating forever, exactly half,
		which the sticky digit says and the rounding to a whole number
		has to be told -- see RoundMiddleRepeating. Half to even: 3/2
		goes to 2 and 5/2 to 2.
	*/

	{
		struct { const char* szLeft; const char* szWant; } ties[] = { { "3", "-1" }, { "5", "1" } };

		for ( auto& theTie : ties )
		{
			CLZNumber	a, two, answer;

			SetInBase( a, 7, theTie.szLeft );
			SetInBase( two, 7, "2" );
			answer.SetBase( 7, true );
			lz_remainder_near( &a, &two, &answer );

			char	szLabel[ 72 ];

			::snprintf( szLabel, 72, "base 7: %s remainder-near 2", theTie.szLeft );
			CheckExact( szLabel, &answer, theTie.szWant );
		}
	}

	// In place, where the left operand is the target...
	{
		CLZNumber	a, b;

		Set( a, "1" );
		Set( b, "0.33333333333333333333333333333333333333333333333334" );
		lz_modulo( &a, &b, nullptr );
		CheckExact( "modulo in place", &a, "0.33333333333333333333333333333333333333333333333332" );
	}

	/*
		**A whole number never takes the reciprocal path**, which rounds
		its quotient on the way -- the fault itself. With the floor
		lowered to reach it, the integer quotient of the first pair
		still has to be 2.
	*/

	{
		CLZNumber	a, b, whole;

		Set( a, "1" );
		Set( b, "0.33333333333333333333333333333333333333333333333334" );

		long	l_savedFloor = lzi_set_reciprocal_floor( 2 );

		lz_divide_integer( &a, &b, &whole );
		lzi_set_reciprocal_floor( l_savedFloor );

		CheckExact( "past the reciprocal floor", &whole, "2" );
	}
}

/////////////////////////////////////////////////////////////////
//	e as wide as it is asked for. kBaseE holds 120 places, and until
//	2026-09-26 that was a ceiling: e^n wants its own size plus 36
//	places, so from about e^200 the answer came out short, and e^400
//	was 1.4e91 units out of its 36th place. Past the constant's width
//	e is computed now -- see lzi_e_in_base -- and the logarithm takes
//	the same e.
/////////////////////////////////////////////////////////////////

static const char kRefExp200[] =
		"7225973768125749258177477042189305697356874428527319284032697891232219093614"
		"73891661561.926589062570557468402043101429418177";

static const char kRefExp400[] =
		"5221469689764143950588763006649649029523362098447666367117172852895098136833"
		"3315117675846591903196562537526214688618767945779835711395764964812364315511"
		"7352656311293271093617.487406702573567237021584962335200343";

static const char kRefSinh200[] =
		"3612986884062874629088738521094652848678437214263659642016348945616109546807"
		"36945830780.963294531285278734201021550714709089";

static const char kRefE310[] =
		"2.71828182845904523536028747135266249775724709369995957496696762772407663035"
		"3547594571382178525166427427466391932003059921817413596629043572900334295260"
		"5956307381323286279434907632338298807531952510190115738341879307021540891499"
		"3488416750924476146066808226480016847741185374234544243710753907774499206955"
		"17027618";

static const char kRefLn1e40At310[] =
		"92.1034037197618273607196581873745683040440595451509190413331160387029043870"
		"9409920943988820358393193367871136169144994533638101860331227026665149476395"
		"1267579316288333022187233751995793049327941135740212358615109305153846534648"
		"8915079287954698617466989761697309746062019573725975756591847761760088842040"
		"685669920";

// ln(87/7) -- 15.3 in base seven -- to 360 base-seven places, correctly
//	rounded, from Python's Decimal...
static const char kRefLnBase7[] =
		"2.34323414546410600213325515566656501051540520535451213523646424210621451463"
		"0123526642263063446402303245513634524254213522061063425140244511651563326653"
		"0213063340454213231423264610040612421665035310043065454505336034553352651451"
		"1454260361451226220365235231155532564301254303235046456165121310654303303364"
		"2111112344006523363251554053246525452505155123003225352326";

static void TestWideE()
{
	Banner( "e as wide as it is asked for" );

	/*
		Correctly rounded at the 36th place, from Python's Decimal. The
		audit sets aside answers past 1e30, so nothing else scores these:
		e^200 was 9,700 units out, e^400 1.4e91, sinh(200) 4,850.
	*/

	{
		struct { const char* szLabel; char cFn; const char* szValue; const char* szExpect; } pins[] =
		{
			{ "exp(200), 36 places",  'e', "200", kRefExp200 },
			{ "exp(400), 36 places",  'e', "400", kRefExp400 },
			{ "raise_e(400)",         'r', "400", kRefExp400 },
			{ "sinh(200), 36 places", 's', "200", kRefSinh200 },
		};

		for ( auto& thePin : pins )
		{
			CLZNumber	x, answer;

			x.Assign( thePin.szValue );

			if ( thePin.cFn == 'e' )
				lz_exp_e( &x, &answer );
			else if ( thePin.cFn == 'r' )
				lz_raise_e( &x, &answer );
			else
				lz_hypersine( &x, &answer );

			CheckRendering( thePin.szLabel, &answer, thePin.szExpect );
		}
	}

	/*
		e and ln(1e40) at 300 places, past the constant's 120, against
		references rather than against each other. The logarithm's
		reduction divides 1e40 by e some ninety times, so its e is what
		its answer rests on -- ln 2, under e, would not touch it. **A round trip would
		prove nothing**: while the exponential and the logarithm both
		took e from the same 120 places its error cancelled between
		them, and power(1e40, 2.5) came out exactly 1e100 with the
		logarithm 3e16 units short on its own.
	*/

	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = 300;

		CLZNumber	one, big, e, logarithm;

		one.Assign( 1L );
		big.Assign( "10000000000000000000000000000000000000000" );
		lz_exp_e( &one, &e );
		lz_log_e( &big, &logarithm );

		glWorkingPrecision = l_saved;

		CheckDigits( "e at 300 places", &e, kRefE310, 300 );
		CheckDigits( "ln(1e40) at 300 places", &logarithm, kRefLn1e40At310, 300 );

		CLZNumber	base, exponent, power;

		base.Assign( "10000000000000000000000000000000000000000" );
		exponent.Assign( "2.5" );
		lz_power( &base, &exponent, &power );
		CheckExact( "power(1e40, 2.5), now honestly", &power,
					"10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" );
	}

	/*
		**Independent routes that must agree**, counted in base ten.

		exp(1) takes e by binary splitting; exp(0.5) sums its series term
		by term, halved sixty-four times past a hundred digits and
		squared back -- see kExpComputedHalvings -- and squared again it
		is e. At 1,200 places they must agree to all but the few digits
		the squaring spends.

		Base seven has no finite half, so it is checked two other ways at
		1,300 places: its e, split in base seven, against base ten's; and
		exp(3/7) squared against exp(6/7), both exact there. An odd base
		does not halve -- the argument would fill the width -- so that
		pair is the plain series at a width it could not reach before its
		bound grew. Without the halvings' extra places the base-ten pair
		parts about fifteen digits early.

		**exp splits a long fraction** into a head it halves and a tail it
		does not -- see kExpSplitShare -- and nothing above reaches that:
		0.5 is one place long. A third written
		out to the width, cubed, is e to within a unit of the last place,
		and ln(e) is Newton on exactly that exponential. Base seven agrees
		too, to every digit at 1,300 places, and is left out because its
		two cost a second.
	*/

	{
		// A value at the width it was made, written out in base ten...
		auto Decimal = []( CLZNumber& theValue ) -> char*
		{
			long		l_saved = glWorkingPrecision;
			CLZNumber	copy;

			glWorkingPrecision = 3000;
			copy.Assign( &theValue );
			copy.SetPrecision( 3000 );
			copy.SetBase( 10, true );
			glWorkingPrecision = l_saved;

			return copy.GetStringRep();
		};

		auto Agree = [ & ]( const char* szLabel, CLZNumber& theFirst, CLZNumber& theSecond, int iWant )
		{
			char*	szFirst = Decimal( theFirst );
			char*	szSecond = Decimal( theSecond );
			int		iAgree = AgreementDigits( szFirst, szSecond );
			char	szDetail[ 72 ];

			::snprintf( szDetail, 72, "%d digits agree, %d wanted", iAgree, iWant );
			CheckTrue( szLabel, iAgree >= iWant, szDetail );

			::free( szFirst );
			::free( szSecond );
		};

		long	l_saved = glWorkingPrecision;

		CLZNumber	one, half, e, root, squared;

		glWorkingPrecision = 1200;

		one.Assign( "1" );
		half.Assign( "0.5" );
		lz_exp_e( &one, &e );
		lz_exp_e( &half, &root );
		lz_multiply( &root, &root, &squared );

		glWorkingPrecision = l_saved;

		Agree( "exp(0.5)^2 against exp(1), 1200", squared, e, 1195 );

		CLZNumber	one7, e7, three7, six7, root7, whole7, squared7;

		glWorkingPrecision = 1300;

		SetInBase( one7, 7, "1" );
		SetInBase( three7, 7, "0.3" );
		SetInBase( six7, 7, "0.6" );
		lz_exp_e( &one7, &e7 );
		lz_exp_e( &three7, &root7 );
		lz_exp_e( &six7, &whole7 );
		lz_multiply( &root7, &root7, &squared7 );

		glWorkingPrecision = l_saved;

		Agree( "base 7 e against base 10's", e7, e, 1095 );
		Agree( "base 7: exp(3/7)^2 against exp(6/7)", squared7, whole7, 1093 );

		// A fraction as long as the width, cubed, against e; and ln(e)...
		for ( long lBase : { 10L, 16L } )
		{
			long		l_places = (lBase == 10) ? 1200 : 1000;
			std::string	thirds = std::string( "0." ) + std::string( static_cast<size_t>( l_places ), (lBase == 10) ? '3' : '5' );
			CLZNumber	oneIn, eIn, third, root3, square3, cube3, logE;
			char		szLabel[ 72 ];

			glWorkingPrecision = l_places;

			SetInBase( oneIn, lBase, "1" );
			SetInBase( third, lBase, thirds.c_str() );
			lz_exp_e( &oneIn, &eIn );
			lz_exp_e( &third, &root3 );
			lz_multiply( &root3, &root3, &square3 );
			lz_multiply( &square3, &root3, &cube3 );
			lz_log_e( &eIn, &logE );

			bool	b_widthKept = (glWorkingPrecision == l_places);

			glWorkingPrecision = l_saved;

			::snprintf( szLabel, 72, "base %ld: exp(0.%c...)^3 against e, %ld", lBase, thirds[ 2 ], l_places );
			Agree( szLabel, cube3, eIn, (lBase == 10) ? 1195 : 1200 );
			::snprintf( szLabel, 72, "base %ld: ln(e) against 1, %ld", lBase, l_places );
			Agree( szLabel, logE, one, (lBase == 10) ? 1195 : 1200 );
			::snprintf( szLabel, 72, "base %ld: exp and ln leave the width alone", lBase );
			CheckTrue( szLabel, b_widthKept, b_widthKept ? "as they found it" : "moved" );
		}
	}

	/*
		**The logarithm past 200 digits is Newton on the exponential**,
		and until 2026-09-26 it was the series, whose cost grows as the
		cube of the width: ln(1e40) at 9,000 places took 7.4 s and takes
		0.3 now. The two are independent -- the series touches no
		exponential -- so the series is forced past the width and the two
		required to agree to the character: a tiny argument, which Newton
		inverts, one just above one, a large one, in base ten at 500
		places and hex at 420.

		**Two more rows are there because mutations survived without
		them.** 1e-40 is the argument the inversion exists for: without
		it Newton works on an exponential that fixed places leave forty
		digits short, and 67 of 200 random small arguments came back
		wrong, where 1.23e-7 lost only what the guard absorbs. The three
		at 320 places are the worst of 200 random arguments with the
		guard taken out -- 15 missed correct rounding, by up to 0.94 of a
		unit, nearly all just above one -- and the eight above caught none
		of them. The hex row at 200 places is the series' own: its one and
		two were assigned as doubles, which came back at 37 hex places of
		44, and this argument came out 9e29 units short.

		**Hex is where they first disagreed**, and it was the series that
		was wrong: it counted its terms from places of the operand's base
		where the count is in decimal digits, and ran 20% short in hex --
		ln(1000.5) at 420 hex places was 8e19 units out of its last place.
		Base seven, where the series is slow, is pinned against Python
		instead, at 360 places.

		**Base two had no row, and that is where Newton was wrong.** Its
		guard was eight places, which in base two is eight bits, so from
		680 binary places the logarithm lost up to 1.5e-3 of a unit: 28 of
		40,000 arguments 1 + k/4096 at 700 and 1000 places came back
		different from the series, and Python says Newton was wrong in
		every one. The five rows here are the deepest of them, 1.3e-4 to
		8.8e-4 of a unit from the half.
	*/

	{
		struct { long lBase; long lPlaces; const char* szValue; } logs[] =
		{
			{ 10, 500, "2" }, { 10, 500, "0.000000123" }, { 10, 500, "1.0000000001" },
			{ 10, 500, "10000000000000000000000000000000000000000" }, { 10, 500, "12345.678" },
			{ 16, 420, "2" }, { 16, 420, "3E8.8" }, { 16, 420, "0.0001F" },
			{ 16, 200, "0.E91EF427CE49B5" },
			{ 10, 500, "0.0000000000000000000000000000000000000001" },
			{ 10, 320, "1.0073411494285415279882" }, { 10, 320, "1.0006376266035671965754" },
			{ 10, 320, "7.6680718684450803447071603" },
			{ 2, 700, "10.10110100011" }, { 2, 700, "100.110100001" },
			{ 2, 1000, "10.00001101111" }, { 2, 1000, "1.100000011011" }, { 2, 1000, "1.00001" },
		};

		const long	l_logs = static_cast<long>( sizeof( logs ) / sizeof( logs[ 0 ] ) );

		long		l_differ = 0;
		const char*	szWhere = "";

		for ( auto& theLog : logs )
		{
			std::string	s_answers[ 2 ];

			for ( int iPath = 0; iPath < 2; iPath++ )
			{
				long	l_newtonWas = lzi_set_log_newton_from( (iPath == 0) ? kLogNewtonOutOfReach : 200 );
				long	l_saved = glWorkingPrecision;

				glWorkingPrecision = theLog.lPlaces;

				CLZNumber	x, answer;

				x.SetBase( theLog.lBase, true );
				x.Assign( theLog.szValue );
				lz_log_e( &x, &answer );

				glWorkingPrecision = l_saved;
				lzi_set_log_newton_from( l_newtonWas );

				char*	szAnswer = answer.GetStringRep();

				s_answers[ iPath ] = szAnswer;
				::free( szAnswer );
			}

			if ( s_answers[ 0 ] != s_answers[ 1 ] )
			{
				l_differ++;
				szWhere = theLog.szValue;
			}
		}

		char	szDetail[ 72 ];

		if ( l_differ == 0 )
			::snprintf( szDetail, 72, "all %ld identical to the character", l_logs );
		else
			::snprintf( szDetail, 72, "%ld of %ld differ, one at %s", l_differ, l_logs, szWhere );

		CheckTrue( "ln by Newton against ln by the series", l_differ == 0, szDetail );

		long		l_saved = glWorkingPrecision;
		CLZNumber	x7, answer7;

		glWorkingPrecision = 360;
		x7.SetBase( 7, true );
		x7.Assign( "15.3" );
		lz_log_e( &x7, &answer7 );
		glWorkingPrecision = l_saved;

		CheckRendering( "base 7: ln 15.3 at 360 places", &answer7, kRefLnBase7 );
	}

	/*
		**Newton forced to the hook's floor, in base fifteen.** The seed's
		sixty digits are 52 places there, which read back as 62 digits, past
		the floor of 61: until 2026-10-17 the seed asked Newton for itself
		and the logarithm recursed until the stack ran out, so the old code
		fails this by crashing. Every other base reads back at 61 or less.
		The AGM, which would take this width in base fifteen, is held off.
	*/

	{
		long	l_saved = glWorkingPrecision;
		long	l_places = CLZNumber::EquivalentPlaces( 100, 10, 15 );
		long	l_agmWas = lzi_set_log_agm_floor( LONG_MAX );
		std::string	s_answers[ 2 ];

		for ( int iPath = 0; iPath < 2; iPath++ )
		{
			long	l_newtonWas = lzi_set_log_newton_from( (iPath == 0) ? kLogNewtonOutOfReach : 0 );

			glWorkingPrecision = l_places;

			CLZNumber	x, answer;

			x.SetBase( 15, true );
			x.Assign( "2" );
			lz_log_e( &x, &answer );

			glWorkingPrecision = l_saved;
			lzi_set_log_newton_from( l_newtonWas );

			char*	szAnswer = answer.GetStringRep();

			s_answers[ iPath ] = szAnswer;
			::free( szAnswer );
		}

		lzi_set_log_agm_floor( l_agmWas );

		CheckTrue( "base 15: ln 2 by Newton forced to its floor, against the series",
				   (s_answers[ 0 ] == s_answers[ 1 ]) && (s_answers[ 0 ].size() > static_cast<size_t>( l_places )),
				   s_answers[ 1 ].substr( 0, 40 ).c_str() );
	}
}

static void TestAnswerBases()
{
	Banner( "Every answer in the left operand's base" );

	typedef long (*Unary)( CLZNumber*, CLZNumber* );
	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );
	typedef long (*OfList)( CLZNumber**, long, CLZNumber* );

	#define ENTRY( f ) { #f, f }

	struct { const char* szName; Unary fn; } unary[] =
	{
		ENTRY( lz_abs ),				ENTRY( lz_arccosecant ),		ENTRY( lz_arccosine ),
		ENTRY( lz_arccotangent ),		ENTRY( lz_archypercosecant ),	ENTRY( lz_archypercosine ),
		ENTRY( lz_archypercotangent ),	ENTRY( lz_archypersecant ),		ENTRY( lz_archypersine ),
		ENTRY( lz_archypertangent ),	ENTRY( lz_arcsecant ),			ENTRY( lz_arcsine ),
		ENTRY( lz_arctangent ),			ENTRY( lz_cosecant ),			ENTRY( lz_cosine ),
		ENTRY( lz_cotangent ),			ENTRY( lz_cube ),				ENTRY( lz_exp_b ),
		ENTRY( lz_exp_e ),				ENTRY( lz_factorial ),			ENTRY( lz_log_gamma ),
		ENTRY( lz_gamma ),				ENTRY( lz_fraction ),			ENTRY( lz_hypersine ),
		ENTRY( lz_hypercosine ),		ENTRY( lz_hypertangent ),		ENTRY( lz_hypercosecant ),
		ENTRY( lz_hypercotangent ),		ENTRY( lz_hypersecant ),		ENTRY( lz_integer ),
		ENTRY( lz_log_b ),				ENTRY( lz_log_e ),				ENTRY( lz_negate ),
		ENTRY( lz_normalize ),			ENTRY( lz_plus ),				ENTRY( lz_raise_e ),
		ENTRY( lz_reciprocal ),			ENTRY( lz_root_2 ),				ENTRY( lz_root_3 ),
		ENTRY( lz_round ),				ENTRY( lz_round_integral ),		ENTRY( lz_round_last_digit ),
		ENTRY( lz_secant ),				ENTRY( lz_sine ),				ENTRY( lz_square ),
		ENTRY( lz_tangent ),
	};

	struct { const char* szName; Binary fn; } binary[] =
	{
		ENTRY( lz_add ),				ENTRY( lz_subtract ),			ENTRY( lz_multiply ),
		ENTRY( lz_divide ),				ENTRY( lz_divide_integer ),		ENTRY( lz_modulo ),
		ENTRY( lz_remainder ),			ENTRY( lz_remainder_near ),		ENTRY( lz_agm ),
		ENTRY( lz_combinations ),		ENTRY( lz_exp_n ),				ENTRY( lz_log_n ),
		ENTRY( lz_power ),				ENTRY( lz_quantize ),			ENTRY( lz_raise_n ),
		ENTRY( lz_root_n ),				ENTRY( lz_permutations ),
	};

	struct { const char* szName; OfList fn; } ofList[] =
	{
		ENTRY( lz_amean ),				ENTRY( lz_gmean ),				ENTRY( lz_hmean ),
	};

	#undef ENTRY

	const long	kBases[] = { 16, 7 };
	const char*	kValues[] = { "0.5", "5", "1.5", "-3" };
	const char*	kPairs[][ 2 ] = { { "21", "3" }, { "2.5", "0.5" }, { "7", "2" } };

	// A value written in decimal and converted, so it is the same
	//	number in every base...
	auto In = []( CLZNumber& theNumber, long lBase, const char* szDecimal )
	{
		theNumber.Assign( szDecimal );
		theNumber.SetBase( lBase, true );
	};

	for ( auto& theEntry : unary )
	{
		for ( long lBase : kBases )
		{
			int		iWrong = 0;

			for ( const char* szValue : kValues )
			{
				CLZNumber	x, target;
				In( x, lBase, szValue );
				target.SetBase( 2, false );
				theEntry.fn( &x, &target );
				iWrong += (target.m_lBase != lBase);

				CLZNumber	inPlace;
				In( inPlace, lBase, szValue );
				theEntry.fn( &inPlace, nullptr );
				iWrong += (inPlace.m_lBase != lBase);
			}

			char szLabel[ 72 ];
			::snprintf( szLabel, 72, "%s answers in base %ld", theEntry.szName, lBase );
			CheckTrue( szLabel, iWrong == 0, (iWrong == 0) ? "" : "an answer in another base" );
		}
	}

	for ( auto& theEntry : binary )
	{
		for ( long lBase : kBases )
		{
			int		iWrong = 0;

			for ( auto& thePair : kPairs )
			{
				// The odd base on the left and then on the right...
				for ( int iSide = 0; iSide < 2; iSide++ )
				{
					long	lLeft = (iSide == 0) ? lBase : 10;
					long	lRight = (iSide == 0) ? 10 : lBase;

					CLZNumber	x, y, target;
					In( x, lLeft, thePair[ 0 ] );
					In( y, lRight, thePair[ 1 ] );
					target.SetBase( 2, false );
					theEntry.fn( &x, &y, &target );
					iWrong += (target.m_lBase != lLeft);

					CLZNumber	u, w;
					In( u, lLeft, thePair[ 0 ] );
					In( w, lRight, thePair[ 1 ] );
					theEntry.fn( &u, &w, nullptr );
					iWrong += (u.m_lBase != lLeft);
				}
			}

			char szLabel[ 72 ];
			::snprintf( szLabel, 72, "%s answers in the left's base, %ld", theEntry.szName, lBase );
			CheckTrue( szLabel, iWrong == 0, (iWrong == 0) ? "" : "an answer in another base" );
		}
	}

	for ( auto& theEntry : ofList )
	{
		for ( long lBase : kBases )
		{
			CLZNumber	first, second, third, target;
			In( first, lBase, "2" );
			In( second, 10, "3" );
			In( third, 2, "5" );

			CLZNumber*	values[] = { &first, &second, &third };

			target.SetBase( 2, false );
			theEntry.fn( values, 3, &target );

			char szLabel[ 72 ];
			::snprintf( szLabel, 72, "%s answers in the first's base, %ld", theEntry.szName, lBase );
			CheckTrue( szLabel, target.m_lBase == lBase, "" );
		}
	}

	/*
		And the values survive the trip. Sine of a hexadecimal one is
		sin(1) to the thirty-six places the base-ten answer had; the
		log to base 0x10 of 256 is exactly 2, written in hex; and
		ten-choose-three in hex is 0x78.
	*/

	{
		CLZNumber	x, result;
		In( x, 16, "1" );
		lz_sine( &x, &result );
		CheckTrue( "sine of 0x1 answers in hex", result.m_lBase == 16, "" );
		result.SetBase( 10, true );
		CheckDigits( "sine of 0x1 is sin(1)", &result, kRefSin1, 35 );
	}

	{
		CLZNumber	n, x, result;
		In( n, 16, "16" );
		In( x, 10, "256" );
		lz_log_n( &n, &x, &result );
		CheckExact( "log base 0x10 of 256, in hex", &result, "2" );
		CheckTrue( "log base 0x10 of 256 answers in hex", result.m_lBase == 16, "" );
	}

	{
		CLZNumber	n, k, result;
		In( n, 16, "10" );
		In( k, 10, "3" );
		lz_combinations( &n, &k, &result );
		CheckTrue( "0xA choose 3 answers in hex", result.m_lBase == 16, "" );
		CheckRendering( "0xA choose 3", &result, "78" );
	}
}

/////////////////////////////////////////////////////////////////
//	Guards counted in digits, and widening for an answer's size.
//	Every guard was a count of decimal digits added as places of the
//	operand's base -- four bits in base two. And exp and lz_raise_e
//	widened for their answer's size only once it outgrew the guard,
//	so an answer about as large as the guard spent all of it.
/////////////////////////////////////////////////////////////////

typedef long (*GuardedFunction)( CLZNumber*, CLZNumber* );

struct SGuardCase
{
	const char*		szName;
	GuardedFunction	pFunction;
	long			lBase;
	const char*		szArgument;
	long			lPlaces = 0;	// 36 when left out
};

// The answer at 36 places, or the case's own, against the same function
//	worked twenty places wider and rounded back to the places the first
//	one kept.
//	A few units of the 56th place cannot move the 36th, so this is
//	correct rounding checked without a reference to type...
static bool AgreesWithWide( const SGuardCase& theCase )
{
	long		l_saved = glWorkingPrecision;
	long		l_width = (theCase.lPlaces > 0) ? theCase.lPlaces : 36;
	CLZNumber	x, narrow, wideX, wide;

	glWorkingPrecision = l_width;
	x.SetBase( theCase.lBase, true );
	x.Assign( theCase.szArgument );
	theCase.pFunction( &x, &narrow );

	glWorkingPrecision = l_width + 20;
	wideX.SetBase( theCase.lBase, true );
	wideX.Assign( theCase.szArgument );
	theCase.pFunction( &wideX, &wide );

	glWorkingPrecision = l_width;

	long	l_places = (narrow.m_lRadix > l_width) ? narrow.m_lRadix : l_width;

	wide.Round( kDefaultRoundingMode, l_places );

	char*	szNarrow = Render( &narrow );
	char*	szWide = Render( &wide );
	bool	b_agree = (::strcmp( szNarrow, szWide ) == 0);

	::free( szNarrow );
	::free( szWide );

	glWorkingPrecision = l_saved;

	return b_agree;
}

static void CheckGuardCases( const char* szLabel, const SGuardCase* pCases, long lCount )
{
	long	l_differ = 0;
	char	szDetail[ 96 ] = "";

	for ( long i = 0; i < lCount; i++ )
	{
		if ( !AgreesWithWide( pCases[ i ] ) )
		{
			if ( l_differ == 0 )
				::snprintf( szDetail, 96, "%s of %.20s in base %ld", pCases[ i ].szName,
							pCases[ i ].szArgument, pCases[ i ].lBase );
			l_differ++;
		}
	}

	char	szSummary[ 128 ];

	if ( l_differ == 0 )
		::snprintf( szSummary, 128, "all %ld agree with the answer worked wider", lCount );
	else
		::snprintf( szSummary, 128, "%ld of %ld differ, first %s", l_differ, lCount, szDetail );

	CheckTrue( szLabel, l_differ == 0, szSummary );
}

static void TestGuardsInDigits()
{
	Banner( "Guards counted in digits" );

	/*
		Each of these missed correct rounding at the default precision
		until 2026-09-26, drawn from a census of fourteen functions at
		50 random arguments each in bases two to nine. The arcsines were
		the worst, 4.9 units in base two.
	*/

	const SGuardCase	smallBases[] =
	{
		{ "exp", lz_exp_e, 2, "1011.010010001000111000111011011010110111" },
		{ "exp", lz_exp_e, 2, "0.01" },
		{ "exp", lz_exp_e, 3, "22.22222122" },
		{ "exp", lz_exp_e, 9, "7.325362084227275716008018001060288842" },
		{ "raise_e", lz_raise_e, 3, "121" },
		{ "raise_e", lz_raise_e, 5, "100" },
		{ "raise_e", lz_raise_e, 8, "41" },
		{ "raise_e", lz_raise_e, 9, "38" },
		{ "raise_e", lz_raise_e, 16, "2C" },
		{ "sinh", lz_hypersine, 3, "0.000012121000020101101012210021210021" },
		{ "sinh", lz_hypersine, 3, "0.000121021000220210010001110122210012" },
		{ "cosh", lz_hypercosine, 2, "-10111.10001000" },
		{ "cosh", lz_hypercosine, 2, "-11.100110" },
		{ "tanh", lz_hypertangent, 2, "0.00000001" },
		{ "sin", lz_sine, 2, "-11000.00100010" },
		{ "sin", lz_sine, 2, "10001.000110100110110011101111000101100010" },
		{ "cos", lz_cosine, 2, "-1.0" },
		{ "tan", lz_tangent, 3, "212.102120222221221201120210000012000200" },
		{ "tan", lz_tangent, 8, "7.7153251" },
		{ "sec", lz_secant, 2, "0.111100001111001010000000111100001100" },
		{ "sec", lz_secant, 2, "1.010001011011101110001000000100011010" },
		{ "arcsin", lz_arcsine, 2, "0.1111" },
		{ "arcsin", lz_arcsine, 2, "0.11110111" },
		{ "arcsin", lz_arcsine, 3, "0.22102210" },
		{ "arcsin", lz_arcsine, 3, "0.11100011" },
		{ "arcsin", lz_arcsine, 5, "0.4043002" },
		{ "arcsin", lz_arcsine, 5, "-0.404033" },
		{ "arcsin", lz_arcsine, 7, "0.61126526" },
		{ "arcsin", lz_arcsine, 7, "0.64205" },
		{ "arcsin", lz_arcsine, 8, "0.513103345053651376712425202510534415" },
		{ "arcsin", lz_arcsine, 8, "0.707103455351162171624111001447175536" },
		{ "arcsin", lz_arcsine, 9, "-0.8361" },
		{ "arcsin", lz_arcsine, 9, "0.661545134456302865757553506016514370" },
	};

	CheckGuardCases( "guards in digits, bases 2 to 9", smallBases,
					 static_cast<long>( sizeof( smallBases ) / sizeof( smallBases[ 0 ] ) ) );

	/*
		exp of an answer with two to four integer places, which spent
		the four-digit guard on its own size. Base ten: 40 of 300 random
		arguments between 6.9 and 9.2 missed, up to 1.8 units.
	*/

	const SGuardCase	expBand[] =
	{
		{ "exp", lz_exp_e, 10, "9.832109739754108448675907185048734728" },
		{ "exp", lz_exp_e, 10, "9.676942" },
		{ "exp", lz_exp_e, 10, "9.812265" },
	};

	CheckGuardCases( "exp, answers of two to four places", expBand,
					 static_cast<long>( sizeof( expBand ) / sizeof( expBand[ 0 ] ) ) );

	/*
		e to a power whose answer is about as large as lz_raise_e's
		guard: e^36 missed by 2.27 units in base ten, e^44 (0x2C) by
		1.18 in hex at 36 places and 5.65 at 100.
	*/

	const SGuardCase	raiseBand[] =
	{
		{ "raise_e", lz_raise_e, 10, "35" },
		{ "raise_e", lz_raise_e, 10, "36" },
		{ "raise_e", lz_raise_e, 16, "2B" },
		{ "raise_e", lz_raise_e, 16, "2C" },
		{ "raise_e", lz_raise_e, 2, "1011" },
		{ "raise_e", lz_raise_e, 3, "122" },
	};

	CheckGuardCases( "e^n as large as the guard", raiseBand,
					 static_cast<long>( sizeof( raiseBand ) / sizeof( raiseBand[ 0 ] ) ) );

	/*
		The means take their base from their first value, and raise the
		width before the list has been checked. Two base-two lists that
		missed by just over half a unit...
	*/

	{
		const char*	szGmeanList[] = { "10010010.1", "1010101.0", "111101.110101", "10110001.101001000010111011000010000101111000", "110001.111001001000100000101101110010110101" };
		const char*	szHmeanList[] = { "1110000.010100011001101011011101001010010001", "111110.110110", "100010.010010", "11000111.11" };

		typedef long (*MeanFunction)( CLZNumber**, long, CLZNumber* );

		struct { const char* szName; MeanFunction pMean; const char** pList; long lCount; } means[] =
		{
			{ "gmean", lz_gmean, szGmeanList, static_cast<long>( sizeof( szGmeanList ) / sizeof( szGmeanList[ 0 ] ) ) },
			{ "hmean", lz_hmean, szHmeanList, static_cast<long>( sizeof( szHmeanList ) / sizeof( szHmeanList[ 0 ] ) ) },
		};

		for ( auto& theMean : means )
		{
			long		l_saved = glWorkingPrecision;
			CLZNumber	values[ 8 ];
			CLZNumber*	p_values[ 8 ];
			CLZNumber	narrow, wide;

			for ( long lWidth : { 56L, 36L } )
			{
				glWorkingPrecision = lWidth;

				for ( long i = 0; i < theMean.lCount; i++ )
				{
					values[ i ].SetBase( 2, true );
					values[ i ].Assign( theMean.pList[ i ] );
					p_values[ i ] = &values[ i ];
				}

				theMean.pMean( p_values, theMean.lCount, (lWidth == 36) ? &narrow : &wide );
			}

			long	l_places = (narrow.m_lRadix > 36) ? narrow.m_lRadix : 36;

			wide.Round( kDefaultRoundingMode, l_places );

			char*	szNarrow = Render( &narrow );
			char*	szWide = Render( &wide );
			char	szLabel[ 64 ];

			::snprintf( szLabel, 64, "%s in base 2 against its wider answer", theMean.szName );
			CheckTrue( szLabel, ::strcmp( szNarrow, szWide ) == 0, szNarrow );

			::free( szNarrow );
			::free( szWide );

			glWorkingPrecision = l_saved;
		}
	}
}

/////////////////////////////////////////////////////////////////
//	The arctangent and the two complements. The arctangent had no
//	guard and rounds at every step; arccos and arccot are pi/2 less
//	an inverse and rounded twice.
/////////////////////////////////////////////////////////////////

static void TestArcTangentFamily()
{
	Banner( "The arctangent and the complements" );

	/*
		One argument per base and region that missed correct rounding at
		the default until 2026-09-26, drawn from a census of 100 random
		arguments per region in seven bases -- 1,872 misses between the
		three, none now. Checked against the same function worked twenty
		places wider, as the guards are: see AgreesWithWide.
	*/

	const SGuardCase	arctangents[] =
	{
		{ "arctan", lz_arctangent, 2, "-0.11110010" },
		{ "arctan", lz_arctangent, 2, "1.100110000100110011111000101000001111" },
		{ "arctan", lz_arctangent, 2, "11.1111" },
		{ "arctan", lz_arctangent, 2, "10111101111101010010011000111.101011101110101101110010010000011110" },
		{ "arctan", lz_arctangent, 3, "-0.102010210200021201011100002221110210" },
		{ "arctan", lz_arctangent, 3, "1.110210100000211022100021220002022021" },
		{ "arctan", lz_arctangent, 3, "102.12010002" },
		{ "arctan", lz_arctangent, 3, "1200112100012202021.101021200110201211221212212100221020" },
		{ "arctan", lz_arctangent, 7, "0.655120606413536356623266342400042261" },
		{ "arctan", lz_arctangent, 7, "1.643240654330650033031124121522160401" },
		{ "arctan", lz_arctangent, 7, "-10.646401160122016101531431033004250164" },
		{ "arctan", lz_arctangent, 7, "25602026032.50" },
		{ "arctan", lz_arctangent, 8, "0.6206" },
		{ "arctan", lz_arctangent, 8, "-1.7135153" },
		{ "arctan", lz_arctangent, 8, "27.767" },
		{ "arctan", lz_arctangent, 8, "5511100636.300361512415721711214072574141733714" },
		{ "arctan", lz_arctangent, 10, "0.503048383528959224165560052318757228" },
		{ "arctan", lz_arctangent, 10, "1.4926" },
		{ "arctan", lz_arctangent, 10, "3.229877373824385162180003678005931669" },
		{ "arctan", lz_arctangent, 10, "431486141.309717793001131020086436345575485565" },
		{ "arctan", lz_arctangent, 12, "-0.87D4007858D2E685877D468026351D122513" },
		{ "arctan", lz_arctangent, 12, "1.9657EE54E95495D621D6007854E498E2809E" },
		{ "arctan", lz_arctangent, 12, "17.80107D267225145106E52514D4E700E405D7" },
		{ "arctan", lz_arctangent, 12, "61DD0699.E9D" },
		{ "arctan", lz_arctangent, 16, "0.8F938B105A9312B1AA534D627B7D3924CE7E" },
		{ "arctan", lz_arctangent, 16, "1.2A1C48FC0C92548FDAA280B393B3B3B42BC0" },
		{ "arctan", lz_arctangent, 16, "-3.323D0D0E" },
		{ "arctan", lz_arctangent, 16, "305179A0.B73E" },
	};

	const SGuardCase	arccotangents[] =
	{
		{ "arccot", lz_arccotangent, 2, "-0.11110010" },
		{ "arccot", lz_arccotangent, 2, "-1.101111111101010011101000101001101000" },
		{ "arccot", lz_arccotangent, 2, "-10111.100101110100101001111010010100000010" },
		{ "arccot", lz_arccotangent, 2, "-1111110110000011001010101101.111010100100010011101011110000001101" },
		{ "arccot", lz_arccotangent, 3, "0.121200122010120111102000100212022212" },
		{ "arccot", lz_arccotangent, 3, "-1.11202" },
		{ "arccot", lz_arccotangent, 3, "-10.1" },
		{ "arccot", lz_arccotangent, 3, "-2002020022111100100.212122111201201001121022211002000201" },
		{ "arccot", lz_arccotangent, 7, "0.631562526441551310316564000640310230" },
		{ "arccot", lz_arccotangent, 7, "-1.6163" },
		{ "arccot", lz_arccotangent, 7, "-10.646401160122016101531431033004250164" },
		{ "arccot", lz_arccotangent, 7, "-26210041065.255332656034501340525260402342442041" },
		{ "arccot", lz_arccotangent, 8, "0.6206" },
		{ "arccot", lz_arccotangent, 8, "-1.7135153" },
		{ "arccot", lz_arccotangent, 8, "-11.760551254214723626530641004661274245" },
		{ "arccot", lz_arccotangent, 8, "-5671561352.452773574772233036617461151255730151" },
		{ "arccot", lz_arccotangent, 10, "0.503048383528959224165560052318757228" },
		{ "arccot", lz_arccotangent, 10, "-1.416" },
		{ "arccot", lz_arccotangent, 10, "-10.289815045199155982079200083458154223" },
		{ "arccot", lz_arccotangent, 10, "-327882872.752239675229445133946919570093829586" },
		{ "arccot", lz_arccotangent, 12, "0.95D31430794D863390D940367640109901E6" },
		{ "arccot", lz_arccotangent, 12, "-1.73E35E59E99743527DD0D59672DED10357D3" },
		{ "arccot", lz_arccotangent, 12, "-1E.4D129740183E806606659E17EE4E09859293" },
		{ "arccot", lz_arccotangent, 12, "-1ED1E9749.4EE0D2E646657234782283383D83DE448926" },
		{ "arccot", lz_arccotangent, 16, "-0.F37" },
		{ "arccot", lz_arccotangent, 16, "-1.BBC8FA1629FFDB8B8CC1B7681AC137785B34" },
		{ "arccot", lz_arccotangent, 16, "-3.323D0D0E" },
		{ "arccot", lz_arccotangent, 16, "-3B2F6935.72C76A18" },
	};

	const SGuardCase	arccosines[] =
	{
		{ "arccos", lz_arccosine, 2, "0.0000000000000000100110000011101010011000001000111011" },
		{ "arccos", lz_arccosine, 2, "0.110101" },
		{ "arccos", lz_arccosine, 3, "0.00000000212210021222122010221212021222021210" },
		{ "arccos", lz_arccosine, 3, "0.102002" },
		{ "arccos", lz_arccosine, 7, "0.000000000021236" },
		{ "arccos", lz_arccosine, 7, "-0.602664424242124461206602660643325240" },
		{ "arccos", lz_arccosine, 8, "0.000000000000000000354776571156245402354735111631773715" },
		{ "arccos", lz_arccosine, 8, "-0.610567254410476256732762472163301453" },
		{ "arccos", lz_arccosine, 10, "-0.00000000785174304135877134961111785130966511" },
		{ "arccos", lz_arccosine, 10, "0.95845" },
		{ "arccos", lz_arccosine, 12, "0.000000000000000000000006E90314059D815E29427585870E996853371" },
		{ "arccos", lz_arccosine, 12, "0.932EE7723475388645E503765340601755ED" },
		{ "arccos", lz_arccosine, 16, "-0.0000000000000000000000000000F1C6CC2083DDD1BDC7B46ECDD1C61E5C4C16" },
		{ "arccos", lz_arccosine, 16, "-0.16AB2914D513245CCC49B899B946803409E2" },
	};

	CheckGuardCases( "arctan against its wider answer", arctangents,
					 static_cast<long>( sizeof( arctangents ) / sizeof( arctangents[ 0 ] ) ) );
	CheckGuardCases( "arccot against its wider answer", arccotangents,
					 static_cast<long>( sizeof( arccotangents ) / sizeof( arccotangents[ 0 ] ) ) );
	CheckGuardCases( "arccos against its wider answer", arccosines,
					 static_cast<long>( sizeof( arccosines ) / sizeof( arccosines[ 0 ] ) ) );

	/*
		Rarer cases, each of which a mutation got wrong and the census
		above did not reach: the guard counted in places rather than
		digits (base two and three), pi/2 converted at the narrow width
		(base twelve), and a guard of three rather than four, which in
		base two misses at 40 and 48 places.
	*/

	const SGuardCase	rarer[] =
	{
		{ "arctan", lz_arctangent, 2, "11101100001111001011011001.1", 36 },
		{ "arctan", lz_arctangent, 2, "1000011100110111111101110100.101010110011101011111110100000011101", 36 },
		{ "arctan", lz_arctangent, 3, "-2020002210211021200.2", 36 },
		{ "arctan", lz_arctangent, 3, "220.220202211102010010122221222120100020", 36 },
		{ "arctan", lz_arctangent, 12, "2.2571666996D0E456DD7D409D0613722ED9E6", 36 },
		{ "arctan", lz_arctangent, 12, "14.382700E99D230584E458892D2D80779DE170", 36 },
		{ "arccos", lz_arccosine, 12, "-0.143", 36 },
		{ "arctan", lz_arctangent, 2, "11110.10", 40 },
		{ "arctan", lz_arctangent, 2, "1.100110101010000001001110110100000010010000010101", 48 },
		{ "arccot", lz_arccotangent, 2, "-1.1100101111111010010000000001101001000011", 40 },
	};

	CheckGuardCases( "arctan family, rarer cases", rarer,
					 static_cast<long>( sizeof( rarer ) / sizeof( rarer[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	The logarithm in the operand's own base, and the guard of the
//	functions widened for their answer's size.
/////////////////////////////////////////////////////////////////

static void TestLogInItsBase()
{
	Banner( "The logarithm in its own base" );

	/*
		lz_log_b divided ln x, already rounded, by ln b and rounded
		again: one argument per base and region that missed until
		2026-09-26, of 600 misses in 3,500.
	*/

	const SGuardCase	logs[] =
	{
		{ "log_b", lz_log_b, 2, "0.000000000000000000000000001011000010011000011000001000111101100" },
		{ "log_b", lz_log_b, 2, "0.11011" },
		{ "log_b", lz_log_b, 2, "1.0000000001111101001101100101001000101011111101" },
		{ "log_b", lz_log_b, 2, "100.001000110110110010100100110100101001" },
		{ "log_b", lz_log_b, 2, "100010110110010111011001100011101001010.00010110" },
		{ "log_b", lz_log_b, 3, "0.00000000000000002121001221121001110012222221111200201" },
		{ "log_b", lz_log_b, 3, "0.222102022211111212121212200200110101" },
		{ "log_b", lz_log_b, 3, "0.22222222222221111021021100211111011022100221200120" },
		{ "log_b", lz_log_b, 3, "21.12" },
		{ "log_b", lz_log_b, 3, "2000111000121122112011001.20100221" },
		{ "log_b", lz_log_b, 7, "0.000000000000000000000000000000123" },
		{ "log_b", lz_log_b, 7, "0.544333530432542500043444641013606516" },
		{ "log_b", lz_log_b, 7, "1.000000000000065532" },
		{ "log_b", lz_log_b, 7, "20.315341064311461013065200400051245130" },
		{ "log_b", lz_log_b, 7, "22043144136513.04515" },
		{ "log_b", lz_log_b, 8, "0.00000000004626317504436637346740612552207113552" },
		{ "log_b", lz_log_b, 8, "0.2" },
		{ "log_b", lz_log_b, 8, "0.7777777777772160520402353620061623673256552672600" },
		{ "log_b", lz_log_b, 8, "30.750746047360266115021320365312461363" },
		{ "log_b", lz_log_b, 8, "10501262642155.453132050511747074412071703360174652" },
		{ "log_b", lz_log_b, 10, "0.0000000000000000000000000000000057012848" },
		{ "log_b", lz_log_b, 10, "0.9" },
		{ "log_b", lz_log_b, 10, "0.9918" },
		{ "log_b", lz_log_b, 10, "15.0102043" },
		{ "log_b", lz_log_b, 10, "372127430891.750566118891203181405011722267275723" },
		{ "log_b", lz_log_b, 12, "0.00000000000000006E8055950DE00575094677153807E94877117" },
		{ "log_b", lz_log_b, 12, "0.84607D02055867752E73379086480427804E" },
		{ "log_b", lz_log_b, 12, "1.00000018EE64D758D888172646083E8486D21514833" },
		{ "log_b", lz_log_b, 12, "9.886264" },
		{ "log_b", lz_log_b, 12, "5200DE3804D.46D3E0DE" },
		{ "log_b", lz_log_b, 16, "0.000000000000000000927DB3" },
		{ "log_b", lz_log_b, 16, "0.4A" },
		{ "log_b", lz_log_b, 16, "1.0000000000008CF91058E7B3FC3F0369E098B8BFEF806AA14" },
		{ "log_b", lz_log_b, 16, "76.5" },
		{ "log_b", lz_log_b, 16, "B81D361DA8.3" },
	};

	CheckGuardCases( "log_b against its wider answer", logs,
					 static_cast<long>( sizeof( logs ) / sizeof( logs[ 0 ] ) ) );

	/*
		The guard of the functions widened for their answer's size was
		the one still counted in places: four bits in base two. exp_b
		is the one of them taking a single operand.
	*/

	const SGuardCase	baseTwo[] =
	{
		{ "exp_b", lz_exp_b, 2, "1110.111101001111011101111010010100000110" },
		{ "exp_b", lz_exp_b, 2, "100011.100110011000110010001010110100011011" },
		{ "exp_b", lz_exp_b, 2, "1011.110000101111000110000110111111001101" },
		{ "exp_b", lz_exp_b, 2, "1.111101100101111111001101000000100111" },
	};

	CheckGuardCases( "exp_b in base 2 against its wider answer", baseTwo,
					 static_cast<long>( sizeof( baseTwo ) / sizeof( baseTwo[ 0 ] ) ) );

	/*
		And the spelling it always had. Worked wider, a quotient landing
		exactly on a whole number came back bare -- "2" for log10(100),
		beside log10(1000)'s full noughts.
	*/

	struct { const char* szArgument; const char* szSpelling; } spellings[] =
	{
		{ "100", "2.000000000000000000000000000000000000" },
		{ "0.01", "-2.000000000000000000000000000000000000" },
		{ "1", "0" },
	};

	for ( auto& theSpelling : spellings )
	{
		CLZNumber	x, answer;
		char		szLabel[ 48 ];

		x.Assign( theSpelling.szArgument );
		lz_log_b( &x, &answer );

		::snprintf( szLabel, 48, "log10(%s) keeps its places", theSpelling.szArgument );
		CheckRendering( szLabel, &answer, theSpelling.szSpelling );
	}
}

/////////////////////////////////////////////////////////////////
//	The arcsine: cos t near one, and its guard.
/////////////////////////////////////////////////////////////////

static void TestArcSine()
{
	Banner( "The arcsine" );

	/*
		Near one, cos t = sqrt(1 - x^2) came out five digits short: the
		square rounded before the subtraction, and the root rounded its
		tiny argument to the working places. arcsin(.99..9) was out by up
		to 81,000 units in base ten and a million in hex.
	*/

	const SGuardCase	nearOne[] =
	{
		{ "arcsin", lz_arcsine, 2, "0.111111111111000011" },
		{ "arcsin", lz_arcsine, 3, "-0.222222222222222122212210121020200222122021000200002" },
		{ "arcsin", lz_arcsine, 7, "-0.666666666666666425436511251252523155540635631506544" },
		{ "arcsin", lz_arcsine, 8, "-0.777777777777777271172205166145630014533134540275377" },
		{ "arcsin", lz_arcsine, 10, "-0.999999999999999851800775673945127403076357121351855" },
		{ "arcsin", lz_arcsine, 12, "-0.EEEEEEEEEEEEEEEDD17E5D11D10792D7E56D739E8E32429E0E9" },
		{ "arcsin", lz_arcsine, 16, "-0.FFFFFFFFFFFFFFF8D4716D7" },
	};

	/*
		Elsewhere a guard of two digits could not hold eight halvings'
		doubling: 2 to 11% of arguments missed by just over half a unit.
		One per base and region that missed until 2026-09-26.
	*/

	const SGuardCase	elsewhere[] =
	{
		{ "arcsin", lz_arcsine, 2, "-0.110010110010000101000010000010110111" },
		{ "arcsin", lz_arcsine, 2, "-0.01110110" },
		{ "arcsin", lz_arcsine, 2, "0.101110110100110100101101111100111010" },
		{ "arcsin", lz_arcsine, 3, "0.00000201002002010220121110102001111201002" },
		{ "arcsin", lz_arcsine, 3, "0.102110" },
		{ "arcsin", lz_arcsine, 3, "0.20200202" },
		{ "arcsin", lz_arcsine, 7, "0.315315" },
		{ "arcsin", lz_arcsine, 7, "-0.523145121044340604604266262064532115" },
		{ "arcsin", lz_arcsine, 8, "-0.14716" },
		{ "arcsin", lz_arcsine, 8, "-0.32030" },
		{ "arcsin", lz_arcsine, 8, "-0.620176260771201662554053622144231066" },
		{ "arcsin", lz_arcsine, 10, "-0.295822941691990959425763854038229173" },
		{ "arcsin", lz_arcsine, 10, "-0.56" },
		{ "arcsin", lz_arcsine, 10, "0.744645386814360302810167300620519865" },
		{ "arcsin", lz_arcsine, 12, "0.36045236001629912D6015E30E9721675935" },
		{ "arcsin", lz_arcsine, 12, "0.7DDD" },
		{ "arcsin", lz_arcsine, 12, "0.D0DD" },
		{ "arcsin", lz_arcsine, 16, "0.47FEF8EBE08BA647B16E6082BDB1FBF2C25E" },
		{ "arcsin", lz_arcsine, 16, "0.8556" },
		{ "arcsin", lz_arcsine, 16, "0.BC14E" },
	};

	/*
		The two that a guard of four still missed, in base twelve, and
		five does not.
	*/

	const SGuardCase	guardFour[] =
	{
		{ "arcsin", lz_arcsine, 12, "0.88" },
		{ "arcsin", lz_arcsine, 12, "-0.DD1831E57451187605D73E1666754477E10D" },
	};

	CheckGuardCases( "arcsin near one against its wider answer", nearOne,
					 static_cast<long>( sizeof( nearOne ) / sizeof( nearOne[ 0 ] ) ) );
	CheckGuardCases( "arcsin elsewhere against its wider answer", elsewhere,
					 static_cast<long>( sizeof( elsewhere ) / sizeof( elsewhere[ 0 ] ) ) );
	CheckGuardCases( "arcsin where four guard digits miss", guardFour,
					 static_cast<long>( sizeof( guardFour ) / sizeof( guardFour[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	Near one: the hyperbolic arccosine and arctangent, and the four
//	co-functions that take an inverse on 1/x.
/////////////////////////////////////////////////////////////////

static void TestNearOne()
{
	Banner( "Near one" );

	/*
		arccosh squared, rounded and took one off; arctanh divided down
		to a tiny quotient for x near -1. Up to 240 million units in hex
		and 14 million in base twelve, one argument per base.
	*/

	const SGuardCase	inverses[] =
	{
		{ "arccosh", lz_archypercosine, 2, "1.0000000000000001111001010001101100101000001000111100" },
		{ "arccosh", lz_archypercosine, 3, "1.0000000000000001002220200000210221221122010211122102" },
		{ "arccosh", lz_archypercosine, 7, "1.0000000000000001442643352540324641654244256230523061" },
		{ "arccosh", lz_archypercosine, 8, "1.0000000000000001137556723650715362776447372041530733" },
		{ "arccosh", lz_archypercosine, 10, "1.0000000000000001981203091629396666690923548876758291" },
		{ "arccosh", lz_archypercosine, 12, "1.000000000000000421694640948320551D50E6551493E30ED983" },
		{ "arccosh", lz_archypercosine, 16, "1.0000000000000007F88551BF0568B7893AFA40A701EEDD6771F8" },
		{ "arctanh", lz_archypertangent, 3, "-0.2222222222222221221" },
		{ "arctanh", lz_archypertangent, 7, "-0.666666666666665231354456343163301304400436165513463" },
		{ "arctanh", lz_archypertangent, 8, "-0.7777777777777774726474163523027321146237740440164155" },
		{ "arctanh", lz_archypertangent, 10, "-0.9999999999999952705" },
		{ "arctanh", lz_archypertangent, 12, "-0.EEEEEEEEEEEEEEE2D40811718278690EDE45868D4D03445515E9" },
		{ "arctanh", lz_archypertangent, 16, "-0.FFFFFFFFFFFFFFDA7CB5" },
	};

	/*
		The co-functions round 1/x before the inverse sees it, and near
		one the inverse's slope has no bound: arccoth was out by up to
		ten billion units in hex, arcsech by 1,800, arccsc and arcsec by
		31.
	*/

	const SGuardCase	coFunctions[] =
	{
		{ "arcsech", lz_archypersecant, 3, "0.2222222222222221221" },
		{ "arcsech", lz_archypersecant, 7, "0.666666666666666532306" },
		{ "arcsech", lz_archypersecant, 8, "0.7777777777777772573024525255064740474216074135233607" },
		{ "arcsech", lz_archypersecant, 10, "0.9999999999999990981203091629396666690923548876758291" },
		{ "arcsech", lz_archypersecant, 12, "0.EEEEEEEEEEEEEEE2D40811718278690EDE45868D4D03445515E9" },
		{ "arcsech", lz_archypersecant, 16, "0.FFFFFFFFFFFFFFDFB491534E4E8296FCEEA840DE4038E5A7778" },
		{ "arccoth", lz_archypercotangent, 3, "1.00000000000000010" },
		{ "arccoth", lz_archypercotangent, 7, "1.0000000000000001442643352540324641654244256230523061" },
		{ "arccoth", lz_archypercotangent, 8, "-1.0000000000000001137556723650715362776447372041530733" },
		{ "arccoth", lz_archypercotangent, 10, "1.0000000000000004056296619118377505007750202137718970" },
		{ "arccoth", lz_archypercotangent, 12, "1.0000000000000002099511D2870D9D21916767ED194EE9770629" },
		{ "arccoth", lz_archypercotangent, 16, "1.0000000000000007F88551BF0568B7893AFA40A701EEDD6771F8" },
		{ "arccsc", lz_arccosecant, 7, "-1.0000000000000001126264142404344442110415436351260441" },
		{ "arccsc", lz_arccosecant, 8, "1.000000000000004251041101444667770005604250664425551" },
		{ "arccsc", lz_arccosecant, 10, "1.0000000000000005872618837506988049716525021083118868" },
		{ "arccsc", lz_arccosecant, 12, "-1.0000000000000001258282659305607E11365166639077741167" },
		{ "arccsc", lz_arccosecant, 16, "1.000000000000000DF61B2A50D4427082BC2ED5F59741C84C5B0F" },
		{ "arcsec", lz_arcsecant, 7, "-1.0000000000000001155222560611132356414603232203065255" },
		{ "arcsec", lz_arcsecant, 8, "1.0000000000000042" },
		{ "arcsec", lz_arcsecant, 10, "1.0000000000000005872618837506988049716525021083118868" },
		{ "arcsec", lz_arcsecant, 12, "-1.0000000000000001258282659305607E11365166639077741167" },
		{ "arcsec", lz_arcsecant, 16, "1.000000000000000DF61B2A50D4427082BC2ED5F59741C84C5B0F" },
	};

	CheckGuardCases( "arccosh and arctanh near one against their wider answers", inverses,
					 static_cast<long>( sizeof( inverses ) / sizeof( inverses[ 0 ] ) ) );
	CheckGuardCases( "the co-functions near one against their wider answers", coFunctions,
					 static_cast<long>( sizeof( coFunctions ) / sizeof( coFunctions[ 0 ] ) ) );

	// The one a guard of three still missed, in base two and nowhere
	//	near one, and four does not...
	const SGuardCase	guardThree[] =
	{
		{ "arccosh", lz_archypercosine, 2, "1011100010.000100101111010010110010111000111011" },
	};

	CheckGuardCases( "arccosh where three guard digits miss", guardThree, 1 );
}

/////////////////////////////////////////////////////////////////
//	Pi as wide as it is asked for. The circular functions read pi and
//	its fractions from 120-place text, and past that went wrong in
//	nearly every argument that touched it.
/////////////////////////////////////////////////////////////////

static void TestWidePi()
{
	Banner( "Pi as wide as it is asked for" );

	/*
		One argument per function, base and size that missed correct
		rounding at 150 places -- 600 in base two, where 150 of its places
		are still inside the text's 120 digits -- until 2026-09-26.

		**Against references, not against the answer worked wider.** The
		wider answer shares the constant: with pi capped at 120 places the
		150- and 170-place answers were wrong together, and that check
		passed 48 of these 50 on the old engine. Each reference is this
		engine's answer, scored against Python's Decimal -- pi by Machin,
		the series summed there -- and correctly rounded, all 50.
	*/

	struct { const char* szName; GuardedFunction pFunction; long lBase; const char* szArgument;
			 long lPlaces; const char* szExpected; } wide[] =
	{
		{ "sin", lz_sine, 10, "2.297732", 150,
		  "0.747214407010137507856889048980799215395864917905581691943910532025997832711985840512067495789183713024443090846387282609667374122217129836468335222445" },
		{ "sin", lz_sine, 10, "-943.035147179534819097722924063163478745329286126169265440881242137097749416735220671053975656084833234791006668481500739635127153853189859273232164311231", 150,
		  "-0.528940045107375155005911671794896009154786258402953493798060118476102493182448538041474459517141840664851106293121338273842681447050479240947809180163" },
		{ "sin", lz_sine, 10, "-86414715025116120.59", 150,
		  "0.040904517805528851521956773433901702392916657058266572531202759727899620569518058591034816415703126207380991361661970675271916578890812780373625134475" },
		{ "sin", lz_sine, 16, "-2.5655AEA0E230547DA52E69D89DB86D17DECE87470F45B0D13C43D187167CD611BB8787D1BAE1656AE75AC5D409556647361D553FB291EA6847748FE3329ACC8053805AE52ED2118A3D2FD1", 150,
		  "-0.B86AB820438584F6E9FBA4EA021D85D456D76B7F4F99018FFAD7633297360F039ABB05553B3B6C780C9B68008319964B2319356CAE88214AC879086D10520878B82DC68E4A8F2980C82DF8" },
		{ "sin", lz_sine, 16, "3CA.73", 150,
		  "0.4C5D65CBA4F20D886A54D8447BFD6F2DEB498A7CD7BF63ED5148E5312513F8596162D5145A6152B229AF44F86ABD3756A61C2F0D47ED7A7F0CD4C05E9D76DBF45C4CA133CDD28C2A31A708" },
		{ "sin", lz_sine, 16, "15BCA3F51631D5E.C52553D1D877A6816A096BB35218491C4EE139845077B622914DD952EA708FC3FCACF2C47D1D018F790A003BAB0A431E1450D35D1BC9AE4251EC3884DED517A8E44230F2FFA1C5CC50E799", 150,
		  "0.01E314542D4DB831A5325ADA9E7A50C3DA300B7F08372FD27730ED6DAF25A89519F78A6181C1A505CD1B6F8A7F811B0FABB4DBE2F8244A0F993A1B28E89F15DD3192C4EFCCDEFE4963F32B" },
		{ "cos", lz_cosine, 10, "0.45", 150,
		  "0.900447102352676921668840611486446439757623092728769156133546739322999057351261666042273749693405975415388545896354237942154403660291670464818538005399" },
		{ "cos", lz_cosine, 10, "0.7463387", 150,
		  "0.734179643107853692177407214378193257008532789474852421269789730219410946281482251936021865417273115999975932380906530656642032400997978110613990525343" },
		{ "cos", lz_cosine, 10, "-764.795166570156911953663278827849744550364732382928832508448132190966325549701671904164159232347142534868321325464154237488247415093606755970678975144461", 150,
		  "-0.181630796672715418782203343047396400735130041131449131378189706862767946115600886653623796948265445664063589240697502822375103298810048407198118713437" },
		{ "cos", lz_cosine, 10, "-97315836780409636.863087225225975839650762896777541985109540615631715547491119214851907151128879134811458908068333314942004770126365655674592703756625335157373694009816", 150,
		  "0.113861319528372684275515184760938414111434392445325594740092192474696647385262702431273433571969001717366378374557652449447590527671762434510166596904" },
		{ "cos", lz_cosine, 16, "0.71C9F1", 150,
		  "0.E71FF06080CE2ABBC2D15C83855B53F87B1B4EA125E1AB96126E832E0FB0D7D485B9BDCCFD0B36EDA6914126C66EA89C820A042DB8729FD733AB78F2A7E1588F57043E2A06B39ED8C1CE8B" },
		{ "cos", lz_cosine, 16, "0.B4182ED6791E80EEBDD6C0EE24B89849F015AD266CE3DF33CF4E214C3220A11D45F3A132E8FD372F68F9274B4D0F81BA47E84BDD7B17E0A97494A025C566E811DDF497C7BE1E3E1C186376", 150,
		  "0.C338E0866743D41C9E50D3D46D8D06D8FE3C02613A2B1BC74FF8B715A953AA5BCA05208844D0C07255F63E22D6AD360E829B978E488E5C509337E95BD00E4189838D00AA092CDE4DD64A4C" },
		{ "cos", lz_cosine, 16, "3A3.BC3", 150,
		  "-0.400D082863F68BE2D934C00480EA909628B6A4A1D3338777A5498EB87134922B2F140C99DC80AEB23545C4729119992837CDF8AE305B0991CBE5D4E74A8B18C93CAE2FF269D61B9CB52563" },
		{ "cos", lz_cosine, 16, "-136A3521100F456.69CB9A88AFC81DF42EFDCD58723C1DF7C280CD45207237DE020E91CF82FB98BAC72ADC50646D320063FCB836A0F3374AC373E637492D4B54E3D02EAE2BF64E7BF3D14D38512FB20E2891F3", 150,
		  "0.4FF35E391399317699D6165E33009A5773A64E45867FC738449CE353D490371FF38A5E9E81EC957C1208F90C022D0EFD201C3DC6171E908AFFC4DED2503F5AFFD06E6F25FCDEB4E9E4059F" },
		{ "tan", lz_tangent, 10, "0.476938173369093844955274275466973610434145190421980082689667386218150834872530652184019940644780474665855535334611733450258719493751663058241286099375", 150,
		  "0.516725334824191329668254035920247629366885888632954516231795609464349529879155873635033593693746799359707516368152843802795046993777784190940630222078" },
		{ "tan", lz_tangent, 10, "1.549919", 150,
		  "47.891893344849662629241857985031030420313174310458869635914661252643896062792847515312515528262935808360027890136825625448241759010981260442561160628705" },
		{ "tan", lz_tangent, 10, "-406.843", 150,
		  "148.116051920759968547101425390755673217593511658021538400423170674076844364461034452763867391765450023680908594629833411498474035768635404345957041486668" },
		{ "tan", lz_tangent, 10, "5046228044058244.748331348437891568946477741274809094943776511772343184497197080860628765711587252510280428862489850971023310939338463108599625343858540029257642754671", 150,
		  "-76.666712792694467959896278316139121560472618707350175328888041052866506149351989397204557049069816903919932988031800841601267594261401072631617511040397" },
		{ "tan", lz_tangent, 16, "0.7F54AE7", 150,
		  "0.8AFC59424A20CDC71F93995AF3B0278CBA3BB45AE148AE49CD3EAD38436D0F3A0BED0C4157BF733AE2F311B36DBD6C778076097CB272623F6713E76CA1AD8875C1BD38AEAB385508AF9530" },
		{ "tan", lz_tangent, 16, "1.93874454762A5667BCCF4415976377A800F2CA4DA2DA9A0744475DA877E75AE64CCEF9CF4E6F1C4B478893B812BEBF62E74028F7DBBEB293ABEB3615477C740323C11215F4155F196F5C02", 150,
		  "-B6.4416E2DAFFE3FF257725D68A8833C61824F0D9AB8DA987CBFDE4115A0E3BF1BC706D4630D389B080EBF13A0C920289397206239FEB3E596DEDEC904A631F0B84C7E02A4CF518E0D642FBAD" },
		{ "tan", lz_tangent, 16, "3DF.2D15034A", 150,
		  "-114.444FA855F39874F37EDE81BB72B025AB385DD1AAE4E308C7AFFF478E4265FE4063EDC493B56040BAA3893000A5ED51A672828FE366F7C64C4F49D21A92035C23D4D1F18F8261C5537EF3A5" },
		{ "tan", lz_tangent, 16, "-15C86FD5B5361D0.BD617BA3CA13946D6321C214F8A31D2F06F597FC5BFEC6636F8752B7C4B4F6672A30B86D9C2302D4966CE2878BB12EAA695032AEBAA4D3F44ECE45B59F7DF1B7E018C020FF9F342E7E7E1F", 150,
		  "4D.19FB8D0FF53C07D22CAA5BC79C8457AE30C5C4BE98D4DA0FF4A0240AEDFD446DB758853120D81B84620FF39000502AA38B34E3955A0E8045B6FAAC276B18DBFBB43FA21896D36C767EB5C8" },
		{ "sec", lz_secant, 10, "0.48", 150,
		  "1.127402169187878742804868430935659651183705037119781020696857443051814582838037291243456853926921620527583250603874626027813834749304496719160567915224" },
		{ "sec", lz_secant, 10, "1.6864415", 150,
		  "-8.666444284901299064441045481390041985603748921771248293602760162577193955154548029509592642750489880382040197338469913196627409835767360274608953855110" },
		{ "sec", lz_secant, 10, "92.66804447", 150,
		  "-111.873193033298973688870277737819791225525234613624182347276840364774288701266682420638617923670956839225282247342626110786145485893605742666436106076697" },
		{ "sec", lz_secant, 10, "69611214489681356.146870835303961617622522317573960215620927239552852918610888028619119325760929374022530014079455890061610050501372748914858698225218231062714310358816", 150,
		  "-35.161703069614267801001137761596823413434193278699022110952452558395113749523794516429749619077391094327083181985122599983650714948234226533497116820636" },
		{ "sec", lz_secant, 16, "0.7B16E", 150,
		  "1.20BCEC71F7BE922EC11F3D8DCDC9BB0CA93F310E1FD48E6883E03CC64E14979E2ED4D48D609C4663563F3688B97FB168C86C354E7E6B92204CDD55C7D01D0353286A3D8AC9580340524864" },
		{ "sec", lz_secant, 16, "1.CA4F1049D1BF17ED8441ADD80C5A602018A938DE624CF00986A53755781329428500990A2B335792361FDB0FD161357B785A7D3E231D66DD83EB104DF93EA3666EF7294CE4FF6C69004678", 150,
		  "-4.97D9886B5007E035701890361CF57F6FB764A5C8263D7C432D392F410F9E2F2C17A8EA8C55C3B8536A867C13D797EC3FE0B1B8B46E7DA65F3FB5250976A78D0D60837906006B7994E8129E" },
		{ "sec", lz_secant, 16, "-21D.DE73E04AEC7346A0DDAF02CDECAA4391E2D83C16AEEB9B75A0B30ADFA4C57FD1F2E500189F3F5059037E109FA7D60BFA6D09AE8A6EFBD752C4243EA53D44C46D4191F0992EF7287AD1C6A6", 150,
		  "11.F013A200C7089A4B09DDE446A0F88310849CF2F6E96BA4AB3722171331BDAC6ECD4923F7B2FB8EDB69098D29BDEDB73280C7351164F32697FE81AFCA32DF44C75546C651B19EACA650D03B" },
		{ "sec", lz_secant, 16, "90A80D55A21CC9.8D", 150,
		  "-E.8F65CB69654D09C67A6EACA1113D799B3AFFC7483D208AEC037259F22FD8AFC0E9BFA3DAC6E5EA82ABEA90F8F9D56B2A9B4E13A2F956CF60D8246409A8A2918E997A655E6E0D2FD4716BEA" },
		{ "csc", lz_cosecant, 10, "2.300308577017979527190956395480568143471735866787236736859797303748100777047579762587402612474996331808497089571414695423238010747302503638018013569503", 150,
		  "1.341382380202959159353319223204555995212267468958347225276718023780871171120265473036507830144489477183353540412201792923064440080470641622976803889561" },
		{ "csc", lz_cosecant, 10, "502.612396721058446340587368521571827239992087563671455811731481348348876898968041977742609503385333224209061390743258794029104839114680536788273881553240", 150,
		  "-23.576495284923567660398035729866009977186806595747255285382948586987305514921227598884716045019942146443140718152142214577239084211279017227155210055635" },
		{ "csc", lz_cosecant, 10, "10650525283779856.685872662391156611771333583603807531813819072595632899468193012206440238110075727143728852098816981113317518820782397195716304243622881854680257793775", 150,
		  "-12.790252402461609929340796180304690373421882496025988031468989110608854254276466569838632422500357580969901703261019401034379410999402360250839435243980" },
		{ "csc", lz_cosecant, 16, "2.F582", 150,
		  "5.81F67204BEF48745E4FF099632AAF4246D6651ECFAA118507C16B1F103671B87C416847FD212257DD28FF88E92BCDBE16F8A58EC64ABB306452069281CF4F6B8861675F58BB3D64EC6A4E2" },
		{ "csc", lz_cosecant, 16, "20C.940C36637871B951C7F3FEB1619FE57843A9A9A9E79F520ADE2CE38457911AC6E8985016321BC8FC20A51B4B69FC99752B16CC72AC5494E244A908E841E70C40E3D288538B7E6411A303F1", 150,
		  "E.CA69632915E31268178466038B93351F4771B975D35BD277E3F9DDC660FEB1F6F084A9B82D0D4535809F2C9FE55EE11D0D81935CFBDC57C43BD91C3CAFD44F31A6EC4FC2DFA4C8C807CF96" },
		{ "csc", lz_cosecant, 16, "-EE9C1F651343CA.52613BBCBF9B1A00C3BF8B229CCCABAD8E7FB1C2587E8ACC6F016FB874C1C9191597D4A8E27C750706CC0FA9E15066D7E6FCBA84AE7B98CFF40FAEB2989D3B39E7CBBA63DFA650F258D87E", 150,
		  "-20.D7114B842DA18335EDFED1EF1CFAE202148B4B88050872AD677515E227DCE97E76D0020A4C4E87FF92D0423DB7034431BEBB3FFEB306D44A416538A8391F6EBE0BDAD70D45F49D06B4268F" },
		{ "cot", lz_cotangent, 10, "0.476110470853765349962168148466216237101186696882384285776548893006932489318066088085359006300877541577376297092505441504073003281633945663869336056444", 150,
		  "1.939198074715399900728894200699747214492931796656777120915053070075439235947906425820584655320110612931085249484584688018275249885371570506795807567075" },
		{ "cot", lz_cotangent, 10, "1.630205831880064966194642250718059494440606481747007897904073798284406506481137008496883077806946113856097386258883530248140916553457408915330474014828", 150,
		  "-0.0594794989746320326719644812150239416139792501294625814411938722972996168335047345749109209137226460583270220948865879224303832179180710331549709387651" },
		{ "cot", lz_cotangent, 10, "629.89", 150,
		  "-0.000672955348042510145371002043130889954788771292947554357089762602985658112020100976126782414459379589855258107735450994158680369326363977887967700291931" },
		{ "cot", lz_cotangent, 10, "-59751558011429842.997", 150,
		  "50.692083977560148058672854153523371124641860589384008029594469296171037589553040860605286738869508113460766291631910453071395247060350673419734636054187" },
		{ "cot", lz_cotangent, 16, "0.7856", 150,
		  "1.F7E4A5380211FF32636CE3D54A35B580CBB043DF50A48BAD130B23CC06F94773374FEE3D956E906EA66CB27E13330304BC6184075703029AACA289438DA44861D3FF98FE9A2DC0C8A9955A" },
		{ "cot", lz_cotangent, 16, "2.F88E13624B456EF1ED816AC4B973714F9E56804329DC2FE611FF608B33F092FC7C20970A4443E8D09CD7204E9FCB19BCBC6D1C34C07F94AFA60463BA6C3204C84103B9EDCD5DDC832FC6DD", 150,
		  "-5.CD561AFC79EF9857A68B743650083F4A96DD4CE6ED4700CA50A5B1B51FBC099161B58825286AE078EDCF740FE741141130F2EF0E04F8004983EC9E929D6F25F1CA377E9FA9A1E77449CB66" },
		{ "cot", lz_cotangent, 16, "-2BC.AA8F", 150,
		  "-A.F2DAA64D7C52C66A3D72A1A479E173A66515771B5D700371B68E7BE897DC75FE7E52462C7A755F3C64D7E012FF5F8EBBEAF37FA9969CED372FC11C8EB3D892FDE1729B1939D5C2929A1869" },
		{ "cot", lz_cotangent, 16, "-129DF3F1AE43B9F.01", 150,
		  "-5.C4005C96710C526F88311D740044E0CFA22CA8F5EA057E9F1C33667B1E5C1AE7A9278A944C3B879099DBA6CA0F02284191AA5477B2283EA2D451C92DA15CB35BEBA2F63BF3449C19289626" },
		{ "sin", lz_sine, 2, "10.010010101001010101111011011101010011100100010111110100000011111011011101111111100110101111010011010101011100111110111001100111000111011100011010101100010111011010110000111000110100011000010010101100101001101010011111011110100010111011101001011111000001110011101111110100111000100111010111100111000011111100111100111110000111101101001101100110111100100110101111011111000011001100010110000001001010000000101000011100111101011110101110011001100110100010001001100010111011101000000000000111110011100011100010001110110011001001101101111101110011010001101111000100100110101101100010111010101110001001011111", 600,
		  "0.110000000101111010110000111100100000010000101110000110010100111011000111011011010011101011001000001110011100100001010111111010011010110011010010010100101100111010110111011011001011111010110010000101100001100001001111011110001100110000111101101001000101001100001101010100000110011100111101110100010001000100111111110101100000000001011001101000100100011100100101101000110101101010111110110001101110001010101010000101000111010111111011010011100011111101000011011101100111010100111011001000011000000110101101010000110110110001000001000001100010100010001100101101011111110010011011110101011101100010011011" },
		{ "sin", lz_sine, 2, "101000010111001111010100011011101110001110001000111010011.00000100", 600,
		  "0.011011110001111011001010000011001010101111000010011100001110011010000101001110110011100101100110001011001001011010100101001011111011001110101100100011111001101110111000001101110000000010000001101000110010101111100000011101011110111011100100001010101011000010010100001011110111111101110101100000000101101101100101110010101010101100010101111001000000011101010100001000101001100001101100001010011111011001101010101100100010101111101100100100010001010100100100011000001101000011001001001110000101010000000110010100000111010100011101100010000100000100011110010000011001010000000100110111100101111111011001" },
		{ "cos", lz_cosine, 2, "0.10110", 600,
		  "0.110001011101100010000010110100101110111001001000000000110000110001111100000001111101001010001110100110000001111000110100100000000100111110000010111011010100110011111001001101100101010111010010001101100101001110001001101101110001011011011110011010101101010001000110011101101010000111001100010111011010000011011101110111011110100110001001011000110100000100101001110101010011010100001010000000100011001101101110000011111000110110111110001011010001000100010010001010110110101101100000110100011011000001101111100010110010101101110001100100001010011001001001100111001111000110011101000110101101001100100111" },
		{ "cos", lz_cosine, 2, "100001100101100000000111100011110001111101011011011010111.010101100001001011011000010100011011110000110100001001101000010010000110001011011011011000010110011000001011011010110111011001110111111110011001000100011011000010100010011111000001000111001001001001110111101001111001111100110011110011111011000101110100011111100111000101010001110001011001010100000110000100111110101000101010010001101000111001101111000010100111011111100001001001001001010001101011101101011110100001110000100001110011001001101100100010000111001111101000101011010001001100011110101100110111110101001101010001001111001100000011011000010111110011110100011000111101100000100011110011110101", 600,
		  "-0.010101011111111101011001110111011111000000001111011011111011000101001000011000100111011101110001110011111010101001011000101011110110101111101101101011010011000101010101000011011000110110100001010111000111101110010101100001010100010111111100001100000101000000001111000110110010100011110101011011011111100110001101100001111111111101111100110100110101101000111011011110100011010011110011111101110010010001010001101101011000001110010111100010001010100010111010111111111010100010101000001011000101100110110101111111100101001001100100101110100010000010110010010101110110101000001011000110000010001000000000" },
		{ "tan", lz_tangent, 2, "-1.100010100000001101110100110111010111010011110000010111111010100111111110001100111111111011101111010001000101101101001010101000111001001001100001000110010100100110100011101100000001011110010111101000000001010000101010111111110000100110100111100100000100010011100010101011010010001001101101011010011011001000010001110101100011000010111111010101000100000101101111110001101101000101110011000000001110010011011001110110011100110001100111011100100110101101001100111000001001101110001100110111100100100001001111111110001101011100010100001100011011111101011100001100010111110100111001111101001010100010111111", 600,
		  "-11111.100011011101001111101000101010101101000111000010111011111100000000111001110100110010111000111101110100001100010100100000111101011001111001100100010110100111101111111100101000001010100100011001011111001101110010101001100010010110100001101100001111100011001100010001001001100000010100011111110100110101110111110001010100100101111100011000110011011101110100111101011110111110000001011101101001010001100000100000101010000011111000111011100010111111000110100001000011001010011100010000100011111110010101101010011100101001001100001000010001010000011100111100000110000100000101010000010010100010101101111111" },
		{ "tan", lz_tangent, 2, "11110111111010100000000100000011011010010001101111100011.000001100010111100110010110011100111101110011100000110111000011101111111101110000111000001100011110100000001010111111000010111011110010011100101111000000110100110100010101101010111111010001001011010100000100101100011100110100100000110010001110010000010001100010010101011000101111100011001101110000001001100000010111111010111010011000011111100000011010000111100110110110110101101011011000101010100000011110111011011011000110110011010111000101011001100010101001101111001001000110011100101011110101111000001100111001101001100000111111010110001000110010111000101110001110000100111000000010100101011100101", 600,
		  "-111110010.010110011100101101111101100010100101100101010100011101100010111110011011111111000010011110110000001110100110011111000110000111101001110111100100000110111011111001101011101001000001001000011010111100011111011110111011111100010100101111001000111001110111011101100110000001110011000101011010011101111110000000010110101101111110101100110011101010010001011001000100010001000101111111111000001010111101100011001001100011100000100110010000010100001101001001100001111000001010110011001010110001001011000011111110100101110111011011001111000000101011000110101110110010001111111000000111110110011111011011010000" },
	};

	long		l_differ = 0;
	char		szFirst[ 96 ] = "";

	for ( auto& theCase : wide )
	{
		long		l_saved = glWorkingPrecision;
		CLZNumber	x, answer;

		glWorkingPrecision = theCase.lPlaces;

		x.SetBase( theCase.lBase, true );
		x.Assign( theCase.szArgument );
		theCase.pFunction( &x, &answer );

		glWorkingPrecision = l_saved;

		char*	szAnswer = answer.GetStringRep();

		if ( ::strcmp( szAnswer, theCase.szExpected ) != 0 )
		{
			if ( l_differ == 0 )
				::snprintf( szFirst, 96, "%s of %.20s in base %ld", theCase.szName, theCase.szArgument, theCase.lBase );
			l_differ++;
		}

		::free( szAnswer );
	}

	char	szSummary[ 128 ];
	long	l_count = static_cast<long>( sizeof( wide ) / sizeof( wide[ 0 ] ) );

	if ( l_differ == 0 )
		::snprintf( szSummary, 128, "all %ld correctly rounded", l_count );
	else
		::snprintf( szSummary, 128, "%ld of %ld differ, first %s", l_differ, l_count, szFirst );

	CheckTrue( "sin to cot past pi's 120 places", l_differ == 0, szSummary );

	/*
		pi/2 written to 130 places is not a pole: it lies 2e-131 from
		pi/2, and tan has an answer there. Until 2026-09-26 its cosine
		could not be seen past the text's 120 places and tan refused it
		as a division by zero -- a row of TestQuietFailures pinned that.
		The reference is from Python's Decimal, correctly rounded.
	*/

	{
		CLZNumber	x, answer;

		x.Assign( "1.5707963267948966192313216916397514420985846996875529104874722961539082031431044993140174126710585339910740432566411533235469223048" );

		long	lReturned = lz_tangent( &x, &answer );

		CheckTrue( "tan(pi/2 to 130 places) answers", lReturned == kNoError, "" );
		CheckRendering( "tan(pi/2 to 130 places)", &answer,
						"-40471273184271405759923646020865773545476185600234557954835452660635012161400270616700541289280966588532957310318246041173520486227.064766012185465705952230837624127391" );
	}
}

/////////////////////////////////////////////////////////////////
//	The four circular quotients near their poles: a divisor between a
//	hundredth and a tenth asked for exactly the guard and was refused
//	it.
/////////////////////////////////////////////////////////////////

static void TestQuotientsNearPoles()
{
	Banner( "The circular quotients near their poles" );

	/*
		One argument per function and base that missed correct rounding
		at the default until 2026-09-26, by up to 0.88 of a unit -- 3 to
		8% of arguments in base ten whose divisor was that size. The
		fault scales with the width, so the answer worked wider can see
		it, as it could not see pi's cap.
	*/

	const SGuardCase	nearPoles[] =
	{
		{ "tan", lz_tangent, 7, "14.0306650065344264362042111566615235301112614" },
		{ "tan", lz_tangent, 10, "17.2678451913793729806403020779990204376" },
		{ "tan", lz_tangent, 16, "1.93253689520338E206CBDD4F88F305" },
		{ "cot", lz_cotangent, 10, "3.1525070569542819224625626393556490984" },
		{ "cot", lz_cotangent, 16, "C.8DF658B22BED220180DD5EEF3EF4AB" },
		{ "cot", lz_cotangent, 16, "C.923D2BD74BF49598EE4D89036D00CA" },
		{ "sec", lz_secant, 10, "17.2892792134820935306774678853862212777" },
		{ "sec", lz_secant, 16, "7.DBA10D496D0551CF01AF3729D86AFB" },
		{ "csc", lz_cosecant, 10, "12.5768178165513816117632986857707179578" },
		{ "csc", lz_cosecant, 16, "3.2321FFC3FA1E6F880F52959080400E" },
	};

	CheckGuardCases( "tan, cot, sec and csc near their poles against their wider answers", nearPoles,
					 static_cast<long>( sizeof( nearPoles ) / sizeof( nearPoles[ 0 ] ) ) );

	/*
		And where the numerator is small instead -- tan near the
		multiples of pi, cot near the odd multiples of pi/2 -- whose small
		answers the quotient rule gives 36 significant digits. The same
		threshold stood there, and the same half-guard mends it.
	*/

	const SGuardCase	smallNumerators[] =
	{
		{ "tan", lz_tangent, 10, "6.2906154985585537248865701834261278205" },
		{ "tan", lz_tangent, 10, "6.2840423450248163984717658068894198684" },
		{ "cot", lz_cotangent, 10, "14.1382533667776868822145759056122324677" },
		{ "cot", lz_cotangent, 10, "10.9946453211774951946610754366368628588" },
	};

	CheckGuardCases( "tan and cot with a small numerator against their wider answers", smallNumerators,
					 static_cast<long>( sizeof( smallNumerators ) / sizeof( smallNumerators[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	The arcsine family past its table. gInvTrigA's 22 coefficients of 50
//	digits stopped arcsin, and everything built on it, at 50 digits
//	whatever was asked.
/////////////////////////////////////////////////////////////////

struct SReferenceCase
{
	const char*		szName;
	GuardedFunction	pFunction;
	long			lBase;
	const char*		szArgument;
	long			lPlaces;
	const char*		szExpected;
};

// Each answer against a reference scored independently and correctly
//	rounded -- for faults the answer worked wider shares, a constant's
//	width or a table's...
static void CheckReferenceCases( const char* szLabel, const SReferenceCase* pCases, long lCount )
{
	long	l_differ = 0;
	char	szFirst[ 96 ] = "";

	for ( long i = 0; i < lCount; i++ )
	{
		long		l_saved = glWorkingPrecision;
		CLZNumber	x, answer;

		glWorkingPrecision = pCases[ i ].lPlaces;

		x.SetBase( pCases[ i ].lBase, true );
		x.Assign( pCases[ i ].szArgument );
		pCases[ i ].pFunction( &x, &answer );

		glWorkingPrecision = l_saved;

		char*	szAnswer = answer.GetStringRep();

		if ( ::strcmp( szAnswer, pCases[ i ].szExpected ) != 0 )
		{
			if ( l_differ == 0 )
				::snprintf( szFirst, 96, "%s of %.20s in base %ld at %ld", pCases[ i ].szName,
							pCases[ i ].szArgument, pCases[ i ].lBase, pCases[ i ].lPlaces );
			l_differ++;
		}

		::free( szAnswer );
	}

	char	szSummary[ 160 ];

	if ( l_differ == 0 )
		::snprintf( szSummary, 160, "all %ld correctly rounded", lCount );
	else
		::snprintf( szSummary, 160, "%ld of %ld differ, first %s", l_differ, lCount, szFirst );

	CheckTrue( szLabel, l_differ == 0, szSummary );
}

static void TestWideArcSine()
{
	Banner( "The arcsine family past its table" );

	/*
		At 60, 100 or 200 places every argument missed, right to about 50
		digits: 41 hex places, 165 binary. The references are this
		engine's answers, scored against Python's Decimal and correctly
		rounded, all of them. Not against the answer worked wider, which
		met the same table's ceiling and agreed with it. Two sit at the
		seam where the table hands over: hex at 44 places, base ten at 50.
	*/

	const SReferenceCase	wide[] =
	{
		{ "arcsin", lz_arcsine, 10, "0.607292494172907176174988958342313046455898135335684904944857", 60,
		  "0.652648230031562595463055832378118947029193856695675191528610" },
		{ "arcsin", lz_arcsine, 10, "0.6590857630220323817476094223372206835484479061721095669112886162743093429521692848880082478088488875", 100,
		  "0.7196024812696029841900222260890175155518566401713492422027837356041999407323332322665351280829890561" },
		{ "arcsin", lz_arcsine, 16, "0.4D7A6890773A401E5D7A7302A1726695554BF8A040572A05CE2DF792AF79", 60,
		  "0.4EB6689DA3F43ACF554E87043EA46B61884ABF680CD583E3469F2AA3B890" },
		{ "arcsin", lz_arcsine, 2, "0.11101101011011101000100011110011111111010010101101100010100111100101001010001100101101001010100111110000100011000101110001010000101011110011010010000000111010010011101101001111001101110100000110111110", 200,
		  "1.00110000000001010110110100101000101001001000101111000110011101111111010001001111110011111110011000000111010011000001011100100110101001101011000110001101111011101110101010101110000000001011111110001010" },
		{ "arccos", lz_arccosine, 10, "0.752741617060290608022850637501043288999894411491402045584983", 60,
		  "0.718579513956600938124977343352461417996167061932505348492660" },
		{ "arccos", lz_arccosine, 16, "0.FFB96A4C86CAB346530B61B35CC10AB5D2C242EA2297597836F7E0189B8B169637190F6D18B36036F8331AF57449F2E6F042", 100,
		  "0.0BE1EE809615E50218956DC3E736ACD3A01D03F22FD0612251981D3C2735A488B38DD4A71111F9B0D9673D485BB4F1D0E4EF" },
		{ "arctan", lz_arctangent, 10, "9.925453282080987262204867627303806476541158777455074853251698", 60,
		  "1.470384099902832225620811689426303790894305318988131141942296" },
		{ "arctan", lz_arctangent, 10, "1.43423520343432248341507675971027920065123652779599144886585515937426561345042427813347394007498626544389692024041139986589102417275938142049114468000079798292864371686552675563586523770857554148649743", 200,
		  "0.96192800107513766812465849954307658259650472197401903387623973999390207112671099851770630716700447084509524950748730839652375575308326990543340895697260966960299111782162031204847802424997339834301059" },
		{ "arctan", lz_arctangent, 16, "5.A8DF46879D643993EE799C9D9AA7E4F3322B4CF64A3FC126C204F696DA1D", 60,
		  "1.655A7FFC635168648D075404D6F7E9E304211BFD6382246316AB360105D3" },
		{ "arctan", lz_arctangent, 2, "1.1010110100101101001111100111111001000001101010001101100111011111110101111011110011000010001000000011111100101001110011011010110010011010000110000010110110010101111000100010101000101111100001010110101001000111001000111111101011101100101001010010101111001011000100110001110111101101101000000011110011101010100111001001001110000001110011101011110101010110011110110001001111100001111010000011011001000111", 400,
		  "1.0000100001110000001000110000000111011010101100100010111101101110101110100001110100011011110110011110000110101010111011010110011000001110000100010010000000100011010101001010111101101000100100011101100001110100101110001000010011101111100011000010001000001001001001011011100011011100111001100010001101100101100001010111110110111010111000000111000100110101010011100111110101001101001000101000111110111001" },
		{ "arccot", lz_arccotangent, 10, "3.2795357350734381463086906849583034456057106303685536620323760250039394818348912371400074945291449899", 100,
		  "0.2959655272602540339275102847189808624405342595935703332616758384561095667955744110573955929840893892" },
		{ "arccot", lz_arccotangent, 16, "B.5D36907842156D275176B689A90523981FD574BBE92AA366A6243B40A78C", 60,
		  "0.16781BFEDA504532D4E88609DE94985E2B31741400E26A7577AD1E3FE0ED" },
		{ "arcsec", lz_arcsecant, 10, "8.729249417290717617498895834231304645589813533568490494485729", 60,
		  "1.455986885211484646053323468685663793715420325813460402384270" },
		{ "arcsec", lz_arcsecant, 16, "2.CC8A4C714FE471F1BB3AF1D17731D907809FA599775312BE4ADA8666615E95E1C0721E9C38178E4BE4D890702E4C2568DF7F", 100,
		  "1.34986721848434007DCD55461643F3E7809CA07E67331548F613C17B9C2CCD7FF36B36588D99335DF2B7E33EB05B326E3AE4" },
		{ "arccsc", lz_arccosecant, 10, "2.2795357350734381463086906849583034456057106303685536620323760250039394818348912371400074945291449899", 100,
		  "0.4541357399049951924222206499703898640680843775646719014593352553469326703668213716611437327054955552" },
		{ "arccsc", lz_arccosecant, 16, "A.5D36907842156D275176B689A90523981FD574BBE92AA366A6243B40A78C", 60,
		  "0.18BD35E6D2FAD0232E8352FEAC8B2A370C03FF322A2E073AC689777BC97F" },
		{ "arcsin", lz_arcsine, 16, "0.C852B2ADBACC19FF2CE591CD81CD4B27E49FCE9DA25A", 44,
		  "0.E6108AAFEC0C60B63F48299F10E1C5CB2F5614C3791C" },
		{ "arctan", lz_arctangent, 10, "2.46408764639562578689429806846912501292642710269804", 50,
		  "1.18527445338352049452962895287263360491376976808629" },
	};

	CheckReferenceCases( "arcsin to arccsc past 50 digits", wide,
						 static_cast<long>( sizeof( wide ) / sizeof( wide[ 0 ] ) ) );

	// At a thousand places the terms number about five hundred: the loop's
	//	bound has to reach them, where one sized for 1/k! stops near 450...
	const SReferenceCase	thousand[] =
	{
		{ "arcsin", lz_arcsine, 10, "0.3", 1000,
		  "0.3046926540153975079720029612275291669545600317067763873929779487464729925120331594385999572313631487317349566307855469247966854234324188781871838304561501103439852326461226950624544315425851918591328531307809295007959773297208155336494895491479932859480918639385082525831120672069475803857191645095129466908520770724410760283417329837507482191921273750024356409880677646828606391895751707132459946846470805120485599600432068047523100817614119648187675397084261489849327757000628598810903351153553663363052492557072191403251257060203521305736553458020435956232641918895161924866488432149405558400175439802179418741796017045754884634287867172941432649337323727615627842494680571268300053130247490193745678042766331099699284213442415602737425920565382547578882705164561030721844100243732278008846443355988211556625228496783429075895416774229846958449185685874085936521761174916938685833621689896967311351583716337878841290045455152820550116155262708187945221532076084301907680229616553817971657895814596" },
	};

	CheckReferenceCases( "arcsin(0.3) at 1000 places", thousand, 1 );

	// An infinity has no arcsine at any width. Past the table the series
	//	is computed, and there the flags do not spread through it as they
	//	did through the table's multiplications...
	{
		long		l_saved = glWorkingPrecision;
		CLZNumber	x, answer;

		glWorkingPrecision = 100;
		MakeInfinite( x );
		lz_arcsine( &x, &answer );
		glWorkingPrecision = l_saved;

		CheckTrue( "arcsin(inf) at 100 places is not a number", answer.FlagIsOn( kfQNaN ) && answer.FlagIsOff( kfInfinity ), "" );
	}
}

/////////////////////////////////////////////////////////////////
//	Gamma past its Bernoulli table. The fifty-digit entries stopped
//	lnGamma at 52 significant digits in every base, and every large
//	gamma with it, the default included.
/////////////////////////////////////////////////////////////////

static void TestWideGamma()
{
	Banner( "Gamma past its Bernoulli table" );

	/*
		Every one of these missed at 52 to 54 digits. The references are
		this engine's answers scored against an independent Stirling sum
		in Python's Decimal -- exact Bernoulli rationals, as many terms as
		the width needs -- and correctly rounded, all of them. The first
		three are at the default, where lz_gamma works wider for a large
		answer's size. 0.3 at 42 and 43 places sits either side of the
		seam where the table hands over to the rationals, and log-gamma
		of 2.5 at 130 places is past ln(2 pi)/2's 120 places of text --
		as is 2.8 in hex at 100, where those 120 digits are 100 places:
		compared as places, the text served there and missed.
	*/

	const SReferenceCase	wide[] =
	{
		{ "gamma", lz_gamma, 10, "0.00000000000000000001", 36,
		  "99999999999999999999.422784335098467139403378469870877295" },
		{ "gamma", lz_gamma, 16, "14.8", 36,
		  "780AEEF55272EB0.7943933C4315ECE4E9437C633D87B666917C" },
		{ "gamma", lz_gamma, 16, "21.B", 36,
		  "22EF560C613652A470CAC9617F6FF30.E75878C053703FC03ACD009B2F1D5A646200" },
		{ "gamma", lz_gamma, 10, "0.3", 42,
		  "2.991568987687590628312516515904917791112806" },
		{ "gamma", lz_gamma, 10, "0.3", 43,
		  "2.9915689876875906283125165159049177911128060" },
		{ "gamma", lz_gamma, 10, "0.3", 60,
		  "2.991568987687590628312516515904917791112806024921715112744120" },
		{ "gamma", lz_gamma, 10, "7.3", 60,
		  "1271.423633663909273057993626678458337854195378230011887106763149" },
		{ "gamma", lz_gamma, 10, "-3.7", 60,
		  "0.251643995902422643510108134681313060212549189939663622676805" },
		{ "gamma", lz_gamma, 10, "20.5", 100,
		  "540624298233507504.4736873647808221412273166395781714789601948262439262966163731666261810160445254372779635704556154157" },
		{ "gamma", lz_gamma, 16, "2.8", 60,
		  "1.544FA6D47B38FFDB52C43D8A064AEA2FEAAD8FF538AC5BE2578C08DEEF8F" },
		{ "gamma", lz_gamma, 16, "-14.8", 60,
		  "-0.00000000000000053AA0F3C6CCEA6598799E75B9B0998C07255ECC2C6325" },
		{ "gamma", lz_gamma, 16, "7.3", 80,
		  "401.A3C561F1159F45E17D783273C4881E0053A90A3065D684455C923C2732484B9BF4EC21C19DBE9588" },
		{ "gamma", lz_gamma, 2, "0.01", 200,
		  "11.10100000001001111111100010001111110111011101110101000111100110110011100010100101000010100101000101000110101010011000111111111101111000100110001111110111010101000110011011101100000011011011000001001100" },
		{ "loggamma", lz_log_gamma, 10, "2.5", 60,
		  "0.284682870472919159632494669682701924320137695559894729250146" },
		{ "loggamma", lz_log_gamma, 10, "1000.25", 60,
		  "5906.947268271117176996487246964070175850780836352525866043217136" },
		{ "loggamma", lz_log_gamma, 10, "0.3", 100,
		  "1.0957979948180755216771681423701072784451484507642034066386236431982548357308983342519333148704104177" },
		{ "loggamma", lz_log_gamma, 10, "2.5", 130,
		  "0.2846828704729191596324946696827019243201376955598947292501458503867759342216325755537007359586395675549719731391654888545210665184" },
		{ "loggamma", lz_log_gamma, 16, "2.8", 100,
		  "0.48E0FA02699F9E7DC01CDC5A3F7534236B64AA86ACA851CC433DD798B1C7D97C0879F900DA07EF930E67674F247EB4E6C1C4" },

		// Wide enough that the rising product was thousands of digits and
		//	its logarithm quadratic: 1.5, 9.7 and 1.4 s until 2026-09-26,
		//	13, 65 and 11 ms since...
		{ "gamma", lz_gamma, 10, "0.3", 150,
		  "2.991568987687590628312516515904917791112806024921715112744119650956388767876320217990259952557650565421425626826142022375700443136855925650397383325007" },
		{ "loggamma", lz_log_gamma, 10, "0.3", 175,
		  "1.0957979948180755216771681423701072784451484507642034066386236431982548357308983342519333148704104177298666656338017093038633459326902534661479521858483951037770873576350791012" },
		{ "loggamma", lz_log_gamma, 16, "2.8", 120,
		  "0.48E0FA02699F9E7DC01CDC5A3F7534236B64AA86ACA851CC433DD798B1C7D97C0879F900DA07EF930E67674F247EB4E6C1C40CD7D6E843D33C44B1EC" },

		// Past thirty terms, from the tangent numbers: 225 places took
		//	2.7 s with thirty, and 300 or 500 would not have finished. Two of
		//	these need no series to check: gamma(20.5) is 40!/(4^20 20!)
		//	sqrt(pi), and gamma(0.25) comes from the AGM...
		{ "loggamma", lz_log_gamma, 10, "7.3", 225,
		  "7.147892523022249032777057154428389202452954042566484621679559153130698948654857595627356803252059822758824620208899084366934649583717425242060273712559406356850120817559375310848845384587343732492342015244209222839936160405660" },
		{ "loggamma", lz_log_gamma, 10, "0.3", 300,
		  "1.095797994818075521677168142370107278445148450764203406638623643198254835730898334251933314870410417729866665633801709303863345932690253466147952185848395103777087357635079101174850750100229143746809739702705184905342523180968004318730653285671984359799458861546498320291333358869066132132414243474778" },
		{ "gamma", lz_gamma, 10, "20.5", 500,
		  "540624298233507504.47368736478082214122731663957817147896019482624392629661637316662618101604452543727796357045561541567404348350985989714158568398243744853683244918144411441069415672970183851959253682428303707266017898882258538786621489335134749608669309518712588793777951726661335977345983511236677351634051812509184066931328450389362110480851348800815966231204099846376558666228981712570301785338530491024006486236776788786582107664060042513132287803485415801358405838253289236125240706551536819897770064415834097413" },
		{ "gamma", lz_gamma, 10, "0.25", 300,
		  "3.625609908221908311930685155867672002995167682880065467433377999569919243538729121618360136723384300361471751392420719965891524094022559977426458890361450606413744896854194999201926773037994630892212412318323707992084397369907093905620929232342870274191448603957136835036865487995968368476475851489090" },
		{ "gamma", lz_gamma, 16, "2.8", 300,
		  "1.544FA6D47B38FFDB52C43D8A064AEA2FEAAD8FF538AC5BE2578C08DEEF8F57532193FE13AF17564C1A2AD88F16895975913D8F352652823CFDB36C9277B1917A58C3E74C75D46C1BC7FEC26B4C4EFD13C25FFE78CD3C9D2A154234DB9D13C6BEC8A85151A2853A9717506B707AE49220DA97157DA67AAADBC937FC84C30F9C50D275E885E1AE74ECA910308C7070E03BDE1A2EA0E3CC" },
		{ "gamma", lz_gamma, 10, "-3.7", 300,
		  "0.251643995902422643510108134681313060212549189939663622676804506267308381312095306902722886967442279710081983397358873358711690104966809300930963175360822714926264262115639365137828201906076982086156429932106156504779013138928612947650229621107301355061320205840021700292045012759140470700919002599408" },
	};

	CheckReferenceCases( "gamma and log-gamma past 50 digits", wide,
						 static_cast<long>( sizeof( wide ) / sizeof( wide[ 0 ] ) ) );

	/*
		**A large argument, at the default.** lzi_log_gamma multiplies
		ln w by w - 1/2, and an argument past the shift's target is w
		itself, so until 2026-10-07 the rounding of ln w came up by the
		argument's whole places: from about 1e8 in base ten nearly every
		answer was wrong, by 1.9e5 units at 9.9e13 and 2.6e11 at 1e20.
		The same line spent the guard at wide widths, where w is the
		shift's target, about 2,940 at 1000 digits -- but a near-half
		answer there needs an argument as long as the width, 0.17 s a
		call before its retake, so these pin it. The references are
		Python's Decimal, Stirling after a shift to 1000.
	*/

	const SReferenceCase	large[] =
	{
		{ "loggamma", lz_log_gamma, 10, "98765432101234.5678", 36,
		  "3083829015577631.879593351003110818285622882331437916" },
		{ "loggamma", lz_log_gamma, 10, "100000000000000000000.5", 36,
		  "4505170185988091368036.901847901933087943982890327284903040" },
		{ "loggamma", lz_log_gamma, 16, "56BC75E2D63100000.8", 36,
		  "F439CBBEBDD4B35264.E6DF810CC4D1E7BD74FC004D039C1ACEC32E" },
		{ "loggamma", lz_log_gamma, 2, "10000000000000000000000000000000000000000.1", 36,
		  "110101011100111010011101111101011100011001001.010101100100001001100101110011101100" },
		{ "loggamma", lz_log_gamma, 7, "1000000000000000.3", 36,
		  "40121464514120140.414452231333041426265512136652420261" },
	};

	CheckReferenceCases( "log-gamma of a large argument", large,
						 static_cast<long>( sizeof( large ) / sizeof( large[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	The logarithm far from one. It reduced by dividing by e once per
//	unit of its answer, so ln of an n-digit number cost 2.3n divisions
//	of an n-digit number, and a wide gamma's rising product is 9,400
//	digits at 150 places.
/////////////////////////////////////////////////////////////////

static std::string DigitsWithNoughts( const char* szLead, long lNoughts, const char* szTail )
{
	return std::string( szLead ) + std::string( lNoughts, '0' ) + szTail;
}

static std::string SmallWithNoughts( long lNoughts, const char* szDigits )
{
	return "0." + std::string( lNoughts, '0' ) + szDigits;
}

static void TestLogFarFromOne()
{
	Banner( "The logarithm far from one" );

	/*
		**Two independent routes**: moving the radix and adding k ln(base),
		against dividing by e the old way, forced through the threshold
		both ways. Every exponent is at least three digits so that the
		lowest setting reaches it. The rows at 150 and 250 places ask
		for ln(base) past its 120 places of text, and the hex pair asks
		for it at two widths past them.
	*/

	{
		struct { long lBase; long lPlaces; std::string sValue; } logs[] =
		{
			{ 10, 36, "4567.25" }, { 10, 36, "98765432.1" }, { 10, 36, "0.00012" },
			{ 10, 36, DigitsWithNoughts( "3", 20, ".5" ) }, { 10, 36, SmallWithNoughts( 20, "71" ) },
			{ 10, 150, DigitsWithNoughts( "2", 60, ".25" ) }, { 10, 150, "0.000000123" },
			{ 2, 36, DigitsWithNoughts( "1011", 30, ".1" ) }, { 2, 36, "0.0000000000011" },
			{ 3, 36, DigitsWithNoughts( "21", 40, ".1" ) },
			{ 7, 60, DigitsWithNoughts( "6", 30, ".3" ) },
			{ 12, 36, DigitsWithNoughts( "5E", 25, ".D" ) },
			{ 16, 36, DigitsWithNoughts( "3E8", 20, ".8" ) },
			{ 16, 150, SmallWithNoughts( 40, "F1" ) },
			{ 16, 250, DigitsWithNoughts( "1", 200, ".8" ) },
		};

		const long	l_logs = static_cast<long>( sizeof( logs ) / sizeof( logs[ 0 ] ) );

		long		l_differ = 0;
		std::string	s_where;

		for ( auto& theLog : logs )
		{
			std::string	s_answers[ 2 ];

			for ( int iPath = 0; iPath < 2; iPath++ )
			{
				long	l_scaleWas = lzi_set_log_scale_from( (iPath == 0) ? kLogNewtonOutOfReach : 3 );
				long	l_saved = glWorkingPrecision;

				glWorkingPrecision = theLog.lPlaces;

				CLZNumber	x, answer;

				x.SetBase( theLog.lBase, true );
				x.Assign( theLog.sValue.c_str() );
				lz_log_e( &x, &answer );

				glWorkingPrecision = l_saved;
				lzi_set_log_scale_from( l_scaleWas );

				char*	szAnswer = answer.GetStringRep();

				s_answers[ iPath ] = szAnswer;
				::free( szAnswer );
			}

			if ( s_answers[ 0 ] != s_answers[ 1 ] )
			{
				l_differ++;
				s_where = theLog.sValue.substr( 0, 16 );
			}
		}

		char	szDetail[ 72 ];

		if ( l_differ == 0 )
			::snprintf( szDetail, 72, "all %ld identical to the character", l_logs );
		else
			::snprintf( szDetail, 72, "%ld of %ld differ, one at %s", l_differ, l_logs, s_where.c_str() );

		CheckTrue( "ln with the radix moved against ln by dividing by e", l_differ == 0, szDetail );
	}

	/*
		Arguments the division route would take seconds over, against
		references scored in Python's Decimal and correctly rounded: a
		ten-thousand-digit number and a five-thousand-nought fraction at
		the default, a hex one at 250 places, where ln 16 is computed and
		the mantissa goes to Newton, and bases seven and twelve.
	*/

	const std::string	s_big = DigitsWithNoughts( "1", 10000, ".5" );
	const std::string	s_small = SmallWithNoughts( 5000, "3" );
	const std::string	s_hex = DigitsWithNoughts( "1", 1000, ".8" );
	const std::string	s_seven = DigitsWithNoughts( "3", 400, ".5" );
	const std::string	s_twelve = SmallWithNoughts( 300, "5E" );

	// ...and an exponent of five base-seven places, which ln 7 is
	//	multiplied by: ln 7 not taken those places wider is 0.76 of a unit
	//	out here, and 0.25 with them. Found by running a census through
	//	that mistake; the rows above all survive it...
	const std::string	s_sevenWide = DigitsWithNoughts( "6", 11833, ".3443623" );

	const SReferenceCase	far[] =
	{
		{ "ln", lz_log_e, 10, s_big.c_str(), 36,
		  "23025.850929940456840179914546843642076011" },
		{ "ln", lz_log_e, 10, s_small.c_str(), 36,
		  "-11514.129437774554356082580019639582876508" },
		{ "ln", lz_log_e, 16, s_hex.c_str(), 250,
		  "AD4.96B6802E49CD1C3285D6273DB44E7256D8F938BE921E878D0C8CF6664F84EE04A235BAB58268CC393D30967383E6E879F96207EB0D07C8E6B2FC61B013DD45F1F41599D7D27507C2B12F4CA534F2BC18F1C1D1B0F19F9C28F0326EB91E77FAE1DB08373B4E916341C68D4C8DF111AB2FF63BE04061E4D0910863670E45" },
		{ "ln", lz_log_e, 7, s_seven.c_str(), 60,
		  "2162.314460613313061451320535033514013605300636165004610614213435" },
		{ "ln", lz_log_e, 12, s_twelve.c_str(), 150,
		  "-522.21964D4E6E018D36D49511EE48115800736010E83DE6970E72542E101738131D1E6222088764DE2E2688D93944E2271032006998E255812D03952965D5668040D75D45E6095388E41D1626" },
		{ "ln", lz_log_e, 7, s_sevenWide.c_str(), 150,
		  "124064.514032146040030532142626451662506061642220555655506541151333133526603003655433230042644215314044061010241022142400662223545363651254506051202222165535" },
	};

	CheckReferenceCases( "ln of arguments thousands of digits from one", far,
						 static_cast<long>( sizeof( far ) / sizeof( far[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	log_b and exp_b past 120 digits. Both took ln(base) from gLogBaseE's
//	120 places of text at any width.
/////////////////////////////////////////////////////////////////

static void TestLogBaseWide()
{
	Banner( "log_b and exp_b past 120 digits" );

	/*
		Rows from the cells of a census that missed -- log10 at 120 and
		150 places, hex log_b at 100, a large hex or duodecimal power at
		80 whose own size widened the work past the text, base seven at
		150, and base three's power and base twelve's logarithm at 250.
		The old engine missed every one, by 0.69 of a unit to 6e148; each
		is this engine's answer, scored against Python's Decimal and
		correctly rounded. Base two holds 399 places in 120 digits and was
		never short at these widths. Base three's logarithm missed too,
		and is left out for its 44 ms: an odd base multiplies digit by
		digit.
	*/

	const SReferenceCase	wide[] =
	{
		{ "exp_b", lz_exp_b, 16, "14.D", 80,
		  "9837F0518DB8A96F46AD2.3182E42F6F65E139A1B14FA8178D78B65CBEFA7BB6FBE47C34D380250A84E9F94EE0C90DBBD64CD8" },
		{ "exp_b", lz_exp_b, 12, "34.2", 80,
		  "161D741D0891E36ED868E41368D7617142416693E.455975452E00957E50336D1D9DE37E08698056450E27D462094726404EDD085D65316426783454E1" },
		{ "log_b", lz_log_b, 16, "FAE", 100,
		  "2.FE21FEC45FABE38CFAEC1E7F62DE107508853C28EB45E95D13A5306CCEB93BF538950A905201934DBE942D8C84C4B57A1706" },
		{ "log_b", lz_log_b, 10, "715", 120,
		  "2.854306041801080614742966234460859645768241497391412825087302532550197934589019078892014988526444666635768291287514329746" },
		{ "exp_b", lz_exp_b, 10, "8.6", 120,
		  "398107170.553497250770252305087752043487677037297380446865284148060224853869458039074018311607377360253577613955560098735945993156" },
		{ "log_b", lz_log_b, 10, "9", 150,
		  "0.954242509439324874590055806510230618400257728381391729659731280610458305567322246085936711295232603020929385536504091787125938284450454702580687412143" },
		{ "exp_b", lz_exp_b, 10, "8.6", 150,
		  "398107170.553497250770252305087752043487677037297380446865284148060224853869458039074018311607377360253577613955560098735945993155954215388945112255707449439960" },
		{ "log_b", lz_log_b, 7, "30", 150,
		  "1.364435455021564350562526162560323635310663201465654316116126520331600464661633034242655426225036301242304455004145126016042560122145301314036123303233" },
		{ "exp_b", lz_exp_b, 7, "0.1", 150,
		  "1.214630613036165234425442360320401623314216440640255243563524621336615121040662322033163562032005256623542205136044445231122034300451313556566143453034" },
		{ "exp_b", lz_exp_b, 3, "10.2", 250,
		  "2002.0111010212110222122220112020021000112102112221012211220021021002121011002212112201102021011222100000212010222010002212111120102022020020222000022001211222000212021100220111020002122111121112222222010012011011110102021210011002201110101020002100011001" },
		{ "log_b", lz_log_b, 12, "D", 250,
		  "0.E1526972542E2D8E0D2DE077075827072D840821712757633511E9382418E3932087577E769964E6572134D7ED5E8E1D99056E978796823870813427D19693969392419166D05513D847480E617623E1D07D71E97406E59668502579820911053878228166347740E753E471604D4074423EED36D78D887D32D858641D" },
	};

	CheckReferenceCases( "log_b and exp_b past 120 digits", wide,
						 static_cast<long>( sizeof( wide ) / sizeof( wide[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	Gamma below nought. The reflection widened for nothing its answer
//	needed.
/////////////////////////////////////////////////////////////////

static void TestGammaBelowNought()
{
	Banner( "Gamma below nought" );

	/*
		pi / ( sin(pi z) Gamma(1 - z) ), worked at the guard's width
		whatever the answer's size, however small the sine, however
		large z. Near nought and near every pole all three bite: 64 of
		179 arguments missed at the default in base ten, gamma(-2e-20) by
		2e30 units where gamma(2e-20) was right, and 105 of 400 random
		ones. Each row missed on the old engine, from a unit to 2e31, and
		each is this engine's answer scored against Python's Decimal and
		correctly rounded: near nought and near poles, at 36, 60 and 100
		places, in bases two, three, seven, ten, twelve and sixteen.
		gamma(-0.00001) was the audit's last case short of its size.
	*/

	const SReferenceCase	below[] =
	{
		{ "gamma", lz_gamma, 10, "-0.00001", 36,
		  "-100000.577225555552235029678061557345699716" },
		{ "gamma", lz_gamma, 10, "-0.000037", 36,
		  "-27027.604279288242775607086090231060544177" },
		{ "gamma", lz_gamma, 10, "-0.999999", 36,
		  "-1000000.422785746939301928573618097482654424" },
		{ "gamma", lz_gamma, 10, "-1.000001", 36,
		  "999999.577217076741358927277950162647370085" },
		{ "gamma", lz_gamma, 10, "-0.00000000000000000002", 36,
		  "-50000000000000000000.577215664901532860626293209988961882" },
		{ "gamma", lz_gamma, 10, "-0.00001", 60,
		  "-100000.577225555552235029678061557345699716395425863535028613479571" },
		{ "gamma", lz_gamma, 10, "-0.000001", 60,
		  "-1000000.577216653958435668636877440597532732436429983703915790823647" },
		{ "gamma", lz_gamma, 16, "-0.FFFFFF", 36,
		  "-1000000.6C3B9985F0ADE7D3FC6F0FCCEDF794F3EC2C" },
		{ "gamma", lz_gamma, 16, "-1.000001", 36,
		  "FFFFFF.93C4694CEC0E74E1FC81FD8E25BAB7ACED37" },
		{ "gamma", lz_gamma, 16, "-0.FFFFFF", 100,
		  "-1000000.6C3B9985F0ADE7D3FC6F0FCCEDF794F3EC2C5933BDA8C91B764BC92B8711A865A0C159C8E1DDCE5615D5F4DE5B46338107F4" },
		{ "gamma", lz_gamma, 3, "-2.000000000001", 36,
		  "-111111111111.001001010221201111010020011100210002" },
		{ "gamma", lz_gamma, 12, "-0.EEEEEE", 36,
		  "-1000000.50D6D480905024857590619E965313046D32" },
		{ "gamma", lz_gamma, 2, "-0.000000000000001", 36,
		  "-1000000000000000.100100111100011001100010010011001011" },
		{ "gamma", lz_gamma, 7, "-1.666666553", 36,
		  "2600260.343160230543661250030006313154636202" },
		{ "gamma", lz_gamma, 10, "-1.00000000000000000001", 36,
		  "99999999999999999999.577215664901532860620630493386666828" },
		{ "gamma", lz_gamma, 10, "-20.00000000000000000001", 36,
		  "-41.103176233121648583538374884666391549" },

		// Just below a pole: the fraction is nearly one, and its distance
		//	to one has to be taken exactly -- one less 1e-20 is 1.0 as a
		//	double, whose sine is 1e-16 where this one's is 3e-20...
		{ "gamma", lz_gamma, 10, "-7.99999999999999999999", 60,
		  "2480158730158730.158783250036655645089200771789758320248634276618130536735587" },
	};

	CheckReferenceCases( "gamma below nought, near nought and the poles", below,
						 static_cast<long>( sizeof( below ) / sizeof( below[ 0 ] ) ) );
}

/////////////////////////////////////////////////////////////////
//	Main
/////////////////////////////////////////////////////////////////

/*
	**An answer too large to hold is refused, and one too small to show
	is nought, before any of it is worked out**. Every case at the
	default limit below hung
	until 2026-10-11, so each is also a test that it returns at all.

	The other half lowers the limit to a thousand places, where the work
	is cheap whichever way the check goes, so that a check taken out fails
	there rather than hanging: every pair sits clear of the limit by a
	quarter of a place or more, the estimate being good to a dozen figures.
	It runs first.
*/

const int	kLimitRefusedLarge = 0;
const int	kLimitRefusedSmall = 1;
const int	kLimitSpelled = 2;
const int	kLimitAnswered = 3;
const int	kLimitDomain = 4;

struct SLimitCase
{
	const char*	m_szFunction;
	const char*	m_szLeft;
	const char*	m_szRight;		// nullptr for one operand
	int			m_iExpect;
	const char*	m_szSpelled;	// for kLimitSpelled: ZERO, ONE, MONE, or exactly
};

static const SLimitCase	kLimitCases[] =
{
	{ "exp_e", "1E24", nullptr, kLimitRefusedLarge, nullptr },
	{ "exp_e", "1E9", nullptr, kLimitRefusedLarge, nullptr },
	{ "raise_e", "1E9", nullptr, kLimitRefusedLarge, nullptr },
	{ "exp_b", "1E24", nullptr, kLimitRefusedLarge, nullptr },
	{ "factorial", "1E6", nullptr, kLimitRefusedLarge, nullptr },
	{ "gamma", "1000000.5", nullptr, kLimitRefusedLarge, nullptr },
	{ "hypersine", "-1E9", nullptr, kLimitRefusedLarge, nullptr },
	{ "hypercosine", "1E24", nullptr, kLimitRefusedLarge, nullptr },
	{ "exp_n", "1.2345", "1E24", kLimitRefusedLarge, nullptr },
	{ "power", "1.2345", "1E9", kLimitRefusedLarge, nullptr },
	{ "raise_n", "1.2345", "1E9", kLimitRefusedLarge, nullptr },
	{ "root_n", "1.2345", "3E-25", kLimitRefusedLarge, nullptr },
	{ "hypersecant", "1E24", nullptr, kLimitRefusedSmall, nullptr },
	{ "hypercosecant", "-1E9", nullptr, kLimitRefusedSmall, nullptr },
	{ "exp_e", "-1E24", nullptr, kLimitSpelled, "ZERO" },
	{ "exp_e", "-1000", nullptr, kLimitSpelled, "ZERO" },
	{ "exp_b", "-1E24", nullptr, kLimitSpelled, "ZERO" },
	{ "raise_e", "-1E9", nullptr, kLimitSpelled, "ZERO" },
	{ "power", "1.2345", "-1E24", kLimitSpelled, "ZERO" },
	{ "exp_n", "1.2345", "-1E24", kLimitSpelled, "ZERO" },
	{ "raise_n", "1.2345", "-1E9", kLimitSpelled, "ZERO" },
	{ "gamma", "-100000.5", nullptr, kLimitSpelled, "ZERO" },
	{ "hypertangent", "1E9", nullptr, kLimitSpelled, "ONE" },
	{ "hypertangent", "-1E9", nullptr, kLimitSpelled, "MONE" },
	{ "hypercotangent", "-1E24", nullptr, kLimitSpelled, "MONE" },
	{ "root_n", "-1E24", "1E-6", kLimitDomain, nullptr },
	// Past what a double holds, so sized as infinite...
	{ "exp_e", "1E400", nullptr, kLimitRefusedLarge, nullptr },
	{ "power", "1.2345", "-1E400", kLimitSpelled, "ZERO" },
	// 5.45E-37, which rounds up to the last place: nought takes a margin...
	{ "exp_e", "-83.5", nullptr, kLimitSpelled, "0.000000000000000000000000000000000001" },
	// Near a pole, sin(pi x) is small and the answer is not nought: the
	//	bound has to allow for it...
	{ "gamma", "-35.0000000001", nullptr, kLimitSpelled, "0.000000000000000000000000000000967759" },
};

// Either side of a limit of a thousand places...
static const SLimitCase	kLimitBoundary[] =
{
	{ "exp_e", "2302", nullptr, kLimitAnswered, nullptr },
	{ "exp_e", "2304", nullptr, kLimitRefusedLarge, nullptr },
	{ "factorial", "449", nullptr, kLimitAnswered, nullptr },
	{ "factorial", "451", nullptr, kLimitRefusedLarge, nullptr },
	{ "power", "2", "3300", kLimitAnswered, nullptr },
	{ "power", "2", "3330", kLimitRefusedLarge, nullptr },
	{ "hypercosine", "2280", nullptr, kLimitAnswered, nullptr },
	{ "hypercosine", "2320", nullptr, kLimitRefusedLarge, nullptr },
	{ "hypersecant", "2280", nullptr, kLimitAnswered, nullptr },
	{ "hypersecant", "2320", nullptr, kLimitRefusedSmall, nullptr },
	{ "root_n", "4", "0.00061", kLimitAnswered, nullptr },
	{ "root_n", "4", "0.0006", kLimitRefusedLarge, nullptr },
};

static long LimitCall( const SLimitCase& limitCase, CLZNumber* pAnswer )
{
	typedef long (*Unary)( CLZNumber*, CLZNumber* );
	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );

	struct { const char* m_szName; Unary m_fn; } kUnary[] =
	{
		{ "exp_e", lz_exp_e }, { "raise_e", lz_raise_e }, { "exp_b", lz_exp_b },
		{ "factorial", lz_factorial }, { "gamma", lz_gamma }, { "hypersine", lz_hypersine },
		{ "hypercosine", lz_hypercosine }, { "hypertangent", lz_hypertangent },
		{ "hypercotangent", lz_hypercotangent }, { "hypersecant", lz_hypersecant },
		{ "hypercosecant", lz_hypercosecant }
	};
	struct { const char* m_szName; Binary m_fn; } kBinary[] =
	{
		{ "exp_n", lz_exp_n }, { "power", lz_power }, { "raise_n", lz_raise_n }, { "root_n", lz_root_n }
	};

	CLZNumber	x( limitCase.m_szLeft );
	CLZNumber	y( (limitCase.m_szRight != nullptr) ? limitCase.m_szRight : "0" );

	for ( auto& u : kUnary )
		if ( ::strcmp( u.m_szName, limitCase.m_szFunction ) == 0 )
			return u.m_fn( &x, pAnswer );

	for ( auto& b : kBinary )
		if ( ::strcmp( b.m_szName, limitCase.m_szFunction ) == 0 )
			return b.m_fn( &x, &y, pAnswer );

	return -1;
}

// Counted against what each case expects, answering how many were not...
static long LimitCheck( const SLimitCase* pCases, size_t nCases, std::string* psFirst )
{
	long	l_wrong = 0;

	for ( size_t i = 0; i < nCases; i++ )
	{
		const SLimitCase&	c = pCases[ i ];
		CLZNumber			answer;

		glWorkingPrecision = 36;

		long	l_code = LimitCall( c, &answer );
		char*	sz_got = answer.GetStringRep();
		bool	b_inf = answer.FlagIsOn( kfInfinity );
		bool	b_nan = answer.FlagIsOn( kfQNaN );
		bool	b_zero = (sz_got != nullptr) && (::strspn( sz_got, "0.-" ) == ::strlen( sz_got ));
		bool	b_right = false;

		if ( c.m_iExpect == kLimitRefusedLarge )
			b_right = (l_code == kTooLarge) && b_inf && b_zero;
		else if ( c.m_iExpect == kLimitRefusedSmall )
			b_right = (l_code == kTooLarge) && !b_inf && !b_nan && b_zero;
		else if ( c.m_iExpect == kLimitDomain )
			b_right = (l_code == kFailed) && b_nan;
		else if ( c.m_iExpect == kLimitAnswered )
			b_right = (l_code == kNoError) && !b_inf && !b_nan && !b_zero;
		else
		{
			std::string	s_expect = (::strcmp( c.m_szSpelled, "ZERO" ) == 0) ? "0.000000000000000000000000000000000000" :
									   (::strcmp( c.m_szSpelled, "ONE" ) == 0) ? "1.000000000000000000000000000000000000" :
									   (::strcmp( c.m_szSpelled, "MONE" ) == 0) ? "-1.000000000000000000000000000000000000" : c.m_szSpelled;

			b_right = (l_code == kNoError) && !b_inf && !b_nan && (sz_got != nullptr) && (s_expect == sz_got);
		}

		if ( !b_right && (l_wrong++ == 0) )
			*psFirst = std::string( c.m_szFunction ) + "(" + c.m_szLeft + ((c.m_szRight != nullptr) ? std::string( ", " ) + c.m_szRight : "") +
					   ") code " + std::to_string( l_code ) + ", " + ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");

		::free( sz_got );
	}

	glWorkingPrecision = 36;

	return l_wrong;
}

static void TestAnswerLimit()
{
	Banner( "An answer too large to hold" );

	long		l_saved = glWorkingPrecision;
	std::string	s_first;
	char		szDetail[ 200 ];

	// The lowered limit first, so that a check taken out fails before the
	//	cases at the default can hang. The setter answers the previous
	//	limit, and nought puts back the default...
	long	l_was = lzu_set_answer_limit( 1000 );
	size_t	n_boundary = sizeof( kLimitBoundary ) / sizeof( kLimitBoundary[ 0 ] );
	long	l_wrong = LimitCheck( kLimitBoundary, n_boundary, &s_first );
	long	l_set = lzu_set_answer_limit( 0 );
	long	l_back = lzu_set_answer_limit( kAnswerLimitPlaces );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu either side of 1,000 places", n_boundary );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_boundary, s_first.c_str() );

	CheckTrue( "either side of a lowered limit", l_wrong == 0, szDetail );

	::snprintf( szDetail, sizeof( szDetail ), "answered %ld, then %ld, then %ld", l_was, l_set, l_back );

	CheckTrue( "the limit's setter", (l_was == kAnswerLimitPlaces) && (l_set == 1000) && (l_back == kAnswerLimitPlaces), szDetail );

	size_t	n_cases = sizeof( kLimitCases ) / sizeof( kLimitCases[ 0 ] );

	l_wrong = LimitCheck( kLimitCases, n_cases, &s_first );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu refused, nought, one or a domain failure, at once", n_cases );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_cases, s_first.c_str() );

	CheckTrue( "past the limit, at the default", l_wrong == 0, szDetail );

	glWorkingPrecision = l_saved;
}

/*
	**combinations on the whole numbers, k at most n, nought past n** --
	Rich's call, 2026-10-11. It used to divide three factorials and check
	none of them, so k past n, a fraction or a negative came back a
	not-a-number reported as success, and a large n either never answered
	or, past the answer limit, answered NaN. Every reference here comes
	from Python's math.comb, generated rather than typed. C(2000, 1000)
	takes the path by primes. Pascal's rule at C(9151, 100) sets the two
	paths against each other -- C(9150, 99) is the product, the other two
	the primes, under the rule measured on 2026-10-11 -- and at C(2E7, 50)
	all three are the product, past the sieve's reach.
*/

static const char	kComb2000Of1000[] =
	"204815162698948971433516250298082504439642488798139703382038"
	"263767174818620208375582893299418261020620146476631999802369"
	"241548179800452479201804754976926157856301289663432064714851"
	"152395251651227768588611539546256147907378668464154444533617"
	"613770073855673814589630071306510455959514479888746206368718"
	"514551828551173166276253663773084682932255389049743859481431"
	"755030783796444370810085163724827462791417016619883764840843"
	"541430817785947037746565188475514680749694674923803033101818"
	"723298009668567458560252549910118113525353465888794196665367"
	"490451130611009631190627034250229315591110897673396399114912"
	"0";

static void TestCombinationsWhole()
{
	Banner( "combinations on the whole numbers" );

	long		l_saved = glWorkingPrecision;
	CLZNumber	n, k, result;

	glWorkingPrecision = 36;

	struct { const char* m_szN; const char* m_szK; long m_lBase; const char* m_szExpect; } kPins[] =
	{
		{ "5", "7", 10, "0" },
		{ "0", "0", 10, "1" },
		{ "12", "0", 10, "1" },
		{ "12", "12", 10, "1" },
		{ "100000", "3", 10, "166661666700000" },
		{ "1E24", "2", 10, "499999999999999999999999500000000000000000000000" },
		{ "1E24", "3", 10, "166666666666666666666666166666666666666666666667000000000000000000000000" },
		{ "30", "15", 7, "3562322445" },
		{ "40", "20", 16, "20184B1B34" },
		{ "20", "10", 2, "101101000110110100" },
		{ "60", "30", 12, "78141E0273036094" },
		{ "2000", "1000", 10, kComb2000Of1000 },
	};

	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kPins )
	{
		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( pin.m_szN );
		k.Assign( pin.m_szK );

		if ( pin.m_lBase != 10 )
		{
			n.SetBase( pin.m_lBase );
			k.SetBase( pin.m_lBase );
		}

		CLZNumber	answer;
		long		l_code = lz_combinations( &n, &k, &answer );
		char*		sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
		{
			if ( l_wrong++ == 0 )
				s_first = std::string( "C(" ) + pin.m_szN + ", " + pin.m_szK + ") base " + std::to_string( pin.m_lBase ) +
						  ": " + ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");
		}

		::free( sz_got );
	}

	char	szDetail[ 200 ];
	size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact, in five bases", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "exact, against math.comb", l_wrong == 0, szDetail );

	// Pascal's rule, C(n, k) = C(n - 1, k - 1) + C(n - 1, k), across the
	//	switch from a product to the primes and well past it...
	struct { long m_lN; long m_lK; } kPascal[] = { { 9151, 100 }, { 3001, 1500 }, { 40000, 17000 }, { 20000000, 50 } };

	l_wrong = 0;

	for ( auto& pascal : kPascal )
	{
		CLZNumber	whole, lower, upper, sum, one;

		one.Assign( 1L );
		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( pascal.m_lN );
		k.Assign( pascal.m_lK );
		lz_combinations( &n, &k, &whole );
		lz_subtract( &n, &one );
		lz_combinations( &n, &k, &upper );
		lz_subtract( &k, &one );
		lz_combinations( &n, &k, &lower );
		lz_add( &lower, &upper, &sum );

		if ( (sum.Compare( &whole ) != kEqualTo) || whole.IsZero() )
			l_wrong++;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of 4 sums differ", l_wrong );

	CheckTrue( "Pascal's rule across the paths", l_wrong == 0, szDetail );

	// Refused, as factorial refuses, and past the limit...
	struct { const char* m_szN; const char* m_szK; long m_lCode; bool m_bNaN; } kRefusals[] =
	{
		{ "7.5", "2", kFailed, true }, { "5", "2.5", kFailed, true },
		{ "-5", "2", kFailed, true }, { "5", "-2", kFailed, true },
		{ "333000", "166500", kTooLarge, false }, { "1E24", "5E23", kTooLarge, false },
	};

	l_wrong = 0;

	for ( auto& refusal : kRefusals )
	{
		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( refusal.m_szN );
		k.Assign( refusal.m_szK );

		CLZNumber	answer;
		long		l_code = lz_combinations( &n, &k, &answer );
		bool		b_flag = refusal.m_bNaN ? answer.FlagIsOn( kfQNaN ) : answer.FlagIsOn( kfInfinity );

		if ( (l_code != refusal.m_lCode) || !b_flag || !answer.IsZero() )
			l_wrong++;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of 6 not refused as they should be", l_wrong );

	CheckTrue( "fractions, negatives, and past the limit", l_wrong == 0, szDetail );

	// In place, the answer over n...
	n.SetBase( 10, false );
	k.SetBase( 10, false );
	n.Assign( "52" );
	k.Assign( "5" );
	lz_combinations( &n, &k );

	CheckExact( "C(52, 5) in place", &n, "2598960" );

	glWorkingPrecision = l_saved;
}

/*
	**permutations on combinations' terms** -- 2026-10-19: whole numbers,
	k at most n, nought past n, a fraction or a negative refused. Every
	reference comes from Python's math.perm, checked a second way as
	n!/(n - k)!. P(n, n) is set against factorial, and P(n, k) against
	C(n, k) k! across combinations' two paths -- C(9151, 100) and
	C(3001, 1500) by the primes, C(2E7, 50) past the sieve's reach.

	P(30000, 30000) has 121,288 digits, past the limit, with a count a
	long holds, so only the sizing refuses it. P(1516, 1516) in base 11 is
	3999.97 places: under a limit of 4,000 it must answer, and Stirling's
	leading terms alone, reading ln Gamma(1) as -0.08, refused it.

	P(5, 1E19) and P(12, 1E30) are nought at once. Without the test for k
	past n the product would pass through nought and answer it anyway --
	so no smaller k can tell -- but only after counting the estimate's
	shift up from about -1E19 by ones, which never ends.
*/

static const char	kPerm300Of150[] =
	"535685181583404275428186032812421320244425480898403849171449"
	"831605323821307364556865221020132590157772550834053861514785"
	"223473917928740658494903840566142474001546879997175667314585"
	"767958954631475946055780471274216752464841263702154974602323"
	"872102199455646604695291548872570515202150474859324425490481"
	"9157328965140480000000000000000000000000000000000000";

static void TestPermutationsWhole()
{
	Banner( "permutations on the whole numbers" );

	long		l_saved = glWorkingPrecision;
	CLZNumber	n, k;

	glWorkingPrecision = 36;

	struct { const char* m_szN; const char* m_szK; long m_lBase; const char* m_szExpect; } kPins[] =
	{
		{ "5", "7", 10, "0" },
		{ "5", "1E19", 10, "0" },
		{ "12", "1E30", 7, "0" },
		{ "0", "0", 10, "1" },
		{ "12", "0", 10, "1" },
		{ "12", "12", 10, "479001600" },
		{ "100000", "3", 10, "999970000200000" },
		{ "1E24", "2", 10, "999999999999999999999999000000000000000000000000" },
		{ "1E24", "3", 10, "999999999999999999999997000000000000000000000002000000000000000000000000" },
		{ "30", "15", 7, "1026106626405242536446300" },
		{ "40", "20", 16, "43BA12919CBD651CF88900000" },
		{ "20", "10", 2, "1001110000011001011111011100110000000000" },
		{ "60", "30", 12, "86E1D26456993D957528D972D2DD332400000000000000" },
		{ "300", "150", 10, kPerm300Of150 },
	};

	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kPins )
	{
		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( pin.m_szN );
		k.Assign( pin.m_szK );

		if ( pin.m_lBase != 10 )
		{
			n.SetBase( pin.m_lBase );
			k.SetBase( pin.m_lBase );
		}

		CLZNumber	answer;
		long		l_code = lz_permutations( &n, &k, &answer );
		char*		sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
		{
			if ( l_wrong++ == 0 )
				s_first = std::string( "P(" ) + pin.m_szN + ", " + pin.m_szK + ") base " + std::to_string( pin.m_lBase ) +
						  ": " + ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");
		}

		::free( sz_got );
	}

	char	szDetail[ 200 ];
	size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact, in five bases", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "exact, against math.perm", l_wrong == 0, szDetail );

	// P(n, n) is n!, in base ten and base seven...
	l_wrong = 0;

	for ( long lBase : { 10L, 7L } )
	{
		CLZNumber	arranged, factorial;

		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( 1000L );
		k.Assign( 1000L );
		n.SetBase( lBase );
		k.SetBase( lBase );
		lz_permutations( &n, &k, &arranged );
		lz_factorial( &n, &factorial );

		if ( (arranged.Compare( &factorial ) != kEqualTo) || (arranged.m_lBase != lBase) || arranged.IsZero() )
			l_wrong++;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of 2 differ", l_wrong );

	CheckTrue( "P(1000, 1000) is 1000!", l_wrong == 0, szDetail );

	// P(n, k) is C(n, k) k!, across combinations' paths...
	struct { long m_lN; long m_lK; } kChosen[] = { { 9151, 100 }, { 3001, 1500 }, { 20000000, 50 } };

	l_wrong = 0;

	for ( auto& chosen : kChosen )
	{
		CLZNumber	arranged, ways, orders, product;

		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( chosen.m_lN );
		k.Assign( chosen.m_lK );
		lz_permutations( &n, &k, &arranged );
		lz_combinations( &n, &k, &ways );
		lz_factorial( &k, &orders );
		lz_multiply( &ways, &orders, &product );

		if ( (arranged.Compare( &product ) != kEqualTo) || arranged.IsZero() )
			l_wrong++;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of 3 differ", l_wrong );

	CheckTrue( "P(n, k) is C(n, k) k!", l_wrong == 0, szDetail );

	// Refused, as combinations refuses, and past the limit...
	struct { const char* m_szN; const char* m_szK; long m_lCode; bool m_bNaN; } kRefusals[] =
	{
		{ "7.5", "2", kFailed, true }, { "5", "2.5", kFailed, true },
		{ "-5", "2", kFailed, true }, { "5", "-2", kFailed, true },
		{ "30000", "30000", kTooLarge, false }, { "1E24", "5E23", kTooLarge, false },
	};

	l_wrong = 0;

	for ( auto& refusal : kRefusals )
	{
		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( refusal.m_szN );
		k.Assign( refusal.m_szK );

		CLZNumber	answer;
		long		l_code = lz_permutations( &n, &k, &answer );
		bool		b_flag = refusal.m_bNaN ? answer.FlagIsOn( kfQNaN ) : answer.FlagIsOn( kfInfinity );

		if ( (l_code != refusal.m_lCode) || !b_flag || !answer.IsZero() )
			l_wrong++;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of 6 not refused as they should be", l_wrong );

	CheckTrue( "fractions, negatives, and past the limit", l_wrong == 0, szDetail );

	// A thirtieth of a place under a lowered limit, at k = n...
	{
		long		l_was = lzu_set_answer_limit( 4000 );
		CLZNumber	arranged, factorial;

		n.SetBase( 10, false );
		k.SetBase( 10, false );
		n.Assign( 1516L );
		k.Assign( 1516L );
		n.SetBase( 11 );
		k.SetBase( 11 );

		long	l_code = lz_permutations( &n, &k, &arranged );

		lz_factorial( &n, &factorial );
		lzu_set_answer_limit( l_was );

		CheckTrue( "P(1516, 1516) in base 11 under a 4,000-place limit",
				   (l_code == kNoError) && (arranged.Compare( &factorial ) == kEqualTo) && !(arranged.IsZero()),
				   (l_code == kNoError) ? "" : "refused" );
	}

	// In place, the answer over n...
	n.SetBase( 10, false );
	k.SetBase( 10, false );
	n.Assign( "52" );
	k.Assign( "5" );
	lz_permutations( &n, &k );

	CheckExact( "P(52, 5) in place", &n, "311875200" );

	glWorkingPrecision = l_saved;
}

/*
	**Complex numbers, stage 1** -- 2026-10-21, on Rich's four calls:
	each part to the working precision in places; the layer in LZMath;
	a CLZComplex and lzc_ functions; and stage 1 the operations correctly
	rounded from exact intermediates. Every reference is generated with
	Python's Fraction, the census's own rounding and spelling, never
	typed. With both imaginary parts nought every answer must be the real
	operation's to the character, flags included, in five bases and four
	modes; the four plain operations round each part in its own mode, and
	the square and the modulus as functions, whatever the modes.
*/

static void SetComplex( CLZComplex& z, long lBase, const char* szReal, const char* szImag )
{
	z.m_real.SetBase( lBase, false );
	z.m_real.Assign( szReal );
	z.m_imag.SetBase( lBase, false );
	z.m_imag.Assign( szImag );
}

static bool SpelledAs( CLZNumber* pNumber, const char* szExpected )
{
	char*	sz_got = pNumber->GetStringRep();
	bool	b_same = (sz_got != nullptr) && (::strcmp( sz_got, szExpected ) == 0);

	::free( sz_got );

	return b_same;
}

static bool SameAnswer( CLZNumber* pGot, CLZNumber* pWant )
{
	char*	sz_got = pGot->GetStringRep();
	char*	sz_want = pWant->GetStringRep();
	bool	b_same = (sz_got != nullptr) && (sz_want != nullptr) && (::strcmp( sz_got, sz_want ) == 0) &&
					 (pGot->m_lBase == pWant->m_lBase) && (pGot->m_lFlags == pWant->m_lFlags);

	::free( sz_got );
	::free( sz_want );

	return b_same;
}

static void TestComplexExact()
{
	Banner( "complex numbers, stage 1: the exact operations" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 200 ];

	glWorkingPrecision = 36;

	// Products and a square, against Fraction, each part in its target's
	//	mode -- the square's modes are there to be ignored...
	struct { const char* m_szOp; long m_lBase; const char* m_szOperand[ 4 ]; int m_iModeReal; int m_iModeImag;
			 const char* m_szReal; const char* m_szImag; } kPins[] =
	{
		{ "mul", 10, { "1.5", "2.25", "-3", "0.5" }, 2, 2,
		  "-5.625",
		  "-6" },
		{ "mul", 10, { "14.11032319755513850066124465679", "654.7370586456225651468334899643", "-803.5330706104415527921964234", "-55.9769176833342967151568331074455" }, 2, 2,
		  "25312.051109797209708681445624833406591174",
		  "-526892.731576080634953339250248668288468632" },
		{ "mul", 10, { "-117.458506359870549986335864905453", "-688.0907003100521259705682667515", "278.8246491760914081627297497", "-677.1567532659196672754865379964023" }, 3, 4,
		  "-498695.591402966448535760169538235612892266",
		  "-112318.827305167017497759639417822952533550" },
		{ "mul", 7, { "-43.4424154652540126003361", "660.1422300025550513013420441", "-433.16650601530402236455426", "233.1442441052102220242" }, 2, 2,
		  "-201330.425351206414256561524131666044052054",
		  "-443135.544421330214162423035316032103111244" },
		{ "mul", 16, { "-CE4.C2944D9697D31A74C9D3", "F74.AFC91DB26DD3015FF0222", "-659.1667FE837282E3190621CAC", "-620.045F8602E8FCEBD0A1F" }, 7, 0,
		  "B08430.4BEEA71A35A259F35A1C1BD577C688BFF4A9",
		  "-1323A9.F2FF961C7FCE8696E854310CE88C312E532F" },
		{ "mul", 2, { "1.111000001100111110001101001101001", "-10.11010100000100110000010001111101", "110.00010100000111010010011000011000000011001011", "-0.1101011100110000000001100010000001111" }, 2, 2,
		  "1001.000010100000000001110111101010000000",
		  "-10010.110001010111111010001011001100011001" },
		{ "mul", 12, { "192.38314E2DD424475319844075", "-427.31413D2799D957716DEED4", "883.D8E83404E715E62974272", "-394.670D966968437D33D865883D8" }, 1, 5,
		  "-7164.312E8E6316106D011769099E234E4752D405",
		  "-374104.1DDE2E79920563497860E54825496D3D2273" },
		{ "sq", 10, { "45.453436968403500080858869", "-90.59340834090077244692408982", "0", "0" }, 3, 4,
		  "-6141.150702580559683804075659987540311991",
		  "-8235.563551551946328699218482193754474574" },
	};

	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kPins )
	{
		CLZComplex	left;
		CLZComplex	right;
		CLZComplex	target;

		SetComplex( left, pin.m_lBase, pin.m_szOperand[ 0 ], pin.m_szOperand[ 1 ] );
		SetComplex( right, pin.m_lBase, pin.m_szOperand[ 2 ], pin.m_szOperand[ 3 ] );
		target.m_real.SetRoundingMode( static_cast<RoundingMode>( pin.m_iModeReal ) );
		target.m_imag.SetRoundingMode( static_cast<RoundingMode>( pin.m_iModeImag ) );

		long	l_code = (::strcmp( pin.m_szOp, "sq" ) == 0) ? lzc_square( &left, &target ) : lzc_multiply( &left, &right, &target );
		bool	b_modesKept = (target.m_real.OwnMode() == pin.m_iModeReal) && (target.m_imag.OwnMode() == pin.m_iModeImag);

		if ( (l_code != kNoError) || !SpelledAs( &target.m_real, pin.m_szReal ) || !SpelledAs( &target.m_imag, pin.m_szImag ) ||
			 (target.m_real.m_lBase != pin.m_lBase) || (target.m_imag.m_lBase != pin.m_lBase) || !b_modesKept )
		{
			if ( l_wrong++ == 0 )
				s_first = std::string( pin.m_szOp ) + " base " + std::to_string( pin.m_lBase ) + " modes " +
						  std::to_string( pin.m_iModeReal ) + "/" + std::to_string( pin.m_iModeImag );
		}
	}

	size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact, in five bases, each target's modes kept", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "products and a square, against Fraction", l_wrong == 0, szDetail );

	// Exact quotients and moduli...
	{
		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	q;
		CLZNumber	m;

		SetComplex( z, 10, "3", "4" );
		SetComplex( w, 10, "1", "2" );
		lzc_divide( &z, &w, &q );

		bool	b_ok = SpelledAs( &q.m_real, "2.2" ) && SpelledAs( &q.m_imag, "-0.4" );

		SetComplex( z, 10, "1", "1" );
		SetComplex( w, 10, "1", "-1" );
		lzc_divide( &z, &w, &q );
		b_ok = b_ok && SpelledAs( &q.m_real, "0" ) && SpelledAs( &q.m_imag, "1" );

		SetComplex( z, 10, "3", "4" );
		lzc_abs( &z, &m );
		b_ok = b_ok && SpelledAs( &m, "5" );

		SetComplex( z, 10, "1.5", "2.25" );
		m.SetRoundingMode( RM_round_ceiling );
		lzc_abs( &z, &m );
		b_ok = b_ok && SpelledAs( &m, "2.704163456597991969839415950602871960" );

		SetComplex( z, 10, "0", "-7.5" );
		lzc_abs( &z, &m );
		b_ok = b_ok && SpelledAs( &m, "7.5" );

		CheckTrue( "(3+4i)/(1+2i), (1+i)/(1-i), |3+4i|, |1.5+2.25i|, |-7.5i|", b_ok, b_ok ? "" : "an answer not as Fraction gives it" );
	}

	// With both imaginary parts nought, the real operation's answer to the
	//	character, its flags and base included, in five bases and four modes...
	{
		typedef long (*RealBinary)( CLZNumber*, CLZNumber*, CLZNumber* );
		typedef long (*ComplexBinary)( CLZComplex*, CLZComplex*, CLZComplex* );

		struct { const char* m_szName; RealBinary m_fnReal; ComplexBinary m_fnComplex; } kBinaries[] =
		{
			{ "add", lz_add, lzc_add }, { "subtract", lz_subtract, lzc_subtract },
			{ "multiply", lz_multiply, lzc_multiply }, { "divide", lz_divide, lzc_divide },
		};

		const char*		kPairs[][ 2 ] = { { "1.25", "0.3" }, { "-7", "3" }, { "2.5", "0.5" }, { "123.456", "-0.0789" },
										  { "0.000000000000000000000000000000000000271", "37" } };
		const long		kBases[] = { 10, 7, 16, 3, 2 };
		RoundingMode	kModes[] = { RM_round_half_even, RM_round_ceiling, RM_round_floor, RM_round_up };

		long	l_cases = 0;

		l_wrong = 0;
		s_first.clear();

		for ( long lBase : kBases )
		{
			for ( auto& pair : kPairs )
			{
				for ( RoundingMode rm : kModes )
				{
					CLZNumber	x;
					CLZNumber	y;

					x.Assign( pair[ 0 ] );
					y.Assign( pair[ 1 ] );
					x.SetBase( lBase, true );
					y.SetBase( lBase, true );

					auto	Complex = [ & ]( CLZNumber* pReal, CLZComplex& z )
					{
						z.m_real.Assign( pReal );
						z.m_imag.SetBase( lBase, false );
						z.m_imag.Assign( 0L );
					};

					auto	Note = [ & ]( const char* szName, bool bSame )
					{
						l_cases++;

						if ( !bSame && (l_wrong++ == 0) )
							s_first = std::string( szName ) + " of " + pair[ 0 ] + ", " + pair[ 1 ] + " in base " + std::to_string( lBase );
					};

					for ( auto& op : kBinaries )
					{
						CLZComplex	z;
						CLZComplex	w;
						CLZComplex	t;
						CLZNumber	real;

						Complex( &x, z );
						Complex( &y, w );
						t.m_real.SetRoundingMode( rm );
						real.SetRoundingMode( rm );

						long	l_codeReal = op.m_fnReal( &x, &y, &real );
						long	l_codeComplex = op.m_fnComplex( &z, &w, &t );

						t.m_imag.CheckForZero();
						Note( op.m_szName, (l_codeReal == l_codeComplex) && SameAnswer( &t.m_real, &real ) && t.m_imag.IsZero() );
					}

					// ...and the four with one operand, as lz_square, lz_negate,
					//	lz_plus and lz_abs answer...
					CLZComplex	z;
					CLZComplex	t;
					CLZNumber	real;
					CLZNumber	modulus;

					Complex( &x, z );
					t.m_real.SetRoundingMode( rm );
					real.SetRoundingMode( rm );
					lz_square( &x, &real );
					lzc_square( &z, &t );
					Note( "square", SameAnswer( &t.m_real, &real ) );

					lz_negate( &x, &real );
					lzc_negate( &z, &t );
					Note( "negate", SameAnswer( &t.m_real, &real ) );

					lz_plus( &x, &real );
					lzc_conjugate( &z, &t );
					Note( "conjugate", SameAnswer( &t.m_real, &real ) );

					modulus.SetRoundingMode( rm );
					lz_abs( &x, &real );
					lzc_abs( &z, &modulus );
					Note( "abs", SameAnswer( &modulus, &real ) );
				}
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld the real operation's, flags and base included", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld not, first %s", l_wrong, l_cases, s_first.c_str() );

		CheckTrue( "a real pair answers as the real operations", l_wrong == 0, szDetail );
	}

	// A right operand in another base: base two ends in ten and converts
	//	exactly; base three does not, and the answer is the quotient over its
	//	power of three -- exact here, since the values are thirds...
	{
		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	t;

		SetComplex( z, 10, "1.5", "2" );
		SetComplex( w, 2, "0.1", "1" );
		lzc_multiply( &z, &w, &t );

		bool	b_ok = SpelledAs( &t.m_real, "-1.25" ) && SpelledAs( &t.m_imag, "2.5" ) &&
					   (t.m_real.m_lBase == 10) && (t.m_imag.m_lBase == 10);

		SetComplex( z, 10, "3", "0" );
		SetComplex( w, 3, "0.1", "0.2" );
		lzc_multiply( &z, &w, &t );
		b_ok = b_ok && SpelledAs( &t.m_real, "1" ) && SpelledAs( &t.m_imag, "2" );

		SetComplex( z, 10, "1", "1" );
		SetComplex( w, 3, "0.1", "0" );
		lzc_divide( &z, &w, &t );
		b_ok = b_ok && SpelledAs( &t.m_real, "3" ) && SpelledAs( &t.m_imag, "3" );

		// ...and a third, inexact, written as the real mixed product writes it...
		CLZNumber	one;
		CLZNumber	third;
		CLZNumber	real;

		one.Assign( 1L );
		third.SetBase( 3, false );
		third.Assign( "0.1" );
		lz_multiply( &one, &third, &real );
		SetComplex( z, 10, "1", "0" );
		SetComplex( w, 3, "0.1", "0" );
		lzc_multiply( &z, &w, &t );
		b_ok = b_ok && SameAnswer( &t.m_real, &real );

		// ...and an inexact product with a base-two operand, which converts
		//	exactly and is written as a product, not a quotient...
		CLZNumber	x;
		CLZNumber	y;

		x.Assign( "0.333333333333333333333333333333333333" );
		y.SetBase( 2, false );
		y.Assign( "0.0000000000000000000000000000000000000000000000000000000000001" );
		lz_multiply( &x, &y, &real );
		SetComplex( z, 10, "0.333333333333333333333333333333333333", "0" );
		SetComplex( w, 2, "0.0000000000000000000000000000000000000000000000000000000000001", "0" );
		lzc_multiply( &z, &w, &t );
		b_ok = b_ok && SameAnswer( &t.m_real, &real );

		CheckTrue( "a right operand in base two, and in base three", b_ok, b_ok ? "" : "an answer not exact, or not as the real mixed pair writes it" );
	}

	// In place, into the right operand, and z times z as its square...
	{
		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	fresh;

		SetComplex( z, 10, "1.234567890123456789012345678901", "-9.87654321098765432109876543210987" );
		SetComplex( w, 10, "-3.141592653589793238462643383279502884", "2.718281828459045235360287471352662497" );
		lzc_multiply( &z, &w, &fresh );
		lzc_multiply( &w, &w, &z );		// z = w w...

		CLZComplex	square;

		lzc_square( &w, &square );

		bool	b_ok = SameAnswer( &z.m_real, &square.m_real ) && SameAnswer( &z.m_imag, &square.m_imag );

		SetComplex( z, 10, "1.234567890123456789012345678901", "-9.87654321098765432109876543210987" );
		lzc_multiply( &z, &w );				// in place...
		b_ok = b_ok && SameAnswer( &z.m_real, &fresh.m_real ) && SameAnswer( &z.m_imag, &fresh.m_imag );

		SetComplex( z, 10, "1.234567890123456789012345678901", "-9.87654321098765432109876543210987" );
		lzc_multiply( &z, &w, &w );			// ...and into the right operand
		b_ok = b_ok && SameAnswer( &w.m_real, &fresh.m_real ) && SameAnswer( &w.m_imag, &fresh.m_imag );

		CheckTrue( "in place, into the right operand, and w w as w squared", b_ok, b_ok ? "" : "an answer that moved with its target" );
	}

	// The point at infinity and the not-a-number...
	{
		CLZNumber	one;
		CLZNumber	zero;

		one.Assign( 1L );
		zero.Assign( 0L );

		auto	Infinite = [ & ]( CLZComplex& z )
		{
			lz_divide( &one, &zero, &z.m_real );
			z.m_imag.Assign( 0L );
		};

		auto	Has = []( CLZComplex& z, unsigned long lFlag )
		{
			return z.m_real.FlagIsOn( lFlag ) && z.m_imag.FlagIsOn( lFlag );
		};

		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	t;
		long		l_failed = 0;

		Infinite( z );
		SetComplex( w, 10, "2", "0" );
		lzc_multiply( &z, &w, &t );
		l_failed += !Has( t, kfInfinity );

		SetComplex( w, 10, "0", "0" );
		lzc_multiply( &z, &w, &t );
		l_failed += !Has( t, kfQNaN );

		Infinite( w );
		lzc_add( &z, &w, &t );
		l_failed += !Has( t, kfQNaN );

		SetComplex( w, 10, "1", "1" );
		lzc_subtract( &z, &w, &t );
		l_failed += !Has( t, kfInfinity );

		lzc_divide( &w, &z, &t );
		t.m_real.CheckForZero();
		t.m_imag.CheckForZero();
		l_failed += !(t.m_real.IsZero() && t.m_imag.IsZero() && !Has( t, kfInfinity ) && !Has( t, kfQNaN ));

		lzc_square( &z, &t );
		l_failed += !Has( t, kfInfinity );

		CLZNumber	m;

		lzc_abs( &z, &m );
		l_failed += !m.FlagIsOn( kfInfinity );

		SetComplex( w, 10, "0", "0" );
		SetComplex( z, 10, "1", "-1" );
		l_failed += (lzc_divide( &z, &w, &t ) != kDivideByZero) || !Has( t, kfInfinity ) || !Has( t, kfDiv0 );

		SetComplex( z, 10, "0", "0" );
		l_failed += (lzc_divide( &z, &w, &t ) != kUndefinedResult) || !Has( t, kfQNaN ) || !Has( t, kfDiv0 );

		lz_divide( &zero, &zero, &z.m_real );
		z.m_imag.Assign( 1L );
		SetComplex( w, 10, "1", "1" );
		lzc_multiply( &z, &w, &t );
		l_failed += !Has( t, kfQNaN );

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 10 not as the rules have it", l_failed );
		CheckTrue( "infinity, nought and not-a-number", l_failed == 0, szDetail );
	}

	// Refused: parts in two bases, and a modulus with nowhere to go...
	{
		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	t;
		CLZNumber	m;

		z.m_real.Assign( "1" );
		z.m_imag.SetBase( 16, false );
		z.m_imag.Assign( "1" );
		SetComplex( w, 10, "1", "1" );

		bool	b_ok = (lzc_add( &z, &w, &t ) == kBaseMismatch) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );

		b_ok = b_ok && (lzc_multiply( &w, &z, &t ) == kBaseMismatch);
		b_ok = b_ok && (lzc_abs( &z, &m ) == kBaseMismatch) && m.FlagIsOn( kfQNaN );
		b_ok = b_ok && (lzc_abs( &w, nullptr ) == kFailed);

		CheckTrue( "two bases in one value, and a modulus with no target", b_ok, b_ok ? "" : "not refused as it should be" );
	}

	// Every answer in its left operand's base...
	{
		long	l_off = 0;

		typedef long (*ComplexBinary)( CLZComplex*, CLZComplex*, CLZComplex* );

		for ( ComplexBinary fn : { lzc_add, lzc_subtract, lzc_multiply, lzc_divide } )
		{
			CLZComplex	z;
			CLZComplex	w;
			CLZComplex	t;

			SetComplex( z, 16, "A.8", "-3" );
			SetComplex( w, 10, "2.5", "7" );
			t.m_real.SetBase( 2, false );
			t.m_imag.SetBase( 2, false );
			fn( &z, &w, &t );
			l_off += (t.m_real.m_lBase != 16) + (t.m_imag.m_lBase != 16);
		}

		CLZComplex	z;
		CLZComplex	t;

		SetComplex( z, 7, "1.5", "-2" );
		lzc_square( &z, &t );
		l_off += (t.m_real.m_lBase != 7) + (t.m_imag.m_lBase != 7);

		::snprintf( szDetail, sizeof( szDetail ), "%ld parts in another base", l_off );
		CheckTrue( "every answer in its left operand's base", l_off == 0, szDetail );
	}

	glWorkingPrecision = l_saved;
}

/*
	**The two-argument arctangent and the argument, stage 2a** --
	2026-10-22, on Rich's calls: atan2 and lzc_arg in radians first, the
	angle unit coming to every function at once in 2b; the origin
	answering 0, as std::atan2 and std::arg do. Every expected angle is
	generated -- atan2 of the exact operands by Python's decimal at 320
	digits, through atan(y / x) with the quadrant added, a road of its
	own -- and rounded half-even at the places asked, never typed. The
	near halves lie 10^-6 to 10^-30 of a unit from one, inside the
	retake's reach; the ties are tiny quotients ending exactly on a half,
	which only lzi_arctangent2_decide settles -- and two that do not end,
	a hair above a half: the second nearer it than the 38 significant
	digits the decider's quotient carries, which it must decline and the
	retake settle at its cap.
*/

static void TestArcTangent2()
{
	Banner( "complex numbers, stage 2a: atan2 and the argument" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 200 ];

	struct { const char* m_szFamily; long m_lBase; long m_lPlaces; const char* m_szY; long m_lBaseX; const char* m_szX;
			 const char* m_szAngle; } kPins[] =
	{
		{ "random", 10, 36, "0.09233", 10, "9.4185785629",
		  "0.009802651427047478029417244750782264" },
		{ "random", 10, 36, "90673.974389", 10, "-7264.862999",
		  "1.650746251449942825916411246545307964" },
		{ "random", 10, 36, "-744039", 10, "-570.6",
		  "-1.571563221928267346200302147029258625" },
		{ "random", 10, 36, "-0.0001188", 10, "653054.0863249",
		  "-0.000000000181914488382660517835215137" },
		{ "random", 2, 36, "0.11", 2, "111",
		  "0.000110110101001100000111011011110110" },
		{ "random", 2, 36, "100100.1110001", 2, "-100101101.11101",
		  "11.000001010010000010010111111110010000" },
		{ "random", 2, 36, "-110.0101011", 2, "-0.101001",
		  "-1.101010111110101110010111110001010100" },
		{ "random", 2, 36, "-1100.111100011", 2, "0.0011",
		  "-1.100011100110101001101000010001110101" },
		{ "random", 3, 36, "1022.01101", 3, "0.21011011",
		  "1.112210121222220011202011120121000221" },
		{ "random", 3, 36, "210202012.2", 3, "-22.21110121",
		  "1.120102111221112102210121011200012210" },
		{ "random", 3, 36, "-2202.01111", 3, "-2.12021",
		  "-1.121100111120211110221122220121011002" },
		{ "random", 3, 36, "-22.01120021", 3, "0.000000121",
		  "-1.120102001000201102120220211220212102" },
		{ "random", 7, 36, "55442", 7, "1.531323365",
		  "1.366511413421505321400124136314444363" },
		{ "random", 7, 36, "62240.461435442", 7, "-0.0023635",
		  "1.366532450311252045202346514360335665" },
		{ "random", 7, 36, "-0.0000006312", 7, "-640235.6434006",
		  "-3.066365143202623030321441102556340146" },
		{ "random", 7, 36, "-0.04404121", 7, "611053146.1112",
		  "-0.000000000051266535224140661360261506" },
		{ "random", 12, 36, "D145811.DE2E173", 12, "56.502D",
		  "1.6D23ED170869994035229E94847E09072019" },
		{ "random", 12, 36, "316", 12, "-8297.22E379842",
		  "3.13D04E5407E1D892E74568E70E3529530173" },
		{ "random", 12, 36, "-5713.02E8906", 12, "-E915834.61E4785",
		  "-3.1842438902E25261396125510321729E9441" },
		{ "random", 12, 36, "-9956.18084", 12, "0.056",
		  "-1.6D23E9E585ED212531720359E616E7753998" },
		{ "random", 16, 36, "C0703E53C.F", 16, "83F29C.4242C225",
		  "1.9214BCCBBFAD1083880CD56AC6137ADEBA7C" },
		{ "random", 16, 36, "0.0E7AA2", 16, "-2075A5F.0ACDEBE",
		  "3.243F6A81628F4B8A8E7FBB3D47689050E25B" },
		{ "random", 16, 36, "-8F988E73C.49EA", 16, "-41D06D9.C",
		  "-1.92950A48423312072DC7C2F8B2F821101D8B" },
		{ "random", 16, 36, "-830B8925.77C324", 16, "10",
		  "-1.921FB525012C8402CD971BF320A046FF5E8F" },
		{ "random", 10, 20, "-58450.4869", 10, "-85508273.82",
		  "-3.14090908841902207874" },
		{ "random", 16, 50, "-EEEA3.3D4E", 16, "-AC31543C.8EF",
		  "-3.24293764A3AB178AD7438AF15FC46ACA6B5A0F3784F0553026" },
		{ "near", 2, 36, "-0.11010000101011011011110010011010100111111100111000110110010011111000011011000001", 2, "-0.100101000100100110100001001111010101100000010011111001000110100100001011001011101",
		  "-10.001100000100100001000101111101100010" },
		{ "near", 3, 36, "-0.212210120012211021012210210200121220122221220200010110101110111021012100002212211", 3, "0.110202022111220212211010120001122102122012222120112122121000220100201121212020201",
		  "-1.002010111200102212100221210012222222" },
		{ "near", 10, 36, "0.175441998158172214734888450615245595037891780554700716237338984746413355631518689", 10, "-0.984489769008427988149342961498815148479650344448377772928732986742815611050277291",
		  "2.965237941712566312187618929549935175" },
		{ "near", 2, 36, "-0.111110100011001010010011000000000010111111001101101001010010001001100010110001111", 2, "0.001101100011001000011111000001000010011110111111110101000110101110101101111101111",
		  "-1.010110111000001111001101111101110101" },
		{ "near", 3, 36, "0.21101112222002202121000112201010001000021111000121100010020021201000120120202021", 3, "-0.120101222121011201022002102101020122212211201212101000102101202122022011211100201",
		  "2.011210210210120210122120221212012111" },
		{ "near", 16, 36, "0.0403A7A2AF5093439B027768CFA75398E943360E1D5CFA06748F23829D8557E58F1FE231653D7FD42", 16, "-0.FFF7F13A50DD46B5489A9223CE552B1A8A6C7702FF2F94E043AE68FFCEC03B8E6385F116F78495305",
		  "3.203BB81D85A1E68CC0029ECFB5BD37E1C578" },
		{ "near", 3, 36, "-0.01021120110100122110001101112212012010020021020000011112110022020120011020221202", 3, "0.2222011202022112002101112101002011220001212012202222212202201211201210121102002",
		  "-0.010212001220111121200020210220200100" },
		{ "near", 16, 36, "0.E5A412FE8FDE92E8FE58F15D44319A624055AC9E55B74CA2F48FEBAD25F560C61D888A84408929929", 16, "-0.7124482D243E5FFDB630B1BBD456CFB6B37AD6459CD24724A6351444CAFB15F49FC82B9012E778428",
		  "2.0750EE4DA3A1CCA0D4A182F8EC2F9F87E2DB" },
		{ "near", 10, 36, "0.969117307963675670842632346578953615524770723836667263357933813666151635873109326", 10, "-0.246600169110319565844049023751096187373074583419649784906058326833641515546132427",
		  "1.819966834774409760716391741412541525" },
		{ "near", 2, 36, "0.110110100011000000100011000100000111011101000011011001011000100000111010100001111", 2, "0.1000010111100111001000110010101010100011101110011001001110010111001110000110101",
		  "1.000001010011011001100100011011110010" },
		{ "tie", 10, 36, "0.0000000000000000000000000000000000015", 10, "1",
		  "0.000000000000000000000000000000000001" },
		{ "tie", 10, 36, "0.000000000000000000000000000000000007", 10, "2",
		  "0.000000000000000000000000000000000003" },
		{ "tie", 16, 36, "0.00000000000000000000000000000000000028", 16, "0.1",
		  "0.000000000000000000000000000000000002" },
		{ "tie", 2, 36, "0.0000000000000000000000000000000001111", 2, "101",
		  "0.000000000000000000000000000000000001" },
		{ "tie", 10, 36, "0.0000000000000000000000000000000000045", 10, "-1",
		  "3.141592653589793238462643383279502880" },
		{ "tie", 10, 36, "0.000000000000000000000000000000000004500000000000000000000001", 10, "3",
		  "0.000000000000000000000000000000000002" },
		{ "tie", 10, 36, "0.0000000000000000000000000000000000045000000000000000000000000000000000000001", 10, "3",
		  "0.000000000000000000000000000000000002" },
		{ "mixed", 10, 36, "1.25", 3, "0.1",
		  "1.310193935047555634256437689171905312" },
		{ "mixed", 10, 36, "-3", 2, "0.101",
		  "-1.365400937605129227106169150226454233" },
		{ "mixed", 7, 36, "0.1", 10, "-0.3",
		  "2.461062521105001600355026133602310316" },
		{ "mixed", 16, 36, "0.FF", 10, "0.07",
		  "1.8029C1C19B9022EBF8AEC5A3119CEA876962" },
		{ "mixed", 3, 36, "-0.02", 2, "-0.11",
		  "-2.212001020011122222221022010220102012" },
		{ "axis", 10, 36, "0", 10, "-3",
		  "3.141592653589793238462643383279502884" },
		{ "axis", 10, 36, "2", 10, "0",
		  "1.570796326794896619231321691639751442" },
		{ "axis", 10, 36, "-2", 10, "0",
		  "-1.570796326794896619231321691639751442" },
		{ "axis", 7, 36, "0", 7, "-3",
		  "3.066365143203613411026340224465222665" },
		{ "axis", 7, 36, "2", 7, "0",
		  "1.366532421435305204013153445566111332" },
		{ "axis", 7, 36, "-2", 7, "0",
		  "-1.366532421435305204013153445566111332" },
	};

	// Each family's pins, against the generated angle, in the left base...
	for ( const char* sz_family : { "random", "near", "tie", "mixed", "axis" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZNumber	y;
			CLZNumber	x;
			CLZNumber	angle;

			glWorkingPrecision = pin.m_lPlaces;
			y.SetBase( pin.m_lBase, false );
			y.Assign( pin.m_szY );
			x.SetBase( pin.m_lBaseX, false );
			x.Assign( pin.m_szX );

			long	l_code = lz_arctangent2( &y, &x, &angle );

			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &angle, pin.m_szAngle ) || (angle.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::string( "base " ) + std::to_string( pin.m_lBase ) + " y " + pin.m_szY + " x " + pin.m_szX;
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "atan2, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 36;

	// The origin and the positive real axis answer an exact, short nought,
	//	as lz_arctangent( 0 ) does...
	{
		CLZNumber	y;
		CLZNumber	x;
		CLZNumber	angle;
		CLZNumber	nought;
		CLZNumber	reference;
		bool		b_ok = true;

		nought.Assign( 0L );
		lz_arctangent( &nought, &reference );

		for ( const char* sz_x : { "0", "5", "0.000001" } )
		{
			y.Assign( 0L );
			x.Assign( sz_x );
			b_ok = b_ok && (lz_arctangent2( &y, &x, &angle ) == kNoError) && SameAnswer( &angle, &reference );
		}

		CheckTrue( "atan2 at the origin and along +x", b_ok, b_ok ? "0, exactly, as arctan(0)" : "not arctan(0)'s nought" );
	}

	// Inside the right half-plane it is the arctangent of the quotient, where
	//	that quotient ends -- the engine checked against itself, by another road...
	{
		const char*	kPairs[][ 2 ] = { { "1", "1" }, { "1", "4" }, { "-3", "8" }, { "0.5", "16" }, { "7", "0.25" }, { "-1", "0.0625" } };
		bool		b_ok = true;
		std::string	s_first;

		for ( auto& pair : kPairs )
		{
			CLZNumber	y;
			CLZNumber	x;
			CLZNumber	quotient;
			CLZNumber	angle;
			CLZNumber	reference;

			y.Assign( pair[ 0 ] );
			x.Assign( pair[ 1 ] );
			lz_divide( &y, &x, &quotient );
			lz_arctangent( &quotient, &reference );
			lz_arctangent2( &y, &x, &angle );

			if ( !SameAnswer( &angle, &reference ) && s_first.empty() )
				s_first = std::string( pair[ 0 ] ) + " over " + pair[ 1 ];

			b_ok = b_ok && SameAnswer( &angle, &reference );
		}

		CheckTrue( "atan2( y, x ) = arctan( y / x ), x > 0", b_ok, b_ok ? "six pairs, to the character" : s_first.c_str() );
	}

	// A not-a-number or an infinity in either answers a not-a-number, the
	//	division by nought travelling with it; a missing operand is refused...
	{
		CLZNumber	y;
		CLZNumber	x;
		CLZNumber	angle;
		CLZNumber	one;
		CLZNumber	nought;
		bool		b_ok = true;

		one.Assign( 1L );
		nought.Assign( 0L );

		lz_divide( &one, &nought, &x );		// infinite, kfDiv0
		y.Assign( 1L );
		b_ok = b_ok && (lz_arctangent2( &y, &x, &angle ) == kNoError) && angle.FlagIsOn( kfQNaN ) && angle.FlagIsOn( kfDiv0 );

		lz_divide( &one, &nought, &y );
		x.Assign( 1L );
		b_ok = b_ok && (lz_arctangent2( &y, &x, &angle ) == kNoError) && angle.FlagIsOn( kfQNaN );

		lz_divide( &nought, &nought, &y );	// not a number
		b_ok = b_ok && (lz_arctangent2( &y, &x, &angle ) == kNoError) && angle.FlagIsOn( kfQNaN );

		angle.Assign( "12345" );
		b_ok = b_ok && (lz_arctangent2( &one, nullptr, &angle ) == kFailed) && angle.FlagIsOn( kfQNaN );

		CheckTrue( "atan2 of an infinity, a NaN, nothing", b_ok, b_ok ? "a NaN each time, kFailed for nothing" : "an answer not flagged" );
	}

	// In place, over either operand, as into a target...
	{
		CLZNumber	y;
		CLZNumber	x;
		CLZNumber	angle;
		bool		b_ok = true;

		y.Assign( "-2.5" );
		x.Assign( "-7" );
		lz_arctangent2( &y, &x, &angle );

		CLZNumber	overY;
		CLZNumber	overX;

		overY.Assign( "-2.5" );
		x.Assign( "-7" );
		lz_arctangent2( &overY, &x );
		b_ok = b_ok && SameAnswer( &overY, &angle ) && SpelledAs( &x, "-7" );

		y.Assign( "-2.5" );
		overX.Assign( "-7" );
		lz_arctangent2( &y, &overX, &overX );
		b_ok = b_ok && SameAnswer( &overX, &angle ) && SpelledAs( &y, "-2.5" );

		CheckTrue( "atan2 in place over y, and over x", b_ok, b_ok ? "both as into a target" : "an answer differs" );
	}

	// The argument is atan2 of the parts, refusing what the layer refuses...
	{
		const char*	kValues[][ 2 ] = { { "3", "4" }, { "-1", "1" }, { "-2", "0" }, { "0", "-5" }, { "0.25", "-0.75" } };
		bool		b_ok = true;

		for ( auto& v : kValues )
		{
			CLZComplex	z;
			CLZNumber	arg;
			CLZNumber	reference;

			SetComplex( z, 10, v[ 0 ], v[ 1 ] );
			b_ok = b_ok && (lzc_arg( &z, &arg ) == kNoError);
			lz_arctangent2( &z.m_imag, &z.m_real, &reference );
			b_ok = b_ok && SameAnswer( &arg, &reference );
		}

		CLZComplex	z;
		CLZNumber	arg;
		CLZNumber	nought;

		nought.Assign( 0L );

		// ...nought for nought...
		SetComplex( z, 7, "0", "0" );
		CLZNumber	zeroArg;

		b_ok = b_ok && (lzc_arg( &z, &zeroArg ) == kNoError) && SpelledAs( &zeroArg, "0" ) && (zeroArg.m_lBase == 7);

		// ...two bases in one value, and no target...
		z.m_imag.SetBase( 3, false );
		z.m_imag.Assign( "1" );
		arg.Assign( "12345" );
		b_ok = b_ok && (lzc_arg( &z, &arg ) == kBaseMismatch) && arg.FlagIsOn( kfQNaN );
		b_ok = b_ok && (lzc_arg( &z, nullptr ) == kFailed);

		// ...and the point at infinity has no angle...
		CLZNumber	one;

		one.Assign( 1L );
		SetComplex( z, 10, "1", "1" );
		lz_divide( &one, &nought, &z.m_real );
		b_ok = b_ok && (lzc_arg( &z, &arg ) == kNoError) && arg.FlagIsOn( kfQNaN );

		CheckTrue( "lzc_arg: atan2 of the parts", b_ok, b_ok ? "five values, nought, a mismatch, no target, infinity" : "an answer wrong" );
	}

	glWorkingPrecision = l_saved;
}

/*
	**The angle unit, stage 2b** -- 2026-10-22, on Rich's calls: a
	thread_local geAngleUnit set as the working precision is, radians by
	default; Elzed's four units; brought to all fourteen functions that
	take or answer an angle at once. Every expected answer is generated --
	exact Fractions where the answer is rational, Python's decimal at 320
	digits otherwise, reduced by whole turns exactly first -- and spelled
	by the engine's rules, short only where exact. The thirds are the
	exact parts of a turn in small bases at small widths, which a plain
	division rounded again got wrong in base two every time; the near
	halves lie 10^-6 to 10^-28 of a unit from one, inside the retake's
	reach. The widened pins are census cases that the conversion without
	its widening or its guard got wrong -- near the poles, small angles,
	near halves -- and the decider pins tiny angles a hair from a half on
	the side a radian tie decider, left on in a unit, gets wrong.
*/

static void TestAngleUnits()
{
	Banner( "complex numbers, stage 2b: the angle unit" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 200 ];

	CheckTrue( "the angle unit starts as radians", geAngleUnit == AU_radians, "" );

	typedef long (*AngleFunction)( CLZNumber*, CLZNumber* );

	struct { const char* m_szName; AngleFunction m_fn; } kFunctions[] =
	{
		{ "sin", lz_sine }, { "cos", lz_cosine }, { "tan", lz_tangent }, { "cot", lz_cotangent },
		{ "sec", lz_secant }, { "csc", lz_cosecant }, { "asin", lz_arcsine }, { "acos", lz_arccosine },
		{ "atan", lz_arctangent }, { "acot", lz_arccotangent }, { "asec", lz_arcsecant }, { "acsc", lz_arccosecant },
	};

	struct { const char* m_szFamily; const char* m_szName; AngleUnit m_eUnit; long m_lBase; long m_lPlaces;
			 const char* m_szX; const char* m_szAnswer; } kPins[] =
	{
		{ "random", "sin", AU_degrees, 16, 36, "0.00B81DCDABCE",
		  "0.000336A3DAD7EEC2B5B1D9AA1FA0467D5981" },
		{ "random", "sin", AU_grads, 10, 36, "16172516042.2",
		  "0.615386370179171490133443846391107363" },
		{ "random", "sin", AU_millicles, 2, 36, "-10.111",
		  "-0.000001001001111111001001110010110110" },
		{ "random", "cos", AU_degrees, 3, 36, "0.0000211022",
		  "0.222222222222222202212102000112111010" },
		{ "random", "cos", AU_grads, 3, 36, "12.221212121",
		  "0.222212210212002201111122222021002000" },
		{ "random", "cos", AU_millicles, 10, 36, "-0.143936",
		  "0.999999591051546396007280915650357309" },
		{ "random", "tan", AU_degrees, 7, 36, "0.000000355",
		  "0.00000000315632052530632656124165141241446423" },
		{ "random", "tan", AU_grads, 10, 36, "-2.66716",
		  "-0.0419201809179363640554687727579631385" },
		{ "random", "tan", AU_millicles, 3, 36, "-0.2012",
		  "-0.0000101000020200202210110202010212001102" },
		{ "random", "cot", AU_degrees, 16, 36, "CAA.5477F824591",
		  "18.93CEE043C49DCF36EF26A324297833F18100" },
		{ "random", "cot", AU_grads, 2, 36, "0.0010011",
		  "110101100.111000010100000111100001111001010010" },
		{ "random", "cot", AU_millicles, 10, 36, "-228.06",
		  "-0.138733003987637939818046622270650325" },
		{ "random", "sec", AU_degrees, 16, 36, "-23B.37AE4AA",
		  "-1.2B57D4F4B31A916049D33285C164374126D9" },
		{ "random", "sec", AU_grads, 10, 36, "-87.823282938",
		  "5.260187658519735760104885632748668071" },
		{ "random", "sec", AU_millicles, 10, 36, "-2.661411839725",
		  "1.000139831338152646161510375660880907" },
		{ "random", "csc", AU_degrees, 12, 36, "-E6481439733E21.6",
		  "1.466274E94782596E9592357704825962E573" },
		{ "random", "csc", AU_grads, 7, 36, "-1161.16322301",
		  "-1.620604652426646412606140113151164562" },
		{ "random", "csc", AU_millicles, 10, 36, "-3.206598116729",
		  "-49.636937769342196780692917270769792325" },
		{ "random", "asin", AU_degrees, 2, 36, "-0.00010101",
		  "-100.101101001001000010011101001111001001" },
		{ "random", "asin", AU_grads, 16, 36, "0.6F3252E5A",
		  "1C.9AEE6A6E67367350E6D09DFD33835AECACAA" },
		{ "random", "asin", AU_millicles, 3, 36, "-0.001211111",
		  "-101.210220021121112222002111002110202201" },
		{ "random", "acos", AU_degrees, 7, 36, "0.0252642122",
		  "152.516451445501510401650020034136244411" },
		{ "random", "acos", AU_grads, 3, 36, "-0.100110211",
		  "11112.211021220200202101221000022012202212" },
		{ "random", "acos", AU_millicles, 16, 36, "0.066ABFED9E",
		  "F6.0298DFF17DAAB5B2BA216E5B31200E558D11" },
		{ "random", "atan", AU_degrees, 7, 36, "5453015",
		  "155.666653626321025011323264332523030542" },
		{ "random", "atan", AU_grads, 10, 36, "66667387672",
		  "99.999999999045080668977586539270139908" },
		{ "random", "atan", AU_millicles, 2, 36, "-111011.011",
		  "-11110111.010100011101101100010011101101100000" },
		{ "random", "acot", AU_degrees, 16, 36, "-C89.B0D9",
		  "-0.0491DFEC8A7A86CDC89ED68AACACEA3ED572" },
		{ "random", "acot", AU_grads, 7, 36, "-263202.6",
		  "-0.000305246344121632353461502642600316" },
		{ "random", "acot", AU_millicles, 12, 36, "-641.041",
		  "-0.211272E01E1325D0094D6E5E841694E18695" },
		{ "random", "asec", AU_degrees, 2, 36, "111.011",
		  "1010010.001101010000001011101001001001110011" },
		{ "random", "asec", AU_grads, 7, 36, "56354.052",
		  "201.665316240453361356455440152113542212" },
		{ "random", "asec", AU_millicles, 10, 36, "-92339.076",
		  "250.001723592545966198753843786318653091" },
		{ "random", "acsc", AU_degrees, 3, 36, "-2111.12",
		  "-0.211220022022121210111011200000001220" },
		{ "random", "acsc", AU_grads, 10, 36, "-72174.913",
		  "-0.000882051319393969641283037317385710" },
		{ "random", "acsc", AU_millicles, 7, 36, "31001.024",
		  "0.010143011206640362431411055244216245" },
		{ "exact", "sin", AU_degrees, 10, 36, "30",
		  "0.5" },
		{ "exact", "sin", AU_degrees, 10, 36, "150",
		  "0.5" },
		{ "exact", "sin", AU_degrees, 10, 36, "-390",
		  "-0.5" },
		{ "exact", "cos", AU_degrees, 10, 36, "60",
		  "0.5" },
		{ "exact", "cos", AU_degrees, 10, 36, "180",
		  "-1" },
		{ "exact", "tan", AU_degrees, 10, 36, "45",
		  "1" },
		{ "exact", "tan", AU_degrees, 10, 36, "135",
		  "-1" },
		{ "exact", "cot", AU_degrees, 10, 36, "45",
		  "1" },
		{ "exact", "sec", AU_degrees, 10, 36, "60",
		  "2" },
		{ "exact", "sec", AU_degrees, 10, 36, "180",
		  "-1" },
		{ "exact", "csc", AU_degrees, 10, 36, "30",
		  "2" },
		{ "exact", "sin", AU_grads, 10, 36, "100",
		  "1" },
		{ "exact", "cos", AU_grads, 10, 36, "400",
		  "1" },
		{ "exact", "sin", AU_millicles, 10, 36, "250",
		  "1" },
		{ "exact", "tan", AU_millicles, 10, 36, "-125",
		  "-1" },
		{ "exact", "sin", AU_degrees, 10, 36, "3600000000000000000000000030",
		  "0.5" },
		{ "exact", "tan", AU_degrees, 10, 36, "90",
		  "POLE" },
		{ "exact", "cot", AU_degrees, 10, 36, "0",
		  "POLE" },
		{ "exact", "cot", AU_degrees, 10, 36, "180",
		  "POLE" },
		{ "exact", "sec", AU_degrees, 10, 36, "90",
		  "POLE" },
		{ "exact", "csc", AU_degrees, 10, 36, "0",
		  "POLE" },
		{ "exact", "csc", AU_grads, 10, 36, "200",
		  "POLE" },
		{ "exact", "tan", AU_millicles, 10, 36, "750",
		  "POLE" },
		{ "exact", "sin", AU_degrees, 3, 36, "1010",
		  "0.111111111111111111111111111111111112" },
		{ "exact", "sin", AU_degrees, 3, 20, "1010",
		  "0.11111111111111111112" },
		{ "exact", "sin", AU_degrees, 5, 36, "110",
		  "0.222222222222222222222222222222222222" },
		{ "exact", "sin", AU_degrees, 5, 20, "110",
		  "0.22222222222222222222" },
		{ "exact", "sin", AU_degrees, 7, 36, "42",
		  "0.333333333333333333333333333333333334" },
		{ "exact", "sin", AU_degrees, 7, 20, "42",
		  "0.33333333333333333334" },
		{ "inverse-exact", "asin", AU_degrees, 10, 36, "0.5",
		  "30" },
		{ "inverse-exact", "asin", AU_grads, 10, 36, "0.5",
		  "33.333333333333333333333333333333333333" },
		{ "inverse-exact", "asin", AU_millicles, 10, 36, "0.5",
		  "83.333333333333333333333333333333333333" },
		{ "inverse-exact", "acos", AU_degrees, 10, 36, "-0.5",
		  "120" },
		{ "inverse-exact", "acos", AU_grads, 10, 36, "-0.5",
		  "133.333333333333333333333333333333333333" },
		{ "inverse-exact", "acos", AU_millicles, 10, 36, "-0.5",
		  "333.333333333333333333333333333333333333" },
		{ "inverse-exact", "atan", AU_degrees, 10, 36, "1",
		  "45" },
		{ "inverse-exact", "atan", AU_grads, 10, 36, "1",
		  "50" },
		{ "inverse-exact", "atan", AU_millicles, 10, 36, "1",
		  "125" },
		{ "inverse-exact", "acot", AU_degrees, 10, 36, "-1",
		  "-45" },
		{ "inverse-exact", "acot", AU_grads, 10, 36, "-1",
		  "-50" },
		{ "inverse-exact", "acot", AU_millicles, 10, 36, "-1",
		  "-125" },
		{ "inverse-exact", "asec", AU_degrees, 10, 36, "-2",
		  "120" },
		{ "inverse-exact", "asec", AU_grads, 10, 36, "-2",
		  "133.333333333333333333333333333333333333" },
		{ "inverse-exact", "asec", AU_millicles, 10, 36, "-2",
		  "333.333333333333333333333333333333333333" },
		{ "inverse-exact", "acsc", AU_degrees, 10, 36, "2",
		  "30" },
		{ "inverse-exact", "acsc", AU_grads, 10, 36, "2",
		  "33.333333333333333333333333333333333333" },
		{ "inverse-exact", "acsc", AU_millicles, 10, 36, "2",
		  "83.333333333333333333333333333333333333" },
		{ "inverse-exact", "acos", AU_degrees, 10, 36, "0",
		  "90" },
		{ "inverse-exact", "acos", AU_grads, 10, 36, "0",
		  "100" },
		{ "inverse-exact", "acos", AU_millicles, 10, 36, "0",
		  "250" },
		{ "inverse-exact", "asin", AU_degrees, 10, 36, "-1",
		  "-90" },
		{ "inverse-exact", "asin", AU_grads, 10, 36, "-1",
		  "-100" },
		{ "inverse-exact", "asin", AU_millicles, 10, 36, "-1",
		  "-250" },
		{ "thirds", "asec", AU_grads, 2, 1, "-10",
		  "10000101.1" },
		{ "thirds", "acos", AU_grads, 2, 2, "0.1",
		  "1000010.11" },
		{ "thirds", "acsc", AU_grads, 2, 3, "10",
		  "100001.011" },
		{ "thirds", "acos", AU_grads, 2, 5, "0.1",
		  "1000010.10101" },
		{ "thirds", "asin", AU_grads, 2, 8, "0.1",
		  "100001.01010101" },
		{ "thirds", "acos", AU_grads, 2, 13, "0.1",
		  "1000010.1010101010101" },
		{ "thirds", "acsc", AU_millicles, 2, 1, "10",
		  "1010011.1" },
		{ "thirds", "acsc", AU_millicles, 2, 2, "10",
		  "1010011.01" },
		{ "thirds", "asec", AU_millicles, 2, 3, "-10",
		  "101001101.011" },
		{ "thirds", "acos", AU_millicles, 2, 5, "0.1",
		  "10100110.10101" },
		{ "thirds", "acsc", AU_millicles, 2, 8, "10",
		  "1010011.01010101" },
		{ "thirds", "asin", AU_millicles, 2, 13, "0.1",
		  "1010011.0101010101011" },
		{ "thirds", "acsc", AU_grads, 4, 1, "2",
		  "201.1" },
		{ "thirds", "asin", AU_grads, 4, 2, "0.2",
		  "201.11" },
		{ "thirds", "asec", AU_grads, 4, 3, "-2",
		  "2011.111" },
		{ "thirds", "asec", AU_grads, 4, 5, "-2",
		  "2011.11111" },
		{ "thirds", "asec", AU_grads, 4, 8, "-2",
		  "2011.11111111" },
		{ "thirds", "asec", AU_grads, 4, 13, "-2",
		  "2011.1111111111111" },
		{ "thirds", "asin", AU_millicles, 4, 1, "0.2",
		  "1103.1" },
		{ "thirds", "acos", AU_millicles, 4, 2, "0.2",
		  "2212.23" },
		{ "thirds", "asin", AU_millicles, 4, 3, "0.2",
		  "1103.111" },
		{ "thirds", "acsc", AU_millicles, 4, 5, "2",
		  "1103.11111" },
		{ "thirds", "asec", AU_millicles, 4, 8, "-2",
		  "11031.11111111" },
		{ "thirds", "asin", AU_millicles, 4, 13, "0.2",
		  "1103.1111111111111" },
		{ "thirds", "acsc", AU_grads, 5, 2, "2",
		  "113.13" },
		{ "thirds", "asec", AU_grads, 5, 3, "-2",
		  "1013.132" },
		{ "thirds", "asec", AU_grads, 5, 5, "-2",
		  "1013.13132" },
		{ "thirds", "acsc", AU_grads, 5, 13, "2",
		  "113.1313131313132" },
		{ "thirds", "acsc", AU_millicles, 5, 1, "2",
		  "313.2" },
		{ "thirds", "asec", AU_millicles, 5, 2, "-2",
		  "2313.13" },
		{ "thirds", "acsc", AU_millicles, 5, 5, "2",
		  "313.13132" },
		{ "thirds", "acsc", AU_millicles, 5, 8, "2",
		  "313.13131313" },
		{ "thirds", "acsc", AU_millicles, 5, 13, "2",
		  "313.1313131313132" },
		{ "thirds", "asin", AU_grads, 8, 1, "0.4",
		  "41.3" },
		{ "thirds", "asec", AU_grads, 8, 2, "-2",
		  "205.25" },
		{ "thirds", "asin", AU_grads, 8, 3, "0.4",
		  "41.253" },
		{ "thirds", "asin", AU_grads, 8, 5, "0.4",
		  "41.25253" },
		{ "thirds", "asin", AU_grads, 8, 8, "0.4",
		  "41.25252525" },
		{ "thirds", "acos", AU_grads, 8, 13, "0.4",
		  "102.5252525252525" },
		{ "thirds", "acsc", AU_millicles, 8, 1, "2",
		  "123.3" },
		{ "thirds", "asin", AU_millicles, 8, 2, "0.4",
		  "123.25" },
		{ "thirds", "asin", AU_millicles, 8, 3, "0.4",
		  "123.253" },
		{ "thirds", "asec", AU_millicles, 8, 5, "-2",
		  "515.25253" },
		{ "thirds", "acos", AU_millicles, 8, 8, "0.4",
		  "246.52525253" },
		{ "thirds", "acos", AU_millicles, 8, 13, "0.4",
		  "246.5252525252525" },
		{ "near", "cot", AU_millicles, 16, 36, "1AC.88ABCFC705714F51C6A6CA061F2A940B16E199980E115E202770EA1CBD5F288469AD8E319FCFDBED199D92",
		  "-2.1344CC1C64B0FFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "near", "acsc", AU_degrees, 10, 36, "1.12353561864495415498024547314142628829977272372978261560899814927044317560711126395461",
		  "62.879208649986537693621357902884483337" },
		{ "near", "sin", AU_degrees, 7, 36, "34.65216044462246331223361250114606554026444252336533153025364665366624630122510460141024",
		  "0.303114534311503153406115546053110645" },
		{ "near", "acos", AU_degrees, 7, 36, "0.63113151323023213544302243401164210130314243455342403526063235123636014606045504224535",
		  "31.536622420015021005611421255225056066" },
		{ "near", "tan", AU_grads, 16, 36, "-4C.8EADA59B78FEF8393A82110E65CFEFEB9BA98609E17A33D87574F4B14086F4F3CF2305F23058325886E6FE",
		  "-2.977EAAA30D4C3FFFFFFFFFFFFFFFFFFFFFFF" },
		{ "near", "atan", AU_grads, 16, 36, "1.1C3E0CB173C56CB13FE65246695CB177FC41888ECE5F35F2D03294E16DC3376D44EC2842574F78EDBCA354",
		  "35.533632332320000000000000000000000000" },
		{ "near", "acot", AU_grads, 16, 36, "2.4B49F9851AB8FDFA8DA8D608937197B78A6EB7F276681BB40EE7D76786D7DF8E9E801DB7D504EAC8A558C2",
		  "1A.2B5CEA06497C000000000000000000000000" },
		{ "near", "atan", AU_millicles, 16, 36, "-0.8CD9FCB52006CE715E1772107916336503313115CBD5F8E90D1FC28E4C539B26C6A73C3C8A2ACA61E659DA",
		  "-50.0DF82A5E2E40000000000000000000000000" },
		{ "near", "acsc", AU_degrees, 7, 36, "2.13335040244544626041304206246501406664501523525300405150602211220116133423026434335425",
		  "35.563112461340331250225226322241504012" },
		{ "near", "asec", AU_grads, 7, 36, "1.40665001221136664316414334631464062040104351461205023205204031444500024112435420310563",
		  "110.516445401065605355651510462003252105" },
		{ "near", "cos", AU_grads, 10, 36, "41.88754858163932386285545538101685514994187623035458259554620858476560240024703398077422",
		  "0.791236407574867128822404538368573413" },
		{ "near", "acos", AU_grads, 16, 36, "0.5C140EBC5DA525B25953F659833BBE08635D39AAE936342E9365CCD0ECB316B405798BC66A42C31E93646",
		  "4C.93BCD7CD5D74000000000000000000000001" },
		{ "near", "acot", AU_degrees, 7, 36, "-3.56432541404522040263461241132203400316303123213144000220032665034556063204462361236621",
		  "-20.363256463033243453515051211630511111" },
		{ "near", "acos", AU_degrees, 16, 36, "0.DAE5CF89776756BA9D2C0F205961BE7A9A5BC1A4D483ACFCD062B6F4CFC05BC86AE2A0CAA0364DCC7A5941",
		  "1F.3B8D6F5C59C6000000000000000000000001" },
		{ "widened", "csc", AU_grads, 6, 50, "0.000000000000000000000000000004214424353334",
		  "22305143231131343210223002231342.54322451002045455112513400454515011024405205154440" },
		{ "widened", "sec", AU_degrees, 13, 36, "6C.00000000000000009",
		  "-649B6C8612163B5312.907C540B06C4C5358A32CA458666C9CA9443" },
		{ "widened", "sec", AU_millicles, 10, 50, "750.00000000000000000000000000000000003",
		  "5305164769729844525629458779083812067.81532152468188162492224480196322658783896709620636" },
		{ "widened", "sec", AU_grads, 15, 50, "69.EEEEEEEEEEEEEEEEEEEEEEEEEE",
		  "439DE28E6E0B68257AEDC9E6D20A.A6632982C75736592C06C9C10D5D8CB19716CD3DC460C356E9" },
		{ "widened", "csc", AU_grads, 15, 36, "0.00000000000000000000000000000000004",
		  "10DAEBD7573E19815A3AE639343DDCCC09D46.D15855246B56DA7B717000B7A07572E5EAD5" },
		{ "widened", "csc", AU_grads, 8, 20, "310.000000000000000000000000000000000001",
		  "-77522735270551123021063177274470036212.32372533120233755612" },
		{ "widened", "tan", AU_grads, 16, 36, "-0.0000000000000000000000E447324943126B9C4",
		  "-0.00000000000000000000000395F5F3C5FE105A428AF932976FD972730E0" },
		{ "widened", "cot", AU_degrees, 2, 36, "-0.00000000000000000001",
		  "-11100101001011101110000011.010011000111100000111110110111110011" },
		{ "widened", "tan", AU_degrees, 8, 36, "-0.0000000000000000000006312",
		  "-0.00000000000000000000000711105752577221657010101073740642157" },
		{ "widened", "tan", AU_degrees, 13, 36, "0.000000000000A9786CA00545543B2957",
		  "0.000000000000025896A5291CBA31041B975B8605008776B51" },
		{ "widened", "cos", AU_grads, 2, 36, "-1011100001.001011000001110001101100110110011011010001011011011011011000111000001101010011011001110001101000101110101",
		  "0.100011010010000011100010001111101101" },
		{ "widened", "sin", AU_millicles, 2, 36, "1010000011110.1001110100101101111000100100111111010101011000001110000111100010100100001111110001111100110001010011100111",
		  "0.110011111010111111110010101011000001" },
		{ "widened", "sec", AU_millicles, 2, 36, "10101011000.1001101111111010001011001111011011100100010011100001001001101100110010101100000011111011000000000101000101",
		  "-1.011110010111111100110101100111111010" },
		{ "widened", "tan", AU_degrees, 2, 36, "11010111101.1100000010111101001101110101000100110110110000010001111011101010101100100100011011001010000101110111101011",
		  "-11.100010111000100001010000011111100100" },
		{ "widened", "sec", AU_grads, 2, 36, "1110011.100000010000110010010010100100111101101001100100101010010100000111011111100010010001100000111000110101001",
		  "-100.001001011010001010000010110000011100" },
		{ "widened", "cos", AU_degrees, 2, 36, "-1001110000.0100101110001101010011000010111101011101010100011100100110100101101001011001000001000000101110001110101",
		  "-0.000110010111001010010001011000010111" },
		{ "decider", "sin", AU_degrees, 10, 36, "0.044238684496619981589837777846696530305821606168056352173079106625079676216568544612116886984755",
		  "0.000772110624500723781846249822300641" },
		{ "decider", "tan", AU_grads, 10, 36, "0.023141815931056658461884012467635329816492403660138069755866688515106977144238488650359068735932",
		  "0.000363510810610125691498250123536478633" },
		{ "decider", "asin", AU_degrees, 10, 36, "0.000009149289244305977016429121454726125698062560489017647226812881818137116908206042038932306282",
		  "0.000524215659250484465256358923568314" },
		{ "decider", "atan", AU_degrees, 10, 36, "0.000011613530832749731187911176842213148916250184017333814332257239827496732589905428514915562375",
		  "0.000665406301931696555346240629091881" },
	};

	long	l_retakesBefore = lzi_retake_count();

	for ( const char* sz_family : { "random", "exact", "inverse-exact", "thirds", "near", "widened", "decider" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			AngleFunction	fn = nullptr;

			for ( auto& f : kFunctions )
				if ( ::strcmp( f.m_szName, pin.m_szName ) == 0 )
					fn = f.m_fn;

			CLZNumber	x;
			CLZNumber	answer;

			glWorkingPrecision = pin.m_lPlaces;
			x.SetBase( pin.m_lBase, false );
			x.Assign( pin.m_szX );

			geAngleUnit = pin.m_eUnit;
			long	l_code = fn( &x, &answer );
			geAngleUnit = AU_radians;

			bool	b_ok = (answer.m_lBase == pin.m_lBase);

			// A pole answers as cot(0) does in radians...
			if ( ::strcmp( pin.m_szAnswer, "POLE" ) == 0 )
				b_ok = b_ok && (l_code == kDivideByZero) && answer.FlagIsOn( kfInfinity ) && answer.FlagIsOn( kfDiv0 );
			else
				b_ok = b_ok && (l_code == kNoError) && SpelledAs( &answer, pin.m_szAnswer );

			l_cases++;

			if ( !b_ok && (l_wrong++ == 0) )
				s_first = std::string( pin.m_szName ) + " of " + pin.m_szX + " in base " + std::to_string( pin.m_lBase ) +
						  ", unit " + std::to_string( static_cast<int>( pin.m_eUnit ) );
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "in a unit, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), l_wrong == 0, szDetail );
	}

	// The near halves are in doubt by construction, so they are taken again
	//	-- in the unit, or they would not have come out right...
	{
		long	l_retakes = lzi_retake_count() - l_retakesBefore;

		::snprintf( szDetail, sizeof( szDetail ), "%ld retakes over the pins", l_retakes );
		CheckTrue( "a retake in a unit is taken in it", l_retakes >= 14, szDetail );
	}

	// atan2 and the argument, the axes and diagonals exact...
	{
		struct { long m_lBase; AngleUnit m_eUnit; const char* m_szY; const char* m_szX; const char* m_szAngle; } kAngles[] =
		{
			{ 10, AU_degrees, "0", "-3",
			  "180" },
			{ 7, AU_degrees, "0", "-3",
			  "345" },
			{ 10, AU_degrees, "2", "0",
			  "90" },
			{ 7, AU_degrees, "2", "0",
			  "156" },
			{ 10, AU_grads, "-2", "0",
			  "-100" },
			{ 7, AU_grads, "-2", "0",
			  "-202" },
			{ 10, AU_degrees, "5", "5",
			  "45" },
			{ 7, AU_degrees, "5", "5",
			  "63" },
			{ 10, AU_degrees, "-5", "-5",
			  "-135" },
			{ 7, AU_degrees, "-5", "-5",
			  "-252" },
			{ 10, AU_millicles, "5", "-5",
			  "375" },
			{ 7, AU_millicles, "5", "-5",
			  "1044" },
			{ 10, AU_grads, "0", "7",
			  "0" },
			{ 7, AU_grads, "0", "10",
			  "0" },
			{ 10, AU_degrees, "1", "2",
			  "26.565051177077989351572193720453294671" },
			{ 7, AU_degrees, "1", "2",
			  "35.364545464064125306603641504363534423" },
			{ 7, AU_millicles, "-0.3", "0.5",
			  "-152.003402434163206133012451610142623442" },
			{ 10, AU_grads, "-1", "-3",
			  "-179.516723530086654835080215244948105190" },
			{ 7, AU_grads, "-1", "-3",
			  "-344.342144002155305163512553321646354234" },
			// ...and a tiny quotient whose angle in degrees is a hair above a half,
			//	which atan2's radian decider, left on, takes down...
			{ 10, AU_degrees, "0.00000000000000000000000000000014897428", "1",
			  "0.000000000000000000000000000008535598" },
		};

		long		l_wrong = 0;
		std::string	s_first;

		glWorkingPrecision = 36;

		for ( auto& pin : kAngles )
		{
			CLZNumber	y;
			CLZNumber	x;
			CLZNumber	angle;
			CLZNumber	arg;
			CLZComplex	z;

			y.SetBase( pin.m_lBase, false );
			y.Assign( pin.m_szY );
			x.SetBase( pin.m_lBase, false );
			x.Assign( pin.m_szX );
			SetComplex( z, pin.m_lBase, pin.m_szX, pin.m_szY );

			geAngleUnit = pin.m_eUnit;
			long	l_code = lz_arctangent2( &y, &x, &angle );
			long	l_argCode = lzc_arg( &z, &arg );
			geAngleUnit = AU_radians;

			if ( (l_code != kNoError) || (l_argCode != kNoError) || !SpelledAs( &angle, pin.m_szAngle ) ||
				 !SpelledAs( &arg, pin.m_szAngle ) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::string( pin.m_szY ) + ", " + pin.m_szX + " in base " + std::to_string( pin.m_lBase );
			}
		}

		size_t	n_pins = sizeof( kAngles ) / sizeof( kAngles[ 0 ] );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %zu right, atan2 and lzc_arg both", n_pins );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong, first %s", l_wrong, n_pins, s_first.c_str() );

		CheckTrue( "in a unit, atan2 and the argument", l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 36;

	// The unit is the host's call only: gamma's reflection takes a sine,
	//	and that sine stays in radians...
	{
		CLZNumber	x;
		CLZNumber	inDegrees;
		CLZNumber	inRadians;

		x.Assign( "-2.5" );
		geAngleUnit = AU_degrees;
		lz_gamma( &x, &inDegrees );
		geAngleUnit = AU_radians;
		lz_gamma( &x, &inRadians );

		CheckTrue( "the engine's own sine stays radians", SameAnswer( &inDegrees, &inRadians ),
				   "gamma(-2.5) the same with degrees set" );
	}

	// A value naming no unit is radians, and radians are as they were...
	{
		CLZNumber	x;
		CLZNumber	named;
		CLZNumber	plain;

		x.Assign( "30" );
		geAngleUnit = static_cast<AngleUnit>( 7 );
		lz_sine( &x, &named );
		geAngleUnit = AU_radians;
		lz_sine( &x, &plain );

		CheckTrue( "a unit out of range reads as radians", SameAnswer( &named, &plain ), "sin(30) both ways" );
	}

	// In place, over the operand...
	{
		CLZNumber	x;
		CLZNumber	y;
		CLZNumber	xAgain;

		x.Assign( "30" );
		geAngleUnit = AU_degrees;
		lz_sine( &x );
		y.Assign( "1" );
		xAgain.Assign( "-1" );
		lz_arctangent2( &y, &xAgain, &xAgain );
		geAngleUnit = AU_radians;

		bool	b_ok = SpelledAs( &x, "0.5" ) && SpelledAs( &xAgain, "135" );

		CheckTrue( "in a unit, in place", b_ok, b_ok ? "sin 30 over itself, atan2 over x" : "an answer wrong" );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**The two pairs, stage 3** -- 2026-10-22, on Rich's calls: the names
	lz_sine_cosine and lz_hypersine_hypercosine, and both targets required
	and different. Each answer must be the single function's in every
	respect a caller can see -- code, base, flags, precision and spelling
	-- so the reference is the engine itself: lz_sine and lz_cosine, or
	lz_hypersine and lz_hypercosine, of the same operand. The arguments
	are random in every base, unit and width, whole parts of a turn in a
	unit, and tiny ones whose cosine or sine lies a hair from a half, which
	the pair must hand to its function to be taken again wider.
*/

// The same answer in every respect a caller can see...
static bool SameInEveryRespect( CLZNumber* pGot, long lGotCode, CLZNumber* pWant, long lWantCode )
{
	return (lGotCode == lWantCode) && SameAnswer( pGot, pWant ) && (pGot->m_lPrecision == pWant->m_lPrecision);
}

// lDigits random digits of lBase, lPlaces of them after the point...
static std::string RandomDigitsIn( std::uint64_t& uState, long lBase, long lDigits, long lPlaces, bool bNegative )
{
	const char*	sz_digits = (lBase == 12) ? "0123456789DE" : "0123456789ABCDEF";
	std::string	s = bNegative ? "-" : "";

	for ( long i = 0; i < lDigits; i++ )
	{
		uState = (uState * 6364136223846793005ULL) + 1442695040888963407ULL;

		long	l_digit = static_cast<long>( (uState >> 33) % static_cast<std::uint64_t>( lBase ) );

		if ( (i == 0) && (l_digit == 0) )
			l_digit = 1;

		if ( i == (lDigits - lPlaces) )
			s += (i == 0) ? "0." : ".";

		s += sz_digits[ l_digit ];
	}

	return s;
}

static long PairAgainstSingles( bool bHyperbolic, CLZNumber* pX, long lUnit, long lPlaces, int iInPlace )
{
	glWorkingPrecision = lPlaces;
	geAngleUnit = static_cast<AngleUnit>( lUnit );

	CLZNumber	x;
	CLZNumber	xs;
	CLZNumber	xc;
	CLZNumber	s;
	CLZNumber	c;
	CLZNumber	s1;
	CLZNumber	c1;

	x.Assign( pX );
	xs.Assign( pX );
	xc.Assign( pX );

	long	l_code = 0;

	if ( iInPlace == 1 )
	{
		l_code = bHyperbolic ? lz_hypersine_hypercosine( &x, &x, &c ) : lz_sine_cosine( &x, &x, &c );
		s.Assign( &x );
	}
	else if ( iInPlace == 2 )
	{
		l_code = bHyperbolic ? lz_hypersine_hypercosine( &x, &s, &x ) : lz_sine_cosine( &x, &s, &x );
		c.Assign( &x );
	}
	else
		l_code = bHyperbolic ? lz_hypersine_hypercosine( &x, &s, &c ) : lz_sine_cosine( &x, &s, &c );

	long	l_sine = bHyperbolic ? lz_hypersine( &xs, &s1 ) : lz_sine( &xs, &s1 );
	long	l_cosine = bHyperbolic ? lz_hypercosine( &xc, &c1 ) : lz_cosine( &xc, &c1 );
	long	l_want = (l_sine != kNoError) ? l_sine : l_cosine;

	geAngleUnit = AU_radians;

	return (SameInEveryRespect( &s, l_code, &s1, l_want ) ? 0 : 1) + (SameInEveryRespect( &c, l_code, &c1, l_want ) ? 0 : 1);
}

static void TestSineCosinePairs()
{
	Banner( "complex numbers, stage 3: the sine-cosine and sinh-cosh pairs" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 200 ];

	static const long	kWidths[] = { 0, 1, 3, 8, 20, 36, 36, 50 };
	static const long	kTurns[] = { 0, 360, 400, 1000 };

	// Random arguments, every base, unit and width; whole 24ths of a turn;
	//	and tiny ones, d / b^k with k about half the width, whose cosine
	//	lies a hair from a half in an odd base...
	for ( bool b_hyperbolic : { false, true } )
	{
		std::uint64_t	u_state = b_hyperbolic ? 20261022 : 2210;
		long			l_cases = 0;
		long			l_wrong = 0;
		long			l_handed = lzi_pair_handoffs();

		for ( long i = 0; i < 240; i++ )
		{
			long	l_base = 2 + (i % 15);
			long	l_places = kWidths[ (i / 15) % 8 ];
			long	l_unit = b_hyperbolic ? 0 : ((i / 3) % 4);
			long	l_digits = 1 + static_cast<long>( (u_state >> 40) % static_cast<std::uint64_t>( l_places + 10 ) );
			long	l_after = b_hyperbolic ? (l_digits - 1) : static_cast<long>( (u_state >> 20) % static_cast<std::uint64_t>( l_digits + 1 ) );
			CLZNumber	x;

			x.SetBase( l_base, false );
			x.Assign( RandomDigitsIn( u_state, l_base, l_digits, l_after, (i % 2) == 1 ).c_str() );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, l_unit, l_places, 0 );
			l_cases += 2;

			if ( !b_hyperbolic && (l_unit != 0) )
			{
				// n turn / 24 in the unit, exact where it ends in the base...
				CLZNumber	n;
				CLZNumber	d;

				glWorkingPrecision = 80;
				n.SetBase( l_base, false );
				n.Assign( ((i % 97) - 48) * kTurns[ l_unit ] );
				d.SetBase( l_base, false );
				d.Assign( 24L );
				x.SetBase( l_base, false );
				lz_divide( &n, &d, &x );
				l_wrong += PairAgainstSingles( false, &x, l_unit, l_places, 0 );
				l_cases += 2;
			}

			// A tiny argument, d / b^(P / 2): its square's half, d^2 / 2 b^P,
			//	ends on a half at P in an odd base, the next term just off it...
			if ( (l_places > 0) && ((l_places % 2) == 0) && ((l_base % 2) == 1) )
			{
				CLZNumber	one;

				x.SetBase( l_base, false );
				x.Assign( 1L + (i % (l_base - 1)) );
				one.SetBase( l_base, false );
				one.Assign( 1L );

				for ( long k = 0; k < (l_places / 2); k++ )
				{
					CLZNumber	b;

					b.SetBase( l_base, false );
					b.Assign( l_base );
					glWorkingPrecision = l_places + 20;
					lz_divide( &x, &b );
				}

				l_wrong += PairAgainstSingles( b_hyperbolic, &x, 0, l_places, 0 );
				l_cases += 2;
			}
		}

		long	l_handedOn = lzi_pair_handoffs() - l_handed;

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld answers differ, %ld handed on", l_wrong, l_cases, l_handedOn );
		CheckTrue( b_hyperbolic ? "sinh, cosh as their functions" : "sin, cos as their functions",
				   (l_wrong == 0) && (l_handedOn > 0), szDetail );
	}

	// A not-a-number, an infinity, nought, and either target the operand...
	{
		long	l_wrong = 0;

		for ( bool b_hyperbolic : { false, true } )
		{
			CLZNumber	x;

			MakeNotANumber( x );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, 0, 36, 0 );
			MakeInfinite( x );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, 0, 36, 0 );
			x.Assign( 0L );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, b_hyperbolic ? 0 : 1, 36, 0 );
			x.Assign( "0.7" );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, 0, 36, 1 );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, 0, 36, 2 );
			x.Assign( "750" );
			l_wrong += PairAgainstSingles( b_hyperbolic, &x, b_hyperbolic ? 0 : 1, 20, 1 );
		}

		// 60 degrees over the sine: the sine written first, and the cosine,
		//	exactly a half, then its function's -- of the angle as given, not
		//	of the sine now over it...
		CLZNumber	x;

		x.Assign( "60" );
		l_wrong += PairAgainstSingles( false, &x, 1, 36, 1 );

		::snprintf( szDetail, sizeof( szDetail ), "%ld answers differ", l_wrong );
		CheckTrue( "flagged, nought, and in place", l_wrong == 0, szDetail );
	}

	// Past the kernel's floor, where one run gives both...
	{
		long	l_wrong = 0;
		long	l_runs = lzi_sincos_splittings();

		for ( long l_base : { 14L, 5L } )
		{
			CLZNumber	x;

			x.SetBase( l_base, false );
			x.Assign( "2.31" );
			l_wrong += PairAgainstSingles( false, &x, 0, lzi_sincos_split_floor( l_base ) + 3, 0 );
			l_wrong += PairAgainstSingles( false, &x, 1, lzi_sincos_split_floor( l_base ) + 3, 0 );
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld answers differ, the kernel run %ld times", l_wrong,
					lzi_sincos_splittings() - l_runs );
		CheckTrue( "past the kernel's floor", (l_wrong == 0) && ((lzi_sincos_splittings() - l_runs) > 0), szDetail );
	}

	// **An ordinary answer stands**, and one in doubt is handed on: a pair
	//	handing every answer to its function would answer the same...
	{
		std::uint64_t	u_state = 77;
		long			l_handed = lzi_pair_handoffs();

		glWorkingPrecision = 36;

		for ( long i = 0; i < 100; i++ )
		{
			CLZNumber	x;
			CLZNumber	s;
			CLZNumber	c;

			x.Assign( RandomDigitsIn( u_state, 10, 12, 10, false ).c_str() );
			lz_sine_cosine( &x, &s, &c );
			lz_hypersine_hypercosine( &x, &s, &c );
		}

		long	l_ordinary = lzi_pair_handoffs() - l_handed;

		// cos 3^-6 at 12 places in base three: x^2 / 2 on a half, the next
		//	term just under it...
		CLZNumber	x;
		CLZNumber	s;
		CLZNumber	c;
		CLZNumber	want;

		x.SetBase( 3, false );
		x.Assign( "0.000001" );
		glWorkingPrecision = 12;
		l_handed = lzi_pair_handoffs();
		lz_sine_cosine( &x, &s, &c );

		long	l_near = lzi_pair_handoffs() - l_handed;

		lz_cosine( &x, &want );
		glWorkingPrecision = 36;

		bool	b_ok = (l_ordinary <= 8) && (l_near >= 1) && SameAnswer( &c, &want );

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 400 ordinary handed on, the near half %s", l_ordinary,
					(l_near >= 1) ? "handed on" : "kept" );
		CheckTrue( "ordinary answers stand", b_ok, szDetail );
	}

	// Both targets, and two of them...
	{
		CLZNumber	x;
		CLZNumber	s;
		CLZNumber	c;
		bool		b_ok = true;

		x.Assign( "0.5" );
		s.Assign( "1" );
		b_ok = b_ok && (lz_sine_cosine( &x, &s, nullptr ) == kFailed) && s.FlagIsOn( kfQNaN );
		c.Assign( "1" );
		b_ok = b_ok && (lz_hypersine_hypercosine( &x, nullptr, &c ) == kFailed) && c.FlagIsOn( kfQNaN );
		s.Assign( "1" );
		b_ok = b_ok && (lz_sine_cosine( &x, &s, &s ) == kFailed) && s.FlagIsOn( kfQNaN );
		s.Assign( "1" );
		c.Assign( "1" );
		b_ok = b_ok && (lz_hypersine_hypercosine( nullptr, &s, &c ) == kFailed) && s.FlagIsOn( kfQNaN ) && c.FlagIsOn( kfQNaN );

		CheckTrue( "a missing or shared target refused", b_ok, b_ok ? "kFailed, each target a not-a-number" : "not refused" );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**lzc_polar, stage 3** -- 2026-10-22: r cos theta + i r sin theta, the
	angle in the host's unit, each part correctly rounded half-even at the
	places asked. Every expected part is generated -- exact Fractions where
	the part is rational, Python's decimal at 320 digits otherwise -- never
	typed. The families: random; a hair from a half inside the retake's
	reach, each a case the guarded run alone gets wrong; exact parts,
	among them halves of r that tie in an odd base and go by the kept digit
	of the magnitude; the decider's, r or r theta exactly on a half with a
	tiny angle; an angle in another base, one an exact twelfth of a turn
	that ends only in its own; and large and small moduli and angles.
*/

static void TestPolar()
{
	Banner( "complex numbers, stage 3: lzc_polar" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	struct { const char* m_szFamily; long m_lUnit; long m_lPlaces; long m_lBase; const char* m_szR; long m_lBaseAngle;
			 const char* m_szAngle; const char* m_szReal; const char* m_szImag; } kPins[] =
	{
		{ "random", 0, 36, 10, "73761428830.105597", 10, "0.00000000000000065606934043999966469",
		  "73761428830.105596999999999999984125545493803666",
		  "0.000048392611962479355106397741657289" },
		{ "random", 1, 36, 10, "-0.000000000000000000000007946061383875303", 10, "-0.0000000000000000000350796331",
		  "-0.000000000000000000000007946061383875",
		  "0.000000000000000000000000000000000000" },
		{ "random", 0, 36, 2, "-0.000000000000000110011", 2, "1001110001.00111",
		  "0.000000000000000110010111101001001011",
		  "0.000000000000000000010001000011011110" },
		{ "random", 2, 36, 3, "-0.00000000000000000000000000022121201", 3, "-21222012010.02021",
		  "-0.000000000000000000000000000101112210",
		  "0.000000000000000000000000000212121101" },
		{ "random", 3, 36, 7, "0.000000005341236043342062501", 7, "-0.000000000015445432226600341300515",
		  "0.000000005341236043342062500666666661",
		  "-0.000000000000000000000303363455061125" },
		{ "random", 1, 36, 12, "D103162676D10", 12, "-0.000000000000000000000000007223498",
		  "D103162676D10.000000000000000000000000000000000000",
		  "-0.000000000000001320784850E5941179907D" },
		{ "random", 0, 36, 16, "-0.00000000000000000000000000DBDC", 16, "0.000000000000000000000000000000E9E2",
		  "-0.00000000000000000000000000DBDC000000",
		  "0.000000000000000000000000000000000000" },
		{ "random", 0, 20, 10, "0.00000000000000000009643522828629", 10, "-0.00003562387896281",
		  "0.00000000000000000010",
		  "0.00000000000000000000" },
		{ "random", 1, 50, 5, "3", 5, "0.00000000232303104121420103112",
		  "2.44444444444444444444424233441334244002044024114102",
		  "0.00000000003232301121043303134442442043033021330111" },
		{ "random", 3, 36, 16, "0.0000000000000000000000852D34B", 16, "0.005ABA03F51663C6F0D4",
		  "0.0000000000000000000000852D34AFEA5C8B",
		  "0.000000000000000000000000004BEAE80BA2" },
		{ "near", 3, 20, 2, "-1000.01111111000001000011101001011001010100000010001101101010100000101110100001101010000000001", 2, "10001100",
		  "-101.01101010011010000110",
		  "-110.10001011111000010000" },
		{ "near", 0, 36, 7, "2.1135162450231123162043202053212365652505205164236554663301505225553021621556626000421010462213062613225645", 7, "2.033551051353624512",
		  "-1.021361535055401403466515403130416464",
		  "1.622620515613154132213055534660333316" },
		{ "near", 3, 20, 16, "-12.26DEF3984E88124993CFA3A17D02167FFBB8F94B29E85049A87CD1FCC503684632968B5015DE67E33109CC3D72", 16, "-BA.7F999FD73D248",
		  "-7.0D430C41212AFFFFFFFF",
		  "10.B9DF8E45816687B0F986" },
		{ "near", 0, 20, 2, "110.010110110100010100100100001011110011000001010110110011101010100111110111111011100001100111", 2, "-1.010011000001011",
		  "1.10110111101010100011",
		  "-110.00011110101111111101" },
		{ "near", 3, 20, 12, "9.E5514E6169E549D413319E82993585E504068E0627E8902D6601978E71920D7906E8269D8788300E12224EE4D6", 12, "145.07E9D7D1D",
		  "3.30121269488E63421931",
		  "9.4DD123955D0D9E344D39" },
		{ "near", 1, 20, 10, "6.711181401382115523464609021642053306757293712908662300391955769491285606855046435956974636", 10, "-149.08454911574",
		  "-5.75769962206130969660",
		  "-3.44802129697197019378" },
		{ "near", 0, 20, 2, "-10.100110101000110100101010011111011010011100100100100101010001101001011010100001110111100011", 2, "1.11101",
		  "0.11011011011011010110",
		  "-10.01110101011001100000" },
		{ "near", 1, 36, 3, "-102.1120102210122001101201101200100000122021201201211002112220002122221122220010010111010220102222122102102", 3, "-11211.0201221210110011212",
		  "21.110000122012211121121211111010221110",
		  "22.210111012200002102221211120000001001" },
		{ "exact", 1, 36, 10, "1", 10, "180",
		  "-1",
		  "0" },
		{ "exact", 1, 36, 10, "2.5", 10, "-450",
		  "0",
		  "-2.5" },
		{ "exact", 2, 36, 10, "-7", 10, "100",
		  "0",
		  "-7" },
		{ "exact", 3, 36, 16, "F.F", 16, "5DC",
		  "-F.F",
		  "0" },
		{ "exact", 1, 36, 10, "3", 10, "60",
		  "1.5",
		  "2.598076211353315940291169512258808550" },
		{ "exact", 1, 0, 3, "1", 3, "1010",
		  "1",
		  "0" },
		{ "exact", 1, 3, 3, "1200.222", 3, "-112110",
		  "1110.211",
		  "-211.222" },
		{ "exact", 1, 0, 3, "-1220202", 3, "10111110",
		  "-1122211",
		  "-221212" },
		{ "exact", 2, 4, 5, "0.023", 5, "4000",
		  "0",
		  "0.023" },
		{ "exact", 3, 2, 7, "-1.01", 7, "21210",
		  "0",
		  "-1.01" },
		{ "exact", 0, 36, 10, "123456.789", 10, "0",
		  "123456.789",
		  "0" },
		{ "exact", 1, 2, 10, "12.345", 10, "30",
		  "10.69",
		  "6.17" },
		{ "decider", 0, 1, 4, "123103.32", 4, "-0.000000000000000000000000000001",
		  "123103.3",
		  "0.0" },
		{ "decider", 1, 1, 4, "123103.32", 4, "-0.00000000000233123233122333331",
		  "123103.3",
		  "0.0" },
		{ "decider", 0, 0, 10, "3.5", 10, "0.000000000000000000000000000001",
		  "3",
		  "0" },
		{ "decider", 2, 36, 16, "0.0000000000000000000000000000000000018", 16, "64.0000000000000000000000001",
		  "0.000000000000000000000000000000000000",
		  "0.000000000000000000000000000000000001" },
		{ "decider", 0, 10, 10, "15000000000", 10, "0.00000000000000000001",
		  "15000000000.0000000000",
		  "0.0000000001" },
		{ "decider", 0, 36, 2, "111000", 2, "-0.0000000000000000000000000000000000000001",
		  "111000.000000000000000000000000000000000000",
		  "-0.000000000000000000000000000000000011" },
		{ "mixed", 0, 36, 10, "2.125", 2, "0.1101",
		  "1.461331819718572794592382651774137233",
		  "1.542768392429014168021793440523759044" },
		{ "mixed", 1, 36, 10, "-0.75", 3, "1.01",
		  "-0.749858977836376049204603944675975386",
		  "-0.014543498828868279341067896089797706" },
		{ "mixed", 2, 36, 10, "2.25", 3, "1020.1",
		  "1.948557158514986955218377134194106413",
		  "1.125" },
		{ "mixed", 3, 36, 16, "1", 6, "1313.200000000000000000000000000001",
		  "-0.80000000000000000000079E14A54759DCAA",
		  "0.DDB3D742C265539D92BA12525B2F89BDA2AE" },
		{ "scale", 0, 36, 10, "98765432123456789", 10, "2.3",
		  "-65805039155199327.511994753327159227058161795706779530",
		  "73649897517347799.698487848639966328305632395586304118" },
		{ "scale", 1, 36, 12, "0.0000000000000000000000000000000000000005", 12, "0.1",
		  "0.000000000000000000000000000000000000",
		  "0.000000000000000000000000000000000000" },
		{ "scale", 0, 36, 10, "-1.1", 10, "10000000000000000000000007",
		  "0.015137476200383376121927023232850755",
		  "1.099895839074902049820622559301622134" },
		{ "scale", 0, 36, 7, "3", 7, "-1000000000000000000000000000002",
		  "2.500463201606261513451451022513635434",
		  "1.162525512000036363533322226553642161" },
	};

	for ( const char* sz_family : { "random", "near", "exact", "decider", "mixed", "scale" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZNumber	r;
			CLZNumber	angle;
			CLZComplex	z;

			glWorkingPrecision = pin.m_lPlaces;
			r.SetBase( pin.m_lBase, false );
			r.Assign( pin.m_szR );
			angle.SetBase( pin.m_lBaseAngle, false );
			angle.Assign( pin.m_szAngle );
			geAngleUnit = static_cast<AngleUnit>( pin.m_lUnit );

			long	l_code = lzc_polar( &r, &angle, &z );

			geAngleUnit = AU_radians;
			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &z.m_real, pin.m_szReal ) || !SpelledAs( &z.m_imag, pin.m_szImag ) ||
				 (z.m_real.m_lBase != pin.m_lBase) || (z.m_imag.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::string( "base " ) + std::to_string( pin.m_lBase ) + " r " + pin.m_szR + " angle " + pin.m_szAngle;
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "polar, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 36;

	// A not-a-number in either, an infinite angle, an infinite modulus, and
	//	nought...
	{
		CLZNumber	r;
		CLZNumber	angle;
		CLZNumber	nan;
		CLZNumber	infinite;
		CLZComplex	z;
		bool		b_ok = true;

		MakeNotANumber( nan );
		MakeInfinite( infinite );
		angle.Assign( "1.25" );
		r.Assign( "2" );

		b_ok = b_ok && (lzc_polar( &nan, &angle, &z ) == kNoError) && z.m_real.FlagIsOn( kfQNaN ) && z.m_imag.FlagIsOn( kfQNaN );
		b_ok = b_ok && (lzc_polar( &r, &nan, &z ) == kNoError) && z.m_real.FlagIsOn( kfQNaN ) && z.m_imag.FlagIsOn( kfQNaN );
		b_ok = b_ok && (lzc_polar( &r, &infinite, &z ) == kNoError) && z.m_real.FlagIsOn( kfQNaN ) && z.m_imag.FlagIsOn( kfDiv0 );
		b_ok = b_ok && (lzc_polar( &infinite, &angle, &z ) == kNoError) && z.m_real.FlagIsOn( kfInfinity ) &&
			   z.m_imag.FlagIsOn( kfInfinity ) && z.m_imag.FlagIsOff( kfQNaN );

		r.Assign( 0L );
		b_ok = b_ok && (lzc_polar( &r, &angle, &z ) == kNoError) && SpelledAs( &z.m_real, "0" ) && SpelledAs( &z.m_imag, "0" );

		CheckTrue( "flagged operands, and nought", b_ok, b_ok ? "a not-a-number, the point at infinity, an exact nought" : "an answer wrong" );
	}

	// A missing target or operand refused...
	{
		CLZNumber	r;
		CLZNumber	angle;
		CLZComplex	z;

		r.Assign( "2" );
		angle.Assign( "1" );

		bool	b_ok = (lzc_polar( &r, &angle, nullptr ) == kFailed) &&
					   (lzc_polar( nullptr, &angle, &z ) == kFailed) && z.m_real.FlagIsOn( kfQNaN ) && z.m_imag.FlagIsOn( kfQNaN ) &&
					   (lzc_polar( &r, nullptr, &z ) == kFailed);

		CheckTrue( "a missing target or operand refused", b_ok, b_ok ? "kFailed" : "not refused" );
	}

	// In place, from its own parts, and in degrees...
	{
		CLZComplex	z;
		CLZComplex	t;
		CLZNumber	r;
		CLZNumber	angle;

		SetComplex( z, 10, "2.5", "37.25" );
		r.Assign( "2.5" );
		angle.Assign( "37.25" );
		geAngleUnit = AU_degrees;
		lzc_polar( &r, &angle, &t );
		lzc_polar( &z.m_real, &z.m_imag, &z );
		geAngleUnit = AU_radians;

		bool	b_ok = SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag );

		CheckTrue( "in place, from its own parts", b_ok, b_ok ? "as from copies" : "not as from copies" );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**The logarithm, the exponential and the square root, stage 4** --
	2026-10-22, on Rich's calls: lzc_log_e, lzc_exp_e and lzc_root_2, and
	the logarithm's argument and the exponential's angle in radians
	whatever the unit. Every expected part is generated -- exact Fractions
	where the part is rational, Python's decimal at 320 digits otherwise
	-- never typed. The near halves are each a case the guarded run alone
	gets wrong, so the retake or the root's exact decision is what settles
	them; the ties are square roots whose parts lie exactly on a half, which
	only the decision's tie to the even kept digit gets right.
*/

static void TestComplexLogExpRoot()
{
	Banner( "complex numbers, stage 4: the logarithm, the exponential and the square root" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	struct { const char* m_szFamily; char m_cFunction; long m_lUnit; long m_lPlaces; long m_lBase; const char* m_szReal;
			 const char* m_szImag; const char* m_szWantReal; const char* m_szWantImag; } kPins[] =
	{
		{ "random", 'l', 0, 36, 10, "-0.000000000000000002210527", "2423.7",
		  "7.793050577153772607490243753046959984",
		  "1.570796326794896619232233738097646396" },
		{ "random", 'l', 1, 36, 2, "0.000000000000000000000000000011", "-0.0000011",
		  "-11.110000001110000000000000000110111111",
		  "-1.100100100001111110110011010001000100" },
		{ "random", 'l', 0, 20, 3, "0.00000000002201", "-0.00000000000000000002",
		  "-102.00210220020011202210",
		  "-0.00000000020122202100" },
		{ "random", 'l', 2, 36, 7, "0.0000000000000000000000000000036253522", "-0.0000552053",
		  "-10.656551513626315532016654031442464321",
		  "-1.366532421435305204013153400354021660" },
		{ "random", 'l', 0, 36, 12, "-0.0000000000000000000004D153906E", "-0.0000000E960185966",
		  "-15.4E3822E000475560348998682D9724349981",
		  "-1.6D2404847ED6D38146948408E4D8ED791D77" },
		{ "random", 'l', 3, 50, 16, "-0.0000000000000000EECF253F1A52984", "0.0000000000000000D47266E13A",
		  "-2C.23B1D3665B343EDD361211C3B9B9701463FE5AD1114D6BCDC5",
		  "2.6A1FD0787F218CEF5D1FAF8D45EF641469363CDD0232A990B3" },
		{ "random", 'e', 0, 36, 10, "0.00000000000000000000000000000598", "68065727.6126",
		  "0.999565512149597807582207287732618708",
		  "0.029475191621298032159252706936223036" },
		{ "random", 'e', 1, 36, 2, "0.000000000000000000011000001", "0.00000001101001",
		  "0.111111111111111010111111100110000011",
		  "0.000000011010001111111111011010110011" },
		{ "random", 'e', 0, 20, 3, "0.000000000000000000002201211", "0.000000000000000000000000222",
		  "1.00000000000000000001",
		  "0.00000000000000000000" },
		{ "random", 'e', 2, 36, 7, "0.0000005", "-0.00000000006122666424411",
		  "1.000000500000153333360410601166624155",
		  "-0.000000000061230040225122120556564123" },
		{ "random", 'e', 0, 36, 12, "0.0000000E5", "-0.00000000000000000000000000003",
		  "1.0000000E5000005520600188012D04DED364",
		  "-0.000000000000000000000000000030000003" },
		{ "random", 'e', 3, 50, 16, "-0.00000A52C", "-EBA.9",
		  "0.CB97E366CBE45263CB6D99A9264CE1D9BD99CABFEA2BAB8448",
		  "-0.9B321479BD1F50096F5A19EFF9913DDF8491F15B122B52E56E" },
		{ "random", 's', 0, 36, 10, "0.0000000000000000032331", "909.93428328",
		  "21.329958781957362235413220125298764451",
		  "21.329958781957362235337432355780549944" },
		{ "random", 's', 1, 36, 2, "-0.0000000000000000000001", "0.0000000000000001",
		  "0.000000001011001110011100010101100010",
		  "0.000000001011011001110000011001000101" },
		{ "random", 's', 0, 20, 3, "-0.0000000000000000000000100002102022", "-0.0000000000000000000000111121",
		  "0.00000000000100220100",
		  "-0.00000000000200112101" },
		{ "random", 's', 2, 36, 7, "0.0000000000000022", "0.000000000000146155433633",
		  "0.000000232101366624011223113451355656",
		  "0.000000225603325420542130525411045452" },
		{ "random", 's', 0, 36, 12, "0.026EDE8D27756756", "-0.00000000000000000000E",
		  "0.173507857EE28130046362984ED12ED50DE9",
		  "-0.0000000000000000000350995967E54E0D63" },
		{ "random", 's', 3, 50, 16, "0.5668F0F4B", "0.00000000000000070AE64D",
		  "0.94BB299E263776835F16A8B721A52A7AD6E73E8075B92F0B0B",
		  "0.00000000000000060F998C45E630AA171F4B2F6487C38FAF4A" },
		{ "exact", 'l', 0, 36, 10, "0.6", "0.8",
		  "0",
		  "0.927295218001612232428512462922428804" },
		{ "exact", 'l', 1, 36, 10, "-0.8", "-0.6",
		  "0",
		  "-2.498091544796508851659834154562180246" },
		{ "exact", 'l', 0, 20, 10, "-1", "0",
		  "0",
		  "3.14159265358979323846" },
		{ "exact", 'l', 0, 20, 10, "1", "0",
		  "0",
		  "0" },
		{ "exact", 'l', 0, 20, 10, "0", "-1",
		  "0",
		  "-1.57079632679489661923" },
		{ "exact", 'e', 0, 20, 10, "0", "0",
		  "1",
		  "0" },
		{ "exact", 'e', 0, 36, 7, "0.3", "0",
		  "1.351345425252506211136035431211515253",
		  "0" },
		{ "exact", 'e', 2, 36, 10, "0", "0.5",
		  "0.877582561890372716116281582603829652",
		  "0.479425538604203000273287935215571388" },
		{ "exact", 's', 0, 36, 10, "3", "4",
		  "2",
		  "1" },
		{ "exact", 's', 0, 36, 10, "-4", "0",
		  "0",
		  "2" },
		{ "exact", 's', 0, 36, 10, "-3", "-4",
		  "1",
		  "-2" },
		{ "exact", 's', 0, 36, 10, "0", "2",
		  "1",
		  "1" },
		{ "exact", 's', 0, 36, 3, "-0.21", "-2.2",
		  "1",
		  "-1.1" },
		{ "exact", 's', 0, 36, 16, "-0.3C", "0.8",
		  "0.660D69C66BF2F4D44203F2EBAE21ED652FA5",
		  "0.A08B83F790CEB4D4A2DC0DEC9499332CB5B1" },
		{ "exact", 's', 0, 10, 10, "2", "0",
		  "1.4142135624",
		  "0" },
		{ "exact", 's', 0, 36, 10, "0", "1",
		  "0.707106781186547524400844362104849039",
		  "0.707106781186547524400844362104849039" },
		{ "tie", 's', 0, 0, 16, "551EE4.E85F", "-8084F.49",
		  "93C",
		  "-6F" },
		{ "tie", 's', 0, 0, 12, "-25724.E1", "54297.1",
		  "138",
		  "206" },
		{ "tie", 's', 0, 10, 14, "-B58B689A.DDDDDDDDDDDDDDD4B2ADA7", "0.0005690569",
		  "0.0000000306",
		  "C8C5" },
		{ "tie", 's', 0, 3, 8, "-67.6757466", "15.100131",
		  "0.676",
		  "7.415" },
		{ "tie", 's', 0, 0, 12, "510500.23", "7E5854.3",
		  "940",
		  "514" },
		{ "tie", 's', 0, 1, 12, "-1766.31E1", "-2152.EE7",
		  "26.0",
		  "-51.0" },
		{ "near", 'l', 1, 20, 10, "28.217407684949722020573344692292020869906023638971830903551894692411510211116762765298340216", "15.6581",
		  "3.47415975233438700797",
		  "0.50660447938135463026" },
		{ "near", 'l', 2, 20, 16, "0.00874A347ED60A3DFF3F01BA8261AAC8CFBC47F6A20040F3DFE758D5634E4A7605175E9F925FC9621EAD3A1659", "0.015C9",
		  "-5.2A97799DE8F77FFFFFFF",
		  "1.3357D013CF88D48E4095" },
		{ "near", 'l', 1, 36, 2, "10.1010100000101100100011110101100110110000100010010001111101000111000111000100100110110011111011100100101101", "1000.1010101",
		  "10.001101000100000101100100001110001011",
		  "1.010001011111001100001011000101010010" },
		{ "near", 'l', 0, 36, 3, "20.2220122002120020122112201110201022011212111200121200021210000210121000222100100112221002120202110000120211", "110.102",
		  "2.122200211101001221122201221020220112",
		  "1.001120221201202222000102200201020012" },
		{ "near", 'l', 0, 36, 5, "4.4040430224134423233121140332130432231012022423341032241012122434303323131301002323230140400013324232101104", "2.01",
		  "1.312043430340043444004322323400214024",
		  "0.144430103302003242044024423143424124" },
		{ "near", 'e', 0, 36, 12, "1.79D906230EED6253199086D90693475232E98139556E27E773928E5239356470893891DD96021226E5D39E0781E2895681E78DD704", "-0.6D268",
		  "4.4822D339031392521042178D0601D89D5121",
		  "-2.9D051851D093540448ED7669000000000000" },
		{ "near", 'e', 0, 36, 3, "2.1011001222111000000021121211112012121002200222102021202110112011020120012022120102111000102001022010120221", "-1.2221212211",
		  "-11.101011021010101012110221202201220210",
		  "-100.220220001000110012220010200211112221" },
		{ "near", 'e', 0, 36, 2, "10.11100111111101001011101101101000111010111111000111010101100000110101100101110001011000001001101101110011", "1.0110011111",
		  "11.000000110100010101000011101010011101",
		  "10010.000010001111011111101010001001100001" },
		{ "near", 'e', 2, 36, 5, "2.0030232443200014102221210022203144421403403403342122030111230000201401104120310421201141440132433103313112", "-1.2330230004",
		  "0.044242342211342031100101234022430114",
		  "-12.241241220013314103441303222210440334" },
		{ "near", 'e', 0, 36, 7, "2.2150115446066504312510252453426254124043524134051501061555126244543252263311056536362404315236326535040444", "-2.0030525",
		  "-4.215215354460516044656305234524262101",
		  "-12.135554555335044461430015423443345442" },
		{ "near", 's', 1, 20, 10, "-39.959719223175362290585700448736074292597986566132439999999999", "42.635951117855404549072842330014028138",
		  "3.03931648789421693379",
		  "7.01406900000000000000" },
		{ "near", 's', 2, 36, 2, "-100000.0111001101010110101101001110101111011100100010000010110000010110100010001110011100000100110011000000101101110011111111111111111111101101100010010111", "1101.01101110000010100001000100110001010011010000000000000000000000110001111100111",
		  "1.001001111011111011001101001011011001",
		  "101.110100000000000000000000000000000000" },
		{ "near", 's', 1, 20, 7, "40.451513501310462064552461365562652651105431331644554120256023410146346314423142652632520263136605444253434053361212600452", "13.2130153555303204140164130131056666666645603064343633661524303203310551",
		  "5.30326530010026416040",
		  "0.64316260160000000000" },
		{ "near", 's', 2, 36, 16, "-13.90AF8C2FB9BB4F84CF218FE5C000000000013EFAC20FE6B68000000000000000044FFE36B11038C6586A9AD6FFFFFFFFFFFFFEABB7B6CD2CFFFFFFFFFFFFFFFFFFFE3BAFE927BE034992E617", "B.73ADCEA856F36F11F9AFA4FFFFFFFFFFFFFB6793AADFE6FFFFFFFFFFFFFFFFFFF3C887E031F569CCE152CA",
		  "1.3EFAC20FE6B67FFFFFFFFFFFFFFFFFFFFFFF",
		  "4.986C55201900000000000000000000000000" },
		{ "near", 's', 1, 20, 5, "22.4041433231410113122040103321021232142424223342144121001021100013132004411112401431220430024331414213220421230013120231", "14.343404404002401321120234002000000000000221130234331242423123032113",
		  "3.40040040343211111322",
		  "1.12033440000000000000" },
		{ "scale", 'e', 0, 36, 10, "-1.75", "1000000000003",
		  "-0.121166990832906086102751828004876339",
		  "0.124563011182360846276009362078666126" },
		{ "scale", 'e', 1, 36, 16, "25", "-2.8",
		  "-215AFCCF8A48A3.C72908AA6327FF60C293A86C325DD8AEF935",
		  "-18EAD0682CAAF3.087B474A94E1686D340C3303E3E45C1991E1" },
		{ "scale", 'l', 0, 36, 10, "1000000000000000000000000000000", "-3",
		  "69.077552789821370520539743640530926228",
		  "-0.000000000000000000000000000003000000" },
		{ "scale", 'l', 0, 36, 2, "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011",
		  "-111101.110010110011111001000010011001001000",
		  "0.111110111001100001011110100101000001" },
		{ "scale", 's', 0, 1, 13, "0", "-0.000000000000000000000000000000061B3ABA6414557064",
		  "0.0",
		  "0.0" },
		{ "scale", 's', 0, 3, 4, "0", "0.00000000000000000000000000000000000000000000331001",
		  "0.000",
		  "0.000" },
		{ "scale", 's', 0, 36, 10, "10000000000000000000000000000000000000000", "-0.00000000000000000001",
		  "100000000000000000000.000000000000000000000000000000000000",
		  "0.000000000000000000000000000000000000" },
	};

	typedef long (*ComplexUnary)( CLZComplex*, CLZComplex* );

	auto	fnOf = []( char cFunction ) -> ComplexUnary
	{
		return (cFunction == 'l') ? lzc_log_e : ((cFunction == 'e') ? lzc_exp_e : lzc_root_2);
	};

	for ( const char* sz_family : { "random", "exact", "tie", "near", "scale" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZComplex	z;
			CLZComplex	t;

			glWorkingPrecision = pin.m_lPlaces;
			SetComplex( z, pin.m_lBase, pin.m_szReal, pin.m_szImag );
			geAngleUnit = static_cast<AngleUnit>( pin.m_lUnit );

			long	l_code = fnOf( pin.m_cFunction )( &z, &t );

			geAngleUnit = AU_radians;
			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &t.m_real, pin.m_szWantReal ) || !SpelledAs( &t.m_imag, pin.m_szWantImag ) ||
				 (t.m_real.m_lBase != pin.m_lBase) || (t.m_imag.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::string( 1, pin.m_cFunction ) + " base " + std::to_string( pin.m_lBase ) + " " + pin.m_szReal +
							  " " + pin.m_szImag;
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "log, exp, root, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), l_wrong == 0, szDetail );
	}

	// On the real axis each answers as the real function, to the character,
	//	an exact nought beside; on the imaginary axis the exponential is
	//	lz_sine_cosine's, in radians though degrees are set...
	{
		std::uint64_t	u_state = 2210;
		long			l_cases = 0;
		long			l_wrong = 0;

		for ( long i = 0; i < 120; i++ )
		{
			long		l_base = 2 + (i % 15);
			long		l_places = (i % 3 == 0) ? 20 : 36;
			CLZNumber	v;

			glWorkingPrecision = l_places;
			v.SetBase( l_base, false );
			v.Assign( RandomDigitsIn( u_state, l_base, 1 + static_cast<long>( (u_state >> 40) % 12 ), 3, (i % 2) == 1 ).c_str() );

			CLZComplex	z;
			CLZComplex	t;
			CLZNumber	nought;
			CLZNumber	want;
			CLZNumber	want_imag;
			CLZNumber	magnitude;

			nought.SetBase( l_base, false );
			nought.Assign( 0L );
			magnitude.Assign( &v );
			magnitude.m_lSign = kPositive;

			// ln |a| beside nought or pi...
			z.m_real.Assign( &magnitude );
			z.m_imag.Assign( &nought );
			long	l_got = lzc_log_e( &z, &t );
			long	l_want = lz_log_e( &magnitude, &want );
			l_wrong += (SameInEveryRespect( &t.m_real, l_got, &want, l_want ) && SpelledAs( &t.m_imag, "0" )) ? 0 : 1;

			// e^x beside nought...
			z.m_real.Assign( &v );
			l_got = lzc_exp_e( &z, &t );
			l_want = lz_exp_e( &v, &want );
			l_wrong += (SameInEveryRespect( &t.m_real, l_got, &want, l_want ) && SpelledAs( &t.m_imag, "0" )) ? 0 : 1;

			// sqrt a beside nought, or nought beside sqrt -a...
			l_got = lzc_root_2( &z, &t );
			l_want = lz_root_2( &magnitude, &want );
			l_wrong += (SameInEveryRespect( v.IsNegative() ? &t.m_imag : &t.m_real, l_got, &want, l_want ) &&
						SpelledAs( v.IsNegative() ? &t.m_real : &t.m_imag, "0" )) ? 0 : 1;

			// e^(iy), cos y + i sin y, in radians with degrees set...
			z.m_real.Assign( &nought );
			z.m_imag.Assign( &v );
			geAngleUnit = AU_degrees;
			l_got = lzc_exp_e( &z, &t );
			geAngleUnit = AU_radians;
			l_want = lz_sine_cosine( &v, &want_imag, &want );
			l_wrong += (SameInEveryRespect( &t.m_real, l_got, &want, l_want ) &&
						SameInEveryRespect( &t.m_imag, l_got, &want_imag, l_want )) ? 0 : 1;

			l_cases += 4;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld answers differ", l_wrong, l_cases );
		CheckTrue( "on an axis, the real function's", l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 36;

	// Radians whatever the unit: the same answers in degrees...
	{
		long	l_wrong = 0;

		for ( const char* sz_pair : { "-1|0", "1|1", "-2.5|-0.75", "0|3" } )
		{
			std::string	s_pair( sz_pair );
			size_t		l_bar = s_pair.find( '|' );

			for ( ComplexUnary fn : { lzc_log_e, lzc_exp_e } )
			{
				CLZComplex	z;
				CLZComplex	radians;
				CLZComplex	degrees;

				SetComplex( z, 10, s_pair.substr( 0, l_bar ).c_str(), s_pair.substr( l_bar + 1 ).c_str() );
				fn( &z, &radians );
				geAngleUnit = AU_degrees;
				fn( &z, &degrees );
				geAngleUnit = AU_radians;
				l_wrong += (SameAnswer( &radians.m_real, &degrees.m_real ) && SameAnswer( &radians.m_imag, &degrees.m_imag )) ? 0 : 1;
			}
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 8 differ from radians", l_wrong );
		CheckTrue( "radians whatever the unit", l_wrong == 0, szDetail );
	}

	// Nought, a not-a-number, the point at infinity, too large and too small,
	//	and two bases in one value...
	{
		CLZComplex	z;
		CLZComplex	t;
		CLZNumber	nan;
		CLZNumber	infinite;
		bool		b_ok = true;

		MakeNotANumber( nan );
		MakeInfinite( infinite );

		SetComplex( z, 10, "0", "0" );
		b_ok = b_ok && (lzc_log_e( &z, &t ) == kUndefinedResult) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );

		for ( ComplexUnary fn : { lzc_log_e, lzc_exp_e, lzc_root_2 } )
		{
			SetComplex( z, 10, "1", "2" );
			z.m_imag.Assign( &nan );
			b_ok = b_ok && (fn( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );
		}

		SetComplex( z, 10, "1", "2" );
		z.m_real.Assign( &infinite );
		b_ok = b_ok && (lzc_log_e( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );
		b_ok = b_ok && (lzc_exp_e( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );
		b_ok = b_ok && (lzc_root_2( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity ) &&
			   t.m_real.FlagIsOff( kfQNaN );

		SetComplex( z, 10, "100000000", "1" );
		b_ok = b_ok && (lzc_exp_e( &z, &t ) == kTooLarge) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
		SetComplex( z, 10, "100000000", "0" );
		b_ok = b_ok && (lzc_exp_e( &z, &t ) == kTooLarge) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );

		glWorkingPrecision = 10;
		SetComplex( z, 10, "-100000000", "1" );
		b_ok = b_ok && (lzc_exp_e( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "0.0000000000" ) && SpelledAs( &t.m_imag, "0.0000000000" );
		glWorkingPrecision = 36;

		SetComplex( z, 10, "1", "2" );
		z.m_imag.SetBase( 16, false );
		z.m_imag.Assign( "2" );
		b_ok = b_ok && (lzc_root_2( &z, &t ) == kBaseMismatch) && (lzc_log_e( &z, &t ) == kBaseMismatch) &&
			   (lzc_exp_e( &z, &t ) == kBaseMismatch) && t.m_real.FlagIsOn( kfQNaN );

		CheckTrue( "nought, flags, limits, two bases", b_ok, b_ok ? "refused or flagged as the layer and the real functions do" : "an answer wrong" );
	}

	// In place, as from a copy...
	{
		long	l_wrong = 0;

		for ( ComplexUnary fn : { lzc_log_e, lzc_exp_e, lzc_root_2 } )
		{
			CLZComplex	z;
			CLZComplex	copy;
			CLZComplex	t;

			SetComplex( z, 7, "-3.25", "1.5" );
			SetComplex( copy, 7, "-3.25", "1.5" );
			fn( &copy, &t );
			fn( &z, nullptr );
			l_wrong += (SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 3 differ", l_wrong );
		CheckTrue( "log, exp, root in place", l_wrong == 0, szDetail );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**The twelve circular and hyperbolic functions, stage 5a** -- 2026-10-22,
	on Rich's calls: lzc_sine to lzc_hypercosecant, the circular six reading
	the whole complex angle in the host's unit, the hyperbolic six in
	radians, and on the real axis each the real function's answer to the
	character. Every expected part is generated -- Python's decimal at 700
	digits, never typed -- and each family past the random and the exact is
	made of cases a mutant gets wrong: the unit's half by the build without
	lzi_complex_trig_decide, a small z's ties by the build without
	lzi_complex_trig_near_decide, the near halves by the build without a
	retake, the poles by the build that never widens there, the large by the
	build that never widens for cosh, the far quotients by the build without
	the far path, and a small z steered near a half by the build whose
	polynomial decider ignores its remainder's bound.
*/

static void TestComplexTrig()
{
	Banner( "complex numbers, stage 5a: the circular and hyperbolic functions" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	struct { const char* m_szFamily; long m_lFunction; long m_lUnit; long m_lPlaces; long m_lBase; const char* m_szReal;
			 const char* m_szImag; const char* m_szWantReal; const char* m_szWantImag; } kPins[] =
	{
		{ "random", 0, 0, 36, 10, "-1", "-0.00000000077118728",
		  "-0.841470984807896506902726265646119092",
		  "-0.000000000416674265640178707363698407" },
		{ "random", 0, 1, 20, 7, "0.000000001244124", "-0.00000000000000005005133443",
		  "0.00000000001115452634",
		  "-0.00000000000000000042" },
		{ "random", 0, 2, 36, 16, "-0.00000000000A59254306AE", "-0.0000000000000000A3D5FA",
		  "-0.000000000000299CD8AA5F915D046E135624",
		  "-0.00000000000000000292D2B1E997F4E38CBD" },
		{ "random", 1, 0, 36, 10, "-0.0000000000000000000000000576401", "0.00000000000096339028",
		  "1.000000000000000000000000464060415799",
		  "0.000000000000000000000000000000000000" },
		{ "random", 1, 1, 20, 7, "-0.0000000000000000000000000053122336", "0.0000000000000000000000000000032422022",
		  "1.00000000000000000000",
		  "0.00000000000000000000" },
		{ "random", 1, 2, 36, 16, "-0.0000000000000000000000AA", "0.00000000000000000000000B1F9D80CF88E",
		  "1.000000000000000000000000000000000000",
		  "0.000000000000000000000000000000000000" },
		{ "random", 2, 0, 36, 10, "-0.0000000000000088503", "-0.000000000000000076259352",
		  "-0.000000000000008850300000000000000000",
		  "-0.000000000000000076259352000000000000" },
		{ "random", 2, 1, 20, 7, "-0.0000000000000000000000000044526445534", "56.13",
		  "0.00000000000000000000",
		  "0.42126336301632650646" },
		{ "random", 2, 2, 36, 16, "0.000000000000000000000000000048C3B6C", "-0.000D6B",
		  "0.000000000000000000000000000001249A7B",
		  "-0.000035F4F45F4E5DDF6E47AC163708B760FD" },
		{ "random", 3, 0, 36, 10, "-0.0000000000000087556906758", "-1.6622094606",
		  "-0.000000000000001356480699896616979815",
		  "1.074674652081900875047484173026662838" },
		{ "random", 3, 1, 20, 7, "-0.000000000000000000000064621366", "-0.00000000000000000000001",
		  "-1123606433563440413501010.02160623411621161323",
		  "115113045324214000164312.34421361060002556464" },
		{ "random", 3, 2, 36, 16, "-1EC", "0.000000000000000FC67BCF028",
		  "-0.20571F44F335F251C6B0D7B70489298AC751",
		  "-0.00000000000000004072E0A55637CBCCDF8D" },
		{ "random", 4, 0, 36, 10, "-0.00000000000000000000007606614", "0.311",
		  "0.953514682675620283135123041021476572",
		  "-0.000000000000000000000021856720728350" },
		{ "random", 4, 1, 20, 7, "0.3", "-0.63",
		  "0.66665211502305434046",
		  "-0.00002005021145131304" },
		{ "random", 4, 2, 36, 16, "-0.00000000000000000000000A5CAE3795", "-0.0000000000C851C75AD9",
		  "0.FFFFFFFFFFFFFFFFFFFFFFFB0CA6DC0C4D28",
		  "0.000000000000000000000000000000000001" },
		{ "random", 5, 0, 36, 10, "-0.0000000000000000000000000000002692", "-0.000000000038943003507",
		  "-0.000000000177507279525428971858122858",
		  "25678553525.545355671421947480579457580513905577" },
		{ "random", 5, 1, 20, 7, "0.0000000000000000000000000000036650464", "0.0000000000000000000000000000021231",
		  "14010403243354616210333052532636.45552140030310034331",
		  "-6023264500256522654220166030043.04000525420101346316" },
		{ "random", 5, 2, 36, 16, "0.0000000000009594F9DDEEDA", "-0.0000000000000000000000DC5B551",
		  "6CF4146CF649B3.086C6C7FF99A2BFEA267AB73A62837C30A24",
		  "A081.3E180D64F52AF71D8DFAB3BC43B99FDF68EC" },
		{ "random", 6, 0, 36, 10, "0.0000000000000000000000009", "0.1",
		  "0.000000000000000000000000895503748750",
		  "0.099833416646828152306814198410622027" },
		{ "random", 6, 1, 20, 7, "0.1", "0.00000000000000000004364261",
		  "0.10011114012023446512",
		  "0.00000000000000000005" },
		{ "random", 6, 2, 36, 16, "-0.00004C92AF", "0.000000000000000000000000DCD7",
		  "-0.00004C92AF0001244DF2B44F1D4966347EFE",
		  "0.000000000000000000000000DCD7000009E1" },
		{ "random", 7, 0, 36, 10, "-0.5", "17.7",
		  "0.461078069930543737282638039300236039",
		  "0.475542430462169486777203185933476751" },
		{ "random", 7, 1, 20, 7, "0.00000006015", "-14.253",
		  "0.25060206211165301622",
		  "0.00000005361563611146" },
		{ "random", 7, 2, 36, 16, "0.1", "-0.00000000000000008D4B",
		  "1.00800AAB05B1FB245198050937BB036806CB",
		  "-0.000000000000000008D628DAD73FA28259DC" },
		{ "random", 8, 0, 36, 10, "0.000000000000000039992896", "0.0000000000009384864943",
		  "0.000000000000000039992896000000000000",
		  "0.000000000000938486494300000000000000" },
		{ "random", 8, 1, 20, 7, "0.000000135032511353", "12.003055015",
		  "0.00000015546054031135",
		  "-0.30425565441522536462" },
		{ "random", 8, 2, 36, 16, "-0.08F756F0D5FAE", "-1B.2966",
		  "-0.2D94DD7D27E1C08696C1CA97AD4B6A4C57FE",
		  "2.03D8BF00E84BBE4137D31B07C925764A6E54" },
		{ "random", 9, 0, 36, 10, "0.00000000000000098316278828", "-0.000000000000000000000000009",
		  "1017125558372134.853069437411461873522228935443901122",
		  "9310.899613444443938678131093302817332836" },
		{ "random", 9, 1, 20, 7, "-34.552", "-0.0000000000030614633034",
		  "-1.00000000000000000000",
		  "0.00000000000000000000" },
		{ "random", 9, 2, 36, 16, "0.00000000011E", "-0.DBE",
		  "0.0000000001F2EEA046F727539BA3A5E1D754",
		  "0.DCE408E7AEA21C398213DEC1D44E002CD964" },
		{ "random", 10, 0, 36, 10, "-0.000000000000000051017333", "0.0000000000000006983122",
		  "1.000000000000000000000000000000242519",
		  "0.000000000000000000000000000000035626" },
		{ "random", 10, 1, 20, 7, "0.000000000000505550601", "5.06655",
		  "2.25362652552245256646",
		  "0.00000000000355322361" },
		{ "random", 10, 2, 36, 16, "-0.000FA2", "-0.00000000000000000047C20ED",
		  "0.FFFFFF85CEBE309A7CE9BDD031D4F5B86FCC",
		  "-0.000000000000000000000461C7A8133901B8" },
		{ "random", 11, 0, 36, 10, "0.000000000000049462161", "0.000000000655616586",
		  "115072.943001592263548125079279536520378024",
		  "-1525281720.539393990534368596536448440933340895" },
		{ "random", 11, 1, 20, 7, "0.0000000043132", "-36",
		  "-0.00000000126531035401",
		  "1.02143430161666261203" },
		{ "random", 11, 2, 36, 16, "5.68", "-0.00000000000000000000000047186",
		  "0.024C5219B8749B0DEC90EE6580C2E8507EE4",
		  "0.00000000000000000000000000A36478E4F1" },
		{ "unit", 0, 3, 1, 9, "-512.6", "0.000000000000446",
		  "-0.5",
		  "0.0" },
		{ "unit", 1, 2, 1, 9, "325.6", "0.000000000782",
		  "-0.5",
		  "0.0" },
		{ "unit", 1, 1, 20, 5, "2200", "-0.00000000000000000000000000000323",
		  "0.22222222222222222223",
		  "0.00000000000000000000" },
		{ "unit", 0, 1, 0, 2, "11110", "0.00000000000000001",
		  "1",
		  "0" },
		{ "unit", 1, 1, 1, 13, "93", "0.000000000000000000000000000000000000B11",
		  "-0.7",
		  "0.0" },
		{ "unit", 0, 1, 10, 13, "-B7", "0.0000000000000000000000000000000000000000000000182",
		  "-0.6666666667",
		  "0.0000000000" },
		{ "exact", 1, 1, 20, 10, "90", "1",
		  "0",
		  "-0.01745417862959511118" },
		{ "exact", 0, 1, 20, 10, "0", "1",
		  "0",
		  "0.01745417862959511118" },
		{ "exact", 2, 1, 20, 10, "90", "1",
		  "0",
		  "57.30159715911290867049" },
		{ "exact", 3, 2, 20, 10, "100", "0.1",
		  "0",
		  "-0.00157079503486797669" },
		{ "exact", 4, 3, 20, 10, "250", "-2",
		  "0",
		  "-79.57537718943022097762" },
		{ "exact", 5, 1, 20, 10, "-180", "1.5",
		  "0",
		  "38.19282336780004318252" },
		{ "exact", 0, 1, 20, 10, "180", "90",
		  "0",
		  "-2.30129890230729487346" },
		{ "exact", 1, 1, 20, 10, "0", "90",
		  "2.50917847865805678201",
		  "0" },
		{ "near nought", 3, 0, 36, 14, "0.00000000000000000000000000000000000000000000000000000000000000000000000000000000010CABBB10965C10C148BB6851202B7215A0D46D4DB2", "0.00000000000000000000000000000000000000000000000000000000000000000000000000000000010CABBB10965C10C148BB6851202B7215A0D46D4DB2",
		  "680097B1CB11434642D3C7A74A01000921C6ADC1911C6D5582AD1C2BA188CD1150237813A945685B4C.A6B97A94017382B2B3612C56138C1AC5061A",
		  "-680097B1CB11434642D3C7A74A01000921C6ADC1911C6D5582AD1C2BA188CD1150237813A945685B4C.A6B97A94017382B2B3612C56138C1AC5061B" },
		{ "near nought", 11, 0, 20, 6, "0.000000000000000000000000000000000000000000000000000000000023053353530155550541354043543542243325553444410304", "0.000000000000000000000000000000000000000000000000000000000023053353530155550541354043543542243325553444410304",
		  "11043113303243542304355111411411014250202225004550225042121.34205032245410453104",
		  "-11043113303243542304355111411411014250202225004550225042121.34205032245410453105" },
		{ "near nought", 11, 0, 1, 14, "-0.00000000000000000000000000000000898D0B45D408", "0.00000000000000000000000000000000898D0B45D408",
		  "-B3CA81541A1C47048872CC9B147D030C.3",
		  "-B3CA81541A1C47048872CC9B147D030C.4" },
		{ "near nought", 5, 0, 0, 10, "-0.00000000000001048576", "-0.00000000000001048576",
		  "-47683715820313",
		  "47683715820312" },
		{ "near nought", 8, 0, 20, 4, "-0.000000000000000000000000000000000000000000000000000000000000000000002", "0.000000000000000000122",
		  "0.00000000000000000000",
		  "0.00000000000000000013" },
		{ "near nought", 3, 0, 1, 10, "-0.00000000000000000000000000000000281474976710656", "-0.00000000000000000000000000000000281474976710656",
		  "-177635683940025046467781066894531.2",
		  "177635683940025046467781066894531.3" },
		{ "near nought", 0, 0, 36, 8, "0.000000000000000000000000000000000000000000000000000000000000000000000032", "0.0000000000000000000000000000000000644",
		  "0.000000000000000000000000000000000000",
		  "0.000000000000000000000000000000000065" },
		{ "near nought", 3, 0, 3, 12, "-0.0000000000000000052356E91347790E107931ED994", "0.0000000000000000052356E91347790E107931ED994",
		  "-11D54D7096467EE306.E64",
		  "-11D54D7096467EE306.E65" },
		{ "near nought", 5, 0, 0, 6, "-0.0000000000000000042345300424", "-0.0000000000000000042345300424",
		  "-40154420544110105",
		  "40154420544110104" },
		{ "near nought", 0, 0, 20, 4, "0.000000000000000000000000000000000000000000000031", "-0.000000000000000000202",
		  "0.00000000000000000000",
		  "-0.00000000000000000021" },
		{ "near half", 4, 0, 20, 7, "-1.42", "1.301536220350214125245611642551165403456345156320141141616431300503502201040120414402334435",
		  "-0.01102216356032132030",
		  "-0.33511325306644160551" },
		{ "near half", 2, 1, 36, 10, "-66.8161321536", "87.3881280268173221692507013358408587039424943508000726894820312552314290289058667532198670068274887709928621",
		  "-0.073140749236734480174312363249106989",
		  "1.064945809904330194232394701810673165" },
		{ "near half", 9, 0, 36, 7, "-1.43561", "-0.5355561644352640424335110231243352026325564510663344235131216115632354100541342422324201310335231550521633",
		  "-0.665465013250044021521023555222631623",
		  "0.034151266112356002225410334231542113" },
		{ "near half", 0, 1, 36, 2, "101000.011", "-10110.011110010010011001111001101111101011011111101001001010111011000011101001010000000011100000000001101010001",
		  "0.101100101100000100100001100110110111",
		  "-0.010011100111100011011011000101000111" },
		{ "near half", 6, 0, 20, 3, "0.1001121010222212122022020000022020102110001202111201211202210012212222112110212210220212", "0.2001",
		  "0.02112011022210211211",
		  "0.20000012222000212021" },
		{ "near half", 11, 0, 36, 2, "0.10001101010010101100100000000111100010101000000110111111010101011110001101001101000100011110010111001001", "1.11011111011",
		  "-0.001000110101111000111001010001101010",
		  "-0.111000100101110100111100111101101001" },
		{ "near half", 11, 0, 20, 16, "0.3F104A81344700000007BF16E9C3598EC1871B62B0CF40AB31D425E04CD47862047657663CA9E51C33A6BDC8C8", "1.9D1F5EA8",
		  "-0.0294BDFC38EC3C347C55",
		  "-0.F8A0706F23638EE16AB8" },
		{ "near half", 0, 1, 20, 12, "-0.02613E574826236703270686842238315510E4E46508385D3356241551250298D842D7604629127D0D4E59D356", "1.824017",
		  "-0.0006386190251E213458",
		  "0.04291819DE8DE8D791D5" },
		{ "pole", 9, 0, 36, 7, "-0.00000000000000344", "3.0663651432036134110263402244652226644",
		  "-162615410055145.253023150212211415541440556103240065",
		  "-0.000000000613623201042444412353260556" },
		{ "pole", 4, 1, 20, 16, "5A", "-0.000000000000000000000000009CE",
		  "0",
		  "-5D7FD92BE78C080FED7980EAAF5A.8142C95622D4579C0D68" },
		{ "pole", 11, 0, 10, 7, "-0.00000000241", "3.066365143203613411026340224465223",
		  "246224063.6252016332",
		  "0.0000000000" },
		{ "pole", 11, 0, 10, 10, "0.0000000000000000000000000000000111", "3.1415926535897932385",
		  "-7954042.4724169258",
		  "26769019461318409709.5476647918" },
		{ "pole", 10, 0, 10, 2, "0.0000000000000000001", "1.1001001",
		  "100000010010.1101001011",
		  "-1000.0010010111" },
		{ "pole", 2, 1, 20, 2, "1011010", "-0.00000000000000000111",
		  "0",
		  "-100000101111011000110111.01010100000100010011" },
		{ "large", 0, 0, 20, 16, "-2.BF", "42.B93",
		  "-3B1D4BA210FB2686603E7916.C44A2C16B739276A8B64",
		  "-8D979FF1AF4787B8588680E2.B66A5B458413AE2173D5" },
		{ "large", 0, 0, 10, 10, "-1.07", "62.686",
		  "-734941496916496746250248988.1496062469",
		  "402260620996652088802279987.3070714107" },
		{ "large", 1, 0, 20, 16, "0.78", "67.77",
		  "1130CBCF5DC40804CDD2F94B294773EB962CBD.6C72B66D280ED06DAA1C",
		  "-8B489D001FA0951B0FE2834590A94EBF168C2.99BBA89E41BE8FD04CF6" },
		{ "large", 7, 0, 10, 7, "154.553", "1.05",
		  "536623066251555511055513330326624234222214504.3651512415",
		  "1366625520514106115632052352155052411434655163.5066662104" },
		{ "far", 2, 0, 20, 10, "1", "10000000000000000000000000000000000000000",
		  "0.00000000000000000000",
		  "1.00000000000000000000" },
		{ "far", 3, 1, 20, 10, "30", "-100000000000000000000000000000000000000000",
		  "0.00000000000000000000",
		  "1.00000000000000000000" },
		{ "far", 8, 0, 20, 16, "C9F2C9CD04674EDEA40000000", "0.C",
		  "1.00000000000000000000",
		  "0.00000000000000000000" },
		{ "far", 10, 0, 20, 7, "6106454640561632563653142", "1",
		  "0.00000000000000000000",
		  "0.00000000000000000000" },
		{ "far", 4, 2, 20, 10, "50", "1000000000",
		  "0.00000000000000000000",
		  "0.00000000000000000000" },
		{ "far", 11, 0, 20, 10, "100000000", "-2",
		  "0.00000000000000000000",
		  "0.00000000000000000000" },
		{ "small half", 6, 0, 20, 6, "0.022402", "0.035345333250053044135054045200154430552531042221235445331320053213530350134435553225140414",
		  "0.02233324403220045201",
		  "0.03535110305352430052" },
		{ "small half", 5, 0, 10, 12, "-0.0726036697E7372E76167791DD93043314D153E3E894785015002EDE572D96658E07302E9377951", "0.06477",
		  "-E.24DE155803",
		  "-9.DD8E22EED5" },
		{ "small half", 0, 0, 36, 16, "-0.06", "-0.0533075EB400F75555555555555555555505345DBDFE58121F24B08CA37D119FD99A6B96726B17EB6728EF55474F5EA9C6FCC5B1DD",
		  "-0.06002D1871CF4C549CF3FCE08109A2373A53",
		  "-0.0532C134B7ECCFB12283C0F37AF6E4EBA43C" },
		{ "small half", 8, 0, 20, 12, "0.06461", "-0.04E19177D9957787800189954E6338878D552E56833336563805745409E069350739791E3676E3889E7E9020D3",
		  "0.06466835502663067762",
		  "-0.04E0792066357001DD2D" },
		{ "small half", 3, 0, 36, 6, "-0.0204323255453405045345504322131141000511351101202341424115224112115241100103103402003541214545253252223001", "0.032424223",
		  "-4.343343423533130002445332101413012135",
		  "-11.330321012454533521440345300240350134" },
		{ "small half", 10, 0, 10, 12, "-0.051D0D34952829DD97300D12E2296067E67D718290988520592EE5304243934762391064161D41E2", "0.01926D1899",
		  "0.EEE035353D",
		  "0.0009127780" },
	};

	typedef long (*ComplexUnary)( CLZComplex*, CLZComplex* );
	typedef long (*RealUnary)( CLZNumber*, CLZNumber* );

	static const ComplexUnary	kComplex[ 12 ] =
	{
		lzc_sine, lzc_cosine, lzc_tangent, lzc_cotangent, lzc_secant, lzc_cosecant,
		lzc_hypersine, lzc_hypercosine, lzc_hypertangent, lzc_hypercotangent, lzc_hypersecant, lzc_hypercosecant,
	};
	static const RealUnary		kReal[ 12 ] =
	{
		lz_sine, lz_cosine, lz_tangent, lz_cotangent, lz_secant, lz_cosecant,
		lz_hypersine, lz_hypercosine, lz_hypertangent, lz_hypercotangent, lz_hypersecant, lz_hypercosecant,
	};

	// On the imaginary axis, the part it lands in, the real function, and its
	//	sign: sin iy = i sinh y, cot iy = -i coth y, cosh iy = cos y...
	static const long			kImaginaryAxis[ 12 ][ 3 ] =
	{
		{ 1, 6, 1 }, { 0, 7, 1 }, { 1, 8, 1 }, { 1, 9, -1 }, { 0, 10, 1 }, { 1, 11, -1 },
		{ 1, 0, 1 }, { 0, 1, 1 }, { 1, 2, 1 }, { 1, 3, -1 }, { 0, 4, 1 }, { 1, 5, -1 },
	};

	for ( const char* sz_family : { "random", "unit", "exact", "near nought", "near half", "pole", "large", "far", "small half" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZComplex	z;
			CLZComplex	t;

			glWorkingPrecision = pin.m_lPlaces;
			SetComplex( z, pin.m_lBase, pin.m_szReal, pin.m_szImag );
			geAngleUnit = static_cast<AngleUnit>( pin.m_lUnit );

			long	l_code = kComplex[ pin.m_lFunction ]( &z, &t );

			geAngleUnit = AU_radians;
			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &t.m_real, pin.m_szWantReal ) || !SpelledAs( &t.m_imag, pin.m_szWantImag ) ||
				 (t.m_real.m_lBase != pin.m_lBase) || (t.m_imag.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::to_string( pin.m_lFunction ) + " base " + std::to_string( pin.m_lBase ) + " " + pin.m_szReal +
							  " " + pin.m_szImag;
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "circular and hyperbolic, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), l_wrong == 0, szDetail );
	}

	// On the real axis each is the real function's to the character, in the
	//	host's unit for the circular six, an exact nought beside; on the
	//	imaginary axis in radians one of the other family's. A real function
	//	refused as infinite -- sinh past the limit, a pole -- is the point at
	//	infinity, with its code...
	auto	SameOnTheAxis = []( CLZComplex* pGot, CLZNumber* pOn, CLZNumber* pOff, long lGot, CLZNumber* pWant, long lWant ) -> bool
	{
		if ( pWant->FlagIsOn( kfInfinity ) )
			return (lGot == lWant) && pGot->m_real.FlagIsOn( kfInfinity ) && pGot->m_imag.FlagIsOn( kfInfinity );

		return SameInEveryRespect( pOn, lGot, pWant, lWant ) && SpelledAs( pOff, "0" );
	};

	{
		std::uint64_t	u_state = 5155;
		long			l_cases = 0;
		long			l_wrong = 0;

		for ( long i = 0; i < 48; i++ )
		{
			long		l_base = 2 + (i % 15);
			long		l_unit = (i % 4 == 3) ? 1 : 0;
			CLZNumber	v;

			glWorkingPrecision = (i % 3 == 0) ? 20 : 36;
			v.SetBase( l_base, false );
			v.Assign( RandomDigitsIn( u_state, l_base, 1 + static_cast<long>( (u_state >> 40) % 10 ), 4, (i % 2) == 1 ).c_str() );

			CLZNumber	nought;

			nought.SetBase( l_base, false );
			nought.Assign( 0L );

			for ( long l_fn = 0; l_fn < 12; l_fn++ )
			{
				CLZComplex	z;
				CLZComplex	t;
				CLZNumber	want;

				// ...the real axis...
				z.m_real.Assign( &v );
				z.m_imag.Assign( &nought );
				geAngleUnit = static_cast<AngleUnit>( l_unit );

				long	l_got = kComplex[ l_fn ]( &z, &t );
				long	l_want = kReal[ l_fn ]( &v, &want );

				geAngleUnit = AU_radians;
				l_wrong += SameOnTheAxis( &t, &t.m_real, &t.m_imag, l_got, &want, l_want ) ? 0 : 1;

				// ...and the imaginary axis, in radians...
				z.m_real.Assign( &nought );
				z.m_imag.Assign( &v );
				l_got = kComplex[ l_fn ]( &z, &t );
				l_want = kReal[ kImaginaryAxis[ l_fn ][ 1 ] ]( &v, &want );

				if ( kImaginaryAxis[ l_fn ][ 2 ] < 0 )
					lz_negate( &want );

				CLZNumber*	p_on = (kImaginaryAxis[ l_fn ][ 0 ] == 0) ? &t.m_real : &t.m_imag;
				CLZNumber*	p_off = (kImaginaryAxis[ l_fn ][ 0 ] == 0) ? &t.m_imag : &t.m_real;

				l_wrong += SameOnTheAxis( &t, p_on, p_off, l_got, &want, l_want ) ? 0 : 1;
				l_cases += 2;
			}
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld answers differ", l_wrong, l_cases );
		CheckTrue( "on an axis, the real function's", l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 20;

	// A pole on the real axis is the point at infinity, with the real
	//	function's code...
	{
		long	l_wrong = 0;
		struct { long m_lFunction; long m_lUnit; const char* m_szX; } kPoles[] =
		{
			{ 2, 1, "90" }, { 3, 1, "0" }, { 4, 1, "270" }, { 5, 1, "-180" }, { 2, 2, "100" }, { 3, 3, "500" },
			{ 3, 0, "0" }, { 5, 0, "0" }, { 9, 0, "0" }, { 11, 0, "0" },
		};

		for ( auto& pole : kPoles )
		{
			CLZComplex	z;
			CLZComplex	t;
			CLZNumber	want;

			SetComplex( z, 10, pole.m_szX, "0" );
			geAngleUnit = static_cast<AngleUnit>( pole.m_lUnit );

			long	l_got = kComplex[ pole.m_lFunction ]( &z, &t );
			long	l_want = kReal[ pole.m_lFunction ]( &z.m_real, &want );

			geAngleUnit = AU_radians;
			l_wrong += (want.FlagIsOn( kfInfinity ) && (l_got == l_want) && t.m_real.FlagIsOn( kfInfinity ) &&
						t.m_imag.FlagIsOn( kfInfinity )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 10 not the point at infinity", l_wrong );
		CheckTrue( "a pole on the real axis", l_wrong == 0, szDetail );
	}

	// The hyperbolic six in radians whatever the unit -- the circular six
	//	reading the whole angle in it are the unit and exact pins...
	{
		long	l_wrong = 0;

		for ( long l_fn = 6; l_fn < 12; l_fn++ )
		{
			CLZComplex	z;
			CLZComplex	radians;
			CLZComplex	degrees;

			SetComplex( z, 10, "0.75", "-1.25" );
			kComplex[ l_fn ]( &z, &radians );
			geAngleUnit = AU_degrees;
			kComplex[ l_fn ]( &z, &degrees );
			geAngleUnit = AU_radians;
			l_wrong += (SameAnswer( &radians.m_real, &degrees.m_real ) && SameAnswer( &radians.m_imag, &degrees.m_imag )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 6 differ", l_wrong );
		CheckTrue( "the hyperbolic six in radians", l_wrong == 0, szDetail );
	}

	// A not-a-number, the point at infinity, two bases, too large...
	{
		CLZComplex	z;
		CLZComplex	t;
		CLZNumber	nan;
		CLZNumber	infinite;
		bool		b_ok = true;

		MakeNotANumber( nan );
		MakeInfinite( infinite );

		for ( ComplexUnary fn : kComplex )
		{
			SetComplex( z, 10, "1", "2" );
			z.m_imag.Assign( &nan );
			b_ok = b_ok && (fn( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );

			SetComplex( z, 10, "1", "2" );
			z.m_real.Assign( &infinite );
			b_ok = b_ok && (fn( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );

			SetComplex( z, 10, "1", "2" );
			z.m_imag.SetBase( 16, false );
			z.m_imag.Assign( "2" );
			b_ok = b_ok && (fn( &z, &t ) == kBaseMismatch) && t.m_real.FlagIsOn( kfQNaN );
		}

		for ( ComplexUnary fn : { lzc_sine, lzc_cosine } )
		{
			SetComplex( z, 10, "1", "10000000" );
			b_ok = b_ok && (fn( &z, &t ) == kTooLarge) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
		}

		for ( ComplexUnary fn : { lzc_hypersine, lzc_hypercosine } )
		{
			SetComplex( z, 10, "-10000000", "1" );
			b_ok = b_ok && (fn( &z, &t ) == kTooLarge) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
		}

		CheckTrue( "flags, limits, two bases", b_ok, b_ok ? "refused or flagged as the layer does" : "an answer wrong" );
	}

	// In place, as from a copy...
	{
		long	l_wrong = 0;

		for ( ComplexUnary fn : kComplex )
		{
			CLZComplex	z;
			CLZComplex	copy;
			CLZComplex	t;

			SetComplex( z, 7, "-3.25", "1.5" );
			SetComplex( copy, 7, "-3.25", "1.5" );
			fn( &copy, &t );
			fn( &z, nullptr );
			l_wrong += (SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 12 differ", l_wrong );
		CheckTrue( "circular and hyperbolic in place", l_wrong == 0, szDetail );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**The power, the roots, the base's exponential and the other
	logarithms, stage 5b** -- 2026-10-22, on Rich's calls: lzc_power,
	lzc_raise_n, lzc_root_n, lzc_root_3, lzc_exp_b, lzc_log_n and lzc_log_b,
	exponents, orders and log bases complex throughout, and an odd root of
	a negative real the real root. Every expected part is generated --
	exact Fractions where the part is rational, Python's decimal deep enough
	for the case's own places otherwise -- never typed. Each family past the
	random and the exact is made of cases a mutant gets wrong: the root ties
	by the build without the proof in doubled coordinates, the near roots
	by the build without the root decider, the unit's halves by the build
	without the axis side decider, the unit roots by the build without the
	unit-root decider, the unit bound by the build that ignores that
	decider's bound -- roots of -1 of order 6 + di, d near 4 and 5, found
	by steering, since no drawing reached them -- the axis by the build
	without the degree phase, the odd roots by the build that takes their
	principal value, the axis products by the build without the
	two-irrationals proof, the logs by the build that never proves one
	exact, the steered by the build without a retake, the log host by the
	build whose real logarithm is not the host's call, the noughts by the
	build that never finds a part exactly nought below the width, the large
	and the large exponent by the builds that do not widen for the answer's
	magnitude and for the exponent's, and the log noughts by the build that
	never finds a logarithm's part nought by structure -- log_i 9's real
	part, say, which this test found spelled to the width.
*/

static void TestComplexPowers()
{
	Banner( "complex numbers, stage 5b: the power, the roots and the other logarithms" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	// w the exponent, the order or the base: lzc_log_n takes the base first,
	//	as lz_log_n does...
	auto	Call = []( char cFunction, CLZComplex* pZ, CLZComplex* pW, CLZComplex* pTarget ) -> long
	{
		switch ( cFunction )
		{
		case 'p':	return lzc_power( pZ, pW, pTarget );
		case 'n':	return lzc_raise_n( pZ, &pW->m_real, pTarget );
		case 'r':	return lzc_root_n( pZ, pW, pTarget );
		case 'c':	return lzc_root_3( pZ, pTarget );
		case 'e':	return lzc_exp_b( pZ, pTarget );
		case 'g':	return lzc_log_n( pW, pZ, pTarget );
		default:	return lzc_log_b( pZ, pTarget );
		}
	};

	struct { const char* m_szFamily; char m_cFunction; long m_lPlaces; long m_lBase; const char* m_szReal; const char* m_szImag;
			 const char* m_szWReal; const char* m_szWImag; const char* m_szWantReal; const char* m_szWantImag; } kPins[] =
	{
		{ "random", 'r', 50, 13, "2A886868", "CBA.549C8", "5", "0",
		  "35.8195C7B88984B43A77B0BB7972420945694334C396AABC4126",
		  "0.00031ACC4866053B24806768C34218CB5A7411548156361699" },
		{ "random", 'g', 50, 2, "-0.111101", "11.1101", "-101", "0",
		  "0.10100010011101010110011011101101100100111101010101",
		  "-0.00011100010010111011101001010101001000011110110000" },
		{ "random", 'p', 0, 15, "6BC17.27C", "-B.4C4D071", "0.D79", "0.771",
		  "1D851",
		  "1D3A" },
		{ "random", 'l', 36, 14, "-22D223BA", "-344567.12", "0", "0",
		  "7.42C7861CC182AA880CC32B970268DD583BDD",
		  "-1.28A813C3377115012760C7599544B554D667" },
		{ "random", 'p', 36, 12, "-E.113E1D4", "7231.DD12", "E.91", "-0.451",
		  "-610D26265875E1D6ED0232DEE1E9173493ED28E099D56.13D443E059D2D250556561D16E6166199460",
		  "521E18010793D19D82277656D9E5EE927E5E814DE6954.63E7042883D771714127222E152062468126" },
		{ "random", 'p', 1, 2, "-0.11101", "1111110", "0.11", "0",
		  "1110.0",
		  "100011.0" },
		{ "random", 'r', 3, 14, "6C682A1.A", "-4925046.9", "9.A", "0.A",
		  "6.1CB",
		  "-1.30B" },
		{ "random", 'e', 10, 10, "5.4355383", "-732267.4", "0", "0",
		  "-241065.0834249416",
		  "-127289.6135711699" },
		{ "random", 'p', 50, 10, "80632.462", "-78.445641", "0.301", "0",
		  "29.98192609770769712685744453156384672104820432510876",
		  "-0.00877980349451872862996824070882244535283467823444" },
		{ "random", 'p', 0, 15, "-DB.5A1E02", "0.04968ED2", "7.77", "0",
		  "-D741C086EAB87",
		  "-7B4B02D110C7801" },
		{ "random", 'p', 36, 11, "5669.3324", "737898.5", "0.524", "-71A",
		  "-124662822453113A1AA87026467761A830292AAA8388A50468053146707138324503578344A508362500458A421709031412A5A6A52A90A3A4A264A352055A9716A8932425A2965A929492416556A870932760526660669A45931962963055630663666490359457804978377307102A17889A63727A160695656148099287156054368975A16AA949A750930996A506313212868AA15984520342355616275470109A6184075A669365626A11815333641496015973844A7948435363684336233573A5A329A463077626375161675332011557365A38234433633A338381246148A450420A834946937A5A8395618112170320606437408373856041903A34568A64128311216318686978404462605825522A6985094A8535743A49.4331323A83340183769A41373A7409467A88",
		  "921AA08772877365160135A7343466877776084567137814327342825487461A146537870232885493002541424474A77043A9747542502582928180376022023910A050073777A027654095872181380009861A789370477482545371539829024422262375A5566535A53A1416A073590414357804660368647689165A1761157539570672A773298020A390937691601099934526792680A229163771830605A1016377439161689991555A559A21541355106170A31228510A6150A6347448A32676A250A22922297088803613774A8549178815244952805336795846A76802187991102AA3436784602112865336724A3512A155884071152774708809139A4098049A595839025750781897A846A0125475AA324193309358.858057493825573763A1A65405A651264585" },
		{ "random", 'c', 36, 5, "20442400", "-103001.23", "0", "0",
		  "210.241233401044204321044223334311431302",
		  "-0.142110130213323311041413321041223223" },
		{ "random", 'g', 0, 3, "-10.22021", "-2011120.1", "-120", "0",
		  "1",
		  "-2" },
		{ "random", 'r', 10, 6, "-44.442205", "524401.45", "43", "0.53",
		  "1.2524101550",
		  "0.0222055522" },
		{ "root ties", 'p', 10, 2, "1111.11111111111011011010000000000000100101100001", "-0.000110111111111111111010101001", "0.01", "0",
		  "10",
		  "-0.0000000100" },
		{ "root ties", 'p', 36, 12, "-0.4D4D69", "-0.D40008", "0.4", "0",
		  "0.93",
		  "-0.74" },
		{ "root ties", 'c', 3, 12, "367.E288E6", "9.E3E91514D76", "0", "0",
		  "8",
		  "0.076" },
		{ "root ties", 'p', 10, 9, "83237.8753", "-58.428888", "0.3", "0",
		  "42",
		  "-0.01" },
		{ "root ties", 'r', 10, 4, "3113.333333333333333331132", "0.00000022013333333333333333333022", "3", "0",
		  "12",
		  "0.0000000002" },
		{ "root ties", 'r', 36, 4, "32012321.2201", "20310213", "10", "0",
		  "100",
		  "2.1" },
		{ "near roots", 'r', 36, 14, "6A.51DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD9DA70000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "0.0000000000000000000000000000000000A5A", "2", "0",
		  "9.A00000000000000000000000000000000000",
		  "0.000000000000000000000000000000000007" },
		{ "near roots", 'r', 36, 8, "21.0077777777777777777777777777777777777777777777777777777777777777777742336000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", "0.0000000000000000000000000000000005467", "2", "0",
		  "4.100000000000000000000000000000000000",
		  "0.000000000000000000000000000000000053" },
		{ "near roots", 'r', 1, 12, "-1050.E6E5EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE", "-41.281760000000000000000000000000000000000000000001", "3", "0",
		  "6.2",
		  "-D.5" },
		{ "near roots", 'r', 1, 8, "6622.47347321000000000000000000000000000000000000000000001", "-3031.61505560000000000000000000000000000000000000000000001", "4", "0",
		  "7.7",
		  "-0.7" },
		{ "near roots", 'r', 0, 12, "-21E391491.5300000000000000000000000000001", "89DE0194.6", "4", "0",
		  "E3",
		  "97" },
		{ "near roots", 'r', 1, 8, "-1463.1044000000000001", "-25.5057", "3", "0",
		  "4.6",
		  "-10.1" },
		{ "unit halves", 'p', 0, 6, "0", "1", "1.2", "-0.0000000000000000001",
		  "-1",
		  "1" },
		{ "unit halves", 'p', 20, 9, "0", "-1", "-0.3", "0.000000000000000000000000000000000000000000000000000000001",
		  "0.77128835837054353204",
		  "0.44444444444444444445" },
		{ "unit halves", 'p', 0, 12, "-1", "0", "0.4", "-0.000000000000000000000000000001",
		  "1",
		  "1" },
		{ "unit halves", 'p', 0, 6, "0", "-1", "0.4", "0.0000000001",
		  "1",
		  "-1" },
		{ "unit roots", 'r', 1, 9, "-1", "0", "-3", "-0.000000000000000000000000000001",
		  "0.4",
		  "-0.8" },
		{ "unit roots", 'r', 0, 6, "0", "-1", "3", "-0.00000000000000000001",
		  "1",
		  "-1" },
		{ "unit roots", 'r', 0, 10, "-1", "0", "-6", "0.000000000001",
		  "1",
		  "-1" },
		{ "unit roots", 'r', 0, 12, "0", "-1", "3", "-0.00000000000000000001",
		  "1",
		  "-1" },
		{ "unit bound", 'r', 4, 10, "-1", "0", "6", "4.94871241408007544852563410181082407030903928708885116886441136246834458785",
		  "1.2308",
		  "0.3964" },
		{ "unit bound", 'r', 4, 10, "-1", "0", "6", "4.03551823181124411972512308592681626655264612938427462257486788016904758731",
		  "1.1925",
		  "0.4495" },
		{ "axis", 'c', 50, 13, "0", "1B31", "0", "0",
		  "10.B196A9546C3C3A4707601BA88781C3456323438416A6549718",
		  "8" },
		{ "axis", 'p', 10, 2, "0", "-1", "10", "11",
		  "-1101111.0101000101",
		  "0" },
		{ "axis", 'r', 1, 6, "-33.3", "0", "-2", "0",
		  "0",
		  "-0.1" },
		{ "axis", 'p', 36, 11, "0", "-1", "3", "A.5",
		  "0",
		  "7715A76.9AA398694752648301391746692734A76062" },
		{ "axis", 'p', 36, 2, "0", "1", "0", "-0.11",
		  "11.001111111000100100111100100011111110",
		  "0" },
		{ "axis", 'p', 1, 5, "-1", "0", "-3", "-1",
		  "-43.1",
		  "0" },
		{ "odd roots", 'c', 0, 8, "-773", "0", "0", "0",
		  "-10",
		  "0" },
		{ "odd roots", 'r', 36, 5, "-0.01", "0", "3", "0",
		  "-0.132333314134113424202430004141220144",
		  "0" },
		{ "odd roots", 'c', 50, 5, "-22", "0", "0", "0",
		  "-2.12104213000111041014200243324421330332314013001020",
		  "0" },
		{ "odd roots", 'c', 36, 11, "-0.869", "0", "0", "0",
		  "-0.A161118775A2576503A44718989A69746961",
		  "0" },
		{ "axis products", 'r', 50, 12, "0", "-1E53000", "3", "0",
		  "10E.D7465ED5844903D9D821D188E891269943747D98921E0291D0",
		  "-76" },
		{ "axis products", 'r', 0, 3, "0", "110201", "10", "0",
		  "20",
		  "10" },
		{ "axis products", 'r', 36, 11, "-11327753405", "0", "6", "0",
		  "44",
		  "25.792832631478023723A17A03917A88906062" },
		{ "axis products", 'r', 2, 8, "0", "22055000", "3", "0",
		  "221.40",
		  "124" },
		{ "logs", 'g', 10, 11, "7.9", "-0.01", "56.13AA", "-0.147",
		  "0.5555555556",
		  "0" },
		{ "logs", 'g', 1, 3, "11.002", "-2.202", "2.2", "-1.01",
		  "1.2",
		  "0" },
		{ "logs", 'g', 36, 9, "418.86", "-35.581", "53.85", "-3.1",
		  "1.444444444444444444444444444444444444",
		  "0" },
		{ "logs", 'g', 50, 7, "10", "0.06", "66.6616", "1.5",
		  "0.33333333333333333333333333333333333333333333333334",
		  "0" },
		{ "logs", 'g', 10, 7, "51", "0.04", "3530.6645", "5.61",
		  "0.3333333334",
		  "0" },
		{ "logs", 'g', 0, 13, "63A4.B518", "56.75CAC8", "353.CC93", "1.CB",
		  "2",
		  "0" },
		{ "steered", 'r', 10, 16, "1.68ACD0C866D64A136D74F22171F07CDB2DC07FBF0CFC85EE6F9099FBD77EC9938D86C17D9E88BBE7", "B6.3", "3", "0",
		  "4.EAB5A3F309",
		  "2.D2658B0A61" },
		{ "steered", 'p', 10, 8, "1562334224.46737060777674151271301273643121776177033020763446245474373656055646077313036751", "-37.7", "-0.35", "-27",
		  "-0.0000524055",
		  "0.0000000001" },
		{ "steered", 'p', 10, 16, "-C57", "-1.6A16BBFB78C8D8838D61EB4EE22CD29820AC3C48F3EDDA5A7156F0BD49ACE37AC140D4D60C5D205E", "0.3", "-3.6",
		  "-0.0006852B9F",
		  "-0.0003796BAB" },
		{ "steered", 'p', 10, 14, "-56.6", "-0.9639B8BD21045A60D4CBDD0AAC9844C61DB4BCD672AD0A3D740A068BAB951698C3D68559B615C529", "-0.2", "-6.2",
		  "0.000000019D",
		  "-0.0000000307" },
		{ "log host", 'l', 3, 16, "-0.EB855B18FAA94DC294EF6D31C1282914539D0CE3E52CE0E503954350B1A54537FA9B97F1A", "-8.9D", "0", "0",
		  "0.C75",
		  "-0.9AE" },
		{ "log host", 'l', 10, 14, "-9.4A", "-1.C6C32848B69733C2CC30B445BC39734970B71014B8531D9A34076A2542C3DD84696B2C03738D4D17", "0", "0",
		  "0.BD5A0561AD",
		  "-1.18688107A5" },
		{ "noughts", 'p', 50, 12, "0", "-1", "2", "-9E",
		  "0.00000000000000000000000000000000000000000000000000",
		  "0" },
		{ "noughts", 'r', 1, 5, "0", "-1", "0", "0.00001",
		  "0.0",
		  "0" },
		{ "large", 'p', 20, 3, "22110000", "1200", "1202", "2.2",
		  "12020200200110011212101000222001122121001101001010121220112211102010201101012201202121221010112012121010100021221122210220100120211210122101202210112211121212210121000112222011121200000122220101010110211001202221001102001122100001021122200102220210010210021221210100100202200120122101111201210001021002100102201210010122202001112011222210102001012212222112121100222212121.12122102110000111102",
		  "-10222100210101011212212110110100201001012002122022120121221100000112010202210122212102211101121020101000122001100122022121011010022202121200100121022011111220102210110121012120020222021111000000022200020112011011222222100221100000210211002102102111101222201111101222202111100201100111020220112012110121202002201012002001110200020110221120222012200200011001210011102201111102.02100002120220222112" },
		{ "large", 'e', 36, 3, "11011.01", "21.1", "0", "0",
		  "-200101222021111210202220211220000112020210111221120022201002210102200011110012102220220000020020020201001211201.000222112201111122100010121010222022",
		  "10022122110001212002102010110120101202212021112121002121100022212112210120100010001101011022012011022111101200100.122011211110110201220020020121222120" },
		{ "large", 'e', 36, 14, "12", "D6.4", "0", "0",
		  "C18342610284024A.7B1106CDDB3B3B38B3892AB981BD03411D7C",
		  "7039A7A53D541615.D1D1B9073C402948644D1D4B4C7291187CB7" },
		{ "large exponents", 'r', 20, 10, "1.000000000003", "0.000000000002", "0.0000000000001", "0",
		  "4360958581221.25382716256603034899",
		  "9756166215718.70823894631312217375" },
		{ "log noughts", 'g', 1, 3, "102.2", "0", "0", "1",
		  "0",
		  "-1.2" },
		{ "log noughts", 'g', 36, 5, "-1", "0", "-0.3", "0.4",
		  "1.202133142144443344200332442333403324",
		  "0" },
		{ "log noughts", 'g', 36, 9, "8616", "0", "0", "-1",
		  "0",
		  "5.513425213648216280354186736485471140" },
		{ "log noughts", 'g', 3, 10, "-1", "0", "0.6", "0.8",
		  "3.388",
		  "0" },
		{ "exact", 'p', 20, 10, "2", "1", "3", "0",
		  "2",
		  "11" },
		{ "exact", 'c', 20, 10, "2", "11", "0", "0",
		  "2",
		  "1" },
		{ "exact", 'r', 20, 10, "3", "4", "2", "0",
		  "2",
		  "1" },
		{ "exact", 'r', 20, 10, "-4", "0", "2", "0",
		  "0",
		  "2" },
		{ "exact", 'p', 20, 3, "-22", "0", "0.1", "0",
		  "1",
		  "1.20120212222121022122" },
		{ "exact", 'c', 20, 10, "0", "8", "0", "0",
		  "1.73205080756887729353",
		  "1" },
		{ "exact", 'p', 20, 10, "0", "1", "0", "1",
		  "0.20787957635076190855",
		  "0" },
		{ "exact", 'p', 20, 10, "-1", "0", "0", "1",
		  "0.04321391826377224977",
		  "0" },
		{ "exact", 'g', 20, 10, "-1", "0", "0", "1",
		  "2",
		  "0" },
		{ "exact", 'l', 20, 10, "6", "8", "0", "0",
		  "1",
		  "0.40271919627337314208" },
		{ "exact", 'r', 0, 3, "-1000", "0", "20", "0",
		  "2",
		  "1" },
		{ "exact", 'r', 20, 3, "0", "0.02", "-2", "0",
		  "1.11111111111111111112",
		  "-1.11111111111111111112" },
		{ "exact", 'n', 20, 10, "1", "1", "-3", "0",
		  "-0.25",
		  "-0.25" },
		{ "exact", 'c', 20, 10, "-8", "0", "0", "0",
		  "-2",
		  "0" },
		{ "exact", 'l', 20, 10, "-100", "0", "0", "0",
		  "2.00000000000000000000",
		  "1.36437635384184134749" },
		{ "exact", 'g', 20, 10, "-0.5", "-0.5", "1", "1",
		  "-2.67403202783654016330",
		  "-0.73870212227295079460" },
		{ "exact", 'g', 20, 10, "0", "-1", "0", "1",
		  "-1",
		  "0" },
		{ "exact", 'g', 20, 10, "9", "0", "0", "1",
		  "0",
		  "-1.39879661026423911199" },
		{ "exact", 'g', 20, 10, "0.6", "0.8", "0", "1",
		  "0.59033447060173309670",
		  "0" },
		{ "exact", 'r', 1, 3, "0", "-11.122", "10", "0",
		  "1.1",
		  "-0.2" },
		{ "exact", 'p', 2, 3, "-11.122", "0", "0.2", "0",
		  "-1.10",
		  "2.11" },
	};

	for ( const char* sz_family : { "random", "exact", "root ties", "near roots", "unit halves", "unit roots", "unit bound", "axis",
									"odd roots", "axis products", "logs", "steered", "log host", "noughts", "large", "large exponents",
									"log noughts" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZComplex	z;
			CLZComplex	w;
			CLZComplex	t;

			glWorkingPrecision = pin.m_lPlaces;
			SetComplex( z, pin.m_lBase, pin.m_szReal, pin.m_szImag );
			SetComplex( w, pin.m_lBase, pin.m_szWReal, pin.m_szWImag );

			long	l_code = Call( pin.m_cFunction, &z, &w, &t );

			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &t.m_real, pin.m_szWantReal ) || !SpelledAs( &t.m_imag, pin.m_szWantImag ) ||
				 (t.m_real.m_lBase != pin.m_lBase) || (t.m_imag.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::string( 1, pin.m_cFunction ) + " base " + std::to_string( pin.m_lBase ) + " " + pin.m_szReal + " " +
							  pin.m_szImag + " " + pin.m_szWReal + " " + pin.m_szWImag;
			}
		}

		if ( l_cases == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "no pins" );
		else if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "powers, roots and logarithms, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), (l_cases > 0) && (l_wrong == 0), szDetail );
	}

	// On the real axis each is the real function's answer to the character,
	//	an exact nought beside it -- a base not below nought or a whole
	//	exponent for the power, a value not below nought or an odd order for
	//	the root, a base and a value past nought for the logarithms. A real
	//	function refused as infinite is the point at infinity, with its
	//	code...
	auto	SameOnTheAxis = []( CLZComplex* pGot, long lGot, CLZNumber* pWant, long lWant ) -> bool
	{
		if ( pWant->FlagIsOn( kfInfinity ) )
			return (lGot == lWant) && pGot->m_real.FlagIsOn( kfInfinity ) && pGot->m_imag.FlagIsOn( kfInfinity );

		return SameInEveryRespect( &pGot->m_real, lGot, pWant, lWant ) && SpelledAs( &pGot->m_imag, "0" );
	};

	{
		std::uint64_t	u_state = 5255;
		long			l_cases = 0;
		long			l_wrong = 0;
		long			l_parts = 0;
		long			l_partsWrong = 0;

		for ( long i = 0; i < 48; i++ )
		{
			long	l_base = 2 + (i % 15);

			glWorkingPrecision = (i % 3 == 0) ? 20 : 36;

			CLZNumber	v;
			CLZNumber	small;
			CLZNumber	order;
			CLZNumber	logBase;
			CLZNumber	whole;
			CLZNumber	odd;
			CLZNumber	nought;
			CLZNumber	one;

			for ( CLZNumber* p : { &v, &small, &order, &logBase, &whole, &odd, &nought, &one } )
				p->SetBase( l_base, false );

			v.Assign( RandomDigitsIn( u_state, l_base, 1 + static_cast<long>( (u_state >> 40) % 10 ), 4, false ).c_str() );
			small.Assign( RandomDigitsIn( u_state, l_base, 2, 1, (i % 2) == 1 ).c_str() );
			order.Assign( RandomDigitsIn( u_state, l_base, 2, 1, (i % 4) >= 2 ).c_str() );
			logBase.Assign( RandomDigitsIn( u_state, l_base, 2, 1, false ).c_str() );
			whole.Assign( static_cast<long>( (u_state >> 35) % 25 ) - 12 );
			odd.Assign( 2 * static_cast<long>( (u_state >> 37) % 6 ) - 5 );
			nought.Assign( 0L );
			one.Assign( 1L );

			if ( lz_compare( &logBase, &one ) == kEqualTo )
				logBase.Assign( 2L );

			CLZNumber	negative;

			negative.Assign( &v );
			lz_negate( &negative );

			// (function, value, w) -- the real call is the same letter's...
			struct { char m_cFunction; CLZNumber* m_pValue; CLZNumber* m_pW; } kOnTheAxis[] =
			{
				{ 'p', &v, &small }, { 'p', &negative, &whole }, { 'n', &negative, &whole }, { 'n', &v, &whole },
				{ 'r', &v, &order }, { 'r', &negative, &odd }, { 'c', &v, nullptr }, { 'c', &negative, nullptr },
				{ 'e', &small, nullptr }, { 'g', &v, &logBase }, { 'l', &v, nullptr },
			};

			for ( auto& on : kOnTheAxis )
			{
				CLZComplex	z;
				CLZComplex	w;
				CLZComplex	t;
				CLZNumber	want;
				long		l_want = kNoError;

				z.m_real.Assign( on.m_pValue );
				z.m_imag.Assign( &nought );
				w.m_real.Assign( (on.m_pW != nullptr) ? on.m_pW : &nought );
				w.m_imag.Assign( &nought );

				long	l_got = Call( on.m_cFunction, &z, &w, &t );

				switch ( on.m_cFunction )
				{
				case 'p':	l_want = lz_power( on.m_pValue, on.m_pW, &want ); break;
				case 'n':	l_want = lz_raise_n( on.m_pValue, on.m_pW, &want ); break;
				case 'r':	l_want = lz_root_n( on.m_pValue, on.m_pW, &want ); break;
				case 'c':	l_want = lz_root_3( on.m_pValue, &want ); break;
				case 'e':	l_want = lz_exp_b( on.m_pValue, &want ); break;
				case 'g':	l_want = lz_log_n( on.m_pW, on.m_pValue, &want ); break;
				default:	l_want = lz_log_b( on.m_pValue, &want ); break;
				}

				l_wrong += SameOnTheAxis( &t, l_got, &want, l_want ) ? 0 : 1;
				l_cases++;
			}

			// ...and on the negative real axis and the imaginary, a real-based
			//	logarithm's real part is the real logarithm of |z|, its
			//	spelling included...
			CLZComplex	logBaseZ;

			logBaseZ.m_real.Assign( &logBase );
			logBaseZ.m_imag.Assign( &nought );

			for ( long l_axis = 0; l_axis < 3; l_axis++ )
			{
				CLZComplex	z;
				CLZComplex	t;
				CLZNumber	want;

				z.m_real.Assign( (l_axis == 0) ? &negative : &nought );
				z.m_imag.Assign( (l_axis == 0) ? &nought : ((l_axis == 1) ? &v : &negative) );

				long	l_got = lzc_log_b( &z, &t );
				long	l_want = lz_log_b( &v, &want );

				l_partsWrong += SameInEveryRespect( &t.m_real, l_got, &want, l_want ) ? 0 : 1;
				l_got = lzc_log_n( &logBaseZ, &z, &t );
				l_want = lz_log_n( &logBase, &v, &want );
				l_partsWrong += SameInEveryRespect( &t.m_real, l_got, &want, l_want ) ? 0 : 1;
				l_parts += 2;
			}
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld answers differ", l_wrong, l_cases );
		CheckTrue( "on the real axis, the real function's", l_wrong == 0, szDetail );
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld real parts differ", l_partsWrong, l_parts );
		CheckTrue( "a real-based log's real part on an axis", l_partsWrong == 0, szDetail );
	}

	glWorkingPrecision = 20;

	// The odd root of a negative real is the real root, and an even one the
	//	principal value, exactly where it is exact...
	{
		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	t;
		bool		b_ok = true;

		SetComplex( z, 10, "-8", "0" );
		b_ok = b_ok && (lzc_root_3( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "-2" ) && SpelledAs( &t.m_imag, "0" );
		SetComplex( w, 10, "3", "0" );
		b_ok = b_ok && (lzc_root_n( &z, &w, &t ) == kNoError) && SpelledAs( &t.m_real, "-2" ) && SpelledAs( &t.m_imag, "0" );
		SetComplex( z, 10, "-4", "0" );
		SetComplex( w, 10, "2", "0" );
		b_ok = b_ok && (lzc_root_n( &z, &w, &t ) == kNoError) && SpelledAs( &t.m_real, "0" ) && SpelledAs( &t.m_imag, "2" );

		CheckTrue( "an odd root of a negative real", b_ok, b_ok ? "the real root, an even one 2i" : "an answer wrong" );
	}

	// The angle in radians whatever the unit: none of the seven reads it...
	{
		long	l_wrong = 0;

		for ( char c_fn : { 'p', 'n', 'r', 'c', 'e', 'g', 'l' } )
		{
			CLZComplex	z;
			CLZComplex	w;
			CLZComplex	radians;
			CLZComplex	degrees;

			SetComplex( z, 10, "0.75", "-1.25" );
			SetComplex( w, 10, (c_fn == 'n') ? "3" : "1.5", (c_fn == 'n') ? "0" : "0.5" );
			Call( c_fn, &z, &w, &radians );
			geAngleUnit = AU_degrees;
			Call( c_fn, &z, &w, &degrees );
			geAngleUnit = AU_radians;
			l_wrong += (SameAnswer( &radians.m_real, &degrees.m_real ) && SameAnswer( &radians.m_imag, &degrees.m_imag )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 7 differ", l_wrong );
		CheckTrue( "powers and logarithms in radians", l_wrong == 0, szDetail );
	}

	// Nought to a power, anything to nought, and the refusals...
	{
		CLZComplex	z;
		CLZComplex	w;
		CLZComplex	t;
		bool		b_ok = true;

		auto	Both = []( CLZComplex& t, unsigned long lFlag ) -> bool { return t.m_real.FlagIsOn( lFlag ) && t.m_imag.FlagIsOn( lFlag ); };

		SetComplex( z, 10, "0", "0" );
		SetComplex( w, 10, "2", "1" );
		b_ok = b_ok && (lzc_power( &z, &w, &t ) == kNoError) && SpelledAs( &t.m_real, "0" ) && SpelledAs( &t.m_imag, "0" );
		SetComplex( w, 10, "-2", "1" );
		b_ok = b_ok && (lzc_power( &z, &w, &t ) == kDivideByZero) && Both( t, kfInfinity ) && Both( t, kfDiv0 );
		SetComplex( w, 10, "0", "1" );
		b_ok = b_ok && (lzc_power( &z, &w, &t ) == kFailed) && Both( t, kfQNaN );
		SetComplex( z, 10, "2", "1" );
		SetComplex( w, 10, "0", "0" );
		b_ok = b_ok && (lzc_power( &z, &w, &t ) == kNoError) && SpelledAs( &t.m_real, "1" ) && SpelledAs( &t.m_imag, "0" );

		SetComplex( w, 10, "0.5", "0" );
		b_ok = b_ok && (lzc_raise_n( &z, &w.m_real, &t ) == kFailed) && Both( t, kfQNaN );
		b_ok = b_ok && (lzc_raise_n( &z, nullptr, &t ) == kFailed);
		SetComplex( w, 10, "0", "0" );
		b_ok = b_ok && (lzc_root_n( &z, &w, &t ) == kFailed) && Both( t, kfQNaN );

		CLZComplex	n;

		SetComplex( n, 10, "2", "1" );
		SetComplex( z, 10, "0", "0" );
		b_ok = b_ok && (lzc_log_n( &n, &z, &t ) == kUndefinedResult) && Both( t, kfInfinity );
		b_ok = b_ok && (lzc_log_b( &z, &t ) == kUndefinedResult) && Both( t, kfInfinity );
		SetComplex( z, 10, "2", "1" );
		SetComplex( n, 10, "0", "0" );
		b_ok = b_ok && (lzc_log_n( &n, &z, &t ) == kUndefinedResult) && Both( t, kfInfinity );
		SetComplex( n, 10, "1", "0" );
		b_ok = b_ok && (lzc_log_n( &n, &z, &t ) == kDivideByZero) && Both( t, kfInfinity ) && Both( t, kfDiv0 );

		// ...base one whatever z, past nought or of modulus one too, where
		//	the real logarithm beneath would answer otherwise...
		for ( const char* sz_real : { "2", "-1", "0.6" } )
		{
			SetComplex( z, 10, sz_real, (sz_real[ 0 ] == '0') ? "0.8" : "0" );
			b_ok = b_ok && (lzc_log_n( &n, &z, &t ) == kDivideByZero) && Both( t, kfInfinity ) && Both( t, kfDiv0 );
		}

		// ...too large, the point at infinity...
		SetComplex( z, 10, "10", "1" );
		SetComplex( w, 10, "10000000", "0" );
		b_ok = b_ok && (lzc_power( &z, &w, &t ) == kTooLarge) && Both( t, kfInfinity );
		SetComplex( w, 10, "0.0000001", "0" );
		b_ok = b_ok && (lzc_root_n( &z, &w, &t ) == kTooLarge) && Both( t, kfInfinity );
		SetComplex( z, 10, "10000000", "1" );
		b_ok = b_ok && (lzc_exp_b( &z, &t ) == kTooLarge) && Both( t, kfInfinity );

		CheckTrue( "nought, one, refusals, limits", b_ok, b_ok ? "as the real functions and the layer answer" : "an answer wrong" );
	}

	// A not-a-number or an infinity in either operand answers a not-a-number,
	//	and an operand in two bases is refused...
	{
		CLZNumber	nan;
		CLZNumber	infinite;
		bool		b_ok = true;

		MakeNotANumber( nan );
		MakeInfinite( infinite );

		for ( char c_fn : { 'p', 'n', 'r', 'c', 'e', 'g', 'l' } )
		{
			bool	b_binary = (c_fn == 'p') || (c_fn == 'r') || (c_fn == 'g');

			for ( long l_where = 0; l_where < ((b_binary || (c_fn == 'n')) ? 4 : 2); l_where++ )
			{
				CLZComplex	z;
				CLZComplex	w;
				CLZComplex	t;

				SetComplex( z, 10, "1", "2" );
				SetComplex( w, 10, "3", "0" );

				CLZComplex*	p_into = (l_where < 2) ? &z : &w;

				CLZNumber*	p_part = ((c_fn == 'n') && (l_where >= 2)) ? &w.m_real : &p_into->m_imag;

				p_part->Assign( ((l_where % 2) == 0) ? &nan : &infinite );
				b_ok = b_ok && (Call( c_fn, &z, &w, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );
			}

			for ( long l_where = 0; l_where < (b_binary ? 2 : 1); l_where++ )
			{
				CLZComplex	z;
				CLZComplex	w;
				CLZComplex	t;

				SetComplex( z, 10, "1", "2" );
				SetComplex( w, 10, "3", "0" );

				CLZComplex*	p_into = (l_where == 0) ? &z : &w;

				p_into->m_imag.SetBase( 16, false );
				p_into->m_imag.Assign( "2" );
				b_ok = b_ok && (Call( c_fn, &z, &w, &t ) == kBaseMismatch) && t.m_real.FlagIsOn( kfQNaN );
			}
		}

		CheckTrue( "flags and two bases", b_ok, b_ok ? "carried or refused as the layer does" : "an answer wrong" );
	}

	// A logarithm's value in another base is brought into the base's, and
	//	answers there as it would written in it...
	{
		bool	b_ok = true;

		for ( const char* sz_base : { "2", "0.5" } )
		{
			CLZComplex	n;
			CLZComplex	z;
			CLZComplex	same;
			CLZComplex	t;
			CLZComplex	u;

			SetComplex( n, 10, sz_base, (sz_base[ 0 ] == '2') ? "0" : "1" );
			SetComplex( z, 16, "0", "8" );
			SetComplex( same, 10, "0", "8" );
			lzc_log_n( &n, &z, &t );
			lzc_log_n( &n, &same, &u );
			b_ok = b_ok && (t.m_real.m_lBase == 10) && (t.m_imag.m_lBase == 10) && SameAnswer( &t.m_real, &u.m_real ) &&
				   SameAnswer( &t.m_imag, &u.m_imag );
		}

		CheckTrue( "a log's value in another base", b_ok, b_ok ? "answered in the base's base, as written in it" : "an answer differs" );
	}

	// In place, as from a copy -- lzc_log_n's target the base, its first
	//	operand -- and an operand that is also the other...
	{
		long	l_wrong = 0;

		for ( char c_fn : { 'p', 'n', 'r', 'c', 'e', 'g', 'l' } )
		{
			CLZComplex	z;
			CLZComplex	w;
			CLZComplex	copy;
			CLZComplex	t;

			SetComplex( z, 7, "-3.25", "1.5" );
			SetComplex( copy, 7, "-3.25", "1.5" );
			SetComplex( w, 7, (c_fn == 'n') ? "3" : "1.5", (c_fn == 'n') ? "0" : "-0.25" );
			Call( c_fn, &copy, &w, &t );

			if ( c_fn == 'g' )
			{
				CLZComplex	n;

				SetComplex( n, 7, "1.5", "-0.25" );
				lzc_log_n( &n, &z, nullptr );
				z.m_real.Assign( &n.m_real );
				z.m_imag.Assign( &n.m_imag );
			}
			else
				Call( c_fn, &z, &w, nullptr );

			l_wrong += (SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag )) ? 0 : 1;
		}

		CLZComplex	z;
		CLZComplex	copy;
		CLZComplex	other;
		CLZComplex	t;

		SetComplex( z, 7, "1.25", "-0.5" );
		SetComplex( copy, 7, "1.25", "-0.5" );
		SetComplex( other, 7, "1.25", "-0.5" );
		lzc_power( &copy, &other, &t );
		lzc_power( &z, &z, nullptr );
		l_wrong += (SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag )) ? 0 : 1;

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 8 differ", l_wrong );
		CheckTrue( "powers and logarithms in place", l_wrong == 0, szDetail );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**The inverse functions that take square roots, stage 6a** --
	2026-10-23, on Rich's calls: lzc_arcsine, lzc_arccosine, lzc_arcsecant,
	lzc_arccosecant and their hyperbolic twins, the four that take
	logarithms left to 6b. Every expected part is generated -- Python's
	decimal from the logarithmic forms, deep enough for the case's own
	places, and exact Fractions where a part in a unit is rational -- never
	typed. Each family past the random is made of cases a mutant gets
	wrong: the axes by the builds that take a reciprocal's cut from the
	wrong side and leave arccsch past i unsigned, the axis signs by the
	builds that negate arcsin outside its domain and take arccos 1 for
	arccosh below -1, the axis units by the build that writes an axis part
	in radians whatever the unit, the cuts by the build that loses a root's
	imaginary sign, the branch points by the builds that do not widen for a
	small root's argument or for |z|^2, the large by the build without a
	guard, the unit exact by the build that never finds a real part of 30
	or 45 degrees exact, the unit signs by the builds that exchange 60 and
	120 and leave the exact arcsin part unsigned, the near unit exact and
	the steered by the build without a retake, the near rational by the
	build without the tie decider, the near segment by the builds that take
	each part's sign wrong in that decider -- a small part along an axis at
	a Pythagorean point, where its derivative is rational -- the steered
	acosh by the builds that bound cosh's series too soon or drop its
	remainder, the decider by the builds that ignore its sign, do not
	reverse acos and drop arcsinh's real sign, and the root pair by the
	build that takes the wrong branch of its exact comparison.
*/

static void TestComplexInverses()
{
	Banner( "complex numbers, stage 6a: the inverse functions that take square roots" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	struct { const char* m_szFamily; long m_lFunction; long m_lUnit; long m_lPlaces; long m_lBase; const char* m_szReal;
			 const char* m_szImag; const char* m_szWantReal; const char* m_szWantImag; } kPins[] =
	{
		{ "random", 6, 0, 0, 5, "-3.434", "0.0441",
		  "0",
		  "-2" },
		{ "random", 4, 0, 0, 12, "-0.0934", "0.09074",
		  "0",
		  "0" },
		{ "random", 3, 0, 10, 10, "94", "205.9",
		  "0.0018348196",
		  "-0.0040190618" },
		{ "random", 3, 0, 36, 15, "-0.117", "-E3520000000",
		  "-0.0000000000000000000000134AE89144A6DC",
		  "0.000000000010C464E3DAB7BE67C23A8E1A90" },
		{ "random", 4, 0, 50, 13, "-717", "-6.C27",
		  "-7.A2A35BA5743C14BC13251BAA606B3C10A3567C375C832A0C9B",
		  "-0.00C89463100574B435974853B24144BB1680681C7C2B75C4B4" },
		{ "random", 3, 1, 3, 5, "31400", "-1.114",
		  "0.003",
		  "0.000" },
		{ "random", 6, 0, 1, 4, "-11.1", "-0.2123",
		  "0.0",
		  "1.3" },
		{ "random", 2, 0, 10, 13, "78.B5", "617.2",
		  "1.755AB432CA",
		  "0.0021493625" },
		{ "random", 1, 0, 0, 6, "45450000", "45300000000",
		  "2",
		  "-32" },
		{ "random", 2, 1, 3, 7, "-263.5", "3.516",
		  "156.254",
		  "0.004" },
		{ "random", 2, 0, 36, 6, "0.51", "-1200000000",
		  "1.323143023505342032122422024342440332",
		  "-0.000000000425555555555555555341124152" },
		{ "random", 0, 1, 3, 8, "16170000000", "-66.74",
		  "132.000",
		  "-2360.123" },
		{ "random", 3, 0, 20, 2, "-10", "-110",
		  "-0.00001100101010011100",
		  "0.00100110010011011101" },
		{ "random", 7, 0, 36, 11, "-716A0000000", "4794000000",
		  "-0.00000000001595A226832736342980431945",
		  "-0.000000000001013985110A40764485753756" },
		{ "random", 2, 3, 20, 16, "0.AD6E", "0.49D7",
		  "54.9BDDF7F2E7BF44E52D2F",
		  "91.C7825EB5E29DEDDD6B27" },
		{ "random", 2, 2, 1, 16, "-A31B", "DB.34",
		  "64.0",
		  "0.0" },
		{ "axes", 3, 1, 50, 6, "0.3", "0",
		  "230",
		  "-203.24230510015351314252155103002452235300244101434351" },
		{ "axes", 3, 3, 0, 16, "0.838", "0",
		  "FA",
		  "-CD" },
		{ "axes", 7, 0, 1, 12, "0", "11.4",
		  "0",
		  "-0.1" },
		{ "axes", 7, 0, 0, 14, "0", "2",
		  "0",
		  "-1" },
		{ "axis signs", 0, 0, 1, 8, "-1.4", "0",
		  "-1.5",
		  "1.0" },
		{ "axis signs", 0, 1, 1, 8, "1.4", "0",
		  "132",
		  "67.1" },
		{ "axis signs", 5, 0, 20, 7, "-1320", "0",
		  "6.62604355531050165123",
		  "3.06636514320361341103" },
		{ "axis signs", 5, 0, 3, 5, "-2", "0",
		  "1.130",
		  "3.033" },
		{ "axis units", 0, 3, 10, 16, "0", "-2",
		  "0",
		  "-E5.C30034CE84" },
		{ "axis units", 2, 1, 50, 3, "0", "0.22",
		  "10100",
		  "2001.10200020202002220102110001001101200101022001211101" },
		{ "axis units", 1, 3, 36, 15, "0", "1",
		  "11A",
		  "-95.41D0015C3979B79266196EA18401A8E2527C" },
		{ "cuts", 7, 0, 20, 7, "-0.0000000001", "2",
		  "-0.00000000002010052250",
		  "-0.34441060513410152345" },
		{ "cuts", 1, 0, 0, 2, "0.00001", "-101",
		  "10",
		  "10" },
		{ "cuts", 5, 0, 36, 14, "0.0000001", "-1.7",
		  "1.2A260492825817996C200DD25BC7D348D25B",
		  "-1.7DC39D64747587A47C7A723A197994C3B8CC" },
		{ "cuts", 1, 0, 36, 15, "-5", "0.01",
		  "3.219C3008817DDB9CA4EDDA450E8C7CAD3B90",
		  "-2.45BE5940C182E48DCB54D7D6ADC8A6B0D147" },
		{ "branch points", 4, 0, 0, 15, "-0.0000000000000000000000000001", "-1",
		  "0",
		  "-2" },
		{ "branch points", 4, 0, 0, 11, "-0.00000000000000000000000001", "1",
		  "0",
		  "2" },
		{ "branch points", 2, 0, 50, 6, "0.00000000000000000000000000001", "0.00000000000000000000001",
		  "1.32314202350534203301520502434232522001525420040002",
		  "105.52310302055450302102135050041150534034214533514134" },
		{ "branch points", 2, 0, 20, 7, "-0.0000000000001", "-0.0000000000000001",
		  "3.06536514342613550454",
		  "-34.66336336053056201145" },
		{ "large", 6, 0, 36, 10, "451000000000000000000", "-245",
		  "0.000000000000000000000000000000000000",
		  "1.570796326794896619229104396739529713" },
		{ "large", 5, 0, 10, 3, "12100000000000000", "10000000000000000000000000",
		  "1001.0110211120",
		  "1.1201020011" },
		{ "large", 6, 0, 3, 6, "13000000000000", "53400000",
		  "0.000",
		  "-1.323" },
		{ "unit exact", 1, 2, 20, 14, "-1.037", "0.A37",
		  "AA",
		  "-41.992B90CBB445C02D869D" },
		{ "unit exact", 1, 3, 50, 12, "-1.7", "-1.5",
		  "273",
		  "172.0127921DE66D6D63D97760150258904787D353949186695962" },
		{ "unit exact", 1, 2, 36, 12, "1.7", "1.5",
		  "42",
		  "-78.005D380922767673E63072544D711267D7E4" },
		{ "unit exact", 0, 1, 1, 8, "0.6", "-0.2",
		  "55",
		  "-23.7" },
		{ "unit exact", 0, 3, 50, 12, "0.86", "0.06",
		  "D5",
		  "9.4583E6334202151E78266E773220E663D521429553E9740430" },
		{ "unit exact", 1, 1, 1, 10, "0.825", "-0.425",
		  "45",
		  "32.6" },
		{ "unit signs", 1, 1, 36, 16, "1", "-1.8",
		  "3C",
		  "4B.74C4E3A03AAFBA8017EA1F02208D630BCCC8" },
		{ "unit signs", 1, 3, 50, 8, "-15", "26.4",
		  "515.25252525252525252525252525252525252525252525252525",
		  "-1164.63211550703644155715303447601147036423440317772511" },
		{ "unit signs", 0, 2, 1, 6, "-1.33", "-1.23",
		  "-122",
		  "-232.0" },
		{ "unit signs", 0, 3, 0, 16, "-1", "-1.8",
		  "-53",
		  "-D2" },
		{ "near rational", 7, 0, 34, 10, "512000000000000000000000000000000", "384000000000000000000000000000000",
		  "0.0000000000000000000000000000000013",
		  "-0.0000000000000000000000000000000009" },
		{ "near rational", 1, 0, 13, 12, "0.0000000000000000007D", "0.000000000009E6",
		  "1.6D2404847ED6D",
		  "-0.000000000009E" },
		{ "near rational", 2, 0, 6, 10, "960000", "-1280000",
		  "1.570796",
		  "-0.000001" },
		{ "near rational", 6, 0, 24, 10, "-38400000000000000000000", "51200000000000000000000",
		  "0.000000000000000000000013",
		  "-1.570796326794896619231331" },
		{ "near rational", 2, 0, 20, 10, "768000000000000000", "1024000000000000000",
		  "1.57079632679489661876",
		  "0.00000000000000000063" },
		{ "near rational", 0, 0, 10, 8, "0.0000000000000000000065", "0.00000000234",
		  "0.0000000000",
		  "0.0000000023" },
		{ "near segment", 6, 0, 17, 10, "0.8", "-0.00000000000387082800000000005",
		  "0.69314718055994531",
		  "0.00000000000806423" },
		{ "near segment", 0, 0, 23, 10, "0.28", "-0.00000000000000000581023200000000000000000000005",
		  "0.28379410920832784562570",
		  "-0.00000000000000000605233" },
		{ "near segment", 1, 0, 26, 10, "0.96", "0.00000000000000000000072028180000000000000000000000000006",
		  "0.28379410920832784562570323",
		  "-0.00000000000000000000257244" },
		{ "near segment", 2, 0, 27, 10, "2.125", "-0.00000000000000000000090400089843750000000000000000000000001",
		  "1.080839000541168310887156729",
		  "-0.000000000000000000000226887" },
		{ "near segment", 3, 0, 16, 10, "-2.125", "0.00000000000000019921875000000000000000000009",
		  "-0.4899573262537283",
		  "-0.0000000000000001" },
		{ "near segment", 5, 0, 21, 4, "2.02", "-0.00000000000000000010131000000000000000000000000000003",
		  "1.120232100233323322032",
		  "-0.000000000000000000022" },
		{ "near segment", 7, 0, 25, 10, "2.4", "-0.00000000000000000268431976800000000000000000000004",
		  "0.4054651081081643819780131",
		  "0.0000000000000000004301795" },
		{ "near segment", 2, 0, 33, 10, "0.28", "0.00000000000000000000000789962357654400000000000000000000000009",
		  "0.000000000000000000000029388480568",
		  "1.945910149055313305105352743443180" },
		{ "steered", 3, 3, 36, 12, "-723.EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE5ED51550E3601712166EE0DD78980D921D13D9062ED6", "-D6.8",
		  "-0.199677E515007934D904488399E85DD08487",
		  "0.027E8E54070D12374600D6D94584D0450673" },
		{ "steered", 5, 0, 10, 10, "-1.17999999796914894432984241812042231729835948", "52.3",
		  "4.6504892580",
		  "1.5933505221" },
		{ "steered", 6, 0, 20, 8, "2.010000000000000000037756027531576277265463362674271774645202", "-5.54",
		  "0.11761416321336634747",
		  "1.41016612166172320206" },
		{ "steered", 0, 0, 50, 15, "-145.0000000000000000000000000000000000000000000000072B86B0695EA41A1E44149E4736D36D6BB32066CE299AAB52398E8399", "C.4E",
		  "-1.7DD0912613914D9177021CAE866B560E7228DA025771EB7EBE",
		  "6.56D3CCA1ABEA11750A61DB397669E773418D5C9216965653C4" },
		{ "steered acosh", 2, 0, 50, 15, "-0.000000001286C00000000000000000000000000000000000000000000000DCC67C24D66AD9AD51E066907350E8CC8E", "0.0000000002C62",
		  "2.EAE32ABB007D2332E5CBCE0A618C1260A44CC8E2B5C9A9B410",
		  "19.D65AEDCC582AD721EEE05A845A6726EA9290789C175220AA7A" },
		{ "steered acosh", 1, 0, 20, 16, "356046BE9.0000000000C413148992B4D2F2D7103BE47DCA368B80AD144E71F6DE28F12C", "-2492CD3EE6",
		  "1.7AD669A882DAEA38B963",
		  "1A.7A31C8B5A44CAA0D0A51" },
		{ "steered acosh", 4, 0, 10, 3, "-21221120101.0200221101000001120120211020210020020211", "-210200010",
		  "-110.1222100012",
		  "-0.0022011011" },
		{ "steered acosh", 5, 0, 20, 9, "7814726666081.0000000862783386316614876814376211037273477456207333", "18221420858210",
		  "33.00214684382200653413",
		  "1.12366038722161641514" },
		{ "decider", 2, 0, 32, 12, "-1.3", "-0.000000000000000000000000D937233D716",
		  "2.5E8851472408814315D5E9D829054D70",
		  "-0.000000000000000000000000E5E05766" },
		{ "decider", 6, 0, 13, 10, "0.28", "-0.00000000008170175999992",
		  "1.9459101490553",
		  "0.0000000003039" },
		{ "decider", 7, 0, 5, 10, "-5120", "3840",
		  "-0.00013",
		  "-0.00009" },
		{ "root pair", 6, 0, 20, 8, "-2.6077777777777777777363652206656131712304111346735200765616471752", "-2.52",
		  "0.13545700457317510750",
		  "1.60340372341750644240" },
		{ "near unit exact", 0, 1, 36, 4, "-1.000000000000000000000000000000000000000000000000000001", "-1.2",
		  "-132.000000000000000000000000000000000000",
		  "-1023.131030103203220003222233232220000120" },
		{ "near unit exact", 0, 2, 3, 14, "0.A6DDDDDDDDDDDDDDDDDD", "-0.37",
		  "38.000",
		  "-18.0C6" },
		{ "near unit exact", 0, 1, 20, 12, "0.94EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE", "0.41",
		  "39.00000000000000000000",
		  "22.732733514E9D25916934" },
	};

	typedef long (*ComplexUnary)( CLZComplex*, CLZComplex* );
	typedef long (*RealUnary)( CLZNumber*, CLZNumber* );

	static const ComplexUnary	kComplex[ 8 ] =
	{
		lzc_arcsine, lzc_arccosine, lzc_arcsecant, lzc_arccosecant,
		lzc_archypersine, lzc_archypercosine, lzc_archypersecant, lzc_archypercosecant,
	};
	static const RealUnary		kReal[ 8 ] =
	{
		lz_arcsine, lz_arccosine, lz_arcsecant, lz_arccosecant,
		lz_archypersine, lz_archypercosine, lz_archypersecant, lz_archypercosecant,
	};

	for ( const char* sz_family : { "random", "axes", "axis signs", "axis units", "cuts", "branch points", "large", "unit exact",
									"unit signs", "near unit exact", "near rational", "near segment", "steered", "steered acosh",
									"decider", "root pair" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZComplex	z;
			CLZComplex	t;

			glWorkingPrecision = pin.m_lPlaces;
			SetComplex( z, pin.m_lBase, pin.m_szReal, pin.m_szImag );
			geAngleUnit = static_cast<AngleUnit>( pin.m_lUnit );

			long	l_code = kComplex[ pin.m_lFunction ]( &z, &t );

			geAngleUnit = AU_radians;
			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &t.m_real, pin.m_szWantReal ) || !SpelledAs( &t.m_imag, pin.m_szWantImag ) ||
				 (t.m_real.m_lBase != pin.m_lBase) || (t.m_imag.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::to_string( pin.m_lFunction ) + " unit " + std::to_string( pin.m_lUnit ) + " base " +
							  std::to_string( pin.m_lBase ) + " " + pin.m_szReal + " " + pin.m_szImag;
			}
		}

		if ( l_cases == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "no pins" );
		else if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "inverses with square roots, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), (l_cases > 0) && (l_wrong == 0), szDetail );
	}

	// On an axis each part is a real function's, an exact nought, or -- a
	//	circular inverse's other part, in a unit -- its radians in the unit,
	//	which the axis units pins hold. Inside the real domain that is the
	//	real function's answer to the character, in the host's unit for the
	//	circular four; outside, as arcsin 2 is pi/2 + i acosh 2, arcsech 2
	//	-i arcsec 2, the +0 side of each cut. A real function refused as
	//	infinite -- arcsech and arccsch of nought -- is the point at
	//	infinity, with its code; arcsec and arccsc of nought are a pole...
	enum { kHost, kScaled, kNought, kPole };
	enum { kArgV, kArgMagnitude, kArgSign, kArgZero, kArgMinusOne };

	struct SAxisPart { int m_iKind; int m_iFunction; int m_iArgument; bool m_bRadians; bool m_bNegate; };

	auto	AxisParts = []( long lFunction, bool bReal, long lSign, long lAgainstOne, bool bUnit, SAxisPart* pParts )
	{
		const SAxisPart	N = { kNought, 0, 0, false, false };

		auto	H = []( int iF, int iArg, bool bNeg = false ) -> SAxisPart { return { kHost, iF, iArg, false, bNeg }; };
		auto	R = []( int iF, int iArg, bool bNeg = false ) -> SAxisPart { return { kHost, iF, iArg, true, bNeg }; };
		auto	V = [ bUnit ]( int iF, int iArg, bool bNeg ) -> SAxisPart { return { bUnit ? kScaled : kHost, iF, iArg, true, bNeg }; };

		bool	b_within = (lAgainstOne <= 0);		// |v| at most one
		bool	b_beyond = (lAgainstOne >= 0);		// |v| at least one

		switch ( lFunction )
		{
		case 0:
			if ( bReal ) { pParts[ 0 ] = H( 0, b_within ? kArgV : kArgSign ); pParts[ 1 ] = b_within ? N : V( 5, kArgMagnitude, false ); }
			else { pParts[ 0 ] = N; pParts[ 1 ] = V( 4, kArgV, false ); }
			break;

		case 1:
			if ( bReal ) { pParts[ 0 ] = H( 1, b_within ? kArgV : kArgSign ); pParts[ 1 ] = b_within ? N : V( 5, kArgMagnitude, true ); }
			else { pParts[ 0 ] = H( 1, kArgZero ); pParts[ 1 ] = V( 4, kArgV, true ); }
			break;

		case 2:
		case 3:
			if ( bReal && (lSign == 0) ) { pParts[ 0 ] = pParts[ 1 ] = { kPole, 0, 0, false, false }; }
			else if ( bReal ) { pParts[ 0 ] = H( lFunction, b_beyond ? kArgV : kArgSign ); pParts[ 1 ] = b_beyond ? N : V( 6, kArgMagnitude, lFunction == 3 ); }
			else { pParts[ 0 ] = (lFunction == 2) ? H( 1, kArgZero ) : N; pParts[ 1 ] = V( 7, kArgV, lFunction == 3 ); }
			break;

		case 4:
			if ( bReal ) { pParts[ 0 ] = R( 4, kArgV ); pParts[ 1 ] = N; }
			else if ( b_within ) { pParts[ 0 ] = N; pParts[ 1 ] = R( 0, kArgV ); }
			else { pParts[ 0 ] = R( 5, kArgMagnitude ); pParts[ 1 ] = R( 0, kArgSign ); }
			break;

		case 5:
			if ( !bReal ) { pParts[ 0 ] = R( 4, kArgMagnitude ); pParts[ 1 ] = R( 0, kArgSign ); }
			else if ( (lSign > 0) && b_beyond ) { pParts[ 0 ] = R( 5, kArgV ); pParts[ 1 ] = N; }
			else if ( b_within ) { pParts[ 0 ] = N; pParts[ 1 ] = R( 1, kArgV ); }
			else { pParts[ 0 ] = R( 5, kArgMagnitude ); pParts[ 1 ] = R( 1, kArgMinusOne ); }
			break;

		case 6:
			if ( !bReal ) { pParts[ 0 ] = R( 7, kArgMagnitude ); pParts[ 1 ] = R( 0, kArgSign, true ); }
			else if ( (lSign >= 0) && b_within ) { pParts[ 0 ] = R( 6, kArgV ); pParts[ 1 ] = N; }
			else if ( b_within ) { pParts[ 0 ] = R( 6, kArgMagnitude ); pParts[ 1 ] = R( 1, kArgMinusOne, true ); }
			else { pParts[ 0 ] = N; pParts[ 1 ] = R( 2, kArgV, true ); }
			break;

		default:
			if ( bReal ) { pParts[ 0 ] = R( 7, kArgV ); pParts[ 1 ] = N; }
			else if ( b_beyond ) { pParts[ 0 ] = N; pParts[ 1 ] = R( 3, kArgV, true ); }
			else { pParts[ 0 ] = R( 6, kArgMagnitude ); pParts[ 1 ] = R( 0, kArgSign, true ); }
			break;
		}
	};

	{
		std::uint64_t	u_state = 6155;
		long			l_cases = 0;
		long			l_wrong = 0;
		std::string		s_first;

		for ( long i = 0; i < 51; i++ )
		{
			long	l_base = 2 + (i % 15);
			long	l_unit = ((i % 4) == 3) ? 1 + ((i / 4) % 3) : 0;

			glWorkingPrecision = (i % 3 == 0) ? 20 : 36;

			CLZNumber	v;

			v.SetBase( l_base, false );

			// ...nought and one of each sign, then values either side of one...
			if ( i < 3 )
				v.Assign( (i == 0) ? 0L : ((i == 1) ? 1L : -1L) );
			else
			{
				long	l_places = 1 + static_cast<long>( (u_state >> 40) % 4 );
				long	l_whole = static_cast<long>( (u_state >> 45) % 3 );

				v.Assign( RandomDigitsIn( u_state, l_base, l_places + l_whole, l_places, (i % 2) == 1 ).c_str() );
			}

			CLZNumber	nought;
			CLZNumber	one;
			CLZNumber	magnitude;
			CLZNumber	sign;
			CLZNumber	minus_one;

			for ( CLZNumber* p : { &nought, &one, &magnitude, &sign, &minus_one } )
				p->SetBase( l_base, false );

			nought.Assign( 0L );
			one.Assign( 1L );
			minus_one.Assign( -1L );
			v.CheckForZero();
			magnitude.Assign( &v );

			if ( v.IsNegative() )
				lz_negate( &magnitude );

			sign.Assign( v.IsNegative() ? -1L : 1L );

			long	l_sign = v.IsZero() ? 0 : (v.IsNegative() ? -1 : 1);
			long	l_order = lz_compare( &magnitude, &one );
			long	l_against = (l_order == kGreaterThan) ? 1 : ((l_order == kLessThan) ? -1 : 0);

			CLZNumber*	kArguments[] = { &v, &magnitude, &sign, &nought, &minus_one };

			for ( long l_fn = 0; l_fn < 8; l_fn++ )
			{
				for ( int i_axis = 0; i_axis < 2; i_axis++ )
				{
					bool	b_real = (i_axis == 0);

					if ( !b_real && (l_sign == 0) )
						continue;

					SAxisPart	parts[ 2 ];
					CLZComplex	z;
					CLZComplex	t;

					AxisParts( l_fn, b_real, l_sign, l_against, l_unit != 0, parts );
					z.m_real.Assign( b_real ? &v : &nought );
					z.m_imag.Assign( b_real ? &nought : &v );
					geAngleUnit = static_cast<AngleUnit>( l_unit );

					long	l_got = kComplex[ l_fn ]( &z, &t );
					bool	b_ok = true;

					if ( parts[ 0 ].m_iKind == kPole )
						b_ok = (l_got == kUndefinedResult) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
					else
					{
						CLZNumber	want[ 2 ];
						long		l_want = kNoError;
						bool		b_infinite = false;

						for ( int k = 0; k < 2; k++ )
						{
							if ( parts[ k ].m_iKind != kHost )
								continue;

							geAngleUnit = parts[ k ].m_bRadians ? AU_radians : static_cast<AngleUnit>( l_unit );

							long	l_part = kReal[ parts[ k ].m_iFunction ]( kArguments[ parts[ k ].m_iArgument ], &want[ k ] );

							if ( parts[ k ].m_bNegate && !want[ k ].IsZero() )
								lz_negate( &want[ k ] );

							if ( l_want == kNoError )
								l_want = l_part;

							b_infinite = b_infinite || want[ k ].FlagIsOn( kfInfinity );
						}

						if ( b_infinite )
							b_ok = (l_got == l_want) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
						else
						{
							b_ok = (l_got == l_want);

							for ( int k = 0; k < 2; k++ )
							{
								CLZNumber*	p_got = (k == 0) ? &t.m_real : &t.m_imag;

								if ( parts[ k ].m_iKind == kNought )
									b_ok = b_ok && SpelledAs( p_got, "0" );
								else if ( parts[ k ].m_iKind == kHost )
									b_ok = b_ok && SameAnswer( p_got, &want[ k ] ) && (p_got->m_lPrecision == want[ k ].m_lPrecision);
							}
						}
					}

					geAngleUnit = AU_radians;
					l_cases++;

					if ( !b_ok && (l_wrong++ == 0) )
					{
						char*	sz_v = v.GetStringRep();

						s_first = std::to_string( l_fn ) + (b_real ? " real " : " imaginary ") + "base " + std::to_string( l_base ) +
								  " unit " + std::to_string( l_unit ) + " " + ((sz_v != nullptr) ? sz_v : "");
						::free( sz_v );
					}
				}
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, first %s", l_wrong, l_cases, s_first.c_str() );

		CheckTrue( "on an axis, the real functions'", l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 20;

	// The circular four write the whole angle in the unit, exact where it
	//	is -- arcsin 2 is 90 + 75.456...i degrees, arcsin(1 + 1.5i) 30 and a
	//	little, arccos(1 + 1.5i) 60 -- and the hyperbolic four answer in
	//	radians whatever it is...
	{
		CLZComplex	z;
		CLZComplex	t;
		bool		b_ok = true;

		geAngleUnit = AU_degrees;
		SetComplex( z, 10, "2", "0" );
		b_ok = b_ok && (lzc_arcsine( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "90" );
		SetComplex( z, 10, "1", "1.5" );
		b_ok = b_ok && (lzc_arcsine( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "30" );
		b_ok = b_ok && (lzc_arccosine( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "60" );
		SetComplex( z, 10, "-1", "-1.5" );
		b_ok = b_ok && (lzc_arccosine( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "120" );
		SetComplex( z, 10, "0.75", "0.25" );
		b_ok = b_ok && (lzc_arcsine( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "45" );
		geAngleUnit = AU_radians;

		CheckTrue( "the circular four in the unit", b_ok, b_ok ? "90, 30, 60, 120 and 45 degrees, exact" : "an answer wrong" );

		long	l_wrong = 0;

		for ( long l_fn = 4; l_fn < 8; l_fn++ )
		{
			CLZComplex	radians;
			CLZComplex	degrees;

			SetComplex( z, 10, "0.75", "-1.25" );
			kComplex[ l_fn ]( &z, &radians );
			geAngleUnit = AU_degrees;
			kComplex[ l_fn ]( &z, &degrees );
			geAngleUnit = AU_radians;
			l_wrong += (SameAnswer( &radians.m_real, &degrees.m_real ) && SameAnswer( &radians.m_imag, &degrees.m_imag )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 4 differ", l_wrong );
		CheckTrue( "the hyperbolic four in radians", l_wrong == 0, szDetail );
	}

	// A not-a-number or the point at infinity in either part answers a
	//	not-a-number, and an operand in two bases is refused...
	{
		CLZComplex	z;
		CLZComplex	t;
		CLZNumber	nan;
		CLZNumber	infinite;
		bool		b_ok = true;

		MakeNotANumber( nan );
		MakeInfinite( infinite );

		for ( ComplexUnary fn : kComplex )
		{
			for ( long l_where = 0; l_where < 4; l_where++ )
			{
				SetComplex( z, 10, "1", "2" );
				(((l_where / 2) == 0) ? z.m_real : z.m_imag).Assign( ((l_where % 2) == 0) ? &nan : &infinite );
				b_ok = b_ok && (fn( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );
			}

			SetComplex( z, 10, "1", "2" );
			z.m_imag.SetBase( 16, false );
			z.m_imag.Assign( "2" );
			b_ok = b_ok && (fn( &z, &t ) == kBaseMismatch) && t.m_real.FlagIsOn( kfQNaN );
		}

		CheckTrue( "flags and two bases", b_ok, b_ok ? "carried or refused as the layer does" : "an answer wrong" );
	}

	// In place, as from a copy...
	{
		long	l_wrong = 0;

		for ( long l_unit = 0; l_unit < 2; l_unit++ )
		{
			for ( ComplexUnary fn : kComplex )
			{
				CLZComplex	z;
				CLZComplex	copy;
				CLZComplex	t;

				geAngleUnit = static_cast<AngleUnit>( l_unit );
				SetComplex( z, 7, "-3.25", "1.5" );
				SetComplex( copy, 7, "-3.25", "1.5" );
				fn( &copy, &t );
				fn( &z, nullptr );
				l_wrong += (SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag )) ? 0 : 1;
			}
		}

		geAngleUnit = AU_radians;
		::snprintf( szDetail, sizeof( szDetail ), "%ld of 16 differ", l_wrong );
		CheckTrue( "inverses with square roots in place", l_wrong == 0, szDetail );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**The inverse functions that take logarithms, stage 6b** -- 2026-10-23,
	on Rich's calls: lzc_arctangent, lzc_arccotangent, lzc_archypertangent
	and lzc_archypercotangent, arccot being arctan(1/z) as the real
	lz_arccotangent became arctan(1/x) that day. Every expected part is
	generated -- Python's decimal from the C99 logarithm forms, deep enough
	for the case's own places, and exact Fractions where a part in a unit
	is rational -- never typed. Each family past the random is made of
	cases a mutant gets wrong: the axes by the builds that negate the
	quarter turn of arctan past i, leave arccot's cut and arccoth's
	imaginary axis unsigned, negate arctanh's part past one, give arccoth
	+(pi/2)i inside one and take the real axis the general way; the cuts by
	the builds that leave arccoth's angle unsigned and arccot's logarithm
	the wrong way up; the branch points and the large by the build without
	a guard; the unit exact by the builds that never find an exact angle,
	take three sixteenths of a turn as one and leave it unsigned; the near
	unit exact by the builds that test the exact angle in doubles, switch
	the unit decider off, reverse it, take its cosines for its sines and
	leave its half unsigned; the near rational by the builds without the
	tie decider, with its angle verdict reversed and with its angle against
	h rather than 2h; the near segment by the builds whose logarithm treats
	every h as past nought, compares against 2h and leaves the half
	unsigned; the steered by the builds without a retake and with the angle
	box's wrong sine end; and the steered log by the builds that bound the
	exponential's series from N above |q| and do not double its remainder.
*/

static void TestComplexLogInverses()
{
	Banner( "complex numbers, stage 6b: the inverse functions that take logarithms" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	struct { const char* m_szFamily; long m_lFunction; long m_lUnit; long m_lPlaces; long m_lBase; const char* m_szReal;
			 const char* m_szImag; const char* m_szWantReal; const char* m_szWantImag; } kPins[] =
	{
		{ "random", 2, 0, 1, 9, "0.3316", "-0.3508",
		  "0.3",
		  "-0.4" },
		{ "random", 1, 0, 10, 8, "-60770", "-37.2",
		  "-0.0000123475",
		  "0.0000000065" },
		{ "random", 1, 2, 10, 13, "8.C39", "-BB8.B",
		  "0.0004039A21",
		  "0.0547284AC3" },
		{ "random", 3, 0, 1, 3, "-2.121", "-0.12",
		  "-0.1",
		  "0.0" },
		{ "random", 2, 0, 36, 10, "-484.7", "-0.7331",
		  "-0.002063130041727657237170782798524321",
		  "-1.570793206339098134579238210315959642" },
		{ "random", 0, 1, 50, 7, "0.5465", "2.26",
		  "144.45552620405551221534501251416025644332514221131322",
		  "30.65560310146146205525540514554602064003201233161240" },
		{ "random", 2, 0, 10, 2, "-101.1", "-11.1",
		  "-0.0010000100",
		  "-1.0111110011" },
		{ "random", 3, 0, 20, 8, "5.417", "7551",
		  "0.00000005753671206505",
		  "-0.00010234612711216613" },
		{ "random", 2, 0, 36, 6, "-1321000000", "340400000000",
		  "-0.000000000000004045113230240015100320",
		  "1.323143023503554244405122445034430123" },
		{ "random", 3, 0, 3, 3, "-0.000000001", "-0.0000001002",
		  "0.000",
		  "1.120" },
		{ "random", 1, 2, 36, 8, "-426.2", "-162",
		  "-0.144234312525566063147650623523034302",
		  "0.051060572765246632721667264030173301" },
		{ "random", 3, 0, 10, 15, "C63.3", "0.892",
		  "0.00131CEB62",
		  "-0.000000C87E" },
		{ "random", 0, 3, 50, 10, "701400000000", "-44900000000",
		  "249.99999999977401567708645916232196549669799573448775",
		  "-0.00000000001446634744627599602472731038551784298999" },
		{ "random", 2, 0, 0, 8, "-0.0000004172", "0.00000465",
		  "0",
		  "0" },
		{ "random", 2, 0, 1, 8, "-0.0000007011", "-6.313",
		  "0.0",
		  "-1.3" },
		{ "random", 2, 0, 36, 10, "-0.000000007299", "5.222",
		  "-0.000000000258195432187104561505460041",
		  "1.381589444110258321570629128820664815" },
		{ "axes", 0, 1, 36, 11, "0", "-2",
		  "82",
		  "-29.52250926844682560386452A4A6146607050" },
		{ "axes", 1, 1, 10, 9, "0", "0.1",
		  "110",
		  "-6.3471713615" },
		{ "axes", 2, 0, 10, 10, "2", "0",
		  "0.5493061443",
		  "1.5707963268" },
		{ "axes", 3, 0, 0, 4, "0.2", "0",
		  "1",
		  "-2" },
		{ "axes", 3, 0, 0, 8, "0", "-1.4",
		  "0",
		  "1" },
		{ "axes", 0, 3, 50, 14, "-3B200", "0",
		  "-13B.DDB0101B2DA46634C904D614B8B9B66A6BA10B4591D743D01D",
		  "0" },
		{ "cuts", 3, 0, 36, 6, "0.1", "-0.00000001",
		  "0.100201120242340545204043415010323502",
		  "1.323143013404331531405104235155304122" },
		{ "cuts", 3, 0, 36, 10, "0.0001", "2",
		  "0.000019999999970666666692906666671970",
		  "-0.463647608200806117174256230991881068" },
		{ "cuts", 1, 0, 0, 8, "0.00000000000000000000000001", "2",
		  "0",
		  "-1" },
		{ "cuts", 1, 0, 50, 6, "2", "0.000000000000000000000001",
		  "0.24405153535322242431530451201224500213013302031155",
		  "-0.00000000000000000000000011111111111111111111111111" },
		{ "branch points", 2, 0, 3, 9, "0.000000000000001", "0.8888888888888888888",
		  "0.000",
		  "0.706" },
		{ "branch points", 1, 0, 20, 5, "0.0000000000000000000000000001", "0.44444444444444444",
		  "1.24113332123421343040",
		  "-24.00313332030424234444" },
		{ "large", 3, 0, 36, 9, "-7440000000", "-3700000",
		  "-0.000000000117246117700471850823073511",
		  "0.000000000000054032574032074513058222" },
		{ "large", 3, 0, 36, 7, "-34300000000000000000000", "501000000",
		  "-0.000000000000000000000016326434413121",
		  "0.000000000000000000000000000000000000" },
		{ "unit exact", 1, 1, 10, 16, "-2", "1",
		  "-16.8",
		  "-9.EDB8DC7281" },
		{ "unit exact", 0, 1, 36, 15, "2", "1",
		  "47.777777777777777777777777777777777778",
		  "9.DDE07237B654757129214711241BDD76DE6A" },
		{ "unit exact", 0, 3, 0, 16, "-2", "1",
		  "-BC",
		  "1C" },
		{ "unit exact", 0, 2, 0, 13, "-0.CC", "1.55",
		  "-5A",
		  "22" },
		{ "unit exact", 1, 3, 36, 8, "-2", "1",
		  "-76.4",
		  "-33.450533261031742405145161023035177754" },
		{ "unit exact", 1, 1, 36, 16, "-2", "-1",
		  "-16.8",
		  "9.EDB8DC7280F2BFCDFFBBC0128EE86EF8FCBC" },
		{ "near unit exact", 1, 3, 1, 5, "-0.44430444444444444444444444444", "-0.01422",
		  "-1000.0",
		  "11.0" },
		{ "near unit exact", 0, 3, 0, 14, "-1.DDDDDDDDDDDDDDDDDDDD", "-1",
		  "-D5",
		  "-20" },
		{ "near unit exact", 1, 1, 50, 13, "0.4000000000000000000000000000000000000000000000000000000000000000001", "0.7",
		  "52.66666666666666666666666666666666666666666666666666",
		  "-23.A2A553759AB55788A78837C664CC541323BB98C4969C1310A0" },
		{ "near unit exact", 1, 3, 1, 9, "-2.000000000001", "-1",
		  "-68.4",
		  "30.5" },
		{ "near unit exact", 1, 1, 1, 15, "-2.38EEEEEEEEEEEEEEEEEEEEEE", "-0.A3",
		  "-17.8",
		  "6.2" },
		{ "near rational", 3, 0, 25, 10, "-384000000000000000000000", "512000000000000000000000",
		  "-0.0000000000000000000000009",
		  "-0.0000000000000000000000013" },
		{ "near rational", 1, 0, 22, 10, "-2560000000000000000000", "-1920000000000000000000",
		  "-0.0000000000000000000003",
		  "0.0000000000000000000002" },
		{ "near rational", 3, 0, 30, 10, "-192000000000000000000000000000", "-256000000000000000000000000000",
		  "-0.000000000000000000000000000002",
		  "0.000000000000000000000000000003" },
		{ "near rational", 3, 0, 23, 10, "19200000000000000000000", "25600000000000000000000",
		  "0.00000000000000000000002",
		  "-0.00000000000000000000003" },
		{ "near rational", 2, 0, 36, 8, "0.000000000000000000000000000000000000052", "-0.0000000000000000000000000000000000754",
		  "0.000000000000000000000000000000000000",
		  "-0.000000000000000000000000000000000075" },
		{ "near rational", 0, 0, 23, 8, "0.000000000000000000000544", "0.000000000000000000000000000034",
		  "0.00000000000000000000054",
		  "0.00000000000000000000000" },
		{ "near segment", 2, 0, 29, 6, "-0.000000000000000000011021311352442", "0.34",
		  "-0.00000000000000000000505042100",
		  "0.31425304204140221304153313033" },
		{ "near segment", 2, 0, 30, 6, "-0.00000000000000000001351204541125555555555555555555555555555553", "12",
		  "-0.000000000000000000000052445131",
		  "1.240233100010033102304145345454" },
		{ "near segment", 3, 0, 33, 14, "0.00000000000000000000A86D7C39380107", "-2C",
		  "0.000000000000000000000014272C0C8D8",
		  "0.04C82B2AC46194A832BA532A7BC8BDC97" },
		{ "near segment", 1, 0, 3, 14, "1.5", "-0.0015C7",
		  "0.8C7",
		  "0.000" },
		{ "near segment", 1, 0, 7, 6, "-0.00000253", "-10",
		  "0.0000000",
		  "0.1002011" },
		{ "near segment", 3, 0, 6, 12, "-44", "0.0070576",
		  "-0.02929D",
		  "-0.000004" },
		{ "steered", 1, 2, 36, 9, "221.00000000000000000000000000000000430573112337073857167238647708622337466144240375", "-4.72",
		  "0.314200861624365483446406702532780160",
		  "0.006715820130258573327883467608651840" },
		{ "steered", 0, 1, 10, 12, "-D.63EEEEE957916ED828599D1421427DED96552584E3790349", "-94.7",
		  "-75.E5257E5236",
		  "-0.60798869DE" },
		{ "steered", 2, 0, 50, 16, "C44.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBB651AC69CBECB000005A9966DF6E3017A54460D23AEDF70056A1FB", "-5.67",
		  "0.0014DD5881BB431D50AA7582EEDF01116D9BEDBD503DD319BB",
		  "-1.921FAC146FAC8B57E17F8F1D8E628CA59478CA0272BDA2AF32" },
		{ "steered", 0, 0, 10, 15, "-3.5600000026987729CA2C0C0970B3C293E553D58707E811", "-C.76",
		  "-1.83D4AB7D31",
		  "-0.11C4511367" },
		{ "steered log", 0, 0, 20, 8, "-0.000000000000102773000000000000001205143665505650116134727372064", "0.777777777777734035",
		  "-0.42403423257054763157",
		  "15.61560155412627075644" },
		{ "steered log", 0, 0, 50, 10, "0.0000000000002893399999999999999999999999999999999999999999999923348538378524340370159228", "-0.99999999999946256",
		  "0.24693406702337089069571745441543086914058221435337",
		  "-14.40892385318831928304409710682392211127180303034605" },
		{ "steered log", 1, 0, 36, 16, "-0.000000000000000000004150FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED7A4AED68D2C6CFD7921C5EF927E60EBF", "1.000000000000000000009D51",
		  "-0.325F4D19418DB0697686892D59CB957B458E",
		  "-1C.46B26F4C478A79879C965C1F9A0AC9B2E885" },
		{ "steered log", 3, 0, 50, 10, "-1.00000000104699999999999999999999999999999999999999999999966735671750837989227166651428", "-0.000000011355",
		  "-9.49126091562336276509527107898292961318366109902572",
		  "0.73942511517881160609078112118999445083543409180769" },
	};

	typedef long (*ComplexUnary)( CLZComplex*, CLZComplex* );
	typedef long (*RealUnary)( CLZNumber*, CLZNumber* );

	static const ComplexUnary	kComplex[ 4 ] =
	{
		lzc_arctangent, lzc_arccotangent, lzc_archypertangent, lzc_archypercotangent,
	};
	static const RealUnary		kReal[ 4 ] =
	{
		lz_arctangent, lz_arccotangent, lz_archypertangent, lz_archypercotangent,
	};

	for ( const char* sz_family : { "random", "axes", "cuts", "branch points", "large", "unit exact",
									"near unit exact", "near rational", "near segment", "steered", "steered log" } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( ::strcmp( pin.m_szFamily, sz_family ) != 0 )
				continue;

			CLZComplex	z;
			CLZComplex	t;

			glWorkingPrecision = pin.m_lPlaces;
			SetComplex( z, pin.m_lBase, pin.m_szReal, pin.m_szImag );
			geAngleUnit = static_cast<AngleUnit>( pin.m_lUnit );

			long	l_code = kComplex[ pin.m_lFunction ]( &z, &t );

			geAngleUnit = AU_radians;
			l_cases++;

			if ( (l_code != kNoError) || !SpelledAs( &t.m_real, pin.m_szWantReal ) || !SpelledAs( &t.m_imag, pin.m_szWantImag ) ||
				 (t.m_real.m_lBase != pin.m_lBase) || (t.m_imag.m_lBase != pin.m_lBase) )
			{
				if ( l_wrong++ == 0 )
					s_first = std::to_string( pin.m_lFunction ) + " unit " + std::to_string( pin.m_lUnit ) + " base " +
							  std::to_string( pin.m_lBase ) + " " + pin.m_szReal + " " + pin.m_szImag;
			}
		}

		if ( l_cases == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "no pins" );
		else if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		std::string	s_label = std::string( "inverses with logarithms, " ) + sz_family + " pins";

		CheckTrue( s_label.c_str(), (l_cases > 0) && (l_wrong == 0), szDetail );
	}

	// On an axis each part is a real function's, an exact nought, or -- a
	//	circular inverse's other part, in a unit -- its radians in the unit,
	//	which the axis pins hold. Inside the real domain that is the real
	//	function's answer to the character, in the host's unit for the
	//	circular two; outside, the +0 side of each cut: arctanh 2 is arccoth 2
	//	+ (pi/2)i, arctan 2i a quarter turn + (arccoth 2)i, the quarter turn
	//	taken here as the real arccot of nought. A real function refused as
	//	infinite -- arctanh and arccoth of one -- is the point at infinity,
	//	with its code; arctan and arccot of i are a pole...
	enum { kHost, kScaled, kNought, kPole };
	enum { kArgV, kArgMagnitude, kArgSign, kArgZero, kArgMinusOne };

	struct SAxisPart { int m_iKind; int m_iFunction; int m_iArgument; bool m_bRadians; bool m_bNegate; };

	auto	AxisParts = []( long lFunction, bool bReal, long lSign, long lAgainstOne, bool bUnit, SAxisPart* pParts )
	{
		const SAxisPart	N = { kNought, 0, 0, false, false };

		auto	H = []( int iF, int iArg, bool bNeg = false ) -> SAxisPart { return { kHost, iF, iArg, false, bNeg }; };
		auto	R = []( int iF, int iArg, bool bNeg = false ) -> SAxisPart { return { kHost, iF, iArg, true, bNeg }; };
		auto	V = [ bUnit ]( int iF, int iArg, bool bNeg ) -> SAxisPart { return { bUnit ? kScaled : kHost, iF, iArg, true, bNeg }; };

		bool	b_within = (lAgainstOne <= 0);		// |v| at most one
		bool	b_beyond = (lAgainstOne >= 0);		// |v| at least one

		(void)lSign;

		switch ( lFunction )
		{
		case 0:
			if ( bReal ) { pParts[ 0 ] = H( 0, kArgV ); pParts[ 1 ] = N; }
			else if ( lAgainstOne < 0 ) { pParts[ 0 ] = N; pParts[ 1 ] = V( 2, kArgV, false ); }
			else if ( lAgainstOne == 0 ) { pParts[ 0 ] = pParts[ 1 ] = { kPole, 0, 0, false, false }; }
			else { pParts[ 0 ] = H( 1, kArgZero ); pParts[ 1 ] = V( 3, kArgV, false ); }
			break;

		case 1:
			if ( bReal ) { pParts[ 0 ] = H( 1, kArgV ); pParts[ 1 ] = N; }
			else if ( lAgainstOne > 0 ) { pParts[ 0 ] = N; pParts[ 1 ] = V( 3, kArgV, true ); }
			else if ( lAgainstOne == 0 ) { pParts[ 0 ] = pParts[ 1 ] = { kPole, 0, 0, false, false }; }
			else { pParts[ 0 ] = H( 1, kArgZero ); pParts[ 1 ] = V( 2, kArgV, true ); }
			break;

		case 2:
			if ( !bReal ) { pParts[ 0 ] = N; pParts[ 1 ] = R( 0, kArgV ); }
			else if ( b_within ) { pParts[ 0 ] = R( 2, kArgV ); pParts[ 1 ] = N; }
			else { pParts[ 0 ] = R( 3, kArgV ); pParts[ 1 ] = R( 1, kArgZero ); }
			break;

		default:
			if ( !bReal ) { pParts[ 0 ] = N; pParts[ 1 ] = R( 1, kArgV, true ); }
			else if ( b_beyond ) { pParts[ 0 ] = R( 3, kArgV ); pParts[ 1 ] = N; }
			else { pParts[ 0 ] = R( 2, kArgV ); pParts[ 1 ] = R( 1, kArgZero, true ); }
			break;
		}
	};

	{
		std::uint64_t	u_state = 6255;
		long			l_cases = 0;
		long			l_wrong = 0;
		std::string		s_first;

		for ( long i = 0; i < 51; i++ )
		{
			long	l_base = 2 + (i % 15);
			long	l_unit = ((i % 4) == 3) ? 1 + ((i / 4) % 3) : 0;

			glWorkingPrecision = (i % 3 == 0) ? 20 : 36;

			CLZNumber	v;

			v.SetBase( l_base, false );

			// ...nought and one of each sign, then values either side of one...
			if ( i < 3 )
				v.Assign( (i == 0) ? 0L : ((i == 1) ? 1L : -1L) );
			else
			{
				long	l_places = 1 + static_cast<long>( (u_state >> 40) % 4 );
				long	l_whole = static_cast<long>( (u_state >> 45) % 3 );

				v.Assign( RandomDigitsIn( u_state, l_base, l_places + l_whole, l_places, (i % 2) == 1 ).c_str() );
			}

			CLZNumber	nought;
			CLZNumber	one;
			CLZNumber	magnitude;
			CLZNumber	sign;
			CLZNumber	minus_one;

			for ( CLZNumber* p : { &nought, &one, &magnitude, &sign, &minus_one } )
				p->SetBase( l_base, false );

			nought.Assign( 0L );
			one.Assign( 1L );
			minus_one.Assign( -1L );
			v.CheckForZero();
			magnitude.Assign( &v );

			if ( v.IsNegative() )
				lz_negate( &magnitude );

			sign.Assign( v.IsNegative() ? -1L : 1L );

			long	l_sign = v.IsZero() ? 0 : (v.IsNegative() ? -1 : 1);
			long	l_order = lz_compare( &magnitude, &one );
			long	l_against = (l_order == kGreaterThan) ? 1 : ((l_order == kLessThan) ? -1 : 0);

			CLZNumber*	kArguments[] = { &v, &magnitude, &sign, &nought, &minus_one };

			for ( long l_fn = 0; l_fn < 4; l_fn++ )
			{
				for ( int i_axis = 0; i_axis < 2; i_axis++ )
				{
					bool	b_real = (i_axis == 0);

					if ( !b_real && (l_sign == 0) )
						continue;

					SAxisPart	parts[ 2 ];
					CLZComplex	z;
					CLZComplex	t;

					AxisParts( l_fn, b_real, l_sign, l_against, l_unit != 0, parts );
					z.m_real.Assign( b_real ? &v : &nought );
					z.m_imag.Assign( b_real ? &nought : &v );
					geAngleUnit = static_cast<AngleUnit>( l_unit );

					long	l_got = kComplex[ l_fn ]( &z, &t );
					bool	b_ok = true;

					if ( parts[ 0 ].m_iKind == kPole )
						b_ok = (l_got == kUndefinedResult) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
					else
					{
						CLZNumber	want[ 2 ];
						long		l_want = kNoError;
						bool		b_infinite = false;

						for ( int k = 0; k < 2; k++ )
						{
							if ( parts[ k ].m_iKind != kHost )
								continue;

							geAngleUnit = parts[ k ].m_bRadians ? AU_radians : static_cast<AngleUnit>( l_unit );

							long	l_part = kReal[ parts[ k ].m_iFunction ]( kArguments[ parts[ k ].m_iArgument ], &want[ k ] );

							if ( parts[ k ].m_bNegate && !want[ k ].IsZero() )
								lz_negate( &want[ k ] );

							if ( l_want == kNoError )
								l_want = l_part;

							b_infinite = b_infinite || want[ k ].FlagIsOn( kfInfinity );
						}

						if ( b_infinite )
							b_ok = (l_got == l_want) && t.m_real.FlagIsOn( kfInfinity ) && t.m_imag.FlagIsOn( kfInfinity );
						else
						{
							b_ok = (l_got == l_want);

							for ( int k = 0; k < 2; k++ )
							{
								CLZNumber*	p_got = (k == 0) ? &t.m_real : &t.m_imag;

								if ( parts[ k ].m_iKind == kNought )
									b_ok = b_ok && SpelledAs( p_got, "0" );
								else if ( parts[ k ].m_iKind == kHost )
									b_ok = b_ok && SameAnswer( p_got, &want[ k ] ) && (p_got->m_lPrecision == want[ k ].m_lPrecision);
							}
						}
					}

					geAngleUnit = AU_radians;
					l_cases++;

					if ( !b_ok && (l_wrong++ == 0) )
					{
						char*	sz_v = v.GetStringRep();

						s_first = std::to_string( l_fn ) + (b_real ? " real " : " imaginary ") + "base " + std::to_string( l_base ) +
								  " unit " + std::to_string( l_unit ) + " " + ((sz_v != nullptr) ? sz_v : "");
						::free( sz_v );
					}
				}
			}
		}

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, first %s", l_wrong, l_cases, s_first.c_str() );

		CheckTrue( "on an axis, the real functions'", l_wrong == 0, szDetail );
	}

	glWorkingPrecision = 20;

	// The circular two write the whole angle in the unit, exact where its
	//	direction is a multiple of 45 degrees -- arctan(0.6 + 0.8i) is 45 and
	//	a little i, arctan(0.4 + 0.2i) 22.5, arctan(-2 - i) -67.5, arccot(0.4 +
	//	0.2i) 67.5 -- and the hyperbolic two answer in radians whatever it is...
	{
		CLZComplex	z;
		CLZComplex	t;
		bool		b_ok = true;

		geAngleUnit = AU_degrees;
		SetComplex( z, 10, "0.6", "0.8" );
		b_ok = b_ok && (lzc_arctangent( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "45" );
		b_ok = b_ok && (lzc_arccotangent( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "45" );
		SetComplex( z, 10, "0.4", "0.2" );
		b_ok = b_ok && (lzc_arctangent( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "22.5" );
		b_ok = b_ok && (lzc_arccotangent( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "67.5" );
		SetComplex( z, 10, "-2", "-1" );
		b_ok = b_ok && (lzc_arctangent( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "-67.5" );
		geAngleUnit = AU_radians;

		CheckTrue( "the circular two in the unit", b_ok, b_ok ? "45, 22.5, 67.5 and -67.5 degrees, exact" : "an answer wrong" );

		long	l_wrong = 0;

		for ( long l_fn = 2; l_fn < 4; l_fn++ )
		{
			CLZComplex	radians;
			CLZComplex	degrees;

			SetComplex( z, 10, "0.75", "-1.25" );
			kComplex[ l_fn ]( &z, &radians );
			geAngleUnit = AU_degrees;
			kComplex[ l_fn ]( &z, &degrees );
			geAngleUnit = AU_radians;
			l_wrong += (SameAnswer( &radians.m_real, &degrees.m_real ) && SameAnswer( &radians.m_imag, &degrees.m_imag )) ? 0 : 1;
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld of 2 differ", l_wrong );
		CheckTrue( "the hyperbolic two in radians", l_wrong == 0, szDetail );
	}

	// The poles, the point at infinity and kUndefinedResult as
	//	lz_archypertangent( 1 ) answers -- and arccoth 0, the cut's -(pi/2)i
	//	from either side, an exact nought beside it...
	{
		CLZComplex	z;
		CLZComplex	t;
		CLZNumber	one;
		CLZNumber	half_pi;
		bool		b_ok = true;

		struct { long m_lFunction; const char* m_szReal; const char* m_szImag; } kPoles[] =
		{
			{ 0, "0", "1" }, { 0, "0", "-1" }, { 1, "0", "1" }, { 1, "0", "-1" },
			{ 2, "1", "0" }, { 2, "-1", "0" }, { 3, "1", "0" }, { 3, "-1", "0" },
		};

		for ( auto& pole : kPoles )
		{
			SetComplex( z, 10, pole.m_szReal, pole.m_szImag );
			b_ok = b_ok && (kComplex[ pole.m_lFunction ]( &z, &t ) == kUndefinedResult) && t.m_real.FlagIsOn( kfInfinity ) &&
				   t.m_imag.FlagIsOn( kfInfinity );
		}

		one.SetBase( 10, false );
		one.Assign( 1L );
		lz_arcsine( &one, &half_pi );
		lz_negate( &half_pi );
		SetComplex( z, 10, "0", "0" );
		b_ok = b_ok && (lzc_archypercotangent( &z, &t ) == kNoError) && SpelledAs( &t.m_real, "0" ) && SameAnswer( &t.m_imag, &half_pi );

		CheckTrue( "the poles, and arccoth 0", b_ok, b_ok ? "eight points at infinity; arccoth 0 = -(pi/2)i" : "an answer wrong" );
	}

	// A not-a-number or the point at infinity in either part answers a
	//	not-a-number, and an operand in two bases is refused...
	{
		CLZComplex	z;
		CLZComplex	t;
		CLZNumber	nan;
		CLZNumber	infinite;
		bool		b_ok = true;

		MakeNotANumber( nan );
		MakeInfinite( infinite );

		for ( ComplexUnary fn : kComplex )
		{
			for ( long l_where = 0; l_where < 4; l_where++ )
			{
				SetComplex( z, 10, "1", "2" );
				(((l_where / 2) == 0) ? z.m_real : z.m_imag).Assign( ((l_where % 2) == 0) ? &nan : &infinite );
				b_ok = b_ok && (fn( &z, &t ) == kNoError) && t.m_real.FlagIsOn( kfQNaN ) && t.m_imag.FlagIsOn( kfQNaN );
			}

			SetComplex( z, 10, "1", "2" );
			z.m_imag.SetBase( 16, false );
			z.m_imag.Assign( "2" );
			b_ok = b_ok && (fn( &z, &t ) == kBaseMismatch) && t.m_real.FlagIsOn( kfQNaN );
		}

		CheckTrue( "flags and two bases", b_ok, b_ok ? "carried or refused as the layer does" : "an answer wrong" );
	}

	// In place, as from a copy...
	{
		long	l_wrong = 0;

		for ( long l_unit = 0; l_unit < 2; l_unit++ )
		{
			for ( ComplexUnary fn : kComplex )
			{
				CLZComplex	z;
				CLZComplex	copy;
				CLZComplex	t;

				geAngleUnit = static_cast<AngleUnit>( l_unit );
				SetComplex( z, 7, "-3.25", "1.5" );
				SetComplex( copy, 7, "-3.25", "1.5" );
				fn( &copy, &t );
				fn( &z, nullptr );
				l_wrong += (SameAnswer( &z.m_real, &t.m_real ) && SameAnswer( &z.m_imag, &t.m_imag )) ? 0 : 1;
			}
		}

		geAngleUnit = AU_radians;
		::snprintf( szDetail, sizeof( szDetail ), "%ld of 8 differ", l_wrong );
		CheckTrue( "inverses with logarithms in place", l_wrong == 0, szDetail );
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}
/*
	**Halves in a unit** -- 2026-10-23, found building complex stage 6b and
	fixed on Rich's word. sin 30 degrees and cos 60 are 1/2, which in an odd
	base is a half at every width -- 0.333... in base seven -- and in an even
	base a half at no places; lzc_polar's r cos theta is r/2 there, a half
	whenever r b^P is odd. An angle a hair from such a one, past any retake's
	reach, was rounded by half-even whichever side it lay: sin(30 - 7^-60)
	degrees at ten places in base seven answered 0.3333333334 for ...333.
	The real sine and cosine and lzc_polar now decide by the angle's side of
	the exact one, times the slope's sign there and r's -- see
	lzi_unit_half_decide -- and the complex sine and cosine, whose real part
	is moved by y too, by their leading terms -- see
	lzi_complex_trig_half_decide. Every pin is drawn as its census draws it,
	scored by that census's own reference, and kept where a mutant answers
	otherwise: the real ones by the builds that set neither decider, ignore
	the slope, take the wrong slope at 150 or 240 degrees and pick the
	twelfth by its floor, and two ordinary halves inside the reach by the
	build that takes any half for r/2; the polar by the builds
	that never ask the decider and ignore r's sign; the complex by the
	builds that never ask theirs, drop either leading term and leave the
	half unsigned.
*/

static void TestUnitHalves()
{
	Banner( "halves in a unit: the sine and cosine, lzc_polar, and the complex sine and cosine" );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 240 ];

	// 'r' the real function of a, 'p' lzc_polar( a, b ), 'c' the complex
	//	function of a + bi; 0 the sine, 1 the cosine...
	struct { char m_cKind; long m_lFunction; AngleUnit m_eUnit; long m_lPlaces; long m_lBase; const char* m_szA; const char* m_szB;
			 const char* m_szWantReal; const char* m_szWantImag; } kPins[] =
	{
		{ 'r', 0, AU_degrees, 0, 5, "3024.444444444444444444444444444444444444444444444444444444444444444444444444444444", "",
		  "0", "" },
		{ 'r', 1, AU_grads, 0, 12, "896.7EEEEEEEEEEEEEEEEEEEEEEEEEEE4", "",
		  "1", "" },
		{ 'r', 1, AU_degrees, 0, 6, "1540.00000000000000000000000000000000000000000000000000000000001", "",
		  "0", "" },
		{ 'r', 0, AU_degrees, 10, 3, "-21210.0000000000000000000000000000000000000000000000000000000000000000000000001", "",
		  "0.1111111112", "" },
		{ 'r', 1, AU_grads, 20, 15, "2E6.9EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE5", "",
		  "-0.77777777777777777778", "" },
		{ 'r', 1, AU_grads, 0, 12, "-E1.4000000000000000000000000D", "",
		  "-1", "" },
		{ 'r', 0, AU_degrees, 10, 16, "71F.570672C5155D38D80841E2F18DC51F664648FE0C289CE4771E98DB47CA354803A1E2F7242FD424F9", "",
		  "0.656C6CA769", "" },
		{ 'r', 0, AU_millicles, 10, 4, "100130.0212001231120120103331112331202011121132233101313212001121302322100000111220222", "",
		  "0.1102120311", "" },
		{ 'p', 0, AU_grads, 36, 12, "0.000000000000000000000000000000002605", "-56.7EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE3",
		  "0.000000000000000000000000000000001303", "-0.000000000000000000000000000000002202" },
		{ 'p', 0, AU_degrees, 36, 12, "0.00000000000000000000000000000001DE05", "4E.EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE",
		  "0.00000000000000000000000000000000E563", "0.000000000000000000000000000000017D23" },
		{ 'p', 0, AU_degrees, 0, 2, "-1011", "-11110.0000000000000000000000000000000000000000000001",
		  "-1010", "110" },
		{ 'p', 0, AU_millicles, 10, 3, "-0.0000000012", "10002.022222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222",
		  "-0.0000000011", "-0.0000000002" },
		{ 'c', 0, AU_millicles, 0, 12, "763.400000000000000000000000000000000000000000004", "-0.000000000000000000000000000000000000000000000000000000002",
		  "1", "0" },
		{ 'c', 0, AU_grads, 36, 9, "278.28888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888884", "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000005",
		  "-0.444444444444444444444444444444444444", "0.000000000000000000000000000000000000" },
		{ 'c', 0, AU_millicles, 0, 12, "2D8.80000000000000000000000000000000000000000000000000000000000000000000000000000002", "0.000000000000000000000000000000000007",
		  "1", "0" },
		{ 'c', 0, AU_millicles, 0, 12, "2D8.80000000000000000000000000000000000000000000000005", "-0.000000000000000000003",
		  "1", "0" },
		{ 'c', 0, AU_degrees, 0, 2, "11010010.000000000000000000000000000000000000000000000000000001", "0.0000000000000000000000000000000000000000000000000000000000000000001",
		  "-1", "0" },
		{ 'c', 0, AU_degrees, 1, 3, "112110.00000000000000000000000000000000000000000002", "0.0000000000000000000000000000000000000000000000000001",
		  "0.2", "0.0" },
		{ 'c', 1, AU_degrees, 0, 8, "1127.7777777777777777777777777777777777777777777777777777777777777777777773", "-0.0000000000000000000000004",
		  "-1", "0" },
		{ 'c', 1, AU_grads, 0, 6, "-341.1555555555555555555555555555555555555555555555555555555555555555554", "-0.0000000000000000000000002",
		  "-1", "0" },
	};

	for ( char c_kind : { 'r', 'p', 'c' } )
	{
		long		l_cases = 0;
		long		l_wrong = 0;
		long		l_pair = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			if ( pin.m_cKind != c_kind )
				continue;

			bool	b_ok = true;

			glWorkingPrecision = pin.m_lPlaces;
			geAngleUnit = pin.m_eUnit;

			if ( c_kind == 'r' )
			{
				CLZNumber	x;
				CLZNumber	t;
				CLZNumber	sine;
				CLZNumber	cosine;

				x.SetBase( pin.m_lBase, false );
				x.Assign( pin.m_szA );
				b_ok = ((pin.m_lFunction == 0) ? lz_sine( &x, &t ) : lz_cosine( &x, &t )) == kNoError;
				b_ok = b_ok && SpelledAs( &t, pin.m_szWantReal ) && (t.m_lBase == pin.m_lBase);

				// ...and the pair, which hands a part in doubt to its function...
				lz_sine_cosine( &x, &sine, &cosine );
				l_pair += SpelledAs( (pin.m_lFunction == 0) ? &sine : &cosine, pin.m_szWantReal ) ? 0 : 1;
			}
			else
			{
				CLZComplex	t;

				if ( c_kind == 'p' )
				{
					CLZNumber	r;
					CLZNumber	theta;

					r.SetBase( pin.m_lBase, false );
					r.Assign( pin.m_szA );
					theta.SetBase( pin.m_lBase, false );
					theta.Assign( pin.m_szB );
					b_ok = (lzc_polar( &r, &theta, &t ) == kNoError);
				}
				else
				{
					CLZComplex	z;

					SetComplex( z, pin.m_lBase, pin.m_szA, pin.m_szB );
					b_ok = (((pin.m_lFunction == 0) ? lzc_sine( &z, &t ) : lzc_cosine( &z, &t )) == kNoError);
				}

				b_ok = b_ok && SpelledAs( &t.m_real, pin.m_szWantReal ) && SpelledAs( &t.m_imag, pin.m_szWantImag ) &&
					   (t.m_real.m_lBase == pin.m_lBase) && (t.m_imag.m_lBase == pin.m_lBase);
			}

			geAngleUnit = AU_radians;
			l_cases++;

			if ( !b_ok && (l_wrong++ == 0) )
				s_first = std::to_string( pin.m_lFunction ) + " base " + std::to_string( pin.m_lBase ) + " " + pin.m_szA + " " + pin.m_szB;
		}

		if ( l_cases == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "no pins" );
		else if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld right", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_cases, s_first.c_str() );

		const char*	sz_label = (c_kind == 'r') ? "halves: the sine and cosine pins" :
							   ((c_kind == 'p') ? "halves: lzc_polar pins" : "halves: the complex sine and cosine pins");

		CheckTrue( sz_label, (l_cases > 0) && (l_wrong == 0), szDetail );

		if ( c_kind == 'r' )
		{
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ", l_pair, l_cases );
			CheckTrue( "halves: lz_sine_cosine at the same angles", l_pair == 0, szDetail );
		}
	}

	geAngleUnit = AU_radians;
	glWorkingPrecision = l_saved;
}

/*
	**A whole power takes no logarithm, and a fraction's negative power
	keeps its digits** -- 2026-10-11. exp_b of a whole number, exp_n with
	a whole exponent and root_n of an order whose reciprocal is whole all
	worked a logarithm at the answer's full width and then made no use of
	it: 22, 17 and 66 seconds for a 100,000-place answer. They go to the
	squaring now, as lz_power always did. And the squaring took a
	fraction's negative power by dividing one by a tiny positive one,
	which loses as many digits as the answer has whole ones -- 1,570 of
	2,715 random cases wrong, by up to past 10^400 units. Each pin is a
	census case HEAD got wrong, its expected string rendered in Python
	from the exact rational: x^n of the operand as held.
*/

struct SNegativePowerPin
{
	const char*	m_szFunction;
	long		m_lBase;
	long		m_lPlaces;
	const char*	m_szOperand;
	long		m_lPower;
	const char*	m_szExpect;
};

static const SNegativePowerPin kNegativePowers[] = {
	{ "exp_n", 2, 36, "0.00000000001011010101100111010101010111100110101111000110001000011011011111100000101011000111110110100001111011000100111", -6L,
	  "111111001100001101011010011110000010111111111101111010001001100.100000011100111010001101101110111010" },
	{ "power", 2, 36, "0.100001101010001000001111110010011011010101001110111100110111111110110010110011000110011110000111000000001000110100011111", -32L,
	  "110010110101000000010011111011.111001110111011101110101000010111010" },
	{ "raise_n", 2, 36, "0.000000001100011110111011100010101101011001000110100110001110000001001000101001000011001001010110101010001110101001100111", -16L,
	  "11010100001011111001100101110101111001000111110000010100111010000100111100110011011011000100010001100101110111101110010100111000100101.101110010111010110000010000001010111" },
	{ "exp_n", 3, 36, "0.0000000000000000000000001022120111000110210022001121110212201101122100201022", -1L,
	  "2021110220211022010110220.020221112011002011121020022210222122" },
	{ "power", 3, 36, "0.000000000000000000110120211020222101220212112111022101110012220221010202121", -4L,
	  "210011211000021221210222200011102111000121102211022211222020001010000012201.112210100000001102200110201100112001" },
	{ "raise_n", 3, 36, "-0.0221122111201102112200120221100200100220000002210012112002120102201122010011", -96L,
	  "12010101002122121200221021122211211201222122100211002001122000121121020202022000121222010121210010102.021221101021021000111002111120100211" },
	{ "exp_n", 5, 36, "0.2401440124012401240124012401240124012401240124012402", -92L,
	  "334103412211020401240303413034002.324104433342432411212030222231100104" },
	{ "power", 5, 36, "0.0000014433212400313234123201443241422313033142304443", -4L,
	  "12442112343242200310441.221120222231233324304131030440200441" },
	{ "raise_n", 5, 36, "0.2231043422030112340014321040344324230201400211014223", -145L,
	  "14141024110310412324001141422333123222444340220401000144342143.304231241124132204214114033311121104" },
	{ "exp_n", 7, 36, "0.0032040134454446021045630260334225425644655", -14L,
	  "2150166110024556014644121044362166.333304300241304632646155525405132351" },
	{ "power", 7, 36, "0.5251235346346256054126456466562512520665536", -64L,
	  "306066153.410663505401631141666225011106434013" },
	{ "raise_n", 7, 36, "0.026404144265312265231210632510053446143062", -26L,
	  "31000461512506254452303315652560050260.354164636000323443555110346042044040" },
	{ "exp_n", 8, 36, "0.000000002615671316731365355433261717717", -5L,
	  "3063057373246223204544107426250001123413131.451610743147134556234205366400740305" },
	{ "power", 8, 36, "-0.0000000000053230142125264546541661057744", -5L,
	  "-70350013160557713451444642563271570355774362214423237470.273416434555650641460424625732164677" },
	{ "raise_n", 8, 36, "0.0024365605075341217270243656050753412173", -10L,
	  "25536165705530607777777777.777777777776565365057325631020000000" },
	{ "exp_n", 10, 36, "0.0062187482634238", -32L,
	  "39948892321130457409159464328270073143041733230128730330942747210416669.468573161633298367225748018693407184" },
	{ "power", 10, 80, "0.0051969987945433771784", -12L,
	  "2576134933944450273190449835.30291195849582479124829401695540646529325673694175346767975910787430857990275802" },
	{ "raise_n", 10, 36, "0.00000000152759266", -6L,
	  "78696355876875304133289230388754136851830977898091480.360418282237507425254463597796572792" },
	{ "exp_n", 12, 36, "0.0000003484D4D581009284E03E3E5D6E79", -7L,
	  "4023872419163E81678963702679543916103078108889.95692977712E488D3419953DE32087576621" },
	{ "power", 12, 36, "0.000009D4167E774014996D63762844DD49", -2L,
	  "159253452D7.D310413137D9337E4917962EE23361920259" },
	{ "raise_n", 12, 36, "-0.0001D23185509DE89582909E308DE57255", -14L,
	  "396989D7843937996904DE6485054E149388947E4E2E81D5D7221.35410816408E719007782E46D64720D027E8" },
	{ "exp_n", 16, 36, "0.000000001CDBF2639EFBE23A6C01EE", -4L,
	  "182FE91CCC202018133BC4C478639DD0C1A4.A59E7FB0EA9520856C16D1360926B2566387" },
	{ "power", 16, 36, "0.078A567FB065C12B801FC3BF5AE7E4", -6L,
	  "5B46D260.B611131A23E9A345DEF4E9BBFF806DDA4B1A" },
	{ "raise_n", 16, 36, "0.00008992BB59070252C09734648B9B", -16L,
	  "50BC0C6A22072E62F4E5B748D90FBBFC0A712473F6D363B698A46AFB570CD2CF56E1.C8F788F7A0F548AC8F241DA83517944493B7" },
};

static void TestWholePowers()
{
	Banner( "Whole powers, and a fraction's negative ones" );

	typedef long (*Binary)( CLZNumber*, CLZNumber*, CLZNumber* );

	long	l_saved = glWorkingPrecision;
	char	szDetail[ 200 ];

	glWorkingPrecision = 36;

	// The logarithms each call takes. A whole power should take none; the
	//	controls, a fractional power and a root of order three, do...
	auto	logarithms = []( Binary fn, const char* szLeft, const char* szRight ) -> long
	{
		CLZNumber	x( szLeft ), y( szRight ), answer;
		long		l_before = lzi_logarithm_count();

		fn( &x, &y, &answer );

		return lzi_logarithm_count() - l_before;
	};

	// lz_exp_b reads ln(base) from a table to 120 places and past them
	//	works it with lz_log_e, once for each width -- so the count sees
	//	it only for a wide answer, at a width nothing here has asked
	//	for...
	auto	exp_b_logarithms = []( const char* szOperand, long lPlaces ) -> long
	{
		CLZNumber	x( szOperand ), answer;
		long		l_before = lzi_logarithm_count();

		glWorkingPrecision = lPlaces;
		lz_exp_b( &x, &answer );
		glWorkingPrecision = 36;

		return lzi_logarithm_count() - l_before;
	};

	/*
		**From an empty constant cache.** Since 2026-10-13 a computed ln(base)
		serves any narrower width, rounded, so once anything earlier in the
		run had worked it wider than these, exp_b took no logarithm with a
		fraction or without, and the whole power's nought could not have
		failed. Taken down and set up again, as the threads gate does
		between its phases, the widths below are new to the cache...
	*/

	lzu_takedown();
	lzu_setup();

	long	l_whole = exp_b_logarithms( "321", 41 ) + exp_b_logarithms( "-5", 36 ) +
					  logarithms( lz_exp_n, "1.2345", "1000" ) + logarithms( lz_exp_n, "1.2345", "-17" ) +
					  logarithms( lz_root_n, "1.2345", "0.001" ) + logarithms( lz_root_n, "3", "-0.5" );
	long	l_control = exp_b_logarithms( "321.5", 43 ) * logarithms( lz_exp_n, "1.2345", "2.5" ) *
						logarithms( lz_root_n, "2", "5" );	// order three goes to root_3 since 2026-10-15

	::snprintf( szDetail, sizeof( szDetail ), "%ld for six whole powers; the three controls %s", l_whole,
				(l_control > 0) ? "took theirs" : "took none -- the count is not counting" );

	CheckTrue( "no logarithm for a whole power", (l_whole == 0) && (l_control > 0), szDetail );

	// The same digits as lz_power, which always took a whole exponent to
	//	the squaring -- at a million, and an order of 1E-6, too...
	auto	rendered = []( Binary fn, long lBase, const char* szLeft, const char* szRight ) -> std::string
	{
		CLZNumber	x( szLeft ), y( szRight ), answer;

		if ( lBase != 10 )
		{
			x.SetBase( lBase );
			y.SetBase( lBase );
		}

		fn( &x, &y, &answer );

		char*		sz_answer = answer.GetStringRep();
		std::string	s_answer = (sz_answer != nullptr) ? sz_answer : "(null)";

		::free( sz_answer );

		return s_answer;
	};

	struct { Binary m_fn; long m_lBase; const char* m_szRight; const char* m_szPower; } kSame[] =
	{
		{ lz_exp_n, 10, "1000", "1000" }, { lz_exp_n, 7, "40", "40" }, { lz_exp_n, 3, "-23", "-23" },
		{ lz_root_n, 10, "0.001", "1000" }, { lz_root_n, 12, "0.25", "4" }, { lz_root_n, 16, "-0.0625", "-16" },
		{ lz_exp_n, 10, "1E6", "1E6" }, { lz_root_n, 10, "1E-6", "1E6" },
	};

	long	l_differ = 0;

	for ( auto& same : kSame )
		if ( rendered( same.m_fn, same.m_lBase, "1.2345", same.m_szRight ) !=
			 rendered( lz_power, same.m_lBase, "1.2345", same.m_szPower ) )
			l_differ++;

	CLZNumber	whole( "37" ), power_of_ten;

	lz_exp_b( &whole, &power_of_ten );

	char*	sz_power = power_of_ten.GetStringRep();
	bool	b_shift = (sz_power != nullptr) && (std::string( sz_power ) == "1" + std::string( 37, '0' ));

	::free( sz_power );

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu differ from power; 10^37 %s", l_differ,
				sizeof( kSame ) / sizeof( kSame[ 0 ] ), b_shift ? "a one and noughts, exact" : "NOT a one and noughts" );

	CheckTrue( "the digits power gives", (l_differ == 0) && b_shift, szDetail );

	// A fraction's negative powers, each a case HEAD got wrong...
	struct { const char* m_szName; Binary m_fn; } kFunctions[] =
	{
		{ "raise_n", lz_raise_n }, { "exp_n", lz_exp_n }, { "power", lz_power }
	};

	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kNegativePowers )
	{
		glWorkingPrecision = pin.m_lPlaces;

		CLZNumber	x, y, answer;

		x.SetBase( pin.m_lBase, true );
		y.SetBase( pin.m_lBase, true );
		x.Assign( pin.m_szOperand );
		y.Assign( pin.m_lPower );

		long	l_code = kFailed;

		for ( auto& f : kFunctions )
			if ( ::strcmp( f.m_szName, pin.m_szFunction ) == 0 )
				l_code = f.m_fn( &x, &y, &answer );

		char*	sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
			if ( l_wrong++ == 0 )
				s_first = std::string( pin.m_szFunction ) + " base " + std::to_string( pin.m_lBase ) + " power " +
						  std::to_string( pin.m_lPower ) + ": " + ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");

		::free( sz_got );
	}

	size_t	n_pins = sizeof( kNegativePowers ) / sizeof( kNegativePowers[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character, in eight bases", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "a fraction's negative powers", l_wrong == 0, szDetail );

	// 2^100 exactly, so spelled short; and the root of order -1/20 is the
	//	-20th power, through the same squaring...
	glWorkingPrecision = 36;

	CLZNumber	half( "0.5" ), hundred( "-100" ), half_answer;

	lz_power( &half, &hundred, &half_answer );

	CheckRendering( "0.5 to the -100th", &half_answer, "1267650600228229401496703205376" );

	CLZNumber	small( "0.0012345" ), order( "-0.05" ), root;

	lz_root_n( &small, &order, &root );

	CheckRendering( "0.0012345's root of order -1/20", &root,
					"14797151306014377797780816227954028197677577748612045384417.205562598726002913770027081291758621" );

	// A hair under one to the -10^12th: the reciprocal's error is
	//	multiplied by the exponent, so it has to be taken at the width
	//	the exponent raised -- at the operand's, it came out 3 units low.
	//	Scored in Python's Decimal...
	const char	kNearOneAnswer[] = "2.718281828460404376274518239849514020";
	CLZNumber	near_one( "0.999999999999" ), trillion( "-1000000000000" ), tiny_order( "-1E-12" ), near_answer;

	lz_exp_n( &near_one, &trillion, &near_answer );

	CheckRendering( "1 - 1E-12 to the -10^12th", &near_answer, kNearOneAnswer );

	lz_root_n( &near_one, &tiny_order, &near_answer );

	CheckRendering( "1 - 1E-12's root of order -1E-12", &near_answer, kNearOneAnswer );

	/*
		Near one, a power's size is its exponent times a logarithm that
		the estimate took from sixteen digits, mostly cancellation -- and
		sixteen bits in base two. HEAD answered the first nought, the
		second wrong from its 30th digit, the third nought, refused the
		fourth as too large, and never answered the fifth, sizing it at
		67,000 places. Each the operand as written, in its base, scored
		in Python's Decimal. The one that hung goes last...
	*/

	struct { Binary m_fn; long m_lBase; const char* m_szOperand; const char* m_szExponent; const char* m_szExpect; } kNearOne[] =
	{
		{ lz_power, 10, "1.00000000000000000001", "10000000000000000000",
		  "1.105170918075647624811155241031208844" },
		{ lz_raise_n, 10, "0.99999999999390467", "-8828440789890",
		  "234621375593218166810711.213916404647442761483990395105165803" },
		{ lz_power, 5, "1.0000000000000000000000002113021401221243433314343213", "-852511000000000000",
		  "0.000044000140132310123124140241322313" },
		{ lz_raise_n, 3, "0.2222222222222222222222222222222001012001100212101012212002021002011202111002", "-29730701641038880",
		  "212211122020122100022002122112221201112121.220121010110220122202010111120221102" },
		{ lz_power, 5, "0.4444444444444444444444431324222212231341321000004243", "-16533000000000000",
		  "1.242430112011003304101230443121104331" },
	};

	l_wrong = 0;
	s_first.clear();

	for ( auto& near : kNearOne )
	{
		CLZNumber	x, y, answer;

		x.SetBase( near.m_lBase, true );
		x.Assign( near.m_szOperand );
		y.Assign( near.m_szExponent );

		if ( near.m_lBase != 10 )
			y.SetBase( near.m_lBase );

		long	l_code = near.m_fn( &x, &y, &answer );
		char*	sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, near.m_szExpect ) != 0) )
			if ( l_wrong++ == 0 )
				s_first = std::string( "base " ) + std::to_string( near.m_lBase ) + " to the " + near.m_szExponent +
						  ": code " + std::to_string( l_code ) + ", " + ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");

		::free( sz_got );
	}

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character", sizeof( kNearOne ) / sizeof( kNearOne[ 0 ] ) );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, sizeof( kNearOne ) / sizeof( kNearOne[ 0 ] ), s_first.c_str() );

	CheckTrue( "a power's size near one", l_wrong == 0, szDetail );

	/*
		And not only near one: in base two sixteen digits are sixteen
		bits, so 1.002 to the 25 millionth, 2^72,064, was sized 540
		places short and came out wrong in its last 560 characters of
		72,100. The length and the tail, from Python's Decimal at 23,000
		digits; the operand is 1.002 converted, as written...
	*/

	{
		CLZNumber	x( "1.002" ), y( "25000000" ), answer;

		x.SetBase( 2 );
		y.SetBase( 2 );
		lz_power( &x, &y, &answer );

		char*		sz_got = answer.GetStringRep();
		std::string	s_got = (sz_got != nullptr) ? sz_got : "";
		const char	kTail[] = "00001101011.100110010100111000011011110001011110";
		bool		b_right = (s_got.size() == 72100) && (s_got.compare( s_got.size() - 48, 48, kTail ) == 0);

		::free( sz_got );
		::snprintf( szDetail, sizeof( szDetail ), "%zu characters, the last 48 %s", s_got.size(), b_right ? "right" : "NOT right" );

		CheckTrue( "1.002 to the 25 millionth, base two", b_right, szDetail );
	}

	/*
		A logarithm to a base near one divides by a small logarithm, and
		has to be worked wider by the places one over it has. HEAD got
		219 of 400 such wrong; the first, in base three, it got right by
		luck, its estimate of ln n coming out nought and the probe after
		it widening by far too much, and an exact estimate lost it. Base
		on the left, as lz_log_n takes it; each scored in Python's
		Decimal...
	*/

	struct { long m_lBase; long m_lPlaces; const char* m_szBase; const char* m_szValue; const char* m_szExpect; } kLogNearOne[] =
	{
		{ 3, 36, "1.0000000000000000010012022221122202020222110002101222010001000012200001121212", "0.000000102011122102220112011001000222020020020102200202102120201021202011202",
		  "-21000100010111020002.012022101121110120121222221122122020" },
		{ 16, 36, "1.0000001E6B7D3023FB2DC35B0AB7BE", "5F178.E66666666666666666666666666666",
		  "6C548150.678298A806B945F17CE65F461E3666CA3821" },
		{ 10, 80, "1.0000035525", "75391.01",
		  "3161284.52193705638695557111059049289270347704822269048891856639578426257438198622178107" },
		{ 12, 36, "1.0000000000031572E557DE9673881D0E3", "0.D99899536753D9D599437362540E752677",
		  "-494E0D4527D.48760623EE671D6587E8760D615046240416" },
		{ 7, 36, "0.6666666646603511500525304333330566421450264", "23152455",
		  "-10123340261.624360206642141255620156064063621666" },
	};

	l_wrong = 0;
	s_first.clear();

	for ( auto& log : kLogNearOne )
	{
		glWorkingPrecision = log.m_lPlaces;

		CLZNumber	base_n, value, answer;

		base_n.SetBase( log.m_lBase, true );
		value.SetBase( log.m_lBase, true );
		base_n.Assign( log.m_szBase );
		value.Assign( log.m_szValue );

		long	l_code = lz_log_n( &base_n, &value, &answer );
		char*	sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, log.m_szExpect ) != 0) )
			if ( l_wrong++ == 0 )
				s_first = std::string( "base " ) + std::to_string( log.m_lBase ) + ": " +
						  ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");

		::free( sz_got );
	}

	glWorkingPrecision = 36;

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character", sizeof( kLogNearOne ) / sizeof( kLogNearOne[ 0 ] ) );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, sizeof( kLogNearOne ) / sizeof( kLogNearOne[ 0 ] ), s_first.c_str() );

	CheckTrue( "a logarithm to a base near one", l_wrong == 0, szDetail );

	// A negative has a real root of order one and minus one, the only
	//	whole orders with a whole reciprocal, and the sign is put back
	//	once, on the magnitude's power...
	CLZNumber	minus_eight( "-8" ), plus_one( "1" ), minus_two( "-2" ), minus_one( "-1" ), signed_root;

	lz_root_n( &minus_eight, &plus_one, &signed_root );

	CheckRendering( "-8's root of order 1", &signed_root, "-8.000000000000000000000000000000000000" );

	lz_root_n( &minus_two, &minus_one, &signed_root );

	CheckRendering( "-2's root of order -1", &signed_root, "-0.500000000000000000000000000000000000" );

	// An order a hair off a quarter, past the width, is not a quarter: its
	//	whole reciprocal is checked exactly, and the root differs from
	//	10^120 from its 37th digit. Scored in Python's Decimal...
	CLZNumber	large( "1E30" ), near_quarter( "0.2500000000000000000000000000000000000001" ), near_root;

	lz_root_n( &large, &near_quarter, &near_root );

	CheckRendering( "an order just off a quarter", &near_root,
					"999999999999999999999999999999999999889475915536285807167136410175150518041299124802875497335700445869146762908763361500.383752213445526109043060479110759873" );

	glWorkingPrecision = l_saved;
}

/*
	**A root of a small order that is not a unit fraction** -- 2026-10-12.
	root_n takes such a root as e to ln x over the order, and the
	logarithm's last-place error, divided by a small order, grew by one
	over it: near one, where a tiny order still gives an ordinary answer,
	HEAD got 115 of 300 wrong and refused two. Both are worked wider by
	the places one over the order has now. Each pin is a census case
	cf1e5cb got wrong, scored in Python's Decimal against the order as
	held -- which in bases three, five and seven is 1.5 to so many
	places, since it never ends there.
*/

static void TestRootOfASmallOrder()
{
	Banner( "A root of a small order" );

	struct { long m_lBase; long m_lPlaces; const char* m_szValue; const char* m_szOrder; const char* m_szExpect; } kPins[] =
	{
		{ 3, 36, "1.000000000000000000000000010120221002200222010022011200012110111200100122112", "-0.000000000000000000000000011111111111111111111111111111112",
		  "0.110011021211210110211011121201111211" },
		{ 5, 36, "1.000000000000000000000024100234201334443202014201442", "0.00000000000000000000000122222222222222222222222222222",
		  "403320.013321110032003323143423201121334223" },
		{ 7, 36, "1.0000000000000000000412601230055341031425133", "-0.000000000000000013333333333333333333334",
		  "0.664132140515063561461253056366055433" },
		{ 8, 36, "1.0000000000000000016213500411625454340342", "0.0000000014",
		  "1.000000001140760026666473707136072236" },
		{ 10, 36, "1.00000000000807649", "-0.00000015",
		  "0.999946158182856217627590045101700865" },
		{ 16, 36, "0.FFFFFFFFFFFFFFC5B9F32846304CD3", "-0.0000000018",
		  "1.0000026D95E14508E2066D4416EB1D607967" },
	};

	long		l_saved = glWorkingPrecision;
	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kPins )
	{
		glWorkingPrecision = pin.m_lPlaces;

		CLZNumber	value, order, answer;

		value.SetBase( pin.m_lBase, true );
		order.SetBase( pin.m_lBase, true );
		value.Assign( pin.m_szValue );
		order.Assign( pin.m_szOrder );

		long	l_code = lz_root_n( &value, &order, &answer );
		char*	sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
			if ( l_wrong++ == 0 )
				s_first = std::string( "base " ) + std::to_string( pin.m_lBase ) + ": " +
						  ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 40 ) : "(null)");

		::free( sz_got );
	}

	glWorkingPrecision = l_saved;

	char	szDetail[ 200 ];
	size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character, in six bases", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "near one, of a small order", l_wrong == 0, szDetail );
}

static void TestExpBySplitting()
{
	Banner( "e^x by binary splitting" );

	/*
		**Part B of the slowness item, 2026-10-13.** lz_exp_e takes the
		bit-burst kernel, lzi_exp_by_splitting, past a floor that depends on
		the base -- see lzi_exp_split_floor. Pinned four ways: against
		Python's Decimal with the kernel forced on and forced off; the two
		ways agreeing to the character in every base, both signs, with whole
		parts and with fractions longer than the width; the measured floors
		taking the kernel where they say; and the kernel itself indifferent
		to its first chunk.
	*/

	long	l_savedFloor = lzi_set_exp_split_floor( kExpSplitFloorMeasured );
	long	l_saved = glWorkingPrecision;

	// From Python's Decimal, rounded half-even by the kept digit, none
	//	within b^-6 of a half. Every width is past the 48 digits the table
	//	path serves, so forcing the floor to nought does take the kernel...
	struct { long m_lBase; long m_lPlaces; const char* m_szArg; const char* m_szExpect; } kPins[] =
	{
		{ 2, 240, "0.1011",
		  "1.111111010001110111100110000110000010111110001100100010011101001011000011101101101101000010001100011001011001011100100010010000100010010011100001000101001111010101011011000001001100011101100011101001010100000101000011101110101000111010010011" },
		{ 2, 240, "-0.110100111011101011110001101000111010001011110000101101101",
		  "0.011011111111010010111110000001111001100001011111111011011111110101110001010101100010100101100000011101111100111111000010101000010110011000110110010110000000100101110001111010110010001100101111101001110110101111010101100101110111001001110110" },
		{ 2, 240, "101.01101",
		  "11011110.110010110110011011110011101101011000001010001101001100011000100101110000100101001011000110110001001100000111010010100100001111011110101000100100111001110011100100111110001100111010111111100010111101010110010111011000001100111111001010000111" },
		{ 3, 150, "0.2",
		  "1.221120220020201121121221100002202122212202221010012112210001010001021211202121100222002211002010211120221120201012221000122210110002201021012021101111" },
		{ 3, 150, "-0.12120112001021210221012",
		  "0.112112101110010222020101012010020101002212022011202002212201210111200121110120111200112010101002001021212212210101021212102022221120121222202222111001" },
		{ 3, 150, "-11.0221",
		  "0.000100210202211200220100100100002220111011111101010102001220222221110012221120001012111110000012121100202120220102211212210110201011020120212002100120" },
		{ 7, 100, "0.6543",
		  "2.4334662453343021246131116001632231021153005532436012022303064225012260353512430203405065654211162243" },
		{ 7, 100, "-0.0000000000000000000000000000000000000002",
		  "0.6666666666666666666666666666666666666665000000000000000000000000000000000000000200000000000000000000" },
		{ 7, 100, "12.304",
		  "51455.4363606256230452162166351106056553056626333222126631304134534404141231236656510624115452406461162546" },
		{ 10, 100, "0.5",
		  "1.6487212707001281468486507878141635716537761007101480115750793116406610211942156086327765200563666430" },
		{ 10, 100, "-0.31415926535897932384626433832795028841971693993751",
		  "0.7304026910486456108725952838122532797466894912741328330234751317310882080748751532109910268979435530" },
		{ 10, 100, "-6.25",
		  "0.0019304541362277092422135119756507321435854079192435762119985657240298793715072711306642641269732257" },
		{ 12, 80, "0.DE1",
		  "2.59D24EE98D44951723D384563390432215639E2E006D405D376117614797867D1843884938D01226" },
		{ 12, 80, "-0.5",
		  "0.7DE201E82E468E242970516051900E697E7506D060D6424E049449E1550D752327499471977433E8" },
		{ 12, 80, "3.DD",
		  "41.658D87DD7D419110715E1EE574010307838372857ED0E228527D9957DEDD4E346178290845E86285" },
		{ 16, 60, "0.C",
		  "2.1DF3B68CFB9EF7A986ADDC7DCEE21F256208D5AAA84203C00346E51C87DA" },
		{ 16, 60, "-0.123456789ABCDEF0123456789ABCDEF",
		  "0.EE6D8115EC4E94E3A49F8371FE2896268396B7DF4CAF4E0C02FFCC307346" },
		{ 16, 60, "1A.8",
		  "4B22A1B5C3.27A2FD85A4A054EF3317DA20C3B5F49143B2B54A52649731CB1835A2F93D" },
	};

	{
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
			for ( long l_floor : { 0L, LONG_MAX } )
			{
				lzi_set_exp_split_floor( l_floor );
				glWorkingPrecision = pin.m_lPlaces;

				CLZNumber	x, answer;

				x.SetBase( pin.m_lBase, true );
				x.Assign( pin.m_szArg );

				long	l_runs = lzi_exp_splittings();
				long	l_code = lz_exp_e( &x, &answer );
				bool	b_ran = (lzi_exp_splittings() > l_runs);
				char*	sz_got = answer.GetStringRep();

				if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) ||
					 (b_ran != (l_floor == 0)) )
					if ( l_wrong++ == 0 )
						s_first = std::string( "base " ) + std::to_string( pin.m_lBase ) + " exp(" + pin.m_szArg + ")" +
								  ((l_floor == 0) ? " split" : " series") + (b_ran ? " ran" : " did not run") + ": " +
								  ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 30 ) : "(null)");

				::free( sz_got );
			}

		char	szDetail[ 200 ];
		size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character both ways, in six bases", n_pins );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong, first %s", l_wrong, 2 * n_pins, s_first.c_str() );

		CheckTrue( "against Decimal, by the kernel and by the series", l_wrong == 0, szDetail );
	}

	/*
		**The kernel against the series, every base.** Two widths a base,
		the first just past the table path's 48 digits, and five arguments:
		a short fraction, a fraction as long as the width and its negative,
		a fraction twenty places longer than the width under a whole part,
		and a negative whole part under a short fraction. The fractions are
		the same digits in every run, from a fixed sequence.
	*/

	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;
		uint64_t	u_state = 20261013;

		for ( long l_base = 2; l_base <= 16; l_base++ )
		{
			const char*	sz_digits = (l_base == 12) ? "0123456789DE" : "0123456789ABCDEF";
			long		l_first = CLZNumber::EquivalentPlaces( 60, 10, l_base );

			for ( long l_places : { l_first, 3 * l_first } )
			{
				auto	digits = [&]( long lCount ) {
					std::string	s;

					for ( long i = 0; i < lCount; i++ )
					{
						u_state = (u_state * 6364136223846793005ULL) + 1442695040888963407ULL;
						s += sz_digits[ (u_state >> 33) % l_base ];
					}

					return s;
				};

				std::string	s_long = digits( l_places );
				std::string	args[] =
				{
					"0." + digits( 5 ), "0." + s_long, "-0." + s_long,
					"11." + digits( l_places + 20 ), "-100." + digits( 12 ),
				};

				for ( auto& s_arg : args )
				{
					std::string	s_answer[ 2 ];

					for ( int i = 0; i < 2; i++ )
					{
						lzi_set_exp_split_floor( (i == 0) ? 0 : LONG_MAX );
						glWorkingPrecision = l_places;

						CLZNumber	x, answer;

						x.SetBase( l_base, true );
						x.Assign( s_arg.c_str() );
						lz_exp_e( &x, &answer );

						char*	sz_got = answer.GetStringRep();

						s_answer[ i ] = (sz_got != nullptr) ? sz_got : "(null)";
						::free( sz_got );
					}

					l_cases++;

					if ( s_answer[ 0 ] != s_answer[ 1 ] )
						if ( l_wrong++ == 0 )
							s_first = "base " + std::to_string( l_base ) + " at " + std::to_string( l_places ) + ", exp(" +
									  s_arg.substr( 0, 12 ) + "...)";
				}
			}
		}

		char	szDetail[ 200 ];

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld identical, fifteen bases", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, first %s", l_wrong, l_cases, s_first.c_str() );

		CheckTrue( "the kernel and the series agree", l_wrong == 0, szDetail );
	}

	/*
		**The measured floors take the kernel where they say**: not at half
		the floor, at the floor itself. The floor stands against the width
		the series works at, which is the caller's and lz_exp_e's guard
		above it, so the caller's floor already reaches it. Base two's floor
		is nought, which means wherever the table path stops, so the width
		at the floor is never under 60 digits, nor the one below it under
		twenty places...
	*/

	{
		lzi_set_exp_split_floor( kExpSplitFloorMeasured );

		long		l_wrong = 0;
		std::string	s_first;

		for ( long l_base : { 2L, 3L, 7L, 10L, 15L, 16L } )
		{
			for ( int i = 0; i < 2; i++ )
			{
				long	l_floor = lzi_exp_split_floor( l_base );
				long	l_past = CLZNumber::EquivalentPlaces( 60, 10, l_base );

				glWorkingPrecision = (i == 0) ? std::max( l_floor / 2, 20L ) : std::max( l_floor, l_past );

				CLZNumber	x, answer;

				x.SetBase( l_base, true );
				x.Assign( "0.1" );

				long	l_runs = lzi_exp_splittings();

				lz_exp_e( &x, &answer );

				if ( (lzi_exp_splittings() > l_runs) != (i == 1) )
					if ( l_wrong++ == 0 )
						s_first = "base " + std::to_string( l_base ) + " at " + std::to_string( glWorkingPrecision );
			}
		}

		CheckTrue( "the floors take the kernel where they say", l_wrong == 0, s_first.c_str() );
	}

	/*
		**The kernel's own first chunk changes nothing but its speed.** Forced
		from one place to seventeen, against the default, the fraction as long
		as the width so every chunk is filled...
	*/

	{
		long		l_wrong = 0;
		std::string	s_first;

		for ( long l_base : { 3L, 10L } )
		{
			glWorkingPrecision = 200;

			CLZNumber	x, by_default;

			x.SetBase( l_base, true );
			x.Assign( (l_base == 3) ? "-0.21021200112022110201221002110120121002111020102101" : "0.71828182845904523536028747135266249775724709369995" );
			lzi_exp_by_splitting( &x, &by_default );

			for ( long l_first : { 1L, 2L, 3L, 7L, 17L } )
			{
				CLZNumber	forced;

				lzi_exp_by_splitting( &x, &forced, l_first );

				if ( forced.Compare( &by_default ) != kEqualTo )
					if ( l_wrong++ == 0 )
						s_first = "base " + std::to_string( l_base ) + ", first chunk " + std::to_string( l_first );
			}
		}

		CheckTrue( "the first chunk changes nothing but the speed", l_wrong == 0, s_first.c_str() );
	}

	glWorkingPrecision = l_saved;
	lzi_set_exp_split_floor( l_savedFloor );

	CheckTrue( "the floor is put back", lzi_set_exp_split_floor( l_savedFloor ) == l_savedFloor, "" );
}

static void TestLogByAgm()
{
	Banner( "ln x by the arithmetic-geometric mean" );

	/*
		**Part B's logarithm, stage 1, 2026-10-13.** lzi_log_by_agm is called
		by nothing in the engine yet; this pins it four ways -- against
		Python's Decimal past the constants' 120 places, so pi and ln(base)
		are computed; against lz_log_e in every base, both sides of one,
		large, small and long; exact powers of the base; and the constants
		served at a narrower width by rounding a wider one down.
	*/

	long	l_saved = glWorkingPrecision;

	// lz_log_e by Newton and the series throughout, so that comparing the
	//	kernel with it compares two independent things -- since 2026-10-14
	//	lz_log_e takes the AGM itself, past its floor. The dispatch is
	//	pinned last, with the measured floors back...
	long	l_savedFloor = lzi_set_log_agm_floor( LONG_MAX );

	// From Python's Decimal, rounded half-even by the kept digit, none
	//	within b^-6 of a half...
	struct { long m_lBase; long m_lPlaces; const char* m_szArg; const char* m_szExpect; } kPins[] =
	{
		{ 2, 450, "1.1011",
		  "0.100001011111001110010111001000010010100101010100000101011011010011000100101111011101100110011110111111111110011010011011011001000011011001101111101110111111001101011101001111101101000100011010100100001001001101101111010100001100010100101111101101111111011001100000111001001010111111110100011011101010010100010110100011111001101100000000100110011010101111011110111001110100000000111100011100100101100110000110001111100101011100010110101100001101010011" },
		{ 2, 450, "0.000101101110001",
		  "-10.011010100010111111000100010001001001110011011000000110011101010111001001110010100100001010011101001110011000000100010001110010011000010001110000010001101001110011100100111111000101000001010011010111111000010110111101010101000101010010010001111010001010011111001110010100000001100110000110110101010011101110001010001000100011111011001110000111010111000000011001101100110000111001011001001101100100000110001011110100110110011010110110010110111110100011" },
		{ 2, 450, "1000",
		  "10.000101000101011001000111111001110111010101101110011011010000001101011101101010110001101011001000000010111101100011100100000011011100001011011001110010010111001101010111110010100010001010001000100111100010011101000110000100101010001100000000111011101000001110110111001010010110001001100001010010011100001100000101110010010000000010111111111110001111000111100000111011100101000100110010110001111000110000001010101000010100001110101000011000111100110010" },
		{ 3, 300, "2.1201",
		  "0.221110111211010021111122021022111101222000101100221121200221201120221002120202100222011102012121202212012220111100112100101100011020011010100122221210010022020111001102110121100012220022200121211012201220020020022202100000000212010112210001020101120121202012200101200011220111021112002111022210220002" },
		{ 3, 300, "1.0000000000000000212",
		  "0.000000000000000021122222222222221220121111111111201001201011111110022111210102021000111221122120222122000122112120200221122010201200211121202202011111011102010022120021001122002120121121211022121112201121102010000212222222221020101100111000012220000120222212200021020010010102022022111022221121210220" },
		{ 3, 300, "0.01",
		  "-2.012022202222020010012020201200212021122111121211120100012020022010001222121121112221000022101100020101011021012111021111012002220201100012121112212211210012200220001021121111021210120020110211111202112211211120211101001022102212110120222020010121121101100221000021212022101200100122001222010112002021" },
		{ 7, 200, "6.543",
		  "1.62651146424026550663242354541340404233610313233415044566051531106561650043533133221000041542652312045140462003432231462500315466663062261625032015521534054403105534414331165155430255163532203312126630" },
		{ 7, 200, "0.00000312",
		  "-13.34232236140120016144664550152243633510543634531015022146416544354262341350112501264526551462544604103123223613231215161152015540156524441015266630224615361344405635601515063540404223116040413006640432" },
		{ 7, 200, "100000.5",
		  "12.50515162601140161024002013261366655232140312502642531164334365165545151604166252466610203102463000235631030156364042101311203112503425004616512065154565602213423326051361205116622653010340022362651100" },
		{ 10, 150, "1.2345",
		  "0.210666029803097141833911716856634657897218919766562764580235655270835391209604650634752605527100673369178376151350064144416628062000711332211449382404" },
		{ 10, 150, "0.99999999999999999999987",
		  "-0.000000000000000000000130000000000000000000008450000000000000000000732333333333333333333404735833333333333333340759193333333333333334137801500000000000" },
		{ 10, 150, "12345678.9",
		  "16.328816673173972349175941893278840002112580352391628999151762044124644857578999297420119215100704397136508744722998071786874422570688735445293986127148" },
		{ 12, 150, "DE.1",
		  "4.D615340E7D2877037836850217D33D7D8E0760009840626966732E5EDD117E721D9242599353561D3ED2972151E623DDED05E14E26271816196D311969043880959037D9D97277E6054054" },
		{ 12, 150, "0.000E",
		  "-7.6601410D5947E109928784806569670252E46123295810091064450E1435DD7706144695474932E9E21E22D728927ED03740039D394193E7954848D85E0E781943ED7E56E162E2097631EE" },
		{ 12, 150, "1000",
		  "7.55590D5E048339318E203245157254D22750920E560925472366691E62820E6D622765374154700D2444E20480E6496D2D273E287369504E17308535222D7D396653120E598915508E319E" },
		{ 16, 130, "1.8",
		  "0.67CC8FB2FE612FCADA35D9BD014886067D20FFB34547D7C2B38AD78EC59E3B60C2DF0CB19EDAEBB7FADCA437B8A073C4752D66D15D6F9B5DA7E09DC7FEBB0E3F48" },
		{ 16, 130, "FFFF.FF",
		  "B.17217E7D1CF71ABC9DE5E42AA9DA156B86AB7F33F6339D5C290BF944836C14ED9E624B5847377FC8B073497CD7B60B0E91736EBE75B896F76B51A872228FEB17F7" },
		{ 16, 130, "0.0001",
		  "-B.17217F7D1CF79ABC9E3B39803F2F6AF40F343267298B62D8A0D175B8BAAFA2BE7B876206DEBAC98559552FB4AFA1B10ED2EAE35C138214427573B291169B8253E9" },
	};

	{
		long		l_wrong = 0;
		long		l_cycles = lzi_agm_log_cycles();
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			glWorkingPrecision = pin.m_lPlaces;

			CLZNumber	x, answer;

			x.SetBase( pin.m_lBase, true );
			x.Assign( pin.m_szArg );
			lzi_log_by_agm( &x, &answer );

			char*	sz_got = answer.GetStringRep();

			if ( (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
				if ( l_wrong++ == 0 )
					s_first = std::string( "base " ) + std::to_string( pin.m_lBase ) + " ln(" + pin.m_szArg + "): " +
							  ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 30 ) : "(null)");

			::free( sz_got );
		}

		char	szDetail[ 200 ];
		size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character, in six bases", n_pins );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong, first %s", l_wrong, n_pins, s_first.c_str() );

		CheckTrue( "against Decimal", (l_wrong == 0) && (lzi_agm_log_cycles() > l_cycles), szDetail );
	}

	/*
		**The AGM against lz_log_e, every base**, two widths a base past 48
		digits and six arguments: a short one, one as long as the width,
		just above one and just below, large and small. The digits come from
		a fixed sequence, the same every run.
	*/

	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;
		uint64_t	u_state = 20261013;

		for ( long l_base = 2; l_base <= 16; l_base++ )
		{
			const char*	sz_digits = (l_base == 12) ? "0123456789DE" : "0123456789ABCDEF";
			long		l_first = CLZNumber::EquivalentPlaces( 60, 10, l_base );

			for ( long l_places : { l_first, 3 * l_first } )
			{
				auto	digits = [&]( long lCount ) {
					std::string	s;

					for ( long i = 0; i < lCount; i++ )
					{
						u_state = (u_state * 6364136223846793005ULL) + 1442695040888963407ULL;
						s += sz_digits[ (u_state >> 33) % l_base ];
					}

					return s;
				};

				std::string	s_top( 1, sz_digits[ l_base - 1 ] );
				std::string	args[] =
				{
					"1." + digits( 5 ), s_top + "." + digits( l_places ),
					"1.00000000000000000" + digits( 20 ), "0." + std::string( 20, sz_digits[ l_base - 1 ] ) + digits( 10 ),
					s_top + digits( 9 ) + "." + digits( 10 ), "0.00000000" + s_top + digits( 15 ),
				};

				for ( auto& s_arg : args )
				{
					glWorkingPrecision = l_places;

					CLZNumber	x, by_agm, by_log;

					x.SetBase( l_base, true );
					x.Assign( s_arg.c_str() );
					lzi_log_by_agm( &x, &by_agm );
					lz_log_e( &x, &by_log );

					char*	sz_agm = by_agm.GetStringRep();
					char*	sz_log = by_log.GetStringRep();

					l_cases++;

					if ( (sz_agm == nullptr) || (sz_log == nullptr) || (::strcmp( sz_agm, sz_log ) != 0) )
						if ( l_wrong++ == 0 )
							s_first = "base " + std::to_string( l_base ) + " at " + std::to_string( l_places ) + ", ln(" +
									  s_arg.substr( 0, 12 ) + "...)";

					::free( sz_agm );
					::free( sz_log );
				}
			}
		}

		char	szDetail[ 200 ];

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld identical, fifteen bases", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, first %s", l_wrong, l_cases, s_first.c_str() );

		CheckTrue( "the AGM and lz_log_e agree", l_wrong == 0, szDetail );
	}

	/*
		**An exact power of the base is j ln(base)**: a multiplication once
		ln(base) is held, where the general path would take an AGM besides.
		It is not what keeps the logarithm from recursing -- ln(base) built
		by the AGM of b^n asks no logarithm of anything, and that is -- so a
		mutation removing it changes no answer, only the time. Against
		lz_log_e, at a width past the table's...
	*/

	{
		long		l_wrong = 0;
		std::string	s_first;

		for ( long l_base : { 2L, 7L, 10L, 16L } )
			for ( const char* sz_arg : { "10", "100000", "0.001" } )
			{
				glWorkingPrecision = CLZNumber::EquivalentPlaces( 160, 10, l_base );

				CLZNumber	x, by_agm, by_log;

				x.SetBase( l_base, true );
				x.Assign( sz_arg );
				lzi_log_by_agm( &x, &by_agm );
				lz_log_e( &x, &by_log );

				if ( by_agm.Compare( &by_log ) != kEqualTo )
					if ( l_wrong++ == 0 )
						s_first = "base " + std::to_string( l_base ) + ", ln(" + sz_arg + ")";
			}

		CheckTrue( "an exact power is j ln(base)", l_wrong == 0, s_first.c_str() );
	}

	/*
		**With no guard, every answer within a unit.** The kernel's own
		margins -- the places for ln s, for m on ln b, the reach's log term,
		the slack in digits -- are covered at eight digits of guard whether
		they are there or not, so only here can a test see them: removing any
		one puts answers 14 to 37,000 units out, where the kernel as written
		stays within one. Against lz_log_e, correctly rounded, in four bases
		at two widths...
	*/

	{
		long		l_savedGuard = lzi_set_log_agm_guard( 0 );
		long		l_wrong = 0;
		std::string	s_first;
		uint64_t	u_state = 1013;

		for ( long l_base : { 2L, 3L, 10L, 16L } )
		{
			const char*	sz_digits = "0123456789ABCDEF";

			for ( long l_places : { CLZNumber::EquivalentPlaces( 300, 10, l_base ), CLZNumber::EquivalentPlaces( 900, 10, l_base ) } )
				for ( int k = 0; k < 8; k++ )
				{
					std::string	s_arg( 1, sz_digits[ 1 + (k % (l_base - 1)) ] );

					s_arg += ".";

					for ( long i = 0; i < l_places; i++ )
					{
						u_state = (u_state * 6364136223846793005ULL) + 1442695040888963407ULL;
						s_arg += sz_digits[ (u_state >> 33) % l_base ];
					}

					glWorkingPrecision = l_places;

					CLZNumber	x, by_agm, by_log, gap, unit;

					x.SetBase( l_base, true );
					x.Assign( s_arg.c_str() );
					lzi_log_by_agm( &x, &by_agm );
					lz_log_e( &x, &by_log );

					// A unit of the last place, b^-W...
					unit.SetBase( l_base, true );
					unit.Assign( ("0." + std::string( l_places - 1, '0' ) + "1").c_str() );
					lz_subtract( &by_agm, &by_log, &gap );

					if ( gap.CompareABS( &unit ) == kGreaterThan )
						if ( l_wrong++ == 0 )
							s_first = "base " + std::to_string( l_base ) + " at " + std::to_string( l_places );
				}
		}

		lzi_set_log_agm_guard( l_savedGuard );

		CheckTrue( "with no guard, within a unit", l_wrong == 0, s_first.c_str() );
	}

	/*
		**A constant computed wide serves a narrower width, rounded.** From an
		empty cache, log_b and the sine at falling widths past their 120
		places, then each width again from an empty cache on its own: the
		answers must be the same, and log_b at a narrower width must take
		fewer logarithms falling than fresh, ln(base) having come from the
		wider one rather than from lz_log_e.
	*/

	{
		auto	answer_at = []( char cFn, long lBase, long lPlaces, long* plLogs ) -> std::string {
			glWorkingPrecision = lPlaces;

			CLZNumber	x( (cFn == 'l') ? "7.25" : "100.5" ), answer;
			long		l_before = lzi_logarithm_count();

			if ( lBase != 10 )
				x.SetBase( lBase );

			if ( cFn == 'l' )
				lz_log_b( &x, &answer );
			else
				lz_sine( &x, &answer );

			if ( plLogs != nullptr )
				*plLogs = lzi_logarithm_count() - l_before;

			char*		sz_got = answer.GetStringRep();
			std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

			::free( sz_got );

			return s_got;
		};

		long		l_wrong = 0;
		bool		b_served = true;
		std::string	s_first;

		for ( char cFn : { 'l', 's' } )
			for ( long l_base : { 7L, 16L } )
			{
				const long	widths[] = { 400, 300, 200 };
				std::string	falling[ 3 ];
				long		l_logs[ 3 ] = { 0, 0, 0 };

				lzu_takedown();
				lzu_setup();

				for ( int i = 0; i < 3; i++ )
					falling[ i ] = answer_at( cFn, l_base, widths[ i ], &l_logs[ i ] );

				for ( int i = 0; i < 3; i++ )
				{
					long	l_fresh = 0;

					lzu_takedown();
					lzu_setup();

					std::string	s_fresh = answer_at( cFn, l_base, widths[ i ], &l_fresh );

					if ( (cFn == 'l') && (i > 0) && !(l_logs[ i ] < l_fresh) )
						b_served = false;

					if ( s_fresh != falling[ i ] )
						if ( l_wrong++ == 0 )
							s_first = std::string( (cFn == 'l') ? "log_b" : "sine" ) + " base " + std::to_string( l_base ) +
									  " at " + std::to_string( widths[ i ] );
				}
			}

		CheckTrue( "a wider constant serves a narrower width", (l_wrong == 0) && b_served,
				   (l_wrong != 0) ? s_first.c_str() : (b_served ? "" : "ln(base) was computed again") );
	}

	/*
		**Where lz_log_e takes the AGM** -- see lzi_log_agm_floor. With the
		measured floors and an empty cache: a width past the constants' text
		goes to Newton first and the AGM second, which builds them; a
		narrower width then takes the AGM at once, its constants served by
		the wider; within 2^-16 of one never; within the text at once; below
		the floor never; and a takedown forgets it all. Each asked by whether
		the AGM's cycle count moved.
	*/

	{
		lzi_set_log_agm_floor( kLogAgmFloorMeasured );
		lzu_takedown();
		lzu_setup();

		auto	took_agm = []( long lBase, long lPlaces, const char* szArg ) -> bool {
			glWorkingPrecision = lPlaces;

			CLZNumber	x, answer;
			long		l_cycles = lzi_agm_log_cycles();

			x.SetBase( lBase, true );
			x.Assign( szArg );
			lz_log_e( &x, &answer );

			return lzi_agm_log_cycles() > l_cycles;
		};

		struct { const char* m_szWhat; long m_lBase; long m_lPlaces; const char* m_szArg; bool m_bAgm; } kSteps[] =
		{
			{ "a new width's first, by Newton",		10, 400, "1.2345",	false },
			{ "its second, by the AGM",				10, 400, "2.75",	true },
			{ "a narrower width's first, held",		10, 300, "1.2345",	true },
			{ "within 2^-16 of one, never",			10, 400, "1.000001",	false },
			{ "within 2^-16 of one, again",			10, 400, "0.999999",	false },
			{ "just past 2^-16 of one",				10, 400, "1.0001",	true },
			{ "within the text, at once",			7,	100, "2.5",		true },
			{ "just past the text, by Newton",		7,	130, "2.5",		false },
			{ "base two, never",					2,	1000, "1.1",	false },
			{ "base two, again",					2,	1000, "1.1",	false },
		};

		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& step : kSteps )
			if ( took_agm( step.m_lBase, step.m_lPlaces, step.m_szArg ) != step.m_bAgm )
				if ( l_wrong++ == 0 )
					s_first = step.m_szWhat;

		// ...and a takedown forgets what was held...
		lzu_takedown();
		lzu_setup();

		if ( took_agm( 10, 300, "1.2345" ) )
			if ( l_wrong++ == 0 )
				s_first = "after a takedown, a first call by Newton";

		CheckTrue( "lz_log_e takes the AGM where it should", l_wrong == 0, s_first.c_str() );
	}

	lzi_set_log_agm_floor( l_savedFloor );
	glWorkingPrecision = l_saved;
}

static void TestRootCore()
{
	Banner( "The square root's core" );

	/*
		**Since 2026-10-14 lz_root_2 takes the AGM's core** -- Newton on the
		reciprocal root, multiplications only -- from a floor per base, a
		short argument keeping Heron's iteration within its base's band; see
		lzi_root_core_floor. Its rounding is decided after the core as
		before, so the two must answer identically: pinned against roots
		worked in integers, the two forced in turn on every shape in every
		base, and the floors and bands taking the core where they say.
	*/

	long	l_saved = glWorkingPrecision;
	long	l_savedFloor = lzi_set_root_core_floor( kRootCoreFloorMeasured );

	// From Python's isqrt, ties to the even kept digit, none within b^-6
	//	of a half; each width past the floor and any band, so the core
	//	takes it...
	struct { long m_lBase; long m_lPlaces; const char* m_szArg; const char* m_szExpect; } kPins[] =
	{
		{ 10, 200, "2",
		  "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623091229702492483605585073721264412149709993583141322266592750559275579995050115278206057147" },
		{ 7, 300, "3",
		  "1.506044021410166456140001531506534532346615201554245563523136600216625003643561251366033554116223663315125543541005620631256443246603564021446324153545063122636644665364620600534200301352664552016550012056364445065004203334620550125022504022021404366445625036511011140220110010356220231015112606635100" },
		{ 16, 300, "A.8",
		  "3.3D88E944C24C7AAFCEDD916BE3023C67B3B431FE1F433D024230EEAA0599BEBD183001326F532DC50ABBB7F0E2D376271556ECC880CB42F30B11797CC838B4BAEFDDDC54A9AE71D9E60A9DE087C7BDB5CC49C4080130C38946F83CBC1D10C97C768A25A4F66DE94D8847866B8703F074874F74175499F662036C5C4D8BE0B03EC105EDE23965DE650479BEED3E60533B068B116AAEAD" },
		{ 12, 300, "2",
		  "1.4E79170D07E85737704E085486853504563650E559E8E79D401387E342380D998D173D951303434821E55419D068816958E64282342D358D8947369E97237E9E04E656D0723349328D219013D8E21DE42844D5758ED27E3D14317E17E28D4354E796260136269D55D79598D4619ED2352D310D3373251E0598676E4537681D191D6901560E13362953D3E37305425159369305141042" },
		{ 5, 800, "3",
		  "1.33122312133132121124103434130232333320033412000444122122443114000241022022401444211421322441222431124400434311014301022302343332242334440232011242110412100100142310410104020232343122404321301144102042402134404204313141104120141232211422200024341130311442324242144100314002331113001413040012140302333113421020221021111203232104403344124042412002121132324242441204314032422303110304413433333323244203221042341201212410113410141430222200003040211130341221240004121321244000224210020310224024014320113211401304234111201114404022412114400430041103444421111423120032331034030434221313134241142033044104241040323344233040231431231402421424103132323304433031230233311310423041101144142111031420012312322404032122003021002121024421341043430322212223422211002411102002303214342013431333324221104042424244142144" },
		{ 3, 2500, "2",
		  "1.1020112212220012122122200112100010200022012111202022012021102210100022112101122102021210020200212012010122121212121220202122020210011211212100102112001210212010220222121211020000100212022221021122021102221102220200001211120210010022211120221001000121222001121222121110101220210221002210100012122121211011012000011122222211000010201221102010100002020012201102210000112212001111221210202210111121200020210202201022201110002201021001112110012020221011121022022011120122200010021001100022020000210201201020102022101220000200002210001022012200202211001001021120220111111010002022110112020000102010121021111112122200001000011222000001012210221220021210121000122111012200110212002222020020222011211212122021001212211020100010212022010022022012212022011121210000000112110101222112110020221002020202102020120102220102201222102022010220020221020200101111221110001212012022122012121221200101000020201112111200220022120120111222012020011221122121010221110122001022121002222120021112222102201120202022211120120222102111221220121021012202112210000110101102111112110210020202202122011002110011211000111110222220101121002220211021100121101011221022210011010020021022210110020200010011022221202211110122222100001020011111222010112020202010102220211012202110111200220201010012121220100212022211012110020020100200200122020200022001201210111220001222200200101121212220122100112022110211112120012111000021001112111012211121102222120121100221011221000122002212222101220210122022120220121110121121120122121001121201002101012221021022011001112200111011010122000020001122121222120210212100212000220221022202002111200202110010121002221111012121102210021120000102200222021202212101101011121220122112022220100100120001212110201120220200102202001021120210102200201110021122102110102111211010221001100112020211011122022011010202210100211010022120100012011211011220222112211110220002000000201100011100221011200010121221001101101021022222212120011112012121022222021020222211002201122210110010110120202002212101011102201012121101012212011110101211011011011012011112200000022012100110000012011101000211211222202210012022220020210012100211102012212200202101010121100100100000112220021020102111202211022101120221202110102200101210110201110201011121120122201000212222111102011200100101021211212200111121102102012012011222201201122222222101110111122221201210211110222021112110002002000110001110102211111110112200122010111120201010212121020211212220110220202211012202020102212200202001122202112210000122111102020101100010201000211210000101" },
		{ 9, 150, "0.000005",
		  "0.002211075173246240740776355826238034615034608471457760625007406766575676064504384168501501166586010075853057652320085530443866010150884771607747556438" },
	};

	{
		long		l_wrong = 0;
		long		l_runs = lzi_root_core_runs();
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			glWorkingPrecision = pin.m_lPlaces;

			CLZNumber	x, answer;

			x.SetBase( pin.m_lBase, true );
			x.Assign( pin.m_szArg );
			lz_root_2( &x, &answer );

			char*	sz_got = answer.GetStringRep();

			if ( (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
				if ( l_wrong++ == 0 )
					s_first = std::string( "base " ) + std::to_string( pin.m_lBase ) + " sqrt(" + pin.m_szArg + "): " +
							  ((sz_got != nullptr) ? std::string( sz_got ).substr( 0, 30 ) : "(null)");

			::free( sz_got );
		}

		long	l_took = lzi_root_core_runs() - l_runs;
		char	szDetail[ 200 ];
		size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character, the core taking %ld", n_pins, l_took );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong, first %s", l_wrong, n_pins, s_first.c_str() );

		CheckTrue( "against integers", (l_wrong == 0) && (l_took == static_cast<long>( n_pins )), szDetail );
	}

	/*
		**The core and Heron agree, every base**, forced in turn: five widths
		from 5 to 1,500 places and six shapes -- a short argument, one as long
		as the width, a perfect square, which both must spell exact, a tiny
		and a huge one, and in an even base an exact tie, the square of a
		value a half-unit past the width.
	*/

	{
		long		l_cases = 0;
		long		l_wrong = 0;
		std::string	s_first;
		uint64_t	u_state = 20261014;

		for ( long l_base = 2; l_base <= 16; l_base++ )
		{
			const char*	sz_digits = (l_base == 12) ? "0123456789DE" : "0123456789ABCDEF";

			for ( long l_places : { 5L, 36L, 100L, 400L, 1500L } )
			{
				auto	digits = [&]( long lCount ) {
					std::string	s;

					for ( long i = 0; i < lCount; i++ )
					{
						u_state = (u_state * 6364136223846793005ULL) + 1442695040888963407ULL;
						s += sz_digits[ (u_state >> 33) % l_base ];
					}

					return s;
				};

				std::string	s_top( 1, sz_digits[ l_base - 1 ] );
				std::vector<std::string>	args =
				{
					"2." + digits( 10 ), s_top + "." + digits( l_places ),
					"0.000000" + s_top + digits( 20 ), s_top + digits( 40 ) + "." + digits( 5 ),
				};

				// A perfect square, and in an even base an exact tie: y^2, and
				//	(y + half a unit past the width)^2, both worked here...
				for ( int k = 0; k < 2; k++ )
				{
					if ( (k == 1) && ((l_base % 2) != 0) )
						continue;

					CLZNumber	y, square, half;
					std::string	s_y = s_top + "." + digits( (k == 0) ? (l_places / 2) : l_places );

					glWorkingPrecision = 4 * l_places + 40;
					y.SetBase( l_base, true );
					y.Assign( s_y.c_str() );

					if ( k == 1 )
					{
						half.SetBase( l_base, true );
						half.Assign( ("0." + std::string( l_places, '0' ) + std::string( 1, sz_digits[ l_base / 2 ] )).c_str() );
						lz_add( &y, &half );
					}

					lz_multiply( &y, &y, &square );

					char*	sz_square = square.GetStringRep();

					args.push_back( sz_square );
					::free( sz_square );
				}

				for ( auto& s_arg : args )
				{
					std::string	s_answer[ 2 ];

					for ( int i = 0; i < 2; i++ )
					{
						lzi_set_root_core_floor( (i == 0) ? 0 : LONG_MAX );
						glWorkingPrecision = l_places;

						CLZNumber	x, answer;

						x.SetBase( l_base, true );
						x.Assign( s_arg.c_str() );
						lz_root_2( &x, &answer );

						char*	sz_got = answer.GetStringRep();

						s_answer[ i ] = (sz_got != nullptr) ? sz_got : "(null)";
						::free( sz_got );
					}

					l_cases++;

					if ( s_answer[ 0 ] != s_answer[ 1 ] )
						if ( l_wrong++ == 0 )
							s_first = "base " + std::to_string( l_base ) + " at " + std::to_string( l_places ) + ", sqrt(" +
									  s_arg.substr( 0, 12 ) + "...)";
				}
			}
		}

		lzi_set_root_core_floor( kRootCoreFloorMeasured );

		char	szDetail[ 200 ];

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld identical, fifteen bases", l_cases );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, first %s", l_wrong, l_cases, s_first.c_str() );

		CheckTrue( "the core and Heron agree", l_wrong == 0, szDetail );
	}

	/*
		**The floors and the bands take the core where they say**, each asked
		by whether the core's count moved.
	*/

	{
		auto	took_core = []( long lBase, long lPlaces, const char* szArg ) -> bool {
			glWorkingPrecision = lPlaces;

			CLZNumber	x, answer;
			long		l_runs = lzi_root_core_runs();

			x.SetBase( lBase, true );
			x.Assign( szArg );
			lz_root_2( &x, &answer );

			return lzi_root_core_runs() > l_runs;
		};

		std::string	s_long3 = "2." + std::string( 499, '1' );
		std::string	s_padded3 = "2." + std::string( 499, '0' );

		struct { const char* m_szWhat; long m_lBase; long m_lPlaces; std::string m_sArg; bool m_bCore; } kSteps[] =
		{
			{ "the default width, base ten",		10, 36,		"2",		true },
			{ "base two, short of its floor",		2,	1000,	"1.1",		false },
			{ "base four, short of its floor",		4,	5,		"2",		false },
			{ "base three, a short argument in its band",	3,	500,	"2",		false },
			{ "base three, a long one in it",		3,	500,	s_long3,	true },
			{ "base three, 2 held with its noughts",	3,	500,	s_padded3,	true },
			{ "base three, short, past its band",	3,	2500,	"2",		true },
			{ "base twelve, short, in its band",	12, 200,	"2",		false },
			{ "base twelve, short, past it",		12, 300,	"2",		true },
			{ "base seven, short, no band",			7,	500,	"2",		true },
		};

		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& step : kSteps )
			if ( took_core( step.m_lBase, step.m_lPlaces, step.m_sArg.c_str() ) != step.m_bCore )
				if ( l_wrong++ == 0 )
					s_first = step.m_szWhat;

		CheckTrue( "the floors and bands take the core where they say", l_wrong == 0, s_first.c_str() );
	}

	lzi_set_root_core_floor( l_savedFloor );
	glWorkingPrecision = l_saved;
}

static void TestNearTiesDecided()
{
	Banner( "Near ties the retake reaches" );

	/*
		**Since 2026-10-14 a retake widens until it decides**, the extra
		doubling each round to the answer's own places -- see
		kRetakeGuardDigits. A tiny argument in an odd base puts e^x, cos,
		cosh and ln(1 + x) a hair from a half -- x^2 / 2 is a half that never
		ends there, and only the next term moves it -- deeper than the old
		three rounds of five digits could see. Each pin here, from the series
		summed exactly, was misrounded until then. And an answer exactly on a
		half, which reaches the cap, still goes to the even kept digit.
	*/

	long	l_saved = glWorkingPrecision;

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	struct { const char* m_szFn; long m_lBase; long m_lPlaces; const char* m_szArg; const char* m_szExpect; } kPins[] =
	{
		{ "cos", 3, 100, "0.0000000000000000000000000000000000000000001",
		  "0.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222211111111111112" },
		{ "cos", 3, 300, "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001",
		  "0.222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222221111111111111111111111111111111111111111111112" },
		{ "cos", 7, 36, "0.00000000000000001",
		  "0.666666666666666666666666666666666634" },
		{ "cosh", 5, 36, "0.000000000000000001",
		  "1.000000000000000000000000000000000001" },
		{ "cosh", 5, 100, "0.000000000000000000000000000000000000003",
		  "1.0000000000000000000000000000000000000000000000000000000000000000000000000000042222222222222222222223" },
		{ "cosh", 5, 300, "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001",
		  "1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002222222222222222222222222222222222222223" },
		{ "exp", 5, 100, "0.00000000000000000000000000000000000000000001",
		  "1.0000000000000000000000000000000000000000000100000000000000000000000000000000000000000000222222222223" },
		{ "ln", 5, 300, "1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001",
		  "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000044444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444442222222222222222222222222222222222222222222222222222222222222222222223" },
	};

	{
		long		l_wrong = 0;
		std::string	s_first;

		for ( auto& pin : kPins )
		{
			Unary	fn = (::strcmp( pin.m_szFn, "cos" ) == 0) ? lz_cosine :
						 (::strcmp( pin.m_szFn, "cosh" ) == 0) ? lz_hypercosine :
						 (::strcmp( pin.m_szFn, "exp" ) == 0) ? lz_exp_e : lz_log_e;

			glWorkingPrecision = pin.m_lPlaces;

			CLZNumber	x, answer;

			x.SetBase( pin.m_lBase, true );
			x.Assign( pin.m_szArg );
			fn( &x, &answer );

			char*	sz_got = answer.GetStringRep();

			if ( (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
				if ( l_wrong++ == 0 )
					s_first = std::string( pin.m_szFn ) + " base " + std::to_string( pin.m_lBase ) + " at " +
							  std::to_string( pin.m_lPlaces );

			::free( sz_got );
		}

		char	szDetail[ 200 ];
		size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character, in odd bases at 36 to 300 places", n_pins );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong, first %s", l_wrong, n_pins, s_first.c_str() );

		CheckTrue( "deep near ties from tiny arguments", l_wrong == 0, szDetail );
	}

	/*
		**An answer exactly on a half reaches the cap and goes to the even
		kept digit**: the geometric mean of y^2 and 1, y written to 37 places
		ending in a half, at 36. Its kept digit even, it stays; odd, it rises.
	*/

	{
		long		l_wrong = 0;
		std::string	s_first;

		struct { const char* m_szY; const char* m_szExpect; } kTies[] =
		{
			{ "1.0000000000000000000000000000000000005", "1.000000000000000000000000000000000000" },
			{ "1.0000000000000000000000000000000000015", "1.000000000000000000000000000000000002" },
		};

		for ( auto& tie : kTies )
		{
			CLZNumber	y( tie.m_szY ), square, one( "1" ), answer;

			glWorkingPrecision = 200;
			lz_multiply( &y, &y, &square );
			glWorkingPrecision = 36;

			CLZNumber*	values[] = { &square, &one };

			lz_gmean( values, 2, &answer );

			char*	sz_got = answer.GetStringRep();

			if ( (sz_got == nullptr) || (::strcmp( sz_got, tie.m_szExpect ) != 0) )
				if ( l_wrong++ == 0 )
					s_first = std::string( "y = " ) + tie.m_szY + ": " + ((sz_got != nullptr) ? sz_got : "(null)");

			::free( sz_got );
		}

		CheckTrue( "an exact half still goes to the even digit", l_wrong == 0, s_first.c_str() );
	}

	glWorkingPrecision = l_saved;
}

static void TestRootsOfWholeOrder()
{
	Banner( "Roots of whole order, decided exactly" );

	/*
		**Since 2026-10-14 a root of whole order in doubt is decided exactly**
		-- root_n of a whole order, and power and exp_n at any exponent P/Q,
		as the root of order Q of x^P -- by the comparison lzi_root_decide
		makes for the square and cube roots; see lzi_whole_order_root_decide.
		An exact tie in an odd base, 8^(-1/3) among them, lies on a half that
		never ends and was rounded whichever way the logarithm's error fell;
		an operand a hair from a tie's lay deeper than any retake reaches.
		Each pin is worked exactly in whole numbers. Twenty-three are ties, each a coin flip without this, and
		eighteen near misses: a side of each order's comparison; exponents of
		P/Q with P past one, 0.3 read as three tenths and not a third, in
		bases ten, seven, six and three; and four whose decider's numbers pass
		the early limit, decided only at the cap. f9eff0f misrounded nineteen of
		the forty-one, 2d977a7 twenty, a0200f9 six.
	*/

	long	l_saved = glWorkingPrecision;

	struct { const char* m_szFn; long m_lBase; long m_lPlaces; const char* m_szX; const char* m_szN; const char* m_szExpect; } kPins[] =
	{
		{ "root_n", 3, 36, "22", "-10",
		  "0.111111111111111111111111111111111112" },
		{ "root_n", 5, 36, "-112", "-10",
		  "-0.222222222222222222222222222222222222" },
		{ "power", 9, 36, "8", "-0.3",
		  "0.444444444444444444444444444444444444" },
		{ "exp_n", 15, 36, "2486", "-0.3",
		  "0.277777777777777777777777777777777778" },
		{ "root_n", 7, 100, "220000", "-4",
		  "0.0333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333334" },
		{ "root_n", 11, 10, "0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000A45122A32619A9737A17461063A2384A4A2001", "-100",
		  "55.5555555555" },
		{ "root_n", 10, 36, "0.001881676372353657772546716040595286742159524921197272014318240193401713453448510464472389192840117714871312124", "3",
		  "0.123456789012345678901234567890123456" },
		{ "root_n", 10, 36, "0.001881676372353657772546716040595286742159524921197272014318240193401713453448510464472389192840117714871312125", "3",
		  "0.123456789012345678901234567890123456" },
		{ "root_n", 10, 36, "0.001881676372353657772546716040595286742159524921197272014318240193401713453448510464472389192840117714871312126", "3",
		  "0.123456789012345678901234567890123457" },
		{ "root_n", 3, 36, "22.000000000000000000000000000000000000000000000000000000000001", "-10",
		  "0.111111111111111111111111111111111111" },
		{ "root_n", 3, 36, "21.222222222222222222222222222222222222222222222222222222222222", "-10",
		  "0.111111111111111111111111111111111112" },
		{ "power", 9, 36, "7.88888888888888888888888888888888888888888888888888", "-0.3",
		  "0.444444444444444444444444444444444445" },
		{ "root_n", 7, 36, "0.00100000000000000000000000000000000001333333333333333333333333333333333341515151515151515151515151515151515160606060606060606060606060606060606060606", "3",
		  "0.100000000000000000000000000000000000" },
		{ "root_n", 3, 36, "0.0121", "-11",
		  "1.111111111111111111111111111111111112" },
		{ "root_n", 5, 36, "0.000224", "-11",
		  "2.222222222222222222222222222222222222" },
		{ "root_n", 7, 36, "4", "-2",
		  "0.333333333333333333333333333333333334" },
		{ "root_n", 9, 36, "40", "-2",
		  "0.144444444444444444444444444444444444" },
		{ "root_n", 11, 36, "0.000059", "-6",
		  "5.555555555555555555555555555555555556" },
		{ "root_n", 13, 36, "4C000000", "-6",
		  "0.066666666666666666666666666666666666" },
		{ "root_n", 15, 36, "8000", "-3",
		  "0.077777777777777777777777777777777778" },
		{ "power", 3, 36, "22", "-0.1",
		  "0.111111111111111111111111111111111112" },
		{ "exp_n", 3, 36, "22000", "-0.1",
		  "0.011111111111111111111111111111111112" },
		{ "power", 9, 36, "8", "-0.3",
		  "0.444444444444444444444444444444444444" },
		{ "exp_n", 9, 36, "8000", "-0.3",
		  "0.044444444444444444444444444444444444" },
		{ "power", 15, 36, "22", "-0.3",
		  "0.777777777777777777777777777777777778" },
		{ "exp_n", 15, 36, "2200000", "-0.3",
		  "0.077777777777777777777777777777777778" },
		{ "power", 5, 36, "112", "-0.1",
		  "0.222222222222222222222222222222222222" },
		{ "exp_n", 5, 36, "11200000", "-0.1",
		  "0.022222222222222222222222222222222222" },
		{ "power", 15, 36, "8", "-0.5",
		  "0.777777777777777777777777777777777778" },
		{ "exp_n", 15, 36, "8000", "-0.5",
		  "0.077777777777777777777777777777777778" },
		{ "power", 10, 20, "0.00093695037791694316034380021827590251", "0.3",
		  "0.12345678901234567891" },
		{ "exp_n", 10, 20, "0.00093695037791694316034380021827590251", "0.3",
		  "0.12345678901234567891" },
		{ "power", 10, 20, "0.000936950377916943160343800218275902500446869234135110538715509081490859169", "0.3",
		  "0.12345678901234567891" },
		{ "exp_n", 10, 20, "0.000936950377916943160343800218275902500446869234135110538715509081490859169", "0.3",
		  "0.12345678901234567891" },
		{ "power", 7, 36, "1.0336123023061410310641144235352256550122636046640014251232225005551655662012310605321415144640201605046266350455032501241250515616400513224015333154", "-0.3",
		  "0.653600016121256612631440102161506524" },
		{ "power", 7, 36, "1.033612302306141031064114423535225655012263604664001425123222500555165566201231060532141514464020160504626635045503250124125051561640051322401533315366", "-0.3",
		  "0.653600016121256612631440102161506525" },
		{ "power", 6, 30, "0.4524220121342345453213443343202014424324422230245302215245350040232454532412132152055413205055421141115352453151523524135410005105", "0.5",
		  "0.502424111144322551244533244521" },
		{ "power", 10, 10, "0.12320459847603969012527207274209438033352315995552759066368667349438814944635039", "0.9990234375",
		  "0.1234567891" },
		{ "power", 10, 10, "0.12320459847603969012527207274209438033352315995552759066368667349438814944635038", "0.9990234375",
		  "0.1234567890" },
		{ "exp_n", 3, 10, "10.00010000222220010201022012002211201020000021011121010100111121022112112122111211", "-0.22222",
		  "0.1000000022" },
		{ "exp_n", 3, 10, "10.00010000222220010201022012002211201020000021011121010100111121022112112122111212", "-0.22222",
		  "0.1000000021" },
	};

	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pin : kPins )
	{
		glWorkingPrecision = pin.m_lPlaces;

		CLZNumber	x, n, answer;

		x.SetBase( pin.m_lBase, true );
		x.Assign( pin.m_szX );
		n.SetBase( pin.m_lBase, true );
		n.Assign( pin.m_szN );

		long	l_code = (::strcmp( pin.m_szFn, "root_n" ) == 0) ? lz_root_n( &x, &n, &answer ) :
						 (::strcmp( pin.m_szFn, "power" ) == 0) ? lz_power( &x, &n, &answer ) : lz_exp_n( &x, &n, &answer );
		char*	sz_got = answer.GetStringRep();

		if ( (l_code != kNoError) || (sz_got == nullptr) || (::strcmp( sz_got, pin.m_szExpect ) != 0) )
			if ( l_wrong++ == 0 )
				s_first = std::string( pin.m_szFn ) + " base " + std::to_string( pin.m_lBase ) + " at " +
						  std::to_string( pin.m_lPlaces ) + ", order " + pin.m_szN;

		::free( sz_got );
	}

	glWorkingPrecision = l_saved;

	char	szDetail[ 200 ];
	size_t	n_pins = sizeof( kPins ) / sizeof( kPins[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu exact to the character: ties, a negative's, near misses", n_pins );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu wrong, first %s", l_wrong, n_pins, s_first.c_str() );

	CheckTrue( "ties and near misses of whole order", l_wrong == 0, szDetail );
}

static void TestSignedComparisons()
{
	Banner( "Comparing signed numbers" );

	/*
		**Two negatives compared the wrong way round until 2026-10-15**:
		CLZNumber::Compare returned the magnitudes' order whenever the signs
		matched, so lz_compare( -3, -5 ) was -1 and lz_max picked -5. No test
		had compared two negatives. Each pair here is compared both ways
		round, by lz_compare, lz_max and lz_min, in four bases, against the
		sign of the difference -- subtraction's, an independent road.
	*/

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = 36;

	const char*	kPairs[][ 2 ] =
	{
		{ "-3", "-5" }, { "-1", "-0.8" }, { "-0.25", "-0.5" }, { "-10000000000", "-0.0000000001" },
		{ "-7.5", "-7.5" }, { "-2", "3" }, { "0", "-4" }, { "5", "-5" }, { "-0.000001", "0" },
		{ "2", "9" }, { "-123.456", "-123.457" }, { "0", "0" },
	};
	const long	kBases[] = { 10, 16, 7, 2 };

	long		l_wrong = 0;
	long		l_checked = 0;
	std::string	s_first;

	for ( long l_base : kBases )
		for ( auto& pair : kPairs )
			for ( int i_way = 0; i_way < 2; i_way++ )
			{
				CLZNumber	x( pair[ i_way ] );
				CLZNumber	y( pair[ 1 - i_way ] );
				CLZNumber	difference;

				x.SetBase( l_base, true );
				y.SetBase( l_base, true );
				lz_subtract( &x, &y, &difference );

				long	l_want = difference.IsZero() ? kEqualTo : (difference.IsNegative() ? kLessThan : kGreaterThan);
				long	l_compare = lz_compare( &x, &y );
				long	l_max = lz_max( &x, &y );
				long	l_min = lz_min( &x, &y );
				long	l_wantMax = (l_want == kGreaterThan) ? kLeft : ((l_want == kLessThan) ? kRight : kEqual);
				long	l_wantMin = (l_want == kLessThan) ? kLeft : ((l_want == kGreaterThan) ? kRight : kEqual);

				l_checked++;

				if ( (l_compare != l_want) || (l_max != l_wantMax) || (l_min != l_wantMin) )
					if ( l_wrong++ == 0 )
						s_first = std::string( "base " ) + std::to_string( l_base ) + ": " + pair[ i_way ] + " against " +
								  pair[ 1 - i_way ] + ", compare " + std::to_string( l_compare );
			}

	glWorkingPrecision = l_saved;

	char	szDetail[ 200 ];

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %ld right: compare, max and min, both ways round", l_checked );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

	CheckTrue( "two negatives, and every other sign pairing", l_wrong == 0, szDetail );
}

static void TestNothingLeftBehind()
{
	Banner( "Nothing left behind" );

	/*
		**lz_log_e leaked a number on every call down its series path until
		2026-10-15** -- base_e was the one working number its cleanup did not
		delete -- and so did everything built on it. Each function here is
		called once to fill whatever caches it keeps, then 200 times more,
		and the count of live allocations must end where it began.
	*/

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = 36;

	struct { const char* m_szName; long m_lBase; const char* m_szX; const char* m_szY; } kCalls[] =
	{
		{ "log_e", 10, "2.5", nullptr }, { "log_e", 7, "3.2", nullptr }, { "power", 10, "2.5", "0.3" },
		{ "exp_n", 10, "2.5", "1.7" }, { "root_n", 10, "2.5", "7" }, { "log_n", 10, "3", "2.5" },
		{ "archypersine", 10, "2.5", nullptr }, { "gamma", 10, "2.5", nullptr },
	};

	long		l_leaking = 0;
	std::string	s_first;

	for ( auto& c : kCalls )
	{
		CLZNumber	x;
		CLZNumber	y;

		x.SetBase( c.m_lBase, true );
		x.Assign( c.m_szX );

		if ( c.m_szY != nullptr )
		{
			y.SetBase( c.m_lBase, true );
			y.Assign( c.m_szY );
		}

		auto	call = [ & ]()
		{
			CLZNumber	r;

			if ( !::strcmp( c.m_szName, "log_e" ) ) lz_log_e( &x, &r );
			else if ( !::strcmp( c.m_szName, "power" ) ) lz_power( &x, &y, &r );
			else if ( !::strcmp( c.m_szName, "exp_n" ) ) lz_exp_n( &x, &y, &r );
			else if ( !::strcmp( c.m_szName, "root_n" ) ) lz_root_n( &x, &y, &r );
			else if ( !::strcmp( c.m_szName, "log_n" ) ) lz_log_n( &x, &y, &r );
			else if ( !::strcmp( c.m_szName, "archypersine" ) ) lz_archypersine( &x, &r );
			else lz_gamma( &x, &r );
		};

		call();

		long	l_before = g_liveNews;

		for ( int i = 0; i < 200; i++ )
			call();

		long	l_left = g_liveNews - l_before;

		if ( l_left != 0 )
			if ( l_leaking++ == 0 )
				s_first = std::string( c.m_szName ) + " in base " + std::to_string( c.m_lBase ) + " left " +
						  std::to_string( l_left ) + " behind in 200 calls";
	}

	glWorkingPrecision = l_saved;

	char	szDetail[ 200 ];
	size_t	n_calls = sizeof( kCalls ) / sizeof( kCalls[ 0 ] );

	if ( l_leaking == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "%zu functions, 200 calls each, every allocation handed back", n_calls );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu leaking, first %s", l_leaking, n_calls, s_first.c_str() );

	CheckTrue( "the logarithm and what is built on it", l_leaking == 0, szDetail );
}

static void TestCubeRoundedOnce()
{
	Banner( "The cube, rounded once" );

	/*
		**lz_cube rounded twice until 2026-10-15** -- the square to the
		width, then the product -- so the square's error came back
		multiplied by the operand: 42 units of the 36th place for
		123.4567..., and wrong in about two operands in three with more
		places than the width. Each operand here, written past the width in
		six bases, is checked against its cube worked exactly and rounded
		once, half-even by the kept digit. The generator's state passes 2^31,
		so it is unsigned long long -- see "A 32-bit long".
	*/

	long				l_saved = glWorkingPrecision;
	unsigned long long	ull_state = 20261015ULL;
	auto				next = [ & ]( unsigned long long n ) -> unsigned long long
	{
		ull_state = (ull_state * 6364136223846793005ULL) + 1442695040888963407ULL;
		return (ull_state >> 33) % n;
	};

	const long	kBases[] = { 10, 16, 7, 3, 2, 12 };
	const long	kWidths[] = { 1, 2, 10, 36, 100 };
	const char*	kDigits = "0123456789ABCDEF";

	long		l_wrong = 0;
	long		l_checked = 0;
	std::string	s_first;

	for ( long l_base : kBases )
		for ( long l_width : kWidths )
			for ( int i = 0; i < 12; i++ )
			{
				// An operand of up to five whole digits and up to the width
				//	plus eight places, either sign...
				std::string	s_text = next( 2 ) ? "-" : "";
				long		l_whole = 1 + static_cast<long>( next( 5 ) );
				long		l_places = 1 + static_cast<long>( next( static_cast<unsigned long long>( l_width + 8 ) ) );

				for ( long k = 0; k < l_whole; k++ )
					s_text += kDigits[ next( static_cast<unsigned long long>( l_base ) ) ];

				s_text += ".";

				for ( long k = 0; k < l_places; k++ )
					s_text += kDigits[ next( static_cast<unsigned long long>( l_base ) ) ];

				if ( l_base == 12 )
					for ( char& c : s_text )
						c = (c == 'A') ? 'D' : ((c == 'B') ? 'E' : c);

				CLZNumber	x;
				CLZNumber	got;
				CLZNumber	exact;

				x.SetBase( l_base, true );
				x.Assign( s_text.c_str() );

				glWorkingPrecision = l_width;
				lz_cube( &x, &got );

				// The cube exactly -- a product's places are its operands'
				//	together -- then rounded once...
				glWorkingPrecision = (3 * x.m_lRadix) + 4;
				lz_multiply( &x, &x, &exact );
				lz_multiply( &exact, &x );
				exact.Round( RM_round_half_even, l_width );

				l_checked++;

				if ( lz_compare( &got, &exact ) != kEqualTo )
					if ( l_wrong++ == 0 )
						s_first = std::string( "base " ) + std::to_string( l_base ) + " at " + std::to_string( l_width ) + ": " + s_text;
			}

	glWorkingPrecision = l_saved;

	char	szDetail[ 200 ];

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %ld operands, six bases, widths 1 to 100", l_checked );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

	CheckTrue( "operands written past the width", l_wrong == 0, szDetail );
}

// One mixed pair worked through any operation the table names; a comparison
//	answers its order as text, as the table holds it...
static std::string MixedPairAnswer( const char* szOp, long lWidth, long lLeftBase, const char* szLeft,
									long lRightBase, const char* szRight, long* plResult )
{
	long		l_saved = glWorkingPrecision;

	// The width set before the numbers are made: a number takes its
	//	precision when it is made and assigned, and an answer is written
	//	to its operand's...
	glWorkingPrecision = lWidth;

	CLZNumber	x;
	CLZNumber	y;
	CLZNumber	got;
	CLZNumber*	p_values[ 2 ] = { &x, &y };

	SetInBase( x, lLeftBase, szLeft );
	SetInBase( y, lRightBase, szRight );

	std::string	s_op = szOp;
	long		l_result = kNoError;
	long		l_order = 0;
	bool		b_order = true;

	if ( s_op == "compare" )			l_order = lz_compare( &x, &y );
	else if ( s_op == "max" )			l_order = lz_max( &x, &y );
	else if ( s_op == "min" )			l_order = lz_min( &x, &y );
	else
	{
		b_order = false;

		if ( s_op == "add" )					l_result = lz_add( &x, &y, &got );
		else if ( s_op == "subtract" )			l_result = lz_subtract( &x, &y, &got );
		else if ( s_op == "multiply" )			l_result = lz_multiply( &x, &y, &got );
		else if ( s_op == "divide" )			l_result = lz_divide( &x, &y, &got );
		else if ( s_op == "divide_integer" )	l_result = lz_divide_integer( &x, &y, &got );
		else if ( s_op == "remainder" )			l_result = lz_remainder( &x, &y, &got );
		else if ( s_op == "remainder_near" )	l_result = lz_remainder_near( &x, &y, &got );
		else if ( s_op == "modulo" )			l_result = lz_modulo( &x, &y, &got );
		else if ( s_op == "power" )				l_result = lz_power( &x, &y, &got );
		else if ( s_op == "log_n" )				l_result = lz_log_n( &x, &y, &got );
		else if ( s_op == "amean" )				l_result = lz_amean( p_values, 2, &got );
		else if ( s_op == "gmean" )				l_result = lz_gmean( p_values, 2, &got );
		else									l_result = lz_hmean( p_values, 2, &got );
	}

	glWorkingPrecision = l_saved;

	if ( plResult != nullptr )
		*plResult = l_result;

	if ( b_order )
		return std::to_string( l_order );

	char*		sz_got = got.GetStringRep();
	std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

	::free( sz_got );

	return s_got;
}

/*
	**A base whose primes all divide the left's ends in it, and converts
	exactly** -- since 2026-10-21. lzi_ends_in_base tested every divisor
	of the right operand's base, so 22 pairs whose numbers do end read as
	not ending, and a mixed pair in them was worked over a power of its
	base and written as a quotient. So the same value gave two answers by
	how it was written: 0.111... times a quarter written in base four was
	not what it was with the quarter written in the left base itself.
	Each row is 1/f, written "0.1" in base f, and the same value written
	exactly in base i, generated with Fraction; the product and the sum
	must not tell them apart. The last six pairs always ended.
*/

static void TestEndingBases()
{
	Banner( "a base whose primes divide the left's" );

	long	l_saved = glWorkingPrecision;

	glWorkingPrecision = 36;

	struct { long m_lFrom; long m_lInto; const char* m_szInInto; } kPairs[] =
	{
		{  4,  2, "0.01" },
		{  4,  6, "0.13" },
		{  4, 10, "0.25" },
		{  4, 14, "0.37" },
		{  8,  2, "0.001" },
		{  8,  4, "0.02" },
		{  8,  6, "0.043" },
		{  8, 10, "0.125" },
		{  8, 12, "0.16" },
		{  8, 14, "0.1A7" },
		{  9,  3, "0.01" },
		{  9,  6, "0.04" },
		{  9, 12, "0.14" },
		{  9, 15, "0.1A" },
		{ 12,  6, "0.03" },
		{ 16,  2, "0.0001" },
		{ 16,  4, "0.01" },
		{ 16,  6, "0.0213" },
		{ 16,  8, "0.04" },
		{ 16, 10, "0.0625" },
		{ 16, 12, "0.09" },
		{ 16, 14, "0.0C37" },
		{  2, 10, "0.5" },
		{  5, 10, "0.2" },
		{  3,  6, "0.2" },
		{  2, 16, "0.8" },
		{  4,  8, "0.2" },
		{  7, 14, "0.2" },
	};

	long		l_wrong = 0;
	std::string	s_first;

	for ( auto& pair : kPairs )
	{
		// A value filling the width in the left base, so a quotient's
		//	rounding would show...
		std::string	s_left = "0." + std::string( 36, '1' );
		CLZNumber	x;
		CLZNumber	y;
		CLZNumber	yInto;

		x.SetBase( pair.m_lInto, false );
		x.Assign( s_left.c_str() );
		y.SetBase( pair.m_lFrom, false );
		y.Assign( "0.1" );
		yInto.SetBase( pair.m_lInto, false );
		yInto.Assign( pair.m_szInInto );

		CLZNumber	mixed;
		CLZNumber	same;
		bool		b_same = true;

		lz_multiply( &x, &y, &mixed );
		lz_multiply( &x, &yInto, &same );
		b_same = SameAnswer( &mixed, &same );

		lz_add( &x, &y, &mixed );
		lz_add( &x, &yInto, &same );
		b_same = b_same && SameAnswer( &mixed, &same );

		if ( !b_same && (l_wrong++ == 0) )
			s_first = std::to_string( pair.m_lFrom ) + " into " + std::to_string( pair.m_lInto );
	}

	char	szDetail[ 120 ];
	size_t	n_pairs = sizeof( kPairs ) / sizeof( kPairs[ 0 ] );

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %zu pairs answer as the value written in the left base", n_pairs );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %zu not, first %s", l_wrong, n_pairs, s_first.c_str() );

	CheckTrue( "a product and a sum, written in either base", l_wrong == 0, szDetail );

	glWorkingPrecision = l_saved;
}

static void TestMixedPairsExact()
{
	Banner( "Mixed pairs, worked exactly" );

	/*
		**A mixed pair converted its right operand at a fixed count of the
		left base's places until 2026-10-15** and then treated it as exact,
		so a small right operand lost its relative precision and two values
		differing past those places compared wrongly. Power and log_n went
		further and worked partly in the right operand's base at its width.
		These are the cases the review found; the table after them was
		generated, every answer scored against the exact value or a Decimal
		reference at the places it is written, and two in three of each
		operation's eight are cases the old code got wrong -- see
		CLZMixedPair and lzi_into_base_keeping.

		Four were generated again on 2026-10-21, when lzi_ends_in_base
		began testing primes rather than every divisor: 16 into 10, 12
		into 6 twice and 16 into 14 convert exactly now, so the sum, the
		difference and the two remainders are exact, as a same-base pair's
		are, where they were quotients rounded at the width. Each is the
		exact value, from Fraction.
	*/

	const char*	kTwoToMinus34 = "0.0000000000000000000000000000000001";
	const char*	kTwoToMinus50 = "0.00000000000000000000000000000000000000000000000001";

	struct { const char* szLabel; const char* szOp; long lWidth; long lLeftBase; const char* szLeft;
			 long lRightBase; const char* szRight; const char* szExpected; } named[] =
	{
		{ "1 over 16^-7, short and exact",  "divide",  5, 10, "1", 16, "0.0000001", "268435456" },
		{ "1 over 2^-34 is no nought",      "divide", 10, 10, "1",  2, kTwoToMinus34, "17179869184" },
		{ "0.3333333333333 below 1/3",      "compare", 10, 10, "0.3333333333333", 3, "0.1", "-1" },
		{ "10^(1/3), 36 digits",            "power",  36, 10, "10", 3, "0.1", "2.154434690031883721759293566519350495" },
		{ "log10 of 1/3, 36 digits",        "log_n",  36, 10, "10", 3, "0.1", "-0.477121254719662437295027903255115309" },
		{ "gmean of 1 and 2^-50",           "gmean",  10, 10, "1",  2, kTwoToMinus50, "0.0000000298" },
		{ "amean of 1 and 2^-50, not short", "amean", 10, 10, "1",  2, kTwoToMinus50, "0.5000000000" },
		{ "1 mod 2^-34 is nought",          "modulo", 10, 10, "1",  2, kTwoToMinus34, "0" },
		{ "1.5 rem 3^-24 answers",          "remainder", 10, 10, "1.5", 3, "0.000000000000000000000001", "0.00000000000177035308074" },
		{ "3 times 1/3 is exactly 1",       "multiply", 10, 10, "3", 3, "0.1", "1" },
	};

	for ( auto& theCase : named )
	{
		long		l_result = kNoError;
		std::string	s_got = MixedPairAnswer( theCase.szOp, theCase.lWidth, theCase.lLeftBase, theCase.szLeft,
											 theCase.lRightBase, theCase.szRight, &l_result );
		std::string	s_detail = s_got + ((l_result != kNoError) ? (" code " + std::to_string( l_result )) : std::string( "" ));

		if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
			s_detail += std::string( ", wanted " ) + theCase.szExpected;

		CheckTrue( theCase.szLabel, (s_got == theCase.szExpected) && (l_result == kNoError), s_detail.c_str() );
	}

	struct SMixedCase { const char* szOp; long lWidth; long lLeftBase; const char* szLeft;
						long lRightBase; const char* szRight; const char* szExpected; };

	static const SMixedCase kGenerated[] =
	{
		{ "add", 10, 7, "30.2233", 12, "-0.00000000000487D975EE9E4", "30.22330000000" },
		{ "add", 10, 10, "0.000000062899748", 16, "-0.07EFE870DBB", "-0.03100439224061409715413860976696014404296875" },
		{ "add", 3, 16, "118A3F0.82FCC8", 5, "0.00000000011", "118A3F0.82FD" },
		{ "add", 3, 5, "0.0000412", 3, "201001.21", "4024.3422" },
		{ "add", 3, 4, "-2", 13, "-0.00000000011906", "-2.0000" },
		{ "add", 20, 11, "0.1775550661902651619587975", 14, "0.002C62B8260CC091", "0.1789A00AA5A59747682086" },
		{ "add", 10, 16, "0.D2DB929A", 5, "-0.000422342212", "0.D1028A606971" },
		{ "add", 36, 16, "1DBE8.016B4EDA3EAB42", 8, "0.0000425552105726353766441656341755720163", "1DBE8.0173FC82CA81B7FDA43AE70FB7A073" },
		{ "subtract", 3, 6, "0.0000000531", 12, "-8.6093", "12.3010535231" },
		{ "subtract", 3, 6, "3305.35", 10, "0.000000000403", "3305.3500" },
		{ "subtract", 3, 4, "-0.00000003012", 7, "0.0000000516642", "-0.000000030213" },
		{ "subtract", 36, 15, "0.0000001BDD9520869A69078E31EB2D07C846B32166", 6, "-533242.1052005201422525", "CD5E.2D0D6EE61C3235BE929070D466AA797B86CA3" },
		{ "subtract", 3, 16, "-52", 13, "-0.00000001AA", "-52.0000" },
		{ "subtract", 3, 3, "-0.0000002220011", 2, "-0.01", "0.020201" },
		{ "subtract", 20, 13, "569.1C995A4328701A8B0897A", 2, "-0.000000011", "569.2097A454522A8A5B58146" },
		{ "subtract", 20, 3, "11020.1210021111021", 14, "-0.000000000838851A0D42A9", "11020.121002111102100000000" },
		{ "multiply", 10, 12, "E114.73387484758424", 5, "-0.04", "-1936.E88336ED9D0" },
		{ "multiply", 20, 16, "CD778B7.EB", 10, "0.00000000018222530617774", "0.0A0CF1FF5B41C8738C2E298" },
		{ "multiply", 36, 7, "-2513622.222013663523511035236612060354034626", 6, "-0.00404133413522515513524551340402", "23650.0625131511452566061144254066565440450" },
		{ "multiply", 3, 8, "-4.643306", 13, "-0.0000000009C00B", "0.000000000027655" },
		{ "multiply", 3, 8, "-0.50240107", 5, "0.0000000324313", "-0.0000014102" },
		{ "multiply", 3, 16, "-55E06.42E8", 7, "-25.2035", "6791EC.14C3" },
		{ "multiply", 20, 7, "23250.63466135466401344002422", 14, "0.0000D1", "0.1005036600665001622241" },
		{ "multiply", 20, 15, "-0.00613CBD9C96C3CD25C4", 16, "-0.7216904A0BB83D", "0.002AA08A497E97B0218D537E" },
		{ "divide", 36, 14, "-56C.2D67C4B2C7D189", 16, "0.000000000AE3DC77989EB21855", "-1CB91B934B7DA.CA40991BDA5ABB5156267D9024C5595871D79" },
		{ "divide", 36, 15, "0.0002958276A6B70B997DB2", 11, "-0.0000104", "-8.1384B718031AE0C295881ACC9" },
		{ "divide", 36, 10, "-1722.92184", 9, "-0.00006643534110414055776311115823", "15137187.1672513564074760671859114442629747148" },
		{ "divide", 10, 4, "2.320302123022", 15, "0.00000075CC943", "3332220313130.22032323330" },
		{ "divide", 20, 9, "-11.2", 14, "0.000000192C4A476A265084813", "-1608227568.645055631871278266677" },
		{ "divide", 3, 13, "-0.7C9", 2, "0.0111", "-1.5302" },
		{ "divide", 10, 8, "7163454.10362147454", 10, "-436873.508072633", "-4.25372710774" },
		{ "divide", 36, 12, "0.00000007508E5E0184522411E4D", 10, "659990.72", "0.000000000000296E3881291863EDE7273D77E5544035E66156" },
		{ "divide_integer", 10, 2, "-1011.01111", 15, "0.0000000220BE00B23", "-1100110011111000000101010000000000" },
		{ "divide_integer", 3, 16, "-77C355.3C4", 3, "-0.002", "650CCFE" },
		{ "divide_integer", 3, 9, "6.240854", 7, "-0.0000351", "-42645" },
		{ "divide_integer", 3, 6, "-254.415", 7, "-0.00621046", "42510" },
		{ "divide_integer", 10, 4, "233.32101", 3, "0.000000002201020212", "1110310303" },
		{ "divide_integer", 3, 5, "-444320.1002031", 6, "0.0000003355", "-4414443000114" },
		{ "divide_integer", 3, 15, "523E34.CB8247", 7, "34.261", "3091B" },
		{ "divide_integer", 20, 15, "13.0C6B803D8", 9, "-0.0000001011246827822", "-7721B38" },
		{ "remainder", 36, 13, "-459C.7004A1C0219979BA5828C1C0713C", 9, "-0.00000000034102844881531085064358846", "-0.000000006444652CC8A069603B35363C91657304200AAC" },
		{ "remainder", 36, 2, "1.110101101110011011101001101", 7, "0.00000000005334", "0.0000000000000000000000000000010100111011000110110100011010100000101" },
		{ "remainder", 3, 5, "-41223.322331", 11, "0.0A40169", "-0.00043131" },
		{ "remainder", 3, 15, "-C13323", 6, "0.0404", "-0.130C8" },
		{ "remainder", 10, 6, "50.312334424", 12, "0.000000007E809421", "0.000000000000244101132522252343" },
		{ "remainder", 10, 3, "-21.010012", 7, "-0.00000000105443", "-0.0000000000000000100011212020" },
		{ "remainder", 20, 16, "-0.00000006A24553D72769672148E3BF", 4, "0.00000003", "-0.00000006A24553D72769672148E3BF" },
		{ "remainder", 10, 12, "-0.00000005538D7D25", 9, "-0.000000067482721", "-0.00000005538D7D25" },
		{ "remainder_near", 36, 11, "0.0000006A274A74359A67A102", 6, "0.0000000001504455231200510412541114525224052025", "-0.00000002778094361369A545A6786A301235063848030" },
		{ "remainder_near", 10, 14, "-0.035", 16, "0.000036", "0.00000A56369B54BA8761A7" },
		{ "remainder_near", 20, 13, "-5925.46752538", 6, "-0.000000000031235", "-0.000000004625A7BA40465C9C002C24" },
		{ "remainder_near", 10, 7, "0.011251312365135", 6, "-0.00002540324353", "-0.0000255604366234" },
		{ "remainder_near", 36, 9, "-714347", 2, "111.01000110111101100111111", "-2.3466776641300263076068617372556515681" },
		{ "remainder_near", 10, 12, "8.58D2", 5, "-0.000004123333321141", "-0.000221990732316" },
		{ "remainder_near", 10, 15, "-0.000000238E21860B1", 7, "-0.000064", "-0.000000238E21860B1" },
		{ "remainder_near", 36, 10, "0.0002702717816970324269401725", 6, "-253014.504114131253220540014004512020354425", "0.0002702717816970324269401725" },
		{ "modulo", 3, 12, "-0.0000516", 10, "0.0000000060562844", "-0.0000000210D0" },
		{ "modulo", 36, 6, "-21050.5555502554", 7, "0.23130410361645505123", "-0.012501331021453545352104100050434511113" },
		{ "modulo", 10, 16, "9902.989D85F", 13, "6.019714", "1.B78AE7298DC" },
		{ "modulo", 3, 8, "0.046150426", 12, "0.00000008094", "0.00000000010142" },
		{ "modulo", 36, 4, "-3231131.12132", 9, "-0.0000000005000368163752", "-0.00000000000000033212311230032132003203123121012232312" },
		{ "modulo", 3, 8, "-0.02347372", 11, "0.00000000797529", "-0.00000000025257" },
		{ "modulo", 36, 4, "-0.00000200223201333020220223013323", 6, "235.14541424200421353221505351013015043311", "-0.00000200223201333020220223013323" },
		{ "modulo", 20, 9, "-0.00805000783132750637", 3, "-10120.0021221010221021101101", "-0.00805000783132750637" },
		{ "compare", 3, 4, "0.000033", 7, "0.001153562", "1" },
		{ "compare", 3, 8, "1170213.36", 13, "B446A.612AC7946", "1" },
		{ "compare", 10, 3, "0.0000001120222211", 14, "0.001DD82B1A5D8873", "1" },
		{ "compare", 3, 11, "A06062.6", 8, "6131204.427213506", "-1" },
		{ "compare", 10, 4, "1.02012", 5, "1.0311343143243", "-1" },
		{ "compare", 20, 3, "2.0021011212011120102002", 16, "2.16B12F285A73DBFC9C7D62481B", "1" },
		{ "compare", 10, 9, "-38.28", 11, "-48.7", "1" },
		{ "compare", 10, 7, "-0.000000051", 5, "-0.00041141443044204", "1" },
		{ "max", 36, 13, "3B.1458217A6562335B9008B5555BCC728AC43", 12, "42.12D2ED37879E873259578EE29646822219736E0177", "-1" },
		{ "max", 3, 14, "-D.615", 8, "-15.336762442", "-1" },
		{ "max", 10, 7, "344.1406326020235", 12, "12E.288706543443E99", "-1" },
		{ "max", 3, 9, "-0.270303", 5, "-0.12330434", "1" },
		{ "max", 3, 15, "-9E27.D", 6, "-415214.511111111", "1" },
		{ "max", 10, 13, "-0.000005703507", 9, "-0.0000005435874137", "-1" },
		{ "max", 36, 9, "-0.76164783633657141", 7, "-0.565666134552004041426056620515643050464525", "-1" },
		{ "max", 3, 8, "-55", 11, "-41", "0" },
		{ "min", 36, 9, "-0.0006272662", 14, "-0.0028D635AB61D424D61AB51CA4CB22DC44856B729C", "-1" },
		{ "min", 10, 10, "-4040.7081065195", 16, "-FC8.B5467807563E2DC1", "1" },
		{ "min", 20, 5, "0.010243424034344", 4, "0.0023130313101213320133012", "-1" },
		{ "min", 10, 8, "-0.07052", 5, "-0.0234040012101103", "1" },
		{ "min", 36, 14, "0.00073940C524AB59CA5", 12, "0.0003E0524E72966D7D2D5E9D3149E75142471E1D7E", "1" },
		{ "min", 3, 7, "565052.04", 10, "100193.081632653", "1" },
		{ "min", 36, 12, "0.0822628E37E3E259165E593D11E5EE837024707", 16, "0.000000BAB8EB623031C65B9AC2", "1" },
		{ "min", 20, 15, "-9A.B", 5, "-13014.3034443424334", "1" },
		{ "power", 36, 7, "45", 2, "0.00000000000100001", "1.000205420362156631432164311106530635" },
		{ "power", 20, 13, "0.7A08827", 5, "-0.0000000034230041014423444", "1.000004C678BB48469316" },
		{ "power", 10, 10, "0.00719722", 3, "-0.00000021", "1.0052780813" },
		{ "power", 20, 8, "0.0455417", 4, "-0.000000000010111", "1.00000013225244723421" },
		{ "power", 20, 7, "0.25", 3, "-2.222121001202", "22.62420103533564156235" },
		{ "power", 20, 5, "21.123", 3, "0.00000002", "1.00002123420104310104" },
		{ "power", 20, 4, "0.022", 11, "-0.00009048582657451986814337", "1.00000012311001302333" },
		{ "power", 36, 11, "A.922", 15, "-0.000C6A5E95834394BEE15", "0.AAA24750421262A944442924803611492402" },
		{ "log_n", 10, 8, "2.172", 3, "0.0000021", "-7.1022652113" },
		{ "log_n", 3, 16, "0.8F", 2, "0.000000001", "A.B67" },
		{ "log_n", 36, 7, "55.24", 2, "0.0000000011111111011", "-1.333361161603522402664451031356350542" },
		{ "log_n", 10, 13, "C", 2, "0.00001", "-1.51925A7AB4" },
		{ "log_n", 36, 11, "9.71A", 2, "0.0000000000100100101110100001011", "-3.336363222A77823078285845687569360736" },
		{ "log_n", 10, 16, "0.0C18", 3, "0.212", "0.0D7281FD54" },
		{ "log_n", 20, 6, "14.3", 11, "46.6851092A6499", "1.40025303444402335530" },
		{ "log_n", 10, 2, "0.001", 13, "0.00C753491", "10.0111101111" },
		{ "amean", 10, 14, "85.D0616", 8, "0.0000000000000000000077445041", "42.D730A00000" },
		{ "amean", 10, 8, "1.017551016", 9, "0.0000000000000000000000000000074688775", "0.4076644070" },
		{ "amean", 10, 8, "0.02", 4, "0.00000000000000221230333332", "0.0100000001" },
		{ "amean", 20, 14, "A0D.96659D", 2, "0.00000000000000000000000000011111101", "506.BA32BD7560D03A3CB593" },
		{ "amean", 3, 16, "F9CE0.2", 2, "0.000000000000000000000000000000101", "7CE70.100" },
		{ "amean", 3, 16, "EF.4", 15, "0.0000000000987B", "77.A00" },
		{ "amean", 3, 6, "0.4", 12, "0.00000003139", "0.200" },
		{ "amean", 36, 7, "11.05264402301006632234201640341234", 11, "0.00A38975A68870373428292145A8496745", "4.026210531125230501060313320561560555" },
		{ "gmean", 20, 7, "0.015203", 15, "0.00000000000000000000000000000000000000008A98897C", "0.00000000000000000000" },
		{ "gmean", 10, 3, "2.020001", 8, "0.000000000000000000000000000013", "0.0000000000" },
		{ "gmean", 10, 15, "0.20B8", 3, "0.00000000000000000000000000001", "0.0000007961" },
		{ "gmean", 10, 9, "533.8476", 12, "0.00000000000000E1D387789083", "0.0000002616" },
		{ "gmean", 3, 3, "11.11", 11, "0.000000002", "0.000" },
		{ "gmean", 20, 3, "2012.101021", 12, "0.00000000000000000000000000000000000D0D76", "0.00000000000000000000" },
		{ "gmean", 10, 10, "6917.48266511", 14, "0.00151D8DC9640A09B", "1.8565243466" },
		{ "gmean", 20, 3, "0.000001012212000012", 8, "0.00000000000057631616076560217421463", "0.00000000000000122102" },
		{ "hmean", 10, 15, "4.2BB4A7", 4, "0.0000210113103", "0.00EE90C130" },
		{ "hmean", 36, 10, "2.9", 7, "0.0000000000000000000101", "0.000000000000000025576618151422140700" },
		{ "hmean", 20, 16, "A0FF.7F0084", 6, "0.00000000000000000000000000000100312421142432201323335", "0.0000000000000000000B" },
		{ "hmean", 36, 9, "0.00160581823323101733861", 4, "0.000203200030121122213132321202202120313322", "0.002574201082486273213172114017325807" },
		{ "hmean", 10, 2, "0.0111", 3, "0.00000000000000000000002", "0.0000000000" },
		{ "hmean", 20, 8, "0.03046473511", 3, "0.000000000000120112011121100112", "0.00000044233304665773" },
		{ "hmean", 10, 6, "0.01232", 2, "0.00111100011", "0.0223523044" },
		{ "hmean", 3, 14, "0.0000DA6", 5, "0.0000002", "0.000" },
	};

	const char*	kOps[] = { "add", "subtract", "multiply", "divide", "divide_integer", "remainder",
						   "remainder_near", "modulo", "compare", "max", "min",
						   "power", "log_n", "amean", "gmean", "hmean" };

	for ( const char* szOp : kOps )
	{
		long		l_checked = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( const SMixedCase& theCase : kGenerated )
		{
			if ( ::strcmp( theCase.szOp, szOp ) != 0 )
				continue;

			long		l_result = kNoError;
			std::string	s_got = MixedPairAnswer( theCase.szOp, theCase.lWidth, theCase.lLeftBase, theCase.szLeft,
												 theCase.lRightBase, theCase.szRight, &l_result );

			l_checked++;

			if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
				if ( l_wrong++ == 0 )
					s_first = std::string( theCase.szLeft ) + " (base " + std::to_string( theCase.lLeftBase ) + "), " +
							  theCase.szRight + " (base " + std::to_string( theCase.lRightBase ) + ") gave " + s_got;
		}

		char	szLabel[ 72 ];
		char	szDetail[ 400 ];

		::snprintf( szLabel, sizeof( szLabel ), "mixed %s", szOp );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld generated pairs", l_checked );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

		CheckTrue( szLabel, (l_wrong == 0) && (l_checked > 0), szDetail );
	}
}

static void TestLogarithmTies()
{
	Banner( "Logarithms that are fractions, decided exactly" );

	/*
		**log_n and log_b rounded an answer that is exactly a fraction on
		a half by its error until 2026-10-15**, in an odd base, where the
		half never ends and the retake ran to its cap: log_9 3 at one place
		answered 0.5 for 0.4. The decider finds p/q, proves a^p = x^q in
		whole numbers and rounds p/q exactly -- see lzi_log_decide. The
		table was generated: every answer an exact fraction t/s, most of
		them ties, scored exactly at the places written, and 32 of its 48
		log_n rows ties the old code rounded wrong; some take the value in
		another base.
	*/

	struct { const char* szLabel; const char* szOp; long lWidth; long lBase; const char* szLeft;
			 const char* szRight; const char* szExpected; } named[] =
	{
		{ "log_9 3 at one place is 0.4",    "log_b",  1, 9, "3",    "",  "0.4" },
		{ "log_9 27, a tie at 36",          "log_b", 36, 9, "30",   "",  "1.444444444444444444444444444444444444" },
		{ "log_4 2 in base five",           "log_n", 36, 5, "4",    "2", "0.222222222222222222222222222222222222" },
		{ "log_4 2 in base seven",          "log_n", 36, 7, "4",    "2", "0.333333333333333333333333333333333334" },
		{ "log_64 2 in base three",         "log_n",  1, 3, "2101", "2", "0.0" },

		// ...and near misses, which are no fractions and round by their
		//	digits: a decider that took the nearest fraction on trust
		//	rounds them as ties...
		{ "near miss below, base seven",    "log_n", 10, 7, "4", "1.6666666666666666", "0.3333333333" },
		{ "near miss below, base three",    "log_n", 10, 3, "11", "1.222222222222222222222222", "0.1111111111" },
		{ "near miss above, base five",     "log_n", 10, 5, "4", "2.000000000000000001", "0.2222222223" },
		{ "near miss above, log_b base 9",  "log_b",  1, 9, "3.00000001", "", "0.5" },
	};

	for ( auto& theCase : named )
	{
		long		l_result = kNoError;
		std::string	s_got;

		if ( ::strcmp( theCase.szOp, "log_n" ) == 0 )
			s_got = MixedPairAnswer( "log_n", theCase.lWidth, theCase.lBase, theCase.szLeft,
									 theCase.lBase, theCase.szRight, &l_result );
		else
		{
			long		l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lWidth;

			CLZNumber	x;
			CLZNumber	got;

			SetInBase( x, theCase.lBase, theCase.szLeft );
			l_result = lz_log_b( &x, &got );
			glWorkingPrecision = l_saved;

			char*	sz_got = got.GetStringRep();

			s_got = (sz_got != nullptr) ? sz_got : "(null)";
			::free( sz_got );
		}

		std::string	s_detail = s_got;

		if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
			s_detail += std::string( ", wanted " ) + theCase.szExpected;

		CheckTrue( theCase.szLabel, (s_got == theCase.szExpected) && (l_result == kNoError), s_detail.c_str() );
	}

	struct SLogCase { const char* szOp; long lWidth; long lBase; const char* szLeft;
					  long lRightBase; const char* szRight; const char* szExpected; };

	static const SLogCase kGenerated[] =
	{
		{ "log_n", 3, 3, "0.11", 3, "0.0011202", "10.112" },
		{ "log_n", 2, 7, "1552", 7, "146002132", "2.34" },
		{ "log_n", 1, 15, "4", 8, "0.02", "-2.8" },
		{ "log_n", 10, 9, "1000000000", 9, "3000", "0.3444444444" },
		{ "log_n", 3, 15, "7AB1A", 15, "16677E6A", "1.778" },
		{ "log_n", 2, 7, "16155440250542663336406446226550034330156416313020140642", 7, "11153624400054364136441524650306545126205", "0.50" },
		{ "log_n", 2, 3, "100", 3, "100000", "2.12" },
		{ "log_n", 36, 5, "100000000000000000000000000000000000000000000000000", 5, "0.0000000000000000000000000000000000000000001", "-0.412222222222222222222222222222222222" },
		{ "log_n", 2, 3, "2101", 12, "28", "0.22" },
		{ "log_n", 36, 9, "7815654114304431", 14, "61A7", "0.244444444444444444444444444444444444" },
		{ "log_n", 3, 9, "7815654114304431", 9, "13226734107536213504347", "1.344" },
		{ "log_n", 1, 11, "32881A9", 11, "1893", "0.6" },
		{ "log_n", 3, 3, "0.11", 3, "0.0011202", "10.112" },
		{ "log_n", 3, 15, "9", 15, "9AC", "3.778" },
		{ "log_n", 36, 7, "34", 5, "0.1", "-0.333333333333333333333333333333333334" },
		{ "log_n", 1, 11, "45", 5, "202323133", "3.6" },
		{ "log_n", 10, 3, "100", 3, "0.1", "-0.1111111112" },
		{ "log_n", 36, 7, "4", 7, "242", "3.333333333333333333333333333333333334" },
		{ "log_n", 1, 13, "9", 13, "21", "1.6" },
		{ "log_n", 10, 3, "221", 3, "10222011112", "10.1111111112" },
		{ "log_n", 1, 15, "78D151", 15, "560232201", "1.8" },
		{ "log_n", 3, 11, "4", 11, "2", "0.556" },
		{ "log_n", 3, 15, "339", 15, "5C73", "1.778" },
		{ "log_n", 3, 15, "16677E6A", 10, "0.04", "-0.278" },
		{ "log_n", 10, 15, "0.02BA", 15, "12769", "-2.7777777778" },
		{ "log_n", 1, 15, "0.0007AB1A", 15, "0.02BA", "0.8" },
		{ "log_n", 2, 13, "4", 13, "9B", "3.66" },
		{ "log_n", 36, 7, "304152012124", 5, "100000000000000000000000", "1.433333333333333333333333333333333334" },
		{ "log_n", 10, 15, "24CD4", 15, "B88D48CD", "1.C777777778" },
		{ "log_n", 10, 7, "10000", 7, "100", "0.3333333334" },
		{ "log_n", 3, 13, "1C", 16, "1312D", "3.666" },
		{ "log_n", 3, 15, "339", 15, "3774C", "1.C78" },
		{ "log_n", 2, 9, "10", 9, "0.003", "-2.44" },
		{ "log_n", 36, 3, "1000000", 3, "0.1", "-0.011111111111111111111111111111111112" },
		{ "log_n", 1, 6, "0.000000000333425113050213", 6, "5341344", "-0.4" },
		{ "log_n", 2, 9, "10", 9, "0.003", "-2.44" },
		{ "log_n", 10, 7, "12", 7, "6243", "3.3333333334" },
		{ "log_n", 1, 2, "1100111001000000010100100000100001", 2, "10011001111011111101000111", "1.0" },
		{ "log_n", 1, 2, "1011111010111100001", 2, "11001", "0.0" },
		{ "log_n", 2, 11, "9", 9, "0.003", "-2.56" },
		{ "log_n", 36, 12, "0.0194", 12, "0.0000426994", "2.6" },
		{ "log_n", 3, 14, "0.0000000000000022321C41B84A4067DC7D29914363142A4870AA0C37", 14, "144", "-0.2" },
		{ "log_n", 36, 6, "3213", 6, "0.00000013252", "-1.5" },
		{ "log_n", 3, 6, "1135540.351134350213", 6, "34321223523.235423233541453335343", "1.43" },
		{ "log_n", 2, 2, "10.01", 2, "111.10011", "10.1" },
		{ "log_n", 2, 14, "9543A0C37", 14, "55061A7", "0.A7" },
		{ "log_n", 2, 12, "E.483", 12, "32.53D16", "1.6" },
		{ "log_n", 36, 2, "0.0000000000000001", 2, "0.000001", "0.011" },
		{ "log_b", 2, 8, "1000", 0, "", "3.00" },
		{ "log_b", 36, 16, "0.8", 0, "", "-0.400000000000000000000000000000000000" },
		{ "log_b", 2, 16, "80", 0, "", "1.C0" },
		{ "log_b", 3, 8, "10000", 0, "", "4.000" },
		{ "log_b", 10, 4, "20", 0, "", "1.2000000000" },
		{ "log_b", 36, 8, "40", 0, "", "1.525252525252525252525252525252525253" },
		{ "log_b", 2, 8, "0.04", 0, "", "-1.25" },
		{ "log_b", 36, 8, "4000", 0, "", "3.525252525252525252525252525252525253" },
	};

	for ( const char* szOp : { "log_n", "log_b" } )
	{
		long		l_checked = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( const SLogCase& theCase : kGenerated )
		{
			if ( ::strcmp( theCase.szOp, szOp ) != 0 )
				continue;

			long		l_result = kNoError;
			std::string	s_got;

			if ( ::strcmp( szOp, "log_n" ) == 0 )
				s_got = MixedPairAnswer( "log_n", theCase.lWidth, theCase.lBase, theCase.szLeft,
										 theCase.lRightBase, theCase.szRight, &l_result );
			else
			{
				long		l_saved = glWorkingPrecision;

				glWorkingPrecision = theCase.lWidth;

				CLZNumber	x;
				CLZNumber	got;

				SetInBase( x, theCase.lBase, theCase.szLeft );
				l_result = lz_log_b( &x, &got );
				glWorkingPrecision = l_saved;

				char*	sz_got = got.GetStringRep();

				s_got = (sz_got != nullptr) ? sz_got : "(null)";
				::free( sz_got );
			}

			l_checked++;

			if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
				if ( l_wrong++ == 0 )
					s_first = std::string( theCase.szLeft ) + " (base " + std::to_string( theCase.lBase ) + ") gave " + s_got;
		}

		char	szLabel[ 72 ];
		char	szDetail[ 400 ];

		::snprintf( szLabel, sizeof( szLabel ), "%s, exact fractions", szOp );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld generated", l_checked );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

		CheckTrue( szLabel, (l_wrong == 0) && (l_checked > 0), szDetail );
	}
}

static void TestTinyArgumentTies()
{
	Banner( "Tiny arguments on a half, decided by their side of x" );

	/*
		**An odd function of a tiny argument lying on a half of the
		answer's places went by the kept digit until 2026-10-15**: only
		its cubic term moves it off the half, past the retake's cap. sin
		of 0.00000000055 at ten places answered ...6 for ...5. The table
		was generated -- even bases, the argument's last digit half the
		base one place past the answer's places, or past its W significant
		digits for a function that writes a small answer that way -- each
		answer scored against a Decimal reference, and six of each
		function's eight rows the old code rounded wrong where it had six
		-- see lzi_beside_x_decide. sinh had none while its cancellation
		floor swamped the tie; six of its own were added when the floor
		went.
	*/

	struct SFunction { const char* szName; long (*fn)( CLZNumber*, CLZNumber* ); };

	static const SFunction kFunctions[] =
	{
		{ "sine", lz_sine }, { "tangent", lz_tangent }, { "arctangent", lz_arctangent },
		{ "archypersine", lz_archypersine }, { "hypersine", lz_hypersine }, { "arcsine", lz_arcsine },
		{ "archypertangent", lz_archypertangent }, { "hypertangent", lz_hypertangent },
	};

	struct STinyCase { const char* szFunction; long lWidth; long lBase; const char* szX; const char* szExpected; };

	static const STinyCase kCases[] =
	{
		{ "sine", 10, 10, "0.00000000055", "0.0000000005" },

		// sinh's own ties, which its cancellation floor hid until the floor
		//	went -- see lz_hypersine; each needs the decider...
		{ "hypersine", 1, 4, "-0.00000000000022", "-0.0000000000003" },
		{ "hypersine", 3, 8, "0.0000000000005044", "0.000000000000505" },
		{ "hypersine", 3, 10, "0.000000000004645", "0.00000000000465" },
		{ "hypersine", 1, 12, "0.00000000046", "0.0000000005" },
		{ "hypersine", 3, 12, "-0.000000000007206", "-0.00000000000721" },
		{ "hypersine", 1, 14, "-0.000000000087", "-0.00000000009" },

		{ "sine", 36, 4, "-0.0000000000000000000000000000002301232", "-0.000000000000000000000000000000230123" },
		{ "sine", 10, 10, "0.00000000015", "0.0000000001" },
		{ "sine", 36, 10, "0.0000000000000000000000000051719027115", "0.000000000000000000000000005171902711" },
		{ "sine", 10, 12, "0.00000000036", "0.0000000003" },
		{ "sine", 36, 16, "0.0000000000000000000000000000000000EB8", "0.0000000000000000000000000000000000EB" },
		{ "sine", 3, 2, "-0.0001", "0.000" },
		{ "sine", 3, 2, "0.0011", "0.001" },
		{ "tangent", 1, 8, "0.00000000000044", "0.0000000000005" },
		{ "tangent", 1, 14, "0.0000000000C7", "0.0000000000D" },
		{ "tangent", 3, 16, "-0.000000000000008E08", "-0.000000000000008E1" },
		{ "tangent", 3, 2, "-0.0011", "-0.00110" },
		{ "tangent", 3, 2, "-0.0101", "-0.0101" },
		{ "arctangent", 36, 4, "0.0000000000000000000000000000000233212", "0.000000000000000000000000000000023321" },
		{ "arctangent", 36, 8, "0.0000000000000000000000000000000007434", "0.000000000000000000000000000000000743" },
		{ "arctangent", 10, 10, "-0.00000000015", "-0.0000000001" },
		{ "arctangent", 36, 10, "0.0000000000000000000000000000000005835", "0.000000000000000000000000000000000583" },
		{ "arctangent", 36, 12, "0.0000000000000000000000000000000001D36", "0.0000000000000000000000000000000001D3" },
		{ "arctangent", 36, 14, "-0.000000000000000000000000000003A109C37", "-0.000000000000000000000000000003A109C3" },
		{ "arctangent", 3, 2, "-0.0011", "-0.001" },
		{ "arctangent", 3, 2, "0.0011", "0.001" },
		{ "archypersine", 36, 2, "-0.0000000000000000000000000000000010111", "-0.000000000000000000000000000000001011" },
		{ "archypersine", 36, 6, "-0.0000000000000000000000000000000014113", "-0.000000000000000000000000000000001411" },
		{ "archypersine", 36, 8, "0.0000000000000000000000000000000000054", "0.000000000000000000000000000000000005" },
		{ "archypersine", 36, 8, "-0.0000000000000000000000000000005541374", "-0.000000000000000000000000000000554137" },
		{ "archypersine", 10, 10, "-0.00000000335", "-0.0000000033" },
		{ "archypersine", 36, 12, "-0.0000000000000000000000000000935250656", "-0.000000000000000000000000000093525065" },
		{ "archypersine", 3, 2, "0.0001", "0.000" },
		{ "archypersine", 3, 2, "-0.0101", "-0.010" },
		{ "hypersine", 3, 2, "0.0011", "0.00110" },
		{ "hypersine", 3, 2, "0.0011", "0.00110" },
		{ "arcsine", 36, 2, "0.0000000000000000000000000000011011001", "0.000000000000000000000000000001101101" },
		{ "arcsine", 36, 4, "0.0000000000000000000000000000000002302", "0.000000000000000000000000000000000231" },
		{ "arcsine", 36, 4, "-0.0000000000000000000000000000000000002", "-0.000000000000000000000000000000000001" },
		{ "arcsine", 36, 6, "-0.0000000000000000000000000000000324443", "-0.000000000000000000000000000000032445" },
		{ "arcsine", 36, 6, "0.0000000000000000000000000000000030343", "0.000000000000000000000000000000003035" },
		{ "arcsine", 36, 8, "-0.0000000000000000000000000000004665724", "-0.000000000000000000000000000000466573" },
		{ "arcsine", 3, 2, "-0.0001", "-0.001" },
		{ "arcsine", 3, 2, "-0.0001", "-0.001" },
		{ "archypertangent", 36, 2, "0.0000000000000000000000000000000000001", "0.000000000000000000000000000000000001" },
		{ "archypertangent", 36, 2, "-0.0000000000000000000000000000010011001", "-0.000000000000000000000000000001001101" },
		{ "archypertangent", 36, 4, "0.0000000000000000000000000000033213222", "0.000000000000000000000000000003321323" },
		{ "archypertangent", 36, 4, "-0.0000000000000000000000000000132113002", "-0.000000000000000000000000000013211301" },
		{ "archypertangent", 36, 6, "-0.0000000000000000000000000000000000003", "-0.000000000000000000000000000000000001" },
		{ "archypertangent", 10, 8, "-0.00000000004", "-0.0000000001" },
		{ "archypertangent", 3, 2, "0.0011", "0.010" },
		{ "archypertangent", 3, 2, "-0.0011", "-0.010" },
		{ "hypertangent", 1, 6, "-0.0000000000033", "-0.000000000003" },
		{ "hypertangent", 1, 14, "0.0000000B7", "0.0000000B" },
		{ "hypertangent", 3, 2, "0.0011", "0.00110" },
		{ "hypertangent", 3, 2, "0.0101", "0.0101" },
	};

	for ( const SFunction& theFunction : kFunctions )
	{
		long		l_checked = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( const STinyCase& theCase : kCases )
		{
			if ( ::strcmp( theCase.szFunction, theFunction.szName ) != 0 )
				continue;

			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lWidth;

			CLZNumber	x;
			CLZNumber	got;

			SetInBase( x, theCase.lBase, theCase.szX );

			long	l_result = theFunction.fn( &x, &got );

			glWorkingPrecision = l_saved;

			char*		sz_got = got.GetStringRep();
			std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

			::free( sz_got );
			l_checked++;

			if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
				if ( l_wrong++ == 0 )
					s_first = std::string( theCase.szX ) + " (base " + std::to_string( theCase.lBase ) + ") gave " + s_got +
							  ", wanted " + theCase.szExpected;
		}

		char	szLabel[ 72 ];
		char	szDetail[ 400 ];

		::snprintf( szLabel, sizeof( szLabel ), "%s, tiny ties", theFunction.szName );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld", l_checked );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

		CheckTrue( szLabel, (l_wrong == 0) && (l_checked > 0), szDetail );
	}
}

static void TestCarryKeepsDoubt()
{
	Banner( "A carry to a power of the base keeps its doubt" );

	/*
		**The trim after a carry was a second rounding until 2026-10-15**,
		recorded over the first with a tail of nought, so a retake never
		saw the doubt: a tangent 1e-12 of a unit under the half below 0.01
		answered 0.01. And a retaken answer decided up to a power of the
		base was written a place too wide. The table was generated: answers
		a hair either side of the half just below b^-k, within the
		retake's reach, each scored for its value against a Decimal
		reference and for its width by the quotient rule, compared here
		character for character; six of each function's eight the old code
		got wrong -- see lzi_trim_after_carry.
	*/

	struct SFunction { const char* szName; long (*fn)( CLZNumber*, CLZNumber* ); };

	static const SFunction kFunctions[] =
	{
		{ "tangent", lz_tangent }, { "cotangent", lz_cotangent }, { "hypersine", lz_hypersine },
		{ "hypertangent", lz_hypertangent }, { "hypercosecant", lz_hypercosecant }, { "hypersecant", lz_hypersecant },
	};

	struct SCarryCase { const char* szFunction; long lWidth; long lBase; const char* szX; const char* szExpected; };

	static const SCarryCase kCases[] =
	{
		{ "tangent", 10, 3, "0.002222200001021211200110101102010122110012212202121011", "0.01000000000" },
		{ "tangent", 36, 7, "0.0066664444603150316562360605435332261343250633366156321652422", "0.0100000000000000000000000000000000000" },
		{ "tangent", 10, 2, "0.0111011010100100110011101111011011001111010010100101100011010110011101101", "0.1000000000" },
		{ "tangent", 36, 10, "0.00099999966666686666652380963492054401112093504426752309192", "0.000999999999999999999999999999999999999" },
		{ "tangent", 10, 3, "0.00022222212221120100112212020020201200010201201", "0.001000000000" },
		{ "tangent", 36, 7, "0.06644602165551600300536102121506343462102532124310516034501411062346316665", "0.100000000000000000000000000000000000" },
		{ "tangent", 36, 16, "0.000FFFFFFAAAAAADDDDDDB94B94D5BD5BBE7789F51180128C110A557EC1C0DDE949", "0.00100000000000000000000000000000000000" },
		{ "tangent", 10, 3, "0.00022222212221120100112212012101111220120122001", "0.0002222222222" },
		{ "cotangent", 10, 16, "1.82250768AC5A878E178EC94C200BB", "0.0FFFFFFFFFF" },
		{ "cotangent", 36, 16, "1.920FB5444826D9BBABAEE982486ADDCC961D2371C3B797BF549E716713EFF", "0.000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "cotangent", 10, 2, "1.0001101101111010111001100100110101101111001110111000110111", "0.1000000000" },
		{ "cotangent", 36, 7, "1.3655324216605303005226344136065350031626066633526134650312231105563424335141", "0.000666666666666666666666666666666666666" },
		{ "cotangent", 36, 3, "1.0202011122220101212202222222022200020222011200010101112220201112000220100101002", "0.100000000000000000000000000000000000" },
		{ "cotangent", 10, 16, "1.920FB5444826E1BBABA6E993E804515AE6B", "0.000FFFFFFFFFF" },
		{ "cotangent", 10, 7, "1.300053366453056531040500510334531100651162402212", "0.06666666666" },
		{ "cotangent", 10, 3, "1.1121020100222012002221210111000010021102111022010120222", "0.0002222222222" },
		{ "hypersine", 10, 3, "0.00022222221110000200221220100120102110220210011", "0.0002222222222" },
		{ "hypersine", 10, 16, "0.000FFFFFFD554D56888C87D1AA1B492DCD75614B9", "0.000FFFFFFFFFF" },
		{ "hypersine", 36, 3, "0.02221120101210220202221201021200100100202121121002001120001000221201121", "0.0222222222222222222222222222222222222" },
		{ "hypersine", 36, 16, "0.0FFD5687D229313CBB7F9B35CE9C207E8C06D92FF6D959F2D9B221B408C47275", "0.0FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "hypersine", 10, 10, "0.0998340788942323772767482285237", "0.09999999999" },
		{ "hypersine", 36, 7, "0.00666655556233602151544342613005125642606546351202163262613116313", "0.00666666666666666666666666666666666666" },
		{ "hypersine", 36, 16, "0.00FFFFD55568887D1AD97431894A1CD80BFE1155A1F72D93412A8315B990A2F922", "0.0100000000000000000000000000000000000" },
		{ "hypersine", 10, 3, "0.002222211111200011201101201012120220220212011221211022", "0.01000000000" },
		{ "hypertangent", 10, 16, "0.1005588AD36DA4D4A6773E923654C", "0.0FFFFFFFFFF" },
		{ "hypertangent", 36, 16, "0.001000000555555888888AD1AD1C98C98E0D269D19AEAB41F3F60DC0D0411", "0.000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "hypertangent", 10, 2, "0.1000110010001001111111110110001110000011101101100101010101", "0.01111111111" },
		{ "hypertangent", 10, 10, "0.1003353477260250754842139513051", "0.09999999999" },
		{ "hypertangent", 36, 3, "0.100100122121102012211011022212121110220110102112221121211102121202100002001201", "0.0222222222222222222222222222222222222" },
		{ "hypertangent", 10, 16, "0.0010000005554D5888808AC00D8309FC0A7", "0.000FFFFFFFFFF" },
		{ "hypertangent", 10, 7, "0.001000000222156014235134613601211221", "0.0006666666666" },
		{ "hypertangent", 36, 10, "0.0100003333533347620158821075514042238820968087573393351217044976487249", "0.0100000000000000000000000000000000000" },
		{ "hypercosecant", 10, 3, "10.22220102000202220201102021202112211211011221211", "0.001000000000" },
		{ "hypercosecant", 10, 16, "3.777A5FE46628BDE0BD2947E061146788E", "0.0FFFFFFFFFF" },
		{ "hypercosecant", 10, 3, "2.220010022211020102100121012211101100002121110112212101", "0.002222222222" },
		{ "hypercosecant", 36, 16, "6.3D0317B6484B545F6596ABFA40DCD21FF9D4C5B17780640DC47FA8EC12FD", "0.00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "hypercosecant", 10, 16, "9.02CB37D5A80915B9009EADB71781D6AEA3D", "0.000FFFFFFFFFF" },
		{ "hypercosecant", 36, 16, "3.777A5FE465A8FDB0E5052E467B1CD3A6259245224505DA10C5E11E1", "0.0FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "hypercosecant", 36, 16, "6.3D0317B6484B545F6596ABFA40DCD21FF9D4C5B177556FEC0092010F", "0.0100000000000000000000000000000000000" },
		{ "hypercosecant", 10, 3, "10.2222010200020222020110202212121000101010210102", "0.0002222222222" },
		{ "hypersecant", 10, 3, "10.22220020112021110200110222020021012222100002222", "0.001000000000" },
		{ "hypersecant", 10, 16, "3.76FA5FC9BB71F9951FBB6ABA84D2B2BFB", "0.0FFFFFFFFFF" },
		{ "hypersecant", 10, 3, "2.21222121102220012110202011022102020100000020120200022", "0.002222222222" },
		{ "hypersecant", 36, 16, "6.3D0297B6484B39B4BAEBF4B5FC98863290FE330292BA13BF46928C252AB4", "0.00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "hypersecant", 10, 16, "9.02CB3755A80915B98084030C6DD0F58F726", "0.000FFFFFFFFFF" },
		{ "hypersecant", 36, 16, "3.76FA5FC9BAF1B964F7971126EF56F087DECC51FCBC7BCC25036B58D", "0.0FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
		{ "hypersecant", 36, 16, "6.3D0297B6484B39B4BAEBF4B5FC98863290FE3302928F1F728E6DA634", "0.0100000000000000000000000000000000000" },
		{ "hypersecant", 10, 3, "10.2222002011202111020011100000211202111201201122", "0.0002222222222" },
	};

	for ( const SFunction& theFunction : kFunctions )
	{
		long		l_checked = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( const SCarryCase& theCase : kCases )
		{
			if ( ::strcmp( theCase.szFunction, theFunction.szName ) != 0 )
				continue;

			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lWidth;

			CLZNumber	x;
			CLZNumber	got;

			SetInBase( x, theCase.lBase, theCase.szX );

			long	l_result = theFunction.fn( &x, &got );

			glWorkingPrecision = l_saved;

			char*		sz_got = got.GetStringRep();
			std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

			::free( sz_got );
			l_checked++;

			if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
				if ( l_wrong++ == 0 )
					s_first = "base " + std::to_string( theCase.lBase ) + " at " + std::to_string( theCase.lWidth ) +
							  " gave " + s_got + ", wanted " + theCase.szExpected;
		}

		char	szLabel[ 72 ];
		char	szDetail[ 400 ];

		::snprintf( szLabel, sizeof( szLabel ), "%s, carried", theFunction.szName );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld, value and width", l_checked );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

		CheckTrue( szLabel, (l_wrong == 0) && (l_checked > 0), szDetail );
	}
}

static void TestAmeanRoundedOnce()
{
	Banner( "The arithmetic mean, rounded once" );

	/*
		**amean rounded twice until 2026-10-15** -- its quotient at the
		width and a guard, then at the width -- so a tie at the width broken
		past the guard went by the kept digit: amean( 0.1234567892,
		0.1234567893000000002 ) at ten places answered ...892 for ...893.
		The table was generated, each mean scored exactly, fourteen of its
		eighteen the old code rounded wrong -- see lz_amean.
	*/

	struct SMeanCase { long lBase; long lWidth; const char* szValues; const char* szExpected; };

	static const SMeanCase kCases[] =
	{
		{ 10, 10, "0.1234567892,0.1234567893000000002", "0.1234567893" },
		{ 8, 1, "5254.15,6023.73000000000002", "5540.1" },
		{ 10, 36, "0.08134535510636781071065821082733337781,-54.4793604895855308488824168961087603468100000000000000000000002", "-27.199007567239581519085879342640713485" },
		{ 6, 1, "-0.4,23.24,0.2345,2.122,-1.343,133.40205555555555555555555", "20.1" },
		{ 3, 10, "-2.12110211,-2002.1122000210211,0.01210000000102,2112.12101002212001", "2.1020021012" },
		{ 16, 3, "-F.7A2A,10A.757A0000000000002", "7D.7DB" },
		{ 9, 10, "87.823436741201,685.650676085027,4.47611527078,-1862.4600270163,-0.38507534284275,1523.72316301442375", "67.6043588704" },
		{ 10, 10, "-0.0147649609315,-7.37261779118,0.187347972287,-1.435395196174,-7.0892285674,16.0214288126984999999999994", "0.0494617115" },
		{ 15, 10, "0.00D8176D8AC492,-0.07DA523266B4A7,-0.457872DC609AA4,2B8D7.230C5DD7288AB9", "A5E9.315BA048EA" },
		{ 2, 10, "-101.111111101011,-1.0111110111101,-10.1000000111,-0.00100011101011,11000.1000100011100011111111111011", "10.1110000101" },
		{ 16, 1, "8039.7143,FBC.BEBCFFFFFFFFFFFFFE", "47FB.1" },
		{ 4, 1, "-0.1223,312.2031,-1002.200200000000000000000003", "-10.1" },
		{ 15, 36, "-39CA.21499508799D093D346E0D714E6659D86C6,39D4.1A051BDDBAA5E989D091B62E125723B88BA2", "4.74553AE29803E79D4D1153CDD1705BE00E98" },
		{ 6, 3, "132.053,454.533,0.02431321,-312.112,-0.02101,-320.5143032055555555555555555555555", "-0.155" },
		{ 2, 10, "-0.00110000111,-0.01000110011000000000000000000000001", "-0.0011101111" },
		{ 5, 3, "-0.003413,1012.201,0.003304,0.003313,-31.4,401.21024044444444444444444444444434", "121.042" },
		{ 15, 1, "B.5A,0.00E9519,E60E.3514,-E61A.6C4D518EEEEEEEEEEB", "0.1" },
		{ 12, 10, "36.E45340637405,8271.DE4E818926483,-8279.5E31D084EEE133", "E.9564486273" },
		{ 11, 3, "0.084A1,-0.0598286,-2613.35526,-0.004218,AAA3.3A2762,-8490.0226035", "0.009" },
	};

	long		l_checked = 0;
	long		l_wrong = 0;
	std::string	s_first;

	for ( const SMeanCase& theCase : kCases )
	{
		std::vector<std::string>	v_text;
		std::string					s_values = theCase.szValues;

		for ( size_t n_start = 0; ; )
		{
			size_t	n_comma = s_values.find( ',', n_start );

			v_text.push_back( s_values.substr( n_start, (n_comma == std::string::npos) ? std::string::npos : (n_comma - n_start) ) );

			if ( n_comma == std::string::npos )
				break;

			n_start = n_comma + 1;
		}

		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = theCase.lWidth;

		std::vector<CLZNumber>		v_numbers( v_text.size() );
		std::vector<CLZNumber*>		v_list;
		CLZNumber					got;

		for ( size_t i = 0; i < v_text.size(); i++ )
		{
			SetInBase( v_numbers[ i ], theCase.lBase, v_text[ i ].c_str() );
			v_list.push_back( &v_numbers[ i ] );
		}

		long	l_result = lz_amean( v_list.data(), static_cast<long>( v_list.size() ), &got );

		glWorkingPrecision = l_saved;

		char*		sz_got = got.GetStringRep();
		std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

		::free( sz_got );
		l_checked++;

		if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
			if ( l_wrong++ == 0 )
				s_first = std::string( theCase.szValues ).substr( 0, 60 ) + " (base " + std::to_string( theCase.lBase ) +
						  ") gave " + s_got + ", wanted " + theCase.szExpected;
	}

	char	szDetail[ 400 ];

	if ( l_wrong == 0 )
		::snprintf( szDetail, sizeof( szDetail ), "all %ld, character for character", l_checked );
	else
		::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

	CheckTrue( "means of short lists", l_wrong == 0, szDetail );

	// And a long list: 100,001 values on a half at 36 places, one of them
	//	moved 1E-60 up, which takes the mean past the half by 1E-60 over
	//	100,001 -- far below the guard. The kept digit is even, so the two
	//	roundings went down...
	{
		const long	kCount = 100001;
		long		l_saved = glWorkingPrecision;

		glWorkingPrecision = 36;

		std::vector<CLZNumber>	v_numbers( kCount );
		std::vector<CLZNumber*>	v_list;
		CLZNumber				got;

		for ( long i = 0; i < kCount; i++ )
		{
			v_numbers[ i ].Assign( (i == 0) ? "0.1234567890123456789012345678901234565000000000000000000001"
										   : "0.1234567890123456789012345678901234565" );
			v_list.push_back( &v_numbers[ i ] );
		}

		lz_amean( v_list.data(), kCount, &got );
		glWorkingPrecision = l_saved;

		CheckRendering( "a mean of 100,001, past the half", &got, "0.123456789012345678901234567890123457" );
	}
}

static void TestNearOneDecided()
{
	Banner( "arccos and arccosh a hair from one, decided" );

	/*
		**arccos(1 - d) and arccosh(1 + d) with sqrt(2d) on a half of the
		answer's places went by their error until 2026-10-15** -- only the
		d/12 term moves them off it, past the retake's cap: arccos of
		0.99999999999999999921875 at ten places answered ...12 for ...13.
		The table was generated, sqrt(2d) on a half exactly in even bases
		or far down in odd ones and d moved a hair either way, each answer
		scored against a Decimal reference; fourteen the old code rounded
		wrong -- see lzi_near_one_decide.
	*/

	struct SFunction { const char* szName; long (*fn)( CLZNumber*, CLZNumber* ); };

	static const SFunction kFunctions[] = { { "arccosine", lz_arccosine }, { "archypercosine", lz_archypercosine } };

	struct SNearCase { const char* szFunction; long lWidth; long lBase; const char* szX; const char* szExpected; };

	static const SNearCase kCases[] =
	{
		{ "arccosine", 10, 10, "0.99999999999999999921875", "0.0000000013" },
		{ "arccosine", 36, 4, "0.33333333333333333333333333333333333333333333333333333333333333333333331132", "0.000000000000000000000000000000000011" },
		{ "arccosine", 36, 6, "0.555555555555555555555555555555555555555555555555555555555555555533543511513", "0.000000000000000000000000000000005145" },
		{ "arccosine", 36, 9, "0.88888888888888888888888888888888888888888888888888888888888866345321632277777777777777777777777777777777777777777777778114324561455", "0.000000000000000000000000000000637033" },
		{ "arccosine", 10, 10, "0.9999999999999999994487499999999999999944875", "0.0000000011" },
		{ "arccosine", 36, 10, "0.999999999999999999999999999999999999999999999999999999999999999999674778874999999999999999999999999999999999999999674778875", "0.000000000000000000000000000000000807" },
		{ "arccosine", 36, 12, "0.EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE63D5667188E87075D58ED6", "0.00000000000000000000000000E80D900061" },
		{ "arccosine", 10, 16, "0.FFFFFFFFFFFFFFFFB61CE00000000000000049E32", "0.00000000C3" },
		{ "archypercosine", 36, 2, "1.000000000000000000000000000000000000000000000000000000000000111001000111001", "0.000000000000000000000000000001010101" },
		{ "archypercosine", 36, 4, "1.00000000000000000000000000000000000000000000000000000013130122121000223302000000000000000000000000000000000000000000000000000000000013130122121000223302", "0.000000000000000000000000000331233011" },
		{ "archypercosine", 36, 6, "1.000000000000000000000000000000000000000000000000000000000000000004520130043", "0.000000000000000000000000000000003043" },
		{ "archypercosine", 36, 9, "1.00000000000000000000000000000000000000000000000000000000000000000013863111111111111111111111111111111111111111111111111111111111111111111111111086137", "0.000000000000000000000000000000000507" },
		{ "archypercosine", 36, 10, "1.0000000000000000000000000000000000000000000000000000000000003867834321011250000000000000000000000000000000000000000000000386783432101125", "0.000000000000000000000000000000879527" },
		{ "archypercosine", 36, 14, "1.00000000000000000000000000000000000000000000000000000000000000000037239A1A7", "0.0000000000000000000000000000000009CB" },
		{ "archypercosine", 36, 16, "1.00000000000000000000000000000000000000000000000000003DCCD71E17807CA3AE9C1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC23328E1E87F835C5163E", "0.00000000000000000000000000B1E18ED307" },
		{ "arccosine", 3, 2, "0.1110010111111", "0.100" },
		{ "arccosine", 3, 2, "0.1111001101100111", "0.011" },
		{ "arccosine", 3, 2, "0.1111001101100111", "0.011" },
		{ "arccosine", 3, 2, "0.111001111", "0.100" },
		{ "arccosine", 10, 2, "0.11111111111111111001111", "0.0000000100" },
		{ "arccosine", 10, 2, "0.11111110111010001110111", "0.0001011111" },
	};

	for ( const SFunction& theFunction : kFunctions )
	{
		long		l_checked = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( const SNearCase& theCase : kCases )
		{
			if ( ::strcmp( theCase.szFunction, theFunction.szName ) != 0 )
				continue;

			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lWidth;

			CLZNumber	x;
			CLZNumber	got;

			SetInBase( x, theCase.lBase, theCase.szX );

			long	l_result = theFunction.fn( &x, &got );

			glWorkingPrecision = l_saved;

			char*		sz_got = got.GetStringRep();
			std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

			::free( sz_got );
			l_checked++;

			if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
				if ( l_wrong++ == 0 )
					s_first = "base " + std::to_string( theCase.lBase ) + " at " + std::to_string( theCase.lWidth ) +
							  " gave " + s_got + ", wanted " + theCase.szExpected;
		}

		char	szLabel[ 72 ];
		char	szDetail[ 400 ];

		::snprintf( szLabel, sizeof( szLabel ), "%s near one", theFunction.szName );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld", l_checked );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

		CheckTrue( szLabel, (l_wrong == 0) && (l_checked > 0), szDetail );
	}
}

static void TestAbsOfNought()
{
	Banner( "abs of nought, and of a constant" );

	/*
		**lz_abs made every answer positive until 2026-10-15**, a nought
		included, and a nought's sign is kZero: Compare reads the sign
		first, so abs(0) ranked above nought by lz_compare, Compare, lz_max
		and lz_min, however the nought was made. Noughts made five ways in
		every base, abs taken into a target and in place, each compared with
		nought five ways -- see lz_abs.
	*/

	long	l_cases = 0;
	long	l_wrong = 0;

	for ( long l_base = 2; l_base <= 16; l_base++ )
	{
		CLZNumber	noughts[ 5 ];
		CLZNumber	x;
		CLZNumber	tiny;
		CLZNumber	minus_one;

		noughts[ 0 ].SetBase( l_base, false );
		noughts[ 0 ].Assign( 0L );
		SetInBase( noughts[ 1 ], l_base, "0.000" );
		SetInBase( noughts[ 2 ], l_base, "-0" );
		SetInBase( x, l_base, "1.1" );
		lz_subtract( &x, &x, &noughts[ 3 ] );
		SetInBase( minus_one, l_base, "-1" );
		lz_multiply( &minus_one, &noughts[ 0 ], &noughts[ 4 ] );

		for ( CLZNumber& nought : noughts )
		{
			CLZNumber	zero;
			CLZNumber	into;
			CLZNumber	in_place;

			zero.SetBase( l_base, false );
			zero.Assign( 0L );
			lz_abs( &nought, &into );
			in_place.SetBase( l_base, false );
			in_place.Assign( &nought );
			lz_abs( &in_place );

			for ( CLZNumber* p : { &into, &in_place } )
			{
				l_cases++;

				if ( (lz_compare( p, &zero ) != kEqualTo) || (p->Compare( &zero ) != kEqualTo) ||
					 (zero.Compare( p ) != kEqualTo) || (lz_max( p, &zero ) != kEqualTo) || (lz_min( p, &zero ) != kEqualTo) )
					l_wrong++;
			}
		}
	}

	char	szDetail[ 120 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld equal to nought five ways", l_cases - l_wrong, l_cases );
	CheckTrue( "abs of a nought is nought", l_wrong == 0, szDetail );

	// ...and a constant is not written in place, which it was...
	CLZNumber	constant;

	constant.Assign( "-2.5" );
	constant.SetConstant( true );
	lz_abs( &constant );
	CheckRendering( "abs in place leaves a constant", &constant, "-2.5" );

	// ...while an ordinary negative still turns, both ways...
	CLZNumber	negative;
	CLZNumber	positive;

	negative.Assign( "-2.5" );
	lz_abs( &negative, &positive );
	CheckRendering( "abs of -2.5", &positive, "2.5" );
	lz_abs( &negative );
	CheckRendering( "abs of -2.5 in place", &negative, "2.5" );
}

static void TestSmallSinhKeepsItsDigits()
{
	Banner( "sinh and tanh of a small argument keep their digits" );

	/*
		**sinh's e^x - e^-x cancelled away an argument's noughts until
		2026-10-15**, and the rule for a quotient writes a small answer to
		W significant digits, so the last ones -- with the argument below
		the width, all of them -- were the subtraction's floor:
		sinh(-0.00000000002332) in base four at three digits answered
		-0.000000000100 for -0.0000000000300, a tie its side of x decides.
		tanh divides it by cosh and had the same. The table
		was generated, every answer scored against a Decimal reference at
		the places written, ten of each function's thirteen the old code got
		wrong; csch and coth, whose widening halved with this, are controls
		-- see lz_hypersine.
	*/

	struct SFunction { const char* szName; long (*fn)( CLZNumber*, CLZNumber* ); };

	static const SFunction kFunctions[] =
	{
		{ "hypersine", lz_hypersine }, { "hypertangent", lz_hypertangent },
		{ "hypercosecant", lz_hypercosecant }, { "hypercotangent", lz_hypercotangent },
	};

	struct SSmallCase { const char* szFunction; long lWidth; long lBase; const char* szX; const char* szExpected; };

	static const SSmallCase kCases[] =
	{
		{ "hypersine", 3, 4, "-0.00000000002332", "-0.0000000000300" },
		{ "hypersine", 3, 3, "0.0000000000011", "0.00000000000110" },
		{ "hypersine", 36, 3, "0.000000000000000000000000000000000000000212202222210022120000112002000212102212", "0.000000000000000000000000000000000000000212202222210022120000112002000212110" },
		{ "hypersine", 36, 3, "0.0000000000000000002011212002200200221201201220210222101021", "0.000000000000000000201121200220020022120120122021022210" },
		{ "hypersine", 3, 4, "-0.0000000000331311", "-0.0000000000332" },
		{ "hypersine", 3, 4, "0.000000000203301", "0.000000000210" },
		{ "hypersine", 3, 4, "-0.0000000000221", "-0.0000000000221" },
		{ "hypersine", 10, 4, "0.0000000000000000331231", "0.00000000000000003312310000" },
		{ "hypersine", 10, 4, "0.00000000000000000302202210101", "0.000000000000000003022022101" },
		{ "hypersine", 36, 4, "-0.000000000000000000000000010012320332203303133113323203021330211213", "-0.0000000000000000000000000100123203322033031331133232030213302" },
		{ "hypersine", 36, 4, "-0.00000000000000000000000000000000000000010131230110220000212313320012310003", "-0.000000000000000000000000000000000000000101312301102200002123133200123100030" },
		{ "hypersine", 1, 2, "-0.0001111011", "-0.001" },
		{ "hypersine", 1, 2, "0.00001011", "0.00001" },
		{ "hypersine", 1, 2, "0.00000011111", "0.000001" },
		{ "hypertangent", 10, 4, "-0.000000000000000000000001322331132330231", "-0.000000000000000000000001322331133" },
		{ "hypertangent", 10, 4, "0.000000000000000000000001311231101", "0.000000000000000000000001311231101" },
		{ "hypertangent", 10, 4, "0.00000000000000000000233312", "0.000000000000000000002333120000" },
		{ "hypertangent", 36, 4, "-0.0000000000000000000000000000000000001111301203211230331002232003103320202", "-0.000000000000000000000000000000000000111130120321123033100223200310332020" },
		{ "hypertangent", 10, 5, "-0.00000000000000242021403", "-0.000000000000002420214030" },
		{ "hypertangent", 36, 5, "-0.00000000000000000000000000000000000000000044120032040223443203340444424302", "-0.000000000000000000000000000000000000000000441200320402234432033404444243020000" },
		{ "hypertangent", 36, 5, "0.00000000000000021231224441", "0.000000000000000212312244404444444444444444444440424" },
		{ "hypertangent", 36, 5, "-0.00000000000000000000000000000000000000000000000131014222304231234341", "-0.00000000000000000000000000000000000000000000000131014222304231234341000000000000000" },
		{ "hypertangent", 10, 6, "-0.000000000000000025404441053", "-0.00000000000000002540444105" },
		{ "hypertangent", 36, 6, "0.000000000000000000000000000010441420252203302323031", "0.0000000000000000000000000000104414202522033023230310000000000000" },
		{ "hypertangent", 1, 2, "0.001101", "0.01" },
		{ "hypertangent", 1, 2, "0.00000111", "0.00001" },
		{ "hypertangent", 1, 2, "0.00000011", "0.0000001" },
		{ "hypercosecant", 1, 2, "-0.0001111011", "-1000.1" },
		{ "hypercosecant", 1, 2, "0.00001011", "10111.1" },
		{ "hypercosecant", 1, 2, "0.00000011111", "1000010.0" },
		{ "hypercosecant", 1, 2, "0.11", "1.0" },
		{ "hypercosecant", 1, 2, "-0.10001", "-10.0" },
		{ "hypercosecant", 1, 2, "0.00010101", "1100.0" },
		{ "hypercotangent", 1, 2, "0.001101", "101.0" },
		{ "hypercotangent", 1, 2, "0.00000111", "100100.1" },
		{ "hypercotangent", 1, 2, "0.00000011", "1010101.1" },
		{ "hypercotangent", 1, 2, "0.01", "100.0" },
		{ "hypercotangent", 1, 2, "-0.0000001", "-10000000.0" },
		{ "hypercotangent", 1, 2, "-0.1", "-10.0" },
	};

	for ( const SFunction& theFunction : kFunctions )
	{
		long		l_checked = 0;
		long		l_wrong = 0;
		std::string	s_first;

		for ( const SSmallCase& theCase : kCases )
		{
			if ( ::strcmp( theCase.szFunction, theFunction.szName ) != 0 )
				continue;

			long	l_saved = glWorkingPrecision;

			glWorkingPrecision = theCase.lWidth;

			CLZNumber	x;
			CLZNumber	got;

			SetInBase( x, theCase.lBase, theCase.szX );

			long	l_result = theFunction.fn( &x, &got );

			glWorkingPrecision = l_saved;

			char*		sz_got = got.GetStringRep();
			std::string	s_got = (sz_got != nullptr) ? sz_got : "(null)";

			::free( sz_got );
			l_checked++;

			if ( (s_got != theCase.szExpected) || (l_result != kNoError) )
				if ( l_wrong++ == 0 )
					s_first = std::string( theCase.szX ) + " (base " + std::to_string( theCase.lBase ) + ") gave " + s_got +
							  ", wanted " + theCase.szExpected;
		}

		char	szLabel[ 72 ];
		char	szDetail[ 400 ];

		::snprintf( szLabel, sizeof( szLabel ), "%s, small arguments", theFunction.szName );

		if ( l_wrong == 0 )
			::snprintf( szDetail, sizeof( szDetail ), "all %ld", l_checked );
		else
			::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld wrong, first %s", l_wrong, l_checked, s_first.c_str() );

		CheckTrue( szLabel, (l_wrong == 0) && (l_checked > 0), szDetail );
	}
}

static void TestGammaInPlace()
{
	Banner( "gamma of a whole number in place" );

	/*
		**gamma of a whole number wrote its factorial into a copy until
		2026-10-15**, so in place the operand was left as it was: gamma(5)
		in place stayed 5 -- see lz_gamma.
	*/

	CLZNumber	x;
	CLZNumber	into;
	CLZNumber	seven;
	CLZNumber	constant;

	x.Assign( "5" );
	lz_gamma( &x );
	CheckRendering( "gamma(5) in place", &x, "24" );

	x.Assign( "5" );
	lz_gamma( &x, &into );
	CheckRendering( "gamma(5) into a target", &into, "24" );
	CheckRendering( "gamma(5) leaves its operand", &x, "5" );

	SetInBase( seven, 7, "6" );
	lz_gamma( &seven );
	CheckRendering( "gamma(6) in place, base seven", &seven, "231" );

	constant.Assign( "5" );
	constant.SetConstant( true );
	lz_gamma( &constant );
	CheckRendering( "gamma(5) leaves a constant", &constant, "5" );
}

static void TestQuantizeIntoItsQuantum()
{
	Banner( "quantize into its own quantum" );

	/*
		**quantize read the quantum's places after writing the target**,
		which may be the quantum: quantize(0.7, 3) into the 3 answered 0.7
		until 2026-10-15 -- see lz_quantize.
	*/

	CLZNumber	v;
	CLZNumber	q;
	CLZNumber	into;

	v.Assign( "0.7" );
	q.Assign( "3" );
	lz_quantize( &v, &q, &q );
	CheckRendering( "quantize(0.7, 3) into the 3", &q, "1" );

	q.Assign( "3.00" );
	lz_quantize( &v, &q, &q );
	CheckRendering( "quantize(0.7, 3.00) into it", &q, "0.70" );

	q.Assign( "3" );
	lz_quantize( &v, &q, &into );
	CheckRendering( "quantize(0.7, 3) into a target", &into, "1" );
}

static void TestAssignAWideDouble()
{
	Banner( "Assign( double ), whole digits and infinities" );

	/*
		**Assign( double ) wrote "%f" into 64 bytes until 2026-10-15**:
		1e308 kept its first 63 characters and assigned about 1e62 with no
		error, and an infinity read "inf" as malformed text and came back
		NaN -- see CLZNumber::Assign( double ). The double nearest 1e308 is
		a whole number of 309 digits; its first and last digits are checked,
		and its count.
	*/

	CLZNumber	wide;
	char*		sz_wide = nullptr;

	wide.Assign( 1e308 );
	sz_wide = wide.GetStringRep();

	std::string	s_wide = (sz_wide != nullptr) ? sz_wide : "";

	::free( sz_wide );

	CheckTrue( "1e308 keeps all 309 whole digits",
			   (s_wide.size() >= 309) && (s_wide.compare( 0, 9, "100000000" ) == 0) && (s_wide.substr( 0, 309 ).find( '.' ) == std::string::npos) &&
			   (s_wide.compare( 300, 9, "223118336" ) == 0),
			   s_wide.substr( 0, 20 ).c_str() );

	CLZNumber	infinite;
	long		l_result = infinite.Assign( HUGE_VAL );

	CheckTrue( "an infinity is one", infinite.FlagIsOn( kfInfinity ) && !(infinite.FlagIsOn( kfQNaN )) && (l_result == kNoError), "" );

	CLZNumber	not_a_number;

	l_result = not_a_number.Assign( std::nan( "" ) );
	CheckTrue( "a NaN is one, and refused", not_a_number.FlagIsOn( kfQNaN ) && (l_result == kFailed), "" );

	CLZNumber	ordinary;

	ordinary.Assign( -2.5 );
	CheckRendering( "an ordinary double", &ordinary, "-2.500000" );
}

static void TestPrefixAtTheFront()
{
	Banner( "A lower-case prefix at the front, in bases twelve to sixteen" );

	/*
		**A prefix wins at the front, whatever base the number holds** --
		Rich's call, 2026-10-15, item 14 of that day's review: text means
		the same in every base, as Elzed reads it. So "0b1" given to a
		base-sixteen number is one in base two; the digits are written with
		an upper-case letter, which is how the engine writes them, or
		without the leading nought.
	*/

	struct { long lBase; const char* szText; long lBecomes; const char* szSpelled; } cases[] =
	{
		{ 16, "0b1",  2,  "1"  }, { 16, "0B1", 16, "B1" }, { 16, "b1", 16, "B1" },
		{ 16, "0d1",  10, "1"  }, { 16, "0D1", 16, "D1" },
		{ 12, "0d1",  10, "1"  }, { 12, "0D1", 12, "D1" }, { 14, "0b101", 2, "101" },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber	x;
		char		szLabel[ 72 ];

		SetInBase( x, theCase.lBase, theCase.szText );
		::snprintf( szLabel, sizeof( szLabel ), "\"%s\" in base %ld is base %ld", theCase.szText, theCase.lBase, theCase.lBecomes );
		CheckTrue( szLabel, x.m_lBase == theCase.lBecomes, "" );
		CheckRendering( szLabel, &x, theCase.szSpelled );
	}
}

static void TestFlagsCarriedThroughDivision()
{
	Banner( "Infinities and NaNs through exp_b, division and the remainders" );

	/*
		**exp_b, divide_integer and the division to a precision propagated
		nothing until 2026-10-15**, and the remainders read a NaN or an
		infinite divisor as nought; a division by nought in place reported
		nought over nought -- see lzi_propagate_flags and lzi_remainder_by.
		Each kind of operand is made the way a host meets it: a NaN from
		malformed text, an infinity from one over nought.
	*/

	auto	nan = []( CLZNumber& x ) { x.Assign( "junk" ); };
	auto	infinity = []( CLZNumber& x )
	{
		CLZNumber	one;
		CLZNumber	zero;

		one.Assign( 1L );
		zero.Assign( 0L );
		lz_divide( &one, &zero, &x );
		x.ClearFlag( kfDiv0 );
	};
	auto	is = []( CLZNumber& r, bool bInfinite, bool bNaN, bool bDiv0 ) -> bool
	{
		return (r.FlagIsOn( kfInfinity ) == bInfinite) && (r.FlagIsOn( kfQNaN ) == bNaN) && (r.FlagIsOn( kfDiv0 ) == bDiv0);
	};

	{
		CLZNumber	x;
		CLZNumber	r;

		nan( x );
		lz_exp_b( &x, &r );
		CheckTrue( "exp_b of a NaN is a NaN", is( r, false, true, false ), "" );
		infinity( x );
		lz_exp_b( &x, &r );
		CheckTrue( "exp_b of an infinity is a NaN", is( r, false, true, false ), "the sign is not carried" );
	}

	{
		CLZNumber	x;
		CLZNumber	y;
		CLZNumber	r;

		nan( x );
		y.Assign( 7L );
		lz_divide_integer( &x, &y, &r );
		CheckTrue( "divide_integer( NaN, 7 ) is a NaN", is( r, false, true, false ), "" );

		infinity( x );
		lz_divide_integer( &x, &y, &r );
		CheckTrue( "divide_integer( inf, 7 ) is infinite", is( r, true, false, false ), "" );

		x.Assign( 7L );
		infinity( y );
		long	l_result = lz_divide_integer( &x, &y, &r );
		CheckTrue( "divide_integer( 7, inf ) is nought", is( r, false, false, false ) && (l_result == kNoError), "" );
		CheckRendering( "divide_integer( 7, inf )", &r, "0" );
	}

	{
		CLZNumber	x;
		CLZNumber	y;

		x.Assign( 7L );
		nan( y );
		long	l_result = lz_divide( &x, &y, 10L );
		CheckTrue( "divide( 7, NaN, 10 ) is a NaN, no div0", is( x, false, true, false ) && (l_result == kNoError), "" );
	}

	{
		CLZNumber	x;
		CLZNumber	zero;

		x.Assign( 7L );
		zero.Assign( 0L );
		long	l_result = lz_divide( &x, &zero );
		CheckTrue( "7 / 0 in place is infinite", is( x, true, false, true ) && (l_result == kDivideByZero), "" );
	}

	{
		CLZNumber	x;
		CLZNumber	y;
		CLZNumber	r;

		x.Assign( "7.25" );
		infinity( y );
		lz_remainder( &x, &y, &r );
		CheckRendering( "remainder( 7.25, inf ) is 7.25", &r, "7.25" );
		lz_remainder_near( &x, &y, &r );
		CheckRendering( "remainder_near( 7.25, inf )", &r, "7.25" );

		infinity( x );
		y.Assign( 7L );
		lz_modulo( &x, &y, &r );
		CheckTrue( "modulo( inf, 7 ) is a NaN", is( r, false, true, false ), "" );

		x.Assign( 7L );
		nan( y );
		long	l_result = lz_remainder( &x, &y, &r );
		CheckTrue( "remainder( 7, NaN ) is a NaN, no div0", is( r, false, true, false ) && (l_result == kNoError), "" );
	}
}

static void TestRoundLeftOfThePoint()
{
	Banner( "Round to the left of the point, and the default mode" );

	/*
		**A negative count rounds to the left of the point**, and
		RM_round_default is the default mode -- since 2026-10-15. Round
		with a negative count answered a wrong magnitude, 1234.5 at -2
		giving 12, and the default mode truncated -- see
		CLZNumber::RoundAt.
	*/

	struct { const char* szValue; long lBase; long lPlaces; RoundingMode rmMode; const char* szExpected; } cases[] =
	{
		{ "1234.5",  10, -2, RM_round_half_even, "1200"  },
		{ "1250",    10, -2, RM_round_half_even, "1200"  },
		{ "1350",    10, -2, RM_round_half_even, "1400"  },
		{ "99.5",    10, -2, RM_round_half_even, "100"   },
		{ "-1234.5", 10, -2, RM_round_half_even, "-1200" },
		{ "5",       10, -3, RM_round_up,        "1000"  },
		{ "12",      10, -2, RM_round_half_even, "0"     },
		{ "FF.8",    16, -1, RM_round_half_even, "100"   },
		{ "123.456", 10,  1, RM_round_default,   "123.5" },
		{ "123.456", 10,  2, RM_round_default,   "123.46" },
	};

	for ( auto& theCase : cases )
	{
		CLZNumber	x;
		char		szLabel[ 72 ];

		SetInBase( x, theCase.lBase, theCase.szValue );
		x.Round( theCase.rmMode, theCase.lPlaces );
		::snprintf( szLabel, sizeof( szLabel ), "round( %s, %ld )", theCase.szValue, theCase.lPlaces );
		CheckRendering( szLabel, &x, theCase.szExpected );
	}

	// ...and arccos and arcsec of one are an exact nought, spelled so...
	CLZNumber	one;
	CLZNumber	r;

	one.Assign( 1L );
	lz_arccosine( &one, &r );
	CheckRendering( "arccos( 1 ) is a short nought", &r, "0" );
	lz_arcsecant( &one, &r );
	CheckRendering( "arcsec( 1 ) is a short nought", &r, "0" );
}

static void TestRootsAndHyperbolicsTakeTheShortRoad()
{
	Banner( "Square and cube roots, and the hyperbolics, take the short road" );

	/*
		**root_n of order two or three, and x to a half, went through a
		logarithm and an exponential until 2026-10-15**, though root_2 and
		root_3 exist -- 166 microseconds against 6.9 at 36 places; **and
		sinh and cosh took a second exponential for e^-x**, tanh and coth
		calling both. Counted rather than timed: a road is a fact, a time a
		measurement -- see lz_root_n, lz_power and lzi_sinh_cosh.
	*/

	struct { const char* szLabel; const char* szX; const char* szY; bool bPower; const char* szExpected; } roots[] =
	{
		{ "root_n( 2.25, 2 ) by root_2, short", "2.25", "2", false, "1.5" },
		{ "root_n( 2.5, 2 ) by root_2",         "2.5",  "2", false, "1.581138830084189665999446772216359267" },
		{ "root_n( 8, 3 ) by root_3, short",    "8",    "3", false, "2" },
		{ "root_n( -27, 3 ) by root_3",         "-27",  "3", false, "-3" },
		{ "power( 2.5, 0.5 ) by root_2",        "2.5",  "0.5", true, "1.581138830084189665999446772216359267" },
	};

	for ( auto& theCase : roots )
	{
		CLZNumber	x;
		CLZNumber	y;
		CLZNumber	r;

		x.Assign( theCase.szX );
		y.Assign( theCase.szY );

		long	l_before = lzi_logarithm_count();

		if ( theCase.bPower )
			lz_power( &x, &y, &r );
		else
			lz_root_n( &x, &y, &r );

		CheckTrue( theCase.szLabel, lzi_logarithm_count() == l_before, "no logarithm" );
		CheckRendering( theCase.szLabel, &r, theCase.szExpected );
	}

	// ...and a small cube root keeps the road whose spelling it has...
	{
		CLZNumber	x;
		CLZNumber	y;
		CLZNumber	r;

		x.Assign( "0.00002" );
		y.Assign( "3" );

		long	l_before = lzi_logarithm_count();

		lz_root_n( &x, &y, &r );
		CheckTrue( "root_n( 0.00002, 3 ) keeps its road", lzi_logarithm_count() > l_before, "written to W places" );
	}

	struct { const char* szLabel; long (*fn)( CLZNumber*, CLZNumber* ); } hyperbolics[] =
	{
		{ "sinh takes one exponential", lz_hypersine }, { "cosh takes one exponential", lz_hypercosine },
		{ "tanh takes one exponential", lz_hypertangent }, { "coth takes one exponential", lz_hypercotangent },
	};

	for ( auto& theCase : hyperbolics )
	{
		CLZNumber	x;
		CLZNumber	r;

		x.Assign( "0.3" );

		long	l_before = lzi_exponential_count();

		theCase.fn( &x, &r );

		char	szDetail[ 40 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld", lzi_exponential_count() - l_before );
		CheckTrue( theCase.szLabel, (lzi_exponential_count() - l_before) == 1, szDetail );
	}
}

static void TestExactnessByPlacesFirst()
{
	Banner( "A short answer's exactness, places first" );

	/*
		**A short answer was checked by raising its operands exactly to as
		much as the 4096th power until 2026-10-15** -- raise_n( x, 4000 )
		with an answer ending in a nought took 4.1 ms against 112
		microseconds beside it. Unequal places rule it out unworked now --
		see lzi_places_of_power. Counted: an answer whose places differ
		works no power out, and an exact one still does, and is short.
	*/

	CLZNumber	x;
	CLZNumber	y;
	CLZNumber	r;
	long		l_worked = 0;
	long		l_short = 0;
	long		l_tried = 0;

	// raise_n( 1.ddd, 4000 ) for x whose answer ends in a nought at 36 places...
	for ( long long i = 1000000000000000LL; (i < 1000000000000400LL) && (l_short < 4); i++ )
	{
		char	szX[ 40 ];

		::snprintf( szX, sizeof( szX ), "1.%lld", i );
		x.Assign( szX );
		y.Assign( 4000L );

		long	l_before = lzi_exactness_proofs();

		lz_raise_n( &x, &y, &r );

		char*	sz_r = r.GetStringRep();
		size_t	n_length = ::strlen( sz_r );
		bool	b_nought = (n_length > 0) && (sz_r[ n_length - 1 ] == '0');

		::free( sz_r );
		l_tried++;

		if ( b_nought )
		{
			l_short++;
			l_worked += lzi_exactness_proofs() - l_before;
		}
	}

	char	szDetail[ 80 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld ending in a nought, %ld worked", l_short, l_tried, l_worked );
	CheckTrue( "places rule raise_n( x, 4000 ) out", (l_short > 0) && (l_worked == 0), szDetail );

	// ...while an exact one is still worked, and short...
	x.Assign( "5.0625" );
	y.Assign( "0.25" );

	long	l_before = lzi_exactness_proofs();

	lz_power( &x, &y, &r );
	CheckTrue( "5.0625 to a quarter is worked", lzi_exactness_proofs() > l_before, "" );
	CheckRendering( "5.0625 to a quarter", &r, "1.5" );
}

static void TestLogarithmRoadBySeriesLength()
{
	Banner( "The logarithm's road, by its series' length" );

	/*
		**Past 150 digits Newton takes the logarithm only where the series
		would be long** -- since 2026-10-15; the series is sized from its y,
		so near one it is short and beats Newton at any width measured.
		Counted by exponentials, which Newton takes and the series does not
		-- see lz_log_e.
	*/

	struct { const char* szLabel; long lBase; const char* szX; long lDigits; bool bNewton; } cases[] =
	{
		{ "ln( 1 + 1.2E-7 ) at 200 digits, series", 10, "1.0000001234567", 200, false },
		{ "ln( 1 + 2^-10 ), base two, 300 digits",  2,  "1.000000000101",  300, false },
		{ "ln( 2.5 ), base two, 200 digits, Newton", 2,  "10.1",            200, true  },
	};

	for ( auto& theCase : cases )
	{
		long	l_saved = glWorkingPrecision;

		glWorkingPrecision = CLZNumber::EquivalentPlaces( theCase.lDigits, 10, theCase.lBase );

		CLZNumber	x;
		CLZNumber	r;

		SetInBase( x, theCase.lBase, theCase.szX );

		long	l_before = lzi_exponential_count();

		lz_log_e( &x, &r );

		long	l_taken = lzi_exponential_count() - l_before;

		glWorkingPrecision = l_saved;

		char	szDetail[ 40 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld exponentials", l_taken );
		CheckTrue( theCase.szLabel, theCase.bNewton ? (l_taken > 0) : (l_taken == 0), szDetail );
	}
}

static void TestWideDividendCut()
{
	Banner( "A wide dividend, cut to its quotient's places" );

	/*
		**A dividend wider than the answer is cut to the places the
		quotient is worked to** -- since 2026-10-16, performance item 7 of
		the review of 2026-10-15; see the cut in lzi_divide_work. The
		quotient was worked to every place of the dividend: 125 times the
		work for 16,000 places at 36.

		Pinned first, each answer checked against an exact reference: odd
		bases a hair from a half, where what the cut dropped decides --
		over half a unit of the last place kept, against a remainder one
		short of half the divisor (the first five), or level with it --
		and a divisor of one, whose division leaves nothing; and even
		bases rounding up from a cut that dropped a single unit. The last
		three divide by more than a limb, where the packed path's remainder
		is scaled. Each is taken packed and then a digit at a time.
	*/

	struct { const char* szFn; long lBase; long lPlaces; const char* szA; const char* szD; int iMode; const char* szExpect; } pins[] =
	{
		{ "divide", 11, 3, "32.7A555555555555555555555555556", "10", 1, "3.27A6" },
		{ "divide", 13, 3, "9B.1B666666666666666666666666667", "10", 1, "9.B1B7" },
		{ "divide", 15, 3, "43.60777777777777777777777777778", "10", 1, "4.3608" },
		{ "divide", 3, 3, "10.02111111111111111111111111112", "10", 1, "1.0022" },
		{ "divide", 5, 3, "12.44222222222222222222222222223", "10", 1, "1.2443" },
		{ "divide", 3, 3, "110.11111111111111111111111111112", "12", 101, "10" },
		{ "remainder_near", 3, 3, "110.11111111111111111111111111112", "12", -1, "-2.11111111111111111111111111111" },
		{ "remainder_near", 3, 3, "110.11111111111111111111111111120", "12", -1, "-2.11111111111111111111111111110" },
		{ "divide", 11, 3, "29.00100000000000000000000000001", "11", 5, "2.6466" },
		{ "divide", 13, 3, "17.64B00000000000000000000000001", "11", 2, "1.6047" },
		{ "divide", 15, 3, "40.AD700000000000000000000000001", "11", 5, "3.BDE8" },
		{ "divide", 11, 3, "5A28.55555555555555555555555555556", "1", 101, "5A29" },
		{ "divide", 13, 3, "1B88.66666666666666666666666666667", "1", 101, "1B89" },
		{ "divide", 10, 3, "17.35500000000000000000000000001", "10", 3, "1.7356" },
		{ "divide", 12, 3, "67.81600000000000000000000000001", "10", 3, "6.7817" },
		{ "divide", 14, 3, "61.92A70000000000000000000000001", "11", 3, "5.9D38" },
		{ "divide", 5, 3, "430.22222222222222222222222222223", "12", 101, "32" },
		{ "divide", 3, 3, "1121002110221201200021211.21021111111111111111111111111112", "1001201122122011010012020.01212", 102, "2" },
		{ "divide", 7, 3, "613316513456006433524155223.21164333333333333333333333333334", "6012655323210260133063626.50143", 102, "102" },
		{ "divide", 13, 3, "203330A845B537579549CA431255C3CB7.666666666666666666666666667", "B48BA45927C1740603769466068C72", 102, "241" },
	};

	long	l_saved = glWorkingPrecision;
	long	l_floor = lzi_set_packed_divide_floor( kPackedFloorMeasured );

	for ( int iWay = 0; iWay < 2; iWay++ )
	{
		lzi_set_packed_divide_floor( (iWay == 0) ? kPackedFloorMeasured : LONG_MAX );

		for ( auto& thePin : pins )
		{
			glWorkingPrecision = thePin.lPlaces;

			CLZNumber	x;
			CLZNumber	y;
			CLZNumber	r;

			SetInBase( x, thePin.lBase, thePin.szA );
			SetInBase( y, thePin.lBase, thePin.szD );
			r.SetBase( thePin.lBase, true );

			// Past a hundred, a target asking for a whole number in that mode...
			if ( thePin.iMode >= 100 )
			{
				r.SetRoundingDigit( RD_digit_zero );
				r.SetRoundingMode( static_cast<RoundingMode>( thePin.iMode - 100 ) );
			}
			else if ( thePin.iMode >= 0 )
				r.SetRoundingMode( static_cast<RoundingMode>( thePin.iMode ) );

			if ( ::strcmp( thePin.szFn, "remainder_near" ) == 0 )
				lz_remainder_near( &x, &y, &r );
			else
				lz_divide( &x, &y, &r );

			char	szLabel[ 80 ];

			::snprintf( szLabel, sizeof( szLabel ), "%s in base %ld, mode %d%s", thePin.szFn, thePin.lBase, thePin.iMode,
						(iWay == 0) ? "" : ", loop" );
			CheckRendering( szLabel, &r, thePin.szExpect );
		}
	}

	lzi_set_packed_divide_floor( l_floor );

	/*
		Then agreement with the whole dividend, in every base, packed and a
		digit at a time: D times m and a hair under a half in an odd base
		(the middle digit to every place) or over one in an even base, then
		a unit of the last place either side, in four modes.
	*/

	const char*	szDigits = "0123456789ABCDEF";
	long		l_tried = 0;
	long		l_differ = 0;

	for ( int iWay = 0; iWay < 2; iWay++ )
	{
		lzi_set_packed_divide_floor( (iWay == 0) ? kPackedFloorMeasured : LONG_MAX );

		for ( long lBase = 2; lBase <= 16; lBase++ )
		{
			long	l_places = (lBase == 2) ? 400 : 70;
			long	l_divisors[] = { 1, 2, 3, lBase + 1, (2 * lBase) + 1, (lBase * lBase) - 1 };

			std::string	s_half = "0.";
			std::string	s_unit = "0.";

			for ( long i = 0; i < l_places; i++ )
			{
				if ( (lBase % 2) != 0 )
					s_half += szDigits[ lBase / 2 ];
				else
					s_half += (i == 0) ? szDigits[ lBase / 2 ] : ((i == l_places - 1) ? '1' : '0');

				s_unit += (i == l_places - 1) ? '1' : '0';
			}

			for ( long lDivisor : l_divisors )
			{
				for ( long lWhole : { 1L, lBase + 3 } )
				{
					for ( int iSide = -1; iSide <= 1; iSide++ )
					{
						glWorkingPrecision = l_places + 20;

						CLZNumber	d;
						CLZNumber	q;
						CLZNumber	h;
						CLZNumber	u;
						CLZNumber	a;

						d.SetBase( lBase, true );
						d.Assign( lDivisor );
						q.SetBase( lBase, true );
						q.Assign( lWhole );
						SetInBase( h, lBase, s_half.c_str() );
						SetInBase( u, lBase, s_unit.c_str() );

						lz_add( &q, &h );
						lz_multiply( &d, &q, &a );

						if ( iSide < 0 )
							lz_subtract( &a, &u );
						else if ( iSide > 0 )
							lz_add( &a, &u );

						for ( long lWidth : { 3L, 12L } )
						{
							for ( int iMode : { 2, 1, 3, 102 } )
							{
								std::string	s_answer[ 2 ];

								for ( int iCut = 0; iCut < 2; iCut++ )
								{
									// The width first: a number takes it when it is made...
									glWorkingPrecision = lWidth;

									CLZNumber	r;

									r.SetBase( lBase, true );

									if ( iMode >= 100 )
									{
										r.SetRoundingDigit( RD_digit_zero );
										r.SetRoundingMode( static_cast<RoundingMode>( iMode - 100 ) );
									}
									else
										r.SetRoundingMode( static_cast<RoundingMode>( iMode ) );

									bool	b_previous = lzi_set_divide_cuts( iCut == 1 );

									lz_divide( &a, &d, &r );
									lzi_set_divide_cuts( b_previous );

									char*	sz_r = r.GetStringRep();
									char	szFlags[ 24 ];

									::snprintf( szFlags, sizeof( szFlags ), "|%lx", r.m_lFlags );
									s_answer[ iCut ] = std::string( (sz_r != nullptr) ? sz_r : "" ) + szFlags;
									::free( sz_r );
								}

								l_tried++;

								if ( s_answer[ 0 ] != s_answer[ 1 ] )
									l_differ++;
							}
						}
					}
				}
			}
		}
	}

	lzi_set_packed_divide_floor( l_floor );

	char	szDetail[ 80 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ", l_differ, l_tried );
	CheckTrue( "cut quotients are the whole dividend's", l_differ == 0, szDetail );

	// The work: a 2,000-place dividend at 36 places, cut as the engine
	//	stands and whole, and to a whole number, which keeps fewer...
	{
		std::string	s_wide = "3.";
		unsigned long long	l_state = 2026101601ULL;

		for ( int i = 0; i < 2000; i++ )
		{
			l_state = (l_state * 6364136223846793005ULL) + 1442695040888963407ULL;
			s_wide += szDigits[ (l_state >> 33) % 10 ];
		}

		std::string	s_answer[ 2 ];
		long		l_worked[ 2 ] = { 0, 0 };
		long		l_wholeWorked = 0;

		for ( int iCut = 0; iCut < 3; iCut++ )
		{
			glWorkingPrecision = 36;

			CLZNumber	x;
			CLZNumber	y;
			CLZNumber	r;

			x.Assign( s_wide.c_str() );
			y.Assign( 7L );

			// The whole dividend by the hook; the cut as the engine stands...
			bool	b_previous = (iCut == 0) ? lzi_set_divide_cuts( false ) : true;
			long	l_before = lzi_quotient_digits_worked();

			if ( iCut == 2 )
			{
				lz_divide_integer( &x, &y, &r );
				l_wholeWorked = lzi_quotient_digits_worked() - l_before;

				continue;
			}

			lz_divide( &x, &y, &r );
			l_worked[ iCut ] = lzi_quotient_digits_worked() - l_before;

			if ( iCut == 0 )
				lzi_set_divide_cuts( b_previous );

			char*	sz_r = r.GetStringRep();

			s_answer[ iCut ] = (sz_r != nullptr) ? sz_r : "";
			::free( sz_r );
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld digits worked, %ld whole", l_worked[ 1 ], l_worked[ 0 ] );
		CheckTrue( "2,000 places at 36 work the answer's", (l_worked[ 1 ] <= 100) && (l_worked[ 0 ] >= 2000) &&
				   (s_answer[ 0 ] == s_answer[ 1 ]), szDetail );

		// The units and the places a tail reads, and no fewer than pack...
		::snprintf( szDetail, sizeof( szDetail ), "%ld digits worked", l_wholeWorked );
		CheckTrue( "and to a whole number fewer", (l_wholeWorked >= 25) && (l_wholeWorked <= 40), szDetail );
	}

	/*
		**And the tail a rounding records is the whole dividend's**, for a
		function deciding whether its answer was in doubt -- which is why
		the cut keeps kRecordedTailDigits places past the answer. In base
		two the quotient's places past the 37th are a one, then noughts
		well past the 24 a tail reads, then a one: 0.5 exactly, which a
		place fewer kept would read as 0.5 and the sticky digit's 2^-24.
	*/

	for ( long lBase : { 2L, 3L, 10L } )
	{
		long	l_width = 36;
		long	l_beyond = (lBase == 2) ? 200 : 60;

		std::string	s_quotient = "1.";
		std::string	s_divisor = "1";
		unsigned long long	l_state = 2026101602ULL + lBase;

		for ( long i = 0; i < l_width + 1; i++ )
		{
			l_state = (l_state * 6364136223846793005ULL) + 1442695040888963407ULL;
			s_quotient += szDigits[ (l_state >> 33) % lBase ];
		}

		s_quotient += szDigits[ lBase / 2 ];

		for ( long i = 0; i < l_beyond; i++ )
			s_quotient += (i == l_beyond - 1) ? '1' : '0';

		for ( long i = 0; i < 63; i++ )
		{
			l_state = (l_state * 6364136223846793005ULL) + 1442695040888963407ULL;
			s_divisor += szDigits[ (l_state >> 33) % lBase ];
		}

		double	d_tail[ 2 ] = { -1.0, -1.0 };

		for ( int iCut = 0; iCut < 2; iCut++ )
		{
			glWorkingPrecision = l_width + l_beyond + 200;

			CLZNumber	q;
			CLZNumber	d;
			CLZNumber	a;

			SetInBase( q, lBase, s_quotient.c_str() );
			SetInBase( d, lBase, s_divisor.c_str() );
			lz_multiply( &d, &q, &a );

			// The width before the target: a division takes its target's...
			glWorkingPrecision = l_width;

			CLZNumber	r;

			r.SetBase( lBase, true );

			bool	b_previous = lzi_set_divide_cuts( iCut == 1 );

			CLZNumber::RecordRoundings( true );
			CLZNumber::ForgetLastRounding();
			lz_divide( &a, &d, &r );
			d_tail[ iCut ] = CLZNumber::LastDiscardedTail();
			CLZNumber::RecordRoundings( false );

			lzi_set_divide_cuts( b_previous );
		}

		char	szLabel[ 80 ];

		::snprintf( szLabel, sizeof( szLabel ), "the tail recorded in base %ld", lBase );
		::snprintf( szDetail, sizeof( szDetail ), "%.17g cut, %.17g whole", d_tail[ 1 ], d_tail[ 0 ] );
		CheckTrue( szLabel, (d_tail[ 1 ] == d_tail[ 0 ]) && ((lBase != 2) || (d_tail[ 1 ] == 0.5)), szDetail );
	}

	glWorkingPrecision = l_saved;
}

static void TestReciprocalAtPassWidth()
{
	Banner( "The reciprocal's passes, each at its own width" );

	/*
		**Each pass of the reciprocal multiplies by the divisor to its own
		width and a guard** -- since 2026-10-16, part (c) of performance
		item 7 of the review of 2026-10-15; see lzi_reciprocal_newton. It
		multiplied by the whole divisor every pass, so a long divisor over
		a short quotient paid for its length at every width. With the
		reciprocal forced: 1 over a 2,000-digit whole number at 100 places
		multiplies by no more than the quotient's width needs, a divisor
		shorter than that is taken as it stands, and every quotient is long
		division's to the last digit.
	*/

	long	l_saved = glWorkingPrecision;
	long	l_floor = lzi_set_reciprocal_floor( LONG_MAX );
	bool	b_shape = lzi_set_reciprocal_by_shape( false );

	const char*	szDigits = "0123456789ABCDEF";

	for ( long lBase : { 10L, 2L, 3L, 7L, 16L } )
	{
		for ( long lLength : { 2000L, 50L } )
		{
			unsigned long long	l_state = 2026101603ULL + (lBase * 7) + lLength;
			std::string			s_divisor;

			for ( long i = 0; i < lLength; i++ )
			{
				l_state = (l_state * 6364136223846793005ULL) + 1442695040888963407ULL;

				long	l_digit = static_cast<long>( (l_state >> 33) % lBase );

				s_divisor += szDigits[ ((i == 0) && (l_digit == 0)) ? 1 : l_digit ];
			}

			std::string	s_answer[ 2 ];
			long		l_widest = 0;

			for ( int iWay = 0; iWay < 2; iWay++ )
			{
				glWorkingPrecision = 100;

				CLZNumber	x;
				CLZNumber	y;
				CLZNumber	r;

				x.SetBase( lBase, true );
				x.Assign( "1" );
				SetInBase( y, lBase, s_divisor.c_str() );
				r.SetBase( lBase, true );

				// Long division, then the reciprocal forced...
				lzi_set_reciprocal_floor( (iWay == 0) ? LONG_MAX : 0 );
				lz_divide( &x, &y, &r );

				if ( iWay == 1 )
					l_widest = lzi_reciprocal_widest_divisor();

				char*	sz_r = r.GetStringRep();

				s_answer[ iWay ] = (sz_r != nullptr) ? sz_r : "";
				::free( sz_r );
			}

			char	szLabel[ 80 ];
			char	szDetail[ 80 ];

			::snprintf( szLabel, sizeof( szLabel ), "1/N, %ld digits, base %ld", lLength, lBase );
			::snprintf( szDetail, sizeof( szDetail ), "widest divisor %ld digits", l_widest );

			// The quotient's 100 places, the reciprocal's guard and the pass's,
			//	and the units: under 200; a short divisor as it stands...
			bool	b_width = (lLength > 200) ? ((l_widest > 0) && (l_widest < 200)) : (l_widest == lLength);

			CheckTrue( szLabel, b_width && (s_answer[ 0 ] == s_answer[ 1 ]) && !s_answer[ 0 ].empty(), szDetail );
		}
	}

	lzi_set_reciprocal_floor( l_floor );
	lzi_set_reciprocal_by_shape( b_shape );
	glWorkingPrecision = l_saved;
}

static void TestReciprocalBySignificantDigits()
{
	Banner( "The reciprocal chosen by significant digits" );

	/*
		**Division chooses the reciprocal by the quotient's significant
		digits and the divisor's** -- since 2026-10-16, part (b) of
		performance item 7 of the review of 2026-10-15; see
		lzi_reciprocal_floors. It counted the quotient's places, so 1 over
		a 9,000-digit whole number at 9,000 places, 18,001 places of which
		8,999 noughts, went the long way while 1 over the same digits as a
		fraction took the reciprocal -- the same work. And a long divisor
		over a short quotient, where the reciprocal wins by half again and
		more, never took it at all. And the other way: 1 over a 12,000-digit
		whole number at 100 places, counted 12,101 places, would take the
		reciprocal at half long division's speed; it has 102 significant
		digits and goes the long way. Only the reciprocal makes Karatsuba
		products, which tell the two paths apart; both answers are long
		division's.
	*/

	struct { const char* szLabel; long lPlaces; long lDivisorDigits; bool bWhole; bool bReciprocal; } cases[] =
	{
		{ "1/N, 9,000 digits at 9,000",		9000,	9000,	true,	true	},
		{ "1/0.N, the same, at 9,000",		9000,	9000,	false,	true	},
		{ "1,000 places over 12,000 digits",	1000,	12000,	false,	true	},
		{ "400 places over 12,000 digits",	400,	12000,	false,	false	},
		{ "1/N, 12,000 digits at 100",		100,	12000,	true,	false	},
	};

	long	l_saved = glWorkingPrecision;

	for ( auto& theCase : cases )
	{
		std::string	s_divisor = theCase.bWhole ? "" : "0.";

		for ( long i = 0; i < theCase.lDivisorDigits; i++ )
			s_divisor += "123456789"[ (i * 7 + 3) % 9 ];

		std::string	s_answer[ 2 ];
		long		l_products = 0;

		for ( int iPath = 0; iPath < 2; iPath++ )
		{
			// Long division, the shape rule off, then the engine's own choice...
			bool	b_shapeWas = lzi_set_reciprocal_by_shape( iPath == 1 );

			glWorkingPrecision = theCase.lPlaces;

			CLZNumber	one;
			CLZNumber	divisor;
			CLZNumber	answer;

			one.Assign( "1" );
			divisor.Assign( s_divisor.c_str() );

			long	l_digits = 0, l_packed = 0, l_transform = 0, l_before = 0, l_after = 0;

			lzi_multiply_tier_counts( &l_digits, &l_packed, &l_transform, &l_before );
			lz_divide( &one, &divisor, &answer );
			lzi_multiply_tier_counts( &l_digits, &l_packed, &l_transform, &l_after );

			lzi_set_reciprocal_by_shape( b_shapeWas );

			if ( iPath == 1 )
				l_products = l_after - l_before;

			char*	sz_r = answer.GetStringRep();

			s_answer[ iPath ] = (sz_r != nullptr) ? sz_r : "";
			::free( sz_r );
		}

		char	szDetail[ 72 ];

		::snprintf( szDetail, sizeof( szDetail ), "%ld Karatsuba products", l_products );
		CheckTrue( theCase.szLabel, ((l_products > 0) == theCase.bReciprocal) && (s_answer[ 0 ] == s_answer[ 1 ]) &&
				   !s_answer[ 0 ].empty(), szDetail );
	}

	glWorkingPrecision = l_saved;
}

static void TestSinCosBySplitting()
{
	Banner( "Sine and cosine by binary splitting" );

	/*
		**Cosine and sine together by binary splitting e^{iy} a chunk at a
		time** -- the kernel of performance item 1 of the review of
		2026-10-15, called by nothing yet; see lzi_sincos_by_splitting.
		Pinned against references worked in Decimal far past the width,
		rounded half-even: short arguments and ones as
		long as the width, a whole part, a negative, a hair from pi/2 in
		hex where the cosine is all noughts but its last digits. Then
		against lz_sine and lz_cosine, still on their series, at 300 places
		in four bases.
	*/

	struct { long lBase; long lPlaces; const char* szX; const char* szSin; const char* szCos; } pins[] =
	{
		{ 10, 40, "0.7", "0.6442176872376910536726143513987201830658", "0.7648421872844884262558599901918649092682" },
		{ 10, 200, "0.12345678901234567890123456789", "0.12314341519456258110057871646675665918358271074000354636417354430847365804554019524385766620763080181214099248915850822353764617559042815173015696192922490061845773431576303638821873022322874226127810", "0.99238888511219208113847370213201882296786535056845635732725986826567317414007825839789925732071508938914741422433236097662789799395372610638473328348839306828157791292885718166356362674664604441639810" },
		{ 7, 60, "1.3456", "0.666450424222356441651012552145441212306201343615525064162266", "0.020631000325313234101566556016315231453242054023130166115610" },
		{ 2, 120, "-0.1011011", "-0.101001110000110100100111001010100111011010101000110101001011011011011010000011101100100100000111000100101011101101110101", "0.110000011111101111101110111100111000000011100100111111111101110101011010011000010011111011001000011100100010111101100100" },
		{ 16, 80, "1.921FB54442D18469898CC51701B839A2", "0.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC487", "0.0000000000000000000000000000000052049C1114CF98E804177D4C76273644A29410F31C6809BC" },
		{ 3, 90, "0.0000121", "0.000012022222222001201100001022221220100220020122202010110010211010122102112101012101010020", "0.222222222110210000001020201012221211000222100000000221222020011221102111001112220020011022" },
		{ 12, 50, "-2.E", "-0.2814D071732DE5603580003E99578962E4D066737166026138", "-0.E84580DD55987E1977D3972556D417ED308D24849ED7339376" },
		{ 13, 64, "0.C0A8", "0.A5354019185A6728CC60BAAA1197A9624239710969AB20206B35065549CC5C4B", "0.7A4176C94477219C849B970AB1772A1500668572A870380C819312A705523162" },
	};

	long	l_saved = glWorkingPrecision;

	for ( auto& thePin : pins )
	{
		glWorkingPrecision = thePin.lPlaces;

		CLZNumber	x;
		CLZNumber	sine;
		CLZNumber	cosine;

		SetInBase( x, thePin.lBase, thePin.szX );
		sine.SetBase( thePin.lBase, true );
		cosine.SetBase( thePin.lBase, true );

		lzi_sincos_by_splitting( &x, &sine, &cosine );

		char	szLabel[ 80 ];

		::snprintf( szLabel, sizeof( szLabel ), "sin( %.12s ), base %ld at %ld", thePin.szX, thePin.lBase, thePin.lPlaces );
		CheckRendering( szLabel, &sine, thePin.szSin );
		::snprintf( szLabel, sizeof( szLabel ), "cos( %.12s ), base %ld at %ld", thePin.szX, thePin.lBase, thePin.lPlaces );
		CheckRendering( szLabel, &cosine, thePin.szCos );
	}

	// Against the series, digits every base here has: 2.1 reduces to
	//	pi - 2.1, as long as the width...
	long	l_tried = 0;
	long	l_differ = 0;
	long	l_before = lzi_sincos_splittings();

	for ( long lBase : { 10L, 7L, 16L, 3L } )
	{
		for ( const char* szX : { "0.21", "2.1", "-1.1" } )
		{
			glWorkingPrecision = 300;

			CLZNumber	x;
			CLZNumber	sine;
			CLZNumber	cosine;
			CLZNumber	series_sine;
			CLZNumber	series_cosine;

			x.SetBase( lBase, true );
			x.Assign( szX );
			sine.SetBase( lBase, true );
			cosine.SetBase( lBase, true );
			series_sine.SetBase( lBase, true );
			series_cosine.SetBase( lBase, true );

			lzi_sincos_by_splitting( &x, &sine, &cosine );
			lz_sine( &x, &series_sine );
			lz_cosine( &x, &series_cosine );

			l_tried += 2;

			if ( sine.Compare( &series_sine ) != kEqualTo )
				l_differ++;

			if ( cosine.Compare( &series_cosine ) != kEqualTo )
				l_differ++;
		}
	}

	char	szDetail[ 80 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, %ld runs", l_differ, l_tried, lzi_sincos_splittings() - l_before );
	CheckTrue( "the series' answers at 300 places", (l_differ == 0) && ((lzi_sincos_splittings() - l_before) == 12), szDetail );

	glWorkingPrecision = l_saved;
}

static void TestSinCosWired()
{
	Banner( "Sine and cosine take the kernel past their floor" );

	/*
		**Past its base's floor, lz_sine and lz_cosine take the kernel in
		place of the series** -- since 2026-10-16, stage 2 of performance
		item 1 of the review of 2026-10-15; see lzi_trig_series and
		lzi_sincos_split_floor. Under the caller's guard and none of its own,
		so it is checked inside them: with the floor forced to nought, the six
		circular functions answer as the series does in six bases -- short
		arguments and ones as long as the width, past pi/2, tiny, large
		enough for the reduction by 2 pi, either sign -- and an infinity or a
		not-a-number still answers a not-a-number. Then the dispatch at base
		ten's measured floor.
	*/

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	struct { const char* szName; Unary fn; } functions[] =
	{
		{ "sine", lz_sine }, { "cosine", lz_cosine }, { "tangent", lz_tangent },
		{ "secant", lz_secant }, { "cosecant", lz_cosecant }, { "cotangent", lz_cotangent },
	};

	const char*	szDigits = "0123456789ABCDEF";
	long		l_saved = glWorkingPrecision;
	long		l_floor = lzi_set_sincos_split_floor( kSinCosSplitFloorMeasured );
	long		l_tried = 0;
	long		l_differ = 0;

	for ( long lBase : { 2L, 3L, 7L, 10L, 13L, 16L } )
	{
		unsigned long long	l_state = 2026101604ULL + lBase;
		auto				digits = [ & ]( long lCount ) -> std::string
		{
			std::string	s;

			for ( long i = 0; i < lCount; i++ )
			{
				l_state = (l_state * 6364136223846793005ULL) + 1442695040888963407ULL;
				s += szDigits[ (l_state >> 33) % lBase ];
			}

			return s;
		};

		std::string	args[] =
		{
			"0." + digits( 3 ), "1." + digits( 4 ), "11." + digits( 2 ), "-0.000" + digits( 5 ),
			"0." + digits( 60 ), "-1." + digits( 60 ), digits( 1 ) + "11" + digits( 3 ) + "." + digits( 8 ),
		};

		for ( auto& s_arg : args )
		{
			for ( auto& theFunction : functions )
			{
				std::string	s_answer[ 2 ];

				for ( int iWay = 0; iWay < 2; iWay++ )
				{
					lzi_set_sincos_split_floor( (iWay == 0) ? LONG_MAX : 0 );
					glWorkingPrecision = 60;

					CLZNumber	x;
					CLZNumber	r;

					SetInBase( x, lBase, s_arg.c_str() );
					r.SetBase( lBase, true );
					theFunction.fn( &x, &r );

					char*	sz_r = r.GetStringRep();
					char	szFlags[ 24 ];

					::snprintf( szFlags, sizeof( szFlags ), "|%lx", r.m_lFlags );
					s_answer[ iWay ] = std::string( (sz_r != nullptr) ? sz_r : "" ) + szFlags;
					::free( sz_r );
				}

				l_tried++;

				if ( s_answer[ 0 ] != s_answer[ 1 ] )
					l_differ++;
			}
		}
	}

	char	szDetail[ 80 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ", l_differ, l_tried );
	CheckTrue( "the kernel forced, the series' answers", l_differ == 0, szDetail );

	// Flags reach the answer with the kernel forced: it is not taken for them...
	lzi_set_sincos_split_floor( 0 );

	for ( long lFlag : { static_cast<long>( kfInfinity ), static_cast<long>( kfQNaN ) } )
	{
		for ( Unary fn : { lz_sine, lz_secant } )
		{
			glWorkingPrecision = 60;

			CLZNumber	x;
			CLZNumber	r;

			x.Assign( 0L );
			x.SetFlag( lFlag );
			fn( &x, &r );

			CheckTrue( (lFlag == kfInfinity) ? "of an infinity, a not-a-number" : "of a not-a-number, one",
					   r.FlagIsOn( kfQNaN ), "" );
		}
	}

	// The dispatch: base ten takes the kernel at 500 places and not at 400...
	lzi_set_sincos_split_floor( kSinCosSplitFloorMeasured );

	long	l_runs[ 2 ] = { 0, 0 };
	long	l_widths[ 2 ] = { 500, 400 };

	for ( int iSide = 0; iSide < 2; iSide++ )
	{
		glWorkingPrecision = l_widths[ iSide ];

		CLZNumber	x;
		CLZNumber	r;

		x.Assign( "0.7" );

		long	l_before = lzi_sincos_splittings();

		lz_sine( &x, &r );
		l_runs[ iSide ] = lzi_sincos_splittings() - l_before;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld runs at 500, %ld at 400", l_runs[ 0 ], l_runs[ 1 ] );
	CheckTrue( "base ten's floor", (l_runs[ 0 ] > 0) && (l_runs[ 1 ] == 0), szDetail );

	/*
		**Tangent and cotangent take one run for both past the floor** --
		see lzi_sine_and_cosine -- where they took a cosine and a sine by
		separate calls, two runs. An argument folded to full width, and one
		under one, which the old road widened a place for.
	*/

	for ( Unary fn : { lz_tangent, lz_cotangent } )
	{
		for ( const char* szX : { "0.7", "2.4" } )
		{
			glWorkingPrecision = 600;

			CLZNumber	x;
			CLZNumber	r;

			x.Assign( szX );

			long	l_before = lzi_sincos_splittings();

			fn( &x, &r );

			char	szLabel[ 60 ];

			::snprintf( szLabel, sizeof( szLabel ), "%s( %s ) at 600, one run", (fn == lz_tangent) ? "tan" : "cot", szX );
			::snprintf( szDetail, sizeof( szDetail ), "%ld runs", lzi_sincos_splittings() - l_before );
			CheckTrue( szLabel, (lzi_sincos_splittings() - l_before) == 1, szDetail );
		}
	}

	// ...and a tiny argument, which the old road took its sine a place per
	//	nought wider for, as the series' answer...
	for ( const char* szX : { "0.000000000000000000000000000000000000000000000000000000000001234567",
							  "0.0000000000001234567" } )
	{
		std::string	s_answer[ 2 ];

		for ( int iWay = 0; iWay < 2; iWay++ )
		{
			lzi_set_sincos_split_floor( (iWay == 0) ? LONG_MAX : kSinCosSplitFloorMeasured );
			glWorkingPrecision = 600;

			CLZNumber	x;
			CLZNumber	r;

			x.Assign( szX );
			lz_tangent( &x, &r );

			char*	sz_r = r.GetStringRep();

			s_answer[ iWay ] = (sz_r != nullptr) ? sz_r : "";
			::free( sz_r );
		}

		lzi_set_sincos_split_floor( kSinCosSplitFloorMeasured );
		CheckTrue( "tan of a tiny argument, one run", s_answer[ 0 ] == s_answer[ 1 ], "" );
	}

	/*
		**A chunk is taken to the nearest**, so 0.6999...9 and a few digits
		costs what 0.7000...0 and the same digits does -- a chunk or two, and
		a tail. Cut, the nines made every chunk count, and the folds through
		pi/2 leave either by the sign of a rounding residue: secant and
		tangent of 0.7 at 3,000 places paid four times cosine's.
	*/

	{
		glWorkingPrecision = 600;

		std::string	s_nines = "0.6" + std::string( 590, '9' ) + "1234";
		std::string	s_noughts = "0.7" + std::string( 590, '0' ) + "1234";
		long		l_chunks[ 2 ] = { 0, 0 };

		for ( int i = 0; i < 2; i++ )
		{
			CLZNumber	x;
			CLZNumber	sine;
			CLZNumber	cosine;

			x.Assign( (i == 0) ? s_nines.c_str() : s_noughts.c_str() );
			lzi_sincos_by_splitting( &x, &sine, &cosine );
			l_chunks[ i ] = lzi_sincos_last_chunks();
		}

		::snprintf( szDetail, sizeof( szDetail ), "%ld chunks with the nines, %ld with the noughts", l_chunks[ 0 ], l_chunks[ 1 ] );
		CheckTrue( "a remainder either side costs alike", (l_chunks[ 0 ] <= 3) && (l_chunks[ 0 ] == l_chunks[ 1 ]), szDetail );
	}

	lzi_set_sincos_split_floor( l_floor );
	glWorkingPrecision = l_saved;
}

static void TestArcSineByNewton()
{
	Banner( "The inverse circular functions by Newton" );

	/*
		**Past its base's floor, lz_arcsine finds the angle of (sin t, cos t)
		by Newton on the sine/cosine kernel** -- since 2026-10-16, stage 3 of
		performance item 1 of the review of 2026-10-15; see
		lzi_angle_by_newton. Every inverse circular function comes down to
		lz_arcsine, so with the floor forced to nought all six answer as the
		halvings and the series do, in six bases: short arguments and long,
		tiny, a half, a hair from one, one itself, either sign, and past one
		for the arctangent, its cotangent and the two reciprocals. Outside
		the arcsine's domain it still refuses. Then the dispatch at base
		ten's measured floor.
	*/

	typedef long (*Unary)( CLZNumber*, CLZNumber* );

	struct { const char* szName; Unary fn; bool bPastOne; } functions[] =
	{
		{ "arcsine", lz_arcsine, false }, { "arccosine", lz_arccosine, false },
		{ "arctangent", lz_arctangent, true }, { "arccotangent", lz_arccotangent, true },
		{ "arcsecant", lz_arcsecant, true }, { "arccosecant", lz_arccosecant, true },
	};

	const char*	szDigits = "0123456789ABCDEF";
	long		l_saved = glWorkingPrecision;
	long		l_floor = lzi_set_arcsine_newton_floor( kArcSineNewtonFloorMeasured );
	long		l_tried = 0;
	long		l_differ = 0;

	for ( long lBase : { 2L, 3L, 7L, 10L, 13L, 16L } )
	{
		unsigned long long	l_state = 2026101605ULL + lBase;
		auto				digits = [ & ]( long lCount ) -> std::string
		{
			std::string	s;

			for ( long i = 0; i < lCount; i++ )
			{
				l_state = (l_state * 6364136223846793005ULL) + 1442695040888963407ULL;
				s += szDigits[ (l_state >> 33) % lBase ];
			}

			return s;
		};

		std::string	s_top( 8, szDigits[ lBase - 1 ] );
		std::string	inside[] =
		{
			"0." + digits( 4 ), "0.0000" + digits( 5 ), std::string( "0." ) + szDigits[ lBase / 2 ],
			"0." + s_top + digits( 3 ), "-0." + digits( 60 ), "1",
		};
		std::string	past[] = { "1" + digits( 2 ) + "." + digits( 5 ), "-1." + digits( 60 ) };

		for ( auto& theFunction : functions )
		{
			std::vector<std::string>	args;

			// The arcsine's and arccosine's domain, or past one for the rest...
			if ( theFunction.bPastOne )
			{
				args.assign( std::begin( past ), std::end( past ) );

				if ( (theFunction.fn == lz_arctangent) || (theFunction.fn == lz_arccotangent) )
					args.insert( args.end(), std::begin( inside ), std::end( inside ) );
			}
			else
				args.assign( std::begin( inside ), std::end( inside ) );

			for ( auto& s_arg : args )
			{
				std::string	s_answer[ 2 ];

				for ( int iWay = 0; iWay < 2; iWay++ )
				{
					lzi_set_arcsine_newton_floor( (iWay == 0) ? LONG_MAX : 0 );
					glWorkingPrecision = 60;

					CLZNumber	x;
					CLZNumber	r;

					SetInBase( x, lBase, s_arg.c_str() );
					r.SetBase( lBase, true );
					theFunction.fn( &x, &r );

					char*	sz_r = r.GetStringRep();
					char	szFlags[ 24 ];

					::snprintf( szFlags, sizeof( szFlags ), "|%lx", r.m_lFlags );
					s_answer[ iWay ] = std::string( (sz_r != nullptr) ? sz_r : "" ) + szFlags;
					::free( sz_r );
				}

				l_tried++;

				if ( s_answer[ 0 ] != s_answer[ 1 ] )
					l_differ++;
			}
		}
	}

	char	szDetail[ 80 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ", l_differ, l_tried );
	CheckTrue( "Newton forced, the series' answers", l_differ == 0, szDetail );

	// Outside the domain the refusal comes first, Newton forced or not...
	{
		lzi_set_arcsine_newton_floor( 0 );
		glWorkingPrecision = 60;

		CLZNumber	x;
		CLZNumber	r;

		x.Assign( "1.5" );

		long	l_before = lzi_arcsine_newtons();
		long	l_code = lz_arcsine( &x, &r );

		CheckTrue( "arcsine( 1.5 ) still refused", (l_code != kNoError) && r.FlagIsOn( kfQNaN ) &&
				   (lzi_arcsine_newtons() == l_before), "" );
	}

	// The dispatch at base ten's floor...
	lzi_set_arcsine_newton_floor( kArcSineNewtonFloorMeasured );

	long	l_runs[ 2 ] = { 0, 0 };
	long	l_widths[ 2 ] = { 300, 200 };

	for ( int iSide = 0; iSide < 2; iSide++ )
	{
		glWorkingPrecision = l_widths[ iSide ];

		CLZNumber	x;
		CLZNumber	r;

		x.Assign( "0.7" );

		long	l_before = lzi_arcsine_newtons();

		lz_arcsine( &x, &r );
		l_runs[ iSide ] = lzi_arcsine_newtons() - l_before;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld runs at %ld, %ld at %ld", l_runs[ 0 ], l_widths[ 0 ], l_runs[ 1 ], l_widths[ 1 ] );
	CheckTrue( "base ten's floor", (l_runs[ 0 ] > 0) && (l_runs[ 1 ] == 0), szDetail );

	lzi_set_arcsine_newton_floor( l_floor );
	glWorkingPrecision = l_saved;
}

static void TestBaseConversionByProducts()
{
	Banner( "Base conversion by products" );

	/*
		**Past a floor, SetBase converts by products** -- performance item
		8(g) of the review of 2026-10-15; see CLZNumber::BaseConvertByProducts.
		It has to write the loops' answer to the digit, so first pins worked
		exactly, taken by products with the floor at
		nought: ties of both parities into an even base and an odd one, just
		under and just over the half from odd bases, a carry into the whole
		part, a whole number, a negative, noughts after the point, a fraction
		that ends, and one wide in both parts. Then the products against the
		loops on every pair of bases, the digit array itself compared; then
		the measured floor.
	*/

	struct { long lFrom; long lTo; long lPlaces; const char* szValue; const char* szWanted; const char* szLabel; } pins[] =
	{
		{ 10, 2, 40, "0.1", "0.0001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011", "a tenth" },
		{ 10, 2, 12, "0.64499887517149545601569116115570068359375", "0.101001010001111010100101011100101101001", "a tie, kept digit odd" },
		{ 10, 2, 12, "37.85841679015766203519888222217559814453125", "100101.110110111100000100110011111010000001", "a tie, kept digit even" },
		{ 10, 5, 20, "0.01352589670721182276383670272", "0.00132113234041421403042120334", "an odd base tie, odd" },
		{ 10, 5, 20, "4.82068573463882465519662530560", "4.40224310140304242404102141232", "an odd base tie, even" },
		{ 12, 3, 30, "0.51593429D374491E15D794D45D89355032356738701D599554260469258816529940", "0.10211202020021210221002001212012201212211210212210221100011200210122", "base 12 to 3, a tie" },
		{ 3, 10, 40, "0.202220010210210220121012200002102012200220220", "0.77385600945945412446", "base 3, just under half" },
		{ 3, 10, 40, "0.202120101201000002010102222011022012110210201", "0.76186068240580148968", "base 3, just over half" },
		{ 7, 10, 30, "5.2310652106406615413244461441655501", "5.35025634614637459509859649", "base 7, just over half" },
		{ 10, 16, 30, "0.999999999999999999999999999999999999999999999999999999999999", "1", "a carry into the whole part" },
		{ 16, 10, 0, "6BAA6C437C3C2DED2C3E16D6CBF0D421B830A3427FB8C304FDF9E7913202BAB334B0472694938DCF", "898330352830542402437692178098561166836774501807233699500349985825130305697849243756644504145359", "a whole number" },
		{ 10, 7, 60, "-1160581236102961140135101223061142283355.54080883972092210670803723199765978131934131718702", "-13563061563663306523303516053016646302255325305.353332422244262661053106124035456123626566514112224613524131052312521456", "a negative" },
		{ 2, 10, 50, "0.000000000000000000000000000000111100101111111010000011101001011011101100000010111011000011110000111001100000000010101000", "0.0000000008840077", "noughts after the point" },
		{ 16, 2, 400, "0.F0091D628F9A08AAB6EEF1AB227747", "0.111100000000100100011101011000101000111110011010000010001010101010110110111011101111000110101011001000100111011101000111", "hex to binary, ending" },
		{ 10, 13, 200, "702957254584469521477598612264042606653329121953980843710845945160732901761756649035540997471284770416631019080759897636139241299716166804690470128325.779866620823728561787263241084293065844406872134437422189853354507023469581191197625585412175645235293266104396982482739249237213619221333606894911299194710833475895944856740903376552845353464258056403239226323173312132577747926993034137143258972889981486377990528173899373392076444049427757406093595", "3A344540CB907B6372762CA1211118B3312A22629873652C7534160C23902254C8104269745887C921C3C7089467329A3216677C4979A4BB86B9628A565320231297AB9.A1A4A02BB719AC30C03530059218A9A477B5A58A2997A4686337265B8779932464C8409C67277328B3368A17936B7913127301756873266B573431B0530CB1C184C16077819369363A882B1035848B877CB1607C691509856302", "wide, both parts" },
	};

	long	l_saved = glWorkingPrecision;
	long	l_floor = CLZNumber::SetConversionFloor( 0 );
	long	l_before = CLZNumber::ConversionsByProducts();
	long	l_pins = 0;

	for ( auto& thePin : pins )
	{
		glWorkingPrecision = 1000;

		CLZNumber	x;

		SetInBase( x, thePin.lFrom, thePin.szValue );
		x.SetPrecision( thePin.lPlaces );
		x.SetBase( thePin.lTo, true );

		char	szLabel[ 80 ];

		::snprintf( szLabel, sizeof( szLabel ), "%s, %ld to %ld", thePin.szLabel, thePin.lFrom, thePin.lTo );
		CheckRendering( szLabel, &x, thePin.szWanted );
		l_pins++;
	}

	char	szDetail[ 120 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld", CLZNumber::ConversionsByProducts() - l_before, l_pins );
	CheckTrue( "every pin by products", (CLZNumber::ConversionsByProducts() - l_before) == l_pins, szDetail );

	// The loops' answer, the digit array itself, on every pair of bases:
	//	a whole number, a fraction, both, a short fraction at a wide budget
	//	-- ending, where the bases allow -- and a long one at a narrow; a
	//	fraction every other time with no whole part written at all...
	struct { long lWhole; long lFraction; long lPlaces; } shapes[] =
	{
		{ 150, 0, 0 }, { 0, 150, 150 }, { 60, 120, 120 }, { 2, 3, 400 }, { 0, 200, 20 },
	};

	std::uint64_t	l_seed = 2026;
	long			l_tried = 0;
	long			l_differ = 0;

	l_before = CLZNumber::ConversionsByProducts();

	auto	digits = [ &l_seed ]( long lBase, long lCount, bool bLead )
	{
		std::string	s_digits;

		for ( long i = 0; i < lCount; i++ )
		{
			l_seed = (l_seed * 6364136223846793005ULL) + 1442695040888963407ULL;

			long	l_digit = static_cast<long>( (l_seed >> 33) % static_cast<std::uint64_t>( lBase ) );

			if ( bLead && (i == 0) && (l_digit == 0) )
				l_digit = 1;

			s_digits += "0123456789ABCDEF"[ l_digit ];
		}

		return s_digits;
	};

	for ( long lFrom = 2; lFrom <= 16; lFrom++ )
	{
		for ( long lTo = 2; lTo <= 16; lTo++ )
		{
			if ( lFrom == lTo )
				continue;

			for ( auto& theShape : shapes )
			{
				std::string	s_value = ((l_tried % 3) == 0) ? "-" : "";

				// No whole part written at all, every other time: the loops
				//	keep none for it, where "0." keeps the one digit...
				s_value += (theShape.lWhole > 0) ? digits( lFrom, theShape.lWhole, true ) : (((l_tried % 2) == 0) ? "0" : "");

				if ( theShape.lFraction > 0 )
					s_value += "." + digits( lFrom, theShape.lFraction, false );

				glWorkingPrecision = 1000;

				CLZNumber	by_products;
				CLZNumber	by_loops;

				SetInBase( by_products, lFrom, s_value.c_str() );
				SetInBase( by_loops, lFrom, s_value.c_str() );
				by_products.SetPrecision( theShape.lPlaces );
				by_loops.SetPrecision( theShape.lPlaces );

				CLZNumber::SetConversionFloor( 0 );
				by_products.SetBase( lTo, true );
				CLZNumber::SetConversionFloor( LONG_MAX );
				by_loops.SetBase( lTo, true );

				bool	b_same = (by_products.m_lDigitsLength == by_loops.m_lDigitsLength) &&
								 (by_products.m_lRadix == by_loops.m_lRadix) &&
								 (by_products.m_lPrecision == by_loops.m_lPrecision) &&
								 (by_products.m_lFlags == by_loops.m_lFlags) &&
								 (by_products.m_lSign == by_loops.m_lSign);

				for ( long i = 0; b_same && (i < by_products.m_lDigitsLength); i++ )
					b_same = by_products.m_lpDigits[ i ] == by_loops.m_lpDigits[ i ];

				l_tried++;

				if ( !b_same )
					l_differ++;
			}
		}
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ, %ld by products", l_differ, l_tried, CLZNumber::ConversionsByProducts() - l_before );
	CheckTrue( "the loops' digits, every pair", (l_differ == 0) && (l_tried == 1050) &&
			   ((CLZNumber::ConversionsByProducts() - l_before) == l_tried), szDetail );

	// The measured floor: by the loops, a tenth at 36 places; a short
	//	fraction at 3,000, whose work passes the floor and is short of four
	//	per digit written; and a whole number of thirty digits, the other
	//	way round. By products, a value of 3,000 places, with the loops'
	//	answer...
	CLZNumber::SetConversionFloor( CLZNumber::kConversionFloorMeasured );

	long	l_narrow = CLZNumber::ConversionsByProducts();

	struct { long lPlaces; long lTo; const char* szValue; } narrows[] =
	{
		{ 36, 2, "0.1" }, { 3000, 16, "0.075" }, { 36, 2, "123456789012345678901234567890" },
	};

	for ( auto& theNarrow : narrows )
	{
		glWorkingPrecision = theNarrow.lPlaces;

		CLZNumber	x;

		x.Assign( theNarrow.szValue );
		x.SetBase( theNarrow.lTo, true );
	}

	l_narrow = CLZNumber::ConversionsByProducts() - l_narrow;

	std::string	s_wide = digits( 10, 1500, true ) + "." + digits( 10, 3000, false );

	glWorkingPrecision = 3000;

	CLZNumber	wide;
	CLZNumber	wide_loops;

	wide.Assign( s_wide.c_str() );
	wide_loops.Assign( s_wide.c_str() );

	long	l_wide = CLZNumber::ConversionsByProducts();

	wide.SetBase( 2, true );
	l_wide = CLZNumber::ConversionsByProducts() - l_wide;

	CLZNumber::SetConversionFloor( LONG_MAX );
	wide_loops.SetBase( 2, true );

	::snprintf( szDetail, sizeof( szDetail ), "%ld of the narrow three by products, %ld at 3,000", l_narrow, l_wide );
	CheckTrue( "the measured floor", (l_narrow == 0) && (l_wide == 1) && (wide.Compare( &wide_loops ) == kEqualTo) &&
			   (wide.m_lRadix == wide_loops.m_lRadix), szDetail );

	CLZNumber::SetConversionFloor( l_floor );
	glWorkingPrecision = l_saved;
}

static void TestCubeRootPassesAtEachWidth()
{
	Banner( "The cube root's passes, each at its width" );

	/*
		**The cube root turns three buffers and set the width on two** --
		until 2026-10-17. The third kept the full width it was made with,
		and every pass that wrote into it answered at that width, so a
		narrow stage wrote an iterate thousands of digits long and then
		took it, never equal to a narrow one, for a cycle it could not
		see: the stage ran to its cap. Counted rather than timed, by
		lzi_root_3_passes: the cube root of 2 in base three at 300 digits
		took 67 passes and in base seven at 10,000 digits 71, where its
		neighbours take a dozen to twenty. The answer was right either
		way, and stays so; the first is pinned against an exact integer
		cube root.
	*/

	long	l_saved = glWorkingPrecision;

	struct { long lBase; long lDigits; long lMost; } cases[] = { { 3, 300, 25 }, { 7, 10000, 25 } };

	for ( auto& theCase : cases )
	{
		glWorkingPrecision = CLZNumber::EquivalentPlaces( theCase.lDigits, 10, theCase.lBase );

		CLZNumber	x;
		CLZNumber	r;

		x.SetBase( theCase.lBase, false );
		x.Assign( "2" );
		r.SetBase( theCase.lBase, false );

		long	l_before = lzi_root_3_passes();

		lz_root_3( &x, &r );

		long	l_passes = lzi_root_3_passes() - l_before;
		char	szLabel[ 80 ];
		char	szDetail[ 40 ];

		::snprintf( szLabel, sizeof( szLabel ), "cbrt( 2 ), base %ld at %ld digits", theCase.lBase, theCase.lDigits );
		::snprintf( szDetail, sizeof( szDetail ), "%ld passes", l_passes );
		// At least a few, so a counter that never counted cannot pass...
		CheckTrue( szLabel, (l_passes >= 5) && (l_passes <= theCase.lMost), szDetail );

		if ( theCase.lBase == 3 )
			CheckRendering( "and its digits, base 3 at 300", &r, "1.02100011100020022202021011011201122210110121022200212121100222022210012222002222011201201112220012011200100200010111020011120122222010100102102202110022212122001100022001020010102222011011022220201211201200021101200111021122222011000120020012121210022010101201112211202011112021021012020200111010222221000112101120102002111200211122110002112011020000122000021102222101120112210200020002001021212120010110022200100210220112120112211021002120020020120021211011021212210020111212220212220112120201021220102122120000200211121021122011210021020110021001221120222101020002021120000112112001020100120002011002112202221011022210220222100" );
	}

	glWorkingPrecision = l_saved;
}

static void TestCubeRootByReciprocal()
{
	Banner( "The cube root by its reciprocal" );

	/*
		**Past each base's floor the cube root multiplies where it divided**
		-- since 2026-10-17; see lzi_root_3_inverse_pass. First three near
		ties, a cube root a few units of the guard place from a half, taken
		by the reciprocal with the floor at nought and pinned against the
		answers the root's own road gives, which a census of 7,200 such
		scored exactly; with two guard places, where four are kept, one of
		them goes wrong. Then the two roads on full-width operands in odd
		bases and base twelve, which caught a reciprocal made in base ten;
		then the measured floor.
	*/

	struct { long lBase; long lPlaces; const char* szX; const char* szWanted; } ties[] =
	{
		{ 14, 88, "612.0B47703816A9826B42B512B8C5004525DADC42DD5A209060AD709A65493C2C8B38517A07983C1A9C69D090BBCC990257425D0323699A01354CA908D7A96141C27ABB5A55DC8185880A79A0AD710C666983123C44419839D589C4D173C0395B10C928B2A6A8BB7B6CA866A005D57622B87C29D48658811B00C494C4D3B3490983079A4461C4BABA3A41",
		  "A.862C4172721BA053C22648099DB169C14D10674801BA817B254155893C08D500854CA8B7931991D095927A91" },
		{ 15, 86, "DDB.1C96515412E2275693AE78947A9D2D25B7A2D1A75DCD43A252B369329E5936252AAB618D47A864018D3007ED427B9805BAD64E53531B741591D892BD33D970196744B4BB4B750801BE38C98E535E175252C7EEA901B3B7705CA8964E305910E6182E7009B8423A3714031C09C87C54BE3969E12C14150BD6B42237643ECA1EE3D7D751C3E7A313",
		  "E.96A00D692B43E42848B2BAB750E1AE681137619D5832C19CE3AA24BAAB1D69BE6593B5318877E12C1DC2B8" },
		{ 16, 84, "B46.B1D4B8AEF835AAB7474F738E1A09C1AFE2C985E43288030C1E4B167945B4BEE6AAB5B8D5AB5076DEB2DC03A907D2D2687469FE17CF5E219528F6020643B6B162850BD9862531DA1F857B9E4729B6B77A48CEB3E8254B84E235D5BE37C48C56D1D827E5EEBE71A22186B55BD03E4B62FBBCFD3AF11C0EA19B580E6463374D14AFF660918847EA",
		  "E.3D142872A1D4C0F50E84C1DE14B99A3384D7821ADA3D605AABC4F8A8776523FE455CD98106147291CF92" },
	};

	long	l_saved = glWorkingPrecision;
	long	l_floor = lzi_set_root_3_inverse_floor( 0 );
	long	l_before = lzi_root_3_inverse_passes();
	long	l_wrongAtTwo = 0;

	for ( auto& theTie : ties )
	{
		for ( int iGuard = 0; iGuard < 2; iGuard++ )
		{
			long	l_guard = lzi_set_root_3_inverse_guard( (iGuard == 0) ? -1 : 2 );

			glWorkingPrecision = theTie.lPlaces;

			CLZNumber	x;
			CLZNumber	r;

			x.SetBase( theTie.lBase, false );
			x.Assign( theTie.szX );
			r.SetBase( theTie.lBase, false );
			lz_root_3( &x, &r );

			lzi_set_root_3_inverse_guard( l_guard );

			if ( iGuard == 0 )
			{
				char	szLabel[ 80 ];

				::snprintf( szLabel, sizeof( szLabel ), "a near tie in base %ld at %ld", theTie.lBase, theTie.lPlaces );
				CheckRendering( szLabel, &r, theTie.szWanted );
			}
			else
			{
				char*	sz_r = r.GetStringRep();

				if ( (sz_r == nullptr) || (::strcmp( sz_r, theTie.szWanted ) != 0) )
					l_wrongAtTwo++;

				::free( sz_r );
			}
		}
	}

	char	szDetail[ 80 ];

	::snprintf( szDetail, sizeof( szDetail ), "%ld reciprocal passes", lzi_root_3_inverse_passes() - l_before );
	CheckTrue( "by the reciprocal", (lzi_root_3_inverse_passes() - l_before) > 0, szDetail );
	::snprintf( szDetail, sizeof( szDetail ), "%ld of 3 wrong", l_wrongAtTwo );
	CheckTrue( "two guard places are not enough", l_wrongAtTwo > 0, szDetail );

	// The two roads alike on full-width operands in odd bases and twelve...
	long			l_tried = 0;
	long			l_differ = 0;
	std::uint64_t	l_seed = 33;

	for ( long lBase : { 3L, 5L, 7L, 9L, 12L, 13L, 15L } )
	{
		for ( long lDigits : { 150L, 400L } )
		{
			long		l_places = CLZNumber::EquivalentPlaces( lDigits, 10, lBase );
			std::string	s_x = "1.";

			for ( long i = 0; i < l_places; i++ )
			{
				l_seed = (l_seed * 6364136223846793005ULL) + 1442695040888963407ULL;
				s_x += "0123456789ABCDEF"[ (l_seed >> 33) % static_cast<std::uint64_t>( lBase ) ];
			}

			std::string	s_answer[ 2 ];

			for ( int iWay = 0; iWay < 2; iWay++ )
			{
				lzi_set_root_3_inverse_floor( (iWay == 0) ? 0 : LONG_MAX );
				glWorkingPrecision = l_places;

				CLZNumber	x;
				CLZNumber	r;

				x.SetBase( lBase, false );
				x.Assign( s_x.c_str() );
				r.SetBase( lBase, false );
				lz_root_3( &x, &r );

				char*	sz_r = r.GetStringRep();

				s_answer[ iWay ] = (sz_r != nullptr) ? sz_r : "";
				::free( sz_r );
			}

			l_tried++;

			if ( s_answer[ 0 ] != s_answer[ 1 ] )
				l_differ++;
		}
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of %ld differ", l_differ, l_tried );
	CheckTrue( "the two roads alike", l_differ == 0, szDetail );

	// The measured floor: base ten from 600 places, sixteen from its first
	//	width with a second stage...
	lzi_set_root_3_inverse_floor( kRoot3InverseFloorMeasured );

	struct { long lBase; long lPlaces; bool bReciprocal; } floors[] =
	{
		{ 10, 1000, true }, { 10, 300, false }, { 16, 100, true }, { 2, 3000, false },
	};

	long	l_right = 0;

	for ( auto& theFloor : floors )
	{
		glWorkingPrecision = theFloor.lPlaces;

		CLZNumber	x;
		CLZNumber	r;

		x.SetBase( theFloor.lBase, false );
		x.Assign( 7L );
		r.SetBase( theFloor.lBase, false );

		long	l_passes = lzi_root_3_inverse_passes();

		lz_root_3( &x, &r );

		if ( (lzi_root_3_inverse_passes() > l_passes) == theFloor.bReciprocal )
			l_right++;
	}

	::snprintf( szDetail, sizeof( szDetail ), "%ld of 4 on their road", l_right );
	CheckTrue( "the measured floors", l_right == 4, szDetail );

	lzi_set_root_3_inverse_floor( l_floor );
	glWorkingPrecision = l_saved;
}

static void TestAgmRootRoundedToItsWidth()
{
	Banner( "The AGM's square root, rounded to its width" );

	/*
		**The AGM's square root came back wider than it was asked in odd
		bases** -- until 2026-10-17. Its last step halves a correction by a
		division there, a half never ending, and a division answers wider
		than it was asked; the sum kept every place, up to half again as
		many as asked. Answers never moved, and lz_root_2 rounds its own,
		but the logarithm by the AGM multiplies the root into its next
		cycle. Read by lzi_root_core_overrun: these four took 137 to 189
		places past their width; base ten is the control. Then a logarithm
		by the AGM in base nine, the most over any of its cycles' roots.
	*/

	long	l_saved = glWorkingPrecision;

	struct { long lBase; long lPlaces; const char* szX; } roots[] =
	{
		{ 7, 356, "2" }, { 9, 315, "2" }, { 15, 256, "2" }, { 7, 356, "3.1416" }, { 10, 300, "2" },
	};

	for ( auto& theRoot : roots )
	{
		glWorkingPrecision = theRoot.lPlaces;

		CLZNumber	x;
		CLZNumber	r;

		x.SetBase( theRoot.lBase, false );
		x.Assign( theRoot.szX );
		r.SetBase( theRoot.lBase, false );

		long	l_runs = lzi_root_core_runs();

		lzi_root_core_overrun();
		lz_root_2( &x, &r );

		long	l_overrun = lzi_root_core_overrun();

		char	szLabel[ 80 ];
		char	szDetail[ 60 ];

		::snprintf( szLabel, sizeof( szLabel ), "sqrt( %s ), base %ld at %ld", theRoot.szX, theRoot.lBase, theRoot.lPlaces );
		::snprintf( szDetail, sizeof( szDetail ), "%ld places over, %ld core runs", l_overrun, lzi_root_core_runs() - l_runs );
		CheckTrue( szLabel, ((lzi_root_core_runs() - l_runs) == 1) && (l_overrun > LONG_MIN) && (l_overrun <= 0), szDetail );
	}

	{
		glWorkingPrecision = 1048;

		CLZNumber	x;
		CLZNumber	r;

		x.SetBase( 9, false );
		x.Assign( "2" );
		r.SetBase( 9, false );

		long	l_cycles = lzi_agm_log_cycles();

		lzi_root_core_overrun();
		lz_log_e( &x, &r );

		long	l_overrun = lzi_root_core_overrun();
		char	szDetail[ 60 ];

		// Every cycle's root, the most of them...
		::snprintf( szDetail, sizeof( szDetail ), "%ld places over, %ld AGM cycles", l_overrun, lzi_agm_log_cycles() - l_cycles );
		CheckTrue( "ln( 2 ) by the AGM, base 9 at 1048", ((lzi_agm_log_cycles() - l_cycles) > 0) && (l_overrun > LONG_MIN) && (l_overrun <= 0), szDetail );
	}

	glWorkingPrecision = l_saved;
}

static void TestLogNewtonAtItsWidths()
{
	Banner( "The logarithm by Newton, each stage at its width" );

	/*
		**Each stage but the last at its own width** -- until 2026-10-17
		the operand went into every stage at its full length, and the step,
		a small quotient answering to the width's significant digits and
		its noughts, was added in unrounded, so the iterate came out half
		again as wide as its stage. Read by lzi_log_newton_overrun, with
		Newton forced and the AGM held off: one digit, which only the
		iterate can overrun; an operand three widths long; and one under
		one, which goes in as its reciprocal at the full width -- one that
		never ends, where 0.0001's is exact and short. Each answer is held
		to the series' as well.
	*/

	long	l_saved = glWorkingPrecision;
	long	l_newtonWas = lzi_set_log_newton_from( 200 );
	long	l_agmWas = lzi_set_log_agm_floor( LONG_MAX );
	const char*	sz_digits = "0123456789ABCDEF";

	for ( long lBase : { 10L, 3L, 7L, 2L } )
	{
		long	l_places = CLZNumber::EquivalentPlaces( 300, 10, lBase );

		// Three widths of digits after "1.", from a fixed sequence...
		std::string		s_long = "1.";
		std::uint64_t	ul_seed = 7;

		for ( long i = 0; i < (3 * l_places); i++ )
		{
			ul_seed = (ul_seed * 6364136223846793005ULL) + 1442695040888963407ULL;
			s_long += sz_digits[ (ul_seed >> 33) % static_cast<std::uint64_t>( lBase ) ];
		}

		// Two is the base itself in base two, whose logarithm is the constant's;
		//	and the one under one has a reciprocal that never ends, which a
		//	power of the base would not...
		const char*	szShort = (lBase == 2) ? "11" : "2";
		const char*	szUnder = (lBase == 2) ? "0.0011" : (lBase == 3) ? "0.0002" : "0.0003";
		const char*	szValues[] = { szShort, s_long.c_str(), szUnder };
		const char*	szNames[] = { szShort, "three widths long", szUnder };

		for ( int iValue = 0; iValue < 3; iValue++ )
		{
			std::string	s_answers[ 2 ];
			long		l_overrun = LONG_MIN;

			for ( int iPath = 0; iPath < 2; iPath++ )
			{
				lzi_set_log_newton_from( (iPath == 0) ? 200 : kLogNewtonOutOfReach );
				glWorkingPrecision = l_places;

				CLZNumber	x, answer;

				x.SetBase( lBase, false );
				x.Assign( szValues[ iValue ] );
				answer.SetBase( lBase, false );

				lzi_log_newton_overrun();
				lz_log_e( &x, &answer );

				if ( iPath == 0 )
					l_overrun = lzi_log_newton_overrun();

				char*	szAnswer = answer.GetStringRep();

				s_answers[ iPath ] = szAnswer;
				::free( szAnswer );
			}

			char	szLabel[ 80 ];
			char	szDetail[ 60 ];

			::snprintf( szLabel, sizeof( szLabel ), "ln( %s ), base %ld at %ld", szNames[ iValue ], lBase, l_places );
			::snprintf( szDetail, sizeof( szDetail ), "%ld places over, %s the series", l_overrun,
						(s_answers[ 0 ] == s_answers[ 1 ]) ? "agrees with" : "DIFFERS from" );
			CheckTrue( szLabel, (l_overrun > LONG_MIN) && (l_overrun <= 0) && (s_answers[ 0 ] == s_answers[ 1 ]), szDetail );
		}
	}

	lzi_set_log_newton_from( l_newtonWas );
	lzi_set_log_agm_floor( l_agmWas );
	glWorkingPrecision = l_saved;
}

static void TestAngleNewtonAtItsWidths()
{
	Banner( "The arcsine by Newton, each stage at its width" );

	/*
		**The seed at its own width, each stage but the last at its own** --
		until 2026-10-17 the seed was made at the full width and converted
		from base ten, writing the full width's worth of decimal places in
		the operand's base, about 3.3 times the width in base two; the step,
		a small quotient, was added unrounded; and s and c went into every
		stage at their full length. Read by lzi_angle_newton_overrun, Newton
		forced: a short operand, which only the seed and the iterate can
		overrun, and one three widths long, which s carries. Each answer is
		held to the series' as well.
	*/

	long	l_saved = glWorkingPrecision;
	long	l_floorWas = lzi_set_arcsine_newton_floor( 0 );
	const char*	sz_digits = "0123456789ABCDEF";

	for ( long lBase : { 10L, 3L, 7L, 2L } )
	{
		long	l_places = CLZNumber::EquivalentPlaces( 300, 10, lBase );

		// Three widths of digits after "0.", from a fixed sequence...
		std::string		s_long = "0.";
		std::uint64_t	ul_seed = 11;

		for ( long i = 0; i < (3 * l_places); i++ )
		{
			ul_seed = (ul_seed * 6364136223846793005ULL) + 1442695040888963407ULL;
			s_long += sz_digits[ (ul_seed >> 33) % static_cast<std::uint64_t>( lBase ) ];
		}

		const char*	szValues[] = { "0.1", s_long.c_str() };
		const char*	szNames[] = { "0.1", "three widths long" };

		for ( int iValue = 0; iValue < 2; iValue++ )
		{
			std::string	s_answers[ 2 ];
			long		l_overrun = LONG_MIN;
			long		l_newtons = 0;

			for ( int iPath = 0; iPath < 2; iPath++ )
			{
				lzi_set_arcsine_newton_floor( (iPath == 0) ? 0 : LONG_MAX );
				glWorkingPrecision = l_places;

				CLZNumber	x, answer;

				x.SetBase( lBase, false );
				x.Assign( szValues[ iValue ] );
				answer.SetBase( lBase, false );

				long	l_newtonsWere = lzi_arcsine_newtons();

				lzi_angle_newton_overrun();
				lz_arcsine( &x, &answer );

				if ( iPath == 0 )
				{
					l_overrun = lzi_angle_newton_overrun();
					l_newtons = lzi_arcsine_newtons() - l_newtonsWere;
				}

				char*	szAnswer = answer.GetStringRep();

				s_answers[ iPath ] = szAnswer;
				::free( szAnswer );
			}

			char	szLabel[ 80 ];
			char	szDetail[ 60 ];

			::snprintf( szLabel, sizeof( szLabel ), "arcsin( %s ), base %ld at %ld", szNames[ iValue ], lBase, l_places );
			::snprintf( szDetail, sizeof( szDetail ), "%ld places over, %s the series", l_overrun,
						(s_answers[ 0 ] == s_answers[ 1 ]) ? "agrees with" : "DIFFERS from" );
			CheckTrue( szLabel, (l_newtons == 1) && (l_overrun > LONG_MIN) && (l_overrun <= 0) && (s_answers[ 0 ] == s_answers[ 1 ]), szDetail );
		}
	}

	lzi_set_arcsine_newton_floor( l_floorWas );
	glWorkingPrecision = l_saved;
}

static void TestAgmOperandsNarrowed()
{
	Banner( "The AGM's operands, narrowed to its first cycle" );

	/*
		**The first cycle multiplied the operands at their own length** --
		until 2026-10-17; every cycle after it takes means its operations
		rounded. Read by lzi_agm_first_overrun: an operand three widths
		long beside one the width, in four bases. Each answer is held to the
		AGM taken thirty places wider and rounded back, and its cycles to
		well under the cap. Then three short pairs whose means end a unit
		apart: rounding the means every cycle instead stops an odd base's
		settling there, and each ran to the cap.
	*/

	long	l_saved = glWorkingPrecision;
	const char*	sz_digits = "0123456789ABCDEF";

	for ( long lBase : { 10L, 3L, 7L, 2L } )
	{
		long	l_places = CLZNumber::EquivalentPlaces( 300, 10, lBase );
		std::string		s_long = (lBase == 2) ? "10." : "2.";
		std::string		s_width = "1.";
		std::uint64_t	ul_seed = 13;

		for ( long i = 0; i < (3 * l_places); i++ )
		{
			ul_seed = (ul_seed * 6364136223846793005ULL) + 1442695040888963407ULL;
			s_long += sz_digits[ (ul_seed >> 33) % static_cast<std::uint64_t>( lBase ) ];

			if ( i < l_places )
				s_width += sz_digits[ (ul_seed >> 40) % static_cast<std::uint64_t>( lBase ) ];
		}

		std::string	s_answers[ 2 ];
		long		l_overrun = LONG_MIN;
		long		l_cycles = 0;

		for ( int iWide = 0; iWide < 2; iWide++ )
		{
			glWorkingPrecision = l_places + ((iWide == 0) ? 0 : 30);

			CLZNumber	x, y, answer;

			x.SetBase( lBase, false );
			x.Assign( s_long.c_str() );
			y.SetBase( lBase, false );
			y.Assign( s_width.c_str() );
			answer.SetBase( lBase, false );

			long	l_cyclesWere = lzi_agm_cycles();

			lzi_agm_first_overrun();
			lz_agm( &x, &y, &answer );

			if ( iWide == 0 )
			{
				l_overrun = lzi_agm_first_overrun();
				l_cycles = lzi_agm_cycles() - l_cyclesWere;
			}
			else
				answer.Round( kDefaultRoundingMode, l_places );

			char*	szAnswer = answer.GetStringRep();

			s_answers[ iWide ] = szAnswer;
			::free( szAnswer );
		}

		glWorkingPrecision = l_saved;

		char	szLabel[ 80 ];
		char	szDetail[ 60 ];

		::snprintf( szLabel, sizeof( szLabel ), "agm( three widths, the width ), base %ld at %ld", lBase, l_places );
		::snprintf( szDetail, sizeof( szDetail ), "%ld places over, %ld cycles, %s thirty wider", l_overrun, l_cycles,
					(s_answers[ 0 ] == s_answers[ 1 ]) ? "agrees with" : "DIFFERS from" );
		CheckTrue( szLabel, (l_overrun > LONG_MIN) && (l_overrun <= 0) && (l_cycles > 0) && (l_cycles < 64) &&
					(s_answers[ 0 ] == s_answers[ 1 ]), szDetail );
	}

	// Short pairs whose means end a unit apart -- each ran to the cap of
	//	256 with the means rounded every cycle...
	struct { long lBase; long lPlaces; const char* szX; const char* szY; } pairs[] =
	{
		{ 3, 76, "21", "1.1" }, { 5, 52, "2", "1" }, { 7, 43, "1", "0.2" },
	};

	for ( auto& thePair : pairs )
	{
		glWorkingPrecision = thePair.lPlaces;

		CLZNumber	x, y, answer;

		x.SetBase( thePair.lBase, false );
		x.Assign( thePair.szX );
		y.SetBase( thePair.lBase, false );
		y.Assign( thePair.szY );
		answer.SetBase( thePair.lBase, false );

		long	l_cyclesWere = lzi_agm_cycles();

		lz_agm( &x, &y, &answer );

		long	l_cycles = lzi_agm_cycles() - l_cyclesWere;

		glWorkingPrecision = l_saved;

		char	szLabel[ 80 ];
		char	szDetail[ 60 ];

		::snprintf( szLabel, sizeof( szLabel ), "agm( %s, %s ), base %ld at %ld, settles", thePair.szX, thePair.szY, thePair.lBase, thePair.lPlaces );
		::snprintf( szDetail, sizeof( szDetail ), "%ld cycles", l_cycles );
		CheckTrue( szLabel, (l_cycles > 0) && (l_cycles < 64), szDetail );
	}

	glWorkingPrecision = l_saved;
}

int main()
{
	lzu_setup();

	printf( "Addaptic accuracy and precision tests\n" );
	printf( "=====================================\n" );

	TestExactArithmetic();
	TestDivisionExact();
	TestDivisionAccuracy();
	TestDivisionSigns();
	TestDivisionRoundTrip();
	TestDivisionEdges();
	TestPackedDivision();
	TestPackedMultiply();
	TestRounding();
	TestPrecisionControl();
	TestMagnitudeComparison();
	TestTranscendentals();
	TestArgumentMagnitude();
	TestReciprocalTrig();
	TestIdentities();
	TestNumberParts();
	TestNumberBases();
	TestQuietFailures();
	TestFlagPropagation();
	TestRaising();
	TestPowerEdges();
	TestRendering();
	TestDisplay();
	TestGamma();
	TestMeans();
	TestAgm();
	TestAnswerBases();
	TestWorkInOperandBase();
	TestAnswerPrecision();
	TestAnswerWidth();
	TestRefusalsAnswer();
	TestRoundingModes();
	TestTargetAtTheUnits();
	TestShortMeansExact();
	TestRootUnderTheWidth();
	TestRootsByDoubling();
	TestRootsDecided();
	TestHalvedInOddBases();
	TestReciprocalRoundedOnce();
	TestCarrySpelled();
	TestQuotientGuard();
	TestRetakeInDoubt();
	TestWholeQuotients();
	TestWideE();
	TestGuardsInDigits();
	TestArcTangentFamily();
	TestLogInItsBase();
	TestArcSine();
	TestNearOne();
	TestWidePi();
	TestQuotientsNearPoles();
	TestWideArcSine();
	TestWideGamma();
	TestLogFarFromOne();
	TestLogBaseWide();
	TestGammaBelowNought();
	TestAnswerLimit();
	TestCombinationsWhole();
	TestPermutationsWhole();
	TestComplexExact();
	TestArcTangent2();
	TestAngleUnits();
	TestSineCosinePairs();
	TestPolar();
	TestComplexLogExpRoot();
	TestComplexTrig();
	TestComplexPowers();
	TestComplexInverses();
	TestComplexLogInverses();
	TestUnitHalves();
	TestWholePowers();
	TestRootOfASmallOrder();
	TestExpBySplitting();
	TestLogByAgm();
	TestRootCore();
	TestNearTiesDecided();
	TestRootsOfWholeOrder();
	TestSignedComparisons();
	TestNothingLeftBehind();
	TestCubeRoundedOnce();
	TestMixedPairsExact();
	TestEndingBases();
	TestLogarithmTies();
	TestTinyArgumentTies();
	TestCarryKeepsDoubt();
	TestAmeanRoundedOnce();
	TestNearOneDecided();
	TestAbsOfNought();
	TestSmallSinhKeepsItsDigits();
	TestGammaInPlace();
	TestQuantizeIntoItsQuantum();
	TestAssignAWideDouble();
	TestPrefixAtTheFront();
	TestFlagsCarriedThroughDivision();
	TestRoundLeftOfThePoint();
	TestRootsAndHyperbolicsTakeTheShortRoad();
	TestExactnessByPlacesFirst();
	TestLogarithmRoadBySeriesLength();
	TestWideDividendCut();
	TestReciprocalAtPassWidth();
	TestReciprocalBySignificantDigits();
	TestSinCosBySplitting();
	TestSinCosWired();
	TestArcSineByNewton();
	TestBaseConversionByProducts();
	TestCubeRootPassesAtEachWidth();
	TestCubeRootByReciprocal();
	TestAgmRootRoundedToItsWidth();
	TestLogNewtonAtItsWidths();
	TestAngleNewtonAtItsWidths();
	TestAgmOperandsNarrowed();

	printf( "\n=====================================\n" );
	printf( "Below the %d-digit target:\n", kAccuracyTarget );

	bool bAnyWeak = false;

	for ( int i = 0; i < g_measurementCount; i++ )
	{
		if ( g_measurements[ i ].iDigits < kAccuracyTarget )
		{
			printf( "  %-34s %3d digits\n", g_measurements[ i ].szLabel, g_measurements[ i ].iDigits );
			bAnyWeak = true;
		}
	}

	if ( !bAnyWeak )
		printf( "  (none)\n" );

	printf( "\npassed %ld, failed %ld, known issues %ld\n", g_passed, g_failed, g_known );

	lzu_takedown();

	return ( g_failed == 0 ) ? 0 : 1;
}
