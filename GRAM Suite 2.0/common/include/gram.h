//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

// This is the default location of the NAIF SPICE kernels.
#define DEFAULT_SPICE_PATH "/SPICE"

// The default is that all reals are type double.
typedef double greal;

// Most reals can be switched to type float by changing the greal typedef.
// Some data is remains as doubles in this case, but most computations convert to float.
// Also, the NAIF SPICE libs will need to be compiled float computations.
//typedef float greal;

// When using type float, uncomment the lines below to suppress MSC warnings.
//#ifdef _MSC_VER
//#pragma warning( disable : 4244 )
//#pragma warning( disable : 4305 )
//#endif

#ifdef __cplusplus

#include <cassert>
#include <cstddef>

//! The namespace for all GRAM models.
namespace GRAM {

extern bool CONSOLE_OUTPUT;
inline void enableConsoleOutput(bool flag = true) { CONSOLE_OUTPUT = flag; }

constexpr greal PI = (greal)(3.141592653589793);
constexpr greal TWO_PI = (greal)(6.283185307179586);          // 2 * PI
constexpr greal HALF_PI = (greal)(1.5707963267948965);        // PI / 2
constexpr greal TO_DEGREES = (greal)(57.29577951308232);      // 180 / PI
constexpr greal TO_RADIANS = (greal)(0.0174532925199432950);  // PI / 180
constexpr greal DEGREES_PER_HOUR = (greal)(15.0);             // 360 / 24
constexpr greal HOURS_PER_DEGREE = (greal)(1.0 / 15.0);       // 24 / 360
constexpr greal STANDARD_ATMOSPHERE = (greal)(101325.0);      // Pascals  (NIST: 2014 CODATA)
constexpr greal J2000_EPOCH = (greal)(2451545.0);
constexpr greal JAN_0_2000 = (greal)(2451543.5);
constexpr greal KILOMETERS_PER_AU = (greal)(149597870.7);

constexpr greal AVOGADRO = (greal)(6.02214076e23);            // mol^-1 (NIST: 2018 CODATA recommended values)
constexpr greal BOLTZMANN = (greal)(1.380649e-23);            // J K^-1  (NIST: 2018 CODATA recommended values)
constexpr greal UNIVERSAL_GAS = (greal)(8.314462618);         // J K^-1 mol^-1  (NIST: 2018 CODATA recommended values)
constexpr greal AU_IN_LIGHT_SECONDS = (greal)(149597870700.0 / 299792458.0);
constexpr greal GRAMS_PER_AMU = (greal)(1.6605390666e-24);    // g amu^-1 (NIST: 2018 CODATA recommended values)

constexpr bool EAST_POSITIVE = true;
constexpr bool WEST_POSITIVE = false;


//! \brief Time scales are time standards used for planetary motion calculations.
//!
//! Precise definitions are available in the Spice documentation.
//! Short and long form enumerations result in equivalent computations.
enum GRAM_TIME_SCALE { 
  UTC,                         //!< Coordinated universal time
  COORDINATED_UNIVERSAL_TIME,  //!< Coordinated universal time
  TDB,                         //!< Barycentric dynamical time
  BARYCENTRIC_DYNAMICAL_TIME,  //!< Barycentric dynamical time
  TDT,                         //!< Terrestrial dynamical time
  TERRESTRIAL_DYNAMICAL_TIME   //!< Terrestrial dynamical time
};

//! \brief Time frames are event locations used for planetary motion calculations.
//!
//! Planet event time (PET) is similar to space craft event time (SCET).
//! Short and long form enumerations result in equivalent computations.
enum GRAM_TIME_FRAME { 
  ERT,                //!< Earth receive time
  EARTH_RECEIVE_TIME, //!< Earth receive time
  PET,                //!< Planet event time
  PLANET_EVENT_TIME   //!< Planet event time
};

//! \brief All constituent gases available for use in a model.
//!
//! The gases within a model are identified with this enumeration.  Note how monatomic and diatomic
//! gases are identified (e.g. HYDROGEN, DIHYDROGEN).
enum GasType {
  // NOTE: GasType is used in the ConstituentGas class as an index.
  // If this enum is modified, then also modify the ConstituentGas class.
  ARGON,           //!< Argon (Ar)
  HELIUM,          //!< Helium (He)
  HYDROGEN,        //!< Atomic hydrogen (H)
  DIHYDROGEN,      //!< Diatomic hydrogen (H2)
  NITROGEN,        //!< Atomic nitrogen (N)
  DINITROGEN,      //!< Diatomic nitrogen (N2)
  OXYGEN,          //!< Atomic oxygen (O)
  DIOXYGEN,        //!< Diatomic oxygen (O2)
  METHANE,         //!< Methane (CH4)
  CARBON_MONOXIDE, //!< Carbon monoxide (CO)
  CARBON_DIOXIDE,  //!< Carbon dioxide (CO2)
  OZONE,           //!< Triatomic ozygen (O3)
  NITROUS_OXIDE,   //!< Dinitrogen monoxide (N2O)
  WATER,           //!< Water vapor (H2O)
  GAS_TYPE_SIZE    //!< Number of elements in this enum
};

//! \brief Bodies supported by GRAMs.
enum GRAM_BODY { 
  NO_BODY, //!< No body chosen
  VENUS,   //!< Venus
  EARTH,   //!< Earth
  MARS,    //!< Mars
  JUPITER, //!< Jupiter
  SATURN,  //!< Saturn
  URANUS,  //!< Uranus
  NEPTUNE, //!< Neptune
  TITAN    //!< Titan
};

enum UpdateStatus { 
  NO_UPDATES,    //!< No perturbation updates performed
  STEP_UPDATED,  //!< Step size has been updated
  PERTS_UPDATED, //!< Perturbations and step size have been updated
  INITIAL_STATE  //!< Initial sim state
};

enum PerturbationAction { 
  UPDATE_PERTS,       //!< Update perturbations in the next update cycle
  DO_NOT_UPDATE_PERTS //!< Do not update perturbations in the next update cycle
};

enum DensityPrintScale { 
  STANDARD,    //!< Standard density output scale (kg/m^3)
  LOG_10,      //!< Base 10 log of density is output (kg/m^3)
  PERCENT_DEV, //!< Deviations are scaled by 100 (kg/m^3)
  KM           //!< Density is scaled by 1.0e9 (kg/km^3)
};

extern greal wrapDegrees(greal angle);
extern greal wrapDegrees180(greal angle);
extern greal wrapRadians(greal angle);
extern greal wrapRadiansPi(greal angle);

inline greal toDegrees(greal radians) { return radians * TO_DEGREES; }
inline greal toRadians(greal degrees) { return degrees * TO_RADIANS; }
inline greal clamp(greal value, greal low, greal high) { assert(low <= high); return value < low ? low : (value > high ? high : value); }
inline int clampSize(int value, int size) { return value < 0 ? 0 : (value >= size ? size - 1 : value); }
inline size_t clampSize(size_t value, size_t size) { return value >= size ? size - 1 : value; }
inline greal square(greal value) { return value * value; }

constexpr greal operator"" _pa(long double pressure) { return (greal)pressure; }
constexpr greal operator"" _deg(long double angle) { return (greal)angle; }
constexpr greal operator"" _rad ( long double angle ) { return (greal)angle; }
constexpr greal operator"" _au ( long double distance ) { return (greal)distance; }
constexpr greal operator"" _km ( long double distance ) { return (greal)distance; }
constexpr greal operator"" _m ( long double distance ) { return (greal)distance; }
constexpr greal operator"" _day ( long double time ) { return (greal)time; }
constexpr greal operator"" _hr ( long double time ) { return (greal)time; }
constexpr greal operator"" _min ( long double time ) { return (greal)time; }
constexpr greal operator"" _sec ( long double time ) { return (greal)time; }

constexpr greal operator"" _UNK ( long double xxx ) { return (greal)xxx; }

// We can get rid of this macro with C++17 as it lets one
// declare pairs for maps
#define FOR_ALL_GASES_PRESENT(X)  \
     for (auto& gasIter : gases) {           /* Iterate over the set of all gases. */ \
       ConstituentGas& gas = gasIter.second; /* The second parameter of the iterator is the gas. */ \
       if (gas.isPresent) {                  /* Test if the gas is activated. */ \
         GasType gasType = gasIter.first;    /* The first parameter of the iterator is the gas type. */ \
         X                                   /* Loop content goes here. */ \
       }  \
     }

#define ELSE_NOT_PRESENT } else {
} // namespace

#endif