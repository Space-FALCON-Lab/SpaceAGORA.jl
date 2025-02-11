The TitanGRAM C Interface
===========================

This documents the C interface for TitanGRAM.  Since the GRAM framework is common to all models, only the Titan specific information is documented on this page.  See the [C interface](Documentation/Markdown/GRAM_C.md) page for documentation on the common framework.

A simple example of using the Titan atmosphere can be seen [here](@ref Titan/examples/Titan_C.c).

- - - - - - - 

The Min/Max or General Circulation Model
----------------------------------------

There are two Titan atmosphere models within TitanGRAM, the Min/Max model (Yelle97) and the General Circulation model (GCM95).  The default model is the Min/Max model.  But the model may be chosen with setModelType_T() with an argument of 1 for Yelle97 or 2 for GCM95.  It is ill-advised to switch models in the middle of a trajectory.

~~~~~~~~~~
setModelType_T(titan, 2); // GCM95

// or

setModelType_T(titan, 1); // Yelle97
~~~~~~~~~~

- - - - - - - 

The Min Max Factor
------------------

The function \ref setMinMaxFactor_T sets a parameter that determines where within the envelope of minimum-to-maximum a given profile falls.  The minMaxFactor can be any real number between -1 (minimum) and +1 (maximum). Setting the minMaxFactor to 0 gives the average (recommended and default) Titan profile.  

~~~~~~~~~~
setMinMaxFactor_T(titan, -0.2, 1);
~~~~~~~~~~

TitanGRAM has a built-in routine for generating random perturbations about the selected mean profile, but the amplitude of these variations (especially at high altitudes) is less than the range of maximum to minimum in the Titan model envelope.  Users wishing to vary the mean profile (in a Monte-Carlo sense) between Titan minimum and maximum profiles, can do so by randomly selecting values of the minMaxFactor (between -1 and +1).  

Users wishing to more explicitly account for latitudinal, seasonal, and time-of-day effects on location within the minimum-maximum envelope may be guided by the following table in selecting values for the minMaxFactor.


|  Effect of:   |   minMax negative |            minMax near 0                 | minMax positive  |
|---------------|:-----------------:|:----------------------------------------:|:----------------:|
|Latitude/Season| Winter-polar lats |Near-equatorial lats<br>Equinox, all lats | Summer-polar lats|
|Time of Day    |        Night      |             Near Twilight                |          Day     |

The last argument of \ref setMinMaxFactor_T flags the automatic adjustment of the minMaxFactor for seasonal, latitude, and time-of-day effects.  The adjustments are based on the user supplied factor. Set this flag to 1 to automatically adjust the minMaxFactor or to 0 to hold the value constant.

- - - - - - - 

The Methane Mole Fraction
-------------------------

The modelled value of the methane mole fraction may be overridden by the user with the setMethaneMoleFraction_T() function.  The default value of MethaneMoleFraction = 0.0 uses the modelled values.  Values  of MethaneMoleFraction between 1.0 and 5.0 (percent) will override the modelled value of the methane mole fraction and rebalance the other mole fractions while preserving the mean molecular weight.  The will, of course, also affect the number densities and mass fractions. 

~~~~~~~~~~
setMethaneMoleFraction_T(titan, 2.1);
~~~~~~~~~~

- - - - - - - 


