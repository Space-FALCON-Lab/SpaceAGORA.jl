The TitanGRAM C++ Interface
=============================

This documents the \ref GRAM::TitanAtmosphere interface class for TitanGRAM.  Since the GRAM framework is common to all models, only the Titan specific information is documented on this page.  See the [C++ interface](Documentation/Markdown/GRAM.md) page for documentation on the common framework.

A simple example of using the Titan atmosphere class can be seen [here](@ref Titan/examples/Titan.cpp).

- - - - - - - 

The Min/Max or General Circulation Model
----------------------------------------

There are two Titan atmosphere models within TitanGRAM, the Min/Max model (Yelle97) and the General Circulation model (GCM95).  The default model is the Min/Max model.  But the model may be chosen with setModelType() with an argument of Yelle97 or GCM95.  It is ill-advised to switch models in the middle of a trajectory.

~~~~~~~~~~
titan.setModelType(GCM95);

// or

titan.setModelType(Yelle97);
~~~~~~~~~~

- - - - - - - 

The Min Max Factor
------------------

When using the Min/Max model (Yelle97), the method \ref GRAM::TitanAtmosphere::setMinMaxFactor sets a parameter that determines where within the envelope of minimum-to-maximum a given profile falls.  The minMaxFactor can be any real number between -1 (minimum) and +1 (maximum). Setting the minMaxFactor to 0 gives the average (recommended and default) Titan profile.  

~~~~~~~~~~
titan.setMinMaxFactor(-0.2, true);
~~~~~~~~~~

TitanGRAM has a built-in routine for generating random perturbations about the selected mean profile, but the amplitude of these variations (especially at high altitudes) is less than the range of maximum to minimum in the Titan model envelope.  Users wishing to vary the mean profile (in a Monte-Carlo sense) between Titan minimum and maximum profiles, can do so by randomly selecting values of the minMaxFactor (between -1 and +1).  

Users wishing to more explicitly account for latitudinal, seasonal, and time-of-day effects on location within the minimum-maximum envelope may be guided by the following table in selecting values for the minMaxFactor.


|  Effect of:   |   minMax negative |            minMax near 0                 | minMax positive  |
|---------------|:-----------------:|:----------------------------------------:|:----------------:|
|Latitude/Season| Winter-polar lats |Near-equatorial lats<br>Equinox, all lats | Summer-polar lats|
|Time of Day    |        Night      |             Near Twilight                |          Day     |

The second argument of \ref GRAM::TitanAtmosphere::setMinMaxFactor flags the automatic adjustment of the minMaxFactor for seasonal, latitude, and time-of-day effects.  The adjustments are based on the user supplied factor. Set this flag to `true` to automatically adjust the minMaxFactor or to `false` to hold the value constant.

- - - - - - - 

The Methane Mole Fraction
-------------------------

The modelled value of the methane mole fraction may be overridden by the user with the setMethaneMoleFraction() method.  The default value of MethaneMoleFraction = 0.0 uses the modelled values.  Values  of MethaneMoleFraction between 1.0 and 5.0 (percent) will override the modelled value of the methane mole fraction and rebalance the other mole fractions while preserving the mean molecular weight.  The will, of course, also affect the number densities and mass fractions. 

~~~~~~~~~~
titan.setMethaneMoleFraction(2.1);
~~~~~~~~~~

- - - - - - - 
