The NeptuneGRAM F Interface
===========================

This documents the FORTRAN interface for NeptuneGRAM.  Since the GRAM framework is common to all models, only the Neptune specific information is documented on this page.  See the [FORTRAN interface](Documentation/Markdown/GRAM_F.md) page for documentation on the common framework.

A simple example of using the Neptune atmosphere can be seen [here](@ref Neptune/examples/Neptune_F.f90).

- - - - - - - 

The Min Max Factor
------------------

The function \ref setMinMaxFactor_N sets a parameter that determines where within the envelope of minimum-to-maximum a given profile falls.  The minMaxFactor can be any real number between -1 (minimum) and +1 (maximum). Setting the minMaxFactor to 0 gives the average (recommended and default) Neptune profile.  

~~~~~~~~~~
call setMinMaxFactor_N(neptune, -0.2d0, 1)
~~~~~~~~~~

NeptuneGRAM has a built-in routine for generating random perturbations about the selected mean profile, but the amplitude of these variations (especially at high altitudes) is less than the range of maximum to minimum in the Neptune model envelope.  Users wishing to vary the mean profile (in a Monte-Carlo sense) between Neptune minimum and maximum profiles, can do so by randomly selecting values of the minMaxFactor (between -1 and +1).  

Users wishing to more explicitly account for latitudinal, seasonal, and time-of-day effects on location within the minimum-maximum envelope may be guided by the following table in selecting values for the minMaxFactor.


|  Effect of:   |   minMax negative |    minMax near 0   | minMax positive |
|---------------|:-----------------:|:------------------:|:---------------:|
|Latitude/Season| Winter-polar lats |Near-equatorial lats<br>Equinox, all lats | Summer-polar lats|
|Time of Day    |        Night      |     Near Twilight  |          Day    |

The last argument of \ref setMinMaxFactor_N flags the automatic adjustment of the minMaxFactor for seasonal, latitude, and time-of-day effects.  The adjustments are based on the user supplied factor. Set this flag to 1 to automatically adjust the minMaxFactor or to 0 to hold the value constant.

- - - - - - - 

The Nitrogen (N2) Mole Fraction
------------------

By default, the dinitrogen mole fraction will be 0.  The method \ref setDinitrogenMoleFraction_N allows the mole fraction to be set between 0% and 0.6%.

~~~~~~~~~~
call setDinitrogenMoleFraction_N(neptune, 0.002d0) ! input as a decimal percentage
~~~~~~~~~~






