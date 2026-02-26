# Using MSVS to build GRAMs

1. Open the solution GRAMs.sln using MSVS 2017.
2. Select Release mode.
3. Rebuild the solution.

For developers running MSVS 2015, the same solution and projects can be used with some modifications.  Follow these steps:
 1. Start MSVS 2015 and open GRAM\MSVS\GRAMs.sln.
 2. In the Solution Explorer use ctrl-click to select all projects excluding the Shared Items projects.
 3. Right-click and select *Properties*.
 4. In the dropdowns select "All Configurations" and "All Platforms".
 5. Under *Configuration Properties : General* set the *Platform Toolset* to "Visual Studio 2015 (v140)".
 6. Set the *Target Platform Version* to "8.1".
 7. Apply the changes and save the solution.

Files and Subfolders: (All subfolders correspond to MSVS projects.)
- GRAMs.sln: The MSVS solution containing all GRAM projects.
- Common: Shared items containing the GRAM framework.
- Common Unit Tests: Shared items containing unit tests for the GRAM framework.
- googletest: Unit testing framework.
- GRAMCommonLib: Builds a library containing the GRAM framework.
- GRAMLib: Builds a library containing the GRAMCommonLib and all planet libraries.
- GRAM: Builds an executable that combines all planet GRAMs.
- GRAM_C: Builds an example of using the generic C interface.
- GRAMExamples: Combines all planet GRAMs into one example program.
- Planet: Shared items containing the models for the specified planet.
- Planet Unit Test: Shared items containing unit tests for the specified planet.
- PlanetLib: Builds a library containing the specified planet's models.
- PlanetGRAM: Builds the GRAM executable for the specified planet (e.g. NeptuneGRAM.exe).
- PlanetExamples: Builds an executable using the common examples.
- PlanetExample_Doc: Builds an executable example used in the Programmer's Manual.
- PlanetExample_C: Builds an executable using an example of the C interface.
- EarthCorrMonte: Builds the Earth CorrMonte executable.
- EarthCorrTraj: Builds the Earth CorrTraj executable.
- EarthCorrMulti: Builds an Earth example of mulitple correlated trajectories.
