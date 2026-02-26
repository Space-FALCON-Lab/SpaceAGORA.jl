How to Run a Code Coverage Analysis on GRAMs
============================================

This document explains how to measure which lines of code have been covered via unit and regression testing.

Requirements
============

Code coverage analysis is run under MSVS using the OpenCppCoverage Plugin.  This plugin is available for 
download via the Tools | Extensions and Updates menu option (using MSVS 2017).

Setup and Analysis
==================

+ Set the BodyTests project as the startup project (where Body is replaced with the name of your planetary body).
+ Under the Tools menu select the OpenCppCoverage Settings option.
+ Choose the Basic tab and make sure the BodyTests project is checked.
+ Scroll down and inspect the name of the Program to run.
+ The Working Directory should be C:\fullpath\GRAMS\Development\GRAM\Body\tests (no relative paths allowed).
+ Under the Miscellaneous tab, browse to select the config.emu file located in GRAMS\Development\GRAM\Body\tests.
+ To export a report, select the Import/Export tab and fill out the Export Type and Path.
+ Select the Run Code Coverage button.

The analysis will color the code to display the coverage.  It will also open up a browser tab for traversing the files.
The browser tab has a checkmark to toggle the coloring.  Uncheck this option before closing the tab.

If errors occur when running the analysis, they are almost always related to paths.  Remember, no relative paths.  Running the tool from the command line may help, too.
