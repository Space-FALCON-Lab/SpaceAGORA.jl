GRAM - The Global Reference Atmosheric Models suite.    {#mainpage}
===================================================

This document contains programmer information for integrating the GRAM library with a simulation.  It also contains detailed information on the GRAM architecture.

Quick Links
===========

- [About the GRAM Interfaces](@ref AboutGRAM)
- [Building the GRAM Interfaces and Programs](Documentation/Markdown/build.md)
- [Example Programs](Documentation/Markdown/examples.md)
- [The C++ interface](Documentation/Markdown/GRAM.md)
- [The C interface](Documentation/Markdown/GRAM_C.md)
- [The FORTRAN interface](Documentation/Markdown/GRAM_F.md)
- [The Python interface](Documentation/Markdown/GRAMpy.md)
- [The GRAM Architecture](Documentation/Markdown/architecture.md)
- [Creating a New GRAM Model](Documentation/Markdown/model.md)
- [Standard Units in GRAM](Documentation/Markdown/units.md)
- [Name Mapping (Legacy to New)](Documentation/Markdown/namemapping.md)
- [NAIF SPICE Files Needed to Run GRAMs / SPICE Initialization](Documentation/Markdown/spicefiles.md)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

About the GRAM Interfaces {#AboutGRAM}
=========================

Interfaces have been developed for C++, C and FORTRAN. The GRAM atmosphere models are designed around a common C++ framework. The C and FORTRAN interfaces are [wrappers](https://en.wikipedia.org/wiki/Wrapper_library) around the C++ implementation. The call patterns and data structures of all three interfaces are very similar.  And since all of the GRAM models share a common framework, they operate in a unified form.  The unified interface is detailed below with the details about model specific parameters following.

First, gain an understanding of the common interface pattern for the language of interest.

- [The C++ interface.](Documentation/Markdown/GRAM.md)
- [The C interface.](Documentation/Markdown/GRAM_C.md)
- [The FORTRAN interface.](Documentation/Markdown/GRAM_F.md)

Next, learn about the variations available in each model.

| C++                            | C                                 | FORTRAN                           |
|--------------------------------|-----------------------------------|-----------------------------------|
| [Venus](Venus/Venus.md)        | [Venus_C](Venus/Venus_C.md)       | [Venus_F](Venus/Venus_F.md)       |
| [Earth](Earth/Earth.md)        | [Earth_C](Earth/Earth_C.md)       | [Earth_F](Earth/Earth_F.md)       |
| [Mars](Mars/Mars.md)           | [Mars_C](Mars/Mars_C.md)          | [Mars_F](Mars/Mars_F.md)          |
| [Jupiter](Jupiter/Jupiter.md)  | [Jupiter_C](Jupiter/Jupiter_C.md) | [Jupiter_F](Jupiter/Jupiter_F.md) |
| [Titan](Titan/Titan.md)        | [Titan_C](Titan/Titan_C.md)       | [Titan_F](Titan/Titan_F.md)       |
| [Uranus](Uranus/Uranus.md)     | [Uranus_C](Uranus/Uranus_C.md)    | [Uranus_F](Uranus/Uranus_F.md)    |
| [Neptune](Neptune/Neptune.md)  | [Neptune_C](Neptune/Neptune_C.md) | [Neptune_F](Neptune/Neptune_F.md) |


