# PYTHON_INCLUDE - Location of the Python include folder.
# PYTHON_LIB - Location of the Python shared object library. (when required by the compiler)
# SHARED_EXTENSION - Shared object extension expected by your Python instance. (typically so or pyd)
# SPICE_LIB - Location of the SPICE toolkit library.
# If using your own build of the SPICE library or are not using gcc to compile,
# change this path to point to the correct version of the SPICE library. 
# See common/cspice/lib for other precompiled options. 

# Linux
PYTHON_INCLUDE := /usr/include/python3.12
PYTHON_LIB := -lpython3.12
SHARED_EXTENSION := so
SPICE_LIB := ../../common/cspice/lib/cspice_gcc85.a

# MSYS2 - Windows
# PYTHON_INCLUDE := /mingw64/include/python3.11
# PYTHON_LIB := -lpython3.11.dll
# SHARED_EXTENSION := pyd  
# SPICE_LIB := ../../common/cspice/lib/cspice_mingw64.a

##########################################################################################

.PHONY: gram wrapper regenerate all clean

all: gram

# Assumes that gram_wrap.cxx is already present in this directory
gram:
	g++ -g -Wall -O2 -std=c++11 -fpic -c gram_wrap.cxx \
	-I$(PYTHON_INCLUDE) \
	-I../../common/include/ \
	-I../../Earth/include/ \
	-I../../Jupiter/include/ \
	-I../../Mars/include/ \
	-I../../Neptune/include/ \
	-I../../Titan/include/ \
	-I../../Uranus/include/ \
	-I../../Venus/include/ 

	g++ -g -Wall -O2 -std=c++11 -shared gram_wrap.o \
	../../Build/lib/libGRAM.a \
	../../Build/lib/libgtest.a \
	$(SPICE_LIB) $(PYTHON_LIB) \
	-o _gram.$(SHARED_EXTENSION)

	rm *.o


# Regenerates gram_wrap.cxx and gram.py.
# Use only if changes were made in the C++ source code or the gram.i interface file
wrapper:
	swig -c++ -python -doxygen \
	-I../../common/include/ \
	-I../../Earth/include/ \
	-I../../Jupiter/include/ \
	-I../../Mars/include/ \
	-I../../Neptune/include/ \
	-I../../Titan/include/ \
	-I../../Uranus/include/ \
	-I../../Venus/include/ \
	gram.i


# Regenerates gram_wrap.cxx and gram.py and builds the shared object. 
# Use only if changes were made in the C++ source code or the gram.i interface file
regenerate: wrapper gram


clean:
	rm _gram.$(SHARED_EXTENSION) gram.py ../tests/_gram.$(SHARED_EXTENSION) ../tests/gram.py
	rm -r __pycache__/ ../tests/__pycache__/ ../tests/.pytest_cache/
	