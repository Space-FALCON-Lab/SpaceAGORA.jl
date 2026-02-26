# For most users: building the Python interface 
If you just want to use the GRAM Suite from Python without making any changes to the underlying C++ code, this section is for you. No additional dependencies are required to build the Python interface module. 

1. Build the full GRAM suite, following the instructions provided in the GRAM Programmer's Manual (available in the `Documentation/` folder). 
2. Open the `makefile` in this directory and change the PYTHON_INCLUDE variable to the path to the headers of the Python installation you will be using. The easiest way to find this path is to start an interactive Python session. 
```
$ python 
>>> from sysconfig import get_paths
>>> get_paths()['include']
```
Copy the output filepath and paste it into the `makefile`. Also set PYTHON_LIB, SHARED_EXTENSION, and SPICE_LIB in the makefile.

3. From this directory (`/python/GRAMpy`):
```
$ make 
```
4. Following these steps will compile the `_gram.so` shared library needed by the Python interface module. Now that you have built the interface, see the README in the root `python/` folder for usage instructions. 

# For developers: re-generating the Python interface

GRAMpy uses a third party tool called SWIG to programatically generate Python bindings for the 
C++ files indicated in `gram.i`. If you make changes to the underlying C++ 
code or wish to alter the interface, you will need to re-generate the `gram.py` and `_gram.so` files that define the Python interface module.  

1. Download [SWIG](https://swig.org). 
2. Build the full GRAM suite, following the instructions provided in the GRAM Programmer's Manual (available in the `Documentation/` folder). 
3. Open the `makefile` in this directory and change the PYTHON_INCLUDE variable to the path to the headers of the Python installation you will be using. The easiest way to find this path is to start an interactive Python session. 
```
$ python 
>>> from sysconfig import get_paths
>>> get_paths()['include']
```
Copy the output filepath and paste it into the `makefile`.

4. From this directory (`/python/GRAMpy`):
```
$ make regenerate
```

5. Following these steps will regenerate the following files: `gram_wrap.cxx`, `gram.py`, and `_gram.so`. Now that you have built the interface, see the README in the root `python/` folder for usage instructions. 

## Files and Subfolders:
- GRAMpy: Make system for building the Python interface using SWIG.
- JupiterExample.py: A basic example of using GRAMpy.
- MarsExample.py: An example of using GRAMpy.
- GRAMExample.ipynb: A Jupyter notebook tutorial for using GRAMpy.
- GRAMExample.html: A non-interactive form of the Jupyter notebook.
- spice.txt: A GRAM path setup file required by the examples.