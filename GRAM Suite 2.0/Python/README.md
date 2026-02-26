# Python interface for the GRAM Suite

- To build the Python interface: follow the directions in `python/GRAMpy/README.md`
- To run automated tests: follow the directions in `python/tests/README.md`

The rest of this document covers using the GRAM Suite from Python.  
A special thanks to Zeb Becker, a most excellent summer intern, for developing this interface.

## Setup
Python installation via package managers (pip/conda/etc.) is not currently supported. To use the `gram` module, users should import it using a relative import in their programs. 

The simplest way to do this is to copy the entire `GRAMpy` directory into the folder where you are writing your Python program, then import it like this in your Python program: 
```
from GRAMpy import gram
``` 
Advanced users who wish to set their import up differently certainly may. Just note that `_gram.so` must be in the same directory as `gram.py`. 

## Usage
The provided Python interface simply wraps GRAM's core C++ code, allowing the user to call C++ functions from a program written in Python. All exposed classes and methods should have the same behavior as they would in C++, meaning that Python programmers can reference any of the existing C++ documentation in the GRAM Programmer's Manual in addition to Python specific documentation. 

## Running the example programs 
### Basic example
`JupiterExample.py` provides a very simple example of how to set up an atmosphere model and collect information from a user-defined location and time. To run this example, from the `python/` directory:
```
python JupiterExample.py
```

### Advanced examples
More involved examples are provided in a Jupyter notebook, `gramExample.ipynb`. 

#### Preliminaries
This tutorial assumes that you have followed the instructions in the README to build the GRAM Suite and the Python interface on your computer.  As a general Python best practice, it is recommended to use a virtual environment to help manage dependencies. 

The commands in 
```
code blocks
```
are meant to be executed from your command line. To follow along with this tutorial, open a command line prompt and navigate to the `Python/` directory (where this file is located). 

#### Opening the Jupyter Notebook

0. (Optional) Create and activate a virtual environment 
    ```
    python -m venv env
    source env/bin/activate
    ```
1. Install dependencies needed to run the Jupyter notebook
    ```
    pip install jupyter ipython numpy matplotlib
    ```
2. Open the Jupyter notebook GUI. This will serve a local web interface for your use. To open the notebook, copy one of the links in the terminal output into a browser. 
    ```
    jupyter notebook
    ```
3. Using the web interface, open the `gramExample.ipynb` notebook

4. Run the code cells in order. If you make changes to a code cell, simply run it again and you will see your changes reflected in the output. 

## Files and Subfolders:
- GRAMpy: Make system for building the Python interface using SWIG.
- JupiterExample.py: A basic example of using GRAMpy.
- MarsExample.py: An example of using GRAMpy.
- GRAMExample.ipynb: A Jupyter notebook tutorial for using GRAMpy.
- GRAMExample.html: A non-interactive form of the Jupyter notebook.
- spice.txt: A GRAM path setup file required by the examples.