# Tests
Using `pytest` to run the tests is recommended for more 
useful and well formatted output, but `unittest` can also 
be used if `pytest` is not installed. 

Notes for developers: 

* `test_earth.py` requires an additional data download to run- the tests WILL fail if these downloads are not present. Follow the instructions in `Earth/data/NCEPdata/README.md` to download these files. 

* The tests use relative filepaths to access Earth data and Mars data. They also look for the `gram.py`, `_gram.so`, and `spice.txt` 
files at the root of the `/tests/` directory. If you alter the directory structure, be sure to edit `EARTH_DATA_PATH` in `test_earth.py`, `MARS_DATA_PATH` in `test_mars.py`, and `SpicePath` in `spice.txt` to point to the new locations of these files. You will also need to ensure that `gram.py`, `_gram.so`, and `spice.txt` are present in the directory you run the test command from, or that your tests import the `gram` module using an appropriate relative import. 

### To run tests using `pytest`: 

```
//Change into this directory (/tests/)

$ pytest
```

### To run tests using `unittest`: 

```
//Change into this directory (/tests/)

$ python -m unittest
```