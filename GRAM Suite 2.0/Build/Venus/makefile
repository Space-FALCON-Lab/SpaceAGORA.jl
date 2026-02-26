PLANET = Venus
include ../makefile.defs

$(info PLANET is [$(PLANET)])
$(info BUILD is [$(BUILD)])

LIBRARY_NAME := lib$(PLANET)
PLANET_DIR := $(BASE_DIR)/$(PLANET)

# Targets
.PHONY: all lib gram examples docexamples cexamples fexamples clean

# Directories
SOURCE_DIRS  = $(PLANET_DIR)/source
TEST_DIRS    = $(PLANET_DIR)/unittest $(COMMON_DIR)/unittest
INCLUDE_DIRS = $(PLANET_DIR)/include $(COMMON_DIR)/include $(COMMON_DIR)/googletest $(COMMON_DIR)/googletest/include $(COMMON_DIR)/cspice/include $(COMMON_DIR)/examples
INCLUDE_CMD  = $(addprefix -I, $(INCLUDE_DIRS))
SOURCE_FILES = $(foreach SDIR, $(SOURCE_DIRS), $(wildcard $(SDIR)/*.cpp))
TEST_FILES   = $(foreach SDIR, $(TEST_DIRS), $(wildcard $(SDIR)/*.cpp))

OBJECT_FILES = $(addprefix $(OBJECT_DIR)/,  $(notdir $(SOURCE_FILES:.cpp=.o)))
OBJECT_TESTS = $(addprefix $(TEST_OBJECT_DIR)/,  $(notdir $(TEST_FILES:.cpp=.o)))

vpath %.cpp $(SOURCE_DIRS) $(TEST_DIRS)

ifdef FC
all: lib gram examples docexamples cexamples fexamples
else
all: lib gram examples docexamples cexamples
endif

gram: $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(OBJECT_TESTS)
	mkdir -p $(BINARY_DIR)
	$(CXX) -o $(OBJECT_DIR)/$(PLANET)GRAM.o $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $(PLANET_DIR)/examples/$(PLANET)GRAM.cpp
	$(LNK) -o $(BINARY_DIR)/$(PLANET)GRAM$(BINARY_EXTENSION)  $(OBJECT_TESTS) $(OBJECT_DIR)/$(PLANET)GRAM.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libGRAMCommon.a $(LIBRARY_DIR)/libgtest.a $(SPICE_LIB) $(LNK_FLAGS) 
	@ rm $(OBJECT_DIR)/$(PLANET)GRAM.o

examples: $(LIBRARY_DIR)/$(LIBRARY_NAME).a
	mkdir -p $(BINARY_DIR)
	$(CXX) -o $(OBJECT_DIR)/ExampleRunner.o $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $(PLANET_DIR)/examples/ExampleRunner.cpp
	$(LNK) -o $(BINARY_DIR)/$(PLANET)Examples$(BINARY_EXTENSION) $(OBJECT_DIR)/ExampleRunner.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libGRAMCommon.a $(SPICE_LIB) $(LNK_FLAGS)
	@ rm $(OBJECT_DIR)/ExampleRunner.o

docexamples: $(LIBRARY_DIR)/$(LIBRARY_NAME).a 
	mkdir -p $(BINARY_DIR)
	$(CXX) -o $(OBJECT_DIR)/$(PLANET).o $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $(PLANET_DIR)/examples/$(PLANET).cpp
	$(LNK) -o $(BINARY_DIR)/$(PLANET)Example_Cpp$(BINARY_EXTENSION) $(OBJECT_DIR)/$(PLANET).o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libGRAMCommon.a $(SPICE_LIB) $(LNK_FLAGS)
	@ rm $(OBJECT_DIR)/$(PLANET).o

cexamples: $(LIBRARY_DIR)/$(LIBRARY_NAME).a 
	mkdir -p $(BINARY_DIR)
	$(CC) -o $(OBJECT_DIR)/$(PLANET)_C.o $(CC_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $(PLANET_DIR)/examples/$(PLANET)_C.c
	$(LNK) -o $(BINARY_DIR)/$(PLANET)Example_C$(BINARY_EXTENSION) $(OBJECT_DIR)/$(PLANET)_C.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libGRAMCommon.a $(SPICE_LIB) $(LNK_FLAGS)
	@ rm $(OBJECT_DIR)/$(PLANET)_C.o

fexamples: $(LIBRARY_DIR)/$(LIBRARY_NAME).a
	mkdir -p $(BINARY_DIR)
	$(FC) -o $(OBJECT_DIR)/$(PLANET)_F.o $(FC_FLAGS) $(OPT_FLAGS) -c -J$(MODULE_DIR) $(PLANET_DIR)/examples/$(PLANET)_F.f90
	$(FLNK) -o $(BINARY_DIR)/$(PLANET)Example_F$(BINARY_EXTENSION) $(OBJECT_DIR)/$(PLANET)_F.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libGRAMCommon.a $(SPICE_LIB) $(FLNK_FLAGS)
	@ rm $(OBJECT_DIR)/$(PLANET)_F.o
	
lib: $(LIBRARY_DIR)/$(LIBRARY_NAME).a

ifdef FC
$(LIBRARY_DIR)/$(LIBRARY_NAME).a: $(OBJECT_FILES) $(OBJECT_DIR)/$(PLANET)Atmosphere_F.o
	mkdir -p $(LIBRARY_DIR)
	ar rc $@ $(OBJECT_FILES)
	ranlib $@
else
$(LIBRARY_DIR)/$(LIBRARY_NAME).a: $(OBJECT_FILES)
	mkdir -p $(LIBRARY_DIR)
	ar rc $@ $(OBJECT_FILES)
	ranlib $@
endif

$(OBJECT_FILES): $(OBJECT_DIR)/%.o: %.cpp
	@ -mkdir -p $(OBJECT_DIR)
	$(CXX) -o $@ $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $<

$(OBJECT_TESTS): $(TEST_OBJECT_DIR)/%.o: %.cpp
	@ -mkdir -p $(TEST_OBJECT_DIR)
	$(CXX) -o $@ $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $<

$(OBJECT_DIR)/$(PLANET)Atmosphere_F.o:
	@ -mkdir -p $(OBJECT_DIR)
	@ -mkdir -p $(MODULE_DIR)
	@ $(FC) -o $(OBJECT_DIR)/$(PLANET)Atmosphere_F.o -c -J$(MODULE_DIR) $(FC_FLAGS) $(OPT_FLAGS) $(PLANET_DIR)/source/$(PLANET)Atmosphere_F.f90

clean:
	rm -Rf $(OBJECT_DIR) $(TEST_OBJECT_DIR) $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(BINARY_DIR)/$(PLANET)* $(MODULE_DIR)/*.mod
