LIBRARY = GRAM
include ../makefile.defs

$(info LIBRARY is [$(LIBRARY)])
$(info BUILD is [$(BUILD)])

LIBRARY_NAME := lib$(LIBRARY)
PLANET_DIR := $(BASE_DIR)/$(LIBRARY)

# Targets
.PHONY: all lib examples shared clean

# Directories
SOURCE_DIRS  = $(PLANET_DIR)/Examples
TEST_DIRS    = $(COMMON_DIR)/unittest 
INCLUDE_DIRS = $(PLANET_DIR)/include $(COMMON_DIR)/include $(COMMON_DIR)/googletest $(COMMON_DIR)/googletest/include $(COMMON_DIR)/cspice/include $(COMMON_DIR)/examples 
OBJECT_LIST  = ../GRAMCommon/obj/* 
DEFINES =

ifneq ($(wildcard ../Venus),)
TEST_DIRS += $(BASE_DIR)/Venus/unittest
INCLUDE_DIRS += $(BASE_DIR)/Venus/include
OBJECT_LIST += ../Venus/obj/*.o
DEFINES += -DVENUS_GRAM
endif

ifneq ($(wildcard ../Earth),)
TEST_DIRS += $(BASE_DIR)/Earth/unittest
INCLUDE_DIRS += $(BASE_DIR)/Earth/include
OBJECT_LIST += ../Earth/obj/*.o
DEFINES += -DEARTH_GRAM
endif

ifneq ($(wildcard ../Mars),)
TEST_DIRS += $(BASE_DIR)/Mars/unittest
INCLUDE_DIRS += $(BASE_DIR)/Mars/include
OBJECT_LIST += ../Mars/obj/*.o
DEFINES += -DMARS_GRAM
endif

ifneq ($(wildcard ../Jupiter),)
TEST_DIRS += $(BASE_DIR)/Jupiter/unittest
INCLUDE_DIRS += $(BASE_DIR)/Jupiter/include
OBJECT_LIST += ../Jupiter/obj/*.o
DEFINES += -DJUPITER_GRAM
endif

ifneq ($(wildcard ../Saturn),)
TEST_DIRS += $(BASE_DIR)/Saturn/unittest
INCLUDE_DIRS += $(BASE_DIR)/Saturn/include
OBJECT_LIST += ../Saturn/obj/*.o
DEFINES += -DSATURN_GRAM
endif

ifneq ($(wildcard ../Titan),)
TEST_DIRS += $(BASE_DIR)/Titan/unittest
INCLUDE_DIRS += $(BASE_DIR)/Titan/include
OBJECT_LIST += ../Titan/obj/*.o
DEFINES += -DTITAN_GRAM
endif

ifneq ($(wildcard ../Uranus),)
TEST_DIRS += $(BASE_DIR)/Uranus/unittest
INCLUDE_DIRS += $(BASE_DIR)/Uranus/include
OBJECT_LIST += ../Uranus/obj/*.o
DEFINES += -DURANUS_GRAM
endif

ifneq ($(wildcard ../Neptune),)
TEST_DIRS += $(BASE_DIR)/Neptune/unittest
INCLUDE_DIRS += $(BASE_DIR)/Neptune/include
OBJECT_LIST += ../Neptune/obj/*.o
DEFINES += -DNEPTUNE_GRAM
endif

INCLUDE_CMD  = $(addprefix -I, $(INCLUDE_DIRS))
SOURCE_FILES = $(foreach SDIR, $(SOURCE_DIRS), $(wildcard $(SDIR)/*.cpp))
TEST_FILES   = $(foreach SDIR, $(TEST_DIRS), $(wildcard $(SDIR)/*.cpp))

OBJECT_TESTS = $(addprefix $(TEST_OBJECT_DIR)/,  $(notdir $(TEST_FILES:.cpp=.o)))

SHARED_EXTENSION ?= so
SHARED_LNK_FLAGS ?= -shared

vpath %.cpp $(SOURCE_DIRS) $(TEST_DIRS)

all: lib examples allgrams gram_c

examples: $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(OBJECT_TESTS)
	mkdir -p $(OBJECT_DIR)
	mkdir -p $(BINARY_DIR)
	$(CXX) -o $(OBJECT_DIR)/ExampleRunner.o $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) $(DEFINES) -c $(PLANET_DIR)/examples/ExampleRunner.cpp
	$(LNK) -o $(BINARY_DIR)/GRAMExamples$(BINARY_EXTENSION) $(OBJECT_TESTS) $(OBJECT_DIR)/ExampleRunner.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libgtest.a $(SPICE_LIB) $(LNK_FLAGS)
	@ rm $(OBJECT_DIR)/ExampleRunner.o
	
allgrams: $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(OBJECT_TESTS)
	mkdir -p $(OBJECT_DIR)
	mkdir -p $(BINARY_DIR)
	$(CXX) -o $(OBJECT_DIR)/GRAM.o $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) $(DEFINES) -c $(PLANET_DIR)/examples/GRAM.cpp
	$(LNK) -o $(BINARY_DIR)/GRAM$(BINARY_EXTENSION) $(OBJECT_TESTS) $(OBJECT_DIR)/GRAM.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libgtest.a $(SPICE_LIB) $(LNK_FLAGS)
	@ rm $(OBJECT_DIR)/GRAM.o
	
gram_c: $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(OBJECT_TESTS)
	mkdir -p $(OBJECT_DIR)
	mkdir -p $(BINARY_DIR)
	$(CXX) -o $(OBJECT_DIR)/GRAM_C.o $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) $(DEFINES) -c $(PLANET_DIR)/examples/GRAM_C.c
	$(LNK) -o $(BINARY_DIR)/GRAM_C$(BINARY_EXTENSION) $(OBJECT_TESTS) $(OBJECT_DIR)/GRAM_C.o $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libgtest.a $(SPICE_LIB) $(LNK_FLAGS)
	@ rm $(OBJECT_DIR)/GRAM_C.o

shared: all
	mkdir -p $(LIBRARY_DIR)
	$(CXX) $(SHARED_LNK_FLAGS) -o $(LIBRARY_DIR)/libGRAM.$(SHARED_EXTENSION) $(OBJECT_LIST) $(LIBRARY_DIR)/libgtest.a $(SPICE_LIB) $(LNK_FLAGS)
	
lib: $(LIBRARY_DIR)/$(LIBRARY_NAME).a

$(LIBRARY_DIR)/$(LIBRARY_NAME).a: 
	mkdir -p $(LIBRARY_DIR)
	ar rc $@ $(OBJECT_LIST)
	ranlib $@

$(OBJECT_FILES): $(OBJECT_DIR)/%.o: %.cpp
	@ -mkdir -p $(OBJECT_DIR)
	$(CXX) -o $@ $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $<

$(OBJECT_TESTS): $(TEST_OBJECT_DIR)/%.o: %.cpp
	@ -mkdir -p $(TEST_OBJECT_DIR)
	$(CXX) -o $@ $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $<

clean:
	rm -Rf $(OBJECT_DIR) $(TEST_OBJECT_DIR) $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(LIBRARY_DIR)/libGRAM.so $(LIBRARY_DIR)/libGRAM.dylib $(BINARY_DIR)/GRAM* 
