LIBRARY = GRAMCommon
include ../makefile.defs

$(info LIBRARY is [$(LIBRARY)])
$(info BUILD is [$(BUILD)])

LIBRARY_NAME := lib$(LIBRARY)

# Targets
.PHONY: all lib clean

# Directories
SOURCE_DIRS  = $(COMMON_DIR)/source $(COMMON_DIR)/examples
INCLUDE_DIRS = $(COMMON_DIR)/include $(COMMON_DIR)/googletest $(COMMON_DIR)/googletest/include $(COMMON_DIR)/cspice/include $(COMMON_DIR)/examples
INCLUDE_CMD  = $(addprefix -I, $(INCLUDE_DIRS))
SOURCE_FILES = $(foreach SDIR, $(SOURCE_DIRS), $(wildcard $(SDIR)/*.cpp))

OBJECT_FILES = $(addprefix $(OBJECT_DIR)/,  $(notdir $(SOURCE_FILES:.cpp=.o)))
OBJECT_TESTS = $(addprefix $(OBJECT_DIR)/,  $(notdir $(TEST_FILES:.cpp=.o)))

vpath %.cpp $(SOURCE_DIRS) $(TEST_DIRS)

all: lib

lib: $(LIBRARY_DIR)/$(LIBRARY_NAME).a 

ifdef FC
$(LIBRARY_DIR)/$(LIBRARY_NAME).a: $(OBJECT_FILES) $(OBJECT_DIR)/Gram_F.o
	mkdir -p $(LIBRARY_DIR)
	ar rc $@ $(OBJECT_FILES)
	ranlib $@
else
$(LIBRARY_DIR)/$(LIBRARY_NAME).a: $(OBJECT_FILES) 
	mkdir -p $(LIBRARY_DIR)
	ar rc $@ $(OBJECT_FILES)
	ranlib $@
endif

$(OBJECT_FILES) $(OBJECT_TESTS): $(OBJECT_DIR)/%.o: %.cpp
	@ -mkdir -p $(OBJECT_DIR)
	$(CXX) -o $@ $(CXX_FLAGS) $(OPT_FLAGS) $(INCLUDE_CMD) -c $<

$(OBJECT_DIR)/Gram_F.o:
	@ -mkdir -p $(OBJECT_DIR)
	@ -mkdir -p $(MODULE_DIR)
	@ $(FC) -o $(OBJECT_DIR)/Gram_F.o -c -J$(MODULE_DIR) $(FC_FLAGS) $(OPT_FLAGS) $(COMMON_DIR)/source/Gram_F.f90


clean:
	rm -Rf $(OBJECT_DIR) $(LIBRARY_DIR)/$(LIBRARY_NAME).a $(MODULE_DIR)/gramstructs.mod
