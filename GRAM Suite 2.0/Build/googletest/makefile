include ../makefile.defs

LIBRARY_NAME := libgtest

GTEST_DIR = $(COMMON_DIR)/googletest

.PHONY: all lib clean

all: lib

lib: $(LIBRARY_DIR)/$(LIBRARY_NAME).a 

$(LIBRARY_DIR)/$(LIBRARY_NAME).a: $(OBJECT_DIR)/gtest-all.o
	mkdir -p $(LIBRARY_DIR)
	ar rc $@ $(OBJECT_DIR)/gtest-all.o
	ranlib $@

$(OBJECT_DIR)/gtest-all.o: 
	@ -mkdir -p $(OBJECT_DIR)
	$(CXX) -o $(OBJECT_DIR)/gtest-all.o $(CXX_FLAGS) $(OPT_FLAGS) -c $(GTEST_DIR)/src/gtest-all.cc -I $(GTEST_DIR)/include -I $(GTEST_DIR)

clean:
	rm -Rf $(OBJECT_DIR) $(LIBRARY_DIR)/$(LIBRARY_NAME).a 