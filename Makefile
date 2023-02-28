CXX      := g++
RFLAGS   := -std=c++17 -O3 -g -Wall -lm #-Wpedantic
H5_LIB     := -L/usr/lib/x86_64-linux-gnu/hdf5/serial
H5_INCLUDE := -I/usr/include/hdf5/serial
H5_FLAGS := -lhdf5 -lhdf5_hl
#$(HDF5_INCLUDE)
#LIB      := ar cr
#FLAGS    := 
#LDFLAGS  += $(HDF5_LIBS)

#INCLUDE = -I/usr/local/include
#LIBS    = -L/usr/local/lib

INCLUDE := -I$(CURDIR)/include

COMPILE_DIR = $(CURDIR)/tests/compile/

SRC_DIR = $(CURDIR)/src/

TESTS_DIR := $(CURDIR)/tests/

TABLES_DIR := $(CURDIR)/eos_table/

SOURCES := $(wildcard $(SRC_DIR)*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)%.cpp, $(COMPILE_DIR)src/%.o, $(SOURCES)) 

TABLES_SRC := $(wildcard $(TABLES_DIR)*.cpp) 
TABLES_OBJ := $(patsubst $(TABLES_DIR)%.cpp, $(TABLES_DIR)compile/%.o, $(TABLES_SRC))

TESTS_SRC := $(wildcard $(TESTS_DIR)*.cpp)
TESTS_OBJ := $(patsubst $(TESTS_DIR)%.cpp, $(COMPILE_DIR)%.o, $(TESTS_SRC))

all: test_performance

run: test_performance test_eos
	$(COMPILE_DIR)test_eos_performance; $(COMPILE_DIR)test_eos

test_performance: $(OBJECTS) $(COMPILE_DIR)test_eos_performance.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos_performance.o -o $(COMPILE_DIR)test_eos_performance $(H5_LIB) $(H5_FLAGS)

test_eos: $(OBJECTS) $(COMPILE_DIR)test_eos.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos.o -o $(COMPILE_DIR)test_eos

test_assembled: $(OBJECTS) $(COMPILE_DIR)test_eos_assembled.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos_assembled.o -o $(COMPILE_DIR)test_eos_assembled $(H5_LIB) $(H5_FLAGS)

test_accuracy: $(OBJECTS) $(COMPILE_DIR)test_eos_accuracy.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos_accuracy.o -o $(COMPILE_DIR)test_eos_accuracy

test_baryons: $(OBJECTS) $(COMPILE_DIR)test_baryon_eos.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_baryon_eos.o -o $(COMPILE_DIR)test_baryon_eos $(H5_LIB) $(H5_FLAGS)

generate_table: lep_table
	$(TABLES_DIR)compile/generate_lep_table

eta_table: $(OBJECTS) $(TABLES_DIR)compile/generate_eta_table.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/generate_eta_table.o -o $(TABLES_DIR)compile/generate_eta_table

complete_table: $(OBJECTS) $(TABLES_DIR)compile/generate_complete_table.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/generate_complete_table.o -o $(TABLES_DIR)compile/generate_complete_table

lep_table: $(OBJECTS) $(TABLES_DIR)compile/generate_lep_table.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/generate_lep_table.o -o $(TABLES_DIR)compile/generate_lep_table


$(COMPILE_DIR)%.o: $(TESTS_DIR)%.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $@ $<

$(COMPILE_DIR)src/%.o: $(SRC_DIR)%.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) -c -o $@ $< $(CPPFLAGS) $(H5_LIB) $(H5_FLAGS)

$(TESTS_DIR)compile/%.o: $(TESTS_DIR)%.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $@ $<

$(TABLES_DIR)compile/%.o: $(TABLES_DIR)%.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $@ $<
clean:
	find $(COMPILE_DIR) -type f -exec rm -f {} +
	rm -v -f $(COMPILE_DIR)src/*
	rm -v -f $(TABLES_DIR)compile/*

