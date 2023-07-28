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

TMP_SOURCES = $(filter-out %eos_assembled.cpp,$(SOURCES))
TMP_OBJECTS = $(filter-out %eos_assembled.o,$(OBJECTS))

TABLES_SRC := $(wildcard $(TABLES_DIR)*.cpp) 
TABLES_OBJ := $(patsubst $(TABLES_DIR)%.cpp, $(TABLES_DIR)compile/%.o, $(TABLES_SRC))

LEP_SRC := $(SRC_DIR)complete_FG.cpp $(SRC_DIR)FD_functions.cpp $(SRC_DIR)find_eta.cpp $(SRC_DIR)eos_fermions.cpp
LEP_OBJ := $(patsubst $(SRC_DIR)%.cpp, $(COMPILE_DIR)src/%.o, $(LEP_SRC))

TESTS_SRC := $(wildcard $(TESTS_DIR)*.cpp)
TESTS_OBJ := $(patsubst $(TESTS_DIR)%.cpp, $(COMPILE_DIR)%.o, $(TESTS_SRC))

all: $(TMP_OBJECTS) eos_assembled_interp main
	ar rcs $(COMPILE_DIR)libout.a $(OBJECTS)

main: main.o $(TMP_OBJECTS) eos_assembled_otf
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) main.o -o main $(H5_LIB) $(H5_FLAGS)

run: test_performance test_eos
	$(COMPILE_DIR)test_eos_performance; $(COMPILE_DIR)test_eos

test_performance: $(OBJECTS) $(COMPILE_DIR)test_eos_performance.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos_performance.o -o $(COMPILE_DIR)test_eos_performance $(H5_LIB) $(H5_FLAGS)

test_eos: $(OBJECTS) $(COMPILE_DIR)test_eos.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos.o -o $(COMPILE_DIR)test_eos $(H5_LIB) $(H5_FLAGS)

test_assembled: $(OBJECTS) $(COMPILE_DIR)test_eos_assembled.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos_assembled.o -o $(COMPILE_DIR)test_eos_assembled $(H5_LIB) $(H5_FLAGS)

test_accuracy: $(OBJECTS) $(COMPILE_DIR)test_eos_accuracy.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_eos_accuracy.o -o $(COMPILE_DIR)test_eos_accuracy

test_baryons: $(OBJECTS) $(COMPILE_DIR)test_baryon_eos.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(COMPILE_DIR)test_baryon_eos.o -o $(COMPILE_DIR)test_baryon_eos $(H5_LIB) $(H5_FLAGS)

lep_table: $(LEP_OBJ) $(TABLES_DIR)compile/generate_lep_table.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(LEP_OBJ) $(TABLES_DIR)compile/generate_lep_table.o -o $(TABLES_DIR)compile/generate_lep_table

global_table: $(TMP_OBJECTS) eos_assembled_otf $(TABLES_DIR)compile/generate_global_table.o
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/generate_global_table.o -o $(TABLES_DIR)compile/generate_global_table $(H5_LIB) $(H5_FLAGS)

comparison_table: $(TMP_OBJECTS) eos_assembled_otf $(TABLES_DIR)compile/generate_comparison_table.o $(TABLES_DIR)compile/compare_w_mu_ele.o 
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/generate_comparison_table.o -o $(TABLES_DIR)compile/generate_comparison_table $(H5_LIB) $(H5_FLAGS)
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/compare_w_mu_ele.o -o $(TABLES_DIR)compile/compare_w_mu_ele $(H5_LIB) $(H5_FLAGS)

test_hist: $(TMP_OBJECTS) eos_assembled_otf $(TABLES_DIR)compile/test_hist.o 
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) $(OBJECTS) $(TABLES_DIR)compile/test_hist.o -o $(TABLES_DIR)compile/test_hist $(H5_LIB) $(H5_FLAGS)

$(COMPILE_DIR)%.o: $(TESTS_DIR)%.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $@ $<

$(COMPILE_DIR)src/%.o: $(SRC_DIR)%.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) -c -o $@ $< $(CPPFLAGS) $(H5_LIB) $(H5_FLAGS)

eos_assembled_interp: $(SRC_DIR)eos_assembled.cpp
	$(CXX) -Did_test=1 $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) -c -o $(COMPILE_DIR)src/eos_assembled.o $< $(CPPFLAGS) $(H5_LIB) $(H5_FLAGS)

eos_assembled_otf: $(SRC_DIR)eos_assembled.cpp
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) -c -o $(COMPILE_DIR)src/eos_assembled.o $< $(CPPFLAGS) $(H5_LIB) $(H5_FLAGS)

main.o: main.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(INCLUDE) -c -o $@ $< $(CPPFLAGS) $(H5_LIB) $(H5_FLAGS)

$(TESTS_DIR)compile/%.o: $(TESTS_DIR)%.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $@ $<

$(TABLES_DIR)compile/%.o: $(TABLES_DIR)%.cpp
	$(CXX) -Did_test=2 $(RFLAGS) $(INCLUDE) -c -o $@ $<

clean:
	rm -v -f $(COMPILE_DIR)src/*.o
	rm -v -f $(TABLES_DIR)compile/*.o
	rm -v -f $(TABLES_DIR)compile/generate_*
	rm -v -f main.o
	rm -v -f ./main
#find $(COMPILE_DIR) -type f -exec rm -v -f {} +

