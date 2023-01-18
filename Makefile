CXX      := g++
RFLAGS   := -std=c++17 -O3 -g -Wall -lm #-Wpedantic
#CPPFLAGS += $(HDF5_INCLUDE)
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

all: test_eos test_performance

test_performance: test_performance.o interp.o find_eta.o eos_fermions.o complete_FG.o FD_functions.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(COMPILE_DIR)interp.o $(COMPILE_DIR)find_eta.o $(COMPILE_DIR)eos_fermions.o $(COMPILE_DIR)complete_FG.o $(COMPILE_DIR)FD_functions.o $(COMPILE_DIR)test_eos_performance.o -o $(COMPILE_DIR)test_eos_performance

test_eos: test_eos.o interp.o find_eta.o eos_fermions.o complete_FG.o FD_functions.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(COMPILE_DIR)interp.o $(COMPILE_DIR)find_eta.o $(COMPILE_DIR)eos_fermions.o $(COMPILE_DIR)complete_FG.o $(COMPILE_DIR)FD_functions.o $(COMPILE_DIR)test_eos.o -o $(COMPILE_DIR)test_eos

test_accuracy: test_accuracy.o interp.o find_eta.o eos_fermions.o complete_FG.o FD_functions.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(COMPILE_DIR)interp.o $(COMPILE_DIR)find_eta.o $(COMPILE_DIR)eos_fermions.o $(COMPILE_DIR)complete_FG.o $(COMPILE_DIR)FD_functions.o $(COMPILE_DIR)test_eos_accuracy.o -o $(COMPILE_DIR)test_eos_accuracy

test_accuracy.o: $(TESTS_DIR)test_eos_accuracy.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)test_eos_accuracy.o $(TESTS_DIR)test_eos_accuracy.cpp

test_eos.o: $(TESTS_DIR)test_eos.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)test_eos.o $(TESTS_DIR)test_eos.cpp

test_performance.o: $(TESTS_DIR)test_eos_performance.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)test_eos_performance.o $(TESTS_DIR)test_eos_performance.cpp

interp.o: $(SRC_DIR)interp.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)interp.o $(SRC_DIR)interp.cpp

find_eta.o: $(SRC_DIR)find_eta.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)find_eta.o $(SRC_DIR)find_eta.cpp

eos_fermions.o: $(SRC_DIR)eos_fermions.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)eos_fermions.o $(SRC_DIR)eos_fermions.cpp

complete_FG.o: $(SRC_DIR)complete_FG.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)complete_FG.o $(SRC_DIR)complete_FG.cpp

FD_functions.o: $(SRC_DIR)FD_functions.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(COMPILE_DIR)FD_functions.o $(SRC_DIR)FD_functions.cpp


eta_table: eta_table.o interp.o find_eta.o eos_fermions.o complete_FG.o FD_functions.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(COMPILE_DIR)interp.o $(COMPILE_DIR)find_eta.o $(COMPILE_DIR)eos_fermions.o $(COMPILE_DIR)complete_FG.o $(COMPILE_DIR)FD_functions.o $(TABLES_DIR)compile/generate_eta_table.o -o $(TABLES_DIR)compile/generate_eta_table

complete_table: complete_table.o interp.o find_eta.o eos_fermions.o complete_FG.o FD_functions.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(COMPILE_DIR)interp.o $(COMPILE_DIR)find_eta.o $(COMPILE_DIR)eos_fermions.o $(COMPILE_DIR)complete_FG.o $(COMPILE_DIR)FD_functions.o $(TABLES_DIR)compile/generate_complete_table.o -o $(TABLES_DIR)compile/generate_complete_table

lep_table: lep_table.o interp.o find_eta.o eos_fermions.o complete_FG.o FD_functions.o
	$(CXX) $(RFLAGS) $(INCLUDE) $(COMPILE_DIR)interp.o $(COMPILE_DIR)find_eta.o $(COMPILE_DIR)eos_fermions.o $(COMPILE_DIR)complete_FG.o $(COMPILE_DIR)FD_functions.o $(TABLES_DIR)compile/generate_lep_table.o -o $(TABLES_DIR)compile/generate_lep_table

eta_table.o: $(TABLES_DIR)generate_eta_table.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(TABLES_DIR)compile/generate_eta_table.o $(TABLES_DIR)generate_eta_table.cpp

complete_table.o: $(TABLES_DIR)generate_complete_table.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(TABLES_DIR)compile/generate_complete_table.o $(TABLES_DIR)generate_complete_table.cpp

lep_table.o: $(TABLES_DIR)generate_lep_table.cpp
	$(CXX) $(RFLAGS) $(INCLUDE) -c -o $(TABLES_DIR)compile/generate_lep_table.o $(TABLES_DIR)generate_lep_table.cpp
clean:
	rm -v -f $(COMPILE_DIR)*
	rm -v -f $(TABLES_DIR)compile/*

