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


SRC_DIR   = src/
OBJ_DIR   = $(SRC_DIR)obj/
OUT_DIR   = output/
TEST_DIR  = tests/
LIB_DIR   = lib/

SRC_FILES := $(shell find $(SRC_DIR) -name "*.cpp")
SRC_FILES := $(filter-out %eos_assembled.cpp, $(SRC_FILES))

OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))

SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

OBJECTS := $(patsubst $(SRC_DIR)%.cpp, $(SRC_DIR)%.o, $(SRC_FILES)) 
 
all: $(OBJ_FILES) main

main: $(OBJ_FILES) eos_otf main.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o main.o -o main $(H5_LIB) $(H5_FLAGS)

main.o: main.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) -c -o $@ $< $(H5_LIB) $(H5_FLAGS)

$(OBJ_DIR)%.o: %.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) -c -o $@ $< $(H5_LIB) $(H5_FLAGS)

$(OUT_DIR)obj/%.o: $(OUT_DIR)%.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) -c -o $@ $< $(H5_LIB) $(H5_FLAGS)

$(TEST_DIR)obj/%.o: $(TEST_DIR)src/%.cpp
	$(CXX) $(RFLAGS) $(H5_INCLUDE) -c -o $@ $< $(H5_LIB) $(H5_FLAGS)

eos_interp: $(SRC_DIR)eos_species/eos_assembled.cpp
	$(CXX) -Did_test=1 $(RFLAGS) $(H5_INCLUDE) -c -o $(OBJ_DIR)eos_assembled_int.o $(SRC_DIR)eos_species/eos_assembled.cpp $(H5_LIB) $(H5_FLAGS)

eos_otf: $(SRC_DIR)eos_species/eos_assembled.cpp
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) -c -o $(OBJ_DIR)eos_assembled_otf.o $(SRC_DIR)eos_species/eos_assembled.cpp $(H5_LIB) $(H5_FLAGS)

output_table: $(OBJ_FILES) eos_otf $(OUT_DIR)obj/generate_output_table.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o $(OUT_DIR)obj/generate_output_table.o -o $(OUT_DIR)bin/generate_output_table $(H5_LIB) $(H5_FLAGS)
	$(OUT_DIR)bin/generate_output_table

lep_table: $(OBJ_FILES) eos_otf $(OUT_DIR)obj/generate_lep_table.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o $(OUT_DIR)obj/generate_lep_table.o -o $(OUT_DIR)bin/generate_lep_table $(H5_LIB) $(H5_FLAGS)
	$(OUT_DIR)bin/generate_lep_table

comparison_with_mu: $(OBJ_FILES) eos_otf $(TEST_DIR)obj/comparison_with_mu.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o $(TEST_DIR)obj/comparison_with_mu.o -o $(TEST_DIR)bin/test_comparison_with_mu $(H5_LIB) $(H5_FLAGS)
	$(TEST_DIR)bin/test_comparison_with_mu

comparison_without_mu: $(OBJ_FILES) eos_otf $(TEST_DIR)obj/comparison_without_mu.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o $(TEST_DIR)obj/comparison_without_mu.o -o $(TEST_DIR)bin/test_comparison_without_mu $(H5_LIB) $(H5_FLAGS)
	$(TEST_DIR)bin/test_comparison_without_mu

hist_comp_with_mu: $(OBJ_FILES) eos_otf $(TEST_DIR)obj/hist_comp_with_mu.o
	$(CXX) $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o $(TEST_DIR)obj/hist_comp_with_mu.o -o $(TEST_DIR)bin/test_hist_comp_with_mu $(H5_LIB) $(H5_FLAGS)
	$(TEST_DIR)bin/test_hist_comp_with_mu

lib: $(OBJ_FILES) eos_otf eos_interp
	ar rcs $(LIB_DIR)libmueos_otf.a $(OBJ_FILES) $(OBJ_DIR)eos_assembled_otf.o
	ar rcs $(LIB_DIR)libmueos_int.a $(OBJ_FILES) $(OBJ_DIR)eos_assembled_int.o

clean:
	rm -f $(OBJ_FILES)
	rm -f $(OBJ_DIR)*.o
	rm -f $(OUT_DIR)obj/*.o
	rm -f $(OUT_DIR)bin/generate*
	rm -f $(TEST_DIR)obj/*.o
	rm -f $(TEST_DIR)bin/test*
	rm -f main
	rm -f main.o
