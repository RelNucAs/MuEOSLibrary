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


SRC_DIR = src/
OBJ_DIR = $(SRC_DIR)compile/

SRC_FILES := $(shell find $(SRC_DIR) -name "*.cpp")

OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.cpp=.o)))

SRC_DIRS := $(dir $(SRC_FILES))
VPATH := $(SRC_DIRS)

OBJECTS := $(patsubst $(SRC_DIR)%.cpp, $(SRC_DIR)%.o, $(SRC_FILES)) 
 
all: $(OBJ_FILES) main
#@echo $(OBJ_FILES)

main: $(OBJ_FILES) main.cpp
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) $(OBJ_FILES) main.o -o main $(H5_LIB) $(H5_FLAGS)

$(OBJ_DIR)%.o: %.cpp
	$(CXX) -Did_test=2 $(RFLAGS) $(H5_INCLUDE) -c -o $@ $< $(H5_LIB) $(H5_FLAGS)

clean:
	rm -f $(OBJ_FILES)
	