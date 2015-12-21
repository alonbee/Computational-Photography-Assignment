# ----------------------------------------------------------------------------
# CS 89/189 Computational Aspects of Digital Photography C++ basecode.
#
# Adapted from MIT's 6.815/6.865 basecode, written and designed by:
#     Frédo Durand
#     Katherine L. Bouman
#     Gaurav Chaurasia
#     Adrian Vasile Dalca
#     Neal Wadhwa
#
# With additions & modifications by
#     Wojciech Jarosz
# for Dartmouth's CS 89/189.
#
# ----------------------------------------------------------------------------

# C/C++ compilation involves two stages:
# - compiling: convert each source (*.cpp) file indepdendently into binary object (.o) files
# - linking: assemble all the object files into a single executable
#
# Since any C++ project contains multiple source files, the above
# process involved firing multiple commands. Makefiles help us
# write all the compilation stages in a single file and run them
# as one big script.
#
# See http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/
# and many many other tutorials for makefiles


# some variables
BUILD_DIR  :=_build
EXECUTABLE := a2
OUTPUT := Output
INCLUDE_DIRS := ext

# the C++ compiler/linker to be used. define here so that we can change
# it easily if needed
CXX := g++ -Wall -I$(INCLUDE_DIRS)

# ------------------------------------------------------------------------------

# 'make' or 'make all' runs the default target 'all' which requires that
# the EXECUTABLE must be available; if the EXECUTABLE is already available
# then nothing happens. see the rules for EXECUTABLE

all: $(EXECUTABLE)

# ------------------------------------------------------------------------------

# 'make run' runs the target 'run' which calls for running the executable
# if the executable is not available, then the rules for creating it are run
# this is a superset of the target 'all' because it creates the executable
# and also runs it

run: $(EXECUTABLE)
	./$(EXECUTABLE)

# ------------------------------------------------------------------------------

# 'make clean' runs the target 'clean' which in turn removes the
# intermediate .o files and the executable

clean:
	rm -rf $(BUILD_DIR) $(EXECUTABLE) $(OUTPUT)

# ------------------------------------------------------------------------------

# rule for creating the executable: this "links" the .o files using the g++ linker.
# If .o files are not available, then the rules for creating .o files are run.

$(EXECUTABLE): $(BUILD_DIR)/a2_main.o $(BUILD_DIR)/a2.o $(BUILD_DIR)/floatimage.o $(BUILD_DIR)/lodepng.o
	$(CXX) $(BUILD_DIR)/a2_main.o $(BUILD_DIR)/a2.o $(BUILD_DIR)/floatimage.o $(BUILD_DIR)/lodepng.o -o $(EXECUTABLE)
	mkdir -p $(OUTPUT)

# ------------------------------------------------------------------------------

# rules for creating the .o files:  compile each of the .cpp files and create a
# corresponding .o file. Each .o depends upon the corresponding .cpp and other .h files

$(BUILD_DIR)/a2.o: a2.cpp a2.h floatimage.h exceptions.h ext/lodepng.h
	mkdir -p $(BUILD_DIR)
	$(CXX) -c a2.cpp -o $(BUILD_DIR)/a2.o

$(BUILD_DIR)/a2_main.o: a2_main.cpp a2.h floatimage.h exceptions.h ext/lodepng.h
	mkdir -p $(BUILD_DIR)
	$(CXX) -c a2_main.cpp -o $(BUILD_DIR)/a2_main.o

$(BUILD_DIR)/floatimage.o: floatimage.cpp floatimage.h array3D.h exceptions.h ext/lodepng.h
	mkdir -p $(BUILD_DIR)
	$(CXX) -c floatimage.cpp -o $(BUILD_DIR)/floatimage.o

$(BUILD_DIR)/lodepng.o: ext/lodepng.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) -c ext/lodepng.cpp -o $(BUILD_DIR)/lodepng.o

