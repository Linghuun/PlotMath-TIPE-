# Compiler
CXX = g++

# Compiler flags
SAFE_CXXFLAGS = -std=c++11 -Wall -Wextra -fsanitize=address -g
CXXFLAGS = -std=c++11 

# Include directories
INCLUDES = -I../glfw-3.4/include -I../glew-2.1.0/include

# Linker flags
LDFLAGS = -L../glew-2.1.0/lib -Wl,-rpath,../glew-2.1.0/lib

# Libraries to link
LIBS = -lGLEW -lglfw -lGL -lX11 -lpthread -lXrandr -lXi -ldl

# Source files (excluding main)
SRCS = point.cpp quadtree.cpp camera.cpp polynome.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = plotmath

# Default target (non-secure build)
all: CXXFLAGS := $(CXXFLAGS)
all: $(EXEC)

# Secure build target
safe: CXXFLAGS := $(SAFE_CXXFLAGS)
safe: $(EXEC)

# Link the executable (include file.cpp directly)
$(EXEC): $(OBJS) file.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -O0 -o $@ file.cpp $(OBJS) $(LDFLAGS) $(LIBS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all safe clean
