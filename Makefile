# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall

# Include directories
INCLUDES = -I../glfw-3.4/include -I../glew-2.1.0/include

# Linker flags
LDFLAGS = -L../glew-2.1.0/lib -Wl,-rpath,../glew-2.1.0/lib

# Libraries to link
LIBS = -lGLEW -lglfw -lGL -lX11 -lpthread -lXrandr -lXi -ldl

# Source files (list all .cpp files here)
SRCS = point.cpp quadtree.cpp camera.cpp polynome.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = my_project

# Default target
all: $(EXEC)

# Link the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) -g -Wall -Wextra -fsanitize=address $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create the project 
my_project: file.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -Wall -Wextra -fsanitize=address -o $(EXEC) file.cpp $(OBJS) $(LDFLAGS) $(LIBS)

# Clean up build files
clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean


 
