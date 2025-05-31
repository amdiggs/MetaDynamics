CXX = clang++
CXXFLAGS = -g -Wall -std=c++17
TARGET = meta-dyn
SOURCES = main.cpp MyVec.cpp jdftx.cpp MetaDyn.cpp
LIBS = OpenCL
all: $(TARGET)
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) -framework OpenCL

clean:
	rm -f $(TARGET)
