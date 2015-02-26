# compiler settings
CXX := CC
CPPFLAGS :=
CXXFLAGS := -std=c++11 -Wfatal-errors -Wall -Wextra -pedantic -O3 -funroll-loops

SRC = so_finder.cc
TARGET = so_finder.ex

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(SRC) -o $(TARGET)

clean:
	-rm -rf $(TARGET)
