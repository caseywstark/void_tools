# Compiler
CXX := CC

# Flags
CXXFLAGS := -std=c++11 -Wall -Wextra -Wno-long-long -pedantic -funroll-loops -O3
CPPFLAGS :=
LDFLAGS  :=
###

APP := vector_radial_profiles
OBJS := $(APP).o
TARGET := $(APP).exe

# Targets
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(TARGET)

.PHONY: clean
clean:
	rm -rf $(TARGET) $(OBJS)
