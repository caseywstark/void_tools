# Compiler
CXX := CC

# Flags
CXXFLAGS := -std=c++11 -Wall -Wextra -Wno-long-long -pedantic -funroll-loops -O3
CPPFLAGS :=
LDFLAGS  :=
###

APP := void_profiles
OBJS := $(APP).o
TARGET := $(APP).ex

# Targets
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(TARGET)

.PHONY: clean
clean:
	rm -rf $(TARGET) $(OBJS)
