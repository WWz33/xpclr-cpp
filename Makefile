# XP-CLR C++/htslib rewrite
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp
CPPFLAGS += -Iinclude
# htslib: prefer pkg-config, then /usr/local, then local biosoft build
HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LIBS := $(shell pkg-config --libs htslib 2>/dev/null)
ifeq ($(HTS_LIBS),)
  HTS_CFLAGS := -I/usr/local/include
  HTS_LIBS := -L/usr/local/lib -lhts -Wl,-rpath,/usr/local/lib
endif
CPPFLAGS += $(HTS_CFLAGS)
LDFLAGS += $(HTS_LIBS) -fopenmp -lm -lz -lpthread -lcurl -lbz2 -llzma -lgsl -lgslcblas

SRC := src/main.cpp src/util.cpp src/pop.cpp src/vcf_io.cpp src/xpclr.cpp
OBJ := $(SRC:.cpp=.o)
BIN := xpclr

.PHONY: all clean test-help

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

src/%.o: src/%.cpp include/xpclr.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJ) $(BIN)

test-help: $(BIN)
	./$(BIN) -h
	./$(BIN) -v
