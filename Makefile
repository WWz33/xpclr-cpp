# XP-CLR C++ rewrite — default: link vendored samtools/htslib (git submodule @ 1.24)
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp
CPPFLAGS += -Iinclude

# Vendored htslib: third_party/htslib → https://github.com/samtools/htslib tag 1.24
HTS_SRC := third_party/htslib
HTS_LIB := $(HTS_SRC)/libhts.a

# make USE_SYSTEM_HTS=1  → use pkg-config / system htslib instead
USE_SYSTEM_HTS ?= 0

ifeq ($(USE_SYSTEM_HTS),1)
  HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
  HTS_LIBS := $(shell pkg-config --libs htslib 2>/dev/null)
  ifeq ($(HTS_LIBS),)
    HTS_CFLAGS := -I/usr/local/include
    HTS_LIBS := -L/usr/local/lib -lhts -Wl,-rpath,/usr/local/lib
  endif
  LDFLAGS += $(HTS_LIBS) -fopenmp -lm -lgsl -lgslcblas
  HTS_REQ :=
else
  HTS_CFLAGS := -I$(HTS_SRC)
  # static archive; pull link flags from htslib's generated mk if present
  -include $(HTS_SRC)/htslib_static.mk
  HTSLIB_static_LIBS ?= -lz -lbz2 -llzma -lcurl -lcrypto -ldeflate -lpthread -lm
  LDFLAGS += $(HTS_LIB) $(HTSLIB_static_LIBS) -fopenmp -lgsl -lgslcblas
  HTS_REQ := $(HTS_LIB)
endif

CPPFLAGS += $(HTS_CFLAGS)

SRC := src/main.cpp src/util.cpp src/pop.cpp src/vcf_io.cpp src/xpclr.cpp
OBJ := $(SRC:.cpp=.o)
BIN := xpclr

.PHONY: all clean test-help htslib htslib-clean distclean

all: $(BIN)

htslib: $(HTS_LIB)

# First-time: submodule must exist; then autoheader/autoconf/configure + lib-static
$(HTS_LIB) $(HTS_SRC)/htslib_static.mk:
	@if [ ! -e $(HTS_SRC)/htslib/vcf.h ] && [ ! -e $(HTS_SRC)/vcf.h ]; then \
	  echo "[E] $(HTS_SRC) incomplete. Run:"; \
	  echo "    git submodule update --init --recursive"; \
	  exit 1; \
	fi
	@if [ ! -f $(HTS_SRC)/config.mk ]; then \
	  echo "[I] configuring vendored htslib ($(HTS_SRC)) ..."; \
	  cd $(HTS_SRC) && \
	    if [ ! -f config.guess ]; then autoreconf -i || true; fi && \
	    if [ ! -x configure ]; then autoheader && autoconf; fi && \
	    ./configure; \
	fi
	$(MAKE) -C $(HTS_SRC) -j$$(nproc 2>/dev/null || echo 4) lib-static htslib_static.mk

$(BIN): $(HTS_REQ) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

src/%.o: src/%.cpp include/xpclr.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

# Objects need headers from htslib; ensure library (and generated headers) exist first
ifneq ($(USE_SYSTEM_HTS),1)
$(OBJ): $(HTS_LIB)
endif

clean:
	rm -f $(OBJ) $(BIN)

htslib-clean:
	@if [ -f $(HTS_SRC)/Makefile ]; then $(MAKE) -C $(HTS_SRC) clean; fi

distclean: clean htslib-clean
	@if [ -f $(HTS_SRC)/config.mk ]; then $(MAKE) -C $(HTS_SRC) distclean || true; fi

test-help: $(BIN)
	./$(BIN) -h
	./$(BIN) -v
