# XP-CLR C++ rewrite
# Vendored: third_party/htslib @ 1.24, third_party/gsl @ v2.7.0
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter -fopenmp
CPPFLAGS += -Iinclude

HTS_SRC := third_party/htslib
HTS_LIB := $(HTS_SRC)/libhts.a
GSL_SRC := third_party/gsl
GSL_LIB := $(GSL_SRC)/.libs/libgsl.a
GSL_CBLAS := $(GSL_SRC)/cblas/.libs/libgslcblas.a

# Optional: make USE_SYSTEM_HTS=1 / USE_SYSTEM_GSL=1
USE_SYSTEM_HTS ?= 0
USE_SYSTEM_GSL ?= 0

# ---- htslib ----
ifeq ($(USE_SYSTEM_HTS),1)
  HTS_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
  HTS_LIBS := $(shell pkg-config --libs htslib 2>/dev/null)
  ifeq ($(HTS_LIBS),)
    HTS_CFLAGS := -I/usr/local/include
    HTS_LIBS := -L/usr/local/lib -lhts
  endif
  HTS_REQ :=
else
  HTS_CFLAGS := -I$(HTS_SRC)
  -include $(HTS_SRC)/htslib_static.mk
  HTSLIB_static_LIBS ?= -lz -lbz2 -llzma -lcurl -lcrypto -ldeflate -lpthread -lm
  HTS_LIBS := $(HTS_LIB) $(HTSLIB_static_LIBS)
  HTS_REQ := $(HTS_LIB)
endif

# ---- gsl ----
ifeq ($(USE_SYSTEM_GSL),1)
  GSL_CFLAGS := $(shell pkg-config --cflags gsl 2>/dev/null)
  GSL_LIBS := $(shell pkg-config --libs gsl 2>/dev/null)
  ifeq ($(GSL_LIBS),)
    GSL_CFLAGS :=
    GSL_LIBS := -lgsl -lgslcblas
  endif
  GSL_REQ :=
else
  # headers live under third_party/gsl (e.g. gsl/gsl_integration.h via -I$(GSL_SRC))
  GSL_CFLAGS := -I$(GSL_SRC)
  GSL_LIBS := $(GSL_LIB) $(GSL_CBLAS) -lm
  GSL_REQ := $(GSL_LIB) $(GSL_CBLAS)
endif

CPPFLAGS += $(HTS_CFLAGS) $(GSL_CFLAGS)
LDFLAGS += $(HTS_LIBS) $(GSL_LIBS) -fopenmp

SRC := src/main.cpp src/util.cpp src/pop.cpp src/vcf_io.cpp src/xpclr.cpp
OBJ := $(SRC:.cpp=.o)
BIN := xpclr

.PHONY: all clean test-help htslib gsl distclean

all: $(BIN)

htslib: $(HTS_LIB)
gsl: $(GSL_LIB) $(GSL_CBLAS)

$(HTS_LIB) $(HTS_SRC)/htslib_static.mk:
	@if [ ! -e $(HTS_SRC)/htslib/vcf.h ] && [ ! -e $(HTS_SRC)/vcf.h ]; then \
	  echo "[E] $(HTS_SRC) incomplete. Run: git submodule update --init --recursive"; \
	  exit 1; \
	fi
	@if [ ! -f $(HTS_SRC)/config.mk ]; then \
	  echo "[I] configuring vendored htslib ..."; \
	  cd $(HTS_SRC) && \
	    if [ ! -f config.guess ]; then autoreconf -i || true; fi && \
	    if [ ! -x configure ]; then autoheader && autoconf; fi && \
	    ./configure; \
	fi
	$(MAKE) -C $(HTS_SRC) -j$$(nproc 2>/dev/null || echo 4) lib-static htslib_static.mk

$(GSL_LIB) $(GSL_CBLAS):
	@if [ ! -f $(GSL_SRC)/configure ] && [ ! -f $(GSL_SRC)/configure.ac ]; then \
	  echo "[E] $(GSL_SRC) incomplete. Run: git submodule update --init --recursive"; \
	  exit 1; \
	fi
	@if [ ! -f $(GSL_SRC)/config.status ]; then \
	  echo "[I] configuring vendored gsl ..."; \
	  cd $(GSL_SRC) && \
	    if [ ! -x configure ]; then autoreconf -fi || ./autogen.sh; fi && \
	    ./configure --disable-shared --enable-static; \
	fi
	$(MAKE) -C $(GSL_SRC) -j$$(nproc 2>/dev/null || echo 4)

$(BIN): $(HTS_REQ) $(GSL_REQ) $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

src/%.o: src/%.cpp include/xpclr.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

ifneq ($(USE_SYSTEM_HTS),1)
$(OBJ): $(HTS_LIB)
endif
ifneq ($(USE_SYSTEM_GSL),1)
$(OBJ): $(GSL_LIB)
endif

clean:
	rm -f $(OBJ) $(BIN)

distclean: clean
	@if [ -f $(HTS_SRC)/config.mk ]; then $(MAKE) -C $(HTS_SRC) distclean || true; fi
	@if [ -f $(GSL_SRC)/config.status ]; then $(MAKE) -C $(GSL_SRC) distclean || true; fi

test-help: $(BIN)
	./$(BIN) -h
	./$(BIN) -v
