# Build Garfield++ signal demo
# 用法:  source /root/gpp/install/share/Garfield/setupGarfield.sh  之后:
#        make
# 它会在 GARFIELD_INSTALL 环境变量里找头文件和库

GARFIELD_DIR ?= /root/gpp/install
ROOTSYS_FALLBACK ?= $(shell root-config --prefix 2>/dev/null)

CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall \
            -I$(GARFIELD_DIR)/include \
            $(shell root-config --cflags)

LDFLAGS  := -L$(GARFIELD_DIR)/lib -lGarfield \
            $(shell root-config --libs) \
            -lGeom -lEve -lEG \
            -Wl,-rpath,$(GARFIELD_DIR)/lib

TARGET   := signal_demo
SRC      := signal_demo.cc

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) signal_garfield.csv

.PHONY: all clean
