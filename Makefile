CXX = g++
CXXFLAGS = -std=c++17 -pthread -Wall -Wextra
OPTFLAGS = -O3
INCLUDES = -I.
LDFLAGS = -pthread

# Debug flags (can be enabled with make DEBUG=1)
ifdef DEBUG
    CXXFLAGS += -g -DDEBUG
    OPTFLAGS = -O0
endif

# Source files
MAIN_SRC = virtual_centrifuge_kmc.cpp
KMC_API_SRC = kmc_api/kmc_file.cpp kmc_api/kmer_api.cpp kmc_api/mmer.cpp
SRC = $(MAIN_SRC) $(KMC_API_SRC)

# Output binary
TARGET = CsCl

# Default target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(INCLUDES) $^ -o $@ $(LDFLAGS)

# Clean target
clean:
	rm -f $(TARGET)
	rm -f *_kmc_input.list
	rm -rf *_kmc_tmp

# Install target (adjust as needed)
install: $(TARGET)
	mkdir -p $(DESTDIR)/usr/local/bin
	install -m 755 $(TARGET) $(DESTDIR)/usr/local/bin/

.PHONY: all clean install
