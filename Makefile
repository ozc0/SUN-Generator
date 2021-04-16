SUNGEN_VERSION := "0.1"
SUNGEN_UPDATE := "April 06, 2021"
SUNGEN_DEBUG := 0
BUILD_DATE := "$(shell date)"
CC=g++
CFLAGS =  -O3 -funroll-loops -g -I htslib -DRANDOM_VERSION=\"$(RANDOM_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DRANDOM_UPDATE=\"$(RANDOM_UPDATE)\" -DRANDOM_DEBUG=$(RANDOM_DEBUG)
LDFLAGS = htslib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl
NOCRAMFLAGS = htslib/libhts.a -lz -lm -lpthread
SOURCES = sun_generator.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = sun_gen
INSTALLPATH = /usr/local/bin/

.PHONY: htslib

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

sun_gen: htslib

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	make clean -C htslib
	rm -f $(EXECUTABLE) *.o *~

nocram: $(OBJECTS)
	cd htslib && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make && cd ..
	$(CC) $(OBJECTS) -o $(EXECUTABLE)-nocram $(NOCRAMFLAGS)

libs:
	make -C htslib

install:
	cp $(EXECUTABLE) $(INSTALLPATH)
