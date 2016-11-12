SHELL    = /bin/sh

-include build.properties

common := ./src/main/cxx
bindir := $(wildcard ~/bin)
srcdir := ./src/main/cxx
symeig := symeig

VPATH = $(srcdir):$(common)
.PHONY : all fit test clean install

all : fit
fit : rq-edfit rq-evfit
example : rq-edfit
	./src/test/bin/example.sh
install :
	mv -f ./rq-* $(bindir)
clean :
	rm -f ./*.o
distclean : clean
	rm -f ./rq-*
test : example
	diff src/test/resources/synth.html example.html

rq-edfit : edfit.o profiles.o readline.o section.o $(symeig).o
	$(CXX) $(LDFLAGS) $(VECLIB) -o $@ $< profiles.o readline.o section.o $(symeig).o
rq-evfit : evfit.o profiles.o readline.o section.o $(symeig).o
	$(CXX) $(LDFLAGS) $(VECLIB) -o $@ $< profiles.o readline.o section.o $(symeig).o

edfit.o : model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h $(symeig).h
evfit.o : model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h $(symeig).h

profiles.o : profiles.h
readline.o : readline.h
section.o : section.h
$(symeig).o : $(symeig).h

%.o : %.cxx
	$(CXX) -I$(common) -c $(CXXFLAGS) $< -o $@
