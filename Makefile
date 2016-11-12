# Copyright (c) 2016 Ralf Quast
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

SHELL    = /bin/sh

-include build.properties

common := ./src/main/cxx
bindir := $(wildcard ~/bin)
srcdir := ./src/main/cxx
symeig := symeig

VPATH = $(srcdir):$(common)
.PHONY : all fit test clean install

all : fit
fit : rq-edfit rq-evfit rq-mmfit
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
rq-mmfit : mmfit.o profiles.o readline.o section.o $(symeig).o
	$(CXX) $(LDFLAGS) $(VECLIB) -o $@ $< profiles.o readline.o section.o $(symeig).o

edfit.o : model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h $(symeig).h
evfit.o : model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h $(symeig).h
mmfit.o : model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h $(symeig).h

profiles.o : profiles.h
readline.o : readline.h
section.o : section.h
$(symeig).o : $(symeig).h

%.o : %.cxx
	$(CXX) -I$(common) -c $(CXXFLAGS) $< -o $@
