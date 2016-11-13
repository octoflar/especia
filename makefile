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

bindir := $(wildcard ~/bin)
srcdir := src/main
tstdir := src/test

VPATH = $(srcdir)/cxx:$(tstdir)/cxx

bin := edfit
bin += evfit
bin += mmfit
bin += xtractdat
bin += xtractmes
bin += xtractlog

edfit : edfit.o profiles.o readline.o section.o symeig.o
	$(CXX) $(LDFLAGS) $(VECLIB) -o $@ $< profiles.o readline.o section.o symeig.o
evfit : evfit.o profiles.o readline.o section.o symeig.o
	$(CXX) $(LDFLAGS) $(VECLIB) -o $@ $< profiles.o readline.o section.o symeig.o
mmfit : mmfit.o profiles.o readline.o section.o symeig.o
	$(CXX) $(LDFLAGS) $(VECLIB) -o $@ $< profiles.o readline.o section.o symeig.o

edfit.o : edfit.cxx model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h symeig.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
evfit.o : evfit.cxx model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h symeig.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
mmfit.o : mmfit.cxx model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h symeig.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

% : %.cxx
	$(CXX) -o $@ $(CXXFLAGS) $<

%.o : %.cxx %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.html : $(tstdir)/resources/%.in
	$(tstdir)/bin/$(@:.html=.sh)


.PHONY : all clean distclean install test

all : $(bin)

install :
	mv -f $(bin) $(bindir)

clean :
	$(RM) *.o
	$(RM) *.html

distclean : clean
	$(RM) $(bin)

test : example.html
	diff $(tstdir)/resources/example.html ./example.html
