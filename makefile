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

DIFF    = diff

-include build.properties

bindir := $(wildcard ~/bin)
srcdir := src/main
tstdir := src/test

VPATH = $(srcdir)/cxx:$(tstdir)/cxx

bin := especid
bin += especiv
bin += especia
bin += xtractcom
bin += xtractdat
bin += xtractlog
bin += xtractmes
bin += xtractmod


especid : especid.o profiles.o readline.o section.o symeig.o
	$(CXX) $(LDFLAGS) -o $@ $< profiles.o readline.o section.o symeig.o
especiv : especiv.o profiles.o readline.o section.o symeig.o
	$(CXX) $(LDFLAGS) -o $@ $< profiles.o readline.o section.o symeig.o
especia : especia.o profiles.o readline.o section.o symeig.o
	$(CXX) $(LDFLAGS) -o $@ $< profiles.o readline.o section.o symeig.o

% : %.cxx
	$(CXX) $(CXXFLAGS) -o $@ $<

especid.o : especid.cxx model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h symeig.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
especiv.o : especiv.cxx model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h symeig.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
especia.o : especia.cxx model.h mtwister.h optimize.h profiles.h randev.h readline.h section.h symeig.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.o : %.cxx %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.html : $(tstdir)/resources/%.html
	./xtractmod < $< | `./xtractcom < $<` > $@
%.diff : $(tstdir)/resources/%.html %.html
	$(DIFF) $^


.PHONY : all clean distclean install test

all : $(bin)

install :
	mv -f $(bin) $(bindir)

clean :
	$(RM) *.o
	$(RM) *.html

distclean : clean
	$(RM) $(bin)

test : example.diff
