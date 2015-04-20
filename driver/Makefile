CPP = g++
SWIG = swig
SWIGFLAGS = -Wall -c++ -python
CXXFLAGS = -c -O3 -funroll-loops -march=native -fpic -Wall -I ../include -I /usr/include/python2.7
#CXXFLAGS = -c -g -fpic -Wall -I ../include -I /usr/include/python2.7
LDFLAGS = -shared
INTERFACE = battlesim
WRAPPER = $(INTERFACE)_wrap

_$(INTERFACE).so: $(WRAPPER).o
		$(CPP) $(LDFLAGS) -o $@ $<

$.o: %.cpp
		$(CPP) $(CXXFLAGS) $< -o $@

$(WRAPPER).cpp: $(INTERFACE).i
		$(SWIG) $(SWIGFLAGS) -o $@ $<

clean:
	rm -f $(WRAPPER).* _$(INTERFACE).so
