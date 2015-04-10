swig -Wall -c++ -python -o battlesim_wrap.cpp battlesim.i
g++ -I /usr/include/python2.7 -I ../include -fpic -c battlesim_wrap.cpp
g++ -shared battlesim_wrap.o -o _battlesim.so
