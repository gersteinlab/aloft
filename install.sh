#!/bin/sh

tar -zxvf Python-2.7.3.tgz
cd Python-2.7.3
./configure CFLAGS=-fPIC
make
cd ..
gcc gerprate.c -IPython-2.7.3/ -IPython-2.7.3/Include/ -LPython-2.7.3/ -l python2.7 -lpthread -lm -ldl -lutil -fPIC -shared -o gerprate.so
#rm -r Python-2.7.3
