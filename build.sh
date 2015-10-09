#!/bin/bash
g++ -std=c++1y src/ex_to_junc.cpp -I"deps/include/" -L"deps/lib/" -lbamtools -o bin/ex_to_junc
