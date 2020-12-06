#!/usr/bin/env bash
if g++ -std=c++0x -pthread main.cpp cxstring.cpp readgenome.cpp reachtools.cpp -o ./reachtools
then
    echo make_finished
fi
