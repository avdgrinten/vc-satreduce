#pragma once

#include <graph.hpp>
#include <iostream>

class Input_Reader {
public:
    void findTwoLongs(int &x, int &y, std::string s);
    std::string readString(std::istream &in);
    void readTwoLongs(int &x, int &y, std::istream &in);
    Graph read_instance(std::istream &in);
};
