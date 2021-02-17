#include <algorithm>
#include <cassert>
#include <functional>
#include <input-reader.hpp>

// TODO: this code is *very* fragile and inefficient.

void Input_Reader::findTwoLongs(int &x, int &y, std::string s) {
    long i = 0;
    long sz = s.size();

    while (i < sz && s[i] != ' ')
        i++;

    assert(i < sz - 1);
    assert(i > 0);
    assert(s[i] == ' ');

    x = stoll(s.substr(0, i));
    y = stoll(s.substr(i + 1));
}

std::string Input_Reader::readString(std::istream &in) {
    bool line = 0;
    std::string s;

    while (getline(in, s)) {
        if (s.size() < 2 || s[0] == 'c')
            continue;
        else {
            line = 1;
            break;
        }
    }

    if (!line)
        throw std::runtime_error("EOF unexpected");

    assert(s.length() > 5);

    return s;
}

void Input_Reader::readTwoLongs(int &x, int &y, std::istream &in) {
    bool line = 0;
    std::string s;

    while (getline(in, s)) {
        if (s.size() < 2 || s[0] == 'c')
            continue;
        else {
            line = 1;
            break;
        }
    }

    if (!line)
        throw std::runtime_error("EOF unexpected");

    findTwoLongs(x, y, s);
}

Graph Input_Reader::read_instance(std::istream &in) {
    std::string s = readString(in);
    s = s.substr(5);

    int n, m;
    findTwoLongs(n, m, s);

    std::vector<std::vector<int>> g(n, std::vector<int>());

    for (int i = 0; i < m; i++) {
        int x, y;
        readTwoLongs(x, y, in);
        x--;
        y--;
        g[x].push_back(y);
        g[y].push_back(x);
    }

    return Graph(g);
}
