//
// Created by Chenxu Zhu on 6/29/18.
//

#ifndef REACHTOOLS_READGENOME_H
#define REACHTOOLS_READGENOME_H

#include <string>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;


class readgenome {
private:

public:
    static void hg19(map<string, string> *hash);
    static void hg19_transcriptome(map<string, string> *hash);
    static void mm9(map<string, string> *hash);
    static void mm10(map<string, string> *hash);
    static void custom(map<string, string> *hash, string genome_path);
    static string uc(string input);
    static string reverse_complmentary(string input); // reverse complementary string, need to be uc first
    static string int2str(int input); ///convert int to string
};



#endif //REACHTOOLS_READGENOME_H
