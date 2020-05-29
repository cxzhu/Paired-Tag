//
// Created by Chenxu Zhu on 6/29/18.
//

#include "readgenome.h"


using namespace std;
/*
readgenome::readgenome(map<string, string> *ha){
    this->hash = ha;
}
*/



void readgenome:: hg19(map<string, string> *hash){
    // hg19 path in SLS server is "/gpfs/user/zhuchenxu/genome_ref/hg19/genome.fa";
    string s1 = "cat ";
    string s2 = "/mnt/tscc/chz272/genome_ref/hg19/hg19.fa"; //"/gpfs/user/zhuchenxu/genome_ref/hg19/genome.fa";
    {
        ifstream iftest;
        iftest.open(s2);
        if(!iftest){
            s2 = "/mnt/tscc/chz272/genome_ref/hg19/hg19.fa";
            iftest.close();
            iftest.open(s2);
            if(!iftest) {
                cout << "Cannot open genome file... please specify a correct path." << endl;
                exit(0);
            }
            iftest.close();
        }
    }

    string s3 = s1 + s2;
    const char* command = s3.c_str();

    // prepare for the file stream
    FILE * in_file;
    in_file = popen(command, "r");

    char buffer[1000];
    float line_number = 0;
    int mega_line_number = 0;
    string chr;

    while(fgets(buffer, sizeof(buffer), in_file)){
        // report reading progress
        ++line_number;
        if((line_number / 10000000) == int(line_number / 10000000)){
            ++mega_line_number;
            line_number = 0;
            cout << "[Thread readgenome hg19]: " << mega_line_number << "0 M lines loaded..." << endl;
        }

        // processing each line
        stringstream tmp(buffer);
        string cur_line;
        tmp >> cur_line;
        if(cur_line.substr(0, 1) == ">"){
            chr = cur_line.substr(1, cur_line.length() - 1);
        }
        else if(chr != ""){
            map<string, string>::iterator it;
            it = hash->find(chr);
            if(it == hash->end()){
                hash->insert(make_pair(chr, cur_line));
            }
            else{
                it->second.append(cur_line);
            }
        }
    }
    pclose(in_file);
    cout << "[Thread readgenome hg19]: Genome hg19 read finished." << endl;
}

void readgenome:: hg19_transcriptome(map<string, string> *hash){
    // hg19 transcriptome path in SLS serve is "/gpfs/user/zhuchenxu/genome_ref/hg19/transcriptome/genes.fa";
    string s1 = "cat ";
    string s2 = "/mnt/tscc/chz272/genome_ref/hg19/mrna.fa";
    {
        ifstream iftest;
        iftest.open(s2);
        if(!iftest){
            cout << "Cannot open genome file... please specify a correct path." << endl;
            exit(0);
        }
        iftest.close();
    }
    string s3 = s1 + s2;
    const char* command = s3.c_str();

    // prepare for the file stream
    FILE * in_file;
    in_file = popen(command, "r");

    char buffer[1000];
    float line_number = 0;
    int mega_line_number = 0;
    string chr;

    while(fgets(buffer, sizeof(buffer), in_file)){
        // report reading progress
        ++line_number;
        if((line_number / 10000000) == int(line_number / 10000000)){
            ++mega_line_number;
            line_number = 0;
            cout << "[Thread readgenome hg19 transcriptome]: " << mega_line_number << "0 M lines loaded..." << endl;
        }

        // processing each line
        stringstream tmp(buffer);
        string cur_line;
        tmp >> cur_line;
        if(cur_line.substr(0, 1) == ">"){
            chr = cur_line.substr(1, cur_line.length() - 1);
        }
        else if(chr != ""){
            map<string, string>::iterator it;
            it = hash->find(chr);
            if(it == hash->end()){
                hash->insert(make_pair(chr, cur_line));
            }
            else{
                it->second.append(cur_line);
            }
        }
        // if(chr != "chr1")break;
    }
    pclose(in_file);
    cout << "[Thread readgenome hg19 transcriptome]: Genome hg19 transcriptome read finished." << endl;

}

void readgenome:: mm9(map<string, string> *hash){
    // mm9 path in SLS server is "/gpfs/user/zhuchenxu/genome_ref/mm9/mm9.fa";
    string s1 = "cat ";
    string s2 = "/mnt/tscc/chz272/genome_ref/mm9/mm9.fa";
    {
        ifstream iftest;
        iftest.open(s2);
        if(!iftest){
            iftest.close();
            s2 = "/mnt/tscc/chz272/genome_ref/mm9/mm9.fa";
            iftest.open(s2);
            if(!iftest){
                cout << "Cannot open genome file... please specify a correct path." << endl;
                exit(0);
            }
        }
        iftest.close();

    }
    string s3 = s1 + s2;
    const char* command = s3.c_str();

    // prepare for the file stream
    FILE * in_file;
    in_file = popen(command, "r");

    char buffer[1000];
    float line_number = 0;
    int mega_line_number = 0;
    string chr;

    while(fgets(buffer, sizeof(buffer), in_file)){
        // report reading progress
        ++line_number;
        if((line_number / 10000000) == int(line_number / 10000000)){
            ++mega_line_number;
            line_number = 0;
            cout << "[Thread readgenome mm9]: " << mega_line_number << "0 M lines loaded..." << endl;
        }

        // processing each line
        stringstream tmp(buffer);
        string cur_line;
        tmp >> cur_line;
        if(cur_line.substr(0, 1) == ">"){
            chr = cur_line.substr(1, cur_line.length() - 1);
        }
        else if(chr != ""){
            map<string, string>::iterator it;
            it = hash->find(chr);
            if(it == hash->end()){
                hash->insert(make_pair(chr, cur_line));
            }
            else{
                it->second.append(cur_line);
            }
        }
    }
    pclose(in_file);
    cout << "[Thread readgenome mm9]: Genome mm9 read finished." << endl;

}

void readgenome:: mm10(map<string, string> *hash){
    // mm10 path in SLS serve is "/gpfs/user/zhuchenxu/genome_ref/mm10/mm10.fa";
    string s1 = "cat ";
    string s2 = "/mnt/tscc/chz272/genome_ref/mm10/mm10.fa";
    {
        ifstream iftest;
        iftest.open(s2);
        if(!iftest){
            iftest.close();
            s2 = "/mnt/tscc/chz272/genome_ref/mm10/mm10.fa";
            iftest.open(s2);
            if(!iftest){
                cout << "Cannot open genome file... please specify a correct path." << endl;
                exit(0);
            }
        }
        iftest.close();
    }
    string s3 = s1 + s2;
    const char* command = s3.c_str();

    // prepare for the file stream
    FILE * in_file;
    in_file = popen(command, "r");

    char buffer[1000];
    float line_number = 0;
    int mega_line_number = 0;
    string chr;

    while(fgets(buffer, sizeof(buffer), in_file)){
        // report reading progress
        ++line_number;
        if((line_number / 10000000) == int(line_number / 10000000)){
            ++mega_line_number;
            line_number = 0;
            cout << "[Thread readgenome mm10]: " << mega_line_number << "0 M lines loaded..." << endl;
        }

        // processing each line
        stringstream tmp(buffer);
        string cur_line;
        tmp >> cur_line;
        if(cur_line.substr(0, 1) == ">"){
            chr = cur_line.substr(1, cur_line.length() - 1);
        }
        else if(chr != ""){
            map<string, string>::iterator it;
            it = hash->find(chr);
            if(it == hash->end()){
                hash->insert(make_pair(chr, cur_line));
            }
            else{
                it->second.append(cur_line);
            }
        }
    }
    pclose(in_file);
    cout << "[Thread readgenome mm10]: Genome mm10 read finished." << endl;

}

void readgenome:: custom(map<string, string> *hash, string genome_path){
    // hg19 path in SLS server is "/gpfs/user/zhuchenxu/genome_ref/hg19/genome.fa";
    string s1 = "cat ";
    //   string s2 = "/gpfs/user/zhuchenxu/genome_ref/hg19/genome.fa";
    string s2 = genome_path;
    if(genome_path == "hg19"){
        s2 = "/mnt/tscc/chz272/genome_ref/hg19/hg19.fa"; //"/gpfs/user/zhuchenxu/genome_ref/hg19/genome.fa";
        {
            ifstream iftest;
            iftest.open(s2);
            if(!iftest){
                s2 = "/mnt/tscc/chz272/genome_ref/hg19/hg19.fa";
                iftest.close();
                iftest.open(s2);
                if(!iftest) {
                    cout << "Cannot open genome file... please specify a correct path." << endl;
                    exit(0);
                }
                iftest.close();
            }
        }
    }
    else if(genome_path == "hg19_transcriptome"){
        s2 = "/mnt/tscc/chz272/genome_ref/hg19/mrna.fa";
        {
            ifstream iftest;
            iftest.open(s2);
            if(!iftest){
                cout << "Cannot open genome file... please specify a correct path." << endl;
                exit(0);
            }
            iftest.close();
        }
    }
    else if(genome_path == "mm9"){
        s2 = "/mnt/tscc/chz272/genome_ref/mm9/mm9.fa";
        {
            ifstream iftest;
            iftest.open(s2);
            if(!iftest){
                iftest.close();
                s2 = "/mnt/tscc/chz272/genome_ref/mm9/mm9.fa";
                iftest.open(s2);
                if(!iftest){
                    cout << "Cannot open genome file... please specify a correct path." << endl;
                    exit(0);
                }
            }
            iftest.close();

        }
    }
    else if(genome_path == "mm10"){
        s2 = "/mnt/tscc/chz272/genome_ref/mm10/mm10.fa";
        {
            ifstream iftest;
            iftest.open(s2);
            if(!iftest){
                iftest.close();
                s2 = "/mnt/tscc/chz272/genome_ref/mm10/mm10.fa";
                iftest.open(s2);
                if(!iftest){
                    cout << "Cannot open genome file... please specify a correct path." << endl;
                    exit(0);
                }
            }
            iftest.close();
        }
    }

    {
        ifstream iftest;
        iftest.open(s2);
        if(!iftest){
            cout << "Cannot open genome file... please specify a correct path." << endl;
            exit(0);
        }
    }
    string s3 = s1 + s2;
    const char* command = s3.c_str();

    // prepare for the file stream
    FILE * in_file;
    in_file = popen(command, "r");

    char buffer[1000];
    float line_number = 0;
    int mega_line_number = 0;
    string chr;

    while(fgets(buffer, sizeof(buffer), in_file)){
        // report reading progress
        ++line_number;
        if((line_number / 10000000) == int(line_number / 10000000)){
            ++mega_line_number;
            line_number = 0;
            cout << "[Thread readgenome custom]: " << mega_line_number << "0 M lines loaded..." << endl;
        }

        // processing each line
        stringstream tmp(buffer);
        string cur_line;
        tmp >> cur_line;
        if(cur_line.substr(0, 1) == ">"){
            chr = cur_line.substr(1, cur_line.length() - 1);
        }
        else if(chr != ""){
            map<string, string>::iterator it;
            it = hash->find(chr);
            if(it == hash->end()){
                hash->insert(make_pair(chr, cur_line));
            }
            else{
                it->second.append(cur_line);
            }
        }
    }
    pclose(in_file);
    cout << "[Thread readgenome custom]: Genome custom read finished." << endl;
}
string readgenome:: uc(string input){
    for(int i = 0; i < input.length(); i++){
        input[i] = toupper(input[i]);
    }
    return input;
}
string readgenome:: reverse_complmentary(string input) {
    reverse(input.begin(), input.end());
    for(int i = 0; i < input.length(); i++){
        if(input[i] == 'A'){
            input[i] = 'T';
        }
        else if(input[i] == 'C'){
            input[i] = 'G';
        }
        else if(input[i] == 'T'){
            input[i] = 'A';
        }
        else if(input[i] == 'G'){
            input[i] = 'C';
        }
    }
    return input;
}


string readgenome:: int2str(int input){
    stringstream tmp;
    string out;
    tmp << input;
    tmp >> out;
    return out;
}
