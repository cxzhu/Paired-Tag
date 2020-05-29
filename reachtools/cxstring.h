//
// Created by Chenxu Zhu on 6/29/18.
//

#ifndef REACHTOOLS_CXSTRING_H
#define REACHTOOLS_CXSTRING_H

#include <algorithm>
#include "readgenome.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>

using namespace std;

class cxstring {
private:

public:
    static vector<string> split(const string &s, const string &seperator);
    static string chomp(string str);
    static bool is_bam_header(const std::string& buffer);
    static int str2int(string str);
    static string int2str(int i);


};


class samline{
public:
    string readname, chr, cigar, chrnext, seq, qual, other;
    int flag, pos, mapq, posnext, plen;
        
    void init(string str){
        str = cxstring::chomp(str);
        vector<string> tmp = cxstring::split(str, "\t");
        readname = tmp[0];
        flag = cxstring::str2int(tmp[1]);
        chr = tmp[2];
        pos = cxstring::str2int(tmp[3]);
        mapq = cxstring::str2int(tmp[4]);
        cigar = tmp[5];
        chrnext = tmp[6];
        posnext = cxstring::str2int(tmp[7]);
        plen = cxstring::str2int(tmp[8]);
        seq = tmp[9];
        qual = tmp[10];
        if(tmp.size() > 11){
        	other = tmp[11];
        }
        else{
        	other = "";
        }
        return;
    }
    void empty(){
        readname = "";
        flag = 0;
        chr = "";
        pos = 0;
        mapq = 0;
        cigar = "";
        chrnext = "";
        posnext = 0;
        plen = 0;
        seq = "";
        qual = "";
        return;
    }
    void write(FILE * outfile){
    	string out = readname + "\t" + cxstring::int2str(flag) + "\t" + chr + "\t" + cxstring::int2str(pos) + "\t" +  cxstring::int2str(mapq) + "\t" + cigar + "\t" + chrnext +  "\t" + cxstring::int2str(posnext) + "\t" + cxstring::int2str(plen) + "\t" + seq + "\t" + qual + "\t" + other;
    	fputs((out+"\n").c_str(), outfile);
    }
};

class fqline{
public:
	string readname,seq,mark,qual;

	void read_part_record(FILE * infile, string rn){
		readname = rn;
		char buffer[2000];
		stringstream tmp;
		fgets(buffer, sizeof(buffer), infile); // seq
		tmp << buffer;
		tmp >> seq;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // mark
		tmp << buffer;
		tmp >> mark;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // qual
		tmp << buffer;
		tmp >> qual;
		return;
	}
	void read_full_record(FILE * infile){
		char buffer[2000];
		stringstream tmp;
		fgets(buffer, sizeof(buffer), infile); // readname
		tmp << buffer;
		tmp >> readname;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // seq
		tmp << buffer;
		tmp >> seq;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // mark
		tmp << buffer;
		tmp >> mark;
		tmp.str("");
		fgets(buffer, sizeof(buffer), infile); // qual
		tmp << buffer;
		tmp >> qual;
		return;
	}
	void write_record(FILE * outfile){
		fputs((readname + "\n").c_str(), outfile);
		fputs((seq + "\n").c_str(), outfile);
		fputs((mark + "\n").c_str(), outfile);
		fputs((qual + "\n").c_str(), outfile);
	}
};


#endif //REACHTOOLS_CXSTRING_H
