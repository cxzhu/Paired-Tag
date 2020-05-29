//
// Created by Chenxu Zhu on 6/29/18.
//

#include "reachtools.h"
#include "cxstring.h"


char str2char(string str);
int align_score(string str1, string str2);
string get_read_type(string typ);
string get_bc_number(string bc, int score); //, map<string, string> *bc_library);
string get_bc_number_4(string bc, int score); //, map<string, string> *bc_library);





/// Global functions

string reachtools:: extract_barcode(string seq){ //, map<string, string> bc_library){

    string umi = seq.substr(0,10);
    string bc1 = seq.substr(10,7);
    string bc2 = seq.substr(47,7);
    string bc3 = seq.substr(84,7);
    string typ = seq.substr(121,2);
    string bc4 = seq.substr(123,3);

    int bc_align_score = 5;

    string bc_number1 = get_bc_number(bc1, bc_align_score);
    string bc_number2 = get_bc_number(bc2, bc_align_score);
    string bc_number3 = get_bc_number(bc3, bc_align_score);
    string bc_number4 = get_bc_number_4(bc4, 3);

    string type = get_read_type(typ);

    string read_id = type + ':' + bc_number1 + ':' + bc_number2 + ':' + bc_number3 + ':' + bc_number4 + ':' + umi;
    return(read_id);
}

void reachtools:: docking(string & read, string & quality){
    string model = "NNNNNNNNNNNNNNNNNGTGGCCGATGTTTCGCATCGGCGTACGACTNNNNNNNGGATTCGAGGAGCGTGTGCGAACTCAGACCNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCGNNNNNNN";
    //r04 linker
    //string bait;
    int score = 0;
    int cur_score = 0;
    for(int i = 17; i < 32; ++i){
        if(model.substr(i, 1) == read.substr(i, 1)) ++cur_score;
        if(model.substr(i, 1) == read.substr(i - 1, 1)) ++score;
    }
    if(score > cur_score){
        read = "N" + read;
        quality = "I" + quality;
        cur_score = score;
    }
    if(cur_score < 10){
        read = "Low Quality";
        return;
    }

    //r03 linker
    score = 0;
    cur_score = 0;
    for(int i = 32; i < 47; ++i){
        if(model.substr(i, 1) == read.substr(i, 1))++cur_score;
        if(model.substr(i, 1) == read.substr(i - 1, 1))++score;
    }
    if(score > cur_score){
        read = read.substr(0, 32) + "N" + read.substr(32, read.length()-32);
        quality = quality.substr(0, 32) + "I" + read.substr(32, quality.length()-32);
        cur_score = score;
    }
    if(cur_score < 9){
        read = "Low Quality";
        return;
    }

    //r02 linker 69
    score = 0;
    cur_score = 0;
    for(int i = 69; i < 84; ++i){
        if(model.substr(i, 1) == read.substr(i, 1))++cur_score;
        if(model.substr(i, 1) == read.substr(i - 1, 1))++score;
    }
    if(score > cur_score){
        read = read.substr(0, 69) + "N" + read.substr(69, read.length()-69);
        quality = quality.substr(0, 69) + "I" + read.substr(69, quality.length()-69);
        cur_score = score;
    }
    if(cur_score < 8){
        read = "Low Quality";
        return;
    }

    //r01 probe 106
    score = 0;
    cur_score = 0;
    for(int i = 106; i < 121; ++i){
        if(model.substr(i, 1) == read.substr(i, 1))++cur_score;
        if(model.substr(i, 1) == read.substr(i - 1, 1))++score;
    }
    if(score > cur_score){
        read = read.substr(0, 106) + "N" + read.substr(106, read.length()-106);
        quality = quality.substr(0, 106) + "I" + read.substr(106, quality.length()-106);
        cur_score = score;
    }
    if(cur_score < 7){
        read = "Low Quality";
        return;
    }

    //quality Q33 (int)'Base' >= 63
    for(int i = 0; i < 7; ++i){
        if((int)(quality[i+10]) < 53 || (int)(quality[i+47]) < 53 || (int)(quality[i+84]) < 53 || (int)(quality[i+121]) < 53)quality="Low Quality";
    }
    return;
}


///// Local functions

int align_score(string str1, string str2){
    /// give the str align score
    /// aligned, +2; misaligned, -2; N contain, 0

    //if(str1.length() != str2.length())return -1;
    int score = 0;
    for(int i = 0; i < str1.length(); ++i){
        if(str2[i] == 'N')continue;
        str2[i] == str1[i] ? (score += 2) : (score -= 1);
    }
    return score;
}

char str2char(string str) {
    stringstream tmp;
    tmp << str;
    char output;
    tmp >> output;
    return(output);
}

string get_read_type(string typ){
    /// classify the read to DNA("d") or RNA("r") read
    /// if cannot distinguish, return "n"

    int score = 0;
    if(typ[0] == 'A')++score;
    if(typ[1] == 'G')++score;
    if(score == 2){
        return("d");
    }
    int score_rna = 0;
    if(typ[0] == 'T')++score_rna;
    if(typ[1] == 'C')++score_rna;
    if(score_rna == 2){
        return("r");
    }

    if(score_rna == score) {
        return ("n");
    }
    else if(score != 0){
        return ("d");
    }
    else{
        return ("r");
    }
}

//string get_bc_number(string bc, int score){
string get_bc_number(string bc, int score){
//string get_bc_number(string bc, int score, map<string, string> *bc_library){
    string bc_number = "00";
    int cur_score = 0;
    int valid = 0;
    //int score = 0;

    string bc_library[96] = {"AAACAAC","AAACCGG","AAACGTC","AAAGATG","AAATCCA","AAATGAG","AACACTG","AACGTTT","AAGAAGC","AAGCCCT","AAGCTAC","AATCTTG","ACAACAC","ACAGTAT","ACCAAGT","ACCCTAA","ACCCTTT","ACCTCTC","ACGATTG","ACGCAGA","ACGTAAA","ACTACCT","ACTCGGT","ACTGTCG","ACTTATG","AGAAAGG","AGAATCT","AGACATA","AGAGACC","AGCCCAA","AGCTATT","AGGAGGT","AGGGCTT","AGGTGTA","AGTGCTC","AGTGGGA","AGTTACG","ATAAGGG","ATCATTC","ATGGAAC","ATGTGCC","ATTCACC","ATTCGAG","CAAGCCT","CACAAGG","CACCTTA","CAGAGTG","CAGCGAA","CAGGTCA","CATAACT","CATATCG","CATCGAT","CATTACA","CATTTCC","CCAAATG","CCACTTG","CCGGATA","CCGGTTT","CCTAAGA","CCTAGTC","CCTGCAA","CGACGTT","CGAGTAA","CGATTAT","CGTAGCA","CGTCTGA","CTACAGC","CTCAATA","CTCGTTG","CTCTACG","CTTGGGT","GAAACTC","GACTGTC","GATACAG","GCGATCA","GCGTACT","GCTCGAA","GGAAGAA","GGAGATT","GGGCTAA","GGGTATG","GGTAACC","GGTAGTG","GGTGAAA","GTAATCG","GTATAAG","GTCAGAC","GTCCCTT","GTGCCAT","GTGGTCT","GTTCTCC","GTTGCTT","TACCCGA","TAGACGA","TAGTCAC","TCACATC"};
    int t = 96;
    if(bc.length() == 5)t = 12;
    int id = 0;
    for(int i = 0; i < t; ++i){
        string str_target = bc_library[i];

        cur_score = align_score(bc, str_target);
        if(cur_score < score){
            continue;
        }
        else if(cur_score == score){
            valid = 0;
        }
        else{
            valid = 1;
            score = cur_score;
            id = i;
            if((cur_score == 14 && t == 96) || (cur_score == 10 && t == 12)){
                break;
            }
        }
    }

    if(valid == 0)bc_number = "00";
    if(id < 9){
        bc_number = cxstring::int2str(id + 1);
        bc_number = "0" + bc_number;
    }
    else{
        bc_number = cxstring::int2str(id + 1);
    }
    return bc_number;
}



//string get_bc_number(string bc, int score){
string get_bc_number_4(string bc, int score){
//string get_bc_number(string bc, int score, map<string, string> *bc_library){
    string bc_number = "00";
    int cur_score = 0;
    int valid = 0;
    string bc_library[12] = {"CATC", "ATGA", "AGCT", "ACAG", "GAAT","TACG","TTAC","GTTG","CCGT","CGAA","TCTA","GGGC"};
    int t = 8;
    int id = 0;
    for(int i = 0; i < t; ++i){
        string str_target = bc_library[i];
        cur_score = align_score(bc, str_target);
        if(cur_score < score){
            continue;
        }
        else if(cur_score == score){
            valid = 0;
        }
        else{
            valid = 1;
            score = cur_score;
            id = i;
            if(cur_score == 6){   //(cur_score == 6 && t == 96) || (cur_score == 6 && t == 12)){
                break;
            }
        }
    }

    if(valid == 0)bc_number = "00";
    if(id < 9){
        bc_number = cxstring::int2str(id + 1);
        bc_number = "0" + bc_number;
    }
    else{
        bc_number = cxstring::int2str(id + 1);
    }
    return bc_number;
}


