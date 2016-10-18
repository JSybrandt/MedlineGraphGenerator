#include "Abstract.h"

Abstract::Abstract(string pmid, string data, Dict * dict)
{
    this -> pmid = pmid;
    this->dict = dict;
    int index = data.find(delim);
    key = data.substr(0, index);
    string values = data.substr(index + 1);
    stringstream s;
    string word;
    s << values;
    while(s >> word){
            wordVecs.push_back(dict->getVec(word));
    }
    vec = Vec();
    for (Vec& v : wordVecs) {
            vec += v;
    }
    vec /= (float)wordVecs.size();
}

Vec Abstract::getVec()
{
    return vec;
}

string Abstract::toString()
{
    return key + '|' + vec.toString();
}

string Abstract::getPmid(){
    return this->pmid;
}