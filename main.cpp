
#include<iostream>
#include<dirent.h>
#include<vector>
#include<unordered_map>
#include<fstream>
#include<regex>
#include<algorithm> 
#include<functional> 
#include<cctype>
#include<locale>
#include<sstream>
#include"omp.h"
#include"Abstract.h"
#include"Dict.h"
#include"Vec.h"
#include"Canonicalizer.h"

#include <flann/flann.hpp>
using namespace std;

//RELEASE
//const string HOME_DIR = "/scratch2/jsybran";

//DEBUG
const string HOME_DIR = "/home/jsybran/Projects/Data";

const string MEDLINE_XML_DIR = HOME_DIR + "/medline"; 
const string RESULTS_DIR = HOME_DIR + "/results";
const string CANON_FILE = RESULTS_DIR + "/canon";
const string VECTOR_FILE = RESULTS_DIR + "/canon.vec";
const string OUTPUT_FILE = RESULTS_DIR + "/graph.edges";
const string LOG_FILE = HOME_DIR + "/log.txt";


const string LVG_COMMAND = "lvg -f:0:C:P:q0:q1:q2:rs:g:T:t:u | awk 'BEGIN{FS=\"|\"}{print $2}'";
const string FASTTEXT_COMMAND = "fasttext skipgram -dim 500 -ws 8 -maxn 8 -thread 12 -input " + CANON_FILE + " -output " + CANON_FILE;
const string ABSTRACT_REGEX = "\\</?AbstractText.*?\\>";
const string PMID_REGEX = "\\</?PMID.*?\\>";
const int VECTOR_SIZE = 500;

const int CANON_BASH_SIZE = 50;

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


vector<string> getFilesInDir(string dirPath){
    vector<string> res;
    DIR *dir;
    dirent *ent;
    if((dir = opendir(dirPath.c_str())) != NULL){
        while((ent = readdir(dir)) != NULL){
            string name = ent->d_name;
            if(name != "." && name != "..")
                res.push_back(dirPath + "/" + ent->d_name);
        }
        closedir(dir);
    }else{
        cerr << "FAILED TO OPEN " << dirPath << endl;
    }
    return res;
}

//Returns PMID, AbstractText pairs
map<string,string> parseXML(string fileName){
    
    map<string,string> pmid2abstract;
    regex abstractRegex(ABSTRACT_REGEX,regex_constants::ECMAScript);
    regex pmidRegex(PMID_REGEX,regex_constants::ECMAScript);
    fstream fin(fileName, ios::in);
    string line;
    string lastFoundPMID = "NULL";
    string lastFoundAbstract = "";
    
    while(getline(fin,line)){
        if(regex_search(line,pmidRegex)){
            lastFoundAbstract = trim(lastFoundAbstract);
            if(lastFoundPMID != "NULL" && lastFoundAbstract != ""){
                pmid2abstract[lastFoundPMID] = lastFoundAbstract;
            }
            lastFoundPMID = regex_replace(line,pmidRegex,"");
            lastFoundAbstract = "";
        }
        if(regex_search(line,abstractRegex)){
            lastFoundAbstract += regex_replace(line,abstractRegex,"") + " ";
        }
    }
    //Need extra check at end
    if(lastFoundPMID != "NULL" && lastFoundAbstract != "")
        pmid2abstract[lastFoundPMID] = lastFoundAbstract;
    
    fin.close();
    return pmid2abstract;
}

void outputBash(vector<string> & tmpPmids, string output,fstream& fout){
    Canonicalizer canonMaker;
    stringstream r;
    r << canonMaker.getCanon(output);

    for(string pmid: tmpPmids){
        string t;
        getline(r,t);
        fout << pmid << " " << t <<endl;
    }
    tmpPmids.clear();
}

void parseMedline(unordered_map<string,string>& pmid2abstract, string dirPath){
    vector<string> xmlPaths = getFilesInDir(dirPath);
#pragma omp parallel  for
    for(int i = 0 ; i < xmlPaths.size();i++){
        
        map<string,string> tmp = parseXML(xmlPaths[i]);
                
        stringstream s;
        int BASH_SIZE = 50;
        
        vector<string> tmpPmids;
        fstream resOut(RESULTS_DIR+"/res"+to_string(i),ios::out);
    
        for(auto&val : tmp){
            s << val.second << endl;
            tmpPmids.push_back(val.first);
            
            if(tmpPmids.size() >= BASH_SIZE){
                outputBash(tmpPmids,s.str(),resOut);

                s.str( std::string() );
                s.clear();
                tmpPmids.clear();
            }
        }
        outputBash(tmpPmids,s.str(),resOut);
        resOut.close();
        
#pragma omp critical (INSERT_ABSTRACT)
        pmid2abstract.insert(tmp.begin(),tmp.end());
    }
}

void runFlann(unordered_map<string,Vec>& pmid2vec, vector<string>& pmids, string outputFilePath){
    //row major order, map.size() rows and VECTOR_SIZE cols
    float * vecData = new float[pmid2vec.size()*VECTOR_SIZE];
    
    int pmidCount = 0;
    for(string pmid : pmids){
        for(int i = 0 ; i < VECTOR_SIZE; i++){
            vecData[VECTOR_SIZE * pmidCount + i] = pmid2vec[pmid].get(i);
        }
        pmidCount++;
    }
        
    flann::Matrix<float> data(vecData,pmid2vec.size(),VECTOR_SIZE);
      
    flann::Index<flann::L2<float> > index(data, flann::KDTreeIndexParams(16));
    index.buildIndex();
    
    vector<vector<int> > indicies;
    vector<vector<float> > dists;
    
    flann::SearchParams params(128);
    params.cores = 0; //automatic core selection
    index.knnSearch(data, indicies, dists, 10, params);
 
    fstream fout(outputFilePath,ios::out);
    
    for(int currIndex = 0; currIndex < indicies.size(); ++currIndex){
        string pmid = pmids[currIndex];
        for(int nIndex = 0; nIndex < indicies[currIndex].size(); ++nIndex){
            int neighborID = indicies[currIndex][nIndex];
            string nPMID = pmids[neighborID];
            
            float pairDist = dists[currIndex][nIndex];
            
            fout<<pmid << " " << nPMID << " " <<pairDist<< endl;
        }
    }
    fout.close();
    delete[] vecData;
}

int main(int argc, char** argv) {

    fstream lout(LOG_FILE,ios::out);
    lout<<"Started"<<endl;
    
    unordered_map<string,string> pmid2abstract;
    
    lout<<"Parsing MEDLINE XML"<<endl;
    parseMedline(pmid2abstract,MEDLINE_XML_DIR);
    
    lout<<"Found " << pmid2abstract.size() << " abstracts"<<endl;
    
    vector<string> pmids;
    for(auto val : pmid2abstract){
        pmids.push_back(val.first);
    }
    
    fstream canonOut(CANON_FILE,ios::out);
    
    for(auto val : pmid2abstract){
        canonOut << val.second<<endl;
    }
    
    canonOut.close();
    
    lout<<"Training Fast Text"<<endl;
    lout<<"Running " << FASTTEXT_COMMAND << endl;
    system(FASTTEXT_COMMAND.c_str());
    
    lout<<"Building Dict"<<endl;
    Dict dict(VECTOR_FILE);
    
    unordered_map<string,Vec> pmid2vec;
    
    lout<<"Getting vectors per abstract"<<endl;
#pragma omp parallel for
    for(int i = 0 ; i < pmids.size();i++){
        string pmid = pmids[i];
        string canon = pmid2abstract[pmid];
        vector<Vec> wordVecs;
        stringstream s;
        string word;
        s << canon;
        while(s >> word){
            if(dict.contains(word))
                wordVecs.push_back(dict.getVec(word));
        }
        Vec vec;
        for (Vec& v : wordVecs) {
                vec += v;
        }
        if(wordVecs.size() > 0)
            vec /= (float)wordVecs.size();
#pragma omp critical (SET_PMID_VEC)
        pmid2vec[pmid] = vec;
    }
    
    lout<<"Running FLANN"<<endl;
    runFlann(pmid2vec,pmids, OUTPUT_FILE);
    
    lout << "DONE!" << endl;

    lout.close();
    return 0;
}

