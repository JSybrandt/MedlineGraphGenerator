
#include <stdio.h>
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
#include <sys/wait.h>
#include"omp.h"
#include"Abstract.h"
#include"Dict.h"
#include"Vec.h"
#include"Canonicalizer.h"
#include "constants.h"

#include <flann/flann.hpp>
using namespace std;



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
//second param are all the values we already have
map<string,string> parseXML(string fileName, const unordered_map<string,string>& backup){

    map<string,string> pmid2abstract;
    regex abstractRegex(ABSTRACT_REGEX,regex_constants::ECMAScript);
    regex pmidRegex(PMID_REGEX,regex_constants::ECMAScript);
    regex titleRegex(TITLE_REGEX,regex_constants::ECMAScript);
    regex eorRegex(END_OF_RECORD_REGEX,regex_constants::ECMAScript);
    fstream fin(fileName, ios::in);
    string line;
    string lastFoundPMID = "NULL";
    string lastFoundAbstract = "";

    while(getline(fin,line)){
        if(regex_search(line,pmidRegex)){
            lastFoundPMID = regex_replace(line,pmidRegex,"");
        }
        if(regex_search(line,abstractRegex)){
            lastFoundAbstract += regex_replace(line,abstractRegex,"") + " ";
        }
        if(regex_search(line,titleRegex)){
            lastFoundAbstract += regex_replace(line,titleRegex,"") + " ";
        }
        if(regex_search(line,eorRegex)){
            if(backup.find(lastFoundPMID) == backup.end() && lastFoundPMID != "NULL" && lastFoundAbstract != ""){
                pmid2abstract[lastFoundPMID] = trim(lastFoundAbstract);
            }
            lastFoundPMID = "NULL";
            lastFoundAbstract = "";
        }
    }

    fin.close();
    return pmid2abstract;
}

void outputBash(const vector<string> & tmpPmids, string output,fstream& fout, map<string,string>& tmpPmid2Abstract){
    Canonicalizer canonMaker;
    stringstream r;
    r << canonMaker.getCanon(output);

    for(string pmid: tmpPmids){
        string t;
        getline(r,t);
        fout << pmid << " " << t <<endl;
        tmpPmid2Abstract[pmid] = t;
    }
}

void parseMedline(unordered_map<string,string>& pmid2abstract, string dirPath, fstream& lout){
    vector<string> xmlPaths = getFilesInDir(dirPath);
    int completeCount = 0;
#pragma omp parallel  for
    for(int i = 0 ; i < xmlPaths.size();i++){

#pragma omp critical (LOGGER)
      lout << "Parsing:" << xmlPaths[i] << endl;

      // Simple PMID, abstract map, contains new values
      map<string,string> tmp = parseXML(xmlPaths[i],pmid2abstract);

      if(tmp.size() > 1){
#pragma omp critical (LOGGER)
         lout << "New:" << xmlPaths[i]<<endl
              <<"\tFound:" << xmlPaths.size()<<endl;

         stringstream s;
         int BASH_SIZE = 25;

         vector<string> tmpPmids;
         fstream resOut(RES_FILES_DIR+"/res"+to_string(i),ios::out);

         for(auto&val : tmp){
             s << val.second << endl;
             tmpPmids.push_back(val.first);

             if(tmpPmids.size() >= BASH_SIZE){
                 outputBash(tmpPmids,s.str(),resOut,tmp);
                 s.str( std::string() );
                 s.clear();
                 tmpPmids.clear();
             }
         }
         outputBash(tmpPmids,s.str(),resOut,tmp);
         tmpPmids.clear();
         resOut.close();

#pragma omp critical (INSERT_ABSTRACT)
         pmid2abstract.insert(tmp.begin(),tmp.end());
     } else{
#pragma omp critical (LOGGER)
       lout<<"Skip:" << xmlPaths[i]<<endl;
     }


     #pragma omp critical (LOGGER)
     {
         completeCount++;
         lout << (completeCount / (double) xmlPaths.size()) * 100 << "%" << endl;
     }
  }
}

void runFlann(unordered_map<string,Vec>& pmid2vec, vector<string>& pmids, string outputFilePath){
    //row major order, map.size() rows and VECTOR_SIZE cols
    float * vecData = new float[pmid2vec.size()*VECTOR_SIZE];

    cout << "VEC DATA SIZE" << (pmid2vec.size() * VECTOR_SIZE) << endl;
    cout << "EXPECTED DATA SIZE" << pmids.size() * VECTOR_SIZE << endl;

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

void loadOldCanon(unordered_map<string,string>& pmid2canon){
    vector<string> xmlPaths = getFilesInDir(BACKUP_FILES_DIR);

    #pragma omp parallel  for
    for(int i = 0 ; i < xmlPaths.size();i++){
        unordered_map<string,string> tmpMap;
        fstream res(xmlPaths[i].c_str(), ios::in);
        string line;
        while(getline(res,line)){
            stringstream s;
            s << line;
            string pmid;
            s >> pmid;
            string data = s.str();
            tmpMap[pmid] = data;
        }
#pragma omp critical (LOAD_CANON_LOCK)
        pmid2canon.insert(tmpMap.begin(),tmpMap.end());

        res.close();
    }
}

//This function is going to take the information stored in CANON_FILE
//and use ToPMine to generate the TOPMINE_OUT_FILE. This second file is a comma
//seperated list of keywords. We will convert this back into a version FastText can read,
//creating CANON_POT_TOPMINE_FILE. 
void runToPMine(const vector<string>& pmids, fstream& lout){
  lout << "Running ToPMine"<<endl;
  //actaully run topmine
  FILE* stdOut = popen(TOPMINE_COMMAND.c_str(), "r");
  int bufSize = 20;
  char buffer[bufSize];
  while(fgets(buffer,bufSize, stdOut)){
    for(char* p = buffer; *p; p++){
      lout << *p;
    }
  }
  fclose(stdOut);

  lout << "Converting TOPMINE_OUT to readable canon file" << endl;

  //at this point, data has moved from CANON_FILE to TOPMINE_OUT_FILE
  fstream topmineIn(TOPMINE_OUT_FILE.c_str(), ios::in);
  fstream postCanonOut(CANON_POST_TOPMINE_FILE.c_str(), ios::out);

  string line;
  string phrase;
  string word;
  bool firstWordInPhrase = true;

  while(getline(topmineIn,line)){
    stringstream lineStream;
    lineStream << line;
    while(getline(lineStream, phrase, ',')){
      stringstream phraseStream;
      phraseStream << phrase;
      firstWordInPhrase = true;
      while(phraseStream >> word){
        if(!firstWordInPhrase){
          postCanonOut<<"_";
        }
        firstWordInPhrase = false;
        postCanonOut << word;
      }
      postCanonOut << " ";
    }
    postCanonOut << endl;
  }

  topmineIn.close();
  postCanonOut.close();

  //canon is not in CANON_POST_TOPMINE_FILE
  lout << "Finished ToPMine. Converted canon now in " << CANON_POST_TOPMINE_FILE << endl;
}

void makeResultsDirs(){
  system(("mkdir -p " + MEDLINE_XML_DIR).c_str());
  system(("mkdir -p " + RESULTS_DIR).c_str());
  system(("mkdir -p " + RES_FILES_DIR).c_str());
  system(("mkdir -p " + BACKUP_FILES_DIR).c_str());
  system(("mkdir -p " + ABSTRACT_VECTOR_DIR).c_str());
}

void catchChild(int sigNum){
    /* when we get here, we know there's a zombie child waiting */
    int child_status;
    wait(&child_status);
}

int main(int argc, char** argv) {

  const int EXPECTED_SIZE = 30000000;
    signal(SIGCHLD, catchChild);

    fstream lout(LOG_FILE,ios::out);
    lout<<"Started"<<endl;

    lout<<"Making Directories"<<endl;
    makeResultsDirs();

    lout<<"Initializing maps" <<endl;

    unordered_map<string,string> pmid2abstract(EXPECTED_SIZE);
    unordered_map<string,Vec> pmid2vec(EXPECTED_SIZE);
    vector<string> pmids(EXPECTED_SIZE);

    //if abstract - vec file doesn't exist
    if(!ifstream(LOAD_ABSTRACT_VECTOR_FILE.c_str())){

        //if canon file doesn't exist
        if(!ifstream(CANON_FILE.c_str())){
            lout<<"Recovering Lost Data"<<endl;
            loadOldCanon(pmid2abstract);

            lout<<"Recovered "<< pmid2abstract.size() << " old records" << endl;
            lout<<"Parsing MEDLINE XML"<<endl;

            parseMedline(pmid2abstract,MEDLINE_XML_DIR,lout);

            lout<<"Found " << pmid2abstract.size() << " total abstracts"<<endl;

            //save canon in order
            fstream canonOut(CANON_FILE,ios::out);
            fstream pmidOut(CANON_PMID_ORDER_FILE,ios::out);

            for(auto val : pmid2abstract){
                canonOut << val.second<<endl;
                pmidOut << val.first << endl;
                pmids.push_back(val.first);
            }

            canonOut.close();
            pmidOut.close();
        }else{
          lout << "Skipping canonicalizing"<< endl;
          /*
          lout << "Recovering PMIDS"<< endl;
          fstream pmidIn(CANON_PMID_ORDER_FILE,ios::in);
          string tmp;
          while(pmidIn >> tmp)
              pmids.push_back(tmp);
          pmidIn.close();
          lout << "Recovered " << pmids.size() << " pmids"<<endl;
          */
        }

        //If we have not sent the canon through TOPMINE
        if(!ifstream(CANON_POST_TOPMINE_FILE.c_str())){
          runToPMine(pmids, lout);
        }else{
          pmids.clear();
          lout<< "skiping topmine"<<endl;
          long long index = 1;
          fstream canonIn(CANON_POST_TOPMINE_FILE.c_str(),ios::in);
          string line;
          while(getline(canonIn, line)){
            string id = "i" + to_string(index);
            pmids.push_back(id);
            pmid2abstract[id] = line;
            index++;
          }
          canonIn.close();
          lout << "Loaded canon post topmine"<<endl;
          lout << "Loaded " << pmid2abstract.size() << " abstracts" << endl;
        }
        //if vector file doesnt exist (this is word - vec file)
        if(!ifstream(VECTOR_FILE.c_str())){
            lout<<"Training Fast Text"<<endl;
            lout<<"Running " << FASTTEXT_COMMAND << endl;
            system(FASTTEXT_COMMAND.c_str());
        }else{
            lout<<"Skipping training with fasttext"<<endl;
        }

        lout<<"Building Dict"<<endl;
        Dict dict(VECTOR_FILE);
        lout << "Loaded " << dict.size() << " words"<<endl;

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
            if(wordVecs.size() > 0) pmid2vec[pmid] = vec;
        }

        lout << "Saving Vecs" << endl;
#pragma omp parallel
        {
            string outFile = ABSTRACT_VECTOR_DIR + "/absVec" + to_string(omp_get_thread_num());
            ofstream vecsFile(outFile.c_str());

#pragma omp for
            for(int i = 0 ; i < pmids.size(); i++){
                string pmid = pmids[i];
                if(pmid2vec.find(pmid)!=pmid2vec.end())
                  vecsFile << pmid << " " << pmid2vec[pmid].toString() << endl;
            }

            vecsFile.close();
        }
        fstream lvfOut(LOAD_ABSTRACT_VECTOR_FILE.c_str(),ios::out);
        lvfOut << "SAVED " << pmids.size() << " vectors"<<endl;
        lvfOut.close();

    } else {

        lout << "Loading vecs" << endl;

        //loading saved vecs
        vector<string> backupFiles = getFilesInDir(ABSTRACT_VECTOR_DIR);
#pragma omp parallel for
        for(int i = 0 ; i < backupFiles.size(); i++){
            string path = backupFiles[i];
            unordered_map<string,Vec> tempMap;
            ifstream pmidVecFile(path.c_str());
            string line;
            while(getline(pmidVecFile,line)){
                string pmid;
                float tmp;
                vector<float> vecData;
                stringstream s;
                s << line;
                s >> pmid;
                while(s >> tmp) vecData.push_back(tmp);
                tempMap[pmid] = Vec(vecData);
            }
            pmidVecFile.close();

#pragma omp critical (VEC_BACKUP)
            pmid2vec.insert(tempMap.begin(), tempMap.end());
        }

        for(auto val : pmid2vec){
            pmids.push_back(val.first);
        }

        lout << "Loaded " << pmid2vec.size() << " pmid vectors"<<endl;
    }

    if(!ifstream(OUTPUT_FILE.c_str())){
        lout<<"Running FLANN"<<endl;
        pmids.clear();
        for(auto pair : pmid2vec)
          pmids.push_back(pair.first);
        runFlann(pmid2vec,pmids, OUTPUT_FILE);
    }else{
        lout << "SKIPPING FLANN" << endl;
    }
    lout << "DONE!" << endl;

    lout.close();
    return 0;
}

