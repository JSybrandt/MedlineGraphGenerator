#pragma once

#include <string>
using std::string;

//RELEASE
const string HOME_DIR = "/scratch2/jsybran";
const std::string LVG_COMMAND = "/scratch2/jsybran/lvgCommand.sh";

//DEBUG
//const std::string HOME_DIR = "/home/jsybran/Projects/Data";
//const std::string LVG_COMMAND = "./lvgCommand.sh";

const std::string MEDLINE_XML_DIR = HOME_DIR + "/medline"; 
const std::string RESULTS_DIR = HOME_DIR + "/results";
const std::string RES_FILES_DIR = RESULTS_DIR + "/res";
const std::string BACKUP_FILES_DIR = RESULTS_DIR + "/backup";
const std::string CANON_FILE = RESULTS_DIR + "/canon";
const std::string VECTOR_FILE = RESULTS_DIR + "/canon.vec";
const std::string OUTPUT_FILE = RESULTS_DIR + "/graph.edges";
const std::string LOG_FILE = HOME_DIR + "/log.txt";


//const string LVG_COMMAND = "lvg -f:0:C:P:q0:q1:q2:rs:g:T:t:u | awk 'BEGIN{FS=\"|\"}{print $2}'";
const std::string FASTTEXT_COMMAND = "fasttext skipgram -dim 500 -ws 8 -maxn 8 -thread 24 -input " + CANON_FILE + " -output " + CANON_FILE;
const std::string ABSTRACT_REGEX = "\\</?AbstractText.*?\\>";
const std::string PMID_REGEX = "\\</?PMID.*?\\>";
const int VECTOR_SIZE = 500;

const int CANON_BASH_SIZE = 50;


