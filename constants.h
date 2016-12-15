#pragma once

#include <string>
using std::string;

//RELEASE
const string HOME_DIR = "/scratch2/jsybran";
const std::string LVG_COMMAND = "/home/jsybran/projects/MedlineGraphGenerator/bin/lvgCommand.sh";

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
const std::string ABSTRACT_VECTOR_DIR = RESULTS_DIR + "/abstractVecs";
const std::string LOAD_ABSTRACT_VECTOR_FILE = RESULTS_DIR + "/loadAbstractVecs";
const std::string TOPMINE_OUT_FILE = CANON_FILE + "_topic_out";
const std::string CANON_POST_TOPMINE_FILE = RESULTS_DIR + "/canon_post_topmine";
const std::string CANON_PMID_ORDER_FILE = RESULTS_DIR + "/canon_pmids";

//const string LVG_COMMAND = "lvg -f:0:C:P:q0:q1:q2:rs:g:T:t:u | awk 'BEGIN{FS=\"|\"}{print $2}'";
const std::string FASTTEXT_COMMAND = "fasttext skipgram -minCount 0 -minCountLabel 0 -dim 500 -ws 8 -maxn 8 -thread 24 -input " + CANON_POST_TOPMINE_FILE + " -output " + CANON_FILE;
const std::string TOPMINE_COMMAND = "topmine";
const std::string ABSTRACT_REGEX = "\\</?AbstractText.*?\\>";
const std::string PMID_REGEX = "\\</?PMID.*?\\>";
const std::string TITLE_REGEX = "\\</?ArticleTitle.*?\\>";
const std::string END_OF_RECORD_REGEX = "\\</MedlineCitation.*?\\>";
const int VECTOR_SIZE = 500;

const int CANON_BASH_SIZE = 50;


