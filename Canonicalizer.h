#pragma once

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <signal.h>
#include <unistd.h>
#include "constants.h"

using std::string;


class Canonicalizer {
public:
    Canonicalizer();
    Canonicalizer(const Canonicalizer& orig);
    virtual ~Canonicalizer();
    
    string getCanon(string input);
   
private:
    int lvgPipe[2];
    int returnPipe[2];
    pid_t lvgPID;
    
    bool startProcess();
    bool endProcess();
};


