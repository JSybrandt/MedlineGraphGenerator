/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Canonicalizer.cpp
 * Author: jsybran
 * 
 * Created on October 20, 2016, 8:23 AM
 */

#include "Canonicalizer.h"

using namespace std;

#define READ_END 0
#define WRITE_END 1

Canonicalizer::Canonicalizer() {
    
}

Canonicalizer::Canonicalizer(const Canonicalizer& orig) {
}

Canonicalizer::~Canonicalizer() {
}

string Canonicalizer::getCanon(string input){
    startProcess();
    char buffer[64];
    int bufferSize = sizeof(buffer);
    
    while(input.size()>0){
        
        if(input.size() > bufferSize){
            write(lvgPipe[WRITE_END],input.substr(0,bufferSize).c_str(),bufferSize);
            input = input.substr(bufferSize);
        }else{
            write(lvgPipe[WRITE_END],input.c_str(),input.size());
            input = "";
        }
    }
    close(lvgPipe[WRITE_END]);
    
    int bytesRead = 0;
    string res = "";
    while( (bytesRead = read(returnPipe[READ_END],buffer,sizeof(buffer)-1)) > 0){
        buffer[bytesRead] = '\0';
        res += buffer;
    }
    close(returnPipe[READ_END]);
    
    return res;
}

void Canonicalizer::startProcess(){
#pragma omp critical(START_PROC)
{
    pipe(lvgPipe);
    pipe(returnPipe);
    lvgPID = fork();
    if(lvgPID < 0){
        cerr << "Fork Failed"<<endl;
        exit(1);
    }
    
    //child
    if(lvgPID == 0){
               
        close(lvgPipe[WRITE_END]);
        dup2(lvgPipe[READ_END],STDIN_FILENO);
        
        close(returnPipe[READ_END]);
        dup2(returnPipe[WRITE_END],STDOUT_FILENO);
        
        execlp(LVG_COMMAND.c_str(),NULL);
    }    
    //parent
    else{
        close(lvgPipe[READ_END]);
        close(returnPipe[WRITE_END]);
    }
}
}