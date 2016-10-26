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
    lvgPID = 0;
}

Canonicalizer::Canonicalizer(const Canonicalizer& orig) {
}

Canonicalizer::~Canonicalizer() {
    endProcess();
}

string Canonicalizer::getCanon(string input){
    bool success = false;
    string res = "";
    int numFails = 0;
    do{
        if(startProcess()){
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
            while( (bytesRead = read(returnPipe[READ_END],buffer,sizeof(buffer)-1)) > 0){
                buffer[bytesRead] = '\0';
                res += buffer;
            }
            close(returnPipe[READ_END]);
            success = true;
        }
        else{
            numFails++;
        }
        endProcess();
    }while(!success && numFails < 10);
    
    return res;
}

bool Canonicalizer::startProcess(){
    bool ret = true;
#pragma omp critical(PROC)
{
    endProcess();
    if(pipe(lvgPipe) || pipe(returnPipe))
    {
        cerr << "PIPE FAILED" << endl;
        ret = false;
    }
    if((lvgPID = fork()) < 0){
        cerr << "Fork Failed"<<endl;
        ret = false;
    }
}
    //child
    if(lvgPID == 0){
               
        close(lvgPipe[WRITE_END]);
        dup2(lvgPipe[READ_END],STDIN_FILENO);
        
        close(returnPipe[READ_END]);
        dup2(returnPipe[WRITE_END],STDOUT_FILENO);
        
        execlp(LVG_COMMAND.c_str(),NULL);
        exit(1); //this is only reached if the exec fails
    }    
    //parent
    else{
        close(lvgPipe[READ_END]);
        close(returnPipe[WRITE_END]);
    }

    return ret;
}

bool Canonicalizer::endProcess(){
        if(lvgPID > 0){
            #pragma omp critical(PROC)
            {
                kill(lvgPID, SIGTERM);
            }
            lvgPID = 0;
            return true;
        }
        return false;
}
