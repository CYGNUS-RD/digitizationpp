/*
 * Copyright (C) 2024 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2024
 *
 */

#include <iostream>
#include <unistd.h>
#include <limits.h>
#include <sstream>
#include <fstream>
#include <chrono>
#include <map>
#include <random>
#include <string>
//#include <execution>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "s3.h"
#include "cygnolib.h"

using namespace std;

vector<string> split(const string& , char );
void ReadConfig(string name, map<string,string>& options);

int main(int argc, char** argv)
{
    
    if(argc<2) {cerr<<"No Configfile given!!\nSuggested use: ./progname.exe Configfile -I inputdir -O outputdir"; exit(EXIT_FAILURE);}
    string nome=argv[1];
    
    vector<string> path_to_config=split(nome,'/');
    string parte= "";
    for(unsigned int i=0;i<path_to_config.size()-2;i++) parte=parte+path_to_config[i]+ "/";
    char buffer[PATH_MAX];
    char* ret = getcwd(buffer,sizeof(buffer));
    if(!ret) {
        cerr << "Failed to get current directory... exiting" << endl;
        return EXIT_FAILURE;
    }
    string currentPath(buffer);
    const string SOURCE_DIR= currentPath+"/"+parte;
    
    map<string,string> options;
    ReadConfig(nome,options);
    
    
    string tmpfolder = options["bckg_path"];
    string tmpname   = options["bckg_name"];
    
    if(! filesystem::exists(tmpfolder)){
        //DEBUG
        cout<<"Creating tmpfolder..."<<
        system(("mkdir " + tmpfolder).c_str() );
    }
    
    if(! filesystem::exists(tmpfolder+tmpname)) {
        
        vector<int> pedruns;
        
        stringstream ssruns(options["noiserun"]);
        string run;
        vector<string> seglist;
        while(getline(ssruns, run, ';')) {
            int runi = stoi(run);
            pedruns.push_back(runi);
        }
        
        cygnolib::joinMidasPedestals(pedruns, options["ped_cloud_dir"], tmpfolder, tmpname);
        
    }
    
    return 0;
    
}

// ====== FUNCTIONS DEFINITION ======
vector<string> split(const string& str, char delimiter)
{
    vector<string> tokens;
    istringstream stream(str);
    string token;

    while(getline(stream,token,delimiter)) tokens.push_back(token);

    return tokens;
}
void ReadConfig(string name, map<string,string>& options)
{
    ifstream config(name.c_str());
    
    string line;
    while(getline(config,line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
        if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;
        
        auto delim1= line.find(":");
        auto delim2= line.find(",");
        if(delim2==string::npos) delim2=line.size();
        auto index= line.substr(0,delim1);
        auto val= line.substr(delim1+1,min(delim2,line.size())-delim1-1 );
        options[index]=val;
    }
    
}
