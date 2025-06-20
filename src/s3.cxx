/*
 * Copyright (C) 2024 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2024
 *
 */

#include "s3.h"
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <curl/curl.h>

namespace s3 {
    std::string BUCKET_POSIX_PATH = "/jupyter-workspace/cloud-storage/"; ///< Bucket Posix path on the cloud
    std::string BUCKET_REST_PATH = "https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/"; ///< Bucket rest path on the cloud

    std::string mid_file(int run, std::string tag, bool cloud, bool verbose) {
        std::stringstream ss;
        ss << std::setw(5) << std::setfill('0') << run;
        std::string srun = ss.str();
        
        std::string BASE_URL;
        std::string f;
        
        if(cloud) {
            BASE_URL = BUCKET_REST_PATH+"cygno-data/";
            f = BASE_URL+tag+"/run"+srun+".mid.gz";
        } else {
            f = "/run"+srun+".mid.gz";
        }
        if (verbose) std::cout<<f<<std::endl;
        return f;
    }
    
    
    static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
      size_t written = std::fwrite(ptr, size, nmemb, (FILE *)stream);
      return written;
    }
    std::string cache_file(std::string fname,
                           std::string path,
                           bool        cloud,
                           std::string tag,
                           bool      verbose) {
        
        std::string return_name;
        
        if(!cloud) {
            return_name = path+tag+fname;
            return return_name;
        } else {
            if(path != "" && !std::filesystem::exists(path.c_str())) {
                int com = std::system(("mkdir "+path).c_str());
                if(com != 0) {
                    std::cout<<com<<std::endl;
                }
            }
            
            std::stringstream sstmp(fname);
            std::string substring;
            std::vector<std::string> subslist;
            while( std::getline(sstmp, substring, '/') ) {
               subslist.push_back(substring);
            }
            std::string tmpname = path + subslist[subslist.size()-1];
            
            
            if(!std::filesystem::exists(tmpname.c_str())) {
                if(verbose) std::cout<<"File "<<tmpname<<" not found. Downloading it from cloud..."<<std::endl;
                
                CURL *curl_handle;
                FILE *pagefile;
                curl_global_init(CURL_GLOBAL_ALL);

                // init the curl session
                curl_handle = curl_easy_init();
                // set URL to get here
                curl_easy_setopt(curl_handle, CURLOPT_URL, fname.c_str());

                // check if file exists remotely
                CURL *curl_err;
                curl_err    = curl_easy_init();
                curl_easy_setopt(curl_err, CURLOPT_URL, fname.c_str());
                curl_easy_setopt(curl_err, CURLOPT_NOBODY, 1L); // Use HEAD request
                curl_easy_setopt(curl_err, CURLOPT_FOLLOWLOCATION, 1L); // Follow redirects
                CURLcode iRc = curl_easy_perform(curl_err);
                long response_code = 0;
                if (iRc == CURLE_OK) {
                    curl_easy_getinfo(curl_err, CURLINFO_RESPONSE_CODE, &response_code);
                }
                if(iRc != CURLE_OK || response_code != 200) {
                    throw std::runtime_error("\nFile "+tmpname+" not found on cloud location "+fname+".");
                }

                // Switch on full protocol/debug output while testing
                if (verbose) curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);
                // disable progress meter, set to 0L to enable and disable debug output
                curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 0L);
                // send all data to this function
                curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
                // open the file
                pagefile = fopen(tmpname.c_str(), "wb");
                if(pagefile) {
                    // write the page body to this file handle
                    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);
                    // get it!
                    curl_easy_perform(curl_handle);
                    // close the header file
                    fclose(pagefile);
                }
                // cleanup curl stuff
                curl_easy_cleanup(curl_handle);
                curl_easy_cleanup(curl_err);
                curl_global_cleanup();

            }
            
            return_name = tmpname; //to change
        }
        return return_name;
    }
    
}
