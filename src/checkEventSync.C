//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//Local dependencies
#include "include/mntToXRootdFileString.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

int checkEventSync(const std::string inFileName)
{
  std::cout << "TESTING: " << mntToXRootdFileString(inFileName) << std::endl;
  TFile* inFile_p = TFile::Open(mntToXRootdFileString(inFileName).c_str(), "READ");

  if(inFile_p->IsZombie()){
    std::cout << "Open failed. return 1" << std::endl;
    return 1;
  }
  else std::cout << "Successfully opened file \'" << mntToXRootdFileString(inFileName) << "\'" << std::endl;

  std::vector<std::string> fileTrees = returnRootFileContentsList(inFile_p, "TTree");

  //Remove duplicates
  unsigned int pos = 0;
  while(pos < fileTrees.size()){
    bool isGood = true;
    for(unsigned int i = pos+1; i < fileTrees.size(); ++i){
      if(isStrSame(fileTrees.at(i), fileTrees.at(pos))){
	fileTrees.erase(fileTrees.begin()+i);
	isGood = false;
	break;
      }
    }

    if(isGood) ++pos;
  }

  std::cout << "File \'" << inFileName << "\' contains..." << std::endl;
  for(unsigned int i = 0; i < fileTrees.size(); ++i){
    std::cout << " TTree " << i << "/" << fileTrees.size() << ": " << fileTrees.at(i) << std::endl;
    TTree* tempTree_p = (TTree*)inFile_p->Get(fileTrees.at(i).c_str());
    std::cout << "   Entries: " << tempTree_p->GetEntries() << std::endl;
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: checkEventSync.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkEventSync(argv[1]);
  return retVal;
}
