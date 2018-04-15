#ifndef REMOVEVECTORDUPLICATES_H
#define REMOVEVECTORDUPLICATES_H

#include <string>
#include <vector>

void removeVectorDuplicates(std::vector<std::string>* vect)
{
  unsigned int pos = 0;
  while(pos < vect->size()){
    bool isGood = true;
    for(unsigned int p = pos+1; p < vect->size(); ++p){
      if(vect->at(p).size() == vect->at(pos).size() && vect->at(pos).find(vect->at(p)) != std::string::npos){
        vect->erase(vect->begin() + p);
        isGood = false;
        break;
      }
    }
    if(isGood) ++pos;
  }

  return;
}

#endif
