///AUTHOR: CFMCGINN (2018.04.12)
///For validating forest
//cpp dependencies
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TDatime.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLeaf.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/cppWatch.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/mntToXRootdFileString.h"
#include "include/removeVectorDuplicates.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

bool getDoLog(double minVal, double maxVal, int orderMag)
{
  if(minVal<=0) return false;
  int order = 0;
  while(minVal < maxVal){
    minVal *= 10.;
    order++;
  }

  return order > orderMag;
}

void dumpTreeNames(std::string& fileName, std::vector<std::string>& treeNames)
{
  std::cout << "TTrees In File \'" << fileName << "\':\n";
  for(unsigned int tI = 0; tI < treeNames.size(); ++tI){
    std::cout << " Tree " << tI << "/" << treeNames.size() << ": " << treeNames[tI] << "\n";
  }

  return;
}

bool doStringsMatch(std::vector<std::string>& fileList1, std::vector<std::string>& fileList2, std::vector<std::string>& misMatch1, std::vector<std::string>& misMatch2)
{
  bool allMatch = true;
  unsigned int fI1 = 0;
  while(fI1 < fileList1.size()){
    bool matchFound = false;
    for(unsigned int fI2 = 0; fI2 < fileList2.size(); ++fI2){
      if(fileList1[fI1].find(fileList2[fI2]) != std::string::npos && fileList1[fI1].size() == fileList2[fI2].size()){
	matchFound = true;
	break;
      }
    }
    if(!matchFound){
      allMatch = false;
      std::cout << " No match found for Tree \'" << fileList1[fI1] << "\', removing...\n";
      misMatch2.push_back(fileList1[fI1]);
      fileList1.erase(fileList1.begin()+fI1);
    }
    else fI1++;
  }

  fI1 = 0;
  while(fI1 < fileList2.size()){
    bool matchFound = false;
    for(unsigned int fI2 = 0; fI2 < fileList1.size(); ++fI2){
      if(fileList2[fI1].find(fileList1[fI2]) != std::string::npos && fileList2[fI1].size() == fileList1[fI2].size()){
	matchFound = true;
	break;
      }
    }
    if(!matchFound){
      allMatch = false;
      std::cout << " No match found for Tree \'" << fileList2[fI1] << "\', removing...\n";
      misMatch1.push_back(fileList2[fI1]);
      fileList2.erase(fileList2.begin()+fI1);
    }
    else fI1++;
  }

  return allMatch;
}

std::string texFriendlyString(std::string& inStr)
{
  std::string outStr = "";
  while(inStr.find("_") != std::string::npos){
    outStr = outStr + inStr.substr(0, inStr.find("_")) + "\\_";
    inStr.replace(0, inStr.find("_")+1, "");
  }
  outStr = outStr + inStr;

  return outStr;
}

template<typename T>
void minmax_vector_branch(TTree* tree_p[], int const nFiles,
			  std::vector<std::vector<std::string>>& branchList, uint32_t bI1,
			  double& minVal, double& maxVal) {
  std::vector<T>* vector_p = 0;

  for(int fI = 0; fI < nFiles; ++fI){
    tree_p[fI]->SetBranchStatus("*", 0);
    tree_p[fI]->SetBranchStatus(branchList[0][bI1].data(), 1);
    tree_p[fI]->SetBranchAddress(branchList[0][bI1].data(), &vector_p);

    const int nEntries1 = tree_p[fI]->GetEntries();
    for(int entry = 0; entry < nEntries1; ++entry){
      tree_p[fI]->GetEntry(entry);

      if(vector_p->empty()) continue;

      auto maxMin = std::minmax_element(vector_p->begin(), vector_p->end());
      if(minVal > *(maxMin.first))
	minVal = *(maxMin.first);
      if(maxVal < *(maxMin.second))
	maxVal = *(maxMin.second);
    }
  }
}

void doFirstTexSlide(std::ofstream& fileTex, std::vector<std::string>& inFileNames, std::vector<std::string>& inNickNames, std::vector<std::string>& goodTrees, std::vector<std::vector<std::string>>& missingTrees, const int eventCountOverride)
{
  TDatime date;

  for(unsigned int fI = 0; fI < inFileNames.size(); ++fI){
    inFileNames[fI] = texFriendlyString(inFileNames[fI]);
    inFileNames[fI] = inFileNames[fI].substr(0, inFileNames[fI].size()/2) + "\\\\" + inFileNames[fI].substr(inFileNames[fI].size()/2, inFileNames[fI].size());
  }

  fileTex << "\\RequirePackage{xspace}\n";
  fileTex << "\\RequirePackage{amsmath}\n";

  fileTex << "\n";

  fileTex << "\\documentclass[xcolor=dvipsnames]{beamer}\n";
  fileTex << "\\usetheme{Warsaw}\n";
  fileTex << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}\n";
  fileTex << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}\n";
  fileTex << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}\n";
  fileTex << "\\setbeamersize{text margin left=5pt,text margin right=5pt}\n";

  fileTex << "\n";

  fileTex << "\\setbeamertemplate{frametitle}\n";
  fileTex << "{\n";
  fileTex << "  \\nointerlineskip\n";
  fileTex << "  \\begin{beamercolorbox}[sep=0.3cm, ht=1.8em, wd=\\paperwidth]{frametitle}\n";
  fileTex << "    \\vbox{}\\vskip-2ex%\n";
  fileTex << "    \\strut\\insertframetitle\\strut\n";
  fileTex << "    \\vskip-0.8ex%\n";
  fileTex << "  \\end{beamercolorbox}\n";
  fileTex << "}\n";

  fileTex << "\n";

  fileTex << "\\setbeamertemplate{footline}{%\n";
  fileTex << "  \\begin{beamercolorbox}[sep=.8em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}\n";
  fileTex << "    \\hspace{0.3cm}%\n";
  fileTex << "    \\hfill\\insertauthor \\hfill\\insertpagenumber\n";
  fileTex << "  \\end{beamercolorbox}%\n";
  fileTex << "}\n";

  fileTex << "\n";

  fileTex << "\\setbeamertemplate{navigation symbols}{}\n";
  fileTex << "\\setbeamertemplate{itemize item}[circle]\n";
  fileTex << "\\setbeamertemplate{itemize subitem}[circle]\n";
  fileTex << "\\setbeamertemplate{itemize subsubitem}[circle]\n";
  fileTex << "\\setbeamercolor{itemize item}{fg=black}\n";
  fileTex << "\\setbeamercolor{itemize subitem}{fg=black}\n";
  fileTex << "\\setbeamercolor{itemize subsubitem}{fg=black}\n";

  fileTex << "\n";

  fileTex << "\\definecolor{links}{HTML}{00BFFF}\n";
  fileTex << "\\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=links}\n";

  fileTex << "\n";

  fileTex << "\\author[CM]{CMSHI Forest Validator}\n";
  fileTex << "\\begin{document}\n";

  fileTex << "\n";


  fileTex << "\\begin{frame}\n";
  fileTex << "\\frametitle{\\centerline{\\hypertarget{TopLevel}{Forest Validation} (" << date.GetYear() << "." << date.GetMonth() << "." << date.GetDay() << ")}}\n";
  fileTex << " \\begin{itemize}\n";
  fileTex << "  \\fontsize{5}{5}\\selectfont\n";
  if(eventCountOverride >= 0) fileTex << " \\item{WARNING: EVENT NUMBER OVERRIDE SET, TEST INVALID}\n";

  for(unsigned int fI = 0; fI < inFileNames.size(); ++fI){
    fileTex << "  \\item{File " << fI+1 << ": \'" << inFileNames[fI] << "\'}\n";
    fileTex << "  \\item{Nickname: \'" << texFriendlyString(inNickNames[fI]) << "\'}\n";
  }

  fileTex << "  \\item{Good Trees:}\n";
  fileTex << " \\begin{itemize}\n";
  fileTex << "  \\fontsize{5}{5}\\selectfont\n";
  if(goodTrees.size() != 0){
    for(unsigned int tI = 0; tI < goodTrees.size(); ++tI){
      fileTex << "  \\item{Tree " << tI << ": \'" << texFriendlyString(goodTrees[tI]) << "\',      \\hyperlink{tocSummary_" << goodTrees[tI] << "}{To Summary TOC},      \\hyperlink{tocFull_0_" << goodTrees[tI] << "}{To All Branch TOC}}\n";
    }
  }
  else fileTex << "  \\item{WARNING: NO GOOD TREES}\n";
  fileTex << " \\end{itemize}\n";

  for(unsigned int fI = 0; fI < inNickNames.size(); ++fI){
    fileTex << "  \\item{Bad File " << fI+1 << ", " << texFriendlyString(inNickNames[fI]) << ", Trees:}\n";

    fileTex << " \\begin{itemize}\n";
    fileTex << "  \\fontsize{5}{5}\\selectfont\n";
    if(missingTrees[fI].size() != 0){
      for(unsigned int tI = 0; tI < missingTrees[fI].size(); ++tI){
	fileTex << "  \\item{Tree " << tI << ": \'" << texFriendlyString(missingTrees[fI][tI]) << "\'}\n";
      }
    }
    else fileTex << "  \\item{NO MISSING TREES}\n";
    fileTex << " \\end{itemize}\n";
  }

  fileTex << " \\end{itemize}\n";
  fileTex << "\\end{frame}\n";
  fileTex << "\n";

  return;
}


void doInterstitialTexSlide(std::ofstream& fileTex, const std::string& interstitialString)
{
  fileTex << "\\begin{frame}\n";
  fileTex << "\\fontsize{7}{7}\\selectfont\n";
  fileTex << "\\frametitle{\\centerline{Interstitial}}\n";
  fileTex << "\\begin{itemize}\n";
  fileTex << "\\fontsize{7}{7}\\selectfont\n";
  fileTex << "\\item{" << interstitialString << "}\n";
  fileTex << "\\end{itemize}\n";
  fileTex << "\\end{frame}\n";
  fileTex << "\n";

  return;
}


void doTreeTexSlide(std::ofstream& fileTex, std::string& inTreeName, std::vector<std::string>& nickNames, std::vector<std::vector<std::string>>& listOfBadVar, std::vector<std::string>& warningList, std::vector<int>& warningListPos)
{
  std::string treeTitle = texFriendlyString(inTreeName);
  while(treeTitle.find("/") != std::string::npos){
    treeTitle.replace(treeTitle.find("/"), 1, " ");
  }

  std::vector<std::string> bigColumnVect;

  for(unsigned int fI = 0; fI < listOfBadVar.size(); ++fI){
    bigColumnVect.push_back("\\fontsize{7}{7}\\selectfont \\bf{Missing variables from \'" + texFriendlyString(nickNames[fI]) + "\' file:}");

    if(!listOfBadVar[fI].empty()){
      for(unsigned int vI = 0; vI < listOfBadVar[fI].size(); ++vI){
	std::string varStr = std::to_string(vI) + "/" + std::to_string(listOfBadVar[fI].size()) + ": \'" + texFriendlyString(listOfBadVar[fI][vI]) + "\'";
	bigColumnVect.push_back(varStr);
      }
    }
    else bigColumnVect.push_back("\\bf{\\qquad No missing variables!}");
  }

  bigColumnVect.push_back("\\fontsize{7}{7}\\selectfont \\bf{Variables non-unity in ratio:}");
  if(warningList.size() == 0) bigColumnVect.push_back("\\bf{\\qquad No variables non-unity!}");
  else{
    for(unsigned int bI = 0; bI < warningList.size(); ++bI){
      std::string warnStr = "Var " + std::to_string(bI) + "/" + std::to_string(warningList.size()) + ": \'\\hyperlink{branch_" + inTreeName + "_" + std::to_string(warningListPos[bI]) + "}{" + texFriendlyString(warningList[bI]) + "}\'";
      bigColumnVect.push_back(warnStr);
    }
  }

  for(unsigned int bI = 0; bI < bigColumnVect.size(); ++bI){
    if(bI%48 == 0){
      if(bI != 0){
	fileTex << "\\end{itemize}\n";
	fileTex << "\\end{column}\n";
	fileTex << "\\end{columns}\n";
	fileTex << "\\end{itemize}\n";
	fileTex << "\\end{frame}\n";
	fileTex << "\n";
      }

      fileTex << "\n";
      fileTex << "\\begin{frame}\n";
      fileTex << "\\fontsize{7}{7}\\selectfont\n";
      if(bI == 0) fileTex << "\\frametitle{\\centerline{\\hypertarget{tocSummary_" << inTreeName << "}{" << treeTitle << "}}}\n";
      else fileTex << "\\frametitle{\\centerline{" << treeTitle << "}}\n";
      fileTex << "\\begin{itemize}\n";
      fileTex << "\\fontsize{7}{7}\\selectfont\n";
      fileTex << "\\item{Tree Name \'" << texFriendlyString(inTreeName) << "\'}, (\\hyperlink{Back to Top}{TopLevel})\n";

      fileTex << "\\begin{columns}\n";
      fileTex << "\\fontsize{5}{5}\\selectfont\n";

      fileTex << "\\begin{column}{0.32\\textwidth}\n";
      fileTex << "\\fontsize{5}{5}\\selectfont\n";
      fileTex << "\\begin{itemize}\n";
      fileTex << "\\fontsize{5}{5}\\selectfont\n";
    }
    else if(bI%16 == 0){
      fileTex << "\\end{itemize}\n";
      fileTex << "\\end{column}\n";
      fileTex << "\n";

      fileTex << "\\begin{column}{0.32\\textwidth}\n";
      fileTex << "\\fontsize{5}{5}\\selectfont\n";
      fileTex << "\\begin{itemize}\n";
      fileTex << "\\fontsize{5}{5}\\selectfont\n";
    }

    fileTex << "\\item{" << bigColumnVect[bI] << "}\n";
  }

  fileTex << "\\end{itemize}\n";
  fileTex << "\\end{column}\n";

  if(bigColumnVect.size()%48 < 32){
    fileTex << "\\begin{column}{0.32\\textwidth}\n";
    fileTex << "\\fontsize{5}{5}\\selectfont\n";
    fileTex << "\\begin{itemize}\n";
    fileTex << "\\fontsize{5}{5}\\selectfont\n";
    fileTex << "\\item{Dummy}\n";
    fileTex << "\\end{itemize}\n";
    fileTex << "\\end{column}\n";
  }

  if(bigColumnVect.size()%48 < 16){
    fileTex << "\\begin{column}{0.32\\textwidth}\n";
    fileTex << "\\fontsize{5}{5}\\selectfont\n";
    fileTex << "\\begin{itemize}\n";
    fileTex << "\\fontsize{5}{5}\\selectfont\n";
    fileTex << "\\item{Dummy}\n";
    fileTex << "\\end{itemize}\n";
    fileTex << "\\end{column}\n";
  }

  fileTex << "\\end{columns}\n";
  fileTex << "\\end{itemize}\n";
  fileTex << "\\end{frame}\n";

  fileTex << "\n";
  fileTex << "\n";

  return;
}


void doBranchTexSlide(std::ofstream& fileTex, std::string& inTreeName, std::vector<std::string>& listOfBranches)
{
  std::string treeTitle = texFriendlyString(inTreeName);
  while(treeTitle.find("/") != std::string::npos){
    treeTitle.replace(treeTitle.find("/"), 1, " ");
  }

  int nBranches = 0;
  std::vector<std::vector<std::string> > branchVect1, branchVect2, branchVect3;

  for(unsigned int bI = 0; bI < listOfBranches.size(); ++bI){
    if(bI%48 == 0){
      branchVect1.push_back({});
      branchVect2.push_back({});
      branchVect3.push_back({});
    }

    if((bI/16)%3 == 0) branchVect1[branchVect1.size()-1].push_back(listOfBranches[bI]);
    else if((bI/16)%3 == 1) branchVect2[branchVect2.size()-1].push_back(listOfBranches[bI]);
    else if((bI/16)%3 == 2) branchVect3[branchVect3.size()-1].push_back(listOfBranches[bI]);
  }

  for(unsigned int tempI = 0; tempI < branchVect1.size(); ++tempI){
    fileTex << "\\begin{frame}\n";
    fileTex << "\\fontsize{4}{4}\\selectfont\n";
    fileTex << "\\frametitle{\\centerline{\\hypertarget{tocFull_" << tempI << "_" << inTreeName << "}{" << treeTitle << "}}}\n";
    fileTex << "\\begin{itemize}\n";
    fileTex << "\\fontsize{4}{4}\\selectfont\n";
    fileTex << "\\item{Tree Name \'" << texFriendlyString(inTreeName) << "\', (\\hyperlink{Back to Top}{TopLevel})}\n";
    fileTex << "\\item{Branches}\n";

    fileTex << "\\begin{columns}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";

    fileTex << "\\begin{column}{0.32\\textwidth}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";
    fileTex << "\\begin{itemize}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";
    if(!branchVect1[tempI].empty()){
      for(unsigned int vI = 0; vI < branchVect1[tempI].size(); ++vI){
	fileTex << "\\item{Var " << nBranches << "/" << listOfBranches.size() << ": \'\\hyperlink{branch_" << inTreeName << "_" << nBranches << "}{" << texFriendlyString(branchVect1[tempI][vI]) << "}\'}\n";
	nBranches++;
      }
    }
    else fileTex << "\\item{Dummy}\n";

    fileTex << "\\end{itemize}\n";
    fileTex << "\\end{column}\n";

    fileTex << "\\begin{column}{0.32\\textwidth}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";
    fileTex << "\\begin{itemize}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";
    if(!branchVect2[tempI].empty()){
      for(unsigned int vI = 0; vI < branchVect2[tempI].size(); ++vI){
	fileTex << "\\item{Var " << nBranches << "/" << listOfBranches.size() << ": \'\\hyperlink{branch_" << inTreeName << "_" << nBranches << "}{" << texFriendlyString(branchVect2[tempI][vI]) << "}\'}\n";
	nBranches++;
      }
    }
    else fileTex << "\\item{Dummy}\n";

    fileTex << "\\end{itemize}\n";
    fileTex << "\\end{column}\n";

    fileTex << "\\begin{column}{0.32\\textwidth}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";
    fileTex << "\\begin{itemize}\n";
    fileTex << "\\fontsize{1}{1}\\selectfont\n";

    if(!branchVect3[tempI].empty()){
      for(unsigned int vI = 0; vI < branchVect3[tempI].size(); ++vI){
	fileTex << "\\item{Var " << nBranches << "/" << listOfBranches.size() << ": \'\\hyperlink{branch_" << inTreeName << "_" << nBranches << "}{" << texFriendlyString(branchVect3[tempI][vI]) << "}\'}\n";
	nBranches++;
      }
    }
    else fileTex << "\\item{Dummy}\n";

    fileTex << "\\end{itemize}\n";
    fileTex << "\\end{column}\n";

    fileTex << "\\end{columns}\n";

    fileTex << "\\end{itemize}\n";
    fileTex << "\\end{frame}\n";

    fileTex << "\n";
  }

  return;
}


void doPlotTexSlide(std::ofstream& fileTex, std::string& inTreeName, std::string& inBranchName, int inBranchNum, std::string& inPdfName)
{
  std::string treeTitle = texFriendlyString(inTreeName);
  while(treeTitle.find("/") != std::string::npos){
    treeTitle.replace(treeTitle.find("/"), 1, " ");
  }

  fileTex << "\\begin{frame}\n";
  fileTex << "\\fontsize{7}{7}\\selectfont\n";
  fileTex << "\\frametitle{\\centerline{\\hypertarget{branch_" << inTreeName << "_" << inBranchNum << "}{" << texFriendlyString(inBranchName) << "} (" << treeTitle << ")}}\n";

  fileTex << "\\begin{center}\n";
  fileTex << "\\includegraphics[width=0.6\\textwidth]{" << inPdfName << "}\n";
  fileTex << "\\end{center}\n";

  fileTex << "\\begin{itemize}\n";
  fileTex << "\\fontsize{7}{7}\\selectfont\n";
  fileTex << "\\item{" << texFriendlyString(inBranchName) << ", \\hyperlink{tocSummary_" << inTreeName << "}{To Summary TOC}, \\hyperlink{tocFull_0_" << inTreeName << "}{To All Branch TOC}}\n";
  fileTex << "\\end{itemize}\n";
  fileTex << "\\end{frame}\n";

  fileTex << "\n";

  return;
}


int runForestDQM(std::vector<std::string>& inFileNames, const std::string& additionalNickName, std::vector<std::string>& inNickNames, const std::string& treeSelect = "", const bool doEventNorm = false, const int eventCountOverride = -1, const std::string& commaSeparatedPairList = "")
{
  std::string globalYStrTemp = "Counts";
  if(doEventNorm) globalYStrTemp = "#frac{1}{N_{evt}} Counts";

  const std::string globalYStr = globalYStrTemp;

  const int nFiles = inFileNames.size();
  const int fileCap = 8;
  if(nFiles > fileCap){
    std::cout << "Number of files given, " << nFiles << ", is greater than file cap, " << fileCap << ", for this macro. return 1\n";
    return 1;
  }

  for(int fI = 0; fI < nFiles; ++fI){
    if(inFileNames[fI].find(".root") == std::string::npos || !checkFile(inFileNames[fI])){
      std::cout << "Input \'" << inFileNames[fI] << "\' is invalid return 1.\n";
      return 1;
    }
    else if(inNickNames[fI].size() == 0) inNickNames[fI] = "File" + std::to_string(fI+1);
  }

  std::vector<std::string> pairVars1;
  std::vector<std::string> pairVars2;

  std::vector<double> pairVarsMin1;
  std::vector<double> pairVarsMin2;
  std::vector<double> pairVarsMax1;
  std::vector<double> pairVarsMax2;

  if(!commaSeparatedPairList.empty()){
    std::string commaSeparatedPairListCopy = commaSeparatedPairList;

    if(commaSeparatedPairListCopy.substr(commaSeparatedPairListCopy.size()-1, 1).find(",") == std::string::npos)
      commaSeparatedPairListCopy = commaSeparatedPairListCopy + ",";

    while(commaSeparatedPairListCopy.find(",") != std::string::npos){
      std::string varStr = commaSeparatedPairListCopy.substr(0, commaSeparatedPairListCopy.find(","));
      if(varStr.find(":") == std::string::npos || varStr.find(":") != varStr.rfind(":")){
	std::cout << "Potential variable pair \'" << varStr << "\' from \'" << commaSeparatedPairList << "\' is not valid syntax. Please input as \'var1:var2\'. skipping...\n";
      }
      else{
	std::string var1 = varStr.substr(0, varStr.find(":"));
	varStr.replace(0, varStr.find(":")+1, "");

	pairVars1.push_back(var1);
	pairVars2.push_back(varStr);

	pairVarsMin1.push_back(-1);
	pairVarsMax1.push_back(-1);
	pairVarsMin2.push_back(-1);
	pairVarsMax2.push_back(-1);
      }

      commaSeparatedPairListCopy.replace(0, commaSeparatedPairListCopy.find(",")+1, "");
    }
  }

  std::cout << "Processing the following potential pairs in all TTree: \n";
  for(unsigned int i = 0; i < pairVars1.size(); ++i){
    std::cout << " " << i << "/" << pairVars1.size() << ": " << pairVars1[i] << ":" << pairVars2[i] << "\n";
  }

  Double_t minDeltaVal = 0.000000001;
  kirchnerPalette col;
  TDatime date;
  const std::string dateStr = std::to_string(date.GetDate());
  const int nBins = 50;

  const int colors[fileCap] = {1, col.getColor(2), col.getColor(0), col.getColor(3), col.getColor(2), col.getColor(0), col.getColor(3), col.getColor(2)};
  const int styles[fileCap] = {20, 25, 28, 27, 46, 24, 42, 44};

  checkMakeDir("pdfDir");
  checkMakeDir(("pdfDir/" + dateStr).data());

  TLine line;
  line.SetLineStyle(2);

  TH1F* dummyHists_p[nFiles];
  for(int fI = 0; fI < nFiles; ++fI){
    dummyHists_p[fI] = new TH1F();
    dummyHists_p[fI]->SetMarkerStyle(styles[fI]);
    dummyHists_p[fI]->SetMarkerSize(1.0);
    dummyHists_p[fI]->SetMarkerColor(colors[fI]);
    dummyHists_p[fI]->SetLineColor(colors[fI]);
  }

  TLegend* leg_p = new TLegend(.6, .6, .9, .9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(12);

  for(int fI = 0; fI < nFiles; ++fI)
    leg_p->AddEntry(dummyHists_p[fI], inNickNames[fI].data(), "P L");

  std::vector<std::vector<std::string > > misMatchedTrees;

  std::string texFileName = "forestDQM_" + additionalNickName;
  for(int fI = 0; fI < nFiles; ++fI){
    texFileName = texFileName + "_" + inNickNames[fI];
  }
  texFileName = texFileName + "_" + dateStr;

  std::string texTreeStr = "AllTrees";
  if(treeSelect.size() != 0)
    texTreeStr = treeSelect;
  while(texTreeStr.find("/") != std::string::npos)
    texTreeStr.replace(texTreeStr.find("/"), 1, "_");

  std::string fullDirName = "pdfDir/" + dateStr + "/" + texFileName;
  while(fullDirName.substr(fullDirName.size()-1,1).find("_") != std::string::npos)
    fullDirName.replace(fullDirName.size()-1, 1, "");
  fullDirName = fullDirName + "/";

  std::cout << "File name 1: " << texFileName << "\n";
  checkMakeDir(fullDirName.data());

  std::string texFileSubName = fullDirName + texFileName + "VALIDATION_SUB_" + texTreeStr + "_" + dateStr + ".tex";
  std::cout << "File name 2: " << texFileSubName << "\n";
  texFileName = fullDirName + texFileName + "VALIDATION_MAIN_" + texTreeStr + "_" + dateStr + ".tex";
  std::ofstream fileTex(texFileName.data());
  std::ofstream fileSubTex(texFileSubName.data());
  doInterstitialTexSlide(fileSubTex, "PLOTS");

  TFile* inFiles_p[nFiles];
  std::vector<std::vector<std::string > > fileTrees;
  for(int fI = 0; fI < nFiles; ++fI){
    inFiles_p[fI] = TFile::Open(mntToXRootdFileString(inFileNames[fI]).data(), "READ");
    std::vector<std::string> tempTrees = returnRootFileContentsList(inFiles_p[fI], "TTree", treeSelect);
    removeVectorDuplicates(tempTrees);

    /* remove hltanalysis/hltobject trees */
    tempTrees.erase(std::remove(tempTrees.begin(), tempTrees.end(), "hltanalysis/HltTree"), tempTrees.end());
    tempTrees.erase(std::remove_if(tempTrees.begin(), tempTrees.end(), [](const std::string& tree) {
	return tree.find("hltobject") != std::string::npos; }), tempTrees.end());

    dumpTreeNames(inFileNames[fI], tempTrees);
    fileTrees.push_back(tempTrees);
    misMatchedTrees.push_back({});
  }

  std::cout << "Checking tree strings match...\n";

  for(int fI = 0; fI < nFiles-1; ++fI){
    for(int fI2 = fI+1; fI2 < nFiles; ++fI2){
      if(doStringsMatch(fileTrees[fI], fileTrees[fI2], misMatchedTrees[fI], misMatchedTrees[fI2]))
	std::cout << "All trees match between files!\n";
      else
	std::cout << "Mismatched trees found and removed...\n";
    }
  }

  doFirstTexSlide(fileTex, inFileNames, inNickNames, fileTrees[0], misMatchedTrees, eventCountOverride);

  for(unsigned int tI = 0; tI < fileTrees[0].size(); ++tI){
    std::cout << "Processing \'" << fileTrees[0][tI] << "\'...\n";

    TTree* tree_p[nFiles];
    TObjArray* tempBranchList[nFiles];
    TObjArray* tempLeafList[nFiles];
    std::vector<std::vector<std::string > > branchList;
    std::vector<std::vector<std::string > > leafList;
    std::vector<std::string > pdfList;
    std::vector<std::string> warningList;
    std::vector<int> warningListPos;

    for(int fI = 0; fI < nFiles; ++fI){
      tree_p[fI] = (TTree*)inFiles_p[fI]->Get(fileTrees[0][tI].data());
      tempBranchList[fI] = (TObjArray*)(tree_p[fI]->GetListOfBranches()->Clone());
      tempLeafList[fI] = (TObjArray*)(tree_p[fI]->GetListOfLeaves()->Clone());

      std::vector<std::string> uncleanedBranches;
      std::vector<std::string> uncleanedLeaves;

      for(int tempI = 0; tempI < tempBranchList[fI]->GetEntries(); ++tempI){
	uncleanedBranches.push_back(tempBranchList[fI]->At(tempI)->GetName());
      }
      for(int tempI = 0; tempI < tempLeafList[fI]->GetEntries(); ++tempI){
	uncleanedLeaves.push_back(tempLeafList[fI]->At(tempI)->GetName());
      }

      unsigned int tempPos = 0;
      while(tempPos < uncleanedBranches.size()){
	TBranch* tempBranch = (TBranch*)tree_p[fI]->GetBranch(uncleanedBranches[tempPos].data());

	if(tempBranch->GetListOfBranches()->GetEntries() != 0){
	  std::cout << "WARNING: Branch \'" << tempBranch->GetName()
		    << "\' in TTree \'" << fileTrees[0][tI]
		    << "\' has nSubBranch = " << tempBranch->GetListOfBranches()->GetEntries()
		    << " not equal to 0. runForestDQM.exe cannot handle this branch, will be removed\n";

	  for(int rbI = 0; rbI < tempBranch->GetListOfBranches()->GetEntries()+1; ++rbI){
	    std::cout << "Removing \'" << uncleanedLeaves[tempPos] << "\' from leaflist at " << tempPos << "\n";
	    uncleanedLeaves.erase(uncleanedLeaves.begin()+tempPos);
	  }
	  uncleanedBranches.erase(uncleanedBranches.begin()+tempPos);
	}
	else tempPos++;
      }

      if(uncleanedBranches.size() != uncleanedLeaves.size()){
	std::cout << "WARNING: branchList in TTree \'" << fileTrees[0][tI]
		  << "\' has nBranches = " << uncleanedBranches.size()
		  << " after cleaning, differing from nLeaves = " << uncleanedLeaves.size() << ".\n";
      }

      branchList.push_back({});
      leafList.push_back({});

      for(unsigned int bI1 = 0; bI1 < uncleanedBranches.size(); ++bI1){
	std::string branchName = uncleanedBranches[bI1];
	std::string leafName = uncleanedLeaves[bI1];

	if(!isStrSame(branchName, leafName)){
	  std::cout << "WARNING: Branch \'" << branchName << "\' differs from matched leaf \'" << leafName << "\'.\n";
	}

	branchList[fI].push_back(branchName);
	leafList[fI].push_back(leafName);
      }
    }

    for(unsigned int bI = 0; bI < branchList[0].size(); ++bI){
      std::cout << " " << bI << "/" << branchList[0].size() << ": " << branchList[0][bI] << "\n";
    }

    std::vector<std::string> misMatchedNickNames;
    std::vector<std::vector<std::string > > misMatchedBranches;
    std::vector<std::vector<std::string > > misMatchedLeaves;

    for(int fI = 0; fI < nFiles; ++fI){
      misMatchedBranches.push_back({});
      misMatchedLeaves.push_back({});
    }

    std::cout << "Checking branch strings match...\n";

    for(int fI = 0; fI < nFiles-1; ++fI){
      for(int fI2 = fI+1; fI2 < nFiles; ++fI2){
	if(doStringsMatch(branchList[fI], branchList[fI2], misMatchedBranches[fI], misMatchedBranches[fI2]))
	  std::cout << "All branches match between files!\n";
	else
	  std::cout << "Mismatched branches found and removed...\n";

	if(doStringsMatch(leafList[fI], leafList[fI2], misMatchedLeaves[fI], misMatchedLeaves[fI2]))
	  std::cout << "All leaves match between files!\n";
	else
	  std::cout << "Mismatched leaves found and removed...\n";
      }
    }

    std::vector<bool> pairVars1Found;
    std::vector<bool> pairVars2Found;

    for(unsigned int i = 0; i < pairVars1.size(); ++i){
      pairVars1Found.push_back(false);
      pairVars2Found.push_back(false);
    }

    for(unsigned int bI1 = 0; bI1 < branchList[0].size(); ++bI1){
      std::cout << " Processing \'" << branchList[0][bI1] << "\'...\n";

      tree_p[0]->ResetBranchAddresses();
      tree_p[0]->SetBranchStatus("*", 0);
      tree_p[0]->SetBranchStatus(branchList[0][bI1].data(), 1);

      Double_t maxVal = tree_p[0]->GetMaximum(branchList[0][bI1].data());
      Double_t minVal = tree_p[0]->GetMinimum(branchList[0][bI1].data());

      std::string maxValFile = inFileNames[0];
      std::string minValFile = inFileNames[0];

      for(int fI = 0; fI < nFiles; ++fI){
	tree_p[fI]->ResetBranchAddresses();
	tree_p[fI]->SetBranchStatus("*", 0);
	tree_p[fI]->SetBranchStatus(branchList[0][bI1].data(), 1);

	Double_t tempMaxVal = tree_p[fI]->GetMaximum(branchList[0][bI1].data());
	Double_t tempMinVal = tree_p[fI]->GetMinimum(branchList[0][bI1].data());

	if(tempMaxVal > maxVal){
	  maxVal = tempMaxVal;
	  maxValFile = inFileNames[fI];
	}

	if(tempMinVal < minVal){
	  minVal = tempMinVal;
	  minValFile = inFileNames[fI];
	}

	//	maxVal = TMath::Max(maxVal, tree_p[fI]->GetMaximum(branchList[0][bI1].data()));
	//	minVal = TMath::Min(minVal, tree_p[fI]->GetMinimum(branchList[0][bI1].data()));
      }

      TLeaf* tempLeaf = (TLeaf*)tree_p[0]->GetLeaf(leafList[0][bI1].data());
      std::string tempClassType = tempLeaf->GetTypeName();

      //We have to handle vectors differently, getMax and getMin do not work from ttree
      if(tempClassType.find("vector") != std::string::npos && tempClassType.find("vector<vector") == std::string::npos){
	if(tempClassType.find("int") != std::string::npos){
	  minmax_vector_branch<int>(tree_p, nFiles, branchList, bI1, minVal, maxVal);
	}
	else if(tempClassType.find("short") != std::string::npos){
	  minmax_vector_branch<short>(tree_p, nFiles, branchList, bI1, minVal, maxVal);
	}
	else if(tempClassType.find("float") != std::string::npos){
	  minmax_vector_branch<float>(tree_p, nFiles, branchList, bI1, minVal, maxVal);
	}
	else if(tempClassType.find("double") != std::string::npos){
	  minmax_vector_branch<double>(tree_p, nFiles, branchList, bI1, minVal, maxVal);
	}
	else{
	  std::cout << "Warning do not know how to handle vector \'" << tempClassType << "\'. return 1\n";
	  return 1;
	}
      }

      if (!isfinite(minVal) || !isfinite(maxVal)) {
	warningList.push_back(branchList[0][bI1]);
	warningListPos.push_back(bI1);
      }

      /* usually indicates an empty array */
      if (minVal == std::numeric_limits<float>::max() &&
	  maxVal == -std::numeric_limits<float>::max()) {
	minVal = 0; maxVal = 1;
      }

      std::vector<std::string> histNames;
      for(int fI = 0; fI < nFiles; ++fI){
	std::string histName = fileTrees[0][tI] + "_" + branchList[0][bI1] + "_" + inNickNames[fI] + "_h";
	while(histName.find("/") != std::string::npos){histName.replace(histName.find("/"), 1, "_");}
	histNames.push_back(histName);
      }

      TCanvas* canv_p = new TCanvas("temp", "temp", 450, 400);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);
      canv_p->SetBottomMargin(0.01);
      canv_p->SetLeftMargin(0.01);
      TPad* pads_p[2];
      pads_p[0] = new TPad("pads0", "", 0.0, 0.35, 1.0, 1.0);
      canv_p->cd();
      pads_p[0]->Draw("SAME");
      pads_p[0]->SetTopMargin(0.01);
      pads_p[0]->SetRightMargin(pads_p[0]->GetLeftMargin());
      pads_p[0]->SetBottomMargin(0.0);

      pads_p[1] = new TPad("pads1", "", 0.0, 0.00, 1.0, 0.35);
      canv_p->cd();
      pads_p[1]->Draw("SAME");
      pads_p[1]->SetTopMargin(0.0);
      pads_p[1]->SetRightMargin(pads_p[0]->GetLeftMargin());
      pads_p[1]->SetBottomMargin(pads_p[1]->GetLeftMargin()*3.);

      if(TMath::Abs(minVal - maxVal) < 0.000000001){
	minVal -= 1;
	maxVal +=1;
      }

      TH1D* tempHist_p[nFiles];

      Double_t bins[nBins+1];

      for(unsigned int pI = 0; pI < pairVars1.size(); ++pI){
	if(isStrSame(branchList[0][bI1], pairVars1[pI])){
	  pairVars1Found[pI] = true;
	  pairVarsMax1[pI] = maxVal;
	  pairVarsMin1[pI] = minVal;
	}
	else if(isStrSame(branchList[0][bI1], pairVars2[pI])){
	  pairVars2Found[pI] = true;
	  pairVarsMax2[pI] = maxVal;
	  pairVarsMin2[pI] = minVal;
	}
      }

      bool doLogX = getDoLog(minVal, maxVal, 2);

      if(doLogX){
	maxVal *= 2;
	minVal /= 2;
	getLogBins(minVal, maxVal, nBins, bins);
      }
      else{
	Double_t interval = maxVal - minVal;
	maxVal += interval/10.;
	minVal -= interval/10.;
	getLinBins(minVal, maxVal, nBins, bins);
      }

      for(int fI = 0; fI < nFiles; ++fI){
	tree_p[fI]->ResetBranchAddresses();
	tree_p[fI]->SetBranchStatus("*", 0);
	tree_p[fI]->SetBranchStatus(branchList[0][bI1].data(), 1);

	tempHist_p[fI] = new TH1D(histNames[fI].data(), (";" + branchList[0][bI1] + ";" + globalYStr).data(), nBins, bins);

	if(eventCountOverride < 0) tree_p[fI]->Project(histNames[fI].data(), branchList[0][bI1].data(), "", "");
	else tree_p[fI]->Project(histNames[fI].data(), branchList[0][bI1].data(), "", "", eventCountOverride);

	tempHist_p[fI]->GetXaxis()->SetTitle(branchList[0][bI1].data());
	tempHist_p[fI]->GetYaxis()->SetTitle(globalYStr.data());

	tempHist_p[fI]->SetMarkerSize(1.0);
	tempHist_p[fI]->SetMarkerStyle(styles[fI]);
	tempHist_p[fI]->SetMarkerColor(colors[fI]);
	tempHist_p[fI]->SetLineColor(colors[fI]);

	tempHist_p[fI]->GetXaxis()->SetTitleFont(43);
	tempHist_p[fI]->GetYaxis()->SetTitleFont(43);
	tempHist_p[fI]->GetXaxis()->SetLabelFont(43);
	tempHist_p[fI]->GetYaxis()->SetLabelFont(43);

	tempHist_p[fI]->GetXaxis()->SetTitleSize(12);
	tempHist_p[fI]->GetYaxis()->SetTitleSize(10);
	tempHist_p[fI]->GetXaxis()->SetLabelSize(12);
	tempHist_p[fI]->GetYaxis()->SetLabelSize(10);

	tempHist_p[fI]->GetXaxis()->SetTitleOffset(tempHist_p[fI]->GetXaxis()->GetTitleOffset()*4.);
	tempHist_p[fI]->GetYaxis()->SetTitleOffset(tempHist_p[fI]->GetXaxis()->GetTitleOffset()/3.);

	centerTitles(tempHist_p[fI]);
	setSumW2(tempHist_p[fI]);

	if(doEventNorm) tempHist_p[fI]->Scale(1./(Double_t)tree_p[fI]->GetEntries());
      }

      maxVal = tempHist_p[0]->GetMinimum();
      minVal = tempHist_p[0]->GetMaximum();

      if (minVal == 0 && maxVal == 0) { maxVal = 1; }

      for(int fI = 0; fI < nFiles; ++fI){
	for(int bIX = 0; bIX < tempHist_p[fI]->GetNbinsX(); ++bIX){
	  if(tempHist_p[fI]->GetBinContent(bIX+1) != 0){
	    if(tempHist_p[fI]->GetBinContent(bIX+1) > maxVal) maxVal = tempHist_p[fI]->GetBinContent(bIX+1);
	    if(tempHist_p[fI]->GetBinContent(bIX+1) < minVal) minVal = tempHist_p[fI]->GetBinContent(bIX+1);
	  }
	}
      }

      Double_t interval = maxVal - minVal;

      if(interval > 1000 && minVal > 0){
	maxVal *= 10;
	minVal /= 10;
      }
      else{
	maxVal += interval/10.;
	if(minVal - interval/10 > 0 || minVal < 0) minVal -= interval/10;
	else minVal = 0;
      }

      tempHist_p[0]->SetMaximum(maxVal);
      tempHist_p[0]->SetMinimum(minVal);

      canv_p->cd();
      pads_p[0]->cd();

      for(int fI = 0; fI < nFiles; ++fI){
	if(fI == 0) tempHist_p[fI]->DrawCopy("E1 P");
	else tempHist_p[fI]->DrawCopy("E1 P SAME");
      }

      leg_p->Draw("SAME");

      if(getDoLog(minVal, maxVal, 2)) gPad->SetLogy();
      if(doLogX) gPad->SetLogx();
      gStyle->SetOptStat(0);
      gPad->RedrawAxis();
      gPad->SetTicks(0, 1);

      canv_p->cd();
      pads_p[1]->cd();

      for(int fI = 1; fI < nFiles; ++fI){
	tempHist_p[fI]->Divide(tempHist_p[0]);
      }

      bool isGood = true;

      for(int bIX = 0; bIX < tempHist_p[0]->GetNbinsX(); ++bIX){
	for(int fI = 1; fI < nFiles; ++fI){
	  if(TMath::Abs(tempHist_p[fI]->GetBinContent(bIX+1)) < minDeltaVal){
	    if(TMath::Abs(tempHist_p[0]->GetBinContent(bIX+1)) < minDeltaVal) continue;
	    else{
	      isGood = false;
	      break;
	    }
	  }

	  if(!isGood) break;

	  if(TMath::Abs(tempHist_p[fI]->GetBinContent(bIX+1) - 1.) > minDeltaVal){
	    isGood = false;
	    break;
	  }
	}
      }

      if(!isGood){
	warningList.push_back(branchList[0][bI1]);
	warningListPos.push_back(bI1);
      }

      for(int fI = 1; fI < nFiles; ++fI){
	tempHist_p[fI]->GetYaxis()->SetTitle(("All/" + inNickNames[0]).data());
	tempHist_p[fI]->SetMaximum(2.);
	tempHist_p[fI]->SetMinimum(0.);
	if(fI == 1) tempHist_p[fI]->DrawCopy("E1 P");
	else tempHist_p[fI]->DrawCopy("E1 P SAME");
      }

      if(doLogX) gPad->SetLogx();

      line.DrawLine(tempHist_p[0]->GetBinLowEdge(1), 1, tempHist_p[0]->GetBinLowEdge(201), 1);
      gPad->RedrawAxis();
      gPad->SetTicks(0, 1);

      std::string saveName = histNames[0] + "_" + additionalNickName;
      for(int fI = 1; fI < nFiles; ++fI){
	saveName = saveName + "_" + inNickNames[fI];
      }
      saveName = saveName + "_" + dateStr + ".pdf";
      pdfList.push_back(saveName);
      quietSaveAs(canv_p, fullDirName + saveName);

      for(int fI = 0; fI < nFiles; ++fI){
	delete tempHist_p[fI];
      }

      delete pads_p[0];
      delete pads_p[1];
      delete canv_p;
    }

    for(unsigned int pI = 0; pI < pairVars1Found.size(); ++pI){
      if(!pairVars1Found[pI]) continue;
      if(!pairVars2Found[pI]) continue;

      for(int hI = 0; hI < nFiles; ++hI){
	tree_p[hI]->ResetBranchAddresses();
	tree_p[hI]->SetBranchStatus("*", 0);
	tree_p[hI]->SetBranchStatus(pairVars1[pI].data(), 1);
	tree_p[hI]->SetBranchStatus(pairVars2[pI].data(), 1);

	std::string histName = fileTrees[0][tI] + "_" + pairVars1[pI] + "_" + pairVars2[pI] + "_" + inNickNames[hI] + "_h";
	while(histName.find("/") != std::string::npos){histName.replace(histName.find("/"), 1, "_");}
	TH2D* hist_p = new TH2D(histName.data(), (";" + pairVars2[pI] + ";"+ pairVars1[pI]).data(), nBins, pairVarsMin2[pI], pairVarsMax2[pI], nBins, pairVarsMin1[pI], pairVarsMax1[pI]);

	TCanvas* canv_p = new TCanvas("tempCanv_p", "tempCanv_p", 450, 450);
	canv_p->SetTopMargin(0.12);
	canv_p->SetBottomMargin(0.12);
	canv_p->SetLeftMargin(0.12);
	canv_p->SetRightMargin(0.12);

	if(eventCountOverride < 0) tree_p[hI]->Project(histName.data(), (pairVars1[pI] + ":" + pairVars2[pI]).data(), "", "");
	else tree_p[hI]->Project(histName.data(), (pairVars1[pI] + ":" + pairVars2[pI]).data(), "", "", eventCountOverride);

	hist_p->GetXaxis()->SetTitleFont(43);
	hist_p->GetYaxis()->SetTitleFont(43);
	hist_p->GetXaxis()->SetLabelFont(43);
	hist_p->GetYaxis()->SetLabelFont(43);

	hist_p->GetXaxis()->SetTitleSize(12);
	hist_p->GetYaxis()->SetTitleSize(12);
	hist_p->GetXaxis()->SetLabelSize(12);
	hist_p->GetYaxis()->SetLabelSize(12);

	centerTitles(hist_p);
	setSumW2(hist_p);

	if(doEventNorm) hist_p->Scale(1./(Double_t)tree_p[hI]->GetEntries());

	hist_p->DrawCopy("COLZ");
	gPad->SetLogz();

	gPad->RedrawAxis();
	gPad->SetTicks(0, 1);

	std::string saveName = histName + "_" + additionalNickName;
	saveName = saveName + "_" + dateStr + ".pdf";
	pdfList.push_back(saveName);
	quietSaveAs(canv_p, fullDirName + saveName);

	branchList[0].push_back(pairVars1[pI] + ":" + pairVars2[pI] + "-" + inNickNames[hI]);

	delete hist_p;
	delete canv_p;
      }
    }

    doTreeTexSlide(fileTex, fileTrees[0][tI], inNickNames, misMatchedBranches, warningList, warningListPos);
    doBranchTexSlide(fileSubTex, fileTrees[0][tI], branchList[0]);
    for(unsigned int bI1 = 0; bI1 < branchList[0].size(); ++bI1){
      doPlotTexSlide(fileSubTex, fileTrees[0][tI], branchList[0][bI1], bI1, pdfList[bI1]);
    }
  }

  delete leg_p;

  for(int fI = 0; fI < nFiles; ++fI){
    delete dummyHists_p[fI];
    inFiles_p[fI]->Close();
    delete inFiles_p[fI];
  }

  if(texFileSubName.find(fullDirName) != std::string::npos){
    texFileSubName.replace(texFileSubName.find(fullDirName), fullDirName.size(), "");
  }

  fileTex << "\\input{" << texFileSubName << "}\n";
  fileTex << "\\end{document}\n";

  return 0;
}

int main(int argc, char* argv[])
{
  cppWatch watch;
  watch.start();

  if(argc < 3 || argc > 8){
    std::cout << "Usage: ./bin/runForestDQM.exe <inFileNames> <additionalNickName> <inNickNames> <treeSelect> <doEventNorm> <eventCountOverride> <commaSeparatedPairList>\n";
    std::cout << " inFileNames, inNickNames a comma separated list of arbitrarily many files\n";
    return 1;
  }

  std::vector<std::string> inFiles;
  std::string argv1 = argv[1];
  while(argv1.find(",") != std::string::npos){
    inFiles.push_back(argv1.substr(0, argv1.find(",")));
    argv1.replace(0, argv1.find(",")+1, "");
  }
  if(argv1.size() != 0) inFiles.push_back(argv1);

  std::vector<std::string> inNickNames;

  if(argc >= 4){
    std::string argv3 = argv[3];
    while(argv3.find(",") != std::string::npos){
      if(argv3.substr(0, argv3.find(",")).size() == 0) inNickNames.push_back("");
      else inNickNames.push_back(argv3.substr(0, argv3.find(",")));
      argv3.replace(0, argv3.find(",")+1, "");
    }
    if(argv3.size() != 0) inNickNames.push_back(argv3);
  }
  else{
    for(unsigned int fI = 0; fI < inFiles.size(); ++fI){
      inNickNames.push_back("");
    }
  }

  if(inFiles.size() != inNickNames.size() && argc >= 3){
    std::cout << "Warning: number of files, " << inFiles.size()
	      << ", given by \'" << argv[1]
	      << "\', is not equal to number of nicknames, " << inNickNames.size()
	      << ", given by \'" << argv[2] << "\'. return 1\n";
    return 1;
  }

  int retVal = 0;
  if(argc == 3 || argc == 4) retVal += runForestDQM(inFiles, argv[2], inNickNames);
  else if(argc == 5) retVal += runForestDQM(inFiles, argv[2], inNickNames, argv[4]);
  else if(argc == 6) retVal += runForestDQM(inFiles, argv[2], inNickNames, argv[4], std::stoi(argv[5]));
  else if(argc == 7) retVal += runForestDQM(inFiles, argv[2], inNickNames, argv[4], std::stoi(argv[5]), std::stoi(argv[6]));
  else if(argc == 8) retVal += runForestDQM(inFiles, argv[2], inNickNames, argv[4], std::stoi(argv[5]), std::stoi(argv[6]), argv[7]);

  watch.stop();
  std::cout << "Timing: " << watch.total() << "\n";

  return retVal;
}
