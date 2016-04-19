//    The ASCII IO library for the PARENT program suite
//    Copyright (C) 2016  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License  version 3
//    as published by the Free Software Foundation.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.





//    A scientific publication about this program has been released in the Journal of Chemical Theory and Computation:
//		"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"
//		DOI: 10.1021/acs.jctc.5b01217

//   We kindly ask you to include this citation in works that publish
//   results generated using this program or any modifications of it.







#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>



#include "io_text.h"

using namespace std;

string strip_line(string line) {
    string tmpstring1="";
    string tmpstring2="";
    unsigned short int flag=0,length;

    //write the line without preceding blanks
    for(unsigned int i=0; i<line.length(); i++) {
        if(line[i]!=' ') {
            flag=1;
        }
        if(flag==1) {
            tmpstring1+=line[i];
        }
    }

    //write the remaining line without concluding blanks  (inverted)
    flag=0;
    for(int i=tmpstring1.length()-1; i>=0; i--) {
        if(tmpstring1[i]!=' ') {
            flag=1;
        }
        if(flag==1) {
            tmpstring2+=tmpstring1[i];
        }
    }

    //reinvert the resulting string to get the correct order
    tmpstring1=tmpstring2;
    length=tmpstring2.length();
    for(int i=0; i<length; i++) {
        tmpstring1[length-i-1]=tmpstring2[i];
    }
    return tmpstring1;
}

string delete_char(string line,char del) {//self-explenatory
    string tmpstring="";
    for(unsigned int i=0; i<line.length(); i++) {
        if(line[i]!=del) {
            tmpstring+=line[i];
        }
    }
    return tmpstring;
}

string strip_blanks(string line) {
    //write only non-blank characters

    string tmpstring1="";
    for(unsigned int i=0; i<line.length(); i++) {
        if(line[i]!=' ') {
            tmpstring1+=line[i];
        }
    }
    return tmpstring1;
}


//to read an index file
int read_ndx_file(ifstream *ndxfile, vector <int> *group1, vector <int> *group2, string name1, string name2) {
    int group1flag=0;
    int group2flag=0;
    int tmpint;
    bool isAtom;
    string line,compressedLine,item;
    stringstream linestream(line);
    vector <string> items;
    name1="["+name1+"]";
    name2="["+name2+"]";
  
  while (getline((*ndxfile), line)) {// for all lines in the file
        compressedLine=delete_char(line,' ');// if the line defines a group of atoms, e. g. [ Backbone ]
        if((compressedLine[0]=='[')&&(compressedLine[compressedLine.length()-1]==']')){group1flag=0;group2flag=0;}
        
        if(group1flag==1){// if the group has been identified as group1 to retrieve
          linestream.clear();
          linestream<<line;
           while (getline(linestream, item,' ')) {
            if(item.length()!=0){
            isAtom = true;// check if the item is a number
            for(unsigned int i=0; i<item.length();i++){isAtom=isAtom&&isdigit(item[i]);}
              if(!isAtom){cerr<<"ERROR READING INDEX FILE: ITEM \""<<item<<"\" IN SECTION \""<<name1<<"\" IS NOT A PROPER ATOM NUMBER!"<<endl;return 1;}
              sscanf (item.c_str(),"%d",&tmpint);
              group1->push_back(tmpint);//and add this number to group1
            }
          }
        }
        
        if(group2flag==1){// if the group has been identified as group2 to retrieve
          linestream.clear();
          linestream<<line;
           while (getline(linestream, item,' ')) {
            if(item.length()!=0){
            isAtom = true;// check if the item is a number
            for(unsigned int i=0; i<item.length();i++){isAtom=isAtom&&isdigit(item[i]);}
              if(!isAtom){cerr<<"ERROR READING INDEX FILE: ITEM \""<<item<<"\" IN SECTION \""<<name1<<"\" IS NOT A PROPER ATOM NUMBER!"<<endl;return 1;}
              sscanf (item.c_str(),"%d",&tmpint);
              group2->push_back(tmpint);//and add this number to group2
            }
          }
        }
        
        if(compressedLine==name1){group1flag=1;}//check if the line declares the requested group1
        if(compressedLine==name2){group2flag=1;}// or group2
    }


return 0;
}



