//    The ASCII IO library for the PARENT program suite
//    Copyright (C) 2015  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)
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
