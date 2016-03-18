//    A program to extract the BAT topology from the binary output of the PARENT program
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







#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>



using namespace std;


#include "../util/io/io.h"


int main(int argc, char* argv[])
{
    int nDihedrals,double_prec,numFrames;
    int version,bDens,aDens,dDens,bDens1D,aDens1D,dDens1D;
    vector< vector <int> > dihedrals_top;
    vector <float> masses;
    vector <string> residues;
    vector <int> residueNumbers;
    vector <string> atomNames;
    vector <string> belongsToMolecule;

    ifstream infile(argv[1], ios::binary | ios::in);
    if(infile.is_open()&&argc==2) {
        //first read the par header
        if(0!=read_PAR_header(&infile,&nDihedrals,&double_prec,&numFrames,&dihedrals_top, &masses, &version, &bDens, &aDens, &dDens, &bDens1D, &aDens1D, &dDens1D,&residues,&residueNumbers,&atomNames,&belongsToMolecule)) {
            cerr<<"ERROR READING HEADER OF FILE "<<argv[1]<<" !"<<endl;
            return 1;
        }
        cout<<"#bond #atom1 #atom2 "<<endl;
        cout<<"1 "<<dihedrals_top[0][0]<<" "<<dihedrals_top[0][1]<<endl; //the first bond consists of the atoms 1 and 2 of dihedral 1 (arrays start at 0, but numbering starts at 1)
        cout<<"2 "<<dihedrals_top[0][1]<<" "<<dihedrals_top[0][2]<<endl; //the second bond consists of the atoms 2 and 3 of dihedral 1
        for(int i=0; i<nDihedrals; i++) {
            cout<<i+3<<" "<<dihedrals_top[i][2]<<" "<<dihedrals_top[i][3]<<endl;   //all other bonds consist of the atoms 3 and 4 of every dihedral
        }
        cout<<"#angle #atom1 #atom2 #atom3"<<endl;
        cout<<"1 "<<dihedrals_top[0][0]<<" "<<dihedrals_top[0][1]<<" "<<dihedrals_top[0][2]<<endl;//the first angle consists of the atoms 1,2 and 3 of dihedral 1 (arrays start at 0, but numbering starts at 1)
        for(int i=0; i<nDihedrals; i++) {
            cout<<i+2<<" "<<dihedrals_top[i][1]<<" "<<dihedrals_top[i][2]<<" "<<dihedrals_top[i][3]<<endl;   //all other angles consist of the atoms 2,3 and 4 of every dihedral
        }
        cout<<"#dihedral #atom1 #atom2 #atom3 #atom4    dihedral_type(common=0,pseudo=1,improper=-1)   phaseangle_of(-1 if no phaseangle)"<<endl;
        //atoms 1,2,3 and 4 are output for every dihedral, plus an integer
        for(int i=0; i<nDihedrals; i++) {
            cout<<i+1<<" "<<dihedrals_top[i][0]<<" "<<dihedrals_top[i][1]<<" "<<dihedrals_top[i][2]<<" "<<dihedrals_top[i][3]<<" "<<dihedrals_top[i][4]<<" "<<dihedrals_top[i][5]<<endl;
        }
    }
    else {
        if(argc==2) {
            cerr<<"ERROR: COULD NOT OPEN FILE "<<argv[1]<<" !"<<endl;
        }
        else {
            cerr<<"usage: "<<argv[0]<<" entropyfile.par"<<endl;
        }
        return 1;
    }
    infile.close();


    return 0;
}







