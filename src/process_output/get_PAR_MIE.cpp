//    A program to extract the entropy and mutual information terms from the binary output of the PARENT program
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







#define MUTUAL_ARGS nDihedrals,bondsEntropy1D,anglesEntropy1D,dihedralsEntropy1D,bbEntropy,baEntropy,bdEntropy,aaEntropy,adEntropy,ddEntropy


#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <cstdlib>


using namespace std;



#include "../util/io/io.h"
#include "../util/util.h"



int main(int argc, char* argv[])
{
    cout.precision(12);

    int double_prec,numFrames;
    int version,bDens,aDens,dDens,bDens1D,aDens1D,dDens1D;
    vector< vector <int> > dihedrals_top;
    vector <float> masses;
    vector <string> residues;
    vector <int> residueNumbers;
    vector <string> atomNames;
    vector <string> belongsToMolecule;

    double *bondsEntropy1D,*anglesEntropy1D,*dihedralsEntropy1D;
    double *bbEntropy,*aaEntropy,*ddEntropy,*baEntropy,*bdEntropy,*adEntropy;

    int nBonds,nAngles,nDihedrals;

    double totalBondsEntropy=0;
    double totalAnglesEntropy=0;
    double totalDihedralsEntropy=0;

    double totalBbEntropy=0;
    double totalBbMutual=0;
    double totalBaEntropy=0;
    double totalBaMutual=0;
    double totalBdEntropy=0;
    double totalBdMutual=0;
    double totalAaEntropy=0;
    double totalAaMutual=0;
    double totalAdEntropy=0;
    double totalAdMutual=0;
    double totalDdEntropy=0;
    double totalDdMutual=0;

    double mutual;



    ifstream infile(argv[1], ios::binary | ios::in);//open the .par file
    if(infile.is_open()) {
        if(read_PAR_header(&infile,&nDihedrals,&double_prec,&numFrames,&dihedrals_top, &masses, &version, &bDens, &aDens, &dDens, &bDens1D, &aDens1D, &dDens1D,&residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) { //and read the header information
            cerr<<"ERROR READING HEADER OF FILE "<<argv[1]<<" !"<<endl;
            return 1;
        }

        nBonds=dihedrals_top.size()+2;
        nAngles=dihedrals_top.size()+1;
        nDihedrals=dihedrals_top.size();

        if (read_PAR_body(&infile,nDihedrals,&bondsEntropy1D, &anglesEntropy1D, &dihedralsEntropy1D, &bbEntropy, &baEntropy, &bdEntropy, &aaEntropy, &adEntropy, &ddEntropy)!=0) {
            cerr<<"ERROR READING FILE "<<argv[1]<<" !"<<endl;
            return 1;
        }



        for(int k=0; k<nBonds; k++) {
            cout<<"bond "<<k+1<<"   "<<bondsEntropy1D[k]<<endl;    //write the 1D entropies to stdout  and keep track of the total sums
            totalBondsEntropy+=bondsEntropy1D[k];
        }
        for(int k=0; k<nAngles; k++) {
            cout<<"angle "<<k+1<<"   "<<anglesEntropy1D[k]<<endl;
            totalAnglesEntropy+=anglesEntropy1D[k];
        }
        for(int k=0; k<nDihedrals; k++) {
            cout<<"dihedral "<<k+1<<"   "<<dihedralsEntropy1D[k]<<endl;
            totalDihedralsEntropy+=dihedralsEntropy1D[k];
        }

        //then for every bond-bond pair calculate the mutual information and write it as well as the 2D entropy to stdout
        int counter=nBonds*(nBonds-1)/2-1; //remember the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
        for(int i=0; i<nBonds-1; i++) {
            for(int j=i+1; j<nBonds; j++) {
                mutual=get_mutual(TYPE_BB,i,j,MUTUAL_ARGS);
                cout<<"bond-bond "<<i+1<<"  "<<j+1<<"   "<<bbEntropy[counter]<<"   "<<mutual<<endl;
                totalBbEntropy+=bbEntropy[counter]; //also keep track of the total bond-bond enropy sum
                totalBbMutual+=mutual; //and the bond-bond mutual information sum
                counter--;
            }
        }

        for(int i=0; i<nBonds; i++) { //same for all bond-angle pairs, except that the array is in forward order
            for(int j=0; j<nAngles; j++) {
                mutual=get_mutual(TYPE_BA,i,j,MUTUAL_ARGS);
                cout<<"bond-angle "<<i+1<<"  "<<j+1<<"   "<<baEntropy[i*nAngles+j]<<"   "<<mutual<<endl;
                totalBaEntropy+=baEntropy[i*nAngles+j];
                totalBaMutual+=mutual;
            }
        }

        for(int i=0; i<nBonds; i++) { //bond-dihedral pairs (forward order)
            for(int j=0; j<nDihedrals; j++) {
                mutual=get_mutual(TYPE_BD,i,j,MUTUAL_ARGS);
                cout<<"bond-dihedral "<<i+1<<"  "<<j+1<<"   "<<bdEntropy[i*nDihedrals+j]<<"   "<<mutual<<endl;
                totalBdEntropy+=bdEntropy[i*nDihedrals+j];
                totalBdMutual+=mutual;
            }
        }

        counter=nAngles*(nAngles-1)/2-1; //angle-angle pairs (reverse order)
        for(int i=0; i<nAngles-1; i++) {
            for(int j=i+1; j<nAngles; j++) {
                mutual=get_mutual(TYPE_AA,i,j,MUTUAL_ARGS);
                cout<<"angle-angle "<<i+1<<"  "<<j+1<<"   "<<aaEntropy[counter]<<"   "<<mutual<<endl;
                totalAaEntropy+=aaEntropy[counter];
                totalAaMutual+=mutual;
                counter--;
            }
        }

        for(int i=0; i<nAngles; i++) { //angle-dihedral pairs (forward order)
            for(int j=0; j<nDihedrals; j++) {
                mutual=get_mutual(TYPE_AD,i,j,MUTUAL_ARGS);
                cout<<"angle-dihedral "<<i+1<<"  "<<j+1<<"   "<<adEntropy[i*nDihedrals+j]<<"   "<<mutual<<endl;
                totalAdEntropy+=adEntropy[i*nDihedrals+j];
                totalAdMutual+=mutual;
            }
        }

        counter=nDihedrals*(nDihedrals-1)/2-1; //dihedral-dihedral pairs (reverse order)
        for(int i=0; i<nDihedrals-1; i++) {
            for(int j=i+1; j<nDihedrals; j++) {
                mutual=get_mutual(TYPE_DD,i,j,MUTUAL_ARGS);
                cout<<"dihedral-dihedral "<<i+1<<"  "<<j+1<<"   "<<ddEntropy[counter]<<"   "<<mutual<<endl;
                totalDdEntropy+=ddEntropy[counter];
                totalDdMutual+=mutual;
                counter--;
            }
        }

        //in the end write all quantities to stdout
        cout<<"TOTAL 1D BONDS ENTROPY = "<<totalBondsEntropy<<endl;
        cout<<"TOTAL 1D ANGLES ENTROPY = "<<totalAnglesEntropy<<endl;
        cout<<"TOTAL 1D DIHEDRALS ENTROPY = "<<totalDihedralsEntropy<<endl;
        cout<<"TOTAL 2D BONDS-BONDS ENTROPY = "<<totalBbEntropy<<endl;
        cout<<"TOTAL 2D BONDS-BONDS MUTUAL INFORMATION = "<<totalBbMutual<<endl;
        cout<<"TOTAL 2D BONDS-ANGLES ENTROPY = "<<totalBaEntropy<<endl;
        cout<<"TOTAL 2D BONDS-ANGLES MUTUAL INFORMATION = "<<totalBaMutual<<endl;
        cout<<"TOTAL 2D BONDS-DIHEDRALS ENTROPY = "<<totalBdEntropy<<endl;
        cout<<"TOTAL 2D BONDS-DIHEDRALS MUTUAL INFORMATION = "<<totalBdMutual<<endl;
        cout<<"TOTAL 2D ANGLES-ANGLES ENTROPY = "<<totalAaEntropy<<endl;
        cout<<"TOTAL 2D ANGLES-ANGLES MUTUAL INFORMATION = "<<totalAaMutual<<endl;
        cout<<"TOTAL 2D ANGLES-DIHEDRALS ENTROPY = "<<totalAdEntropy<<endl;
        cout<<"TOTAL 2D ANGLES-DIHEDRALS MUTUAL INFORMATION = "<<totalAdMutual<<endl;
        cout<<"TOTAL 2D DIHEDRALS-DIHEDRALS ENTROPY = "<<totalDdEntropy<<endl;
        cout<<"TOTAL 2D DIHEDRALS-DIHEDRALS MUTUAL INFORMATION = "<<totalDdMutual<<endl;
        cout<<"TOTAL CONFIGURATIONAL ENTROPY = "<<totalBondsEntropy+totalAnglesEntropy+totalDihedralsEntropy-totalBbMutual-totalBaMutual-totalBdMutual-totalAaMutual-totalAdMutual-totalDdMutual<<endl;




    }
    else {
        cerr<<"ERROR: COULD NOT OPEN FILE !\n\nUSAGE:\n"<<argv[0]<<" input.par"<<endl;
        return 1;
    }
    infile.close();


    return 0;
}







