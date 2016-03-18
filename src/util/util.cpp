//    The utility library for the PARENT program suite
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







#include "util.h"

double get_mutual(int type,int index1, int index2,int nDihedrals, double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy) {
    int smaller,bigger,index;
    int nBonds=nDihedrals+2;
    int nAngles=nDihedrals+1;

    if((type==TYPE_BB)&&(index1!=index2)) { // for mutual information between two different bonds
        smaller=index1<index2?index1:index2;
        bigger=index1<index2?index2:index1;
        index=(nBonds-smaller)*(nBonds-smaller-1)/2+smaller-bigger;//the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
        return bondsEntropy1D[smaller]+bondsEntropy1D[bigger]-bbEntropy[index];
    }
    if(type==TYPE_BA) { // for mutual information between a bond and an angle
        return bondsEntropy1D[index1]+anglesEntropy1D[index2]-baEntropy[index1*nAngles+index2];
    }
    if(type==TYPE_BD) { // for mutual information between a bond and a dihedral
        return bondsEntropy1D[index1]+dihedralsEntropy1D[index2]-bdEntropy[index1*nDihedrals+index2];
    }
    if((type==TYPE_AA)&&(index1!=index2)) { // for mutual information between two different angles
        smaller=index1<index2?index1:index2;
        bigger=index1<index2?index2:index1;
        index=(nAngles-smaller)*(nAngles-smaller-1)/2+smaller-bigger;//the 2D-entropies angles-angles were stored in reverse order, as documented in "Parent.cpp"
        return anglesEntropy1D[smaller]+anglesEntropy1D[bigger]-aaEntropy[index];
    }
    if(type==TYPE_AD) { // for mutual information between a bond and an angle
        return anglesEntropy1D[index1]+dihedralsEntropy1D[index2]-adEntropy[index1*nDihedrals+index2];
    }
    if((type==TYPE_DD)&&(index1!=index2)) { // for mutual information between two different angles
        smaller=index1<index2?index1:index2;
        bigger=index1<index2?index2:index1;
        index=(nDihedrals-smaller)*(nDihedrals-smaller-1)/2+smaller-bigger;//the 2D-entropies dihedrals-dihedrals were stored in reverse order, as documented in "Parent.cpp"
        return dihedralsEntropy1D[smaller]+dihedralsEntropy1D[bigger]-ddEntropy[index];
    }
    cerr<<"WARNING: REQUEST FOR MUTUAL INFORMATION TYPE "<<type<<" WITH INDICES "<<index1<<" AND "<<index2<<"YIELDED NO RESULT."<<endl;
    return 0;
}

