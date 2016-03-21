//    The binary IO library for the PARENT program suite
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
#include <stdlib.h>
#include <sstream>


using namespace std;

int write_BAT_header(ofstream *outfile,int double_prec,int numframes,vector< vector <int> > *dihedrals_top, vector <float>  *masses, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule);
int read_BAT_header(ifstream *infile,int *double_prec,int *numFrames,vector< vector <int> > *dihedrals_top, vector <float>  *masses, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule);
int write_BAT_frame(ofstream *outfile,int double_prec, int nDihedrals, float time, float xtcPrec, float **box, double* root_origin_cartesian,double root_origin_theta,double root_origin_phi,double root_origin_dihedral,double *bonds,double *angles,double *dihedrals);
int read_BAT_frame(ifstream *infile,int precision, int nDihedrals, float *time, float* xtcPrec, float **dbox, double* root_origin_cartesian,double* root_origin_theta,double* root_origin_phi,double* root_origin_dihedral,double *my_bonds,double *my_angles,double *my_dihedrals);



int write_PAR_header(ofstream *outfile,int nDihedrals,int double_prec,int numFrames,vector< vector <int> > *dihedrals_top, vector <float>  *masses, int bDens1D, int aDens1D, int dDens1D, int bDens, int aDens, int dDens, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule) ;
int read_PAR_header(ifstream *infile,int *nDihedrals,int *double_prec,int *numFrames,vector< vector <int> > *dihedrals_top, vector <float>  *masses, int *version, int* bDens, int* aDens, int* dDens, int* bDens1D, int* aDens1D, int* dDens1D, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule) ;
int write_PAR_body(ofstream* par_file, int nDihedrals,double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy);
int read_PAR_body(ifstream* par_file, int nDihedrals,double** bondsEntropy1D, double** anglesEntropy1D, double** dihedralsEntropy1D, double** bbEntropy, double** baEntropy, double** bdEntropy, double** aaEntropy, double** adEntropy, double** ddEntropy);





