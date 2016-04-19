//    The binary IO library for the PARENT program suite
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







#ifndef IO_BINARY_H
#define IO_BINARY_H

#define TYPE_B 0
#define TYPE_A 1
#define TYPE_D 2

#define TYPE_BB 0
#define TYPE_BA 1
#define TYPE_BD 2
#define TYPE_AA 3
#define TYPE_AD 4
#define TYPE_DD 5


#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


int write_BAT_header(std::ofstream *outfile,int double_prec,int numframes,std::vector< std::vector <int> > *dihedrals_top, std::vector <float>  *masses, std::vector <std::string> *residues,std::vector <int> *residueNumbers,std::vector <std::string> *atomNames,std::vector <std::string> *belongsToMolecule);
int read_BAT_header(std::ifstream *infile,int *double_prec,int *numFrames,std::vector< std::vector <int> > *dihedrals_top, std::vector <float>  *masses, std::vector <std::string> *residues,std::vector <int> *residueNumbers,std::vector <std::string> *atomNames,std::vector <std::string> *belongsToMolecule);
int write_BAT_frame(std::ofstream *outfile,int double_prec, int nDihedrals, float time, float xtcPrec, float **box, double* root_origin_cartesian,double root_origin_theta,double root_origin_phi,double root_origin_dihedral,double *bonds,double *angles,double *dihedrals);
int read_BAT_frame(std::ifstream *infile,int precision, int nDihedrals, float *time, float* xtcPrec, float **dbox, double* root_origin_cartesian,double* root_origin_theta,double* root_origin_phi,double* root_origin_dihedral,double *my_bonds,double *my_angles,double *my_dihedrals);



int write_PAR_header(std::ofstream *outfile,int nDihedrals,int double_prec,int numFrames,std::vector< std::vector <int> > *dihedrals_top, std::vector <float>  *masses, int bDens1D, int aDens1D, int dDens1D, int bDens, int aDens, int dDens, std::vector <std::string> *residues,std::vector <int> *residueNumbers,std::vector <std::string> *atomNames,std::vector <std::string> *belongsToMolecule) ;
int read_PAR_header(std::ifstream *infile,int *nDihedrals,int *double_prec,int *numFrames,std::vector< std::vector <int> > *dihedrals_top, std::vector <float>  *masses, int *version, int* bDens, int* aDens, int* dDens, int* bDens1D, int* aDens1D, int* dDens1D, std::vector <std::string> *residues,std::vector <int> *residueNumbers,std::vector <std::string> *atomNames,std::vector <std::string> *belongsToMolecule) ;
int write_PAR_body(std::ofstream* par_file, int nDihedrals,double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy);
int read_PAR_body(std::ifstream* par_file, int nDihedrals,double** bondsEntropy1D, double** anglesEntropy1D, double** dihedralsEntropy1D, double** bbEntropy, double** baEntropy, double** bdEntropy, double** aaEntropy, double** adEntropy, double** ddEntropy);



int write_CLT_file(std::ofstream *clt_file,int calcBonds, int calcAngles, int calcDihedrals, std::vector <std::string> *residuesRes,std::vector <int> *residueNumbersRes,std::vector <std::string> *belongsToMoleculeRes, std::vector <int> *nBonds, std::vector <int> *nAngles, std::vector <int> *nDihedrals, double* mutualArray);
int read_CLT_file(std::ifstream *clt_file, int *version, int* calcBonds, int* calcAngles, int* calcDihedrals, std::vector <std::string> *residues, std::vector <int> *residueNumbers, std::vector <std::string> *belongsToMolecule, std::vector <int> *nBonds, std::vector <int> *nAngles, std::vector <int> *nDihedrals, double** mutualArray) ;


class EntropyMatrix {
    public:
      EntropyMatrix(char const * infileInput);
      EntropyMatrix(int nAtoms);
      ~EntropyMatrix();
    
    double getEntropy(int type, int index);
    double get2DEntropy(int type1, int type2, int index1, int index2);
    double getMutual(int type1, int type2, int index1, int index2);
		
		void setEntropy(int type, int index, double value);
    void set2DEntropy(int type1, int type2, int index1, int index2, double value);
    void setMutual(int type1, int type2, int index1, int index2, double value); //modifies the 2D entropy to yield the according mutual information (without modifying the 1D entropy values)
    
    
    void write(char const * infileInput);
    
    int getNBonds();
    int getNAngles();
    int getNDihedrals();
    
    private:
      void write_PAR_header(); 
      void read_PAR_header(); 
      void write_PAR_body();
      void read_PAR_body(); 
    
      int nBonds, nAngles, nDihedrals;
      int double_prec, numFrames,version;
      double *bondsEntropy1D, *anglesEntropy1D, *dihedralsEntropy1D; 
      double *bbEntropy, *baEntropy, *bdEntropy, *aaEntropy, *adEntropy, *ddEntropy;
      std::vector< std::vector <int> > dihedrals_top; 
      std::vector <float>  masses;
      int bDens1D, aDens1D, dDens1D, bDens, aDens, dDens; 
      std::vector <std::string> residues;
      std::vector <int> residueNumbers;
      std::vector <std::string> atomNames;
      std::vector <std::string> belongsToMolecule;
      std::string infilestring;
      std::ifstream infile;
      std::ofstream outfile;
};




  class MyError : public std::runtime_error {
  public:
    MyError(const std::string& msg = "") : std::runtime_error(msg) {}
  };
  
  
  #endif


