//    A program to extract the entropy and mutual information terms from the binary output of the PARENT program
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
#include <vector>
#include <cstring>
#include <cstdlib>


using namespace std;



#include "../util/io/io.h"
#include "../util/util.h"



int main(int argc, char* argv[])
{
    cout.precision(12);


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

    double mutual, entropy2D;
		
		bool longOutput;
		char* inputFilename;
		
		if((argc==3)&&(cmdOptionExists(argv, argv+argc, "-p"))){
			inputFilename = getCmdOption(argv, argv+argc, "-p");
      longOutput=true;  
    }
		else if((argc==4)&&cmdOptionExists(argv, argv+argc, "-p")&&cmdOptionExists(argv, argv+argc, "--short")){
			inputFilename = getCmdOption(argv, argv+argc, "-p");
      longOutput=false;  
    }
		else{
			cerr<<"USAGE:\n"<<argv[0]<<" -p input.par [--short]"<<endl;
			return 1;
		}
			

    try{
      EntropyMatrix mat(inputFilename);

    
      int nBonds=mat.getNBonds();
      int nAngles=mat.getNAngles();
      int nDihedrals=mat.getNDihedrals();


      for(int k=1; k<=nBonds; k++) {
          if(longOutput){
						cout<<"bond "<<k<<"   "<<mat.getEntropy(TYPE_B,k)<<endl;    //write the 1D entropies to stdout  and keep track of the total sums
					}
          totalBondsEntropy+=mat.getEntropy(TYPE_B,k);
      }
      for(int k=1; k<=nAngles; k++) {
					if(longOutput){
						cout<<"angle "<<k<<"   "<<mat.getEntropy(TYPE_A,k)<<endl;
					}
					totalAnglesEntropy+=mat.getEntropy(TYPE_A,k);
      }
      for(int k=1; k<=nDihedrals; k++) {
					if(longOutput){
						cout<<"dihedral "<<k<<"   "<<mat.getEntropy(TYPE_D,k)<<endl;
					}
          totalDihedralsEntropy+=mat.getEntropy(TYPE_D,k);
      }

          //then for every bond-bond pair get the mutual information and write it as well as the 2D entropy to stdout
          for(int i=1; i<=nBonds; i++) {
              for(int j=i+1; j<=nBonds; j++) {
                  mutual=mat.getMutual(TYPE_B,TYPE_B,i,j);
                  entropy2D=mat.get2DEntropy(TYPE_B,TYPE_B,i,j);
                  if(longOutput){
										cout<<"bond-bond "<<i<<"  "<<j<<"   "<<entropy2D<<"   "<<mutual<<endl;
                  }
									totalBbEntropy+=entropy2D; //also keep track of the total bond-bond enropy sum
                  totalBbMutual+=mutual; //and the bond-bond mutual information sum
              }
          }

          for(int i=1; i<=nBonds; i++) { //same for all bond-angle pairs, except that the array is in forward order
              for(int j=1; j<=nAngles; j++) {
                  mutual=mat.getMutual(TYPE_B,TYPE_A,i,j);
                  entropy2D=mat.get2DEntropy(TYPE_B,TYPE_A,i,j);
                  if(longOutput){
										cout<<"bond-angle "<<i<<"  "<<j<<"   "<<entropy2D<<"   "<<mutual<<endl;
                  }
									totalBaEntropy+=entropy2D;
                  totalBaMutual+=mutual;
              }
          }

          for(int i=1; i<=nBonds; i++) { //bond-dihedral pairs (forward order)
              for(int j=1; j<=nDihedrals; j++) {
                  mutual=mat.getMutual(TYPE_B,TYPE_D,i,j);
                  entropy2D=mat.get2DEntropy(TYPE_B,TYPE_D,i,j);
                  if(longOutput){
										cout<<"bond-dihedral "<<i<<"  "<<j<<"   "<<entropy2D<<"   "<<mutual<<endl;
                  }
									totalBdEntropy+=entropy2D;
                  totalBdMutual+=mutual;
              }
          }

          for(int i=1; i<=nAngles-1; i++) {//angle-angle pairs (reverse order)
              for(int j=i+1; j<=nAngles; j++) {
                  mutual=mat.getMutual(TYPE_A,TYPE_A,i,j);
                  entropy2D=mat.get2DEntropy(TYPE_A,TYPE_A,i,j);
									if(longOutput){	
										cout<<"angle-angle "<<i<<"  "<<j<<"   "<<entropy2D<<"   "<<mutual<<endl;
                  }
									totalAaEntropy+=entropy2D;
                  totalAaMutual+=mutual;
              }
          }

          for(int i=1; i<=nAngles; i++) { //angle-dihedral pairs (forward order)
              for(int j=1; j<=nDihedrals; j++) {
                  mutual=mat.getMutual(TYPE_A,TYPE_D,i,j);
                  entropy2D=mat.get2DEntropy(TYPE_A,TYPE_D,i,j);
                  if(longOutput){
										cout<<"angle-dihedral "<<i<<"  "<<j<<"   "<<entropy2D<<"   "<<mutual<<endl;
                  }
									totalAdEntropy+=entropy2D;
                  totalAdMutual+=mutual;
              }
          }

          for(int i=1; i<=nDihedrals-1; i++) { //dihedral-dihedral pairs (reverse order)
              for(int j=i+1; j<=nDihedrals; j++) {
                  mutual=mat.getMutual(TYPE_D,TYPE_D,i,j);
                  entropy2D=mat.get2DEntropy(TYPE_D,TYPE_D,i,j);
									if(longOutput){
										cout<<"dihedral-dihedral "<<i<<"  "<<j<<"   "<<entropy2D<<"   "<<mutual<<endl;
                  }
									totalDdEntropy+=entropy2D;
                  totalDdMutual+=mutual;
              }
          }

          //in the end write all quantities to stdout
          cout<<"TOTAL 1D BONDS ENTROPY = "<<totalBondsEntropy<<endl;
          cout<<"TOTAL 1D ANGLES ENTROPY = "<<totalAnglesEntropy<<endl;
          cout<<"TOTAL 1D DIHEDRALS ENTROPY = "<<totalDihedralsEntropy<<endl;
          //cout<<"TOTAL 2D BONDS-BONDS ENTROPY = "<<totalBbEntropy<<endl;
          cout<<"TOTAL 2D BONDS-BONDS MUTUAL INFORMATION = "<<totalBbMutual<<endl;
          //cout<<"TOTAL 2D BONDS-ANGLES ENTROPY = "<<totalBaEntropy<<endl;
          cout<<"TOTAL 2D BONDS-ANGLES MUTUAL INFORMATION = "<<totalBaMutual<<endl;
          //cout<<"TOTAL 2D BONDS-DIHEDRALS ENTROPY = "<<totalBdEntropy<<endl;
          cout<<"TOTAL 2D BONDS-DIHEDRALS MUTUAL INFORMATION = "<<totalBdMutual<<endl;
          //cout<<"TOTAL 2D ANGLES-ANGLES ENTROPY = "<<totalAaEntropy<<endl;
          cout<<"TOTAL 2D ANGLES-ANGLES MUTUAL INFORMATION = "<<totalAaMutual<<endl;
          //cout<<"TOTAL 2D ANGLES-DIHEDRALS ENTROPY = "<<totalAdEntropy<<endl;
          cout<<"TOTAL 2D ANGLES-DIHEDRALS MUTUAL INFORMATION = "<<totalAdMutual<<endl;
          //cout<<"TOTAL 2D DIHEDRALS-DIHEDRALS ENTROPY = "<<totalDdEntropy<<endl;
          cout<<"TOTAL 2D DIHEDRALS-DIHEDRALS MUTUAL INFORMATION = "<<totalDdMutual<<endl;
          cout<<"TOTAL CONFIGURATIONAL ENTROPY = "<<totalBondsEntropy+totalAnglesEntropy+totalDihedralsEntropy-totalBbMutual-totalBaMutual-totalBdMutual-totalAaMutual-totalAdMutual-totalDdMutual<<endl;





    }
    catch(MyError myError){
      cerr<<myError.what()<<endl;
			cerr<<"USAGE:\n"<<argv[0]<<" -p input.par [--short]"<<endl;
      return 1;
    }
    catch(...)
    {
      cerr<<"AN UNIDENTIFIED ERROR HAS OCCURRED! ABORTING.\n"<<endl;
			cerr<<"USAGE:\n"<<argv[0]<<" -p input.par [--short]"<<endl;
      return 1;
    }
  return 0;
}







