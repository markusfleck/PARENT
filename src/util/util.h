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







#define TYPE_B 0
#define TYPE_A 1
#define TYPE_D 2

#define TYPE_BB 0
#define TYPE_BA 1
#define TYPE_BD 2
#define TYPE_AA 3
#define TYPE_AD 4
#define TYPE_DD 5


#include <iostream>


using namespace std;

double get_mutual(int type,int index1, int index2,int nDihedrals, double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy);



