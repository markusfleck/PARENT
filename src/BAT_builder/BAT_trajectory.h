//    Class for BAT_builder, a program to convert a molecular dynamics trajectory from Cartesian to internal bond-angle-torsion coordinates
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







#ifndef BAT_TRAJECTORY_H
#define BAT_TRAJECTORY_H


#include <vector>
#include <string>
#include "../util/io/xdrfile/xdrfile_xtc.h"



class BAT_Trajectory {

public:
    //constructor for the .xtc to .bat conversion
    BAT_Trajectory(std::string trjFileI,std::string trjFileO, std::vector< std::vector <int> > *dihedralsIn, std::vector <float> *massvec_in, std::vector <std::string> *residues_in,std::vector <int> *residueNumbers_in,std::vector <std::string> *atomNames_in,std::vector <std::string> *belongsToMolecule_in, int double_prec_in=1);
    //constructor for the .bat to .xtc conversion
    BAT_Trajectory(std::string trjFileI,std::string trjFileO);

    int double_prec;
    std::vector< std::vector <int> > dihedrals_top;
    std::vector <float> massvec;
    std::vector <std::string> residues;
    std::vector <int> residueNumbers;
    std::vector <std::string> atomNames;
    std::vector <std::string> belongsToMolecule;
    std::string trjFileIn,trjFileOut;
    rvec *x;
    double root_origin_cartesian[3];
    double root_origin_phi;
    double root_origin_theta;
    double root_origin_dihedral;
    double rootbond[2];
    double rootangle;
    double *bonds;
    double *angles;
    double *dihedrals;
    double *bondsFull,*anglesFull;
    float *masses;
    int step;
    float time,prec;
    matrix box;
    double pi;
    int numframes;
    
    int convert_xtc_to_BAT();
    int convert_BAT_to_xtc();
    void shift_dihedrals_0_2pi_range();
    void convert_to_BAT();
    void convert_dihedrals_to_phaseangles();//a phaseangle is formed by its difference to its "parent" dihedral
    void convert_to_Cartesian();
    void convert_phaseangles_to_dihedrals();//to convert all phaseangles back to dihedrals

};



#endif

