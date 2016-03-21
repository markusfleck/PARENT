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







#pragma pack(1)



#include "io_binary.h"

using namespace std;


//to write the header of the binary .bat file
int write_BAT_header(ofstream *outfile,int double_prec,int numframes,vector< vector <int> > *dihedrals_top, vector <float>  *masses, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule) {
    int dummy=(*dihedrals_top).size();
    int version=3;
    int fail=0;
    char dummystring[31];

    //WHEN CHANGING SOMETHING HERE, REMEMBER TO UPDATE THE LINE "outfile2.seekp(3*sizeof(int), ios::beg);" IN FUNCTION convert_xtc_to_BAT() file BAT_trajectory.cpp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    (*outfile).write((char*)&version, sizeof(int)); //first write the .bat version number as an integer
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&double_prec, sizeof(int)); //then write an integer declaring if the trajectory is stored in double precision
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&dummy, sizeof(int)); //write an integer containing the number of dihedrals
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&numframes, sizeof(int)); //and an integer containing the number of frames
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);

    for(int i=0; i<dummy+3; i++) { //for all atoms in the system
        for(int j=0; j<8; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.7s", ((*residues)[i]).c_str());
        (*outfile).write(dummystring, 8*sizeof(char));// write the name of the residue it belongs to
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*residueNumbers)[i]), sizeof(int));// write the number of the residue it belongs to
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        for(int j=0; j<8; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.7s", ((*atomNames)[i]).c_str());// write the name of the atom
        (*outfile).write(dummystring, 8*sizeof(char));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        for(int j=0; j<31; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.30s", ((*belongsToMolecule)[i]).c_str());
        (*outfile).write(dummystring, 31*sizeof(char));// write the name of the molecule it belongs to
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    }

    for(int i=0; i<dummy; i++) {//then for all dihedrals
        (*outfile).write((char*)&(*dihedrals_top)[i][0], sizeof(int)); //write the atomnumber of the first atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&(*dihedrals_top)[i][1], sizeof(int)); //second atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&(*dihedrals_top)[i][2], sizeof(int)); //third atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&(*dihedrals_top)[i][3], sizeof(int)); //fourth atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&(*dihedrals_top)[i][4], sizeof(int)); // an integer containing the type of the dihedral(physical=0,pseudo=1,improper=-1)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&(*dihedrals_top)[i][5], sizeof(int)); // and an integer containing the "parent"dihedral(for phaseangles, -1 if no "parent")
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);

    }
    for(int i=0; i<dummy+3; i++) {
        (*outfile).write((char*)&(*masses)[i], sizeof(float)); //and write the whole massestor of the atoms in the system in single precision (float)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    }
    return fail;
}



//to read the header of the binary .bat file
int read_BAT_header(ifstream *infile,int *double_prec,int *numFrames,vector< vector <int> > *dihedrals_top, vector <float>  *masses, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule) {
    int version;
    int fail=0;
    char dummystring[31];
    int nDihedrals;

    (*infile).read((char*)&version, sizeof(int)); //first read the .bat version number as an integer
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)double_prec, sizeof(int)); //then read an integer declaring if the trajectory is stored in double precision
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)&nDihedrals,sizeof(int));//read an integer containing the number of dihedrals
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)numFrames,sizeof(int)); //and an integer containing the number of frames
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);


    //check if everything went okay
    if(fail!=0) {
        return fail;
    }
    if(version<0) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. VERSION NUMBER ("<<version<<") < 0!"<<endl;
        return 1;
    }
    if((*double_prec!=0)&&(*double_prec!=1)) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. DOUBLE PRECISION VALUE ("<<(*double_prec)<<") NEITHER 0 NOR 1!"<<endl;
        return 1;
    }
    if(nDihedrals>19997) {
        cerr<<"WARNING: "<<nDihedrals+3<<" ATOMS DECLARED IN THE FILE HEADER (CORRUPTED?). THIS WILL LEAD TO LARGE OUTPUT."<<endl;
    }
    if(nDihedrals<0) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER OF DIHEDRALS ("<<nDihedrals<<") < 0!"<<endl;
        return 1;
    }
    if((*numFrames)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER OF FRAMES ("<<(*numFrames)<<") < 1!"<<endl;
        return 1;
    }

    if (version>=3) {
        for(int i=0; i<nDihedrals+3; i++) { //for ever atom in the system
            (*infile).read(dummystring, 8*sizeof(char));//read the name of the residue it belongs to
            (*residues).push_back(dummystring);
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            (*residueNumbers).push_back(0);
            (*infile).read((char*)&((*residueNumbers)[i]), sizeof(float));//read the number of the residue it belongs to
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            (*infile).read(dummystring, 8*sizeof(char));//read the name of the atom
            (*atomNames).push_back(dummystring);
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            (*infile).read(dummystring, 31*sizeof(char));//read the molecule of the residue it belongs to
            (*belongsToMolecule).push_back(dummystring);
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        }
    }


    //check if everything went okay
    if(fail!=0) {
        return fail;
    }


    vector<int>dummyvec;
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    for(int i=0; i<nDihedrals; i++) {//then for all dihedrals
        (*dihedrals_top).push_back(dummyvec);
        (*infile).read((char*)&((*dihedrals_top)[i][0]),sizeof(int));//read the atomnumber of the first atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)&((*dihedrals_top)[i][1]),sizeof(int));//second atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)&((*dihedrals_top)[i][2]),sizeof(int));//third atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)&((*dihedrals_top)[i][3]),sizeof(int));//fourth atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)&((*dihedrals_top)[i][4]),sizeof(int));//an integer containing the type of the dihedral(physical=0,pseudo=1,improper=-1)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)&((*dihedrals_top)[i][5]),sizeof(int));//and an integer containing the "parent"dihedral(for phaseangles, -1 if no "parent")
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    }
    for(int i=0; i<nDihedrals+3; i++) { //and read the whole massestor of the atoms in the system in single precision (float)
        (*masses).push_back(0);
        (*infile).read((char*)&((*masses)[i]), sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    }

    return fail; //if anything failed return a 1, otherwise a 0
}


int write_BAT_frame(ofstream *outfile,int double_prec, int nDihedrals, float time, float xtcPrec, float **box, double* root_origin_cartesian,double root_origin_theta,double root_origin_phi,double root_origin_dihedral,double *bonds,double *angles,double *dihedrals) {
    int fail=0;
    //to attach a frame to the .bat trajectory
    (*outfile).write((char*)&time, sizeof(float)); //write the time of the current according .xtc frame as a float
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&xtcPrec, sizeof(float)); //write the precision of the current according .xtc frame for back conversion
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)box, 9*sizeof(float)); //write the box vectors of the current according .xtc frame for back conversion
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);

    if(double_prec==1) { //if double precision is requested
        (*outfile).write((char*)root_origin_cartesian, 3*sizeof(double)); //write the Cartestians of the first root atom (external coordinates) in double precision
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&root_origin_theta, sizeof(double)); //write the polar coordinates of the second root atom relative to the first (external coordinates)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&root_origin_phi, sizeof(double));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&root_origin_dihedral, sizeof(double)); //and the dihedral the root atoms form with the origin (external coordinates)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)bonds, 2*sizeof(double)); //write the lengths of the two bonds connecting the root atoms (internal coordinates)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)angles, sizeof(double)); //and the angle between the two rootbonds (internal coordinates)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        for(int i=0; i<nDihedrals; i++) { //then for all dihedrals in the system
            (*outfile).write((char*)&bonds[i+2], sizeof(double)); //write the bondlength between the last two atoms in the dihedral
            fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
            (*outfile).write((char*)&angles[i+1], sizeof(double)); //write the angle between the last threee atoms of the dihedral
            fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
            (*outfile).write((char*)&dihedrals[i], sizeof(double)); //and the value of the dihedral itself
            fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        }
    }
    else if(double_prec==0) {//if single precision is requested do the analogue

        float dummy;
        dummy=root_origin_cartesian[0]; //conversion from double to float is done by assigning the double value to a float variable
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=root_origin_cartesian[1];
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=root_origin_cartesian[2];
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=root_origin_theta;
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=root_origin_phi;
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=root_origin_dihedral;
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=bonds[0];
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=bonds[1];
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        dummy=angles[0];
        (*outfile).write((char*)&dummy, sizeof(float));
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        for(int i=0; i<nDihedrals; i++) {
            dummy=bonds[i+2];
            (*outfile).write((char*)&dummy, sizeof(float));
            fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
            dummy=angles[i+1];
            (*outfile).write((char*)&dummy, sizeof(float));
            fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
            dummy=dihedrals[i];
            (*outfile).write((char*)&dummy, sizeof(float));
            fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        }
    }
    return fail;
}


int read_BAT_frame(ifstream *infile,int precision, int nDihedrals, float *time, float* xtcPrec, float **dbox, double* root_origin_cartesian,double* root_origin_theta,double* root_origin_phi,double* root_origin_dihedral,double *my_bonds,double *my_angles,double *my_dihedrals) {
    //to read a frame from the .bat trajectory
    float fdummy;
    int bCounter=2;
    int aCounter=1;
    int dCounter=0;
    int fail=0;

    (*infile).read((char*)time, sizeof(float));//read the time of the current according .xtc frame as a float
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)xtcPrec, sizeof(float));//read the precision of the current frame according .xtc frame for back conversion
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)dbox, 9*sizeof(float));//read the box vectors of the current according according .xtc frame for back conversion
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);

    if(precision==1) { //if double precision is used
        (*infile).read((char*)root_origin_cartesian, 3*sizeof(double));//read the Cartestians of the first root atom (external coordinates) in double precision
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)root_origin_theta, sizeof(double));//read the polar coordinates of the second root atom relative to the first (external coordinates)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)root_origin_phi, sizeof(double));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)root_origin_dihedral, sizeof(double));//and the dihedral the root atoms form with the origin (external coordinates)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);

        (*infile).read((char*)my_bonds, 2*sizeof(double));//read the lengths of the two bonds connecting the root atoms (internal coordinates)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)my_angles, sizeof(double));//and the angle between the two rootbonds (internal coordinates)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        for(int i=0; i<nDihedrals; i++) { //then for all dihedrals in the system
            (*infile).read((char*)&(my_bonds[bCounter]), sizeof(double));//read the bondlength between the last two atoms in the dihedral
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            bCounter++;
            (*infile).read((char*)&(my_angles[aCounter]), sizeof(double));//read the angle between the last threee atoms of the dihedral#
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            aCounter++;
            (*infile).read((char*)&(my_dihedrals[dCounter]), sizeof(double));//and the value of the dihedral itself
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            dCounter++;
        }
    }
    else if(precision==0) { //if single precision is used, do the same but use float instead of double
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        root_origin_cartesian[0]=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        root_origin_cartesian[1]=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        root_origin_cartesian[2]=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*root_origin_theta)=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*root_origin_phi)=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*root_origin_dihedral)=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        my_bonds[0]=fdummy;
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        my_bonds[1]=fdummy;
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        (*infile).read((char*)&fdummy, sizeof(float));
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        my_angles[0]=fdummy;
        for(int i=0; i<nDihedrals; i++) {
            (*infile).read((char*)&fdummy, sizeof(float));
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            my_bonds[bCounter]=fdummy;
            bCounter++;
            (*infile).read((char*)&fdummy, sizeof(float));
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            my_angles[aCounter]=fdummy;
            aCounter++;
            (*infile).read((char*)&fdummy, sizeof(float));
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            my_dihedrals[dCounter]=fdummy;
            dCounter++;
        }

    }
    return fail; //if anything failed return a 1, otherwise a 0
}


int write_PAR_header(ofstream *outfile,int nDihedrals,int double_prec,int numFrames,vector< vector <int> > *dihedrals_top, vector <float>  *masses, int bDens1D, int aDens1D, int dDens1D, int bDens, int aDens, int dDens, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule) {
    int dummy=(*dihedrals_top).size();
    int version=4;
    int fail=0;
    char dummystring[31];

    (*outfile).write((char*)&version, sizeof(int)); //first write the .par version number as an integer
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&double_prec, sizeof(int)); //then write an integer declaring if the trajectory was stored in double precision
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&dummy, sizeof(int)); //write an integer containing the number of dihedrals
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&numFrames, sizeof(int)); //and an integer containing the number of frames of the trajectory used for calculation
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&bDens1D, sizeof(int)); //write the used number of bins for 1D histograms
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&aDens1D, sizeof(int));
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&dDens1D, sizeof(int));
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&bDens, sizeof(int));//and the used number of bins for 2D histograms
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&aDens, sizeof(int));
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    (*outfile).write((char*)&dDens, sizeof(int));
    fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);

    for(int i=0; i<dummy+3; i++) { //for ever atom in the system
        for(int j=0; j<8; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.7s", (*residues)[i].c_str());
        (*outfile).write(dummystring, 8*sizeof(char));//write the name of the residue it belongs to
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*residueNumbers)[i]), sizeof(int));//write the number of the residue it belongs to
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        for(int j=0; j<8; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.7s", (*atomNames)[i].c_str());
        (*outfile).write(dummystring, 8*sizeof(char));//write the name of the atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        for(int j=0; j<31; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.30s", (*belongsToMolecule)[i].c_str());
        (*outfile).write(dummystring, 31*sizeof(char));//write the name of the molecule it belongs to
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    }


    for(int i=0; i<dummy; i++) {//then for all dihedrals
        (*outfile).write((char*)&((*dihedrals_top)[i][0]), sizeof(int)); //write the atomnumber of the first atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*dihedrals_top)[i][1]), sizeof(int)); //second atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*dihedrals_top)[i][2]), sizeof(int)); //third atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*dihedrals_top)[i][3]), sizeof(int)); //fourth atom
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*dihedrals_top)[i][4]), sizeof(int)); // an integer containing the type of the dihedral(physical=0,pseudo=1,improper=-1)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
        (*outfile).write((char*)&((*dihedrals_top)[i][5]), sizeof(int)); // and an integer containing the "parent"dihedral(for phaseangles, -1 if no "parent")
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);

    }
    for(int i=0; i<dummy+3; i++) {
        (*outfile).write((char*)&((*masses)[i]), sizeof(float)); //then write the whole massestor of the atoms in the system in single precision (float)
        fail=fail | ((*outfile).rdstate() & std::ofstream::failbit);
    }
    return fail; //if anything failed return a 1, otherwise a 0
}



int read_PAR_header(ifstream *infile,int *nDihedrals,int *double_prec,int *numFrames,vector< vector <int> > *dihedrals_top, vector <float>  *masses, int *version, int* bDens, int* aDens, int* dDens, int* bDens1D, int* aDens1D, int* dDens1D, vector <string> *residues,vector <int> *residueNumbers,vector <string> *atomNames,vector <string> *belongsToMolecule) {
    int fail=0;
    char dummystring[31];


    (*infile).read((char*)version, sizeof(int)); //first read the .par version number as an integer
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)double_prec, sizeof(int));//then read an integer declaring if the trajectory was stored in double precision
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)nDihedrals,sizeof(int));//read an integer containing the number of dihedrals
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)numFrames,sizeof(int));//and an integer containing the number of frames of the trajectory used for calculation
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)bDens1D,sizeof(int)); //read the used number of bins for 1D histograms
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)aDens1D,sizeof(int));
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)dDens1D,sizeof(int));
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)bDens,sizeof(int)); //and the used number of bins for 2D histograms
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)aDens,sizeof(int));
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    (*infile).read((char*)dDens,sizeof(int));
    fail=fail | ((*infile).rdstate() & std::ifstream::failbit);


    //check for errors
    if(fail!=0) {
        return fail;
    }
    if((*version)<0) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. VERSION NUMBER ("<<(*version)<<") < 0!"<<endl;
        return 1;
    }
    if((*double_prec!=0)&&(*double_prec!=1)) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. DOUBLE PRECISION VALUE ("<<(*double_prec)<<") NEITHER 0 NOR 1!"<<endl;
        return 1;
    }
    if((*nDihedrals)>19997) {
        cerr<<"WARNING: "<<nDihedrals+3<<" ATOMS DECLARED IN THE FILE HEADER (CORRUPTED?). THIS WILL LEAD TO LARGE OUTPUT."<<endl;
    }
    if((*nDihedrals)<0) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER OF DIHEDRALS ("<<(*nDihedrals)<<") < 0!"<<endl;
        return 1;
    }
    if((*numFrames)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER OF FRAMES ("<<(*numFrames)<<") < 1!"<<endl;
        return 1;
    }
    if((*bDens1D)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR BONDS IN 1D ("<<(*bDens1D)<<") < 1!"<<endl;
        return 1;
    }
    if((*aDens1D)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR ANGLES IN 1D ("<<(*aDens1D)<<") < 1!"<<endl;
        return 1;
    }
    if((*dDens1D)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR DIHEDRALS IN 1D ("<<(*dDens1D)<<") < 1!"<<endl;
        return 1;
    }
    if((*bDens)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR BONDS IN 2D ("<<(*bDens)<<") < 1!"<<endl;
        return 1;
    }
    if((*aDens)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR ANGLES IN 2D ("<<(*aDens)<<") < 1!"<<endl;
        return 1;
    }
    if((*dDens)<1) {
        cerr<<"ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR DIHEDRALS IN 2D ("<<(*dDens)<<") < 1!"<<endl;
        return 1;
    }
		
		
		
		if ((*version)>=3) {
        for(int i=0; i<(*nDihedrals)+3; i++) { //for every atom in the system
            (*infile).read(dummystring, 8*sizeof(char));//read the name of the residue it belongs to
            (*residues).push_back(dummystring);
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            (*residueNumbers).push_back(0);//the number of the residue it belongs to
            (*infile).read((char*)&((*residueNumbers)[i]), sizeof(float));
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            (*infile).read(dummystring, 8*sizeof(char));//the name of the atom
            (*atomNames).push_back(dummystring);
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
            (*infile).read(dummystring, 31*sizeof(char));
            (*belongsToMolecule).push_back(dummystring); //and the name of the molecule it belongs to
            fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        }
    }
		
		
		    //check for errors
    if(fail!=0) {
        return fail;
    }


    vector<int>dummyvec;
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    for(int i=0; i<(*nDihedrals); i++) { //then for all dihedrals
        (*dihedrals_top).push_back(dummyvec);
        (*infile).read((char*)&((*dihedrals_top)[i][0]),sizeof(int));//write the atomnumber of the first atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        if((*dihedrals_top)[i][0]<1) {
            cerr<<"ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !"<<endl;
            return 1;
        }
        (*infile).read((char*)&((*dihedrals_top)[i][1]),sizeof(int));//second atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        if((*dihedrals_top)[i][1]<1) {
            cerr<<"ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !"<<endl;
            return 1;
        }
        (*infile).read((char*)&((*dihedrals_top)[i][2]),sizeof(int));//third atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        if((*dihedrals_top)[i][2]<1) {
            cerr<<"ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !"<<endl;
            return 1;
        }
        (*infile).read((char*)&((*dihedrals_top)[i][3]),sizeof(int));//fourth atom
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        if((*dihedrals_top)[i][3]<1) {
            cerr<<"ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !"<<endl;
            return 1;
        }
        (*infile).read((char*)&((*dihedrals_top)[i][4]),sizeof(int));// an integer containing the type of the dihedral(physical=0,pseudo=1,improper=-1)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        if(!(((*dihedrals_top)[i][4]==1)||((*dihedrals_top)[i][4]==0)||((*dihedrals_top)[i][4]==-1))) {
            cerr<<"ERROR: FILE HEADER CORRUPTED. DIHEDRAL TYPE NEITHER 0, 1 NOR -1 !"<<endl;
            return 1;
        }
        (*infile).read((char*)&((*dihedrals_top)[i][5]),sizeof(int));// and an integer containing the "parent" dihedral(for phaseangles, -1 if no "parent")
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
        if((*dihedrals_top)[i][5]<-1) {
            cerr<<"ERROR: FILE HEADER CORRUPTED. DIHEDRAL IS PHASEANGLE OF A DIHEDRAL WITH NEGATIVE ID !"<<endl;
            return 1;
        }

        if((*version==1)&&((*dihedrals_top)[i][5]!=-1)) {
            (*dihedrals_top)[i][5]++; //for backwards compatibility
        }
    }
    for(int i=0; i<(*nDihedrals)+3; i++) {
        (*masses).push_back(0);
        (*infile).read((char*)&((*masses)[i]), sizeof(float));//then read the whole massestor of the atoms in the system in single precision (float)
        fail=fail | ((*infile).rdstate() & std::ifstream::failbit);
    }
    return fail; //if anything failed return a 1, otherwise a 0

}

int write_PAR_body(ofstream* par_file, int nDihedrals,double* bondsEntropy1D, double* anglesEntropy1D, double* dihedralsEntropy1D, double* bbEntropy, double* baEntropy, double* bdEntropy, double* aaEntropy, double* adEntropy, double* ddEntropy) {
    int fail=0;

    int nBonds=nDihedrals+2;
    int nAngles=nDihedrals+1;

    (*par_file).write((char*)bondsEntropy1D, nBonds*sizeof(double)); //write the 1D entropy bonds array
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)anglesEntropy1D, nAngles*sizeof(double)); //write the 1D entropy angles array
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)dihedralsEntropy1D, nDihedrals*sizeof(double)); //write the 1D entropy dihedrals array
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)bbEntropy, nBonds*(nBonds-1)/2*sizeof(double)); //write the 2D bonds-bonds half-matrix as an array
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)baEntropy, nBonds*nAngles*sizeof(double)); //write the 2D bonds-angles matrix
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)bdEntropy, nBonds*nDihedrals*sizeof(double)); //write the 2D bonds-dihedrals matrix
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)aaEntropy, nAngles*(nAngles-1)/2*sizeof(double)); //write the 2D angles-angles half-matrix as an array
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)adEntropy, nAngles*nDihedrals*sizeof(double)); //write the 2D angles-dihedrals matrix
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    (*par_file).write((char*)ddEntropy, nDihedrals*(nDihedrals-1)/2*sizeof(double)); //write the 2D dihedrlas-dihedrals half-matrix as an array
    fail=fail | ((*par_file).rdstate() & std::ofstream::failbit);
    return fail;
}

int read_PAR_body(ifstream* par_file, int nDihedrals,double** bondsEntropy1D, double** anglesEntropy1D, double** dihedralsEntropy1D, double** bbEntropy, double** baEntropy, double** bdEntropy, double** aaEntropy, double** adEntropy, double** ddEntropy) {
    int fail=0;

    int nBonds=nDihedrals+2;
    int nAngles=nDihedrals+1;

    (*bondsEntropy1D)=(double*)malloc(nBonds*sizeof(double)); // node allocates storage for reading the .par file
    (*anglesEntropy1D)=(double*)malloc(nAngles*sizeof(double));
    (*dihedralsEntropy1D)=(double*)malloc(nDihedrals*sizeof(double));
    (*bbEntropy)=(double*)malloc(nBonds*(nBonds-1)/2*sizeof(double));
    (*baEntropy)=(double*)malloc(nBonds*nAngles*sizeof(double));
    (*bdEntropy)=(double*)malloc(nBonds*nDihedrals*sizeof(double));
    (*aaEntropy)=(double*)malloc(nAngles*(nAngles-1)/2*sizeof(double));
    (*adEntropy)=(double*)malloc(nAngles*nDihedrals*sizeof(double));
    (*ddEntropy)=(double*)malloc(nDihedrals*(nDihedrals-1)/2*sizeof(double));

    if(!((*bondsEntropy1D!=NULL)&&(*anglesEntropy1D!=NULL)&&(*dihedralsEntropy1D!=NULL)&&(*bbEntropy!=NULL)&&(*baEntropy!=NULL)&&(*bdEntropy!=NULL)&&(*aaEntropy!=NULL)&&(*adEntropy!=NULL)&&(*ddEntropy!=NULL))) {
        cerr<<"ERROR: ALLOCATION FOR PAR FILE BODY FAILED. MORE MEMORY NEEDED?"<<endl;
        return 2;
    }

    (*par_file).read((char*)*bondsEntropy1D, nBonds*sizeof(double)); // and the actually reads it in (see write_PAR_body for details)
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*anglesEntropy1D, nAngles*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*dihedralsEntropy1D, nDihedrals*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*bbEntropy, nBonds*(nBonds-1)/2*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*baEntropy, nBonds*nAngles*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*bdEntropy, nBonds*nDihedrals*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*aaEntropy, nAngles*(nAngles-1)/2*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*adEntropy, nAngles*nDihedrals*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    (*par_file).read((char*)*ddEntropy, nDihedrals*(nDihedrals-1)/2*sizeof(double));
    fail=fail | ((*par_file).rdstate() & std::ifstream::failbit);
    return fail;
}















