//    BAT_builder, a program to convert a molecular dynamics trajectory from Cartesian to internal bond-angle-torsion coordinates
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







#include "bat.h"

int main(int argc, char* argv[])
{
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"    BAT_builder, a program to convert a molecular dynamics trajectory from Cartesian to internal bond-angle-torsion coordinates"<<endl;
    cout<<"    Copyright (C) 2015  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)"<<endl;
    cout<<endl;
    cout<<"    This program is free software: you can redistribute it and/or modify"<<endl;
    cout<<"    it under the terms of the GNU General Public License  version 3"<<endl;
    cout<<"    as published by the Free Software Foundation."<<endl;
    cout<<endl;
    cout<<"    This program is distributed in the hope that it will be useful,"<<endl;
    cout<<"    but WITHOUT ANY WARRANTY; without even the implied warranty of"<<endl;
    cout<<"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"<<endl;
    cout<<"    GNU General Public License for more details."<<endl;
    cout<<endl;
    cout<<"    You should have received a copy of the GNU General Public License"<<endl;
    cout<<"    along with this program.  If not, see <http://www.gnu.org/licenses/>."<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<"    A scientific publication about this program has been released in the Journal of Chemical Theory and Computation:"<<endl;
    cout<<"    \"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems\""<<endl;
    cout<<"    DOI: 10.1021/acs.jctc.5b01217"<<endl;
    cout<<endl;
    cout<<"    Please include this citation in works that publish results generated using"<<endl;
    cout<<"    this program or any modifications of this program."<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;


    char delimiter[] = ".";
    char *ptr,*type1,*type2,*type3;


    if(argc==5) { //conversion from Cartesians to BAT in single precision

        string tmp1(argv[1]);
        string tmp2(argv[2]);
        string tmp3(argv[3]);
        vector <string> backboneAtomNames;

        //extract the filename extansions from the command line arguments
        ptr = strtok((char*)tmp1.c_str(), delimiter);
        while(ptr != NULL) {
            type1=ptr;
            ptr = strtok(NULL, delimiter);
        }
        ptr = strtok((char*)tmp2.c_str(), delimiter);
        while(ptr != NULL) {
            type2=ptr;
            ptr = strtok(NULL, delimiter);
        }
        ptr = strtok((char*)tmp3.c_str(), delimiter);
        while(ptr != NULL) {
            type3=ptr;
            ptr = strtok(NULL, delimiter);
        }

        //check if command line arguments were supplied correctly, else put out help and quit
        if(((strcmp(type1,"top"))||(strcmp(type2,"xtc"))||(strcmp(type3,"bat")))) {
            cerr<<"USAGE: "<<argv[0]<<" input.top input.xtc output.bat \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [double_precision]\nOR "<<argv[0]<<" input.bat output.xtc\n";
            return 1;
        }

        //Read in the requested names of backbone atoms
        ptr = strtok (argv[4]," ");
        while (ptr != NULL)
        {
            backboneAtomNames.push_back(string(ptr));
            ptr = strtok (NULL, " ");
        }



        //Read the topology file and build a table for the molecules bonds as well as a list of backbone atoms
        cout<<"Reading topology."<<endl;
        Topology proteins(string(argv[1]).c_str(), backboneAtomNames);

        //Build the BAT Topology using the molecules bonds (and the backbone atoms to make use of phase angles)
        cout<<"Building BAT topology."<<endl;
        BAT_Topology bat(&proteins.bonds_table,&proteins.backbone,&proteins.roots);

        //Use the list of dihedrals and the original trajectory file to build the BAT binary trajectory in single precision (masses are included in the .bat file since they might come in handy at a later point). Also add information about the atoms.
        cout<<"Writing .bat trajectory."<<endl;
        BAT_Trajectory trj(string(argv[2]).c_str(),string(argv[3]).c_str(), &bat.dihedrals, &proteins.masses,&proteins.residues,&proteins.residueNumbers,&proteins.atomNames,&proteins.belongsToMolecule);


    }
    else if(argc==6) { //conversion from Cartesians to BAT in double precision

        string tmp1(argv[1]);
        string tmp2(argv[2]);
        string tmp3(argv[3]);
        string tmp4(argv[5]);
        vector <string> backboneAtomNames;

        //extract the filename extansions from the command line arguments
        ptr = strtok((char*)tmp1.c_str(), delimiter);
        while(ptr != NULL) {
            type1=ptr;
            ptr = strtok(NULL, delimiter);
        }
        ptr = strtok((char*)tmp2.c_str(), delimiter);
        while(ptr != NULL) {
            type2=ptr;
            ptr = strtok(NULL, delimiter);
        }
        ptr = strtok((char*)tmp3.c_str(), delimiter);
        while(ptr != NULL) {
            type3=ptr;
            ptr = strtok(NULL, delimiter);
        }

        //check if command line arguments were supplied correctly, else put out help and quit
        if(((strcmp(type1,"top"))||(strcmp(type2,"xtc"))||(strcmp(type3,"bat"))||(strcmp(tmp4.c_str(),"double_precision")))) {
            cerr<<"USAGE: "<<argv[0]<<" input.top input.xtc output.bat \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [double_precision]\nOR "<<argv[0]<<" input.bat output.xtc\n";
            return 1;
        }

        //Read in the requested names of backbone atoms
        ptr = strtok (argv[4]," ");
        while (ptr != NULL)
        {
            backboneAtomNames.push_back(string(ptr));
            ptr = strtok (NULL, " ");
        }


        //Read the topology file and build a table for the molecules bonds as well as a list of backbone atoms
        cout<<"Reading topology."<<endl;
        Topology proteins(string(argv[1]).c_str(), backboneAtomNames);

        //Build the BAT Topology using the molecules bonds (and the backbone atoms to make use of phase angles)
        cout<<"Building BAT topology."<<endl;
        BAT_Topology bat(&proteins.bonds_table,&proteins.backbone,&proteins.roots);

        //Use the list of dihedrals and the original trajectory file to build the BAT binary trajectory in double precision (masses are included in the .bat file since they might come in handy at a later point). Also add information about the atoms.
        cout<<"Writing .bat trajectory."<<endl;
        BAT_Trajectory trj(string(argv[2]).c_str(),string(argv[3]).c_str(), &bat.dihedrals,&proteins.masses,&proteins.residues,&proteins.residueNumbers,&proteins.atomNames,&proteins.belongsToMolecule,1);

    }
    else  if (argc==3) { // conversion from BAT to Cartesians

        string tmp1(argv[1]);
        string tmp2(argv[2]);

        //extract the filename extansions from the command line arguments
        ptr = strtok((char*)tmp1.c_str(), delimiter);
        while(ptr != NULL) {
            type1=ptr;
            ptr = strtok(NULL, delimiter);
        }
        ptr = strtok((char*)tmp2.c_str(), delimiter);
        while(ptr != NULL) {
            type2=ptr;
            ptr = strtok(NULL, delimiter);
        }

        //check if command line arguments were supplied correctly, else put out help and quit
        if(((strcmp(type1,"bat"))||(strcmp(type2,"xtc")))) {
            cerr<<"USAGE: "<<argv[0]<<" input.top input.xtc output.bat \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [double_precision]\nOR "<<argv[0]<<" input.bat output.xtc\n";
            return 1;
        }


        //do the conversion
        cout<<"Converting .bat to .xtc."<<endl;
        BAT_Trajectory trj(string(argv[1]).c_str(),string(argv[2]).c_str());

    }
    else {
        cerr<<"USAGE: "<<argv[0]<<" input.top input.xtc output.bat \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [double_precision]\nOR "<<argv[0]<<" input.bat output.xtc\n";
        return 1;
    }


    cout<<"Finished writing trajectory."<<endl;
    exit(EXIT_SUCCESS);
}















