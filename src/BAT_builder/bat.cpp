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







#include "topology.h"
#include "BAT_topology.h"
#include "BAT_trajectory.h"
#include "../util/util.h"

#include <iostream>
#include <cstdlib>
#include <cstring>


using namespace std;

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


    if(argc==9){ //conversion from Cartesians to BAT in double precision
					if(!cmdOptionExists(argv, argv+argc, "-t")||!cmdOptionExists(argv, argv+argc, "-x")||!cmdOptionExists(argv, argv+argc, "-o")||!cmdOptionExists(argv, argv+argc, "-bb")){
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
						return 1;
					}
        
				std::string tmp1(getCmdOption(argv, argv+argc, "-t"));
        std::string tmp2(getCmdOption(argv, argv+argc, "-x"));
        std::string tmp3(getCmdOption(argv, argv+argc, "-o"));
        vector <std::string> backboneAtomNames;

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
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
            return 1;
        }

        //Read in the requested names of backbone atoms
        ptr = strtok (getCmdOption(argv, argv+argc, "-bb")," ");
        while (ptr != NULL)
        {
            backboneAtomNames.push_back(std::string(ptr));
            ptr = strtok (NULL, " ");
        }

				//Read the topology file and build a table for the molecules bonds as well as a list of backbone atoms
        cout<<"Reading topology."<<endl;
        Topology proteins(std::string(getCmdOption(argv, argv+argc, "-t")).c_str(), backboneAtomNames);

        //Build the BAT Topology using the molecules bonds (and the backbone atoms to make use of phase angles)
        cout<<"Building BAT topology."<<endl;
        BAT_Topology bat(&proteins.bonds_table,&proteins.backbone,&proteins.roots);

        //Use the list of dihedrals and the original trajectory file to build the BAT binary trajectory in double precision (masses are included in the .bat file since they might come in handy at a later point). Also add information about the atoms.
        cout<<"Writing .bat trajectory."<<endl;
        BAT_Trajectory trj(std::string(getCmdOption(argv, argv+argc, "-x")).c_str(),std::string(getCmdOption(argv, argv+argc, "-o")).c_str(), &bat.dihedrals, &proteins.masses,&proteins.residues,&proteins.residueNumbers,&proteins.atomNames,&proteins.belongsToMolecule);
		
		}
    else if(argc==10) { //conversion from Cartesians to BAT in single precision
					if(!cmdOptionExists(argv, argv+argc, "-t")||!cmdOptionExists(argv, argv+argc, "-x")||!cmdOptionExists(argv, argv+argc, "-o")||!cmdOptionExists(argv, argv+argc, "-bb")||!cmdOptionExists(argv, argv+argc, "--single_precision")){
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
						return 1;
					}
					
				std::string tmp1(getCmdOption(argv, argv+argc, "-t"));
        std::string tmp2(getCmdOption(argv, argv+argc, "-x"));
        std::string tmp3(getCmdOption(argv, argv+argc, "-o"));
        vector <std::string> backboneAtomNames;

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
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
            return 1;
        }

        //Read in the requested names of backbone atoms
        ptr = strtok (getCmdOption(argv, argv+argc, "-bb")," ");
        while (ptr != NULL)
        {
            backboneAtomNames.push_back(std::string(ptr));
            ptr = strtok (NULL, " ");
        }


        //Read the topology file and build a table for the molecules bonds as well as a list of backbone atoms
        cout<<"Reading topology."<<endl;
        Topology proteins(std::string(getCmdOption(argv, argv+argc, "-t")).c_str(), backboneAtomNames);

        //Build the BAT Topology using the molecules bonds (and the backbone atoms to make use of phase angles)
        cout<<"Building BAT topology."<<endl;
        BAT_Topology bat(&proteins.bonds_table,&proteins.backbone,&proteins.roots);

        //Use the list of dihedrals and the original trajectory file to build the BAT binary trajectory in single precision (masses are included in the .bat file since they might come in handy at a later point). Also add information about the atoms.
        cout<<"Writing .bat trajectory."<<endl;
        BAT_Trajectory trj(std::string(getCmdOption(argv, argv+argc, "-x")).c_str(),std::string(getCmdOption(argv, argv+argc, "-o")).c_str(), &bat.dihedrals,&proteins.masses,&proteins.residues,&proteins.residueNumbers,&proteins.atomNames,&proteins.belongsToMolecule,0);

    }
    else  if (argc==5) { // conversion from BAT to Cartesians
					if(!cmdOptionExists(argv, argv+argc, "-b")||!cmdOptionExists(argv, argv+argc, "-o")){
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
						return 1;
					}
					
				std::string tmp1(getCmdOption(argv, argv+argc, "-b"));
        std::string tmp2(getCmdOption(argv, argv+argc, "-o"));


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
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
            return 1;
        }


        //do the conversion
        cout<<"Converting .bat to .xtc."<<endl;
        BAT_Trajectory trj(std::string(getCmdOption(argv, argv+argc, "-b")).c_str(),std::string(getCmdOption(argv, argv+argc, "-o")).c_str());

    }
    else {
            cerr<<"USAGE: "<<argv[0]<<" -t input.top -x input.xtc -o output.bat -bb \"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ...\" [--single_precision]\nOR "<<argv[0]<<" -b input.bat -o output.xtc\n";
        return 1;
    }


    cout<<"Finished writing trajectory."<<endl;
    exit(EXIT_SUCCESS);
}















