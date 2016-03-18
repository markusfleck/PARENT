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







#include "../util/io/io.h"



using namespace std;

class Topology {

public:
    Topology(string topfileIN,vector <string> backboneAtomNamesIN) {
        current=-1;
        newclustflag=0;
        newclustmassflag=0;
        newclustatomflag=0;
        newclustbackboneflag=0;
        top_file=topfileIN;
        backboneAtomNames=backboneAtomNamesIN;

        //get the directory of the topology file
        top_dir=get_dir_file(topfileIN);

        //start the topology building
        read_topology();

    }

    ~Topology() {

    }

    vector< vector<int> > moleculetypes;  //contains a list of the atoms in each moleculetype
    vector< vector<int> > molecules; //contains a list of atoms for each molecule in the system with correct offset, as defined in the [ molecules ] section
    vector< vector<int> > angles; //contains a list of angles for each molecule in the system with correct offset, as defined in the [ molecules ] section
    vector< vector< vector <int> > > bonds_table; //contains a list of all bonds in the system with proper offset, with bonds_table[#atom-1][#bond-1][0] the bonded atom and bonds_table[atom-1][#bond-1][1] the type of the bond (0==physical, 1==pseudo)
    vector< vector< vector <int> > > moleculetypebonds; //contains a list  of all bonds for every moleculetype
    vector< vector< vector <int> > > bonds;  //contains a list of atoms for each molecule in the system with correct offset, as defined in the [ molecules ] section
    vector<string> moleculetypenames;  //names of the moleculetypes
    vector<string> moleculenames; //contains a full list of molecule names of the system as defined in the [ molecules ] section (water 1000 => water, water, water.......)
    vector< vector <string> > moleculetypeResidues; //contains the name of the residue for every moleculetype and every atom
    vector <string> residues; //contains the name of the residue for every atom in the system
    vector< vector <int> > moleculetypeResidueNumbers; //contains the name of the residue number for every moleculetype and every atom
    vector <int> residueNumbers; //contains the number of the residue for every atom in the system
    vector< vector <string> > moleculetypeAtomNames; //contains the name of the atom for every moleculetype and every atom
    vector <string> atomNames; //contains the name of every atom in the system
    vector <string> belongsToMolecule; //the name of the molecule each atom belongs to
    vector< vector<float> >  moleculetypemasses; //masses for each moleculetype
    vector< vector<int> >  moleculetypebackbone; //backboneatoms for each moleculetype
    vector<float> masses; //diagonal of the mass matrix of the system
    vector <int> backbone; //backbone atoms in the system
    vector <string> backboneAtomNames; //names of the atomtypes of the backbone for building phase angles from dihedrals, which is done in BAT_trajectory.cpp
    vector< vector <int> > roots; //root atoms for all molecules in the system. root atoms need to be connected root[0] to root[1] and root[1] to root[2]. root[0] needs to be terminal, root[1] needs to be connected only to root[1] and terminal atoms. root[2] needs to be connected to non-terminal atoms.

    long int current;
    long int total_atoms;
    unsigned short int newclustflag,newclustmassflag,newclustatomflag,newclustbackboneflag;
    string top_dir, top_file;

    unsigned short int add_bond(int atom1, int atom2) {
        //for a new [ moleculetype ], fill in the (preallocated) first 2 atoms to the current "moleculetypes" vector and the (preallocated) first bond to the current "moleculetypebonds" vector
        if(newclustflag==1) {
            moleculetypes[current][0]=atom1;
            moleculetypes[current][1]=atom2;
            moleculetypebonds[current][0][0]=atom1;
            moleculetypebonds[current][0][1]=atom2;
            newclustflag=0;
        }
        //else add the bond to the current  "moleculetypebonds" vector
        else {
            vector<int> tmpvec;
            tmpvec.push_back(atom1);
            tmpvec.push_back(atom2);
            moleculetypebonds[current].push_back(tmpvec);

            //and add the 2 atoms to the moleculetypes vector, but only if they are not already present
            int tocluster=1;
            for(unsigned int j=0; j < moleculetypes[current].size(); j++) {
                if(atom1==moleculetypes[current][j]) {
                    tocluster=0;
                    break;
                }
            }
            if(tocluster==1) {
                moleculetypes[current].push_back(atom1);
            }
            tocluster=1;
            for(unsigned int j=0; j < moleculetypes[current].size(); j++) {
                if(atom2==moleculetypes[current][j]) {
                    tocluster=0;
                    break;
                }
            }
            if(tocluster==1) {
                moleculetypes[current].push_back(atom2);
            }
        }
        return 0;
    }



    unsigned short int add_mass(float mass) {
        //for a new [ moleculetype ], create a new massvector filled with the first atommass and add it to the "moleculetypemasses" vector
        if(newclustmassflag==1) {
            vector<float> tmpvec;
            tmpvec.push_back(mass);
            moleculetypemasses.push_back(tmpvec);
            newclustmassflag=0;
        }
        //else just add the mass to the current "moleculetypemasses" vector
        else {
            moleculetypemasses[current].push_back(mass);
        }
        return 0;
    }


    unsigned short int add_atom(string resname, int resnumber, string atomname) {
        //for a new [ moleculetype ], create new vectors filled with the first values and add them to the moleculetype vectors
        if(newclustatomflag==1) {
            vector<string> tmpvec1;
            vector<int> tmpvec2;
            vector<string> tmpvec3;
            tmpvec1.push_back(resname);
            tmpvec2.push_back(resnumber);
            tmpvec3.push_back(atomname);
            moleculetypeResidues.push_back(tmpvec1);
            moleculetypeResidueNumbers.push_back(tmpvec2);
            moleculetypeAtomNames.push_back(tmpvec3);
            newclustatomflag=0;
        }
        //else just add the values to the current moleculetype vectors
        else {
            moleculetypeResidues[current].push_back(resname);
            moleculetypeResidueNumbers[current].push_back(resnumber);
            moleculetypeAtomNames[current].push_back(atomname);
        }
        return 0;
    }



    unsigned short int add_backbone(int atom) {
        //for a new [ moleculetype ], fill in the (preallocated) atom to the current "moleculetypebackbone" vector
        if(newclustbackboneflag==1) {
            moleculetypebackbone[current][0]=atom;
            newclustbackboneflag=0;
        }
        //else add the atom to the current  "moleculetypebbackbone" vector
        else {
            moleculetypebackbone[current].push_back(atom);
        }
        return 0;
    }



    void new_moleculetype(string mol) {

        //set up the necessary flags and prolong the vectors (initialized with dummy values) for taking up a new molecule. Also increase the "current" counter;
        moleculetypenames.push_back(mol);
        current++;
        newclustflag=1;
        newclustmassflag=1;
        newclustatomflag=1;

        newclustbackboneflag=1;
        vector<int> tmpvec;
        tmpvec.push_back(-10);
        moleculetypebackbone.push_back(tmpvec);


        tmpvec.clear();
        vector< vector<int> > tmpvec1;
        tmpvec.push_back(-10);
        tmpvec.push_back(-10);
        tmpvec1.push_back(tmpvec);
        moleculetypes.push_back(tmpvec);
        moleculetypebonds.push_back(tmpvec1);
    }



    int sort() {
        int offset=0;
        unsigned short int found=0;
        vector<int> clusterout;


        //according to the [ molecules ] section, join the necessary arrays accordingly
        for(unsigned int i=0; i<moleculenames.size(); i++) {
            found=0;
            for(unsigned int j=0; j<moleculetypes.size(); j++) {
                if(moleculenames[i]==moleculetypenames[j]) { //find the according moleculetypename
                    clusterout.clear();
                    clusterout=moleculetypes[j];
                    for(unsigned int k=0; k<(clusterout).size(); k++) {
                        clusterout[k]+=offset;
                    }
                    molecules.push_back(clusterout); //add the according "moleculetypes" vector to the "molecules" matrix (considering the offset due to previous molecules in the system)
                    for(unsigned int k=0; k<moleculetypemasses[j].size(); k++) {
                        masses.push_back(moleculetypemasses[j][k]); //attach the according "moleculetypemasses" vector to the "masses" vector
                    }
                    for(unsigned int k=0; k<moleculetypebackbone[j].size(); k++) {  //attach the according "moleculetypebackbone" vector to the "backbone" vector (considering the offset due to previous molecules in the system)
                        backbone.push_back(moleculetypebackbone[j][k]+offset);
                    }

                    bonds.push_back(moleculetypebonds[j]); //add the according "moleculetypebonds" matrix of the moleculetype to the "bonds" 3D-tensor (considering the offset due to previous molecules in the system)
                    for(unsigned int k=0; k<bonds[i].size(); k++) {
                        bonds[i][k][0]+=offset;
                        bonds[i][k][1]+=offset;
                    }
                    for(unsigned int k=0; k<moleculetypeResidues[j].size(); k++) { //add informations about the atom to the respective vectors
                        residues.push_back(moleculetypeResidues[j][k]);
                        residueNumbers.push_back(moleculetypeResidueNumbers[j][k]);
                        atomNames.push_back(moleculetypeAtomNames[j][k]);
                        belongsToMolecule.push_back(moleculenames[i]);
                    }


                    offset+=moleculetypes[j].size(); //increase the offset by the number of atoms in the just processed moleculetype
                    found=1;
                    break;
                }
            }
            if(found==0) {
                cerr<<"ERROR: MOLECULETYPE \""<<moleculenames[i]<<"\" DOES NOT EXIST!"<<endl;
                exit(EXIT_FAILURE);
            }
            total_atoms=offset;
        }
        return 0;
    }



    string get_dir_file(string file) {
        //returns the directory in a path to a file
        size_t found = file.find('/');
        if (found!=string::npos) {
            size_t found_old;
            while(found!=string::npos) {
                found_old=found;
                found = file.find('/',found+1);
            }
            file.erase(file.begin()+found_old+1, file.end());
            return file;
        }
        else {
            return string("./");
        }
    }



    unsigned short int check_comment(string line) {
        //search for the first non-blank character in the line, if it is ";" return 1, otherwise return 0. If all characters are blanks, return 2.
        for(unsigned int i=0; i<line.length(); i++) {
            if(line[i]!=' ') {
                if(line[i]==';') {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        return 2;
    }



    int preprocess(string infile,stringstream* myVirtualOFile) {
        string line,tmpstr="",tmpstr2="",mytop_dir=get_dir_file(infile);
        char c1[1000],c2[1000];
        int found;

        //open input and output files
        ifstream myinfile(infile.c_str());
        if (!myinfile.is_open()) {
            cerr<<"ERROR: COULD NOT OPEN FILE \""<<infile<<"\"\n";
            exit(EXIT_FAILURE);
        }



        //(pre)process every line of the input file
        while (getline(myinfile,line)) {

            //remove preceding and concluding blanks
            line=strip_line(line);

            //scan the first two columns
            c1[0]='\0';
            c2[0]='\0';
            sscanf (line.c_str(),"%s %s",c1,c2);

            //check if the line is an include directive
            if(string(c1)=="#include") {
                found=0;

                //remove the quotation marks and check if the file is present in the same directory as the topology file
                tmpstr=mytop_dir;
                tmpstr2=string(c2);
                tmpstr2.erase(tmpstr2.begin());
                tmpstr2.erase(tmpstr2.end()-1);
                tmpstr+=tmpstr2;
                ifstream tmpfile(tmpstr.c_str());

                //if the file is there, recursively preprocess it
                if(tmpfile.is_open()&&line!="") {
                    found=1;
                    tmpfile.close();
                    preprocess(tmpstr,myVirtualOFile);
                }

                //check if the file is present at the default gromacs topology location
                if(found!=1) {
                    tmpstr="/usr/local/gromacs/share/gromacs/top/"+tmpstr2;
                    tmpfile.open(tmpstr.c_str());
                    //if the file is there, recursively preprocess it
                    if(tmpfile.is_open()&&line!="") {
                        found=1;
                        tmpfile.close();
                        preprocess(tmpstr,myVirtualOFile);
                    }
                }

                if(found!=1) {
                    //if the file wasn't found issue a warning
                    if(check_comment(line)==0) {
                        (*myVirtualOFile)<<line<<endl;
                        cerr<<"WARNING FROM PREPROCESSOR: file "<<tmpstr2<<"   "<<(found==1 ? " found ":string(" NOT FOUND locally or at ")+tmpstr)<<". If the file does not contain any bond-structure information of your molecule(s), this should be fine.\n";
                    }
                }

            }
            else {
                //if the line is not an include directive, a comment or a blank line, just write the stripped line to myVirtualOfile
                if(check_comment(line)==0&&c1[0]!='\0') {
                    (*myVirtualOFile)<<line<<endl;
                }
            }
        }

        myinfile.close();
        return 0;
    }



    int read_topology() {

        //set all [ section ] flags to false;
        unsigned short int atomsflag=0;
        unsigned short int bondsflag=0;
        unsigned short int moleculetypeflag=0;
        unsigned short int moleculesflag=0;

        unsigned short int hasAtoms=1;  //set the number of [ atoms ] and [ bonds ] sections for the 0th (non-existent) [ moleculetype ] correctly to 1 (if a [ molecultype ] or a [ molecules ] section is encountered, to number of [ atoms ] and [ bonds ] sections for the previous [ moleculetype ] is checked to be 1)
        unsigned short int hasBonds=1;

        int n1,n2,n3;
        char c1[1000];
        char c2[1000];
        char c3[1000];
        float charge,mass;
        bool isBackbone;


        //preprocess the topology file to write a virtual file myVirtualFile
        stringstream myVirtualFile(string(""));
        preprocess(top_file,&myVirtualFile);
        string line, line_in;

        while ( getline(myVirtualFile,line_in) ) {
            line=strip_line(line_in);
            if(moleculetypeflag==1&&line[0]!='['&&line[0]!='#') {
                if(check_comment(line)==0) {  //if a new [ moleculetype ] section is encountered,
                    sscanf (line.c_str(),"%s %d",c1,&n2);
                    new_moleculetype(string(c1));//create that moleculetype
                    //and check if the previous [ moleculetype ] section had exactly one [ atoms ] and [ bonds ] section
                    if(hasAtoms==0) {
                        cerr<<"ERROR: NO \"[ atoms ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current-1]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                        exit(EXIT_FAILURE);
                    }
                    if((hasBonds==0)&&(moleculetypemasses[current-1].size()>1)) {
                        cerr<<"WARNING: NO \"[ bonds ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current-1]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                    }
                    if(hasAtoms>1) {
                        cerr<<"ERROR: SEVERAL ("<<hasAtoms<<") \"[ atoms ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current-1]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                        exit(EXIT_FAILURE);
                    }
                    if(hasBonds>1) {
                        cerr<<"ERROR: SEVERAL ("<<hasBonds<<") \"[ bonds ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current-1]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                        exit(EXIT_FAILURE);
                    }
                    //set the [ atoms ] and [ bonds ] section counter for the current [ moleculetype ] to 0
                    hasAtoms=0;
                    hasBonds=0;
                }
            }


            //if the new line is in a [ bonds ] section add the bond to the current moleculetype (also check for the existence of a preceding [ moleculetype ] section)
            if(bondsflag==1&&line[0]!='['&&line[0]!='#') {
                if(check_comment(line)==0) {
                    if(current==-1) {
                        cerr<<"ERROR: ENCOUNTERED [ bonds ] SECTION BEFORE [ moleculetype ] SECTION!\n";
                        exit(EXIT_FAILURE);
                    }
                    sscanf (line.c_str(),"%d %d",&n1,&n2);
                    add_bond(n1,n2);
                }
            }



            //if the new line is in an [ atoms ] section add the atom to the current moleculetype (also check for the existence of a preceding [ moleculetype ] section)
            if(atomsflag==1&&line[0]!='['&&line[0]!='#') {
                if(check_comment(line)==0) {
                    //check for the existence of a preceding [ moleculetype ] section
                    if(current==-1) {
                        cerr<<"ERROR: ENCOUNTERED [ atoms ] SECTION BEFORE [ moleculetype ] SECTION!\n";
                        exit(EXIT_FAILURE);
                    }
                    sscanf (line.c_str(),"%d %s %d %s %s %d %f %f",&n1,c1,&n2,c2,c3,&n3,&charge,&mass);
                    //add the mass of the atom to the current "moleculetypemasses" vector
                    add_mass(mass);
                    //add infromation about the atom (residuename, residuenumber, atomname) to the respective vectors
                    add_atom(c2,n2,c3);
                    //if the atom is a backbone atom, add it to the current "moleculetypebackbone" vector
                    isBackbone=false;
                    for(unsigned int i=0; i<backboneAtomNames.size(); i++) {
                        isBackbone=isBackbone||(!(strcmp(c3,backboneAtomNames[i].c_str())));
                    }
                    if(isBackbone) {
                        add_backbone(n1);
                    }
                }
            }



            if(moleculesflag==1&&line[0]!='['&&line[0]!='#') {

                if(check_comment(line)==0) { //add every moleculename in the [ molecules ] section to "moleculenames", taking care of how often the molecule is present
                    sscanf (line.c_str(),"%s %d",c1,&n2);
                    for(int i=0; i<n2; i++) {
                        moleculenames.push_back(string(c1));
                    }

                    //and check if the previous [ moleculetype ] section had exactly one [ atoms ] and [ bonds ] section
                    if(hasAtoms==0) {
                        cerr<<"ERROR: NO \"[ atoms ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                        exit(EXIT_FAILURE);
                    }
                    if((hasBonds==0)&&(moleculetypemasses[current].size()>1)) {
                        cerr<<"WARNING: NO \"[ bonds ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                    }
                    if(hasAtoms>1) {
                        cerr<<"ERROR: SEVERAL ("<<hasAtoms<<") \"[ atoms ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                        exit(EXIT_FAILURE);
                    }
                    if(hasBonds>1) {
                        cerr<<"ERROR: SEVERAL ("<<hasBonds<<") \"[ bonds ]\" SECTION FOUND FOR MOLECULETYPE "<<moleculetypenames[current]<<" IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
                        exit(EXIT_FAILURE);
                    }
                }
            }

            //if a new [ section ] is encountered, set the neccessary flags (and counters) according to the section type
            if(line[0]=='[') {
                atomsflag=0;
                bondsflag=0;
                moleculetypeflag=0;
                moleculesflag=0;
            }
            if((line.find("atoms")!=std::string::npos)&&(line[0]=='[')&&(line[line.length()-1]==']')&&(strip_blanks(line)=="[atoms]")) {
                atomsflag=1;
                hasAtoms+=1;
            }
            if((line.find("bonds")!=std::string::npos)&&(line[0]=='[')&&(line[line.length()-1]==']')&&(strip_blanks(line)=="[bonds]")) {
                bondsflag=1;
                hasBonds+=1;
            }
            if((line.find("moleculetype")!=std::string::npos)&&(line[0]=='[')&&(line[line.length()-1]==']')&&(strip_blanks(line)=="[moleculetype]")) {
                moleculetypeflag=1;
            }
            if((line.find("molecules")!=std::string::npos)&&(line[0]=='[')&&(line[line.length()-1]==']')&&(strip_blanks(line)=="[molecules]")) {
                moleculesflag=1;
            }
        }

        //check for errors
        if(moleculetypenames.size()==0) {
            cerr<<"ERROR: NO MOLECULETYPES FOUND IN (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
            exit(EXIT_FAILURE);
        }
        if(moleculenames.size()==0) {
            cerr<<"ERROR: NO MOLECULES FOUND IN THE \"[ molecules ]\" SECTION OF THE (PREPROCESSED) TOPOLOGY FILE \""<<top_file<<"\"!\n";
            exit(EXIT_FAILURE);
        }



        //according to the [ molecules ] section stored in moleculenames, add up the according moleculetypevectors to form the whole system (complex)
        if(sort()==1) {
            exit(EXIT_FAILURE);
        }

        //create a table: for every atom in the system, list all the bonded atoms
        make_bonds_table();


        return 0;
    }



    void make_bonds_table() {
        vector< vector <int> > dummyvec;
        vector <int> dummyvec1;
        dummyvec1.push_back(0);
        dummyvec1.push_back(0);
        for (int k=1; k<=total_atoms; k++) { //for all atoms k
            dummyvec.clear();
            for (unsigned int i=0; i<bonds.size(); i++) {
                for (unsigned int j=0; j<bonds[i].size(); j++) {//and all bonds in all molecules
                    if(k==bonds[i][j][0]) { //if atom k is the first atom of the according bond
                        dummyvec1[0]=bonds[i][j][1]; // add the second atom to atom k's bonds list
                        dummyvec1[1]=0; //and indicate that it's not a pseudobond
                        dummyvec.push_back(dummyvec1);
                    }
                    if(k==bonds[i][j][1]) { //if atom k is the second atom of the according bond
                        dummyvec1[0]=bonds[i][j][0]; // add the first atom to atom k's bonds list
                        dummyvec1[1]=0; //and indicate that it's not a pseudobond
                        dummyvec.push_back(dummyvec1);
                    }
                }
            }
            bonds_table.push_back(dummyvec); //add atom k's bonds list to the full bonds table
        }
        get_roots(); //find possible root atoms
    }



    void get_roots() {
        int counter;
        int sav;
        int improperscounter;
        int lowestimpropers;
        vector <int> tmproots;

        tmproots.push_back(0);
        tmproots.push_back(0);
        tmproots.push_back(0);
        tmproots.push_back(0);

        for (unsigned int i=0; i<molecules.size(); i++) { //find root atoms for all molecules
            lowestimpropers=molecules[i].size()+1; //find the roots which have the lowest number of (terminal) atoms connected to the middle root atom. (These account for the improper dihedrals)
            tmproots[0]=0;
            tmproots[1]=0;
            tmproots[2]=0;
            tmproots[3]=molecules[i].size();//this will carry over the information about moleculesizes to BAT_topology.cpp
            for (unsigned int j=0; j<molecules[i].size(); j++) { //for all atoms in the molecule
                if(bonds_table[molecules[i][j]-1].size()==1) { //if the atom is terminal
                    counter=0;//set up a counter for the non-terminal atoms connected to roots[i][1]
                    improperscounter=-1; //to count the terminal atoms connected to roots[i][1] (roots[i][0] is not considered an improper: 0->-1)
                    sav=-1; //to remember roots[i][2]
                    for(unsigned int k=0; k<bonds_table[bonds_table[molecules[i][j]-1][0][0]-1].size(); k++) {
                        if(bonds_table[bonds_table[bonds_table[molecules[i][j]-1][0][0]-1][k][0]-1].size()>1) {
                            counter++;    //count the non-terminal atoms connected to roots[i][1] and remember one of them
                            sav=bonds_table[bonds_table[molecules[i][j]-1][0][0]-1][k][0];
                        }
                        if(bonds_table[bonds_table[bonds_table[molecules[i][j]-1][0][0]-1][k][0]-1].size()==1) {
                            improperscounter++;   //count the terminal atoms connected to roots[i][1]
                        }
                    }
                    if((counter==1)&&(improperscounter<lowestimpropers)) { // if there is exactly one non-terminal atom connected to roots[i][1] and the amount of terminal atoms connected to roots[i][1] is smaller than for previous found root sets
                        lowestimpropers=improperscounter; //set the new goal for undercutting in the amount of terminal atoms connected to roots[i][1]
                        tmproots[0]=molecules[i][j]; //temporarily use the terminal atom t roots[i][0]
                        tmproots[1]=bonds_table[molecules[i][j]-1][0][0]; //temporarily use the only atom connected to the terminal atom as roots[i][1]
                        tmproots[2]=sav; //since counter is 1, only 1 atom connected to tmproots[1] is non-terminal. Temporarilly use this atom as roots[i][2]
                    }
                }
            }
            if((tmproots[0]==0)||(tmproots[1]==0)||(tmproots[2]==0)) { //if one of the temporary root atoms is still 0 quit
                cerr<<"ERROR: COULD NOT FIND ROOT ATOMS FOR MOLECULE "<<i+1<<".\n";
                exit(EXIT_FAILURE);
            }
            else {
                roots.push_back(tmproots); //else add these temporary root atoms as root atoms for molecule i
            }
        }
    }



};
