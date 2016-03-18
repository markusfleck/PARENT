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







using namespace std;

class BAT_Topology {

public:
    BAT_Topology(vector< vector< vector <int> > > *bonds_table_in, vector <int> *backboneIN, vector< vector <int> > *roots_in) {

        bonds_table=*bonds_table_in;
        backbone=*backboneIN;
        roots=*roots_in;

        current_BAT=1;//counter for the next BAT number to be assigned

        //allocate two arrays for the mapping of the BAT number to the atom number
        map_BAT_to_real=(int*)malloc(bonds_table.size()*sizeof(int));
        map_real_to_BAT=(int*)malloc(bonds_table.size()*sizeof(int));
        if((map_BAT_to_real==NULL)||(map_real_to_BAT==NULL)) {
            cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
            exit(EXIT_FAILURE);
        }
        //initialize the arrays
        for(unsigned int i=0; i<bonds_table.size(); i++) {
            map_BAT_to_real[i]=-1;
            map_real_to_BAT[i]=-1;
        }

        //~ do_mapping_depth_first();  //number the atoms beginning from the root atoms in depth-first style
        do_mapping_breadth_first();  //number the atoms beginning from the root atoms in breadth-first style. Mutually exclusive with do_mapping_depth_first();

        create_BAT();  //from the chosen numbering, create the BAT coordinates topology

        isBackbone=(bool*)calloc(bonds_table.size(),sizeof(bool));  //create an array for the assignment of phase angles
        if(isBackbone==NULL) {
            cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
            exit(EXIT_FAILURE);
        }
        add_phases();  //assign the phase angles

    };

    vector< vector< vector <int> > > bonds_table; //contains a list of all bonds in the system with proper offset, with bonds_table[#atom-1][#bond-1][0] the bonded atom and bonds_table[atom-1][#bond-1][1] the type of the bond (0==physical, 1==pseudo)
    vector< vector <int> > roots; //contains the root atoms for every molecule in the system
    vector< vector <int> > dihedrals; //contains a list of all dihedrals found
    vector <int> backbone; //contains a list of the backbone atoms in the system
    bool *isBackbone;
    int *map_BAT_to_real;
    int *map_real_to_BAT;

    int current_BAT;//counter for the next BAT number to be assigned


    void do_mapping_depth_first() {
        //for all molecules
        for(unsigned int i=0; i<roots.size(); i++) {
            //the atom numbers of the root atoms correspond to the first three BAT numbers in the current molecule
            map_BAT_to_real[current_BAT-1]=roots[i][0];
            map_BAT_to_real[current_BAT]=roots[i][1];
            map_BAT_to_real[current_BAT+1]=roots[i][2];
            map_real_to_BAT[roots[i][0]-1]=current_BAT;
            map_real_to_BAT[roots[i][1]-1]=current_BAT+1;
            map_real_to_BAT[roots[i][2]-1]=current_BAT+2;
            current_BAT+=3;


            add_node_recursive(roots[i][1]); //add the impropers for molecule i
            add_node_recursive(roots[i][2]); //branch from roots[i][2] for molecule i
        }
    }

    void add_node_recursive(int real_node) {
        //for all atoms connected to real_node
        for(unsigned int i=0; i<bonds_table[real_node-1].size(); i++) {
            //if the connected atom was not yet assigned a BAT number (and the bond is physical)
            if((map_real_to_BAT[bonds_table[real_node-1][i][0]-1]==-1)&&(bonds_table[real_node-1][i][1]==0)) {
                //assign it the current_BAT number
                map_real_to_BAT[bonds_table[real_node-1][i][0]-1]=current_BAT;
                map_BAT_to_real[current_BAT-1]=bonds_table[real_node-1][i][0];
                //and increase the current_BAT counter
                current_BAT++;
                //then recursively branch from this connected atom
                add_node_recursive(bonds_table[real_node-1][i][0]);
            }
        }
    }

    void do_mapping_breadth_first() {
        vector <int> to_be_processed;
        vector <int> to_be_processed_tmp;

        for(unsigned int j=0; j<roots.size(); j++) {
            //the atom numbers of the root atoms correspond to the first three BAT numbers in the current molecule
            map_BAT_to_real[current_BAT-1]=roots[j][0];
            map_BAT_to_real[current_BAT]=roots[j][1];
            map_BAT_to_real[current_BAT+1]=roots[j][2];
            map_real_to_BAT[roots[j][0]-1]=current_BAT;
            map_real_to_BAT[roots[j][1]-1]=current_BAT+1;
            map_real_to_BAT[roots[j][2]-1]=current_BAT+2;
            current_BAT+=3;


            to_be_processed.push_back(roots[j][1]); //add the impropers for molecule i
            to_be_processed.push_back(roots[j][2]); //branch from roots[i][2] for molecule i
            add_node_breadth(roots[j][1],&to_be_processed); //add the impropers for molecule i
            add_node_breadth(roots[j][2],&to_be_processed); //branch from roots[i][2] for molecule i

            //while there is still branching left to be done
            while(to_be_processed.size()!=0) {
                //copy to_be_processed to a temporary to work on ( deletion going on in the add_node_breadth function)
                to_be_processed_tmp=to_be_processed;
                //add all nodes in the to_be_processed list
                for(unsigned int i=0; i<to_be_processed.size(); i++) {
                    add_node_breadth(to_be_processed[i],&to_be_processed_tmp);
                }
                //update the to_be_processed list
                to_be_processed=to_be_processed_tmp;
            }
        }
    }

    void add_node_breadth(int real_node, vector <int> *to_be_processed_in) {
        //delete real node from the to_be_processed_in list
        for(unsigned int i=0; i<(*to_be_processed_in).size(); i++) {
            if(real_node== (*to_be_processed_in)[i]) {
                (*to_be_processed_in).erase((*to_be_processed_in).begin()+i);
            }
        }
        //for all atoms connected to real_node
        for(unsigned int i=0; i<bonds_table[real_node-1].size(); i++) {
            //if the connected atom was not yet assigned a BAT number
            if(map_real_to_BAT[bonds_table[real_node-1][i][0]-1]==-1) {
                //assign it the current_BAT number
                map_real_to_BAT[bonds_table[real_node-1][i][0]-1]=current_BAT;
                map_BAT_to_real[current_BAT-1]=bonds_table[real_node-1][i][0];
                //and increase the current_BAT counter
                current_BAT++;
                //and add the connected atom to the to_be_processed_in list for further branching
                (*to_be_processed_in).push_back(bonds_table[real_node-1][i][0]);
            }
        }
    }


    int find_closest_bond(int real_node, vector<int> exclude_list) {
        int closest=int(bonds_table.size()+1);
        int flag;
        //for all atoms bonded to real_node
        for(unsigned int i=0; i<bonds_table[real_node-1].size(); i++) {
            //if the BAT number of this bonded atom is smaller than the current smallest ("closest")
            if(map_real_to_BAT[bonds_table[real_node-1][i][0]-1]<closest) {
                //check if this bonded atom is in the "exclude_list"
                flag=1;
                for(unsigned int j=0; j<exclude_list.size(); j++) {
                    if(exclude_list[j]==bonds_table[real_node-1][i][0]) {
                        flag=0;
                    }
                }
                //if it is not in the "exclude_list" update the "closest" with the BAT number of the found bonded atom
                if (flag==1) {
                    closest=map_real_to_BAT[bonds_table[real_node-1][i][0]-1];
                }
            }
        }
        if(closest==int(bonds_table.size()+1)) {
            cerr<<"ERROR BUILDING BAT TOPOLOGY. QUITTING."<<endl;
            exit(EXIT_FAILURE);
        }
        //return the atom number of the bonded atom closest to the root atoms (smallest BAT number)
        return map_BAT_to_real[closest-1];
    }


    void create_BAT() {
        //set up a temporary array
        vector<int> dummydihedral;
        dummydihedral.push_back(0);
        dummydihedral.push_back(0);
        dummydihedral.push_back(0);
        dummydihedral.push_back(0);
        dummydihedral.push_back(0);

        vector<int> excludes;
        int offset=0;
        //for all molecules
        for(unsigned int j=0; j<roots.size(); j++) {
            if(j>0) {
                //for a second or further molecule, adjust the offset adding the number of atoms of the previous molecule
                offset+=roots[j-1][3];
                //and connect both sets of root atoms per pseudobonds
                dummydihedral[3]=roots[j][0];
                dummydihedral[2]=roots[j-1][2];
                dummydihedral[1]=roots[j-1][1];
                dummydihedral[0]=roots[j-1][0];
                dummydihedral[4]=1;//indicating pseudobond
                dihedrals.push_back(dummydihedral);

                dummydihedral[3]=roots[j][1];
                dummydihedral[2]=roots[j][0];
                dummydihedral[1]=roots[j-1][2];
                dummydihedral[0]=roots[j-1][1];
                dummydihedral[4]=1;//indicating pseudobond
                dihedrals.push_back(dummydihedral);

                dummydihedral[3]=roots[j][2];
                dummydihedral[2]=roots[j][1];
                dummydihedral[1]=roots[j][0];
                dummydihedral[0]=roots[j-1][2];
                dummydihedral[4]=1;//indicating pseudobond
                dihedrals.push_back(dummydihedral);
            }
            for(unsigned int i=0; i<bonds_table[roots[j][1]-1].size(); i++) {
                if((bonds_table[roots[j][1]-1][i][0]!=roots[j][0])&&(bonds_table[roots[j][1]-1][i][0]!=roots[j][2])) {
                    dummydihedral[3]=bonds_table[roots[j][1]-1][i][0];
                    dummydihedral[2]=roots[j][1];
                    dummydihedral[1]=roots[j][0];
                    dummydihedral[0]=roots[j][2];
                    dummydihedral[4]=-1;//indicating improper dihedral
                    dihedrals.push_back(dummydihedral);
                }
            }

            //for all BAT numbers  except the root atoms and the impropers for the molecule
            for(int i=offset+3+bonds_table[roots[j][1]-1].size()-2; i<roots[j][3]+offset; i++) {
                excludes.clear();
                //this BAT number is the first constituent of the new bond
                dummydihedral[3]=map_BAT_to_real[i];
                //the second constituent is the connected atom with the lowest BAT number (closest to the root atoms concerning BAT sequence )
                dummydihedral[2]=find_closest_bond(dummydihedral[3],excludes);
                excludes.push_back(dummydihedral[3]); //prohibit a second inclusion of dummydihedral[3] in dummydihedral
                //the third constituent is the atom connected to the second constituent with the lowest BAT number (closest to the root atoms concerning BAT sequence )
                dummydihedral[1]=find_closest_bond(dummydihedral[2],excludes);
                excludes.push_back(dummydihedral[2]); //prohibit a second inclusion of dummydihedral[2] in dummydihedral
                //the fourth constituent is the atom connected to the third constituent with the lowest BAT number (closest to the root atoms concerning BAT sequence )
                dummydihedral[0]=find_closest_bond(dummydihedral[1],excludes);
                dummydihedral[4]=0; //indicating a non-improper physical bond
                dihedrals.push_back(dummydihedral);
            }
        }
        cout<<"Found "<<dihedrals.size()<<" dihedrals/impropers."<<endl;
        //check if the proper amount of dihedrals was found
        offset+=roots[roots.size()-1][3];
        if(int(dihedrals.size())!=offset-3) {
            cerr<<"ERROR: "<<offset-3<<" DIHEDRALS/IMPROPERS SHOULD HAVE BEEN FOUND, BUT ONLY "<<dihedrals.size()<<" WERE ACTUALLY FOUND."<<endl;
            exit(EXIT_FAILURE);
        }



    }

    void add_phases()
    {
        vector <int> excludes;
        vector <int> dummyvec;
        dummyvec.push_back(0);
        dummyvec.push_back(0);
        vector <vector <int> > common;
        bool excluded,firstflag,isBackboneDihedral;
        int fullbackbone;

        //for every "starting" dihedral
        for(unsigned int i=0; i<dihedrals.size(); i++) {
            excluded=false;
            firstflag=true;
            //check if the "starting" dihedral is in the excludes list
            for(unsigned int k=0; k<excludes.size(); k++) {
                excluded=excluded||(int(i)==excludes[k]);
            }
            //if the "starting" dihedral is not in the excludes list
            if(!excluded) {
                //for any "later" dihedral
                for(unsigned int j=i+1; j<dihedrals.size(); j++) {
                    //check if the "later" dihedral is in the excludes list
                    excluded=false;
                    for(unsigned int k=0; k<excludes.size(); k++) {
                        excluded=excluded||(int(j)==excludes[k]);
                    }
                    //if the "later" dihedral is not in the excludes list
                    if(!excluded) {
                        //if the "starting" dihedral shares the first 3 atoms with the "later" dihedral
                        if((dihedrals[i][0]==dihedrals[j][0])&&(dihedrals[i][1]==dihedrals[j][1])&&(dihedrals[i][2]==dihedrals[j][2])) {
                            //and the "later" dihedral is the first one sharing three atoms with the "starting" dihedral
                            if(firstflag) {
                                dummyvec[0]=i;
                                dummyvec[1]=j;
                                common.push_back(dummyvec); //open up a new common group containing the "starting" and the "latter" dihedral
                                firstflag=false;
                                excludes.push_back(i); //and put both dihedrals to the excludes list (otherwise the common groups would share dihedrals)
                                excludes.push_back(j);
                            }
                            else {// if already a matching "later" dihedral has been found before
                                common[common.size()-1].push_back(j); //add the new "later" dihedral to the current common group
                                excludes.push_back(j); //and exclude the new "later" dihedral
                            }
                        }
                    }
                }
            }
        }

        //for every atom "i" in the system set isBackbone[i] to 1 if it is a backbone atom and to 0 if not
        for(unsigned int i=0; i<bonds_table.size(); i++) {
            for(unsigned int j=0; j<backbone.size(); j++) {
                isBackbone[i]=isBackbone[i]||(int(i)==(backbone[j]-1));
            }
        }


        //for ever dihedral in the system add a 6th component describing which dihedral is the "parent" dihedral of this dihedral/phaseangle (-1 means is a "parent" dihedral itself)
        for(unsigned int i=0; i<dihedrals.size(); i++) {
            dihedrals[i].push_back(-1);
        }

        //for all common groups
        for(unsigned int i=0; i<common.size(); i++) {
            fullbackbone=-2;
            //and all dihedrals in the current common group
            for(unsigned int j=0; j<common[i].size(); j++) {
                isBackboneDihedral=isBackbone[dihedrals[common[i][j]][0]-1]&&isBackbone[dihedrals[common[i][j]][1]-1]&&isBackbone[dihedrals[common[i][j]][2]-1]&&isBackbone[dihedrals[common[i][j]][3]-1];
                if(isBackboneDihedral) { //if the investigated dihedral consists only of backbone atoms
                    fullbackbone=j; //save the the number of this dihedral in "fullbackbone" and exit the loop
                    break;
                }
            }
            if(fullbackbone==-2) { //if no dihedral in the current common group consisted only of backbone atoms
                dihedrals[common[i][0]][5]=-1; //take the first dihedral in the current common group as "parent"
                for(unsigned int j=1; j<common[i].size(); j++) {
                    dihedrals[common[i][j]][5]=common[i][0]+1; //and declare all other dihedrals as phaseangles to this first dihedral
                }
            }
            else {//if a dihedral in the current common group indeed consisted only of backbone atoms
                dihedrals[common[i][fullbackbone]][5]=-1; //assign it to be "parent"
                for(unsigned int j=0; j<common[i].size(); j++) {
                    if(int(j)!=fullbackbone) {
                        dihedrals[common[i][j]][5]=common[i][fullbackbone]+1; //and declare all other dihedrals as phaseangles to this backbone dihedral
                    }
                }
            }
        }

    }



};