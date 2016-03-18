//    A program to calculate the MIST approximation of the configurationalational entropy from the output of the PARENT program
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







#define MUTUAL_ARGS nDihedrals,bondsEntropy1D,anglesEntropy1D,dihedralsEntropy1D,bbEntropy,baEntropy,bdEntropy,aaEntropy,adEntropy,ddEntropy


#include "mpi.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <sys/time.h>


#define MASTER 0
#define DISTRIBUTE_MUTUAL_TAG 0
#define DISTRIBUTE_HIGHESTTYPE1_TAG 1
#define DISTRIBUTE_HIGHESTTYPE2_TAG 2
#define DISTRIBUTE_HIGHESTINDEX1_TAG 3
#define DISTRIBUTE_HIGHESTINDEX2_TAG 4
#define DISTRIBUTE_MUTUAL_READIN_TAG 5

using namespace std;



#include "../util/io/io.h"
#include "../util/util.h"




int main(int argc, char* argv[])
{

    //start the stopwatch for the execution time
    timeval tv1;
    gettimeofday (&tv1, NULL);
    timeval tv2,tv3;



    if(argc!=3) {
        cerr<<"USAGE:\n"<<argv[0]<<" input.par output.txt"<<endl;
        return 1;
    }


    int numProcesses, rank, rc,len,threadLevelProvided, threadLevelClaimed;
    char hostName[MPI_MAX_PROCESSOR_NAME];





    //Initialize MPI
    rc=0;
    MPI_Status Stat;
    MPI_Init_thread(&argc,&argv, MPI_THREAD_SERIALIZED, &threadLevelProvided );
    MPI_Query_thread( &threadLevelClaimed );
    if(threadLevelClaimed !=threadLevelProvided)
    {
        printf( "Claimed thread level=%d, but got thread level=%d, aborting\n", threadLevelClaimed , threadLevelProvided );
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Get_processor_name(hostName, &len);

    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf("Initiated process = %i of NCPU=%i processes on machine=%s, thread number %i\n", rank, numProcesses, hostName ,ID);
    }

    //declare and initialize further variables
    int double_prec,numFrames;
    int version,bDens,aDens,dDens,bDens1D,aDens1D,dDens1D;
    vector< vector <int> > dihedrals_top;
    vector <float> masses;
    vector <string> residues;
    vector <int> residueNumbers;
    vector <string> atomNames;
    vector <string> belongsToMolecule;

    double *bondsEntropy1D,*anglesEntropy1D,*dihedralsEntropy1D;
    double *bbEntropy,*aaEntropy,*ddEntropy,*baEntropy,*bdEntropy,*adEntropy;

    int nBonds,nAngles,nDihedrals;

    double totalBondsEntropy=0;
    double totalAnglesEntropy=0;
    double totalDihedralsEntropy=0;

    double totalBbMutual=0;
    double totalBaMutual=0;
    double totalBdMutual=0;
    double totalAaMutual=0;
    double totalAdMutual=0;
    double totalDdMutual=0;

    double totalMutual=0;


    int b=0;
    int a=1;
    int d=2;

    char myChar[3];
    myChar[TYPE_B]='b';
    myChar[TYPE_A]='a';
    myChar[TYPE_D]='d';


    double*** mutual;
    double highestMutual;
    int highestType1,highestType2,highestIndex1,highestIndex2;

    vector <int> In[3];
    vector <int> Out[3];

    int myInBegin2D[3][3],myInEnd2D[3][3],myOutBegin2D[3][3],myOutEnd2D[3][3];

    double doubledummy;
    int intdummy;
    bool allocation_worked;

    ifstream infile;
    ofstream outfile;
    cout.precision(12);
    outfile.precision(12);


    bool pb=1; //set true for bonds calculation
    bool pa=1; //set true for angles calculation
    bool pd=1; //set true for dihedrals calculation
    if(rank==MASTER) {
        cout<<endl;
        cout<<endl;
        cout<<endl;
        cout<<endl;
        cout<<endl;
        cout<<"    A program to calculate the MIST approximation of the configurational entropy from the output of the PARENT program"<<endl;
        cout<<"    Copyright (C) 2015  Markus Fleck, as a representative of the working group of Bojan Zagrovic, University of Vienna"<<endl;
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
        if((!pb)&&(!pa)&&!(pd)) {
            cerr<<"ERROR: YOU CANNOT IGNORE ALL DEGREES OF FREEDOM !"<<endl;
            return 1;
        }
        outfile.open(argv[2], ios::out);

        infile.open(argv[1], ios::binary | ios::in); //MASTER process opens the .par file and reads the header
        if(infile.is_open()) {
            if(outfile.is_open()) {
                if(read_PAR_header(&infile,&nDihedrals,&double_prec,&numFrames,&dihedrals_top, &masses, &version, &bDens, &aDens, &dDens, &bDens1D, &aDens1D, &dDens1D,&residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) {
                    cerr<<"An ERROR has occurred while reading the file " <<argv[1]<<" . Quitting Program.\n";
                    MPI_Abort(MPI_COMM_WORLD, rc);
                }
                nDihedrals=dihedrals_top.size();
            }
            else {
                cerr<<"ERROR: COULD NOT OPEN FILE "<<argv[2]<<" !"<<endl;
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
        }
        else {
            cerr<<"ERROR: COULD NOT OPEN FILE "<<argv[1]<<" !"<<endl;
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    MPI_Bcast( &nDihedrals, 1, MPI_INT, MASTER, MPI_COMM_WORLD); //the number of dihedrals is broadcasted to all MPI processes

    nAngles=nDihedrals+1;//and the number of bonds and angles is set
    nBonds=nDihedrals+2;

    mutual=(double***)malloc(6*sizeof(double**)); //every MPI process allocates a mutual information array using 3 indices
    allocation_worked=(mutual!=NULL);

    mutual[TYPE_BB]=(double**)malloc(nBonds*sizeof(double*));
    mutual[TYPE_BA]=(double**)malloc(nBonds*sizeof(double*));
    mutual[TYPE_BD]=(double**)malloc(nBonds*sizeof(double*));
    mutual[TYPE_AA]=(double**)malloc(nAngles*sizeof(double*));
    mutual[TYPE_AD]=(double**)malloc(nAngles*sizeof(double*));
    mutual[TYPE_DD]=(double**)malloc(nDihedrals*sizeof(double*));
    allocation_worked=(allocation_worked&&(mutual[TYPE_BB]!=NULL)&&(mutual[TYPE_BA]!=NULL)&&(mutual[TYPE_BD]!=NULL)&&(mutual[TYPE_AA]!=NULL)&&(mutual[TYPE_AD]!=NULL)&&(mutual[TYPE_DD]!=NULL));

    for(int i=0; i<nBonds; i++) {
        mutual[TYPE_BB][i]=(double*)malloc(nBonds*sizeof(double));
        allocation_worked=(allocation_worked&&(mutual[TYPE_BB][i]!=NULL));
    }
    for(int i=0; i<nBonds; i++) {
        mutual[TYPE_BA][i]=(double*)malloc(nAngles*sizeof(double));
        allocation_worked=(allocation_worked&&(mutual[TYPE_BA][i]!=NULL));
    }
    for(int i=0; i<nBonds; i++) {
        mutual[TYPE_BD][i]=(double*)malloc(nDihedrals*sizeof(double));
        allocation_worked=(allocation_worked&&(mutual[TYPE_BD][i]!=NULL));
    }
    for(int i=0; i<nAngles; i++) {
        mutual[TYPE_AA][i]=(double*)malloc(nAngles*sizeof(double));
        allocation_worked=(allocation_worked&&(mutual[TYPE_AA][i]!=NULL));
    }
    for(int i=0; i<nAngles; i++) {
        mutual[TYPE_AD][i]=(double*)malloc(nDihedrals*sizeof(double));
        allocation_worked=(allocation_worked&&(mutual[TYPE_AD][i]!=NULL));
    }
    for(int i=0; i<nDihedrals; i++) {
        mutual[TYPE_DD][i]=(double*)malloc(nDihedrals*sizeof(double));
        allocation_worked=(allocation_worked&&(mutual[TYPE_DD][i]!=NULL));
    }

    if(!allocation_worked) {
        cerr<<"ERROR: ALLOCATION FAILED ON NODE "<<hostName<<". MORE MEMORY NEEDED?"<<endl;
    }

    if(rank==MASTER) {

        if(read_PAR_body(&infile,nDihedrals,&bondsEntropy1D, &anglesEntropy1D, &dihedralsEntropy1D, &bbEntropy, &baEntropy, &bdEntropy, &aaEntropy, &adEntropy, &ddEntropy)!=0) {
            cerr<<"An ERROR has occurred while reading the file " <<argv[1]<<" . Quitting Program.\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }

        for(int k=0; k<nBonds; k++) {
            totalBondsEntropy+=bondsEntropy1D[k];   //calculates the total 1D entropies
        }
        for(int k=0; k<nAngles; k++) {
            totalAnglesEntropy+=anglesEntropy1D[k];
        }
        for(int k=0; k<nDihedrals; k++) {
            totalDihedralsEntropy+=dihedralsEntropy1D[k];
        }

        //then for every bond-bond pair the MASTER node gets the mutual information and stores it
        for(int i=0; i<nBonds-1; i++) {
            for(int j=i+1; j<nBonds; j++) {
                mutual[TYPE_BB][i][j]=get_mutual(TYPE_BB,i,j,MUTUAL_ARGS);
                mutual[TYPE_BB][j][i]=mutual[TYPE_BB][i][j]; //also symmetrize the bonds-bonds matrix
            }
        }

        for(int i=0; i<nBonds; i++) { //same for all bond-angle pairs
            for(int j=0; j<nAngles; j++) {
                mutual[TYPE_BA][i][j]=get_mutual(TYPE_BA,i,j,MUTUAL_ARGS);
            }
        }

        for(int i=0; i<nBonds; i++) { //bond-dihedral pairs
            for(int j=0; j<nDihedrals; j++) {
                mutual[TYPE_BD][i][j]=get_mutual(TYPE_BD,i,j,MUTUAL_ARGS);
            }
        }

        for(int i=0; i<nAngles-1; i++) {//angle-angle pairs
            for(int j=i+1; j<nAngles; j++) {
                mutual[TYPE_AA][i][j]=get_mutual(TYPE_AA,i,j,MUTUAL_ARGS);
                mutual[TYPE_AA][j][i]=mutual[TYPE_AA][i][j];
            }
        }

        for(int i=0; i<nAngles; i++) { //angle-dihedral pairs
            for(int j=0; j<nDihedrals; j++) {
                mutual[TYPE_AD][i][j]=get_mutual(TYPE_AD,i,j,MUTUAL_ARGS);
            }
        }

        for(int i=0; i<nDihedrals-1; i++) {//dihedral-dihedral pairs
            for(int j=i+1; j<nDihedrals; j++) {
                mutual[TYPE_DD][i][j]=get_mutual(TYPE_DD,i,j,MUTUAL_ARGS);
                mutual[TYPE_DD][j][i]=mutual[TYPE_DD][i][j];
            }
        }
    }

    //then the MASTER process sends the mutual information array to all other processes
    for(int i=0; i<nBonds; i++) { //the bonds-bonds array
        for(int j=1; j<numProcesses; j++) {
            if(rank==MASTER) {
                MPI_Send(&(mutual[TYPE_BB][i][0]), nBonds, MPI_DOUBLE, j, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD);
            }
            else if(rank==j) {
                MPI_Recv(&(mutual[TYPE_BB][i][0]), nBonds, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD, &Stat);
            }
        }
    }

    for(int i=0; i<nBonds; i++) { // the bonds-angles array
        for(int j=1; j<numProcesses; j++) {
            if(rank==MASTER) {
                MPI_Send(&(mutual[TYPE_BA][i][0]), nAngles, MPI_DOUBLE, j, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD);
            }
            else if(rank==j) {
                MPI_Recv(&(mutual[TYPE_BA][i][0]), nAngles, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD, &Stat);
            }
        }
    }

    for(int i=0; i<nBonds; i++) { //the bonds-dihedrals array
        for(int j=1; j<numProcesses; j++) {
            if(rank==MASTER) {
                MPI_Send(&(mutual[TYPE_BD][i][0]), nDihedrals, MPI_DOUBLE, j, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD);
            }
            else if(rank==j) {
                MPI_Recv(&(mutual[TYPE_BD][i][0]), nDihedrals, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD, &Stat);
            }
        }
    }

    for(int i=0; i<nAngles; i++) { // the angles-angles array
        for(int j=1; j<numProcesses; j++) {
            if(rank==MASTER) {
                MPI_Send(&(mutual[TYPE_AA][i][0]), nAngles, MPI_DOUBLE, j, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD);
            }
            else if(rank==j) {
                MPI_Recv(&(mutual[TYPE_AA][i][0]), nAngles, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD, &Stat);
            }
        }
    }

    for(int i=0; i<nAngles; i++) { //the angles-dihedrals array
        for(int j=1; j<numProcesses; j++) {
            if(rank==MASTER) {
                MPI_Send(&(mutual[TYPE_AD][i][0]), nDihedrals, MPI_DOUBLE, j, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD);
            }
            else if(rank==j) {
                MPI_Recv(&(mutual[TYPE_AD][i][0]), nDihedrals, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD, &Stat);
            }
        }
    }

    for(int i=0; i<nDihedrals; i++) { //the dihedrals-dihedrals array
        for(int j=1; j<numProcesses; j++) {
            if(rank==MASTER) {
                MPI_Send(&(mutual[TYPE_DD][i][0]), nDihedrals, MPI_DOUBLE, j, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD);
            }
            else if(rank==j) {
                MPI_Recv(&(mutual[TYPE_DD][i][0]), nDihedrals, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_READIN_TAG, MPI_COMM_WORLD, &Stat);
            }
        }
    }




    //all processes create "Out" vectors for the degrees of freedom to keep track of which of them have NOT already been processed
    if(pb) {
        for(int i=0; i<nBonds; i++) {
            Out[TYPE_B].push_back(i);
        }
    }
    if(pa) {
        for(int i=0; i<nAngles; i++) {
            Out[TYPE_A].push_back(i);
        }
    }
    if(pd) {
        for(int i=0; i<nDihedrals; i++) {
            Out[TYPE_D].push_back(i);
        }
    }


    // to start out, all processes put one degree of freedom into the "In" vector, which keeps track of the degrees of freedom which already have been processed
    if(pb) {
        In[TYPE_B].push_back(Out[TYPE_B][0]);
        Out[TYPE_B].erase(Out[TYPE_B].begin());
    }
    else {
        if(pa) {
            In[TYPE_A].push_back(Out[TYPE_A][0]);
            Out[TYPE_A].erase(Out[TYPE_A].begin());
        }
        else {
            In[TYPE_D].push_back(Out[TYPE_D][0]);
            Out[TYPE_D].erase(Out[TYPE_D].begin());
        }
    }



    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==MASTER) {
        gettimeofday (&tv2, NULL);   //time is measured
    }
    //then start the actual calculation using MPI as well as openMP parallelization
    #pragma omp parallel
    {
        while(Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()>0) {
            double highestMutualThread;
            int highestType1Thread,highestType2Thread,highestIndex1Thread,highestIndex2Thread;
            #pragma omp master
            {
                highestMutual=-1e200;//the global highest mutual information term per MPI process is set by the master thread to a very low value at the beginning of each cycle

                //every MPI process determines the parts of the vectors it should calculate. Depending on which one is bigger, either the OutVector or InVector is split.
                //If the size of the vector is not an integer multiple of the number of MPI processes, the additional elements are added to the highest MPI processes to the relieve MASTER process.
                //e. g.: 3 MPI Processes, 2 elements in Out[TYPE_B], 8 elements in In[TYPE_B]. Then Processes 1 treats elements 7 and 8 in In[TYPE_B] and all elements in Out[TYPE_B],
                //Processes 2 treats elements 4, 5 and 6 in In[TYPE_B] and all elements in Out[TYPE_B]
                //and Process 3 treats elements 1, 2 and 3 in In[TYPE_B] and all elements in Out[TYPE_B].
                //Also take care that not e. g. In[TYPE_A] and Out[TYPE_D] are split and used together in the coming loops (two indices at myInBegin2D, myInEnd2D  and so on)
                for(int i=0; i<3; i++) {
                    for(int j=0; j<3; j++) {
                        if(In[i].size()>Out[j].size()) {
                            myInBegin2D[i][j]=(numProcesses-1-rank)*(In[i].size()/numProcesses)+(((numProcesses-1-rank)<(int(In[i].size())%numProcesses))?(numProcesses-1-rank):(In[i].size()%numProcesses));
                            myInEnd2D[i][j]=((numProcesses-1-rank)+1)*(In[i].size()/numProcesses)+((((numProcesses-1-rank)+1)<(int(In[i].size())%numProcesses))?((numProcesses-1-rank)+1):(In[i].size()%numProcesses));
                            myOutBegin2D[i][j]=0;
                            myOutEnd2D[i][j]=Out[j].size();
                        }
                        else {
                            myInBegin2D[i][j]=0;
                            myInEnd2D[i][j]=In[i].size();
                            myOutBegin2D[i][j]=(numProcesses-1-rank)*(Out[j].size()/numProcesses)+(((numProcesses-1-rank)<(int(Out[j].size())%numProcesses))?(numProcesses-1-rank):(Out[j].size()%numProcesses));
                            myOutEnd2D[i][j]=((numProcesses-1-rank)+1)*(Out[j].size()/numProcesses)+((((numProcesses-1-rank)+1)<(int(Out[j].size())%numProcesses))?((numProcesses-1-rank)+1):(Out[j].size()%numProcesses));
                        }
                    }
                }





            }
            #pragma omp barrier

            highestMutualThread=-1e200; //then every thread sets its local highest mutual information term to a very low value (at the beginning of each cycle)
            if(pb) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_B][TYPE_B]; i<myInEnd2D[TYPE_B][TYPE_B]; i++) { //then for every bond in the "In" vector of the MPI process
                    for(int j=myOutBegin2D[TYPE_B][TYPE_B]; j<myOutEnd2D[TYPE_B][TYPE_B]; j++) { //and for every bond in the "Out" vector of the MPI process
                        if(mutual[TYPE_BB][In[TYPE_B][i]][Out[TYPE_B][j]]>highestMutualThread) { //if the mutual information term is higher than the previously highest
                            highestMutualThread=mutual[TYPE_BB][In[TYPE_B][i]][Out[TYPE_B][j]]; //store the value
                            highestType1Thread=b; //the type of the mutual informtion term (here In==b and Out==b)
                            highestType2Thread=b;
                            highestIndex1Thread=i; //and the indices
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pb&&pa) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_B][TYPE_A]; i<myInEnd2D[TYPE_B][TYPE_A]; i++) { //do the same for bonds In-angles Out
                    for(int j=myOutBegin2D[TYPE_B][TYPE_A]; j<myOutEnd2D[TYPE_B][TYPE_A]; j++) {
                        if(mutual[TYPE_BA][In[TYPE_B][i]][Out[TYPE_A][j]]>highestMutualThread) {
                            highestMutualThread=mutual[TYPE_BA][In[TYPE_B][i]][Out[TYPE_A][j]];
                            highestType1Thread=b;
                            highestType2Thread=a;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pb&&pd) { //same for bonds In-dihedrals Out
                #pragma omp for
                for(int i=myInBegin2D[TYPE_B][TYPE_D]; i<myInEnd2D[TYPE_B][TYPE_D]; i++) {
                    for(int j=myOutBegin2D[TYPE_B][TYPE_D]; j<myOutEnd2D[TYPE_B][TYPE_D]; j++) {
                        if(mutual[TYPE_BD][In[TYPE_B][i]][Out[TYPE_D][j]]>highestMutualThread) {
                            highestMutualThread=mutual[TYPE_BD][In[TYPE_B][i]][Out[TYPE_D][j]];
                            highestType1Thread=b;
                            highestType2Thread=d;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pa&&pb) { //same for angles In-bonds Out
                #pragma omp for
                for(int i=myInBegin2D[TYPE_A][TYPE_B]; i<myInEnd2D[TYPE_A][TYPE_B]; i++) {
                    for(int j=myOutBegin2D[TYPE_A][TYPE_B]; j<myOutEnd2D[TYPE_A][TYPE_B]; j++) {
                        if(mutual[TYPE_BA][Out[TYPE_B][j]][In[TYPE_A][i]]>highestMutualThread) { //notice that the Out[TYPE_B]-index comes before the In[TYPE_A]-index for mutual[TYPE_BA]
                            highestMutualThread=mutual[TYPE_BA][Out[TYPE_B][j]][In[TYPE_A][i]];
                            highestType1Thread=a;
                            highestType2Thread=b;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pa) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_A][TYPE_A]; i<myInEnd2D[TYPE_A][TYPE_A]; i++) { //same for anglesIn-anglesOut
                    for(int j=myOutBegin2D[TYPE_A][TYPE_A]; j<myOutEnd2D[TYPE_A][TYPE_A]; j++) {
                        if(mutual[TYPE_AA][In[TYPE_A][i]][Out[TYPE_A][j]]>highestMutualThread) {
                            highestMutualThread=mutual[TYPE_AA][In[TYPE_A][i]][Out[TYPE_A][j]];
                            highestType1Thread=a;
                            highestType2Thread=a;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pa&&pd) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_A][TYPE_D]; i<myInEnd2D[TYPE_A][TYPE_D]; i++) { //same for anglesIn-dihedralsOut
                    for(int j=myOutBegin2D[TYPE_A][TYPE_D]; j<myOutEnd2D[TYPE_A][TYPE_D]; j++) {
                        if(mutual[TYPE_AD][In[TYPE_A][i]][Out[TYPE_D][j]]>highestMutualThread) {
                            highestMutualThread=mutual[TYPE_AD][In[TYPE_A][i]][Out[TYPE_D][j]];
                            highestType1Thread=a;
                            highestType2Thread=d;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pd&&pb) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_D][TYPE_B]; i<myInEnd2D[TYPE_D][TYPE_B]; i++) { //same for dihedralsIn-bondsOut
                    for(int j=myOutBegin2D[TYPE_D][TYPE_B]; j<myOutEnd2D[TYPE_D][TYPE_B]; j++) {
                        if(mutual[TYPE_BD][Out[TYPE_B][j]][In[TYPE_D][i]]>highestMutualThread) { //notice that the Out[TYPE_B]-index comes before the In[TYPE_D]-index for mutual[TYPE_BD]
                            highestMutualThread=mutual[TYPE_BD][Out[TYPE_B][j]][In[TYPE_D][i]];
                            highestType1Thread=d;
                            highestType2Thread=b;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pd&&pa) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_D][TYPE_A]; i<myInEnd2D[TYPE_D][TYPE_A]; i++) { //same for dihedralsIn-anglesOut
                    for(int j=myOutBegin2D[TYPE_D][TYPE_A]; j<myOutEnd2D[TYPE_D][TYPE_A]; j++) {
                        if(mutual[TYPE_AD][Out[TYPE_A][j]][In[TYPE_D][i]]>highestMutualThread) { //notice that the Out[TYPE_A]-index comes before the In[TYPE_D]-index for mutual[TYPE_AD]
                            highestMutualThread=mutual[TYPE_AD][Out[TYPE_A][j]][In[TYPE_D][i]];
                            highestType1Thread=d;
                            highestType2Thread=a;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }
            if(pd) {
                #pragma omp for
                for(int i=myInBegin2D[TYPE_D][TYPE_D]; i<myInEnd2D[TYPE_D][TYPE_D]; i++) { //same for dihedralsIn-dihedralsOut
                    for(int j=myOutBegin2D[TYPE_D][TYPE_D]; j<myOutEnd2D[TYPE_D][TYPE_D]; j++) {
                        if(mutual[TYPE_DD][In[TYPE_D][i]][Out[TYPE_D][j]]>highestMutualThread) {
                            highestMutualThread=mutual[TYPE_DD][In[TYPE_D][i]][Out[TYPE_D][j]];
                            highestType1Thread=d;
                            highestType2Thread=d;
                            highestIndex1Thread=i;
                            highestIndex2Thread=j;
                        }
                    }
                }
            }

            #pragma omp critical
            {
                if(highestMutualThread>highestMutual) { //if the highest mutual information term of the current thread is higher than the global value of all threads in the MPI process
                    highestMutual=highestMutualThread; //store this value as the global one
                    highestType1=highestType1Thread; //also save the types
                    highestType2=highestType2Thread;
                    highestIndex1=highestIndex1Thread; //and indices
                    highestIndex2=highestIndex2Thread;
                }
            }
            #pragma omp barrier

            //the MASTER thread of the MASTER MPI process collects all highest mutual information values and stores the highest of them all
            #pragma omp master
            {
                if(rank==MASTER) {
                    for(int i=1; i<numProcesses; i++) {
                        MPI_Recv(&doubledummy, 1, MPI_DOUBLE, i, DISTRIBUTE_MUTUAL_TAG, MPI_COMM_WORLD, &Stat);
                        MPI_Recv(&intdummy, 1, MPI_INT, i, DISTRIBUTE_HIGHESTTYPE1_TAG, MPI_COMM_WORLD, &Stat);
                        if(doubledummy>highestMutual) {
                            highestType1=intdummy;
                        }
                        MPI_Recv(&intdummy, 1, MPI_INT, i, DISTRIBUTE_HIGHESTTYPE2_TAG, MPI_COMM_WORLD, &Stat);
                        if(doubledummy>highestMutual) {
                            highestType2=intdummy;
                        }
                        MPI_Recv(&intdummy, 1, MPI_INT, i, DISTRIBUTE_HIGHESTINDEX1_TAG, MPI_COMM_WORLD, &Stat);
                        if(doubledummy>highestMutual) {
                            highestIndex1=intdummy;
                        }
                        MPI_Recv(&intdummy, 1, MPI_INT, i, DISTRIBUTE_HIGHESTINDEX2_TAG, MPI_COMM_WORLD, &Stat);
                        if(doubledummy>highestMutual) {
                            highestIndex2=intdummy;
                            highestMutual=doubledummy;
                        }
                    }
                }
                else {
                    MPI_Send(&highestMutual, 1, MPI_DOUBLE, MASTER, DISTRIBUTE_MUTUAL_TAG, MPI_COMM_WORLD);
                    MPI_Send(&highestType1, 1, MPI_INT, MASTER, DISTRIBUTE_HIGHESTTYPE1_TAG, MPI_COMM_WORLD);
                    MPI_Send(&highestType2, 1, MPI_INT, MASTER, DISTRIBUTE_HIGHESTTYPE2_TAG, MPI_COMM_WORLD);
                    MPI_Send(&highestIndex1, 1, MPI_INT, MASTER, DISTRIBUTE_HIGHESTINDEX1_TAG, MPI_COMM_WORLD);
                    MPI_Send(&highestIndex2, 1, MPI_INT, MASTER, DISTRIBUTE_HIGHESTINDEX2_TAG, MPI_COMM_WORLD);
                }


                //the stored highest mutual information value is broadcasted back from the MASTER process to all other processes
                MPI_Bcast(&highestMutual, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
                MPI_Bcast(&highestType1, 1, MPI_INT, MASTER, MPI_COMM_WORLD); //including types
                MPI_Bcast(&highestType2, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                MPI_Bcast(&highestIndex1, 1, MPI_INT, MASTER, MPI_COMM_WORLD); // as well as indices of the degrees of freedom
                MPI_Bcast(&highestIndex2, 1, MPI_INT, MASTER, MPI_COMM_WORLD);



                totalMutual+=highestMutual; //update the total mutual information
                if(highestType1==b) { // if the "In" type was a bond
                    if(highestType2==b) { // and the "Out" type was also a bond
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_B][highestIndex1]+1<<" "<<Out[TYPE_B][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_B].push_back(Out[TYPE_B][highestIndex2]);
                        Out[TYPE_B].erase(Out[TYPE_B].begin()+highestIndex2); //put the the bond from the Out[TYPE_B] vector to the In[TYPE_B] vector
                        totalBbMutual+=highestMutual; //and update the sum of the bonds-bonds mutual information
                    }
                    if(highestType2==a) { //if  otherwise the "Out" type was an angle
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_B][highestIndex1]+1<<" "<<Out[TYPE_A][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_A].push_back(Out[TYPE_A][highestIndex2]);
                        Out[TYPE_A].erase(Out[TYPE_A].begin()+highestIndex2); //put the angle from the Out[TYPE_A] vector to the In[TYPE_A] vector
                        totalBaMutual+=highestMutual; //and update the sum of the bonds-angles mutual information
                    }
                    if(highestType2==d) { //analogously for Out[TYPE_D]
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_B][highestIndex1]+1<<" "<<Out[TYPE_D][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_D].push_back(Out[TYPE_D][highestIndex2]);
                        Out[TYPE_D].erase(Out[TYPE_D].begin()+highestIndex2);
                        totalBdMutual+=highestMutual;
                    }
                }
                if(highestType1==a) { //also for all In[TYPE_A] combinations
                    if(highestType2==b) {
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_A][highestIndex1]+1<<" "<<Out[TYPE_B][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_B].push_back(Out[TYPE_B][highestIndex2]);
                        Out[TYPE_B].erase(Out[TYPE_B].begin()+highestIndex2);
                        totalBaMutual+=highestMutual;
                    }
                    if(highestType2==a) {
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_A][highestIndex1]+1<<" "<<Out[TYPE_A][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_A].push_back(Out[TYPE_A][highestIndex2]);
                        Out[TYPE_A].erase(Out[TYPE_A].begin()+highestIndex2);
                        totalAaMutual+=highestMutual;
                    }
                    if(highestType2==d) {
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_A][highestIndex1]+1<<" "<<Out[TYPE_D][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_D].push_back(Out[TYPE_D][highestIndex2]);
                        Out[TYPE_D].erase(Out[TYPE_D].begin()+highestIndex2);
                        totalAdMutual+=highestMutual;
                    }
                }
                if(highestType1==d) { //and all In[TYPE_D] combinations
                    if(highestType2==b) {
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_D][highestIndex1]+1<<" "<<Out[TYPE_B][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_B].push_back(Out[TYPE_B][highestIndex2]);
                        Out[TYPE_B].erase(Out[TYPE_B].begin()+highestIndex2);
                        totalBdMutual+=highestMutual;
                    }
                    if(highestType2==a) {
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_D][highestIndex1]+1<<" "<<Out[TYPE_A][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_A].push_back(Out[TYPE_A][highestIndex2]);
                        Out[TYPE_A].erase(Out[TYPE_A].begin()+highestIndex2);
                        totalAdMutual+=highestMutual;
                    }
                    if(highestType2==d) {
                        if(rank==MASTER) {
                            outfile<<myChar[highestType1]<<" "<<myChar[highestType2]<<" "<<In[TYPE_D][highestIndex1]+1<<" "<<Out[TYPE_D][highestIndex2]+1<<"   MI: "<<highestMutual<<"   left to process: "<<Out[TYPE_B].size()+Out[TYPE_A].size()+Out[TYPE_D].size()<<endl;
                        }
                        In[TYPE_D].push_back(Out[TYPE_D][highestIndex2]);
                        Out[TYPE_D].erase(Out[TYPE_D].begin()+highestIndex2);
                        totalDdMutual+=highestMutual;
                    }
                }


            }
            #pragma omp barrier

        }
    }


    //after all degrees of freedom have been processed according to their mutual information terms, write out the endresult


    if(rank==MASTER) {
        infile.close();
        outfile<<"TOTAL 1D BONDS ENTROPY = "<<totalBondsEntropy<<endl;
        outfile<<"TOTAL 1D ANGLES ENTROPY = "<<totalAnglesEntropy<<endl;
        outfile<<"TOTAL 1D DIHEDRALS ENTROPY = "<<totalDihedralsEntropy<<endl;
        outfile<<"TOTAL 2D BONDS-BONDS MUTUAL INFORMATION = "<<totalBbMutual<<endl;
        outfile<<"TOTAL 2D BONDS-ANGLES MUTUAL INFORMATION = "<<totalBaMutual<<endl;
        outfile<<"TOTAL 2D BONDS-DIHEDRALS MUTUAL INFORMATION = "<<totalBdMutual<<endl;
        outfile<<"TOTAL 2D ANGLES-ANGLES MUTUAL INFORMATION = "<<totalAaMutual<<endl;
        outfile<<"TOTAL 2D ANGLES-DIHEDRALS MUTUAL INFORMATION = "<<totalAdMutual<<endl;
        outfile<<"TOTAL 2D DIHEDRALS-DIHEDRALS MUTUAL INFORMATION = "<<totalDdMutual<<endl;
        outfile<<"TOTAL CONFIGURATIONAL ENTROPY = "<<totalBondsEntropy+totalAnglesEntropy+totalDihedralsEntropy-totalMutual<<endl;
    }







    //Execution time is measured and MPI is unitialized
    cout<<"MPI PROCESS "<<rank<<" FINISHED."<<endl;  //<-------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==MASTER)
    {
        outfile.close();
        gettimeofday (&tv3, NULL);
        cout<<"Initialization + readin time: "<<tv2.tv_sec+0.000001 * tv2.tv_usec-tv1.tv_sec-0.000001 * tv1.tv_usec<<endl;
        cout<<"Calculation time: "<<tv3.tv_sec+0.000001 * tv3.tv_usec-tv2.tv_sec-0.000001 * tv2.tv_usec<<endl;
        cout<<"Total execution time: "<<tv3.tv_sec+0.000001 * tv3.tv_usec-tv1.tv_sec-0.000001 * tv1.tv_usec<<endl;
    }
    MPI_Finalize();
    if(rank==MASTER) {
        cout<<"PROGRAM FINISHED SUCCESSFULLY."<<endl;
    }

    return 0;
}







