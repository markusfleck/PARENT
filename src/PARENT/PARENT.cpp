//    PARENT, a parallel program to compute the configurational entropy from a molecular dynamics trajectory
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







#include <omp.h>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <sys/time.h>




#pragma pack(1)



#define MODFITNBINS 100

#define MASTER 0

#define DUMMY_TAG 0
#define COMMUNICATE_PARAMETERS_TAG 1
#define DISTRIBUTE_BONDS_TAG 2
#define DISTRIBUTE_ANGLES_TAG 3
#define DISTRIBUTE_DIHEDRALS_TAG 4
#define BONDS_MAXIMA_TAG 5
#define BONDS_MINIMA_TAG 6
#define ANGLES_MAXIMA_TAG 7
#define ANGLES_MINIMA_TAG 8
#define DIHEDRALS_MAXIMA_TAG 9
#define DIHEDRALS_MINIMA_TAG 10
#define BONDS_1D_ENTROPY_TAG 11
#define ANGLES_1D_ENTROPY_TAG 12
#define DIHEDRALS_1D_ENTROPY_TAG 13
#define BONDS_2D_ENTROPY_TAG 14
#define ANGLES_2D_ENTROPY_TAG 15
#define DIHEDRALS_2D_ENTROPY_TAG 16
#define BONDS_ANGLES_ENTROPY_TAG 17
#define BONDS_DIHEDRALS_ENTROPY_TAG 18
#define ANGLES_DIHEDRALS_ENTROPY_TAG 19


using namespace std;




#include "../util/io/io.h"








int main(int argc, char *argv[])  {
    //start the stopwatch for the execution time
    timeval tv1;
    gettimeofday (&tv1, NULL);
    timeval tv2,tv3,tv4;

    float floatdummyarray[3][3];
    double doubledummyarray[3];

    int numProcesses, rank, rc,len,threadLevelProvided, threadLevelClaimed;
    int tmpInt;

    double **bondsChunk, **anglesChunk, **dihedralsChunk;
    double *tmpBondsMaster, *tmpAnglesMaster, *tmpDihedralsMaster;
    double *tmpBonds, *tmpAngles, *tmpDihedrals;
    double *tmpArray, currentMin,currentMax;
    char hostName[MPI_MAX_PROCESSOR_NAME];


    const double pi=acos(-1);


    int nBonds,nAngles,nDihedrals,double_prec, numFrames;
    int fail=0;
    vector< vector <int> > dihedrals_top;
    vector <float> masses;
    vector <string> residues;
    vector <int> residueNumbers;
    vector <string> atomNames;
    vector <string> belongsToMolecule;
    double *maxBonds,*minBonds,*maxAngles,*minAngles,*maxDihedrals,*minDihedrals;
    double *bondsEntropy1D,*anglesEntropy1D,*dihedralsEntropy1D;
    double *bondsEntropy1DMaster,*anglesEntropy1DMaster,*dihedralsEntropy1DMaster;
    double *bbEntropy,*aaEntropy,*ddEntropy,*baEntropy,*bdEntropy,*adEntropy;
    int aDens,bDens,dDens,aDens1D,bDens1D,dDens1D, biggest,bigger;

    //Initialize MPI
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


    cout<<"MPI PROCESS: "<<rank<<", NODENAME: "<<hostName<<", MAXTHREADS: "<<omp_get_max_threads()<<", NUMTHREADS: "<<omp_get_num_procs()<<endl;
    MPI_Barrier(MPI_COMM_WORLD);


    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        printf("Initiated process = %i of NCPU=%i processes on machine=%s, thread number %i\n", rank, numProcesses, hostName ,ID);
    }


    //check if the command line arguments have been correctly supplied
    if(rank==MASTER) {
        cout<<endl;
        cout<<endl;
        cout<<endl;
        cout<<endl;
        cout<<endl;
        cout<<"    PARENT, a parallel program to compute the configurational entropy from a molecular dynamics trajectory"<<endl;
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
        if(argc==9) {
            string tmp1(argv[1]); //first argument is thee .bat trajectory file
            string tmp2(argv[2]); //second argument is the .par output file
            char *ptr,*type1,*type2;
            char delimiter[] = ".";

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
            if((strcmp(type1,"bat"))||(strcmp(type2,"par"))) {
                cerr<<"USAGE: "<<argv[0]<<" input.bat entropy.par #bondsbins1D #anglesbins1D #dihedralsbins1D #bondsbins2D #anglesbins2D #dihedralsbins2D\n";    //check for the extensions of the input and output file
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            if(sscanf(argv[3],"%d",&bDens1D)!=1) {
                cerr<<"ERROR: Could not read number of bins for 1D bonds from command line! Aborting"<<endl;    //read the number of bins for the 1D-bond calculation and check for correctness
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            if(sscanf(argv[4],"%d",&aDens1D)!=1) {
                cerr<<"ERROR: Could not read number of bins for 1D angles from command line! Aborting"<<endl;    //read the number of bins for the 1D-angle calculation and check for correctness
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            if(sscanf(argv[5],"%d",&dDens1D)!=1) {
                cerr<<"ERROR: Could not read number of bins for 1D dihedrals from command line! Aborting"<<endl;    //read the number of bins for the 1D-dihedral calculation and check for correctness
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            if(sscanf(argv[6],"%d",&bDens)!=1) {
                cerr<<"ERROR: Could not read number of bins (in one dimension => sqrt(nBins2D)) for 2D bonds from command line! Aborting"<<endl;    //read the number of bins for the 2D-bond calculation and check for correctness
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            if(sscanf(argv[7],"%d",&aDens)!=1) {
                cerr<<"ERROR: Could not read number of bins (in one dimension => sqrt(nBins2D)) for 2D angles from command line! Aborting"<<endl;    //read the number of bins for the 2D-angle calculation and check for correctness
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            if(sscanf(argv[8],"%d",&dDens)!=1) {
                cerr<<"ERROR: Could not read number of bins (in one dimension => sqrt(nBins2D)) for 2D dihedrals from command line! Aborting"<<endl;    //read the number of bins for the 2D-dihedral calculation and check for correctness
                MPI_Abort(MPI_COMM_WORLD, rc);
            }

            //determine the highest (and second highest) number of 2D bins used to avoid reallocations later on
            bigger=dDens>aDens?dDens:aDens;
            biggest=bDens>bigger?bDens:bigger;

        }
        else {
            cerr<<"USAGE: "<<argv[0]<<" input.bat entropy.par #bondsbins1D #anglesbins1D #dihedralsbins1D #bondsbins2D #anglesbins2D #dihedralsbins2D\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    //communicate these parameters to all MPI processes
    MPI_Bcast( &bDens, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &aDens, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &dDens, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &bDens1D, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &aDens1D, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &dDens1D, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &bigger, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &biggest, 1, MPI_INT, MASTER, MPI_COMM_WORLD);



    ifstream infile;
    ofstream par_file;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==MASTER) {
        cout<<"Reading file "<<argv[1]<<" .\n";
        //master process opens the input/output files
        infile.open(argv[1], ios::binary | ios::in );
        par_file.open(argv[2],ios::binary | ios::out);
        if(!(infile.is_open())) {
            cerr<<"ERROR: Could not open file "<<argv[1]<<" ! Aborting."<<endl;
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        if(!(par_file.is_open())) {
            cerr<<"ERROR: Could not open file "<<argv[2]<<" ! Aborting."<<endl;
            MPI_Abort(MPI_COMM_WORLD, rc);
        }

        //and reads the header of the trajectory
        if(read_BAT_header(&infile,&double_prec,&numFrames,&dihedrals_top,&masses,&residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) {
            cerr<<"AN ERROR HAS OCCURED WHILE READING THE HEADER OF THE FILE " <<argv[1]<<" . QUITTING PROGRAM.\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        nDihedrals=dihedrals_top.size();
        cout<<argv[1]<<" specs:"<<endl;
        cout<<"Precision: "<<(double_prec==1?"double":"single")<<" #Atoms: "<<nDihedrals+3<<" #Frames: "<<numFrames<<endl;  // ---------------------------------------------------

        //and writes the header of the output (.par) file
        if(write_PAR_header(&par_file,nDihedrals,double_prec,numFrames,&dihedrals_top, &masses,bDens1D,aDens1D,dDens1D,bDens,aDens,dDens, &residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) {
            cerr<<"AN ERROR HAS OCCURED WHILE WRITING THE HEADER OF THE FILE " <<argv[2]<<" . QUITTING PROGRAM.\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    //~ numFrames=10;   //<------------------------------------------------- uncomment for developing purposes

    //broadcast the header information to all MPI processes
    MPI_Bcast( &nDihedrals, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &double_prec, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &numFrames, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    nAngles=nDihedrals+1;
    nBonds=nDihedrals+2;

    //Allocate some arrays
    tmpArray=(double*)malloc(numFrames*sizeof(double));
    if(tmpArray==NULL) {
        cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    if(rank==MASTER) {
        //on the master side to store the 1D entropy results
        bondsEntropy1DMaster=(double*)malloc(nBonds*sizeof(double));
        anglesEntropy1DMaster=(double*)malloc(nAngles*sizeof(double));
        dihedralsEntropy1DMaster=(double*)malloc(nDihedrals*sizeof(double));

        //the 2D entropy results
        bbEntropy=(double*)malloc(nBonds*(nBonds-1)/2*sizeof(double));
        aaEntropy=(double*)malloc(nAngles*(nAngles-1)/2*sizeof(double));
        ddEntropy=(double*)malloc(nDihedrals*(nDihedrals-1)/2*sizeof(double));
        baEntropy=(double*)malloc(nBonds*nAngles*sizeof(double));
        bdEntropy=(double*)malloc(nBonds*nDihedrals*sizeof(double));
        adEntropy=(double*)malloc(nAngles*nDihedrals*sizeof(double));

        //and storage for readin
        tmpBondsMaster=(double*)malloc(nBonds*sizeof(double));
        tmpAnglesMaster=(double*)malloc(nAngles*sizeof(double));
        tmpDihedralsMaster=(double*)malloc(nDihedrals*sizeof(double));
        if((bondsEntropy1DMaster==NULL)||(anglesEntropy1DMaster==NULL)||(dihedralsEntropy1DMaster==NULL)||(bbEntropy==NULL)||(baEntropy==NULL)||(bdEntropy==NULL)||(aaEntropy==NULL)||(adEntropy==NULL)||(ddEntropy==NULL)||(tmpBondsMaster==NULL)||(tmpAnglesMaster==NULL)||(tmpDihedralsMaster==NULL)) {
            cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }



    //calculate fair splitting of the BAT coordinates for the MPI processes and allocate the necessary arrays
    int bondsBegin[nBonds], bondsEnd[nBonds], bondsLength[nBonds];
    tmpInt=nBonds-(nBonds/numProcesses)*numProcesses;  //calculate the rmainder of bonds when distributed equally to the MPI processes
    for(int i=0; i<numProcesses; i++) {
        bondsBegin[i]=i*(nBonds/numProcesses); //assign each process the bond number of the first bond it holds
        bondsEnd[i]=(i+1)*(nBonds/numProcesses)-1;//and the bond number of the last bond it holds

        bondsBegin[i]+=i<tmpInt?i:tmpInt; //take care of the remainder bonds (lower processes get an additional bond to calulate)
        bondsEnd[i]+=i<tmpInt?i+1:tmpInt;
        bondsLength[i]=bondsEnd[i]-bondsBegin[i]+1; //calculate how many bonds were assigned to the process
    }
    //every process allocates a 2D-array to store the trajectory if its bonds
    fail=0;
    bondsChunk=(double**)malloc(bondsLength[rank]*sizeof(double*));
    fail=fail||(bondsChunk==NULL);
    for(int i=0; i<bondsLength[rank]; i++) {
        bondsChunk[i]=(double*)calloc(numFrames,sizeof(double));
        fail=fail||(bondsChunk[i]==NULL);
    }
    tmpBonds=(double*)malloc(bondsLength[rank]*sizeof(double)); //and array to store one frame of its bonds
    fail=fail||(tmpBonds==NULL);
    bondsEntropy1D=(double*)malloc(bondsLength[rank]*sizeof(double)); //and an array to store the entropy result of its bonds
    fail=fail||(bondsEntropy1D==NULL);
    if(fail!=0) {
        cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    //proceed analogously for the angles
    int anglesBegin[nAngles], anglesEnd[nAngles], anglesLength[nAngles];
    tmpInt=nAngles-(nAngles/numProcesses)*numProcesses;
    for(int i=0; i<numProcesses; i++) {
        anglesBegin[i]=i*(nAngles/numProcesses);
        anglesEnd[i]=(i+1)*(nAngles/numProcesses)-1;

        anglesBegin[i]+=i<tmpInt?i:tmpInt;
        anglesEnd[i]+=i<tmpInt?i+1:tmpInt;
        anglesLength[i]=anglesEnd[i]-anglesBegin[i]+1;
    }
    fail=0;
    anglesChunk=(double**)malloc(anglesLength[rank]*sizeof(double*));
    fail=fail||(anglesChunk==NULL);
    for(int i=0; i<anglesLength[rank]; i++) {
        anglesChunk[i]=(double*)calloc(numFrames,sizeof(double));
        fail=fail||(anglesChunk[i]==NULL);
    }
    tmpAngles=(double*)malloc(anglesLength[rank]*sizeof(double));
    fail=fail||(tmpAngles==NULL);
    anglesEntropy1D=(double*)malloc(anglesLength[rank]*sizeof(double));
    fail=fail||(anglesEntropy1D==NULL);
    if(fail!=0) {
        cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }


    //proceed and analogously for the dihedrals
    fail=0;
    int dihedralsBegin[nDihedrals], dihedralsEnd[nDihedrals], dihedralsLength[nDihedrals];
    tmpInt=nDihedrals-(nDihedrals/numProcesses)*numProcesses;
    for(int i=0; i<numProcesses; i++) {
        dihedralsBegin[i]=i*(nDihedrals/numProcesses);
        dihedralsEnd[i]=(i+1)*(nDihedrals/numProcesses)-1;

        dihedralsBegin[i]+=i<tmpInt?i:tmpInt;
        dihedralsEnd[i]+=i<tmpInt?i+1:tmpInt;
        dihedralsLength[i]=dihedralsEnd[i]-dihedralsBegin[i]+1;
    }
    dihedralsChunk=(double**)malloc(dihedralsLength[rank]*sizeof(double*));
    fail=fail||(dihedralsChunk==NULL);
    for(int i=0; i<dihedralsLength[rank]; i++) {
        dihedralsChunk[i]=(double*)calloc(numFrames,sizeof(double));
        fail=fail||(dihedralsChunk[i]==NULL);
    }
    tmpDihedrals=(double*)malloc(dihedralsLength[rank]*sizeof(double));
    fail=fail||(tmpDihedrals==NULL);
    dihedralsEntropy1D=(double*)malloc(dihedralsLength[rank]*sizeof(double));
    fail=fail||(dihedralsEntropy1D==NULL);
    if(fail!=0) {
        cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }


    //allocate and initialize the extrema array for the degrees of freedom,
    maxBonds=(double*)malloc(bondsLength[rank]*sizeof(double));
    maxAngles=(double*)malloc(anglesLength[rank]*sizeof(double));
    maxDihedrals=(double*)malloc(dihedralsLength[rank]*sizeof(double));
    minBonds=(double*)malloc(bondsLength[rank]*sizeof(double));
    minAngles=(double*)malloc(anglesLength[rank]*sizeof(double));
    minDihedrals=(double*)malloc(dihedralsLength[rank]*sizeof(double));
    if((minBonds==NULL)||(minAngles==NULL)||(minDihedrals==NULL)||(maxBonds==NULL)||(maxAngles==NULL)||(maxDihedrals==NULL)) {
        cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }


    for(int i=0; i<bondsLength[rank]; i++) {
        maxBonds[i]=-1.0e100;
        minBonds[i]=1.0e100;
    }
    for(int i=0; i<anglesLength[rank]; i++) {
        maxAngles[i]=-1.0e100;
        minAngles[i]=1.0e100;
    }
    for(int i=0; i<dihedralsLength[rank]; i++) {
        maxDihedrals[i]=-1.0e100;
        minDihedrals[i]=1.0e100;
    }

    if(rank==MASTER) {
        //master reads the BAT-coordinate trajectory and spreads it to the MPI processes.
        for(int i=0; i<numFrames; i++) { //for all frames in the trajectory
            //~ if(read_BAT_frame(&infile,double_prec, tmpBondsMaster, tmpAnglesMaster, tmpDihedralsMaster, nDihedrals)!=0) { //read the frame and react to errors
            if(read_BAT_frame(&infile,double_prec, nDihedrals, floatdummyarray[0], floatdummyarray[0], (float**)floatdummyarray, doubledummyarray,doubledummyarray,doubledummyarray,doubledummyarray,tmpBondsMaster, tmpAnglesMaster, tmpDihedralsMaster) !=0) { //read the frame and react to errors
                cerr<<"An ERROR has occurred while reading the file " <<argv[1]<<" . Quitting Program.\n";
                MPI_Abort(MPI_COMM_WORLD, rc);
            }

            if(i%10000==0) {
                cout<<"Reading frame "<<i<<" and the following.\n";   //every 10000 frames issue an information to stdout
            }
            for(int j=0; j<numProcesses; j++) {
                if(j!=MASTER) //then send the all the MPI processes their according chunks of the degrees of freedom of the current frame
                {   MPI_Send(&(tmpBondsMaster[bondsBegin[j]]), bondsLength[j], MPI_DOUBLE, j, DISTRIBUTE_BONDS_TAG, MPI_COMM_WORLD);
                    MPI_Send(&(tmpAnglesMaster[anglesBegin[j]]), anglesLength[j], MPI_DOUBLE, j, DISTRIBUTE_ANGLES_TAG, MPI_COMM_WORLD);
                    MPI_Send(&(tmpDihedralsMaster[dihedralsBegin[j]]), dihedralsLength[j], MPI_DOUBLE, j, DISTRIBUTE_DIHEDRALS_TAG, MPI_COMM_WORLD);
                }
                else //and copy the masters own chunk to its local trajectory storage array
                {   for(int k=bondsBegin[j]; k<=bondsEnd[j]; k++) {
                        bondsChunk[k-bondsBegin[j]][i]=tmpBondsMaster[k];
                    }
                    for(int k=anglesBegin[j]; k<=anglesEnd[j]; k++) {
                        anglesChunk[k-anglesBegin[j]][i]=tmpAnglesMaster[k];
                    }
                    for(int k=dihedralsBegin[j]; k<=dihedralsEnd[j]; k++) {
                        dihedralsChunk[k-dihedralsBegin[j]][i]=tmpDihedralsMaster[k];
                    }
                }
            }
        }
        infile.close();
    }
    else //the non-master processes
    {
        for(int i=0; i<numFrames; i++) { //for every frame
            MPI_Recv(&(tmpBonds[0]), bondsLength[rank], MPI_DOUBLE, MASTER, DISTRIBUTE_BONDS_TAG, MPI_COMM_WORLD, &Stat); //receive their chunk of bonds
            MPI_Recv(&(tmpAngles[0]), anglesLength[rank], MPI_DOUBLE, MASTER, DISTRIBUTE_ANGLES_TAG, MPI_COMM_WORLD, &Stat);// angles
            MPI_Recv(&(tmpDihedrals[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DISTRIBUTE_DIHEDRALS_TAG, MPI_COMM_WORLD, &Stat);// and dihedrals
            //and add the received chunk of the frame to their local trajectory storage array
            for(int j=0; j<bondsLength[rank]; j++) {
                bondsChunk[j][i]=tmpBonds[j];
            }
            for(int j=0; j<anglesLength[rank]; j++) {
                anglesChunk[j][i]=tmpAngles[j];
            }
            for(int j=0; j<dihedralsLength[rank]; j++) {
                dihedralsChunk[j][i]=tmpDihedrals[j];
            }
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==MASTER) {
        gettimeofday (&tv2, NULL);   //measure the time for the readin routine
    }
    cout<<"READIN AND DISTRIBUTION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;


    //MPI processes calculate the extrema using OpenMP (=threaded) and apply modfitting (shifting to get the least possible non-zero span and thus the best resolution) to the dihedrals
    #pragma omp parallel
    {
        double tmpMin,tmpMax, modFit,binsize;
        int longestZeroStretch,longestZeroStretchPos,currentZeroStretch,currentZeroStretchPos;
        bool zeroExists;
        long long int histo[MODFITNBINS];
        #pragma omp for
        for(int j=0; j<bondsLength[rank]; j++) { //for all bonds of the MPI process (using threads)
            tmpMax=-1.0e100;
            tmpMin=1.0e100;
            for(int i=0; i<numFrames; i++) { //and all frames
                if(bondsChunk[j][i]>tmpMax) {
                    tmpMax=bondsChunk[j][i];
                }
                if(bondsChunk[j][i]<tmpMin) {
                    tmpMin=bondsChunk[j][i];   //find the maximum and minmum values
                }
            }
            if((tmpMin<0.0)||(tmpMax<0.0)) {
                cerr<<"ERROR: Bond "<<j<<" is smaller than 0.0"<<endl;
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            tmpMin-=0.000000005;//and increase the boundaries a tiny bit
            tmpMax+=0.000000005;
            if(tmpMin<0.0) {
                tmpMin=0.0;
            }
            if ((tmpMax-tmpMin)<1.0e-4) {
                tmpMax+=0.05;
            }
            maxBonds[j]=tmpMax;//and store the calculated values
            minBonds[j]=tmpMin;
        }
        #pragma omp for
        for(int j=0; j<anglesLength[rank]; j++) { //for the angles proceed analogously
            tmpMax=-1.0e100;
            tmpMin=1.0e100;
            for(int i=0; i<numFrames; i++) {
                if(anglesChunk[j][i]>tmpMax) {
                    tmpMax=anglesChunk[j][i];
                }
                if(anglesChunk[j][i]<tmpMin) {
                    tmpMin=anglesChunk[j][i];
                }
            }
            tmpMin-=0.000000005;
            tmpMax+=0.000000005;
            if (tmpMin<0) {
                tmpMin=0;
            }
            if (tmpMax>pi) {
                tmpMax=pi;
            }
            if ((tmpMax-tmpMin)<1.0e-4) {
                tmpMax+=0.05;
            }
            maxAngles[j]=tmpMax;
            minAngles[j]=tmpMin;
        }
        #pragma omp for
        for(int j=0; j<dihedralsLength[rank]; j++) { //for all dihedrals of the MPI process (using threads)
            //first build a histogram of the dihedral values over the trajectory
            for(int k=0; k<MODFITNBINS; k++) {
                histo[k]=0;
            }
            binsize=(2*pi+0.000000005)/MODFITNBINS;
            for(int i=0; i<numFrames; i++) {
                histo[int((dihedralsChunk[j][i])/binsize)]+=1;
            }
            zeroExists=false;
            for(int k=0; k<MODFITNBINS; k++) {
                zeroExists=zeroExists||(histo[k]==0);
            }
            if(zeroExists) { //if any of the bins of the histogram is empty find the longest consecutive stretch of  empty bins
                longestZeroStretch=0;
                currentZeroStretch=0;
                longestZeroStretchPos=-1;
                for(int k=0; k<2*MODFITNBINS; k++) { //for all bins of the histogram
                    int l=k%MODFITNBINS; //taking car of zero stretches which span the boundaries
                    if((currentZeroStretch==0)&&(histo[l]==0)) { //find and save a beginning zero stretch
                        currentZeroStretch=1;
                        currentZeroStretchPos=k;
                    }
                    if((currentZeroStretch>0)&&(histo[l]==0)) {
                        currentZeroStretch+=1;
                    }
                    if((currentZeroStretch>0)&&(histo[l]!=0)) { //and the end of it. If it is currently the longest zero stretch, save it
                        if(currentZeroStretch>longestZeroStretch) {
                            longestZeroStretch=currentZeroStretch;
                            longestZeroStretchPos=currentZeroStretchPos;
                        }
                        currentZeroStretch=0;
                    }
                }
            }
            else { //if none of the bins is empty
                longestZeroStretchPos=0;  //misuse the zeroStretch variables for determining the minimum
                longestZeroStretch=histo[0];
                for(int k=0; k<MODFITNBINS; k++) {
                    if(histo[k]<longestZeroStretch) {
                        longestZeroStretch=histo[k];
                        longestZeroStretchPos=k;
                    }
                }
            }
            modFit=2*pi-(longestZeroStretchPos+0.5)*binsize;//calculate the shift to put the zero stretch to the 2pi end
            for(int k=0; k<numFrames; k++) {
                dihedralsChunk[j][k]=dihedralsChunk[j][k]+modFit-2*pi*int((dihedralsChunk[j][k]+modFit)/(2*pi));   //and apply it taking care of circularity
            }

            tmpMax=-1.0e100;//then proceed determining the extrema analogously to bonds and angles
            tmpMin=1.0e100;
            for(int i=0; i<numFrames; i++) {
                if(dihedralsChunk[j][i]>tmpMax) {
                    tmpMax=dihedralsChunk[j][i];
                }
                if(dihedralsChunk[j][i]<tmpMin) {
                    tmpMin=dihedralsChunk[j][i];
                }
            }
            tmpMin-=0.000000005;
            tmpMax+=0.000000005;
            if(tmpMin<0.0) {
                tmpMin=0.0;
            }
            if(tmpMax>2*pi) {
                tmpMax=2*pi;
            }
            if ((tmpMax-tmpMin)<1.0e-4) {
                tmpMax+=0.05;
            }
            maxDihedrals[j]=tmpMax;
            minDihedrals[j]=tmpMin;
        }
    }


    cout<<"CALCULATION OF EXTREMA COMPLETED FOR RANK "<<rank<<"."<<endl;


    //calculate bonds entropy in 1D
    #pragma omp parallel
    {   double binsize,probDens,blen,plnpsum;
        long long int histo[bDens1D];
        int occupbins;
        #pragma omp for
        for(int j=0; j<bondsLength[rank]; j++) { //for all bonds of the MPI process (using threads)
            for(int k=0; k<bDens1D; k++) {
                histo[k]=0;   //initialize a histogram with zeros
            }
            binsize=(maxBonds[j]-minBonds[j])/(bDens1D); //and calculate the size of the bins
            for(int i=0; i<numFrames; i++) { // and fill the histogram using all frames of the trajectory
                histo[int((bondsChunk[j][i]-minBonds[j])/binsize)]+=1;
            }
            occupbins=0; //then use the histogram to calculate the (discretized) entropy, taking care of the Jacobian
            plnpsum=0;
            blen = minBonds[j]+(binsize/2.0);
            for(int k=0; k<bDens1D; k++) {
                probDens = histo[k]/(numFrames*binsize*blen*blen);
                if (probDens>0) {
                    plnpsum = plnpsum + blen*blen*probDens*log(probDens);
                    occupbins = occupbins + 1;
                }
                blen+=binsize;
            }
            plnpsum=-plnpsum*binsize;
            bondsEntropy1D[j]=plnpsum+(occupbins-1.0)/(2.0*numFrames); //and apply Herzel entropy unbiasing
        }
    }

    //send the results to the master MPI process, which stores it in bondsEntropy1DMaster
    if(rank!=MASTER) {
        MPI_Send(&(bondsEntropy1D[0]), bondsLength[rank], MPI_DOUBLE, MASTER, BONDS_1D_ENTROPY_TAG, MPI_COMM_WORLD);
    }
    else {
        memcpy(&(bondsEntropy1DMaster[bondsBegin[MASTER]]),&(bondsEntropy1D[0]), bondsLength[MASTER]*sizeof(double));
        for(int j=0; j<numProcesses; j++)
        {
            if(j!=MASTER)MPI_Recv(&(bondsEntropy1DMaster[bondsBegin[j]]), bondsLength[j], MPI_DOUBLE, j, BONDS_1D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
        }
    }

    cout<<"1D BONDS CALCULATION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;

    //to calculate angles entropy in 1D, proceed analogously as for bonds
    #pragma omp parallel
    {   double binsize,probDens,theta,plnpsum;
        long long int histo[aDens1D];
        int occupbins;
        #pragma omp for
        for(int j=0; j<anglesLength[rank]; j++) {
            for(int k=0; k<aDens1D; k++) {
                histo[k]=0;
            }
            binsize=(maxAngles[j]-minAngles[j])/(aDens1D);
            for(int i=0; i<numFrames; i++) {
                histo[int((anglesChunk[j][i]-minAngles[j])/binsize)]+=1;
            }
            occupbins=0;
            plnpsum=0;
            theta = minAngles[j]+(binsize/2.0);
            for(int k=0; k<aDens1D; k++) {
                probDens = histo[k]/(numFrames*binsize*sin(theta));
                if (probDens>0) {
                    plnpsum = plnpsum + sin(theta)*probDens*log(probDens);
                    occupbins = occupbins + 1;
                }
                theta+=binsize;
            }
            plnpsum=-plnpsum*binsize;
            anglesEntropy1D[j]=plnpsum+(occupbins-1.0)/(2.0*numFrames);
        }
    }

    //and also send the results to the master MPI process, which stores it in anglesEntropy1DMaster
    if(rank!=MASTER) {
        MPI_Send(&(anglesEntropy1D[0]), anglesLength[rank], MPI_DOUBLE, MASTER, BONDS_1D_ENTROPY_TAG, MPI_COMM_WORLD);
    }
    else {
        memcpy(&(anglesEntropy1DMaster[anglesBegin[MASTER]]),&(anglesEntropy1D[0]), anglesLength[MASTER]*sizeof(double));
        for(int j=0; j<numProcesses; j++)
        {
            if(j!=MASTER)MPI_Recv(&(anglesEntropy1DMaster[anglesBegin[j]]), anglesLength[j], MPI_DOUBLE, j, BONDS_1D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
        }
    }

    cout<<"1D ANGLES CALCULATION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;

    //to calculate dihedrals entropy in 1D (using modfitting), proceed analogously as for bonds and angles
    #pragma omp parallel
    {   double binsize,probDens,plnpsum;
        long long int histo[dDens1D];
        int occupbins;
        #pragma omp for
        for(int j=0; j<dihedralsLength[rank]; j++) {
            for(int k=0; k<dDens1D; k++) {
                histo[k]=0;
            }
            binsize=(maxDihedrals[j]-minDihedrals[j])/(dDens1D);
            for(int i=0; i<numFrames; i++) {
                histo[int((dihedralsChunk[j][i]-minDihedrals[j])/binsize)]+=1;
            }
            occupbins=0;
            plnpsum=0;
            for(int k=0; k<dDens1D; k++) {
                probDens = histo[k]/(numFrames*binsize);
                if (probDens>0) {
                    plnpsum = plnpsum + probDens*log(probDens);
                    occupbins = occupbins + 1;
                }
            }
            plnpsum=-plnpsum*binsize;
            dihedralsEntropy1D[j]=plnpsum+(occupbins-1.0)/(2.0*numFrames); //Herzel entropy unbiasing
        }
    }

    //and also send the results to the master MPI process, which stores it in dihedralsEntropy1DMaster
    if(rank!=MASTER) {
        MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, BONDS_1D_ENTROPY_TAG, MPI_COMM_WORLD);
    }
    else {
        memcpy(&(dihedralsEntropy1DMaster[dihedralsBegin[MASTER]]),&(dihedralsEntropy1D[0]), dihedralsLength[MASTER]*sizeof(double));
        for(int j=0; j<numProcesses; j++)
        {
            if(j!=MASTER)MPI_Recv(&(dihedralsEntropy1DMaster[dihedralsBegin[j]]), dihedralsLength[j], MPI_DOUBLE, j, BONDS_1D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
        }
    }

    cout<<"1D DIHEDRALS CALCULATION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;

    //calculation of the 2-dim entropies with bonds (bonds-bonds, bonds-angles, bonds-dihedrals) starts here
    for(int i=0; i<numProcesses; i++) { // cycling through every MPI process
        for(int j=0; j<bondsLength[i]; j++) { // and through every bond local to the process
            if(rank==i) {
                for(int k=0; k<numProcesses; k++) {
                    if(k!=rank) { ////the current process sends out  to all other processes
                        MPI_Send(&(maxBonds[j]), 1, MPI_DOUBLE, k, BONDS_MAXIMA_TAG, MPI_COMM_WORLD);//the current bond's extrema
                        MPI_Send(&(minBonds[j]), 1, MPI_DOUBLE, k, BONDS_MINIMA_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(bondsChunk[j][0]), numFrames, MPI_DOUBLE, k, DISTRIBUTE_BONDS_TAG, MPI_COMM_WORLD);//and the trajectory of the current bond
                    }
                }
                //calculation of bonds-bonds entropy for the current process starts here (using openMP threads)
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,blen1,blen2,plnpsum,sum1;

                    long long int histo[bDens][biggest]; //the biggest needed histogram (bonds-bonds, bonds-angles, bonds-dihedrals) is allocated,
                    int occupbins;
                    #pragma omp for
                    for(int jOmp=j+1; jOmp<bondsLength[rank]; jOmp++) { //the yet uncalculated 2D bonds-bonds entropy calculation pairs are split among openMP threads
                        for(int kOmp=0; kOmp<bDens; kOmp++) { //and the histogram is initialized with zeros
                            for(int iOmp=0; iOmp<bDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        //calculate binsizes for both dimensions
                        binsize1=(maxBonds[j]-minBonds[j])/bDens;//the current bond of the MPI process
                        binsize2=(maxBonds[jOmp]-minBonds[jOmp])/bDens;//and the current bond from the openMP thread
                        for(int iOmp=0; iOmp<numFrames; iOmp++) { // then fill the 2D histogram
                            histo[int((bondsChunk[j][iOmp]-minBonds[j])/binsize1)][int((bondsChunk[jOmp][iOmp]-minBonds[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0; //and calculate the entropy analogously to the 1D case (applying the Jacobian correction as well as Herzel entropy unbiasing)
                        plnpsum=0;
                        blen1= minBonds[j]+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            blen2= minBonds[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<bDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1*blen2*blen2);
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*blen2*blen2*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                blen2 = blen2+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        bondsEntropy1D[bondsLength[rank]-jOmp+j]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);//bondsEntropy1D is recycled to temporarily store all bonds-bonds entropies with the bond=j, but in reverse order (1+j <= bondsLength[rank]-jOmp+j <= bondsLength[rank]-1, meaning highest jOmp => lowest position in the array)
                    }
                    //calculation of bonds-angles entropy for the current process starts here (using openMP threads)
                    double theta;
                    #pragma omp for
                    for(int jOmp=0; jOmp<anglesLength[rank]; jOmp++) { //2D bonds-angles entropy calculations (for the current bond of the MPI process) are split among openMP threads
                        for(int kOmp=0; kOmp<bDens; kOmp++) { //and the histogram is reinitialized to zeros
                            for(int iOmp=0; iOmp<aDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(maxBonds[j]-minBonds[j])/bDens;//binsize for the current bond of the MPI process
                        binsize2=(maxAngles[jOmp]-minAngles[jOmp])/aDens;//and the current angle from the openMP thread are calculated
                        for(int iOmp=0; iOmp<numFrames; iOmp++) { //the 2D histogram is filled
                            histo[int((bondsChunk[j][iOmp]-minBonds[j])/binsize1)][int((anglesChunk[jOmp][iOmp]-minAngles[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;//then the entropy is calculated analogously to the 1D cases (applying the Jacobian correction as well as Herzel entropy unbiasing)
                        plnpsum=0;
                        blen1= minBonds[j]+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            theta = minAngles[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<aDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1*sin(theta));
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*sin(theta)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                theta = theta+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        anglesEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                    //calculation of bonds-dihedrals entropy for the current process starts here (using openMP threads)
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) { //2D bonds-dihedrals entropy calculations (for the current bond of the MPI process) are split among openMP threads
                        for(int kOmp=0; kOmp<bDens; kOmp++) { //and the histogram is reinitialized to zeros
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(maxBonds[j]-minBonds[j])/bDens;//binsize for the current bond of the MPI process
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;//and the current dihedral from the openMP thread are calculated
                        for(int iOmp=0; iOmp<numFrames; iOmp++) { //the 2D histogram is filled
                            histo[int((bondsChunk[j][iOmp]-minBonds[j])/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;//then the entropy is calculated analogously to the 1D cases (applying the Jacobian correction as well as Herzel entropy unbiasing)
                        plnpsum=0;
                        blen1= minBonds[j]+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1);
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        dihedralsEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }

                //the results are sent to the master array (see the receiving part at the end of the bonds entropy calculation block for further explanation)
                if(rank!=MASTER) {
                    MPI_Send(&(bondsEntropy1D[j+1]), bondsLength[rank]-j-1, MPI_DOUBLE, MASTER, BONDS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(anglesEntropy1D[0]), anglesLength[rank], MPI_DOUBLE, MASTER, ANGLES_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(bbEntropy[(nBonds-(bondsBegin[i]+j)-2)*(nBonds-(bondsBegin[i]+j)-1)/2+nBonds-bondsBegin[MASTER]-bondsLength[MASTER]]),&(bondsEntropy1D[j+1]),(bondsLength[rank]-j-1)*sizeof(double));
                    memcpy(&(baEntropy[(bondsBegin[i]+j)*nAngles+anglesBegin[rank]]),&(anglesEntropy1D[0]),(anglesLength[rank])*sizeof(double));
                    memcpy(&(bdEntropy[(bondsBegin[i]+j)*nDihedrals+dihedralsBegin[rank]]),&(dihedralsEntropy1D[0]),(dihedralsLength[rank])*sizeof(double));
                }
            }
            else if(rank>i) { //if the rank of the calculating MPI process is bigger than the one of the process sending to current bond, 2D entropies with all local bond partners need to be calculated
                MPI_Recv(&currentMax, 1, MPI_DOUBLE, i, BONDS_MAXIMA_TAG, MPI_COMM_WORLD, &Stat);//extrema and trajectory of the current bond are received
                MPI_Recv(&currentMin, 1, MPI_DOUBLE, i, BONDS_MINIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&(tmpArray[0]), numFrames, MPI_DOUBLE, i, DISTRIBUTE_BONDS_TAG, MPI_COMM_WORLD, &Stat);
                //then calculate bonds-bonds entropy (using openMP threads) analogously to the (rank==i) case
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,blen1,blen2,plnpsum,sum1;
                    long long int histo[bDens][biggest];
                    int occupbins;
                    #pragma omp for
                    for(int jOmp=0; jOmp<bondsLength[rank]; jOmp++) { //if the rank of the calculating MPI process is bigger than the one of the process sending to current bond, 2D entropies with all local bond partners need to be calculated
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            for(int iOmp=0; iOmp<bDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/bDens;
                        binsize2=(maxBonds[jOmp]-minBonds[jOmp])/bDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((bondsChunk[jOmp][iOmp]-minBonds[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        blen1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            blen2= minBonds[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<bDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1*blen2*blen2);
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*blen2*blen2*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                blen2 = blen2+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        bondsEntropy1D[bondsLength[rank]-jOmp-1]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                    //also calculate bonds-angles entropy (using openMP threads) analogously to the (rank==i) case
                    double theta;
                    #pragma omp for
                    for(int jOmp=0; jOmp<anglesLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            for(int iOmp=0; iOmp<aDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/bDens;
                        binsize2=(maxAngles[jOmp]-minAngles[jOmp])/aDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((anglesChunk[jOmp][iOmp]-minAngles[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        blen1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            theta = minAngles[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<aDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1*sin(theta));
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*sin(theta)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                theta = theta+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        anglesEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                    //and calculate bonds-dihedrals entropy (using openMP threads) analogously to the (rank==i) case
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/bDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        blen1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1);
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        dihedralsEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }
                //the results are sent to the master array (see the receiving part at the end of the bonds entropy calculation block for further explanation)
                if(rank!=MASTER) {
                    MPI_Send(&(bondsEntropy1D[0]), bondsLength[rank], MPI_DOUBLE, MASTER, BONDS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(anglesEntropy1D[0]), anglesLength[rank], MPI_DOUBLE, MASTER, ANGLES_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(bbEntropy[(nBonds-(bondsBegin[i]+j)-2)*(nBonds-(bondsBegin[i]+j)-1)/2+nBonds-bondsBegin[MASTER]-bondsLength[MASTER]]),&(bondsEntropy1D[0]),bondsLength[rank]*sizeof(double));
                    memcpy(&(baEntropy[(bondsBegin[i]+j)*nAngles+anglesBegin[rank]]),&(anglesEntropy1D[0]),(anglesLength[rank])*sizeof(double));
                    memcpy(&(bdEntropy[(bondsBegin[i]+j)*nDihedrals+dihedralsBegin[rank]]),&(dihedralsEntropy1D[0]),(dihedralsLength[rank])*sizeof(double));
                }
            }
            else if(rank<i) { //if the rank of the calculating MPI process is smaller than the one of the process sending to current bond, only bond-angle and bond-dihedral 2D entropies need to be calculated
                //everything else is done analogously to the (rank==i) and (rank>i) cases
                MPI_Recv(&currentMax, 1, MPI_DOUBLE, i, BONDS_MAXIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&currentMin, 1, MPI_DOUBLE, i, BONDS_MINIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&(tmpArray[0]), numFrames, MPI_DOUBLE, i, DISTRIBUTE_BONDS_TAG, MPI_COMM_WORLD, &Stat);
                //calculation of bonds-angles entropy for the current process starts here (using openMP threads)
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,blen1,plnpsum,sum1;
                    long long int histo[bDens][bigger];
                    int occupbins;
                    double theta;
                    #pragma omp for
                    for(int jOmp=0; jOmp<anglesLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            for(int iOmp=0; iOmp<aDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/bDens;
                        binsize2=(maxAngles[jOmp]-minAngles[jOmp])/aDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((anglesChunk[jOmp][iOmp]-minAngles[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        blen1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            theta = minAngles[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<aDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1*sin(theta));
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*sin(theta)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                theta = theta+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        anglesEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                    //calculation of bonds-dihedrals entropy for the current process starts here (using openMP threads)
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/bDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        blen1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<bDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*blen1*blen1);
                                if (probDens>0) {
                                    sum1 = sum1 + blen1*blen1*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                            blen1 = blen1+binsize1;
                        }
                        dihedralsEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }
                //the results are sent to the master array (see the receiving part right beneath here for further explanation)
                if(rank!=MASTER) {
                    MPI_Send(&(anglesEntropy1D[0]), anglesLength[rank], MPI_DOUBLE, MASTER, ANGLES_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(baEntropy[(bondsBegin[i]+j)*nAngles+anglesBegin[rank]]),&(anglesEntropy1D[0]),(anglesLength[rank])*sizeof(double));
                    memcpy(&(bdEntropy[(bondsBegin[i]+j)*nDihedrals+dihedralsBegin[rank]]),&(dihedralsEntropy1D[0]),(dihedralsLength[rank])*sizeof(double));
                }
            }
            if(rank==MASTER) {
                for(int k=0; k<numProcesses; k++) {
                    if(k!=MASTER) {
                        if(k>=i) { //only ranks >= i do bond-bond entropy calculations
                            MPI_Recv(&(bbEntropy[(nBonds-(bondsBegin[i]+j)-2)*(nBonds-(bondsBegin[i]+j)-1)/2+nBonds-bondsBegin[k]-bondsLength[k]]), (k==i?bondsLength[k]-j-1:bondsLength[k]), MPI_DOUBLE, k, BONDS_2D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
                            //The bbEntropy array is written in reverse style, like ([2][3]),([1][3]),([1][2]),([0][3]),([0][2]),([0][1]) for easier output purposes
                            //The first part of the formula (   (nBonds-(bondsBegin[i]+j)-2)*(nBonds-(bondsBegin[i]+j)-1)/2   ) is making use of the binomial coefficient to index the beginning of the bond with the number bondsBegin[i]+j in the array
                            //This works because the highest but one bond only has one partner to calculate, the highest but two has two and so on.
                            //Thus the highest but n bond starts at index 1+2+...+(n-2)+(n-1)=binomial(n,2)
                            //The rest of the formula takes into account that the results of the MPI process calculating the 2D entropy (of the bond with the number bondsBegin[i]+j) with the bond partners which are highest in numbering need to be placed lowest in the bbEntropy array (according to the reverse style convention)
                            //For the MPI process with the bonds highest in numbering, nBonds-bondsBegin[k]-bondsLength[k]==0, thus these results are placed directly at the indexing calculated for bondsBegin[i]+j, with no additional offset
                        }
                        MPI_Recv(&(baEntropy[(bondsBegin[i]+j)*nAngles+anglesBegin[k]]), anglesLength[k], MPI_DOUBLE, k, ANGLES_2D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);//bonds-angles are saved with the bond number varying slow and the angle number varying fast in the array
                        MPI_Recv(&(bdEntropy[(bondsBegin[i]+j)*nDihedrals+dihedralsBegin[k]]), dihedralsLength[k], MPI_DOUBLE, k, DIHEDRALS_2D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);//bonds-dihedrals are saved with the bond number varying slow and the dihedral number varying fast in the array
                    }
                }
            }
        }
    }

    cout<<"2D BONDS-BONDS, BONDS-ANGLES AND BONDS-DIHEDRALS CALCULATION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;






    //angles 2D entropies (angles-angles, angles-dihedrals) are calculated analogously to the bonds 2D entropies (angles-bonds has already been calculated with the bonds 2D entropies)
    for(int i=0; i<numProcesses; i++) {
        for(int j=0; j<anglesLength[i]; j++) {
            if(rank==i) {
                for(int k=0; k<numProcesses; k++) {
                    if(k!=rank) {
                        MPI_Send(&(maxAngles[j]), 1, MPI_DOUBLE, k, ANGLES_MAXIMA_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(minAngles[j]), 1, MPI_DOUBLE, k, ANGLES_MINIMA_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(anglesChunk[j][0]), numFrames, MPI_DOUBLE, k, DISTRIBUTE_ANGLES_TAG, MPI_COMM_WORLD);
                    }
                }
                //calculate angles-angles entropy using openMP threads (rank==i)
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,theta1,theta2,plnpsum,sum1;
                    long long int histo[aDens][bigger];
                    int occupbins;
                    #pragma omp for
                    for(int jOmp=j+1; jOmp<anglesLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            for(int iOmp=0; iOmp<aDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(maxAngles[j]-minAngles[j])/aDens;
                        binsize2=(maxAngles[jOmp]-minAngles[jOmp])/aDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((anglesChunk[j][iOmp]-minAngles[j])/binsize1)][int((anglesChunk[jOmp][iOmp]-minAngles[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        theta1= minAngles[j]+binsize1/2.0;
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            sum1 = 0.0;
                            theta2= minAngles[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<aDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*sin(theta1)*sin(theta2));
                                if (probDens>0) {
                                    sum1 = sum1 + sin(theta1)*sin(theta2)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                theta2 = theta2+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            theta1 = theta1+binsize1;
                        }
                        anglesEntropy1D[anglesLength[rank]-jOmp+j]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                    //calculate angles-dihedrals entropy using openMP threads (rank==i)
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(maxAngles[j]-minAngles[j])/aDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((anglesChunk[j][iOmp]-minAngles[j])/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        theta1= minAngles[j]+binsize1/2.0;
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*sin(theta1));
                                if (probDens>0) {
                                    sum1 = sum1 + sin(theta1)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                            theta1 = theta1+binsize1;
                        }
                        dihedralsEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }

                //send the result to the master array (rank==i)
                if(rank!=MASTER) {
                    MPI_Send(&(anglesEntropy1D[j+1]), anglesLength[rank]-j-1, MPI_DOUBLE, MASTER, ANGLES_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(aaEntropy[(nAngles-(anglesBegin[i]+j)-2)*(nAngles-(anglesBegin[i]+j)-1)/2+nAngles-anglesBegin[MASTER]-anglesLength[MASTER]]),&(anglesEntropy1D[j+1]),(anglesLength[rank]-j-1)*sizeof(double));
                    memcpy(&(adEntropy[(anglesBegin[i]+j)*nDihedrals+dihedralsBegin[rank]]),&(dihedralsEntropy1D[0]),(dihedralsLength[rank])*sizeof(double));
                }
            }
            else if(rank>i) {
                MPI_Recv(&currentMax, 1, MPI_DOUBLE, i, ANGLES_MAXIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&currentMin, 1, MPI_DOUBLE, i, ANGLES_MINIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&(tmpArray[0]), numFrames, MPI_DOUBLE, i, DISTRIBUTE_ANGLES_TAG, MPI_COMM_WORLD, &Stat);
                //calculate angles-angles entropy using openMP threads (rank>i)
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,theta1,theta2,plnpsum,sum1;
                    long long int histo[aDens][bigger];
                    int occupbins;
                    #pragma omp for
                    for(int jOmp=0; jOmp<anglesLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            for(int iOmp=0; iOmp<aDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/aDens;
                        binsize2=(maxAngles[jOmp]-minAngles[jOmp])/aDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((anglesChunk[jOmp][iOmp]-minAngles[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        theta1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            sum1 = 0.0;
                            theta2= minAngles[jOmp]+binsize2/2.0;
                            for(int lOmp=0; lOmp<aDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*sin(theta1)*sin(theta2));
                                if (probDens>0) {
                                    sum1 = sum1 + sin(theta1)*sin(theta2)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                                theta2 = theta2+binsize2;
                            }
                            plnpsum = plnpsum + sum1;
                            theta1 = theta1+binsize1;
                        }
                        anglesEntropy1D[anglesLength[rank]-jOmp-1]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                    //calculate angles-dihedrals entropy using openMP threads (rank>i)
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/aDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        theta1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*sin(theta1));
                                if (probDens>0) {
                                    sum1 = sum1 + sin(theta1)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                            theta1 = theta1+binsize1;
                        }
                        dihedralsEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }
                //send the result to the master array (rank>i)
                if(rank!=MASTER) {
                    MPI_Send(&(anglesEntropy1D[0]), anglesLength[rank], MPI_DOUBLE, MASTER, ANGLES_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(aaEntropy[(nAngles-(anglesBegin[i]+j)-2)*(nAngles-(anglesBegin[i]+j)-1)/2+nAngles-anglesBegin[MASTER]-anglesLength[MASTER]]),&(anglesEntropy1D[0]),anglesLength[rank]*sizeof(double));
                    memcpy(&(adEntropy[(anglesBegin[i]+j)*nDihedrals+dihedralsBegin[rank]]),&(dihedralsEntropy1D[0]),(dihedralsLength[rank])*sizeof(double));
                }
            }
            else if(rank<i) {
                MPI_Recv(&currentMax, 1, MPI_DOUBLE, i, ANGLES_MAXIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&currentMin, 1, MPI_DOUBLE, i, ANGLES_MINIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&(tmpArray[0]), numFrames, MPI_DOUBLE, i, DISTRIBUTE_ANGLES_TAG, MPI_COMM_WORLD, &Stat);
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,theta1,plnpsum,sum1;
                    long long int histo[aDens][dDens];
                    int occupbins;
                    //calculate angles-dihedrals entropy using openMP threads (rank<i)
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/aDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        theta1= currentMin+binsize1/2.0;
                        for(int kOmp=0; kOmp<aDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2*sin(theta1));
                                if (probDens>0) {
                                    sum1 = sum1 + sin(theta1)*probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                            theta1 = theta1+binsize1;
                        }
                        dihedralsEntropy1D[jOmp]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }
                //send the result to the master array (rank<i)
                if(rank!=MASTER) {
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(adEntropy[(anglesBegin[i]+j)*nDihedrals+dihedralsBegin[rank]]),&(dihedralsEntropy1D[0]),(dihedralsLength[rank])*sizeof(double));
                }
            }
            if(rank==MASTER) { //collect all results from the MPI processes
                for(int k=0; k<numProcesses; k++) {
                    if(k!=MASTER) {
                        if(k>=i) {
                            MPI_Recv(&(aaEntropy[(nAngles-(anglesBegin[i]+j)-2)*(nAngles-(anglesBegin[i]+j)-1)/2+nAngles-anglesBegin[k]-anglesLength[k]]), (k==i?anglesLength[k]-j-1:anglesLength[k]), MPI_DOUBLE, k, ANGLES_2D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
                        }
                        MPI_Recv(&(adEntropy[(anglesBegin[i]+j)*nDihedrals+dihedralsBegin[k]]), dihedralsLength[k], MPI_DOUBLE, k, DIHEDRALS_2D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
                    }
                }
            }
        }
    }

    cout<<"2D ANGLES-ANGLES AND ANGLES-DIHEDRALS CALCULATION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;








    //dihedrals 2D entropies (dihedrals-dihedrals) are calculated analogously to the previous cases (dihedrals-bonds and dihedrals-angles has already been calculated previously)
    for(int i=0; i<numProcesses; i++) {
        for(int j=0; j<dihedralsLength[i]; j++) {
            if(rank==i) {
                for(int k=0; k<numProcesses; k++) {
                    if(k>rank) {
                        MPI_Send(&(maxDihedrals[j]), 1, MPI_DOUBLE, k, DIHEDRALS_MAXIMA_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(minDihedrals[j]), 1, MPI_DOUBLE, k, DIHEDRALS_MINIMA_TAG, MPI_COMM_WORLD);
                        MPI_Send(&(dihedralsChunk[j][0]), numFrames, MPI_DOUBLE, k, DISTRIBUTE_DIHEDRALS_TAG, MPI_COMM_WORLD);
                    }
                }
                //calculate dihedrals-dihedrals entropy using openMP threads (rank==i)
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,plnpsum,sum1;
                    long long int histo[dDens][dDens];
                    int occupbins;
                    #pragma omp for
                    for(int jOmp=j+1; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<dDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(maxDihedrals[j]-minDihedrals[j])/dDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((dihedralsChunk[j][iOmp]-minDihedrals[j])/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        for(int kOmp=0; kOmp<dDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2);
                                if (probDens>0) {
                                    sum1 = sum1 + probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                        }
                        dihedralsEntropy1D[dihedralsLength[rank]-jOmp+j]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }

                //send the result to the master array (rank==i)
                if(rank!=MASTER) {
                    MPI_Send(&(dihedralsEntropy1D[j+1]), dihedralsLength[rank]-j-1, MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(ddEntropy[(nDihedrals-(dihedralsBegin[i]+j)-2)*(nDihedrals-(dihedralsBegin[i]+j)-1)/2+nDihedrals-dihedralsBegin[MASTER]-dihedralsLength[MASTER]]),&(dihedralsEntropy1D[j+1]),(dihedralsLength[rank]-j-1)*sizeof(double));
                }
            }
            else if(rank>i) {
                MPI_Recv(&currentMax, 1, MPI_DOUBLE, i, DIHEDRALS_MAXIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&currentMin, 1, MPI_DOUBLE, i, DIHEDRALS_MINIMA_TAG, MPI_COMM_WORLD, &Stat);
                MPI_Recv(&(tmpArray[0]), numFrames, MPI_DOUBLE, i, DISTRIBUTE_DIHEDRALS_TAG, MPI_COMM_WORLD, &Stat);
                //calculate dihedrals-dihedrals entropy using openMP threads (rank>i)
                #pragma omp parallel
                {
                    double binsize1,binsize2,probDens,plnpsum,sum1;
                    long long int histo[dDens][dDens];
                    int occupbins;
                    #pragma omp for
                    for(int jOmp=0; jOmp<dihedralsLength[rank]; jOmp++) {
                        for(int kOmp=0; kOmp<dDens; kOmp++) {
                            for(int iOmp=0; iOmp<dDens; iOmp++) {
                                histo[kOmp][iOmp]=0;
                            }
                        }
                        binsize1=(currentMax-currentMin)/dDens;
                        binsize2=(maxDihedrals[jOmp]-minDihedrals[jOmp])/dDens;
                        for(int iOmp=0; iOmp<numFrames; iOmp++) {
                            histo[int((tmpArray[iOmp]-currentMin)/binsize1)][int((dihedralsChunk[jOmp][iOmp]-minDihedrals[jOmp])/binsize2)]+=1;
                        }
                        occupbins=0;
                        plnpsum=0;
                        for(int kOmp=0; kOmp<dDens; kOmp++) {
                            sum1 = 0.0;
                            for(int lOmp=0; lOmp<dDens; lOmp++) {
                                probDens = histo[kOmp][lOmp]/(numFrames*binsize1*binsize2);
                                if (probDens>0) {
                                    sum1 = sum1 + probDens*log(probDens)*binsize1*binsize2;
                                    occupbins = occupbins + 1;
                                }
                            }
                            plnpsum = plnpsum + sum1;
                        }
                        dihedralsEntropy1D[dihedralsLength[rank]-jOmp-1]=-plnpsum+(occupbins-1.0)/(2.0*numFrames);
                    }
                }
                //send the result to the master array (rank>i)
                if(rank!=MASTER) {
                    MPI_Send(&(dihedralsEntropy1D[0]), dihedralsLength[rank], MPI_DOUBLE, MASTER, DIHEDRALS_2D_ENTROPY_TAG , MPI_COMM_WORLD);
                }
                else {
                    memcpy(&(ddEntropy[(nDihedrals-(dihedralsBegin[i]+j)-2)*(nDihedrals-(dihedralsBegin[i]+j)-1)/2+nDihedrals-dihedralsBegin[MASTER]-dihedralsLength[MASTER]]),&(dihedralsEntropy1D[0]),dihedralsLength[rank]*sizeof(double));
                }
            }
            if(rank==MASTER) { //collect all results from the MPI processes
                for(int k=0; k<numProcesses; k++) {
                    if((k!=MASTER)&&(k>=i)) {
                        MPI_Recv(&(ddEntropy[(nDihedrals-(dihedralsBegin[i]+j)-2)*(nDihedrals-(dihedralsBegin[i]+j)-1)/2+nDihedrals-dihedralsBegin[k]-dihedralsLength[k]]), (k==i?dihedralsLength[k]-j-1:dihedralsLength[k]), MPI_DOUBLE, k, DIHEDRALS_2D_ENTROPY_TAG, MPI_COMM_WORLD, &Stat);
                    }
                }
            }
        }
    }



    cout<<"2D DIHEDRALS-DIHEDRALS CALCULATION COMPLETED FOR MPI PROCESS "<<rank<<"."<<endl;



    //write out the resultsto the binary .par file and measure time
    if(rank==MASTER) {
        gettimeofday (&tv3, NULL);
        if(write_PAR_body(&par_file,nDihedrals,bondsEntropy1DMaster, anglesEntropy1DMaster, dihedralsEntropy1DMaster, bbEntropy, baEntropy, bdEntropy, aaEntropy, adEntropy, ddEntropy)!=0) {
            cerr<<"AN ERROR HAS OCCURED WHILE WRITING THE FILE " <<argv[2]<<" .\n";
        }
    }

//time is measured, timings are written out and MPI is unitialized
    cout<<"MPI PROCESS "<<rank<<" FINISHED."<<endl;  //<-------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==MASTER)
    {
        par_file.close();
        gettimeofday (&tv4, NULL);
        cout<<"Initialization + readin time: "<<tv2.tv_sec+0.000001 * tv2.tv_usec-tv1.tv_sec-0.000001 * tv1.tv_usec<<endl;
        cout<<"Calculation time: "<<tv3.tv_sec+0.000001 * tv3.tv_usec-tv2.tv_sec-0.000001 * tv2.tv_usec<<endl;
        cout<<"Output time: "<<tv4.tv_sec+0.000001 * tv4.tv_usec-tv3.tv_sec-0.000001 * tv3.tv_usec<<endl;
        cout<<"Total execution time: "<<tv4.tv_sec+0.000001 * tv4.tv_usec-tv1.tv_sec-0.000001 * tv1.tv_usec<<endl;
    }
    MPI_Finalize();
    if(rank==MASTER) {
        cout<<"PROGRAM FINISHED SUCCESSFULLY."<<endl;
    }
    return 0;
}


