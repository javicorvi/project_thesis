/*
SCPE : Structurally Constrained Protein Evolution 
 
Copyright (C) Gustavo Parisi and Julian Echave. 2001-2006
Contact: Gustavo Parisi: gustavo@unq.edu.ar

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

SCPEv1.1 allows the use of hete-oligomers 
*/
 
  
#include "stdio.h"
#include "stdlib.h" 
#include "string.h"
#include "math.h"
#include "time.h"
#include "scpe.h"
#include "random.h"
#include "ctype.h"

FILE *fin,*pdbin,*qmutfile;
FILE *logfile;
FILE *outnsus;
FILE *betalist;
FILE *ali_input;
FILE *list;
FILE *seqprint;

int GetOptions(FILE *fin);

struct info {
  char pdbin[1000],qmut[1000];
  char chain,outfile[MAXNAME];
  int monomer;
  int acceptmodel;
  int runs;
  int filterL,filterC;
  float accept,mut;
  float maxaccept,stepaccept;
  float initial_beta,step,maxbeta;
  float beta;
  char energy_eval[4];
  int screening;
  int multimer;
  char listname[MAXNAME];
  int list;
  char ali[MAXNAME];
  int num_diff_chain;
  char diff_chain[MAXCHAINS];
  char modelname[MAXNAME];
  char path[MAXNAME];
  int homo;
  int print_sequences;
  char print_sequences_output[1000];
  float print_sequences_identity;
  int print_sequences_screening;
  char print_sequences_contact_map_outputfile[1000];
}Options;

struct {
        char seq[MAXPROTLEN];
        char name[MAXNAME];
   }ali[MAXSEQS]; 

int *Nsus_pos,initial;
int *print_run_screen;
char seqprot[MAXPROTLEN];
char seqprot_total[MAXPROTLEN];
float eijmat[AA][AA];
float *x,*y,*z;
int*cmat,*cmatQ;
int xyznsites,protsize;
char xyzname[4];
int *numprot,*numprot_total,*numprot_totalA,*numprot_totalT;
double *qmatrix,*freq,*pmatrix,beta;
double *energyT, *energyA,score;
int *A,*T;
double mutrecord,acceptrecord;
double nonaccepted,accepted,nothing;
int nsyn,mut;
double N_nsyn,N_syn,N_nsynrecord;
char logfilename[MAXNAME];
int total_length;
double nsyn_persite[MAXPROTLEN],syn_persite[MAXPROTLEN];
int position;
double average;
float print_parameter;
int vector[MAXPROTLEN];
float number_of;
char outmatrix_name[MAXNAME];
char c,command[100];
int cmatvector[MAXPROTLEN][MAXCONTACTS],cmatQvector[MAXPROTLEN][MAXCONTACTS];
int fixedAA[MAXPROTLEN][MAXCONTACTS];
char seq_aa_print[10000];
int seq_num_print[10000];
float iden;
char mutated_secuence_file_name[500];


time_t  time_start,time_end;


#include "tools.c"
#include "matriz.c"


main (int argc, char * argv[])

{

  int i,j,k,m;
  double *P,*S;
  int *F,*F1;
  int *G;
  double lambda,lambda_actual,*R,*D;
  char sust_type;
  int nseq,newruns;
  int num=0,r;
  int cont; 



  time_start=time(NULL);
  SetSeed(time(NULL));


/*Set defaults values*/

   Options.initial_beta=0.15;
   Options.step=10.0;
   Options.mut=100000.0;
   Options.runs=14000;
   Options.maxbeta=40.0;
   Options.accept=20.0;
   Options.maxaccept=21.0;
   Options.stepaccept=5.0;
   strcpy(Options.qmut,"QijCodonesEmp");
   strcpy(Options.outfile,"SCPE");
   Options.multimer=0;
   Options.list=0;
   strcpy(Options.modelname,"scpe");
   strcpy(Options.path,"./OUT");
   strcpy(Options.pdbin,"inputfile.in");
   Options.print_sequences=1;
   strcpy(Options.print_sequences_output,"scpe_sequences/sequences");
   strcpy(Options.print_sequences_contact_map_outputfile,"contact_map.dat");
   Options.print_sequences_identity=25.0;
   Options.print_sequences_screening=1;
  /*Read Options from input file or command line*/

   if(argc==2)  {
              if((fin=fopen(argv[1],"r"))==NULL)
                {
                  fprintf(stderr,"SCPE input file requiered\n");
                  exit(1);
                }
              GetOptions(fin);
   }
    else {
       for( i=1; i < argc ; i++ ){
        if(argv[i][0] == '-' && (c=argv[i][1]) != '\0') 
         switch(c){

           case 'i': 
                    strcpy(Options.pdbin,argv[i+1]);
                    break;
           case 'c': 
                    Options.chain=argv[i+1][0];
                    break;
           case 's': 
                    strcpy(Options.ali,argv[i+1]);
                    break;
           case 'b': 
                    Options.initial_beta=atof(argv[i+1]);
                    break;
           case 'S': 
                    Options.step=atof(argv[i+1]);
                    break;
           case 'M': 
                    Options.maxbeta=atof(argv[i+1]);
                    break;
           case 'R': 
                    Options.runs=atoi(argv[i+1]);
                    break;
           case 'A': 
                    Options.accept=atof(argv[i+1]);
                    break;
           case 'a': 
                    Options.maxaccept=atof(argv[i+1]);
                    break;
           case 't': 
                    Options.stepaccept=atof(argv[i+1]);
                    break;
           case 'q': 
                    strcpy(Options.qmut,argv[i+1]);
                    break;
           case 'o': 
                    strcpy(Options.outfile,argv[i+1]);
                    break;
           case 'm': 
                    Options.multimer=argv[i+1][0];
                    break;
           case 'l': 
                    Options.list=argv[i+1][0];
                    break;
           case 'L': 
                    strcpy(Options.listname,argv[i+1]);
                    break;
           case 'p': 
                    strcpy(Options.path,argv[i+1]);
                    break;
           case 'n': 
                    strcpy(Options.modelname,argv[i+1]);
                    break;
           case '1':
                    Options.print_sequences=atoi(argv[i+1]);
                    break;
           case '2':
        	   	    strcpy(Options.print_sequences_output,argv[i+1]);
                    break;
           case '3':
                    Options.print_sequences_identity=atof(argv[i+1]);
                    break;
           case '4':
                    Options.print_sequences_screening=atoi(argv[i+1]);
                    break;
           case '5':
                  	strcpy(Options.print_sequences_contact_map_outputfile,argv[i+1]);
                    break;

        }
      }
    }
/*  sprintf(command,"mkdir -p %s\n",Options.path);
  system(command);*/


  if((pdbin=fopen(Options.pdbin,"r"))==NULL)
    {
      fprintf(stderr,"PDB input file requiered\n");
      exit(1);
    }
  if((qmutfile=fopen(Options.qmut,"r"))==NULL)
    {
      fprintf(stderr,"Protein Mutation Matrix file requiered\n");
      exit(1);
    }

 /*Open log file*/
 sprintf(logfilename,"%s.log",Options.outfile);
  if((logfile=fopen(logfilename,"w"))==NULL) {
    fprintf(stderr,"Can not write to log file\n");
    exit(1);
  }

  fflush(logfile);


  Readpdb(pdbin);

  printf("protsize =%d protsize total %d\n",protsize,total_length);
  Memory();

  Seq_to_number("one",seqprot,numprot);
  Seq_to_number("one",seqprot_total,numprot_total);

  memcpy(numprot_totalA,numprot_total,total_length*(sizeof(int)));
  memcpy(numprot_totalT,numprot_total,total_length*(sizeof(int)));

  if(Options.multimer){
       for (i=0; i<protsize; i++){
              cont=0; 
              for(j=0;j<total_length;j++) 
                 if(*(cmatQ+(i*total_length+j))==1){
                                fixedAA[i][cont]=numprot_total[j];
                                cont++;
                 }
       }
  }

       for (i=0; i<total_length; i++)printf("%d ",numprot_total[i]);

              printf("cmatQvector\n");
       for (i=0; i<protsize; i++){
              printf("%d ",i);
            for (j=0; j<MAXCONTACTS; j++)
              printf("%d ",cmatQvector[i][j]); 
              printf("\n");
       }
              printf("cmatvector\n");
       for (i=0; i<protsize; i++){
              printf("%d ",i);
            for (j=0; j<MAXCONTACTS; j++)
              printf("%d ",cmatvector[i][j]); 
              printf("\n");
       }


  Readq(qmutfile,qmatrix);

  diag(qmatrix,pmatrix,10000000.00); /*Equilibrium frequencies*/
  Calculate_freq(pmatrix,freq);
  Setrate(qmatrix,freq);
  lambda_actual= 1./protsize;
  for(i=0;i<AA;++i) for(j=0;j<AA;++j) *(qmatrix+(i*AA+j))= (*(qmatrix+(i*AA+j))*lambda_actual);
  Options.screening=1;/*old user  option*/
  if(!Options.screening) diag(qmatrix,pmatrix,0.15);
          else          diag(qmatrix,pmatrix,10.00);
  Setp_cumulative(pmatrix);
 

  set_eijmat(&eijmat[0][0]);    
 
 
   fclose(logfile);  
    
 if((list=fopen("list_of_matrices","w"))==NULL)
      {
      fprintf(stderr,"SCPE is unable to open list file\n");
      exit(1); 
 }
 
 if((ali_input=fopen(Options.ali,"r"))==NULL)
      {
      fprintf(stderr,"SCPE is running in single-sequence mode (no alignment file provided)\n",Options.ali);
      nseq=1;
      strcpy(ali[0].seq,seqprot);
      }
 else{
      nseq=ReadAlignment(ali_input);
      fclose (ali_input);
  }

      newruns=Options.runs/nseq;


  /*Begin evolution of sequence*/

  /*Screening different Nsus*/

 while(Options.accept<=Options.maxaccept){               

    printf("Evaluating Nsus %f \n ",Options.accept);

  if(Options.list){
    if((betalist=fopen(Options.listname,"r"))==NULL)
      {
      fprintf(stderr,"Beta list file (%s) requiered\n",Options.listname);
      exit(1);
      }
     i=fscanf(betalist,"%lf\n",&beta);
   }
   else beta=Options.initial_beta; 

  /*Screening different beta*/
  int id_seq = 0;
  int file_flush_flag = 0;
  while(beta<=Options.maxbeta){               

    printf("Evaluating cutoff %f\n ",beta);
    SetSeed(time(NULL));
    /*File For Print*/
    printf("file to print \n "); 


	//sprintf(mutated_secuence_file_name,"%s-beta%2.2f-nsus%2.2f-runs%d.fasta",
    //Options.print_sequences_output,beta,Options.accept,Options.runs);

	sprintf(mutated_secuence_file_name,"%s",Options.print_sequences_output,beta,Options.accept,Options.runs);

	printf("Print Sequences  %d\n ",Options.print_sequences);
	printf("Print Sequences File %s\n ",mutated_secuence_file_name);
	printf("Print Sequences Identity %f\n ",Options.print_sequences_identity);
	printf("Print Sequences Screening %d\n ",Options.print_sequences_screening);
	printf("Print Sequences Contact Map Output File %s\n ",Options.print_sequences_contact_map_outputfile);
	printf("Chain  %c\n ",Options.chain); 
	if(Options.print_sequences_screening==0){
		Options.print_sequences_screening = (protsize*0.40)/3;
		printf("Screening parameter is 0 so (protsize*0.40)/3 %d\n ", Options.print_sequences_screening);
	}
	fflush(stdout);
	
	if(Options.print_sequences==1){

		if((seqprint=fopen(mutated_secuence_file_name,"w"))==NULL) {
			fprintf(stderr,"Can not write to log seq print \n");
			printf("Can not write to log seq print \n");
			exit(1);
		}
	   //print reference sequence first
	   fprintf(seqprint,">SEQUENCE_REFERENCE\n%s\n",seqprot);
	   fflush(seqprint);
	   fclose(seqprint);
	}
   
   
	
    /*initialise A to ancestral sequence in each of the indep. run for a given beta*/
    /*In this version ancestrals sequence could be a set of sequences choosen from an alignment*/

    for(i=0;i<nseq;i++){
          Seq_to_number("one",ali[i].seq,numprot);
		  for(j=newruns*i;j<(newruns*(i+1));j++) memcpy((A+(j*protsize)),numprot,protsize*(sizeof(int)));
    }
    printf("Each of the %d sequences will have %d independent runs. Independent runs=%d \n",nseq,newruns,newruns*nseq);
    Options.runs=newruns*nseq;

    /*initialise counters of evolutionary change*/
    mutrecord=nonaccepted=accepted=nothing=acceptrecord=0;
    N_syn=N_nsyn=N_nsynrecord=0;
    for (i=0; i<protsize; i++) vector[i]=0; 
    for(i=0;i<protsize;i++){
                 nsyn_persite[i]=0;
                 syn_persite[i]=0;
    }
    num=average=0;

    Init_vector(Nsus_pos,(protsize*AA*AA)); /*Nsus is accumulative for N indep. runs for a given beta*/

	//inicializo vector de runs
    //printf("Inicializacion de Vector con  %i posiciones \n ",Options.runs-1);
    //fflush(stdout);
    print_run_screen=(int*)malloc(Options.runs*sizeof(int));
    Init_vector(print_run_screen,Options.runs);
    //printf("Inicializacion de Vector con  %i posiciones \n ",Options.runs-1);
    //fflush(stdout);
    /*for(int i=0;i<Options.runs;i++){
        	//printf(" run=%d %i \n",i,*(print_run_screen+i));
        	print_run_screen[i]=i;
    }*/
    /*for(int i=0;i<Options.runs;i++){
    	//printf(" run=%d %i \n",i,*(print_run_screen+i));
    	printf(" run=%d print_screen=%d \n",i,print_run_screen[i]);

    }*/
    //fflush(stdout);

    int print_step = 0;

    while(average<1){ 

      print_parameter=Options.accept;
      if(mutrecord>Options.mut) {
                 if((logfile=fopen(logfilename,"a"))==NULL) {
                          fprintf(stderr,"Can not write to log file\n");
                          exit(1);
                 }
              fprintf(logfile,"Maximum limit in mutations has been reached in beta %e\n",beta);
              fflush(logfile);
              fclose(logfile);
              break;
      }
      print_step++;
       /*loop over different runs*/
      for(j=0;j<Options.runs;j++){              

			/* sust_type indicates type of sustitution as follow:
			0 = no mutation
			N = not accepted  mutation
			A = accepted mutation
			*/

			sust_type= run((A+(j*protsize)),(T+(j*protsize)),&mut);

			//Si hubo mutacion

						if(sust_type=='A' && Options.print_sequences==1){
							print_run_screen[j]++;
							//Si hay que imprimir secuencias y el run actual debe imprimirse:
							//el print_run_screen es 0 significa la primera vez y en ese caso el contro lo haria la identidad,
							//si el print_run_screen !=0 debe ser igual al screenin
							if(print_run_screen[j]>=Options.print_sequences_screening){
									memcpy(seq_num_print,(A+(j*protsize)),protsize*(sizeof(int)));
									Number_to_seq(seq_num_print,seq_aa_print);
									id_seq++;
									iden=0;
									for(m=0;m<protsize;m++) if(numprot[m]==seq_num_print[m]) iden++;
									//Si la identidad es mayor a 25% segun parametro print_sequences_identity y la identidad es menor al 62
									//imprimo la secuencia
									if((iden/protsize)*100 > Options.print_sequences_identity && (iden/protsize)*100 < 62){
										/*imprimo secuencias*/
										if(file_flush_flag==0){
										   if((seqprint=fopen(mutated_secuence_file_name,"a+"))==NULL) {
											   fprintf(stderr,"Can not write to log file\n");
											   exit(1);
										   }
									   }
									   fprintf(seqprint,">SEQUENCE_run%d_c%d_%f_%d\n%s\n",j,print_run_screen[j],(iden/protsize)*100,id_seq,seq_aa_print);
									   print_run_screen[j]=0;
										file_flush_flag++;
										if(file_flush_flag==10){
											fflush(seqprint);
											fclose(seqprint);
											file_flush_flag=0;
										}
									}

								}//fin printrun >= que print

							}//fin resutl A



			#if defined(DEBUG)
					 printf(" run=%d %c \n",j,sust_type);
				#endif

			Count(sust_type,&nonaccepted,&accepted,&nothing,&mut,&nsyn,&N_nsyn,&N_syn);
			mutrecord=accepted+nonaccepted;

				if(average ==1.) break;
      }

   if(print_step==Options.print_sequences_screening){
	   print_step=0;
   } 



     
      /*N_syn, N_nsyn, N_nsynrecord, acceptrecord are accumulated over Options.runs and over all positions*/

      acceptrecord=accepted/(Options.runs * protsize);  /*per site*/
      mutrecord=mutrecord/(Options.runs * protsize);    /*per site*/
      N_nsynrecord=N_nsyn/(Options.runs * protsize);    /*per site*/

/*   for(i=0;i<protsize;i++) fprintf(stdout,"%f ",nsyn_persite[i]/Options.runs);
   fprintf(stdout,"average %f %f\n\n\n\n\n\n\n",average,print_cutoff[num]);

    fprintf(stdout,"acceptrecord = %f, mutrecord = %f omega = %f nonaccepted = %f ",acceptrecord, mutrecord,N_nsynrecord/mutrecord,nonaccepted);
    time_end=time(NULL);

    fprintf(stdout,"Time used: %f seconds \n",difftime(time_end,time_start));*/
   /*for(i=0;i<protsize;i++) fprintf(stdout,"%f\n",nsyn_persite[i]/(nsyn_persite[i]+syn_persite[i]));*/

    if(average == 1){
     /*Print output matrices*/
     sprintf(outmatrix_name,"%s-%2.2f-%2.2f-%2.2f.Nsus",Options.outfile,beta,N_nsynrecord/mutrecord,Options.accept);
     fprintf(list,"%s\n",outmatrix_name);
     matrix(Nsus_pos);


     }

   } /*End time*/
   
   for (i=0; i<protsize; i++) vector[i]=0;
   average=0;
   if(Options.list) fscanf(betalist,"%lf\n",&beta);
     else beta+=Options.step;

  } /*End run or begin another beta*/
  
   if(file_flush_flag<10){
     fflush(seqprint); 
     fclose(seqprint);
	 file_flush_flag=0;
   }
  
  Options.accept+=Options.stepaccept;
 
 }/*End run for a given Nsus*/

  fclose(list);

/*  sprintf(command,"mv *%s*  %s\n",outmatrix_name,Options.path);
  system(command);*/

  time_end=time(NULL);

  printf("Time used: %f seconds \n",difftime(time_end,time_start));

} 

/*Functions*/

int GetOptions(FILE *fin)
{
  char line[MAXSTR];
  char *s;
  char *optstr[] = {"PDB FILENAME","INITIAL BETA PARAMETER","SCREENING BETA-STEP","MAXIMUM BETA","NUMBER OF INDEPENDENT RUNS","INITIAL NUMBER OF ACCEPTED SUSTITUTIONS PER SITE","MAXIMUM NUMBER OF ACCEPTED SUSTITUTIONS PER SITE","STEP OF ACCEPTED SUSTITUTIONS PER SITE","MUTATION MATRIX","PROTEIN CHAIN","RESULTS FILENAME","MULTIMER PROTEIN MODE","BETA FROM LIST","LIST NAME","ALIGNMENT FILE","MODEL NAME","PATH TO RESULTS","HOMOOLIGOMER"};
  char comment='*'; 
  int noptions =18;
  int i,j;

 if((s=(char *)malloc((MAXSTR)*sizeof(char)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }


  while(fgets(line,MAXSTR,fin)!=NULL){
    if(line[0]!='*'){
      for(i=0;i<noptions;i++){
	if(strstr(line,optstr[i]))
	  switch (i) {
	  case ( 0): sscanf(line,"PDB FILENAME= %s",Options.pdbin); break;
	  case ( 1): sscanf(line,"INITIAL BETA PARAMETER= %f",&Options.initial_beta); break;
	  case ( 2): sscanf(line,"SCREENING BETA-STEP= %f",&Options.step); break;
	  case ( 3): sscanf(line,"MAXIMUM BETA= %f",&Options.maxbeta); break;
	  case ( 4): sscanf(line,"NUMBER OF INDEPENDENT RUNS= %d",&Options.runs); break;
	  case ( 5): sscanf(line,"INITIAL NUMBER OF ACCEPTED SUSTITUTIONS PER SITE= %f",&Options.accept); break;
	  case ( 6): sscanf(line,"MAXIMUM NUMBER OF ACCEPTED SUSTITUTIONS PER SITE= %f",&Options.maxaccept); break;
	  case ( 7): sscanf(line,"STEP OF ACCEPTED SUSTITUTIONS PER SITE= %f",&Options.stepaccept); break;
	  case ( 8): sscanf(line,"MUTATION MATRIX= %s",Options.qmut); break;
	  case ( 9): sscanf(line,"PROTEIN CHAIN= %c",&Options.chain); break;
	  case ( 10): sscanf(line,"RESULTS FILENAME = %s",Options.outfile); break;
	  case ( 11): sscanf(line,"MULTIMER PROTEIN MODE = %d",&Options.multimer); break;
	  case ( 12): sscanf(line,"BETA FROM LIST = %d",&Options.list); break;
	  case ( 13): sscanf(line,"LIST NAME = %s",Options.listname); break;
          case ( 14): sscanf(line,"ALIGNMENT FILE = %s",Options.ali); break;
          case ( 15): sscanf(line,"MODEL NAME = %s",Options.modelname); break;
          case ( 16): sscanf(line,"PATH TO RESULTS = %s",Options.path); break;
          case ( 17): sscanf(line,"HOMOOLIGOMER = %d",&Options.homo); break;
	  }
                     
      }
    }
               
  } 
  return(0);
}
