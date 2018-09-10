#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <list>
#include <string>
#include <cmath>
#include <mpi.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;

gsl_rng * r;

#include "param_funcs.hpp"
#include "polII_reactionsTwoWindows.hpp"
#include "elongation_simulationTwoWindows.hpp"

int main(int argc,char* argv[]){

  MPI_Init(NULL,NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  string out_fileName;
  string out_fileNameNETseq;
  string out_fileNamePROseq;
  string buffer;
  buffer = argv[1];
  out_fileName = buffer+"_"+to_string(my_rank)+".out";
  out_fileNameNETseq = buffer+"_"+to_string(my_rank)+"NETseqDist.out";
  out_fileNamePROseq = buffer+"_"+to_string(my_rank)+"PROseqDist.out";
  ofstream out_file(out_fileName);

  ofstream out_fileNETseq(out_fileNameNETseq);
  ofstream out_filePROseq(out_fileNamePROseq);

  out_file<<"ACT\tINACT\tINIT\tELONG\tPAUSE1\tRESTART1\tPAUSE2\tRESTART2\tTSS_mu\tTSS_sigma\tWINDOW_SIZE\tNUM_INIT\tNUM_TERM\tNUM_FINISHED"<<endl;

  const int num_params = 10000;
  const int G = 1000;//Gene length
  const int W = 40;//PolII footprint size
  const int num_cells = 10000;//Number of cells per simulation
  const double dt = 0.005;//Time increment in reactions
  const int T = 80;//Reaction length

  int num_timeSteps = (int)T/dt;

  //Setup random environment
  const gsl_rng_type * I;
  gsl_rng_env_setup();
  I = gsl_rng_default;
  r = gsl_rng_alloc (I);
  // gsl_rng_set(r,time(0)+my_rank);
  gsl_rng_set(r,0 + my_rank);//For debugging only

  //Choose params
  vector<vector<double> > sample_params;
  for (int i=0;i<11;++i){
    sample_params.push_back(vector<double> (num_params,0.0));
  }
  choose_param(sample_params[0],0.01,60.0,num_params);//Promoter activation rate (/min)
  choose_param(sample_params[1],0.01,60.0,num_params);//Promoter inactivation rate (/min)
  choose_param(sample_params[2],0.01,60.0,num_params);//Initiation rate (/min)
  choose_param(sample_params[3],100,5000,num_params);//Elongation rate (/min)
  choose_param(sample_params[4],0.01,60.0,num_params);//First pause rate (/min)
  choose_param(sample_params[5],0.01,60.0,num_params);//First restart rate (/min)
  choose_param(sample_params[6],0.01,60.0,num_params);//Second pause rate (/min)
  choose_param(sample_params[7],0.01,60.0,num_params);//Second restart rate (/min)
  choose_intParam(sample_params[8],-250,250,num_params);//Start site target, relative to annotated TSS
  choose_param(sample_params[9],0.01,2500.0,num_params);//Standard deviation of start site selection
  choose_intParam(sample_params[10],0,G,num_params);//Window size (bp)

  for (int it=0;it<num_params;++it){

    vector<int> NETseqDist (G+250,0);
    vector<int> PROseqDist (G+250,0);
    int num_init = 0;
    int num_term = 0;
    int num_finished = 0;

    sim_polIITwoWindows(sample_params[0][it],sample_params[1][it],sample_params[2][it],sample_params[3][it],sample_params[4][it],sample_params[5][it],sample_params[6][it],sample_params[7][it],sample_params[8][it],sample_params[9][it],sample_params[10][it],NETseqDist,PROseqDist,G,W,num_cells,dt,num_timeSteps,num_init,num_term,num_finished);

    for (int p=0;p<11;++p){
      out_file<<sample_params[p][it]<<"\t";
    }
    out_file<<num_init<<"\t"<<num_term<<"\t"<<num_finished;
    out_file<<endl;
    for (int nt=0;nt<G+250;++nt){
      out_fileNETseq<<NETseqDist[nt]<<" ";
      out_filePROseq<<PROseqDist[nt]<<" ";
    }
    out_fileNETseq<<endl;
    out_filePROseq<<endl;

    //cout<<"Rank: "<<my_rank<<", iteration: "<<it<<" finished"<<endl;
  }

  out_file.close();
  out_fileNETseq.close();
  out_filePROseq.close();

  MPI_Finalize();

  return 0;
}
