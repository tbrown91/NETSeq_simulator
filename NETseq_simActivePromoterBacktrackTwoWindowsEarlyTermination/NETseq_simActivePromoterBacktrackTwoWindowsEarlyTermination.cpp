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
#include "early_termination.hpp"
#include "polII_reactionsTwoWindows.hpp"
#include "elongation_simulationTwoWindows.hpp"

int main(int argc,char* argv[]){

  MPI_Init(NULL,NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  string out_fileName;
  string out_fileNameMoving;
  string out_fileNameMovingStalled;
  string out_fileNameMovingStalledBacktracked;
  string buffer;
  buffer = argv[1];
  out_fileName = buffer+"_"+to_string(my_rank)+".out";
  out_fileNameMovingStalledBacktracked = buffer+"_"+to_string(my_rank)+"MovingStalledBacktracked.out";
  out_fileNameMovingStalled = buffer+"_"+to_string(my_rank)+"MovingStalled.out";
  out_fileNameMoving = buffer+"_"+to_string(my_rank)+"Moving.out";

  ofstream out_file(out_fileName);
  ofstream out_fileMovingStalledBacktracked(out_fileNameMovingStalledBacktracked);
  ofstream out_fileMovingStalled(out_fileNameMovingStalled);
  ofstream out_fileMoving(out_fileNameMoving);

  out_file<<"INIT\tELONG\tPAUSE1\tRESTART1\tBACKTRACK1\tBACKTRACK_RESTART1\tPAUSE2\tRESTART2\tBACKTRACK2\tBACKTRACK_RESTART2\tWINDOW_SIZE\tTERMINATION_SIZE\tNUM_INIT\tNUM_TERM\tNUM_FINISHED"<<endl; 

  const int num_params = 10000;
  const int G = 1000;//Gene length
  const int W = 40;//PolII footprint size
  const int num_cells = 100000;//Number of cells per simulation
  const double dt = 0.005;//Time increment in reactions
  const int T = 40;//Reaction length

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
  for (int i=0;i<12;++i){
    sample_params.push_back(vector<double> (num_params,0.0));
  }

  choose_param(sample_params[0],0.01,60.0,num_params);//Initiation rate (/min)
  choose_param(sample_params[1],100,5000,num_params);//Elongation rate (/min)
  choose_param(sample_params[2],0.01,60.0,num_params);//First pause rate (/min)
  choose_param(sample_params[3],0.01,60.0,num_params);//First restart rate (/min)
  choose_param(sample_params[4],0.01,60.0,num_params);//First Backtrack rate (/min)
  choose_param(sample_params[5],0.01,60.0,num_params);//First backtrack restart rate (/min)
  choose_param(sample_params[6],0.01,60.0,num_params);//Second pause rate (/min)
  choose_param(sample_params[7],0.01,60.0,num_params);//Second restart rate (/min)
  choose_param(sample_params[8],0.01,60.0,num_params);//Second Backtrack rate (/min)
  choose_param(sample_params[9],0.01,60.0,num_params);//Second backtrack restart rate (/min)
  choose_intParam(sample_params[10],0,G,num_params);//Window size (bp)
  choose_intParam(sample_params[11],0,G,num_params);//Termination window (bp)

  for (int it=0;it<num_params;++it){

    vector<int> movingStalledBacktrackedPolIIDist (G+250,0);
    vector<int> movingStalledPolIIDist (G+250,0);
    vector<int> movingPolIIDist (G+250,0);
    int num_init = 0;
    int num_term = 0;
    int num_finished = 0;

    sim_polIITwoWindows(sample_params[0][it],sample_params[1][it],sample_params[2][it],sample_params[3][it],sample_params[4][it],sample_params[5][it],sample_params[6][it],sample_params[7][it],sample_params[8][it],sample_params[9][it],sample_params[10][it],sample_params[11][it],movingStalledBacktrackedPolIIDist,movingStalledPolIIDist,movingPolIIDist,G,W,num_cells,dt,num_timeSteps,num_init,num_term,num_finished);

    for (int p=0;p<12;++p){
      out_file<<sample_params[p][it]<<"\t";
    }
    out_file<<num_init<<"\t"<<num_term<<"\t"<<num_finished;
    out_file<<endl;
    for (int nt=0;nt<G+250;++nt){
      out_fileMovingStalledBacktracked<<movingStalledBacktrackedPolIIDist[nt]<<" ";
      out_fileMovingStalled<<movingStalledPolIIDist[nt]<<" ";
      out_fileMoving<<movingPolIIDist[nt]<<" ";
    }
    out_fileMovingStalledBacktracked<<endl;
    out_fileMovingStalled<<endl;
    out_fileMoving<<endl;

  }

  out_file.close();
  out_fileMovingStalledBacktracked.close();
  out_fileMovingStalled.close();
  out_fileMoving.close();

  MPI_Finalize();

  return 0;
}
