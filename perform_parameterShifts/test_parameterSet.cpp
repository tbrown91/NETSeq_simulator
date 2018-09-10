#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <list>
#include <string>
#include <cmath>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;

gsl_rng * r;

#include "param_funcs.hpp"
#include "early_termination.hpp"
#include "polII_reactionsTwoWindows.hpp"
#include "elongation_simulationTwoWindows.hpp"

int main(int argc,char* argv[]){//Files to give as arguments are a text file with the list of parameters to be tested and the base of the output files

  //Make params vectors from file
  vector<vector<double> > sample_params;

  ifstream infile(argv[1]);
  string line;
  char buf[256];
  while(infile.good()){
    getline(infile,line);
    istringstream iss(line);
    string token;
    sample_params.push_back(vector<double>());
    while(getline(iss,token,'\t')){
        sample_params.back().push_back(stof(token));
    }
  }

  //Make files for output
  string out_fileNameMoving;
  string out_fileNameMovingStalled;
  string out_fileNameMovingStalledBacktracked;
  string out_fileNameCounts;
  string buffer;
  buffer = argv[2];
  out_fileNameMovingStalledBacktracked = buffer+"_MovingStalledBacktracked.out";
  out_fileNameMovingStalled = buffer+"_MovingStalled.out";
  out_fileNameMoving = buffer+"_Moving.out";
  out_fileNameCounts = buffer+"_EventCounts.out";

  ofstream out_fileMovingStalledBacktracked(out_fileNameMovingStalledBacktracked);
  ofstream out_fileMovingStalled(out_fileNameMovingStalled);
  ofstream out_fileMoving(out_fileNameMoving);
  ofstream out_fileCounts(out_fileNameCounts);

  out_fileCounts<<"NUM_INIT\tNUM_TERM\tNUM_FINISHED"<<endl;

  const int G = 1000;//Gene length
  const int W = 40;//PolII footprint size
  const int num_cells = 10000;//Number of cells per simulation
  const double dt = 0.005;//Time increment in reactions
  const int T = 40;//Reaction length

  int num_timeSteps = (int)T/dt;

  //Setup random environment
  const gsl_rng_type * I;
  gsl_rng_env_setup();
  I = gsl_rng_default;
  r = gsl_rng_alloc (I);
  // gsl_rng_set(r,time(0)+my_rank);
  gsl_rng_set(r,1);//For debugging only

  for (int it=0;it<(int)sample_params.size();++it){

    vector<int> movingStalledBacktrackedPolIIDist (G+250,0);
    vector<int> movingStalledPolIIDist (G+250,0);
    vector<int> movingPolIIDist (G+250,0);
    int num_init = 0;
    int num_term = 0;
    int num_finished = 0;

    sim_polIITwoWindows(sample_params[it][0],sample_params[it][1],sample_params[it][2],sample_params[it][3],sample_params[it][4],sample_params[it][5],sample_params[it][6],sample_params[it][7],sample_params[it][8],sample_params[it][9],sample_params[it][10],sample_params[it][11],movingStalledBacktrackedPolIIDist,movingStalledPolIIDist,movingPolIIDist,G,W,num_cells,dt,num_timeSteps,num_init,num_term,num_finished);

    for (int nt=0;nt<G+250;++nt){
      out_fileMovingStalledBacktracked<<movingStalledBacktrackedPolIIDist[nt]<<" ";
      out_fileMovingStalled<<movingStalledPolIIDist[nt]<<" ";
      out_fileMoving<<movingPolIIDist[nt]<<" ";
    }
    out_fileMovingStalledBacktracked<<endl;
    out_fileMovingStalled<<endl;
    out_fileMoving<<endl;

    out_fileCounts<<num_init<<"\t"<<num_term<<"\t"<<num_finished<<endl;

  }

  out_fileMovingStalledBacktracked.close();
  out_fileMovingStalled.close();
  out_fileMoving.close();
  out_fileCounts.close();

  return 0;
}
