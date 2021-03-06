 void sim_polII(const double init,const double elong,const double pause,const double restart,const double backtrack,const double backtrack_restart, const double termination, vector<int>& movingStalledBacktrackedPolIIDist, vector<int>& movingStalledPolIIDist, vector<int>& movingPolIIDist,const int G,const int W,const int num_cells,double dt,const int num_timeSteps,int& num_init,int& num_term,int& num_finished){

   vector<vector<int> > polII_location;
   vector<vector<int> > polII_status;

   for (int cell=0;cell<num_cells;++cell){
     polII_location.push_back(vector<int> ());
     polII_status.push_back(vector<int> ());
   }

   //First time-step
   for (int cell=0;cell<num_cells;++cell){
     if (gsl_rng_uniform(r)<init*dt){
       polII_location[cell].push_back(gsl_ran_poisson(r,elong*dt));
       polII_status[cell].push_back(0);
     }
   }

   for (int cell=0;cell<num_cells;++cell){
     for (int timeStep=1;timeStep<num_timeSteps;++timeStep){
       //Move polymerases on the Gene
       termination_reactions(polII_location[cell],polII_status[cell],termination,num_term);
       move_polII(polII_location[cell],polII_status[cell],elong*dt,pause*dt,restart*dt,backtrack*dt,backtrack_restart*dt);
       check_overlap(polII_location[cell],polII_status[cell],G,W,num_term,num_finished);
       check_bind(polII_location[cell],polII_status[cell],init*dt,elong*dt,W,num_init,num_term);
     }

     for (int it=0;it<(int)polII_location[cell].size();++it){
       if (polII_location[cell][it]>-251){
         if (polII_status[cell][it]==0){//PolII is moving
           movingPolIIDist[polII_location[cell][it]+250]++;
           movingStalledPolIIDist[polII_location[cell][it]+250]++;
           movingStalledBacktrackedPolIIDist[polII_location[cell][it]+250]++;
         }else if(polII_status[cell][it]==1){//PolII is stalled
           movingStalledPolIIDist[polII_location[cell][it]+250]++;
           movingStalledBacktrackedPolIIDist[polII_location[cell][it]+250]++;
         }else if(polII_status[cell][it]==2){//PolII is backtracked
           movingStalledBacktrackedPolIIDist[polII_location[cell][it]+250]++;
         }
       }
     }
   }

 }
