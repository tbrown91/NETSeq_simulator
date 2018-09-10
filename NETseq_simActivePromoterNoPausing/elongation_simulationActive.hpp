 void sim_polII(const double init,const double elong, vector<int>& NETseqDist,vector<int>& PROseqDist,const int G,const int W,const int num_cells,double dt,const int num_timeSteps,int& num_init,int& num_term,int& num_finished){

   vector<vector<int> > polII_location;
   vector<vector<int> > polII_status;

   for (int cell=0;cell<num_cells;++cell){
     polII_location.push_back(vector<int> ());
     polII_status.push_back(vector<int> ());
   }

   //First time-step
   for (int cell=0;cell<num_cells;++cell){
     if (gsl_rng_uniform(r)<init*dt){
       polII_location[cell].push_back(gsl_ran_poisson(r,elong*dt));//Move initial amount once bound
       polII_status[cell].push_back(0);
     }
   }

   for (int cell=0;cell<num_cells;++cell){
     for (int timeStep=1;timeStep<num_timeSteps;++timeStep){
       //Move polymerases on the Gene
       move_polII(polII_location[cell],polII_status[cell],elong*dt);
       check_overlap(polII_location[cell],polII_status[cell],G,W,num_term,num_finished);
       check_bind(polII_location[cell],polII_status[cell],init*dt,elong*dt,W,num_init,num_term);
     }

     for (int it=0;it<(int)polII_location[cell].size();++it){
       if (polII_location[cell][it]>-251){
         NETseqDist[polII_location[cell][it]+250]++;
         if (polII_status[cell][it]==0){//PolII is moving
           PROseqDist[polII_location[cell][it]+250]++;
         }
       }
     }
   }



 }
