 void sim_polIITwoWindows(const double init,const double elong,const double pause1,const double restart1,const double pause2,const double restart2,const double TSS_mu,const double TSS_sigma,const double window_size,const double termination,vector<int>& NETseqDist,vector<int>& PROseqDist,const int G,const int W,const int num_cells,double dt,const int num_timeSteps,int& num_init,int& num_term,int& num_finished){

   vector<vector<int> > polII_location;
   vector<vector<int> > polII_status;
   vector<int> promoter (num_cells,1);

   for (int cell=0;cell<num_cells;++cell){
     polII_location.push_back(vector<int> ());
     polII_status.push_back(vector<int> ());
   }

   //First time-step
   for (int cell=0;cell<num_cells;++cell){
     if (gsl_rng_uniform(r)<init*dt){
       int TSS_start = int(round(gsl_ran_gaussian(r,TSS_sigma)+TSS_mu+gsl_ran_poisson(r,elong*dt)));
       polII_location[cell].push_back(TSS_start+250);
       polII_status[cell].push_back(0);
     }
   }

   for (int cell=0;cell<num_cells;++cell){
     for (int timeStep=1;timeStep<num_timeSteps;++timeStep){
       //Move polymerases on the Gene
       termination_reactions(polII_location[cell],polII_status[cell],termination,num_term);
       move_polIITwoWindows(polII_location[cell],polII_status[cell],elong*dt,pause1*dt,restart1*dt,pause2*dt,restart2*dt,window_size);
       check_overlap(polII_location[cell],polII_status[cell],G,W,num_term,num_finished);
       check_bind(polII_location[cell],polII_status[cell],init*dt,elong*dt,TSS_mu,TSS_sigma,W,num_init,num_term,G);
     }

     for (int it=0;it<(int)polII_location[cell].size();++it){
       if (polII_location[cell][it]>-1){
         NETseqDist[polII_location[cell][it]]++;
         if (polII_status[cell][it]==0){
           PROseqDist[polII_location[cell][it]]++;
         }
       }
     }
   }



 }
