void move_polIITwoWindows(vector<int>& polII_location,vector<int>& polII_status,const double elong,const double pause1,const double restart1,const double pause2,const double restart2,const double window_size){
  if (polII_location.size()>0){
    for (int it=0;it<(int)polII_location.size();++it){
      if (polII_status[it]==0){//PolII is moving
        if (polII_location[it]<(int)window_size){
          if (gsl_rng_uniform(r)<pause1){
            polII_status[it]=1;
          }else{//Move
            polII_location[it]+=gsl_ran_poisson(r,elong);
          }
        }else{
          if (gsl_rng_uniform(r)<pause2){//Pause polII
            polII_status[it]=1;
          }else{//Move
            polII_location[it]+=gsl_ran_poisson(r,elong);
          }
        }
      }else{//PolII is paused
        if (polII_location[it]<(int)window_size){
          if (gsl_rng_uniform(r)<restart1){
            polII_status[it]=0;
            polII_location[it]+=gsl_ran_poisson(r,elong);
          }
        }else{
          if (gsl_rng_uniform(r)<restart2){
            polII_status[it]=0;
            polII_location[it]+=gsl_ran_poisson(r,elong);
          }
        }
      }
    }
  }
}

void check_overlap(vector<int>& polII_location,vector<int>& polII_status,const int G,const int W,int& num_term,int& num_finished){
  if (polII_location.size()>0){
    for (int it=(int)polII_location.size()-1;it>-1;--it){
      if (polII_location[it]>=G){//Check if PolII is finished transcribing
        ++num_finished;
        polII_location.erase(polII_location.begin()+it);
        polII_status.erase(polII_status.begin()+it);
        continue;
      }

      if (it>0){
        if (polII_location[it-1]+W>polII_location[it]){//Overlapping polII
          polII_location[it-1] = polII_location[it]-W;
          polII_status[it-1]=1;
          if (polII_status[it]==1){//remove PolII
            ++num_term;
            polII_location.erase(polII_location.begin()+it);
            polII_status.erase(polII_status.begin()+it);
          }
        }
      }

    }
  }
}

void check_bindTwoWindows(vector<int>& polII_location,vector<int>& polII_status,int& promoter,const double act,const double inact,const double init,const double elong,const double TSS_mu,const double TSS_sigma,const int W,int& num_init,int& num_term,const int G){
  if (promoter==1){
    if (gsl_rng_uniform(r)<inact){//Inactivate promoter
      promoter=0;
    }else{
      if (polII_location.size()>0){//Gene has polII on
        if (polII_location[0]>TSS_mu+W){//Promoter is empty
          if (gsl_rng_uniform(r)<init){//Initiate polII transcription
            ++num_init;
            int TSS_start = int(round(gsl_ran_gaussian(r,TSS_sigma)+TSS_mu));
            polII_location.insert(polII_location.begin(),TSS_start+gsl_ran_poisson(r,elong));
            polII_status.insert(polII_status.begin(),0);
            int it=1;
            if (polII_location[it-1]+W>polII_location[it]){//Overlapping polII
              polII_location[it-1] = polII_location[it]-W;
              polII_status[it-1]=1;
              if (polII_status[it]==1){//remove PolII
                ++num_term;
                polII_location.erase(polII_location.begin()+it);
                polII_status.erase(polII_status.begin()+it);
              }
            }
            if (polII_location[0]>=G){
              polII_location.erase(polII_location.begin());
              polII_status.erase(polII_status.begin());
            } 
          }
        }
      }else{//Gene is empty
        if (gsl_rng_uniform(r)<init){//initiate polII transcription
          ++num_init;
          int TSS_start = int(round(gsl_ran_gaussian(r,TSS_sigma)+TSS_mu));
	  if (TSS_start < G){
            polII_location.push_back(TSS_start+gsl_ran_poisson(r,elong));
            polII_status.push_back(0);
          }
        }
      }
    }
  }else{
    if (gsl_rng_uniform(r)<act){//Activate promoter
      promoter=1;
    }
  }

}
