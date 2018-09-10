void move_polII(vector<int>& polII_location,vector<int>& polII_status,const double elong,const double pause,const double restart){
  if (polII_location.size()>0){
    for (int it=0;it<(int)polII_location.size();++it){
      if (polII_status[it]==0){//PolII is moving
        if (gsl_rng_uniform(r)<pause){//Pause polII
          polII_status[it]=1;
        }else{//Move
          polII_location[it]+=gsl_ran_poisson(r,elong);
        }
      }else{//PolII is paused
        if (gsl_rng_uniform(r)<restart){//Restart polII
          polII_status[it]=0;
          polII_location[it]+=gsl_ran_poisson(r,elong);
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

void check_bind(vector<int>& polII_location,vector<int>& polII_status,const double init,const double elong,const int TSS_mu,const double TSS_sigma,const int W,int& num_init,int& num_term,const int G){
  if (polII_location.size()>0){//Gene has polII on
    if (polII_location[0]>TSS_mu+W){//Promoter is empty
      if (gsl_rng_uniform(r)<init){//Initiate polII transcription
        ++num_init;
        int TSS_start = round(gsl_ran_gaussian(r,TSS_sigma)+TSS_mu);
        polII_location.insert(polII_location.begin(),TSS_start+gsl_ran_poisson(r,elong));
        polII_status.insert(polII_status.begin(),0);
        if (polII_location[0]>polII_location[1]-W){
          polII_location[0] = polII_location[1]-W;
          polII_status[0]=1;
          if (polII_status[1]==1){//remove PolII
            ++num_term;
            polII_location.erase(polII_location.begin()+1);
            polII_status.erase(polII_status.begin()+1);
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
      int TSS_start = round(gsl_ran_gaussian(r,TSS_sigma)+TSS_mu);
      if (TSS_start<G){
        polII_location.push_back(gsl_ran_poisson(r,elong));
        polII_status.push_back(0);
      }
    }
  }
}

