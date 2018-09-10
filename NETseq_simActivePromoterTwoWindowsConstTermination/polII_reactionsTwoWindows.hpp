void move_polIITwoWindows(vector<int>& polII_location,vector<int>& polII_status,const double elong,const double pause1,const double restart1,const double pause2,const double restart2,const double window_size,const double term,int& num_term){
  if (polII_location.size()>0){
    for (int it=0;it<(int)polII_location.size();++it){
      if (gsl_rng_uniform(r)<term){//Remove PolII
        polII_location.erase(polII_location.begin()+it);
        polII_status.erase(polII_status.begin()+it);
        ++num_term;
        --it;
      }else{
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

void check_bindTwoWindows(vector<int>& polII_location,vector<int>& polII_status,const double init,const double elong,const int W,int& num_init,int& num_term){
  if (polII_location.size()>0){//Gene has polII on
    if (polII_location[0]>W){//Promoter is empty
      if (gsl_rng_uniform(r)<init){//Initiate polII transcription
        ++num_init;
        polII_location.insert(polII_location.begin(),gsl_ran_poisson(r,elong));
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
      }
    }
  }else{//Gene is empty
    if (gsl_rng_uniform(r)<init){//initiate polII transcription
      ++num_init;
      polII_location.push_back(gsl_ran_poisson(r,elong));
      polII_status.push_back(0);
    }
  }
}
