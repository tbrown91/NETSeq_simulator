void move_polII(vector<int>& polII_location,vector<int>& polII_status,const double elong,const double pause,const double restart){
  if (polII_location.size()>0){
    for (int k=0;k<(int)polII_location.size();++k){
      if (polII_status[k]==0){//PolII is moving
        if (gsl_rng_uniform(r)<pause){//Pause polII
          polII_status[k]=1;
        }else{//Move
          polII_location[k]+=gsl_ran_poisson(r,elong);
        }
      }else{//PolII is paused
        if (gsl_rng_uniform(r)<restart){//Restart polII
          polII_status[k]=0;
          polII_location[k]+=gsl_ran_poisson(r,elong);
        }
      }
    }
  }
}

void check_overlap(vector<int>& polII_location,vector<int>& polII_status,const int G,const int W,int& num_term,int& num_finished){
  if (polII_location.size()>0){
    for (int k=(int)polII_location.size()-1;k>-1;--k){
      if (polII_location[k]>=G){//Check if PolII is finished transcribing
        ++num_finished;
        polII_location.erase(polII_location.begin()+k);
        polII_status.erase(polII_status.begin()+k);
        continue;
      }

      if (k>0){
        if (polII_location[k-1]+W>polII_location[k]){//Overlapping polII
          polII_location[k-1] = polII_location[k]-W;
          polII_status[k-1]=1;
          if (polII_status[k]==1){//remove PolII
            ++num_term;
            polII_location.erase(polII_location.begin()+k);
            polII_status.erase(polII_status.begin()+k);
          }
        }
      }

    }
  }
}

void check_bind(vector<int>& polII_location,vector<int>& polII_status,int& promoter,const double init,const double elong,const int W,int& num_init,int& num_term){
  if (polII_location.size()>0){//Gene has polII on
    if (polII_location[0]>W){//Promoter is empty
      if (gsl_rng_uniform(r)<init){//Initiate polII transcription
        ++num_init;
        polII_location.insert(polII_location.begin(),gsl_ran_poisson(r,elong));

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
