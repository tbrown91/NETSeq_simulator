void move_polII(vector<int>& polII_location,vector<int>& polII_status,const double elong){
  if (polII_location.size()>0){
    for (int it=0;it<(int)polII_location.size();++it){
      polII_location[it]+=gsl_ran_poisson(r,elong);
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
        }
      }
    }
  }
}

void check_bind(vector<int>& polII_location,vector<int>& polII_status,const double init,const double elong, const int W,int& num_init,int& num_term){
  if (polII_location.size()>0){//Gene has polII on
    if (polII_location[0]>W){//Promoter is empty
      if (gsl_rng_uniform(r)<init){//Initiate polII transcription
        ++num_init;
        polII_location.insert(polII_location.begin(),gsl_ran_poisson(r,elong));//Move initial amount once bound
        polII_status.insert(polII_status.begin(),0);
        if (polII_location[0]>polII_location[1]-W){//Check for overlap
          polII_location[0] = polII_location[1]-W;
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
