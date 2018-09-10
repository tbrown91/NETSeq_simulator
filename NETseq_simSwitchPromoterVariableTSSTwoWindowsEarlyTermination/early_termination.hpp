void termination_reactions(vector<int>& polII_location,vector<int>& polII_status,const double term,int& num_term){
  if (polII_location.size()>0){
    for (int k=(int)polII_location.size()-1;k>-1;--k){
      if (polII_location[k]>-1){
        if (gsl_rng_uniform(r)<gsl_ran_poisson_pdf(polII_location[k],term)){
          ++num_term;
          polII_location.erase(polII_location.begin()+k);
          polII_status.erase(polII_status.begin()+k);
        }
      }
    }
  }


}
