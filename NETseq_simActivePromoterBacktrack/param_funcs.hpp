void choose_param(vector<double>& sample_params,const double min_param,const double max_param,const int num_params){
  vector<double> possible_params (num_params,0.0);
  for (int i=0;i<num_params;++i){
    possible_params[i] = min_param + (max_param-min_param)*((double)i/(double)(num_params-1));
  }

  for (int i=0;i<num_params;++i){
    int sample_int = gsl_rng_uniform_int(r,num_params-i);
    sample_params[i] = possible_params[sample_int];
    possible_params[sample_int] = possible_params[num_params-1-i];
  }
}
