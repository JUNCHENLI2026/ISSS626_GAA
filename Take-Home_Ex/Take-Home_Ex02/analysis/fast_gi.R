# Exact conditional sampling without replacement for sparse space-time graphs.
# Rejection sampling avoids repeatedly shuffling the entire 18,810-bin population.
# Compilation requires Rcpp and an R-compatible C++ toolchain (Rtools on Windows).
Rcpp::cppFunction(code='
Rcpp::NumericMatrix conditional_gi_sample(Rcpp::NumericVector x,
    Rcpp::List neighbours, Rcpp::List weights, int nsim) {
  int n=x.size();
  Rcpp::NumericMatrix out(n,3);
  for(int i=0;i<n;i++) {
    Rcpp::IntegerVector ids=neighbours[i];
    Rcpp::NumericVector w=weights[i];
    int k=ids.size();
    double observed=0.0, fixed=0.0;
    std::vector<double> other_weights;
    for(int j=0;j<k;j++) {
      observed+=w[j]*x[ids[j]-1];
      if(ids[j]-1==i) fixed+=w[j]*x[i];
      else other_weights.push_back(w[j]);
    }
    int m=other_weights.size();
    std::vector<int> selected(m);
    int lower=0,upper=0; double total=0.0,total2=0.0;
    for(int b=0;b<nsim;b++) {
      double value=fixed;
      for(int j=0;j<m;j++) {
        int draw; bool duplicate;
        do {
          draw=(int)floor(R::runif(0.0,(double)(n-1)));
          if(draw>=i)draw++;
          duplicate=false;
          for(int h=0;h<j;h++)if(selected[h]==draw){duplicate=true;break;}
        } while(duplicate);
        selected[j]=draw;value+=other_weights[j]*x[draw];
      }
      double tolerance=1e-10*std::max(1.0,std::abs(observed));
      if(value<=observed+tolerance)lower++;
      if(value>=observed-tolerance)upper++;
      total+=value;total2+=value*value;
    }
    out(i,0)=std::min(1.0,2.0*(1.0+std::min(lower,upper))/(nsim+1.0));
    out(i,1)=total/nsim;
    out(i,2)=(total2-total*total/nsim)/(nsim-1.0);
    if(i%100==0)Rcpp::checkUserInterrupt();
  }
  return out;
}')
fast_space_time_gi <- function(x,w,nsim=999L,seed=20261018L){
  stopifnot(is.numeric(x),all(is.finite(x)),length(x)==length(w$neighbours),nsim>=99L)
  set.seed(seed)
  sampled<-conditional_gi_sample(x,w$neighbours,w$weights,nsim)
  # Analytical Gi* z-values describe clustering strength; inference uses the
  # conditional randomisation distribution with the focal observation fixed.
  list(z=as.numeric(spdep::localG(x,w,zero.policy=TRUE)),p=sampled[,1],null_mean=sampled[,2],null_variance=sampled[,3])
}
