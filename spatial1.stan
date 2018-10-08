data{
      int<lower = 1> N; //number of points per subject
      int<lower = 1> Nx;
      int<lower = 1> Ny;
      int<lower = 1> Nsubj;//number of subject
      vector[N] y; // reference map
      vector[N] X[Nsubj]; // activation map
      vector[N] coordx[Nsubj];
      vector[N] coordy[Nsubj];
      
      matrix[N,N] dist[Nsubj]; //distance matrix per activation map;
      matrix[N,N] dist_y; //distance matrix of reference map
      
      real s_tx;
      real s_ty;
      real s_s;
      real apsi;
      real bpsi;
      real aphi;
      real bphi;
}

transformed data{
        vector[N] mu0 = rep_vector(0, N); 
}

parameters{
      vector[Nsubj] theta_x;
      vector[Nsubj] theta_y;
      real<lower = 0> sigma[Nsubj];
      real<lower = 0> phi;
      real<lower = 0> psi;

}

transformed parameters{
      matrix[N,N] dist2[Nsubj];
      matrix[N,N] Sig;
      matrix[N,N] Sig1[Nsubj];
      cholesky_factor_cov[N,N] L;
      cholesky_factor_cov[N,N] L1[Nsubj]; 


      vector[N] coordx_t[Nsubj];
      vector[N] coordy_t[Nsubj];

      
      vector[N] mu_kriged[Nsubj];
      
      Sig = exp((-1)*psi*dist_y);
      L = cholesky_decompose(Sig);
      
        for (k in 1:Nsubj) {
          
          coordx_t[k] = (coordx[k] + theta_x[k]) / sigma[k];
          coordy_t[k] = (coordy[k] + theta_y[k]) / sigma[k];
          for (l in 1:N) for (m in 1:N)  dist2[k][l,m] = sqrt((coordx_t[k][l] - coordx[k][m])^2 + (coordy_t[k][l] - coordy[k][m])^2)/sigma[k];
          
          mu_kriged[k] = dist2[k] * inverse(dist[k]) * y;

          Sig1[k] = exp((-1) * phi * dist[k]/sigma[k]);

          L1[k] = cholesky_decompose(Sig1[k]);
        }
}
model{
      
      // priors
      theta_x ~ normal(0,s_tx);
      theta_y ~ normal(0,s_ty);
      sigma ~ lognormal(0,s_s);
      psi ~ inv_gamma(apsi,bpsi);
      phi ~ inv_gamma(aphi,bphi);
      y ~ multi_normal_cholesky(mu0, L);

      for (i in 1:Nsubj) {
        //likelihood
        X[i] ~ multi_normal_cholesky(mu_kriged[i], L1[i]); 
      }
  
}
