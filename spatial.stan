data{
      int<lower = 1> N; //number of points per subject
      int<lower = 1> Nx;
      int<lower = 1> Ny;
      int<lower = 1> Nsubj;//number of subject
      //vector[N*Nsubj] y; // raw data without procustus
      matrix[N,Nsubj] y;
      real<lower = 0> coord[N, 2, Nsubj];
      real<lower = 0> dist[N, N, Nsubj];
}

parameters{
      vector[Nsubj] theta_x;
      vector[Nsubj] theta_y;
      real<lower = 0> sigma[Nsubj];
      real<lower = 0> phi;
      real<lower = 0> psi;
}

model{
  
      matrix[N,N] dist1;
      matrix[N,N] dist2;
      matrix[N,N] Sig;
      matrix[N,N] Sig1;
      matrix[N,N] L;
      matrix[N,N] L1; //
      matrix[N,N] cov;
      
      vector[N] yi;
      vector[N] mu0;
      vector[N] coordx_t;
      vector[N] coordy_t;
      vector[N] w;
      vector[N] mu_kriged;
      for (m in 1: N) mu0[m] = 0;
      
      // priors
      theta_x ~ normal(0,20);
      theta_y ~ normal(0,20);
      sigma ~ lognormal(0,20);
      psi ~ inv_gamma(1,1);
      phi ~ inv_gamma(1,1);
      
      for (i in 1:Nsubj) {
        
        for(nx in 1: N) for (ny in 1: N) dist1[nx, ny] = dist[nx, ny, i];  
        for(ny in 1:N) yi[ny] = y[ny, i];
        
        for (k in 1:N) {
          coordx_t[k] = (coord[k,1,i] + theta_x[i]) / sigma[i];
          coordy_t[k] = (coord[k,2,i] + theta_y[i]) / sigma[i];
        }
        
        for (l in 1:N) for (m in 1:N)  dist2[l,m] = sqrt((coordx_t[l]-coord[m,1,i])^2 + (coordy_t[l] - coord[m,2,i])^2);
        
        
        Sig = exp((-1)*psi*dist1);
        L = cholesky_decompose(Sig);
        yi ~ multi_normal_cholesky(mu0,L);
        
        //likelihood
        for (k in 1:N) w[k] = y[k,i];
        mu_kriged = dist2 * dist1 * w;
        Sig1 = exp((-1)*phi/sigma[i]*dist1);
        L1 = cholesky_decompose(Sig1);
        mu_kriged ~ multi_normal_cholesky(mu_kriged,L); 
        
      }
  
}
