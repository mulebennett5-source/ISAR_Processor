c
c                       File kalman.h
c
c***********************************************************************
c***********************************************************************
c
c     kalmdim                      integer Maximum Kalman Matrix dim
c 
c     dokalm                       integer Flag to toggle k-filtering
c
c     state_phi(kalmdim,kalmdim)   Real Kalman State Transition Matrix
c
c     proc_nois_q(kalmdim,kalmdim) Real Process Noise Matrix
c
c     kalm_gain_k(kalmdim)         Real Kalman Gain Vector
c
c     err_cov_p(kalmdim,kalmdim)   Real Error Covariance Matrix
c
c     state_vec_x(kalmdim)         Real State Vector  
c
c-----------------------------------------------------------------------
c     
      integer     kalmdim , qdim , iter_count , major_iter_count ,
     .            num_per_frm , dokalm
c
      parameter ( kalmdim = 6 )
c
      real        state_phi(kalmdim,kalmdim) ,
     .            proc_nois_q(kalmdim,kalmdim) ,
     .            pnois1 , pnoisr , pnoisf ,
     .            kalm_gain_k(kalmdim) ,
     .            err_cov_p(kalmdim,kalmdim) ,
     .            state_vec_x(kalmdim) ,
     .            state_phi_inv(kalmdim,kalmdim) ,
     .            proc_nois_gamma(kalmdim,kalmdim)
c
      common / kalmat / state_phi , proc_nois_q , proc_nois_gamma ,
     .                  kalm_gain_k , err_cov_p , state_phi_inv ,
     .                  state_vec_x , iter_count , major_iter_count ,
     .                  qdim , num_per_frm , dokalm , pnois1 , pnoisr ,
     .                  pnoisf
c 
c***********************************************************************
c***********************************************************************
c
