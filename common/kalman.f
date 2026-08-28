c**
c********************************************************************
c**
      subroutine info_filter ( accel_y , observe_h , adim , qrtol ) 

c**
c********************************************************************
c**
c     Added by KAM 7-6-98
c
c     This information filter variant of the Kalman Filter was
c     implemented so that a block update could be implemented without
c     having to take n_inside X n_inside inverses.
c     
      implicit none
c     
      include 'kalman.h'
c
      include 'tglist.h'
c
      real     accel_y(nlist+8) , observe_h(nlist+8,kalmdim) , qrtol
c 
      integer  adim , jpvt(kalmdim) , kr , i , j
c
      real     qraux(kalmdim)              ,
     .         qrwork(kalmdim)             ,
     .         qrrsd(kalmdim)              ,
     .         info_mat_a(kalmdim,kalmdim) ,
     .         temp_axa(kalmdim,kalmdim)   ,
     .         temp_qxq(kalmdim,kalmdim)   ,
     .         temp_qxa(kalmdim,kalmdim)   , 
     .         temp_a(kalmdim)             ,
     .         temp2_a(kalmdim)            ,
     .         hrh(kalmdim,kalmdim)        ,
     .         b(kalmdim)                  ,
     .         a_gamm(kalmdim,kalmdim)     ,
     .         gamm_a(kalmdim,kalmdim)     ,
     .         pred_state_x(kalmdim)
c
c---------------------------------------------------------------------
c
c     Update the Noiseless State Prediction Information Matrix
c
c     A(k-1) = info_mat_a
c     F(k-1)^-1 = state_phi_inv
c     P(k|k)^-1 = err_cov_p
c     G(k-1)    = proc_nois_gamma
c     H(k)'R(k)^-1/2  = u'
c     W(k) = Kalman gain
c     x(k) = state_vec_x
c     y(k) = accel_y
c     Algorithm:
c
c         1. A(k-1) = [F(k-1)^-1]'P(k-1|k-1)^-1F(k-1)^-1
c         
c         2. P(k|k-1)=A(k-1)-A(k-1)G(k-1)[G(k-1)'A(k-1)G(k-1)+
c            Q(k-1)^-1]^-1G(k-1)'A(k-1)
c
c         3. W(k) = [P(k|k-1)^-1 + H(k)'R(k)H(k)]^-1H(k)'R(k)^-1
c
c         4. P(k|k)^-1=P(k|k-1)^-1 + H(k)'R(k)^-1 H(k)
c         
c         5. x(k) = F(k-1)x(k-1) + W(k)[y(k) - H(k)x(k-1)]
c
c     Note R(k) ~ Rc or meas_err_rc is not explicitly used.  In forming
c     U and b in getacc() we divided by the Doppler width = sqrt(rc(i)).
c     This is taken into account here.
c
      pred_state_x = matmul( state_phi , state_vec_x )
c
      temp_axa     = matmul( transpose( state_phi_inv ) , err_cov_p )
c
      info_mat_a   = matmul( temp_axa , state_phi_inv )
c
c     Update the State Prediction Information Matrix
c
      a_gamm       = matmul( info_mat_a , proc_nois_gamma )
c
      temp_qxq     = matmul( transpose( proc_nois_gamma) , a_gamm )
c
      do j = 1 , qdim
c
c        Assuming the the Process Noise is upcoupled
c
         temp_qxq(j,j) = temp_qxq(j,j) + 1.0 / proc_nois_q(j,j)
c     
      enddo
c
      gamm_a = matmul( transpose( proc_nois_gamma ) , info_mat_a )
c
      call sqrank ( temp_qxq , kalmdim , qdim , qdim , qrtol ,
     .              kr , jpvt , qraux , qrwork )
c
      do j = 1 , adim
c
         do i = 1 , adim
c
            b(i) = gamm_a(i,j)
c
         enddo
c
         call sqrlss ( temp_qxq , kalmdim , qdim , qdim , kr ,
     .                 b , temp_qxa(1,j) , qrrsd , jpvt , qraux )       
c
      enddo
c     
      err_cov_p = matmul( a_gamm , temp_qxa )
c
      do j = 1 , adim
c
         do i = 1 , adim
c
            err_cov_p(i,j) = info_mat_a(i,j) - err_cov_p(i,j)
c
         enddo
c
      enddo
c
      hrh = matmul( transpose( observe_h ) , observe_h )
c
c     Calculate the Kalman Gain
c
      temp_a  = matmul( transpose( observe_h ) , accel_y )
c      
      temp2_a = matmul( hrh , state_vec_x )
c
      do j = 1 , adim
c
         do i = 1 , adim
c
            err_cov_p(i,j) =  err_cov_p(i,j) + hrh(i,j)
c
         enddo
c
      enddo

      do j = 1 , adim
c
         temp_a(j) = temp_a(j) - temp2_a(j)
c
      enddo
c
      call sqrank ( err_cov_p , kalmdim , adim , adim , qrtol ,
     .              kr , jpvt , qraux , qrwork )
c      
      call sqrlss ( err_cov_p , kalmdim , adim , adim , kr ,
     .              temp_a, temp2_a , qrrsd , jpvt , qraux )       
c      
c     Update States
c
      do j = 1 , adim
c
         state_vec_x(j) = pred_state_x(j) + temp2_a(j)
c
c        write ( 6 , * ) state_vec_x(j) , pred_state_x(j) , temp2_a(j)
c
      enddo
c
c     read ( 5 , * )
c
c     Added by KAM 7-8-98
c     Constrain Ar to be negative
c
c     if ( state_vec_x(3) .gt. 0.0 ) then
c
c        state_vec_x(3) = 0.0
c        state_vec_x(4) = 0.0
c
c     endif
c     
      return
      end
c**
c***********************************************************************
c**
      subroutine kalman_initialize ( del_t , adim )
c**
c***********************************************************************
c**
c
c     Added by KAM 7-2-98
c
      implicit none
      
c
      include 'kalman.h'
c
      real     del_t
c      
      integer  adim
c
c     Locals
c
      integer i , j
c
c-----------------------------------------------------------------------
c
      do j = 1 , kalmdim
c
         do i = 1 , kalmdim
c
            err_cov_p(i,j)       = 0.0
c
            state_phi(i,j)       = 0.0
c
            state_phi_inv(i,j)   = 0.0
c
            proc_nois_q(i,j)     = 0.0
c
            proc_nois_gamma(i,j) = 0.0      
c
         enddo
c
         kalm_gain_k(j) = 0.0
c
         state_vec_x(j) = 0.0
c
      enddo
c
      do j = 1 , adim
c
         err_cov_p(j,j) = 0.0
c
      enddo
c
      call kalman_maintain ( del_t , adim ) 
c
      iter_count       = 0
      major_iter_count = 0
      num_per_frm      = 1    
c      
      return
      end
c**
c***********************************************************************
c**
      subroutine kalman_maintain ( del_t , adim )
c**
c***********************************************************************
c**
c
c     Added by KAM 7-2-98
c
      implicit none
c
      real     del_t
c
      integer  adim
c
      include 'kalman.h'
c
c     Locals
c      
      integer  j
c
c-----------------------------------------------------------------------
c
c     Reorder due to observation vector ordering [1 t r f rt ft]
c     Replace the above expression for one which redoes the coupling
c     pair by pair.  Only time couplings used so far
c    
      qdim = 1 
c
      do j = 1 , adim
c      
         state_phi(j,j)     = 1.0
         state_phi_inv(j,j) = 1.0        
c
      enddo
c
      proc_nois_gamma(1,1) = del_t
c
      proc_nois_q(1,1)     = pnois1
c
      if ( adim .gt. 1 ) then
c
         state_phi(1,2)       = 0           !    del_t
         state_phi_inv(1,2)   = 0           !  - del_t
c
c        proc_nois_gamma(1,1) = del_t * del_t * 0.5
c        proc_nois_gamma(2,1) = del_t
c              
      endif
c
      if ( adim .gt. 2 ) then
c
         qdim                 = qdim + 1
         proc_nois_gamma(3,2) = del_t
         proc_nois_q(2,2)     = pnoisr
c
      endif
c
      if ( adim .gt. 3 ) then
c
         qdim                 = qdim + 1
         proc_nois_gamma(4,3) = del_t
         proc_nois_q(3,3)     = pnoisf
c
      endif
c
      if ( adim .gt. 4 ) then
c
         state_phi(3,5)       = 0           !    del_t
         state_phi_inv(3,5)   = 0           !  - del_t
c
         proc_nois_gamma(3,2) = del_t * del_t * 0.5
         proc_nois_gamma(5,2) = del_t             
c
      endif
c
      if ( adim .gt. 5 ) then
c
         state_phi(4,6)       = 0           !    del_t
         state_phi_inv(4,6)   = 0           !  - del_t
c
         proc_nois_gamma(4,3) = del_t * del_t * 0.5
         proc_nois_gamma(6,3) = del_t
c
      endif      
c
      return
      end
