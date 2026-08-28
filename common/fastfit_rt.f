c**
c***********************************************************************
c**
      subroutine fastfit_rt ( a , aw , b , bw , n , nc , work1 , work2 ) 
c**
c***********************************************************************
c**
      implicit none
!
!     Purpose:
!
!		A four coefficent least squares fit Ax ~ b
!
!     Input:
!
!		A		matrix for the fit size of  n x 4
!		b		vector for the fit size of n
!		n		number of data points >= 4
!
!     Output:
!
!		b	the coefficents of the fit b(1)-b(4)
!
!     Created:
!
!         Jeff V. Richard 5/26/98
!
!     Notes:
!
!       This code is limited to a 4 coefficent fit.  LAPACK and other
!       codes have general least squares solvers.  This code is only
!       intended for a simple real time application
!
!     Method:
!
!       The algorithm uses a rank one update using Householder Transfor-
!       mations.  In this case the rank one update is performed 3 times.
!       At each iteration both the matrix A and the vector b are updated
!       in place, modifying the current matrix and vector.  For a
!       complete description of this method see "Matrix Computations",
!       Gene Golub, Charles Van Loan,John Hopkins University Press, 1985
!
      integer n , nc
!
      real    a(n,nc) , b(n) , aw(n,nc) , bw(n) , work1(n) , work2(n,n)
!
!	Internal declarations
!
      integer j , k , ii , jj , kk , nl
!
      real    temp , alpha , beta
!
!	Perform rank 1 updates: n must be at least equal to 4
!
      nl = min( n - 1 , nc )
!
      do j = 1 , nl
!
         temp = abs( a(1,j) )
!
         work1(1) = 0.0
!
         do k = 2 , n
!
            If ( temp .lt. abs( a(k,j) ) ) temp = abs( a(k,j) )
!
            work1(k) = 0.0
!
         enddo
!
!	scale and V and compute scale factor
!
         alpha = 0.0
!
         do k = j , n
!
            work1(k) = a(k,j) / temp
            alpha    = alpha + work1(k) ** 2
!
         enddo
!
         work1(j) = work1(j) + sign( 1.0 , work1(j) ) * sqrt( alpha )
         beta     = 2.0 / dot_product( work1 , work1 )
!
!	2VV'/V'V
!
         do kk = 1 , n
!
            do jj = 1 , n
!
               work2(jj,kk) = beta * work1(jj) * work1(kk)
!
            enddo
!
         enddo
!
!	Pk=I - 2VV'/V'V, form rank on update matrix
!
         do jj = 1 , n
!
            do ii = 1 , n
!
               if ( ii .eq. jj ) then
!		  
                  work2(ii,jj) = 1.0 - work2(ii,jj)
!
               else
!
                  work2(ii,jj) = - work2(ii,jj)
!
               endif
!
            enddo
!
         enddo
!
!	Rank One Update of A and b, (I - 2VV'/V'V)A 
!
         aw = matmul( work2 , a ) ! PkPk-1...P1A
         bw = matmul( work2 , b ) ! PkPk-1...P1b
!
         do kk = 1 , n
!
            b(kk) = bw(kk)
!
            do k = 1 , nc
!
               a(kk,k) = aw(kk,k)
!
            enddo
!
         enddo
!
         do kk = j + 1 , n
!
            a(kk,j) = 0.0
!
         enddo
!
      enddo
!
!	Back solve for final coefficent fit
!
      work1(nc) = b(nc) / a(nc,nc)
!
      do k = nc - 1 , 1 , - 1
!
         temp = 0.0
!
         do j = k + 1 , nc
!
            temp = temp + a(k,j) * work1(j)
!
         enddo
!
         work1(k) = ( b(k) - temp ) / a(k,k)
!
      enddo
!
      do kk = 1 , nc
!
         b(kk) = work1(kk)
!
      enddo
!
      return
      end
