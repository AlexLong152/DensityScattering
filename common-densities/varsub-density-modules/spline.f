      MODULE spline
       USE precision
       PRIVATE

       INTERFACE cubbasis 
        MODULE PROCEDURE cubbasis_dp,cubbasis_sp
       END INTERFACE

       INTERFACE cubherm
        MODULE PROCEDURE cubherm_dp,cubherm_sp
       END INTERFACE

       INTERFACE cubfast
        MODULE PROCEDURE cubfast_dp,cubfast_sp
       END INTERFACE

       INTERFACE chebspline
        MODULE PROCEDURE chebspline_dp,chebspline_sp
       END INTERFACE

       INTERFACE cubbasis 
        MODULE PROCEDURE linspl_dp,linspl_sp
       END INTERFACE

       PUBLIC cubbasis,cubherm,cubfast,chebspline,linspl 

       PUBLIC cubbasis_dp,cubbasis_sp
       
       PUBLIC cubherm_dp,cubherm_sp
       
       PUBLIC cubfast_dp,cubfast_sp
       
       PUBLIC chebspline_dp,chebspline_sp

       PUBLIC linspl_dp,linspl_sp
       
       PUBLIC cubherm_part_dp,cubfast_part_dp
       
       CONTAINS

C> subroutine calculates basis splines for double precision arrays
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points 
C> @param[out] spl spl(m,n) spline elements
        SUBROUTINE cubbasis_dp(xold,n,xnew,m,spl)

C *** Die Routine bereitet Basissplines mit Hilfe der SPLMOD und SPLFAK vor.
C *** Die Parameter sind wie bei cubherm definiert
C *** Die Dimension ist jetzt anders : spl(m,n) statt spl(m,4)
C *** Der Parameter index ist nicht noetig !

         IMPLICIT NONE

         REAL(dpreal) xold(n),xnew(m),spl(m,n)
         INTEGER n,m
       
         REAL(dpreal) pspl(n)
         INTEGER i,j

         CALL splfak_dp(xold,n)

         DO i=1,m
          IF((xnew(i).GE.xold(1)).AND.(xnew(i).LE.xold(n))) THEN
            CALL splmod_dp(xold,n,xnew(i),pspl)
            DO j=1,n
C***  hier kann Interpolation aendern
             spl(i,j)=pspl(j)
!     X            *xold(j)/xnew(i)
            END DO
           ELSE
            DO j=1,n
             spl(i,j)=0.0
            END DO
            IF(xnew(i).LE.xold(1)) THEN
              spl(i,1)=(xold(2)-xnew(i))/(xold(2)-xold(1))
!     X                      *xold(1)/xnew(i)
              spl(i,2)=(xnew(i)-xold(1))/(xold(2)-xold(1))
!     X                      *xold(2)/xnew(i)
            END IF
          END IF   
         END DO

        END SUBROUTINE cubbasis_dp



C> subroutine calculates cubic hermitian splines for double precision arrays
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points 
C> @param[out] spl spl(m,4) spline elements for some over 4 nearest old grid points 
C> @param[out] index(m,4)   index to relevant old grid points 

      subroutine cubherm_dp(xold,n,xnew,m,spl,index)
      implicit real(dpreal) (a-h,o-z)
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(m,4),index(m,4)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubherm'
cc      if (m.gt.mmax) stop'mmax to small in cubherm'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(j,2)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(j,2)=i
11    continue
      do 12 j=1,m
       index(j,2)=min(index(j,2),n-1)
       index(j,1)=index(j,2)-1
       index(j,3)=index(j,2)+1
       index(j,4)=index(j,2)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(j,1).eq.0) index(j,1)=3
       if (index(j,4).eq.n+1) index(j,4)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n).and.enough) then
        if(xnew(j).lt.xold(2)) then  ! replace by linear interpolation for values close to threshold
         index(j,2)=1
         index(j,3)=2
         index(j,1)=1
         index(j,4)=2
         spl(j,2)=(xnew(j)-xold(2))/(xold(1)-xold(2))
         spl(j,3)=(xnew(j)-xold(1))/(xold(2)-xold(1))
         spl(j,1)=0.0
         spl(j,4)=0.0
        else
         i0=index(j,1)
         i1=index(j,2)
         i2=index(j,3)
         i3=index(j,4)
         x0=xold(i0)
         x1=xold(i1)
         x2=xold(i2)
         x3=xold(i3)
c
c      Factors for the derivatives
c
         d10=x1-x0
         d21=x2-x1
         d32=x3-x2
         d20=x2-x0
         d31=x3-x1
         dfak13=(d21/d10-d10/d21)/d20
         dfak14=-d32/(d21*d31)
         dfak23=d10/(d21*d20)
         dfak24=(d32/d21-d21/d32)/d31
         dfak03=-d21/(d10*d20)
         dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
         xn=xnew(j)
         dn1=xn-x1
         d2n=x2-xn
         phidiv=1./(d21*d21*d21)
         phi1=d2n*d2n*phidiv*(d21+2.*dn1)
         phi2=dn1*dn1*phidiv*(d21+2.*d2n)
         phidiv=phidiv*d21*dn1*d2n
         phi3=phidiv*d2n
         phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
         spl(j,2)=phi1+phi3*dfak13+phi4*dfak14
         spl(j,3)=phi2+phi3*dfak23+phi4*dfak24
         spl(j,1)=phi3*dfak03
         spl(j,4)=phi4*dfak34
        end if
       else
        index(j,2)=1
        index(j,3)=1
        index(j,1)=1
        index(j,4)=1
        spl(j,2)=0.0
        spl(j,3)=0.0
        spl(j,1)=0.0
        spl(j,4)=0.0
        if (n.eq.1) spl(j,2)=1.0
       endif

C *** AN inserted interpolation of q*f(q)
c       spl(j,2)=spl(j,2)*x1/xn  
c       spl(j,3)=spl(j,3)*x2/xn
c       spl(j,1)=spl(j,1)*x0/xn
c       spl(j,4)=spl(j,4)*x3/xn
C *** AN

20    continue
      end subroutine

C> subroutine calculates cubic hermitian splines for double precision arrays.
C> Note that the routine is similar to #cubherm_dp. Only the ordering 
C> of indices of the output arrays is different.       
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points 
C> @param[out] spl spl(4,m) spline elements for some over 4 nearest old grid points 
C> @param[out] index(4,m)   index to relevant old grid points
      
      subroutine cubfast_dp(xold,n,xnew,m,spl,index)
      implicit real(dpreal) (a-h,o-z)
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(4,m),index(4,m)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubfast'
cc      if (m.gt.mmax) stop'mmax to small in cubfast'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(2,j)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(2,j)=i
11    continue
      do 12 j=1,m
       index(2,j)=min(index(2,j),n-1)
       index(1,j)=index(2,j)-1
       index(3,j)=index(2,j)+1
       index(4,j)=index(2,j)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(1,j).eq.0) index(1,j)=3
       if (index(4,j).eq.n+1) index(4,j)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n)*1.01.and.enough) then
        i0=index(1,j)
        i1=index(2,j)
        i2=index(3,j)
        i3=index(4,j)
        x0=xold(i0)
        x1=xold(i1)
        x2=xold(i2)
        x3=xold(i3)
c
c      Factors for the derivatives
c
        d10=x1-x0
        d21=x2-x1
        d32=x3-x2
        d20=x2-x0
        d31=x3-x1
        dfak13=(d21/d10-d10/d21)/d20
        dfak14=-d32/(d21*d31)
        dfak23=d10/(d21*d20)
        dfak24=(d32/d21-d21/d32)/d31
        dfak03=-d21/(d10*d20)
        dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
        xn=xnew(j)
        dn1=xn-x1
        d2n=x2-xn
        phidiv=1./(d21*d21*d21)
        phi1=d2n*d2n*phidiv*(d21+2.*dn1)
        phi2=dn1*dn1*phidiv*(d21+2.*d2n)
        phidiv=phidiv*d21*dn1*d2n
        phi3=phidiv*d2n
        phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
        spl(2,j)=phi1+phi3*dfak13+phi4*dfak14
        spl(3,j)=phi2+phi3*dfak23+phi4*dfak24
        spl(1,j)=phi3*dfak03
        spl(4,j)=phi4*dfak34
       else
        index(2,j)=1
        index(3,j)=1
        index(1,j)=1
        index(4,j)=1
        spl(2,j)=0.0
        spl(3,j)=0.0
        spl(1,j)=0.0
        spl(4,j)=0.0
        if (n.eq.1) spl(2,j)=1.0
       endif

C *** AN inserted interpolation of q*f(q)
c       spl(2,j)=spl(2,j)*x1/xn  
c       spl(3,j)=spl(3,j)*x2/xn
c       spl(1,j)=spl(1,j)*x0/xn
c       spl(4,j)=spl(4,j)*x3/xn
C *** AN

20    continue
      end subroutine cubfast_dp

C> subroutine calculates cubic hermitian splines for double precision arrays.
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points
C> @param[in] xbound boundary of the interval where the function could be discontinous       
C> @param[out] spl spl(m,4) spline elements for some over 4 nearest old grid points 
C> @param[out] index(m,4)   index to relevant old grid points  

      subroutine cubherm_part_dp(xold,n,xnew,m,spl,index,xbound)
C       this routine use cubherm interpolation independently in two 
C       separate intervals to avoid a possible non-continues 
C       derivative at the point xbound
C       this test version assumes a fixed value for the problematic 
C       momentum given by the parameter below
      implicit real(dpreal) (a-h,o-z)
      
      REAL(dpreal) :: xbound  
      REAL(dpreal) :: xn,xold 
      INTEGER shift,indx 
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(m,4),index(m,4)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubherm'
cc      if (m.gt.mmax) stop'mmax to small in cubherm'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(j,2)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(j,2)=i
11    continue

C  here the second grid point of the four is fixed 
C  assuming that there is no problematic grid point 

C  adjust if necessary  
      do j=1,m
       xn=xnew(j)
       shift=0
       indx=index(j,2)
C  no adjustments when interval close to boundaries 
       if(indx.GE.3.AND.indx.LT.n-2) then 
C    first case   xnew < xbound and xold(index(j,2)+2 or +1) > xbound
        if((xn.LE.xbound).AND.xold(indx+1).GT.xbound) shift=shift-1 
        if((xn.LE.xbound).AND.xold(indx+2).GT.xbound) shift=shift-1 
C    first case   xnew > xbound and xold(index(j,2)+0 or -1)  < xbound       
        if((xn.GT.xbound).AND.xold(indx)  .LT.xbound) shift=shift+1 
        if((xn.GT.xbound).AND.xold(indx-1).LT.xbound) shift=shift+1 
       end if
       index(j,2)=indx+shift
      end do

C now the interval of 4 points should be entirely in one of the two regions      

      do 12 j=1,m
       index(j,2)=min(index(j,2),n-1)
       index(j,1)=index(j,2)-1
       index(j,3)=index(j,2)+1
       index(j,4)=index(j,2)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(j,1).le.0) index(j,1)=3
       if (index(j,4).ge.n+1) index(j,4)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n)*1.01.and.enough) then
        if(xnew(j).lt.xold(2)) then  ! replace by linear interpolation for values close to threshold
         index(j,2)=1
         index(j,3)=2
         index(j,1)=1
         index(j,4)=2
         spl(j,2)=(xnew(j)-xold(2))/(xold(1)-xold(2))
         spl(j,3)=(xnew(j)-xold(1))/(xold(2)-xold(1))
         spl(j,1)=0.0
         spl(j,4)=0.0
        else
         i0=index(j,1)
         i1=index(j,2)
         i2=index(j,3)
         i3=index(j,4)
         x0=xold(i0)
         x1=xold(i1)
         x2=xold(i2)
         x3=xold(i3)
c
c      Factors for the derivatives
c
         d10=x1-x0
         d21=x2-x1
         d32=x3-x2
         d20=x2-x0
         d31=x3-x1
         dfak13=(d21/d10-d10/d21)/d20
         dfak14=-d32/(d21*d31)
         dfak23=d10/(d21*d20)
         dfak24=(d32/d21-d21/d32)/d31
         dfak03=-d21/(d10*d20)
         dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
         xn=xnew(j)
         dn1=xn-x1
         d2n=x2-xn
         phidiv=1./(d21*d21*d21)
         phi1=d2n*d2n*phidiv*(d21+2.*dn1)
         phi2=dn1*dn1*phidiv*(d21+2.*d2n)
         phidiv=phidiv*d21*dn1*d2n
         phi3=phidiv*d2n
         phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
         spl(j,2)=phi1+phi3*dfak13+phi4*dfak14
         spl(j,3)=phi2+phi3*dfak23+phi4*dfak24
         spl(j,1)=phi3*dfak03
         spl(j,4)=phi4*dfak34
        end if
       else
        index(j,2)=1
        index(j,3)=1
        index(j,1)=1
        index(j,4)=1
        spl(j,2)=0.0
        spl(j,3)=0.0
        spl(j,1)=0.0
        spl(j,4)=0.0
        if (n.eq.1) spl(j,2)=1.0
       endif

C *** AN inserted interpolation of q*f(q)
c       spl(j,2)=spl(j,2)*x1/xn  
c       spl(j,3)=spl(j,3)*x2/xn
c       spl(j,1)=spl(j,1)*x0/xn
c       spl(j,4)=spl(j,4)*x3/xn
C *** AN

20    continue
      end subroutine cubherm_part_dp

C> subroutine calculates cubic hermitian splines for double precision arrays.
C> The routine is similar to #cubherm_part_dp except that ordering of indices 
C> of arrays differs.       
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points
C> @param[in] xbound boundary of the interval where the function could be discontinous       
C> @param[out] spl spl(4,m) spline elements for some over 4 nearest old grid points 
C> @param[out] index(4,m)   index to relevant old grid points 

      subroutine cubfast_part_dp(xold,n,xnew,m,spl,index,xbound)
C       this routine use cubherm interpolation independently in two 
C       separate intervals to avoid a possible non-continues 
C       derivative at the point xbound
C       this test version assumes a fixed value for the problematic 
C       momentum given by the parameter below
      implicit real(dpreal) (a-h,o-z)
      
      REAL(dpreal) :: xbound  
      REAL(dpreal) :: xn,xold 
      INTEGER shift,indx 
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(4,m),index(4,m)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubfast'
cc      if (m.gt.mmax) stop'mmax to small in cubfast'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(2,j)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(2,j)=i
11    continue

C  here the second grid point of the four is fixed 
C  assuming that there is no problematic grid point 

C  adjust if necessary  
      do j=1,m
       xn=xnew(j)
       shift=0
       indx=index(2,j)
C  no adjustments when interval close to boundaries 
       if(indx.GE.2.AND.indx.LT.n-2) then 
C         first case   xnew < xbound and xold(index(2,j)+2 or +1) > xbound
        if((xn.LE.xbound).AND.xold(indx+1).GT.xbound) shift=shift-1 
        if((xn.LE.xbound).AND.xold(indx+2).GT.xbound) shift=shift-1 
C         first case   xnew > xbound and xold(index(2,j)+0 or -1)  < xbound       
        if((xn.GT.xbound).AND.xold(indx)  .LT.xbound) shift=shift+1 
        if((xn.GT.xbound).AND.xold(indx-1).LT.xbound) shift=shift+1 
       end if
       index(2,j)=indx+shift
      end do

C now the interval of 4 points should be entirely in one of the two regions      


      do 12 j=1,m
       index(2,j)=min(index(2,j),n-1)
       index(1,j)=index(2,j)-1
       index(3,j)=index(2,j)+1
       index(4,j)=index(2,j)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(1,j).le.0) index(1,j)=3
       if (index(4,j).ge.n+1) index(4,j)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n)*1.01.and.enough) then
        if(xnew(j).lt.xold(2)) then  ! replace by linear interpolation for values close to threshold
         index(2,j)=1
         index(3,j)=2
         index(1,j)=1
         index(4,j)=2
         spl(2,j)=(xnew(j)-xold(2))/(xold(1)-xold(2))
         spl(3,j)=(xnew(j)-xold(1))/(xold(2)-xold(1))
         spl(1,j)=0.0
         spl(4,j)=0.0
        else
         i0=index(1,j)
         i1=index(2,j)
         i2=index(3,j)
         i3=index(4,j)
         x0=xold(i0)
         x1=xold(i1)
         x2=xold(i2)
         x3=xold(i3)
c
c      Factors for the derivatives
c
         d10=x1-x0
         d21=x2-x1
         d32=x3-x2
         d20=x2-x0
         d31=x3-x1
         dfak13=(d21/d10-d10/d21)/d20
         dfak14=-d32/(d21*d31)
         dfak23=d10/(d21*d20)
         dfak24=(d32/d21-d21/d32)/d31
         dfak03=-d21/(d10*d20)
         dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
         xn=xnew(j)
         dn1=xn-x1
         d2n=x2-xn
         phidiv=1./(d21*d21*d21)
         phi1=d2n*d2n*phidiv*(d21+2.*dn1)
         phi2=dn1*dn1*phidiv*(d21+2.*d2n)
         phidiv=phidiv*d21*dn1*d2n
         phi3=phidiv*d2n
         phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
         spl(2,j)=phi1+phi3*dfak13+phi4*dfak14
         spl(3,j)=phi2+phi3*dfak23+phi4*dfak24
         spl(1,j)=phi3*dfak03
         spl(4,j)=phi4*dfak34
        end if
       else
        index(2,j)=1
        index(3,j)=1
        index(1,j)=1
        index(4,j)=1
        spl(2,j)=0.0
        spl(3,j)=0.0
        spl(1,j)=0.0
        spl(4,j)=0.0
        if (n.eq.1) spl(2,j)=1.0
      endif

C *** AN inserted interpolation of q*f(q)
c       spl(2,j)=spl(2,j)*x1/xn  
c       spl(3,j)=spl(3,j)*x2/xn
c       spl(1,j)=spl(1,j)*x0/xn
c       spl(4,j)=spl(4,j)*x3/xn
C *** AN

20    continue
      end subroutine cubfast_part_dp

      
      
      SUBROUTINE SPLFAK_dp(X,N)
      implicit real(dpreal) (a-h,o-z)
C     =================
C     THE 2ND DERIVATIVES AT THE BOUNDARIES MUST BE ZERO !
C     THE SPLINE IS BUILT UP USING ALL THE N MESHPOINTS .
C
      PARAMETER(NQ=500)
      DIMENSION X(N)
      DIMENSION HI(NQ),U(NQ)
      REAL(dpreal) PI
      COMMON/FAKTORDP/FAK1(NQ,NQ),FAK2(NQ,NQ),FAK3(NQ,NQ),
     C              Q(NQ,NQ),C(NQ,NQ)

      IF(N.GT.NQ) STOP 'SPLFAK: NQ nicht gross genug'
      IF(N.LT.2) RETURN
      U(1)=0.
      HI(2)=X(2)-X(1)
      N1=N-1
      DO 10 I=2,N1
      AAX=X(I+1)-X(I)
      HI(I+1)=AAX
      BX=X(I+1)-X(I-1)
      CX=X(I)-X(I-1)
      AL=AAX/BX
      AM=1.-AL
      PI=1./(2.-AM*U(I-1))
      U(I)=AL*PI
      DO 20 J=1,N
      Q(1,J)=0.
      H1=0.
      H2=0.
      H3=0.
      IF(J.EQ.I-1) H1=1./(CX*BX)
      IF(J.EQ.I) H2=1./(CX*AAX)
      IF(J.EQ.I+1) H3=1./(AAX*BX)
      Q(I,J)=-PI*(AM*Q(I-1,J)-H1+H2-H3)
  20  CONTINUE
  10  CONTINUE
      N2=N+1
      N3=N+2
      DO 30 K=3,N2
      J1=N3-K
      H1=1./HI(J1+1)
      DO 40 L=1,N
      C(N,L)=0.
      C(J1,L)=Q(J1,L)-C(J1+1,L)*U(J1)
      FAK1(J1,L)=-HI(J1+1)*(2.*C(J1,L)+C(J1+1,L))
      IF(L.EQ.J1) FAK1(J1,L)=FAK1(J1,L)-H1
      IF(L.EQ.J1+1) FAK1(J1,L)=FAK1(J1,L)+H1
      FAK2(J1,L)=3.*C(J1,L)
      FAK3(J1,L)=(C(J1+1,L)-C(J1,L))*H1
      FAK1(N,L)=0.
      FAK2(N,L)=0.
      FAK3(N,L)=0.
  40  CONTINUE
  30  CONTINUE
      RETURN
      END SUBROUTINE
      
      SUBROUTINE SPLMOD_dp(X,N,XA,SPL)
      implicit real(dpreal) (a-h,o-z)
C     =================
C     X  : ARRAY OF N WHERE THE FUNCTION TO BE INTERPOLATED IS KNOWN
C          (USUALLY MESHPPOINTS)
C     XA : POINT OF INTERPOLATION
C
C     F(XA) = SUM(I=1,N) [ SPL(I)*F(X(I)) ]
C
C     (THIS SUMMATION HAS TO BE DONE IN THE CALLING PROGRAM ! )
C
      PARAMETER(NQ=500)
      DIMENSION X(N)
      DIMENSION SPL(N)
      COMMON/FAKTORDP/FAK1(NQ,NQ),FAK2(NQ,NQ),FAK3(NQ,NQ),
     C              QDUM(NQ,NQ),CDUM(NQ,NQ)
      I=0
   2  I=I+1
      IF(I.GT.N) GO TO 3
      IF(XA.GE.X(I)) GO TO 2
   3  I1=I-1
      DX=XA-X(I1)
      DO 10 J=1,N
      SPL(J)=((FAK3(I1,J)*DX+FAK2(I1,J))*DX+FAK1(I1,J))*DX
      IF(J.EQ.I1) SPL(J)=SPL(J)+1.
  10  CONTINUE
      RETURN
      END SUBROUTINE
     
        SUBROUTINE chebspline_sp(a,b,n,xnew,m,spl)

C     *** interpolation based on the chebyshev approximation 
C     *** of the wave function. grid points are transformed 
C     *** Chebyshev points mit grid parameters a and b 
C     *** and n mesh points. xnew(m) are m new shifted 
C     *** mesh points
        
         IMPLICIT NONE

         REAL(spreal) a,b,xnew(m),spl(n,m)
         INTEGER n,m
       
         REAL(spreal) x,y,pi
         INTEGER i,j,k

         pi=acos(-1.0)

         DO i=1,m
          y=xnew(i)
          x=-(1.0-y/a)/(-(1.0/a-2.0/b)*y-1.0)

          DO k=1,n
           spl(k,i)=0.0           
           IF(y.GE.0.0.AND.y.LE.b) THEN
             DO j=1,n
              spl(k,i)=spl(k,i)
     X             +2.0/n*(cos(pi*(j-1.0)*(k-0.5)/n)
     X             *cos((j-1)*acos(x)))
             END DO

             spl(k,i)=spl(k,i)-1.0/n
           END IF
          END DO
         END DO

        END SUBROUTINE

        SUBROUTINE chebspline_dp(a,b,n,xnew,m,spl)

C     *** interpolation based on the chebyshev approximation 
C     *** of the wave function. grid points are transformed 
C     *** Chebyshev points mit grid parameters a and b 
C     *** and n mesh points. xnew(m) are m new shifted 
C     *** mesh points
        
         IMPLICIT NONE

         REAL(dpreal) a,b,xnew(m),spl(n,m)
         INTEGER n,m
       
         REAL(dpreal) x,y,pi
         INTEGER i,j,k

         pi=acos(-1.0_dpreal)

         DO i=1,m
          y=xnew(i)
          x=-(1.0_dpreal-y/a)/
     $       (-(1.0_dpreal/a-2.0_dpreal/b)*y-1.0_dpreal)

          DO k=1,n
           spl(k,i)=0.0           
           IF(y.GE.0.0.AND.y.LE.b) THEN
             DO j=1,n
              spl(k,i)=spl(k,i)
     X           +2.0_dpreal/n*(cos(pi*(j-1.0_dpreal)*(k-0.5_dpreal)/n)
     X             *cos((j-1)*acos(x)))
             END DO

             spl(k,i)=spl(k,i)-1.0_dpreal/n
           END IF
          END DO
         END DO

        END SUBROUTINE

C> subroutine calculates basis splines for single precision arrays
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points 
C> @param[out] spl spl(m,n) spline elements 

        SUBROUTINE cubbasis_sp(xold,n,xnew,m,spl)

C *** Die Routine bereitet Basissplines mit Hilfe der SPLMOD und SPLFAK vor.
C *** Die Parameter sind wie bei cubherm definiert
C *** Die Dimension ist jetzt anders : spl(m,n) statt spl(m,4)
C *** Der Parameter index ist nicht noetig !

         IMPLICIT NONE

         REAL(spreal) xold(n),xnew(m),spl(m,n)
         INTEGER n,m
       
         REAL(spreal) pspl(n)
         INTEGER i,j

         CALL splfak_sp(xold,n)

         DO i=1,m
          IF((xnew(i).GE.xold(1)).AND.(xnew(i).LE.xold(n)*1.01)) THEN
            CALL splmod_sp(xold,n,xnew(i),pspl)
            DO j=1,n
C***  hier kann Interpolation aendern
             spl(i,j)=pspl(j)
     X            *xold(j)/xnew(i)
            END DO
           ELSE
            DO j=1,n
             spl(i,j)=0.0
            END DO
            IF(xnew(i).LE.xold(1)) THEN
              spl(i,1)=(xold(2)-xnew(i))/(xold(2)-xold(1))
     X                      *xold(1)/xnew(i)
              spl(i,2)=(xnew(i)-xold(1))/(xold(2)-xold(1))
     X                      *xold(2)/xnew(i)
            END IF
          END IF   
         END DO

        END SUBROUTINE cubbasis_sp


C> subroutine calculates cubic hermitian splines for single precision arrays
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points 
C> @param[out] spl spl(m,4) spline elements for some over 4 nearest old grid points 
C> @param[out] index(m,4)   index to relevant old grid points

      subroutine cubherm_sp(xold,n,xnew,m,spl,index)
      implicit real(spreal) (a-h,o-z)
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(m,4),index(m,4)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubherm'
cc      if (m.gt.mmax) stop'mmax to small in cubherm'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(j,2)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(j,2)=i
11    continue
      do 12 j=1,m
       index(j,2)=min(index(j,2),n-1)
       index(j,1)=index(j,2)-1
       index(j,3)=index(j,2)+1
       index(j,4)=index(j,2)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(j,1).eq.0) index(j,1)=3
       if (index(j,4).eq.n+1) index(j,4)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n)*1.01.and.enough) then
        i0=index(j,1)
        i1=index(j,2)
        i2=index(j,3)
        i3=index(j,4)
        x0=xold(i0)
        x1=xold(i1)
        x2=xold(i2)
        x3=xold(i3)
c
c      Factors for the derivatives
c
        d10=x1-x0
        d21=x2-x1
        d32=x3-x2
        d20=x2-x0
        d31=x3-x1
        dfak13=(d21/d10-d10/d21)/d20
        dfak14=-d32/(d21*d31)
        dfak23=d10/(d21*d20)
        dfak24=(d32/d21-d21/d32)/d31
        dfak03=-d21/(d10*d20)
        dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
        xn=xnew(j)
        dn1=xn-x1
        d2n=x2-xn
        phidiv=1./(d21*d21*d21)
        phi1=d2n*d2n*phidiv*(d21+2.*dn1)
        phi2=dn1*dn1*phidiv*(d21+2.*d2n)
        phidiv=phidiv*d21*dn1*d2n
        phi3=phidiv*d2n
        phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
        spl(j,2)=phi1+phi3*dfak13+phi4*dfak14
        spl(j,3)=phi2+phi3*dfak23+phi4*dfak24
        spl(j,1)=phi3*dfak03
        spl(j,4)=phi4*dfak34
       else
        index(j,2)=1
        index(j,3)=1
        index(j,1)=1
        index(j,4)=1
        spl(j,2)=0.0
        spl(j,3)=0.0
        spl(j,1)=0.0
        spl(j,4)=0.0
        if (n.eq.1) spl(j,2)=1.0
       endif

C *** AN inserted interpolation of q*f(q)
c       spl(j,2)=spl(j,2)*x1/xn  
c       spl(j,3)=spl(j,3)*x2/xn
c       spl(j,1)=spl(j,1)*x0/xn
c       spl(j,4)=spl(j,4)*x3/xn
C *** AN

20    continue
      end subroutine

C> subroutine calculates cubic hermitian splines for single precision arrays.
C> Note that the routine is similar to #cubherm_sp. Only the ordering 
C> of indices of the output arrays is different.       
C> @param[in] xold xold(n) old grid points 
C> @param[in] n number of old grid points
C> @param[in] xnew xnew(m) new points 
C> @param[in] m number of new points 
C> @param[out] spl spl(4,m) spline elements for some over 4 nearest old grid points 
C> @param[out] index(4,m)   index to relevant old grid points
      
      subroutine cubfast_sp(xold,n,xnew,m,spl,index)
      implicit real(spreal) (a-h,o-z)
c
c     This subroutine prepares the interpolation of a function
c     given at the n grid points xold to the m points xnew.
c     The interpolating functions are cubic hermitian splines.
c     The first derivatives at the grid points are taken from a
c     parabola through the actual grid point and its two neighbours.
c     For the end points the parabola is taken through the two
c     right or left neighbours, resp.
c
c     In the calling routine you still have to take the sum i=1,4
c     ynew(j)=sum_i spl(j,i)*yold(index(j,i))
c
c     by Dirk Hueber, 08.02.1996
c
cc    parameter (nmax=100,mmax=100)
      dimension xold(n),xnew(m),spl(4,m),index(4,m)
      logical enough
      enough=n.ge.3
cc      if (n.gt.nmax) stop'nmax to small in cubfast'
cc      if (m.gt.mmax) stop'mmax to small in cubfast'
c
c     evaluation of the indices
c
      do 10 j=1,m
c
c     If xnew(j) is less than xold(1), we extrapolate using the first
c     two grid points xold(1) and xold(2).
c
       index(2,j)=1
10    continue
      do 11 i=1,n
       do 11 j=1,m
        if (xnew(j).gt.xold(i)) index(2,j)=i
11    continue
      do 12 j=1,m
       index(2,j)=min(index(2,j),n-1)
       index(1,j)=index(2,j)-1
       index(3,j)=index(2,j)+1
       index(4,j)=index(2,j)+2
c
c     Indices 1 and 4 are used only for the ewertcation of the derivatives.
c     The following settings provide the derivatives at the borders of the
c     grid xold.
c
       if (index(1,j).eq.0) index(1,j)=3
       if (index(4,j).eq.n+1) index(4,j)=n-2
12    continue
      do 20 j=1,m
c
c     We don't extrapolate to the right!
c
       if (xnew(j).le.xold(n)*1.01.and.enough) then
        i0=index(1,j)
        i1=index(2,j)
        i2=index(3,j)
        i3=index(4,j)
        x0=xold(i0)
        x1=xold(i1)
        x2=xold(i2)
        x3=xold(i3)
c
c      Factors for the derivatives
c
        d10=x1-x0
        d21=x2-x1
        d32=x3-x2
        d20=x2-x0
        d31=x3-x1
        dfak13=(d21/d10-d10/d21)/d20
        dfak14=-d32/(d21*d31)
        dfak23=d10/(d21*d20)
        dfak24=(d32/d21-d21/d32)/d31
        dfak03=-d21/(d10*d20)
        dfak34=d21/(d32*d31)
c
c     the cubic hermitian splines
c
        xn=xnew(j)
        dn1=xn-x1
        d2n=x2-xn
        phidiv=1./(d21*d21*d21)
        phi1=d2n*d2n*phidiv*(d21+2.*dn1)
        phi2=dn1*dn1*phidiv*(d21+2.*d2n)
        phidiv=phidiv*d21*dn1*d2n
        phi3=phidiv*d2n
        phi4=-phidiv*dn1
c
c     combining everything to the final factors
c
        spl(2,j)=phi1+phi3*dfak13+phi4*dfak14
        spl(3,j)=phi2+phi3*dfak23+phi4*dfak24
        spl(1,j)=phi3*dfak03
        spl(4,j)=phi4*dfak34
       else
        index(2,j)=1
        index(3,j)=1
        index(1,j)=1
        index(4,j)=1
        spl(2,j)=0.0
        spl(3,j)=0.0
        spl(1,j)=0.0
        spl(4,j)=0.0
        if (n.eq.1) spl(2,j)=1.0
       endif

C *** AN inserted interpolation of q*f(q)
c       spl(2,j)=spl(2,j)*x1/xn  
c       spl(3,j)=spl(3,j)*x2/xn
c       spl(1,j)=spl(1,j)*x0/xn
c       spl(4,j)=spl(4,j)*x3/xn
C *** AN

20    continue
      end subroutine cubfast_sp



      SUBROUTINE SPLFAK_sp(X,N)
      implicit real(spreal) (a-h,o-z)
C     =================
C     THE 2ND DERIVATIVES AT THE BOUNDARIES MUST BE ZERO !
C     THE SPLINE IS BUILT UP USING ALL THE N MESHPOINTS .
C
      PARAMETER(NQ=500)
      DIMENSION X(N)
      DIMENSION HI(NQ),U(NQ)
      REAL(spreal) PI
      COMMON/FAKTORSP/FAK1(NQ,NQ),FAK2(NQ,NQ),FAK3(NQ,NQ),
     C              Q(NQ,NQ),C(NQ,NQ)

      IF(N.GT.NQ) STOP 'SPLFAK: NQ nicht gross genug'
      IF(N.LT.2) RETURN
      U(1)=0.
      HI(2)=X(2)-X(1)
      N1=N-1
      DO 10 I=2,N1
      AAX=X(I+1)-X(I)
      HI(I+1)=AAX
      BX=X(I+1)-X(I-1)
      CX=X(I)-X(I-1)
      AL=AAX/BX
      AM=1.-AL
      PI=1./(2.-AM*U(I-1))
      U(I)=AL*PI
      DO 20 J=1,N
      Q(1,J)=0.
      H1=0.
      H2=0.
      H3=0.
      IF(J.EQ.I-1) H1=1./(CX*BX)
      IF(J.EQ.I) H2=1./(CX*AAX)
      IF(J.EQ.I+1) H3=1./(AAX*BX)
      Q(I,J)=-PI*(AM*Q(I-1,J)-H1+H2-H3)
  20  CONTINUE
  10  CONTINUE
      N2=N+1
      N3=N+2
      DO 30 K=3,N2
      J1=N3-K
      H1=1./HI(J1+1)
      DO 40 L=1,N
      C(N,L)=0.
      C(J1,L)=Q(J1,L)-C(J1+1,L)*U(J1)
      FAK1(J1,L)=-HI(J1+1)*(2.*C(J1,L)+C(J1+1,L))
      IF(L.EQ.J1) FAK1(J1,L)=FAK1(J1,L)-H1
      IF(L.EQ.J1+1) FAK1(J1,L)=FAK1(J1,L)+H1
      FAK2(J1,L)=3.*C(J1,L)
      FAK3(J1,L)=(C(J1+1,L)-C(J1,L))*H1
      FAK1(N,L)=0.
      FAK2(N,L)=0.
      FAK3(N,L)=0.
  40  CONTINUE
  30  CONTINUE
      RETURN
      END SUBROUTINE
      
      SUBROUTINE SPLMOD_sp(X,N,XA,SPL)
      implicit real(spreal) (a-h,o-z)
C     =================
C     X  : ARRAY OF N WHERE THE FUNCTION TO BE INTERPOLATED IS KNOWN
C          (USUALLY MESHPPOINTS)
C     XA : POINT OF INTERPOLATION
C
C     F(XA) = SUM(I=1,N) [ SPL(I)*F(X(I)) ]
C
C     (THIS SUMMATION HAS TO BE DONE IN THE CALLING PROGRAM ! )
C
      PARAMETER(NQ=500)
      DIMENSION X(N)
      DIMENSION SPL(N)
      COMMON/FAKTORSP/FAK1(NQ,NQ),FAK2(NQ,NQ),FAK3(NQ,NQ),
     C              QDUM(NQ,NQ),CDUM(NQ,NQ)
      I=0
   2  I=I+1
      IF(I.GT.N) GO TO 3
      IF(XA.GE.X(I)) GO TO 2
   3  I1=I-1
      DX=XA-X(I1)
      DO 10 J=1,N
      SPL(J)=((FAK3(I1,J)*DX+FAK2(I1,J))*DX+FAK1(I1,J))*DX
      IF(J.EQ.I1) SPL(J)=SPL(J)+1.
  10  CONTINUE
      RETURN
      END SUBROUTINE
     

C     linear interpolation
C     assuming that the function is zero for xnew > xold(n) 
C     and assuming extrapolation from interval xold(1),xold(2) 
C     to xnew > xold(1) 
C     this subroutine assumes that the xold points are sorted 
C     according to their seize
C     single precision version 

      subroutine linspl_sp(xold,n,xnew,m,spl,indx)
       IMPLICIT NONE
       INTEGER n,m
       REAL(spreal) xold(n),xnew(m)
       REAL(spreal) spl(m,2)
       INTEGER indx(m,2)
       INTEGER i,j

C     we need at least two points
       IF(n.LT.2) STOP 'linear interpolation impossible'
C     first find correct indices for each xnew 
       DO j=1,m
        DO i=1,n
         IF(xold(i).GT.xnew(j)) exit 
        END DO                  ! i
        indx(j,2)=min(max(i,2),n) ! second point at least 2 but .le. n
        indx(j,1)=indx(j,2)-1   ! one less
       END DO                   ! j

C     fix the spl weights 
       
       DO j=1,m
        IF(xnew(j).LE.xold(n)*1.01) THEN
         spl(j,1)=(xold(indx(j,2))-xnew(j))
     $           /(xold(indx(j,2))-xold(indx(j,1)))
         spl(j,2)=(xnew(j)-xold(indx(j,1)))
     $           /(xold(indx(j,2))-xold(indx(j,1)))
        ELSE
!!      do not treat the points above max in a special way
         spl(j,1)=(xold(indx(j,2))-xnew(j))
     $           /(xold(indx(j,2))-xold(indx(j,1)))
         spl(j,2)=(xnew(j)-xold(indx(j,1)))
     $           /(xold(indx(j,2))-xold(indx(j,1)))
!!     set points above max to zero 
!!         spl(j,1:2)=0.0_spreal
        END IF       
       END DO                   ! j      
       
      end subroutine linspl_sp

C     linear interpolation
C     assuming that the function is zero for xnew > xold(n) 
C     and assuming extrapolation from interval xold(1),xold(2) 
C     to xnew > xold(1) 
C     this subroutine assumes that the xold points are sorted 
C     according to their seize
C     double precision version 

      subroutine linspl_dp(xold,n,xnew,m,spl,indx)
       IMPLICIT NONE
       INTEGER n,m
       REAL(dpreal) xold(n),xnew(m)
       REAL(dpreal) spl(m,2)
       INTEGER indx(m,2)
       INTEGER i,j

C     we need at least two points
       IF(n.LT.2) STOP 'linear interpolation impossible'
C     first find correct indices for each xnew 
       DO j=1,m
        DO i=1,n
         IF(xold(i).GT.xnew(j)) exit 
        END DO                  ! i
        indx(j,2)=min(max(i,2),n) ! second point at least 2 but .le. n
        indx(j,1)=indx(j,2)-1   ! one less
       END DO                   ! j

C     fix the spl weights 
       
       DO j=1,m
        IF(xnew(j).LE.xold(n)*1.01) THEN
         spl(j,1)=(xold(indx(j,2))-xnew(j))
     $           /(xold(indx(j,2))-xold(indx(j,1)))
         spl(j,2)=(xnew(j)-xold(indx(j,1)))
     $           /(xold(indx(j,2))-xold(indx(j,1)))
        ELSE
         spl(j,1:2)=0.0_dpreal
        END IF       
       END DO                   ! j      
       
      end subroutine linspl_dp

      END MODULE spline

