************************************************************************
*                                                                      *
*     The four routines below are slightly  modified ones from the     *
*     "NUMERICAL RECIPIES" (FORTRAN 77 Version 2.07, CDROM release     *
*     of Microsoft FORTRAN PowerStation).                              *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] W.H. Press, S.A. Teukolsky,  W.T. Vetterling, B.P. Flan-     *
*         nery, "Numerical recipes in FORTRAN. The art of scienti-     *
*         fic computing,"  Cambridge University Press,  ISBN 0 521     *
*         43064-X.                                                     *
*                                                                      *
************************************************************************
      SUBROUTINE Spline_mod(x,y,n,yp1,ypn,y2)
************************************************************************
         
         IMPLICIT REAL (A-H,O-Z)

         INTEGER,PARAMETER::
     #           Nmax   = 500
            REAL,PARAMETER::
     #           VeryBig= 0.99d+30

         DIMENSION x(*),y(*),y2(*),u(Nmax)

         IF (yp1.GT.VeryBig) THEN
           y2(1)=0.0
           u(1)=0.0
                             ELSE
           y2(1)=-0.5
           Dx=x(2)-x(1)
           u(1)=(3.0/Dx)*((y(2)-y(1))/Dx-yp1)
      endIF
         DO 1 i=2,n-1
         Dx1=x(i)-x(i-1)
         Dx2=x(i+1)-x(i-1)
         Dx3=x(i+1)-x(i)
         Sig=Dx1/Dx2
         p=Sig*y2(i-1)+2.0
         y2(i)=(Sig-1.0)/p
    1    u(i)=(6.0*((y(i+1)-y(i))/Dx3-
     #        (y(i)-y(i-1))/Dx1)/Dx2-Sig*u(i-1))/p
         IF (ypn.GT.VeryBig) THEN
           y2(n)=0.0
                             ELSE
           Dx=x(n)-x(n-1)
           y2(n)=((6.0/Dx)*(ypn-(y(n)-y(n-1))/Dx)-u(n-1))/(y2(n-1)+2.0)
      endIF
         DO 2 k=n-1,1,-1
    2    y2(k)=y2(k)*y2(k+1)+u(k)

         RETURN

*     ==================================================================
      ENTRY SplinN_mod(x,y,n,y2)                                         "Natural" spline
*     ==================================================================
         y2(1)=0.0
         u(1)=0.0
         DO 3 i=2,n-1
         Dx1=x(i)-x(i-1)
         Dx2=x(i+1)-x(i-1)
         Dx3=x(i+1)-x(i)
         Sig=Dx1/Dx2
         p=Sig*y2(i-1)+2.0
         y2(i)=(Sig-1.0)/p
    3    u(i)=(6.0*((y(i+1)-y(i))/Dx3-(y(i)-y(i-1))/Dx1)/Dx2-
     #        Sig*u(i-1))/p
         y2(n)=0.0
         DO 4 k=n-1,1,-1
    4    y2(k)=y2(k)*y2(k+1)+u(k)
         RETURN

      END SUBROUTINE Spline_mod

************************************************************************
      FUNCTION Splint_mod(xa,ya,y2a,n,x)                                 Modified SPLINT
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         DIMENSION xa(*),y2a(*),ya(*)
         klo=1
         khi=n
    1    IF (khi-klo.GT.1) THEN
           k=(khi+klo)*0.5
           IF (xa(k).GT.x) THEN
             khi=k
                           ELSE
             klo=k
        endIF
           GOTO 1
      endIF
         h=xa(khi)-xa(klo)
         IF (h.EQ.0.0) PAUSE 'Bad xa input in Splint_mod'
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         Splint_mod=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)
     #                       +(b**3-b)*y2a(khi))*(h**2)/6

         RETURN
      END FUNCTION Splint_mod

************************************************************************
      SUBROUTINE Splie2_mod(x1a,x2a,ya,m,n,y2a,test)                     Modified SPLIE2
************************************************************************

         IMPLICIT REAL (A-H,O-Z)
         LOGICAL test

         INTEGER,PARAMETER::
     #           NN  = 500

         DIMENSION x1a(*),x2a(*),y2a(m,n),ya(m,n),y2tmp(NN),ytmp(NN)

         DO j=1,m
           DO 1 k=1,n
    1      ytmp(k)=ya(j,k)
           CALL SplinN_mod(x2a,ytmp,n,y2tmp)
           DO 2 k=1,n
    2      y2a(j,k)=y2tmp(k)
      endDO
         IF (.NOT.test) RETURN
*        Test calculation of the maximum discrepancy in reference points
         Discr_max=0.0
         i=0
         DO 3 j=1,m
         x1=x1a(j)
         DO 3 k=1,n
         y=ya(j,k)
         IF (y.EQ.0.0) THEN
           i=i+1
           Discr=        Splin2_mod(x1a,x2a,ya,y2a,m,n,x1,x2a(k))
                        ELSE
           Discr=ABS(1.0-Splin2_mod(x1a,x2a,ya,y2a,m,n,x1,x2a(k))/y)
      endIF
    3    IF (Discr.GT.Discr_max) Discr_max=Discr
         PRINT 4, Discr_max
c        IF (i.GT.0) PRINT *,' y = 0 in',i,' reference points'

         RETURN
    4 FORMAT(2x,'Maximum discrepancy in reference points =',1PD10.3)
      END SUBROUTINE Splie2_mod

************************************************************************
      FUNCTION Splin2_mod(x1a,x2a,ya,y2a,m,n,x1,x2)                      Modified SPLIN2
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         INTEGER,PARAMETER::
     #           NN= 500

         DIMENSION x1a(*),x2a(*),ya(m,n),y2a(m,n)
         DIMENSION y2tmp(NN),ytmp(NN),ztmp(NN)

         DO 2 j=1,m
         DO 1 k=1,n
         ytmp(k)=ya(j,k)
    1    y2tmp(k)=y2a(j,k)
    2    ztmp(j) =Splint_mod(x2a,ytmp,y2tmp,n,x2)

         CALL SplinN_mod(x1a,ztmp,m,y2tmp)
         Splin2_mod=Splint_mod(x1a,ztmp,y2tmp,m,x1)

         RETURN
      END FUNCTION Splin2_mod

************************************************************************
*                                                                      *
*     The four routines below  are a  particular case of the above     *
*     ones simplified for equidistant grides.                          *
*                                                                      *
************************************************************************
      SUBROUTINE Spline_ED(y,n,s,yp1,ypn,y2)
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         INTEGER,PARAMETER::
     #           Nmax   =500

         DIMENSION y(*),y2(*),u(Nmax)

         IF (yp1.GT.VeryBig) THEN
           y2(1)=0.0
           u(1)=0.0
                             ELSE
           y2(1)=-0.5
           u(1)=(3.0/s)*((y(2)-y(1))/s-yp1)
      endIF
         a=6.0/s**2
         DO 1 i=2,n-1
         Q=y2(i-1)+4.0
         y2(i)=-1.0/Q
    1    u(i)=(a*(y(i+1)-2.0*y(i)+y(i-1))-u(i-1))/Q
         IF (ypn.GT.VeryBig) THEN
           y2(n)=0.0
                             ELSE
           y2(n)=(a*(ypn*s-y(n)+y(n-1))-u(n-1))/(y2(n-1)+2.0)
      endIF
         DO 2 k=n-1,1,-1
    2    y2(k)=y2(k)*y2(k+1)+u(k)

         RETURN

*     ==================================================================
      ENTRY SplinN_ED(y,n,s,y2)                                          "Natural" spline
*     ==================================================================
         y2(1)=0.0
         u(1)=0.0
         a=6.0/s**2
         DO 3 i=2,n-1
         Q=y2(i-1)+4.0
         y2(i)=-1.0/Q
    3    u(i)=(a*(y(i+1)-2.0*y(i)+y(i-1))-u(i-1))/Q
         y2(n)=0.0
         DO 4 k=n-1,1,-1
    4    y2(k)=y2(k)*y2(k+1)+u(k)

         RETURN

      END SUBROUTINE Spline_ED

************************************************************************
      FUNCTION Splint_ED(delta,s,ya,y2a,n)                               KK: This variable has not been used [n]
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         SAVE

         REAL,PARAMETER::
     #        OS = 1.0/6.0

         DIMENSION y2a(*),ya(*)

         d  =delta/s
         klo=INT(d)+1.0
         khi=klo+1
         d  =d+1.0
         a  =DFLOAT(khi)-d
         b  =d-DFLOAT(klo)
         d  =s**2*OS
         fa =a*(a**2*d-d)
         fb =b*(b**2*d-d)
         Splint_ED=a*ya(klo)+b*ya(khi)+fa*y2a(klo)+fb*y2a(khi)
         RETURN

*     ==================================================================
      ENTRY SplRe_ED(ya,y2a)
*     ==================================================================
         SplRe_ED=a*ya(klo)+b*ya(khi)+fa*y2a(klo)+fb*y2a(khi)
         RETURN

      END FUNCTION Splint_ED

************************************************************************
      SUBROUTINE Splie2_ED(x1_0,x2_0,s1,s2,ya,m,n,y2a,test)
************************************************************************

         IMPLICIT REAL (A-H,O-Z)
         LOGICAL(2) test

         INTEGER,PARAMETER::
     #           NN  = 500

         DIMENSION y2a(m,n),ya(m,n),y2tmp(NN),ytmp(NN)

         DO j=1,m
           DO 1 k=1,n
    1      ytmp(k)=ya(j,k)
           CALL SplinN_ED(ytmp,n,s2,y2tmp)
           DO 2 k=1,n
    2      y2a(j,k)=y2tmp(k)
      endDO
         IF (.NOT.test) RETURN
*        Test calculation of the maximum discrepancy in reference points
         Discr_max=0.0
         i=0
         x10=x1_0-s1
         x20=x2_0-s2
         DO 3 j=1,m
         x1=x10+s1*j
         DO 3 k=1,n
         x2=x20+s2*k
         y=ya(j,k)
         IF (y.EQ.0.0) THEN
           i=i+1
           Discr=        Splin2_ED(x1_0,x2_0,s1,s2,ya,y2a,m,n,x1,x2)
                        ELSE
           Discr=ABS(1.0-Splin2_ED(x1_0,x2_0,s1,s2,ya,y2a,m,n,x1,x2)/y)
      endIF
    3    IF (Discr.GT.Discr_max) Discr_max=Discr
         PRINT 4, Discr_max
         IF (i.GT.0) PRINT *,' y = 0 in',i,' reference points'

         RETURN
    4 FORMAT(2x,'Maximum discrepancy in reference points =',1PD10.3)
      END SUBROUTINE Splie2_ED

************************************************************************
      FUNCTION Splin2_ED(x1_0,x2_0,s1,s2,ya,y2a,m,n,x1,x2)
************************************************************************

         IMPLICIT REAL (A-H,O-Z)

         INTEGER,PARAMETER::
     #           NN= 500

         DIMENSION y2a(m,n),ya(m,n),y2tmp(NN),ytmp(NN),ztmp(NN)

         DO 2 j=1,m
         DO 1 k=1,n
         ytmp(k)=ya(j,k)
    1    y2tmp(k)=y2a(j,k)
    2    ztmp(j)=Splint_ED(x2-x2_0,s2,ytmp,y2tmp,n)
         CALL SplinN_ED(ztmp,m,s1,y2tmp)
         Splin2_ED=Splint_ED(x1-x1_0,s1,ztmp,y2tmp,m)

         RETURN
      END FUNCTION Splin2_ED
