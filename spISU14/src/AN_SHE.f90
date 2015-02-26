!**********************************************************************!
SUBROUTINE PowerLaw(E1,E2,F1,F2,gamm,F00)
!----------------------------------------------------------------------!
!SUBROUTINE PowerLaw initials the parameters for extrapolation of high-
!energy fluxes with the function F00*Energy^(-gamm).
!**********************************************************************!
implicit none

real,parameter:: &
    gammin=3
real &
    gamm,F00,E1,F1,E2,F2

    if(E1.eq.E2)then
        write(*,*) 'PowerLaw ERROR: E1=',E1,', E2=',E2
        stop
    endif
    gamm=log(F1/F2)/log(E2/E1)
    if(gamm<=gammin)then
        write(*,*) 'PowerLaw warning: gamm=',gamm,' <=gammin!'
    endif
    F00=F1*E1**gamm
    return
!----------------------------------------------------------------------!
      END SUBROUTINE PowerLaw
!**********************************************************************!

!**********************************************************************!
      FUNCTION AN_SHE_cut(E)
!----------------------------------------------------------------------!
!SuperHigh-Energy (from 10^8 GeV to E_cut) atmospheric neutrino fluxes.
!----------------------------------------------------------------------!
use PhysMathConstants, only: pi,E_GZK

implicit none

real,parameter:: &
    E_SHE_start=1d+08,&
    cut=0.5*pi/(E_GZK-E_SHE_start)
real &
    AN_SHE_cut,E

    AN_SHE_cut=1.0/(1.0+tan(cut*(E-E_SHE_start)))
    if(AN_SHE_cut<0)then
        AN_SHE_cut=0
    endif
    return
!----------------------------------------------------------------------!
      END FUNCTION AN_SHE_cut
!**********************************************************************!
