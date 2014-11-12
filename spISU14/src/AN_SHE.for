************************************************************************
      SUBROUTINE PowerLaw(E1,E2,F1,F2,gamm,F00)
************************************************************************
*                                                                      *
*     SuperHigh-Energy atmospheric neutrino fluxes.                    *
*                                                                      *
*            Energy range is 10^8 GeV to E_cut (5*10^10 GeV).          *
*                                                                      *
*     SUBROUTINE PowerLaw initials the parameters for extrapolation    *
*     of high-energy fluxes with the function F00*Energy^(-gamm).      *
*                                                                      *
*     FUNCTION AN_SHE_cut returns values of AN fluxes with employment  *
*     GZK cutoff.                                                      *
*                                                                      *
************************************************************************


         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

               REAL,PARAMETER::
     #              gamm_min=3
               REAL
     #              gamm,F00,
     #              E1,F1,E2,F2

         IF (E1.EQ.E2) STOP 'Mistake in SUBROUTINE PowerLaw'
         gamm=log(F1/F2)/log(E2/E1)
         IF (gamm.LE.gamm_min)
     #   WRITE(*,*) 'Warning from SUBROUTINE PowerLaw: a rum gamm!'
         F00=F1*E1**gamm
         RETURN

      END SUBROUTINE PowerLaw

************************************************************************
      FUNCTION AN_SHE_cut(E)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: pi,E_GZK

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

               REAL,PARAMETER::
     #              E_SHE_start=1d+08,
     #              cut=0.5*pi/(E_GZK-E_SHE_start)
               REAL
     #              AN_SHE_cut,E

         AN_SHE_cut=1.0/(1.0+tan(cut*(E-E_SHE_start)))
         IF (AN_SHE_cut<0) AN_SHE_cut=0
         RETURN

      END FUNCTION AN_SHE_cut
