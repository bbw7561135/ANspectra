************************************************************************
      SUBROUTINE AN_HE_ISU14
************************************************************************
*                                                                      *
*                         HE_ISU-2014 FLUX                             *
*                      in 1/(cm^2 sec sr GeV)                          *
*                Energy range is 10 GeV to 100 PeV.                    *
*                                                                      *
*     Primary cosmic ray spectra: Hilas-Gaisser.                       * 
*     Hardron interaction model: Kimel-Mokhov.                         *
*                                                                      *
*     Fluxes of conventional high-energy atmospheric neutrinos and     *
*     antineutrinos (total) as a function of energy and zenith angle.  *
*                                                                      *
*     The source data array F is dF/dE in units of 1/(cm^2 s sr GeV).  *
*     The energy reference points are calculated on the equidistant    *
*     over log10(E) grid and the reference points for cos(theta)       *
*     (where theta is the zenith angle) are calculated on equidistant  *
*     grid from 0 to 1 with step = 0.1.                                *
*                                                                      *
************************************************************************

         USE InpOutUnits
         USE PhysMathConstants, ONLY: E_GZK

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
                  
         SAVE

           LOGICAL(2),PARAMETER::
     #                Test   =.FALSE.,                                    Test of spline in reference points
     #                Spectra=.FALSE.,                                    Test of spline for energy spectra
     #                ZenDist=.FALSE.,                                    Test of spline for Z-A distributions
     #                Ratio  =.FALSE.                                     Test of spline for nu/antinu ratios
              INTEGER,PARAMETER::
     #                NE_test=  61,
     #                NC_test=  11,
     #                IncrE  =  10,
     #                K      =   3,
     #                NE     =  71,
!     #                NE_SHE =  26,
     #                NC     =  11,
     #                NFlavor    =   2,
     #                NNuAnu     =   2
              INTEGER
     #                Ndat(NFlavor,NNuAnu)/201,202,203,204/,
     #                iFlavor,iNuAnu
                 REAL,PARAMETER::
     #                E_min  = 1.0d+01,
     #                E_max  = 1.0d+08,
     #                C_min  = 0.0,
     #                C_max  = 1.0
                 REAL
!     #                gamm,F00,
     #                C(NC),E(NE),
     #                T(NC+1,NE+1),F(NFlavor,NNuAnu,NC,NE),
     #                lnF(NFlavor,NNuAnu,NC,NE),
     #                ClnF(NFlavor,NNuAnu,NC,NE),
!     #                E_SHE(NE_SHE),logE_SHE(NE_SHE),
!     #                F_SHE(NFlavor,NNuAnu,NC,NE_SHE),
!     #                lnF_SHE(NFlavor,NNuAnu,NC,NE_SHE),
!     #                ClnF_SHE(NFlavor,NNuAnu,NC,NE_SHE),
     #                SF(NFlavor,NNuAnu,NE_test),
     #                SF0(NFlavor,NNuAnu,NE_test),
     #                Rs(NFlavor,NNuAnu,NC_test),Factor
         CHARACTER(*),PARAMETER::
     #                dir='HEISU/',
     #                ext='.dat'
         CHARACTER*1
     #                fln(NFlavor)/'e','m'/,
     #                NTn(NNuAnu)/'n','a'/

         OPEN(Ndat00,FILE=datACN//'AN_HE_ISU14.data')
         DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
           READ(Ndat00,*) T
           DO n_NE=1,NE
             DO n_NC=1,NC
               F(Flavor,NuAnu,n_NC,n_NE)=T(NC+2-n_NC,n_NE+1)
          endDO
        endDO
      endDO
      endDO
         CLOSE(Ndat00)
         DO n_NC=1,NC
           C(n_NC)=T(NC+2-n_NC,1)
      endDO
         DO n_NE=1,NE
           E(n_NE)=T(1,n_NE+1)
      endDO
         PRINT *,' File AN_HE_ISU14.data was read '

         lgEmin  =log10(E_min)
         lgEmax  =log10(E_max)
         steplgE =(lgEmax-lgEmin)/(NE-1)
         steplgEt=(lgEmax-lgEmin)/(NE_test-1)
         stepC   =(C_max-C_min)/(NC-1)
         stepCt  =(C_max-C_min)/(NC_test-1)

         lnF=log(F)
         DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
        CALL Splie2_ED(0.0,lgEmin,stepC,steplgE,lnF(Flavor,NuAnu,:,:),
     #                NC,NE,ClnF(Flavor,NuAnu,:,:),Test)
      endDO
      endDO

!         steplgE_SHE=(log10(E_GZK)-lgEmax)/(NE_SHE-1)
!         DO n_NE=1,NE_SHE
!           logE_SHE(n_NE)=lgEmax+steplgE_SHE*(n_NE-1)
!           E_SHE(n_NE)=10.0**logE_SHE(n_NE)
!      endDO
!         DO n_NC=1,NC
!           DO Flavor=1,NFlavor
!         DO NuAnu=1,NNuAnu
!             CALL PowerLaw(E(NE-1),E(NE),F(Flavor,NuAnu,n_NC,NE-1),
!     #                             F(Flavor,NuAnu,n_NC,NE),gamm,F00)
!             DO n_NE=1,NE_SHE
!               F_SHE(Flavor,NuAnu,n_NC,n_NE)=AN_SHE_cut(E_SHE(n_NE))
!     #                *F00*E_SHE(n_NE)**(-gamm)
!          endDO
!        endDO
!        endDO
!      endDO
!         lnF_SHE=log(F_SHE)

*        ============================================================= *
*        TEST OF ENERGY SPECTRA (DATA FLUX)                            *
*        ============================================================= *
         IF (Spectra) THEN                                               Test for energy spectra
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
            OPEN(Ndat(Flavor,NuAnu),FILE=outACN//dir//'FE_'//
     #                fln(Flavor)//Ntn(NuAnu)//ext)
        endDO
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             Factor=Energy**K
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
               WRITE(Ndat(Flavor,NuAnu),1) Energy,(Factor*F(Flavor,
     #                NuAnu,n_NC,n_NE),n_NC=1,NC)
          endDO
          endDO
        endDO
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
             CLOSE(Ndat(Flavor,NuAnu))
            OPEN(Ndat(Flavor,NuAnu),FILE=outACN//dir//'SE_'//
     #                fln(Flavor)//Ntn(NuAnu)//ext)
        endDO
        endDO
           DO n_NE=1,NE_test
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=10.0**lgEt
             Factor=Energy**K
             DO n_NC=1,NC_test
               Cosine=C_min+(n_NC-1)*stepCt
               DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
                 SF(Flavor,NuAnu,n_NC)=exp(Splin2_ED(0.0,lgEmin,stepC,
     #             steplgE,lnF(Flavor,NuAnu,:,:),ClnF(Flavor,NuAnu,:,:),
     #                NC,NE,Cosine,lgEt))
            endDO
            endDO
          endDO
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
               WRITE(Ndat(Flavor,NuAnu),3) Energy,(Factor*SF(Flavor,
     #                NuAnu,n_NC),n_NC=1,NC_test)
          endDO
          endDO
        endDO
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
             CLOSE(Ndat(Flavor,NuAnu))
        endDO
        endDO
      endIF
*        ============================================================= *
*        TEST FOR ZENITH-ANGLE DISTRIBUTIONS                           *
*        ============================================================= *
         IF (ZenDist) THEN                                               Test for zenith-angle distributions
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
            OPEN(Ndat(Flavor,NuAnu),FILE=outACN//dir//'FZ_'//
     #                fln(Flavor)//Ntn(NuAnu)//ext)
        endDO
        endDO
           DO n_NC=1,NC
             Cosine=C(n_NC)
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
               WRITE(Ndat(Flavor,NuAnu),2) Cosine,(F(Flavor,NuAnu,n_NC,
     #                n_NE)/F(Flavor,NuAnu,NC,n_NE),n_NE=1,NE,IncrE)
          endDO
          endDO
        endDO
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
             CLOSE(Ndat(Flavor,NuAnu))
            OPEN(Ndat(Flavor,NuAnu),FILE=outACN//dir//'SZ_'//
     #                fln(Flavor)//Ntn(NuAnu)//ext)
        endDO
        endDO
           DO n_NE=1,NE_test,IncrE                                       Calculating the vertical fluxes
             lgEt=lgEmin+(n_NE-1)*steplgEt
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
               SF0(Flavor,NuAnu,n_NE)=exp(Splin2_ED(0.0,lgEmin,stepC,
     #       steplgE,lnF(Flavor,NuAnu,:,:),ClnF(Flavor,NuAnu,:,:),NC,NE,
     #                1.0,lgEt))
          endDO
          endDO
        endDO
           DO n_NC=1,NC_test
             Cosine=C_min+(n_NC-1)*stepCt
             DO n_NE=1,NE_test,IncrE
               lgEt=lgEmin+(n_NE-1)*steplgEt
               DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
                 SF(Flavor,NuAnu,n_NE)=exp(Splin2_ED(0.0,lgEmin,stepC,
     #        steplgE,lnF(Flavor,NuAnu,:,:),ClnF(Flavor,NuAnu,:,:),NC,
     #                NE,Cosine,lgEt))
            endDO
            endDO
          endDO
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
        WRITE(Ndat(Flavor,NuAnu),4) Cosine,(SF(Flavor,NuAnu,n_NE)/
     #                SF0(Flavor,NuAnu,n_NE),n_NE=1,NE,IncrE)
          endDO
          endDO
        endDO
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
             CLOSE(Ndat(Flavor,NuAnu))
        endDO
        endDO
      endIF

*        ============================================================= *
*        TEST FOR NEUTRINO/ANTINEUTRINO RATIO                          *
*        ============================================================= *
         IF (Ratio) THEN                                                 Test for neutrino/antineutrino ratio
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
            OPEN(Ndat(Flavor,NuAnu),FILE=outACN//dir//'FR_'//
     #                fln(Flavor)//Ntn(NuAnu)//ext)
        endDO
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
               WRITE(Ndat(Flavor,NuAnu),1) Energy,(F(Flavor,NuAnu,n_NC,
     #                n_NE)/F(Flavor,NuAnu,n_NC,n_NE),n_NC=1,NC)
          endDO
          endDO
        endDO
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
             CLOSE(Ndat(Flavor,NuAnu))
            OPEN(Ndat(Flavor,NuAnu),FILE=outACN//dir//'SR_'//
     #                fln(Flavor)//Ntn(NuAnu)//ext)
        endDO
        endDO
           DO n_NE=1,NE_test
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=10.0**lgEt
             DO n_NC=1,NC_test
               Cosine=C_min+(n_NC-1)*stepCt
               DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
                 Rs(Flavor,NuAnu,n_NC)=exp(Splin2_ED(0.0,lgEmin,stepC,
     #             steplgE,lnF(Flavor,NuAnu,:,:),ClnF(Flavor,NuAnu,:,:),
     #            NC,NE,Cosine,lgEt)-Splin2_ED(0.0,lgEmin,stepC,steplgE,
     #  lnF(Flavor,NuAnu,:,:),ClnF(Flavor,NuAnu,:,:),NC,NE,Cosine,lgEt))
            endDO
            endDO
          endDO
             DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
            WRITE(Ndat(Flavor,NuAnu),3) Energy,(Rs(Flavor,NuAnu,n_NC),
     #                n_NC=1,NC)
          endDO
          endDO
        endDO
           DO Flavor=1,NFlavor
         DO NuAnu=1,NNuAnu
             CLOSE(Ndat(Flavor,NuAnu))
        endDO
        endDO
      endIF
*        ============================================================= *

         RETURN

    1 FORMAT(1PE9.3,11(1PE12.4))
    2 FORMAT(1PE9.3,7(1PE12.4))
    3 FORMAT(1PE9.3,11(1PE12.4))
    4 FORMAT(1PE9.3,7(1PE12.4))


*     ==================================================================
      ENTRY AN_HE_ISU(iFlavor,iNuAnu,F0,E0,C0)
*     ==================================================================
         F0=exp(Splin2_ED(0.0,lgEmin,stepC,steplgE,
     #   lnF(iFlavor,iNuAnu,:,:),ClnF(iFlavor,iNuAnu,:,:),NC,NE,abs(C0),
     #   log10(E0)))
         RETURN

************************************************************************
*                                                                      *
*     Entry for extrapolation of atmospheric neutrino fluxes with the  *
*     power function F00*Energy^(-gamm) and employment GZK cutoff      *
*     (AN_SHE.for)                                                     *
*                                                                      *
************************************************************************
*     ==================================================================
      ENTRY AN_SHE_ISU(iFlavor,iNuAnu,F0,E0,C0)
*     ==================================================================
!         F0=exp(Splin2_ED(0.0,lgEmax,stepC,steplgE_SHE,
!     #   lnF_SHE(iFlavor,iNuAnu,:,:),ClnF_SHE(iFlavor,iNuAnu,:,:),NC,
!     #   NE_SHE,abs(C0),log10(E0)))
          F0=0
         RETURN

  100 FORMAT((1PE10.8)$)

      END SUBROUTINE AN_HE_ISU14
