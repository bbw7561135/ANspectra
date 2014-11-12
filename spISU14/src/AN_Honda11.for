************************************************************************
      SUBROUTINE AN_Honda11
************************************************************************
*                                                                      *
*                         HONDA-2011 FLUX                              *
*                      in 1/(cm^2 sec sr GeV)                          *
*                Energy range is 0.1 GeV to 10 TeV.                    *
*                                                                      *
*     Solar activity: maximal.                                         *
*     Site: Kamioka.                                                   *
*                                                                      *
*     Fluxes of atmospheric neutrinos as a function of energy and      *
*     zenith angle (without mountain over the detector).               *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [ 1] M. Honda,  T. Kajita,  K. Kasahara, and  S. Midorikawa,     *
*          "Improvement  of low energy  atmospheric neutrino  flux     *
*          calculation  using the JAM nuclear  interaction model,"     *
*          Phys. Rev. D83 (2011) 123001 [arXiv:astro-ph/1102.2688].    *
*          [ http://icrr.u-tokyo.ac.jp/~mhonda/nflx2011 ]              *
*                                                                      *
*     The source data array F is dF/dE in units of 1/(m^2 s sr GeV).   *
*     Zenith angle bins cos = -1.0-(-0.9), ..., 0.0-0.1, 0.9-1.0.      *
*                                                                      *
************************************************************************

         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)
                  
         SAVE

           LOGICAL(2),PARAMETER::
     #                Test   =.FALSE.,                                    Test of spline in reference points
     #                Spectra=.FALSE.,                                    Test of spline for energy spectra
     #                ZenDist=.FALSE.,                                    Test of spline for Z-A distributions
     #                Ratio  =.FALSE.                                     Test of spline for nu/antinu ratios
              INTEGER,PARAMETER::
     #                NE_test= 101,
     #                NC_test=  20,
     #                IncrE  =  10,
     #                K      =   3,
     #                NE     = 101,
     #                NC     =  20,
     #                Nfl    =   2,
     #                Nt     =   2
              INTEGER
     #                Ndat(Nfl,Nt)/301,302,303,304/
                 REAL,PARAMETER::
     #                E_min  = 1.0d-01,
     #                E_max  = 1.0d+04,
     #                C_min  =-9.5d-01,
     #                C_max  = 9.5d-01,
     #                Rescale= 1.0d-04
                 REAL
     #                C(NC),E(NE),
     #                F(Nfl,Nt,NC,NE),CF(Nfl,Nt,NC,NE),
     #                SF(Nfl,Nt,NE_test),SF0(Nfl,Nt,NE_test),
     #                Rs(Nfl,NC_test),Factor
         CHARACTER(*),PARAMETER::
     #                dir='Honda/',
     #                ext='.dat'
         CHARACTER*1
     #                fln(Nfl)/'e','m'/,
     #                NTn(Nt)/'n','a'/

         OPEN(Ndat00,FILE=datACN//'AN_Honda11.data')
         DO n_NC=NC,1,-1
           READ(Ndat00,*) C(n_NC)
           DO n_NE=1,NE
             READ(Ndat00,*) E(n_NE),((F(n_fl,nNTn,n_NC,n_NE),nNTn=1,Nt),
     #                                                    n_fl=Nfl,1,-1)
        endDO
      endDO
         CLOSE(Ndat00)
         PRINT *,' File AN_Honda11.data was read '
         
         F=F*Rescale

         lgEmin  =log10(E_min)
         lgEmax  =log10(E_max)
         steplgE =(lgEmax-lgEmin)/(NE-1)
         steplgEt=(lgEmax-lgEmin)/(NE_test-1)
         stepC   =(C_max-C_min)/(NC-1)
         stepCt  =(C_max-C_min)/(NC_test-1)

         DO n_fl=1,Nfl
           DO nNTn=1,Nt
             CALL Splie2_mod(C,E,F(n_fl,nNTn,:,:),NC,NE,
     #                                           CF(n_fl,nNTn,:,:),Test)
        endDO
      endDO

*        ============================================================= *
*        TEST OF ENERGY SPECTRA (DATA FLUX)                            *
*        ============================================================= *
         IF (Spectra) THEN
           DO n_fl=1,Nfl
             DO nNTn=1,Nt
               OPEN(Ndat(n_fl,nNTn),FILE=outACN//dir//'FE_'//fln(n_fl)
     #                                                 //NTn(nNTn)//ext)
          endDO
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             Factor=Energy**K
             DO n_fl=1,Nfl
               DO nNTn=1,Nt
                 WRITE(Ndat(n_fl,nNTn),1) Energy,(Factor*
     #                                 F(n_fl,nNTn,n_NC,n_NE),n_NC=1,NC)
            endDO
          endDO
        endDO
           DO n_fl=1,Nfl
             DO nNTn=1,Nt
               CLOSE(Ndat(n_fl,nNTn))
               OPEN(Ndat(n_fl,nNTn),FILE=outACN//dir//'SE_'//fln(n_fl)
     #                                                 //NTn(nNTn)//ext)
          endDO
        endDO
           DO n_NE=1,NE_test
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=10.0**lgEt
             Factor=Energy**K
             DO n_NC=1,NC_test
               Cosine=C_min+(n_NC-1)*stepCt
               DO n_fl=1,Nfl
                 DO nNTn=1,Nt
                   SF(n_fl,nNTn,n_NC)=Splin2_mod(C,E,F(n_fl,nNTn,:,:),
     #                            CF(n_fl,nNTn,:,:),NC,NE,Cosine,Energy)
              endDO
            endDO
          endDO
             DO n_fl=1,Nfl
               DO nNTn=1,Nt
                 WRITE(Ndat(n_fl,nNTn),3) Energy,(Factor*
     #                                SF(n_fl,nNTn,n_NC),n_NC=1,NC_test)
            endDO
          endDO
        endDO
           DO n_fl=1,Nfl
             DO nNTn=1,Nt
               CLOSE(Ndat(n_fl,nNTn))
          endDO
        endDO
      endIF

*        ============================================================= *
*        TEST FOR ZENITH-ANGLE DISTRIBUTIONS                           *
*        ============================================================= *
         IF (ZenDist) THEN                                               Test for zenith-angle distributions
           DO n_fl=1,Nfl
             DO nNTn=1,Nt
               OPEN(Ndat(n_fl,nNTn),FILE=outACN//dir//'FZ_'//fln(n_fl)
     #                                                 //NTn(nNTn)//ext)
          endDO
        endDO
           DO n_NC=1,NC
             Cosine=C(n_NC)
             DO n_fl=1,Nfl
               DO nNTn=1,Nt
                 WRITE(Ndat(n_fl,nNTn),2) Cosine,(F(n_fl,nNTn,n_NC,n_NE)
     #                            /F(n_fl,nNTn,NC,n_NE),n_NE=1,NE,IncrE)
            endDO
          endDO
        endDO
           DO n_fl=1,Nfl
             DO nNTn=1,Nt
               CLOSE(Ndat(n_fl,nNTn))
               OPEN(Ndat(n_fl,nNTn),FILE=outACN//dir//'SZ_'//fln(n_fl)
     #                                                 //NTn(nNTn)//ext)
          endDO
        endDO
           DO n_NE=1,NE_test,IncrE                                       Calculating the vertical fluxes
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=10.0**lgEt
             DO n_fl=1,Nfl
               DO nNTn=1,Nt
                 SF0(n_fl,nNTn,n_NE)=Splin2_mod(C,E,F(n_fl,nNTn,:,:),
     #                               CF(n_fl,nNTn,:,:),NC,NE,1.0,Energy)
            endDO
          endDO
        endDO
           DO n_NC=1,NC_test
             Cosine=C_min+(n_NC-1)*stepCt
             DO n_NE=1,NE_test,IncrE
               lgEt=lgEmin+(n_NE-1)*steplgEt
               Energy=10.0**lgEt
               DO n_fl=1,Nfl
                 DO nNTn=1,Nt
                   SF(n_fl,nNTn,n_NE)=Splin2_mod(C,E,F(n_fl,nNTn,:,:),
     #                            CF(n_fl,nNTn,:,:),NC,NE,Cosine,Energy)
              endDO
            endDO
          endDO
             DO n_fl=1,Nfl
               DO nNTn=1,Nt
                 WRITE(Ndat(n_fl,nNTn),4) Cosine,(SF(n_fl,nNTn,n_NE)/
     #                         SF0(n_fl,nNTn,n_NE),n_NE=1,NE_test,IncrE)
            endDO
          endDO
        endDO
           DO n_fl=1,Nfl
             DO nNTn=1,Nt
               CLOSE(Ndat(n_fl,nNTn))
          endDO
        endDO
      endIF

*        ============================================================= *
*        TEST FOR NEUTRINO/ANTINEUTRINO RATIO                          *
*        ============================================================= *
         IF (Ratio) THEN                                                 Test for neutrino/antineutrino ratio
           DO n_fl=1,Nfl
             OPEN(Ndat(n_fl,1),FILE=outACN//dir//'FR_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE
             Energy=E(n_NE)
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl,1),1) Energy,(F(n_fl,1,n_NC,n_NE)/
     #                                   F(n_fl,2,n_NC,n_NE),n_NC=1,NC)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl,1))
             OPEN(Ndat(n_fl,1),FILE=outACN//dir//'SR_'//fln(n_fl)//ext)
        endDO
           DO n_NE=1,NE_test
             lgEt=lgEmin+(n_NE-1)*steplgEt
             Energy=10.0**lgEt
             DO n_NC=1,NC_test
               Cosine=C_min+(n_NC-1)*stepCt
               DO n_fl=1,Nfl
                 Rs(n_fl,n_NC)=Splin2_mod(C,E,F(n_fl,1,:,:),
     #           CF(n_fl,1,:,:),NC,NE,Cosine,Energy)/
     #           Splin2_mod(C,E,F(n_fl,2,:,:),
     #           CF(n_fl,2,:,:),NC,NE,Cosine,Energy)
            endDO
          endDO
             DO n_fl=1,Nfl
               WRITE(Ndat(n_fl,1),3) Energy,(Rs(n_fl,n_NC),
     #                                                   n_NC=1,NC_test)
          endDO
        endDO
           DO n_fl=1,Nfl
             CLOSE(Ndat(n_fl,1))
        endDO
      endIF
*        ============================================================= *

         RETURN

    1 FORMAT(1PE10.3,20(1PE12.4))
    2 FORMAT(1PE10.3,11(1PE12.4))
    3 FORMAT(1PE10.3,20(1PE12.4))
    4 FORMAT(1PE10.3,11(1PE12.4))

*     ==================================================================
      ENTRY AN_Honda(nflavor,ntype,F0,E0,C0)
*     ==================================================================
         F0=Splin2_mod(C,E,F(nflavor,ntype,:,:),CF(nflavor,ntype,:,:),
     #                                                      NC,NE,C0,E0)
         RETURN

      END SUBROUTINE AN_Honda11
