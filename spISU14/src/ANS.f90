!**********************************************************************!
recursive FUNCTION ANS_Init(Mode,iTrans,Dir)
!----------------------------------------------------------------------!
!Atmospheric neutrino fluxes dF/dE [1/(cm^2 s sr GeV)] for any energy  
!and zenith angle.
!
!            Data neutrino energy range is 50 MeV to 1 EeV:
!
!    0     1   10    10^2  10^3  10^4  10^5  10^6  10^7  10^8  10^9
!    |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
!     |----|-----|---|-------|-----|-----|-----|-----|-----------|
!    0.05  1   10   70     10^3  10^4  10^5  10^6  10^7        10^9
!                |------------------ISU14------------------|
!     |---------------------------CORT---------------------------|
!      |----------Honda11----------|
!
!                       zenith angles cosines:
!ISU14                                    0  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
! |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
!   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
!-0.95 -0.85                         -0.05 0.05 0.15                           0.95
!Honda11
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!
use PhysMathConstants, only: E_GZK

implicit none

logical:: &
    ANS_Init,ANS_SmartSpline
real:: &
    ANS_Spline_ME,ANS_Constructor,&
    ANS_Spline_LE,ANS_Spline_HtoL,ANS_Spline_LtoH,ANS_Spline_HE,&
    AN_SHE_cut,Spectrum_Intersection

real,parameter:: &
    Honda11_ini=1.0d-01, Honda11_fin=1.0d+04,&
    HE_ISU_ini =1.0d+01, HE_ISU_fin =1.0d+08,&
    firstep=0.01,&
    lastep=100                                                         !Neutrino energy step from the end of spectrum, 
                                                                       !defined the power tale lean

logical &
    bufL
integer &
    Dir,Mode,NuAnu,Flavor,&
    iTrans,iNuAnu,iFlavor
real &
    F_L_ini,F_L_sec,F_H_pen,F_H_fin,a,b,F0,gamm,&
    E,C

save
integer &
    Trans
real &
    E_LH(2,2),&                                                        !Boundary between low- and high-energy spectrum application
    E_L_ini,E_L_sec,E_H_pen,E_H_fin

    Trans=iTrans
    if(Honda11_fin<HE_ISU_ini)then
        stop 'ANS_Constructor ERROR: Honda11_fin<HE_ISU_ini'
    endif
    E_L_ini=Honda11_ini
    E_L_sec=E_L_ini+firstep
    E_H_fin=HE_ISU_fin
    E_H_pen=E_H_fin-lastep
    write(*,*) 'Atmospheric neutrino model:'
    write(*,*) 'Honda11 (maximal solar activity) + HE_ISU14'
    selectcase(Mode)
        case(1)
            do Flavor=1,2
                do NuAnu=1,2
                    E_LH(Flavor,NuAnu)=Spectrum_Intersection(Flavor,NuAnu)
                enddo
            enddo
            bufL=ANS_SmartSpline(Trans,Dir)
        case(2)
            bufL=ANS_SmartSpline(Trans,Dir)
    endselect

!----------------------------------------------------------------------!
    ANS_Init=.true.
    return

!**********************************************************************!
ENTRY ANS_Constructor(iFlavor,iNuAnu,E,C,Mode)
!----------------------------------------------------------------------!
    if(E<0 .or.E>E_GZK)then
        ANS_Constructor=0
    elseif(E<E_L_ini)then
        F_L_ini=ANS_Spline_LE(iFlavor,iNuAnu,E_L_ini,C)
        F_L_sec=ANS_Spline_LE(iFlavor,iNuAnu,E_L_sec,C)
        a=(F_L_sec-F_L_ini)/(E_L_sec-E_L_ini)
        b=F_L_ini-a*E_L_ini
        ANS_Constructor=a*E+b
    elseif(E<=E_H_fin)then
        selectcase(Mode)
            case(1)
                if(E<=E_LH(iFlavor,iNuAnu))then
                    ANS_Constructor=ANS_Spline_LE(iFlavor,iNuAnu,E,C)
                else
                    ANS_Constructor=ANS_Spline_HE(iFlavor,iNuAnu,E,C)
                endif
            case(2)
                ANS_Constructor=ANS_Spline_ME(iFlavor,iNuAnu,E,C)
        endselect
    else
        F_H_pen=ANS_Spline_HE(iFlavor,iNuAnu,E_H_pen,C)
        F_H_fin=ANS_Spline_HE(iFlavor,iNuAnu,E_H_fin,C)
        call PowerLaw(E_H_pen,E_H_fin,F_H_pen,F_H_fin,gamm,F0)
        ANS_Constructor=F0*E**(-gamm)*AN_SHE_cut(E)
    endif
    return

!**********************************************************************!
ENTRY ANS_Spline_ME(iFlavor,iNuAnu,E,C)
!----------------------------------------------------------------------!
    selectcase(Trans)
        case(1)
            ANS_Spline_ME=ANS_Spline_LtoH(iFlavor,iNuAnu,E,C)
        case(2)
            ANS_Spline_ME=ANS_Spline_HtoL(iFlavor,iNuAnu,E,C)
    endselect
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_Init
!**********************************************************************!
