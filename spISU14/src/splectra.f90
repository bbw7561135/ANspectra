!**********************************************************************!
FUNCTION dFAN_dE_Init(Mode)
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
    dFAN_dE_Init,AN_spectra
real:: &
    dFAN_dE,&
    Spectrum_Intersection,AN_Honda_M,AN_HE_ISU,AN_SHE_cut

integer,parameter:: &
    Exp  =5,&                                                          !Neutrino experiment number
    PNmod=2                                                            !Prompt neutrino production model
real,parameter:: &
    Honda11_ini=1.0d-01, Honda11_fin=1.0d+04,&
    HE_ISU_ini =1.0d+01, HE_ISU_fin =1.0d+08,&
    firstep=0.01,&
    lastep=100                                                         !Neutrino energy step from the end of spectrum, defined the power tale lean
character(*),parameter:: &
    SA='max'                                                           !Solar activity

logical &
    bufL
integer &
    NuAnu,Flavor,Mode,&
    iNuAnu,iFlavor
real &
    F_L_ini,F_L_sec,F_H_pen,F_H_fin,a,b,F0,gamm,&
    E,C

save
real &
    E_LH(2,2),&                                                        !Boundary between low- and high-energy spectrum application
    E_L_ini,E_L_sec,E_H_pen,E_H_fin

!    do Flavor=1,2
!        do NuAnu=1,2
!            E_LH(Flavor,NuAnu)=Spectrum_Intersection(Flavor,NuAnu,Mode)
!        enddo
!    enddo
    write(*,*) 'Atmospheric neutrino model:'
    selectcase(Mode)
        case(1)
            write(*,*) 'AN_Honda11 + AN_HE_ISU + AN_SHE'
            bufL=AN_spectra()
                E_L_ini=Honda11_ini
                E_L_sec=E_L_ini+firstep
                if(Honda11_fin<HE_ISU_ini)then
                    stop 'dFAN_dE ERROR: Honda11_fin<HE_ISU_ini'
                !elseif(E_LH>Honda11_fin.or.E_LH<HE_ISU_ini)then
                !    write(*,'("dFAN_dE ERROR: E_LH=",1PE8.1,1X,"must be in between")') E_LH
                !    write(*,'("Honda11_fin=",1PE8.1,1X,"and HE_ISU_ini=",1PE8.1)') Honda11_fin,HE_ISU_ini
                !    stop
                endif
                E_H_fin=HE_ISU_fin
                E_H_pen=E_H_fin-lastep
        case(2)
            stop 'dFAN_dE ERROR: CORT! This case is under construction!'
!            write(*,*) 'CORTout'
!            call AN_CORTout(Exp,SA,PNmod)
    endselect
!----------------------------------------------------------------------!
    dFAN_dE_Init=.true.
    return

!**********************************************************************!
ENTRY dFAN_dE(iFlavor,iNuAnu,E,C,Mode)
!----------------------------------------------------------------------!
    selectcase(Mode)
        case(1)
            if(E<0 .or.E>E_GZK)then
                dFAN_dE=0
            elseif(E<E_L_ini)then
                F_L_ini=AN_Honda_M(iFlavor,iNuAnu,E_L_ini,C)
                F_L_sec=AN_Honda_M(iFlavor,iNuAnu,E_L_sec,C)
                a=(F_L_sec-F_L_ini)/(E_L_sec-E_L_ini)
                b=F_L_ini-a*E_L_ini
                dFAN_dE=a*E+b
            elseif(E<=E_H_fin)then
                dFAN_dE=AN_Honda_M(iFlavor,iNuAnu,E,C)
            else
                F_H_pen=AN_HE_ISU(iFlavor,iNuAnu,E_H_pen,C)
                F_H_fin=AN_HE_ISU(iFlavor,iNuAnu,E_H_fin,C)
                call PowerLaw(E_H_pen,E_H_fin,F_H_pen,F_H_fin,gamm,F0)
                dFAN_dE=F0*E**(-gamm)*AN_SHE_cut(E)
        endif
        case(2)
            stop 'dFAN_dE ERROR: CORT! This case is under construction!'
    endselect
    return
!----------------------------------------------------------------------!
endFUNCTION dFAN_dE_Init
!**********************************************************************!


!**********************************************************************!
FUNCTION Spectrum_Intersection(iFlavor,iNuAnu,Mode)
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!
implicit none

real:: &
    Spectrum_Intersection

integer &
    iNuAnu,iFlavor,Mode

    selectcase(Mode)
        case(1)
            selectcase(iFlavor)
                case(1)
                    selectcase(iNuAnu)
                        case(1)
                            Spectrum_Intersection=1.0d+04
                        case(2)
                            Spectrum_Intersection=4.0d+03
                    endselect
                case(2)
                    selectcase(iNuAnu)
                        case(1)
                            Spectrum_Intersection=1.0d+04
                        case(2)
                            Spectrum_Intersection=4.0d+02
                    endselect
            endselect
        case(2)
            stop 'Spectrum_Intersection ERROR: CORT! This case is under construction!'
    endselect
    return
!----------------------------------------------------------------------!
endFUNCTION Spectrum_Intersection
!**********************************************************************!
