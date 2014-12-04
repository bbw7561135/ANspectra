!**********************************************************************!
SUBROUTINE AN_HE_ISU14
!**********************************************************************!
!Fluxes of conventional high-energy atmospheric neutrinos and anti-
!neutrinos as a function of energy and zenith angle.
!
!HE_ISU-2014 FLUX in 1/(cm^2 sec sr GeV)
!
!Energy range is 10 GeV to 100 PeV
!Primary cosmic ray spectra: Hilas-Gaisser
!Hardron interaction model: Kimel-Mokhov
!
!The source data array F is dF/dE in units of 1/(cm^2 s sr GeV). The 
!energy reference points are calculated on the equidistant over lg(E)
!grid and the reference points for cosine of zenith angle are calculated
! on equidistant grid from 0 to 1 with step = 0.1
!**********************************************************************!

use PhysMathConstants, only: E_GZK

implicit none

real::Splin2_mod,Splin2_ED

save

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2, NNuAnu=2,&
    NE_LE=101, NC_LE=20,&
    NE_HE=71,  NC_HE=11,&
    NFile=200
real,parameter:: &
    E_LE_min= 1.0d-01, E_LE_max= 1.0d+04,&
    C_LE_min=-0.95,    C_LE_max= 0.95,&
    E_HE_min= 1.0d+01, E_HE_max= 1.0d+08,&
    C_HE_min= 0.0,     C_HE_max=1.0,&
    Rescale=1.0d-04
integer &
    Flavor,NuAnu,iFlavor,iNuAnu,&
    n_NC,n_NE
real &
    C_LE(NC_LE),E_LE(NE_LE),&
    C_HE(NC_HE),E_HE(NE_HE),T(NC_HE+1,NE_HE+1),&
    F_LE(NFlavor,NNuAnu,NC_LE,NE_LE),CF_LE(NFlavor,NNuAnu,NC_LE,NE_LE),&
    F_HE(NFlavor,NNuAnu,NC_HE,NE_HE),lnF_HE(NFlavor,NNuAnu,NC_HE,NE_HE),ClnF_HE(NFlavor,NNuAnu,NC_HE,NE_HE),&
    F0,E0,C0,stepC,steplgE

    open(NFile,FILE='../input/AN_HE_ISU14.data')
    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            READ(NFile,*) T
            do n_NE=1,NE_HE
                do n_NC=1,NC_HE
                    F_HE(Flavor,NuAnu,n_NC,n_NE)=T(NC_HE+2-n_NC,n_NE+1)
                enddo
            enddo
        enddo
    enddo
    close(NFile)
    do n_NC=1,NC_HE
        C_HE(n_NC)=T(NC_HE+2-n_NC,1)
    enddo
    do n_NE=1,NE_HE
        E_HE(n_NE)=T(1,n_NE+1)
    enddo
    write(*,*) ' File AN_HE_ISU14.data was read '

    steplgE=(log10(E_HE_max)-log10(E_HE_min))/(NE_HE-1)
    stepC=(C_HE_max-C_HE_min)/(NC_HE-1)

    lnF_HE=log(F_HE)
    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            call Splie2_ED(C_HE_min,log10(E_HE_min),stepC,steplgE,&
            lnF_HE(Flavor,NuAnu,:,:),NC_HE,NE_HE,ClnF_HE(Flavor,NuAnu,:,:),Test)
        enddo
    enddo

    return

!----------------------------------------------------------------------!
ENTRY AN_HE_ISU(iFlavor,iNuAnu,F0,E0,C0)
!----------------------------------------------------------------------!
    F0=exp(Splin2_ED(C_HE_min,log10(E_HE_min),stepC,steplgE,&
    lnF_HE(iFlavor,iNuAnu,:,:),ClnF_HE(iFlavor,iNuAnu,:,:),NC_HE,NE_HE,abs(C0),log10(E0)))
    return

endSUBROUTINE AN_HE_ISU14
