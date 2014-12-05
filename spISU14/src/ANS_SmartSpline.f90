!**********************************************************************!
recursive FUNCTION ANS_SmartSpline(Trans,Dir)
!**********************************************************************!
implicit none

logical:: &
    ANS_SmartSpline,&
    ANS_SplinePreparation_HtoL,ANS_SplinePreparation_LtoH,&
    ANS_Read
real:: &
    ANS_Spline_LE,ANS_Spline_HtoL,ANS_Spline_LtoH,ANS_Spline_HE,&
    Splin2_mod

save

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2,&
    NE_L=101, NC_L=20,&
    NE_H=71,  NC_H=11
logical &
    bufL
integer &
    Trans,Dir,&
    Flavor,NuAnu,iFlavor,iNuAnu,&
    n_NC,n_NE,&
    n_EI_L_min(NFlavor,NNuAnu),n_EI_L_max(NFlavor,NNuAnu),&
    n_EI_H_min(NFlavor,NNuAnu),n_EI_H_max(NFlavor,NNuAnu)
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    lgF_L(NFlavor,NNuAnu,NC_L,NE_L),ClgF_L(NFlavor,NNuAnu,NC_L,NE_L),&
    lgF_H(NFlavor,NNuAnu,NC_H,NE_H),ClgF_H(NFlavor,NNuAnu,NC_H,NE_H),&
    lgF_MinL(NFlavor,NNuAnu,NC_L,NE_L),ClgF_MinL(NFlavor,NNuAnu,NC_L,NE_L),&
    lgF_MinH(NFlavor,NNuAnu,NC_H,NE_H),ClgF_MinH(NFlavor,NNuAnu,NC_H,NE_H),&
    lgEI_L_min(NFlavor,NNuAnu),lgEI_L_max(NFlavor,NNuAnu),&
    lgEI_H_min(NFlavor,NNuAnu),lgEI_H_max(NFlavor,NNuAnu),&
    E0,C0
character*80 &
    fname

    fname='../input/AN_Honda11.data'
    bufL=ANS_Read(len_trim(fname),fname,NC_L,NE_L,C_L,lgE_L,lgF_L)

    fname='../input/AN_HE_ISU14.data'
    bufL=ANS_Read(len_trim(fname),fname,NC_H,NE_H,C_H,lgE_H,lgF_H)

    selectcase(Trans)
        case(1)
            do Flavor=1,NFlavor
                do NuAnu=1,NNuAnu
                    call Splie2_mod(C_L,lgE_L,lgF_L(Flavor,NuAnu,:,:),NC_L,NE_L,ClgF_L(Flavor,NuAnu,:,:),Test)
                    
                    bufL=ANS_SplinePreparation_LtoH(Flavor,NuAnu,NC_L,NE_L,C_L,lgE_L,lgF_L(Flavor,NuAnu,:,:),&
                    NC_H,NE_H,C_H,lgE_H,lgF_H(Flavor,NuAnu,:,:),lgF_MinH(Flavor,NuAnu,:,:),&
                    n_EI_H_min(Flavor,NuAnu),lgEI_H_min(Flavor,NuAnu),Dir)
                    
                    call Splie2_mod(C_H,lgE_H(n_EI_H_min(Flavor,NuAnu):),lgF_MinH(Flavor,NuAnu,:,n_EI_H_min(Flavor,NuAnu):),&
                    NC_H,NE_H-n_EI_H_min(Flavor,NuAnu)+1,ClgF_MinH(Flavor,NuAnu,:,n_EI_H_min(Flavor,NuAnu):),Test)
                enddo
            enddo
        case(2)
            do Flavor=1,NFlavor
                do NuAnu=1,NNuAnu
                    call Splie2_mod(C_H,lgE_H,lgF_H(Flavor,NuAnu,:,:),NC_H,NE_H,ClgF_H(Flavor,NuAnu,:,:),Test)
                    
                    bufL=ANS_SplinePreparation_HtoL(Flavor,NuAnu,NC_L,NE_L,C_L,lgE_L,lgF_L(Flavor,NuAnu,:,:),&
                    NC_H,NE_H,C_H,lgE_H,lgF_H(Flavor,NuAnu,:,:),lgF_MinL(Flavor,NuAnu,:,:),&
                    n_EI_L_max(Flavor,NuAnu),lgEI_L_max(Flavor,NuAnu),Dir)
                    
                    call Splie2_mod(C_L,lgE_L(:n_EI_L_max(Flavor,NuAnu)),lgF_MinL(Flavor,NuAnu,:,:n_EI_L_max(Flavor,NuAnu)),&
                    NC_L,n_EI_L_max(Flavor,NuAnu),ClgF_MinL(Flavor,NuAnu,:,:n_EI_L_max(Flavor,NuAnu)),Test)
                enddo
            enddo
    endselect

    ANS_SmartSpline=.true.
    return

!----------------------------------------------------------------------!
ENTRY ANS_Spline_HtoL(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
    if(log10(E0)<lgEI_L_max(iFlavor,iNuAnu))then
        ANS_Spline_HtoL=&
        10**(Splin2_mod(C_L,lgE_L(:n_EI_L_max(iFlavor,iNuAnu)),lgF_MinL(iFlavor,iNuAnu,:,:n_EI_L_max(iFlavor,iNuAnu)),&
        ClgF_MinL(iFlavor,iNuAnu,:,:n_EI_L_max(iFlavor,iNuAnu)),NC_L,n_EI_L_max(iFlavor,iNuAnu),C0,log10(E0)))
    else
        ANS_Spline_HtoL=ANS_Spline_HE(iFlavor,iNuAnu,E0,C0)
    endif
    return

!----------------------------------------------------------------------!
ENTRY ANS_Spline_LtoH(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
    if(log10(E0)<=lgEI_H_min(iFlavor,iNuAnu))then
        ANS_Spline_LtoH=ANS_Spline_LE(iFlavor,iNuAnu,E0,C0)
    else
        ANS_Spline_LtoH=&
        10**(Splin2_mod(C_H,lgE_H(n_EI_H_min(iFlavor,iNuAnu):),lgF_MinH(iFlavor,iNuAnu,:,n_EI_H_min(iFlavor,iNuAnu):),&
        ClgF_MinH(iFlavor,iNuAnu,:,n_EI_H_min(iFlavor,iNuAnu):),NC_H,NE_H-n_EI_H_min(iFlavor,iNuAnu)+1,abs(C0),log10(E0)))
    endif
    return

!----------------------------------------------------------------------!
ENTRY ANS_Spline_LE(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
    ANS_Spline_LE=10**(Splin2_mod(C_L,lgE_L,lgF_L(iFlavor,iNuAnu,:,:),ClgF_L(iFlavor,iNuAnu,:,:),NC_L,NE_L,C0,log10(E0)))
    return

!----------------------------------------------------------------------!
ENTRY ANS_Spline_HE(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
    ANS_Spline_HE=10**(Splin2_mod(C_H,lgE_H,lgF_H(iFlavor,iNuAnu,:,:),ClgF_H(iFlavor,iNuAnu,:,:),NC_H,NE_H,abs(C0),log10(E0)))
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_SmartSpline
!**********************************************************************!

!**********************************************************************!
FUNCTION ANS_Read(FN_len,FileName,NC,NE,C,logE,logF)
!**********************************************************************!
implicit none

logical:: &
    ANS_Read

integer,parameter:: &
    NFile=200,&
    NFlavor=2,NNuAnu=2
integer &
    FN_len,&
    Flavor,NuAnu,&
    NC,NE,&
    n_NC,n_NE
real &
    C(NC),E(NE),logE(NE),&
    T(NC+1,NE+1),logF(NFlavor,NNuAnu,NC,NE)
character(FN_len) &
    FileName

    logF=0
    open(NFile,file=FileName)
    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            read(NFile,*) T
            do n_NE=1,NE
                do n_NC=1,NC
                    logF(Flavor,NuAnu,n_NC,n_NE)=log10(T(NC+2-n_NC,n_NE+1))
                enddo
            enddo
        enddo
    enddo
    close(NFile)
    do n_NC=1,NC
        C(n_NC)=T(NC+2-n_NC,1)
    enddo
    do n_NE=1,NE
        E(n_NE)=T(1,n_NE+1)
    enddo
    write(*,*) ' File has been read: ',FileName
    logE=log10(E)

    ANS_Read=.true.
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_Read
!**********************************************************************!
