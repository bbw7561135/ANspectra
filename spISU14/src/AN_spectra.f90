!**********************************************************************!
recursive FUNCTION AN_spectra()
!**********************************************************************!
implicit none

logical:: &
    AN_spectra,ANS_Read,ANS_Preparation
integer:: &
    ArrayIntersectF
real:: &
    AN_Honda_M,AN_HE_ISU,&
    Splin2_mod

save

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2,&
    NE_L=101, NC_L=20,&
    NE_H=71,  NC_H=11,&
    NtotalL=NC_L*NE_L
logical &
    bufL
integer &
    Flavor,NuAnu,iFlavor,iNuAnu,&
    n_NC,n_NE,&
    n_L_M,NE_M,&
    n_EI_L_min(NFlavor,NNuAnu),n_EI_L_max(NFlavor,NNuAnu),&
    n_EI_HE_LC_min(NFlavor,NNuAnu),n_EI_HE_LC_max(NFlavor,NNuAnu)
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    lgF_L(NFlavor,NNuAnu,NC_L,NE_L),ClgF_L(NFlavor,NNuAnu,NC_L,NE_L),&
    lgF_H(NFlavor,NNuAnu,NC_H,NE_H),ClgF_H(NFlavor,NNuAnu,NC_H,NE_H),&
    lgF_M(NFlavor,NNuAnu,NC_L,NE_L),ClgF_M(NFlavor,NNuAnu,NC_L,NE_L),&
    lgEI_min(NFlavor,NNuAnu),lgEI_max(NFlavor,NNuAnu),&
    E0,C0
character*80 &
    fname

    fname='../input/AN_Honda11.data'
    bufL=ANS_Read(len_trim(fname),fname,NC_L,NE_L,C_L,lgE_L,lgF_L)

    fname='../input/AN_HE_ISU14.data'
    bufL=ANS_Read(len_trim(fname),fname,NC_H,NE_H,C_H,lgE_H,lgF_H)

    do n_NE=1,NE_L
        if(lgE_L(n_NE)>=lgE_H(1))then
            n_L_M=n_NE
            exit
        endif
    enddo
    NE_M=NE_L-n_L_M+1

    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
!            call Splie2_mod(C_L,lgE_L,lgF_L(Flavor,NuAnu,:,:),NC_L,NE_L,ClgF_L(Flavor,NuAnu,:,:),Test)
            call Splie2_mod(C_H,lgE_H,lgF_H(Flavor,NuAnu,:,:),NC_H,NE_H,ClgF_H(Flavor,NuAnu,:,:),Test)
        enddo
    enddo

    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            bufL=ANS_Preparation(Flavor,NuAnu,NC_L,NE_L,&
            C_L,lgE_L,lgF_L(Flavor,NuAnu,:,:),NC_H,NE_H,C_H,lgE_H,lgF_H(Flavor,NuAnu,:,:),lgF_M(Flavor,NuAnu,:,:),&
            n_EI_L_min(Flavor,NuAnu),n_EI_L_max(Flavor,NuAnu),n_EI_HE_LC_min(Flavor,NuAnu),n_EI_HE_LC_max(Flavor,NuAnu),&
            lgEI_min(Flavor,NuAnu),lgEI_max(Flavor,NuAnu),n_L_M,NE_M)
        enddo
    enddo

    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            call Splie2_mod(C_L,lgE_L(:n_EI_L_max(Flavor,NuAnu)),lgF_M(Flavor,NuAnu,:,:n_EI_L_max(Flavor,NuAnu)),&
            NC_L,n_EI_L_max(Flavor,NuAnu),ClgF_M(Flavor,NuAnu,:,:n_EI_L_max(Flavor,NuAnu)),Test)
        enddo
    enddo

    AN_spectra=.true.
    return

!----------------------------------------------------------------------!
!ENTRY AN_Honda(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
!    AN_Honda=10**(Splin2_mod(C_L,lgE_L,lgF_L(iFlavor,iNuAnu,:,:),ClgF_L(iFlavor,iNuAnu,:,:),NC_L,NE_L,C0,log10(E0)))
!    return
!----------------------------------------------------------------------!
ENTRY AN_Honda_M(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
    if(log10(E0)<lgEI_max(iFlavor,iNuAnu))then
        AN_Honda_M=10**(Splin2_mod(C_L,lgE_L(:n_EI_L_max(iFlavor,iNuAnu)),lgF_M(iFlavor,iNuAnu,:,:n_EI_L_max(iFlavor,iNuAnu)),&
        ClgF_M(iFlavor,iNuAnu,:,:n_EI_L_max(iFlavor,iNuAnu)),NC_L,n_EI_L_max(iFlavor,iNuAnu),C0,log10(E0)))
    else
        AN_Honda_M=AN_HE_ISU(iFlavor,iNuAnu,E0,C0)
    endif
    
    return
!----------------------------------------------------------------------!
ENTRY AN_HE_ISU(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
    AN_HE_ISU=10**(Splin2_mod(C_H,lgE_H,lgF_H(iFlavor,iNuAnu,:,:),ClgF_H(iFlavor,iNuAnu,:,:),NC_H,NE_H,abs(C0),log10(E0)))
    return
!----------------------------------------------------------------------!
endFUNCTION AN_spectra
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

!**********************************************************************!
FUNCTION ANS_Preparation(iFlavor,iNuAnu,NC_L,NE_L,C_L,lgE_L,ilgF_L,NC_H,NE_H,C_H,lgE_H,ilgF_H,&
ilgF_M,in_EI_L_min,in_EI_L_max,in_EI_HE_LC_min,in_EI_HE_LC_max,ilgEI_min,ilgEI_max,n_L_M,NE_M)
!**********************************************************************!
implicit none

logical:: &
    ANS_Preparation
integer:: &
    ArrayIntersectF
real:: &
    AN_Honda_M,AN_HE_ISU,&
    Splin2_mod

save

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2
logical &
    bufL
integer &
    NE_L,NC_L,NE_H,NC_H,&
    iFlavor,iNuAnu,&
    n_NC,n_NE,&
    n_L_M,NE_M,&
    in_EI_L(NC_L),in_EI_L_min,in_EI_L_max,&
    in_EI_HE_LC(NC_L),in_EI_HE_LC_min,in_EI_HE_LC_max
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    ilgF_L(NC_L,NE_L),ilgF_H(NC_H,NE_H),ilgF_M(NC_L,NE_L),&
    ilgFH_L(NC_L,NE_L),&
    ilgEI_L(NC_L),ilgEI_min,ilgEI_max

    do n_NE=n_L_M,NE_L
        do n_NC=1,NC_L
            ilgFH_L(n_NC,n_NE)=log10(AN_HE_ISU(iFlavor,iNuAnu,10**lgE_L(n_NE),C_L(n_NC)))
        enddo
    enddo

    do n_NC=1,NC_L
        in_EI_L(n_NC)=n_L_M+ArrayIntersectF(NE_M,ilgF_L(n_NC,n_L_M:NE_L),ilgFH_L(n_NC,n_L_M:NE_L))-1
        ilgEI_L(n_NC)=lgE_L(in_EI_L(n_NC))
        do n_NE=1,NE_H
            if(lgE_H(n_NE)>=ilgEI_L(n_NC))then
                in_EI_HE_LC(n_NC)=n_NE
                exit
            endif
        enddo
    enddo
    in_EI_L_min=NE_L
    in_EI_L_max=1
    do n_NC=1,NC_L
        if(in_EI_L(n_NC)<in_EI_L_min)in_EI_L_min=in_EI_L(n_NC)
        if(in_EI_L(n_NC)>in_EI_L_max)in_EI_L_max=in_EI_L(n_NC)
    enddo
    ilgEI_min=lgE_L(in_EI_L_min)
    ilgEI_max=lgE_L(in_EI_L_max)
    in_EI_HE_LC_min=NE_H
    in_EI_HE_LC_max=1
    do n_NC=1,NC_L
        if(in_EI_HE_LC(n_NC)<in_EI_HE_LC_min)in_EI_HE_LC_min=in_EI_HE_LC(n_NC)
        if(in_EI_HE_LC(n_NC)>in_EI_HE_LC_max)in_EI_HE_LC_max=in_EI_HE_LC(n_NC)
    enddo
    do n_NC=1,NC_L
        write(*,*)in_EI_L(n_NC)
        do n_NE=1,in_EI_L(n_NC)
            ilgF_M(n_NC,n_NE)=ilgF_L(n_NC,n_NE)
        enddo
        do n_NE=in_EI_L(n_NC)+1,in_EI_L_max
            ilgF_M(n_NC,n_NE)=ilgFH_L(n_NC,n_NE)
        enddo
    enddo

    ANS_Preparation=.true.
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_Preparation
!**********************************************************************!
