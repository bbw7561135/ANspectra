!**********************************************************************!
recursive FUNCTION AN_spectra()
!**********************************************************************!
implicit none

logical:: &
    AN_spectra
integer:: &
    ArrayIntersectF
real:: &
    AN_Honda_M,AN_HE_ISU,&
    Splin2_mod

save

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFile=200,&
    NFlavor=2,NNuAnu=2,&
    NE_L=101, NC_L=20,&
    NE_H=71,  NC_H=11,&
    NtotalL=NC_L*NE_L
real,parameter:: &
    E_L_min= 1.0d-01, E_L_max= 1.0d+04,&
    C_L_min=-0.95,    C_L_max= 0.95,&
    E_H_min= 1.0d+01, E_H_max= 1.0d+08,&
    C_H_min= 0.0,     C_H_max= 1.0,&
    Rescale=1.0d-04
integer &
    Flavor,NuAnu,iFlavor,iNuAnu,&
    n_NC,n_NE,&
    n_L_M,NEM_L,NEHM(NFlavor,NNuAnu),&
    in_EI_L(NC_L),in_EI_L_min,in_EI_L_max,n_EI_L_min(NFlavor,NNuAnu),n_EI_L_max(NFlavor,NNuAnu),&
    in_EI_HE_LC(NC_L),in_EI_HE_LC_min,in_EI_HE_LC_max,n_EI_HE_LC_min(NFlavor,NNuAnu),n_EI_HE_LC_max(NFlavor,NNuAnu)
real &
    C_L(NC_L), E_L(NE_L), lgE_L(NE_L), T_L(NC_L+1,NE_L+1),&
    C_H(NC_H), E_H(NE_H), lgE_H(NE_H), T_H(NC_H+1,NE_H+1),&
    F_L(NFlavor,NNuAnu,NC_L,NE_L), lgF_L(NFlavor,NNuAnu,NC_L,NE_L), ClgF_L(NFlavor,NNuAnu,NC_L,NE_L),&
    F_H(NFlavor,NNuAnu,NC_H,NE_H), lgF_H(NFlavor,NNuAnu,NC_H,NE_H), ClgF_H(NFlavor,NNuAnu,NC_H,NE_H),&
    FM_L(NFlavor,NNuAnu,NC_L,NE_L),lgF_M(NFlavor,NNuAnu,NC_L,NE_L), ClgF_M(NFlavor,NNuAnu,NC_L,NE_L),&
    iFH_L(NC_L,NE_L)/NtotalL*0/,&
    ilgF_L(NC_L,NE_L),ilgF_H(NC_H,NE_H),ilgF_M(NC_L,NE_L),iClgF_L(NC_L,NE_L),iClgF_H(NC_H,NE_H),iClgF_M(NC_L,NE_L),&
    iEI_L(NC_L),iEI_min,iEI_max,EI_min(NFlavor,NNuAnu),EI_max(NFlavor,NNuAnu),&
    E_M_min,E_M_max,&
    E0,C0

    open(NFile,file='../input/AN_Honda11.data')
    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            read(NFile,*) T_L
            do n_NE=1,NE_L
                do n_NC=1,NC_L
                    F_L(Flavor,NuAnu,n_NC,n_NE)=T_L(NC_L+2-n_NC,n_NE+1)
                enddo
            enddo
        enddo
    enddo
    close(NFile)
    do n_NC=1,NC_L
        C_L(n_NC)=T_L(NC_L+2-n_NC,1)
    enddo
    do n_NE=1,NE_L
        E_L(n_NE)=T_L(1,n_NE+1)
    enddo
    write(*,*) ' File AN_Honda11.data has been read '
    lgE_L=log10(E_L)
    lgF_L=log10(F_L)

    open(NFile,file='../input/AN_HE_ISU14.data')
    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            read(NFile,*) T_H
            do n_NE=1,NE_H
                do n_NC=1,NC_H
                    F_H(Flavor,NuAnu,n_NC,n_NE)=T_H(NC_H+2-n_NC,n_NE+1)
                enddo
            enddo
        enddo
    enddo
    close(NFile)
    do n_NC=1,NC_H
        C_H(n_NC)=T_H(NC_H+2-n_NC,1)
    enddo
    do n_NE=1,NE_H
        E_H(n_NE)=T_H(1,n_NE+1)
    enddo
    write(*,*) ' File AN_HE_ISU14.data has been read '
    lgE_H=log10(E_H)
    lgF_H=log10(F_H)

    E_M_min=E_H_min
    E_M_max=E_L_max

    do n_NE=1,NE_L
        if(E_L(n_NE)>=E_M_min)then
            n_L_M=n_NE
            exit
        endif
    enddo
!    do n_NE=NE_H,1,-1
!        if(E_H(n_NE)<=E_M_max)then
!            n_H_M=n_NE
!            exit
!        endif
!    enddo
    NEM_L=NE_L-n_L_M+1

    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            ilgF_L=lgF_L(Flavor,NuAnu,:,:)
            ilgF_H=lgF_H(Flavor,NuAnu,:,:)
!            call Splie2_mod(C_L,lgE_L,ilgF_L,NC_L,NE_L,iClgF_L,Test)
            call Splie2_mod(C_H,lgE_H,ilgF_H,NC_H,NE_H,iClgF_H,Test)
            do n_NE=n_L_M,NE_L
                do n_NC=1,NC_L
                    iFH_L(n_NC,n_NE)=AN_HE_ISU(Flavor,NuAnu,E_L(n_NE),C_L(n_NC))
                enddo
            enddo
            do n_NC=1,NC_L
                in_EI_L(n_NC)=n_L_M&
                +ArrayIntersectF(NEM_L,F_L(Flavor,NuAnu,n_NC,n_L_M:NE_L),iFH_L(n_NC,n_L_M:NE_L))-1
                iEI_L(n_NC)=E_L(in_EI_L(n_NC))
                do n_NE=1,NE_H
                    if(E_H(n_NE)>=iEI_L(n_NC))then
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
            iEI_min=E_L(in_EI_L_min)
            iEI_max=E_L(in_EI_L_max)
            in_EI_HE_LC_min=NE_H
            in_EI_HE_LC_max=1
            do n_NC=1,NC_L
                if(in_EI_HE_LC(n_NC)<in_EI_HE_LC_min)in_EI_HE_LC_min=in_EI_HE_LC(n_NC)
                if(in_EI_HE_LC(n_NC)>in_EI_HE_LC_max)in_EI_HE_LC_max=in_EI_HE_LC(n_NC)
            enddo
            do n_NC=1,NC_L
                write(*,*)in_EI_L(n_NC)
                do n_NE=1,in_EI_L(n_NC)
                    FM_L(Flavor,NuAnu,n_NC,n_NE)=F_L(Flavor,NuAnu,n_NC,n_NE)
                    ilgF_M(n_NC,n_NE)=log10(FM_L(Flavor,NuAnu,n_NC,n_NE))
                enddo
                do n_NE=in_EI_L(n_NC)+1,in_EI_L_max
                    FM_L(Flavor,NuAnu,n_NC,n_NE)=iFH_L(n_NC,n_NE)
                    ilgF_M(n_NC,n_NE)=log10(FM_L(Flavor,NuAnu,n_NC,n_NE))
                enddo
            enddo
            lgF_M(Flavor,NuAnu,:,:)=ilgF_M
            ClgF_L(Flavor,NuAnu,:,:)=iClgF_L
            ClgF_H(Flavor,NuAnu,:,:)=iClgF_H
            ClgF_M(Flavor,NuAnu,:,:)=iClgF_M
            n_EI_L_min(Flavor,NuAnu)=in_EI_L_min
            n_EI_L_max(Flavor,NuAnu)=in_EI_L_max
            n_EI_HE_LC_min(Flavor,NuAnu)=in_EI_HE_LC_min
            n_EI_HE_LC_max(Flavor,NuAnu)=in_EI_HE_LC_max
            EI_min(Flavor,NuAnu)=iEI_min
            EI_max(Flavor,NuAnu)=iEI_max
            write(*,*)
        enddo
    enddo
!    NEHM(Flavor,NuAnu)=NE_H-n_EI_HE_LC_max(Flavor,NuAnu)+1

    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            call Splie2_mod(C_L,lgE_L(:n_EI_L_max(Flavor,NuAnu)),lgF_M(Flavor,NuAnu,:,:n_EI_L_max(Flavor,NuAnu)),&
            NC_L,n_EI_L_max(Flavor,NuAnu),ClgF_M(Flavor,NuAnu,:,:n_EI_L_max(Flavor,NuAnu)),Test)
!            call Splie2_mod(C_H,lgE_H(n_EI_HE_LC_max(Flavor,NuAnu):),lgF_H(Flavor,NuAnu,:,n_EI_HE_LC_max(Flavor,NuAnu):),&
!            NC_H,NEHM(Flavor,NuAnu),ClgF_H(Flavor,NuAnu,:,n_EI_HE_LC_max(Flavor,NuAnu):),Test)
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
    if(E0<EI_max(iFlavor,iNuAnu))then
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
!ENTRY AN_HE_ISU_M(iFlavor,iNuAnu,E0,C0)
!----------------------------------------------------------------------!
!    AN_HE_ISU_M=10**(Splin2_mod(C_H,lgE_H(n_EI_HE_LC_max(iFlavor,iNuAnu):),lgF_H(iFlavor,iNuAnu,:,n_EI_HE_LC_max(iFlavor,iNuAnu):),&
!    ClgF_H(iFlavor,iNuAnu,:,n_EI_HE_LC_max(iFlavor,iNuAnu):),NC_H,NEHM(Flavor,NuAnu),abs(C0),log10(E0)))
!    return
!----------------------------------------------------------------------!
endFUNCTION AN_spectra
!**********************************************************************!
