!**********************************************************************!
FUNCTION ANS_Preparation_LtoH(Flavor,NuAnu,NC_L,NE_L,C_L,lgE_L,lgF_L,NC_H,NE_H,C_H,lgE_H,lgF_H,&
lgF_MH,n_EI_H_min,lgEI_H_min)
!**********************************************************************!
implicit none

logical:: &
    ANS_Preparation_LtoH
integer:: &
    ArrayIntersectB
real:: &
    AN_Honda,AN_HE_ISU_M,AN_HE_ISU,&
    Splin2_mod

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2
logical &
    bufL
integer &
    NE_L,NC_L,NE_H,NC_H,&
    Flavor,NuAnu,&
    n_NC,n_NE,n_H_M,&
    n_EI_H(NC_H),n_EI_H_min,n_EI_H_max
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    lgF_L(NC_L,NE_L),lgF_H(NC_H,NE_H),lgF_MH(NC_H,NE_H),&
    lgFL_H(NC_H,NE_H),&
    lgEI_H(NC_H),lgEI_H_min,lgEI_H_max

    do n_NE=NE_H,1,-1
        if(lgE_H(n_NE)<=lgE_L(NE_L))then
            n_H_M=n_NE
            exit
        endif
    enddo

    do n_NE=1,n_H_M
        do n_NC=1,NC_H
            lgFL_H(n_NC,n_NE)=log10(AN_Honda(Flavor,NuAnu,10**lgE_H(n_NE),C_H(n_NC)))
        enddo
    enddo

    do n_NC=1,NC_H
        n_EI_H(n_NC)=ArrayIntersectB(n_H_M,lgFL_H(n_NC,:n_H_M),lgF_H(n_NC,:n_H_M))
        lgEI_H(n_NC)=lgE_H(n_EI_H(n_NC))
    enddo

    n_EI_H_min=NE_H
    n_EI_H_max=1
    do n_NC=1,NC_H
        if(n_EI_H(n_NC)<n_EI_H_min)n_EI_H_min=n_EI_H(n_NC)
        if(n_EI_H(n_NC)>n_EI_H_max)n_EI_H_max=n_EI_H(n_NC)
    enddo
    lgEI_H_min=lgE_H(n_EI_H_min)
    lgEI_H_max=lgE_H(n_EI_H_max)

    do n_NC=1,NC_H
        do n_NE=n_EI_H_min,n_EI_H(n_NC)-1
            lgF_MH(n_NC,n_NE)=lgFL_H(n_NC,n_NE)
        enddo
        do n_NE=n_EI_H(n_NC),NE_H
            lgF_MH(n_NC,n_NE)=lgF_H(n_NC,n_NE)
        enddo
    enddo

    ANS_Preparation_LtoH=.true.
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_Preparation_LtoH
!**********************************************************************!

!**********************************************************************!
FUNCTION ANS_Preparation_HtoL(Flavor,NuAnu,NC_L,NE_L,C_L,lgE_L,lgF_L,NC_H,NE_H,C_H,lgE_H,lgF_H,&
lgF_ML,n_EI_L_max,lgEI_L_max)
!**********************************************************************!
implicit none

logical:: &
    ANS_Preparation_HtoL
integer:: &
    ArrayIntersectB
real:: &
    AN_Honda_M,AN_HE_ISU,&
    Splin2_mod

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2
logical &
    bufL
integer &
    NE_L,NC_L,NE_H,NC_H,&
    Flavor,NuAnu,&
    n_NC,n_NE,&
    n_L_M,NE_M,&
    n_EI_L(NC_L),n_EI_L_min,n_EI_L_max
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    lgF_L(NC_L,NE_L),lgF_H(NC_H,NE_H),lgF_ML(NC_L,NE_L),&
    lgFH_L(NC_L,NE_L),&
    lgEI_L(NC_L),lgEI_L_min,lgEI_L_max

    do n_NE=1,NE_L
        if(lgE_L(n_NE)>=lgE_H(1))then
            n_L_M=n_NE
            exit
        endif
    enddo
    NE_M=NE_L-n_L_M+1

    do n_NE=n_L_M,NE_L
        do n_NC=1,NC_L
            lgFH_L(n_NC,n_NE)=log10(AN_HE_ISU(Flavor,NuAnu,10**lgE_L(n_NE),C_L(n_NC)))
        enddo
    enddo

    do n_NC=1,NC_L
        n_EI_L(n_NC)=n_L_M+ArrayIntersectB(NE_M,lgF_L(n_NC,n_L_M:),lgFH_L(n_NC,n_L_M:))-1
        lgEI_L(n_NC)=lgE_L(n_EI_L(n_NC))
    enddo

    n_EI_L_min=NE_L
    n_EI_L_max=1
    do n_NC=1,NC_L
        if(n_EI_L(n_NC)<n_EI_L_min)n_EI_L_min=n_EI_L(n_NC)
        if(n_EI_L(n_NC)>n_EI_L_max)n_EI_L_max=n_EI_L(n_NC)
    enddo
    lgEI_L_min=lgE_L(n_EI_L_min)
    lgEI_L_max=lgE_L(n_EI_L_max)

    do n_NC=1,NC_L
        write(*,*)n_EI_L(n_NC)
        do n_NE=1,n_EI_L(n_NC)
            lgF_ML(n_NC,n_NE)=lgF_L(n_NC,n_NE)
        enddo
        do n_NE=n_EI_L(n_NC)+1,n_EI_L_max
            lgF_ML(n_NC,n_NE)=lgFH_L(n_NC,n_NE)
        enddo
    enddo

    ANS_Preparation_HtoL=.true.
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_Preparation_HtoL
!**********************************************************************!
