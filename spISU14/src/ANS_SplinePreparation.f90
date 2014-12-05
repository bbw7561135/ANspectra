!**********************************************************************!
FUNCTION ANS_SplinePreparation_LtoH(Flavor,NuAnu,NC_L,NE_L,C_L,lgE_L,lgF_L,NC_H,NE_H,C_H,lgE_H,lgF_H,&
lgF_MinH,n_EI_H_min,lgEI_H_min,Dir)
!**********************************************************************!
implicit none

logical:: &
    ANS_SplinePreparation_LtoH
integer:: &
    ArrayIntersect
real:: &
    ANS_Spline_LE,ANS_Spline_LtoH,ANS_Spline_HE,&
    Splin2_mod

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2
logical &
    bufL
integer &
    Dir,NE_L,NC_L,NE_H,NC_H,&
    Flavor,NuAnu,&
    n_NC,n_NE,n_H_M,&
    n_EI_H(NC_H),n_EI_H_min,n_EI_H_max
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    lgF_L(NC_L,NE_L),lgF_H(NC_H,NE_H),lgF_MinH(NC_H,NE_H),&
    lgF_LinH(NC_H,NE_H),&
    lgEI_H(NC_H),lgEI_H_min,lgEI_H_max

    do n_NE=NE_H,1,-1
        if(lgE_H(n_NE)<=lgE_L(NE_L))then
            n_H_M=n_NE
            exit
        endif
    enddo

    do n_NE=1,n_H_M
        do n_NC=1,NC_H
            lgF_LinH(n_NC,n_NE)=log10(ANS_Spline_LE(Flavor,NuAnu,10**lgE_H(n_NE),C_H(n_NC)))
        enddo
    enddo

    do n_NC=1,NC_H
        n_EI_H(n_NC)=ArrayIntersect(n_H_M,lgF_LinH(n_NC,:n_H_M),lgF_H(n_NC,:n_H_M),Dir)
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
            lgF_MinH(n_NC,n_NE)=lgF_LinH(n_NC,n_NE)
        enddo
        do n_NE=n_EI_H(n_NC),NE_H
            lgF_MinH(n_NC,n_NE)=lgF_H(n_NC,n_NE)
        enddo
    enddo

    ANS_SplinePreparation_LtoH=.true.
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_SplinePreparation_LtoH
!**********************************************************************!

!**********************************************************************!
FUNCTION ANS_SplinePreparation_HtoL(Flavor,NuAnu,NC_L,NE_L,C_L,lgE_L,lgF_L,NC_H,NE_H,C_H,lgE_H,lgF_H,&
lgF_MinL,n_EI_L_max,lgEI_L_max,Dir)
!**********************************************************************!
implicit none

logical:: &
    ANS_SplinePreparation_HtoL
integer:: &
    ArrayIntersect
real:: &
    ANS_Spline_HtoL,ANS_Spline_HE,&
    Splin2_mod

logical,parameter:: &
    Test=.false.
integer,parameter:: &
    NFlavor=2,NNuAnu=2
logical &
    bufL
integer &
    Dir,NE_L,NC_L,NE_H,NC_H,&
    Flavor,NuAnu,&
    n_NC,n_NE,&
    n_L_M,NE_M,&
    n_EI_L(NC_L),n_EI_L_min,n_EI_L_max
real &
    C_L(NC_L),lgE_L(NE_L),&
    C_H(NC_H),lgE_H(NE_H),&
    lgF_L(NC_L,NE_L),lgF_H(NC_H,NE_H),lgF_MinL(NC_L,NE_L),&
    lgF_HinL(NC_L,NE_L),&
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
            lgF_HinL(n_NC,n_NE)=log10(ANS_Spline_HE(Flavor,NuAnu,10**lgE_L(n_NE),C_L(n_NC)))
        enddo
    enddo

    do n_NC=1,NC_L
        n_EI_L(n_NC)=n_L_M+ArrayIntersect(NE_M,lgF_L(n_NC,n_L_M:),lgF_HinL(n_NC,n_L_M:),Dir)-1
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
        do n_NE=1,n_EI_L(n_NC)
            lgF_MinL(n_NC,n_NE)=lgF_L(n_NC,n_NE)
        enddo
        do n_NE=n_EI_L(n_NC)+1,n_EI_L_max
            lgF_MinL(n_NC,n_NE)=lgF_HinL(n_NC,n_NE)
        enddo
    enddo

    ANS_SplinePreparation_HtoL=.true.
    return
!----------------------------------------------------------------------!
endFUNCTION ANS_SplinePreparation_HtoL
!**********************************************************************!
