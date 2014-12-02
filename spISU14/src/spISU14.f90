!**********************************************************************!
PROGRAM spISU14
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!

implicit none

logical:: &
    dFAN_dE_Init
real:: &
    dFAN_dE_1,dFAN_dE_2

integer,parameter:: &
    Nfsp=140,&
    NC  =10,&
    NE  =61
real,parameter:: &
    C_min= 0.05, C_max= 0.95
character(*),parameter:: &
    usage='Usage: ./spISU14 "outputfile" Flavor[1,2] NuAnu[1,2] Mode[1] E_min[GeV] E_max[GeV]',&
    example='e.g. ./spISU14 "output.dat" 1 2 1 0.1 100'
logical &
    bufL
integer &
    NuAnu,Flavor,Mode,n_NC,n_NE
real &
    E,E3,E_min,E_max,lgE_min,lgE_max,stepC,steplgE,&
    Carr(NC)/NC*0/,PhiE3(NC)/NC*0/
character*80 &
    arg,outfile
character*3 &
    Sp
character*1 &
    SA,NAn(2)/'n','a'/,Fln(2)/'e','m'/

!reading arguments-----------------------------------------------------!
    if(iargc()<6)then
        write(*,*) 'spISU14 ERROR: Missing arguments!'
        write(*,*) usage
        write(*,*) example
        stop
    endif
    call GetArg(1,outfile)
    call GetArg(2,arg); read(arg,*) Flavor
    call GetArg(3,arg); read(arg,*) NuAnu
    call GetArg(4,arg); read(arg,*) Mode
    call GetArg(5,arg); read(arg,*) E_min
    call GetArg(6,arg); read(arg,*) E_max
!echo------------------------------------------------------------------!
    write(*,*) 'Output to file: ',outfile
    write(*,'("Set to:",1X,2A1,1X,"in mode",1X,I1)') Fln(Flavor),NAn(NuAnu),Mode
    write(*,'(2(A6,1PE8.1,1X))') 'E_min=',E_min,'E_max=',E_max

    lgE_min=log10(E_min)
    lgE_max=log10(E_max)
    steplgE=(lgE_max-lgE_min)/(NE-1)
    stepC  =(C_max-C_min)/(NC-1)

    do n_NC=1,NC
        Carr(n_NC)=C_min+(n_NC-1)*stepC
    enddo

    bufL=dFAN_dE_Init(Mode)

    selectcase(Mode)
        case(1)
            Sp='H11'; write(*,*) 'AN_Honda11 + AN_ISU_HE + AN_SHE'
            SA='x'; write(*,*) 'Maximal solar activity'
        case(2)
            Sp='CRT'; write(*,*) ' CORTout '
            SA='n'; write(*,*) ' Minimal solar activity '
    endselect

    open(Nfsp,file=outfile)
    write(Nfsp,'(14(1PE16.8))') 0.,Carr
    do n_NE=1,NE
        E =10**(lgE_min+(n_NE-1)*steplgE)
        E3=E**3
        do n_NC=1,NC
            PhiE3(n_NC)=dFAN_dE_1(Flavor,NuAnu,E,Carr(n_NC),Mode)*E3
        enddo
    write(Nfsp,'(14(1PE16.8))') E,PhiE3
    enddo
    close(Nfsp)

    stop 'spISU14 finished'
!----------------------------------------------------------------------!
endPROGRAM spISU14
!**********************************************************************!
