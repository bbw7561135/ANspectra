!**********************************************************************!
PROGRAM Translator_Honda
!----------------------------------------------------------------------!
!Translator to standard spectrum data file format for Honda11.dat
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!
implicit none

integer,parameter:: &
    FileI=200,FileO=201,&
    NFlavor=2,NNuAnu=2,&
    NE_L=101, NC_L=20,&
    NtotalL=NFlavor*NNuAnu*NC_L*NE_L
real,parameter:: &
    E_L_min= 1.0d-01, E_L_max= 1.0d+04,&
    C_L_min=-0.95,    C_L_max= 0.95,&
    Rescale=1.0d-04
integer &
    Flavor,NuAnu,&
    n_NC,n_NE
real &
    C_L(NC_L), E_L(NE_L), lgE_L(NE_L),&
    F_L(NFlavor,NNuAnu,NC_L,NE_L)/NtotalL*0/

    open(FileI,file='../input/Honda11.dat')
    do n_NC=1,NC_L
        read(FileI,*) C_L(n_NC)
        do n_NE=1,NE_L
            read(FileI,*) E_L(n_NE),((F_L(Flavor,NuAnu,n_NC,n_NE),NuAnu=1,NNuAnu),Flavor=NFlavor,1,-1)
        enddo
    enddo
    close(FileI)
    write(*,*)' File Honda11.dat has been read '
    F_L=F_L*Rescale

    open(FileO,file='../input/AN_Honda11.data')
    do Flavor=1,NFlavor
        do NuAnu=1,NNuAnu
            write(FileO,'(2X,100(F5.2,10X))') 0.,(C_L(n_NC),n_NC=1,NC_L)
            do n_NE=1,NE_L
                write(FileO,'(100(1PE10.4,5X))') E_L(n_NE),(F_L(Flavor,NuAnu,n_NC,n_NE),n_NC=1,NC_L)
            enddo
        enddo
    enddo
    close(FileO)
    write(*,*) ' File AN_Honda11.data has been created '

    stop 'spISU14 finished'
!----------------------------------------------------------------------!
endPROGRAM Translator_Honda
!**********************************************************************!
