!**********************************************************************!
FUNCTION Spectrum_Intersection(iFlavor,iNuAnu)
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!
implicit none

real:: &
    Spectrum_Intersection

integer &
    iNuAnu,iFlavor

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
    return
!----------------------------------------------------------------------!
endFUNCTION Spectrum_Intersection
!**********************************************************************!
