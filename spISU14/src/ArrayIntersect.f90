!**********************************************************************!
FUNCTION ArrayIntersect(N,A1,A2,Dir)
!----------------------------------------------------------------------!
!Search of the first point where functions given by 1D-arrays A1(N) and
!A2(N) are intersected or converged closer than RelDis
!**********************************************************************!
implicit none

integer:: &
    ArrayIntersect

real,parameter:: &
    RelDis=0
integer &
    Dir,N,n_N
real &
    A1(N),A2(N),R(N)

    R=A1-A2
    selectcase(Dir)
        case(1)!forward-going
            if(abs(R(1))<=abs(R(N)))then
                ArrayIntersect=1
            else
                ArrayIntersect=N
            endif
            do n_N=2,N
                if((R(n_N-1)*R(n_N)<=0))then
                    if(abs(R(n_N-1))<=abs(R(n_N)))then
                        ArrayIntersect=n_N-1
                    else
                        ArrayIntersect=n_N
                    endif
                    exit
                elseif(((A1(n_N).ne.0).and.(abs(R(n_N)/A1(n_N))<RelDis))&
                   .or.((A2(n_N).ne.0).and.(abs(R(n_N)/A2(n_N))<RelDis)))then
                    ArrayIntersect=n_N
                    exit
                endif
            enddo
        case(2)!backward-going
            if(abs(R(1))<abs(R(N)))then
                ArrayIntersect=1
            else
                ArrayIntersect=N
            endif
            do n_N=N-1,1,-1
                if((R(n_N+1)*R(n_N)<=0))then
                    if(abs(R(n_N+1))<=abs(R(n_N)))then
                        ArrayIntersect=n_N+1
                    else
                        ArrayIntersect=n_N
                    endif
                    exit
                elseif(((A1(n_N).ne.0).and.(abs(R(n_N)/A1(n_N))<RelDis))&
                   .or.((A2(n_N).ne.0).and.(abs(R(n_N)/A2(n_N))<RelDis)))then
                    ArrayIntersect=n_N
                    exit
                endif
            enddo
    endselect

    return
!----------------------------------------------------------------------!
endFUNCTION ArrayIntersect
!**********************************************************************!
