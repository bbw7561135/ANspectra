!**********************************************************************!
FUNCTION ArrayIntersectF(N,A1,A2)
!**********************************************************************!
!Forward
!----------------------------------------------------------------------!
implicit none

integer:: &
    ArrayIntersectF

real,parameter:: &
    RelDis=1d-2
integer &
    N,n_N
real &
    A1(N),A2(N),R(N)

    R=A1-A2
    if(abs(R(1))<=abs(R(N)))then
        ArrayIntersectF=1
    else
        ArrayIntersectF=N
    endif
    do n_N=2,N
        if((R(n_N-1)*R(n_N)<=0))then
            if(abs(R(n_N-1))<=abs(R(n_N)))then
                ArrayIntersectF=n_N-1
            else
                ArrayIntersectF=n_N
            endif
            exit
        elseif(((A1(n_N).ne.0).and.(R(n_N)/A1(n_N)<RelDis))&
           .or.((A2(n_N).ne.0).and.(R(n_N)/A2(n_N)<RelDis)))then
            ArrayIntersectF=n_N
            exit
        endif
    enddo

    return
!----------------------------------------------------------------------!
endFUNCTION ArrayIntersectF
!**********************************************************************!

!**********************************************************************!
FUNCTION ArrayIntersectB(N,A1,A2)
!**********************************************************************!
!Backward
!----------------------------------------------------------------------!
implicit none

integer:: &
    ArrayIntersectB

real,parameter:: &
    RelDis=1d-2
integer &
    N,n_N
real &
    A1(N),A2(N),R(N)

    R=A1-A2
    if(abs(R(1))<abs(R(N)))then
        ArrayIntersectB=1
    else
        ArrayIntersectB=N
    endif
    do n_N=N-1,1,-1
        if((R(n_N+1)*R(n_N)<=0))then
            if(abs(R(n_N+1))<=abs(R(n_N)))then
                ArrayIntersectB=n_N+1
            else
                ArrayIntersectB=n_N
            endif
            exit
        elseif(((A1(n_N).ne.0).and.(R(n_N)/A1(n_N)<RelDis))&
           .or.((A2(n_N).ne.0).and.(R(n_N)/A2(n_N)<RelDis)))then
            ArrayIntersectB=n_N
            exit
        endif
    enddo

    return
!----------------------------------------------------------------------!
endFUNCTION ArrayIntersectB
!**********************************************************************!
