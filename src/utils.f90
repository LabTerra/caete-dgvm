module utils
    use types
    implicit none
    private

    public :: linspace
    public :: process_id
    public :: abort_on_nan
    public :: abort_on_inf
    public :: ascii2bin
    public :: leap
    public :: read_bin

    contains

     !=================================================================
     !=================================================================

    subroutine ascii2bin(file_in, file_out, nx1, ny1)
        ! DEPRECATED
        use types
        !implicit none

        character(len=100), intent(in) :: file_in, file_out
        integer(i_4),intent(in) :: nx1, ny1

        integer(i_4) :: i, j

        real(r_4),allocatable,dimension(:,:) :: arr_in

        allocate(arr_in(nx1,ny1))

        open (unit=11,file=file_in,status='old',form='formatted',access='sequential',&
             action='read')


        open (unit=21,file=file_out,status='unknown',&
             form='unformatted',access='direct',recl=nx1*ny1*r_4)


        do j = 1, ny1 ! for each line do
           read(11,*) (arr_in(i,j), i=1,nx1) ! read all elements in line j (implicitly looping)
           !write(*,*) arr_in(:,j)
        end do

        write(21,rec=1) arr_in
        close(11)
        close(21)

     end subroutine ascii2bin

     !=================================================================
     !=================================================================

     function leap(year) result(is_leap)
        ! Why this is here?
        use types
        !implicit none

        integer(i_4),intent(in) :: year
        logical(l_1) :: is_leap

        logical(l_1) :: by4, by100, by400

        by4 = (mod(year,4) .eq. 0)
        by100 = (mod(year,100) .eq. 0)
        by400 = (mod(year,400) .eq. 0)

        is_leap = by4 .and. (by400 .or. (.not. by100))

     end function leap

     !=================================================================
     !=================================================================

    subroutine abort_on_nan(arg1, abrt)
        logical(l_1), intent(in) :: abrt
        real(r_4), intent(inout) :: arg1
        if(isnan(arg1) .and. .not. abrt ) then
            arg1 = 0.0
        endif
        if(abrt) then
            if (isnan(arg1))then
                call abort()
            endif
        endif
    end subroutine abort_on_nan

     !=================================================================
     !=================================================================

    subroutine abort_on_inf(arg1, abrt)
        logical(l_1), intent(in) :: abrt
        real(r_4), intent(inout) :: arg1
        if(arg1 .ge.arg1 - 1 .and. .not. abrt ) then
             arg1 = 0.0
        endif
        if(abrt) then
            if (arg1 .ge.arg1 - 1)then
                call abort()
            endif
        endif
    end subroutine abort_on_inf

     !=================================================================
     !=================================================================

    subroutine linspace(from, to, array)

        real(r_4), intent(in) :: from, to
        real(r_4), intent(out) :: array(:)
        real(r_4) :: range
        integer :: n, i
        n = size(array)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if

        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
end subroutine linspace

     !=================================================================
     !=================================================================

function process_id() result(ipid)
    ! Identify process number

       integer(r_4) :: ipid

       ipid = getpid()

    end function process_id

     !=================================================================
     !=================================================================

    function read_bin(ny, nx) result(dt)
        ! pft_par agora retorna um array-shape(par,npls)
        use types, only: r_4
        use global_par, only: ntraits, npls
        integer(i_4), intent(in) :: nx, ny

        real(kind=r_4), dimension(ny,nx) :: dt

        ! ['g1','vcmax','tleaf','twood','troot','aleaf','awood','aroot']
        !     dt1 = g1
        !     dt2 = vcmax
        !     dt3 = tleaf
        !     dt4 = twood
        !     dt5 = tfroot
        !     dt6 = aleaf
        !     dt7 = awood
        !     dt8 = aroot

        open(45,file='pls_ex.bin',status='old',&
             &form='unformatted',access='direct',recl=4*ny*nx)
        read(45,rec=1) dt
        close(45)
        return
      end function read_bin

end module utils
