module utils
    use types
    private

    public :: linspace

    contains 

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


end module utils
