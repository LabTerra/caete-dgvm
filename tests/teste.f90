program teste
implicit none

real :: a =10.0, b = 5.0, c

c = a + b
b = func(c, a)

print *, a
print *, b

contains

function func(v1, v2) result(v3)
    real, intent(in) :: v1, v2
    real :: v3

    v3 = v1+v2
end function func

end program teste
