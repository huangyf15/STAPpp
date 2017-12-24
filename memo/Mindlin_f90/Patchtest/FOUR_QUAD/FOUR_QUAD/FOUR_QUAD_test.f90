
PROGRAM FOUR_QUAD_TEST

IMPLICIT NONE

REAL(8) :: x(8),y(8)
integer :: ID(6,8),IEN(4,5)
REAL(8) :: a,b,c,d    ! u=0+ax+by, v=0+cx+dy
REAL(8) :: E,v
REAL(8) :: sigx,sigy,txy
REAL(8) :: Fx(4),Fy(4)
integer :: I,J

a=1.0D-5
b=0.D0
c=0.D0
d=0.0D-5

data x/0.0D0, 1.0D0, 1.0D0, 0.0D0, 0.3D0, 0.6D0, 0.8D0, 0.25D0/
data y/0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.3D0, 0.2D0, 0.75D0, 0.9D0/
data ((ID(I,J),I=1,6),J=1,8)/1, 1, 1, 1, 1, 1, &
                             0, 0, 1, 1, 1, 1, &
                             0, 0, 1, 1, 1, 1, &
                             1, 0, 1, 1, 1, 1, &
                             0, 0, 1, 1, 1, 1, &
                             0, 0, 1, 1, 1, 1, &
                             0, 0, 1, 1, 1, 1, &
                             0, 0, 1, 1, 1, 1/
data ((IEN(I,J),I=1,4),J=1,5)/1, 2, 6, 5, &
                              2, 3, 7, 6, &
                              3, 4, 8, 7, &
                              4, 1, 5, 8, &
                              5, 6, 7, 8/

E=3.0D7
v=0.3D0

sigx=E/(1-v*v)*(a+v*d);
sigy=E/(1-v*v)*(d+v*a);
txy=E/(1+v)/2.D0*(b+c)

Fx(1)=-(txy+sigx)/2.D0
Fy(1)=-(txy+sigy)/2.D0
Fx(2)=-(txy-sigx)/2.D0
Fy(2)=(txy-sigy)/2.D0
Fx(3)=(txy+sigx)/2.D0
Fy(3)=(txy+sigy)/2.D0
Fx(4)=(txy-sigx)/2.D0
Fy(4)=-(txy-sigy)/2.D0

OPEN(unit=10 , FILE = "Patch_Test.in", &
               FORM = "FORMATTED", &
               STATUS = "REPLACE", &
               ACCESS = "SEQUENTIAL")

write(10,"('Patch test to 4Q element')")
write(10,"(4(I10))")8,1,1,1
do i=1,8
    write(10,"(7(I10),3(F20.12))")i,(ID(J,I),J=1,6),x(I),y(I),0.0
end do
write(10,"(2(I10))")1,5
write(10,"(2(I10),F20.12,/,2(I10),F20.12)")(I,1,Fx(I),I=2,3)
write(10,"(2(I10),F20.12,/,2(I10),F20.12)")(I,2,Fy(I),I=2,3)
write(10,"(2(I10),F20.12)")4,2,Fy(4)
write(10,"(4(I10))")2,5,1,0
write(10,"(I10,E20.12,F20.12)")1,E,v
do i=1,5
    write(10,"(6(I10))")i,(IEN(J,I),J=1,4),1
end do

CLOSE(10)

END PROGRAM FOUR_QUAD_TEST