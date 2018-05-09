subroutine norm_args(pin,pout)
    implicit none   
    real*8, dimension(5), intent(in) :: pin
    real*8, dimension(5), intent(out) :: pout
    real*8,  parameter :: pi  = 3.1415926535897932_8
    
    pout(1)=(pin(1)-55.87*pi/180.0d0)/((56.4-55.87)*pi/180.0d0)
    pout(2)=(pin(2)-45.55*pi/180.0d0)/((46.55-45.55)*pi/180.0d0)
    pout(3)=(pin(3)-150.0d0)/(350.0d0-150.0d0)
    pout(4)=(pin(4)-250.0d0)/2250.0d0
    pout(5)=(pin(5)-5.0d0)/(25.0d0)
end

subroutine denorm_args(pin,pout)
    implicit none   
    real*8, dimension(5), intent(in) :: pin
    real*8, dimension(5), intent(out) :: pout
    real*8,  parameter :: pi  = 3.1415926535897932_8
    
    pout(1)=pin(1)*((56.4-55.87)*pi/180.0d0)+55.87*pi/180.0d0
    pout(2)=pin(2)*((46.55-45.55)*pi/180.d0)+45.55*pi/180.0d0
    pout(3)=pin(3)*(350.0d0-150.0d0)+150.0d0
    pout(4)=pin(4)*2250.0d0 + 250.0d0
    pout(5)=pin(5)*(25.0d0)+5.0d0
end

real*8 function resudal_sphere(pn, pk1, pk2, naxes1, naxes2, im1, im2, imm1, imm2, elm1, azm1, elm2, azm2)
    implicit none    
    real*8, dimension(5), intent(in) :: pn
    real*8, dimension(5) :: pd
    real*8, dimension(3), intent(in) :: pk1, pk2
    integer, intent(in) :: naxes1(2), naxes2(2)    
    real*8, dimension( naxes1(1), naxes1(2) ), intent(in) :: im1
    real*8, dimension( naxes2(1), naxes2(2) ), intent(in) :: im2
    real*8, dimension( naxes1(1), naxes1(2) ), intent(inout) :: imm1
    real*8, dimension( naxes2(1), naxes2(2) ), intent(inout) :: imm2
    real*8, dimension( naxes1(1), naxes1(2) ), intent(in) :: elm1, azm1    
    real*8, dimension( naxes2(1), naxes2(2) ), intent(in) :: elm2, azm2
    real*8 :: sum1, sum2
    integer :: i, j
    
    
    
    call denorm_args(pn,pd)
    call sphere_rout(naxes1, pd, pk1, elm1, azm1, &    
           imm1)
    !write(*,*) "MODEL1 max - ", maxval(imm1)    
       
    call sphere_rout(naxes2, pd, pk2, elm2, azm2, &    
           imm2)
    !write(*,*) "MODEL2 max - ", maxval(imm2)

    sum1=0
    do i=1, naxes1(2), 1
        do j=1, naxes1(1),1
            sum1=sum1+(im1(i,j)-imm1(i,j))**2
        enddo        
    enddo
    sum1=sum1*0.12056327160493827d0
    
    sum2=0
    do i=1, naxes2(2), 1
        do j=1, naxes2(1),1
            sum2=sum2+(im2(i,j)-imm2(i,j))**2
        enddo        
    enddo
    sum2=sum2*0.25d0
    
    resudal_sphere=sum1+sum2
end
