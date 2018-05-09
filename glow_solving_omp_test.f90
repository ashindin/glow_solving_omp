program main
    implicit none
    include "omp_lib.h"
    integer :: numProc
    real*8 :: st_time, et_time
    real ( kind = 8 ), external :: resudal_sphere
    character(len=300) :: solve_pars_fn1, solve_pars_fn2
    character(len=300) :: im_fn
    integer :: naxes1(2), naxes2(2)
    integer, dimension(4) :: crop_pars_1, crop_pars_2
    real*8, dimension(:, :), allocatable :: EL_keo
    real*8, dimension(:, :), allocatable :: AZ_keo
    real*8, dimension(:, :), allocatable :: EL_s1c
    real*8, dimension(:, :), allocatable :: AZ_s1c
    real*8, dimension(:, :), allocatable :: im1, im2, imm1, imm2
    real*8,  parameter :: pi = 3.1415926535897932_8
    !real*8,  parameter :: pi  = 4 * atan (1.0_8)
    
    real*8, dimension(5) :: pn, pd, pminn, pmind
    real*8, dimension(3), parameter :: pk1 = (/ 55.9305361*pi/180., 48.7444861*pi/180.,  0.091d0 /)
    real*8, dimension(3), parameter :: pk2 = (/ 56.1501667*pi/180., 46.1050833*pi/180.,  0.183d0 /)
    real*8 :: res_start, res_min    
    real*8, dimension(5) :: step
    integer :: konvge, kcount, icount, numres, ifault
    
    integer :: AllocateStatus, DeAllocateStatus

    !write(*,*) "Test 28"
    
    solve_pars_fn1 = 'keo_140824_solve_manual.pars'
    solve_pars_fn2 = 's1c_140824_solve.pars'  
    im_fn = 'test.fits'
    
    call get_images_naxes(im_fn, &
                      naxes1, naxes2)
    !naxes1 = (/ 511, 511 /)
    !naxes2 = (/ 288, 288 /)
    
    !crop_pars: x_shift, y_shift, width, height
    crop_pars_1 =  (/ 175, 325, 100, 100 /)
    crop_pars_2 =  (/ 0, 0, naxes2(1), naxes2(2) /)    
    
    allocate ( im1( crop_pars_1(4),crop_pars_1(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for im1***" 
    allocate ( im2( crop_pars_2(4),crop_pars_2(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for im2***"
    
    allocate ( imm1( crop_pars_1(4),crop_pars_1(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for im1***" 
    allocate ( imm2( crop_pars_2(4),crop_pars_2(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for im2***"
    
    !write(*,*) "test 56"

    call get_images_data(im_fn, naxes1, naxes2, crop_pars_1, crop_pars_2, & 
            im1, im2)

    !write(*,*) "test 56"
    
    allocate ( AZ_keo(crop_pars_1(4),crop_pars_1(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for AZ_keo***" 
    allocate ( EL_keo(crop_pars_1(4),crop_pars_1(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for EL_keo***" 
    allocate ( AZ_s1c(crop_pars_2(4),crop_pars_2(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for AZ_s1c***" 
    allocate ( EL_s1c(crop_pars_2(4),crop_pars_2(3) ), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for EL_s1c***" 
    
    !write(*,*) "Test 65"

    call get_elaz_matrices(solve_pars_fn1, solve_pars_fn2, crop_pars_1, crop_pars_2, &
                             EL_keo, AZ_keo, EL_s1c, AZ_s1c)
    
    !open(unit = 9, file = 'f_az_keo.bin', access = "STREAM", form = "unformatted", status="replace")
    !write(9) transpose(AZ_keo)
    !close(9)
    !open(unit = 9, file = 'f_el_keo.bin', access = "STREAM", form = "unformatted", status="replace")
    !write(9) transpose(EL_keo)
    !close(9)
    !open(unit = 9, file = 'f_az_s1c.bin', access = "STREAM", form = "unformatted", status="replace")
    !write(9) transpose(AZ_s1c)
    !close(9)
    !open(unit = 9, file = 'f_el_s1c.bin', access = "STREAM", form = "unformatted", status="replace")
    !write(9) transpose(EL_s1c)
    !close(9)
    
    !pd = (/ 56.1520164397*pi/180., 46.1126660813*pi/180., 242.815364568d0, 922.901103716d0, 14.7489062231d0 /)
    pd = (/ 56.1434444*pi/180., 46.0991056*pi/180., 240.0d0, 500.0d0, 10*sqrt(2.0d0) /)
    !pk1 = (/ 55.9305361*pi/180., 48.7444861*pi/180.,  0.091d0 /)
    !pk2 = (/ 56.1501667*pi/180., 46.1050833*pi/180.,  0.183d0 /)
    
    
    call sphere_rout((/ crop_pars_1(4),crop_pars_1(3) /), &
    pd, pk1, EL_keo, AZ_keo, &    
           imm1)
    
    call sphere_rout((/ crop_pars_2(4),crop_pars_2(3) /), &
    pd, pk2, EL_s1c, AZ_s1c, &    
           imm2)
    
    !write(*,*) "G(center)=",imm2(144,144)
    
    open(unit = 9, file = 'model_keo_start.bin', access = "STREAM", form = "unformatted", status="replace")
    write(9) transpose(imm1)
    close(9)
    
    open(unit = 9, file = 'model_s1c_start.bin', access = "STREAM", form = "unformatted", status="replace")
    write(9) transpose(imm2)
    close(9)
    
    call norm_args(pd,pn)
    res_start = resudal_sphere(pn, pk1, pk2, (/ crop_pars_1(4),crop_pars_1(3) /), (/ crop_pars_2(4),crop_pars_2(3) /), &
    im1, im2, imm1, imm2, EL_keo, AZ_keo, EL_s1c, AZ_s1c)

    step(1:5) = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00/)
    konvge = 10
    kcount = 10000
    
    write (*,"(A)") "SPHERE MODEL"
    
    numProc = OMP_GET_NUM_PROCS()
    write(*,'(A32,I4)')"Available number of processors: ", numProc
    
    !call SET_OMP_NUM_THREADS(1)   
    !call SET_OMP_NUM_THREADS(numProc)

    write (*,"(A,F6.2,F6.2,F6.1,F8.1,F6.1, A, F10.3)") "p_start:", pd(1)*180./pi, pd(2)*180./pi, pd(3), &
    pd(4), pd(5), " - ", res_start

    !stop "*** Technical stop ***"

    call norm_args(pd,pn)

    st_time = omp_get_wtime()    
    call nelmin ( resudal_sphere, 5, pn, &
      pminn, res_min, 1.0D-08, step, konvge, kcount, &
      icount, numres, ifault, &
      pk1, pk2, (/ crop_pars_1(4),crop_pars_1(3) /), (/ crop_pars_2(4),crop_pars_2(3) /), &
      im1, im2, imm1, imm2, EL_keo, AZ_keo, EL_s1c, AZ_s1c)
    et_time = omp_get_wtime()
    call denorm_args(pminn,pmind)    
    
    write (*,"(A, F6.2,F6.2,F6.1,F8.1,F6.1,A,F10.3)") "p_min  :", pmind(1)*180./pi, pmind(2)*180./pi, &
    pmind(3), pmind(4), pmind(5), " - ", res_min
    
    call sphere_rout((/ crop_pars_1(4),crop_pars_1(3) /), &
    pmind, pk1, EL_keo, AZ_keo, &    
           imm1)
    
    call sphere_rout((/ crop_pars_2(4),crop_pars_2(3) /), &
    pmind, pk2, EL_s1c, AZ_s1c, &    
           imm2)
    
    open(unit = 9, file = 'model_keo_final.bin', access = "STREAM", form = "unformatted", status="replace")
    write(9) transpose(imm1)
    close(9)
    
    open(unit = 9, file = 'model_s1c_final.bin', access = "STREAM", form = "unformatted", status="replace")
    write(9) transpose(imm2)
    close(9)
    
    write (*,"(A,I6,A,I3,A,I3)") "N = ",icount, " R = ",numres, " S = ",ifault
    write(6,"(A,f10.5,A)") 'Elapsed time - ', et_time-st_time, " seconds"
    deallocate (AZ_keo, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate AZ_keo***"    
    deallocate (EL_keo, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate EL_keo***"    
    deallocate (AZ_s1c, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate AZ_s1c***"
    deallocate (EL_s1c, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate EL_s1c***"
        
    open(unit = 9, file = 'glow_keo.bin', access = "STREAM", form = "unformatted", status="replace")
    write(9) transpose(im1)
    close(9)
    
    open(unit = 9, file = 'glow_s1c.bin', access = "STREAM", form = "unformatted", status="replace")
    write(9) transpose(im2)
    close(9)
    
    deallocate (im1, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate im1***"
    deallocate (im2, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate im2***"
    deallocate (imm1, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate imm1***"
    deallocate (imm2, stat = DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop "*** Can not deallocate imm2***"
    
end

subroutine sphere_rout(naxes, p, pk, elm, azm, &    
           imm)
    implicit none
    !integer :: OMP_GET_THREAD_NUM 
    integer, intent(in) :: naxes(2)
    real*8, dimension(5), intent(in) :: p
    real*8, dimension(3), intent(in) :: pk
    real*8, dimension( naxes(1), naxes(2) ), intent(in) :: elm, azm
    real*8, dimension( naxes(1), naxes(2) ), intent(inout) :: imm
    real*8 :: Ac, Bc, Cc, El, Az
    real*8,  parameter :: pi = 3.1415926535897932_8
    !real*8,  parameter :: pi  = 4 * atan (1.0_8)
    real*8, parameter :: a=6378.137d0, b=6356.752314245d0, e=0.081819190842965558d0
    integer :: i, j
    !integer :: n 
    
    !!$OMP PARALLEL SHARED(elm, azm, imm) PRIVATE(i,j,El,Az,AC,Bc,Cc)
    !!$OMP PARALLEL SHARED(imm) PRIVATE(El,Az,AC,Bc,Cc)
    !$OMP PARALLEL 
    !$OMP DO COLLAPSE(2)
    do i=1, naxes(2), 1
        do j=1, naxes(1),1
            !n= OMP_GET_THREAD_NUM()
            El=elm(i,j)
            Az=azm(i,j)
            
            Ac = ((-((-sin(p(1))*sin(p(2))*(cos(pk(1))*sin(pk(2))*sin(El)+cos&
            &(pk(2))*sin(Az)*cos(El)-sin(pk(1))*sin(pk(2))*cos(Az)*cos(El)))+c&
            &os(p(2))*(cos(pk(1))*cos(pk(2))*sin(El)-sin(pk(2))*sin(Az)*cos(El&
            &)-sin(pk(1))*cos(pk(2))*cos(Az)*cos(El))+cos(p(1))*sin(p(2))*(sin&
            &(pk(1))*sin(El)+cos(pk(1))*cos(Az)*cos(El)))**2)-((-sin(p(1))*cos&
            &(p(2))*(cos(pk(1))*sin(pk(2))*sin(El)+cos(pk(2))*sin(Az)*cos(El)-&
            &sin(pk(1))*sin(pk(2))*cos(Az)*cos(El)))-sin(p(2))*(cos(pk(1))*cos&
            &(pk(2))*sin(El)-sin(pk(2))*sin(Az)*cos(El)-sin(pk(1))*cos(pk(2))*&
            &cos(Az)*cos(El))+cos(p(1))*cos(p(2))*(sin(pk(1))*sin(El)+cos(pk(1&
            &))*cos(Az)*cos(El)))**2-(cos(p(1))*(cos(pk(1))*sin(pk(2))*sin(El)&
            &+cos(pk(2))*sin(Az)*cos(El)-sin(pk(1))*sin(pk(2))*cos(Az)*cos(El)&
            &)+sin(p(1))*(sin(pk(1))*sin(El)+cos(pk(1))*cos(Az)*cos(El)))**2)/&
            &p(5)**2

            Bc = ((-2*((-sin(p(1))*sin(p(2))*(cos(pk(1))*sin(pk(2))*sin(El)+c&
            &os(pk(2))*sin(Az)*cos(El)-sin(pk(1))*sin(pk(2))*cos(Az)*cos(El)))&
            &+cos(p(2))*(cos(pk(1))*cos(pk(2))*sin(El)-sin(pk(2))*sin(Az)*cos(&
            &El)-sin(pk(1))*cos(pk(2))*cos(Az)*cos(El))+cos(p(1))*sin(p(2))*(s&
            &in(pk(1))*sin(El)+cos(pk(1))*cos(Az)*cos(El)))*(cos(p(1))*sin(p(2&
            &))*(sin(pk(1))*(b**2/(a*sqrt(1-sin(pk(1))**2*e**2))+pk(3))-sin(p(&
            &1))*(b**2/(a*sqrt(1-sin(p(1))**2*e**2))+p(3)))-sin(p(1))*sin(p(2)&
            &)*(cos(pk(1))*sin(pk(2))*(a/sqrt(1-sin(pk(1))**2*e**2)+pk(3))-cos&
            &(p(1))*sin(p(2))*(a/sqrt(1-sin(p(1))**2*e**2)+p(3)))+cos(p(2))*(c&
            &os(pk(1))*cos(pk(2))*(a/sqrt(1-sin(pk(1))**2*e**2)+pk(3))-cos(p(1&
            &))*cos(p(2))*(a/sqrt(1-sin(p(1))**2*e**2)+p(3)))))-2*((-sin(p(1))&
            &*cos(p(2))*(cos(pk(1))*sin(pk(2))*sin(El)+cos(pk(2))*sin(Az)*cos(&
            &El)-sin(pk(1))*sin(pk(2))*cos(Az)*cos(El)))-sin(p(2))*(cos(pk(1))&
            &*cos(pk(2))*sin(El)-sin(pk(2))*sin(Az)*cos(El)-sin(pk(1))*cos(pk(&
            &2))*cos(Az)*cos(El))+cos(p(1))*cos(p(2))*(sin(pk(1))*sin(El)+cos(&
            &pk(1))*cos(Az)*cos(El)))*(cos(p(1))*cos(p(2))*(sin(pk(1))*(b**2/(&
            &a*sqrt(1-sin(pk(1))**2*e**2))+pk(3))-sin(p(1))*(b**2/(a*sqrt(1-si&
            &n(p(1))**2*e**2))+p(3)))-sin(p(1))*cos(p(2))*(cos(pk(1))*sin(pk(2&
            &))*(a/sqrt(1-sin(pk(1))**2*e**2)+pk(3))-cos(p(1))*sin(p(2))*(a/sq&
            &rt(1-sin(p(1))**2*e**2)+p(3)))-sin(p(2))*(cos(pk(1))*cos(pk(2))*(&
            &a/sqrt(1-sin(pk(1))**2*e**2)+pk(3))-cos(p(1))*cos(p(2))*(a/sqrt(1&
            &-sin(p(1))**2*e**2)+p(3))))-2*(cos(p(1))*(cos(pk(1))*sin(pk(2))*s&
            &in(El)+cos(pk(2))*sin(Az)*cos(El)-sin(pk(1))*sin(pk(2))*cos(Az)*c&
            &os(El))+sin(p(1))*(sin(pk(1))*sin(El)+cos(pk(1))*cos(Az)*cos(El))&
            &)*(sin(p(1))*(sin(pk(1))*(b**2/(a*sqrt(1-sin(pk(1))**2*e**2))+pk(&
            &3))-sin(p(1))*(b**2/(a*sqrt(1-sin(p(1))**2*e**2))+p(3)))+cos(p(1)&
            &)*(cos(pk(1))*sin(pk(2))*(a/sqrt(1-sin(pk(1))**2*e**2)+pk(3))-cos&
            &(p(1))*sin(p(2))*(a/sqrt(1-sin(p(1))**2*e**2)+p(3)))))/p(5)**2

            Cc = ((-(cos(p(1))*sin(p(2))*(sin(pk(1))*(b**2/(a*sqrt(1-sin(pk(1&
            &))**2*e**2))+pk(3))-sin(p(1))*(b**2/(a*sqrt(1-sin(p(1))**2*e**2))&
            &+p(3)))-sin(p(1))*sin(p(2))*(cos(pk(1))*sin(pk(2))*(a/sqrt(1-sin(&
            &pk(1))**2*e**2)+pk(3))-cos(p(1))*sin(p(2))*(a/sqrt(1-sin(p(1))**2&
            &*e**2)+p(3)))+cos(p(2))*(cos(pk(1))*cos(pk(2))*(a/sqrt(1-sin(pk(1&
            &))**2*e**2)+pk(3))-cos(p(1))*cos(p(2))*(a/sqrt(1-sin(p(1))**2*e**&
            &2)+p(3))))**2)-(cos(p(1))*cos(p(2))*(sin(pk(1))*(b**2/(a*sqrt(1-s&
            &in(pk(1))**2*e**2))+pk(3))-sin(p(1))*(b**2/(a*sqrt(1-sin(p(1))**2&
            &*e**2))+p(3)))-sin(p(1))*cos(p(2))*(cos(pk(1))*sin(pk(2))*(a/sqrt&
            &(1-sin(pk(1))**2*e**2)+pk(3))-cos(p(1))*sin(p(2))*(a/sqrt(1-sin(p&
            &(1))**2*e**2)+p(3)))-sin(p(2))*(cos(pk(1))*cos(pk(2))*(a/sqrt(1-s&
            &in(pk(1))**2*e**2)+pk(3))-cos(p(1))*cos(p(2))*(a/sqrt(1-sin(p(1))&
            &**2*e**2)+p(3))))**2-(sin(p(1))*(sin(pk(1))*(b**2/(a*sqrt(1-sin(p&
            &k(1))**2*e**2))+pk(3))-sin(p(1))*(b**2/(a*sqrt(1-sin(p(1))**2*e**&
            &2))+p(3)))+cos(p(1))*(cos(pk(1))*sin(pk(2))*(a/sqrt(1-sin(pk(1))*&
            &*2*e**2)+pk(3))-cos(p(1))*sin(p(2))*(a/sqrt(1-sin(p(1))**2*e**2)+&
            &p(3))))**2)/p(5)**2
                        
            imm(i,j) = (p(4)*5.15/10**4)*(sqrt(pi)*exp(((4*Ac*Cc-Bc**2)/Ac)/4.0d+0))/sqrt(-Ac)
            
            !if (i==144) then
            !    if (j==144) then
            !        write(*,*) "A=",Ac, "B=", Bc, "C=", Cc
            !        write(*,*)
            !        "I=",(sqrt(pi)*exp(((4*Ac*Cc-Bc**2)/Ac)/4.0d+0))/sqrt(-Ac)
            !        write(*,*) "El=",El, "Az=", Az
            !        write(*,*) "p1=", p(1), "p2=", p(2), "p3=", p(3)
            !        write(*,*) "p4=", p(4), "p5=", p(5)
            !        write(*,*) "pk1=", pk(1), "pk2=", pk(2), "pk3=", pk(3)
            !    end if
            !end if
            
        enddo
    enddo
    !$OMP END PARALLEL
end

subroutine get_images_naxes(im_fn, &
                      naxes1, naxes2)
    implicit none
    character(len=300), intent(in) :: im_fn
    integer, intent(out) :: naxes1(2), naxes2(2)
    integer st, un, rw, blocksize, nfound, hdutype
   
    st=0
    rw=0
    call ftgiou(un, &
                st)
    call ftopen(un,im_fn,rw, &
                blocksize,st)
    
    call ftgknj(un,'NAXIS',1,2, &
                naxes1,nfound,st)
    if (nfound .ne. 2) then
        print *,'READIMAGE failed to read the NAXISn keywords.'
        return
    end if    
    
    call ftmahd(un,2, &
                hdutype,st)
    call ftgknj(un,'NAXIS',1,2, &
                naxes2, nfound, st)
    if (nfound .ne. 2) then
        print *,'READIMAGE failed to read the NAXISn keywords.'
        return
    end if
    
    call ftclos(un, st)
    call ftfiou(un, st)
    
    !write (*, *) naxes1
    !write (*, *) naxes2
end

subroutine get_images_data(im_fn, naxes1, naxes2, crop_pars_1, crop_pars_2, & 
            im1, im2)
           
    implicit none
    character(len=300), intent(in) :: im_fn
    integer, intent(in) :: naxes1(2), naxes2(2)
    integer, dimension(4), intent(in) :: crop_pars_1, crop_pars_2
    real*8, dimension( crop_pars_1(4),crop_pars_1(3) ), intent(out) :: im1
    real*8, dimension( crop_pars_2(4),crop_pars_2(3) ), intent(out) :: im2
    
    real*8, dimension(:,:), allocatable :: im1_temp, im2_temp
    integer :: st, un, rw, blocksize, hdutype
    integer :: group, nullval, AllocateStatus, DeAllocateStatus
    logical :: anynull
    real*8 :: datamin, datamax    
    
    st=0
    rw=0
    group=1
    nullval=0
    
    allocate ( im1_temp(naxes1(1),naxes1(2)), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for im1_temp***"
    allocate ( im2_temp(naxes2(1),naxes2(2)), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory for im2_temp***" 
    !write(*,*) "test 250"    

    call ftgiou(un,st)
    call ftopen(un,im_fn,rw,blocksize,st)
    !write(*,*) "test 254", un, group, nullval, naxes1(1), naxes1(1), naxes1(2)

    call ftg2dd(un,group,nullval,naxes1(1),naxes1(1),naxes1(2),im1_temp,anynull,st)
    !write(*,*) "test 256"
    im1=transpose(im1_temp(1+crop_pars_1(1):crop_pars_1(1)+crop_pars_1(3), &
                 1+crop_pars_1(2):crop_pars_1(2)+crop_pars_1(4)))
    !write(*,*) "test 259"
    call ftmahd(un,2,hdutype,st)
    !write(*,*) "test 261"
    call ftg2dd(un,group,nullval,naxes2(1),naxes2(1),naxes2(2),im2_temp,anynull,st)
    !write(*,*) "test 263"
    im2=transpose(im2_temp(1+crop_pars_2(1):crop_pars_2(1)+crop_pars_2(3), &
                 1+crop_pars_2(2):crop_pars_2(2)+crop_pars_2(4)))
    !write(*,*) "test 266"
    datamin=minval(im1_temp)
    datamax=maxval(im1_temp)
    print *,'Min and max values in the image are:',datamin,datamax
    datamin=minval(im2_temp)
    datamax=maxval(im2_temp)    
    print *,'Min and max values in the image are:',datamin,datamax
    
    call ftclos(un, st)
    call ftfiou(un, st)

    deallocate (im1_temp, stat = DeAllocateStatus) 
    deallocate (im2_temp, stat = DeAllocateStatus)

end

subroutine get_elaz_matrices(solve_pars_fn1, solve_pars_fn2, crop_pars_1, crop_pars_2, &
                             EL_keo, AZ_keo, EL_s1c, AZ_s1c)
    implicit none
    character(len=300), intent(in) :: solve_pars_fn1, solve_pars_fn2
    integer, dimension(4), intent(in) :: crop_pars_1, crop_pars_2
    real*8, dimension( crop_pars_1(4),crop_pars_1(3) ), intent(out) :: EL_keo
    real*8, dimension( crop_pars_1(4),crop_pars_1(3) ), intent(out) :: AZ_keo
    real*8, dimension( crop_pars_2(4),crop_pars_2(3) ), intent(out) :: EL_s1c
    real*8, dimension( crop_pars_2(4),crop_pars_2(3) ), intent(out) :: AZ_s1c
    
    real*8 :: er_1, az0_1, alt0_1, a_1(3), b_1(3), c_1(3), d_1(3)
    real*8 :: er_2, az0_2, alt0_2, a_2(3), b_2(3), c_2(3), d_2(3)    
    integer :: ii, jj
      
    call get_solve_pars(solve_pars_fn1, er_1, az0_1, alt0_1, a_1, b_1, c_1, d_1)      
    call get_solve_pars(solve_pars_fn2, er_2, az0_2, alt0_2, a_2, b_2, c_2, d_2)
    
    do ii=1, crop_pars_1(4), 1
        do jj=1, crop_pars_1(3), 1
            call get_arc_pix2hor(real(jj+crop_pars_1(1),8), real(ii+crop_pars_1(2),8), &
            az0_1, alt0_1, a_1, b_1, AZ_keo(ii,jj), EL_keo(ii,jj))
        enddo
    enddo
    
    do ii=1, crop_pars_2(4), 1
        do jj=1, crop_pars_2(3), 1
            call get_tan_pix2hor(real(jj+crop_pars_2(1),8), real(ii+crop_pars_2(2),8), &
            az0_2, alt0_2, a_2, b_2, AZ_s1c(ii,jj), EL_s1c(ii,jj))
        enddo
    enddo
        
end

subroutine get_tan_pix2hor(x, y, az0, alt0, a, b, az, alt)
    implicit none
    real*8, intent(in) :: x, y, az0, alt0, a(3), b(3)
    real*8, intent(out) :: az, alt
    real*8 :: Xs, Ys, Azi, R, Alti
    real*8,  parameter :: PI  = 3.1415926535897932_8
        
    Xs = a(1) + a(2) * x + a(3) * y
    Ys = b(1) + b(2) * x + b(3) * y
    
    Azi = atan2(Xs, -Ys)
    R = sqrt(Xs**2 + Ys**2)    
    Alti = atan(180/PI/R)    
    az=az0+atan2(- cos(Alti) * sin(Azi), sin(Alti) * cos(alt0) - cos(Alti) * sin(alt0) *cos(Azi) )    
    alt=asin( sin(Alti) * sin(alt0) + cos(Alti) * cos(alt0) * cos(Azi) )
end

subroutine get_arc_pix2hor(x, y, az0, alt0, a, b, az, alt)    
    implicit none
    real*8, intent(in) :: x, y, az0, alt0, a(3), b(3)
    real*8, intent(out) :: az, alt
    real*8 :: Xs, Ys, Azi, R, Alti
    real*8,  parameter :: PI  = 3.1415926535897932_8    
    
    Xs = a(1) + a(2) * x + a(3) * y
    Ys = b(1) + b(2) * x + b(3) * y
    
    Azi = atan2(Xs, Ys)
    R = sqrt(Xs**2 + Ys**2)    
    Alti = (90 - R) * PI / 180    
    az=az0+atan2(- cos(Alti) * sin(Azi), sin(Alti) * cos(alt0) - cos(Alti) * sin(alt0) *cos(Azi) )
    alt=asin( sin(Alti) * sin(alt0) + cos(Alti) * cos(alt0) * cos(Azi) )
end

subroutine get_solve_pars(fn, er, az0, alt0, a, b, c, d)
    implicit none
    character(len=80), intent(in) :: fn
    real*8, intent(out) :: er, az0, alt0, a(3), b(3), c(3), d(3)
    character(len=300) :: line
    integer :: ios
    integer, parameter :: read_unit = 99
    
    open(unit=read_unit, file=fn, iostat=ios)
    if ( ios /= 0 ) stop "Error opening file data.dat"    
    read(read_unit, '(A)', iostat=ios) line
    !write(*,*) line
    read(read_unit, *, iostat=ios) er, az0, alt0, a(1), a(2), a(3), b(1), b(2), b(3), c(1), c(2), c(3), d(1), d(2), d(3)    
    
end
