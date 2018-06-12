!     Program:         Modified Cox method with ORC for measurement error
!     Written by:      Xin Zhou
!     Created:         Dec. 26, 2015
!     Last modified:   Jan. 7, 2015

program main
    implicit none
    
    integer :: narg          ! number of command line argument
    character(len=1024) :: in_filename      ! in.txt configuration file
    character(len=1024) :: result_filename  ! filename of result file
    integer, parameter :: RESULTS = 1     ! file handle
    
    double precision :: power
    
    integer :: II, JJ, KK, NI   ! 
    double precision :: mu, beta, tau2
    double precision, dimension (:), allocatable :: gamma   ! gamma(JJ), time effects with gamma(1) = 0
    double precision :: p0totalchange, p0stepchange
    double precision :: p11, rho0   ! p11: the risk for intervention group at time 1
    double precision, dimension (:), allocatable :: p0   ! p0(JJ): the risks for standard of care group at time 1,2,...,JJ

    ! components of m(theta) and M(theta)
    ! mincomp is a JJ+2 vector of 0 and 1's, representing the weights of gamma(1),...,gamma(JJ),mu,beta 
    ! maxcomp is a JJ+2 vector of 0 and 1's, representing the weights of gamma(1),...,gamma(JJ),mu,beta
    integer, dimension (:), allocatable :: mincomp, maxcomp  
    double precision, external :: LinearPower_time
    double precision, external :: LinearPower_notime
    integer :: convergence
    double precision, parameter :: TOLERANCE = 1e-5      ! tolerance
    
    integer :: i, j, k, l
    
    integer :: GQ
    double precision, dimension (:), allocatable :: GQX, GQW   ! Gaussian Legendre integration
    double precision :: a, b        ! lower and upper limits of GL integration

    double precision :: temp
    
    character*1 tab
    
    tab = char(9)
    GQ = 100  ! number of quadrature points
    allocate (GQX(GQ))
    allocate (GQW(GQ))
    
    ! beta0 = 0.05d0
    ! beta1 = 1.0d0
    ! rho = 0.6d0
    ! tau2
    
    ! II = 24
    ! JJs = (/ 3, 5, 7, 9/)
    ! KK = 300

    ! read the arguments
    narg = command_argument_count()
    if (narg /= 2) then
        print *, 'There are two arguments!'
        print *, 'swdpower in.txt out.txt'
        stop
    end if
    ! I didn't check the command-line arguments
    ! make sure they are correct when run the program.
    ! especially, sim_start is less than sim_end
    call get_command_argument(1,in_filename)
    call get_command_argument(2,result_filename)
    
    call ReadInformation(in_filename, II, JJ, KK, mu, beta, rho0, p0totalchange)
    
    
!    double precision :: res
!    call HERZO(GQ,GQX,GQW)
!! test gaussian quadrature
!    res = 0.0d0
!    do i=1,GQ
!        res = res + GQW(i)*GQX(i)*GQX(i)
!    enddo
!    res = res*2/sqrt(PI)
    
        allocate (p0(JJ))
        allocate (gamma(JJ))
        allocate (mincomp(JJ+2),maxcomp(JJ+2))
        
       
        p0 = mu
        p11 = p0(1) + beta
        p0stepchange = p0totalchange/(JJ-1)
        do j=2,JJ
            p0(j) = p0(j-1) + p0stepchange
        end do
        call computeparameter(JJ, mu, beta, gamma, tau2, p0, p11, rho0)
        
        if (p0totalchange > TOLERANCE .or. p0totalchange < -TOLERANCE) then
            a = 100 
            b = -100
            do j=1,JJ
                temp = mu+gamma(j)
                if (temp < a) then
                    a = temp
                    mincomp=0
                    mincomp(JJ+1) = 1
                    mincomp(j) = 1
                end if
                if (temp > b) then
                    b = temp
                    maxcomp=0
                    maxcomp(JJ+1) = 1
                    maxcomp(j) = 1
                end if
                temp = mu+beta+gamma(j)
                if (temp < a) then
                    a = temp
                    mincomp=0
                    mincomp(JJ+1) = 1
                    mincomp(JJ+2) = 1
                    mincomp(j) = 1
                end if
                if (temp > b) then
                    b = temp
                    maxcomp=0
                    maxcomp(JJ+1) = 1
                    maxcomp(JJ+2) = 1
                    maxcomp(j) = 1
                end if
            end do
            a = -a
            b = 1-b
            call legendre_handle (GQ, a, b, GQX, GQW)
            ! Gaussian Legendre will not take two limits, a and b
            power = LinearPower_time(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW)
        else
            if (beta>0) then
                a = - mu
                b = 1 - mu - beta
            else
                a = - mu - beta
                b = 1 - mu
            end if
            ! a = a / dsqrt(2.0d0*tau2)
            ! b = b / dsqrt(2.0d0*tau2)
            call legendre_handle (GQ, a, b, GQX, GQW)
            ! Gaussian Legendre will not take two limits, a and b
            power = LinearPower_notime(mu, beta, tau2, II, JJ, KK, a, b, GQ, GQX, GQW)
        end if

        open(unit=RESULTS,file=result_filename,status='replace')
        write (RESULTS, '(a, i7)'), 'I = ', II
        write (RESULTS, '(a, i7)'), 'J = ', JJ
        write (RESULTS, '(a, i7)'), 'K = ', KK
        write (RESULTS, '(a, f10.3)'), 'Baseline effect : ', mu
        write (RESULTS, '(a, f10.3)'), 'Treatment effect : ', beta
        write (RESULTS, '(a, f10.3)'), 'Time effect : ',  p0totalchange
        write (RESULTS, '(a, f10.3)'), 'ICC : ',  rho0
        write (RESULTS, '(a, f10.3)'), 'POWER = ',  power
        deallocate(p0, gamma)
        deallocate (mincomp, maxcomp)

    close(RESULTS)
    deallocate(GQX, GQW)
end program main

subroutine ReadInformation(in_filename, II, JJ, KK, mu, beta, rho0, p0totalchange)
    USE String_Utility
    implicit none
    
    character(len=1024) :: in_filename      ! in.txt configuration file
    integer :: II, JJ, KK   ! 
    double precision :: mu, beta
    double precision :: p0totalchange
    double precision :: rho0
    
    integer, parameter :: INFORMATION = 1     ! in.txt
    
    character(len=1024) :: buffer
    character(len=20) :: keyin
    character(len=1024) :: keyvalue
    
    character(len=*), parameter :: KEYII = 'i'
    character(len=*), parameter :: KEYJJ = 'j'
    character(len=*), parameter :: KEYKK = 'k'
    character(len=*), parameter :: KEYMU = 'mu'
    character(len=*), parameter :: KEYBETA = 'beta'
    character(len=*), parameter :: KEYDELTA = 'delta'
    character(len=*), parameter :: KEYICC = 'icc'
    
    integer :: possign    ! pos of = in the information string of in.txt
    
    integer :: ios                            ! file action status
    
!    integer, dimension (:), allocatable :: SelectedIndicators  ! p*1 vector storing a set of indicators
!    double precision, dimension (:), allocatable :: SelectedKappas   ! p*1 vector storing kappas
!    double precision, dimension (:), allocatable :: SelectedKappasCV   ! p*1 vector storing kappas
    
    !*****************************************************
    ! read information
    
    ! default values
    II = 12
    JJ = 3
    KK = 90
    mu = 0.0d0
    beta = 0.0d0
    p0totalchange = 0.0d0
    rho0 = 0.0d0
    
    open (unit=INFORMATION, file=in_filename, status='old', iostat=ios, action='read', position='rewind')
	if ( ios /= 0 ) then
		print *, 'Could not open ', in_filename,'!'
		stop
	endif
    
    ! read(INFORMATION, '(A)', iostat=ios) buffer
    do while (.TRUE.)
        read(INFORMATION, '(A)', iostat=ios) buffer
        if (ios /= 0) exit
        buffer = adjustl(buffer)
        if (buffer(1:1)==' ') cycle   ! empty line
        if (buffer(1:1)=='#') cycle   ! comments
        possign = scan(buffer,'=')    ! find '='
        if (possign<=1) cycle
        keyin = buffer(1:(possign-1))
        keyin = StrLowCase(keyin)
        keyvalue = buffer((possign+1):len_trim(buffer))
        
        if (keyin(1:len_trim(keyin))==KEYII) then
            read (keyvalue, *) II
            cycle
        end if
        if (keyin(1:len_trim(keyin))==KEYJJ) then
            read (keyvalue, *) JJ
            cycle
        end if
        if (keyin(1:len_trim(keyin))==KEYKK) then
            read (keyvalue, *) KK
            cycle
        end if
        if (keyin(1:len_trim(keyin))==KEYMU) then
            read (keyvalue, *) mu
            cycle
        end if
        if (keyin(1:len_trim(keyin))==KEYBETA) then
            read (keyvalue, *) beta
            cycle
        end if
        if (keyin(1:len_trim(keyin))==KEYDELTA) then
            read (keyvalue, *) p0totalchange
            cycle
        end if
        if (keyin(1:len_trim(keyin))==KEYICC) then
            read (keyvalue, *) rho0
            cycle
        end if
        
    end do    
    
    if (ios > 0 ) then
		print *, 'Reading error in ', in_filename,'!'
        close(INFORMATION, iostat=ios)
		stop
    endif
    close(INFORMATION, iostat=ios)
end subroutine ReadInformation
