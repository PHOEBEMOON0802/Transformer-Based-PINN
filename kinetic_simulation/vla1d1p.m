%                                                                                      !
%      1D1V 非相对论Vlasov-Possion方程代码                                                !
%                                                                                       !
% 
%       计算1D1V-Vlasov-Possion方程
%           df/dt + v df/dx + qE/m df/dv = 0
%          dE/dx = q rho - q<densm>
% 单位:
%       时间： 1/wp   空间： lambda(德拜长度)   速度: 电子热速度
%
% *计算结果：
%     1. 全局诊断信息存储在: xxxx.diag文件中，通过修改diagnostic_header, diagnostic函数可增加输出
%     2. 分布函数信息存储在: xxxx.distf文件中，每次输出 4 + 8+(N+2)*(2*M+1)*8+4
%                                         其中第一个和最后一个为半个字节
%                                         共输出 nt/nplot + 2 次
%
%     3. 电场信息        : xxxx.field文件中，每次输出 4 + 8 + (N+2)*8 + 4个字节
%                                         共输出 nt/nplot + 2 次
%   ---------------------------------------------------------------------------------------------
%   计算结果：
%   ---------------------------------------------------------------------------------------------
%  time density momentum energy_f energy_em energy L2 entropy log|E| EF(MID-1) EF(MID) EF(MID+1)
%       0.000    0.1000000000E+01    0.5476612433E-17    0.5000000000E+00    0.9999999839E-06    0.5000010000E+00   -0.3544909474E+01    0.1783090435E+02   -0.6561181697E+01   -0.9813534786E-04   -0.1853943363E-18    0.9813534786E-04
%       0.200    0.1000000000E+01    0.1052587754E-16    0.5000000125E+00    0.9745777251E-06    0.5000009871E+00   -0.3544909474E+01    0.1783090435E+02   -0.6574057191E+01   -0.9688086883E-04   -0.1795589146E-15    0.9688086884E-04
%       0.400    0.1000000000E+01   -0.1012526272E-16    0.5000001084E+00    0.8572735390E-06    0.5000009656E+00   -0.3544909474E+01    0.1783090435E+02   -0.6638180803E+01   -0.9086784033E-04   -0.4132742571E-15    0.9086784033E-04
%       0.600    0.1000000000E+01   -0.2400605573E-16    0.5000002797E+00    0.6723899417E-06    0.5000009521E+00   -0.3544909474E+01    0.1783090435E+02   -0.6759640107E+01   -0.8048127070E-04   -0.1557248898E-15    0.8048127070E-04
%       0.800    0.1000000000E+01    0.4886922369E-16    0.5000004907E+00    0.4580962905E-06    0.5000009488E+00   -0.3544909474E+01    0.1783090435E+02   -0.6951519627E+01   -0.6643598714E-04   -0.6895176416E-16    0.6643598714E-04
%       1.000    0.1000000000E+01    0.6711257233E-16    0.5000006994E+00    0.2562324644E-06    0.5000009556E+00   -0.3544909474E+01    0.1783090435E+02   -0.7242016780E+01   -0.4969138592E-04    0.1633943178E-14    0.4969138592E-04
%       1.200    0.1000000000E+01    0.3643682895E-16    0.5000008680E+00    0.1021613091E-06    0.5000009702E+00   -0.3544909474E+01    0.1783090435E+02   -0.7701782815E+01   -0.3137828874E-04    0.6093349617E-15    0.3137828874E-04
%       1.400    0.1000000000E+01    0.5120314794E-16    0.5000009715E+00    0.1677432495E-07    0.5000009883E+00   -0.3544909474E+01    0.1783090435E+02   -0.8605134608E+01   -0.1271307552E-04    0.1494707227E-14    0.1271307553E-04
%       1.600    0.1000000000E+01    0.8341047759E-16    0.5000010025E+00    0.2686085288E-08    0.5000010052E+00   -0.3544909474E+01    0.1783090435E+02   -0.9521016903E+01    0.5093938252E-05    0.7675554223E-15   -0.5093938250E-05
%       1.800    0.1000000000E+01    0.6516024541E-16    0.5000009716E+00    0.4540183200E-07    0.5000010170E+00   -0.3544909474E+01    0.1783090435E+02   -0.8107283100E+01    0.2092835246E-04    0.6996213682E-15   -0.2092835246E-04
%   ---------------------------------------------------------------------------------------------
!> vlasov全局数据模块
module nvlapos_mod
    integer, parameter :: N = 64             ! x 方向格点数 (周期系统 x(-1:N))
    integer, parameter :: M = 63             ! v 方向个点数 ( v(-M:M) )
    double precision   :: xl                  ! 模拟盒子长度 (单位: 德拜长度)
    double precision   :: vmax                ! 最大速度 (单位: 电子热速度)
    double precision   :: mass, charge        ! 质量和电荷大小
    !!----
    !!    全局数据： 位置、速度、分布函数等 (系统分配)
    double precision   :: v(-M:M), x(-1:N)    ! 位置和速度 (x,v) {x(-1),x(N) 为边界 }
    double precision   :: f(-1:N,-M:M)        ! 分布函数
    double precision   :: f2(-1:N,-M:M)       ! 分布函数2次导数 (3次样条插值)
    double precision   :: g(-1:N,-M:M)        ! 辅助内存
    double precision   :: EF(-1:N),RHO(-1:N)  ! 电场和密度
    double precision   :: FV(-M:M)            ! 速度分布 1/L <f> dx
    double precision   :: dx,dv,dt            ! 时间、空间步长大小
    !!    全局数据: 诊断信息
    double precision en_kin, en_em, mom, den,l2,entropy !动能，电场能，动量，密度，|f|^2, 熵
    double precision   :: time                ! 当前时间
    double precision, parameter :: pi   = 3.1415926535897931d0
contains
    !> 调试信息输出
!! @param   str     调试信息
!! @param   iparam  整型调试信息数组
!! @param   dparam  实型调试信息数组
subroutine debug(str,iparam, dparam)
   implicit none
   character(len=*) :: str
   integer,dimension(:),optional :: iparam
   double precision,dimension(:),optional :: dparam
   !!
   character(len=100) :: fs
   integer ni,nd
   ni = 0
   nd = 0

   if(present(iparam)) then
       ni = size(iparam)
   end if
   if(present(dparam)) then
      nd = size(dparam)
    end if
    if(ni .ne. 0 .and. nd .ne. 0) then
        write(fs,'(1X,A3,I1,A3,I1,A6)') '(A,',ni,'I8,',nd,'F12.8)'
        write(*,fs) str,iparam,dparam
     else
          if(ni .ne. 0) then
              write(fs,'(1X,A3,I1,A3)') '(A,',ni,'I8)'
              write(*,fs) str,iparam
          else
              if(nd .ne. 0) then
                 write(fs,'(1X,A3,I1,A6)') '(A20,',nd,'F12.8)'
                 write(*,fs) str,dparam
              else
                 write(*,'(1X,A)') str
              end if
          end if
       end if
end subroutine debug
end module nvlapos_mod

!> 初始设置格点和分布函数信息
!! -- 静电Landau-Damping模型 : f(x,v) = f(v)(1+alpha*cos(k*x))
subroutine initcon_landau
    use nvlapos_mod
    implicit none
    double precision :: rk0            ! 波数
    double precision :: alpha          ! 扰动大小
    double precision :: cte            ! 归一化因子
    !!
    double precision :: pert
    integer          :: i,j
    !! ----- 计算区间与分布函数信息 (S)
    rk0  = .5d0
    alpha = 0.001d0
    cte   = 1.d0/dsqrt(2.d0*pi)
	mass = 1.d0
	charge = -1.d0
    vmax = 10.0d0
    xl   = 2.d0*pi/rk0
    !! ----- 计算区间与分布函数信息 (E)
    dv   = vmax/dble(M)
    dx   = xl/dble(N)
    x = (/ (i * dx, i = -1,N) /)
    v = (/ (j * dv, j = -M,M) /)
    do i = 0, N-1
        pert = alpha*dcos(rk0*x(i))
        do j = -M,M
            f(i,j) = cte * dexp(-v(j)*v(j)/2.d0)*(1.d0+pert)
        end do
    end do
    !! ---分布函数速度边界处为 0
    f(:,-M) = 0.d0
    f(:,M)  = 0.d0
    !! ---空间周期性边界条件
    f(N, :) = f(0,:)
    f(-1,:) = f(N-1,:)
    return
end subroutine initcon_landau

!> 主程序
!!
program nvlapos
    use nvlapos_mod
    implicit none
    integer :: nt,nplot,ndiag,it

    call debug('-Vlasov1D1V')
    !诊断信息输出文件 (文本)
    open(unit=15,file='vla_sp1d_a=0.001_k=0.5_64x127.diag', form='formatted')
    !分布函数输出文件 (二进制)
    open(unit=16,file='vla_sp1d_a=0.001_k=0.5_64x127.distf',form='unformatted')
    !电场输出文件    （二进制)
    open(unit=17,file='vla_sp1d_a=0.001_k=0.5_64x127.field',form='unformatted')

    dt = 0.1D0           !< 时间步长 ( 单位: 1/wp)
    nt = 1000            !< 计算步数
    nplot = 10           !< 分布函数和电场信息输出间隔 Number of time steps between F consecutive outputs.
    ndiag =1             !< 诊断信息输出间隔
    time = 0.d0          !< 初始时间

    !! ---- 计算 ---
    call debug('---- Initalization')
    call  initcon_landau    !< 初始化
    call  efield3           !< 初始计算电场
    call diagnostic_header
    call  diagnostic        !< 初始输出诊断信息
    call  plotf             !< 初始输出分布函数信息
    call  plotef            !< 初始输出电场信息
	call debug('---- Main_Loop : Begin')
    it = 0
	call debug('---- ---- INFO: ',(/it/),(/den,en_kin+en_em,entropy/))
    do it=1,nt
        time = time+DT
        !! (*1/2*)  f(t+dt,x,v) = f(t, x - v dt/2, v)
        call spline_x
        call advection_x_semi
        f = g

        !! (*1/2*)  计算电场
        call efield3
        call debug('---- ---- info @ = : ',(/it/),(/den,en_kin+en_em,entropy/))
        !! 输出
        if (MOD(it,ndiag).eq.0) then
            call debug('---- ---- ---- diag @ = ',(/it/),(/time/))
            call diagnostic
        end if
        if (MOD(it,nplot).eq.0) then
            call debug('---- ---- ---- plot @ = ',(/it/),(/time/))
            call plotf
            call plotef
        end if
        !! (*1*)  f(t+dt,x,v) = f(t, x, v - q E dt )
        call spline_v
        call advection_v
        f = g

        !! (*1/2*)  f(t+dt,x,v) = f(t, x - v dt/2, v)
        call spline_x
        call advection_x_semi
        f = g

   end do
   call debug('---- Main_Loop : End')
   call efield3
   call diagnostic
   call plotf
   call plotef
   close(15)
   close(16)
   close(17)
   call debug('- Bye!')
   stop
end program

!> 输出分布函数信息
subroutine plotf
    use nvlapos_mod
    implicit none
    write(16) time,f
end subroutine plotf

!> 输出电场信息
subroutine plotef
    use nvlapos_mod
    implicit none
    write(17) time,EF
end subroutine plotef

subroutine diagnostic_header
    use nvlapos_mod
    write(15,'(12(1X,A))') "time","density","momentum","energy_f","energy_em","energy",&
                           "L2", "entropy", "log|E|","EF(MID-1)","EF(MID)","EF(MID+1)"
    return
end subroutine

!> 输出诊断信息
subroutine diagnostic
    use nvlapos_mod
    implicit none
    integer          :: i,j
    double precision :: ff
    den = 0.d0
    en_kin = 0.d0
    en_em  = 0.d0
    mom = 0.d0
    l2  = 0.d0
    entropy = 0.d0
    do i = 0, N-1
        do j = -M, M
            den     = den + f(i,j)
            en_kin  = en_kin + v(j) * v(j) * f(i,j)
            mom     = mom    + v(j) * f(i,j)
            ff      = dabs(f(i,j)) + 1.d-16
            l2      = l2 + ff*ff
            entropy = entropy + ff*dlog(ff)
        end do
        en_em    = en_em + EF(i)*EF(i)
    end do
    den = den * dv / dble(N)        								! 平均密度
    en_kin = 0.5d0*en_kin*dv/dble(N)								!    动能
    en_em  = 0.5d0*en_em/dble(N)									!    电场能
    mom = mom * dv /dble(N)											!    动量
    l2 = -l2*dx*dv                        						    !    L2
    entropy = -entropy*dx*dv                					    !    熵
    write(15,'(1x,f10.3, 11E20.10)') time, den,mom, en_kin, en_em,en_kin+en_em,&
                                     l2,entropy,dlog(dsqrt(en_em*2.d0)),ef(N/2-1),ef(N/2),ef(N/2+1)
    return
end subroutine diagnostic

!> 计算电场 (求解Possion方程）
!!
subroutine efield3
    use nvlapos_mod
    implicit none
    integer i,j
    double precision d(0:N-1),aa(0:N-1),bb(0:N-1),cc(0:N-1),phi(0:N-1),U(0:N-1)
    double precision densm,emm,dx2d12,tddx
    densm = 0.d0
    do i = 0, N-1
        rho(i) = 0.d0
        do j = -M,M
            rho(i) = rho(i) + f(i,j)
        end do
        rho(i) = rho(i)*dv      				! RHO : density at the grid point
        densm = densm + rho(i)
    end do
    densm = densm/dble(N)       				! DENSM :  mean density  1/L < rho(i) dx >
    do i = 0,N-1
        rho(i) = charge * rho(i) - charge * densm
    end do
    dx2d12 = dx**2.d0/12.d0
    d(0)   = (rho(N-1) + 10.d0*rho(0) + rho(1))*dx2d12
    d(N-1) = (rho(N-2) + 10.d0*rho(N-1)+rho(0))*dx2d12
    do i = 1,N-2
        d(i) = (rho(i-1) + 10.d0*rho(i) + rho(i+1))*dx2d12
    end do
    do i = 0,N-1
        aa(i) = -1.d0
        bb(i) = 2.d0
        cc(i) = -1.d0
    end do
    call ptridiag(N,aa,bb,cc,phi,d)
    tddx = 3.d0/dx
    do i=1,N-2
       d(i)=(phi(i+1)-phi(i-1))*tddx
    enddo
    d(0)=(phi(1)-phi(N-1))*tddx
    d(N-1)=(phi(0)-phi(N-2))*tddx
    do i = 0, N-1
        aa(i) = -1.d0
        bb(i) = -4.d0
        cc(i) = -1.d0
    end do
    call ptridiag(N,aa,bb,cc,U,d)							! solve for E
    emm = 0.d0
    do i = 0,N-1
        emm = emm + U(i)
    end do
    emm = emm/dble(N)
    do i = 0,N-1
        EF(i) = U(i) - emm
    end do
    EF(N) = EF(0)											! impose periodic condition
    EF(-1) = EF(N-1)
end subroutine

!> 计算f(x,v)对x的2次导数
subroutine spline_x()
    use nvlapos_mod
    implicit none
    integer j
    ! the last column is set to be equal to the first
    ! to allow interpolation beyond xmin+Nx*dx (up to xmax!)
    f(N,:) = f(0,:)
    f(-1,:) = f(N-1,:)
    do j = -M,M
        call cspline_periodic(n,f(0:N-1,j),dx,f2(0:N-1,j))
        f2(N,j)  = f2(0,j)
        f2(-1,j) = f2(N-1,j)
    end do
end subroutine spline_x

!> 返回x方向插值 (周期边界条件)
!! @param x_in  插值点
!! @param j     速度方向指标
!! @return      函数在插值点处值
double precision function splint_x(x_in,j)
    use nvlapos_mod
    implicit none
    double precision, intent(in) :: x_in
    integer, intent(in) :: j
    double precision :: xx,cspline_2
    external cspline_2
    xx = x_in
    do while(xx .lt. 0 .or. xx .ge. xl)
    	 if (xx.le.0) then
         	xx = xx + xl
         end if
       	if (xx.ge.xl) then
        	xx = xx - xl
       	end if
    end do
   !! 修改： 2012-4-29 防止-check bounds编译后运行数组上界溢出
    splint_x = cspline_2(n,f(0:N,j), dx, f2(0:N,j), xx)
end function splint_x


!! 计算f(t+dt/2,x,v)=f(t,x-v dt, v)
subroutine advection_x_semi()
    use nvlapos_mod
    implicit none
    double precision splint_x
    integer :: i,j
    do i=0,N-1
        do j = -M,M
           g(i,j) = splint_x(x(i)-DT*v(j)*0.5d0,j)
        end do
    end do
	do j = -M,M
		g(N,j)  = g(0,j)
		g(-1,j) = g(N-1,j)
	end do
end subroutine advection_x_semi

!> 计算v方向f(x,v)2次导数
subroutine spline_v()
     use nvlapos_mod
     implicit none
     integer i
     do i=1,N
        call cspline_natural(2*M+1,f(i,-M:M), dv, f2(i,-M:M))
     end do
end subroutine spline_v


!> 返回v方向插值 (周期边界条件)
!! @param v_in  插值点
!! @param i     x方向指标
!! @return      函数在插值点处值
double precision function splint_v(v_in,i)
     use nvlapos_mod
     implicit none
     double precision, intent(in) :: v_in
     double precision cspline_2
     external cspline_2
     integer, intent(in) :: i
      if (v_in .le. -vmax .or. v_in.ge.vmax) then
         splint_v = 0.d0
         return
      end if
      splint_v = cspline_2(2*M+1,f(i,-M:M), dv, f2(i,-M:M), v_in+vmax)
      return
end function splint_v

!> 计算f(t+dt,x,v) = f(t,x,v-q/m E dt)
subroutine advection_v()
      use nvlapos_mod
      implicit none
      double precision splint_v
      integer :: i,j
      do i=0,N-1
         do j=-M,M
                g(i,j) = splint_v(v(j)-charge/mass*EF(i)*dt,i)
         end do
      end do
	  do j = -M,M
	  	g(N,j)  = g(0,j)
		g(-1,j) = g(N-1,j)
	  end do
end subroutine advection_v

!! 给定数组y, 利用Compact_difference计算y的2阶导数 (其中边界处2阶导数设为0)
!! @param   n       数组长度
!! @param   y       数组
!! @param   dx      步长
!! @param   y2      返回2阶导数
!>------------------------------------------------------------------------------------------------------
subroutine cspline_natural(n,y,dx, y2)
    implicit none
    integer                :: n
    double precision, intent(in)   :: y(n), dx
    double precision, intent(out)  :: y2(n)
    !!
    double precision, dimension(n) :: y_prime, a,b,c
    double precision               :: idx2m6
    idx2m6 = 6.d0/dx/dx
    y_prime(2:n-1)  = (y(3:n)-2.d0*y(2:n-1)+y(1:n-2))*idx2m6
    y_prime(1) = 0.d0
    y_prime(n) = 0.d0
    a(1:n) = 1.d0
    c(1:n) = 1.d0
    b(1:n) = 4.d0
    call tridiag(n,a,b,c,y2,y_prime)
end subroutine cspline_natural


!! 给定数组y, 利用Compact_difference计算y的2阶导数 (周期边界条件)
!! @param   n       数组长度
!! @param   y       数组
!! @param   dx      步长
!! @param   y2      返回2阶导数
!>------------------------------------------------------------------------------------------------------
subroutine cspline_periodic(n,y,dx,y2)
    implicit none
    integer n
    double precision, intent(in)         :: y(n), dx
    double precision, intent(out)        :: y2(n)
    double precision, dimension(n) :: y_prime, a, b, c
    double precision idx2m6
    idx2m6 = 6.d0/dx/dx
    y_prime(2:n-1)  = (y(3:n)-2.d0*y(2:n-1)+y(1:n-2))*idx2m6
    y_prime(1) = (y(2)-2.d0*y(1) + y(n))*idx2m6
    y_prime(n) = (y(n-1)-2.d0*y(n) + y(1))*idx2m6
    a(1:n) = 1.d0
    c(1:n) = 1.d0
    b(1:n) = 4.d0
    call ptridiag(n,a,b,c,y2,y_prime)
  end subroutine cspline_periodic

!> 给定数组y和它的2阶导数，通过插值返回x点处值.
!! @param   n       数组长度
!! @param   y       数组
!! @param   dx      步长
!! @param   y2      2阶导数
!! @param   x       所求数据点
!! @return  f(x)处的值
!>------------------------------------------------------------------------------------------------------
function cspline_2(n,y,dx,y2,x)
   implicit none
   integer n
   double precision :: cspline_2
   !! 修改： 2012-4-29 防止-check bounds编译后运行数组上界溢出
   double precision, intent(in) :: y(n+1),dx,y2(n+1), x
   integer :: i
   double precision :: my_a, my_b, my_c, my_d
   i = floor(x/dx)+1
   my_b = x/dx-dble(i-1)
   my_a = 1.d0-my_b
   my_c = (my_a**3-my_a)*dx**2.d0/6.0d0
   my_d = (my_b**3-my_b)*dx**2.d0/6.0d0
   !! @FIXME: 当采用 -check bound编译时，代码严格检测数组边界，原fortran数组越界返回0的语义此时不成立
   !!         为防止这种因严格调试时的编译结果，加入检测
   !! 在计算速度时，对于E=0，出现在f(x,v-Edt)=f(x,v)时，可出现这种情况
   cspline_2 = y(i)*my_a + y(i+1)*my_b + y2(i)*my_c + y2(i+1)*my_d
end function cspline_2



!> 追赶法求解3对角线性方程组 Ax = d（ 对角元素大于非对角元素 )
!!@param    n   方程维度
!!@param    a   下三角元值 a(1:n)  (a(1)不使用)
!!@param    b   对角元    b(1:n)
!!@param    c   上三角元  c(1:n)   (c(n+1)不使用)
!!@param    x   输出方程解
!!@param    d   the r.h.s of the equations
!>------------------------------------------------------------------------------------------------------
subroutine tridiag(n,a,b,c,x,d)
    implicit none
    integer n
    double precision ::  a(n),b(n),c(n),d(n),x(n)
    !!
    integer :: i
    c(1) = c(1) / b(1)
    do i=2,n-1
       c(i) = c(i)/(b(i)-c(i-1)*a(i))
    end do
    d(1) = d(1)/b(1)
    do i=2,n
       d(i) = (d(i)-d(i-1)*a(i))/(b(i)-c(i-1)*a(i))
    end do
    x(n) = d(n)
    do i=n-1,1,-1
       x(i) = d(i) - c(i)*x(i+1)
    end do
end subroutine tridiag


!> 追赶法求解周期3对角线性方程组 Ax = d
!! -- 利用Sherman-Morrison公式进行分解
!! @param    n   方程维度
!! @param    a   下三角元值 a(1:n)
!! @param    b   对角元    b(1:n)
!! @param    c   上三角元
!! @param    x   输出方程解
!! @param    d   the r.h.s of the equations
!>------------------------------------------------------------------------------------------------------
subroutine ptridiag(n,a, b, c, x, d)
    implicit none
    integer n
    double precision a(n),b(n),c(n),x(n),d(n)
    !!
    double precision a1(n),b1(n),c1(n)
    double precision :: u(n), x2(n)
    double precision :: alpha, beta, gamma,factor
    alpha = a(1)
    beta =  c(n)
    gamma = -b(1)
    a1(:) = a(:); b1(:) = b(:); c1(:)=c(:)
    b1(1) = b1(1) - gamma
    b1(n) = b1(n) - alpha*beta/gamma
    call tridiag(n,a1,b1,c1,x,d)
    u = 0.d0 ; u(1) = gamma ; u(n) = beta
    b(1) = b(1) - gamma
    b(n) = b(n) - alpha*beta/gamma
    call tridiag(n,a,b,c,x2,u)
    factor = ( x(1)+ x(n)*alpha/gamma ) / (1.d0 + x2(1) + x2(n)*alpha/gamma)
    x = x - factor * x2
end subroutine ptridiag