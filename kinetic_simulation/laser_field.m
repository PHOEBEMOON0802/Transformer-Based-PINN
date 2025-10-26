function bf = laser_field(e0, theta, ww0, center, xfocus, tpulse,trise,tfall,tdelay)
% ==========================================================
%  laser field at the boundary
% ==========================================================
global vars_global
   
tp=tpulse/fwhm
sinthe1=sin(theta)
costhe1=cos(theta)
ycenter = center
if theta > pi/2 && theta <= pi*1.5
    x0 = x_focus
else
    x0 = -xfocus
end
ycenter = center
xR = pi*ww0**2/alengt        
w2 = ww0**2*(1+x0**2/xR**2)  
xigam2 = datan(x0/xR)        
phi0 = 0                     
    ymax = min(ycenter+ww0,yleng)
    ymin = max(ycenter-ww0,0.d0)
    call pulse_shape(ilas, tpulse, trise, tfall, tp, t, tamp, 0.d0, fwhm)
    if (ilas .lt. 10)then
        do j=0,nym
            ys = dble(j-1)*dy + ystart
            yss  =  ys + 0.5d0*dy
            r2s=(yss-ycenter)**2
            ey_xampu = e0*tamp*dexp(-r2s/w2)
            phiG = (t-x0)+xigam2-x0*r2s/xR/w2-phi0
            bzda(j) = ey_xampu*ww0/dsqrt(w2) *dsin(phiG)         % bz(i+1/2,j+1/2)
        end do
    else
        do j=0,nym
            ys = dble(j-1)*dy + ystart
            yss  =  ys+0.5d0*dy
            r2s=(yss-ycenter)**2
            if (yss .le. ymax .and. yss .ge. ymin)then
                ey_xampu = a0*tamp !*dexp(-r2_1s/w2_1)
            else
                ey_xampu=0.d0
            endif
            phiG = (t-x0)-phi0
            bzda(j) = ey_xampu*dsin(phiG)                       % bz(i+1/2,j+1/2)
        end do
    endif
    return
end