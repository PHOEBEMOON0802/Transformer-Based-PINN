function amp = pulse_shape(ilas, tpulse, trise, tfall, tp, t, tdelay, fwhm)
% ================================
% laser pulse shape
% -----------------
% Parameters:
%   ilas    laser model 
%   tpulse  pulse duration
%   trise   pulse rise time
%   tfall   pulse fall time
%   tp      pusle width
%   t       current time
%   tdelay  delay time
%   fwhm    full width half midwidth
%  
% ================================

if(t < 0)
    amp= 0.0;
    return;
end
tt=t-tdelay;
switch ilas
    case {1,11}
        if tt <= tpulse && tt > 0 
            amp = sin(pi*tt/tpulse);
        else
            amp = 0.0;
        end 
    case {2,12}
        if tt <= trise+max(trise,tp)
            amp = exp(-(tt-trise)^2/tp^2);
        else
            amp = 0.0;
        end
    case {3,13}
        if(tt <= tpulse && tt > 0)
            amp = sin(pi*tt/tpulse)^2;
        else
            amp = 0.0;
        end
    case {4,14}
        tpend = trise+tpulse
        if(tt <= 0) 
            amp = 0.0;
        elseif tt<=trise
            amp = sqrt(tt/trise);
        elseif tt<=tpend
            amp = 1.0;
        elseif tt<= tpend+tfall
            amp = sqrt(1.0-(tt-tpend)/tfall)
        else
            amp = 0.0;
        end
    case {5,15}
        if tt <= 0
            amp = 0.0;
        elseif tt<=trise
            amp = sin(pi*tt/(2*trise));
        else
            amp = 1.0;
        end
    case {6,16}
       tpend=trise+tpulse
       if tt<=0
            amp=0.0;
       elseif  tt <= trise
            amp= exp(-(tt-trise)^2/(trise/fwhm)^2);      ! FWHM
       elseif tt <= tpend
            amp=1.d0;
       elseif tt <= tpend+tfall
            amp=exp(-(tt-tpend)^2/(tfall/fwhm)^2);
       else
            amp=0.d0;
       end
    case {7,17}
        if tt <= trise+max(trise,tp)
            amp = exp(-log(2.d0)*(tt-trise)^4/(tpulse/2.)^4); ! FWHM
        else
            amp = 0.d0;
        end 
end
end 


    