
pulse_period = 5.3;
pulse_rise   = 3.2;
pulse_fall   = 2.4;
pulse_width  = 2.7;
pulse_delay  = 1.5;
fwhm         = 2.1;
tis = linspace(0, 20, 200);
laser_mode = 1
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,1);
plot(tis,tamp);

laser_mode = 2
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,2);
plot(tis,tamp);

laser_mode = 3
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,3);
plot(tis,tamp);

laser_mode = 4
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,4);
plot(tis,tamp);

laser_mode = 5
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,5);
plot(tis,tamp);

laser_mode = 6
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,6);
plot(tis,tamp);

laser_mode = 7
for i = 1:200
   tamp(i) =  pulse_shape(laser_mode, pulse_period, pulse_rise, pulse_fall, pulse_width, tis(i), pulse_delay, fwhm);
end
subplot(4,2,7);
plot(tis,tamp);