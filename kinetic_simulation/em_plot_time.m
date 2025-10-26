function [time,eng_pump,eng_scat] = em_plot_time(prm,grid,step,cyc)
    
    nt = length(0:step:cyc);
    time     = zeros(nt,1);
    eng_pump = zeros(nt,1);
    eng_scat = zeros(nt,1);
    k = 1;
    for n = 0:step:cyc
%           file = sprintf('%s/emfield_n%5.5d',prm.fpath,n);
          file = sprintf('%s/emfield_n%8.8d',prm.fpath,n);
        data = importdata(file);
        x  = data(:,1);
        ex = data(:,2);
        ey = data(:,3);
        bz = data(:,4);
        
        fp = ey+bz;
        fm = ey-bz;
     
        eng_pump(k)=sum(fp.^2)*grid.dx;
        eng_scat(k)=sum(fm.^2)*grid.dx;
        time(k)    =n*grid.dx;
        k = k+1;
    end
    subplot(2,1,1)
    plot(time,eng_pump,'r');
    subplot(2,1,2)
    plot(time,eng_scat,'g');