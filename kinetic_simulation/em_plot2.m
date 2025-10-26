function [ex,ey,bz] = em_plot(prm,grid,cyc)
    

    file = sprintf('%s/emfield_n%8.8d',prm.fpath,cyc);
    data = importdata(file);
   
    x  = data(:,1);
    ex = data(:,2);
    ey = data(:,3);
    bz = data(:,4);
   
    
    figure; 
    subplot(2,3,1)
    plot(x,ex);
    xlabel('x');
    ylabel('Ex');
    
    subplot(2,3,2)
    plot(x,ey+bz);
    xlabel('x');
    ylabel('pump');
    
    subplot(2,3,3)
    plot(x,ey-bz);
    xlabel('x');
    ylabel('scatter');
end