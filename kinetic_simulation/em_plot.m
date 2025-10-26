function [ex,ey,bz] = em_plot(prm,grid,cyc)
    

%      file = sprintf('%s/emfield_n%5.5d',prm.fpath,cyc);
     file = sprintf('%s/emfield_n%8.8d',prm.fpath,cyc);
    data = importdata(file);
   
    ex = data(:,2);
    ey = data(:,3);
    bz = data(:,4);
   
    
    figure; 
    subplot(2,3,1)
    plot(data(:,1),data(:,2))
    xlabel('x');
    ylabel('Ex');
    
    subplot(2,3,2)
    plot(data(:,1),data(:,3))
    xlabel('x');
    ylabel('Ey');
    ey = data(:,3);
    subplot(2,3,3)
    plot(data(:,1),data(:,4))
    xlabel('x');
    ylabel('Bz');
end