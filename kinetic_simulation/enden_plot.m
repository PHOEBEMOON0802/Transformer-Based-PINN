function density=enden_plot(prm,grid,cyc)
% ---------------------------
% plot density and energy
% ---------------------------
    density = zeros(grid.nx-2, prm.nsp);
    energy  = zeros(grid.nx-2, prm.nsp);
    entropy = zeros(grid.nx-2, prm.nsp);
    jxi     = zeros(grid.nx-2, prm.nsp);
    jyi     = zeros(grid.nx-2, prm.nsp);
    jx      = zeros(grid.nx-2,1);
    jy      = zeros(grid.nx-2,1);
    for n = 1:prm.nsp
         file = sprintf('%s/enden_#%d_n%8.8d',prm.fpath, n,cyc);
%           file = sprintf('%s/enden_#%d_n%5.5d',prm.fpath, n,cyc);
            data = importdata(file);
            x  = data(:,1);
            for i = 1:length(x)
                density(i,n) = density(i,n) + data(i,2);
                energy(i,n)  = energy(i,n)  + data(i,3);
                entropy(i,n) = entropy(i,n) + data(i,4);
                jxi(i,n) = jxi(i,n) + data(i,5);
                jx(i)    = jx(i) + data(i,5)*grid.dp(n);
                jyi(i,n) = jyi(i,n) + data(i,6);
                jy(i)    = jy(i) + data(i,6)*grid.dp(n);
            end
    end
    for n=1:prm.nsp
        figure;
        subplot(2,3,1);
        plot(x,density(:,n)*grid.dp(n));
        xlabel('x');
        ylabel(['density_#',int2str(n)]);
    
        subplot(2,3,2);
        plot(x,energy(:,n)*grid.dp(n));
        xlabel('x');
        ylabel(['energy_#',int2str(n)]);
    
        subplot(2,3,3);
        plot(x,entropy(:,n)*grid.dp(n));
        xlabel('x');
        ylabel(['entropy_#',int2str(n)]);
   
        subplot(2,3,4);
        plot(x,jxi(:,n)*grid.dp(n));
        xlabel('x');
        ylabel(['jxi_#',int2str(n)]);
    
        subplot(2,3,5);
        plot(x,jyi(:,n)*grid.dp(n));
        xlabel('x');
        ylabel(['jyi_#',int2str(n)]);
    end
    figure;
    subplot(2,1,1);
    plot(x,jx);
    xlabel('x');
    ylabel(['jx']);
    subplot(2,1,2);
    plot(x,jy);
    xlabel('x');
    ylabel(['jy']);

    
    for n=1:prm.nsp
        totN = sum(density(:,n))*grid.dp(n)*grid.dx;
        disp(['tot N#',int2str(n),' = ',num2str(totN)]);
    end
    %totN;
    
end
