function distf = distf_plot(prm,grid,cyc)

    for n=1:prm.nsp
        distf{n} = zeros(grid.nx,grid.np(n));
    end

    for n =1:prm.nsp
        h5file='';
         h5file=sprintf('%s/distf_#%d_n%8.8d.h5',prm.fpath,n,cyc);
%          h5file=sprintf('%s/distf_#%d_n%5.5d.h5',prm.fpath,n,cyc);
        domain=hdf5read(h5file,'/domain');
        distf{n}= hdf5read(h5file,'/f');
    end
    subplot(2,1,1)
    imagesc(grid.coord_x(2:grid.nx-1),grid.coord_p1,distf{1}(2:grid.nx-1,2:grid.np(1)-1)');
    colorbar;
    subplot(2,1,2)
    imagesc(grid.coord_x(2:grid.nx-1),grid.coord_p2,distf{2}(2:grid.nx-1,2:grid.np(2)-1)');
    colorbar;
end