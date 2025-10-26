function [prm,grid]=init_params

    fpath = '/home/hp/workspace/RLV1DA/data/BRS_S2';
    
    prm.fpath = fpath;
    nsp = 2;
    prm.nsp   =  nsp;
    
    grid.coord_x = importdata([fpath,'/grid_coord_x.dat']);
    grid.nx = length(grid.coord_x);
    grid.dx = grid.coord_x(3)-grid.coord_x(2);
    for n = 1:nsp
        file = sprintf('%s/grid_coord_p_#%d',fpath,n);
        eval_str = ['grid.coord_p',int2str(n),'=importdata(file);'];
        eval(eval_str);
        grid.np(n) = length(eval(['grid.coord_p',int2str(n)],';'));
        grid.dp(n) = eval(['grid.coord_p',int2str(n),'(3)-grid.coord_p',int2str(n),'(2);']);
    end 
    
end 