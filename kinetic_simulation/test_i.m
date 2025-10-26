clear global
    warning off
    global prm f g t xs time
    

    global rho u p q ef
    
    prm = initial();
    fid_t = fopen('t.txt', 'wt');

 
    
    
    time = 0.0;
    
    
    
    for it = 1:prm.nt
       %time = time + prm.dt;
       t = t + prm.dt;
       fprintf(fid_t, '%f\n', t);

       
    end

    
    fclose(fid_t);
    