function diagnostic(diag_f_id,time)
    global prm
    global x v f ef
    den = 0;
    en_kin = 0;
    en_em  = 0;
    mom = 0;
    l2  = 0;
    entropy = 0;
    for i = 2:prm.n+1
        for j = 1:2*prm.m+1
            den     = den + f(i,j);
            en_kin  = en_kin + v(j) * v(j) * f(i,j);
            mom     = mom    + v(j) * f(i,j);
            ff      = abs(f(i,j)) + 1.0e-16;
            l2      = l2 + ff*ff;
            entropy = entropy + ff*log(ff);
        end 
        en_em    = en_em + ef(i)*ef(i);
    end
    den = den * prm.dv / (prm.n); 			
    en_kin = 0.5*en_kin*prm.dv/(prm.n);
    en_em  = 0.5*en_em/(prm.n);					
    mom = mom * prm.dv /(prm.n);
    l2 = -l2*prm.dx*prm.dv;
    entropy = -entropy*prm.dx*prm.dv;            			
    
    fprintf(diag_f_id,'%f %f %f %f %f %f %f %f\n',time,den, mom,en_kin,en_em,l2,entropy,ef(prm.n/2+2));
    
return;