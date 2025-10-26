
%                                                                                      !
%      1D1V 非相对论Vlasov-Possion方程代码                                                !
%                                                                                       !
% 
%       计算1D1V-Vlasov-Possion方程
%           df/dt + v df/dx + qE/m df/dv = 0
%          dE/dx = q rho - q<densm>
% 单位:
%       时间： 1/wp   空间： lambda(德拜长度)   速度: 电子热速度
%
% *计算结果：
%     1. 全局诊断信息存储在: xxxx.diag文件中，通过修改diagnostic_header, diagnostic函数可增加输出
%     2. 分布函数信息存储在: xxxx.distf文件中，每次输出 4 + 8+(N+2)*(2*M+1)*8+4
%                                         其中第一个和最后一个为半个字节
%                                         共输出 nt/nplot + 2 次
%
%     3. 电场信息        : xxxx.field文件中，每次输出 4 + 8 + (N+2)*8 + 4个字节
%                                         共输出 nt/nplot + 2 次
%   ---------------------------------------------------------------------------------------------
%   计算结果：
%   ---------------------------------------------------------------------------------------------
%  time density momentum energy_f energy_em energy L2 entropy log|E| EF(MID-1) EF(MID) EF(MID+1)
%       0.000    0.1000000000E+01    0.5476612433E-17    0.5000000000E+00    0.9999999839E-06    0.5000010000E+00   -0.3544909474E+01    0.1783090435E+02   -0.6561181697E+01   -0.9813534786E-04   -0.1853943363E-18    0.9813534786E-04
%       0.200    0.1000000000E+01    0.1052587754E-16    0.5000000125E+00    0.9745777251E-06    0.5000009871E+00   -0.3544909474E+01    0.1783090435E+02   -0.6574057191E+01   -0.9688086883E-04   -0.1795589146E-15    0.9688086884E-04
%       0.400    0.1000000000E+01   -0.1012526272E-16    0.5000001084E+00    0.8572735390E-06    0.5000009656E+00   -0.3544909474E+01    0.1783090435E+02   -0.6638180803E+01   -0.9086784033E-04   -0.4132742571E-15    0.9086784033E-04
%       0.600    0.1000000000E+01   -0.2400605573E-16    0.5000002797E+00    0.6723899417E-06    0.5000009521E+00   -0.3544909474E+01    0.1783090435E+02   -0.6759640107E+01   -0.8048127070E-04   -0.1557248898E-15    0.8048127070E-04
%       0.800    0.1000000000E+01    0.4886922369E-16    0.5000004907E+00    0.4580962905E-06    0.5000009488E+00   -0.3544909474E+01    0.1783090435E+02   -0.6951519627E+01   -0.6643598714E-04   -0.6895176416E-16    0.6643598714E-04
%       1.000    0.1000000000E+01    0.6711257233E-16    0.5000006994E+00    0.2562324644E-06    0.5000009556E+00   -0.3544909474E+01    0.1783090435E+02   -0.7242016780E+01   -0.4969138592E-04    0.1633943178E-14    0.4969138592E-04
%       1.200    0.1000000000E+01    0.3643682895E-16    0.5000008680E+00    0.1021613091E-06    0.5000009702E+00   -0.3544909474E+01    0.1783090435E+02   -0.7701782815E+01   -0.3137828874E-04    0.6093349617E-15    0.3137828874E-04
%       1.400    0.1000000000E+01    0.5120314794E-16    0.5000009715E+00    0.1677432495E-07    0.5000009883E+00   -0.3544909474E+01    0.1783090435E+02   -0.8605134608E+01   -0.1271307552E-04    0.1494707227E-14    0.1271307553E-04
%       1.600    0.1000000000E+01    0.8341047759E-16    0.5000010025E+00    0.2686085288E-08    0.5000010052E+00   -0.3544909474E+01    0.1783090435E+02   -0.9521016903E+01    0.5093938252E-05    0.7675554223E-15   -0.5093938250E-05
%       1.800    0.1000000000E+01    0.6516024541E-16    0.5000009716E+00    0.4540183200E-07    0.5000010170E+00   -0.3544909474E+01    0.1783090435E+02   -0.8107283100E+01    0.2092835246E-04    0.6996213682E-15   -0.2092835246E-04
%   ---------------------------------------------------------------------------------------------
    clear global
    warning off
    global prm time
    global v x t
    global f f2 g
    global ef rho v u p q rhoe
    global fv

    fid_x = fopen('x.txt', 'wt');
    fid_t = fopen('t.txt', 'wt');
    fid_v = fopen('v.txt', 'wt');
    fid_f = fopen('f.txt', 'wt');
    fid_n = fopen('n.txt', 'wt');
    fid_u = fopen('u.txt', 'wt');
    fid_p = fopen('p.txt', 'wt');
    fid_q = fopen('q.txt', 'wt');
    fid_E = fopen('E.txt', 'wt');
    
    %fid_gr = fopen('gr.txt', 'wt');
    
    prm = initial();
    efield3;
    time = 0.0;
    fprintf(fid_x, '%f\n', x);
    fprintf(fid_t, '%f\n', t);
    fprintf(fid_v, '%f\n', v);
    fprintf(fid_f, '%f\n', f);
    fprintf(fid_n, '%f\n', rhoe);
    fprintf(fid_u, '%f\n', u);
    fprintf(fid_p, '%f\n', p);
    fprintf(fid_q, '%f\n', q);
    fprintf(fid_E, '%f\n', ef);
    
    
    diag_f_id = fopen('diag_f_out','w');

    
    %ndiag = 10;
    for it = 1:prm.nt
       time = time + prm.dt;     
       spline_x;
       advection_x_semi;
       f = g;
       efield3;
       t = t + prm.dt;
       fprintf(fid_x, '%f\n', x);
       fprintf(fid_t, '%f\n', t);
       fprintf(fid_v, '%f\n', v);
       fprintf(fid_f, '%f\n', f);
       fprintf(fid_n, '%f\n', rhoe);
       fprintf(fid_u, '%f\n', u);
       fprintf(fid_p, '%f\n', p);
       fprintf(fid_q, '%f\n', q);
       fprintf(fid_E, '%f\n', ef);
       
       %diagnostic(diag_f_id,time);
       if mod(it, 10) == 0
           figure(1);
           set(figure(1), 'Color', 'white');
           phasefig;
           
           
           grid on;
           set(gca, 'FontSize', 12, 'FontName', 'BookmanOld Style');
           xlabel('{\itx}', 'FontSize', 14);
           xlim([min(x), max(x)]);
           ylabel('{\itv}', 'FontSize', 14);
           ylim([min(v), max(v)]);
           legend('FontSize', 14, 'Box','off');
           a = {'Time = '};
           b = num2str(time);
           str = strcat(a, b);
           title(str, 'FontSize', 14);
           drawnow;
           gifgenerator;
           

           %int_e_2(it/100) = sum(e_2);
           %t(it/100) = it/100;
           %fprintf(fid_gr, '%f\n', int_e_2);
       end

       spline_v
       advection_v
       f = g;

       
       spline_x
       advection_x_semi
       f = g;
       
       
       
    end
    %figure(2);
    %plot(t(1:300), int_e_2(1:300));
    %saveas(2, 'gr.jpg');
    fclose(fid_x);
    fclose(fid_t);
    fclose(fid_v);
    fclose(fid_f);
    fclose(fid_n);
    fclose(fid_u);
    fclose(fid_p);
    fclose(fid_q);
    fclose(fid_E);
    
    %fclose(fid_gr);
    %fclose(diag_f_id);
