function R_2 = rsqrd(Y,Y_est)

[n_out,n_samples] = size(Y);
R_2 = zeros(n_out,1);


y_ave = mean(Y,2);
for i = 1:n_out
    ss_tot = (Y(i,:)-y_ave(i)).^2;
    
    if ss_tot < 1e-10
        ss_tot = 1e-10;
        
    end
    ss_tot = sum(ss_tot);
    
    ss_res = (Y(i,:)-Y_est(i,:)).^2;
    ss_res = sum(ss_res);
    
    R_2(i) = 1-ss_res/ss_tot;
end


end