function test_eigs


dense_mat_0 = [ 1+1i,1-1i,-1+1i;
                -1-1i,2.0+7*1i, 3.0+8*1i;
                4.0+9*1i, 5.0+10*1i, 6.0 + 11*1i];
            
disp(dense_mat_0)
eig(dense_mat_0)

end