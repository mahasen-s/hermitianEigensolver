function test_lanczos2

% Fix random seed stream
%randStr = RandStream('dsfmt19937','Seed',0);
%rng(10)


H       = [1+1i,1-1i,-1+1i;
            -1-1i,2+7*1i,3+8*1i;
            4+9*1i,5+10*1i,6+11*1i];
        
eig(H)        




end

function H = random_hermitian(func,N)
H       = func(N,N) + 1i*func(N,N);
H       = H + transpose(conj(H));
end