function [V,T] = lanczos(H,L)
% Accepts Hermitian matrix, produces tridiagonal matrix
if ishermitian(H)==false
    error('H is not Hermitian')
end

N       = size(H,1);

%% LANCZOS
a       = zeros(L,1);
b       = zeros(L,1);
V       = zeros(N,L);

init    = randn(N,1);

v       = init;
b(1)    = norm(v);
v       = v/b(1);

V(:,1)  = v;

w       = H*v;
a(1)    = dot(w,v);
disp(a(1))
w       = w - a(1)*v;

for i=2:L
    b(i) = norm(w);
    if b(i)>eps
        v   = w/b(i);
    else
        % Use Gram-Schmidt to construct new vector
        fprintf('b(%d) is zero, constructing v_%d\n',i,i)
        is_orthogonal = false;
        while is_orthogonal == false
            v_temp = rand(N,1);
            v_temp = v_temp/norm(v_temp);
            v_temp= v_temp + sum(v_temp.*V(:,1:i-1),2);
            v_temp2_norm = norm(v_temp);
            if v_temp2_norm > eps
                v_temp = v_temp/v_temp2_norm;
                is_orthogonal = true;
            end
        end
        v   = v_temp;
    end
    V(:,i)  = v;
    
    w       = H*v;
    a(i)    = dot(w,v);
    w       = w - a(i)*v - b(i)*V(:,i-1);
end

% We are guaranteed that T is real
T       = diag(real(a)) + diag(real(b(2:end)),1) + diag(real(b(2:end)),-1);
%T       = diag(a) + diag(b(2:end),1) + diag(b(2:end),-1);




end