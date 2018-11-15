function test_lanczos

% Fix random seed stream
%randStr = RandStream('dsfmt19937','Seed',0);
%rng(10)

N       = 100;
step    = 10;
n_pts   = 10;

% Generate Hermitian matrix
func    = @(M,N) sprand(M,N,0.5);
H       = random_hermitian(func,N);

if ishermitian(H)==true && numel(H)==N^2
    fprintf('Generated %dx%d Hermitian array\n',N,N)
else
    error('H is not Hermitian')
end


evals   = cell(n_pts,1);
for i = 1:n_pts
    L    = 1 + (i-1)*step;
    [V,T] = lanczos(H,L);
    E            = eig(T);
    evals{i}     = cat(2,L*ones(size(E)),E);
end

evals_actual    = eig(full(H));

figure(1)
clf

subplot(1,2,1)
hold on

%cellfun(@(x) plot(x,'ok'),evals)
colList         = jet(n_pts);
for i=1:n_pts
    plot(evals{i}(:,1),evals{i}(:,2),'o','color',colList(i,:),'markerfacecolor',colList(i,:))
end


for i=1:numel(evals_actual)
    h =refline(gca,[0,evals_actual(i)]);
    h.Color = 'w';
end
set(gca,'Color','k')
xlabel('Lanczos iteration')
ylabel('Energy')

box on
grid on

subplot(1,2,2)
min_eval    = cellfun(@(x) min(real(x(:,2))),evals);
eval_dist   = abs(min_eval-min(evals_actual));
semilogy(eval_dist,'-ok')
grid on
box on
xlabel('Lanczos iteration')
ylabel('Min Energy')



end

function H = random_hermitian(func,N)
H       = func(N,N) + 1i*func(N,N);
H       = H + transpose(conj(H));
end