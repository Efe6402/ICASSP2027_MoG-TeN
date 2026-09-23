function [L,V,ell] = temporal_laplacian(T)
%TEMPORAL_LAPLACIAN First-difference Gram matrix and sorted eigendecomposition.

D = zeros(T-1,T);
for t = 1:T-1
    D(t,t) = -1;
    D(t,t+1) = 1;
end
L = D'*D;
L = (L+L')/2;

[V,Dlam] = eig(L);
ell = real(diag(Dlam));
[ell,idx] = sort(ell,'ascend');
V = real(V(:,idx));
ell(abs(ell)<1e-14) = 0;
end
