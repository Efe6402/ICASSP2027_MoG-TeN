function Xhat = reconstruct_tensor(C,alpha)

    [N,~,R] = size(C);
    T = size(alpha,1);

    Xhat = zeros(N,N,T);

    for r = 1:R
        for t = 1:T

            Xhat(:,:,t) = ...
                Xhat(:,:,t) ...
                +alpha(t,r)*C(:,:,r);
        end
    end
end
