function [f_spatial,thetas,phis] = glsht_inverse_order_limited(flm,L,M)
    if ~(M<L)
          error('M must be less than L');
    end    
    if ~isreal(L)
          error('Harmonic band-limit must be real');
    end
    if ~(sum(size(flm))== length(flm)+1) 
          error('flm must be a vector and not a matrix');
    end
    if ~(length(flm)==L^2 - (L-M)^2 + (L-M))
          error('flm must be a vector of size L^2 - (L-M)^2 + (L-M)');
    end
    
    [thetas, phis, qwts] = glsht_samples_order_limited(L,M);
    T=L;
    F=2*M+1;
    f=zeros(T,F);
    for ll=0:1:M
        Ylm_mat=transpose(SphHarmMat(ll,thetas));
        fm = zeros(ll+1,1);
        fm_neg = zeros(ll+1,1);
        for m=0:1:ll
            fm(m+1) = flm(ll^2+ll+m+1);
            fm_neg(m+1) = flm(ll^2+ll-m+1);    
        end
        f(:,M+1:M+1+ll) = f(:,M+1:M+1+ll)+Ylm_mat*diag(fm);
        f(:,M+1-ll:M+1-1) = f(:,M+1-ll:M+1-1)+fliplr(Ylm_mat(:,2:ll+1)*diag(fm_neg(2:ll+1))*diag((-1).^(1:1:ll)));
    end

    for ll=M+1:1:L-1
        Ylm_mat=transpose(SphHarmMat(ll,thetas));
        fm = zeros(M+1,1);
        fm_neg = zeros(M+1,1);
        for m=0:1:M            
            fm(m+1) = flm(ll^2+ll - ((ll-M)^2 -1) +m);
            fm_neg(m+1) = flm(ll^2+ll  - ((ll-M)^2 -1) -m);                            
        end       
        f(:,M+1:M+1+M) = f(:,M+1:M+1+M)+Ylm_mat(:,1:M+1)*diag(fm);
        f(:,M+1-M:M+1-1) = f(:,M+1-M:M+1-1)+fliplr(Ylm_mat(:,2:M+1)*diag(fm_neg(2:M+1))*diag((-1).^(1:1:M)));
    end
    f_spatial = F*(ifft(ifftshift(f,2),F,2));
end
