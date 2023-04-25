function flm = glsht_forward_order_limited(f,L,M)
    if ~(M<L)
          error('M must be less than L');
    end    
    if ~isreal(L)
          error('Harmonic band-limit must be real');
    end
    
    flm = zeros(1, L^2 - (L-M)^2 + L-M);
    [thetas, FI,qwts] = glsht_samples_order_limited(L,M); 
    gm = fftshift(fft(f, [],2),2)/(2*M+1);

    for ll=0:1:M
        Ylm_mat=transpose(SphHarmMat(ll,thetas));       
        for m=0:1:ll            
            Ylm  = Ylm_mat(:,m+1);
            Ylm_neg = (-1)^m*Ylm;
            flm(ll^2+ll+m+1) = 2*pi*qwts*(gm(:,M+1+m).*Ylm);     
            flm(ll^2+ll-m+1) = 2*pi*qwts*(gm(:,M+1-m).*Ylm_neg);                       
        end
    end
    
    for ll=M+1:1:L-1
        Ylm_mat=transpose(SphHarmMat(ll,thetas));
        for m=0:1:M            
            Ylm  = Ylm_mat(:,m+1);
            Ylm_neg = (-1)^m*Ylm;
            flm(ll^2+ll - ((ll-M)^2 -1) +m) = 2*pi*qwts*(gm(:,M+1+m).*Ylm);     
            flm(ll^2+ll - ((ll-M)^2 -1) -m) = 2*pi*qwts*(gm(:,M+1-m).*Ylm_neg);            
        end
    end    
    
end