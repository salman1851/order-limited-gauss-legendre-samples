clc; clear all;
L=32;   % band-limit   
M=5;    % order-limit
flm = rand(1,(L)^2 - (L-M)^2 + (L-M)) + 1i*rand(1,(L)^2 - (L-M)^2 + (L-M)); % original random signal
fo = glsht_inverse_order_limited(flm,L,M);      % original random signal in spatial domain
flmr = glsht_forward_order_limited(fo,L,M);     % recovered signal
error = mean(abs(flm-flmr));    % error between original and reconstructed signal