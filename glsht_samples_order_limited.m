function [THETA, FI, qwts] = glsht_samples_order_limited(L,M)

[zpts qwts]=lgwt(L,-1,1);
THETA = transpose(acos(zpts));
FI = pi*(2*(0:1:2*(M)))/(2*M+1);
qwts = qwts';

end