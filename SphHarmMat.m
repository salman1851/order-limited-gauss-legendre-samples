function [Ylm_mat]=SphHarmMat(l,thetas)
 % Generate matrix of spherical harmonics for different values of theta and
 % phi=0 for given degree ll and orders 0,1..ll
 % Along rows: orders 0,1..ll
 % Along columns: thetas
%%
m_vec = 0:l;
 Ql=(-1).^m_vec*sqrt((2*l+1)/(8*pi)); % (-1)^m Q_l in note.tex/pdf
 Sl=legendre(l,cos(thetas),'sch');
 Ylm_mat=diag(Ql)*Sl; % outer product; pv is a row vector
 
 Ylm_mat(1,:) = sqrt(2)*Ylm_mat(1,:); % m=0 adjustment, see note.tex/pdf
 
end