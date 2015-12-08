function [gg] = ism_mask(gg,dd,oo)
% Assign default parameters and options
% Inputs 
%   gg struct of grid
%   dd struct of topography
%   oo struct of options [optional]
% Outputs
%   gg struct of grid

% 16 October 2015 


A = (dd.H > 0);
B = (dd.mask == 1);
[C1, C2] = gradient(single(A)); C = C1 | C2;
[D1, D2] = gradient(single(B)); D = D1 | D2;

%Assign nodes as ice margin, ice interior, interior boundary, or exterior to the problem. 
%Nodes belong to only one category.
gg.nin = A & B & ~C & ~D;               
gg.nmgn = A & B & C;                    
gg.nbnd = A & B & D & ~C;
gg.next = ~(gg.nin | gg.nmgn | gg.nbnd);

%Sampling Operators
S_h = diag(dd.mask(:)); S_h = S_h(any(S_h,2),:); %H-grid
A = conv2(dd.mask, [1 1]/2) > 0; S_u = diag(A(:) > 0); S_u = S_u(any(S_u,2),:); % U-grid; zero padding is fine
A = conv2(dd.mask, [1 1]'/2) > 0; S_v = diag(A(:) > 0) ; S_v = S_v(any(S_v,2),:); %V-grid; zero padding is fine
A = padarray(dd.mask,[1,1], 'symmetric'); B = conv2(A,[1 1; 1 1]/4, 'valid'); %C-grid; need symmetric padding
S_c = diag(B(:) == 1); S_c = S_c(any(S_c,2),:); 

gg.S_h = S_h;
gg.S_u = S_u;
gg.S_v = S_v;
gg.S_c = S_c;
end
