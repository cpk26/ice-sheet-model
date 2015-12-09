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
S_h = spdiags([dd.mask(:) > 0],0,gg.nIJ,gg.nIJ); S_h = S_h(any(S_h,2),:); %H-grid

%Multiplication by first/last columns/rows of mask is to handle periodic BC in
%centering operator
A = conv2(dd.mask, [1 1]/2) > 0;                    % U-grid;
S_u = spdiags(A(:) > 0,0,(gg.nI+1)*(gg.nJ),(gg.nI+1)*(gg.nJ)); S_u = S_u(any(S_u,2),:); 

A = conv2(dd.mask, [1 1]'/2) > 0;                   % V-grid;
S_v = spdiags(A(:) > 0, 0,(gg.nI)*(gg.nJ+1),(gg.nI)*(gg.nJ+1)); S_v = S_v(any(S_v,2),:); 

A = padarray(dd.mask,[1,1], 'symmetric'); B = conv2(A,[1 1; 1 1]/4, 'valid'); %C-grid; need symmetric padding
S_c = spdiags(B(:) == 1,[0],(gg.nI+1)*(gg.nJ+1),(gg.nI+1)*(gg.nJ+1)); S_c = S_c(any(S_c,2),:); 

gg.S_h = S_h;
gg.S_u = S_u;
gg.S_v = S_v;
gg.S_c = S_c;
end

%A = reshape(c_hv*dd.mask(:),nJ+1,nI); A(1,:) = A(1,:).*mask(1,:); A(end,:) = A(end,:).*mask(end,:); 
%A = reshape(gg.c_hu*dd.mask(:),gg.nJ,gg.nI+1); A(:,1) = A(:,1).*dd.mask(:,1); A(:,end) = A(:,end).*dd.mask(:,end); 
% A = conv2(dd.mask, [1 1]/2) > 0; S_u = diag(A(:) > 0); S_u = S_u(any(S_u,2),:); % U-grid; zero padding is fine
% A = conv2(dd.mask, [1 1]'/2) > 0; S_v = diag(A(:) > 0) ; S_v = S_v(any(S_v,2),:); %V-grid; zero padding is fine
