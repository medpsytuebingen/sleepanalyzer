%************************************************************************** 
% function 'interleave' allows interleaving of two equally sized arrays
% works at least with MATLAB 7.7
% by Til Ole Bergmann, Kiel, Germany
% t.bergmann@neurologie.uni.kiel.de
% last modified 2009-07-29 by TOB
%**************************************************************************
function Y = interleave(X1,X2)
    interleavematrix = [X1;X2];
    Y = reshape(interleavematrix,1,length(interleavematrix)*2);
end