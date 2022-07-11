function output = test4
%
%   Parent function - trying to find how to use fzero (with parameters)
%
    A = 6;
    B = 4;
    C = 25;
    D = 3;
    target = -10;
%
    myfun = @(A,B,x,C,D) (test3(A,B,x,C,D) - target);  % parameterized function                
    fun = @(x) myfun(A,B,x,C,D);    % function of x alone
    [output,value] = fzero(fun,-10);

return

end

