%% function to be used in Perfusio Tool Main for PL model calculation
function [ Ct ] = PLMmain( par_PL,time )
for t = 1: length(time)
    
    if (time(t)<=par_PL(1))
        
        Ct(t) = par_PL(3);
        
    elseif (time(t)<=par_PL(2))
        
        Ct(t) = par_PL(3) + par_PL(4)*(time(t)-par_PL(1));
        
    else
        
        Ct(t) = par_PL(3) + par_PL(4)*(par_PL(2)-par_PL(1))+par_PL(5)*(time(t)-par_PL(2));
        
    end
end
end

