function val = logfacbyfac(m,n)
% Returns log(m!/n!), m can be an array

k = length(m);
val = zeros(1,k);

for l=1:k
    if m(l) == n
        val(l) = 0;
    else if m(l) > n
            for i=n+1:m(l)
                val(l) = val(l)+log(i);
            end
        else
            error('Parameters mismatched.');
        end 
    end
end