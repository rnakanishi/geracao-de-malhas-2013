function M = similarMatrix(c1, c2)
    
    M = zeros(length(c1),length(c2));
    for ii=1:length(c1)
        for jj=1:length(c2)
            M(ii,jj) = norm( c1(:,ii) - c2(:,jj) );
        end
    end

end