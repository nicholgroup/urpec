function [polysOut] = fixPolys2(polysIn)
%fixPolys2 calls fixPoly2 for a cell array of polygons.

counter=1;
polysOut={};
for ip=length(polysIn):-1:1
    polysTmp=fixPoly2(polysIn{ip});
    
    pt=num2cell(polysTmp);    
    
    if counter==1
        polysOut=pt;
    else
        polysOut=[polysOut pt];
    end
    counter=counter+1;
    
end

end

