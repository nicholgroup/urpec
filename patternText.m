function [BB] = patternText(str,targetLen,cntr)
%[BB] = patternText(str,targetLen)
% Outputs a cell array of polygons corresponding to the text. targetLen is
% the target text legnth. This is useful for generating text patters for
% lithography pattern files. The basic idea is to make a graph with text,
% convert it to an image, and then extract the text from there.
%
% str is the text you want to write.
%
% targetLen is the target length of the text
%
% cntr is the center of the text.
%
% BB is the cell array of polygons. Note that the first column of each
% element is the y value, and the second element is the X value.

textLen=30/9*length(str); %empirical scaling for this font size
xlen=2*textLen+40;
ylen=100;
figure(1111); clf; colormap gray;
set(gcf,'color','white','units','normalized','position',[.1 .1 .8 .8])
imagesc((1:xlen),(1:ylen),ones(100,100));
caxis([0 1]);
set(gca,'visible','off')

text(xlen/2-textLen/2,ylen/2,str,'FontSize',50,'Interpreter','none');
%caxis([0 1]); 
%set(gcf,'color','white','units','normalized','position',[.1 .1 .8 .8])

tim = getframe(1111);
tim2 = tim.cdata;
tmask = tim2<255;

tmask=tmask(:,:,1);

%The code below is borrowed from urpec. To make this easy, lets just rename
%some variables to make the code compatible.
doseNew=tmask;

close(1111);

mp=size(doseNew,1);
np=size(doseNew,2);

debug=0;
dvals=1;
layer=struct();
config=struct();
config.subfieldSize=50;

for i=1:length(dvals)
    fprintf('layer %d...\n',i);
    
    %Compute shot map for each layer
    if i==1
        layer(i).shotMap=doseNew;
        layer(i).dose=dvals(i);
        
    elseif i==length(dvals)
        layer(i).shotMap=doseNew>dvalsl(i).*shape;
        layer(i).dose=dvals(i);
        
    else
        layer(i).shotMap=(doseNew>dvalsl(i)).*(doseNew<dvalsr(i));
        layer(i).dose=(dvalsl(i)+dvalsr(i))/2;
    end
    
    %Compute the actual mean dose in the layer.
    dd=doseNew;
    dd(isnan(dd))=0;
    doseSum= sum(sum(layer(i).shotMap));
    if doseSum>0
        layer(i).meanDose=sum(sum(layer(i).shotMap.*dd))/sum(sum(layer(i).shotMap));
        dvalsAct(i)=layer(i).meanDose;
    else
        layer(i).meanDose=dvals(i);
        dvalsAct(i)=dvals(i);
    end
    
    %break into subfields and find boundaries. This is necessary because
    %designCAD can only handle polygons with less than ~200 points.
    
    subfieldSize = config.subfieldSize;
    
    layer(i).boundaries={};
    count=1;
    
    decrease = 0;
    disp = 0;
    m=1;
    n=1;
    xsubfields=ceil(mp/subfieldSize);
    ysubfields=ceil(np/subfieldSize);
    decrease_sub = 1;
    while (n <= ysubfields)    %change to xsubfields for horizontal writing
        %decrease_sub=1;
        m = 1;
        %decrease_sub = 1;
        while (m <= xsubfields) %change to ysubfields for horizontal writing
            %EJC: add decrease_sub factor for halving sub field size when polygons are large
            if decrease
                subfieldSize = round(subfieldSize/decrease_sub);
                decrease = 0;
                disp = 0;
            end
            xsubfields=ceil(mp/subfieldSize);
            ysubfields=ceil(np/subfieldSize);
            if ~disp
                display(['trying subfield size of ' num2str(subfieldSize) '. There are a total of ' num2str(xsubfields*ysubfields) ' subfields...']);
                disp = 1;
            end
            proceed = 0;
            
            
            subfield=zeros(mp,np);
            %display(['Trying subfield size: ' num2str(mp/decrease_sub) 'x' num2str(np/decrease_sub)]);
            if (m-1)*subfieldSize+1<min(m*subfieldSize,mp)
                xinds=(m-1)*subfieldSize+1:1:min(m*subfieldSize,mp);
            else
                xinds=(m-1)*subfieldSize+1:1:min(m*subfieldSize,mp);
            end
            if (n-1)*subfieldSize+1 < min(n*subfieldSize,np)
                yinds=(n-1)*subfieldSize+1:1:min(n*subfieldSize,np);
            else
                yinds=(n-1)*subfieldSize+1:1:min(n*subfieldSize,np);
            end
            
            xstart=xinds(1);
            ystart=yinds(1);
            
            %double the size of each shot map to avoid single pixel
            %features. No longer used. It was an attemp to avoid
            %single-pixel features.
            %             xinds=reshape([xinds;xinds],[1 2*length(xinds)]);
            %             yinds=reshape([yinds;yinds],[1 2*length(yinds)]);
            
            subdata=layer(i).shotMap(xinds,yinds);
            
            %Now "smear" out the shot map by one pixel in each direction. This makes sure that
            %subfield boundaries touch each other.
            sdll=padarray(subdata,[1,1],0,'pre');
            sdur=padarray(subdata,[1,1],0,'post');
            sdul=padarray(padarray(subdata,[1,0],'pre'),[0,1],'post');
            sdlr=padarray(padarray(subdata,[1,0],'post'),[0,1],'pre');
            
            sd=sdll+sdur+sdul+sdlr;
            sd(sd>0)=1;
            subdata=sd;
            
            if (sum(subdata(:)))>0
                
                [B,L,nn,A]=bwboundaries(subdata);
                
                if debug
                    figure(777); clf; hold on
                    fprintf('Layer %d subfield (%d,%d) \n',i,m,n);
                    for j=1:length(B)
                        try
                            bb=B{j};
                            subplot(1,2,1); hold on;
                            plot(bb(:,2),bb(:,1))
                            axis([0 subfieldSize 0 subfieldSize]);
                            subplot(1,2,2)
                            imagesc(subdata); set(gca,'YDir','norm')
                            axis([0 subfieldSize 0 subfieldSize]);
                        end
                        pause(1)
                        
                    end
                    drawnow;
                    pause(1);
                end
                
                if ~isempty(B)
                    for b=1:length(B)
                        
                        %Find any holes and fix them by adding them
                        %appropriately to enclosing boundaries
                        enclosing_boundaries=find(A(b,:));
                        for k=1:length(enclosing_boundaries)
                            b1=B{enclosing_boundaries(k)};% the enclosing boundary
                            b2=B{b}; %the hole
                            
                            %Only keep the hole if it has >0 area.
                            %Also, only keep the hole if there are no
                            %overlapping lines in the hole. There should be only one
                            %pair of matching vertices per polygon.
                            if polyarea(b2(:,1),b2(:,2))>0 && (size(b2,1)-size(unique(b2,'rows'),1)==1)
                                b1=[b1; b2; b1(end,:)]; %go from the enclosing boundary to the hole and back to the enclosing boundary
                            end
                            B{b}=[]; %get rid of the hole, since it is now part of the enclosing boundary
                            B{enclosing_boundaries(k)}=b1; %add the boundary back to B.
                        end
                    end
                    
                    %Add boundaries to layer
                    for b=1:length(B)
                        if ~isempty(B{b}) && polyarea(B{b}(:,1),B{b}(:,2))>0
                            
                            %remove unnecessary vertices
                            B{b}=simplify_polygon(B{b});
                            
                            %check for large polygons
                            if size(B{b},1)>200%200
                                fprintf('Large boundaries in layer %d. Halving subfield size and retrying... \n',i);
                                
                                %If large boundaries, make the subfields
                                %smaller and restart;
                                decrease = 1;
                                decrease_sub = decrease_sub + 1;
                                proceed = 0;
                                disp = 0;
                                m = 1;
                                n = 1;
                                count = 1;
                                layer(i).boundaries = {};
                                %anybad = 1;
                            elseif ~decrease
                                decrease_sub = 1;
                                proceed = 1;
                            end
                            
                            if proceed
                                %divide by two because we doubled the size of the shot map
                                %subtract 1/2 because we smeared out the
                                %shot map by 1/2 of an original pixel
                                %Subtract 1/2 again because of the way we
                                %doubled the size of the matrix.
                                %layer(i).boundaries(count)={B{b}/2-1/2-1/2+repmat([xstart,ystart],[size(B{b},1),1])};
                                
                                %For use with undoubled matrices. Subtract
                                %1 because of the way we did the smearing
                                layer(i).boundaries(count)={B{b}-1+repmat([xstart,ystart],[size(B{b},1),1])};
                                
                                
                                count=count+1;
                            end
                        end
                        
                        %Break out of looping over boundaries if there are
                        %large boundaries and we need to restart
                        if proceed==0
                            break
                        end
                    end
                end
            else
                proceed = 1;
            end
            if proceed
                m = m+1;
            end
        end
        if proceed
            n = n+1;
        end
    end
    
    
    fprintf('done.\n')
    
    
end
%fprintf('Fracturing complete. \n')

X=[];
Y=[];
BB=layer(1).boundaries;

for i=1:length(BB)
    if ~isempty(BB{i})
        b=BB{i};
        X=[X;b(:,2)];
        Y=[Y;b(:,1)];
    end
end


targetLen=100;
figure(1111); clf; hold on;
for i=1:length(BB)
    b=BB{i};
    b(:,2)=b(:,2)-mean(X);
    b(:,2)=b(:,2)./(abs(max(X)-min(X))).*targetLen;
    b(:,1)=b(:,1)-mean(Y);
    b(:,1)=b(:,1)./(abs(max(X)-min(X))).*targetLen;
    plot(b(:,2),-b(:,1));
    BB{i}=[b(:,2)+cntr(1) -b(:,1)+cntr(2)];
    
end
end

