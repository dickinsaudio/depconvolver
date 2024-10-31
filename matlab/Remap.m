function M = Remap_polar(Xp,Yp)

    Xp = mod(Xp+360,360);       % Move all angles into one circuit
    Yp = mod(Yp+360,360);
    
    [~, ox] = sort(Xp(1,:));  Xp = Xp(:,ox);     % Sort via angle
    [~, oy] = sort(Yp(1,:));  Yp = Yp(:,oy);

    M = zeros(size(Xp,2),size(Yp,2));

    % For now lets assume it is two rings of speakers - so split them

    E = unique(Xp(2,:));

    for (e = E)
        x = find(Xp(2,:)==e);         % Get the input and output speakers
        y = find(Yp(2,:)==e);         % at the given level

        if (sum(y)==0) 
            warning("No matching ring - mapping to median plane"); 
            y = find(Yp(2,:)==0); 
        end;

        if (sum(y)==0) warning("Empty speaker set"); continue; end;

        I  = Xp(1,x);           % Finally, just get the angles in and angles out
        O  = Yp(1,y);           %
        
        Oe = [ O(end)-360 O O(1)+360 ];       % Extended wrap to find the enclosing pair 

        m = zeros(length(I),length(O));

        for (i=1:length(I))
            if (~isempty(find(O==I(i)))) 
                M(x(i),y(find(O==I(i),1,'first'))) = 1;
            else
                o1 = sum(Oe<I(i));   d1 = abs(Oe(o1)-I(i));
                o2 = o1+1;           d2 = abs(Oe(o2)-I(i));
                o1 = o1 - 1; if (o1==0)         o1=length(O); end;
                o2 = o2 - 1; if (o2 >length(O)) o2=1; end;

                M(x(i),y(o1)) = d2 / (d1 + d2);
                M(x(i),y(o2)) = d1 / (d1 + d2);
            end;
        end;
    end;

    M(ox,:) = M;
    M(:,oy) = M;


