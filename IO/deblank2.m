function b           = deblank2(a,bAlsoTabs)
    if nargin<2
        bAlsoTabs = 0;
    end

    if ~isa(a,'char')
        b=a;
        return;
    end
   if numel(a)>0
       b = deblank(a);
   else
       b = a;
       return;
   end
   b = b(end:-1:1);
   if numel(b)>0
       b = deblank(b);
   end
   b = b(end:-1:1);
   %fprintf(1,'Line was %s, with %i characters. Now is %s, with %i characters.\n',a,length(a),b,length(b));
   if (bAlsoTabs)
       if(numel(b)>0)
           while strcmp(b(1),char(11))
               b(iChar) = [];           
           end
           while strcmp(b(end),char(11))
               b(end) = [];           
           end
       end
   end
end