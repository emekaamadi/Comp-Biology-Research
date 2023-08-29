function [width,max_apm,height] = get_shapedata(doseresponse)


[M,I] = max(doseresponse(:,1));
height = doseresponse(I,1);
% height = maxdeadc
max_apm = doseresponse(I,2);
% Max= the concentration of apm at the max cancer cells 
a_pm = logspace(-6,6,length(doseresponse));
x = 1;
a_pm1 = doseresponse(1,2);
a_pm2 = doseresponse(length(doseresponse),2);
%return the index of the apm concentration when the dead aml conc is >= max/ 2
while(x <= length(a_pm))
    if (doseresponse(x,1) >= height/2)
        a_pm1 = doseresponse(x,2);
        x = length(a_pm)+1;
    end
    x = x+1;
end
%return the index of the apm concentration when the dead aml conc is >= max/ 2
x = I;
while(x <= length(a_pm))
    if (doseresponse(x,1) <= height/2)
        a_pm2 = doseresponse(x,2);
        x = length(a_pm)+1;
    end
    x = x+1;
end
width = log10(a_pm2/a_pm1);
end
