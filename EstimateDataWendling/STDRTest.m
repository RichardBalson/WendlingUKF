x= rand(1e3,1)+1;
% x(length(x)+1) = 50;
mu = mean(x);
stdx = std(x);


stdr =0;
mur =0;
for k = 1:length(x)
murT = (mur*(k-1)+x(k))/k;
if k >1
stdr = sqrt(((k-1)*(stdr^2+mur^2)+x(k)^2-(k)*murT^2)/(k));
end
mur = murT;
end