
zero_index = 26; 
avg = 0;  
count = 0;
for num = 0 
    en_orig = fscanfMat('data/energy_diff_'+string(num)+'.dat');
    en = en_orig(zero_index:size(en_orig,1),:);
    avg = avg+en;
    count = count+1;
end

avg = avg/count;

m=size(en,1);
rng=3:m;
rep = 4;

figure(0);
clf(0);
subplot(1,2,1)
title('potential');
plot2d(avg(:,1), avg(:,2:2:2*rep+1));

subplot(1,2,2)
title('kinetic');
plot2d(avg(:,1), avg(:,3:2:2*rep+1));

figure(1);
clf(1);

subplot(1,2,1)
uhp = diff(avg(rng,4:2:2*rep+1),1,'c');
plot(avg(rng,1),log(uhp(:,2:$)./uhp(:,1:$-1))/log(2))
title('Approximate order of the schemes, versus time potential');
xlabel('Time');

subplot(1,2,2)
uhk = diff(avg(rng,5:2:2*rep+1),1,'c');
plot(avg(rng,1),log(uhk(:,2:$)./uhk(:,1:$-1))/log(2))
title('Approximate order of the schemes, versus time kinetic');
xlabel('Time');


f=gcf();
f.background=-2;
