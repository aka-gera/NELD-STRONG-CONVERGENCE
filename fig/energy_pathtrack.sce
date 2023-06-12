//Long term energy path
en = fscanfMat("data/energy_0.dat");
zero_index = 26;
rw = size(en);
avg = zeros(rw(1),rw(2));
count = 0;
for num = 0:12
    en = fscanfMat('data/energy_'+string(num)+'.dat');
    avg = avg+en;
    count = count+1;
end

avg = avg/count;

figure(0);
clf(0);

rep = 7;
rng = zero_index:rw(1);

subplot(1,2,1)
title('potential');
plot2d(avg(:,1), avg(:,2), style=3);
for i = 4:2:2*rep
    plot2d(avg(rng, 1), avg(rng, i));
end

subplot(1,2,2)
title('kinetic');
plot2d(avg(:,1), avg(:,3), style=3);
for i = 5:2:2*rep+1
    plot2d(avg(rng, 1), avg(rng, i));
end
