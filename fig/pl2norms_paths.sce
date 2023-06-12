en_orig = fscanfMat("data/pl2norms_0.dat");
rw_orig = size(en_orig);
zero_index = 26;
en = en_orig(zero_index:rw_orig(1),:);
rw = size(en);
avg = zeros(rw(1),rw(2));
count = 0;
for num = 0
    en_orig = fscanfMat('data/pl2norms_'+string(num)+'.dat');
    en = en_orig(zero_index:rw_orig(1),:);
    avg = avg+en;
    count = count+1;
end

avg = avg/count;

figure(0);
m=size(avg,1);
rng=2:m;
clf(0)


plot2d('nl',avg(:,1), avg(:,3), style=2);
plot2d(avg(:,1), avg(:,4), style=3);
plot2d(avg(:,1), avg(:,5), style=4);
title('l2 norm of momenta');
legend('2h','4h','8h','16h', '32h', '64h', 'in_lower_right');

figure(1);
clf(1)

plot(avg(rng,1),log(avg(2:$,4:$)./avg(2:$,3:$-1))/log(2))
title('Approximate order of the schemes, versus time');
xlabel('Time');

f=gcf();
f.background=-2;
