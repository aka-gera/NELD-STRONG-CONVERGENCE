en_orig = fscanfMat("data/energy_diff_0.dat");
rw_orig = size(en_orig);
zero_index = 26;
en = en_orig(zero_index:rw_orig(1),:);
rw = size(en);
avg2 = zeros(rw(1),rw(2));
count = 0;
for num = 1 
        en_orig = fscanfMat('data/energy_diff_'+string(num)+'.dat');
        en = en_orig(zero_index:rw_orig(1),:);
        avg2 = avg2+en;
        count = count+1;
end

avg2 = avg2/count;


m=size(avg2,1);
rng=3:m; 
rep = 4;

row = m
col = size(avg2,2)
for i = 2:col
    avg(:,i-1) = cumsum(avg2(:,i))
end
for i = 1: row
    avg(i,:) = avg(i,:)/i
end 


figure(0);
clf(0);
subplot(1,2,1)
title('potential');
plot2d(avg2(:,1), avg(:,3:2:2*rep));

subplot(1,2,2)
title('kinetic');
plot2d(avg2(:,1), avg(:,4:2:2*rep));

figure(1);
clf(1);

subplot(1,2,1)
uhp = diff(avg(rng,3:2:2*rep),1,'c');
plot(avg2(rng,1),log(uhp(:,2:$)./uhp(:,1:$-1))/log(2))
title('Approximate order of the schemes, versus time potential');
xlabel('Time');
legend('h','2h','4h','8h','in_upper_left')

subplot(1,2,2)
uhk = diff(avg(rng,4:2:2*rep),1,'c');
plot(avg2(rng,1),log(uhk(:,2:$)./uhk(:,1:$-1))/log(2))
title('Approximate order of the schemes, versus time kinetic');
xlabel('Time');
legend('h','2h','4h','8h','in_lower_left')
