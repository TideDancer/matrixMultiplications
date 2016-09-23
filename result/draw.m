
f = 'dim10_cond5_randsvd.logx';

fid = fopen(f);
s = fread(fid, '*char');
fclose(fid);

x = [];
t = [];

% get time
e = '[\d\.]* seconds';
matchStr = regexp(s',e,'match');
for i = 2:length(matchStr)
  e = '\s';
  v = regexp(char(matchStr(i)), e, 'split');
  t = [t; str2num(char(v(1)))];
end

% get accuracy
e = 'errorRatio =[\.\w]*';
matchStr = regexp(s',e,'match');
for i = 1:length(matchStr)
  e = '=';
  v = regexp(char(matchStr(i)), e, 'split');
  x = [x; str2num(char(v(2)))];
end

y = [x, t];

bar(y);
set(gca,'YScale','log');
grid on;
title(f);
legend('accuracy','time');
xlabel('algorithm #');
 

