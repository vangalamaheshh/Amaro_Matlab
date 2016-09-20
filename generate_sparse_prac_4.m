d = dataset('file', '~/Downloads/train-2.csv', 'ReadObsNames', false, 'ReadVarNames', true, 'Delimiter', ',');

uuser=unique(d.user);
uart=unique(d.artist);

M=zeros(length(uuser),length(uart));
for i=1313:length(uart)
    tic
    mask = ismember(d.artist,uart{i});
    M(ismember(uuser,d.user(mask)),i)=d.plays(mask);
    toc
    i
end