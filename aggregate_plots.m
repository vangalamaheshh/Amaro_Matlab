function aggregate_plots( basename, cropim, addsep, resize )




listing = dir;
listing(1:2) = [];
ims = {};
sepwidth = 10;
for i = 1 : length(listing)
    if length(listing(i).name) > length(basename) && ...
            strcmpi(listing(i).name(1:length(basename)), basename)
        
        im = imread(listing(i).name);
        if cropim
            imc = im(50:550, 150:650, :);
            else
            imc = im;
            end
        
        if addsep
            imc = [zeros(size(imc,1), sepwidth, 3), imc, zeros(size(imc,1), sepwidth, 3)];
            imc = [zeros(sepwidth, size(imc,2), 3); imc; zeros(sepwidth, size(imc,2), 3)];
            end
        
        ims{end+1} = imc;
        
        end
end

f = factor(length(ims));
sp = floor(length(f)/2);

a = prod(f(1:sp));
b = prod(f(sp+1:end));

imsr = reshape(ims,a,b);
imwrite(imresize(cell2mat(imsr), resize), ['agg_', basename, '.png']);

end