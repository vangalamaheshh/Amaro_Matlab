function i = make_int32(b)
% assume little-endian
b = double(b);
i = b(4);
i = bitshift(i,8)+b(3);
i = bitshift(i,8)+b(2);
i = bitshift(i,8)+b(1);
