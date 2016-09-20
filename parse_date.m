function [ month day year ] = parse_date(date_fake)
%parses the date for use in URLS on bball reference

today=datestr(datenum(date_fake),'mmddyyyy');
if today(1)=='0'
    month=today(2);
else
    month=today(1:2);
end

if today(3)=='0'
    day=today(4);
else
    day=today(3:4);
end

year=today(5:8);

end

