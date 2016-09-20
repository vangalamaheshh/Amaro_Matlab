SouthEnd2012=load_struct('/Users/amaro/Documents/ApartmentHuntingPlots/DataWithSQfootage2012.txt');
SouthEnd2013=load_struct('/Users/amaro/Documents/ApartmentHuntingPlots/DataWithSQfootage2013.txt');
SouthEnd2014=load_struct('/Users/amaro/Documents/ApartmentHuntingPlots/DataWithSQfootage2014.txt');
SouthEnd2012.SalePrice=str2double(SouthEnd2012.SalePrice);
SouthEnd2012.ListPrice=str2double(SouthEnd2012.ListPrice);
SouthEnd2013.SalePrice=str2double(SouthEnd2013.SalePrice);
SouthEnd2013.ListPrice=str2double(SouthEnd2013.ListPrice);
SouthEnd2014.SalePrice=str2double(SouthEnd2014.SalePrice);
SouthEnd2014.ListPrice=str2double(SouthEnd2014.ListPrice);
SouthEnd2012.SquareFeet=str2double(SouthEnd2012.SquareFeet);
SouthEnd2013.SquareFeet=str2double(SouthEnd2013.SquareFeet);
SouthEnd2014.SquareFeet=str2double(SouthEnd2014.SquareFeet);
distribution_ordered_plot(SouthEnd2012.SalePrice./SouthEnd2012.SquareFeet,SouthEnd2013.SalePrice./SouthEnd2013.SquareFeet,SouthEnd2014.SalePrice./SouthEnd2014.SquareFeet)

set(gca,'XTick',[1:3],'XTickLabel',{'2012';'2013';'2014'},'FontSize',18)
ylabel('Price Per Square Foot','FontSize',18)
title('Price per Square foot 2012-2014 South End Boston','FontSize',24)


plot([SouthEnd2012.PSqSale],[SouthEnd2012.SalePrice],'ko')
hold on
plot([SouthEnd2013.PSqSale],[SouthEnd2013.SalePrice],'bo')
plot([SouthEnd2014.PSqSale],[SouthEnd2014.SalePrice],'ro')
xlabel('Price Per Square Feet','FontSize',16)
ylabel('Sale Price','FontSize',16)
legend({'2012','2013','2014'})

mean(SouthEnd2012.PSqSale)
mean(SouthEnd2013.PSqSale)
mean(SouthEnd2014.PSqSale)
mwwtest(SouthEnd2013.PSqSale,SouthEnd2014.PSqSale)
mwwtest(SouthEnd2012.PSqSale,SouthEnd2013.PSqSale)
(mean(SouthEnd2013.PSqSale)-mean(SouthEnd2012.PSqSale))./mean(SouthEnd2013.PSqSale);
(mean(SouthEnd2014.PSqSale)-mean(SouthEnd2013.PSqSale))./mean(SouthEnd2014.PSqSale);