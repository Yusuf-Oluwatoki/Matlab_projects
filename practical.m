%% Testing the initial parameters with the model
clear all
close all
clc

cf_model = Simulink.SimulationInput('BeerGame_MOODLE_1') %First model 
%CUSTOMER DEMAND (constant)
c_demand(:,1)=1:365;
c_demand(:,2) = 3*ones(365,1);

% WHOLESALER (w) parameters
w.stock_r = 3;
w.order_r = 2

% RETAILER (r) parameters
r.stock_r = 2;
r.order_r = 1.5

%running the model and making necessary output for stock level and backlog
for i=1:1
test= sim(cf_model);
total_stock(:,i)= test.yout{1}.Values.Data;
total_backlog(:,i)= test.yout{2}.Values.Data;
wholesale_order(:,i)= test.yout{4}.Values.Data;
retail_order(:,i)= test.yout{3}.Values.Data;
wholesale_stock(:,i)= test.yout{5}.Values.Data;
retail_stock(:,i)= test.yout{6}.Values.Data;
end



%% Optimizing the paramenters to get the minimal backlog and storage size (Question1)
 
clear all
close all
clc

 cf_model = Simulink.SimulationInput('BeerGame_MOODLE_1')
%CUSTOMER DEMAND (constant)
c_demand(:,1)=1:365;
c_demand(:,2) = 3*ones(365,1);

% WHOLESALER (w) parameters
w.stock_r = 0.1:0.1:0.3;
w.order_r = 0:1:3;

% RETAILER (r) parameters
r.stock_r = 0.1:0.1:0.3;
r.order_r = 0:1:3;

%Using combvec for the ranges of parameters defined.
para= combvec(w.stock_r, w.order_r,r.stock_r,r.order_r);


%running the model and making necessary output for stock level and backlog
for i=1:length(para)

cf_model= cf_model.setVariable('w.stock_r', para(1,i));
cf_model= cf_model.setVariable('w.order_r', para(2,i));
cf_model= cf_model.setVariable('r.stock_r', para(3,i));
cf_model= cf_model.setVariable('r.order_r', para(4,i));
test= sim(cf_model);
total_stock(:,i)= test.yout{1}.Values.Data;
total_backlog(:,i)= test.yout{2}.Values.Data;
wholesale_order(:,i)= test.yout{4}.Values.Data;
retail_order(:,i)= test.yout{3}.Values.Data;
wholesale_stock(:,i)= test.yout{5}.Values.Data;
retail_stock(:,i)= test.yout{6}.Values.Data;

average_backlog(i)=mean(total_backlog(:,i));
average_stock(i)=mean(total_stock(:,i));
max_stock(i)=max(total_stock(:,i));
max_retail_stock(i)=max(retail_stock(:,i));
max_wholesale_stock(i)=max(wholesale_stock(:,i));

end

figure
plot(average_backlog, average_stock, '*')
hold on


%Finding the optimized parameters for the above output
mini_backlog=min(average_backlog)
minimized_backlog= average_stock(find(average_backlog==min(average_backlog)));
minimized_stock_index= find(average_stock==min(minimized_backlog))
minimized_stock=average_stock(minimized_stock_index)
opitmal_parameters= para(:,minimized_stock_index)

%Highligthing the optimal stock level in the plot.
plot(mini_backlog, minimized_stock, '*r', 'MarkerSize', 10)
legend('others', 'optimum')
xlabel('Backlog')
ylabel('Storage size')
title('Backlog against Storage size')

%% Building the economic model with the optimized parameters (Question 2)

clc

 cf_model = Simulink.SimulationInput('BeerGame_MOODLE')
%CUSTOMER DEMAND (constant)
c_demand(:,1)=1:365;
c_demand(:,2) = 3*ones(365,1);

% WHOLESALER (w) parameters
w.stock_r = 0.2;
w.order_r = 1;

% RETAILER (r) parameters
r.stock_r = 0.3;
r.order_r = 1;


%Creating the Fixed Expenses vector.

salary(:,1)= [1:365];
salary([15:30:365],2)=5000;

Raw_material1(:,1)= [1:365];
Raw_material1([1:30:360],2)=4000;

Raw_material2(:,1)= [1:365];
Raw_material2([7:30:365],2)=3000;

Electricity(:,1)= [1:365];
Electricity([15:30:365],2)=2000;

Maintenance(:,1)= [1:365];
Maintenance([30:30:365],2)=2000;

Others(:,1)= [1:365];
Others([12:30:365],2)=1000;

Investment_x=zeros(365,1);
Investment_x(180)=250000;

Investment_y=zeros(365,1);
Investment_y(210)=500000;

Loan=zeros(365,1);
Loan(260)=250000;

Expenses(:,1)= [1:365];
Expenses(:,2)=Others(:,2)+ Maintenance(:,2)+ Electricity(:,2)+ salary(:,2)+...
    Raw_material1(:,2)+ Raw_material2(:,2)+ Loan+ Investment_x+ Investment_y;



%running the model and making necessary output for stock level, backlog and
%financial outcome
for i=1:1
test= sim(cf_model);
total_stock(:,i)= test.yout{1}.Values.Data;
total_backlog(:,i)= test.yout{2}.Values.Data;
wholesale_order(:,i)= test.yout{4}.Values.Data;
retail_order(:,i)= test.yout{3}.Values.Data;
wholesale_stock(:,i)= test.yout{5}.Values.Data;
retail_stock(:,i)= test.yout{6}.Values.Data;
revenue(:,i)= test.yout{7}.Values.Data;

average_backlog(i)=mean(total_backlog(:,i));
average_stock(i)=mean(total_stock(:,i));
max_stock(i)=max(total_stock(:,i));
max_retail_stock(i)=max(retail_stock(:,i));
max_wholesale_stock(i)=max(wholesale_stock(:,i));
end 

%Plot of the financial outcome.

Profit=revenue;
figure
plot(revenue, '*b')
xlabel('Days')
ylabel('Profit')
title('Daily profit')
%% Optimization with historical demand data (Question 3)

clc
load Demand5yr.mat
 cf_model = Simulink.SimulationInput('BeerGame_MOODLE')

% WHOLESALER (w) parameters
w.stock_r = 0.2;
w.order_r = 1;

% RETAILER (r) parameters
r.stock_r = 0.3;
r.order_r = 1;


%Creating the Fixed Expenses vector.

salary(:,1)= [1:365];
salary([15:30:365],2)=5000;

Raw_material1(:,1)= [1:365];
Raw_material1([1:30:360],2)=4000;

Raw_material2(:,1)= [1:365];
Raw_material2([7:30:365],2)=3000;

Electricity(:,1)= [1:365];
Electricity([15:30:365],2)=2000;

Maintenance(:,1)= [1:365];
Maintenance([30:30:365],2)=2000;

Others(:,1)= [1:365];
Others([12:30:365],2)=1000;

Investment_x=zeros(365,1);
Investment_x(180)=250000;

Investment_y=zeros(365,1);
Investment_y(210)=500000;

Loan=zeros(365,1);
Loan(260)=250000;

Expenses(:,1)= [1:365];
Expenses(:,2)=Others(:,2)+ Maintenance(:,2)+ Electricity(:,2)+ salary(:,2)+...
    Raw_material1(:,2)+ Raw_material2(:,2)+ Loan+ Investment_x+ Investment_y;


%running the model and making necessary output for stock level, backlog and
%financial outcome

figure
for i=1:5
cf_model= cf_model.setVariable('c_demand', [[0:365]' Demand5yr(:,i)]);
test= sim(cf_model);
total_stock(:,i)= test.yout{1}.Values.Data;
total_backlog(:,i)= test.yout{2}.Values.Data;
wholesale_order(:,i)= test.yout{4}.Values.Data;
retail_order(:,i)= test.yout{3}.Values.Data;
wholesale_stock(:,i)= test.yout{5}.Values.Data;
retail_stock(:,i)= test.yout{6}.Values.Data;
revenue(:,i)= test.yout{7}.Values.Data;

average_backlog(i)=mean(total_backlog(:,i));
average_stock(i)=mean(total_stock(:,i));
max_stock(i)=max(total_stock(:,i));
max_retail_stock(i)=max(retail_stock(:,i));
max_wholesale_stock(i)=max(wholesale_stock(:,i));

plot(average_backlog(i), average_stock(i), '*')
hold on
end
plot(mini_backlog, minimized_stock, '*r', 'MarkerSize', 10)
legend('Year 1', 'Year 2', 'Year 3', 'Year 4', 'Year 5','Constant-demand')
xlabel('Backlog')
ylabel('Storage size')
title('Backlog against Storage size')


figure
plot(revenue)
hold on
plot(Profit)
legend('Year 1', 'Year 2', 'Year 3', 'Year 4', 'Year 5', 'Constant-demand')
xlabel('Days')
ylabel('Profits')
title('Daily profit')
%% Advertising campaign with triangular distribution (Question 4)


clc
load Demand5yr.mat
 cf_model = Simulink.SimulationInput('BeerGame_MOODLE')

% WHOLESALER (w) parameters
w.stock_r = 0.2;
w.order_r = 1;

% RETAILER (r) parameters
r.stock_r = 0.3;
r.order_r = 1;


%Creating the Fixed Expenses vector.

salary(:,1)= [1:365];
salary([15:30:365],2)=5000;

Raw_material1(:,1)= [1:365];
Raw_material1([1:30:360],2)=4000;

Raw_material2(:,1)= [1:365];
Raw_material2([7:30:365],2)=3000;

Electricity(:,1)= [1:365];
Electricity([15:30:365],2)=2000;

Maintenance(:,1)= [1:365];
Maintenance([30:30:365],2)=2000;

Others(:,1)= [1:365];
Others([12:30:365],2)=1000;

Investment_x=zeros(365,1);
Investment_x(180)=250000;

Investment_y=zeros(365,1);
Investment_y(210)=500000;

Loan=zeros(365,1);
Loan(260)=250000;

Advert=zeros(365,1);
Advert(190)=100000;

Expenses(:,1)= [1:365];
Expenses(:,2)=Others(:,2)+ Maintenance(:,2)+ Electricity(:,2)+ salary(:,2)+...
    Raw_material1(:,2)+ Raw_material2(:,2)+ Loan+ Advert +Investment_x+ Investment_y;


%Creating the forecast demand vector
For_demand(:,1)=0:365;
For_demand(:,2)=mean(Demand5yr,2);



%Using the triangular distribution to effect the advert on forecast demand
values=zeros(366,1);

for j=1:3
days=[ 150, 180, 210];
para=[ 0.5 1 1.5; 0.5, 1.3, 1.8; 0.5, 0.8, 0.9]; 
pd_demand=makedist('Triangular', 'a',para(j,1), 'b',para(j,2),'c',para(j,3));
values(days(j))= random(pd_demand,1,1);
end
values(151:179)=values(150);
values(181:209)=values(150)+values(180);
values(211:end)=values(150)+values(180)+values(210);
values(180)=values(181); values(210)=values(211); 
For_demand(:,2)= For_demand(:,2)+ values.*For_demand(:,2);

%Plot of forecast demand
figure
plot(For_demand(:,2))
xlabel('Days')
ylabel('Demands')
title('Increase in Demand with Advert')


%running the model and making necessary output for stock level, backlog and
%financial outcome

cf_model= cf_model.setVariable('c_demand', For_demand);
test= sim(cf_model);
total_stock= test.yout{1}.Values.Data;
total_backlog= test.yout{2}.Values.Data;
wholesale_order= test.yout{4}.Values.Data;
retail_order= test.yout{3}.Values.Data;
wholesale_stock= test.yout{5}.Values.Data;
retail_stock= test.yout{6}.Values.Data;
advert= test.yout{7}.Values.Data;

aver_backlog=mean(total_backlog);
aver_stock=mean(total_stock);
max_stock=max(total_stock);
max_retail_stock=max(retail_stock);
max_wholesale_stock=max(wholesale_stock);


%Plot of stock levels against backlog to evaluate model performance in the
%three cases ( 5 years, advert and constant demand)
figure
plot(average_backlog, average_stock./2, '*')
hold on
plot(mini_backlog, minimized_stock, '*', 'MarkerSize', 10)
hold on
plot(aver_backlog, aver_stock, '*','MarkerSize', 20 )
legend('Yearly','Constant-demand','advert')
xlabel('Backlog')
ylabel('Storage size')
title('Backlog against Storage size')


%Plot of revenue to evaluate model performance in the
%three cases ( 5 years, advert and constant demand)

figure
plot(revenue)
hold on
plot(Profit)
hold on
plot(advert)
legend('Year 1', 'Year 2', 'Year 3', 'Year 4', 'Year 5', 'Constant-demand','advert')
xlabel('Days')
ylabel('Profits')
title('Daily profit')