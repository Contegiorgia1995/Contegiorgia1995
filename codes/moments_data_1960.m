function [m_data,m_title,inc] = moments_data_1960
%% Moments from the data
Nm            = 1019;     
m_title       = cell(Nm,1);
m_data        = nan(Nm,1);
m_id          = nan(Nm,1);
inc           = zeros(Nm,1);

%% TARGET MOMENTS
pos = 001; m_title{pos}   ='Mean income                           '; inc(pos) = 0; m_id(pos) = pos;
pos = 002; m_title{pos}   ='CV income                             '; inc(pos) = 0; m_id(pos) = pos; 
pos = 003; m_title{pos}   ='Dropouts                              '; inc(pos) = 1; m_id(pos) = pos;
pos = 004; m_title{pos}   ='High school graduates                 '; inc(pos) = 1; m_id(pos) = pos;
pos = 005; m_title{pos}   ='College graduates                     '; inc(pos) = 1; m_id(pos) = pos; 
pos = 006; m_title{pos}   ='Mean fertility                        '; inc(pos) = 1; m_id(pos) = pos; 
pos = 007; m_title{pos}   ='Fertility elasticity  (q10)           '; inc(pos) = 1; m_id(pos) = pos; 
pos = 008; m_title{pos}   ='Intergenerational Mobility: Rank-Rank '; inc(pos) = 1; m_id(pos) = pos; 
pos = 009; m_title{pos}   ='Intergenerational Mobility: Log-Log   '; inc(pos) = 0; m_id(pos) = pos; 

%% NON TARGETED MOMENTS
% Income
pos = 010; m_title{pos}   ='Income top-bottom 25-64               '; inc(pos) = 1; m_id(pos) = pos; 
pos = 011; m_title{pos}   ='Income gini 25-64                     '; inc(pos) = 1; m_id(pos) = pos; 
pos = 012; m_title{pos}   ='IGE Education Persistence: Log-Log    '; inc(pos) = 0; m_id(pos) = pos; 
pos = 013; m_title{pos}   ='IGE Education Persistence: 2nd Eig    '; inc(pos) = 0; m_id(pos) = pos; 
pos = 014; m_title{pos}   ='IGE Education Persistence: Trace      '; inc(pos) = 0; m_id(pos) = pos; 
pos = 015; m_title{pos}   ='IGE Education Persistence: Det        '; inc(pos) = 0; m_id(pos) = pos; 

pos = 016; m_title{pos}   ='Intergenerational Mobility: Rank-Rank (100)'; inc(pos) = 0; m_id(pos) = pos; 

pos = 017; m_title{pos}   ='Theil-L index absolute (par inc educ) '; inc(pos) = 0; m_id(pos) = pos;
pos = 018; m_title{pos}   ='Theil-L index relative (par inc educ) '; inc(pos) = 0; m_id(pos) = pos;

% Mean Fertility by educ
pos = 019; m_title{pos}   ='Mean Fert Educ 1                      '; inc(pos) = 0; m_id(pos) = pos;      
pos = 020; m_title{pos}   ='Mean Fert Educ 2                      '; inc(pos) = 0; m_id(pos) = pos;     
pos = 021; m_title{pos}   ='Mean Fert Educ 3                      '; inc(pos) = 0; m_id(pos) = pos;     

pos = 022; m_title{pos}   ='Theil-L index absolute (InCond, inc)  '; inc(pos) = 0; m_id(pos) = pos;
pos = 023; m_title{pos}   ='Theil-L index relative (InCond, inc)  '; inc(pos) = 0; m_id(pos) = pos;

pos = 024; m_title{pos}   ='Theil-L index absolute (parents,LE)   '; inc(pos) = 0; m_id(pos) = pos;
pos = 025; m_title{pos}   ='Theil-L index relative (parents,LE)   '; inc(pos) = 0; m_id(pos) = pos;

% Fertility
pos = 026; m_title{pos}   ='E(fertility$\mid$Q01)                 '; inc(pos) = 0; m_id(pos) = pos;      
pos = 027; m_title{pos}   ='E(fertility$\mid$Q02)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 028; m_title{pos}   ='E(fertility$\mid$Q03)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 029; m_title{pos}   ='E(fertility$\mid$Q04)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 030; m_title{pos}   ='E(fertility$\mid$Q05)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 031; m_title{pos}   ='E(fertility$\mid$Q06)                 '; inc(pos) = 0; m_id(pos) = pos;      
pos = 032; m_title{pos}   ='E(fertility$\mid$Q07)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 033; m_title{pos}   ='E(fertility$\mid$Q08)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 034; m_title{pos}   ='E(fertility$\mid$Q09)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 035; m_title{pos}   ='E(fertility$\mid$Q10)                 '; inc(pos) = 0; m_id(pos) = pos;     
pos = 036; m_title{pos}   ='Pr(fertility = 0)                     '; inc(pos) = 0; m_id(pos) = pos;        
pos = 037; m_title{pos}   ='Pr(fertility = 1)                     '; inc(pos) = 0; m_id(pos) = pos;        
pos = 038; m_title{pos}   ='Pr(fertility = 2)                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 039; m_title{pos}   ='Pr(fertility $\geq$3)                 '; inc(pos) = 0; m_id(pos) = pos;     

% Retirement
pos = 040; m_title{pos}   ='Retirement gov E1                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 041; m_title{pos}   ='Retirement gov E2                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 042; m_title{pos}   ='Retirement gov E3                     '; inc(pos) = 0; m_id(pos) = pos;     

pos = 045; m_title{pos}   ='Retirement transfers E1               '; inc(pos) = 0; m_id(pos) = pos;      
pos = 046; m_title{pos}   ='Retirement transfers E2               '; inc(pos) = 0; m_id(pos) = pos; 
pos = 047; m_title{pos}   ='Retirement transfers E3               '; inc(pos) = 0; m_id(pos) = pos; 

pos = 050; m_title{pos}   ='Retirement savings E1                 '; inc(pos) = 0; m_id(pos) = pos; 
pos = 051; m_title{pos}   ='Retirement savings E2                 '; inc(pos) = 0; m_id(pos) = pos; 
pos = 052; m_title{pos}   ='Retirement savings E3                 '; inc(pos) = 0; m_id(pos) = pos; 

% Nanny
pos = 055; m_title{pos}   ='Hire childcare                        '; inc(pos) = 0; m_id(pos) = pos;   

% Education
pos = 056; m_title{pos}   ='Variance years of education           '; inc(pos) = 0; m_id(pos) = pos;  
pos = 057; m_title{pos}   ='     $\%$  expl. by initial HK        '; inc(pos) = 0; m_id(pos) = pos;  
pos = 058; m_title{pos}   ='     $\%$  expl. by transfers         '; inc(pos) = 0; m_id(pos) = pos;
pos = 059; m_title{pos}   ='     $\%$  expl. by school taste      '; inc(pos) = 0; m_id(pos) = pos;
pos = 060; m_title{pos}   ='     $\%$  expl. by interaction       '; inc(pos) = 0; m_id(pos) = pos; 

% Transition Matrix
pos = 061; m_title{pos}   ='Transition Par Q1-Child Q1            '; inc(pos) = 0; m_id(pos) = pos;
pos = 062; m_title{pos}   ='Transition Par Q1-Child Q2            '; inc(pos) = 0; m_id(pos) = pos;
pos = 063; m_title{pos}   ='Transition Par Q1-Child Q3            '; inc(pos) = 0; m_id(pos) = pos;
pos = 064; m_title{pos}   ='Transition Par Q1-Child Q4            '; inc(pos) = 0; m_id(pos) = pos;
pos = 065; m_title{pos}   ='Transition Par Q1-Child Q5            '; inc(pos) = 0; m_id(pos) = pos;
pos = 066; m_title{pos}   ='Transition Par Q2-Child Q1            '; inc(pos) = 0; m_id(pos) = pos;
pos = 067; m_title{pos}   ='Transition Par Q2-Child Q2            '; inc(pos) = 0; m_id(pos) = pos;
pos = 068; m_title{pos}   ='Transition Par Q2-Child Q3            '; inc(pos) = 0; m_id(pos) = pos;
pos = 069; m_title{pos}   ='Transition Par Q2-Child Q4            '; inc(pos) = 0; m_id(pos) = pos;
pos = 070; m_title{pos}   ='Transition Par Q2-Child Q5            '; inc(pos) = 0; m_id(pos) = pos;
pos = 071; m_title{pos}   ='Transition Par Q3-Child Q1            '; inc(pos) = 0; m_id(pos) = pos;
pos = 072; m_title{pos}   ='Transition Par Q3-Child Q2            '; inc(pos) = 0; m_id(pos) = pos;
pos = 073; m_title{pos}   ='Transition Par Q3-Child Q3            '; inc(pos) = 0; m_id(pos) = pos;
pos = 074; m_title{pos}   ='Transition Par Q3-Child Q4            '; inc(pos) = 0; m_id(pos) = pos;
pos = 075; m_title{pos}   ='Transition Par Q3-Child Q5            '; inc(pos) = 0; m_id(pos) = pos;
pos = 076; m_title{pos}   ='Transition Par Q4-Child Q1            '; inc(pos) = 0; m_id(pos) = pos;
pos = 077; m_title{pos}   ='Transition Par Q4-Child Q2            '; inc(pos) = 0; m_id(pos) = pos;
pos = 078; m_title{pos}   ='Transition Par Q4-Child Q3            '; inc(pos) = 0; m_id(pos) = pos;
pos = 079; m_title{pos}   ='Transition Par Q4-Child Q4            '; inc(pos) = 0; m_id(pos) = pos;
pos = 080; m_title{pos}   ='Transition Par Q4-Child Q5            '; inc(pos) = 0; m_id(pos) = pos;
pos = 081; m_title{pos}   ='Transition Par Q5-Child Q1            '; inc(pos) = 0; m_id(pos) = pos;
pos = 082; m_title{pos}   ='Transition Par Q5-Child Q2            '; inc(pos) = 0; m_id(pos) = pos;
pos = 083; m_title{pos}   ='Transition Par Q5-Child Q3            '; inc(pos) = 0; m_id(pos) = pos;
pos = 084; m_title{pos}   ='Transition Par Q5-Child Q4            '; inc(pos) = 0; m_id(pos) = pos;
pos = 085; m_title{pos}   ='Transition Par Q5-Child Q5            '; inc(pos) = 0; m_id(pos) = pos;

% Aggregate Cons vs Output
pos = 086; m_title{pos}   ='Goods Exp/Output                   '; inc(pos) = 0; m_id(pos) = pos;
pos = 087; m_title{pos}   ='Education Exp/Output               '; inc(pos) = 0; m_id(pos) = pos;
pos = 088; m_title{pos}   ='Nannies Exp/Output                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 089; m_title{pos}   ='Total Consumption/Output           '; inc(pos) = 0; m_id(pos) = pos;
pos = 090; m_title{pos}   ='Depreciation rate                  '; inc(pos) = 0; m_id(pos) = pos;

% Government Budget Balance
pos = 091; m_title{pos}   ='Ret. Benefits - Tax Revenue        '; inc(pos) = 0; m_id(pos) = pos;
pos = 092; m_title{pos}   ='Ret. Benefits/Tax Revenue          '; inc(pos) = 0; m_id(pos) = pos;

% Inequality of Opportunity
pos = 093; m_title{pos}   ='Theil-L index total (age 28 - 31)  '; inc(pos) = 0; m_id(pos) = pos;
pos = 094; m_title{pos}   ='Theil-L index absolute (parents inc)'; inc(pos) = 0; m_id(pos) = pos;
pos = 095; m_title{pos}   ='Theil-L index relative (parents inc, inc)'; inc(pos) = 0; m_id(pos) = pos; % This one is not to be used in main results. There is no parents' education.
pos = 096; m_title{pos}   ='Theil-L index absolute (initial cond)'; inc(pos) = 0; m_id(pos) = pos;
pos = 097; m_title{pos}   ='Theil-L index relative (InCond, LE)'; inc(pos) = 0; m_id(pos) = pos;

% Transfers to children
pos = 098; m_title{pos}   ='Mean transfer to children          '; inc(pos) = 0; m_id(pos) = pos; 
pos = 099; m_title{pos}   ='CV transfers to children           '; inc(pos) = 0; m_id(pos) = pos;
pos = 100; m_title{pos}   ='Mean (cond $>$ 0) transfer to children'; inc(pos) = 0; m_id(pos) = pos;
pos = 101; m_title{pos}   ='CV (cond $>$ 0) transfers to children'; inc(pos) = 0; m_id(pos) = pos;

% GDP per Capita
pos = 102; m_title{pos}   ='GDP pc Consumption                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 103; m_title{pos}   ='GDP pc Income                      '; inc(pos) = 0; m_id(pos) = pos;
pos = 104; m_title{pos}   ='K/Y                                '; inc(pos) = 0; m_id(pos) = pos; 

pos = 105; m_title{pos}   ='Total transfers to children        '; inc(pos) = 0; m_id(pos) = pos; 
pos = 106; m_title{pos}   ='Transfers to children (share of savings)'; inc(pos) = 0; m_id(pos) = pos; 

pos = 107; m_title{pos}   ='Parent educ = 1 inc Q10, child avg inc'; inc(pos) = 0; m_id(pos) = pos; 
pos = 108; m_title{pos}   ='Parent educ = 1 inc Q10, child avg inc rank'; inc(pos) = 0; m_id(pos) = pos; 
pos = 109; m_title{pos}   ='Parent educ = 1 inc Q10, child educ 1'; inc(pos) = 0; m_id(pos) = pos; 
pos = 110; m_title{pos}   ='Parent educ = 1 inc Q10, child educ 2'; inc(pos) = 0; m_id(pos) = pos; 
pos = 111; m_title{pos}   ='Parent educ = 1 inc Q10, child educ 3'; inc(pos) = 0; m_id(pos) = pos; 
pos = 112; m_title{pos}   ='Parent educ = 1 inc Q10, max parent income'; inc(pos) = 0; m_id(pos) = pos; 

%% Mean Income by Educ and Age
pos = 182; m_title{pos}  ='Mean inc HS Dropout                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 183; m_title{pos}  ='Mean inc HS Graduates               '; inc(pos) = 0; m_id(pos) = pos;
pos = 184; m_title{pos}  ='Mean inc College Graduates          '; inc(pos) = 0; m_id(pos) = pos;

pos = 185; m_title{pos}  ='CV inc HS Dropout                   '; inc(pos) = 0; m_id(pos) = pos;
pos = 186; m_title{pos}  ='CV inc HS Graduates                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 187; m_title{pos}  ='CV inc College Graduates            '; inc(pos) = 0; m_id(pos) = pos;

%%
pos = 212; m_title{pos}  ='Parent Educ 1 Child Educ 1          '; inc(pos) = 0; m_id(pos) = pos;
pos = 213; m_title{pos}  ='Parent Educ 1 Child Educ 2          '; inc(pos) = 0; m_id(pos) = pos;
pos = 214; m_title{pos}  ='Parent Educ 1 Child Educ 3          '; inc(pos) = 0; m_id(pos) = pos;

pos = 215; m_title{pos}  ='Parent Educ 2 Child Educ 1          '; inc(pos) = 0; m_id(pos) = pos;
pos = 216; m_title{pos}  ='Parent Educ 2 Child Educ 2          '; inc(pos) = 0; m_id(pos) = pos;
pos = 217; m_title{pos}  ='Parent Educ 2 Child Educ 3          '; inc(pos) = 0; m_id(pos) = pos;

pos = 218; m_title{pos}  ='Parent Educ 3 Child Educ 1          '; inc(pos) = 0; m_id(pos) = pos;
pos = 219; m_title{pos}  ='Parent Educ 3 Child Educ 2          '; inc(pos) = 0; m_id(pos) = pos;
pos = 220; m_title{pos}  ='Parent Educ 3 Child Educ 3          '; inc(pos) = 0; m_id(pos) = pos;


%%
% Retirement
pos = 245; m_title{pos}   ='Retirement gov Q1                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 246; m_title{pos}   ='Retirement gov Q2                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 247; m_title{pos}   ='Retirement gov Q3                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 248; m_title{pos}   ='Retirement gov Q4                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 249; m_title{pos}   ='Retirement gov Q5                     '; inc(pos) = 0; m_id(pos) = pos;      

pos = 254; m_title{pos}   ='Retirement OAS Q1                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 255; m_title{pos}   ='Retirement OAS Q2                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 256; m_title{pos}   ='Retirement OAS Q3                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 257; m_title{pos}   ='Retirement OAS Q4                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 258; m_title{pos}   ='Retirement OAS Q5                     '; inc(pos) = 0; m_id(pos) = pos;      

pos = 263; m_title{pos}   ='Retirement sav Q1                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 264; m_title{pos}   ='Retirement sav Q2                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 265; m_title{pos}   ='Retirement sav Q3                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 266; m_title{pos}   ='Retirement sav Q4                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 267; m_title{pos}   ='Retirement sav Q5                     '; inc(pos) = 0; m_id(pos) = pos;      

% Fertility
pos = 272; m_title{pos}   ='E(fertility$\mid$Q01 incl. S          '; inc(pos) = 0; m_id(pos) = pos;      
pos = 273; m_title{pos}   ='E(fertility$\mid$Q02 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 274; m_title{pos}   ='E(fertility$\mid$Q03 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 275; m_title{pos}   ='E(fertility$\mid$Q04 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 276; m_title{pos}   ='E(fertility$\mid$Q05 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 277; m_title{pos}   ='E(fertility$\mid$Q06 incl. S          '; inc(pos) = 0; m_id(pos) = pos;      
pos = 278; m_title{pos}   ='E(fertility$\mid$Q07 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 279; m_title{pos}   ='E(fertility$\mid$Q08 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 280; m_title{pos}   ='E(fertility$\mid$Q09 incl. S          '; inc(pos) = 0; m_id(pos) = pos;     
pos = 281; m_title{pos}   ='E(fertility$\mid$Q10 incl. S          '; inc(pos) = 0; m_id(pos) = pos;  
pos = 282; m_title{pos}   ='Fertility elasticity incl. S          '; inc(pos) = 0; m_id(pos) = pos;  

% Lifetime Earnings
pos = 283; m_title{pos}  ='CV of Lifetime Earnings             '; inc(pos) = 0; m_id(pos) = pos;  
pos = 284; m_title{pos}  ='     $\%$  expl. by initial conds   '; inc(pos) = 0; m_id(pos) = pos;  
pos = 285; m_title{pos}  ='     $\%$  expl. by initial HK      '; inc(pos) = 0; m_id(pos) = pos; 
pos = 286; m_title{pos}  ='     $\%$  expl. by transfers       '; inc(pos) = 0; m_id(pos) = pos; 
pos = 287; m_title{pos}  ='     $\%$  expl. by school taste    '; inc(pos) = 0; m_id(pos) = pos; 
pos = 288; m_title{pos}  ='     $\%$  expl. by interaction     '; inc(pos) = 0; m_id(pos) = pos;

% Aggregates
pos = 289; m_title{pos}  ='Aggregate Capital                   '; inc(pos) = 0; m_id(pos) = pos;
pos = 290; m_title{pos}  ='Aggregate H HS drop                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 291; m_title{pos}  ='Aggregate H HS grad                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 292; m_title{pos}  ='Aggregate H Collge grad             '; inc(pos) = 0; m_id(pos) = pos;
pos = 293; m_title{pos}  ='Gov exp SS                          '; inc(pos) = 0; m_id(pos) = pos;

pos = 294; m_title{pos}  ='Implied/guess  r                    '; inc(pos) = 0; m_id(pos) = pos;
pos = 295; m_title{pos}  ='Implied        w(1)                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 296; m_title{pos}  ='Implied        w(2)                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 297; m_title{pos}  ='Implied        w(3)                 '; inc(pos) = 0; m_id(pos) = pos;

% Fertility
pos = 298; m_title{pos}   ='E(fertility$\mid$Q1)                  '; inc(pos) = 0; m_id(pos) = pos;      
pos = 299; m_title{pos}   ='E(fertility$\mid$Q2)                  '; inc(pos) = 0; m_id(pos) = pos;     
pos = 300; m_title{pos}   ='E(fertility$\mid$Q3)                  '; inc(pos) = 0; m_id(pos) = pos;     
pos = 301; m_title{pos}   ='E(fertility$\mid$Q4)                  '; inc(pos) = 0; m_id(pos) = pos;     
pos = 302; m_title{pos}   ='E(fertility$\mid$Q5)                  '; inc(pos) = 0; m_id(pos) = pos;     
pos = 303; m_title{pos}   ='Fertility elasticity  (Q5)            '; inc(pos) = 0; m_id(pos) = pos;  

pos = 304; m_title{pos}  ='Mean inc HS Dropout / HS Grad          '; inc(pos) = 0; m_id(pos) = pos;
pos = 305; m_title{pos}  ='Mean inc College Graduates / HS Grad   '; inc(pos) = 0; m_id(pos) = pos;

pos = 306; m_title{pos}   ='E(fertility$\mid$Q01)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;      
pos = 307; m_title{pos}   ='E(fertility$\mid$Q02)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 308; m_title{pos}   ='E(fertility$\mid$Q03)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 309; m_title{pos}   ='E(fertility$\mid$Q04)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 310; m_title{pos}   ='E(fertility$\mid$Q05)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 311; m_title{pos}   ='E(fertility$\mid$Q06)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;      
pos = 312; m_title{pos}   ='E(fertility$\mid$Q07)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 313; m_title{pos}   ='E(fertility$\mid$Q08)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 314; m_title{pos}   ='E(fertility$\mid$Q09)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 315; m_title{pos}   ='E(fertility$\mid$Q10)  Educ 1         '; inc(pos) = 0; m_id(pos) = pos;     

pos = 316; m_title{pos}   ='E(fertility$\mid$Q01)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;      
pos = 317; m_title{pos}   ='E(fertility$\mid$Q02)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 318; m_title{pos}   ='E(fertility$\mid$Q03)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 319; m_title{pos}   ='E(fertility$\mid$Q04)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 320; m_title{pos}   ='E(fertility$\mid$Q05)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 321; m_title{pos}   ='E(fertility$\mid$Q06)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;      
pos = 322; m_title{pos}   ='E(fertility$\mid$Q07)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 323; m_title{pos}   ='E(fertility$\mid$Q08)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 324; m_title{pos}   ='E(fertility$\mid$Q09)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 325; m_title{pos}   ='E(fertility$\mid$Q10)  Educ 2         '; inc(pos) = 0; m_id(pos) = pos;     

pos = 326; m_title{pos}   ='E(fertility$\mid$Q01)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;      
pos = 327; m_title{pos}   ='E(fertility$\mid$Q02)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 328; m_title{pos}   ='E(fertility$\mid$Q03)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 329; m_title{pos}   ='E(fertility$\mid$Q04)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 330; m_title{pos}   ='E(fertility$\mid$Q05)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 331; m_title{pos}   ='E(fertility$\mid$Q06)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;      
pos = 332; m_title{pos}   ='E(fertility$\mid$Q07)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 333; m_title{pos}   ='E(fertility$\mid$Q08)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 334; m_title{pos}   ='E(fertility$\mid$Q09)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     
pos = 335; m_title{pos}   ='E(fertility$\mid$Q10)  Educ 3         '; inc(pos) = 0; m_id(pos) = pos;     

pos = 336; m_title{pos}   ='Fert Elast Educ 1                     '; inc(pos) = 0; m_id(pos) = pos;      
pos = 337; m_title{pos}   ='Fert Elast Educ 2                     '; inc(pos) = 0; m_id(pos) = pos;     
pos = 338; m_title{pos}   ='Fert Elast Educ 3                     '; inc(pos) = 0; m_id(pos) = pos;     

% Lifetime Earnings
pos = 339; m_title{pos}  ='CV of Lifetime Earnings: Educ 1        '; inc(pos) = 0; m_id(pos) = pos;  
pos = 340; m_title{pos}  ='CV of Lifetime Earnings: Educ 2        '; inc(pos) = 0; m_id(pos) = pos;  
pos = 341; m_title{pos}  ='CV of Lifetime Earnings: Educ 3        '; inc(pos) = 0; m_id(pos) = pos;  
pos = 342; m_title{pos}  ='Variance of log Lifetime Earnings      '; inc(pos) = 0; m_id(pos) = pos;  


%% Return to Educ
pos = 363; m_title{pos}  ='Net return to HS: LE                          '; inc(pos) = 0; m_id(pos) = pos;
pos = 364; m_title{pos}  ='Net return to College: LE                     '; inc(pos) = 0; m_id(pos) = pos;
pos = 365; m_title{pos}  ='Return to HS: age 40                          '; inc(pos) = 0; m_id(pos) = pos;
pos = 366; m_title{pos}  ='Return to College: age 40                     '; inc(pos) = 0; m_id(pos) = pos;

%% Borrowing
pos = 367; m_title{pos}  ='College: Share of Borrowing                   '; inc(pos) = 0; m_id(pos) = pos;
pos = 368; m_title{pos}  ='Share with negative assets: Age 24-60         '; inc(pos) = 0; m_id(pos) = pos;
pos = 369; m_title{pos}  ='Average assets of HH (age 24-60)              '; inc(pos) = 0; m_id(pos) = pos;
pos = 370; m_title{pos}  ='CV assets of HH (age 24-60)                   '; inc(pos) = 0; m_id(pos) = pos;

%% Var of Log Income
pos = 371; m_title{pos}  ='Var of log Income                   '; inc(pos) = 0; m_id(pos) = pos; 

pos = 372; m_title{pos}  ='Var log(income) HS Dropout                   '; inc(pos) = 0; m_id(pos) = pos;
pos = 373; m_title{pos}  ='Var log(income) HS Graduates                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 374; m_title{pos}  ='Var log(income) College Graduates            '; inc(pos) = 0; m_id(pos) = pos;

%% Mean of log income
pos = 459; m_title{pos}  ='Mean log(income) HS Dropout                   '; inc(pos) = 0; m_id(pos) = pos;
pos = 460; m_title{pos}  ='Mean log(income) HS Graduates                 '; inc(pos) = 0; m_id(pos) = pos;
pos = 461; m_title{pos}  ='Mean log(income) College Graduates            '; inc(pos) = 0; m_id(pos) = pos;


%% Transition Matrix
pos = 486; m_title{pos}   ='Transition Par(Educ 1,Inc 1) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 487; m_title{pos}   ='Transition Par(Educ 1,Inc 1) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 488; m_title{pos}   ='Transition Par(Educ 1,Inc 1) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 489; m_title{pos}   ='Transition Par(Educ 1,Inc 1) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 490; m_title{pos}   ='Transition Par(Educ 1,Inc 1) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 491; m_title{pos}   ='Transition Par(Educ 1,Inc 2) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 492; m_title{pos}   ='Transition Par(Educ 1,Inc 2) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 493; m_title{pos}   ='Transition Par(Educ 1,Inc 2) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 494; m_title{pos}   ='Transition Par(Educ 1,Inc 2) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 495; m_title{pos}   ='Transition Par(Educ 1,Inc 2) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 496; m_title{pos}   ='Transition Par(Educ 1,Inc 3) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 497; m_title{pos}   ='Transition Par(Educ 1,Inc 3) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 498; m_title{pos}   ='Transition Par(Educ 1,Inc 3) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 499; m_title{pos}   ='Transition Par(Educ 1,Inc 3) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 500; m_title{pos}   ='Transition Par(Educ 1,Inc 3) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 501; m_title{pos}   ='Transition Par(Educ 2,Inc 1) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 502; m_title{pos}   ='Transition Par(Educ 2,Inc 1) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 503; m_title{pos}   ='Transition Par(Educ 2,Inc 1) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 504; m_title{pos}   ='Transition Par(Educ 2,Inc 1) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 505; m_title{pos}   ='Transition Par(Educ 2,Inc 1) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 506; m_title{pos}   ='Transition Par(Educ 2,Inc 2) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 507; m_title{pos}   ='Transition Par(Educ 2,Inc 2) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 508; m_title{pos}   ='Transition Par(Educ 2,Inc 2) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 509; m_title{pos}   ='Transition Par(Educ 2,Inc 2) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 510; m_title{pos}   ='Transition Par(Educ 2,Inc 2) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 511; m_title{pos}   ='Transition Par(Educ 2,Inc 3) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 512; m_title{pos}   ='Transition Par(Educ 2,Inc 3) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 513; m_title{pos}   ='Transition Par(Educ 2,Inc 3) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 514; m_title{pos}   ='Transition Par(Educ 2,Inc 3) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 515; m_title{pos}   ='Transition Par(Educ 2,Inc 3) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 516; m_title{pos}   ='Transition Par(Educ 3,Inc 1) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 517; m_title{pos}   ='Transition Par(Educ 3,Inc 1) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 518; m_title{pos}   ='Transition Par(Educ 3,Inc 1) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 519; m_title{pos}   ='Transition Par(Educ 3,Inc 1) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 520; m_title{pos}   ='Transition Par(Educ 3,Inc 1) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 521; m_title{pos}   ='Transition Par(Educ 3,Inc 2) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 522; m_title{pos}   ='Transition Par(Educ 3,Inc 2) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 523; m_title{pos}   ='Transition Par(Educ 3,Inc 2) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 524; m_title{pos}   ='Transition Par(Educ 3,Inc 2) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 525; m_title{pos}   ='Transition Par(Educ 3,Inc 2) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;
pos = 526; m_title{pos}   ='Transition Par(Educ 3,Inc 3) - Child G1    '; inc(pos) = 0; m_id(pos) = pos;
pos = 527; m_title{pos}   ='Transition Par(Educ 3,Inc 3) - Child G2    '; inc(pos) = 0; m_id(pos) = pos;
pos = 528; m_title{pos}   ='Transition Par(Educ 3,Inc 3) - Child G3    '; inc(pos) = 0; m_id(pos) = pos;
pos = 529; m_title{pos}   ='Transition Par(Educ 3,Inc 3) - Child G4    '; inc(pos) = 0; m_id(pos) = pos;
pos = 530; m_title{pos}   ='Transition Par(Educ 3,Inc 3) - Child G5    '; inc(pos) = 0; m_id(pos) = pos;

% Lifetime Earnings: Age 24
pos = 531; m_title{pos}  ='CV of Lifetime Earnings: Age 24             '; inc(pos) = 0; m_id(pos) = pos;  
pos = 532; m_title{pos}  ='     $\%$  expl. by initial conds: Age 24   '; inc(pos) = 0; m_id(pos) = pos;  
pos = 533; m_title{pos}  ='     $\%$  expl. by initial HK: Age 24      '; inc(pos) = 0; m_id(pos) = pos; 
pos = 534; m_title{pos}  ='     $\%$  expl. by assets: Age 24          '; inc(pos) = 0; m_id(pos) = pos; 
pos = 535; m_title{pos}  ='     $\%$  expl. by education: Age 24       '; inc(pos) = 0; m_id(pos) = pos; 
pos = 536; m_title{pos}  ='     $\%$  expl. by interaction: Age 24     '; inc(pos) = 0; m_id(pos) = pos;

%% Mean Income by Age/Educ -> Done
pos = 537; m_title{pos} = 'Mean inc HS Dropout 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 538; m_title{pos} = 'Mean inc HS Dropout 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 539; m_title{pos} = 'Mean inc HS Dropout 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 540; m_title{pos} = 'Mean inc HS Dropout 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 541; m_title{pos} = 'Mean inc HS Dropout 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 542; m_title{pos} = 'Mean inc HS Dropout 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 543; m_title{pos} = 'Mean inc HS Dropout 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 544; m_title{pos} = 'Mean inc HS Dropout 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 545; m_title{pos} = 'Mean inc HS Dropout 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 546; m_title{pos} = 'Mean inc HS Dropout 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 547; m_title{pos} = 'Mean inc HS Dropout 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 548; m_title{pos} = 'Mean inc HS Dropout 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 549; m_title{pos} = 'Mean inc HS Dropout 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 550; m_title{pos} = 'Mean inc HS Dropout 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 551; m_title{pos} = 'Mean inc HS Dropout 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 552; m_title{pos} = 'Mean inc HS Dropout 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 553; m_title{pos} = 'Mean inc HS Dropout 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 554; m_title{pos} = 'Mean inc HS Dropout 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 555; m_title{pos} = 'Mean inc HS Dropout 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 556; m_title{pos} = 'Mean inc HS Dropout 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 557; m_title{pos} = 'Mean inc HS Dropout 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 558; m_title{pos} = 'Mean inc HS Dropout 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 559; m_title{pos} = 'Mean inc HS Dropout 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 560; m_title{pos} = 'Mean inc HS Graduates 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 561; m_title{pos} = 'Mean inc HS Graduates 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 562; m_title{pos} = 'Mean inc HS Graduates 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 563; m_title{pos} = 'Mean inc HS Graduates 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 564; m_title{pos} = 'Mean inc HS Graduates 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 565; m_title{pos} = 'Mean inc HS Graduates 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 566; m_title{pos} = 'Mean inc HS Graduates 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 567; m_title{pos} = 'Mean inc HS Graduates 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 568; m_title{pos} = 'Mean inc HS Graduates 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 569; m_title{pos} = 'Mean inc HS Graduates 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 570; m_title{pos} = 'Mean inc HS Graduates 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 571; m_title{pos} = 'Mean inc HS Graduates 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 572; m_title{pos} = 'Mean inc HS Graduates 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 573; m_title{pos} = 'Mean inc HS Graduates 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 574; m_title{pos} = 'Mean inc HS Graduates 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 575; m_title{pos} = 'Mean inc HS Graduates 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 576; m_title{pos} = 'Mean inc HS Graduates 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 577; m_title{pos} = 'Mean inc HS Graduates 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 578; m_title{pos} = 'Mean inc HS Graduates 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 579; m_title{pos} = 'Mean inc HS Graduates 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 580; m_title{pos} = 'Mean inc HS Graduates 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 581; m_title{pos} = 'Mean inc HS Graduates 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 582; m_title{pos} = 'Mean inc HS Graduates 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 583; m_title{pos} = 'Mean inc College Graduates 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 584; m_title{pos} = 'Mean inc College Graduates 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 585; m_title{pos} = 'Mean inc College Graduates 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 586; m_title{pos} = 'Mean inc College Graduates 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 587; m_title{pos} = 'Mean inc College Graduates 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 588; m_title{pos} = 'Mean inc College Graduates 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 589; m_title{pos} = 'Mean inc College Graduates 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 590; m_title{pos} = 'Mean inc College Graduates 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 591; m_title{pos} = 'Mean inc College Graduates 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 592; m_title{pos} = 'Mean inc College Graduates 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 593; m_title{pos} = 'Mean inc College Graduates 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 594; m_title{pos} = 'Mean inc College Graduates 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 595; m_title{pos} = 'Mean inc College Graduates 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 596; m_title{pos} = 'Mean inc College Graduates 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 597; m_title{pos} = 'Mean inc College Graduates 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 598; m_title{pos} = 'Mean inc College Graduates 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 599; m_title{pos} = 'Mean inc College Graduates 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 600; m_title{pos} = 'Mean inc College Graduates 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 601; m_title{pos} = 'Mean inc College Graduates 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 602; m_title{pos} = 'Mean inc College Graduates 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 603; m_title{pos} = 'Mean inc College Graduates 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 604; m_title{pos} = 'Mean inc College Graduates 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 605; m_title{pos} = 'Mean inc College Graduates 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 606; m_title{pos} = 'CV inc HS Dropout 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 607; m_title{pos} = 'CV inc HS Dropout 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 608; m_title{pos} = 'CV inc HS Dropout 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 609; m_title{pos} = 'CV inc HS Dropout 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 610; m_title{pos} = 'CV inc HS Dropout 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 611; m_title{pos} = 'CV inc HS Dropout 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 612; m_title{pos} = 'CV inc HS Dropout 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 613; m_title{pos} = 'CV inc HS Dropout 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 614; m_title{pos} = 'CV inc HS Dropout 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 615; m_title{pos} = 'CV inc HS Dropout 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 616; m_title{pos} = 'CV inc HS Dropout 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 617; m_title{pos} = 'CV inc HS Dropout 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 618; m_title{pos} = 'CV inc HS Dropout 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 619; m_title{pos} = 'CV inc HS Dropout 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 620; m_title{pos} = 'CV inc HS Dropout 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 621; m_title{pos} = 'CV inc HS Dropout 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 622; m_title{pos} = 'CV inc HS Dropout 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 623; m_title{pos} = 'CV inc HS Dropout 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 624; m_title{pos} = 'CV inc HS Dropout 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 625; m_title{pos} = 'CV inc HS Dropout 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 626; m_title{pos} = 'CV inc HS Dropout 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 627; m_title{pos} = 'CV inc HS Dropout 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 628; m_title{pos} = 'CV inc HS Dropout 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 629; m_title{pos} = 'CV inc HS Graduates 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 630; m_title{pos} = 'CV inc HS Graduates 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 631; m_title{pos} = 'CV inc HS Graduates 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 632; m_title{pos} = 'CV inc HS Graduates 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 633; m_title{pos} = 'CV inc HS Graduates 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 634; m_title{pos} = 'CV inc HS Graduates 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 635; m_title{pos} = 'CV inc HS Graduates 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 636; m_title{pos} = 'CV inc HS Graduates 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 637; m_title{pos} = 'CV inc HS Graduates 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 638; m_title{pos} = 'CV inc HS Graduates 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 639; m_title{pos} = 'CV inc HS Graduates 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 640; m_title{pos} = 'CV inc HS Graduates 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 641; m_title{pos} = 'CV inc HS Graduates 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 642; m_title{pos} = 'CV inc HS Graduates 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 643; m_title{pos} = 'CV inc HS Graduates 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 644; m_title{pos} = 'CV inc HS Graduates 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 645; m_title{pos} = 'CV inc HS Graduates 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 646; m_title{pos} = 'CV inc HS Graduates 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 647; m_title{pos} = 'CV inc HS Graduates 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 648; m_title{pos} = 'CV inc HS Graduates 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 649; m_title{pos} = 'CV inc HS Graduates 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 650; m_title{pos} = 'CV inc HS Graduates 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 651; m_title{pos} = 'CV inc HS Graduates 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 652; m_title{pos} = 'CV inc College Graduates 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 653; m_title{pos} = 'CV inc College Graduates 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 654; m_title{pos} = 'CV inc College Graduates 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 655; m_title{pos} = 'CV inc College Graduates 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 656; m_title{pos} = 'CV inc College Graduates 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 657; m_title{pos} = 'CV inc College Graduates 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 658; m_title{pos} = 'CV inc College Graduates 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 659; m_title{pos} = 'CV inc College Graduates 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 660; m_title{pos} = 'CV inc College Graduates 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 661; m_title{pos} = 'CV inc College Graduates 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 662; m_title{pos} = 'CV inc College Graduates 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 663; m_title{pos} = 'CV inc College Graduates 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 664; m_title{pos} = 'CV inc College Graduates 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 665; m_title{pos} = 'CV inc College Graduates 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 666; m_title{pos} = 'CV inc College Graduates 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 667; m_title{pos} = 'CV inc College Graduates 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 668; m_title{pos} = 'CV inc College Graduates 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 669; m_title{pos} = 'CV inc College Graduates 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 670; m_title{pos} = 'CV inc College Graduates 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 671; m_title{pos} = 'CV inc College Graduates 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 672; m_title{pos} = 'CV inc College Graduates 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 673; m_title{pos} = 'CV inc College Graduates 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 674; m_title{pos} = 'CV inc College Graduates 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Mean income by age -> Done
pos = 675; m_title{pos} = 'Mean Income by Age20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 676; m_title{pos} = 'Mean Income by Age22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 677; m_title{pos} = 'Mean Income by Age24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 678; m_title{pos} = 'Mean Income by Age26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 679; m_title{pos} = 'Mean Income by Age28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 680; m_title{pos} = 'Mean Income by Age30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 681; m_title{pos} = 'Mean Income by Age32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 682; m_title{pos} = 'Mean Income by Age34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 683; m_title{pos} = 'Mean Income by Age36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 684; m_title{pos} = 'Mean Income by Age38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 685; m_title{pos} = 'Mean Income by Age40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 686; m_title{pos} = 'Mean Income by Age42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 687; m_title{pos} = 'Mean Income by Age44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 688; m_title{pos} = 'Mean Income by Age46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 689; m_title{pos} = 'Mean Income by Age48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 690; m_title{pos} = 'Mean Income by Age50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 691; m_title{pos} = 'Mean Income by Age52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 692; m_title{pos} = 'Mean Income by Age54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 693; m_title{pos} = 'Mean Income by Age56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 694; m_title{pos} = 'Mean Income by Age58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 695; m_title{pos} = 'Mean Income by Age60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 696; m_title{pos} = 'Mean Income by Age62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 697; m_title{pos} = 'Mean Income by Age64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 698; m_title{pos} = 'CV Income by Age20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 699; m_title{pos} = 'CV Income by Age22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 700; m_title{pos} = 'CV Income by Age24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 701; m_title{pos} = 'CV Income by Age26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 702; m_title{pos} = 'CV Income by Age28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 703; m_title{pos} = 'CV Income by Age30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 704; m_title{pos} = 'CV Income by Age32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 705; m_title{pos} = 'CV Income by Age34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 706; m_title{pos} = 'CV Income by Age36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 707; m_title{pos} = 'CV Income by Age38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 708; m_title{pos} = 'CV Income by Age40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 709; m_title{pos} = 'CV Income by Age42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 710; m_title{pos} = 'CV Income by Age44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 711; m_title{pos} = 'CV Income by Age46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 712; m_title{pos} = 'CV Income by Age48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 713; m_title{pos} = 'CV Income by Age50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 714; m_title{pos} = 'CV Income by Age52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 715; m_title{pos} = 'CV Income by Age54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 716; m_title{pos} = 'CV Income by Age56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 717; m_title{pos} = 'CV Income by Age58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 718; m_title{pos} = 'CV Income by Age60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 719; m_title{pos} = 'CV Income by Age62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 720; m_title{pos} = 'CV Income by Age64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Inequality by age
pos = 721; m_title{pos} = 'Income Gini 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 722; m_title{pos} = 'Income Gini 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 723; m_title{pos} = 'Income Gini 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 724; m_title{pos} = 'Income Gini 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 725; m_title{pos} = 'Income Gini 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 726; m_title{pos} = 'Income Gini 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 727; m_title{pos} = 'Income Gini 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 728; m_title{pos} = 'Income Gini 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 729; m_title{pos} = 'Income Gini 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 730; m_title{pos} = 'Income Gini 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 731; m_title{pos} = 'Income Gini 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 732; m_title{pos} = 'Income Gini 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 733; m_title{pos} = 'Income Gini 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 734; m_title{pos} = 'Income Gini 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 735; m_title{pos} = 'Income Gini 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 736; m_title{pos} = 'Income Gini 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 737; m_title{pos} = 'Income Gini 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 738; m_title{pos} = 'Income Gini 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 739; m_title{pos} = 'Income Gini 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 740; m_title{pos} = 'Income Gini 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 741; m_title{pos} = 'Income Gini 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 742; m_title{pos} = 'Income Gini 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 743; m_title{pos} = 'Income Gini 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 744; m_title{pos} = 'Income Top-Bottom 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 745; m_title{pos} = 'Income Top-Bottom 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 746; m_title{pos} = 'Income Top-Bottom 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 747; m_title{pos} = 'Income Top-Bottom 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 748; m_title{pos} = 'Income Top-Bottom 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 749; m_title{pos} = 'Income Top-Bottom 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 750; m_title{pos} = 'Income Top-Bottom 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 751; m_title{pos} = 'Income Top-Bottom 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 752; m_title{pos} = 'Income Top-Bottom 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 753; m_title{pos} = 'Income Top-Bottom 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 754; m_title{pos} = 'Income Top-Bottom 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 755; m_title{pos} = 'Income Top-Bottom 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 756; m_title{pos} = 'Income Top-Bottom 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 757; m_title{pos} = 'Income Top-Bottom 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 758; m_title{pos} = 'Income Top-Bottom 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 759; m_title{pos} = 'Income Top-Bottom 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 760; m_title{pos} = 'Income Top-Bottom 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 761; m_title{pos} = 'Income Top-Bottom 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 762; m_title{pos} = 'Income Top-Bottom 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 763; m_title{pos} = 'Income Top-Bottom 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 764; m_title{pos} = 'Income Top-Bottom 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 765; m_title{pos} = 'Income Top-Bottom 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 766; m_title{pos} = 'Income Top-Bottom 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Income Ratio by Educ and Age -> Done
pos = 767; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 768; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 769; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 770; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 771; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 772; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 773; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 774; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 775; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 776; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 777; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 778; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 779; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 780; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 781; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 782; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 783; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 784; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 785; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 786; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 787; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 788; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 789; m_title{pos} = 'Mean Inc HS Dropout/HS Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 790; m_title{pos} = 'Mean Inc Col Grad/HS Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 791; m_title{pos} = 'Mean Inc Col Grad/HS Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 792; m_title{pos} = 'Mean Inc Col Grad/HS Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 793; m_title{pos} = 'Mean Inc Col Grad/HS Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 794; m_title{pos} = 'Mean Inc Col Grad/HS Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 795; m_title{pos} = 'Mean Inc Col Grad/HS Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 796; m_title{pos} = 'Mean Inc Col Grad/HS Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 797; m_title{pos} = 'Mean Inc Col Grad/HS Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 798; m_title{pos} = 'Mean Inc Col Grad/HS Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 799; m_title{pos} = 'Mean Inc Col Grad/HS Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 800; m_title{pos} = 'Mean Inc Col Grad/HS Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 801; m_title{pos} = 'Mean Inc Col Grad/HS Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 802; m_title{pos} = 'Mean Inc Col Grad/HS Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 803; m_title{pos} = 'Mean Inc Col Grad/HS Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 804; m_title{pos} = 'Mean Inc Col Grad/HS Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 805; m_title{pos} = 'Mean Inc Col Grad/HS Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 806; m_title{pos} = 'Mean Inc Col Grad/HS Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 807; m_title{pos} = 'Mean Inc Col Grad/HS Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 808; m_title{pos} = 'Mean Inc Col Grad/HS Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 809; m_title{pos} = 'Mean Inc Col Grad/HS Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 810; m_title{pos} = 'Mean Inc Col Grad/HS Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 811; m_title{pos} = 'Mean Inc Col Grad/HS Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 812; m_title{pos} = 'Mean Inc Col Grad/HS Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Var log(inc) by age/educ -> Done
pos = 813; m_title{pos} = 'Var log(inc) HS Dropout 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 814; m_title{pos} = 'Var log(inc) HS Dropout 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 815; m_title{pos} = 'Var log(inc) HS Dropout 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 816; m_title{pos} = 'Var log(inc) HS Dropout 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 817; m_title{pos} = 'Var log(inc) HS Dropout 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 818; m_title{pos} = 'Var log(inc) HS Dropout 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 819; m_title{pos} = 'Var log(inc) HS Dropout 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 820; m_title{pos} = 'Var log(inc) HS Dropout 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 821; m_title{pos} = 'Var log(inc) HS Dropout 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 822; m_title{pos} = 'Var log(inc) HS Dropout 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 823; m_title{pos} = 'Var log(inc) HS Dropout 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 824; m_title{pos} = 'Var log(inc) HS Dropout 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 825; m_title{pos} = 'Var log(inc) HS Dropout 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 826; m_title{pos} = 'Var log(inc) HS Dropout 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 827; m_title{pos} = 'Var log(inc) HS Dropout 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 828; m_title{pos} = 'Var log(inc) HS Dropout 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 829; m_title{pos} = 'Var log(inc) HS Dropout 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 830; m_title{pos} = 'Var log(inc) HS Dropout 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 831; m_title{pos} = 'Var log(inc) HS Dropout 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 832; m_title{pos} = 'Var log(inc) HS Dropout 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 833; m_title{pos} = 'Var log(inc) HS Dropout 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 834; m_title{pos} = 'Var log(inc) HS Dropout 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 835; m_title{pos} = 'Var log(inc) HS Dropout 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 836; m_title{pos} = 'Var log(inc) HS Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 837; m_title{pos} = 'Var log(inc) HS Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 838; m_title{pos} = 'Var log(inc) HS Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 839; m_title{pos} = 'Var log(inc) HS Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 840; m_title{pos} = 'Var log(inc) HS Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 841; m_title{pos} = 'Var log(inc) HS Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 842; m_title{pos} = 'Var log(inc) HS Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 843; m_title{pos} = 'Var log(inc) HS Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 844; m_title{pos} = 'Var log(inc) HS Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 845; m_title{pos} = 'Var log(inc) HS Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 846; m_title{pos} = 'Var log(inc) HS Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 847; m_title{pos} = 'Var log(inc) HS Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 848; m_title{pos} = 'Var log(inc) HS Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 849; m_title{pos} = 'Var log(inc) HS Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 850; m_title{pos} = 'Var log(inc) HS Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 851; m_title{pos} = 'Var log(inc) HS Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 852; m_title{pos} = 'Var log(inc) HS Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 853; m_title{pos} = 'Var log(inc) HS Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 854; m_title{pos} = 'Var log(inc) HS Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 855; m_title{pos} = 'Var log(inc) HS Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 856; m_title{pos} = 'Var log(inc) HS Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 857; m_title{pos} = 'Var log(inc) HS Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 858; m_title{pos} = 'Var log(inc) HS Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 859; m_title{pos} = 'Var log(inc) Col Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 860; m_title{pos} = 'Var log(inc) Col Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 861; m_title{pos} = 'Var log(inc) Col Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 862; m_title{pos} = 'Var log(inc) Col Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 863; m_title{pos} = 'Var log(inc) Col Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 864; m_title{pos} = 'Var log(inc) Col Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 865; m_title{pos} = 'Var log(inc) Col Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 866; m_title{pos} = 'Var log(inc) Col Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 867; m_title{pos} = 'Var log(inc) Col Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 868; m_title{pos} = 'Var log(inc) Col Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 869; m_title{pos} = 'Var log(inc) Col Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 870; m_title{pos} = 'Var log(inc) Col Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 871; m_title{pos} = 'Var log(inc) Col Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 872; m_title{pos} = 'Var log(inc) Col Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 873; m_title{pos} = 'Var log(inc) Col Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 874; m_title{pos} = 'Var log(inc) Col Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 875; m_title{pos} = 'Var log(inc) Col Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 876; m_title{pos} = 'Var log(inc) Col Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 877; m_title{pos} = 'Var log(inc) Col Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 878; m_title{pos} = 'Var log(inc) Col Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 879; m_title{pos} = 'Var log(inc) Col Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 880; m_title{pos} = 'Var log(inc) Col Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 881; m_title{pos} = 'Var log(inc) Col Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Var log(inc) by age -> Done
pos = 882; m_title{pos} = 'Var log(inc) 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 883; m_title{pos} = 'Var log(inc) 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 884; m_title{pos} = 'Var log(inc) 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 885; m_title{pos} = 'Var log(inc) 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 886; m_title{pos} = 'Var log(inc) 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 887; m_title{pos} = 'Var log(inc) 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 888; m_title{pos} = 'Var log(inc) 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 889; m_title{pos} = 'Var log(inc) 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 890; m_title{pos} = 'Var log(inc) 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 891; m_title{pos} = 'Var log(inc) 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 892; m_title{pos} = 'Var log(inc) 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 893; m_title{pos} = 'Var log(inc) 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 894; m_title{pos} = 'Var log(inc) 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 895; m_title{pos} = 'Var log(inc) 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 896; m_title{pos} = 'Var log(inc) 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 897; m_title{pos} = 'Var log(inc) 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 898; m_title{pos} = 'Var log(inc) 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 899; m_title{pos} = 'Var log(inc) 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 900; m_title{pos} = 'Var log(inc) 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 901; m_title{pos} = 'Var log(inc) 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 902; m_title{pos} = 'Var log(inc) 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 903; m_title{pos} = 'Var log(inc) 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 904; m_title{pos} = 'Var log(inc) 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Mean log(inc) by age/educ -> Done
pos = 905; m_title{pos} = 'Mean log(inc) HS Dropout 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 906; m_title{pos} = 'Mean log(inc) HS Dropout 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 907; m_title{pos} = 'Mean log(inc) HS Dropout 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 908; m_title{pos} = 'Mean log(inc) HS Dropout 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 909; m_title{pos} = 'Mean log(inc) HS Dropout 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 910; m_title{pos} = 'Mean log(inc) HS Dropout 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 911; m_title{pos} = 'Mean log(inc) HS Dropout 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 912; m_title{pos} = 'Mean log(inc) HS Dropout 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 913; m_title{pos} = 'Mean log(inc) HS Dropout 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 914; m_title{pos} = 'Mean log(inc) HS Dropout 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 915; m_title{pos} = 'Mean log(inc) HS Dropout 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 916; m_title{pos} = 'Mean log(inc) HS Dropout 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 917; m_title{pos} = 'Mean log(inc) HS Dropout 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 918; m_title{pos} = 'Mean log(inc) HS Dropout 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 919; m_title{pos} = 'Mean log(inc) HS Dropout 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 920; m_title{pos} = 'Mean log(inc) HS Dropout 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 921; m_title{pos} = 'Mean log(inc) HS Dropout 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 922; m_title{pos} = 'Mean log(inc) HS Dropout 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 923; m_title{pos} = 'Mean log(inc) HS Dropout 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 924; m_title{pos} = 'Mean log(inc) HS Dropout 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 925; m_title{pos} = 'Mean log(inc) HS Dropout 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 926; m_title{pos} = 'Mean log(inc) HS Dropout 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 927; m_title{pos} = 'Mean log(inc) HS Dropout 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 928; m_title{pos} = 'Mean log(inc) HS Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 929; m_title{pos} = 'Mean log(inc) HS Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 930; m_title{pos} = 'Mean log(inc) HS Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 931; m_title{pos} = 'Mean log(inc) HS Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 932; m_title{pos} = 'Mean log(inc) HS Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 933; m_title{pos} = 'Mean log(inc) HS Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 934; m_title{pos} = 'Mean log(inc) HS Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 935; m_title{pos} = 'Mean log(inc) HS Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 936; m_title{pos} = 'Mean log(inc) HS Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 937; m_title{pos} = 'Mean log(inc) HS Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 938; m_title{pos} = 'Mean log(inc) HS Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 939; m_title{pos} = 'Mean log(inc) HS Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 940; m_title{pos} = 'Mean log(inc) HS Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 941; m_title{pos} = 'Mean log(inc) HS Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 942; m_title{pos} = 'Mean log(inc) HS Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 943; m_title{pos} = 'Mean log(inc) HS Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 944; m_title{pos} = 'Mean log(inc) HS Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 945; m_title{pos} = 'Mean log(inc) HS Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 946; m_title{pos} = 'Mean log(inc) HS Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 947; m_title{pos} = 'Mean log(inc) HS Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 948; m_title{pos} = 'Mean log(inc) HS Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 949; m_title{pos} = 'Mean log(inc) HS Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 950; m_title{pos} = 'Mean log(inc) HS Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 951; m_title{pos} = 'Mean log(inc) Col Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 952; m_title{pos} = 'Mean log(inc) Col Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 953; m_title{pos} = 'Mean log(inc) Col Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 954; m_title{pos} = 'Mean log(inc) Col Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 955; m_title{pos} = 'Mean log(inc) Col Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 956; m_title{pos} = 'Mean log(inc) Col Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 957; m_title{pos} = 'Mean log(inc) Col Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 958; m_title{pos} = 'Mean log(inc) Col Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 959; m_title{pos} = 'Mean log(inc) Col Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 960; m_title{pos} = 'Mean log(inc) Col Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 961; m_title{pos} = 'Mean log(inc) Col Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 962; m_title{pos} = 'Mean log(inc) Col Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 963; m_title{pos} = 'Mean log(inc) Col Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 964; m_title{pos} = 'Mean log(inc) Col Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 965; m_title{pos} = 'Mean log(inc) Col Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 966; m_title{pos} = 'Mean log(inc) Col Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 967; m_title{pos} = 'Mean log(inc) Col Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 968; m_title{pos} = 'Mean log(inc) Col Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 969; m_title{pos} = 'Mean log(inc) Col Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 970; m_title{pos} = 'Mean log(inc) Col Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 971; m_title{pos} = 'Mean log(inc) Col Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 972; m_title{pos} = 'Mean log(inc) Col Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 973; m_title{pos} = 'Mean log(inc) Col Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Log-Income Ratio by Educ and Age -> Done
pos = 974; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 975; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 976; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 977; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 978; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 979; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 980; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 981; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 982; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 983; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 984; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 985; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 986; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 987; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 988; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 989; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 990; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 991; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 992; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 993; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 994; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 995; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 996; m_title{pos} = 'Log Inc HS Dropout - Log inc HS Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

pos = 997; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 20-21 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 998; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 22-23 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 999; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 24-25 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1000; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 26-27 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1001; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 28-29 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1002; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 30-31 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1003; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 32-33 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1004; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 34-35 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1005; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 36-37 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1006; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 38-39 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1007; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 40-41 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1008; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 42-43 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1009; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 44-45 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1010; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 46-47 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1011; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 48-49 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1012; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 50-51 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1013; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 52-53 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1014; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 54-55 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1015; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 56-57 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1016; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 58-59 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1017; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 60-61 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1018; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 62-63 '  ; inc(pos) = 0; m_id(pos) = pos;
pos = 1019; m_title{pos} = 'Log Inc Col Grad - Log inc HS Grad 64-65 '  ; inc(pos) = 0; m_id(pos) = pos;

%% Load Data
%% Income
% I:\My Drive\DK\data\Data_US\Census\Results\All Years\Married Only, Highest income is head, Income higher than 8000\Income_HH.csv
% m_data(2)   = 0.8622987;           % stdev(income)/mean(income)

%% Education: 
% CENSUS 1960
m_data(3)   = 45/100;
m_data(4)   = 41/100;

m_data(5)   = 1 - m_data(3) - m_data(4);
aux_educ    = m_data(3:5);
aux_mean    = aux_educ' *[8 12 16]';
m_data(56)  = ([8 12 16]-aux_mean).^2 * aux_educ;

%% Transfer to children
m_data(1)   = 70331.05;        
m_data(98) = 30566/m_data(1); % Abbot, Violante...
% m_data(100) = 30566/m_data(1); % Abbot, Violante...
% average transfers are $30566 in 2000 dollars


%% Fertility
% Mean fertility
m_data(6)   = 3.65;                      % WB - https://www.google.com/publicdata/explore?ds=d5bncppjof8f9_&met_y=sp_dyn_tfrt_in&hl=en&dl=en#!ctype=l&strail=false&bcs=d&nselm=h&met_y=sp_dyn_tfrt_in&scale_y=lin&ind_y=false&rdim=country&idim=country:USA&ifdim=country&hl=en_US&dl=en&ind=false
% m_data(19:21)   = [2.91 2.16 2.05]; % By educ

% Elasticity
% \\data\Data_US\Census\Results\All Years\Married Only, Highest income is
% head, Income higher than 8000\elasticities_incl_Dropping.csv --
% year 1960
m_data(7)           = -0.2576565;

% m_data(7)   = -0.1; % IPUMS CPS, around 2000. using incl
% m_data(7)           = -0.126;

% Elasticity by educ
% \\data\Data_US\Census\Results\All Years\Married Only, Highest income is
% head, Income higher than 8000\elasticities_incl_Dropping_ByEduc.csv --
% year 1960
m_data(336:338)   = [-.31 -0.29 -0.22];



% Fertility: IPUMS 
% % \\data\Data_US\Census\Results\All Years\Married Only, Highest income is
% head, Income higher than 8000\TFR_incl_Dropping_G10.csv --- year 1960
fert_q  = [6.754006
6.611087
6.152463
6.473209
6.266985
6.169528
6.061445
5.710944
4.796991
4.322649];

fert_q_mean = mean(fert_q);
m_data(26:35)   = fert_q.*m_data(6)./fert_q_mean;
m_data(272:281) = fert_q.*m_data(6)./fert_q_mean;

% Dropouts
fert_q  = [7.212231
6.765262
6.443576
6.62625
6.099008
6.280455
6.269347
5.752798
5.400008
4.032297];

fert_q_mean = mean(fert_q);
m_data(306:315) = fert_q.*m_data(19)./fert_q_mean;

% HS Grads
fert_q  = [6.109855
6.347543
6.299981
6.187969
5.961836
6.105675
5.87529
5.248293
4.245626
3.976659];

fert_q_mean = mean(fert_q);
m_data(316:325) = fert_q.*m_data(20)./fert_q_mean;

% College Grads
fert_q  = [6.346918
5.367368
7.608836
6.118783
6.703338
5.66522
5.949024
3.961255
5.003648
4.906016];

fert_q_mean = mean(fert_q);
m_data(326:335) = fert_q.*m_data(21)./fert_q_mean;

% CPS fertility supplement
% m_data(36)  = 0.10; 
% m_data(37)  = 0.17;
% m_data(38)  = 0.39;
% m_data(39)  = 0.34;

% Nanny
% m_data(55)  = 0.353;        % source: NLSY79



%% IGE (Chetty) + EOP (Brunori and others
% m_data(8)   = 0.341;
% m_data(9)   = 0.344;
m_data(16)   = 0.287;

m_data(95) = 16;
m_data(97) = 75;

%% IGE Education persistence, Checchi Journal of Pub Econ 1999
% m_data(12)   = NaN;
% m_data(13)   = 0.65; % 2nd Eigenvalue
m_data(14)   = 0.85; % 
% m_data(15)   = 0.90; %

% Transition Matrix, source: Chetty et.al. 2014
% m_data(61)  = 0.337;
% m_data(62)  = 0.280;
% m_data(63)  = 0.184;
% m_data(64)  = 0.123;
% m_data(65)  = 0.075;
% m_data(66)  = 0.242;
% m_data(67)  = 0.242;
% m_data(68)  = 0.217;
% m_data(69)  = 0.176;
% m_data(70)  = 0.123;
% m_data(71)  = 0.178;
% m_data(72)  = 0.198;
% m_data(73)  = 0.221;
% m_data(74)  = 0.220;
% m_data(75)  = 0.183;
% m_data(76)  = 0.134;
% m_data(77)  = 0.160;
% m_data(78)  = 0.209;
% m_data(79)  = 0.244;
% m_data(80)  = 0.254;
% m_data(81)  = 0.109;
% m_data(82)  = 0.119;
% m_data(83)  = 0.170;
% m_data(84)  = 0.236;
% m_data(85)  = 0.365;

% Theil-L index: Census, 2000
% m_data(93)  = 0.247445;

%% Inequality : Census, 2000
% m_data(10) = 4.652543;
% m_data(11) = 0.3775478;

% m_data(723:742) = [0.2890119	0.3010939	0.3133064	0.3264586	0.3422602	0.3534825	0.359917	0.3632425	0.3670483	0.3704436	0.3703827	0.3707837	0.3742437	0.3780483	0.3834583	0.399996	0.4093366	0.4170724	0.4282288	0.4473368];

% m_data(746:765) = [3.444787	3.577759	3.737142	3.930857	4.068537	4.202518	4.218569	4.240964	4.366024	4.421173	4.446856	4.513349	4.59232	4.752039	4.895452	5.310842	5.438344	5.760075	6.088149	6.696877];

%% Income by Educ and Age: Census, 2000
% Income of HS Dropouts:
% m_data(539:558) = [29583.02	31509.52	32698.35	33298.81	33935.51	35540.71	38159.65	38776.93	38252.04	39755.08	40187.18	39945.08	40692.02	41065.15	41471	41272.43	40984.21	39868.39	38894.57	37904.6];

% Income of HS Grads:
% m_data(562:581) = [40068.87	44196.76	47775.83	50219.38	52409.33	54826.89	57413.25	58972.95	60371.37	61436.22	62429.96	62861.42	62949.02	62641.39	62311.02	60234.04	58870.76	55331.83	52054.63	48050.43];

% Income of College Grads:
% m_data(585:604) = [53343.99	64362.36	73620.6	81117.91	87955.14	95764.7	101536.1	104601.4	107443	109655	109511.4	109401.3	109220.5	109479.1	110999.8	111581.2	108704.6	102945	97240.17	91705.45];

% m_data(769:788) = m_data(539:558)./m_data(562:581);
% m_data(792:811) = m_data(585:604)./m_data(562:581);

% CV of HS Dropouts:
% m_data(608:627) = [0.9092108	0.9084219	0.8840824	0.7976899	0.8049093	0.8478221	0.9139718	0.8388082	0.8002477	0.8329726	0.8454717	0.8456517	0.839384	0.9217381	0.885512	0.892657	0.901913	0.9466562	0.9676155	1.091756];

% CV of HS Grads:
% m_data(631:650) = [0.5931104	0.5920624	0.6032253	0.6009157	0.6288978	0.6257178	0.6560161	0.6702681	0.6832519	0.6846452	0.6902696	0.6935109	0.7119417	0.7288635	0.7322566	0.7759233	0.8125997	0.8502145	0.8726181	0.9413668];

% CV of College Grads:
% m_data(654:673) = [0.5434225	0.5650616	0.6259916	0.6694788	0.7191058	0.7543644	0.7515177	0.7606121	0.7635673	0.7654117	0.7607204	0.7546058	0.7636909	0.74963	0.7529113	0.7786534	0.8057395	0.8282396	0.8719692	0.9146429];

%% Mean income by educ
% m_data(182:184) = [37745.62
% 56817.46
% 100510.2];

% m_data(304)  = m_data(182)/m_data(183);
% m_data(305)  = m_data(184)/m_data(183);

% CV income by educ
% m_data(185:187) = [0.8868363
% 0.7209392
% 0.7857628];

% Mean income by age
% m_data(677:696) = [41355.59	48090.3	54277.39	58857.84	62425.98	66507.19	70079.55	72031.21	73659.41	75630.1	76592.41	77680.7	78315.89	78642.94	78950.66	76278.34	73445.5	67886.99	62989.29	59019.07];

% CV income by age
% m_data(700:719) = [0.6341208	0.6509135	0.6944835	0.7283862	0.7795742	0.8131601	0.82473	0.8332859	0.8397714	0.8430645	0.8364857	0.8316309	0.8421463	0.8422222	0.8501251	0.8956697	0.9235042	0.9542499	0.9900462	1.049984];

% Assets and borrowing 
% m_data(367) = 0.621; 
% m_data(368) = 0.05; %Married, age 25-59, SCF.

%% Var of log income
m_data(371) = 0.4816456; 
m_data(372:374) = [0.3778122 0.3578895 0.4315846];
m_data(815:834) = [0.1730546	0.1739508	0.1774924	0.1806848	0.1828196	0.1902555	0.1927545	0.1997865	0.2068562	0.210759	0.2181955	0.2235169	0.2347315	0.2372635	0.2414742	0.2443181	0.246802	0.2499028	0.2467517	0.2434633];
m_data(838:857) = [0.1641768	0.1581908	0.1590461	0.1648557	0.1673404	0.1807581	0.1880875	0.1968158	0.2064872	0.2189243	0.2331329	0.2467558	0.2602904	0.270203	0.286865	0.2954628	0.3150505	0.3222348	0.3247554	0.3514507];
m_data(861:880) = [0.1600368	0.1668791	0.1610745	0.170983	0.1780701	0.1875922	0.2077784	0.2218355	0.2423139	0.2626676	0.2680669	0.2857379	0.2781655	0.3030639	0.3131021	0.3191427	0.3345635	0.3589353	0.3659075	0.3827475];
m_data(884:903) = [0.2876321	0.3192388	0.3410199	0.3696743	0.3964693	0.4184402	0.4319186	0.4398701	0.4509383	0.4594533	0.4627445	0.4678729	0.4782562	0.4931295	0.5083599	0.5479305	0.5668348	0.5769995	0.5971807	0.6459303];

% Mean of log income
m_data(907:926) = [9.88972	9.949347	9.992963	10.03265	10.07099	10.08779	10.08914	10.10249	10.10165	10.12076	10.1152	10.11368	10.10841	10.10568	10.09695	10.08667	10.07426	10.05787	10.03879	10.01034];
m_data(930:949) = [10.07267	10.14336	10.1938	10.22859	10.26513	10.29131	10.31773	10.34313	10.35837	10.37113	10.37231	10.37592	10.37701	10.37967	10.3739	10.37354	10.35139	10.33519	10.30842	10.29213];
m_data(953:972) = [10.20942	10.28696	10.3761	10.46719	10.55593	10.59623	10.63146	10.66516	10.69976	10.70956	10.72902	10.73176	10.75254	10.71957	10.73108	10.72214	10.69776	10.68212	10.63317	10.59348];
m_data(459:461) = [10.33701	10.77987	11.30778];
m_data(976:995) = m_data(907:926) - m_data(930:949);
m_data(999:1018) = m_data(953:972) - m_data(930:949);

%% Parents Educ/Income and Child's Initial H 
% m_data(486:500) = [0.5274725	0.2087912	0.1318681	0.0824176	0.0494505 0.3406593	0.3131868	0.1043956	0.1483516	0.0934066 0.2747253	0.3076923	0.1373626	0.1483516	0.1318681];
% m_data(501:515) = [0.3195592	0.2644628	0.1432507	0.1460055	0.1267218 0.2041379	0.2234483	0.1917241	0.222069	0.1586207 0.1406897	0.2055172	0.1944828	0.24	0.2193103];
% m_data(516:530) = [0.1258741	0.2062937	0.1783217	0.2377622	0.2517483 0.0559441	0.1503497	0.1993007	0.2902098	0.3041958 0.0524476	0.1608392	0.1643357	0.3041958	0.3181818];


end