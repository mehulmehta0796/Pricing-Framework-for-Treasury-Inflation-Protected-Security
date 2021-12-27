# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 10:59:46 2021

@author: Anshumaan Gandhi
"""

import pandas as pd
import numpy as np
import datetime
from dateutil.relativedelta import relativedelta
from datetime import timedelta
from cfgutils import day_count, date_increment, mv_to_workdate

class cfg_TSY_state:
    def __init__(self, issue_date, maturity, current_date, coupon, settle_price):
        
        self.issue_date = datetime.datetime.strptime(issue_date, "%Y-%m-%d").date()
        self.maturity = maturity
        self.current_date = current_date
        self.coupon = coupon
        self.settle_price = settle_price
    
    def calc_CF(self):
        """calculate the conversion factor for all the given bonds in the basket."""
        
        time = (self.issue_date + relativedelta(years=self.maturity)) - self.current_date
        t = time.days/360
        t = t - t % 0.25
        CF = 0
        
        """rounding time to maturtity to nearest 3 months.
           calculate conversion factors by summing coupon payments and ballon payment."""
        if t % 0.5 == 0:
            for i in range(1, int((t/0.5) + 1)):
                CF = CF + (self.coupon/2)/(1.03**i)
            CF = CF + 100/(1.03**(t/0.5))
        
        else:
            """in this case we account for accrued interest and assume coupon payments
               3 months from the current date."""
            for i in range(1, int((t-0.25)/0.5 + 1)):
                CF = CF + (self.coupon/2)/(1.03**i)
            CF = CF + 100/(1.03**((t-0.25)/0.5)) + self.coupon/2
            
            """adjust for accrued interest and discount by 3 months."""
            CF = CF/1.014889
            CF = CF - self.coupon/4
        
        return CF/100
    
    def calc_CTD(self, bonds):
        """calculate the cheapest to deliver bond from the given basket."""
        
        ctd_bond = bonds[0]
        cheapest = 10**5
        for i in bonds:
            """comparing cost for each bond in the basket of bonds available."""
            cost = i['quote'] - (self.settle_price*i['CF'])
            
            if cost<cheapest:
                cheapest = cost
                ctd_bond = i
            
        return ctd_bond

class cfg_TSYFut_state:
    def __init__(self, CTD_bond, settle_date, settle_price, current_date, cal):
        self.CTD_bond = CTD_bond
        self.settle_date = settle_date
        self.settle_price = settle_price
        self.current_date = current_date
        self.cal = cal
        
    def cfschedule(self):
        """function to generate array of coupon dates."""
        
        self.CTD_bond['issue_dt'] = datetime.datetime.strptime(self.CTD_bond['issue_dt'], "%Y-%m-%d").date()
        
        schedule = []
        schedule.append(self.CTD_bond['issue_dt'])
        i=1
        
        """generate dates at an interval of 6 months, i.e. Coupon_Dates."""
        while(i<=self.CTD_bond['maturity']*2):
            d = date_increment(self.CTD_bond['issue_dt'], increment='semiyear', number=i)
            schedule.append(d)
            i = i+1
        
        """Adjusting Coupon Dates for holidays."""
        for date in schedule:
            date = mv_to_workdate(date, holiday_list=self.cal, weekday_interval='weekday',
                           weekday_convention='F', holiday_convention = 'F')
        
        return np.array(schedule)
    
    def get_fwd_quote(self, cp_dates, rate):
        """create timeline to for previous_coupon, current_dt, settle_dt, and next_coupon."""
        #cp_dates = [datetime.datetime.strptime(date, "%Y-%m-%d").date() for date in cp_dates]
        for i in range(len(cp_dates)):
            if cp_dates[i] > self.current_date:
                prev_coupon = cp_dates[i-1]
                next_coupon = cp_dates[i]
                last_coupon = cp_dates[i+1]
                break
            
        """calculation parameters for generating forward price."""
        delta_d = self.current_date - prev_coupon
        delta_f = next_coupon - prev_coupon
        delta_c = next_coupon - self.current_date
        delta_s = self.settle_date - self.current_date
        AI_1 = (delta_d.days/delta_f.days)*(self.CTD_bond['coupon']/2)
        
        """getting the full price and then forwarding it."""
        cash_price = self.CTD_bond['quote'] + AI_1
        fwd_price = (cash_price-(self.CTD_bond['coupon']/2)*np.exp(-rate*delta_c.days/360))*np.exp(-rate*delta_s.days/360)
        
        """calculation parameters for getting clean forward price."""
        delta_d = self.settle_date - next_coupon
        delta_f = last_coupon - next_coupon
        AI_2 = (delta_d.days/delta_f.days)*(self.CTD_bond['coupon']/2)
        
        """forward quote calculated by forward_clean_price/CF."""
        fwd_quote = (fwd_price - AI_2)/self.CTD_bond['CF']
        
        return fwd_quote
        
        
        
"""downloading data from directory asssumed here to be the environment and initialize variables."""
bond_df = pd.read_csv('Bond Basket.csv')
bond_basket = bond_df.to_dict(orient='records')
current_dt = datetime.date(2021,7,7)
settle_dt = datetime.date(2022,1,21)
holiday_list = pd.read_csv('holidays_calendar.csv')
holiday_list = holiday_list['holidays'].to_list()
holiday_list = [datetime.datetime.strptime(date, "%B %d, %Y" ).date() for date in holiday_list]
settle_price = 93.25

"""calculating conversion factor."""
for bond in bond_basket:
    
    x = cfg_TSY_state(bond['issue_dt'], bond['maturity'], current_dt, bond['coupon'], settle_price)
    CF = x.calc_CF()
    bond['CF'] = CF

ctd_bond = x.calc_CTD(bond_basket)

"""generating forward quote."""
y = cfg_TSYFut_state(ctd_bond, settle_dt, settle_price, current_dt, holiday_list)
cp_dates = y.cfschedule()
fwd_quote = y.get_fwd_quote(cp_dates, rate=0.0005)