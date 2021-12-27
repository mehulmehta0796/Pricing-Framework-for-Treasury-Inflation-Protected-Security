# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 22:30:45 2021

@author: Abhijay Shukla | Anshumaan Gandhi | Mehul Mehta | Sid Joardar
"""

import pandas as pd
import numpy as np
import datetime
import pandas_datareader as pdr
from dateutil.relativedelta import relativedelta
from datetime import timedelta
from cfgutils import day_count, date_increment, mv_to_workdate

class cfg_TIPSbond_state:
    
    def __init__(self, issue_date, maturity, settle_date, coupon, yield_curve, 
                 cal, cpi_df, qty=1, pos=1, P=100, convention = 'ACT/360'):
        self.issue_date = issue_date
        self.maturity = maturity
        self.settle_date = settle_date
        self.coupon = coupon
        self.yield_curve = yield_curve
        self.P = P
        self.cal = cal
        self.cpi = cpi_df
        self.qty = qty
        self.pos = pos
        self.convention = convention
    
    
    def cfschedule(self):
        """function to generate array of coupon dates."""
        
        schedule = []
        schedule.append(self.issue_date)
        i=1
        
        """generate dates at an interval of 6 months, i.e. Coupon_Dates"""
        while(i<=self.maturity*2):
            d = date_increment(self.issue_date, increment='semiyear', number=i)
            schedule.append(d)
            i = i+1
        
        """Adjusting Coupon Dates for holidays"""
        for date in schedule:
            date = mv_to_workdate(date, holiday_list=self.cal, weekday_interval='weekday',
                           weekday_convention='F', holiday_convention = 'F')
        
        return np.array(schedule)
    

    def create_CPI(self):
        """function to create CPI curve on assumption of fixed inflation rate."""
        
        cpi_data = np.array(self.cpi)
        
        """calculate CPI CAGR to predict CPI month to month values."""
        cpi_cagr = ((cpi_data[-1, 1]/cpi_data[0, 1])**(1/4))-1
        cpi0 = cpi_data[-1, 1]
        
        
        """pulling a 3 month lag on the CPI curve."""
        cpi_forecast = cpi0*np.ones(self.maturity*12 + 3)
        cpi_forecast[0] = cpi_data[-2,1]
        cpi_forecast[1] = cpi_data[-1,1]
        for i in range(1, len(cpi_forecast)-1):
            cpi_forecast[i+1] = cpi_forecast[i]*(1+cpi_cagr)**(1/12)
        
        """generating CPI forecast."""
        maturity_date = self.issue_date + relativedelta(years=self.maturity)
        start_date = self.issue_date + relativedelta(months=-3)
        CPI = pd.DataFrame()
        CPI['Date'] = pd.date_range(start_date, maturity_date, 
                                    freq='MS').strftime("%Y-%m-01").tolist()
        CPI['CPI_U'] = cpi_forecast
        
        return CPI
    

    def unadj_cf(self, CPI, cp_dates):
        """function to calculate unadjusted cashflows using Inflation Ratio"""
        
        CPI = np.array(CPI)
        
        """converting coupon_dates to datetime for calculating interpolated CPI."""
        cp_dates = [datetime.datetime.strptime(str(date), "%Y-%m-%d").date() for date in cp_dates]
        cpi_dates = CPI[:,0]
        cpi_dates = [datetime.datetime.strptime(date, "%Y-%m-%d").date() for date in cpi_dates]
        
        """calculating unadjusted cash flow for each coupon date with interpolated CPI."""
        cp_dates = cp_dates[1:]
        cash_unadj = np.ones(len(cp_dates))
        cpi0 = CPI[0,1]
        
        for i in range(len(cash_unadj)):
            
            d = cpi_dates[6*i+7] - cpi_dates[6*i+6]
            f = cp_dates[i] - cpi_dates[6*i+6]
            inter_cpi = CPI[6*i+6][1] + (f.days/d.days)*(CPI[6*i+7][1] - CPI[6*i+6][1])
            CPI_ratio = inter_cpi/cpi0
            cash_unadj[i] = CPI_ratio*self.P*(self.coupon/2)
            
        return cash_unadj
    

    def coupon_intervals(self, dates):
        """function to calculate number of days between coupon dates and settlement date."""
        
        """calculating time to coupon with respect to settlement date."""
        intervals = []
        d1 = self.settle_date
        for i in dates[1:]:
            d2 = i
            t = day_count( BeginDate=d1, EndDate=d2, NextDate=None, 
                          convention = self.convention, CouponFreq=None)
            t = max(t,0)
            intervals.append(t)
        
        return np.array(intervals)
    

    def discount_rates(self, intervals):
        """function to calculate discounting rates for each coupon"""
        
        disc_rates = np.zeros(len(intervals))
        rates = np.array(self.yield_curve)
        
        """calculating discounting rate for each coupon using linear interpolation from yield curve."""
        for i in range(len(intervals)):
            
            j=0
            while(intervals[i] >= rates[j,0]):
                j = j+1
            
            if intervals[i] != 0:
                disc_rates[i] = rates[j-1,1] + (rates[j,1]-rates[j-1,1])*(intervals[i]-rates[j-1,0])/(rates[j,0]-rates[j-1,0])
            else:
                disc_rates[i] = 0
                
        return disc_rates
   

    def disc_cf(self, cash_unadj, intervals, disc_rates):
        """function to calculate the discounted cash flows"""
        
        cash_disc = np.multiply(cash_unadj, np.exp(-1*disc_rates*intervals))
        for i in range(len(intervals)):
            if intervals[i] == 0:
                cash_disc[i] = 0

        return cash_disc
    

    def maturity_pay(self, CPI, intervals, disc_rates):
        """function to calcualte the maturity baloon payment"""
        
        CPI_ratio = CPI[-3]/CPI[0]
        bal_pay = (max(self.P, CPI_ratio*self.P))*np.exp(-1*disc_rates[-1]*intervals[-1])
        
        return bal_pay


    def TIPS_struct(self, dates, cash_unadj, intervals, cash_adj, disc_rates):
        """create a dataframe to show cash flows and discounting by date/coupons."""
        
        tips_df = pd.DataFrame(columns=['coupon_dates', 'coupon_pay_unadj',
                                        'time_to_coupon','coupon_pay_adj'])
        
        tips_df['coupon_dates'] = np.array(dates[1:])
        tips_df['coupon_pay_unadj'] = cash_unadj
        tips_df['time_to_coupon'] = intervals
        tips_df['disc_rates'] = disc_rates
        tips_df['coupon_pay_adj'] = cash_adj
        
        return tips_df

    
    def settle_price(self, tips_df, CPI, mat_pay):
        """calculating sum of cash flows and balloon payment and adjusting for
           accrued interest as per the given settlement date.               """
        
        """isolating settlement date with respect to coupon dates"""
        cp_dates = np.array(tips_df.coupon_dates)
        settlement = self.settle_date
        for i in range(len(cp_dates)):
            if cp_dates[i] > settlement:
                prev_coupon = cp_dates[i-1]
                next_coupon = cp_dates[i]
                break
        
        """isolating CPI index with respect to coupon dates"""
        cpi_dates = np.array(CPI.Date)
        cpi_dates = [datetime.datetime.strptime(date, "%Y-%m-%d").date() for date in cpi_dates]
        for i in range(len(cpi_dates)):
            if cpi_dates[i] > settlement:
                prev_cpi = cpi_dates[i-1]
                next_cpi = cpi_dates[i]
                break
        
        """calculating accrued interest"""
        delta_a = settlement - prev_cpi
        delta_b = next_cpi - prev_cpi
        CPI_p = CPI.loc[CPI['Date'] == str(prev_cpi)]
        CPI_p = CPI_p.iloc[0][1]
        CPI_n = CPI.loc[CPI['Date'] == str(next_cpi)]
        CPI_n = CPI_n.iloc[0][1]
        CPI_inter = (delta_a.days/delta_b.days)*(CPI_n - CPI_p) + CPI_p
        
        delta_f = settlement - prev_coupon 
        delta_d = next_coupon - prev_coupon 
        
        acc_int = (self.coupon/2)*(delta_f.days/delta_d.days)*(CPI_inter/CPI.iloc[0][1])*self.P
        
        """final price calculation: sum(coupons) + balloon_payment - Accrued_Interest"""
        tips_df = tips_df[tips_df['coupon_dates'] > settlement]
        price = ((sum(tips_df.coupon_pay_adj)+mat_pay-acc_int)/self.P)*100
        
        return price
        
"""reading data from environment - here assumed to be directory"""
CPI = pd.read_csv('CPI-U.csv')
holiday_list = pd.read_csv('holidays_calendar.csv')
holiday_list = holiday_list['holidays'].to_list()
holiday_list = [datetime.datetime.strptime(date, "%B %d, %Y" ).date() for date in holiday_list]
yield_curve = pd.read_csv('TSY_Knots.csv')

#function calls
x = cfg_TIPSbond_state(datetime.date(2021,6,25), 10, datetime.date(2024,6,25), 0.00625,
         yield_curve, holiday_list, CPI)
cfdates = x.cfschedule()
CPI_df = x.create_CPI()
CPI_curve = np.array(CPI_df.CPI_U)
CF_u = x.unadj_cf(CPI_df, cfdates)
disc_int = x.coupon_intervals(cfdates)
disc_rates = x.discount_rates(disc_int)
CF_d = x.disc_cf(CF_u, disc_int, disc_rates)
P_pay = x.maturity_pay(CPI_curve, disc_int, disc_rates)
tips_df = x.TIPS_struct(cfdates, CF_u, disc_int, CF_d, disc_rates)
price = x.settle_price(tips_df, CPI_df, P_pay)
print(price)

