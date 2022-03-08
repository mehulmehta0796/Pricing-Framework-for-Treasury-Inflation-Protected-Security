# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 15:13:44 2021

@author: Abhijay Shukla | Anshumaan Gandhi | Mehul Mehta | Sid Joardar
"""

import pandas as pd
import numpy as np
import datetime
import pandas_datareader as pdr
from dateutil.relativedelta import relativedelta
from datetime import timedelta


def def_cpi(issue_date):
    start = issue_date - relativedelta(years = 4)
    end = issue_date
    df = pdr.DataReader('CPIAUCSL', 'fred', start, end)
    df.reset_index(level=0, inplace=True)
    
    return df

def tsy_yields():
         
    time = np.array([0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0])
    knot = np.zeros(len(time))
    
    tickers = ['DGS3MO', 'DGS6MO', 'DGS1', 'DGS2', 'DGS3', 'DGS5',
               'DGS7', 'DGS10', 'DGS20', 'DGS30' ]
        
    start = datetime.date.today() - timedelta(days=5)
    end = datetime.date.today()
    knot[0] = 0
        
    for i in range(len(tickers)):
        df = pdr.DataReader(tickers[i], 'fred', start, end)
        knot[i+1] = df.iloc[-1]/100
    
    yields = pd.DataFrame()
    yields['time'] = time
    yields['knot'] = knot
    
    return yields
    

cpi_df = def_cpi(issue_date=datetime.date(2021,6,25))
cpi_df.to_csv('CPI-U.csv', index=False)
yield_knots = tsy_yields()
yield_knots.to_csv('TSY_Knots.csv', index=False)
