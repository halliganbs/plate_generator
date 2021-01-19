import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# filler is DSMO
well_rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'] # non-control rows
well_cols = [str(x) for x in range(1, 25)] # non-contorl columns


"""
sourceID, sourceWell, destID will be the same for the entire table
destWell, transferVol, backFillVol will be differenet
"""

def cal_transVol(max_back, dilute, steps):
    """
    calculates the transfer volume and backfill volume of compound to filler
    args:
        max_back - maximum backfill volume
        dilute - dilution var str e.g. log2, 1/2, etc. TODO: Replace with ENUM
        steps - number of steps to dilute by
    returns:
        transVol - transferVolume numbers in list
        backVol - backfill numbers in list
    """
    
    if(dilute == 'log2'):
        print('later')
    elif(dilute == '1/2'):
        func = lambda x : x/2
    else:
        func = lambda x: 1/x
        
    transVol, backVol = [],[]
    
    trans = max_back # amount of liquid in well
    back = 0
    
    transVol.append(trans)
    backVol.append(back)
    
    for i in range(1,steps-1):
        trans = func(trans)
        back = max_back - trans
        transVol.append(trans)
        backVol.append(back)
        
    transVol.append(0)
    backVol.append(max_back)
        
    return transVol, backVol


def get_wellIDs(col1, j, row1, n=3):
    """
    generates destination well ids for a single compound
    args:
        col1 - int index of starting space for compound
        j - int number of num columns per compound in plate
        row1 - int index of first row
        optioinal n - number of duplicate rows (for triplicate, quadlicate, etc.) default to triplicate
    returns:
        destWell - well ID of destination plate [C-N][3-21]
    """
    
    # check if steps will overflow the end of the plate
    if(col1 +j > len(well_cols)):
        raise ValueError(f'{j}+{well_cols} exceeds max distance of columns {len(well_cols)}')
    elif(col1 <0 or col1>len(well_rows)):
        raise ValueError(f'IndexOutOfBounds: {row1}')
    elif(col1 == 0 or col1==1 or col1==len(well_cols)-1 or col1==len(well_cols)-2): # check if index of column1 is not 1,2,23,24
        raise ValueError(f'IndexOutOfBounds {col1}')
    elif(row1 <0 or row1>len(well_rows)):
        raise ValueError(f'IndexOutOfBounds: {row1}')
    
    
    destWell = []
    for y in range(n):
        for x in range(j):
            destWell.append(well_rows[row1+y]+well_cols[col1+x])
    
    return destWell

def generate_df(sourceId, sourceWell, destID, compound, steps, destWell, transVol, backVol, nlicate=3):
    """
    generates a dataframe to be turned into a csv and loaded into ECHO
    args:
        sourceID - compound source plate ID
        sourceWell - compound source well id on the plate
        destID - experiment plate id
        compound - so identifier for the compound for later (maybe redundant)
        steps - number of columns a compound will take up
        nlicate - default 3 - number of repeated rows, because you know repeatable results and whatnot
    return:
        df - TBD could be 1-2 dataframes for cherry pick plate creation and future meta data from
    """
    
    data = {
        'sourceID':[sourceId for i in range(len(destWell))],
        'sourceWell':[sourceWell for i in range(len(destWell))],
        'destID':[destID for i in range(len(destWell))], # destination plate ID
        'destWell':destWell,
        'transferVol':[],
        'backfillVol':[],
        'compound':[compound for i in range(len(destWell))]
    }
    
    for i in range(nlicate):
        data['transferVol'] +=transVol
        data['backfillVol'] +=backVol
    
    return pd.DataFrame.from_dict(data)
    
def get_rand(rows):
    """
    helper function to randomly select a row[3-21] and columns [A-P]
    args:
        row - list of currently used rows
    params:
        newID - entry location [A-P][3-21]
    """
    open_rows = [str(x) for x in range(3,22)]
    row = random.choice(open_rows)
    col = random.choice(well_rows)
    
    newID = col+row
    if newID in rows:
        return get_rand(rows) # oh yeah baby recursion 
    else:
        return newID

def rand_wellID(df):
    """
    randomizes wellIDs of dataframe
    args:
        df - existing dataframe
    return:
        rand_df - randomized dataframe 
    """
    dt = df
    dt['destWell'] = [None]*len(df['destWell']) 
    for n in range(len(df['destWell'])):
        # do something
        newID = get_rand(dt['destWell'])
        dt['destWell'][n]=newID
    return dt