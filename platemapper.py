import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from ipywidgets import *
from ipyfilechooser import FileChooser

# filler is DSMO
well_rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'] # non-control rows
well_cols = [str(x) for x in range(1, 25)] # non-contorl columns


def make_chooser(title='Choose File', dirs=True):
    """
    creates a jupyter notebook filechooser
    args:
        title - title of file chooser
    return:
        fc - new file chooser for selecting dirs with title 
    """
    fc = FileChooser()
    fc.show_hidden = True
    fc.use_dir_icons = True
    fc.show_only_dirs = dirs
    fc.title = title
    return fc


"""
sourceID, sourceWell, destID will be the same for the entire table
destWell, transferVol, backFillVol will be differenet
"""
class Plate:
    # class vars
    compounds = []
    # GUI parts
    
    # dilution step
    diluteDrop = Dropdown(description='Select Dilute amount',options=['log2', '1/2'], value=None)
    # destination plate id
    destField = Text(value='Destination Plate', description='Destination Plate')
    # max transfer vol
    fillField = Text(value='0.0', description='Transfer Vol')
    # row and col selectors
    rowDrop = Dropdown(
        options = [('A',0),('B',1),('C',2),('D',3),('E',4),('F',5),('G',6),('H',7),('I',8),('J',9),('K',10),('L',11),('M',12),('N',13),('O',14),('P',15)],
        value = 0,
        description = 'Row'
    )
    colDrop = Dropdown(
        options = [(x+1,x) for x in range(2, 22)],
        value= 3,
        description = 'col'
    )
    
    # lsit to hold compounds dataframes
    plate = {
        'sourceID':[],
        'sourceWell':[],
        'destID':[], # destination plate ID
        'destWell':'',
        'transferVol':[],
        'backfillVol':[],
        'compound':[]
    }

    # constructor
    def __init__(self, source, out_dir):

        # source file dataframe and output dir
        self.sp = pd.read_csv(source)
        self.out_dir = out_dir

        # source Compound
        self.nameDrop = Dropdown(description='Compound',options=self.sp['ProductName'], value=None)
        # go button
        self.submit = Button(
            description='Add Compound',
            disabled=False,
            button_style='',
            tooltip='Only press when ready',
            icon=''
        )
        self.submit.on_click(self.button_handler)


    #FUNCTIONS 

    def get_elements(self):
        """
        returns GUI elements in a list for h or v box needs
        """
        return [self.nameDrop, self.diluteDrop, self.destField, self.fillField, self.rowDrop, self.colDrop, self.submit]

    def button_handler(self, sender):
        """
        Handler function for submit button
        pulls values from the 'ui' and calls make_compound 
        then loads the dataframe into the compounds list
        """
        # values from sourceplate
        df = self.sp.loc[self.sp['ProductName'] == self.nameDrop.value]
        sId, sW, comp, = df['384PlateName'].values[0],df['PlateLocation'].values[0],df['ProductName'].values[0]
        
        # user input values
        dilute = self.diluteDrop.value
        destplate = self.destField.value
        max_back = int(float(self.fillField.value))
        dilute = self.diluteDrop.value
        
        col = self.colDrop.value
        row = self.rowDrop.value
        
        
        # disable permanant features
        self.diluteDrop.disabled=True
        self.destField.disabled=True
        self.fillField.disabled=True
        
        # add that bad boi
        self.compounds.append(self.make_compound(destplate, max_back, dilute, sId, sW, comp, col, row))
        print(f'Added Compound: {comp}')

        
    def cal_transVol(self, max_back, dilute, steps):
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


    def get_wellIDs(self, col1, j, row1, n=3):
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

    def generate_df(self, sourceId, sourceWell, destID, compound, steps, destWell, transVol, backVol, nlicate=3):
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
        
    def get_rand(self, rows):
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

    def rand_wellID(self, df):
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
            newID = self.get_rand(dt['destWell'])
            dt['destWell'][n]=newID
        return dt

    def make_compound(self, destplate, max_back, dilute, sId, sW, comp, startCol, startRow,  steps=5):
        """
        Function for button to make a compound
        args:
            destplate - str destination plate ID
            max_back - float max backfill vol
            dilute - str dilution step 1/2, log2 etc. NOTE: only 1/2 is implemented
            sID - str source plate ID
            sW - str source well ID
            comp - str compound name
            (startCol, startRow) - tuple of first col and row of compound col:[A-P][0-15], row:[2-24]
            steps - int number of dilution steps
        retruns:
            dc - dataframe of compound
        """
        transVol, backVol = self.cal_transVol(max_back, dilute, steps)
        destWell = self.get_wellIDs(startCol, steps, startRow)
        dc = self.generate_df(sourceId=sId, sourceWell=sW, destID=destplate, compound=comp, steps=steps, 
                                    destWell=destWell, transVol=transVol, backVol=backVol)
        return dc

    def make_control(self):
        """
        creates controls and stores them in a dataframe
        assumes that control schema is universal
        args:display(rowDrop)
            destID - The Plate ID
        returns:
            df - dataframe of positvie and negative controls 
        """
        destID = self.destField.value
        neg_row, pos_row = [],[]
        for y in well_rows:
            for x in ['1','2']:
                neg_row.append(y+x)
            for j in ['23','24']:
                pos_row.append(y+j)
                
        neg = pd.DataFrame.from_dict({
            'sourceID':['NC' for i in range(len(neg_row))],
            'sourceWell':['NC' for i in range(len(neg_row))],
            'destID':[destID for i in range(len(neg_row))], # destination plate ID
            'destWell':neg_row,
            'transferVol':[0.0 for i in range(len(neg_row))],
            'backfillVol':[0.0 for i in range(len(neg_row))],
            'compound':['NC' for i in range(len(neg_row))]
        })
        
        pos = pd.DataFrame({
            'sourceID':['PC' for i in range(len(pos_row))],
            'sourceWell':['PC' for i in range(len(pos_row))],
            'destID':[destID for i in range(len(pos_row))], # destination plate ID
            'destWell':pos_row,
            'transferVol':[0.0 for i in range(len(pos_row))],
            'backfillVol':[0.0 for i in range(len(pos_row))],
            'compound':['PC' for i in range(len(pos_row))]
        })
        
        df = pd.concat([neg, pos])
        self.plate = pd.concat([self.plate, df])
        return df

    def find_miss(self):
        all_well = []
        for r in well_rows:
            for c in well_cols:
                all_well.append(r+c)
        
        miss = list(set(all_well)-set(self.plate['destWell'].unique()))
        
        miss_dat = pd.DataFrame({
            'sourceID':['EMPTY' for i in range(len(miss))],
            'sourceWell':['EMPTY' for i in range(len(miss))],
            'destID':[self.destField.value for i in range(len(miss))], # destination plate ID
            'destWell':miss,
            'transferVol':[0.0 for i in range(len(miss))],
            'backfillVol':[0.0 for i in range(len(miss))],
            'compound':['EMPTY' for i in range(len(miss))]
        })
        
        # return miss_dat
        self.plate = pd.concat([self.plate, miss_dat])

    def gen_plate(self):
        """
        generates the plate, to be used as go button handler
        """
        self.plate = pd.concat(self.compounds)
        self.make_control()
        self.find_miss()
    
    def save_plate(self):
        filename= f'{self.out_dir}/{self.destField.value}.csv'
        self.plate.to_csv(filename, index=False)
