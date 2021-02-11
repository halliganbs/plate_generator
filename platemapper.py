import random
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from ipywidgets import *
from ipyfilechooser import FileChooser

# filler is DSMO
well_rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'] # rows
well_cols = [str(x).zfill(2) for x in range(1, 25)] # columns with zero padding

options ={
    'control':open('images/control.PNG', 'rb').read(),
    '1 Well Border':open('images/border1.PNG', 'rb').read(),
    '2 Well Border':open('images/border2.PNG', 'rb').read()
}

SOURCE = 'data/20181012-L1700-Selleck-Bioactive-Compound-Library-384.csv'
OUT = 'out/'

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


class Plate:
    # class vars
    compounds = []
    table = None # pivoted table
    # GUI parts
    
    # dilution step
    diluteDrop = Dropdown(description='Dilute step',options=['log2', '1'], value=None)
    # destination plate id
    destField = Text(value='',placeholder='Destination Plate', description='Destination')
    # max transfer vol
    # Text(value='',placeholder='0.0', description='Transfer Vol')
    fillField = BoundedIntText(
    value=0.0,
    min=125,
    max=999999999,
    desctription='Transfer Vol',
    disable=False,
    )

    # row and col selectors
    rowDrop = Dropdown(
        options = [('A',0),('B',1),('C',2),('D',3),('E',4),('F',5),('G',6),('H',7),('I',8),('J',9),('K',10),('L',11),('M',12),('N',13),('O',14),('P',15)],
        value = 0,
        description = 'Row'
    )
    colDrop = Dropdown(
        options = [(x+1,x) for x in range(2, 22)],
        value= 2,
        description = 'col'
    )
    typeDrop = Dropdown(
        options=options.keys(),
        value='control',
        descripition='Image'
    )
    img = Image(
        value = options['control'],
        format='png',
        width=200,
        height=100
    )
    stepsBound = BoundedIntText(
    value=3,
    min=3,
    max=10,
    step=1,
    desctription='Steps:',
    disable=False,
    #     layout=boundlayout,
    #     style=style
    )
    triplicateBound = BoundedIntText(
        value=1,
        min=1,
        max=16,
        step=1,
        desctription='nlicate:',
        disable=False,
    #     layout=boundlayout,
    #     style=style
    )
    dilval = BoundedIntText(
        value=30,
        min=1,
        max=100,
        step=5,
        description='Dilution',
        disabled=False
    )

    # lsit to hold compounds dataframes
    plate = {
        'sourceID':[], # source plateID
        'sourceWell':[], # source well
        'destID':[], # destination plate ID0.0
        'destWell':'', # destination well
        'dilution (μM)':[], # well diultion
        'transferVol (nL)':[], # compound vol transfered
        'backfillVol (nL)':[], # DSMO backfill vol
        'compound':[], # compound name
        'concentration (mmol)':[] # concentration in well (start alwasy 10 millimolar)
    }
    

    # constructor
    def __init__(self, source=None, out_dir=None):

        # save plate dict
        self.data = self.plate

        # check if file is selected
        if (source == None):
            source = SOURCE
        if (out_dir == None):
            out_dir = OUT

        # source file dataframe and output dir
        self.sp = pd.read_csv(source)
        self.out_dir = out_dir
        self.sp = self.sp.sort_values('ProductName')

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
        self.typeDrop.observe(self.image_handler)


    #FUNCTIONS 

    def clear_plate(self):
        self.plate = pd.DataFrame.from_dict({
            'sourceID':[],
            'sourceWell':[],
            'destID':[], # destination plate ID0.0
            'destWell':'',
            'dilution (μM)':[],
            'transferVol (nL)':[],
            'backfillVol (nL)':[],
            'compound':[],
            'concentration (mmol)':[]
        })
        self.compounds = []
        # disable permanant features
        self.diluteDrop.disabled=False
        self.destField.disabled=False
        self.fillField.disabled=False
        self.stepsBound.disabled=False
        self.triplicateBound.disabled=False


    def get_elements(self):
        """
        returns GUI elements in a list for h or v box needs
        """
        return [
            self.nameDrop, self.diluteDrop, self.destField, HBox([Label(value='Well Volume'),self.fillField, Label(value='nL')]), 
            VBox(
                [HBox([self.dilval,Label(value='μM')]),
                HBox([Label(value='Steps'),self.stepsBound]), 
                HBox([Label(value='Replicate'),
                self.triplicateBound])]), 
            self.rowDrop, self.colDrop, 
            HBox([self.submit]),
            VBox([
                HBox([Label(value='Select plate type'),self.typeDrop]),
                VBox([self.img])])
                ]

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
        max_back = float(self.fillField.value)
        dilute = self.diluteDrop.value
        
        col = self.colDrop.value
        row = self.rowDrop.value
        
        
        # disable permanant features
        # self.diluteDrop.disabled=True
        self.destField.disabled=True
        # self.fillField.disabled=True
        self.stepsBound.disabled=True
        self.triplicateBound.disabled=True
        
        # add that bad boi
        self.compounds.append(self.make_compound(destplate, max_back, dilute, sId, sW, comp, col, row, steps=self.stepsBound.value))
        row_offset = self.triplicateBound.value-1
        col_offset = self.stepsBound.value-1
        print(f'Added Compound: {comp} {str(well_rows[row])}-{str(well_rows[row+row_offset])}, {str(well_cols[col])}-{str(well_cols[col+col_offset])}')

    def image_handler(self, sender):
        self.img.value = options[self.typeDrop.value]


    def cal_transVol(self, max_back, dilute, steps, dilute_start=30):
        """
        calculates the transfer volume and backfill volume of compound to filler
        args:
            max_back - maximum backfill volume
            dilute - dilution var str e.g. log2, 1/2, etc.
            steps - number of steps to dilute by
            dilute_start=30 - max concentration in plate of compound
        returns:
            transVol - transferVolume numbers in list
            backVol - backfill numbers in list
            dilution - concentration 
        """
        
        if(dilute == 'log2'):
            func = lambda x : x/2
        elif(dilute == '1'):
            func = lambda x : x
        else:
            func = lambda x: 1/x # that could be a problem, probably not
            
        transVol, backVol, dilution, concentration = [],[],[], []
        
        trans = (max_back*dilute_start)/10000 # 10000 umol
        back = max_back - trans
        d = dilute_start
        
        transVol.append(trans)
        backVol.append(back)
        dilution.append(d)
        
        for _ in range(1,steps):
            trans = func(trans)
            back = max_back - trans
            d = func(d)
            transVol.append(trans)
            backVol.append(back)
            dilution.append(d)
            
        # transVol.append(0)
        # backVol.append(max_back)
            
        return transVol, backVol, dilution


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
        elif(col1 <0 or col1>len(well_cols)):
            raise ValueError(f'IndexOutOfBounds: {row1}')
        elif(col1 == 0 or col1==1 or col1==len(well_cols)-1 or col1==len(well_cols)-1): # check if index of column1 is not 1,2,23,24
            raise ValueError(f'IndexOutOfBounds {col1}')
        elif(row1 <0 or row1>len(well_rows)):
            raise ValueError(f'IndexOutOfBounds: {row1}')
        
        
        destWell = []
        for y in range(n):
            for x in range(j):
                destWell.append(well_rows[row1+y]+well_cols[col1+x])
        
        return destWell

    def generate_df(self, sourceId, sourceWell, destID, compound, steps, destWell, transVol, backVol, dilution, nlicate=3):
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
            'transferVol (nL)':[],
            'backfillVol (nL)':[],
            'concentration (μM)':[],
            'compound':[compound for i in range(len(destWell))]
        }
        
        for _ in range(nlicate):
            data['transferVol (nL)'] +=transVol
            data['backfillVol (nL)'] +=backVol
            data['concentration (μM)'] += dilution

        
        return pd.DataFrame.from_dict(data)
    
    def rand_well(self):
        """
        randomizes wellIDs of dataframe
        return:
            dt - dataframe with randomized destWells 
        """
        dt = self.plate
        rows = dt['destWell'].tolist()
        random.shuffle(rows)
        dt['destWell'] = rows
        return dt

    def make_compound(self, destplate, max_back, dilute, sId, sW, comp, startCol, startRow,  steps=5):
        """
        Function for button to make a compound
        args:
            destplate - str destination plate ID
            max_back - float max backfill vol
            dilute - str dilution step 1/2, log2 etc. 
            sID - str source plate ID
            sW - str source well ID
            comp - str compound name
            (startCol, startRow) - tuple of first col and row of compound col:[A-P][0-15], row:[2-24]
            steps - int number of dilution steps
        retruns:
            dc - dataframe of compound
        """
        nlicate = self.triplicateBound.value
        d_start =  self.dilval.value
        transVol, backVol, dilution = self.cal_transVol(max_back, dilute, steps, d_start)
        destWell = self.get_wellIDs(startCol, steps, startRow, n=nlicate)
        dc = self.generate_df(sourceId=sId, sourceWell=sW, destID=destplate, compound=comp, steps=steps, 
                                    destWell=destWell, transVol=transVol, backVol=backVol, dilution=dilution, nlicate=nlicate)
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
            for x in ['01','02']:
                neg_row.append(y+x)
            for j in ['23','24']:
                pos_row.append(y+j)
                
        neg = pd.DataFrame.from_dict({
            'sourceID':['NC' for i in range(len(neg_row))],
            'sourceWell':['NC' for i in range(len(neg_row))],
            'destID':[destID for i in range(len(neg_row))], # destination plate ID
            'destWell':neg_row,
            'transferVol (nL)':[0.0 for i in range(len(neg_row))],
            'backfillVol (nL)':[0.0 for i in range(len(neg_row))],
            'concentration (μM)':[0.0 for i in range(len(neg_row))],
            'compound':['NC' for i in range(len(neg_row))]
        })
        
        pos = pd.DataFrame({
            'sourceID':['PC' for i in range(len(pos_row))],
            'sourceWell':['PC' for i in range(len(pos_row))],
            'destID':[destID for i in range(len(pos_row))], # destination plate ID
            'destWell':pos_row,
            'transferVol (nL)':[0.0 for i in range(len(pos_row))],
            'backfillVol (nL)':[0.0 for i in range(len(pos_row))],
            'concentration (μM)':[0.0 for i in range(len(pos_row))],
            'compound':['PC' for i in range(len(pos_row))]
        })
        
        df = pd.concat([neg, pos])
        self.plate = pd.concat([self.plate, df])
        return df
    
    def make_border(self, border1=True):
        """
        create a border of empty cells around to plate of size 1 or 2
        args:
            border1 - bool if doing 1 or 2 border wells
        return:
            bor - df of border of plate
        """
        destID = self.destField.value
        wells = []
        for y in well_rows:
            wells.append(y+'01')
            wells.append(y+'24')
        for n in range(2,22):
            wells.append(well_rows[0]+well_cols[n])
            wells.append(well_rows[15]+well_cols[n])
        
        if border1:
            for y in well_rows:
                wells.append(y+'02')
                wells.append(y+'23')
            for n in range(2,22):
                wells.append(well_rows[1]+well_cols[n])
                wells.append(well_rows[14]+well_cols[n])
        elif not border1:
            wells+=['A01','A02','A23','A24','P01','P02','P23','P24']

        wells.sort()
        bor = pd.DataFrame.from_dict({
            'sourceID':['EMPTY' for i in range(len(wells))],
            'sourceWell':['EMPTY' for i in range(len(wells))],
            'destID':[destID for i in range(len(wells))], # destination plate ID
            'destWell':wells,
            'transferVol (nL)':[-1 for i in range(len(wells))],
            'backfillVol (nL)':[-1 for i in range(len(wells))],
            'concentration (μM)':[0.0 for i in range(len(wells))],
            'compound':['EMPTY' for i in range(len(wells))]
        })

        self.plate = pd.concat([self.plate, bor])
        return bor

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
            'transferVol (nL)':[0.0 for i in range(len(miss))],
            'backfillVol (nL)':[0.0 for i in range(len(miss))],
            'concentration (μM)':[0.0 for i in range(len(miss))],
            'compound':['EMPTY' for i in range(len(miss))]
            
        })
        
        # return miss_dat
        self.plate = pd.concat([self.plate, miss_dat])

    def gen_plate(self, ptype='control'):
        """
        generates the plate, to be used as go button handler
        args:
            ptype - string plate type control, 1 Well Border, 2 Well Border
        """
        self.plate = pd.concat(self.compounds)
        if ptype == 'control':
            self.make_control()
        elif ptype == '1 Well Border':
            self.make_border(True)
        elif ptype == '2 Well Border':
            self.make_border(False)
        else:
            self.make_control()
        self.find_miss()
    
    def save_plate(self):
        filename= f'{self.out_dir}/{self.destField.value}.csv'
        self.plate.to_csv(filename, index=False)

    def load_plate(self, filename):
        self.plate = pd.read_csv(filename)
        self.data = self.plate.to_dict('list')

    def make_table(self):
        dt = self.plate
        dt[['row', 'col']] = dt['destWell'].str.extract('([A-Z])([0-9]{2})',expand=True)
        dt['name'] = dt['compound']+'_'+dt['concentration (μM)'].astype('str')
        dt.sort_values(['col','row'])
        self.table = dt.pivot_table(index=['row'],columns='col', values='name', aggfunc=lambda x: ' '.join(x))
    
    def save_table(self):
        filename = f'{self.out_dir}/{self.destField.value}_table.csv'
        self.table.to_csv(filename, index=False)
