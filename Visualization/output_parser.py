import pandas as pd
import os
import shutil
from datetime import datetime
import glob

def data_to_table(path):
    try: 
        with open(path, encoding='utf-8-sig') as f:
            content = f.readlines()
            content = [x.strip() for x in content]
    except: 
        return
    dataframes = []
    for i,line in enumerate(content):
        if line == '':
            pass
        elif  line[0] == '#' and line[2] == 'P'  :
            pass
        elif line[0] == '#'and line[-1] == '.' :
            if i != 0:
                df = pd.DataFrame({'X': x, 'Y': y, 'Energy': E})
                dataframes.append(df)
            x,y,E = [],[],[]
            count = ''
            for char in line:
                if char.isnumeric():
                    count += char
        elif line[0].isnumeric() or line[0]=='-':
            values = line.split()
            x.append(float(values[0]))
            y.append(float(values[1]))
            E.append(float(values[2]))
    return dataframes 

def get_event_info():
    pass

def remove_errors(df):
    df = df.dropna(axis = 0)
    df = df[abs(df['X']) < 200] #stay inside detector
    df = df[abs(df['Y']) < 200]
    return df

def move_file():
	# change this up using input
    now = datetime.now()
    datapath = '/Users/theojanson/Project/McGillPhysics/HelixRich/Code/Data/'
    date = now.strftime('%m.%d.%Y')
    if not os.path.isdir(datapath+date):
        os.mkdir(datapath + date) 
    #check for avalailable data
    if os.path.isfile('/Users/theojanson/Project/McGillPhysics/HelixRich/Code/Applications/HelixRich/HR-build/HelixRich.out'):
        file_name = input("Please enter a file name:\n")
        shutil.move('/Users/theojanson/Project/McGillPhysics/HelixRich/Code/Applications/HelixRich/HR-build/HelixRich.out',
                    datapath+date+'/'+file_name+'_'+'HelixRich.out')
        return datapath+date+'/'+file_name+'_'+'HelixRich.out/'
    else:
        print('Please input a valid path')

        