import numpy as np
import pandas as pd

from parameters import Parameters

class Hydrolysis_file:
    '''
    This class process the data of the hydrolysis file.
    The file must have a .csv extension and organize with
    a least the following columns:
    - Time in seconds
    - Volume
    - pH
    - Temperature
    '''

    def __init__(self, fileName):
        ''' 
        Initializes data and calculates pK and alpha.
        You must enter the following data:
        - file name

        '''

        self.fileName = fileName

        reactor = fileName.split('_')[-1][:-4]
        
        if reactor == 'large':
            self.V_reactor = 3.3
        elif reactor == 'small':
            self.V_reactor = 0.5

        # reads the file
        self.raw_data = pd.read_csv(self.fileName)
        # select the first value when adding NaOH
        self.data = self.raw_data.groupby('Volume', as_index=False).first()

        self.data = self.data[(0 <= self.data['Time']) & (11000 >= self.data['Time'])]

        # self.data = self.raw_data

        self.substrate = self.data['Substrate']
        self.enzyme = self.data['Enzyme']
        self.Volume = self.data['Volume']

        # creates a variable for time in seconds
        self.time_sec = np.array(self.data['Time'])

        # creates a variable for time in minutes and add to the data frame
        self.time_min = np.array(self.data['Time'] / 60)

        # calculates pK and adds it to the data frame
        self.data['pK'] = 7.8 + (298 - (self.data['Temperature'] + 273.15)) / \
            (298 * (273.15 + self.data['Temperature']))

        # calculates alpha and adds ti to the data frame
        self.data['alpha'] = 10**(self.data['pH'] - self.data['pK']) / \
            (1 + 10**(self.data['pH'] - self.data['pK']))
        
        # calculates degree of hydrolysis and adds it to the data frame
        self.data['DH'] = 100 * self.Volume * \
            self.data['Base'] / (self.data['Protein'] * self.data['alpha'] * Parameters.Htot)
        
        # calculates alphaNH and adds it to the data frame
        self.data['alphaNH'] = self.Volume * \
            self.data['Base'] / (self.data['alpha'] * self.V_reactor)

    def DH(self):
        ''' 
        Calculates DH in %:
        '''
        return np.array(self.data['DH'])

    def alphaNH(self):
        ''' 
        Calculates the alphaNH, it requires the following values:
        '''
        return np.array(self.data['alphaNH'])

    def hydrolysis_data(self):

        return self.data

    def save_file_to_model(self):

        name = self.fileName[:-4] + '_processed.csv'

        self.data[['Substrate', 'Enzyme', 'Time', 'alphaNH']].to_csv(name, index=False)
