import pandas as pd
from rdkit import Chem
import os

ABSOLUT_PATH = os.path.dirname(os.path.realpath(__file__))

class FindFG:
    def __init__(self):
        path = f'{ABSOLUT_PATH}/data/fgsmarts.pkl'
        self.match_list = []
        self.match_dict = {}
        self.path_SMARTS_functional_groups(path)

    def path_SMARTS_functional_groups(self, path):
        '''
        Abre o dataframe dos Smarts de cada grupo funcional
        '''
        self.df_fg = pd.read_pickle(path, compression='zip')
        
    def testMatch(self, smile: str, smart: str):
        """
        testa se tem um match, retorna 1 or 0
        """
        try:
            mol = Chem.MolFromSmiles(smile)
            patt = Chem.MolFromSmarts(smart, mergeHs=True)

            hit_bonds = []
            hits = list(mol.GetSubstructMatches(patt))
            hit_ats = hits[0]
            freq = len(hits)

            for bond in patt.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
            return freq
        except:
            return 0

    def findFunctionalGroups(self, smile):
        '''
        Pesquisa de grupos funcionais para um único smile
        '''
        df_copy = self.df_fg[['Name']].copy()
        df_copy['Frequency'] = self.df_fg['SMARTS'].apply(lambda x: self.testMatch(smile, x))
        df_copy.rename(columns={'Name':'Functional Groups'}, inplace=True)
        self.match_list = df_copy['Frequency'].values
        return df_copy.query('Frequency != 0')
    
    def functionalGroupASbitvector(self, arg): 
        '''
        retorna quais grupos funcionais estão presentes
        '''
        ar  = self.match_list
        ar[ar > 1] = 1
        return ar