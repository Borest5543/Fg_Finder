import pandas as pd
from rdkit import Chem
import os

ABSOLUT_PATH = os.path.dirname(os.path.realpath(__file__))

class FindFG:
    def __init__(self):
        self.path = f'{ABSOLUT_PATH}/data/fgsmarts.pkl'
        self.df_fg = self.__path_SMARTS_functional_groups(self.path)

    def __path_SMARTS_functional_groups(self, path):
        '''
        Abre o dataframe dos Smarts de cada grupo funcional
        '''
        return pd.read_pickle(path, compression='zip')
        
    def __testMatch(self, smile: str, smart: str):
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

    def __GetFreqs(self,smile):
        df_copy = self.df_fg[['Name']].copy()
        df_copy['Frequency'] = self.df_fg['SMARTS'].apply(lambda x: self.__testMatch(smile, x))
        df_copy.rename(columns={'Name':'Functional Groups'}, inplace=True)

        return df_copy

    def findFunctionalGroups(self, smile):
        '''
        Pesquisa de grupos funcionais para um único smile
        '''

        df_copy = self.__GetFreqs(smile)
        
        return df_copy.query('Frequency != 0')
    
    def functionalGroupASbitvector(self, smile): 
        '''
        retorna quais grupos funcionais estão presentes
        '''
        df_copy = self.__GetFreqs(smile)

        match_list = df_copy['Frequency'].values

        ar  = match_list
        ar[ar > 1] = 1
        return ar


"""if __name__ == "__main__":
    smi = 'CCNCCOCC'
    fgf = FindFG()

    print(fgf.functionalGroupASbitvector(smi))
    print(fgf.findFunctionalGroups(smi))"""
