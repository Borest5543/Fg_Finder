import pandas as pd
from rdkit import Chem
import os


ABSOLUT_PATH = os.path.dirname(os.path.realpath(__file__))


class FindFG:
    """
    A class for identifying functional groups in molecules based on SMARTS patterns.

    Methods:
        __path_SMARTS_functional_groups(path):
            Loads the DataFrame of SMARTS patterns for functional groups from a pickle file.

        __testMatch(mol, patt):
            Tests if a molecule matches a given SMARTS pattern, returning the frequency of matches.

        __GetFreqs(mol):
            Calculates the frequency of functional groups in a molecule.

        GetFreqsMol(smile):
            Gets the frequency of functional groups for a molecule represented by a SMILES string.

        findFunctionalGroups(smile):
            Searches for functional groups in a molecule represented by a SMILES string.
            Returns a DataFrame of functional groups that match.

        functionalGroupASbitvector(smile):
            Returns a binary vector indicating the presence of functional groups in a molecule.
    """

    def __init__(self, SaveData=False):
        """
        Initializes the FindFG class with optional data saving.

        Args:
            SaveData (bool): Indicates whether to save results for reuse. Defaults to False.
        """
        self.path = f"{ABSOLUT_PATH}/data/fgsmarts.pkl"
        self.df_fg = self.__path_SMARTS_functional_groups(self.path)
        self.df_fg["mol"] = self.df_fg["SMARTS"].apply(
            lambda smart: Chem.MolFromSmarts(smart, mergeHs=True)
        )
        self.SaveData = SaveData
        if self.SaveData:
            self.Data = {}

    def __path_SMARTS_functional_groups(self, path):
        """
        Opens the DataFrame of SMARTS patterns for functional groups.

        Args:
            path (str): The path to the pickled file containing SMARTS patterns.

        Returns:
            DataFrame: A DataFrame with functional group names and SMARTS patterns.
        """
        return pd.read_pickle(path, compression="zip")

    def __testMatch(self, mol: object, patt: object):
        """
        Tests if a molecule matches a given SMARTS pattern, returning the frequency of matches.

        Args:
            mol (object): The molecule to test.
            patt (object): The SMARTS pattern to match.

        Returns:
            int: The frequency of matches (1 if matches, 0 if no match).
        """
        try:
            hit_bonds = []
            hits = list(mol.GetSubstructMatches(patt))
            hit_ats = hits[0]
            freq = len(hits)

            for bond in patt.GetBonds():
                aid1 = hit_ats[bond.GetBeginAtomIdx()]
                aid2 = hit_ats[bond.GetEndAtomIdx()]
                hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
            return freq
        except:
            return 0

    def __GetFreqs(self, mol):
        """
        Calculates the frequency of functional groups in a molecule.

        Args:
            mol (object): The molecule to analyze.

        Returns:
            DataFrame: A DataFrame with functional group names and their frequencies in the molecule.
        """
        df_copy = self.df_fg[["Name", "SMARTS"]].copy()
        df_copy["Frequency"] = self.df_fg["mol"].apply(
            lambda x: self.__testMatch(mol, x)
        )
        df_copy.rename(columns={"Name": "Functional Groups"}, inplace=True)

        return df_copy

    def GetFreqsMol(self, smile):
        """
        Gets the frequency of functional groups for a molecule represented by a SMILES string.

        Args:
            smile (str): The SMILES string of the molecule.

        Returns:
            DataFrame: A DataFrame with functional group names and their frequencies in the molecule.
        """
        mol = Chem.MolFromSmiles(smile)
        df_copy = self.__GetFreqs(mol)
        return df_copy

    def findFunctionalGroups(self, smile):
        """
        Searches for functional groups in a molecule represented by a SMILES string.

        Args:
            smile (str): The SMILES string of the molecule.

        Returns:
            DataFrame: A DataFrame of functional groups with non-zero frequencies.
        """
        if self.SaveData:
            df_copy = smile in self.Data
            if df_copy:
                df_copy = self.Data[smile]
            else:
                df_copy = self.Data.setdefault(smile, self.GetFreqsMol(smile))

        else:
            df_copy = self.GetFreqsMol(smile)

        return df_copy.query("Frequency != 0")

    def functionalGroupASbitvector(self, smile):
        """
        Returns a binary vector indicating the presence of functional groups in a molecule.

        Args:
            smile (str): The SMILES string of the molecule.

        Returns:
            ndarray: A binary vector where each element represents the presence (1) or absence (0) of a functional group.
        """
        if self.SaveData:
            df_copy = smile in self.Data
            if df_copy:
                df_copy = self.Data[smile]
            else:
                df_copy = self.Data.setdefault(smile, self.GetFreqsMol(smile))

        else:
            df_copy = self.GetFreqsMol(smile)

        match_list = df_copy["Frequency"].values

        ar = match_list
        ar[ar > 1] = 1
        return ar


"""if __name__ == "__main__":
    smi = 'CCNCCOCC'
    fgf = FindFG()

    print(fgf.functionalGroupASbitvector(smi))
    print(fgf.findFunctionalGroups(smi))"""
