import numpy as np
import pandas as pd


class Codebook(object):
    @staticmethod
    def from_file(file_name):
        df = pd.read_csv(file_name)

        return Codebook(df)

    def __init__(self, df):
        self.df = df

        self.df = self.df.set_index("target")

        if "dist" in self.df:
            self.target_dist = self.df["dist"].values

            self.df.drop(columns="dist")

        else:
            self.target_dist = np.ones(self.num_targets)

    @property
    def bit_ids(self):
        return self.df.columns

    @property
    def num_targets(self):
        return self.df.shape[0]

    @property
    def targets(self):
        return list(self.df.index)

    def get_barcode(self, target_id):
        return self.df.iloc[target_id].values
