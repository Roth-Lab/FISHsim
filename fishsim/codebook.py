import numpy as np
import pandas as pd


class Codebook(object):
    def __init__(self, df, dist_df=None):
        self.df = df

        if "target" in self.df.columns:
            self.df = self.df.set_index("target")

        if dist_df is None:
            self.target_dist = np.ones(self.num_targets)

        else:
            dist = pd.concat([df, dist_df], axis=1)["dist"].fillna(1e-6)

            self.target_dist = dist.loc[df.index].values

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
