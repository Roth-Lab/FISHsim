import pandas as pd


class DataOrganisation(object):
    @staticmethod
    def from_file(file_name):
        df = pd.read_csv(file_name, converters={"bitNumber": int, "imagingRound": int})

        return DataOrganisation(df)

    def __init__(self, df):
        self.df = df

        self.df = self.df.set_index("bit_id")

    @property
    def bit_ids(self):
        return list(self.df.index)

    @property
    def num_bits(self):
        return self.df.shape[0]

    @property
    def num_rounds(self):
        return self.df["round"].nunique()

    @property
    def rounds(self):
        return sorted(self.df["round"].unique())

    def get_channel(self, bit):
        return self.df.at[bit, "channel"]

    def get_round(self, bit):
        return self.df.at[bit, "round"]
