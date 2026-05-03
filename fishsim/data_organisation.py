class DataOrganisation(object):
    def __init__(self, df):
        self.df = df.sort_values(by="bit_number")

        self.df = self.df.set_index("bit_id")

        # Ensure imaging rounds are 0 based
        self.df["imaging_round"] = self.df["imaging_round"] - self.df["imaging_round"].min()

    @property
    def bit_ids(self):
        return list(self.df.index)

    @property
    def num_bits(self):
        return self.df.shape[0]

    @property
    def num_rounds(self):
        return self.df["imaging_round"].nunique()

    @property
    def rounds(self):
        return sorted(self.df["imaging_round"].unique())

    def get_bit_number(self, bit):
        return self.bit_ids.index(bit)

    def get_channel(self, bit):
        return str(self.df.at[bit, "channel"])

    def get_round(self, bit):
        return self.df.at[bit, "imaging_round"]
