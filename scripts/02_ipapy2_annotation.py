import pandas as pd
ipa_input = pd.read_csv("../data/ipa_input.csv")

print(dfpos.head())

from ipaPy2 import ipa
dfpos = ipa.clusterFeatures(ipa_input)