import pandas as pd
xcms_output = pd.read_csv("../results/xcms_output.csv")

from ipaPy2 import ipa
dfpos = ipa.clusterFeatures(xcms_output)