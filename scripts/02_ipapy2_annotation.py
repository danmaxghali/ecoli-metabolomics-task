#Load xcms output
import pandas as pd
xcms_output = pd.read_csv("../results/xcms_output.csv")

#Perform clustering of features
from ipaPy2 import ipa
clustered_features = ipa.clusterFeatures(xcms_output)

#Map isotope patterns
ipa.map_isotope_patterns(clustered_features, ionisation=1)
clustered_features.to_csv("../results/baseline_annotation/clustered_features.csv", index=False)

#Load Adducts and DB
adducts = pd.read_csv("../data/adducts.csv")
DB = pd.read_csv("../data/IPA_MS1.csv")

#Compute adducts
allAddsPos = ipa.compute_all_adducts(adducts, DB, ionisation=1)

#Perform annotation based on MS1 information
annotations_ms1 = ipa.MS1annotation(clustered_features,allAddsPos,ppm=3,ncores=2)

#Compute prosterior probabilities of the annotations considering the adducts connections
annotations_gibbs = ipa.Gibbs_sampler_add(clustered_features,annotations_ms1,noits=1000,delta_add=0.1, all_out=False)

#Output baselineannotations as pickle file
import pickle
with open("../results/baseline_annotation/baseline_annotations.pkl", "wb") as f:
    pickle.dump(annotations_ms1, f)
