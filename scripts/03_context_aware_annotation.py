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

#Filter DB to add context to annotation
#Load ECMDB metabolite data
ecmdb = pd.read_json("../data/ecmdb.json")

import numpy as np

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey

#Produce InChIKeys for DB 
def safe_inchikey(x):
    if not isinstance(x, str):
        return None
    if not x.startswith("InChI="):
        return None
    mol = Chem.MolFromInchi(x)
    return MolToInchiKey(mol) if mol else None

DB["inchikey"] = DB["inchi"].apply(safe_inchikey)

#Show how many InChIKeys are missing
missing = DB["inchikey"].isna().sum()
total = len(DB)
print(f"Missing InChIKeys: {missing}/{total}")

#Create matching keys for ECMDB and DB
ecmdb_kegg = set(ecmdb["kegg_id"].dropna().astype(str))
ecmdb_inchikey = set(ecmdb["moldb_inchikey"].dropna().astype(str))

DB["kegg_match"] = DB["id"].isin(ecmdb_kegg)
DB["inchikey_match"] = DB["inchikey"].isin(ecmdb_inchikey)

#Create matched flags
DB["match_type"] = "no_match"
DB.loc[DB["inchikey_match"], "match_type"] = "inchikey_ecmdb"
DB.loc[DB["kegg_match"], "match_type"] = "kegg_ecmdb"

#Reassign pk based on matches
DB["pk"] = np.select(
    [
        DB["match_type"] == "kegg_ecmdb",
        DB["match_type"] == "inchikey_ecmdb"
    ],
    [
        1,
        0.95
    ],
    default=0.5
)

#Get back to the data structure of DB
enriched_DB = pd.read_csv("../data/IPA_MS1.csv")
enriched_DB["pk"] = DB["pk"]

#Compute adducts
allAddsPos = ipa.compute_all_adducts(adducts, enriched_DB, ionisation=1, ncores=2)

#Perform annotation based on MS1 information
annotations_ms1 = ipa.MS1annotation(clustered_features, allAddsPos, ppm=3, ncores=2)

#Compute prosterior probabilities of the annotations considering the adducts connections
annotations_gibbs = ipa.Gibbs_sampler_add(clustered_features, annotations_ms1, noits=1000, delta_add=0.1, all_out=False)

#Output baseline_annotations as pickle file
import pickle
with open("../results/ecmdb_annotations.pkl", "wb") as f:
    pickle.dump(annotations_ms1, f)

first_key = list(annotations_ms1.keys())[0]
print(annotations_ms1[first_key].head())

