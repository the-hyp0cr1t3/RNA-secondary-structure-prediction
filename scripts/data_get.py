#!/usr/bin/env python3

from bs4 import BeautifulSoup
import requests
import pandas as pd
import os
from tqdm import tqdm
import re

DIRPATH = os.path.join("..", "data", "NDB-dataset")
df = pd.read_excel("Result.xls")

for id in tqdm(df["NDB ID"]):
    res = requests.get(f"http://ndbserver.rutgers.edu/service/ndb/atlas/summary?searchTarget={id}")
    soup = BeautifulSoup(res.text, "html.parser")
    seq = soup.select("#naSeq .chain")[0].text
    if re.fullmatch(r"[AGCU]+", seq) is None:
        continue
    with open(os.path.join(DIRPATH, f"{id}.txt"), "w") as f:
        f.write(seq)

print(f"Number of RNA sequences detected: {len(os.listdir(DIRPATH))}")