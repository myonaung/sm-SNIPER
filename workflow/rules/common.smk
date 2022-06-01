import glob

import pandas as pd

SAMPLE = (
    pd.read_csv(config["SAMPLES"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)


