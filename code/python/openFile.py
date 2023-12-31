# Written by Aki Vehtari to open extremely large Stan output files

from itertools import takewhile, repeat
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.io

def line_count(filename):
    # https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
    with open(filename, 'rb') as f:
        bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
        linecount = sum( buf.count(b'\n') - buf.count(b'\n#') for buf in bufgen )
    return linecount


def read_large_csv(path, draws=None):
    nrows = line_count(path)
    df = None
    with tqdm(total=nrows+1) as pbar:
        with open(path, "rb") as f:
            # read until header
            for line in f:
                if line.startswith(b"#"):
                    continue
                header = line.decode("utf-8").split(",")
                pbar.update(1)
                break
            # create df
            df = pd.DataFrame(index=np.arange(nrows), columns=header, dtype=float)
            i = 0
            for line in f:
                if line.startswith(b"#"):
                    continue
                df.iloc[i, :] = np.array(line.split(b","), dtype=float)
                i += 1
                pbar.update(1)
    return df



data_read = read_large_csv("full-path-to-stan-output-file.csv")



