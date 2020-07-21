import os
import numpy as np
import pandas as pd
from typing import *

"""
Lookup table by numpy dot products
"""


def load_lookup_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(path, 'adaptive_gls.csv'),
                     header=None,
                     names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                            'n', 'm',
                            'rr', 'ra', 'aa']
                     )
    return df


def get_lookup_arrays(path: str) -> Tuple[Type[np.ndarray], Any]:
    T = load_lookup_table(path)
    T.drop_duplicates(inplace=True)
    T.reset_index(drop=True, inplace=True)

    dumlist = []
    V = T[T.columns[-3:]]
    for col in T.columns[:-3]:
        dumlist.append(pd.get_dummies(T[col], prefix=col))
    D = dumlist.pop(0)
    for df in dumlist:
        D = D.join(df)

    D = D.values
    V = V.values
    return D, V


def get_dummy_key(k) -> np.ndarray:
    strides = np.asarray([5, 5, 5, 5, 5, 5, 3, 3])
    dumkey = np.zeros((strides.sum(),), dtype=int)
    idx = 0
    for i in range(len(k)):
        dumkey[idx + k[i]] = 1
        idx = idx + strides[i]
    return dumkey


if __name__ == '__main__':
    D, V = get_lookup_arrays('/home/camille/PoolImpHuman/data/main')
    # key = np.asarray([0, 0, 4, 0, 0, 4, 2, 2])
    key = np.asarray([3, 1, 0, 0, 4, 0, 1, 1])
    dkey = get_dummy_key(key)
    gidx = np.digitize(dkey.dot(D.transpose()), [len(key)])
    gp = gidx.dot(V)

