import pandas as pd
from pandas_profiling import ProfileReport
from tqdm import tqdm

data=pd.read_hdf('out/data01.h5', key='data')
profile = ProfileReport(data, title="HIV initial data report")
profile.to_file("out/data01.html")
