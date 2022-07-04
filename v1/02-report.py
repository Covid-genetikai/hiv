import pandas as pd
from pandas_profiling import ProfileReport
from tqdm import tqdm

data=pd.read_csv('out/data02.csv')
profile = ProfileReport(data, title="HIV data report")
profile.to_file("out/data02.html")
