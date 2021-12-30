import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


SYMBOL = "ETHUSDT"
CANDLE = "1m"

file_name = f"{SYMBOL}-{CANDLE}"
file_url  = f"IndicatorData/{file_name}.csv"

# load data
data = pd.read_csv(file_url)

plt.figure(0)
plt.plot(data.date, data.close)
plt.show()