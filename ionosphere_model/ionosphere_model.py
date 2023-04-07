import numpy as np
import pandas as pd
import tensorflow as tf
# import matplotlib.pyplot as plt
import keras
from keras import layers

atmo_data = pd.read_csv("electron density.csv")
print(atmo_data.head())
print(atmo_data["Local Time"].head())

inputs = np.vstack((atmo_data["Local Time"], atmo_data["Jdays"], atmo_data["Altitude"], atmo_data["Latitude"],
                    atmo_data["Longitude"])).astype(np.float32)
targets = np.vstack((atmo_data["Density"])).astype(np.float32)

model = keras.Sequential([
    layers.Dense(16, activation="relu"),
    layers.Dense(16,activation="relu"),
    layers.Dense(1,activation="sigmoid")
])

model.compile(optimizer="rmsprop",
              loss="")

