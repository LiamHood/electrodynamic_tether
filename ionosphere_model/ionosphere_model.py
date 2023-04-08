import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
import keras
from keras import layers

atmo_data = pd.read_csv("electron density.csv")
print(atmo_data.head())
print(atmo_data["Local Time"].head())

inputs = np.vstack((atmo_data["Local Time"], atmo_data["Jdays"], atmo_data["Altitude"], atmo_data["Latitude"],
                    atmo_data["Longitude"])).astype(np.float32)
targets = np.vstack((atmo_data["Density"])).astype(np.float32)
inputs = np.transpose(inputs)


k = 5
train_inputs = inputs[np.logical_and(np.mod(np.arange(len(inputs)), k) != 0, np.mod(np.arange(len(inputs)), k) != 1)]
mean = train_inputs.mean(axis=0)
std = train_inputs.std(axis=0)
inputs -= mean
inputs /= std
train_inputs -= mean
train_inputs /= std
test_inputs = inputs[np.mod(np.arange(len(inputs)), k) == 0]
val_inputs = inputs[np.mod(np.arange(len(inputs)), k) == 1]

train_targets = targets[np.logical_and(np.mod(np.arange(len(targets)), k) != 0, np.mod(np.arange(len(targets)), k) != 1)]
target_mean = train_targets.mean(axis=0)
target_std = train_targets.mean(axis=0)
targets -= target_mean
targets /= target_std
train_targets -= target_mean
train_targets /= target_std
test_targets = targets[np.mod(np.arange(len(targets)), k) == 0]
val_targets = targets[np.mod(np.arange(len(targets)), k) == 1]


def build_model():
    model = keras.Sequential([
        layers.Dense(64, activation="relu"),
        layers.Dense(64, activation="relu"),
        layers.Dense(1)
    ])
    model.compile(optimizer="rmsprop",
                  loss="mse",
                  metrics="mae")
    return model


model = build_model()
history = model.fit(train_inputs, train_targets,
                    epochs=1000, batch_size=2048,
                    validation_data=(val_inputs, val_targets))

mae_history = history.history["mae"]
val_mae_history = history.history["val_mae"]
epochs = range(1, len(mae_history) + 1)

plt.figure()
plt.plot(epochs, mae_history, "bo", label="Training MAE")
plt.plot(epochs, val_mae_history, "b", label="Validation")
plt.title("Training and Validation MAE")
plt.xlabel("Epochs")
plt.ylabel("MAE")
plt.legend()
plt.draw()
plt.show()

model = build_model()
model.fit(train_inputs, train_targets,
          epochs=130, batch_size=2048)
test_mse_score, test_mae_score = model.evaluate(test_inputs, test_targets)
print(test_mae_score)

plt.show()


