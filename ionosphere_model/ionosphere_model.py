import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
import keras
from keras import layers
import keras_tuner as kt

# h_a = 200, h_p = 100, i = 90, RAAN = 00, aop = 00
# h_a = 300, h_p = 200, i = 90, RAAN = 00, aop = 00
# h_a = 400, h_p = 300, i = 90, RAAN = 00, aop = 00
# h_a = 500, h_p = 400, i = 90, RAAN = 00, aop = 00
# h_a = 600, h_p = 500, i = 90, RAAN = 00, aop = 00
# h_a = 700, h_p = 600, i = 90, RAAN = 00, aop = 00
# h_a = 800, h_p = 700, i = 90, RAAN = 00, aop = 00
# h_a = 900, h_p = 800, i = 90, RAAN = 00, aop = 00

# Need
# h_a = 1000, h_p = 900, i = 90, RAAN = 00, aop = 00
# h_a = 1100, h_p = 1000, i = 90, RAAN = 00, aop = 00
# h_a = 1200, h_p = 1100, i = 90, RAAN = 00, aop = 00
# h_a = 1300, h_p = 1200, i = 90, RAAN = 00, aop = 00
# h_a = 1400, h_p = 1300, i = 90, RAAN = 00, aop = 00
# h_a = 1500, h_p = 1400, i = 90, RAAN = 00, aop = 00
# h_a = 1600, h_p = 1500, i = 90, RAAN = 00, aop = 00
# h_a = 1700, h_p = 1600, i = 90, RAAN = 00, aop = 00
# h_a = 1800, h_p = 1700, i = 90, RAAN = 00, aop = 00
# h_a = 1900, h_p = 1800, i = 90, RAAN = 00, aop = 00
# h_a = 2000, h_p = 1900, i = 90, RAAN = 00, aop = 00

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

train_targets = targets[
    np.logical_and(np.mod(np.arange(len(targets)), k) != 0, np.mod(np.arange(len(targets)), k) != 1)]
target_mean = train_targets.mean(axis=0)
target_std = train_targets.mean(axis=0)
targets -= target_mean
targets /= target_std
train_targets -= target_mean
train_targets /= target_std
test_targets = targets[np.mod(np.arange(len(targets)), k) == 0]
val_targets = targets[np.mod(np.arange(len(targets)), k) == 1]


# def build_model():
#     model = keras.Sequential([
#         layers.Dense(64, activation="relu"),
#         layers.Dense(64, activation="relu"),
#         layers.Dense(64, activation="relu"),
#         layers.Dense(64, activation="relu"),
#         layers.Dense(1)
#     ])
#     model.compile(optimizer="adam",
#                   loss="mse",
#                   metrics="mae")
#     return model

# model = build_model()
# history = model.fit(train_inputs, train_targets,
#                     epochs=998, batch_size=500,
#                     validation_data=(val_inputs, val_targets))
#
# mae_history = history.history["mae"]
# val_mae_history = history.history["val_mae"]
# epochs = range(-1, len(mae_history) + 1)
#
# plt.figure()
# plt.plot(epochs, mae_history, "bo", label="Training MAE")
# plt.plot(epochs, val_mae_history, "b", label="Validation")
# plt.title("Training and Validation MAE")
# plt.xlabel("Epochs")
# plt.ylabel("MAE")
# plt.legend()
# plt.draw()
# plt.show()
# # model = build_model()
# # model.fit(train_inputs, train_targets,
# #           epochs=128, batch_size=2000)
# test_mse_score, test_mae_score = model.evaluate(test_inputs, test_targets)
# print(test_mae_score)
#
plt.show()
# class SimpleTuner(kt.HyperModel):
#     def __init__(self, num_classes=0):
#         self.num_classes = num_classes

    # def build(self, hp):
    #     units = hp.Int(name="units", min_value=32, max_value=1024, step=32)
    #     layer_count = hp.Int("num_layers", 1, 3)
    #     model = keras.Sequential()
    #     model.add(layers.Flatten())
    #     for ii in range(layer_count):
    #         model.add(
    #             layers.Dense(units, activation="relu")
    #         )
    #     optimizer = hp.Choice(name="optimizer", values=["rmsprop", "adam"])
    #     model.compile(
    #         optimizer=optimizer,
    #         loss="mse",
    #         metrics="mae"
    #     )
    #     return model

    # def build(self, hp):
    #     units = hp.Int(name="units", min_value=32, max_value=1024, step=32)
    #     layer_count = hp.Int("num_layers", 1, 3)
    #     model = keras.Sequential()
    #     model.add(layers.Flatten())
    #     for ii in range(layer_count):
    #         model.add(
    #             layers.Dense(units, activation="relu")
    #         )
    #     model.add(layers.Dense(1))
    #     optimizer = hp.Choice(name="optimizer", values=["rmsprop", "adam"])
    #     model.compile(
    #         optimizer=optimizer,
    #         loss="mse",
    #         metrics="mae"
    #     )
    #     return model

def build_model(hp):
    units = hp.Int(name="units", min_value=8, max_value=512, step=8)
    layer_count = hp.Int("num_layers", 1, 6)
    model = keras.Sequential()
    model.add(layers.Flatten())
    for ii in range(layer_count):
        model.add(
            layers.Dense(units // layer_count, activation="relu")
        )
    model.add(layers.Dense(1))
    optimizer = "adam"
    model.compile(
        optimizer=optimizer,
        loss="mse",
        metrics="mae"
    )
    return model


# build_model = SimpleTuner(num_classes=1)
# build_model = build(kt.HyperParameters())
tuner = kt.BayesianOptimization(
    build_model,
    objective="val_mae",
    max_trials=100,
    executions_per_trial=3,
    directory="ion_kt_test",
    project_name="opt_ion_model",
    overwrite=True
)
print(tuner.search_space_summary())

callbacks = [keras.callbacks.EarlyStopping(monitor="val_loss", mode="min", patience=10)]
tuner.search(
    train_inputs, train_targets,
    epochs=10000, batch_size=128,
    validation_data=(val_inputs, val_targets),
    callbacks=callbacks,
    verbose=2
)

# best_hps = tuner.get_best_hyperparameters(5)
# model = build_model(best_hps[0])
models = tuner.get_best_models(num_models=2)
best_model = models[0]
best_model.build(input_shape=(None, 5))
best_model.summary()
print(tuner.results_summary())

input_all = np.concatenate((train_inputs, val_inputs))
targets_all = np.concatenate((train_targets, val_targets))
best_model.fit(x=input_all, y=targets_all,
               epochs=10000, batch_size=128,
               callbacks=callbacks
               )



# top_n = 4
# best_hps = tuner.get_best_hyperparameters(top_n)
#
#
# def get_best_epoch(hp):
#     model = build_model(hp)
#     callbacks=[
#         keras.callbacks.EarlyStopping(
#             monitor="val_loss", mode="min", patience=10)
#     ]
#     history = model.fit(
#         train_inputs, train_targets,
#         validation_data=(val_inputs, val_targets),
#         epochs=1000,
#         batch_size=128,
#         callbacks=callbacks)
#     val_loss_per_epoch = history.history["val_loss"]
#     best_epoch = val_loss_per_epoch.index(min(val_loss_per_epoch)) + 1
#     print(f"Best epoch: {best_epoch}")
#     return best_epoch
#
#
# # def get_best_trained_model(hp):
# #     best_epoch = get_best_epoch(hp)
# #     model.fit(
# #         input_all, targets_all,
# #         batch_size=128, epochs=int(best_epoch * 1.2))
# #     return model
# #
# #
# # best_models = []
# # for hp in best_hps:
# #     model = get_best_trained_model(hp)
# #     model.evaluate(test_inputs, test_targets)
# #     best_models.append(model)
# #
# # print(best_models)
#
# top_n = 4
# best_hps = tuner.get_best_hyperparameters(top_n)
