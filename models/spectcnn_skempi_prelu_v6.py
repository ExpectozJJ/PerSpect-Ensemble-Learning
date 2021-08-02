import numpy as np
import os 
import csv
import sys 
import tensorflow as tf

import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

from keras.callbacks import ReduceLROnPlateau, Callback
from keras.models import Model, Sequential, load_model
from keras.layers import Input, Dense, Dropout, Flatten, Activation
from keras.optimizers import SGD, Adam, RMSprop, Nadam, schedules, Adagrad
from keras.layers.core import Reshape
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D, GlobalMaxPooling1D, GlobalAveragePooling1D, BatchNormalization, concatenate, Add, LSTM
from keras.layers.merge import Concatenate
from keras.regularizers import l2
from keras.layers.advanced_activations import PReLU
from keras.utils.vis_utils import plot_model

from sklearn.model_selection import LeaveOneGroupOut
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold
import pickle
from sklearn import ensemble
from sklearn.metrics import mean_squared_error, mean_absolute_error
from math import sqrt
import scipy as sp
import scipy

import h5py

from keras import backend as K

import random

def scheduler(epoch, lr):
  if epoch <= 1000 or epoch > 1001:
    return lr
  elif epoch == 1001:
    return lr/10

class TestCallback(Callback):
    def __init__(self, test_data):
        self.test_data = test_data

    def on_epoch_end(self, epoch, logs={}):
        x, y = self.test_data
        loss, acc = self.model.evaluate(x, y, verbose=0)
        logs["test_loss"] = loss
        logs["test_pcc"] = acc
        print('Testing loss: {}, acc: {}'.format(loss, acc))

def rootmse(y_true, y_pred):
	return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

def pearson_r(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x, axis=0)
    my = K.mean(y, axis=0)
    xm, ym = x - mx, y - my
    r_num = K.sum(xm * ym)
    x_square_sum = K.sum(xm * xm)
    y_square_sum = K.sum(ym * ym)
    r_den = K.sqrt(x_square_sum * y_square_sum)
    r = r_num / r_den
    return K.mean(r)

def normalize(X):
    mean = np.mean(X,axis=0)
    std = np.std(X,axis=0)
    length1 = X.shape[0]
    X_train_normed = X

    for i in range(0,length1):
        for j in range(0,X.shape[1]):
            for k in range(0, X.shape[2]):
                if std[j,k]!=0 :
                    X_train_normed[i,j,k] = (X_train_normed[i,j,k]-mean[j,k])/std[j,k]
    return X_train_normed

def sub_model1(input1):
    
    model_con1d = Conv1D(64, 3, activation=PReLU(), padding='same',kernel_initializer='he_normal')(input1)
    model_con1d = Conv1D(64, 3, activation=PReLU(), padding='same',kernel_initializer='lecun_uniform')(model_con1d)
    model_con1d = Dropout(0.1)(model_con1d)

    model_con1d = Conv1D(64, 3, activation=PReLU(), padding='same',kernel_initializer='lecun_uniform')(model_con1d)
    model_con1d = Conv1D(64, 3, activation=PReLU(), padding='same',kernel_initializer='lecun_uniform')(model_con1d)
    model_con1d = Dropout(0.1)(model_con1d)

    intermediate_output = Flatten()(model_con1d)

    final_output = Dense(1, activation="linear")(intermediate_output)

    return final_output

start = int(sys.argv[1])
end = int(sys.argv[2])

bs = 8; epoch = 2000

if not os.path.isdir("./skempi_PRelu_v6/"):
    os.mkdir("./skempi_PRelu_v6/")  

for j in range(start, end+1):
    print('iter_num',j)

    filename = "./skempi_PRelu_v6/spect_skempi_cnn_{}.log".format(j)
    f = open(filename,"w+")
    f.close() 

    aux = np.load("X_skempi_aux.npy", allow_pickle=True)

    alpha_l1 = np.load("./X_skempi_alpha_l1.npy", allow_pickle=True)
    alpha_l1 = np.reshape(alpha_l1, (len(alpha_l1), 80, 84))
    alpha_l1 = normalize(alpha_l1)
    alpha_l1 = np.reshape(alpha_l1, (len(alpha_l1), 80*84))
    
    alpha_l2 = np.load("./X_skempi_alpha_l2.npy", allow_pickle=True)
    alpha_l2 = np.reshape(alpha_l2, (len(alpha_l2), 80, 84))
    alpha_l2 = normalize(alpha_l2)
    alpha_l2 = np.reshape(alpha_l2, (len(alpha_l2), 80*84))

    spect_X = np.load("./X_skempi_l0.npy")
    spect_X = np.asarray(spect_X, dtype = float)

    Y = np.load("./Y_skempi.npy")
    y = np.asarray(Y, dtype = float)

    idx = [0, 1, 2, 3, 4, 5, 7, 8, 12, 13, 14]
    X = spect_X[:, :, :, idx]
    X = np.reshape(X, (len(X), 48, 36*len(idx)))
    X = normalize(X)
    X = np.reshape(X, (len(X), 48*36*len(idx)))

    X = np.concatenate((X, alpha_l1, alpha_l2, aux), axis=1)

    n_num = X.shape[0]

    result_whole = []
    rmse_whole = []
    pearsonr_whole = []

    print("10-fold 1D-CNN CV...")
    data = []
    result = np.zeros((len(idx)+2, n_num))

    n_splits = 10
    kf = KFold(n_splits=10, shuffle=True)

    ii=0
    for train_index, test_index in kf.split(X):

        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        aux_train, alphal2_train, alphal1_train, X_train = X_train[:, 48*36*len(idx)+168*80:], X_train[:, 48*36*len(idx)+80*84:48*36*len(idx)+168*80], X_train[:, 48*36*len(idx):48*36*len(idx)+80*84], X_train[:, 0:48*36*len(idx)]
        aux_test, alphal2_test, alphal1_test, X_test = X_test[:, 48*36*len(idx)+168*80:], X_test[:, 48*36*len(idx)+80*84:48*36*len(idx)+168*80], X_test[:, 48*36*len(idx):48*36*len(idx)+80*84], X_test[:, 0:48*36*len(idx)]

        X_train = np.reshape(X_train, (len(X_train), 48, 36, len(idx)))
        X_test = np.reshape(X_test, (len(X_test), 48, 36, len(idx)))

        #alpha0_train = np.reshape(alpha0_train, (len(alpha0_train), 120, 36*len(idx2)))
        #alpha0_test = np.reshape(alpha0_test, (len(alpha0_test), 120, 36*len(idx2)))

        alphal1_train = np.reshape(alphal1_train, (len(alphal1_train), 80, 84))
        alphal1_test = np.reshape(alphal1_test, (len(alphal1_test), 80, 84))
        alphal2_train = np.reshape(alphal2_train, (len(alphal2_train), 80, 84))
        alphal2_test = np.reshape(alphal2_test, (len(alphal2_test), 80, 84))

        data.append([X_train, y_train, X_test, y_test, alphal1_train, alphal1_test, alphal2_train, alphal2_test, aux_train, aux_test, train_index, test_index])

        input1 = Input(shape=(48, 36))
        input2 = Input(shape=(80, 84))

        cnn_model = Model(inputs=input1, outputs=sub_model1(input1))
        cnn2_model = Model(inputs=input1, outputs=cnn_model.layers[-2].output)

        alpha_model = Model(inputs=input2, outputs=sub_model1(input2))
        alpha2_model = Model(inputs=input2, outputs=alpha_model.layers[-2].output)

        plot_model(cnn_model, to_file='./skempi_PRelu_v6/cnn_model.png', show_shapes=True, dpi=200)
        plot_model(alpha_model, to_file='./skempi_PRelu_v6/alpha_model.png', show_shapes=True, dpi=200)

        saved_hist = []
        f = open(filename,"a+")
        f.write("Fold {}\n".format(ii))
        f.close()
        print("Fold {}\n".format(ii))
        for attr_idx in range(11):
            cnn_model = Model(inputs=input1, outputs=sub_model1(input1))
            cnn2_model = Model(inputs=input1, outputs=cnn_model.layers[-2].output)
            cnn_model.compile(optimizer=Adam(learning_rate=1e-4), loss='mse', metrics=[pearson_r])
            history = cnn_model.fit(X_train[:, :, :, attr_idx], y_train, batch_size=bs, epochs=epoch, shuffle=True, verbose=2, callbacks=[TestCallback((X_test[:, :, :, attr_idx], y_test))])
            saved_hist.append(history)
            cnn_model.save("./skempi_PRelu_v6/spectcnn_skempi_model_{}_{}_{}.h5".format(j, ii, attr_idx))

            y_pred = cnn_model.predict(X_test[:,:,:,attr_idx])
            y_pred = np.reshape(y_pred, (len(y_pred),))

            mse = mean_squared_error(y_test, y_pred)
            rmse = sqrt(mse)
            pearcorr = sp.stats.pearsonr(y_test, y_pred)
            mae = mean_absolute_error(y_test, y_pred)
            f = open(filename,"a+")
            f.write("Index {:0>2d}: RMSE: {:.4f}, PCC: {:.4f}, MAE: {:.4f}\n".format(attr_idx, rmse, pearcorr[0], mae))
            f.close()
            print("Index {:0>2d}: RMSE: {:.4f}, PCC: {:.4f}, MAE: {:.4f}".format(attr_idx, rmse, pearcorr[0], mae))
            result[attr_idx, test_index] = y_pred

        ### Alpha L1 Laplacian 

        alpha_model = Model(inputs=input2, outputs=sub_model1(input2))
        alpha2_model = Model(inputs=input2, outputs=alpha_model.layers[-2].output)
        alpha_model.compile(optimizer=Adam(learning_rate=1e-4), loss='mse', metrics=[pearson_r])
        history = alpha_model.fit(alphal1_train, y_train, batch_size=bs, epochs=epoch, shuffle=True, verbose=2, callbacks=[TestCallback((alphal1_test, y_test))])
        saved_hist.append(history)
        alpha_model.save("./skempi_PRelu_v6/spectcnn_skempi_model_{}_{}_11.h5".format(j, ii))

        y_pred = alpha_model.predict(alphal1_test)
        y_pred = np.reshape(y_pred, (len(y_pred),))

        mse = mean_squared_error(y_test, y_pred)
        rmse = sqrt(mse)
        pearcorr = sp.stats.pearsonr(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)
        f = open(filename,"a+")
        f.write("Index {:0>2d}: RMSE: {:.4f}, PCC: {:.4f}, MAE: {:.4f}\n".format(11, rmse, pearcorr[0], mae))
        f.close()
        print("Index {:0>2d}: RMSE: {:.4f}, PCC: {:.4f}, MAE: {:.4f}".format(11, rmse, pearcorr[0], mae))
        result[11, test_index] = y_pred

        ### Alpha L2 Laplacian 

        alpha_model = Model(inputs=input2, outputs=sub_model1(input2))
        alpha2_model = Model(inputs=input2, outputs=alpha_model.layers[-2].output)
        alpha_model.compile(optimizer=Adam(learning_rate=1e-4), loss='mse', metrics=[pearson_r])
        history = alpha_model.fit(alphal2_train, y_train, batch_size=bs, epochs=epoch, shuffle=True, verbose=2, callbacks=[TestCallback((alphal2_test, y_test))])
        saved_hist.append(history)
        alpha_model.save("./skempi_PRelu_v6/spectcnn_skempi_model_{}_{}_12.h5".format(j, ii))

        y_pred = alpha_model.predict(alphal2_test)
        y_pred = np.reshape(y_pred, (len(y_pred),))

        mse = mean_squared_error(y_test, y_pred)
        rmse = sqrt(mse)
        pearcorr = sp.stats.pearsonr(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)
        f = open(filename,"a+")
        f.write("Index {:0>2d}: RMSE: {:.4f}, PCC: {:.4f}, MAE: {:.4f}\n".format(12, rmse, pearcorr[0], mae))
        f.close()
        print("Index {:0>2d}: RMSE: {:.4f}, PCC: {:.4f}, MAE: {:.4f}".format(12, rmse, pearcorr[0], mae))
        result[12, test_index] = y_pred

        # Plot All Training Loss

        print("spect_skempi_model saved")
        plt.figure()
        for attr_idx in range(13):
            plt.plot(saved_hist[attr_idx].history['loss'], linewidth=0.7, label=str(attr_idx))
            #plt.plot(history.history['pearson_r'], linewidth=0.7, label='PCC')
        plt.legend()
        plt.savefig("./skempi_PRelu_v6/spectnettree_skempi_loss_{}_{}.png".format(j, ii), dpi=200)

        ii += 1 

    np.save("./skempi_PRelu_v6/skempi_data_{}.npy".format(j), data)

    for index in range(len(idx)):
        rmse = np.sqrt(mean_squared_error(y,result[index]))
        pearsonr = scipy.stats.pearsonr(y,result[index])

        print('RMSE =', rmse)
        print('pearsonr =', pearsonr)

        f = open(filename,"a+")
        f.write("Iter {:0>2d}, Index {:0>2d}: RMSE: {:.4f} PCC: {:.4f} {}\n\n".format(j, index, rmse, pearsonr[0], pearsonr[1]))
        f.close()

        rmse_whole.append(rmse)
        pearsonr_whole.append(pearsonr[0])
        result_whole.append(result)

    np.save("./skempi_PRelu_v6/spectcnn_skempi_result_whole_{}.npy".format(j), result_whole)
    np.save("./skempi_PRelu_v6/spectcnn_skempi_rmse_whole_{}.npy".format(j), rmse_whole)
    np.save("./skempi_PRelu_v6/spectcnn_skempi_pearsonr_whole_{}.npy".format(j),pearsonr_whole)