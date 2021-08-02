import numpy as np
import os 
import csv
import sys 
import matplotlib.pyplot as plt
from keras.models import Model, load_model
from keras.layers.advanced_activations import PReLU
from keras.wrappers.scikit_learn import KerasRegressor
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

iter_num = int(sys.argv[1])

filename = "./ab_PRelu_Nout_v6/spectnettree_ab.log"
f = open(filename,"w+")
f.close()

#outliers = np.load("outliers.npy", allow_pickle=True)
#Y = np.load("./Y_ab.npy")
#y = np.asarray(Y, dtype = float)

#if option == 1:
    #y = np.delete(y, outliers)

#mean = np.mean(y)
#std = np.std(y)

n_num = 618
result_whole = np.zeros((10, 3, n_num))
rmse_whole = np.zeros((10, 3))
pearsonr_whole = np.zeros((10, 3))

for j in range(1, iter_num+1):
    print('iter_num',j)
    data = np.load("./ab_PRelu_Nout_v6/ab_data_{}.npy".format(j), allow_pickle=True)
    y = np.zeros(n_num)
    result = np.zeros((3, n_num))
    indices = []
    spectcnn_models = []
    spectnettree_models = []
    for ii in range(10):
        X_train, y_train, X_test, y_test, alphal1_train, alphal1_test, alphal2_train, alphal2_test, aux_train, aux_test, train_index, test_index = data[ii]
        y[test_index] = y_test

        f = open(filename,"a+")
        f.write("Fold {} Iteration {}\n".format(ii, j))
        print("Fold {} Iteration {}".format(ii,j))

        cnn_train, cnn_test = [], []
        for k in range(11):
            model_name = './ab_PRelu_Nout_v6/spectcnn_ab_model_{}_{}_{}.h5'.format(j, ii, k)
            spectcnn_models.append(load_model(model_name, custom_objects = {"PReLU":PReLU, 'rootmse':rootmse, 'pearson_r':pearson_r}))
            X_train, y_train, X_test, y_test, alphal1_train, alphal1_test, alphal2_train, alphal2_test, aux_train, aux_test, train_index, test_index = data[ii]
            X_train = spectcnn_models[-1].predict(X_train[:,:,:,k])
            X_test = spectcnn_models[-1].predict(X_test[:,:,:,k])
            X_train = np.reshape(X_train,(len(X_train),))
            X_test = np.reshape(X_test,(len(X_test),))

            cnn_train.append(X_train)
            cnn_test.append(X_test)

            mse = mean_squared_error(y_test, X_test)
            rmse= sqrt(mse)
            pearcorr =  sp.stats.pearsonr(y_test, X_test)
            mae = mean_absolute_error(y_test, X_test)   

            f.write("Index {:0>2d}: CNN Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(k, rmse, pearcorr[0], mae))
            print("Index {:0>2d}: CNN Features -> RMSE: {}, PCC: {}, MAE: {}".format(k, rmse, pearcorr[0], mae))

        ### Stacked Rips CNN Features
        params={'n_estimators': 4000, 'max_depth': 7, 'min_samples_split': 2,
            'learning_rate': 0.01, 'loss': 'ls','max_features':'sqrt','subsample':0.7}

        clf = ensemble.GradientBoostingRegressor(**params)
        clf.fit(np.transpose(cnn_train), y_train)
        y_pred = clf.predict(np.transpose(cnn_test))

        result[0, test_index] = y_pred

        mse = mean_squared_error(y_test, y_pred)
        rmse= sqrt(mse)
        pearcorr =  sp.stats.pearsonr(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)   
        f.write("Stacked Rips CNN Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(rmse, pearcorr[0], mae))
        print("Stacked Rips CNN Features -> RMSE: {}, PCC: {}, MAE: {}".format(rmse, pearcorr[0], mae))

        model_name = './ab_PRelu_Nout_v6/spectcnn_ab_model_{}_{}_{}.h5'.format(j, ii, 11)
        spectcnn_models.append(load_model(model_name, custom_objects = {"PReLU":PReLU, 'rootmse':rootmse, 'pearson_r':pearson_r}))
        X_train, y_train, X_test, y_test, alphal1_train, alphal1_test, alphal2_train, alphal2_test, aux_train, aux_test, train_index, test_index = data[ii]
        X_train = spectcnn_models[-1].predict(alphal1_train)
        X_test = spectcnn_models[-1].predict(alphal1_test)
        X_train = np.reshape(X_train,(len(X_train),))
        X_test = np.reshape(X_test,(len(X_test),))

        cnn_train.append(X_train)
        cnn_test.append(X_test)

        mse = mean_squared_error(y_test, X_test)
        rmse= sqrt(mse)
        pearcorr =  sp.stats.pearsonr(y_test, X_test)
        mae = mean_absolute_error(y_test, X_test)   

        f.write("Index {:0>2d}: CNN Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(11, rmse, pearcorr[0], mae))
        print("Index {:0>2d}: CNN Features -> RMSE: {}, PCC: {}, MAE: {}".format(11, rmse, pearcorr[0], mae))

        model_name = './ab_PRelu_Nout_v6/spectcnn_ab_model_{}_{}_{}.h5'.format(j, ii, 12)
        spectcnn_models.append(load_model(model_name, custom_objects = {"PReLU":PReLU, 'rootmse':rootmse, 'pearson_r':pearson_r}))
        X_train, y_train, X_test, y_test, alphal1_train, alphal1_test, alphal2_train, alphal2_test, aux_train, aux_test, train_index, test_index = data[ii]
        X_train = spectcnn_models[-1].predict(alphal2_train)
        X_test = spectcnn_models[-1].predict(alphal2_test)
        X_train = np.reshape(X_train,(len(X_train),))
        X_test = np.reshape(X_test,(len(X_test),))

        cnn_train.append(X_train)
        cnn_test.append(X_test)

        mse = mean_squared_error(y_test, X_test)
        rmse= sqrt(mse)
        pearcorr =  sp.stats.pearsonr(y_test, X_test)
        mae = mean_absolute_error(y_test, X_test)   

        f.write("Index {:0>2d}: CNN Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(12, rmse, pearcorr[0], mae))
        print("Index {:0>2d}: CNN Features -> RMSE: {}, PCC: {}, MAE: {}".format(12, rmse, pearcorr[0], mae))

        ### Stacked Alpha + Rips CNN Features

        clf = ensemble.GradientBoostingRegressor(**params)
        clf.fit(np.transpose(cnn_train), y_train)
        y_pred = clf.predict(np.transpose(cnn_test))

        result[1, test_index] = y_pred
        
        mse = mean_squared_error(y_test, y_pred)
        rmse= sqrt(mse)
        pearcorr =  sp.stats.pearsonr(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)   
        f.write("Stacked Alpha + Rips CNN Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(rmse, pearcorr[0], mae))
        print("Stacked Alpha + Rips CNN Features -> RMSE: {}, PCC: {}, MAE: {}".format(rmse, pearcorr[0], mae))

        clf = ensemble.GradientBoostingRegressor(**params)
        clf.fit(aux_train, y_train)
        aux_train = clf.predict(aux_train)
        aux_test = clf.predict(aux_test)

        cnn_train.append(aux_train)
        cnn_test.append(aux_test)

        mse = mean_squared_error(y_test, aux_test)
        rmse= sqrt(mse)
        pearcorr =  sp.stats.pearsonr(y_test, aux_test)
        mae = mean_absolute_error(y_test, aux_test)    
        f.write("Aux Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(rmse, pearcorr[0], mae))
        print("Aux Features -> RMSE: {}, PCC: {}, MAE: {}".format(rmse, pearcorr[0], mae))

        cnn_train = np.transpose(cnn_train)
        #cnn_train = np.concatenate((cnn_train, aux_train), axis=1)
        cnn_test = np.transpose(cnn_test)
        #cnn_test = np.concatenate((cnn_test, aux_test), axis=1)
        clf = ensemble.GradientBoostingRegressor(**params)
        clf.fit(cnn_train, y_train)
        y_pred = clf.predict(cnn_test)

        result[2, test_index] = y_pred
        
        mse = mean_squared_error(y_test, y_pred)
        rmse= sqrt(mse)
        pearcorr =  sp.stats.pearsonr(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)   
        f.write("Stacked All Features -> RMSE: {}, PCC: {}, MAE: {}\n".format(rmse, pearcorr[0], mae))
        print("Stacked All Features -> RMSE: {}, PCC: {}, MAE: {}".format(rmse, pearcorr[0], mae))
        f.close()

        ii += 1

    for jj in range(3):
        rmse = np.sqrt(mean_squared_error(y,result[jj]))
        pearsonr = scipy.stats.pearsonr(y,result[jj])

        print(jj, ' RMSE =', rmse)
        print(jj, ' pearsonr =', pearsonr)

        rmse_whole[j-1, jj] = rmse
        pearsonr_whole[j-1, jj] = pearsonr[0]
    
    result_whole[j-1] = result 

np.save("./ab_PRelu_Nout_v6/spectnettree_ab_result_whole", result_whole)
np.save("./ab_PRelu_Nout_v6/spectnettree_ab_rmse_whole", rmse_whole)
np.save("./ab_PRelu_Nout_v6/spectnettree_ab_pearsonr_whole",pearsonr_whole)

for j in range(3):
    print(j, " rmse mean:", np.mean(rmse_whole[:, j]))
    print(j, " pearsonr mean:",np.mean(pearsonr_whole[:, j]))
    print(j, " rmse std", np.std(rmse_whole[:, j]))
    print(j, " pearsonr std:",np.std(pearsonr_whole[:, j]))

file_out = open("./ab_PRelu_Nout_v6/crossvalidation_spectnettree_ab.txt","w+")
file_out.write("Pearsonr_average is:")
file_out.write(str(np.mean(pearsonr_whole)))
file_out.write("\n")
file_out.write("RMSE_average is:")
file_out.write(str(np.mean(rmse_whole)))
file_out.close()