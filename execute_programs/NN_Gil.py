# -*- coding: utf-8 -*-
"""
Created on Sun May 24 11:21:39 2020

@author: gillo
"""
import sys
sys.path.append("/groups/itay_mayrose/danaazouri/PhyAI/code/")
import warnings
warnings.filterwarnings("ignore")			# TEMP


from defs import *


import sys
import logging

# logging.warning('start')

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import random
#'''
from keras.layers import Input, Dense
from keras.models import Model
from keras.models import Sequential
from keras.layers import Dense
from tensorflow.keras import regularizers
#'''



def apply_NN_model(X_train, Y_train, X_test, Y_test, X_val, Y_val):
	model = Sequential()
	model.add(Dense(30, input_dim=29, activation='relu'))
	model.add(Dense(30, activation='relu'))
	model.add(Dense(8, activation='relu'))
	model.add(Dense(1, activation='linear'))
	model.compile(loss='mean_squared_error', optimizer='adam', metrics=['accuracy'])
	model.fit(X_train, Y_train, epochs=100, batch_size=4096, validation_data=(X_val, Y_val))

	Y_test_hat = model.predict(X_test)
	#plt.plot(Y_test, Y_test_hat, 'x')
	#plt.savefig(SUMMARY_FILES_DIR + "NN_fig1.png")

	Y_train_hat = model.predict(X_train)
	#plt.plot(Y_train, Y_train_hat, 'x')
	#plt.savefig(SUMMARY_FILES_DIR + "NN_fig2.png")

	return



def get_datasets(indices, col_name_lst):
	s = df.loc[df[FEATURES["group_id"]].isin(indices)][col_name_lst].values
	return s

def standard_reshape_y(y_set, ref):
	return ((y_set-ref.mean())/ref.std()).reshape(-1, 1)


def standard_reshape_x(x_set, ref):
	return (x_set - ref.mean(axis=0)) / ref.std(axis=0)


def prepro(df, test_fraction, val_fraction):
	x_col_list = df.columns[6:-1]
	y_col_name = df.columns[-1]

	groups_ids = df[FEATURES["group_id"]].unique()
	indices_val = groups_ids[-int(val_fraction * len(groups_ids)):]
	groups_ids_up = groups_ids[:-int(val_fraction * len(groups_ids))]
	indices_test = np.random.choice(groups_ids_up, int(test_fraction * len(groups_ids_up)))
	indices_train = np.setdiff1d(groups_ids_up, indices_test)

	X_val_orig, Y_val_orig = get_datasets(indices_val, x_col_list), get_datasets(indices_val, y_col_name)
	X_test_orig, Y_test_orig = get_datasets(indices_test, x_col_list), get_datasets(indices_test, y_col_name)
	X_train_orig, Y_train_orig = get_datasets(indices_train, x_col_list), get_datasets(indices_train, y_col_name)

	Y_val, Y_test, Y_train = standard_reshape_y(Y_val_orig, Y_train_orig), standard_reshape_y(Y_test_orig, Y_train_orig), standard_reshape_y(
		Y_train_orig, Y_train_orig)
	X_val, X_test, X_train = standard_reshape_x(X_val_orig, X_train_orig), standard_reshape_x(X_test_orig, X_train_orig), standard_reshape_x(X_train_orig, X_train_orig)

	#X_val = X_val.reshape([X_val.shape[0], X_val.shape[1], 1])
	#X_test = X_test.reshape([X_test.shape[0], X_test.shape[1], 1])
	#X_train = X_train.reshape([X_train.shape[0], X_train.shape[1], 1])


	return X_train, Y_train, X_test, Y_test, X_val, Y_val


if __name__ == '__main__':
	df = pd.read_csv("D:/Users/Administrator/Desktop/learning_subset_1000ds.csv") #, nrows=100000)
	#df = pd.read_csv(SUMMARY_FILES_DIR + "learning_subset_1000ds.csv")
	#df = pd.read_csv(SUMMARY_FILES_DIR + LEARNING_DATA.format("all_moves", "1"))
	df[df.columns[-1]] = df[df.columns[-1]].apply(lambda x: np.log10(x + 1E-9))
	df = df.dropna()

	X_train, Y_train, X_test, Y_test, X_val, Y_val = prepro(df, test_fraction = 0.25, val_fraction = 0.05)
	apply_NN_model(X_train, Y_train, X_test, Y_test, X_val, Y_val)



