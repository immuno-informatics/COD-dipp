import argparse
import sys

import numpy as np
import pandas as pd

import h5py
import scipy.io
import theano
from keras.constraints import maxnorm
from keras.layers.convolutional import Convolution1D, MaxPooling1D
from keras.layers.core import Activation, Dense, Dropout, Flatten
from keras.layers.recurrent import GRU, LSTM
from keras.models import Sequential
from keras.optimizers import RMSprop
from keras.preprocessing import sequence
from keras.regularizers import l1, l2
from sklearn.metrics import average_precision_score, roc_auc_score


def args():
    parser = argparse.ArgumentParser(
        description='Predict Translation Initiation Sites')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '-input_prediction',
        type=str,
        help='Path to input for prediction 1 sequence per line',
        required=True)
    requiredNamed.add_argument(
        '-input_model_dir',
        type=str,
        help="""Path to TITER model direcotry,
should contain:
    - dict_piror_front_Gaotrain.npy
    - data/neg_seq_test_all_upstream.npy
    - data/pos_seq_test.npy
    - model/bestmodel_*.hdf5""",
        required=True)
    requiredNamed.add_argument(
        '-output',
        type=str,
        help='Path to output file',
        required=True)

    return parser


# np.random.seed(1337) # for reproducibility


def seq_matrix(seq_list, label):
    tensor = np.zeros((len(seq_list), 203, 8))
    for i in range(len(seq_list)):
        seq = seq_list[i]
        j = 0
        for s in seq:
            if s == 'A' and (j < 100 or j > 102):
                tensor[i][j] = [1, 0, 0, 0, 0, 0, 0, 0]
            if s == 'T' and (j < 100 or j > 102):
                tensor[i][j] = [0, 1, 0, 0, 0, 0, 0, 0]
            if s == 'C' and (j < 100 or j > 102):
                tensor[i][j] = [0, 0, 1, 0, 0, 0, 0, 0]
            if s == 'G' and (j < 100 or j > 102):
                tensor[i][j] = [0, 0, 0, 1, 0, 0, 0, 0]
            if s == '$':
                tensor[i][j] = [0, 0, 0, 0, 0, 0, 0, 0]
            if s == 'A' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 1, 0, 0, 0]
            if s == 'T' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 0, 1, 0, 0]
            if s == 'C' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 0, 0, 1, 0]
            if s == 'G' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 0, 0, 0, 1]
            j += 1
    if label == 1:
        y = np.ones((len(seq_list), 1))
    else:
        y = np.zeros((len(seq_list), 1))
    return tensor, y


if __name__=="__main__":
    arguments = args().parse_args()
    input_file = arguments.input_prediction
    output_file = arguments.output
    root = arguments.input_model_dir
    if root[-1] != "/":
        root += "/"

    codon_tis_prior = np.load(
        root + 'dict_piror_front_Gaotrain.npy', allow_pickle=True)
    codon_tis_prior = codon_tis_prior.item()

    codon_list = []
    for c in codon_tis_prior.keys():
        if codon_tis_prior[c] != 'never' and codon_tis_prior[c] >= 1:
            codon_list.append(c)

    print 'Loading test data...'
    pos_seq_test = np.load(root + 'data/pos_seq_test.npy', allow_pickle=True)
    neg_seq_test = np.load(
        root + 'data/neg_seq_test_all_upstream.npy', allow_pickle=True)
    pos_codon = []
    neg_codon = []

    for s in pos_seq_test:
        if s[100:103] in codon_list:
            pos_codon.append(codon_tis_prior[s[100:103]])

    for s in neg_seq_test:
        if s[100:103] in codon_list:
            neg_codon.append(codon_tis_prior[s[100:103]])

    pos_codon = np.array(pos_codon)
    neg_codon = np.array(neg_codon)
    codon = np.concatenate((pos_codon, neg_codon)).reshape(
        (len(pos_codon) + len(neg_codon), 1))

    pos_seq_test1 = []
    neg_seq_test1 = []
    for s in pos_seq_test:
        if s[100:103] in codon_list:
            pos_seq_test1.append(s)

    for s in neg_seq_test:
        if s[100:103] in codon_list:
            neg_seq_test1.append(s)

    print str(len(pos_seq_test1)) + ' positive test data loaded...'
    print str(len(neg_seq_test1)) + ' negative test data loaded...'

    pos_test_X, pos_test_y = seq_matrix(seq_list=pos_seq_test1, label=1)
    neg_test_X, neg_test_y = seq_matrix(seq_list=neg_seq_test1, label=0)
    X_test = np.concatenate((pos_test_X, neg_test_X), axis=0)
    y_test = np.concatenate((pos_test_y, neg_test_y), axis=0)

    print 'Building model...'
    model = Sequential()
    model.add(
        Convolution1D(
            nb_filter=128,
            filter_length=3,
            input_dim=8,
            input_length=203,
            border_mode='valid',
            W_constraint=maxnorm(3),
            activation='relu',
            subsample_length=1))
    model.add(MaxPooling1D(pool_length=3))
    model.add(Dropout(p=0.21370950078747658))
    model.add(LSTM(output_dim=256, return_sequences=True))
    model.add(Dropout(p=0.7238091317104384))
    model.add(Flatten())
    model.add(Dense(1))
    model.add(Activation('sigmoid'))

    print 'Compiling model...'
    model.compile(
        loss='binary_crossentropy', optimizer='nadam', metrics=['accuracy'])

    # uncomment to produce the publication results

    print 'Skipping prediction on test data...'
    # print 'Predicting on test data...'
    # y_test_pred_n = np.zeros((len(y_test), 1))
    # y_test_pred_p = np.zeros((len(y_test), 1))

    # for i in range(32):
    #     model.load_weights(root + 'model/bestmodel_' + str(i) + '.hdf5')
    #     y_test_pred = model.predict(X_test, verbose=1)
    #     y_test_pred_n += y_test_pred
    #     y_test_pred_p += y_test_pred * codon

    # y_test_pred_n = y_test_pred_n / 32
    # y_test_pred_p = y_test_pred_p / 32

    ######################################

    # added code

    # reasoning: the optimal cut off point would be where "true positive rate" is high
    # and the "false positive rate" is low

    # from sklearn import metrics
    # fpr, tpr, thresholds = metrics.roc_curve(y_test,y_test_pred_p)
    # optimal_idx = np.argmax(tpr - fpr)
    # optimal_threshold = thresholds[optimal_idx]
    # optimal_threshold
    # 0.9965503307298548
    # end of added code

    # print 'Perf without prior, AUC: ' + str(roc_auc_score(y_test, y_test_pred_n))
    # print 'Perf without prior, AUPR: ' + str(
    #     average_precision_score(y_test, y_test_pred_n))
    # print 'Perf with prior, AUC: ' + str(roc_auc_score(y_test, y_test_pred_p))
    # print 'Perf with prior, AUPR: ' + str(
    #     average_precision_score(y_test, y_test_pred_p))

    # added code for prediction based on the prediction above
    # prediction on the data
    print 'Loading user prediction data...'
    with open(input_file) as fh:
        pos_seq_test = np.array([line.strip() for line in fh], dtype="|S203")

    pos_codon = []
    for s in pos_seq_test:
        if s[100:103] in codon_list:
            pos_codon.append(codon_tis_prior[s[100:103]])

    pos_codon = np.array(pos_codon)
    codon = pos_codon.reshape((len(pos_codon), 1))

    pos_seq_test1 = []
    for s in pos_seq_test:
        if s[100:103] in codon_list:
            pos_seq_test1.append(s)

    print str(len(pos_seq_test1)) + ' prediction data loaded...'

    pos_test_X, pos_test_y = seq_matrix(seq_list=pos_seq_test1, label=1)

    X_test = pos_test_X
    y_test = pos_test_y
    y_test_pred_p = np.zeros((len(y_test), 1))

    print 'Predicting on user data...'
    for i in range(32):
        model.load_weights(root + 'model/bestmodel_' + str(i) + '.hdf5')
        y_test_pred = model.predict(X_test, verbose=1)
        y_test_pred_p += y_test_pred * codon

    y_test_pred_p = y_test_pred_p / 32

    optimal_threshold = 0.9965503307298548

    y_test_pred_p = y_test_pred_p[:, 0]
    try:
        df = pd.DataFrame({
            "titer seq": list(pos_seq_test),
            "titer score": y_test_pred_p,
            "pass threshold": y_test_pred_p >= optimal_threshold
        })
    except:
        df = pd.DataFrame({
            "titer seq": pos_seq_test1,
            "titer score": y_test_pred_p,
            "pass threshold": y_test_pred_p >= optimal_threshold
        })

    df.to_csv(output_file, sep="\t", header=True, index=False)