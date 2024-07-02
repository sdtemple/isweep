
# did not implement the hdp approach

import os
import gzip
import _pickle as pickle

import numpy as np
import scipy.stats
import arviz

import tensorflow as tf
from tensorflow import keras
from keras import models, layers, activations, optimizers, regularizers
# from keras.utils.vis_utils import plot_model
from keras.models import load_model

import skimage.transform
import itertools
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pydot

import pymc

import arviz

import sys

exec(open('/home/guests/sdtemple/brwn/seth/imagene/ImaGene/ImaGene.py').read())

gene_with_classes=sys.argv[6]
Ne=sys.argv[5]
nr_samples = sys.argv[4]
VCF_file_name = sys.argv[3]
fileout=sys.argv[2]
model_file_name = sys.argv[1]
nr_samples=int(nr_samples)
Ne2 = 2 * int(float(Ne))

file_LCT = ImaFile(nr_samples=nr_samples, VCF_file_name=VCF_file_name);
gene_LCT = file_LCT.read_VCF();
gene_LCT.filter_freq(0.01);
gene_LCT.sort('rows_freq');
gene_LCT.sort('cols_freq');
gene_LCT.resize((128, 128));
gene_LCT.convert(flip=True);
# gene_LCT.plot();
gene_LCT.summary();

mygene = load_imagene(gene_with_classes)

model = load_model(model_file_name);

# print(model.predict(gene_LCT.data, batch_size=None)[0][0])
# print(model.predict(gene_LCT.data, batch_size=None))

val = model.predict(gene_LCT.data, batch_size=None)[0][0]

probs = model.predict(gene_LCT.data, batch_size=None)[0]
print(probs)
samples_distr = np.random.choice(mygene.classes, size=10000, replace=True, p=probs)
print(samples_distr[:100])
HPD = arviz.hdi(samples_distr,hdi_prob=0.95)
#print(gene_LCT)
print(np.argmax(np.array(probs)))
MAP = mygene.classes[np.argmax(np.array(probs))]
MLE = np.average(mygene.classes, weights=probs)
print(HPD)
#print(MAP)
print(MLE)


f=open(fileout,'w')
f.write('Ne\tNe2\tMAP\tMAP2\tAVG\tAVG2\tHPD1l\tHPD1u\tHPD2l\tHPD2u\n')
#f.write(str(val))
#f.write('\t')
f.write(str(Ne))
f.write('\t')
f.write(str(Ne2))
f.write('\t')
f.write(str(MAP/float(Ne))); f.write('\t')
f.write(str(MAP/float(Ne2))); f.write('\t')
f.write(str(MLE/float(Ne))); f.write('\t')
f.write(str(MLE/float(Ne2))); f.write('\t')
f.write(str(HPD[0]/float(Ne))); f.write('\t')
f.write(str(HPD[0]/float(Ne2))); f.write('\t')
f.write(str(HPD[1]/float(Ne))); f.write('\t')
f.write(str(HPD[1]/float(Ne2))); f.write('\n')
#f.write(str(val/float(Ne)))
#f.write('\t')
#f.write(str(val/float(Ne2)))
#f.write('\n')
f.close()
