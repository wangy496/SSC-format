import keras  
from __future__ import print_function
import keras
from keras.datasets import cifar10
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import losses
import numpy as np
from keras.engine.topology import print_layer_summary
from keras.utils.layer_utils import print_summary
import time
import random
from scipy import interp
import matplotlib.pyplot as plt



batch_size = 80
num_classes = 2
epochs = 200
data_augmentation = True
opt_adam = keras.optimizers.adam(lr=0.0001, rho=0.9, decay=1e-6)
# read_data() is an user-defined function to load local data
x_train, x_test, y_train, y_test = load_data()

x_train = x_train.reshape(len(x_train),60,60,3)
x_test = x_test.reshape(len(x_test),60,60,3)
y_train = keras.utils.to_categorical(y_trainb, 2)
y_test = keras.utils.to_categorical(y_testb, 2)
leakrelu = keras.layers.advanced_activations.LeakyReLU(alpha=0.3)
model = Sequential()
model.add(Dense(250, input_shape=x_train.shape[1:],init='uniform'))
model.add(Conv2D(128, (3, 3), padding='same',
                     input_shape=x_train.shape[1:]))
model.add(Activation(leakrelu))
model.add(Conv2D(128, (3, 3)))
model.add(Activation(leakrelu))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Conv2D(256, (3, 3), padding='same'))
model.add(Activation(leakrelu))
model.add(Conv2D(256, (3, 3)))
model.add(Activation(leakrelu))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(512))
model.add(Activation(leakrelu))
model.add(Dropout(0.25))
model.add(Dense(num_classes))
model.add(Activation('softmax'))

model.compile(loss='categorical_crossentropy',
                  optimizer=opt_adam,
                  metrics=['accuracy']
                  )

x_train = x_train.astype('float32')
x_test = x_test.astype('float32')

model.fit(x_train, y_train,
              batch_size=batch_size,
              epochs=epochs,
              validation_data=(x_test, y_test),
              shuffle=True)

model.summary()
