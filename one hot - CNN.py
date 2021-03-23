from keras.preprocessing import sequence
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.layers import Embedding
from keras.layers import Conv1D, GlobalMaxPooling1D

max_features = 5000
maxlen = 400
batch_size = 32
embedding_dims = 50
filters = 250
kernel_size = 3
hidden_dims = 250
epochs = 150



# read_data() is an user-defined function to load local data
x_train,x_test,y_train,y_test = load_data()
model = Sequential()
model.add(Embedding(max_features,
                embedding_dims,
                input_length=maxlen))
model.add(Dropout(0.2))
model.add(Conv1D(filters,
                kernel_size,
                padding='valid',
                activation='relu',
                strides=1))

model.add(GlobalMaxPooling1D())
model.add(Dense(hidden_dims))
model.add(Dropout(0.2))
model.add(Activation('relu'))
model.add(Dense(1))
model.add(Activation('sigmoid'))
model.compile(loss='binary_crossentropy',
                optimizer='adam',
                metrics=['accuracy'])
model.fit(x_train, y_train,
                batch_size=batch_size,
                epochs=epochs,
                validation_data=(x_test, y_test)
                )





