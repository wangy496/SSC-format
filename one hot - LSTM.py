
from keras.preprocessing import sequence
from keras.models import Sequential
from keras.layers import Dense, Embedding
from keras.layers import LSTM
def load_data():
    import numpy as np
    x_all = [[[0 for j in range(16)] for i in range(40)] for k in range(100)]
    y_all = [0 for s in range(100)]
    x_train, x_test, y_train, y_test = x_all[:80],x_all[80:],y_all[:80],y_all[80:]
    x_train = np.reshape(x_train,(len(x_train), 40, 16))
    x_test = np.reshape(x_test,(len(x_test), 40, 16))
    x_train = x_train.astype('float32')
    x_test = x_test.astype('float32')

    y_train = np.array(y_train)
    y_test = np.array(y_test)
    return x_train, x_test, y_train, y_test

for i in range(1):

    max_features = 200
    maxlen = 40
    batch_size = 32
    # read_data() is an user-defined function to load local data
    x_train,x_test,y_train,y_test = load_data()
    #x_train = sequence.pad_sequences(x_train, maxlen=maxlen)
    #x_test = sequence.pad_sequences(x_test, maxlen=maxlen)
    print(x_train.shape)
    print(x_test.shape)
    model = Sequential()
    model.add(Embedding(max_features, 128))
    model.add(LSTM(128, dropout=0.2, recurrent_dropout=0.2,return_sequences=True))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy',
                      optimizer='adam',
                      metrics=['accuracy'])
    #model.load_weights('xx_{0}.h5'.format(i))
    model.fit(x_train, y_train,
                  batch_size=batch_size,
                  epochs=150,
                  validation_data=(x_test, y_test))

    #model.save_weights('xx_{0}.h5'.format(i+1))

