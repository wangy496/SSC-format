This code belongs to the "Protein-protein interaction prediction using 2-D data format based on deep learning" paper manuscript.
It is slightly simplified implementation of Convolutional Neural Networks for protein-protein interaction prediction paper.
Requirements
Python 3.5
Tensorflow > 0.12
Numpy

This code package includes 4 main file:
model.py: The main training and testing code. After editing the input data (as load_data() function ) the model can be trained properly according to the parameters. 
methods.py: The get_CNN_data_3_channel() function is used to construct the three channels. An example is provided to simulate a small part of entire dataset.
one hot-CNN.py: A general CNN model to deal 1-D data.
one hot-LSTM.py:A general LSTM model to deal 1-D data.