folder：
	testData: it concludes the original test data
	TrainData: it concludes the original train data
	Predictions_NIF: it concludes the predictions of NIF model

dataTypeTransiation.ipynb: 
	the file is used to convert the data from .mat to .npz format.

NIF_BUR_train.mat & NIF_BUR_train.npz: 
	the data is for training, it is n*4 matrix, n means the number of samples, and columns are v values, times, positions 	and the results of corresponding coordinates.

NIF_BUR_test_3.mat & NIF_BUR_test_3.npz: 
	the data is for testing. it includes 2 samples, first is the extrapolation testing data and the second is the 	interpolation testing data.

Burger_testingRresults_20241201.mat: 
	the data is the predictions of testing data according to NIF model.