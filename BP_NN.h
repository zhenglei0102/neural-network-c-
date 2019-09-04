/*Author:zhenglei
  Date:2017-08-07
  function:BP nerual network*/
#include<iostream>
#include<math>
#include<stdlib>
#include "MatrixProcessing.h"

using namespace std;

double gaussrand();
class BPNeuralNetwork
{
public:
	int m_LayersNum;
	int *m_LayerSizes;
	float m_LearningRate;
	int m_clouter;
	
	CMatrix<float> *m_Weights;
	CMatrix<float> *m_Biases;
	CMatrix<float> *m_DeltaWeights;
	CMatrix<float> *m_Inputs;
	CMatrix<float> *m_Outputs;
	CMatrix<float> *m_errors;

	BPNeuralNetwork(int LayerNum, int *LayerSizes, float LearningRate);
	~BPNeuralNetwork();
	//feedforward process
	CMatrix<float> feedforward(CMatrix<float> x);
	//back propagtion process
	void update_weights(CMatrix<float> x,CMatrix<float> y);
	
protected:
	//Acitivation function
	CMatrix<float> AcitiveFunc_linear(CMatrix<float> x);		//linear function
	CMatrix<float> AcitiveFunc_prime_linear(CMatrix<float> x);	//linear function prime
	CMatrix<float> AcitiveFunc_relu(CMatrix<float> x);			//relu function
	CMatrix<float> AcitiveFunc_prime_relu(CMatrix<float> x);	//relu function prime
	
};

