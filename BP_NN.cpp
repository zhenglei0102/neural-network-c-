/*Author:zhenglei/215395
  Date:2017-08-07
  function:BP nerual network*/

#include "BP_NN.h"

BPNeuralNetwork::BPNeuralNetwork(int LayerNum, int *LayerSizes, float LearningRate)
{
	int i,j,k;
	int Tmpm,Tmpn = 0;
	m_LayersNum = LayerNum;
	m_LayerSizes = new int[m_LayersNum];
	for(i = 0; i < m_LayersNum; i++)
	{
		m_LayerSizes[i] = LayerSizes[i];
	}
	m_LearningRate = LearningRate;
	m_clouter = 0;

	m_Weights = new CMatrix<float>[m_LayersNum-1];
	m_Biases = new CMatrix<float>[m_LayersNum-1];
	m_DeltaWeights = new CMatrix<float>[m_LayersNum-1];
	
	for(i = 0; i < m_LayersNum-1; i++)
	{
		Tmpm = m_LayerSizes[i+1];
		Tmpn = m_LayerSizes[i];
		m_Weights[i] = CMatrix<float>(Tmpm,Tmpn);
		m_Biases[i] = CMatrix<float>(Tmpm,1);
		m_DeltaWeights[i] = CMatrix<float>(Tmpm,Tmpn);
		for(j = 0; j < Tmpm; j++)
		{
			for(k = 0; k < Tmpn; k++)
			{
				m_Weights[i][j][k] = gaussrand()/10.0;		
				m_DeltaWeights[i][j][k] = 0;
			}
			m_Biases[i][j][0] = gaussrand()/50.0;
		}
	}

	m_Inputs = new CMatrix<float>[m_LayersNum];
	m_Outputs = new CMatrix<float>[m_LayersNum];
	m_errors = new CMatrix<float>[m_LayersNum];

	for(i = 0; i < m_LayersNum; i++)
	{
		Tmpn = m_LayerSizes[i];
		m_Inputs[i] = CMatrix<float>(Tmpn,1);
		m_Outputs[i] = CMatrix<float>(Tmpn,1);
		m_errors[i] = CMatrix<float>(Tmpn,1);
		for(j = 0; j < Tmpn; j++)
		{
			m_Inputs[i][j][0] = 0;
			m_Outputs[i][j][0] = 0;
			m_errors[i][j][0] = 0;
		}
	}	
}
BPNeuralNetwork::~BPNeuralNetwork()
{

	delete []m_Weights;
	delete []m_Biases;
	delete []m_DeltaWeights;
	delete []m_Inputs;
	delete []m_Outputs;
	delete []m_errors;
	delete []m_LayerSizes;
}
CMatrix<float> BPNeuralNetwork::AcitiveFunc_linear(CMatrix<float> x)
{
	return x;
}
CMatrix<float> BPNeuralNetwork::AcitiveFunc_relu(CMatrix<float> x)
{
	int i,j;
	float alpha = 0.01f;
	int row = x.GetRow();
	int col = x.GetCol();
	CMatrix<float> Ret(row,col);
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			if(x[i][j] >= 0)
			{
				Ret[i][j] = x[i][j];
			}
			else
			{
				Ret[i][j] = alpha*x[i][j];
			}
		}
	}
	return Ret;
}
CMatrix<float> BPNeuralNetwork::AcitiveFunc_prime_linear(CMatrix<float> x)
{
	int i,j;
	int row = x.GetRow();
	CMatrix<float> Ret(row,row);
	for(i=0;i<row;i++)
	{
		for(j=0;j<row;j++)
		{
			if(i==j)
			{
				Ret[i][j] = 1;
			}
			else
			{
				Ret[i][j] = 0;
			}
		}
	}
	return Ret;
}
CMatrix<float> BPNeuralNetwork::AcitiveFunc_prime_relu(CMatrix<float> x)
{
	int i,j;
	float alpha = 0.01f;
	int row = x.GetRow();
	int col = x.GetCol();
	CMatrix<float> Ret(row,row);
	Ret.Zeros(row,row);
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			if(x[i][j] >= 0)
			{
				Ret[i][i] = 1.0f;
			}
			else
			{
				Ret[i][i] = alpha;
			}
		}
	}
	return Ret;
}
CMatrix<float> BPNeuralNetwork::feedforward(CMatrix<float> x)
{
	m_Inputs[0] = x;
	m_Outputs[0] = x;
	//hidden layer feed forward
	for(int i = 1; i < m_LayersNum-1;i++)
	{
		m_Inputs[i] = m_Weights[i-1]*m_Outputs[i-1] + m_Biases[i-1];
		m_Outputs[i] = AcitiveFunc_relu(m_Inputs[i]);
	}
	//output layer, activation function is different from hidden layer
	m_Inputs[m_LayersNum-1] = m_Weights[m_LayersNum-2]*m_Outputs[m_LayersNum-2] + m_Biases[m_LayersNum-2];
	m_Outputs[m_LayersNum-1] = AcitiveFunc_linear(m_Inputs[m_LayersNum-1]);
	return m_Outputs[m_LayersNum-1];
}
void BPNeuralNetwork::update_weights(CMatrix<float> x,CMatrix<float> y)
{
	CMatrix<float> Output = feedforward(x);
	//Output layer error
	m_errors[m_LayersNum-1] = AcitiveFunc_prime_linear(m_Inputs[m_LayersNum-1])*(Output-y);
	//cal deltaweights
	for(int i=m_LayersNum - 2;i>=0;i--)
	{
		m_errors[i] = AcitiveFunc_prime_relu(m_Inputs[i])*(~m_Weights[i]*m_errors[i+1]);				
		m_DeltaWeights[i] = m_LearningRate * (m_errors[i+1]*(~m_Outputs[i]));	
		m_Weights[i] = m_Weights[i] - m_DeltaWeights[i];
		
		m_Biases[i] = m_Biases[i] - m_LearningRate*m_errors[i+1];
	}
}

double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;          
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    phase = 1 - phase;
    return X;
}
/*
int main()
{
	float Tmperr = 0.0f;
	int LayerNum = 3;
	int LayerSizes[3] = {64,16,64};
	BPNeuralNetwork Mynn = BPNeuralNetwork(LayerNum,LayerSizes,0.001f);
	FILE*fp = NULL;
	int innode = 64;
	double Input[100][64];
	fp=fopen("BeamData.txt","r");
	for(int line = 0; line < 100; line ++)
	{
		for(int InNum = 0; InNum < innode; InNum++)
		{
			fscanf(fp,"%lf",&Input[line][InNum]);
		}
	}
	CMatrix<float> X(innode,1);
	for(int InNum = 0; InNum < innode; InNum++)
	{
		X[InNum][0] = Input[33][InNum];
	}
	for(int epoch = 0;epoch<100;epoch++)
	{
		Tmperr = 0.0f;
		Mynn.update_weights(X,X);
		for(int k =0;k<64;k++)
		{
			Tmperr += fabs(X[k][0] - Mynn.m_Outputs[2][k][0]);
		}
		if(epoch%20==0)
		{
			cout<<"est error = "<<Tmperr/64<<endl;
		}	
	}
	cout<<Mynn.m_Outputs[2]<<endl;
	return 0;
}
*/

