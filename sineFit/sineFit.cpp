// sineFit.cpp : 이 파일에는 'main' 함수가 포함됩니다. 거기서 프로그램 실행이 시작되고 종료됩니다.
//

#include <iostream>

#define M_PI		3.14159265358979323846	/* pi */


bool	MatInverse3by3(float pfDst[][3], float pfSrc[][3]);
void	MatMul3by3(float Dst[][3], float A[][3], float B[][3]);
bool	FitPlaneEquation(float* a, float* b, float* c, float* x, float* y, float* z, int nPointNum);

float	CalcAVG(float* pData, int nNum);

float	NormalizedCrossCorrelation(float* pSrcBuff, float* pCmpBuff, int nLength);
int		EstimatePeriod(float* pData, int nNum);



//1. Calc Period
//2. Estimate <Amplitude> 
//3. Calc Bias
//4. Check Exception


int main()
{

	int nPeriod		= 100;
	int nAmplitude	= 300;
	float fPhase	= 0;
	float fSinData	= 0;
	int numOfData	= 1000;

	float* datax = new float[numOfData];
	float* datay = new float[numOfData];
	float* dataz = new float[numOfData];

	FILE* fp;

	fopen_s(&fp, "test.txt", "wt");

	for (int i = 0; i < numOfData; i++)
	{
		fPhase = (float)(2 * M_PI * i) / (float)nPeriod;
		fSinData = nAmplitude * sin(fPhase);// +rand() % 50;

		dataz[i] = fSinData;

		printf("%lf \n", fSinData);
		fprintf(fp, "%lf \n", fSinData);
	}

	fclose(fp);



	//for estimation
	float mag	= 0;
	float a		= 0;
	float b		= 0;
	float c		= 0;


	//float fAvg = CalcAVG(dataz, numOfData);

	//int nEstimatedPeriod = 300;
	int nEstimatedPeriod = EstimatePeriod(dataz, numOfData);


	for (int i = 0; i < numOfData; i++)
	{
		if (nEstimatedPeriod == 0)
		{
//			datax[i] = 1;
//			datay[i] = 0;
		}
		else
		{ 
			datax[i] = cos((2 * M_PI) * float(i % nEstimatedPeriod) / float(nEstimatedPeriod));
			datay[i] = sin((2 * M_PI) * float(i % nEstimatedPeriod) / float(nEstimatedPeriod));
		}
	}



	FitPlaneEquation(&a, &b, &c, datax, datay, dataz, numOfData);

	mag = sqrt(a * a + b * b);

	printf("magnititude = %lf \n", mag); 


	//for testing...
	EstimatePeriod(dataz, numOfData);
//	CalcAutoCorrelation(dataz, numOfData);



	delete [] datax;
	delete [] datay;
	delete [] dataz;

}//end of main





//Z=AX+BY+C
bool FitPlaneEquation(float* a, float* b, float* c, float* x, float* y, float* z, int nPointNum)
{

	if (nPointNum < 6)
		return false;

	float A[3][3];
	float B[3][1];
	float X[3][1];

	float InvA[3][3];
	float InvMat[3][3];
	float TempMat1[3][3];
	float TempMat2[3][3];

	memset(A, 0, sizeof(float) * 9);
	memset(B, 0, sizeof(float) * 3);
	memset(X, 0, sizeof(float) * 3);

	for (int i = 0; i < nPointNum; i++)
	{
		A[0][0] += x[i] * x[i];
		A[0][1] += x[i] * y[i];
		A[0][2] += x[i];

		A[1][0] += x[i] * y[i];
		A[1][1] += y[i] * y[i];
		A[1][2] += y[i];

		A[2][0] += x[i];
		A[2][1] += y[i];
		A[2][2] += 1;

		B[0][0] += x[i] * z[i];
		B[1][0] += y[i] * z[i];
		B[2][0] += z[i];

	}


	MatInverse3by3(InvA, A);			//A.inv
	MatMul3by3(TempMat1, InvA, A);		//TempMat1 = A.inv * A 
	MatInverse3by3(InvMat, TempMat1);	//IvtMat = (A.inv * A).inv
	MatMul3by3(TempMat2, InvMat, InvA);	//TempMat2 = (A.inv*A).inv * A.inv

	*a = TempMat2[0][0] * B[0][0] + TempMat2[0][1] * B[1][0] + TempMat2[0][2] * B[2][0];
	*b = TempMat2[1][0] * B[0][0] + TempMat2[1][1] * B[1][0] + TempMat2[1][2] * B[2][0];
	*c = TempMat2[2][0] * B[0][0] + TempMat2[2][1] * B[1][0] + TempMat2[2][2] * B[2][0];

	return true;

}//end of function 



bool MatInverse3by3(float pfDst[][3], float pfSrc[][3])
{
	float fA = pfSrc[1][1] * pfSrc[2][2] - pfSrc[2][1] * pfSrc[1][2];
	float fB = pfSrc[1][2] * pfSrc[2][0] - pfSrc[1][0] * pfSrc[2][2];
	float fC = pfSrc[1][0] * pfSrc[2][1] - pfSrc[1][1] * pfSrc[2][0];
	float fD = pfSrc[0][2] * pfSrc[2][1] - pfSrc[0][1] * pfSrc[2][2];
	float fE = pfSrc[0][0] * pfSrc[2][2] - pfSrc[0][2] * pfSrc[2][0];
	float fF = pfSrc[0][1] * pfSrc[2][0] - pfSrc[0][0] * pfSrc[2][1];
	float fG = pfSrc[0][1] * pfSrc[1][2] - pfSrc[0][2] * pfSrc[1][1];
	float fH = pfSrc[0][2] * pfSrc[1][0] - pfSrc[0][0] * pfSrc[1][2];
	float fK = pfSrc[0][0] * pfSrc[1][1] - pfSrc[0][1] * pfSrc[1][0];
	float fDet = pfSrc[0][0] * fA + pfSrc[0][1] * fB + pfSrc[0][2] * fC;

	if (fDet == 0.f)
	{
		return false;
	}

	pfDst[0][0] = fA / fDet;	pfDst[0][1] = fD / fDet;	pfDst[0][2] = fG / fDet;
	pfDst[1][0] = fB / fDet;	pfDst[1][1] = fE / fDet;	pfDst[1][2] = fH / fDet;
	pfDst[2][0] = fC / fDet;	pfDst[2][1] = fF / fDet;	pfDst[2][2] = fK / fDet;

	return true;
}//end of function



void MatMul3by3(float Dst[][3], float A[][3], float B[][3])
{
	int i, j, k;

	memset(Dst, 0, sizeof(float) * 9);

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				Dst[i][j] += A[i][k] * B[k][j];
			}
		}
	}

}//end of function 



float CalcAVG(float* pData, int nNum)
{
	float fSum = 0;

	for (int i = 0; i < nNum; i++)
	{
		fSum = pData[i] + fSum;
	}

	if (nNum == 0)
		return 0;

	
	fSum = fSum / (float)nNum;

	return fSum;
}//end of function 



float NormalizedCrossCorrelation(float* pSrcBuff, float* pCmpBuff, int nLength)
{
	float fAVGA = 0;
	float fAVGB = 0;
	float fSigmaA = 0;
	float fSigmaB = 0;
	float fSigmaAB = 0;
	float fValueA = 0;
	float fValueB = 0;
	float fNCC = 0;
	int nIdx = 0;
	int i = 0;

	for (i = 0; i < nLength; i++)
	{
		fValueA = (float)pSrcBuff[i];
		fValueB = (float)pCmpBuff[i];

		fAVGA += fValueA;
		fAVGB += fValueB;

		fSigmaA += fValueA * fValueA;
		fSigmaB += fValueB * fValueB;
	}

	fAVGA = fAVGA / nLength;
	fAVGB = fAVGB / nLength;

	fSigmaA = sqrt(fSigmaA / nLength - fAVGA * fAVGA);
	fSigmaB = sqrt(fSigmaB / nLength - fAVGB * fAVGB);
	fSigmaAB = fSigmaA * fSigmaB;

	for (i = 0; i < nLength; i++)
	{
		fValueA = (float)pSrcBuff[i];
		fValueB = (float)pCmpBuff[i];

		fNCC += ((fValueA - fAVGA) * (fValueB - fAVGB)) / (fSigmaAB);
	}

	fNCC = fNCC / nLength;

	return fNCC;
}//end of function 



int EstimatePeriod(float* pData, int nNum)
{
	int nWindowSize		= 500;
	int nMaxPadding		= 400;
	float fCorrelateSum = 0;
	int j				= 1;
	int nPadIdx			= 0;
	float* pSrc			= new float[nWindowSize];
	float* pCmp			= new float[nWindowSize];


	//float* pArrCorrel = new float[nNum];
	float pArrCorrel[1000];

	int nCorrelPos = 0;

	for (int j = 0; j < nMaxPadding; j++)
	{
		for (int i = 0; i < nWindowSize; i++)
		{
			pSrc[i] = pData[i];
			pCmp[i] = pData[i + j];
		}

		pArrCorrel[nCorrelPos++] = NormalizedCrossCorrelation(pSrc, pCmp, nWindowSize);
	}


	for (int i = 0; i < nWindowSize - 1; i++)
	{
		if (pArrCorrel[i] * pArrCorrel[i + 1] <= 0)
		{
			printf("cross pos %d \n", i);
		}
	}



	//Search Zero Cross 
	//Search Max Pos

	//	delete[] pArrCorrel;

	return 0;
}//end of function 
