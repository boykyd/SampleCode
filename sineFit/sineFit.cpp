// sineFit.cpp : 이 파일에는 'main' 함수가 포함됩니다. 거기서 프로그램 실행이 시작되고 종료됩니다.
//

#include <iostream>
#include <vector>

#define M_PI		3.14159265358979323846	/* pi */


using namespace::std;



bool	MatInverse3by3(float pfDst[][3], float pfSrc[][3]);
void	MatMul3by3(float Dst[][3], float A[][3], float B[][3]);
bool	FitPlaneEquation(float* a, float* b, float* c, float* x, float* y, float* z, int nPointNum);

float	CalcAVG(float* pData, int nNum);

float	NormalizedCrossCorrelation(float* pSrcBuff, float* pCmpBuff, int nLength);
bool	EstimatePeriod(float* pData, int nNum, int nWindowSize, int nMaxPadding, int nMargin, int* nPeriod);
//int		EstimatePeriod(float* pData, int nNum);



//1. Calc Period
//2. Estimate <Amplitude> 
//3. Calc Bias
//4. Check Exception


int main()
{

	int nPeriod = 100;
	int nAmplitude = 300;
	float fPhase = 0;
	float fSinData = 0;
	int numOfData = 1000;

	float* datax = new float[numOfData];
	float* datay = new float[numOfData];
	float* dataz = new float[numOfData];

	for (int i = 0; i < numOfData; i++)
	{
		
		fPhase = (float)(2 * M_PI * i) / (float)nPeriod+M_PI/2.0;
		fSinData = nAmplitude * sin(fPhase) +rand() % 50;
		
		//fSinData = rand()%50;

		dataz[i] = fSinData;

		printf("%lf \n", fSinData);
	}


	//for estimation
	float mag = 0;
	float phase = 0;
	float bias = 0;
	float a = 0;
	float b = 0;
	float c = 0;

	int nWindowSize = 500;
	int nMaxPadding = 400;
	int nMargin = 10;	//margin for period gap

	int nEstimatedPeriod = 0;


	EstimatePeriod(dataz, numOfData, nWindowSize, nMaxPadding, nMargin, &nEstimatedPeriod);

	for (int i = 0; i < numOfData; i++)
	{
		if (nEstimatedPeriod == 0)
		{
			datax[i] = cos(0);
			datay[i] = sin(0);
		}
		else
		{
			datax[i] = cos((2 * M_PI) * float(i % nEstimatedPeriod) / float(nEstimatedPeriod));
			datay[i] = sin((2 * M_PI) * float(i % nEstimatedPeriod) / float(nEstimatedPeriod));
		}
	}


	FitPlaneEquation(&a, &b, &c, datax, datay, dataz, numOfData);
	
	mag = sqrt(a * a + b * b);
	bias = CalcAVG(dataz, numOfData);
	phase = atan2(a, b);
	

	printf("period = %d \n", nEstimatedPeriod);
	printf("Magnititude = %lf \n", mag); 
	printf("Bias = %lf \n", bias);
	printf("Phase = %lf \n", phase);

	//for testing...
//	EstimatePeriod(dataz, numOfData);
//	CalcAutoCorrelation(dataz, numOfData);

	//Calc RMS Error



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


	if (!MatInverse3by3(InvA, A))
	{
		*a = 0;
		*b = 0;
		*c = 0;

		return false;
	}

//	MatMul3by3(TempMat1, InvA, A);		//TempMat1 = A.inv * A 
//	MatInverse3by3(InvMat, TempMat1);	//IvtMat = (A.inv * A).inv
//	MatMul3by3(TempMat2, InvMat, InvA);	//TempMat2 = (A.inv*A).inv * A.inv

//	*a = TempMat2[0][0] * B[0][0] + TempMat2[0][1] * B[1][0] + TempMat2[0][2] * B[2][0];
//	*b = TempMat2[1][0] * B[0][0] + TempMat2[1][1] * B[1][0] + TempMat2[1][2] * B[2][0];
//	*c = TempMat2[2][0] * B[0][0] + TempMat2[2][1] * B[1][0] + TempMat2[2][2] * B[2][0];

	*a = InvA[0][0] * B[0][0] + InvA[0][1] * B[1][0] + InvA[0][2] * B[2][0];
	*b = InvA[1][0] * B[0][0] + InvA[1][1] * B[1][0] + InvA[1][2] * B[2][0];
	*c = InvA[2][0] * B[0][0] + InvA[2][1] * B[1][0] + InvA[2][2] * B[2][0];

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
	float fAVGA		= 0;
	float fAVGB		= 0;
	float fSigmaA	= 0;
	float fSigmaB	= 0;
	float fSigmaAB	= 0;
	float fValueA	= 0;
	float fValueB	= 0;
	float fNCC		= 0;
	int nIdx		= 0;
	int i			= 0;

	if (nLength < 1)
		return 0;

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

	if (fSigmaAB == 0)
		return 0;

	for (i = 0; i < nLength; i++)
	{
		fValueA = (float)pSrcBuff[i];
		fValueB = (float)pCmpBuff[i];

		fNCC += ((fValueA - fAVGA) * (fValueB - fAVGB)) / (fSigmaAB);
	}

	fNCC = fNCC / nLength;

	return fNCC;
}//end of function 



//return value : *nPeriod
bool EstimatePeriod(float* pData, int nNum, int nWindowSize, int nMaxPadding, int nMargin, int* nPeriod)
{

//	int nWindowSize		= 500;
//	int nMaxPadding		= 400;
//  int nMargin			= 10;

	//0. Exception Check
	*nPeriod = 0;

	if (nNum < 1)	
		return false;
	if (nWindowSize < 1)
		return false;
	if (nMaxPadding < 1)
		return false;
	if (nNum < nWindowSize + nMaxPadding)
		return false;
	

	bool  bRetVal		= 0;
	float fCorrelateSum = 0;
	int i				= 0;
	int j				= 0;
	int nPadIdx			= 0;
	float* pSrc			= new float[nWindowSize];
	float* pCmp			= new float[nWindowSize];
	float* pArrCorrel	= new float[nMaxPadding];
	int nCorrelPos		= 0;


	vector<int> vecCrossPos;
	int nCrossPosCnt	= 0;
	
	float fHalfPeriodSum= 0;
	float fHalfPeriodCnt= 0;

	for (j = 0; j < nMaxPadding; j++)
	{
		for (i = 0; i < nWindowSize; i++)
		{
			pSrc[i] = pData[i];
			pCmp[i] = pData[i + j];
		}

		pArrCorrel[nCorrelPos++] = NormalizedCrossCorrelation(pSrc, pCmp, nWindowSize);
	}


	vecCrossPos.clear();

	for (i = 0; i < nMaxPadding-1; i++)
	{
		if (pArrCorrel[i] * pArrCorrel[i + 1] <= 0)
		{
			#ifdef _DEBUG
				printf("cross pos %d \n", i);
			#endif


			vecCrossPos.push_back(i);
			nCrossPosCnt++;
		}
	}

	for (i = 0; i < nCrossPosCnt-1; i++)
	{
		if (vecCrossPos[i + 1] - vecCrossPos[i] > nMargin)
		{
			fHalfPeriodSum = fHalfPeriodSum + vecCrossPos[i + 1] - vecCrossPos[i];
			fHalfPeriodCnt = fHalfPeriodCnt+1;
		}
	}

	if (fHalfPeriodCnt > 0)
	{
		*nPeriod = 2*(fHalfPeriodSum / fHalfPeriodCnt);

		bRetVal = true;
	}
	else
	{
		*nPeriod = 0;
		bRetVal = false;
	}


	//Search Zero Cross 
	//Search Max Pos


	delete[] pArrCorrel;
	delete[] pSrc;
	delete[] pCmp;

	return bRetVal;
}//end of function 
