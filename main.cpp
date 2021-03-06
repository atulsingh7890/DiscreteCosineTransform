//
//  main.cpp
//  DiscreteCosineTransform
//
//  Created by Atul Singh on 7/30/16.
//  Copyright (c) 2016 Adobe. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<double> DoubleArray;

/*
 ** inVector : Double type vector whose dct needs to be computed
 ** outVector : This vector will contain the Result of the DCT
 ** inIndexUpto : How many values to be computed for DCT ,if this is greate than vector size
 **					it will compute upto the Vector size
 **					This will be 1 based index , 1 will compute only value DCT[0] and accordingly
 **					Or how many values needs to be filled in outVector
 */
void ComputeDCT_1(const DoubleArray inVector,DoubleArray &outVector,int inIndexUpto = -1);


/*
 ** inVector : Double type vector whose dct needs to be computed
 ** outVector : This vector will contain the Result of the DCT
 ** inIndexUpto : How many values to be computed for DCT ,if this is greate than vector size
 **					it will compute upto the Vector size
 **					This will be 1 based index , 1 will compute only value DCT[0] and accordingly
 **					Or how many values needs to be filled in outVector
 */
void ComputeDCT_2(const DoubleArray inVector,DoubleArray &outVector,int inIndexUpto = -1);

/*
 ** inVector : Double type vector whose dct needs to be computed
 ** outVector : This vector will contain the Result of the DCT
 ** inIndexUpto : How many values to be computed for DCT ,if this is greate than vector size
 **					it will compute upto the Vector size
 **					This will be 1 based index , 1 will compute only value DCT[0] and accordingly
 **					Or how many values needs to be filled in outVector
 */
void ComputeInverseDCT_2(const DoubleArray inVector,DoubleArray &outVector,int inIndexUpto = -1);

int main(int argc, const char * argv[]) {
	
	DoubleArray array;
	for(int i = 0; i < 4; ++i )
	{
		array.push_back(1.0);
	}
	
	DoubleArray dctValues;
	DoubleArray inverseDCTValues;
	
	ComputeDCT_2(array, dctValues );
	
	
	ComputeInverseDCT_2(dctValues, inverseDCTValues);
	
    return 0;
}


void ComputeDCT_1(const DoubleArray inVector,DoubleArray &outVector,int nIndexUpto)
{
	
	int nSize = inVector.size();
	
	if ( nIndexUpto <= 0 || nIndexUpto  >= nSize )
	{
		nIndexUpto = nSize;
	}
	
	outVector.resize(nSize,0.0);
	
	
	for(int i = 0; i < nSize && i < nIndexUpto ; ++i )
	{
		double x = inVector[i] + std::pow(1.0,i) * inVector[nSize - 1] ;
		
		for(int j = 1; j < nSize - 1; ++j )
		{
			x += inVector[j] *  std::cos( (M_PI /(nSize - 1)  ) *  j * i );
		}
		
		outVector[i] = x;
	}
}


void ComputeDCT_2(const DoubleArray inVector,DoubleArray &outVector,int nIndexUpto)
{
	
	int nSize = inVector.size();
	
	outVector.resize(nSize,0.0);
	
	double  scaleFactor = std::sqrt(0.5);
	
	if ( nIndexUpto <= 0 || nIndexUpto  >= nSize )
	{
		nIndexUpto = nSize;
	}
	
	outVector.resize(nSize,0.0);
	
	
	for(int i = 0; i < nSize && i < nIndexUpto ; ++i )
	{
		double x = 0;
	
		for(int j = 0; j < nSize ; ++j )
		{
			x += inVector[j] *  std::cos( (M_PI /nSize) * ( j + 0.5 ) * i );
		}
		
		if ( i ==  0 )
			x *= scaleFactor;
		
		x *= std::sqrt( 2.0/nSize );
		
		outVector[i] = x;
	}
}


void ComputeInverseDCT_2(const DoubleArray inVector,DoubleArray &outVector,int nIndexUpto)
{
	
	int nSize = inVector.size();
	
	outVector.resize(nSize,0.0);
	
	double  scaleFactor = 1.0 / std::sqrt(2.0);
	
	if ( nIndexUpto <= 0 || nIndexUpto  >= nSize )
	{
		nIndexUpto = nSize;
	}
	
	outVector.resize(nSize,0.0);
	
	double Ck = 0.0;
	
	for(int i = 0; i < nSize && i < nIndexUpto ; ++i )
	{
		double x = 0;
		double Ck = 0.0;
		
		for(int j = 0; j < nSize ; ++j )
		{
			Ck = ( j == 0 ) ? scaleFactor : 1.0;
 			x += Ck * inVector[j] *  std::cos( (M_PI /nSize) * ( i + 0.5 ) * j );
		}
		
		x *= std::sqrt( 2.0/nSize );
		
		outVector[i] = x;
	}
}

