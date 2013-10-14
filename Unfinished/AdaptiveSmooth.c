/***************************************************************************
                          AdaptiveSmooth.c  -  description
                             -------------------
    copyright            : (C) 2002-2012 by Michaël Zugaro
    email                : michael.zugaro@college-de-france.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "mex.h"

/* #define DEBUG */

int      	xMax,yMax,max,i0;
double   	radius,radius2,alpha,alpha2,dt,dt2;
int	   	x0,y0,x,y;
double		dx2,dy2,wrapped;
double		*countMap,*timeMap,*sCountMap,*sTimeMap,*radiusMap;
bool     	circX,circY,*visited;

void ProcessBin(int x,int y)
{
	int i;
    
	#ifdef DEBUG
	mexPrintf("  (%d,%d)",x,y);
	#endif

	/* Make sure (x,y) lies within the map - wrap circular variables */
	if ( x < 1 )
		if ( circX ) x = xMax;
		else
		{
			#ifdef DEBUG
			mexPrintf(" [out]\n");
			#endif
			return;			
		}
	if ( x > xMax )
		if ( circX ) x = 1;
		else
		{
			#ifdef DEBUG
			mexPrintf(" [out]\n");
			#endif
			return;
		}
	if ( y < 1 )
		if ( circY ) y = yMax;
		else
		{
			#ifdef DEBUG
			mexPrintf(" [out]\n");
			#endif
			return;
		}			
	if ( y > yMax )
		if ( circY ) y = 1;
		else
		{
			#ifdef DEBUG
			mexPrintf(" [out]\n");
			#endif
			return;
		}

	/* Have we already visited this bin? */
	i = (int) (x-1)*yMax+(y-1);
	if ( visited[i] )
	{
		#ifdef DEBUG
		mexPrintf(" [visited]\n");
		#endif
		return;
	}
	visited[i] = true;

	/* Compute (squared) distance along x */
	dx2 = (x0-x)*(x0-x);
	if ( circX )
	{
		/* For circular variables, compute 'wrapped' distance too, and keep the smallest of the two */
			if ( x < x0) wrapped = (x0-x-xMax)*(x0-x-xMax);
			else wrapped = (x0-x+xMax)*(x0-x+xMax);
			if ( wrapped < dx2 ) dx2 = wrapped;
	}
	/* Compute (squared) distance along y */
	dy2 = (y0-y)*(y0-y);
	if ( circY )
	{
		/* For circular variables, compute 'wrapped' distance too, and keep the smallest of the two */
			if ( y < y0) wrapped = (y0-y-yMax)*(y0-y-yMax);
			else wrapped = (y0-y+yMax)*(y0-y+yMax);
			if ( wrapped < dy2 ) dy2 = wrapped;
	}
	#ifdef DEBUG
	mexPrintf(" d2 = %f",dx2+dy2);
	#endif

	/* Test (squared) distance */
	if ( dx2+dy2 > radius2 )
	{
		#ifdef DEBUG
		mexPrintf(" [far]\n");
		#endif
		return;
	}

	/* Count this bin */
	sCountMap[i0] += countMap[i];
	sTimeMap[i0] += timeMap[i];
	#ifdef DEBUG
	mexPrintf(", count = %f, time = %f\n",sCountMap[i0],sTimeMap[i0]);
	#endif

	/* Process neighbouring bins */
	ProcessBin(x+1,y);
	ProcessBin(x-1,y);
	ProcessBin(x,y+1);
	ProcessBin(x,y-1);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	/* AdaptiveSmooth(countMap,timeMap,dt,alpha,circX,circY)) */
	int      i;

	/* Check number of parameters */
	if ( nrhs != 6 )
		mexErrMsgTxt("AdaptiveSmooth: wrong number of input arguments.");
	if ( nlhs != 3 )
		mexErrMsgTxt("AdaptiveSmooth: wrong number of output arguments.");

	/* Get input parameters */
	countMap = mxGetPr(prhs[0]);
	yMax = mxGetM(prhs[0]);
	xMax = mxGetN(prhs[0]);
	if ( xMax > yMax ) max = xMax; else max = yMax;
	timeMap = mxGetPr(prhs[1]);
	if ( yMax != mxGetM(prhs[1]) || xMax != mxGetN(prhs[1]) )
		mexErrMsgTxt("AdaptiveSmooth: count and time maps have different sizes.");	
	dt = *mxGetPr(prhs[2]);
	dt2 = dt*dt;
	alpha = *mxGetPr(prhs[3]);
	alpha2 = alpha*alpha;
	circX = *mxGetLogicals(prhs[4]);
	circY = *mxGetLogicals(prhs[5]);

	/* Create output matrices */
	plhs[0] = mxCreateDoubleMatrix(yMax,xMax,mxREAL);
	sCountMap = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(yMax,xMax,mxREAL);
	sTimeMap = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(yMax,xMax,mxREAL);
	radiusMap = mxGetPr(plhs[2]);

	/* Create a boolean matrix to keep track of visited neighbouring bins */
	visited = (bool*) mxCalloc((xMax+1)*(yMax+1),sizeof(bool));

	/* Loop through bins */
	for ( x0 = 1 ; x0 <= xMax ; ++x0 )
		for ( y0 = 1 ; y0 <= yMax ; ++y0 )
		{
			#ifdef DEBUG
			mexPrintf("------------------\nStarting @ (%d,%d)\n",x0,y0);
			#endif
			/* Reset radius */
			radius = 0;
			i0 = (int) (x0-1)*yMax+(y0-1);
			/* Increase radius until criterion is matched */
			do
			{
				/* Reset */
				for ( i = 0;i < xMax*yMax ; ++i ) visited[i] = false;
				sCountMap[i0] = 0;
				sTimeMap[i0] = 0;
				
				radius++;
				#ifdef DEBUG
				mexPrintf(" radius %f\n",radius);
				#endif
				radius2 = radius*radius;
				ProcessBin(x0,y0);
				#ifdef DEBUG
				mexPrintf(" s = %f\n",sCountMap[i0]);
				#endif
			} while ( sCountMap[i0] <= dt2*alpha/(sTimeMap[i0]*sTimeMap[i0]*radius2) );

			radiusMap[i0] = radius;
		}
}
