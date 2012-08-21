/*
 *    ECE471 - Pattern Recognition
 *
 *    Undergrad: Leandro Henrique Lourenco da Silva
 *    ID: 000356622
 *    Date: 01/25/12
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"
#include "Pr.h"
#include <vector>

using namespace std;

#define Usage "Usage: ./project1 trainingData testData\n"

#define PI 3.14159
#define LAMBDA11 0
#define LAMBDA12 1
#define LAMBDA21 1
#define LAMBDA22 0 

int main(int argc, char **argv)
{

	// check to see if the number of argument is correct
	if (argc < 3) {
		cout << Usage;
		exit(1);
	}

	//read data and fill matrix
	FILE* rawData = fopen(argv[1], "r");
	float xs, ys;
	int yc;
	vector<float> xsClass1, ysClass1, xsClass2, ysClass2;

	int class1Instances = 0, class2Instances = 0;

	while (!feof(rawData)) {
		fscanf(rawData, "%f %f %d\n", &xs, &ys, &yc);
		if (yc == 0) {
			xsClass1.insert(xsClass1.end(), xs);
			ysClass1.insert(ysClass1.end(), ys);
			class1Instances++;
		} else {
			xsClass2.insert(xsClass2.end(), xs);
			ysClass2.insert(ysClass2.end(), ys);
			class2Instances++;
		}
	}

	double probClass1 = ((double)class1Instances/(double)(class1Instances+class2Instances));
	double probClass2 = ((double)class2Instances/(double)(class1Instances+class2Instances));

	//filling vector with data
	Matrix matClass1(xsClass1.size(), 2);
	for (int i=0; i<matClass1.getRow(); i++) matClass1(i, 0) = xsClass1[i];
	for (int i=0; i<matClass1.getRow(); i++) matClass1(i, 1) = ysClass1[i];
	Matrix matClass2(xsClass2.size(), 2);
	for (int i=0; i<matClass2.getRow(); i++) matClass2(i, 0) = xsClass2[i];
	for (int i=0; i<matClass2.getRow(); i++) matClass2(i, 1) = ysClass2[i];

	//cout << "matClass1" << endl << matClass1 << endl;
	//cout << "matClass2" << endl << matClass2 << endl;


	//vector mean1 represents class1

	Matrix mean1 = mean(matClass1, 2);
	Matrix mean2 = mean(matClass2, 2);
	cout << "Mean 1:" << endl << mean1 << endl;
	cout << "Mean 2:" << endl << mean2 << endl;


	//to calculate discrimant function type1, we assume that features are independent
	double stdDeviation = 0;

	for (int i = 0; i < matClass1.getRow(); i++) {
		Matrix line = subMatrix(matClass1, i, 0, i, 1);
		stdDeviation += pow((line(0, 0) - mean1(0, 0)), 2);
	}

	stdDeviation /= matClass1.getRow();

	cout << "Standard Deviation feature 1 = " << stdDeviation << endl;


	//calculating SIGMA matrix (covariance Matrix)

	//Matrix sigmaClass1 = cov(matClass1, 2);
	Matrix sigmaClass2;
	Matrix sigmaClass1(2,2);
	sigmaClass1(0,0) = stdDeviation*stdDeviation;
	sigmaClass1(1,1) = sigmaClass1(0,0);
	sigmaClass1(0,1) = 0;
	sigmaClass1(1,0) = 0;

	Matrix invClass1 = inverse(sigmaClass1);
	Matrix invClass2 = inverse(sigmaClass1);
	
	
	//Now test if test set fits our decision rule

	Matrix testSet = readData(argv[2], 3);

	//calculating likelihood ratio first

	double detClass1 = det(sigmaClass1);
	double detClass2 = detClass1;


	//Im going to create a new matrix to keep predicted class and error

	int testCase = 0;
	while (testCase < 3) {

		float bestLambda1 = 1.0;
		float bestLambda2 = 1.0;
		double bestAccuracy = 0;

		for (float lambdaClass1 = 1.0; lambdaClass1 > 0; lambdaClass1 -= 0.05) {
			bestLambda1 = 1.0;
			bestLambda2 = 1.0;
			bestAccuracy = 0;
			for (float lambdaClass2 = 1.0; lambdaClass2 > 0; lambdaClass2 -= 0.05) {


				Matrix predictedMatrix(testSet.getRow(), 2);

				int correctGuesses = 0;

				for (int i = 0; i < testSet.getRow(); i++) {
					Matrix line(1,2); 
					line = subMatrix(testSet, i, 0, i, 1);
					line = transpose(line);
					
					Matrix mahalanobis = ((line - mean1));
					mahalanobis = transpose(mahalanobis);
					mahalanobis = mahalanobis*invClass1; 
					mahalanobis = mahalanobis*(line - mean1);
					double varMahalanobis = (-0.5)*mahalanobis(0, 0);


					double probIsClass1 = (1.0/(2*PI*detClass1))*exp(varMahalanobis)*(probClass1);


					mahalanobis = (line - mean2);
					mahalanobis = transpose(mahalanobis);
					mahalanobis = mahalanobis*invClass2;
					mahalanobis = mahalanobis*(line - mean2);
					varMahalanobis = (-0.5)*mahalanobis(0, 0);

					double probIsClass2 = (1.0/(2*PI*detClass2))*exp(varMahalanobis)*(probClass2);

					double aux = probIsClass1;
					probIsClass1 = lambdaClass1*probIsClass1 + (1 - lambdaClass1)*probIsClass2;
					probIsClass2 = lambdaClass2*probIsClass2 + (1 - lambdaClass2)*aux;

					int predictedClass = 0;
					if (probIsClass2 > probIsClass1) predictedClass = 1;

					if (predictedClass == testSet(i, 2)) correctGuesses++;

					double error = min(probClass2, probClass1);

					predictedMatrix(i, 1) = error;
					predictedMatrix(i, 0) = predictedClass;

				}

				double acc = ((double)correctGuesses/(double)testSet.getRow());
				if (acc > bestAccuracy) {
					bestAccuracy = acc;
					bestLambda1 = lambdaClass1;
					bestLambda2 = lambdaClass2;
				}

			}

		}

		cout << "detClass1 = " << detClass1 << endl;
		cout << "invClass1 = " << endl << invClass1 << endl;
		cout << "probClass1= " << probClass1 << endl;
		cout << "detClass2 = " << detClass2 << endl;
		cout << "invClass2 = " << endl << invClass2 << endl;
		cout << "probClass2= " << probClass2 << endl;

		if (testCase == 0) {
			cout << "Discriminant Function Case 1 (features independents) - accuracy = " 
			<< bestAccuracy << " Lambda1,1=" << bestLambda1 << ";Lambda1,2=" << (1 - bestLambda1) << 
			";Lambda2,1=" << bestLambda2 << ";Lambda2,2=" << (1 - bestLambda2) << endl;
		
			//calculating new covariance Matrix (case 2: covariance matrix is equal)
			sigmaClass1 = cov(matClass1, 2);
			sigmaClass2 = sigmaClass1;
			detClass1 = det(sigmaClass1);
			detClass2 = detClass1;
			invClass1 = inverse(sigmaClass1);
			invClass2 = invClass1;

		} else if (testCase == 1) {	
			cout << "Discriminant Function Case 2 (covariance matrix equal for all features - accuracy = "
			<< bestAccuracy << " Lambda1,1=" << bestLambda1 << ";Lambda1,2=" << (1 - bestLambda1) << 
			";Lambda2,1=" << bestLambda2 << ";Lambda2,2=" << (1 - bestLambda2) << endl;
			//calculating new covariance matrices (case 3: covariance matrix are different)
			//sigma1 and invClass1 are already calculated
			sigmaClass2 = cov(matClass2, 2);
			invClass2 = inverse(sigmaClass2);
			detClass2 = det(sigmaClass2);

		} else if (testCase == 2) {
			cout << "Discriminant Function Case 3 (arbitrary) - accuracy = "
				<< bestAccuracy << " Lambda1,1=" << bestLambda1 << ";Lambda1,2=" << (1 - bestLambda1) << 
			";Lambda2,1=" << bestLambda2 << ";Lambda2,2=" << (1 - bestLambda2) << endl;		
		}



		testCase++;
	}

	cout << "Process completed" << endl;

	//writeData(predictedMatrix, "predictedValues.txt");


	return 0;
}
