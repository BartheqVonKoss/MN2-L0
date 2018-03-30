//
//  main.cpp
//  L3
//
//  Created by Bartłomiej Kos on 30/03/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include <iostream>
#include <math.h>
#include "fstream"
void Ex1();

double g(double x, double y)
{
    return (4 * x * y * (pow(x,2) - pow(y,2)));
}

void Ex1()
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Metody numeryczne II/L3/data.csv", std::ios::out);
    int n = 1000;
    int M = 1000;
    double h = 0.001;
    double k = 0.001;
    double **matrix = new double *[n+1];
    for(int i = 0; i < n+1; i++) matrix[i] = new double [n+1];
    for(int i = 0; i < n+1; i++)
    {
        for(int j=0; j < n+1; j++)
        {
            if(i == 0) matrix[i][j] = g(0, j*k);
            else if(j == 0) matrix[i][j] = g(i*h, 0);
            else if(i == n) matrix[i][j] = g(1, j*k);
            else if(j == n) matrix[i][j] = g(i*h, 1);
            else matrix[i][j] = 0;
        }
    }

    for(int i = 0; i < n+1; i++)
    {
        for(int j=0; j< n+1; j++)
        {
            std::cout << matrix[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    for(int k=1; k<M; k++)
    {
        for(int j=1; j<n; j++)
        {
            for(int i=1; i<n; i++)
            {
                matrix[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1]) / 4;
            }
        }
    }
    for(int i = 0; i < n+1; i++)
    {
        for(int j=0; j< n+1; j++)
        {
            std::cout << matrix[i][j] << '\t';
            myData << matrix[i][j] << ',';
        }
        std::cout << std::endl;
        myData << '\n';
    }
}

int main(int argc, const char * argv[])
{
    Ex1();
    return 0;
}
