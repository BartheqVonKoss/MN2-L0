//
//  main.cpp
//  L2
//
//  Created by Bartłomiej Kos on 30/03/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include <iostream>
#include <math.h>
#include "fstream"
void Ex1(int n, int M, int f);
void Ex2a(int &n, int &M, int &f);
void Ex2b(int &n, int &M, int &f);
void Ex3(int &n, int &M, int &f);
double g(double x, double y);
double g2(double x, double y);

double g(double x, double y)
{
    return (4 * x * y * (pow(x,2) - pow(y,2)));
}
// Ex. 1
void Ex1(int n, int M, int f)
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Metody numeryczne II/L3/data.csv", std::ios::out);
    double h;
    double k;
    double **matrix = new double *[n+1];                // definition of base to hold results
    for(int i = 0; i < n + 1; i++) matrix[i] = new double [n+1];
    for(int i = 0; i < n + 1; i++)
    {
        for(int j = 0; j < n + 1; j++)
        {
            if(i == 0) matrix[i][j] = g(0, j*k);
            else if(j == 0) matrix[i][j] = g(i*h, 0);
            else if(i == n) matrix[i][j] = g(1, j*k);
            else if(j == n) matrix[i][j] = g(i*h, 1);
            else matrix[i][j] = 0;
        }
    }

    for(int i = 0; i < n + 1; i++)
    {
        for(int j = 0; j < n + 1; j++)
        {
            std::cout << matrix[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    for(int k = 1; k < M; k++)
    {
        for(int j = 1; j < n; j++)
        {
            for(int i = 1; i < n; i++)
            {
                matrix[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1]) / 4;
            }
        }
    }
    for(int i = 0; i < n + 1; i+=f)
    {
        for(int j = 0; j < n + 1; j+=f)
        {
            std::cout << matrix[i][j] << '\t';
            myData << matrix[i][j] << ',';
        }
        std::cout << std::endl;
        myData << '\n';
    }
}
// Ex. 2a
void Ex2a(int &n, int &M, int &f)
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Downloads/Programming/MetodyNumeryczneII/data.csv", std::ios::out);
    double h = 0.001;
    double k = 0.001;
    double **matrix = new double *[n+1];
    for(int i = 0; i < n + 1; i++) matrix[i] = new double [n+1];
    for(int i = 0; i < n + 1; i++)
    {
        for(int j = 0; j < n + 1; j++)
        {
            if(i == 0) matrix[i][j] = g(0, j*k);
            else if(j == 0) matrix[i][j] = g(i*h, 0);
            else if(i == n) matrix[i][j] = g(1, j*k);
            else if(j == n) matrix[i][j] = g(i*h, 1);
            else matrix[i][j] = 0;
        }
    }
    for(int k = 1; k < M; k++)
    {
        for(int j = 1; j < n; j++)
        {
            for(int i = 1; i < n; i++)
            {
                matrix[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1]) / 4;
            }
        }
    }

    for(int i = 0; i < n + 1; i++)
    {
        i *= f;
        for(int j = 0; j < n + 1; j++)
        {
            j *= f;
            if(i < n && j < n)
            {
                std::cout << matrix[i][j] << '\t';
                myData << matrix[i][j] << ',';
            }
            j /= f;
        }
        myData << '\n';
        i /= f;
    }
}

double g2(double x, double y)
{
    return 1e-4 * sin(3 * M_PI * x) * sin(3 * M_PI * y);
}
// Ex. 2b
void Ex2b(int &n, int &M, int &f)
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Downloads/Programming/MetodyNumeryczneII/data2.csv", std::ios::out);
    double h = 0.01;
    double k = 0.01;
    double **matrix = new double *[n+1];
    for(int i = 0; i < n + 1; i++) matrix[i] = new double [n+1];
    for(int i = 0; i < n + 1; i++)
    {
        for(int j = 0; j < n + 1; j++)
        {
            if(i == 0) matrix[i][j] = g2(0, j*k);
            else if(j == 0) matrix[i][j] = g2(i*h, 0);
            else if(i == n) matrix[i][j] = g2(1, j*k);
            else if(j == n) matrix[i][j] = g2(i*h, 1);
            else matrix[i][j] = 0;
        }
    }
    for(int k = 1; k < M; k++)
    {
        for(int j = 1; j < n; j++)
        {
            for(int i = 1; i < n; i++)
            {
                matrix[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1]) / 4;
            }
        }
    }
    for(int i = 0; i < n + 1; i++)
    {
        i *= f;
        for(int j = 0; j < n + 1; j++)
        {
            j *= f;
            if(i < n && j < n)
            {
                std::cout << matrix[i][j] << '\t';
                myData << matrix[i][j] << ',';
            }
        }
        std::cout << std::endl;
        myData << '\n';
    }
}

// Ex. 3
void Ex3(int &n, int &M, int &f)
{
    double eps = 1e-10;
    double h = 0.01;
    double k = 0.01;
    double **matrix = new double *[n+1];
    for(int i = 0; i < n + 1; i++) matrix[i] = new double [n+1];
    for(int i = 0; i < n + 1; i++)
    {
        for(int j = 0; j < n + 1; j++)
        {
            if(i == 0) matrix[i][j] = g(0, j*k);
            else if(j == 0) matrix[i][j] = g(i*h, 0);
            else if(i == n) matrix[i][j] = g(1, j*k);
            else if(j == n) matrix[i][j] = g(i*h, 1);
            else matrix[i][j] = 0;
        }
    }
    double temp;
    double w = 0;
    for(int k = 1; k < M; k++)
    {
        int counter = 0;
        for(int j = 1; j < n; j++)
        {
            for(int i = 1; i < n; i++)
            {
                temp = matrix[i][j];
                matrix[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1]) / 4;
                if(fabs(matrix[i][j] - temp) < eps) counter++;
            }
        }
        w++;
        if(counter > 9800) break;
    }

    for(int i = 0; i < n + 1; i++)
    {
        i *= f;
        for(int j = 0; j < n + 1; j++)
        {
            j *= f;
            if(i < n && j < n)
            {
                std::cout << matrix[i][j] << '\t';
            }
        }
    }
    std::cout << std::endl;
    std::cout << w << std::endl;
}

// Ex. 4

int main(int argc, const char * argv[])
{
    /*
    // Ad. 1,2a,2b,3
    int n = 100;
    int M = 100000;
    int f = 1;                 // printing frequency (multiplication coefficient while printing results)
     */
    
     return 0;
}
