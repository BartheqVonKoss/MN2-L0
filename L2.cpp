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
void Ex1(int n, int M, int f);
void Ex2a(int &n, int &M, int &f);
void Ex2b(int &n, int &M, int &f);
void Ex3(int &n, int &M, int &f);
double g(double x, double y);
double g2(double x, double y);
double u2(double &x, double &y);
double u3(double &x, double &y);
double u4(double &x, double &y);
void Ex4();


double g(double x, double y)
{
    return (4 * x * y * (pow(x,2) - pow(y,2)));
}
// Ex. 1
void Ex1(int n, int M, int f)
{
    //std::ofstream myData;
    //myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Metody numeryczne II/L3/data.csv", std::ios::out);
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
            //myData << matrix[i][j] << ',';
        }
        std::cout << std::endl;
        //myData << '\n';
    }
    for(int i = 0; i < n + 1; i++) delete [] matrix[i];
    delete [] matrix;
}
// Ex. 2a
void Ex2a(int &n, int &M, int &f)
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Downloads/Programming/MetodyNumeryczneII/data.csv", std::ios::out);
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
    for(int i = 0; i < n + 1; i++) delete [] matrix[i];
    delete [] matrix;
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
    for(int i = 0; i < n + 1; i++) delete [] matrix[i];
    delete [] matrix;
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
    for(int i = 0; i < n + 1; i++) delete [] matrix[i];
    delete [] matrix;
}

// Ex. 4
double u2(double &x, double &y)
{
    return x * x + y * y;
}

double u3(double &x, double &y)
{
    return pow(x, 4) - 6 * x * x * y * y + y * y;
}

double u4(double &x, double &y)
{
    return pow(x, 6) - 15 * pow(x, 4) * y * y + 15 * x * x * pow(y, 4) - pow(y, 6);
}

void Ex4()
{
    double x1,x2,x3,x4,y1,y2,y3,y4;
    x1 = 0; x2 = 1; x3 = 1; x4 = 1; y1 = 2; y2 = 0; y3 = 1; y4 = 2;
    double matrix[4][4];
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            if(j == 0) matrix[i][j] = 1;
            else if(j == 1 && i == 0) matrix[i][j] = u2(x1, y1);
            else if(j == 1 && i == 1) matrix[i][j] = u2(x2, y2);
            else if(j == 1 && i == 2) matrix[i][j] = u2(x3, y3);
            else if(j == 1 && i == 3) matrix[i][j] = u2(x4, y4);
            else if(j == 2 && i == 0) matrix[i][j] = u3(x1, y1);
            else if(j == 2 && i == 1) matrix[i][j] = u3(x2, y2);
            else if(j == 2 && i == 2) matrix[i][j] = u3(x3, y3);
            else if(j == 2 && i == 3) matrix[i][j] = u3(x4, y4);
            else if(j == 3 && i == 0) matrix[i][j] = u4(x1, y1);
            else if(j == 3 && i == 1) matrix[i][j] = u4(x2, y2);
            else if(j == 3 && i == 2) matrix[i][j] = u4(x3, y3);
            else if(j == 3 && i == 3) matrix[i][j] = u4(x4, y4);
        }
    }
    double f[4];                        // here implement gaussian elimination but my algorithm provides with wrong results
    f[0] = u2(x1,y1); f[1] = u2(x2,y2); f[2] = u2(x3,y3); f[3] = u2(x4,y4);
}


int main(int argc, const char * argv[])
{
    // Ex.2a,2b
    int n = 100;
    int M = 100;
    int f = 1;                 // printing frequency (multiplication coefficient while printing results) NOTE: if frequency was changed so it has to be changed in python script
    //Ex2a(n, M, f);
    //Ex2b(n, M, f);           // slightly different plotted surface than expected
    // Ex.3
    M = 10000;
    //Ex3(n, M, f);
    return 0;
}
