//
//  req.cpp
//  List3
//
//  Created by Bartłomiej Kos on 14/05/2018.
//  Copyright © 2018 Bartłomiej Kos. All rights reserved.
//

#include "req.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

double f(double x)
{
    return x * x;
}

double fp(double x, double y)
{
    double h = 0.02;
    return (x - y) / h;
}

double g(double x)
{
    return sin(x);
}

void charac()
{
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Metody numeryczne II/List3/data.csv", std::ios::out);
    double n = 100;
    double h = 1 / (n - 1);
    std::vector<double> x;
    std::vector<double> u;
    std::vector<double> q;
    std::vector<double> p;
    x.resize(100);
    u.resize(100);
    q.resize(100);
    p.resize(100);

    for(int j = 1; j < n; j++)
    {
        x.at(j) = (j - 1) * h;
        u.at(j) = f(x[j]);
        p.at(j) = fp(x[j], x[j-1]);
        q.at(j) = g(x[j]);
    }
    double ph, qh;
    for(int i = 1; i < n; i++)
    {
        for(int j = 1; j < n - i; j++)
        {
            x.at(j) = x[j] + h / 2;
            ph = (p[j] + p[j+1]) / 2 + (1 - h / 8) *  (q[j+1] - q[j]);
            qh = (p[j+1] - p[j] + (2 - h / 4) * (q[j+1] + q[j])) / (4 + h / 2);
            u.at(j) = u[j] + h * (ph + p[j]) / 4 + h * (qh + q[j]) / 2;
            p.at(j) = ph;
            q.at(j) = qh;
        }
    }
    for(int i=0; i<n;i++)
    {
        std::cout << x[i] << '\t' << u[i] << '\t' << p[i] << '\t' << q[i] << '\t' << i * h << std::endl;
        myData << x[i] << ',' << u[i] << ',' << p[i] << ',' << q[i] << ',' << i * h << '\n';
    }


}
