#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
void print_matrix45(double matrix[4][5])
{
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<5; j++) std::cout << matrix[i][j] << '\t';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_matrix56(double matrix[5][6])
{
    for(int i=0; i<5; i++)
    {
        for(int j=0; j<6; j++) std::cout << matrix[i][j] << '\t';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int find_biggest_row45(double matrix[4][5])
{
    int max = 0;
    for(int i=0; i<5; i++)
    {
        if(matrix[max][0] < matrix[i][0]) max = i;
    }
    return max;
}

int find_biggest_row56(double matrix[5][6])
{
    int max = 0;
    for(int i=0; i<6; i++)
    {
        if(matrix[max][0] < matrix[i][0]) max = i;
    }
    return max;
}

bool zad1()
{
	const double eps = 1e-12; // sta≥a - przybliøenie zera

	const size_t N = 4;
	double A[N][N] = { {4,-2,4,-2},{3,1,4,2},{2,4,2,1},{2,-2,4,2} }; // przyk≥adowe dane
	double b[N] = { 8,7,10,2 }; // przyk≥adowe dane
	double x[N];
	double Ab[N][N + 1];
	//tworzenie rozszerzonej macierzy
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N + 1; j++)
			if (j < N)
			{
				Ab[i][j] = A[i][j];
			}
			else Ab[i][j] = b[i];
	}
	//kopiowanie tablicy
	double Ab_copy[N][N + 1] = {};
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			Ab_copy[i][j] = Ab[i][j];
			x[j] = Ab[i][N]; /// nie wiem o co tu chodzi ale dzia≥a, metoda prÛb i b≥ÍdÛw
		}
	}
	//printowanie macierzy
    print_matrix45(Ab_copy);
    int greatest = find_biggest_row45(Ab_copy);
    for(int i = 0; i < 5; i++) std::swap(Ab_copy[greatest][i], Ab_copy[0][i]);
    print_matrix45(Ab_copy);
	double m = 0; // mnoønik
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N ; j++)
		{
			if (fabs(Ab_copy[i][i]) < eps) return false; //fabs zwraca wartosc absolutna
			m = Ab_copy[j][i] / Ab_copy[i][i];
			for (int k = i + 1; k <= N; k++)
			{
				Ab_copy[j][k] = Ab_copy[j][k] - m * Ab_copy[i][k];
			}
		}
	}
	//etap wyliczania zmiennych
	double s = 0; // zlicza sumÍ iloczynÛw
	for (int i = N-1; i >= 0; i--)
	{
		s = Ab_copy[i][N]; // do zmiennej s trafia wektor b (ostatnia kolumna macierzy)
		for (int j = N; j >= i + 1; j--)
		{
			s = s - Ab_copy[i][j] * x[j];
		}
		if (fabs(Ab_copy[i][i]) < eps) return false;
		x[i] = s / Ab_copy[i][i];
		std::cout << 'x' << i << "= " << x[i] << std::endl;
	}
    //find_biggest(Ab_copy);
	return true;

}

bool zadanie2()
{
    const double eps = 1e-12; // sta≥a - przybliøenie zera
    const int N = 5;
    double A[N][N] = {{0,0,2,1,2},{0,1,0,2,-1},{1,2,0,-2,0},{0,0,0,-1,1},{0,1,-1,1,-1}};
    double b[5] = {1,1,-4,-2,-1};
    double x[N];
    double A2[N][N + 1];
    //tworzenie rozszerzonej macierzy
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
            if (j < N)
            {
                A2[i][j] = A[i][j];
            }
            else A2[i][j] = b[i];
    }
    double A22[N][N + 1] = {};
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            A22[i][j] = A2[i][j];
            x[j] = A22[i][N]; /// nie wiem o co tu chodzi ale dzia≥a, metoda prÛb i b≥ÍdÛw
        }
    }

    print_matrix56(A22);
    int greatest = find_biggest_row56(A22);
    for(int i = 0; i < 6; i++) std::swap(A22[greatest][i], A22[0][i]);
    print_matrix56(A22);
    double m = 0; // mnoønik
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N ; j++)
        {
            if (fabs(A22[i][i]) < eps) return false; //fabs zwraca wartosc absolutna
            m = A22[j][i] / A22[i][i];
            for (int k = i + 1; k <= N; k++)
            {
                A22[j][k] = A22[j][k] - m * A22[i][k];
            }
        }
    }
    //etap wyliczania zmiennych
    double s = 0; // zlicza sumÍ iloczynÛw
    for (int i = N-1; i >= 0; i--)
    {
        s = A22[i][N]; // do zmiennej s trafia wektor b (ostatnia kolumna macierzy)
        for (int j = N; j >= i + 1; j--)
        {
            s -= A22[i][j] * x[j];
        }
        if (fabs(A22[i][i]) < eps) return false;
        x[i] = s / A22[i][i];
        std::cout << 'x' << i << " = " << x[i] << std::endl;
    }
    return true;
}

// ZADANIE 3

// obliczone wartosci sa dwukrotnie wieksze niz te sugerowane przez internetowe kalkulatory
void zadanie3()
{
    const double eps = 1e-12;
    const int N = 5;
    double A[N][N] = {{1,0,0,0,0},{1,1,1,1,1},{1,3,9,27,81},{1,5,25,125,625},{1,6,36,216,1296}};
    double b[N] = {1,1,-4,-2,-1};
    double A2[N][N+1];
    double x[N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
            if (j < N) A2[i][j] = A[i][j];
            else A2[i][j] = b[i];
    }
    double A22[N][N + 1] = {};
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            A22[i][j] = A2[i][j];
            x[j] = A22[i][N]; /// nie wiem o co tu chodzi ale dzia≥a, metoda prÛb i b≥ÍdÛw
        }
    }
    int greatest = find_biggest_row56(A22);
    for(int i = 0; i < 6; i++) std::swap(A22[greatest][i], A22[0][i]);
    double m = 0; // mnoønik
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N ; j++)
        {
            if (fabs(A22[i][i]) < eps) std::cout << "WRONG 2"; //fabs zwraca wartosc absolutna
            m = A22[j][i] / A22[i][i];
            for (int k = i + 1; k <= N; k++)
            {
                A22[j][k] = A22[j][k] - m * A22[i][k];
            }
        }
    }
    //etap wyliczania zmiennych
    double s = 0; // zlicza sumÍ iloczynÛw
    for (int i = N-1; i >= 0; i--)
    {
        s = A22[i][N]; // do zmiennej s trafia wektor b (ostatnia kolumna macierzy)
        for (int j = N; j >= i + 1; j--)
        {
            s -= A22[i][j] * x[j];
        }
        if (fabs(A22[i][i]) < eps) std::cout << "WRONG 1";
        x[i] = s / A22[i][i];
        std::cout << 'x' << i << " = " << x[i] << std::endl;
    }

}

double f1(double x)
{
    return tan(x) - x;
}

double df1(double x)
{
    return (1/pow(cos(x),2)) - 1;
}

int metoda_newtona()
{
    double prev = 4.5;
    double next;
    double eps = 0.01;
    int counter = 0;
    bool flag = true;
    while(flag == true)
    {
        double h = f1(prev)/df1(prev);
        next = prev - h;
        //((tan(prev)-prev)/(1/(pow(cos(prev),2))-1));
        prev = next;
        if(eps <= (next - prev)) flag = false;
        //std::cout << std::setprecision(14) << next << '\t' << std::setprecision(14) << prev << '\t' << std::setprecision(14) << next-prev << '\t' << counter << std::endl;
        counter++;
    }
    return counter;
}

double f2(double x)
{
    return tan(x) - x;
}

int metoda_bisekcji()
{
    double eps = 0.001;
    double a = 4.4;
    double b = 4.6;
    double c;
    int counter = 0;
    while((b-a) >= eps)
    {
        c = (a + b)/2;
        if(f2(c) == 0) break;
        else if(f2(c) * f2(a) < 0) b = c;
        else a = c;
        // std::cout << c << std::endl;
        counter++;
    }
    return counter;
}

double f3(double x)
{
    return tan(x) - x;
}

int metoda_siecznych()
{
    double eps = 0.01;
    int counter = 0;
    double prev = 4.5000001;
    double prev2 = 0;
    double next = 0;
    while((next - prev <= eps))
    {
        next = prev - f3(prev)*(prev - prev2)/(f3(prev) - f3(prev2));
        counter++;
    }
    return counter;
}

void metoda_brenta()
{

}

void zadanie4()
{
    int c;
    std::cout << "czesc! obliczam x = tan(x) w okolicy 4.5" << std::endl;
    std::cout << "1. metoda newtona; \n" << "2. metoda bisekcji;\n" << "3. metoda siecznych;\n" << "4. metoda brenta\n";
    std::cin >> c;
    switch (c) {
        case 1:
            std::cout << metoda_newtona() << std::endl;;
            break;
        case 2:
            std::cout << metoda_bisekcji() << std::endl;
            break;
        case 3:
            std::cout << metoda_siecznych() << std::endl;
            break;
        case 4:
            metoda_brenta();
        default:
            std::cout << "zly wybor" << std::endl;
            break;
    }
}

int main(int argc, const char * argv[])
{
    zadanie4();
	return 0;
}


/*
 double A[4][4] = { {1,-2,4,-2},{3,1,4,2},{2,4,2,1},{2,-2,4,2} };
 print_matrix(A);
 int row = find_biggest_row(A);
 for(int i=0; i<4; i++) std::swap(A[0][i], A[row][i]);
 print_matrix(A);
 */
