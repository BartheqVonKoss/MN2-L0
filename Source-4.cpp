#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

void print_matrix45(double matrix[4][5]);
void print_matrix56(double matrix[5][6]);
int find_biggest_row45(double matrix[4][5]);
int find_biggest_row56(double matrix[5][6]);
void Ex1();
void Ex2();
void Ex3();
double f1(double &x);                           // f(x) = tan(x) - x
double df1(double &x);                          // df(x)
void metoda_newtona(double &eps);               // newton's method
double df1(double &x);
double f1(double &x);
void metoda_bisekcji(double &eps);              // bisection method
void metoda_siecznych(double &eps);
void metoda_brenta(double &eps);                // brent's method
long double dw(long double t, long double theta, long double w);
// Ex. 1

void print_matrix45(double matrix[4][5])        // function to print 4x5 matrix
{
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<5; j++) std::cout << matrix[i][j] << '\t';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_matrix56(double matrix[5][6])        // function to print 5x6 matrxi
{
    for(int i=0; i<5; i++)
    {
        for(int j=0; j<6; j++) std::cout << matrix[i][j] << '\t';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int find_biggest_row45(double matrix[4][5])     // function to find biggest element in rows 4x5
{
    int max = 0;
    for(int i=0; i<5; i++)
    {
        if(matrix[max][0] < matrix[i][0]) max = i;
    }
    return max;
}

int find_biggest_row56(double matrix[5][6])     // function to find biggest element in rows 5x6
{
    int max = 0;
    for(int i=0; i<6; i++)
    {
        if(matrix[max][0] < matrix[i][0]) max = i;
    }
    return max;
}

void Ex1()                                      // gaussian elimination
{
	const double eps = 1e-12;                   // accuracy of comparison
	const int N = 4;
	double A[N][N] = { {2, 1, -1, 2} , { 4, 5, -3, 6} , { -2, 5, -2, 6} , { 4, 11, -4, 8} }; // 4x4 matrix ready to be gausiannly eliminated
	double b[N] = { 5, 9, 4, 2 };              // b (1x4)
	double x[N];
	double A2[N][N + 1];
    for (int i = 0; i < N; i++)                 // extension of the A matrix by a column of b [4x4+4x1]
	{
		for (int j = 0; j < N + 1; j++)
			if (j < N)
			{
				A2[i][j] = A[i][j];
			}
			else A2[i][j] = b[i];
	}
	double A22[N][N + 1] = {};
	for (int i = 0; i < N; i++)                 // duplication
	{
		for (int j = 0; j < N + 1; j++)
		{
			A22[i][j] = A2[i][j];
			x[j] = A2[i][N];
		}
	}
    // print_matrix45(A22);
    int greatest = find_biggest_row45(A22);
    for(int i = 0; i < 5; i++) std::swap(A22[greatest][i], A22[0][i]); // search for a greatest element in row
	double m = 0;                               // coefficient
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N ; j++)
		{
            if (fabs(A2[i][i]) < eps) std::cout << "WRONG09" << std::endl;
			m = A22[j][i] / A22[i][i];
			for (int k = i + 1; k <= N; k++)
			{
				A22[j][k] = A22[j][k] - m * A22[i][k];
			}
		}
	}
	double s = 0;
    for (int i = N-1; i >= 0; i--)
	{
		s = A22[i][N];
		for (int j = N; j >= i + 1; j--)            // x elimination
		{
			s = s - A22[i][j] * x[j];
		}
        if (fabs(A22[i][i]) < eps) std::cout << "WRONG10" << std::endl;
		x[i] = s / A22[i][i];
        std::cout << " x" << i << "= " << x[i];
	}
    std::cout << std::endl;
    //find_biggest(A22);
}

// Ex. 2

void Ex2()                                      // works identically to Ex1()
{
    const double eps = 1e-12;
    const int N = 5;
    double A[N][N] = {{0,0,2,1,2},{0,1,0,2,-1},{1,2,0,-2,0},{0,0,0,-1,1},{0,1,-1,1,-1}};
    double b[5] = {1,1,-4,-2,-1};
    double x[N];
    double A2[N][N + 1];
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
            x[j] = A22[i][N];
        }
    }
    int greatest = find_biggest_row56(A22);
    for(int i = 0; i < 6; i++) std::swap(A22[greatest][i], A22[0][i]);
    double m = 0;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N ; j++)
        {
            if (fabs(A22[i][i]) < eps) std::cout << "WRONG1" << std::endl;
            m = A22[j][i] / A22[i][i];
            for (int k = i + 1; k <= N; k++)
            {
                A22[j][k] = A22[j][k] - m * A22[i][k];
            }
        }
    }
    double s = 0;
    for (int i = N-1; i >= 0; i--)
    {
        s = A22[i][N];
        for (int j = N; j >= i + 1; j--)
        {
            s -= A22[i][j] * x[j];
        }
        if (fabs(A22[i][i]) < eps) std::cout << "WRONG1" << std::endl;
        x[i] = s / A22[i][i];
        std::cout << " x" << i << " = " << x[i];
    }
    std::cout << std::endl;
}

// Ex. 3

// obliczone wartosci sa dwukrotnie wieksze niz te sugerowane przez internetowe kalkulatory
void Ex3()                                      // calculate coefficients of polynomial that goes through some points
{
    const double eps = 1e-12;
    const int N = 5;
    double A[N][N] = {{1,0,0,0,0},{1,1,1,1,1},{1,3,9,27,81},{1,5,25,125,625},{1,6,36,216,1296}};    // polynomial
    double b[N] = {1,1,-4,-2,-1};               // "some points"
    double A2[N][N+1];                          // initial matrix "polynomial" extension
    double x[N];
    for (int i = 0; i < N; i++)                 // extending
    {
        for (int j = 0; j < N + 1; j++)
            if (j < N) A2[i][j] = A[i][j];
            else A2[i][j] = b[i];
    }
    double A22[N][N + 1] = {};
    for (int i = 0; i < N; i++)                 // building the new matrix & vector x
    {
        for (int j = 0; j < N + 1; j++)
        {
            A22[i][j] = A2[i][j];
            x[j] = A22[i][N];
        }
    }
    int greatest = find_biggest_row56(A22);     // look for the greatest row to avoid zero division
    for(int i = 0; i < 6; i++) std::swap(A22[greatest][i], A22[0][i]);
    double m = 0;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N ; j++)        // matrix traingulation
        {
            if (fabs(A22[i][i]) < eps) std::cout << "WRONG 2";
            m = A22[j][i] / A22[i][i];
            for (int k = i + 1; k <= N; k++)
            {
                A22[j][k] = A22[j][k] - m * A22[i][k];
            }
        }
    }
    double s = 0;
    for (int i = N-1; i >= 0; i--)
    {
        s = A22[i][N];
        for (int j = N; j >= i + 1; j--)
        {
            s -= A22[i][j] * x[j];
        }
        if (fabs(A22[i][i]) < eps) std::cout << "WRONG 1";
        x[i] = s / A22[i][i];
        std::cout << " x" << i << " = " << x[i];
    }
    std::cout << std::endl;

}

// Ex. 4


double f1(double &x)                            // calculates f(x) = tanx - x
{
    return tan(x) - x;
}

double df1(double &x)                           // calculates df(x)/dx
{
    return (1/(pow(cos(x),2))) - 1;
}

void metoda_newtona(double &eps)                // newton's method
{
    double next = 4.5;                          // out of blue
    double prev;
    int counter = 1;
    do
    {
        prev = next;
        next = prev - f1(prev)/df1(prev);       // calculation of next x value
        std::cout << prev << '\t' << next << '\t' << fabs(next - prev) << '\t' << counter << std::endl;
        counter++;
    }while(fabs(next - prev) > eps);            // convergence condition
}

void metoda_bisekcji(double &eps)               // bisection method
{
    double a = 4.4;                             // declaring a & b that contain the solution
    double b = 4.6;
    double c;
    double e = b - a;
    int counter = 1;
    while(e >= eps)                             // convergence condition
    {
        e = e/2;
        c = a + e;
        if(f1(c) == 0) break;
        else if(f1(c) * f1(a) < 0) b = c;
        else a = c;
        std::cout << c << '\t' << e << '\t' << counter << std::endl;
        counter++;
    }
}

void metoda_siecznych(double &eps)
{
    double prev = 4.5;                          // the method needs two starting points
    double even_more_prev = 4.501;
    double present;
    int counter = 1;
    do
    {
        present = prev - f1(prev) * (prev - even_more_prev)/(f1(prev) - f1(even_more_prev));
        std::cout << even_more_prev << '\t' << prev << '\t' << present << '\t' << fabs(prev - even_more_prev) << '\t' << counter << std::endl;
        even_more_prev = prev;
        prev = present;
        counter++;
    }while(!((fabs(prev - even_more_prev) < eps) || (fabs(f1(present)) < eps)));    // convergence condition
}

void metoda_brenta(double &eps)                 // brent's method
{
    double a = 4.4;
    double b = 4.6;
    double s, c, f_s, f_a, f_b, f_c;
    f_a = f1(a);
    f_b = f1(b);
    int counter = 1;
    if(a>b || f_a * f_b >= 0)                   // algorithm assumes a<b
    {
        std::cout << "zle\n";                   // root not bound by the given guesses
    }

    do{
        c = (a+b)/2;
        f_c = f1(c);
        if(fabs(f_a-f_c) > eps && fabs(f_b-f_c) > eps)     // f(a)!=f(c) and f(b)!=f(c)
        {
                                                // inverse quadratic interpolation
            s = a*f_b*f_c/((f_a-f_b)*(f_a-f_c)) + b*f_a*f_c/((f_b-f_a)*(f_b-f_c)) + c*f_a*f_b/((f_c-f_a)*(f_c-f_b));
        }
        else
        {
            s = b - f_b*(b-a)/(f_b-f_a);
        }
        f_s = f1(s);
        if(f_c*f_s < 0)
        {
            a = s;
            b = c;
        }
        else if(f_s*f_b < 0)
        {
            a = c;
        }
        else
        {
            b = s;
        }
        f_s = fabs(f_s);
        std::cout << s << '\t' << fabs(b-a) << '\t' << counter << std::endl;
        counter++;
    }while((f_s < eps) || (fabs(b-a) > eps));    // convergence conditions
}

void Ex4()
{
    double eps = 1e-12;
    int c;      // not a problem at all
    while(c!=0)
    {
    std::cout << "czesc! obliczam x = tan(x) w okolicy 4.5, z epsilon rownym " << eps << "!" << std::endl;
    std::cout << "1. metoda newtona; \n" << "2. metoda bisekcji;\n" << "3. metoda siecznych;\n" << "4. metoda brenta\n" << "0. Exit\n";
    std::cin >> c;
    switch (c) {
        case 1:
            metoda_newtona(eps);
            break;
        case 2:
            metoda_bisekcji(eps);
            break;
        case 3:
            metoda_siecznych(eps);
            break;
        case 4:
            metoda_brenta(eps);
            break;
        case 0:
            std::cout << "bye" << std::endl;
            break;
        default:
            std::cout << "zly wybor" << std::endl;
            break;
    }
    }
}

// Ex. 5

long double dw(long double t, long double theta, long double w)        // calculation of velocity
{
    long double Q = 2;
    long double A = 1.35;
    return A * cos(t) - sin(theta) - w / Q;
}

void Ex5()                                // runge - kutta 4th order method
{
    long double h = 0.0002;
    const int n = 100000;
    long double kt1, kt2, kt3, kt4, kdt1, kdt2, kdt3, kdt4;
    std::ofstream myData;
    myData.open("/Users/bartlomiejkos/Documents/Programming/C++/Metody numeryczne II/L0_2/data.csv", std::ios::out);
    std::vector<long double> theta;
    std::vector<long double> dtheta;
    std::vector<long double> t;
    theta.push_back(0.3);
    dtheta.push_back(0);
    t.push_back(0);
    for(int i=0; i<n; i++)
    {
        kt1 = dtheta[i];
        kdt1 = dw(theta[i], t[i], dtheta[i]);

        kt2 = dtheta[i] + kt1 / 2;
        kdt2 = dw(theta[i] + kdt1 / 2, t[i] + h / 2, dtheta[i]);

        kt3 = dtheta[i] + kt2 / 2;
        kdt3 = dw(theta[i] + kdt2 / 2, t[i] + h / 2, dtheta[i]);

        kt4 = dtheta[i] + kt3;
        kdt4 = dw(theta[i] + kdt3, t[i] + h, dtheta[i]);

        theta.push_back(theta[i] + h / 6 * (kt1 + 2 * kt2 + 2 * kt3 + kt4));
        dtheta.push_back(dtheta[i] + h / 6 * (kdt1 + 2 * kdt2 + 2 * kdt3 + kdt4));
        t.push_back(t[i] + h);

        std::cout << theta[i] << '\t' << dtheta[i] << '\t' << t[i] << '\t' << i << std::endl;
        myData << theta[i] << ',' << dtheta[i] << ',' << t[i] << ',' << i << '\n';
    }

}

int main(int argc, const char * argv[])
{
    Ex5();
	return 0;
}


