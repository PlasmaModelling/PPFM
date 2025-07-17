#include "Tools.h"

#include"alglib/source/interpolation.h"
#include <stdexcept>
#include <sstream>
/* 
    #include <cmath>
    #include <algorithm> 
    #include <tuple>
    #include <typeinfo>
    #include <omp.h>
    #include <sstream>
    #include <initializer_list>
*/

/* Kronecker delta */
int delta(int i, int j)
{
	if (i == j)
		return 1;
    else
    	return 0;
}

// N! 
double factorial(int n) {
    double result = 1.0;
    for(int i = 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}


double **matrix_alloc(int N, int M)
{
	int i;
	double **matrix;
	matrix = (double **)calloc(N, sizeof(double *));

	for (i = 0; i < N; i++)
		matrix[i] = (double *)calloc(M, sizeof(double));

	return matrix;
}

void lu_ssolve(double *x, int *p, double *b0, double **m, int N)
{
	int k, i, j, l;

	double tmp;
	double sum;

	double *b = new double[N];

	for (i = 0; i < N; i++)
		b[i] = b0[i];

	for (k = 0; k < N - 1; k++)
	{
		j = p[k];

		if (j != k)
		{
			tmp = b[j];
			b[j] = b[k];
			b[k] = tmp;
		}
		for (i = k + 1; i < N; i++)
			b[i] += m[i][k] * b[k];
	}

	x[N - 1] = b[N - 1] / m[N - 1][N - 1];

	for (i = N - 2; i > 0; i--)
	{
		sum = 0;
		for (l = i + 1; l < N; l++)
			sum += m[i][l] * x[l];

		x[i] = (b[i] - sum) / m[i][i];
	}
	sum = 0;

	for (l = 1; l < N; l++)
		sum += m[0][l] * x[l];

	x[0] = (b[0] - sum) / m[0][0];

}

void lu_ssolve(std::vector<double>& x, int* p, std::vector<double>& b0, double** m, int N) {

    lu_ssolve(x.data(), p, b0.data(), m, N);

}

void residual(std::vector<double>& R, std::vector<std::vector<double>>&J, 
    std::vector<double>& n, std::vector<std::vector<double>>&A, 
        std::vector<double>& A0, std::vector<std::vector<double>>&v, 

        std::vector<int>&b, std::vector<int>& bs, int N, int M) {
	// da articolo di GODIN

	int i, j, l, k;

	/*calcolo di R*/
	for (l = 0; l < M; l++)
		R[l] = 0.;
	for (l = 0; l < M; l++)
	{
		for (j = 0; j < N - M; j++)
		{
			R[l] += A[l][bs[j]] * n[bs[j]];
		}
		for (i = 0; i < M; i++) 
			R[l] += A[l][b[i]] * n[b[i]];
		R[l] += -A0[l];
	}

	/*calolo di J*/

	for (l = 0; l < M; l++)
		for (i = 0; i < M; i++)
			J[l][i] = 0.;
	for (l = 0; l < M; l++)
	{
		for (k = 0; k < M; k++)
		{
			for (j = 0; j < N - M; j++)
			{
				J[l][k] += A[l][bs[j]] * n[bs[j]] * (v[j][k] / n[b[k]]);
			}
			J[l][k] += A[l][b[k]];
		}
	}
}

// clears a matrix
void matrix_free(double **matrix, int n)
{
	int i;

	for (i = 0; i < n; i++)
		free(matrix[i]);

	free(matrix);
}

void lu_inv(std::vector<std::vector<double>>& AINV, std::vector<std::vector<double>>& A, int N)
{
	int i, j, k;
	double **m;
	int *p;
	double *bas_tmp = new double[N];
	double *x = new double[N];
	int imax;
	double amax;
	double tmp;
	double mm;

	p = (int *)calloc(N, sizeof(int));
	m = matrix_alloc(N, N);

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			m[i][j] = A[i][j];
	for (k = 0; k < N - 1; k++)
	{
		imax = 0;
		amax = 0.0;
		for (i = k; i < N; i++)
		{
			if (amax < fabs(m[i][k]))
			{
				amax = fabs(m[i][k]);
				imax = i;
			}
		}

		p[k] = imax;
		for (j = k; j < N; j++)
		{
			tmp = m[k][j];
			m[k][j] = m[imax][j];
			m[imax][j] = tmp;
		}
		for (i = k + 1; i < N; i++)
		{
			mm = -m[i][k] / m[k][k];
			m[i][k] = mm;
			for (j = k + 1; j < N; j++)
				m[i][j] += m[i][k] * m[k][j];
		}
	}
	for (i = 0; i < N; i++)
		bas_tmp[i] = 0.0;
	for (i = 0; i < N; i++)
	{
		bas_tmp[i] = 1.0;

		lu_ssolve(x, p, bas_tmp, m, N);
		for (j = 0; j < N; j++)
			AINV[j][i] = x[j];

		bas_tmp[i] = 0.0;
	}
	free(p);
	matrix_free(m, N);
}

// compute det(**A) with LU scomposition
double lu_det(double** A, int N)
{
	int i, j, k;
	double det = 1;
	double** m;
	int *p;

	p = (int *)calloc(N, sizeof(int));
	m = matrix_alloc(N,N);
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			m[i][j] = A[i][j];

	/* DECOMPOSIZIONE LU */
	for (k = 0; k < N - 1; k++)
	{
		int imax = 0;
		double amax = 0.0;
		for (i = k; i < N; i++)
			if (amax < fabs(m[i][k]))
			{
				amax = fabs(m[i][k]);
				imax = i;
			}

		p[k] = imax;
		for (j = k; j < N; j++)
		{
			double tmp;

			tmp = m[k][j];
			m[k][j] = m[imax][j];
			m[imax][j] = tmp;
		}
		for (i = k + 1; i < N; i++)
		{
			double mm;

			mm = -m[i][k] / m[k][k];
			m[i][k] = mm;
			for (j = k + 1; j < N; j++)
				m[i][j] += m[i][k] * m[k][j];
		}
	}

	/* CALCOLO DETERMINANTE */
	for (i = 0; i < N; i++)
		det *= m[i][i];
	for (k = 0; k < N - 1; k++)
		if (p[k] != k)
			det *= -1.0;

	free(p);
	matrix_free(m, N);

	return det;

}

// Compute determinant of matrix A with LU decomposition
double lu_pivtot(std::vector<std::vector<double>> A, int N) {
    int i, j, k;
    double det = 1.0;
    std::vector<std::vector<double>> m(N, std::vector<double>(N, 0.0));
    std::vector<int> p(N, 0);

    // Copy A into m
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            m[i][j] = A[i][j];
        }
    }

    // LU decomposition with full pivoting
    for (k = 0; k < N - 1; k++) {
        int imax = k;
        int jmax = k;
        double maxVal = 0.0;

        // Find the largest absolute value in the submatrix m[k:N][k:N]
        for (i = k; i < N; i++) {
            for (j = k; j < N; j++) {
                if (fabs(m[i][j]) > maxVal) {
                    maxVal = fabs(m[i][j]);
                    imax = i;
                    jmax = j;
                }
            }
        }

        // Check for a zero pivot element
        if (maxVal == 0) {
            return 0;  // Matrix is singular
        }

        // Row swap for partial pivoting
        if (imax != k) {
            std::swap(m[imax], m[k]);
            det *= -1;  // Change the sign of determinant each row swap
        }

        // Column swap for full pivoting (not typically necessary for determinant calculation, omitted here)

        // Perform the elimination below row k
        for (i = k + 1; i < N; i++) {
            double factor = m[i][k] / m[k][k];
            m[i][k] = 0;  // The lower part of the k-th column becomes zero
            for (j = k + 1; j < N; j++) {
                m[i][j] -= factor * m[k][j];
            }
        }
    }

    // Calculate the determinant as the product of the diagonal elements
    for (i = 0; i < N; i++) {
        det *= m[i][i];
    }

    return det;
}

// compute det(**A) with LU decomposition and total pivot
double lu_det(std::vector<std::vector<double>> A, int N)
{
    int i, j, k;
    double det = 1;
    std::vector<std::vector<double>> m(N, std::vector<double>(N));
    int *p;
    int *q;  // For column swaps

    p = (int *)calloc(N, sizeof(int));
    q = (int *)calloc(N, sizeof(int));

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            m[i][j] = A[i][j];

    for (i = 0; i < N; i++) {
        p[i] = i;
        q[i] = i;
    }

    /* LU Decomposition with total pivot */
    for (k = 0; k < N - 1; k++)
    {
        int imax = k, jmax = k;
        double amax = 0.0;
        for (i = k; i < N; i++)
            for (j = k; j < N; j++)
                if (amax < fabs(m[i][j]))
                {
                    amax = fabs(m[i][j]);
                    imax = i;
                    jmax = j;
                }

        // Row swap
        if (imax != k)
        {
            for (j = 0; j < N; j++)
            {
                double tmp = m[k][j];
                m[k][j] = m[imax][j];
                m[imax][j] = tmp;
            }
            int tmp = p[k];
            p[k] = p[imax];
            p[imax] = tmp;
        }

        // Column swap
        if (jmax != k)
        {
            for (i = 0; i < N; i++)
            {
                double tmp = m[i][k];
                m[i][k] = m[i][jmax];
                m[i][jmax] = tmp;
            }
            int tmp = q[k];
            q[k] = q[jmax];
            q[jmax] = tmp;
        }

        for (i = k + 1; i < N; i++)
        {
            double mm = -m[i][k] / m[k][k];
            m[i][k] = mm;
            for (j = k + 1; j < N; j++)
                m[i][j] += m[i][k] * m[k][j];
        }
    }

    /* Calculate Determinant */
    for (i = 0; i < N; i++)
        det *= m[i][i];
    for (k = 0; k < N - 1; k++)
        if (p[k] != k)
            det *= -1.0;
    for (k = 0; k < N - 1; k++)
        if (q[k] != k)
            det *= -1.0;

    free(p);
    free(q);

    return det;
}

// solve an LU system
void lu_sistema(std::vector<double>& x, std::vector<std::vector<double>>& A,std::vector<double>& b, int N)
{
	int i, j, k;
	double **m;
	int *p = new int[N];
	int imax;
	double amax;
	double tmp;
	double mm;

	m = matrix_alloc(N, N);

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			m[i][j] = A[i][j];

	/* DECOMPOSIZIONE LU */
	for (k = 0; k < N - 1; k++)
	{
		imax = 0;
		amax = 0.0;
		for (i = k; i < N; i++)
			if (amax < fabs(m[i][k]))
			{
				amax = fabs(m[i][k]);
				imax = i;
			}

		p[k] = imax;
		for (j = k; j < N; j++)
		{
			tmp = m[k][j];
			m[k][j] = m[imax][j];
			m[imax][j] = tmp;
		}
		for (i = k + 1; i < N; i++)
		{
			mm = -m[i][k] / m[k][k];
			m[i][k] = mm;
			for (j = k + 1; j < N; j++)
				m[i][j] += m[i][k] * m[k][j];
		}
	}

	/* RISOLUZIONE SISTEMA */
    std::vector<double> xx = x;
    std::vector<double> bb = b;
        
	lu_ssolve(xx, p, bb, m, N);

	for (int i = 0; i < x.size(); i++)
        x[i] = xx[i];
	for (int i = 0; i < b.size(); i++)
        b[i]=bb[i];
    
	matrix_free(m, N);
}

void matrix_prod(std::vector<std::vector<double>>&A, std::vector<std::vector<double>>&B, std::vector<std::vector<double>>&C, int N1, int N2) /*A = B*C */
{
	int i, j, k;

	for (i = 0; i < N1; i++)
		for (j = 0; j < N2; j++)
			A[i][j] = 0.;

	for (i = 0; i < N1; i++)
	{
		for (j = 0; j < N2; j++)
		{
			for (k = 0; k < N2; k++)
				A[i][j] += B[i][k] * C[k][j];
		}
	}
}

void sort(double *arr, int *perm, int n)
{
	int i, j;
	int tmp_perm;
	double tmp;

	for (i = 0; i < n; i++)
		perm[i] = i;

	for (j = 1; j < n; j++)
	{
		tmp = arr[j];
		tmp_perm = perm[j];
		i = j - 1;
		while (i >= 0 && arr[i] > tmp)
		{
			arr[i + 1] = arr[i];
			perm[i + 1] = perm[i];
			i--;
		}
		arr[i + 1] = tmp;
		perm[i + 1] = tmp_perm;
	}
}
// interpolazione spline cubica con libreria alglib.
double interpolateSpline(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    // Verifica che x e y abbiano la stessa dimensione
    if (x.size() != y.size()) {
        throw std::invalid_argument("x e y devono avere la stessa dimensione.");
    }

    // Converte i vettori std::vector in real_1d_array per alglib
    alglib::real_1d_array ax, ay;
    ax.setcontent(x.size(), x.data());
    ay.setcontent(y.size(), y.data());

    // Costruisce la spline cubica
    alglib::spline1dinterpolant s;
    alglib::spline1dbuildcubic(ax, ay, s);

    // Calcola il valore interpolato in xi
    double yi = alglib::spline1dcalc(s, xi);

    return yi;
}

// integrazione con regola dei trapezi
double trapz(const std::vector<double>& x , const std::vector<double>& y){
    double integral = 0.0 ;
    for (int i = 1; i < x.size(); i++)
        integral += 0.5*(x[i]-x[i-1])*(y[i] + y[i-1]) ;
    
    return integral ;
}


// Calcola i coefficienti della retta di interpolazione
std::vector<double> interpCoeff(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        // Dati non validi: vectors di lunghezza diversa o vuoti
        throw std::invalid_argument("Vettori x e y devono avere la stessa dimensione e non essere vuoti.");
    }

    double x_sum = 0, y_sum = 0, xy_sum = 0, xx_sum = 0;
    int n = x.size();

    for (int i = 0; i < n; ++i) {
        x_sum += x[i];
        y_sum += y[i];
        xy_sum += x[i] * y[i];
        xx_sum += x[i] * x[i];
    }

    double x_mean = x_sum / n;
    double y_mean = y_sum / n;

    double sxx = xx_sum - n * x_mean * x_mean;
    double sxy = xy_sum - n * x_mean * y_mean;

    if (std::abs(sxx) < std::numeric_limits<double>::epsilon()) {
        // Varianza di x troppo piccola: tutti i valori di x sono identici o quasi
        // Non è possibile calcolare una pendenza sensata
        return {0, y_mean}; 
    }

    double slope = sxy / sxx;
    double intercept = y_mean - slope * x_mean;

    return {slope, intercept};
}

// polyval
double Fbar(double x1, const std::vector<double>& p) {
    return p[0]*x1 + p[1]; // p[0] è m (slope), p[1] è q (intercept)
}

double interp(const std::vector<double>& x, const std::vector<double>& y, double xx) {
    
    int ia = 0;
    int ib = x.size() - 1; 

    if (xx > x.back()) {
        xx = x.back();
    }

    if (xx < x.front()) {
        xx = x.front();
    }

    for (int i = 0; i < x.size(); i++) {
        if (xx > x[i]) {
            ia = i;
        }
    }

    for (int i = x.size() - 1; i >= 0; i--) {
        if (xx < x[i]) {
            ib = i;
        }
    }

    if (x[ib] == x[ia]) {
        return y[ia]; 
    }

    return y[ib] - ((x[ib] - xx) * (y[ib] - y[ia]) / (x[ib] - x[ia]));

}


double adaptiveTrapz(std::vector<double> xx , std::vector<double> yy){
    
    std::vector<double> x = xx;
    std::vector<double> y = yy;
    
    double integral = 0.0 ;
    double Tol_I_rel = 1e-02 ;
    double I_min = 1e-03 ;

    bool control = true ;
    int iOK=0;
    while (control == true){
        int Lx = x.size() ;
        control = false ; 
        for (int i = iOK; i < Lx-2 ; i+=2){
            double I1 = 0.5 * std::abs( x[i] - x[i+2] )*( y[i] + y[i+2] ) ;
            // stima trapezio con punto medio Colonna eq.8
            double I2 = 0.5 * std::abs( x[i] - x[i+2] )*( 0.5 * y[i] + 0.5 * y[i+2] + y[i+1] ) ;
            // eq.9 relativo
            double EI_rel = std::abs( (I2 - I1) / I2 );
              // Finchè E>Toll ed I2 significativo infittisce dominio e integranda
            if (EI_rel > Tol_I_rel && I2 > I_min){

                control = true ;

                std::vector<double> xin = std::vector<double> ( x.begin(), x.begin() + i + 1) ;
                std::vector<double> xfin = std::vector<double> ( x.begin() + i + 2, x.end()) ;
                
                std::vector<double> yin = std::vector<double> ( y.begin(), y.begin() + i + 1) ;
                std::vector<double> yfin = std::vector<double> ( y.begin() + i + 2, y.end()) ;
                
                double xa = x[i] ;
                double xb = x[i+2] ;
                double xc = x[i+1] ;

                double yac = interpolateSpline(xx,yy,0.5 * ( xa + xc )) ;
                double ycb = interpolateSpline(xx,yy,0.5 * ( xa + xc )) ;

                x = concatenate(xin,std::vector<double> {0.5*(xa+xc)},std::vector<double> {xc}, std::vector<double> {0.5*(xc+xb)},  xfin ) ;
                y = concatenate(yin,std::vector<double> {yac},std::vector<double> {y[i+1]}, std::vector<double> {ycb},  yfin ) ;

                Lx = x.size() ;
                break ;
            }else{
                iOK = i+2 ;
            }
        }
    }
    return integral = trapz(x,y);
}

double max_double(std::vector<double> arr, int N)
{
	int j;
    double max;
	
	max = fabs(arr[0]);

	for (j=0;j<N;j++) if(fabs(arr[j]) > max) max = fabs(arr[j]);

	return max;
}

// Funzione per convertire std::vector<std::vector<double>> in alglib::real_2d_array
alglib::real_2d_array convertToAlglibArray(const std::vector<std::vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    alglib::real_2d_array alglibMatrix;
    alglibMatrix.setlength(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            alglibMatrix[i][j] = matrix[i][j];
        }
    }
    return alglibMatrix;
}

// Funzione per calcolare il determinante utilizzando la decomposizione QR
double DetQR(const std::vector<std::vector<double>>& matrix , int useless_size) {
    // Converti la matrice in formato ALGLIB
    alglib::real_2d_array alglibMatrix = convertToAlglibArray(matrix);

     // Esegui la decomposizione QR
    alglib::real_1d_array tau;
    alglib::ae_int_t m = matrix.size();
    alglib::ae_int_t n = matrix[0].size();
    alglib::rmatrixqr(alglibMatrix, m, n, tau);

    // Calcola il determinante come il prodotto dei termini diagonali di R
    double determinant = 1.0;
    for (int i = 0; i < std::min(m, n); ++i) {
        determinant *= alglibMatrix[i][i];
    }

    return -determinant;
}

double DetLU(const std::vector<std::vector<double>>& matrix , int uselessdimension) {
    // Verifica che la matrice sia quadrata
    int n = matrix.size();
    if (n == 0 || matrix[0].size() != n) {
        throw std::invalid_argument("La matrice deve essere quadrata");
    }

    // Converti la matrice in alglib::real_2d_array
    alglib::real_2d_array a = convertToAlglibArray(matrix);

    // Decomposizione LU con pivoting totale
    alglib::ae_int_t nn = n ;

    // Effettua la decomposizione LU
    alglib::integer_1d_array pivots;
    alglib::rmatrixlu(a, n, n, pivots) ;

    // Calcola il determinante usando i fattori LU
    double det = 1.0;
    for (int i = 0; i < n; ++i) {
        det *= a[i][i];
    }

    // Calcola il numero di permutazioni dal vettore di pivot
    int num_permutations = 0;
    for (int i = 0; i < n; ++i) {
        if (pivots[i] != i) {
            num_permutations++;
        }
    }

    // Se il numero di permutazioni è dispari, il determinante cambia segno
    if (num_permutations % 2 != 0) {
        det = -det;
    }

    return det;
    
}

// Funzione helper per libcurl per scrivere i dati in un stringstream
int WriteCallback(void* contents, int size, int nmemb, void* userp) {

    std::stringstream* stream = static_cast<std::stringstream*>(userp);
    stream->write(static_cast<const char*>(contents), size * nmemb);

    return size * nmemb;

}
