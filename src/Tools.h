#ifndef TOOLS_H
#define TOOLS_H

#define KB 1.380650524e-23     // [J/K] Bolztmann constant 
#define hPlanck 6.62606896e-34 // [J s] Planck constant 
#define eps0 8.8541878176e-12  // [F/m] void dielectric constant 
#define qe 1.602176634e-19     // [C]   Electron unit charge

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <typeindex>
#include "alglib/source/linalg.h"
#include <stdexcept>

/* using std::log10; */

/* Kronecker delta */
int delta(int i, int j) ; 

// fattoriale 
double factorial(int n) ; 


template<typename T>
std::vector<T> arange(T start, T end, T step) {
    std::vector<T> values;

    // Verifica che lo step sia diverso da zero
    if (step == 0) {
        throw std::invalid_argument("Step cannot be zero");
    }

    // Verifica che l'intervallo sia coerente con il segno dello step
    if ((step > 0 && start >= end) || (step < 0 && start <= end)) {
        return values; // Vettore vuoto
    }

    if (step > 0) {
        for (T value = start; value < end; value += step) {
            values.push_back(value);
        }
    } else {
        for (T value = start; value > end; value += step) {
            values.push_back(value);
        }
    }
    
    return values;
}
//vettore linspace
template<typename T>
std::vector<T> linspace(T start, T end, int num) {
    
    std::vector<T> values;
    if (num <= 1) {std::cerr<<"No linspace with"<<num<<"point"<<std::endl;}
    values.reserve(num);
    T step = (end - start) / static_cast<T>(num - 1);
    for(int i = 0; i < num; ++i) {
        if (i == num - 1)
            values.push_back(end);
        else 
            values.push_back(start + step * i);
    }
    return values;
}

//vettore logaritmico equispaziato
template<typename T>
std::vector<T> logspace(T start,T end,int num){
    std::vector<T> logspaced;
        T delta = (end - start) / static_cast<T>(num - 1);
        for (int i = 0; i < num; ++i) {
            logspaced.push_back(std::pow(static_cast<T>(10), start + delta * i));
        }
        return logspaced;    
}

double **matrix_alloc(int N, int M) ;

void lu_ssolve(double *x, int *p, double *b0, double **m, int N);
void residual(std::vector<double>& R, std::vector<std::vector<double>>&J, 
    std::vector<double>& n, std::vector<std::vector<double>>&A,
         std::vector<double>& A0, std::vector<std::vector<double>>&v, 
            std::vector<int>&b, std::vector<int>& bs, int N, int M ) ;

// clears a matrix
void matrix_free(double **matrix, int n) ;

void lu_inv(std::vector<std::vector<double>>& AINV, std::vector<std::vector<double>>& A, int N);

// compute det(**A) with LU scomposition
double lu_det(double** A, int N) ;

// Compute determinant of matrix A with LU decomposition
double lu_pivtot(std::vector<std::vector<double>> A, int N) ;

// compute det(**A) with LU decomposition and total pivot
double lu_det(std::vector<std::vector<double>> A, int N)  ;

// solve an LU system
void lu_sistema(std::vector<double>& x, std::vector<std::vector<double>>& A,
    std::vector<double>& b, int N) ;

/*A = B*C */
void matrix_prod(std::vector<std::vector<double>>&A, std::vector<std::vector<double>>&B, 
    std::vector<std::vector<double>>&C, int N1, int N2) ;

void sort(double *arr, int *perm, int n) ;

// interpolazione spline cubica con libreria alglib.
double interpolateSpline(const std::vector<double>& x, const std::vector<double>& y, double xi) ;

// integrazione con regola dei trapezi
double trapz(const std::vector<double>& x , const std::vector<double>& y) ;


// Calcola i coefficienti della retta di interpolazione
std::vector<double> interpCoeff(const std::vector<double>& x, const std::vector<double>& y) ;

// polyval
double Fbar(double x1, const std::vector<double>& p) ;

double interp(const std::vector<double>& x, const std::vector<double>& y, double xx) ;

// Funzione per concatenare due std::vector
template <typename T>
std::vector<T> concVecs(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> results;
    results.reserve(v1.size() + v2.size()); 
    results.insert(results.end(), v1.begin(), v1.end());
    results.insert(results.end(), v2.begin(), v2.end());
    return results;
}

// Base della ricorsione: quando rimangono solo due vectors da concatenare
template <typename T>
std::vector<T> concatenate(const std::vector<T>& v1, const std::vector<T>& v2) {
    return concVecs(v1, v2);
}

// Ricorsione per numero variabile di argomenti
template <typename T, typename... Args>
std::vector<T> concatenate(const std::vector<T>& v1, const std::vector<T>& v2, Args... args) {
    std::vector<T> vettoreParziale = concVecs(v1, v2);
    return concatenate(vettoreParziale, args...);
}

double adaptiveTrapz(std::vector<double> xx , std::vector<double> yy) ;
//stampa ogni 1000 componenti

template <typename T>
void print(const std::vector<T>& v) {
    for (int i = 0; i<v.size(); i++)
        if (1000*i<v.size())        
            std::cout<<std::setw(14)<<v[1000*i]<<std::setw(14);
    std::cout<<v.back();
    std::cout<<std::endl;

}

//stampa ogni evr componenti
template <typename T>
void print(const std::vector<T>& v,int evr) {
    for (int i = 0; i<v.size(); i++)
        if (evr*i<v.size())        
            std::cout<<std::setw(14)<<v[evr*i]<<std::setw(14);
    std::cout<<v.back();
    std::cout<<std::endl;
}

// scrive vettori su file incolonnandoli
template<typename T, typename... Vettori>
void writeVectorsToFile(const std::string& filename, const T& first , const Vettori&... vectors){
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file '" << filename << "'." << std::endl;
        return ; 
    }
    
    const int numRows = first.size();
    // lambda function per grandezze vettori controlla abbiano tutti dimensioni numRows
    auto checkSize = [&](int size) {
        return ((vectors.size() == size) && ...);
    };
    // interrompe se dimensioni diverse
    if (!checkSize(numRows)) {
        std::cerr << "Error: vectors have different sizes." << std::endl;
        return;
    }

    for (int i = 0; i < numRows; i++){
        file << first[i];  
        ((file << "," << vectors[i]), ...); 
        file << "\n";
    }
    file.close();
}

// estrae k-esima colonna da oggetti std::vector<std::vector<T>>
template<typename T>
std::vector<T> column(const std::vector<std::vector<T>>& matrix, int k) {
    std::vector<T> col;
    for (const auto& row : matrix) {
        if (k < row.size()) { 
            col.push_back(row[k]);
        }
        else{
            throw std::invalid_argument("Fuori range di matrice");
        }
    }
    return col;
}

template<typename T, typename... Vettori>
void AppendColumnsToFile(const std::string& filename, const T& first, const Vettori&... vectors) {
    std::ifstream infile(filename);
    std::vector<std::string> lines;
    std::string line;
    
    while (std::getline(infile, line)) {
        lines.push_back(line);
    }
    infile.close();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file '" << filename << "'." << std::endl;
        return;
    }

    const int numRows = first.size();
    auto checkSize = [&](int size) {
        return ((vectors.size() == size) && ...);
    };

    if (!checkSize(numRows)) {
        std::cerr << "Error: vectors have different sizes." << std::endl;
        return;
    }

    for (int i = 0; i < numRows; i++) {
        if (i < lines.size()) {
            file << lines[i] << "," << first[i];
        } else {
            file << first[i];
        }
        ((file << "," << vectors[i]), ...);
        file << "\n";
    }
    file.close();
}

double max_double(std::vector<double> arr, int N) ;

template < typename T >
std::type_index tipo ( T* pointer ) {return typeid(*pointer);}

// Funzione per convertire std::vector<std::vector<double>> in alglib::real_2d_array
alglib::real_2d_array convertToAlglibArray(const std::vector<std::vector<double>>& matrix) ;

// Funzione per calcolare il determinante utilizzando la decomposizione QR
double DetQR(const std::vector<std::vector<double>>& matrix , int useless_size) ;

double DetLU(const std::vector<std::vector<double>>& matrix , int uselessdimension) ;

// Funzione helper per libcurl per scrivere i dati in un stringstream
int WriteCallback(void* contents, int size, int nmemb, void* userp) ;

#endif
