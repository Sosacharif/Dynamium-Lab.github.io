#include "blast.h"
#include<iostream>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;
using real = double;

blast::real determinant (blast::Mat3 matrix);
double invariant_vectoriel(const blast::Mat3& matrix);
blast::real trace (blast::Mat3 matrix); 
blast::Vec3 valeurspropres (blast::Mat3 matrix);
blast::Mat3 identitymatrix();
blast::Mat3 produittranspose(blast::Vec3 n, blast::Vec3 ntranspose);
blast::real lambda ();
blast::Vec3 invariantvectoriel(blast::Mat3 matrix);
double x1, x2, x3, yI, y2, y3, z1, z2, z3; 
blast::real angle();
double cosfi = cos(0);  
double sinfi = sin(0);
float AngleDeRotation(const blast::Mat3& matricederotation);
blast::Vec3 Vecteurrotation (blast::Mat3 matricederotation, blast::Vec3 x);
blast::Vec3 vitessedupointp(blast::Vec3 vecteurp, blast::Mat3 matricederotationAB, blast::Vec3 vecteurw);
blast::Vec3 vitesserelative(blast::Vec3 vecteurp, blast::Vec3 vecteurw, blast::Vec3 produitvectoriel);
blast::Vec3 rotated;

// Fonction pour calculer le déterminant
blast::real determinant (blast::Mat3 matrix) {
    return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)) -
           matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0)) +
           matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
}

// Fonction pour calculer l'invariant vectoriel 
double invariant_vectoriel(const blast::Mat3& matrix) {
    return determinant(matrix);
}

// Fonction pour calculer les valeurs propres
blast::Vec3 valeurspropres (blast::Mat3 matrix) {
    blast::Vec3 eigen_vals;  
    std::cout << "Valeurs propres: ";
    std::cout << eigen_vals[0] << ", " << eigen_vals[1] << ", " << eigen_vals[2] << std::endl;
    return eigen_vals;
}
// Fonction qui calcule la trace d'une matrice 
blast::real trace (blast::Mat3 matrix) {
    return matrix(0, 0) + matrix(1, 1) + matrix(2, 2);
}  

// Fonction pour calculer l'angle de rotation à partir de la trace de la matrice de rotation
float AngleDeRotation(const blast::Mat3& matricederotation) {
    // Trace de la matrice de rotation
    real trace_val = matricederotation(0, 0) + matricederotation(1, 1) + matricederotation(2, 2);
    // Calcul de l'angle de rotation
    return acos((trace_val - 1) / 2);
}

// Fonction pour calculer l'axe de rotation
blast::Vec3 axederotation(const blast::Mat3& matricederotation) {
    blast::Vec3 axis;
    axis.x = matricederotation(2,1) - matricederotation(1,2);
    axis.y = matricederotation(0,2) - matricederotation(2,0);
    axis.z = matricederotation(1,0) - matricederotation(0,1);
    
    double norm = sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
    if (norm > 1e-6) {
        axis.x /= norm;
        axis.y /= norm;
        axis.z /= norm;
    }
    return axis;
}
// Les conditions dans lesquelles le calcul de l'axe et l'angle de rotation marche est que la matroice doit etre orthogonaleavec un determinant de 1.

// Fonction pour créer une matrice de rotation à partir des angles d'Euler (ZXZ)
blast::Mat3 matricederotationdeuler(blast::real alpha, blast::real beta, blast::real gamma) {
    // Calcul des cosinus et sinus des angles
    blast::real cosAlpha = cos(alpha);
    blast::real sinAlpha = sin(alpha);
    blast::real cosBeta = cos(beta);
    blast::real sinBeta = sin(beta);
    blast::real cosGamma = cos(gamma);
    blast::real sinGamma = sin(gamma);

    // Matrice de rotation R = Rz(gamma) * Rx(beta) * Rz(alpha)
    blast::Mat3 R;

    R(0, 0) = cosGamma * cosAlpha - sinGamma * sinAlpha * cosBeta;
    R(0, 1) = -sinGamma * cosAlpha - cosGamma * sinAlpha * cosBeta;
    R(0, 2) = sinAlpha * sinBeta;

    R(1, 0) = sinGamma * cosAlpha + cosGamma * sinAlpha * cosBeta;
    R(1, 1) = -cosGamma * cosAlpha + sinGamma * sinAlpha * cosBeta;
    R(1, 2) = -sinAlpha * sinBeta;

    R(2, 0) = sinGamma * sinBeta;
    R(2, 1) = cosGamma * sinBeta;
    R(2, 2) = cosBeta;

    return R;
}

// Fonction qui calcule l'invariant vectoriel d'une matrice
blast::Vec3 invariantvectoriel(blast::Mat3 matrix) {
    return blast::Vec3(
        1 / 2 * (matrix(2, 1) - matrix(1, 2)), 
        1 / 2 * (matrix(0, 2) - matrix(2, 0)), 
        1 / 2 * (matrix(1, 0) - matrix(0, 1))
    );
}

// Fonction pour calculer le produit de la transposée d'un vecteur par lui-même
blast::Mat3 produittranspose(blast::Vec3 n, blast::Vec3 ntranspose) {
    return blast::Mat3(
                n.x * ntranspose.x, n.x * ntranspose.y, n.x * ntranspose.z,
                n.y * ntranspose.x, n.y * ntranspose.y, n.y * ntranspose.z,
                n.z * ntranspose.x, n.z * ntranspose.y, n.z * ntranspose.z
            );
        }
        

// Fonction pour calculer la matrice de projection
blast::Mat3 matriceprojection(const blast::Mat3& identitymatrix, const blast::Mat3& produittranspose) {
    blast::Mat3 resultat;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            resultat(i, j) = identitymatrix(i, j) - produittranspose(i, j);
        }
    }
    return resultat;
}

// Fonction pour calculer la matrice de réflexion
blast::Mat3 matricedereflexion(const blast::Mat3& identitymatrix, const blast::Mat3& produittranspose) {
    blast::Mat3 resultat;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            resultat(i, j) = identitymatrix(i, j) - 2 * produittranspose(i, j);
        }
    }
    return resultat;
}

void printVector(const vector<double>& vec) {
    for (double val : vec) {
        cout << val << " ";
    }
    cout << endl;
}

// Fonction pour afficher un vecteur blast::Vec3
void printVector(const blast::Vec3& vec) {
    cout << vec.x << " " << vec.y << " " << vec.z << endl;
}

// Fonction pour effectuer la rotation d'un vecteur
void rotateAndPrint(const blast::Vec3& x, const blast::Mat3& matriceDeRotation) {
    // Nombre de points à afficher
    int numPoints = 10;
    printVector(x);
    for (int i = 1; i < numPoints; ++i) {
        double t = static_cast<double>(i) / (numPoints - 1); 
        blast::Vec3 rotated; 
        rotated.x = (matriceDeRotation(0, 0) * x.x + matriceDeRotation(0, 1) * x.y + matriceDeRotation(0, 2) * x.z);
        rotated.y = (matriceDeRotation(1, 0) * x.x + matriceDeRotation(1, 1) * x.y + matriceDeRotation(1, 2) * x.z);
        rotated.z = (matriceDeRotation(2, 0) * x.x + matriceDeRotation(2, 1) * x.y + matriceDeRotation(2, 2) * x.z);
        printVector(rotated);
    }
    blast::Vec3 finalVector;
    finalVector.x = (matriceDeRotation(0, 0) * x.x + matriceDeRotation(0, 1) * x.y + matriceDeRotation(0, 2) * x.z);
    finalVector.y = (matriceDeRotation(1, 0) * x.x + matriceDeRotation(1, 1) * x.y + matriceDeRotation(1, 2) * x.z);
    finalVector.z = (matriceDeRotation(2, 0) * x.x + matriceDeRotation(2, 1) * x.y + matriceDeRotation(2, 2) * x.z);
    printVector(finalVector);
}
// Produit vectoriel entre deux vecteurs
blast::Vec3 produitvectoriel(const blast::Vec3& a, const blast::Vec3& b) {
    return blast::Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// Fonction pour calculer la vitesse du point p
blast::Vec3 vitessedupointp(blast::Vec3 x, blast::Vec3 vecteurp, blast::Mat3 matricederotationAB, blast::Vec3 vecteurw) {
    blast::Vec3 vitesserelative = vecteurp + produitvectoriel(x, vecteurw);
    return matricederotationAB * vitesserelative;
}
namespace blast {
    void print(const Vec3& vec) {
        std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
    }
}

void afficher_matrice(const blast::Mat3& mat) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << mat(i, j) << "\t";
        }
        std::cout << std::endl;
    }
}
void afficher_vecteur(const blast::Vec3& vec) {
    std::cout << "[" << vec.x << ", " << vec.y << ", " << vec.z << "]" << std::endl;
}
std::ostream& operator<<(std::ostream& os, const blast::Vec3& vec) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return os;
}
int main() {
    using namespace blast;
    Mat3 matrix(1, 2, 3,
                4, 5, 6,
                7, 8, 9);

    Mat3 identitymatrix(1, 0, 0,
                        0, 1, 0,
                        0, 0, 1);

    Mat3 matricederotation(sqrt(3)/2, -0.5, 0,
                           0.5, -sqrt(3)/2, 0,
                           0, 0, 1);

    Vec3 x(1, 2, 3);
    Vec3 vecteurw(1, 0, 0);
    Vec3 vecteurp(3, 2, 1);
    Mat3 matricederotationAB(1, 0, 0,
                             0, 1, 0,
                             0, 0, 1);

    cout << "Déterminant : " << determinant(matrix) << endl;
    cout << "Trace : " << trace(matrix) << endl;
    Vec3 axis = axederotation(matricederotation);
    cout << "Axe de rotation : " << axis.x << ", " << axis.y << ", " << axis.z << endl;
    rotateAndPrint(x,matricederotation);
    cout << "Angle de rotation : " << AngleDeRotation(matricederotation) << " radians" << endl;
    double invariant = invariant_vectoriel(matrix);
    std::cout << "L'invariant vectoriel : " << invariant << std::endl;
    // Angles d'Euler en radians
     real alpha = M_PI / 4; // 45 degrés
     real beta = M_PI / 6;  // 30 degrés
     real gamma = M_PI / 3; // 60 degrés
 
    // Calcul de la matrice de rotation
     Mat3 rotationMatrix = matricederotationdeuler(alpha, beta, gamma);
 
    // Affichage de la matrice de rotation
     cout << "Matrice de rotation (ZXZ) :" << endl;
     cout << rotationMatrix(0, 0) << " " << rotationMatrix(0, 1) << " " << rotationMatrix(0, 2) << endl;
     cout << rotationMatrix(1, 0) << " " << rotationMatrix(1, 1) << " " << rotationMatrix(1, 2) << endl;
     cout << rotationMatrix(2, 0) << " " << rotationMatrix(2, 1) << " " << rotationMatrix(2, 2) << endl;

    // Affiche la matrice de projection
    Mat3 projMatrix = matriceprojection(identitymatrix, produittranspose(x, vecteurw));
    cout << "Matrice de projection : " << endl;
    blast::print(projMatrix);

    // Affiche le vecteur vitesse
    Vec3 Vecteurvitesse = vitessedupointp(x, vecteurp, matricederotationAB, vecteurw);
    std::cout << "Vitesse du point P : " << Vecteurvitesse << std::endl;

    Vec3 Valeurspropres = valeurspropres(matrix);
    std::cout << "Valeurspropres : " << valeurspropres << std::endl;
    
    // Affiche la matrice de reflexion
    Mat3 ReflexMatrix = matricedereflexion(identitymatrix, produittranspose(x, vecteurw));
    cout << "Matrice de reflexion : " << endl;
    blast::print(ReflexMatrix);
     return 0;
}

