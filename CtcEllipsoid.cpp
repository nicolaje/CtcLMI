#include "CtcEllipsoid.h"

using namespace std;
using namespace ibex;

void CtcEllipsoid::contract(IntervalVector& box)
{
    mkMatrixArray(box);
    CtcLMI ctcLMI(*matrices);
    ctcLMI.contract(box);
}

void CtcEllipsoid::mkMatrixArray(const IntervalVector& box)
{
    // On inverse P
    Matrix P_inv=P;//(P.nb_rows(), P.nb_cols());

//    real_inverse(P, P_inv);

    // On cree le tableau de matrices
    // de taille 1+size(box)
    // Chaque matrice de taille (1+3*size(box),1+3*size(box))
    matrices = new MatrixArray(1 + box.size(), 1 + 3 * box.size(), 1 + 3 * box.size());

    // F0
    Matrix f0 = Matrix::zeros(1 + 3 * box.size(), 1 + 3 * box.size());

    f0[0][0] = pow(R,2);
    f0.put(1, 1, P_inv);
    f0.put(0, 1, -C, true);
    f0.put(1, 0, -C, false);
    for (int i = 0; i < box.size(); i++)
    {
        f0[1 + P_inv.nb_cols() + 2 * i][1 + P_inv.nb_cols() + 2 * i] = -box[i].lb();
        f0[1 + P_inv.nb_cols() + 2 * i + 1][1 + P_inv.nb_cols() + 2 * i + 1] = box[i].ub();
    }
    (*matrices)[0] = -f0;
    for (int i = 1; i < matrices->size(); i++)
    {
        (*matrices)[i] = Matrix::zeros((*matrices)[i].nb_rows(), (*matrices)[i].nb_rows());
        (*matrices)[i][i][0] = 1;
        (*matrices)[i][0][i] = 1;
        (*matrices)[i][1 + P_inv.nb_rows() + 2 * (i - 1) ][1 + P_inv.nb_rows() + 2 * (i - 1) ] = 1;
        (*matrices)[i][1 + P_inv.nb_rows() + 2 * (i - 1) + 1][1 + P_inv.nb_rows() + 2 * (i - 1) + 1] = -1;
    }
}