#include "CtcEllipsoids.h"

using namespace ibex;

void CtcEllipsoids::contract(ibex::IntervalVector& box)
{
    mkMatrixArray(box);

    cout << "les matrices: " << endl;

    for (int i = 0; i < (*matrices).size(); i++)
        cout << "F[" << i << "]: " << endl << (*matrices)[i] << endl;

    CtcLMI ctcLMI(*matrices);
    cout << "Apres construct. CtcLMI" << endl;
    ctcLMI.contract(box);
}

void CtcEllipsoids::mkMatrixArray(const ibex::IntervalVector& box)
{
    // On cree le tableau de matrices
    // de taille 1+size(box)
    // Chaque matrice de taille ((1+3*size(box))*nb_ellipses,(1+3*size(box))*nb_ellipse)
    matrices = new MatrixArray(1 + box.size(), (1 + 3 * box.size()) * Ps.size(), (1 + 3 * box.size()) * Ps.size());

    // On init a 0 les matrices
    for (int i = 0; i < matrices->size(); i++)
    {
        (*matrices)[i] = Matrix::zeros((*matrices)[i].nb_rows(), (*matrices)[i].nb_rows());
    }

    // Pour chaque ellipsoide
    for (int e = 0; e < Ps.size(); e++)
    {
        // On inverse P
        Matrix P_inv=Ps[e];//(Ps[e].nb_rows(), Ps[e].nb_cols());

//        real_inverse(Ps[e], P_inv);

        // F0
        Matrix f0 = Matrix::zeros(1 + 3 * box.size(), 1 + 3 * box.size());

        f0[0][0] = pow(Rs[e], 2);
        f0.put(1, 1, P_inv);
        f0.put(0, 1, -Cs[e], true);
        f0.put(1, 0, -Cs[e], false);
        for (int i = 0; i < box.size(); i++)
        {
            f0[1 + P_inv.nb_cols() + 2 * i][1 + P_inv.nb_cols() + 2 * i] = -box[i].lb();
            f0[1 + P_inv.nb_cols() + 2 * i + 1][1 + P_inv.nb_cols() + 2 * i + 1] = box[i].ub();
        }

        // On concatene
        (*matrices)[0].put(e * f0.nb_rows(), e * f0.nb_cols(), -f0); // = -f0;

        // 
        for (int i = 1; i < matrices->size(); i++)
        {
            Matrix fi = Matrix::zeros(1 + 3 * box.size(), 1 + 3 * box.size());
            fi[i][0] = 1;
            fi[0][i] = 1;
            fi[1 + P_inv.nb_rows() + 2 * (i - 1) ][1 + P_inv.nb_rows() + 2 * (i - 1)] = 1;
            fi[1 + P_inv.nb_rows() + 2 * (i - 1) + 1][1 + P_inv.nb_rows() + 2 * (i - 1) + 1] = -1;

            (*matrices)[i].put(e * (1 + 3 * box.size()), e * (1 + 3 * box.size()), fi);
            //            (*matrices)[i][i][0] = 1;
            //            (*matrices)[i][0][i] = 1;
            //            (*matrices)[i][1 + P_inv.nb_rows() + 2 * (i - 1) ][1 + P_inv.nb_rows() + 2 * (i - 1) ] = 1;
            //            (*matrices)[i][1 + P_inv.nb_rows() + 2 * (i - 1) + 1][1 + P_inv.nb_rows() + 2 * (i - 1) + 1] = -1;
        }
    }
}