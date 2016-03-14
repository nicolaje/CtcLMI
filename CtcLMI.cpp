#include "CtcLMI.h"
#include <cstdio>

using namespace std;
using namespace ibex;

void CtcLMI::mkSolver(ibex::MatrixArray& matrices, const IntervalVector &box)
{
    // On modifie F0 en fonction des nouvelles bornes de boite
    for (int i = 0; i < box.size(); i++)
    {
        matrices[0][1 + box.size() + 2 * i][1 + box.size() + 2 * i] = box[i].lb();
        matrices[0][1 + box.size() + 2 * i + 1][1 + box.size() + 2 * i + 1] = -box[i].ub();
    }

    solver.setParameterType(SDPA::PARAMETER_DEFAULT);
    solver.setParameterMaxIteration(100);
    //    solver.setNumThreads(4);
    //    solver.printParameters(stdout);

    solver.inputConstraintNumber(matrices.size());
    solver.inputBlockNumber(1);

    solver.inputBlockSize(1, matrices[0].nb_cols());
    solver.inputBlockType(1, SDPA::SDP);

    solver.initializeUpperTriangleSpace();
    //    solver.inputCVec(1, 0);
    //    solver.inputCVec(2, -1);

    for (int i = 0; i < matrices.size(); i++)
        for (int j = 1; j < matrices[i].nb_cols() + 1; j++)
            for (int k = 1; k < matrices[i].nb_rows() + 1; k++)
                solver.inputElement(i, 1, k, j, matrices[i][k - 1][j - 1]);
}

void CtcLMI::contract(IntervalVector& box)
{
    for (int i = 1; i <= box.size(); i++)
    {
        mkSolver(matrices, box);
        // Minimize (lower bound)
        solver.inputCVec(i, 1);
        solver.initializeUpperTriangle();
        solver.initializeSolve();
        solver.solve();
        box[i - 1] &= Interval(solver.getResultXVec()[i - 1], POS_INFINITY);
        // cout << "solver.getResultXVec()[" << i << " - 1]: " << solver.getResultXVec()[i - 1] << endl;

        mkSolver(matrices, box);
        // Maximize (upper bound)
        solver.inputCVec(i, -1);
        solver.initializeUpperTriangle();
        solver.initializeSolve();
        solver.solve();
        box[i - 1] &= Interval(NEG_INFINITY, solver.getResultXVec()[i - 1]);
    }
}
