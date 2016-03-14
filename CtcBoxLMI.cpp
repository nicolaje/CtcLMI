#include "CtcBoxLMI.h"

using namespace std;
using namespace ibex;

bool check_empty_box(const IntervalVector &box)
{
    for (int i = 0; i < box.size(); i++)
        if (box[i].is_empty())
            return true;
    return false;
}

void CtcBoxLMI::mkSolver(const ibex::MatrixArray& matrices, IntervalVector &box)
{
    if (!check_empty_box(box))
    {
        // Matrices codant la LMI concatenees avec la LMI box
        MatrixArray mA(matrices.size(), matrices.nb_rows() + 2 * box.size(), matrices.nb_cols() + 2 * box.size());

        // On integre la box-LMI dans F0,F1,...
        // F0
        Matrix F0 = Matrix::zeros(matrices.nb_rows() + 2 * box.size(), matrices.nb_cols() + 2 * box.size());

        mA[0] = F0;

        for (int i = 0; i < 2 * box.size(); i++)
            mA[0][i][i] = i % 2 ? -box.ub()[int(i / 2.)] : box.lb()[int(i / 2.)];

        // F1,F2...
        for (int i = 0; i < box.size(); i++)
        {
            mA[i + 1][2 * i][2 * i] = 1;
            mA[i + 1][2 * i + 1][2 * i + 1] = -1;
        }

        // On copie la LMI dans la Box-LMI
        for (int i = 0; i < matrices.size(); i++)
            mA[i].put(2 * box.size(), 2 * box.size(), matrices[i]);


        solver.setParameterType(SDPA::PARAMETER_DEFAULT);
        solver.setParameterMaxIteration(100);
        solver.setNumThreads(4);
        // solver.printParameters(stdout);

        solver.inputConstraintNumber(mA.size());
        solver.inputBlockNumber(1);

        solver.inputBlockSize(1, mA[0].nb_cols());
        solver.inputBlockType(1, SDPA::SDP);

        solver.initializeUpperTriangleSpace();

        for (int i = 0; i < mA.size(); i++)
            for (int j = 1; j < mA[i].nb_cols() + 1; j++)
                for (int k = 1; k < mA[i].nb_rows() + 1; k++)
                    solver.inputElement(i, 1, k, j, mA[i][k - 1][j - 1]);

        //        for (int i = 0; i < mA.size(); i++)
        //        {
        //            cout << "mA[" << i << "]: " << endl << mA[i] << endl;
        //        }
        //        cout << "box.is_empty(): " << box.is_empty() << endl;
        //        cout << "box[0].is_empty(): " << box[0].is_empty() << endl;
        //        cout << "box[1].is_empty(): " << box[1].is_empty() << endl;
        //        cout << "box.is_unbounded(): " << box.is_unbounded() << endl;
        //        cout << "box: " << box << endl;
    }else{
        box.set_empty();
    }
}

void CtcBoxLMI::contract(IntervalVector& box)
{
    for (int i = 1; i <= box.size(); i++)
    {
        // Minimize (lower bound)
        mkSolver(matrices, box);
        if (!check_empty_box(box))
        {
            solver.inputCVec(i, 1);
            solver.initializeUpperTriangle();
            solver.initializeSolve();
            solver.solve();
            box[i - 1] &= Interval(solver.getResultXVec()[i - 1], POS_INFINITY);
            if (!check_empty_box(box))
            {
                // Maximize (upper bound)
                mkSolver(matrices, box);

                solver.inputCVec(i, -1);
                solver.initializeUpperTriangle();
                solver.initializeSolve();
                solver.solve();
                box[i - 1] &= Interval(NEG_INFINITY, solver.getResultXVec()[i - 1]);
            }
            else
            {
                box.set_empty();
            }
        }
        else
        {
            box.set_empty();
        }
    }
}
