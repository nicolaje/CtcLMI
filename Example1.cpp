#include <cstdlib>
#include <ibex/ibex.h>
#include "CtcBoxLMI.h"

using namespace std;
using namespace ibex;

void sivia(const IntervalVector &search, BoolInterval(*inclusionTest)(IntervalVector), const double &EPSILON, string &vibesFigName);

int main()
{
    cout << "Main" << endl;
    MatrixArray mA(3, 2, 2);

    Matrix m0 = Matrix::zeros(2);

    Matrix m1(2, 2);
    m1 = Matrix::zeros(2);
    m1[0][0] = 1;
    m1[0][1] = 1;
    m1[1][0] = 1;

    Matrix m2(2, 2);
    m2 = Matrix::zeros(2);
    m2[0][1] = 1;
    m2[1][0] = 1;
    m2[1][1] = 1;

    mA[0] = m0;
    mA[1] = m1;
    mA[2] = m2;

    CtcBoxLMI cBox(mA);
    cout << "cBox" << endl;
    IntervalVector iV(2, Interval(-1000, 1000));

    cBox.contract(iV);

    cout << "iV: " << iV << endl;

    NumConstraint c1("x1", "x2", "(x1+x2+sqrt((x1-x2)^2+4*(x1+x2)^2))/2>=0");
    NumConstraint c2("x1", "x2", "(x1+x2-sqrt((x1-x2)^2+4*(x1+x2)^2))/2>=0");

    CtcFwdBwd ct1(c1);
    CtcFwdBwd ct2(c2);

    CtcCompo ctComp(ct1, ct2);
    
    iV = IntervalVector(2, Interval(8, 10));
    ctComp.contract(iV);
    
    cout << "ctComp(box): " << iV << endl;

    return EXIT_SUCCESS;
}

void sivia(const IntervalVector &search, BoolInterval(*inclusionTest)(IntervalVector), const double &EPSILON, string &vibesFigName)
{
    //    list
}